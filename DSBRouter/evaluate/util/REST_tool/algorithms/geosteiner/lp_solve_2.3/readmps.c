#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "patchlevel.h"
#include "lpkit.h"

/* defines */
#define VARNAME    0
#define CONSTRNAME 1

#define NAME    -1
#define ROWS    0
#define COLUMNS 1
#define RHS     2
#define BOUNDS  3
#define RANGES  4

/* global vars */
short  Column_ready;
short  Int_section;
lprec *Mlp;
REAL  *Last_column;
char   Last_col_name[25];
int    Lineno;
short  Debug;
short  Unconstrained_rows_found;

/*
A:  MPS format was named after an early IBM LP product and has emerged
as a de facto standard ASCII medium among most of the commercial LP
codes.  Essentially all commercial LP codes accept this format, but if
you are using public domain software and have MPS files, you may need
to write your own reader routine for this.  It's not too hard.  The
main things to know about MPS format are that it is column oriented (as
opposed to entering the model as equations), and everything (variables,
rows, etc.) gets a name.  MPS format is described in more detail in
Murtagh's book, referenced in another section. Also,

ftp://softlib.cs.rice.edu/pub/miplib/mps_format

is a nice short introduction.

MPS is an old format, so it is set up as though you were using punch
cards, and is not free format. Fields start in column 1, 5, 15, 25, 40
and 50.  Sections of an MPS file are marked by so-called header cards,
which are distinguished by their starting in column 1.  Although it is
typical to use upper-case throughout the file (like I said, MPS has
long historical roots), many MPS-readers will accept mixed-case for
anything except the header cards, and some allow mixed-case anywhere.
The names that you choose for the individual entities (constraints or
variables) are not important to the solver; you should pick names that
are meaningful to you, or will be easy for a post-processing code to
read.

Here is a little sample model written in MPS format (explained in more
detail below):

NAME          TESTPROB
ROWS
 N  COST
 L  LIM1
 G  LIM2
 E  MYEQN
COLUMNS
    XONE      COST                 1   LIM1                 1
    XONE      LIM2                 1
    YTWO      COST                 4   LIM1                 1
    YTWO      MYEQN               -1
    ZTHREE    COST                 9   LIM2                 1
    ZTHREE    MYEQN                1
RHS
    RHS1      LIM1                 5   LIM2                10
    RHS1      MYEQN                7
BOUNDS
 UP BND1      XONE                 4
 LO BND1      YTWO                -1
 UP BND1      YTWO                 1
ENDATA

means:

Optimize
 COST:    XONE + 4 YTWO + 9 ZTHREE
Subject To
 LIM1:    XONE + YTWO <= 5
 LIM2:    XONE + ZTHREE >= 10
 MYEQN:   - YTWO + ZTHREE  = 7
Bounds
 0 <= XONE <= 4
-1 <= YTWO <= 1
End

*/

/* copy a MPS name, only trailing spaces are removed. In MPS names can have
   embedded spaces! */
static void namecpy(char *into, char *from)
{
  int i;

  /* copy at most 8 characters of from, stop at end of string or newline */
  for(i = 0; (from[i] != '\0') &&  (from[i] != '\n') && (i < 8); i++)
    into[i] = from[i];

  /* end with end of string */
  into[i] = '\0';

  /* remove trailing spaces, if any */
  for(i--; (i >= 0) && (into[i] == ' '); i--)
    into[i] = '\0';
}

/* scan an MPS line, and pick up the information in the fields that are
   present */
static int scan_line(char* line, char *field1, char *field2, char *field3,
		     double *field4, char *field5, double *field6)
{
  int items = 0, line_len;
  char buf[16];

  line_len = strlen(line);
  
  if(line_len >= 1) { /* spaces or N/L/G/E or UP/LO */
    strncpy(buf, line, 4);
    buf[4] = '\0';
    sscanf(buf, "%s", field1);
    items++;
  }
  else
    field1[0] = '\0';

  line += 4;

  if(line_len >= 5) { /* name */
    namecpy(field2, line);
    items++;
  }
  else
    field2[0] = '\0';

  line += 10;

  if(line_len >= 14) { /* name */
    namecpy(field3, line);
    items++;
  }
  else
    field3[0] = '\0';

  line += 10;

  if(line_len >= 25) { /* number */
    strncpy(buf, line, 15);
    buf[15] = '\0';
    *field4 = atof(buf);
    items++;
  }
  else
    *field4 = 0;

  line += 15;

  if(line_len >= 40) { /* name */
    namecpy(field5, line);
    items++;
  }
  else
    field5[0] = '\0';
  line += 10;

  if(line_len >= 50) { /* number */
    strncpy(buf, line, 15);
    buf[15] = '\0';
    *field6 = atof(buf);
    items++;
  }
  else
    *field6 = 0;
    
  return(items);
}

static void addmpscolumn(void)
{
  int i;
  if (Column_ready) {
    add_column(Mlp, Last_column);
    set_col_name(Mlp, Mlp->columns, Last_col_name);
    set_int(Mlp, Mlp->columns, Int_section);
  }
  Column_ready = FALSE;
  for(i = 0; i <= Mlp->rows; i++)
    Last_column[i] = 0;   
}   

static int find_row(char *name)
{
  hashelem *hp;

  hp = findhash(name, Mlp->rowname_hashtab);

  if (hp == NULL) {
      if(Unconstrained_rows_found) { /* just ignore them in this case */
	return(-1);
      }
      else {
	fprintf(stderr, "Unknown row name (%s) on line %d\n", name, Lineno);
	exit(EXIT_FAILURE);
      }
  }
  return(hp->index);
}

static int find_var(char *name)
{
  hashelem *hp;

  hp = findhash(name, Mlp->colname_hashtab);

  if (hp == NULL) {
    fprintf(stderr, "Unknown variable name (%s) on line %d\n", name, Lineno);
    exit(EXIT_FAILURE);
  }
  return(hp->index);
}

lprec *read_mps(FILE *input, short verbose)
{
  char field1[5], field2[10], field3[10], field5[10], line[BUFSIZ], tmp[15];
  double field4, field6; 
  int section;
  int items;
  int row;
  int var;
  int i;
  nstring probname; 
  int OF_found = FALSE;

  Debug = verbose;
  Mlp = make_lp(0, 0);
  strcpy(Last_col_name, "");
  Int_section = FALSE;
  Column_ready = FALSE;
  Lineno = 0;

  /* let's initialize line to all zero's */
  memset(line, '\0', BUFSIZ);

  while(fgets(line, BUFSIZ - 1, input)) {
    Lineno++;

    /* skip lines which start with "*", they are comment */
    if(line[0] == '*') {
      if(Debug)
	fprintf(stderr, "Comment on line %d: %s", Lineno, line);
      continue;
    }

    if(Debug)
      fprintf(stderr, "Line %d: %s", Lineno, line);

    /* first check for "special" lines: NAME, ROWS, BOUNDS .... */
    /* this must start in the first position of line */
    if(line[0] != ' ') {
      sscanf(line, "%s", tmp);
      if(strcmp(tmp, "NAME") == 0) {
	section = NAME;
	sscanf(line, "NAME %s", probname);
	strcpy(Mlp->lp_name, probname);
      }
      else if(strcmp(tmp, "ROWS") == 0) {
	section = ROWS;
	if(Debug)
	  fprintf(stderr, "Switching to ROWS section\n");
      }
      else if(strcmp(tmp, "COLUMNS") == 0) {
	CALLOC(Last_column, Mlp->rows + 1);
	section = COLUMNS;
	if(Debug)
	  fprintf(stderr, "Switching to COLUMNS section\n");
      }
      else if(strcmp(tmp, "RHS") == 0) {
	addmpscolumn();
	section = RHS;
	if(Debug)
	  fprintf(stderr, "Switching to RHS section\n");
      }
      else if(strcmp(tmp, "BOUNDS") == 0) {
	section = BOUNDS;
	if(Debug)
	  fprintf(stderr, "Switching to BOUNDS section\n");
      }
      else if(strcmp(tmp, "RANGES") == 0) {
	section = RANGES;
	if(Debug)
	  fprintf(stderr, "Switching to RANGES section\n");
      }
      else if(strcmp(tmp, "ENDATA") == 0) {
	if(Debug)
	  fprintf(stderr, "Finished reading MPS file\n");
      }
      else { /* line does not start with space and does not match above */
	fprintf(stderr, "Unrecognized line %d: %s\n", Lineno, line);
	exit(EXIT_FAILURE);
      }
    }
    else { /* normal line, process */
      items = scan_line(line, field1, field2, field3, &field4, field5,
			&field6);

      switch(section) {

      case NAME:
	fprintf(stderr, "Error, extra line under NAME line\n");
	exit(EXIT_FAILURE);
	break;

      case ROWS:
	/* field1: rel. operator; field2: name of constraint */

	if(Debug) {
	  fprintf(stderr, "Rows line: ");
	  fprintf(stderr, "%s %s\n", field1, field2);
	}
         
	if(strcmp(field1, "N") == 0) {
	  if(!OF_found) { /* take the first N row as OF, ignore others */
	    set_row_name(Mlp, 0, field2); 
	    OF_found = TRUE;
	  }
	  else if(!Unconstrained_rows_found) {
	    fprintf(stderr, "Unconstrained row %s will be ignored\n", field2);
	    fprintf(stderr,
		    "Further messages of this kind will be surpressed\n");
	    Unconstrained_rows_found = TRUE;
	  }
	}
	else if(strcmp(field1, "L") == 0) {
	  str_add_constraint(Mlp, "" ,REL_LE ,0);
	  set_row_name(Mlp, Mlp->rows, field2); 
	}
	else if(strcmp(field1, "G") == 0) {
	  str_add_constraint(Mlp, "" ,REL_GE ,0);
	  set_row_name(Mlp, Mlp->rows, field2); 
	}
	else if(strcmp(field1, "E") == 0) {
	  str_add_constraint(Mlp, "",REL_EQ ,0);
	  set_row_name(Mlp, Mlp->rows, field2); 
	}
	else {
	  fprintf(stderr, "Unknown relat '%s' on line %d\n", field1, Lineno);
	  exit(EXIT_FAILURE);
	}
	break;

      case COLUMNS:
	/* field2: variable; field3: constraint; field4: coef */
	/* optional: field5: constraint; field6: coef */

	if(Debug) {
	  fprintf(stderr, "Columns line: ");
	  fprintf(stderr, "%s %s %g %s %g\n", field2, field3, field4, field5,
		  field6);
	}

	if((items == 4) || (items == 6)) {
	  if (strcmp(field2, Last_col_name) != 0 && Column_ready) { 
	    addmpscolumn();
	    strcpy(Last_col_name, field2);
	    Column_ready = TRUE;
	  }
	  else {
	    strcpy(Last_col_name, field2);
	    Column_ready = TRUE;
	  }
	  if((row = find_row(field3)) >= 0) {
	    Last_column[row] = (REAL)field4;             
	  }
	}
	if(items == 6) {
	  if((row = find_row(field5)) >= 0) {
	    Last_column[row] = (REAL)field6;
	  }
	}

	if(items == 5) { /* there might be an INTEND or INTORG marker */
	  /* look for "    <name>  'MARKER'                 'INTORG'" */
	  /* or "    <name>  'MARKER'                 'INTEND'" */
	  if(strcmp(field3, "'MARKER'") ==0) {
	    addmpscolumn();
	    if(strcmp(field5, "'INTORG'") == 0) {
	      Int_section = TRUE;
	      if(Debug)
		fprintf(stderr, "Switching to integer section\n");
	    }
	    else if(strcmp(field5, "'INTEND'") == 0) {
	      Int_section = FALSE;
	      if(Debug)
		fprintf(stderr, "Switching to non-integer section\n");
	    }
	    else
	      fprintf(stderr, "Unknown marker (ignored) at line %d: %s\n",
		      Lineno, field5);
	  }
	}

	if((items != 4) && (items != 6) && (items != 5)) { /* Wrong! */
	  fprintf(stderr,
		  "Wrong number of items (%d) in COLUMNS section (line %d)\n",
		  items, Lineno);
	  exit(EXIT_FAILURE);
	}
	break;

      case RHS:
	/* field2: uninteresting name; field3: constraint name */
	/* field4: value */
	/* optional: field5: constraint name; field6: value */

	if(Debug) {
	  fprintf(stderr, "RHS line: ");
	  fprintf(stderr, "%s %s %g %s %g\n", field2, field3, field4, field5,
		  field6);
	}

	if((items != 4) && (items != 6)) {
	  fprintf(stderr,
		  "Wrong number of items (%d) in RHS section line %d\n",
		  items, Lineno);
	  exit(EXIT_FAILURE);
	}

	if((row = find_row(field3)) >= 0) {
	  set_rh(Mlp, row, (REAL)field4);
	}

	if(items == 6) {
	  if((row = find_row(field5)) >= 0) {
	    set_rh(Mlp, row, (REAL)field6);
	  }
	}

	break;

      case BOUNDS:
	/* field1: bound type; field2: uninteresting name; */
	/* field3: variable name; field4: value */

	if(Debug) {
	  fprintf(stderr, "BOUNDS line: %s %s %s %g\n", field1, field2, field3,
		  field4);
	}

	var = find_var(field3);
	
	if(strcmp(field1, "UP") == 0) {
	  /* upper bound */
	  set_upbo(Mlp, var, field4);
	}
	else if(strcmp(field1, "LO") == 0) {
	  /* lower bound */
	  set_lowbo(Mlp, var, field4);
	}
	else if(strcmp(field1, "FX") == 0) {
	  /* fixed, upper _and_ lower  */
	  set_lowbo(Mlp, var, field4);
	  set_upbo(Mlp, var, field4);
	}
	else if(strcmp(field1, "PL") == 0)
	  /* normal, 0 <= var <= inf, do nothing */;
	else if(strcmp(field1, "BV") == 0) { /* binary variable */
	  set_upbo(Mlp, var, 1);
	  set_int(Mlp, var, TRUE);
	}
	/* hack for free and negative variables. Ugly, and does not
	   always work. MB */
	else if(strcmp(field1, "FR") == 0) { /* free variable */
	  fprintf(stderr, 
		  "Free variable %s is split in a positive part %s and a negative part %s_\n",
		  field3, field3, field3);
	  get_column(Mlp, var, Last_column);
	  for (i = 0; i <= Mlp->rows; i++)
	    Last_column[i]*=-1;
	  add_column(Mlp, Last_column);
	  strcat(field3, "_");
	  set_col_name(Mlp, Mlp->columns, field3);
	  /* should lower and upper bounds of both variables be adjusted? What
	     if lower and upper bound are specified later? MB */
	}
	else if(strcmp(field1, "MI") == 0) { /* negative variable */
	  fprintf(stderr,
		  "Negative variable %s will be represented by - %s-\n",
		  field3, field3);
	  get_column(Mlp, var, Last_column);
	  del_column(Mlp, var);
	  for (i = 0; i <= Mlp->rows; i++)
	    Last_column[i] *= -1;
	  add_column(Mlp, Last_column);
	  strcat(field3, "-");
	  set_row_name(Mlp, var, field3);
	  /* should lower and upper bounds of variable be adjusted? What if
	     lower and upper bound are specified later? (does not work!) MB */
	}
	else {
	  fprintf(stderr, "BOUND type %s on line %d is not supported\n",
		  field1, Lineno);
	  exit(EXIT_FAILURE);
	}
	      
	break;
               
      case RANGES:

	/* We have to implement the following semantics:

	  D. The RANGES section is for constraints of the form: h <=
	  constraint <= u .  The range of the constraint is r = u - h .  The
	  value of r is specified in the RANGES section, and the value of u or
	  h is specified in the RHS section.  If b is the value entered in the
	  RHS section, and r is the value entered in the RANGES section, then
	  u and h are thus defined:

          row type       sign of r       h          u
          ----------------------------------------------
             G            + or -         b        b + |r|
             L            + or -       b - |r|      b
             E              +            b        b + |r|
             E              -          b - |r|      b
	*/

	/* field2: uninteresting name; field3: constraint name */
	/* field4: value */
	/* optional: field5: constraint name; field6: value */

	if(Debug) {
	  fprintf(stderr, "RANGES line: ");
	  fprintf(stderr, "%s %s %g %s %g\n", field2, field3, field4, field5,
		  field6);
	}

	if((items != 4) && (items != 6)) {
	  fprintf(stderr,
		  "Wrong number of items (%d) in RANGES section line %d\n",
		  items, Lineno);
	  exit(EXIT_FAILURE);
	}

	if(((row = find_row(field3)) >= 0) && (field4 != 0)) {
	  /* find out constraint type. If ch_sign[row] is TRUE, it is REL_GE. If
	     ch_sign[row] is FALSE, it is an equality constraint if
	     orig_upbo[row] == 0. For a REL_LE constraint, orig_upbo[row] should be
	     +infinity */
	 
	  if(my_abs(field4) >= Mlp->infinite) {
	    fprintf(stderr,
		    "Warning, Range for row %s >= infinity (value %g) on line %d, ignoring\n",
		    field3, field4, Lineno);
	  }
	  else if(Mlp->ch_sign[row]) {
	    /* REL_GE */
	    Mlp->orig_upbo[row] = my_abs(field4);
	  }
	  else if(Mlp->orig_upbo[row] == 0 && field4 >= 0) {
	    /*  REL_EQ with positive sign of r value */
	    set_constr_type(Mlp, row, REL_GE);
	    Mlp->orig_upbo[row] = field4;
	  }
	  else if(Mlp->orig_upbo[row] == Mlp->infinite) {
	    /* REL_LE */
	    Mlp->orig_upbo[row] = my_abs(field4);
	  }
	  else if(Mlp->orig_upbo[row] == 0 && field4 < 0) {
	    /* REL_EQ with negative sign of r value */
	    set_constr_type(Mlp, row, REL_LE);
	    Mlp->orig_upbo[row] = -field4;
	  }
	  else { /* let's be paranoid */
	    fprintf(stderr,
		    "Cannot figure out row type, row = %d, ch_sign = %d, upbo = %g\n",
		    row, Mlp->ch_sign[row], (double)Mlp->orig_upbo[row]);
	  }
	}

	if(items == 6) {
	  if(((row = find_row(field5)) >= 0) && (field6 != 0)) {
	    /* find out constraint type. If ch_sign[row] is TRUE, it is REL_GE. If
	       ch_sign[row] is FALSE, it is an equality constraint if
	       orig_upbo[row] == 0. For a REL_LE constraint, orig_upbo[row] should
	       be +infinity */
	    
	    if(my_abs(field6) >= Mlp->infinite) {
	      fprintf(stderr,
		      "Warning, Range for row %s >= infinity (value %g) on line %d, ignoring\n",
		      field5, field6, Lineno);
	    }
	    else if(Mlp->ch_sign[row]) {
	      /* REL_GE */
	      Mlp->orig_upbo[row] = my_abs(field6);
	    }
	    else if(Mlp->orig_upbo[row] == 0 && field6 >= 0) {
	      /*  REL_EQ with positive sign of r value */
	      set_constr_type(Mlp, row, REL_GE);
	      Mlp->orig_upbo[row] = field6;
	    }
	    else if(Mlp->orig_upbo[row] == Mlp->infinite) {
	      /* REL_LE */
	      Mlp->orig_upbo[row] = my_abs(field6);
	    }
	    else if(Mlp->orig_upbo[row] == 0 && field6 < 0) {
	      /* REL_EQ with negative sign of r value */
	      set_constr_type(Mlp, row, REL_LE);
	      Mlp->orig_upbo[row] = -field6;
	    }
	    else { /* let's be paranoid */
	      fprintf(stderr,
		      "Cannot figure out row type, row = %d, ch_sign = %d, upbo = %g\n",
		      row, Mlp->ch_sign[row], (double)Mlp->orig_upbo[row]);
	    }
	  }
	  
	}
	break;
      }
    }
  }
  return(Mlp);
}
