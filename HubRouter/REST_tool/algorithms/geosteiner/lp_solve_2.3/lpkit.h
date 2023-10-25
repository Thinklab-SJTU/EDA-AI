/*
  Main header file of the LP_SOLVE toolkit.
  
  Original by Jeroen Dirks, 21-2-95
  Maintained by Michel Berkelaar

  include this file in your program and link with liblps.a

04/24/2014: warme: fix 64-bit architecture issues.
*/

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "hash.h"

#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1	/* older <stdlib.h> seem to lack this... */
#endif

#define FALSE   0
#define TRUE    1

#define DEFNUMINV 50
#define INITIAL_MAT_SIZE 10000

/* solve status values */
#define OPTIMAL     	0
#define MILP_FAIL   	1
#define INFEASIBLE  	2
#define UNBOUNDED   	3
#define FAILURE     	4
#define RUNNING     	5

/* test_branch extra status values */
#define STOP_AT_INVERT	6

/* lag_solve extra status values */
#define FEAS_FOUND   	6
#define NO_FEAS_FOUND 	7
#define BREAK_BB	8

#define FIRST_NI	0
#define RAND_NI		1

#define REL_LE      0
#define REL_EQ      1
#define REL_GE      2
#define REL_OF      3

#if 1
/* DMW - gcc does fabs inline without conditional jumps! */
#define my_abs(x)	(fabs (x))
#else
#define	my_abs(x)	((x) < 0 ? -(x) : (x))
#endif
#define my_min(x, y)    ((x) < (y) ? (x) : (y))
#define my_max(x, y)    ((x) > (y) ? (x) : (y))

#define MAX_WARN_COUNT 20

#ifdef CHECK
#define my_round(val, eps) { \
	REAL absv; \
        absv = ((val) < 0 ? -(val) : (val)); \
        if(absv < (eps)) \
          val = 0; \
	if(Warn_count < MAX_WARN_COUNT) \
	  { \
	    if(absv > 0.5 * (eps) && absv < 2 * (eps)) \
	      { \
		Warn_count++; \
		printf("%% Warning Value close to epsilon V: %e E: %e\n", \
		       (double)absv, (double)(eps)); \
		if(Warn_count == MAX_WARN_COUNT) \
		  { \
		    printf("%% *** Surpressing further rounding warnings\n"); \
		  } \
	      } \
	  } \
}

#else
#if 1
/* DMW - gcc does fabs inline without conditional jumps! */
#define my_round(val,eps) {if (fabs(val) < (eps)) val = 0;}
#else
#define my_round(val,eps) if (((val) < 0 ? -(val) : (val)) < (eps)) val = 0;
#endif
#endif


#define DEF_INFINITE  1e24 /* limit for dynamic range */
#define DEF_EPSB      5.01e-7 /* for rounding RHS values to 0 determine	
				 infeasibility basis */
#define DEF_EPSEL     1e-8 /* for rounding other values (vectors) to 0 */
#define DEF_EPSD      1e-6 /* for rounding reduced costs to zero */
#define DEF_EPSILON   1e-3 /* to determine if a float value is integer */
 
#define PREJ          1e-3  /* pivot reject (try others first) */

#ifndef REAL /* to allow -DREAL=<float type> while compiling */
#define REAL double
#endif

#define ETA_START_SIZE 10000 /* start size of array Eta. Realloced if needed */
#define FNAMLEN 64
#define NAMELEN 25
#define MAXSTRL (NAMELEN-1)
#define STD_ROW_NAME_PREFIX "r_"

#define CALLOC(ptr, nr)\
  if(!(ptr = calloc((size_t)(nr), sizeof(*ptr))) && nr) {\
    fprintf(stderr, "calloc of %ld bytes failed on line %d of file %s\n",\
            (long) (nr * sizeof(*ptr)), __LINE__, __FILE__);\
    exit(EXIT_FAILURE);\
  }

#define MALLOC(ptr, nr)\
  if(!(ptr = malloc((size_t)((nr) * sizeof(*ptr)))) && nr) {\
    fprintf(stderr, "malloc of %ld bytes failed on line %d of file %s\n",\
            (long) (nr * sizeof(*ptr)), __LINE__, __FILE__);\
    exit(EXIT_FAILURE);\
  }

#define REALLOC(ptr, nr)\
  if(!(ptr = realloc(ptr, (size_t)((nr) * sizeof(*ptr)))) && nr) {\
    fprintf(stderr, "realloc of %ld bytes failed on line %d of file %s\n",\
            (long) (nr * sizeof(*ptr)), __LINE__, __FILE__);\
    exit(EXIT_FAILURE);\
  }

#define MALLOCCPY(nptr, optr, nr)\
  {MALLOC(nptr, nr); memcpy(nptr, optr, (size_t)((nr) * sizeof(*optr)));}

#define MEMCPY(nptr, optr, nr)\
  memcpy(nptr, optr, (size_t)((nr) * sizeof(*optr)));

typedef char nstring[NAMELEN];

typedef struct _matrec
{
  int row_nr;
  REAL value;
} matrec;

typedef struct _column
{
  int            row;
  REAL           value;
  struct _column *next ;
} column;

typedef struct _constraint_name
{
  char                    name[NAMELEN];
  int                     row;
  struct _constraint_name *next;
} constraint_name;

typedef struct _bound
{
  REAL          upbo;
  REAL          lowbo;
} bound;

typedef struct _tmp_store_struct
{
  nstring name;
  int     row;
  REAL    value;
  REAL    rhs_value;
  short   relat;
} tmp_store_struct;

typedef struct _rside /* contains relational operator and rhs value */
{
  REAL          value;
  struct _rside *next;
  short         relat;
} rside;


/* fields indicated with ## may be modified directly */
/* pointers will have their array size in the comments */

typedef struct _lprec
{
  nstring   lp_name;		/* the name of the lp */

  short     verbose;            /* ## Verbose flag */
  short     print_duals;        /* ## PrintDuals flag for PrintSolution */
  short     print_sol;          /* ## used in lp_solve */
  short     debug;              /* ## Print B&B information */
  short     print_at_invert;    /* ## Print information at every reinversion */
  short     trace;              /* ## Print information on pivot selection */
  short     anti_degen;		/* ## Do perturbations */
  short     do_presolve;        /* perform matrix presolving */

  int	    rows;               /* Nr of constraint rows in the problem */
  int       rows_alloc;      	/* The allocated memory for Rows sized data */
  int       columns;            /* The number of columns (= variables) */
  int       columns_alloc;  
  int       sum;                /* The size of the variables + the slacks */
  int       sum_alloc;

  short     names_used;         /* Flag to indicate if names for rows and
				   columns are used */
  nstring   *row_name;		/* rows_alloc+1 */
  nstring   *col_name;		/* columns_alloc+1 */

 /* Row[0] of the sparce matrix is the objective function */

  int       non_zeros;          /* The number of elements in the sparce matrix*/
  int       mat_alloc;		/* The allocated size for matrix sized 
				   structures */
  matrec    *mat;               /* mat_alloc :The sparse matrix */
  int       *col_end;           /* columns_alloc+1 :Cend[i] is the index of the
		 		   first element after column i.
				   column[i] is stored in elements 
				   col_end[i-1] to col_end[i]-1 */
  int       *col_no;            /* mat_alloc :From Row 1 on, col_no contains the
				   column nr. of the
                                   nonzero elements, row by row */
  short     row_end_valid;	/* true if row_end & col_no are valid */
  int       *row_end;           /* rows_alloc+1 :row_end[i] is the index of the 
				   first element in Colno after row i */
  REAL      *orig_rh;           /* rows_alloc+1 :The RHS after scaling & sign
				  changing, but before `Bound transformation' */
  REAL      *rh;		/* rows_alloc+1 :As orig_rh, but after Bound 
				   transformation */
  REAL      *rhs;		/* rows_alloc+1 :The RHS of the current
				   simplex tableau */
  short     *must_be_int;       /* sum_alloc+1 :TRUE if variable must be 
				   Integer */
  REAL      *orig_upbo;         /* sum_alloc+1 :Bound before transformations */
  REAL      *orig_lowbo;	/*  "       "                   */
  REAL      *upbo;              /*  " " :Upper bound after transformation &
				   B&B work */
  REAL      *lowbo;             /*  "       "  :Lower bound after transformation
				   & B&B work */

  short     basis_valid;        /* TRUE is the basis is still valid */
  int       *bas;               /* rows_alloc+1 :The basis column list */
  short     *basis;             /* sum_alloc+1 : basis[i] is TRUE if the column
				   is in the basis */
  short     *lower;             /*  "       "  :TRUE if the variable is at its 
				   lower bound (or in the basis), it is FALSE
				   if the variable is at its upper bound */
  REAL	    *colsign;		/*  "       "  : colsign[i] = lower[i]?1:-1 */

  short     eta_valid;          /* TRUE if current Eta structures are valid */
  int       eta_alloc;          /* The allocated memory for Eta */
  int       eta_size;           /* The number of Eta columns */
  int       num_inv;            /* The number of real pivots */
  int       max_num_inv;        /* ## The number of real pivots between 
				   reinversions */
  REAL      *eta_value;         /* eta_alloc :The Structure containing the
				   values of Eta */
  int       *eta_row_nr;         /*  "     "  :The Structure containing the Row
				   indexes of Eta */
  int       *eta_col_end;       /* rows_alloc + MaxNumInv : eta_col_end[i] is
				   the start index of the next Eta column */

  short	    bb_rule;		/* what rule for selecting B&B variables */

  short     break_at_int;       /* TRUE if stop at first integer better than
                                   break_value */
  REAL      break_value;        

  REAL      obj_bound;          /* ## Objective function bound for speedup of 
				   B&B */
  int       iter;               /* The number of iterations in the simplex
				   solver (LP) */
  int       total_iter;         /* The total number of iterations (B&B)
				   (ILP) */
  int       max_level;          /* The Deepest B&B level of the last solution */
  int	    total_nodes;	/* total number of nodes processed in b&b */
  REAL      *solution;          /* sum_alloc+1 :The Solution of the last LP, 
				   0 = The Optimal Value, 
                                   1..rows The Slacks, 
				   rows+1..sum The Variables */
  REAL      *best_solution;     /*  "       "  :The Best 'Integer' Solution */
  REAL      *duals;             /* rows_alloc+1 :The dual variables of the
				   last LP */
  
  short     maximise;           /* TRUE if the goal is to maximise the 
				   objective function */
  short     floor_first;        /* TRUE if B&B does floor bound first */
  short     *ch_sign;           /* rows_alloc+1 :TRUE if the Row in the matrix
				   has changed sign 
                                   (a`x > b, x>=0) is translated to 
				   s + -a`x = -b with x>=0, s>=0) */ 

  short     scaling_used;	/* TRUE if scaling is used */
  short     columns_scaled;     /* TRUE is the columns are scaled too, Only use
		 		   if all variables are non-integer */
  REAL      *scale;             /* sum_alloc+1:0..Rows the scaling of the Rows,
				   Rows+1..Sum the scaling of the columns */

  int	    nr_lagrange;	/* Nr. of Langrangian relaxation constraints */
  REAL	    **lag_row;		/* NumLagrange, columns+1:Pointer to pointer of 
				   rows */
  REAL      *lag_rhs;		/* NumLagrange :Pointer to pointer of Rhs */
  REAL      *lambda;		/* NumLagrange :Lambda Values */
  short     *lag_con_type;      /* NumLagrange :TRUE if constraint type REL_EQ */
  REAL      lag_bound;		/* the lagrangian lower bound */

  short     valid;		/* Has this lp pased the 'test' */
  REAL      infinite;           /* ## numerical stuff */
  REAL      epsilon;            /* ## */
  REAL      epsb;               /* ## */
  REAL      epsd;               /* ## */
  REAL      epsel;              /* ## */
  hashtable *rowname_hashtab;   /* hash table to store row names */
  hashtable *colname_hashtab;   /* hash table to store column names */
} lprec;


/* Saved basis status for rapidly testing branch variables... */
struct basis_save {
	int		eta_size;
	int *		bas;
	short *		basis;
	short *		lower;
	REAL *		rhs;

	REAL *		drow;
	REAL *		prow;
	REAL *		Pcol;
};


/* function interface for the user */

lprec *make_lp(int rows, int columns);
/* create and initialise a lprec structure
   defaults:
   Empty (Rows * Columns) matrix,
   Minimise the objective function
   constraints all type <=
   Upperbounds all Infinite
   no integer variables
   floor first in B&B
   no scaling
   default basis */

lprec *read_lp_file(FILE *input, short verbose, nstring lp_name);
/* create and read an .lp file from input (input must be open) */

void delete_lp(lprec *lp);
/* Remove problem from memory */

lprec *copy_lp(lprec *lp);
/* copy a lp structure */

void set_mat(lprec *lp, int row, int column, REAL value);
/* fill in element (Row,Column) of the matrix
   Row in [0..Rows] and Column in [1..Columns] */

void set_obj_fn(lprec *lp, REAL *row);
/* set the objective function (Row 0) of the matrix */
void str_set_obj_fn(lprec *lp, char *row);
/* The same, but with string input */

void add_constraint(lprec *lp, REAL *row, short constr_type, REAL rh);
/* Add a constraint to the problem,
   row is the constraint row,
   rh is the right hand side,
   constr_type is the type of constraint (REL_LE (<=), REL_GE(>=), REL_EQ(=)) */
void str_add_constraint(lprec *lp, char *row_string ,short constr_type, REAL rh);
/* The same, but with string input */

void add_rows(lprec *lp,
	      int ccnt,
	      int rcnt,
	      REAL *rh,
	      short *ctype,
	      int *rmatbeg,
	      int *rmatind,
	      REAL *rmatval);
/* Efficiently add several new rows, possibly referencing several new
   columns.
	ccnt is number of new columns,
	rcnt is number of new rows,
	rh is vector of new right hand sides,
	ctype is vector of new constraint types,
	rmatbeg[i] (0 <= i <= rcnt) row i has coefficient rmatval[j]
		in column rmatind[j], for rmatbeg [i] <= j < rmatbeg[i+1].
		rmatbeg[0] must be zero.
	rmatind is vector of column numbers,
	rmatval is vector of non-zero coefficients. */

void del_constraint(lprec *lp,int del_row);
/* Remove constrain nr del_row from the problem */

void delete_row_set(lprec *lp, int *row_flags);
/* Efficiently delete a specified subset of the rows. */

void add_lag_con(lprec *lp, REAL *row, short con_type, REAL rhs);
/* add a Lagrangian constraint of form Row' x contype Rhs */
void str_add_lag_con(lprec *lp, char *row, short con_type, REAL rhs);
/* The same, but with string input */

void add_column(lprec *lp, REAL *column);
/* Add a Column to the problem */
void str_add_column(lprec *lp, char *col_string);
/* The same, but with string input */

void del_column(lprec *lp, int column);
/* Delete a column */

void set_upbo(lprec *lp, int column, REAL value);
/* Set the upperbound of a variable */

void set_lowbo(lprec *lp, int column, REAL value);
/* Set the lowerbound of a variable */

void set_bounds(lprec *lp, int column, REAL lower, REAL upper);
/* Set both the lowerbound and upperbound of a variable */

void set_int(lprec *lp, int column, short must_be_int);
/* Set the type of variable, if must_be_int = TRUE then the variable must be integer */

void set_rh(lprec *lp, int row, REAL value);
/* Set the right hand side of a constraint row */

void set_rh_vec(lprec *lp, REAL *rh);
/* Set the right hand side vector */
void str_set_rh_vec(lprec *lp, char *rh_string);
/* The same, but with string input */

void set_maxim(lprec *lp);
/* maximise the objective function */

void set_minim(lprec *lp);
/* minimise the objective function */

void set_constr_type(lprec *lp, int row, short con_type);
/* Set the type of constraint in row Row (REL_LE, REL_GE, REL_EQ) */

void set_row_name(lprec *lp, int row, nstring new_name);
/* Set the name of a constraint row, make sure that the name has < 25 characters */

void set_col_name(lprec *lp, int column, nstring new_name);
/* Set the name of a varaible column, make sure that the name has < 25 characters */

void auto_scale(lprec *lp);
/* Automatic scaling of the problem */

void unscale(lprec *lp);
/* Remove all scaling from the problem */

int solve(lprec *lp);
/* Solve the problem */

int lag_solve(lprec *lp, REAL start_bound, int num_iter, short verbose);
/* Do NumIter iterations with Lagrangian relaxation constraints */

void reset_basis(lprec *lp);
/* Reset the basis of a problem, can be usefull in case of degeneracy - JD */

REAL mat_elm(lprec *lp, int row, int column);
/* get a single element from the matrix */

void get_row(lprec *lp, int row_nr, REAL *row);
/* fill row with the row row_nr from the problem */

void get_column(lprec *lp, int col_nr, REAL *column);
/* fill column with the column col_nr from the problem */

void get_reduced_costs(lprec *lp, REAL *rc);
/* get the reduced costs vector */

void get_slack_vars(lprec *lp, REAL *slacks);
/* get the values of the slack variables */

void save_LP_basis (lprec *lp, struct basis_save *basp);
/* save off the basis for rapid testing of branch variables */

REAL try_branch (lprec *lp, int var, int dir, REAL *, REAL, struct basis_save *basp);
/* give a prospective branch-variable a "test run" to see how well is works */

void destroy_LP_basis (struct basis_save *basp);
/* free up the memory for a saved basis */

short is_feasible(lprec *lp, REAL *values);
/* returns TRUE if the vector in values is a feasible solution to the lp */

short column_in_lp(lprec *lp, REAL *column);
/* returns TRUE if column is already present in lp. (Does not look at bounds
   and types, only looks at matrix values */

lprec *read_mps(FILE *input, short verbose);
/* read a MPS file */

void write_MPS(lprec *lp, FILE *output);
/* write a MPS file to output */

void write_LP(lprec *lp, FILE *output);
/* write a LP file to output */

void print_lp(lprec *lp);
/* Print the current problem, only usefull in very small (test) problems. 
  Shows the effect of scaling */

void print_solution(lprec *lp);
/* Print the solution to stdout */

void print_duals(lprec *lp);
/* Print the dual variables of the solution */

void print_scales(lprec *lp);
/* If scaling is used, print the scaling factors */

void dump_lp (lprec *lp, char *filename);
/* Dump LP problem (in binary form) to given filename */

lprec * load_lp (FILE * fp);
/* Load LP problem (in binary form) from given input stream */




/* functions used internaly by the lp toolkit */
void error(char *format, ...);
void inc_mat_space(lprec *lp, int max_extra);
void inc_row_space(lprec *lp);
void inc_col_space(lprec *lp);
void unscale_columns(lprec *lp);
void btran(lprec *lp, REAL *row);
short invert(lprec *lp);
void presolve(lprec *lp);


/* define yyparse() to make compilers happy. There should be some system
   include file for this */
int yyparse(void);
