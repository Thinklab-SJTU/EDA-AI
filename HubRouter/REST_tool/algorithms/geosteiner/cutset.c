/***********************************************************************

	$Id: cutset.c,v 1.14 2016/09/24 17:53:55 warme Exp $

	File:	cutset.c
	Rev:	e-3
	Date:	09/24/2016

	Copyright (c) 1996, 2016 by David M. Warme.  This work is
	licensed under a Creative Commons Attribution 4.0 International
	License.

************************************************************************

	Separation routines for cutset constraints.

************************************************************************

	Modification Log:

	b-1:	11/14/96	warme
		: Split off from bb.c.
	b-2:	02/28/2001	warme
		: Changes for 3.1 release.
		: Split add_cutset_to_list, and its subroutines
		:  off to new cutsubs.c file.
	c-1:	08/05/2002	benny
		: Numerous changes for library release.
		: Uses parameters.
		: Uses channels for trace output.
	e-1:	04/14/2015	warme
		: Changes for 5.0 release.
	e-2:	09/05/2016	warme
		: Change notices for 5.1 release.
	e-3:	09/24/2016	warme
		: Reorganize include files, apply prefixes.
		: Upgrade fatals.

************************************************************************/

#include "cutset.h"

#include "bb.h"
#include "constrnt.h"
#include "fatal.h"
#include "flow.h"
#include "geosteiner.h"
#include "logic.h"
#include "memory.h"
#include "parmblk.h"
#include "sec_heur.h"
#include "steiner.h"
#include <string.h>


/*
 * Global Routines
 */

void		_gst_build_cutset_separation_formulation (
						bitmap_t *	vert_mask,
						bitmap_t *	edge_mask,
						struct bbinfo *	bbip);
struct constraint * _gst_find_cutset_constraints (double *	x,
					     bitmap_t *		vert_mask,
					     bitmap_t *		edge_mask,
					     struct bbinfo *	bbip);
struct constraint * _gst_find_fractional_cutsets (double *	x,
					     bitmap_t *		vert_mask,
					     bitmap_t *		edge_mask,
					     struct bbinfo *	bbip);
void		_gst_free_cutset_separation_formulation (
						struct cs_info * csip);


/*
 * Local Routines
 */

static struct constraint * enumerate_cuts (int			i,
					   int			ntaken,
					   bool *		cstack,
					   struct constraint *	cutlist,
					   bitmap_t *		comps,
					   int			ncomps,
					   bitmap_t *		cut_terms,
					   double *		x,
					   bitmap_t *		vert_mask,
					   bitmap_t *		edge_mask,
					   struct bbinfo *	bbip);
static int		find_comps (bitmap_t *			sol,
				    bitmap_t *			comps,
				    struct gst_hypergraph *	cip);
static struct constraint * simple_cuts (struct constraint *	cutlist,
					bitmap_t *		comps,
					int			ncomps,
					double *		x,
					bitmap_t *		vert_mask,
					bitmap_t *		edge_mask,
					struct bbinfo *		bbip);

/*
 * This routine quickly finds cutsets of zero weight -- totally
 * disconnected solutions.  It first computes the connected components
 * of the solution and then uses either a combinatorially thorough
 * method to generate a complete set of constraints, or a quicker
 * method to generate one constraint per component.  The method used
 * depends upon the number of connected components found.
 */

	struct constraint *
_gst_find_cutset_constraints (

double *		x,		/* IN - LP solution */
bitmap_t *		vert_mask,	/* IN - set of valid vertices */
bitmap_t *		edge_mask,	/* IN - set of valid hyperedges */
struct bbinfo *		bbip		/* IN - branch and bound info */
)
{
int			i;
int			nverts;
int			nedges;
int			kmasks;
int			nmasks;
int			ncomps;
struct gst_hypergraph *	cip;
bitmap_t *		sol;
bitmap_t *		cut_terms;
struct constraint *	cutlist;
bitmap_t *		comps;
bool *			cstack;

	cip	= bbip -> cip;
	nverts	= cip -> num_verts;
	nedges	= cip -> num_edges;
	kmasks	= cip -> num_vert_masks;
	nmasks	= cip -> num_edge_masks;

	comps	= NEWA (nverts * kmasks, bitmap_t);
	sol	= NEWA (nmasks, bitmap_t);

	cutlist = NULL;

	/* Get bit-mask of all full-sets even partially present in soln. */
	memset (sol, 0, nmasks * sizeof (*sol));
	for (i = 0; i < nedges; i++) {
		if (NOT BITON (edge_mask, i)) continue;
		if (x [i] > FUZZ) {
			SETBIT (sol, i);
		}
	}

	/* Determine the number of connected components it has.	*/
	ncomps = find_comps (sol, comps, cip);
	FATAL_ERROR_IF (ncomps <= 0);
	if (ncomps EQ 1) {
		free ((char *) sol);
		free ((char *) comps);
		return (NULL);
	}

#if 1
	gst_channel_printf (bbip -> params -> print_solve_trace,
			    "@cutset: %d connected components.\n", ncomps);
#endif

	if (ncomps > bbip -> params -> max_cutset_enumerate_comps) {
		/* Too many components to bother trying all combinations */
		/* of them -- just cut around each component...		 */
		cutlist = simple_cuts (cutlist,
				       comps,
				       ncomps,
				       x,
				       vert_mask,
				       edge_mask,
				       bbip);
	}
	else {
		/* Enumerate all cuts induced by the connected	*/
		/* components.  For each cut, we may have to	*/
		/* generate a constraint (if it has not been	*/
		/* generated before).				*/
		cut_terms  = NEWA (kmasks, bitmap_t);
		cstack = NEWA (nverts, bool);
		cstack [0] = TRUE;	/* Always take 1st component... */
		cutlist = enumerate_cuts (1,
					  1,
					  cstack,
					  cutlist,
					  comps,
					  ncomps,
					  cut_terms,
					  x,
					  vert_mask,
					  edge_mask,
					  bbip);
		free ((char *) cstack);
		free ((char *) cut_terms);
	}

	free ((char *) sol);
	free ((char *) comps);

	return (cutlist);
}

/*
 * This routine finds the connected components of the given set of
 * full-sets.  The result is a partition of the given terminals that
 * we represent as an array of terminal bit-masks, one per component.
 */

	static
	int
find_comps (

bitmap_t *		sol,	/* IN - solution to get components of. */
bitmap_t *		comps,	/* OUT - components. */
struct gst_hypergraph *	cip	/* IN - compatibility info. */
)
{
int			i;
int			e;
int			v;
int			nverts;
int			nedges;
int			kmasks;
int			nmasks;
int			ncomps;
int *			ep1;
int *			ep2;
int *			vp1;
int *			vp2;
bitmap_t *		visited;
int *			sp;
int *			stack;

	nverts = cip -> num_verts;
	nedges = cip -> num_edges;
	kmasks = cip -> num_vert_masks;
	nmasks = cip -> num_edge_masks;

	visited = NEWA (kmasks, bitmap_t);
	stack	= NEWA (nverts, int);

	for (i = 0; i < kmasks; i++) {
		visited [i] = 0;
	}

	ncomps = 0;
	for (;;) {
		/* Find next remaining full-set... */
		for (e = 0; e < nedges; e++) {
			if (BITON (sol, e)) break;
		}
		if (e >= nedges) break;
		CLRBIT (sol, e);

		/* Identify next component... */
		for (i = 0; i < kmasks; i++) {
			comps [i] = 0;
		}
		sp = &stack [0];
		vp1 = cip -> edge [e];
		vp2 = cip -> edge [e + 1];
		while (vp1 < vp2) {
			v = *vp1++;
			if (BITON (visited, v)) continue;
			SETBIT (visited, v);
			SETBIT (comps, v);
			*sp++ = v;
		}

		while (sp > stack) {
			v = *--sp;
			ep1 = cip -> term_trees [v];
			ep2 = cip -> term_trees [v + 1];
			while (ep1 < ep2) {
				e = *ep1++;
				if (NOT BITON (sol, e)) continue;
				CLRBIT (sol, e);
				vp1 = cip -> edge [e];
				vp2 = cip -> edge [e + 1];
				while (vp1 < vp2) {
					v = *vp1++;
					if (BITON (visited, v)) continue;
					SETBIT (visited, v);
					SETBIT (comps, v);
					*sp++ = v;
				}
			}
		}
		++ncomps;
		comps += kmasks;
	}

	free ((char *) stack);
	free ((char *) visited);

	return (ncomps);
}

/*
 * This routine recursively enumerates all non-trivial partitions of the
 * connected components into two disjoint groups.  Such a partition
 * defines a CUT of the original Steiner tree problem.  For each such
 * cut, we determine the set of all full-sets that span the cut.  If
 * this cutset has not been seen before, we add it to the cutlist.
 */

	static
	struct constraint *
enumerate_cuts (

int			i,		/* IN - current recursion level */
int			ntaken,		/* IN - number of components taken */
bool *			cstack,		/* IN - which components taken */
struct constraint *	cutlist,	/* IN - list of cutsets */
bitmap_t *		comps,		/* IN - connected components */
int			ncomps,		/* IN - number of components */
bitmap_t *		cut_terms,	/* IN - another temporary cut-set */
double *		x,		/* IN - current LP solution */
bitmap_t *		vert_mask,	/* IN - set of valid vertices */
bitmap_t *		edge_mask,	/* IN - set of valid hyperedges */
struct bbinfo *		bbip		/* IN - branch and bound info */
)
{
int			j;
int			kmasks;
struct gst_hypergraph *	cip;
bitmap_t *		bp1;
bitmap_t *		bp2;
bitmap_t *		bp3;

	if (i >= ncomps) {
		/* Base case. */
		if ((ntaken <= 0) OR (ntaken >= ncomps)) {
			/* Ignore trivial partitions... */
			return (cutlist);
		}

		cip = bbip -> cip;
		kmasks = cip -> num_vert_masks;

		/* Compute the set of terminals in the cut... */
		memset (cut_terms, 0, kmasks * sizeof (*cut_terms));
		for (j = 0; j < ncomps; j++) {
			if (NOT cstack [j]) continue;
			bp1 = cut_terms;
			bp2 = bp1 + kmasks;
			bp3 = &comps [j * kmasks];
			while (bp1 < bp2) {
				*bp1++ |= *bp3++;
			}
		}

		/* Add cutset to the list... */
		cutlist = _gst_add_cutset_to_list (cut_terms,
						   cutlist,
						   x,
						   vert_mask,
						   edge_mask,
						   bbip);

		return (cutlist);
	}

	/* Recurse here... */
	cstack [i] = FALSE;
	cutlist = enumerate_cuts (i + 1,
				  ntaken,
				  cstack,
				  cutlist,
				  comps,
				  ncomps,
				  cut_terms,
				  x,
				  vert_mask,
				  edge_mask,
				  bbip);

	cstack [i] = TRUE;
	cutlist = enumerate_cuts (i + 1,
				  ntaken + 1,
				  cstack,
				  cutlist,
				  comps,
				  ncomps,
				  cut_terms,
				  x,
				  vert_mask,
				  edge_mask,
				  bbip);

	return (cutlist);
}

/*
 * This routine generates a single cutset constraint for each
 * component (and merges appropriately similar cuts).  This method
 * can be MUCH faster than the enumeration of all cuts!
 */

	static
	struct constraint *
simple_cuts (

struct constraint *	cutlist,	/* IN - list of cutsets */
bitmap_t *		comps,		/* IN - connected components */
int			ncomps,		/* IN - number of components */
double *		x,		/* IN - current LP solution */
bitmap_t *		vert_mask,	/* IN - set of valid vertices */
bitmap_t *		edge_mask,	/* IN - set of valid hyperedges */
struct bbinfo *		bbip		/* IN - branch and bound info */
)
{
int			i;
int			kmasks;
struct gst_hypergraph *	cip;
bitmap_t *		cut_terms;

	cip = bbip -> cip;
	kmasks = cip -> num_vert_masks;

	for (i = 0; i < ncomps; i++) {
		/* The terminals in this component. */
		cut_terms = &comps [kmasks * i];

		/* Add this cutset to the list. */
		cutlist = _gst_add_cutset_to_list (cut_terms,
						   cutlist,
						   x,
						   vert_mask,
						   edge_mask,
						   bbip);
	}

	return (cutlist);
}

/*
 * This routine formulates the Linear Program used to find violated
 * cutset constraints.
 */

	void
_gst_build_cutset_separation_formulation (

bitmap_t *		vert_mask,	/* IN - set of valid vertices */
bitmap_t *		edge_mask,	/* IN - set of valid hyperedges */
struct bbinfo *		bbip		/* IN - branch and bound info */
)
{
int			i;
int			j;
int			k;
int			nedges;
int			nverts;
int			sumsizes;
int			num_used_trees;
int			num_nodes;
int			num_arcs;
struct gst_hypergraph *	cip;
struct cs_info *	csip;
struct flow_prob *	prob;
int			p0;
int			q0;
int			arc_num;
int *			outlist;
int *			inlist;
int *			counts;
int **			ptrs;
int *			vp1;
int *			vp2;
int *			ip;

	cip	= bbip -> cip;
	nedges	= cip -> num_edges;
	nverts	= cip -> num_verts;

	csip = NEW (struct cs_info);
	memset (csip, 0, sizeof (*csip));
	bbip -> csip = csip;

	/* Compute the sum of all hyperedge cardinalities. */
	sumsizes = 0;
	num_used_trees = 0;
	for (i = 0; i < nedges; i++) {
		if (NOT BITON (edge_mask, i)) continue;
		sumsizes += cip -> edge_size [i];
		++num_used_trees;
	}

#if 0
	gst_channel_printf (bbip -> params -> print_solve_trace,
			    " %% nedges = %d, nverts = %d, sumsizes = %d,"
			    " num_used_trees = %d\n",
			    nedges, nverts, sumsizes, num_used_trees);
#endif

	/* The network we construct contains two additional	*/
	/* nodes Pi and Qi per full set Fi (in addition to the	*/
	/* terminals), and the following set of edges for each	*/
	/* full set Fi:						*/
	/*							*/
	/*	- 1 edge from Pi to Qi.				*/
	/*	- 2 other edges per terminal Tj in Fi:		*/
	/*		- from Tj to Pi				*/
	/*		- from Qi to Tj				*/
	/*							*/
	/* NOTE: we have nodes Pi and Qi regardless of whether	*/
	/*	 or not full set Fi is in the problem -- this	*/
	/*	 makes it easier to compute the proper node	*/
	/*	 numbers for Pi and Qi...  We also have the	*/
	/*	 corresponding Pi -> Qi arc regardless of	*/
	/*	 whether or not full set Fi is in the problem.	*/

	num_nodes = nverts + 2 * nedges;
	num_arcs = nedges + 2 * sumsizes;

	prob = &(csip -> prob);

	prob -> num_nodes	= num_nodes;
	prob -> num_arcs	= num_arcs;
	prob -> source		= 0;		/* specified dynamically */
	prob -> sink		= 0;		/* specified dynamically */
	prob -> arc_src		= NEWA (num_arcs, int);
	prob -> arc_dst		= NEWA (num_arcs, int);
	prob -> capacity	= NEWA (num_arcs, double);

	csip -> arc_to_fset	= NEWA (num_arcs, int);

	/* The terminals become nodes 0 through nverts-1,	*/
	/* The "Pi" nodes start at nverts, and the "Qi" nodes	*/
	/* follow.						*/
	p0	= nverts;
	q0	= p0 + nedges;

	/* Generate the proper graph widget for each full set	*/
	/* that is properly a part of the problem...		*/
	/* The arcs numbered 0..nedges-1 are the Pi -> Qi arcs,	*/
	/* giving the total flow through each full set.		*/
	arc_num = 0;
	for (i = 0; i < nedges; i++) {
		/* Generate Pi -> Qi arc. */
		prob -> arc_src [arc_num] = p0 + i;
		prob -> arc_dst [arc_num] = q0 + i;
		csip -> arc_to_fset [arc_num] = i;
		++arc_num;
	}

	for (i = 0; i < nedges; i++) {
		if (NOT BITON (edge_mask, i)) continue;

		vp1 = cip -> edge [i];
		vp2 = cip -> edge [i + 1];
		while (vp1 < vp2) {
			j = *vp1++;

			/* Generate Tj -> Pi arc. */
			prob -> arc_src [arc_num] = j;
			prob -> arc_dst [arc_num] = p0 + i;
			csip -> arc_to_fset [arc_num] = i;
			++arc_num;

			/* Generate Qi -> Tj arc. */
			prob -> arc_src [arc_num] = q0 + i;
			prob -> arc_dst [arc_num] = j;
			csip -> arc_to_fset [arc_num] = i;
			++arc_num;
		}
	}

	FATAL_ERROR_IF (arc_num NE num_arcs);

	/* We have now specified the directed flow graph as a	*/
	/* list of directed arcs.  Time to construct the	*/
	/* adjacency lists -- for each node we build a list of	*/
	/* outgoing and incoming arc numbers.  Do the outgoing	*/
	/* lists first...					*/

	outlist	= NEWA (num_arcs, int);
	inlist	= NEWA (num_arcs, int);

	prob -> out = NEWA (num_nodes + 1, int *);
	prob -> in  = NEWA (num_nodes + 1, int *);

	counts	= NEWA (num_nodes, int);
	ptrs	= NEWA (num_nodes, int *);

	for (i = 0; i < num_nodes; i++) {
		counts [i] = 0;
	}
	for (i = 0; i < num_arcs; i++) {
		++(counts [prob -> arc_src [i]]);
	}
	ip = outlist;
	for (i = 0; i < num_nodes; i++) {
		ptrs [i] = ip;
		prob -> out [i] = ip;
		ip += counts [i];
	}
	prob -> out [i] = ip;
	for (i = 0; i < num_arcs; i++) {
		j = prob -> arc_src [i];
		ip = ptrs [j]++;
		*ip = i;
	}

	/* Now do the incoming arc lists... */
	for (i = 0; i < num_nodes; i++) {
		counts [i] = 0;
	}
	for (i = 0; i < num_arcs; i++) {
		++(counts [prob -> arc_dst [i]]);
	}
	ip = inlist;
	for (i = 0; i < num_nodes; i++) {
		ptrs [i] = ip;
		prob -> in [i] = ip;
		ip += counts [i];
	}
	prob -> in [i] = ip;
	for (i = 0; i < num_arcs; i++) {
		k = prob -> arc_dst [i];
		ip = ptrs [k]++;
		*ip = i;
	}

	/* Free temporary memory used to build things... */
	free ((char *) counts);
	free ((char *) ptrs);

	/* Initialize the buffers used to hold flow solutions */
	/* and temporary data... */
	_gst_create_flow_solution_data (prob, &(csip -> soln));
	_gst_create_flow_temp_data (prob, &(csip -> temp));
}

/*
 * This routine performs the separation procedure for cutset constraints.
 * Given an LP solution to the main problem, the task is to find one or
 * more cutsets constraints that are violated (have total weight less
 * than 1) -- or to prove that no such cutsets exist.  We do this by picking
 * the first valid terminal, and computing the max-flow/min-cut to each of
 * the other N-1 terminals.  If the flow is less that 1, we have a
 * violated constraint.  The actual cutset of full-sets consists of those
 * full-sets that span the cut.
 */

	struct constraint *
_gst_find_fractional_cutsets (

double *		x,		/* IN - LP solution to separate. */
bitmap_t *		vert_mask,	/* IN - set of valid vertices */
bitmap_t *		edge_mask,	/* IN - set of valid hyperedges */
struct bbinfo *		bbip		/* IN - branch and bound info. */
)
{
int			i;
int			kmasks;
int			nverts;
int			num_arcs;
struct gst_hypergraph *	cip;
struct cs_info *	csip;
double *		capacity;
int *			arc_to_fset;
int			I;
int			J;
struct constraint *	cutlist;
double			z;

	cip	= bbip -> cip;
	csip	= bbip -> csip;

	nverts	= cip -> num_verts;
	kmasks	= cip -> num_vert_masks;

	num_arcs	= csip -> prob.num_arcs;
	capacity	= csip -> prob.capacity;
	arc_to_fset	= csip -> arc_to_fset;

	/* Distribute the LP solution weight for each full set	*/
	/* to every edge in the max-flow network that belongs	*/
	/* to that full set.					*/
	for (i = 0; i < num_arcs; i++) {
		capacity [i] = x [arc_to_fset [i]];
	}

	/* Find first valid terminal. */
	I = 0;
	for (;;) {
		if (I >= nverts) {
			/* No valid terminal found... */
			FATAL_ERROR;
		}
		if (BITON (vert_mask, I)) break;
		++I;
	}

	cutlist = NULL;

	/* Loop through remaining terminals J.  Do a single max-flow/	*/
	/* min-cut problem for each one...				*/
	for (J = I + 1; J < nverts; J++) {
		if (NOT BITON (vert_mask, J)) continue;

		csip -> prob.source = J;
		csip -> prob.sink   = I;

		_gst_compute_max_flow (&(csip -> prob),
				       &(csip -> temp),
				       &(csip -> soln));

		z = csip -> soln.z;

		if ((1.0 - z) < FUZZ) continue;

		/* We have a violated constraint!  Keep	*/
		/* only bits for valid terminals...	*/
		for (i = 0; i < kmasks; i++) {
			csip -> soln.cut [i] &= vert_mask [i];
		}

		/* Add this cutset to the list! */
		cutlist = _gst_add_cutset_to_list (csip -> soln.cut,
						   cutlist,
						   x,
						   vert_mask,
						   edge_mask,
						   bbip);
	}

	return (cutlist);
}

/*
 * Free all memory associated with the cutset separation formulation.
 */

	void
_gst_free_cutset_separation_formulation (

struct cs_info *	csip		/* IN - cutset formulation to free */
)
{
struct flow_prob *	prob;

	if (csip EQ NULL) return;

	_gst_free_flow_temp_data (&(csip -> temp));
	_gst_free_flow_solution_data (&(csip -> soln));

	prob = &(csip -> prob);
	if (prob -> in NE NULL) {
		if (prob -> in [0] NE NULL) {
			free (prob -> in [0]);
		}
		free (prob -> in);
	}
	if (prob -> out NE NULL) {
		if (prob -> out [0] NE NULL) {
			free (prob -> out [0]);
		}
		free (prob -> out);
	}
	if (csip -> arc_to_fset NE NULL) {
		free (csip -> arc_to_fset);
	}
	if (prob -> capacity NE NULL) {
		free (prob -> capacity);
	}
	if (prob -> arc_dst NE NULL) {
		free (prob -> arc_dst);
	}
	if (prob -> arc_src NE NULL) {
		free (prob -> arc_src);
	}
	free (csip);
}
