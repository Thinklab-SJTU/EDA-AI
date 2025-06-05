/***********************************************************************

	$Id: solver.h,v 1.25 2016/09/24 17:05:43 warme Exp $

	File:	solver.h
	Rev:	e-3
	Date:	09/24/2016

	Copyright (c) 2002, 2016 by Pawel Winter, Martin Zachariasen.
	This work is licensed under a Creative Commons Attribution
	4.0 International License.

************************************************************************


************************************************************************

	Modification Log:

	a-1:	04/30/2002	benny
		: Created.
	b-1:	04/24/2014	warme
		: Declare gst_set_hook_root_node_completed().
	e-1:	04/14/2015	warme
		: Changes for 5.0 release.
	e-2:	09/05/2016	warme
		: Change notices for 5.1 release.
	e-3:	09/24/2016	warme
		: Reorganize include files, add prefixes.

************************************************************************/

#ifndef SOLVER_H
#define SOLVER_H

#include "bitmaskmacros.h"
#include "cputime.h"

struct gst_hypergraph;
struct gst_param;
struct gst_solver;

typedef void (*root_node_completed_func)(struct gst_solver *,
					 struct gst_hypergraph *,
					 struct gst_param *,
					 double *);

struct gst_solver {
	struct gst_hypergraph *
			H;
	struct gst_param *
			params;

	struct ubinfo *	ubip;		/* Info for upper bound heuristic */

	struct bbinfo *	bbip;		/* Branch-and-cut info */

	double		upperbound;
	double		lowerbound;

	cpu_time_t	t0;		/* Start time for gst_hg_solve() */
	cpu_time_t	p2time;		/* Total time used in phase 2 */

	int		solsize;	/* Allocated size of solutions[] array */
	int		nsols;		/* Number of valid solutions[] */
	int		feasible_updates;  /* The number of updates done */
	struct gst_hg_solution * solutions; /* Set of K best feasible solutions */

	struct gst_proplist * proplist;

	volatile int	preempt;	/* Breaks the solution process ASAP. */
					/* Its value indicates the reason    */
					/* for exiting. See below.	     */

	/* Latest optimizer used in solution process (b-and-c,backtrack,none)*/
	int		latest_optimizer_run;

	/* Latest hypergraph version used in solution process */
	int		solution_version;

	/* Callback hooks */
	root_node_completed_func	hook_root_node_completed;
};

#define PREEMPT_SOLVER(solver,reason)		\
	if ((solver) -> preempt EQ 0) {		\
		(solver) -> preempt = reason;	\
	}

/*
 * This structure is simply used to remember a solution. The edge list should
 * be sorted to ease any comparison between solutions.
 */

struct gst_hg_solution {
	int		nedges;	/* number of edges in solution */
	int		hash;	/* hash value based on the set of edges */
	int *		edges;	/* set of edges in solution */
	bitmap_t *	edge_mask; /* set of edges in solution */
	double		length;	/* solution length */
};


extern int	gst_set_hook_root_node_completed (struct gst_solver *,
						  root_node_completed_func);
extern bool	_gst_update_best_solution_set (
				struct gst_solver *	solver,
				double *		x,
				int			nedges,
				int *			edges,
				bitmap_t *		smt);
extern void	_gst_update_solver_properties (struct gst_solver * solver);

#endif
