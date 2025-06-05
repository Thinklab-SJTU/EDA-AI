/***********************************************************************

	$Id: hgmst.c,v 1.15 2016/09/24 17:38:05 warme Exp $

	File:	hgmst.c
	Rev:	e-3
	Date:	09/24/2016

	Copyright (c) 2002, 2016 by Pawel Winter, Martin Zachariasen.
	This work is licensed under a Creative Commons Attribution
	4.0 International License.

************************************************************************

	A simple function to solve hypergraph MST problems.

************************************************************************

	Modification Log:

	a-1:	04/30/02	benny
		: Created.
	e-1:	04/14/2015	warme
		: Changes for 5.0 release.
	e-2:	09/05/2016	warme
		: Change notices for 5.1 release.
	e-3:	09/24/2016	warme
		: Reorganize include files, fix -Wall issues.

************************************************************************/

#include "geosteiner.h"
#include "logic.h"
#include "prepostlude.h"
#include "steiner.h"

/*
 * Given hypergraph information this function finds the MST. If any of the
 * pointers length, nmstedges or mstedges is NULL then the corresponding
 * information is not returned.
 */

	int
gst_hgmst (

int		nterms,
int		nedges,
int *		edge_sizes,
int *		edges,
double *	weights,
double *	length,
int *		nmstedges,
int *		mstedges,
int *		soln_status,
gst_param_ptr	params
)
{
int		res;
gst_hg_ptr	H;
gst_solver_ptr	solver;

	GST_PRELUDE

	res = 0;

	/* Setup hypergraph */
	H = gst_create_hg (NULL);
	gst_set_hg_number_of_vertices (H, nterms);
	gst_set_hg_edges (H, nedges, edge_sizes, edges, weights);

	/* Setup solver and find solution */
	solver = gst_create_solver (H, params, NULL);

	res = gst_hg_solve (solver, NULL);

	if (res EQ 0) {
		/* No serious problems encountered */
		/* Update the wanted information */
		gst_hg_solution (solver, nmstedges, mstedges, length, 0);
		gst_get_solver_status (solver, soln_status);
	}

	/* Free the temporary structures */
	gst_free_hg (H);
	gst_free_solver (solver);

	GST_POSTLUDE

	return res;
}
