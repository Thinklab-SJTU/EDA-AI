/***********************************************************************

	$Id: cutsubs.c,v 1.14 2016/09/24 17:53:19 warme Exp $

	File:	cutsubs.c
	Rev:	e-3
	Date:	09/24/2016

	Copyright (c) 1996, 2016 by David M. Warme.  This work is
	licensed under a Creative Commons Attribution 4.0 International
	License.

************************************************************************

	Processing of violated cutsets.

************************************************************************

	Modification Log:

	a-1:	01/16/2001
		: Split off from cutset.c.  Modified to take cut
		:  vertices rather than set of spanning edges.
	c-1:	08/05/2002	benny
		: A few changes for library release.
	e-1:	04/14/2015	warme
		: Changes for 5.0 release.
	e-2:	09/05/2016	warme
		: Change notices for 5.1 release.
	e-3:	09/24/2016	warme
		: Reorganize include files, apply prefixes.

************************************************************************/

#include "bb.h"
#include "constrnt.h"
#include "cutset.h"
#include "logic.h"
#include "steiner.h"

#if ZERO_WEIGHT_CUTSETS_METHOD EQ ZWC_METHOD_SUBTOURS
 #include "sec_heur.h"
#endif


/*
 * Global Routines
 */

struct constraint * _gst_add_cutset_to_list (
					bitmap_t *		cutset_tmp,
					struct constraint *	clist,
					double *		x,
					bitmap_t *		vert_mask,
					bitmap_t *		edge_mask,
					struct bbinfo *		bbip);


/*
 * Local Routines
 */

#if ZERO_WEIGHT_CUTSETS_METHOD EQ ZWC_METHOD_CUTSET
static bool		is_subset (bitmap_t *, bitmap_t *, int);
#endif

/*
 * This routine handles all of the details of adding a new cutset to the
 * given list of constraints.
 *
 * This is the smart version: we actually generate a pair of subtour
 * constraints instead of the weaker cutset constraint!
 */

#if ZERO_WEIGHT_CUTSETS_METHOD EQ ZWC_METHOD_SUBTOURS

	struct constraint *
_gst_add_cutset_to_list (

bitmap_t *		cutset_tmp,	/* IN - new cutset to add */
struct constraint *	clist,		/* IN - existing constraints */
double *		x,		/* IN - current LP solution */
bitmap_t *		vert_mask,	/* IN - set of valid terminals */
bitmap_t *		edge_mask,	/* IN - set of valid hyperedges */
struct bbinfo *		bbip		/* IN - branch and bound info */
)
{
int			i;
int			kmasks;
struct gst_hypergraph *	cip;

	cip = bbip -> cip;
	kmasks = cip -> num_vert_masks;

	clist = _gst_check_subtour (cutset_tmp,
				    clist,
				    x,
				    edge_mask,
				    bbip);

	/* Generate the complementary set of terminals... */
	for (i = 0; i < kmasks; i++) {
		cutset_tmp [i] ^= vert_mask [i];
	}

	/* Generate second subtour... */

	clist = _gst_check_subtour (cutset_tmp,
				    clist,
				    x,
				    edge_mask,
				    bbip);

	return (clist);
}

#endif

/*
 * This routine handles all of the details of adding a new cutset to the
 * given list of cutsets.  Any existing constraints that are looser
 * (supersets) are deleted.  The new cutset is ignored if it is looser
 * than any existing cutset.
 */

#if ZERO_WEIGHT_CUTSETS_METHOD EQ ZWC_METHOD_CUTSET

	struct constraint *
_gst_add_cutset_to_list (

bitmap_t *		cutset,		/* IN - new cutset to add */
struct constraint *	cutlist,	/* IN - list to add to */
double *		x,		/* IN - current LP solution */
bitmap_t *		vert_mask,	/* IN - set of valid terminals */
bitmap_t *		edge_mask,	/* IN - set of valid hyperedges */
struct bbinfo *		bbip		/* IN - branch and bound info */
)
{
int			i;
int			j;
int			nedges;
int			nmasks;
int			num_in_cut;
int			count;
struct gst_hypergraph *	cip;
int *			vp1;
int *			vp2;
struct constraint *	p;
struct constraint **	hookp;
bitmap_t *		cut_edges;
double			z;

	cip	= bbip -> cip;
	nedges	= cip -> num_edges;
	nmasks	= cip -> num_edge_masks;

	cut_edges = NEWA (nmasks, bitmap_t);
	memset (cut_edges, 0, nmasks * sizeof (*cut_edges));

	count = 0;
	z = 0.0;
	for (i = 0; i < cip -> num_edges; i++) {
		if (NOT BITON (edge_mask, i)) continue;
		num_in_cut = 0;
		vp1 = cip -> edge [i];
		vp2 = cip -> edge [i + 1];
		while (vp1 < vp2) {
			j = *vp1++;
			if (BITON (cutset, j)) {
				++num_in_cut;
			}
		}
		if (num_in_cut <= 0) {
			/* this hyperedge resides entirely	*/
			/* outside of the cut...  doesn't span!	*/
			continue;
		}
		if (num_in_cut >= cip -> edge_size [i]) {
			/* this hyperedge resides entirely	*/
			/* within the cut...  doesn't span!	*/
			continue;
		}
		SETBIT (cut_edges, i);
		++count;
		z += x [i];
	}

	/* Check for an all-zero cutset.  These occasionally	*/
	/* happen because of numeric issues...			*/

	if (count <= 0) {
		/* Empty cutset!  OOOPS! */
#if 1
		gst_channel_printf (bbip -> params -> print_solve_trace,
				    "WARNING!  empty cutset!\n");
#endif
		free ((char *) cut_edges);
		return (cutlist);
	}

	if (z >= 1.0 - FUZZ) {
#if 1
		gst_channel_printf (bbip -> params -> print_solve_trace,
				    "WARNING!  bogus cutset!\n");
#endif
		free ((char *) cut_edges);
		return (cutlist);
	}

	/* If this new cutset is a superset of an existing one,	*/
	/* then there is nothing to add, and nothing to delete.	*/
	for (p = cutlist; p NE NULL; p = p -> next) {
		if (is_subset (p -> mask, cut_edges, nmasks)) {
			free (cut_edges);
			return (cutlist);
		}
	}

	/* Delete all current cutsets which have this new one	*/
	/* as a subset.						*/
	hookp = &cutlist;
	while ((p = *hookp) NE NULL) {
		if (p -> type NE CT_CUTSET) {
			hookp = &(p -> next);
		}
		else if (is_subset (cut_edges, p -> mask, nmasks)) {
			*hookp = p -> next;
			free ((char *) (p -> mask));
			free ((char *) p);
		}
		else {
			hookp = &(p -> next);
		}
	}

	p = NEW (struct constraint);
	p -> next	= NULL;
	p -> iteration	= 0;
	p -> type	= CT_CUTSET;
	p -> mask	= cut_edges;
	*hookp = p;

	return (cutlist);
}

#endif

/*
 * This routine returns TRUE if-and-only-if the first bit-mask is
 * a subset of the second.
 */

#if ZERO_WEIGHT_CUTSETS_METHOD EQ ZWC_METHOD_CUTSET

	static
	bool
is_subset (

bitmap_t *	bp1,		/* IN - first set. */
bitmap_t *	bp2,		/* IN - second set. */
int		nmasks		/* IN - number of masks in set. */
)
{
int		i;

	for (i = 0; i < nmasks; i++) {
		if ((*bp1 & *bp2) NE *bp1) return (FALSE);
		++bp1;
		++bp2;
	}

	return (TRUE);
}

#endif
