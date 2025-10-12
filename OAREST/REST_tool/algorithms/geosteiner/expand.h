/***********************************************************************

	$Id: expand.h,v 1.1 2016/09/24 16:54:22 warme Exp $

	File:	expand.h
	Rev:	e-1
	Date:	09/24/2016

	Copyright (c) 1996, 2016 by David M. Warme.  This work is
	licensed under a Creative Commons Attribution 4.0 International
	License.

************************************************************************

	Expand logical constraints into the proper physical encoding
	as used in the LP.

************************************************************************

	Modification Log:

	e-1:	09/24/2016	warme
		: Split off from constrnt.h.

************************************************************************/

#ifndef EXPAND_H
#define	EXPAND_H

#include "bitmaskmacros.h"


/*
 * Function Prototypes
 */

struct constraint;
struct gst_hypergraph;
struct rcoef;

extern struct rcoef *	_gst_expand_constraint (
					struct constraint *	lcp,
					struct rcoef *		cp,
					bitmap_t *		edge_mask,
					struct gst_hypergraph *	cip);

#endif
