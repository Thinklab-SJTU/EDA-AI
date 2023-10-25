/***********************************************************************

	$Id: sec2.h,v 1.7 2016/09/24 17:15:04 warme Exp $

	File:	sec2.h
	Rev:	e-3
	Date:	09/24/2016

	Copyright (c) 1997, 2016 by David M. Warme.  This work is
	licensed under a Creative Commons Attribution 4.0 International
	License.

************************************************************************

	Declarations for the deterministic separation procedure
	for the "generalized SEC's" that uses a reduction to min-cut
	in a bipartite network.

************************************************************************

	Modification Log:

	a-1:	05/16/1997	warme
		: Created.
	e-1:	04/14/2015	warme
		: Changes for 5.0 release.
	e-2:	09/05/2016	warme
		: Change notices for 5.1 release.
	e-3:	09/24/2016	warme
		: Reorganize include files, apply prefixes.

************************************************************************/

#ifndef SEC2_H
#define	SEC2_H

#include "bitmaskmacros.h"

struct bbinfo;
struct comp;
struct constraint;

/*
 * Function Prototypes
 */

extern struct constraint *	_gst_sec_flow_separator (
					struct comp **		comp_hookp,
					double *		x,
					bitmap_t *		edge_mask,
					struct bbinfo *		bbip,
					struct constraint *	cp);

#endif
