/***********************************************************************

	$Id: localcut.h,v 1.8 2016/09/24 17:35:03 warme Exp $

	File:	localcut.h
	Rev:	e-3
	Date:	09/24/2016

	Copyright (c) 1997, 2016 by David M. Warme.  This work is
	licensed under a Creative Commons Attribution 4.0 International
	License.

************************************************************************

	Declarations pertaining to the local cut generator.

************************************************************************

	Modification Log:

	a-1:	12/07/99	warme
		: Created.
	e-1:	04/14/2015	warme
		: Changes for 5.0 release.
	e-2:	09/05/2016	warme
		: Change notices for 5.1 release.
	e-3:	09/24/2016	warme
		: Reorganize include files, apply prefixes.
		: Make features unconditional.

************************************************************************/

#ifndef LOCALCUT_H
#define	LOCALCUT_H

#include "bitmaskmacros.h"

struct bbinfo;
struct comp;
struct constraint;
struct gst_channel;

extern struct constraint *	_gst_find_local_cuts (
					double *		x,
					bitmap_t *		vert_mask,
					bitmap_t *		edge_mask,
					struct bbinfo *		bbip,
					struct constraint *	cp);
extern struct constraint *	_gst_find_local_cuts_in_component (
					struct comp *		comp,
					double *		x,
					bitmap_t *		vert_mask,
					bitmap_t *		edge_mask,
					struct bbinfo *		bbip,
					struct constraint *	cp);
extern void			_gst_print_forests (
					struct comp *		comp,
					bitmap_t *		flist,
					int			n,
					struct gst_channel *	chan);

#endif
