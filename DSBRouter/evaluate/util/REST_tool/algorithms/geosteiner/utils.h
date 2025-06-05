/***********************************************************************

	$Id: utils.h,v 1.1 2016/09/24 16:42:59 warme Exp $

	File:	utils.h
	Rev:	e-1
	Date:	09/24/2016

	Copyright (c) 1997, 2016 by David M. Warme & Martin Zachariasen.
	This work is licensed under a Creative Commons Attribution
	4.0 International License.

************************************************************************

	Various utility routines.

************************************************************************

	Modification Log:

	e-1:	09/24/2016	warme
		: Split off from steiner.h.

************************************************************************/

#ifndef UTILS_H_INCLUDED
#define	UTILS_H_INCLUDED

#include "bitmaskmacros.h"

struct gst_channel;

extern void	_gst_print_mask (struct gst_channel *	chan,
				 char *			msg,
				 bitmap_t *		bp1,
				 int			n);
extern char *	_gst_strdup (const char * s);

#endif
