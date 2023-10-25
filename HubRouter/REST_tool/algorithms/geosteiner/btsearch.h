/***********************************************************************

	$Id: btsearch.h,v 1.10 2016/09/24 18:00:34 warme Exp $

	File:	btsearch.h
	Rev:	e-3
	Date:	09/24/2016

	Copyright (c) 1993, 2016 by David M. Warme.  This work is
	licensed under a Creative Commons Attribution 4.0 International
	License.

************************************************************************

	Declarations for the special backtrack search used
	by the FST pruning code.

************************************************************************

	Modification Log:

	a-1:	11/30/2000	warme
		: Created.
	e-1:	04/14/2015	warme
		: Changes for 5.0 release.
	e-2:	09/05/2016	warme
		: Change notices for 5.1 release.
	e-3:	09/24/2016	warme
		: Reorganize include files, apply prefixes.

************************************************************************/

#ifndef BTSEARCH_H
#define	BTSEARCH_H

struct gst_solver;


/*
 * Function Prototypes
 */

extern void	_gst_backtrack_search (struct gst_solver * solver);

#endif
