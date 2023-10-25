/***********************************************************************

	$Id: dt.h,v 1.1 2016/09/24 16:55:20 warme Exp $

	File:	dt.h
	Rev:	e-1
	Date:	09/24/2016

	Copyright (c) 1997, 2016 by David M. Warme & Martin Zachariasen.
	This work is licensed under a Creative Commons Attribution
	4.0 International License.

************************************************************************

	Delaunay triangulation.

************************************************************************

	Modification Log:

	e-1:	09/24/2016	warme
		: Split off from steiner.h.

************************************************************************/

#ifndef DT_H_INCLUDED
#define	DT_H_INCLUDED

struct pset;

extern void	_gst_delaunay_triangulation (
				struct pset *	pts,
				int *		numberofedges,
				int **		edgelist,
				int *		numberoftriangles,
				int **		trianglelist,
				int **		neighborlist);

#endif
