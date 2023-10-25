/***********************************************************************

	$Id: p1read.h,v 1.1 2016/09/24 16:46:25 warme Exp $

	File:	p1read.h
	Rev:	e-1
	Date:	09/24/2016

	Copyright (c) 1997, 2016 by David M. Warme & Martin Zachariasen.
	This work is licensed under a Creative Commons Attribution
	4.0 International License.

************************************************************************

	Utility functions pertaining to reading of problem instances.

************************************************************************

	Modification Log:

	e-1:	09/24/2016	warme
		: Split off from steiner.h.

************************************************************************/

#ifndef P1READ_H
#define	P1READ_H

struct gst_hypergraph;

extern int **	_gst_compute_basic_incompat (struct gst_hypergraph * cip);
extern void	_gst_init_term_trees (struct gst_hypergraph * cip);

#endif
