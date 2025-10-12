/***********************************************************************

	$Id: fputils.h,v 1.2 2016/09/30 20:03:32 warme Exp $

	File:	fputils.h
	Rev:	e-2
	Date:	09/30/2016

	Copyright (c) 1993, 2016 by David M. Warme.  This work is
	licensed under a Creative Commons Attribution 4.0 International
	License.

************************************************************************

	General declarations for the Steiner Tree program.

************************************************************************

	Modification Log:

	e-1:	09/24/2016	warme
		: Split off from steiner.h.
	e-2:	09/30/2016	warme
		: Added control of long double arithmetic.

************************************************************************/

#ifndef	FPUTILS_H
#define	FPUTILS_H

/*
 * The following structure is used to save floating point information.
 * Throughout most of GeoSteiner, this is just uninterpreted black-box
 * data.  The actual layout and meaning of this data is processor-specific.
 * We use an array of doubles to guarantee proper alignment, which might
 * be required by some implementations.  This structure must be large
 * enough to hold the info for all processors we support.
 */

struct fpsave {
	double	data [4];
};

/*
 * Function Prototypes.
 */

extern int	_gst_enable_long_double_precision ();
extern void	_gst_restore_long_double_precision (int prevState);

extern void	_gst_restore_floating_point_configuration (struct fpsave * sp);
extern void	_gst_save_floating_point_configuration (struct fpsave * sp);
extern void	_gst_set_floating_point_configuration (struct fpsave * sp);
extern void	_gst_store_double (double * dp, double x);

#endif
