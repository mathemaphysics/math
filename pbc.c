#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "math_config.h"

#if HAVE_DOUBLE == 1  /* This is set in the configure script with --enable-double=yes */
	typedef double t_real;
#else
	typedef float t_real;
#endif

/**
 * Calculate the PBC vector distance between
 * two points.
 * @param x_in First vector
 * @param y_in Second vector
 * @param bx_in Box side legnths
 * @param z_out PBC vector distance between x_in and y_in
 */
void pbc_vec_real3( t_real *x_in, t_real *y_in, t_real *bx_in, t_real *z_out )
{
	t_real xx,yy;

	/* Take linear diff */
	z_out[0] = y_in[0] - x_in[0];
	z_out[1] = y_in[1] - x_in[1];
	z_out[2] = y_in[2] - x_in[2];

	/* Calculate anint stuff */
	xx = z_out[0] / bx_in[0];
#if HAVE_DOUBLE == 0
	canintf_( (float*)(&xx), (float*)(&yy) );
#else
	canint_( (double*)(&xx), (double*)(&yy) );
#endif
	z_out[0] = z_out[0] - bx_in[0] * yy;
	xx = z_out[1] / bx_in[1];
#if HAVE_DOUBLE == 0
	canintf_( (float*)(&xx), (float*)(&yy) );
#else
	canint_( (double*)(&xx), (double*)(&yy) );
#endif
	z_out[1] = z_out[1] - bx_in[1] * yy;
	xx = z_out[2] / bx_in[2];
#if HAVE_DOUBLE == 0
	canintf_( (float*)(&xx), (float*)(&yy) );
#else
	canint_( (double*)(&xx), (double*)(&yy) );
#endif
	z_out[2] = z_out[2] - bx_in[2] * yy;
}

/**
 * Calculate the closest distance between two
 * points in a rectangular domain assuming the
 * domain is defined by coordinates in bx_in.
 * @param x_in First vector
 * @param y_in Second vector
 * @param bx_in Box dimensions for PBC
 * @return Distance between x_in and y_in
 */
t_real pbc_dist_real3( t_real *x_in, t_real *y_in, t_real *bx_in )
{
	t_real xx,yy,vec[3];
	t_real dr;

	/* Take linear diff */
	vec[0] = y_in[0] - x_in[0];
	vec[1] = y_in[1] - x_in[1];
	vec[2] = y_in[2] - x_in[2];

	/* Calculate anint stuff */
	xx = vec[0] / bx_in[0];
#if HAVE_DOUBLE == 0
	canintf_( (float*)(&xx), (float*)(&yy) );
#else
	canint_( (double*)(&xx), (double*)(&yy) );
#endif
	vec[0] = vec[0] - bx_in[0] * yy;
	xx = vec[1] / bx_in[1];
#if HAVE_DOUBLE == 0
	canintf_( (float*)(&xx), (float*)(&yy) );
#else
	canint_( (double*)(&xx), (double*)(&yy) );
#endif
	vec[1] = vec[1] - bx_in[1] * yy;
	xx = vec[2] / bx_in[2];
#if HAVE_DOUBLE == 0
	canintf_( (float*)(&xx), (float*)(&yy) );
#else
	canint_( (double*)(&xx), (double*)(&yy) );
#endif
	vec[2] = vec[2] - bx_in[2] * yy;
#if HAVE_DOUBLE == 0
	dr = (t_real) sqrtf( vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2] );
#else
    dr = (t_real) sqrt( vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2] );
#endif

	/* Return some stuff */
	return dr;
}

// vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab:smarttab
