#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "combinadic-cuda.h"
#include "points-cuda.h"
#include "trees-cuda.h"

#define ZERO_THRESHOLD 1e-6

#define CUDA_RAND_MAX 2147483648

typedef struct
{
	int dim;
	int type;
	double *orig;
	double *params;
} shape_t;

__device__ void domain_indicator( double *x, shape_t *sh, int *ret )
{
	int i;
	double sum;

	switch( sh->type )
	{
		case 0:
			*ret = 1;
			for(i=0;i<sh->dim;i++)
			{
				if( x[i] < sh->orig[i] || x[i] > sh->params[0] + sh->orig[i] )
				{
					*ret = 0;
					break;
				}
			}
			break;
		case 1:
			*ret = 1;
			for(i=0;i<sh->dim;i++)
			{
				if( x[i] < sh->orig[i] || x[i] > sh->params[i] + sh->orig[i] )
				{
					*ret = 0;
					break;
				}
			}
			break;
		case 2:
			sum = 0.0;
			for(i=0;i<sh->dim;i++);
				sum += pow( x[i] - sh->orig[i], 2.0 );
			if( sum < sh->params[0] * sh->params[0] )
				*ret = 1;
			else
				*ret = 0;
			break;
		case 3:
			break;
		default:
			*ret = 0;
			break;
	}
}

/**
 * A random number generator for use inside device functions
 */
__device__ unsigned int cuda_rand( unsigned int *m_z, unsigned int *m_w )
{
	*m_z = 36969 * ( (*m_z) & 65535 ) + ( (*m_z) >> 16 );
	*m_w = 18000 * ( (*m_w) & 65535 ) + ( (*m_w) >> 16 );

	return ( ( (*m_z) << 16 ) + (*m_w) ) % CUDA_RAND_MAX;
}

/**
 * Calculate a set of dim_in - 1 axes perpendicular to vec_in
 */
__device__ void local_axes( int dim_in, double *vec_in, double *basis_out )
{
	int i,j,k;
	unsigned int mz,mw;
	double sum;

	sum = 0.0;
	for(i=0;i<dim_in;i++)
		basis_out[i] = vec_in[i], sum += basis_out[i] * basis_out[i];
	sum = sqrt( sum );
	for(i=0;i<dim_in;i++)
		basis_out[i] /= sum;
	mz = 150; mw = 40;
	for(i=1;i<dim_in;i++)
	{
		/* This is a horrible idea */
		for(j=0;j<dim_in;j++)
			basis_out[i*dim_in+j] = 0.5 - ( (double) cuda_rand( &mz, &mw ) / (double) CUDA_RAND_MAX );
		sum = 0.0;
		for(j=0;j<dim_in;j++)
			sum += basis_out[i*dim_in+j] * basis_out[i*dim_in+j];
		sum = sqrt( sum );
		for(j=0;j<dim_in;j++)
			basis_out[i*dim_in+j] /= sum;
		for(j=0;j<i;j++)
		{
			sum = 0.0;
			for(k=0;k<dim_in;k++)
				sum += basis_out[i*dim_in+k] * basis_out[j*dim_in+k];
			for(k=0;k<dim_in;k++)
				basis_out[i*dim_in+k] -= sum * basis_out[j*dim_in+k];
		}
		sum = 0.0;
		for(j=0;j<dim_in;j++)
			sum += basis_out[i*dim_in+j] * basis_out[i*dim_in+j];
		sum = sqrt( sum );
		for(j=0;j<dim_in;j++)
			basis_out[i*dim_in+j] /= sum;
	}
}

/**
 * Calculate a rectangular box about the intersection of two spheres
 * of potentially different radii
 * @param dim_in Dimension in which the spheres live
 * @param ctr1_in Center of the first sphere
 * @param rad1_in Radius of the first sphere
 * @param ctr2_in Center of the second sphere
 * @param rad2_in Radius of the second sphere
 * @param ctr_out Center of the disc containing the intersection
 * @param rad_out Radius of the disc containing the intersection
 * @param qbox_out Contains a local basis covering the intersection of the spheres
 * @return Returns 1 if spheres intersect, 0 if not
 */
__device__ int sphere_intersection( int dim_in, double *ctr1_in, double rad1_in, double *ctr2_in, double rad2_in, double *ctr_out, double *rad_out, double *qbox_out )
{
	int i,j;
	double r,rr,sum;
	double *ax = (double*) malloc( dim_in * sizeof(double) );

	/* Calculate the origin and axis of the cylinder */
	sum = 0.0;
	for(i=0;i<dim_in;i++)
		ax[i] = ctr2_in[i] - ctr1_in[i], sum += ax[i] * ax[i];
	sum = sqrt( sum );
	if( sum > rad1_in + rad2_in )
		return 0;

	/* If circles are identical */
	if( sum < ZERO_THRESHOLD )
	{
		/* Then build a box around the smaller sphere */
		for(i=0;i<dim_in;i++)
			ctr_out[i] = ctr1_in[i];
		*rad_out = ( rad1_in < rad2_in ) ? rad1_in : rad2_in;
		for(i=0;i<dim_in;i++)
			for(j=0;j<dim_in;j++)
				qbox_out[i*dim_in+j] = ( i == j ) ? *rad_out : 0.0;
		return 1;
	}

	/* Otherwise */
	for(i=0;i<dim_in;i++)
		ax[i] /= sum;
	r = rad1_in + rad2_in - sum;
	for(i=0;i<dim_in;i++)
		ctr_out[i] = ctr1_in[i] + ( sum - rad2_in ) * ax[i];
	for(i=0;i<dim_in;i++)
		ax[i] *= r;

	/* Calculate the radius of the cylinder */
	rr = sqrt( fabs( rad1_in * rad1_in - rad2_in * rad2_in ) );
	if( sum < rr )
		*rad_out = ( rad1_in < rad2_in ) ? rad1_in : rad2_in;
	else
		*rad_out = sqrt( ( -sum + rad2_in - rad1_in ) * ( -sum - rad2_in + rad1_in )
			* ( -sum + rad2_in + rad1_in ) * ( sum + rad2_in + rad1_in ) ) / 2.0 / sum;

	/* Generate the local coordinates */
	local_axes( dim_in, ax, qbox_out );
	for(i=0;i<dim_in;i++)
		qbox_out[i] = ax[i];
	for(i=1;i<dim_in;i++)
		for(j=0;j<dim_in;j++)
			qbox_out[i*dim_in+j] *= (*rad_out);

	/* Move ctr_out to the middle of the box */
	for(i=0;i<dim_in;i++)
		qbox_out[i] *= 0.5, ctr_out[i] += qbox_out[i];

	free( ax );

	return 1;
}

/**
 * Generate a quadrature point in the lens of spherical intersection
 */
__device__ void lens_gauss_point( int dim_in,
		double *ctr1_in, double rad1_in,
		double *ctr2_in, double rad2_in,
		double cr_in, double *nqbox_in, long *index_in,
		double *qpts_in, double *qwts_in,
		double *qp_out, double *qw_out )
{
	int i,j;
	double ssum,x1,x2;
	double *dr,*vec,*wec;

	/* Must allocate memory via malloc in __device__ code */
	dr = (double*) malloc( ( dim_in - 1 ) * sizeof(double) );
	vec = (double*) malloc( dim_in * sizeof(double) );
	wec = (double*) malloc( dim_in * sizeof(double) );

	/* Calculate distance between centers */
	ssum = 0.0;
	for(i=0;i<dim_in;i++)
		ssum += pow( ctr2_in[i] - ctr1_in[i], 2.0 );
	ssum = sqrt( ssum );
	
	/* Calculate the limits for each dimension */
	for(i=0;i<dim_in-1;i++)
	{
		dr[i] = cr_in * cr_in;
		for(j=0;j<i;j++)
			dr[i] -= dr[j] * qpts_in[index_in[j]] * dr[j] * qpts_in[index_in[j]];
		dr[i] = sqrt( dr[i] );
		vec[i] = dr[i] * qpts_in[index_in[i]];
	}

	/* Now for the final dimension which is a function of sphere separation */
	x1 = rad1_in * rad1_in;
	for(i=0;i<dim_in-1;i++)
		x1 -= vec[i] * vec[i];
	x1 = sqrt( x1 );
	x2 = rad2_in * rad2_in;
	for(i=0;i<dim_in-1;i++)
		x2 -= vec[i] * vec[i];
	x2 = ssum - sqrt( x2 );

	/* Project the point onto the whole domain */
	for(i=0;i<dim_in;i++)
      		wec[i] = ctr1_in[i]; /* Translate to the origin, ctr1_in */
	for(i=1;i<dim_in;i++)
		for(j=0;j<dim_in;j++)
			wec[j] += vec[i-1] * nqbox_in[i*dim_in+j];

	/* Index along the axis is given by index_in[dim_in-1] */
	for(i=0;i<dim_in;i++)
		qp_out[i] = wec[i] + ( 0.5 * ( x1 + x2 ) + 0.5 * ( x1 - x2 ) * qpts_in[index_in[dim_in-1]] ) * nqbox_in[i];
	*qw_out = qwts_in[index_in[dim_in-1]];
	for(i=1;i<dim_in;i++)
		*qw_out *= qwts_in[index_in[i-1]] * 2.0 * dr[i-1];
	*qw_out *= fabs( x2 - x1 );
}

__device__ void sphere_gauss_point( int dim_in, double *ctr_in, double rad_in, double *nqbox_in, long *index_in, double *qpts_in, double *qwts_in, double *qp_out, double *qw_out )
{
	int i,j;
	double *dr,*vec;

	/* Allocate directly */
	dr = (double*) malloc( dim_in * sizeof(double) );
	vec = (double*) malloc( dim_in * sizeof(double) );

	/* Calculate dimension limits */
	for(i=0;i<dim_in;i++)
	{
		dr[i] = rad_in * rad_in;
		for(j=0;j<i;j++)
			dr[i] -= dr[j] * qpts_in[index_in[j]] * dr[j] * qpts_in[index_in[j]];
		dr[i] = sqrt( dr[i] );
		vec[i] = dr[i] * qpts_in[index_in[i]];
	}

	/* Project the point onto the entire domain by affine transformation */
	for(i=0;i<dim_in;i++)
		qp_out[i] = ctr_in[i];
	for(i=0;i<dim_in;i++)
		for(j=0;j<dim_in;j++)
			qp_out[j] += vec[i] * nqbox_in[i*dim_in+j];
	*qw_out = 1.0;
	for(i=0;i<dim_in;i++)
		(*qw_out) *= 2.0 * dr[i] * qwts_in[index_in[i]];
}

