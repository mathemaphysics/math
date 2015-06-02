#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linalg.h"
#include "config.h"

typedef struct
{
	int np;
	int dim;
	double *pts;
	double *prm;
	double *vals;
	double *coeffs;
	double *vdm;
} mqrbf_t;

double mqf( int dim_in, double prm_in, double *x_in )
{
	int i;
	double sum = 0.0;
	for(i=0;i<dim_in;i++)
		sum += x_in[i] * x_in[i];
	return sqrt( 1.0 + prm_in * sum );
}

int mqrbf_init( mqrbf_t *obj_in, int dim_in, int np_in, double *pts_in, double *prm_in, double *vals_in )
{
	int i,j;
	if( obj_in == NULL )
		return -1;
	obj_in->np = np_in;
	obj_in->dim = dim_in;
	obj_in->pts = (double*) malloc( np_in * dim_in * sizeof(double) );
	if( obj_in->pts == NULL )
		return -2;
	for(i=0;i<np_in;i++)
		for(j=0;j<dim_in;j++)
			obj_in->pts[i*dim_in+j] = pts_in[i*dim_in+j];
	obj_in->prm = (double*) malloc( np_in * sizeof(double) );
	if( obj_in->prm == NULL )
		return -3;
	for(i=0;i<np_in;i++)
		obj_in->prm[i] = prm_in[i];
	obj_in->vals = (double*) malloc( np_in * sizeof(double) );
        if( obj_in->vals == NULL )
                return -4;
	for(i=0;i<np_in;i++)
		obj_in->vals[i] = vals_in[i];
	obj_in->coeffs = (double*) malloc( np_in * sizeof(double) );
        if( obj_in->coeffs == NULL )
                return -5;
	obj_in->vdm = (double*) malloc( np_in * np_in * sizeof(double) );
	if( obj_in->vdm == NULL )
		return -6;

	return 0;
}

int mqrbf_free( mqrbf_t *obj_in )
{
	free( obj_in->pts );
	free( obj_in->prm );
	free( obj_in->vals );
	free( obj_in->coeffs );
	free( obj_in->vdm );

	return 0;
}

int mqrbf_generate( mqrbf_t *obj_in )
{
	int i,j,k;
	double x[obj_in->dim];

	for(i=0;i<obj_in->np;i++)
	{
		for(j=0;j<obj_in->np;j++)
		{
			for(k=0;k<obj_in->dim;k++)
				x[k] = obj_in->pts[i*obj_in->dim+k] - obj_in->pts[j*obj_in->dim+k];
			obj_in->vdm[i*obj_in->np+j] = mqf( obj_in->dim, obj_in->prm[i], x );
		}
	}

	int ret,n = obj_in->np,*ipiv = (int*) malloc( n * sizeof(int) );
	j = 1;
	for(i=0;i<obj_in->np;i++)
		obj_in->coeffs[i] = obj_in->vals[i];
	dgesv_( &n, &j, obj_in->vdm, &n, ipiv, obj_in->coeffs, &n, &ret );

	return 0;
}

