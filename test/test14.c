#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mqrbf.h"

int main()
{
	int i,j,k;
	int dim = 2;
	int gx = 20;
	int gy = 20;
	double dx = 1.0 / (double) ( gx - 1 );
	double dy = 1.0 / (double) ( gy - 1 );
	double sum = 0.0;
	int np = gx * gy;
	double *pts = (double*) malloc( np * dim * sizeof(double) );
	double *prm = (double*) malloc( np * sizeof(double) );
	double *val = (double*) malloc( np * sizeof(double) );
	double x[2];
	for(k=0,i=0;i<gx;i++)
	{
		for(j=0;j<gy;j++)
		{
			pts[k*dim+0] = (double) i * dx;
			pts[k*dim+1] = (double) j * dy;
			++k;
		}
	}
	for(i=0;i<np;i++)
		prm[i] = 10.0;
	for(i=0;i<np;i++)
		val[i] = exp( -1.0 * ( (pts[i*dim+0]-0.5) * (pts[i*dim+0]-0.5) + (pts[i*dim+1]-0.5) * (pts[i*dim+1]-0.5) ) / 0.05 );

	mqrbf_t mq;
	mqrbf_init( &mq, dim, np, pts, prm, val );
	mqrbf_generate( &mq );
	fprintf( stderr, "coeffs = " );
	for(i=0;i<np;i++)
		fprintf( stderr, "%15.7f", mq.coeffs[i] );
	fprintf( stderr, "\n" );

	for(i=0;i<np;i++)
	{
		sum = 0.0;
		for(j=0;j<np;j++)
		{
			x[0] = pts[i*dim+0] - pts[j*dim+0];
			x[1] = pts[i*dim+1] - pts[j*dim+1];
			sum += mq.coeffs[j] * mqf( dim, mq.prm[j], x );
		}
		printf( "%15.7f%15.7f%15.7f\n", pts[i*dim+0], pts[i*dim+1], sum );
	}

	return 0;
}

