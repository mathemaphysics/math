#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mls.h"

double gaussian( int dim_in, double a_in,  double *x_in )
{
	int i;
	double val = 0.0;
	for(i=0;i<dim_in;i++)
		val += x_in[i] * x_in[i] / a_in / a_in;
	return exp( -1.0 * val );
}

int main( int argc, char **argv )
{
	mls_t ml;
	int i,j;
	int pdim = binomial( 3 + 2, 2 );
	double pts[3*2] = { 0.01, 0.0232, 0.332, 0.1132, -0.2283, 0.9187 }; /* Two 2-dimensional points */
	double vals[3] = { 0.2325, -0.1128, -0.2234 };
	double vec[100],mat[100];
	double x[2] = { 0.232, -0.9238 };
	double xt[2] = { -0.2283, 0.9187 };
	double dlts[3] = { 1.0, 1.0, 1.0 };

	mls_init( &ml, 3, 2, 3, pts, dlts, vals, gaussian );
	mls_basis_evaluate( &ml, x, vec );
	mls_matrix_evaluate( &ml, 1.4, x, mat );
	for(i=0;i<pdim;i++)
		printf( "vec[%d] = %f\n", i, vec[i] );
	printf( "\n\n" );
	for(i=0;i<pdim;i++)
	{
		for(j=0;j<pdim;j++)
			printf( "%10.5f", mat[i*pdim+j] );
		printf( "\n" );
	}
	printf( "phi(idx) = %10.5f\n", mls_term_evaluate( &ml, 2, xt ) );

	printf( "total = %10.5f\n", mls_evaluate( &ml, xt ) );

	return 0;
}

