#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mls.h"

double gaussian( int dim_in, double a_in, double *x_in )
{
	int i;
	double val = 0.0;
	for(i=0;i<dim_in;i++)
		val += x_in[i] * x_in[i];
	return exp( -1.0 * val / 0.005 / a_in / a_in );
}

double step( int dim_in, double a_in, double *x_in )
{
	if( fabs( x_in[0] ) < 0.2 )
		return 1.0;
	else
		return 0.0;
}

int main( int argc, char **argv )
{
	mls_t ml;
	int i,j;
	double pts[10*1] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
	double vals[10] = { 0.0, 0.0, 0.0, 0.05, 0.10, 0.12, 0.14, 0.12, 0.10, 0.05 };
	double xt[1] = { 0.511 };
	double dlts[10] = { 1.0, 3.4, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

	mls_init( &ml, 10, 1, 3, pts, dlts, vals, gaussian );
	double step = 0.9 / (double) 100;
	FILE *fp = fopen( "plot.dat", "w" );
	for(i=0;i<100;i++)
	{
		xt[0] = (double) i * step;
		fprintf( fp, "%f %f\n", xt[0], mls_evaluate( &ml, xt ) );
	}
	fclose( fp );

	fp = fopen( "basis.dat", "w" );
	for(i=0;i<10;i++)
	//i = 1;
	for(j=0;j<100;j++)
	{
		xt[0] = (double) j * step;
		fprintf( fp, "%f %f\n", xt[0], mls_term_evaluate( &ml, i, xt ) );
	}
	fclose( fp );

	return 0;
}

