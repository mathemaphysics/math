#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "lagrange.h"

int main( int argc, char **argv )
{
	long i,j,n=5;
	lagrange_t lag;
	double nodes[5] = { 0.0, 0.333, 0.666, 1.0, 1.333 };

	/* Initialize lagrange basis */
	lagrange_init( &lag, 1, n - 1, nodes );
	lagrange_generate( &lag );

	long idx = 2;
	double x;
	double x_i = 0.0;
	double x_f = 1.333;
	long ndiv = 100;
	double stp = ( x_f - x_i ) / (double) ndiv;
	double vals[101],sums[101];
	for(i=0;i<n;i++)
		vals[i] = pow( lag.nodes[i], 5.0 );
	for(i=0;i<=100;i++)
		sums[i] = 0.0;

	for(j=0;j<n;j++)
		for(i=0;i<=100;i++)
	{
		x = x_i + (double) i * stp;
		sums[i] += vals[j] * lagrange_evaluate( &lag, j, &x );
	}

	FILE *fp = fopen( "test9.out", "w" );
	for(i=0;i<=ndiv;i++)
	{
		x = x_i + (double) i * stp;
		fprintf( fp, "%f %f\n", x, sums[i] );
	}
	fclose( fp );

	return 0;
}
