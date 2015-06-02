#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main()
{
	int i;
	int dim,np;
	double *pts;
	double *rad;

	load_points( "cloud.out", &dim, &np, &pts, 0 );
	rad = (double*) malloc( np * sizeof(double) );
	generate_cloud_supports( dim, np, pts, 6, 1.2, rad );
	fprintf( stderr, "rad = " );
	for(i=0;i<np;i++)
		fprintf( stderr, "%15.7f", rad[i] );
	fprintf( stderr, "\n" );

	free( rad );
	free( pts );
}

