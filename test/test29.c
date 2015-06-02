#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "points.h"

int main( int argc, char **argv )
{
	int i,a,n,*id,dim,npts;
	double *pts,x[2] = { 1.0, 1.0 };
	double r = 0.12;
	kdnode_t root;

	load_points( argv[1], &dim, &npts, &pts, 0 );
	kdtree_build_average( &root, dim, npts, pts, 0, 0 );
	kdtree_walk( &root );

	/* Range query test */
	id = NULL, a = 0, n = 0;
	kdtree_range_query( &root, x, r, &id, &a, &n );
	fprintf( stderr, "n = %5d\n", n );
	for(i=0;i<n;i++)
		fprintf( stderr, "%5d: %5d x = %15.7f%15.7f\n", i, id[i], pts[id[i]*dim+0], pts[id[i]*dim+1] );

	return 0;
}

