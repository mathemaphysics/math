#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "lagrange.h"
#include "cgrad.h"

int main()
{
	int i,j,m,ret,n = 100,nv = 2,p = 50;
	double val,*amat,*bv,*vv;
	double *xv,*alpha,*beta,*gamma,a,b;

	amat = (double*) malloc( n * n * sizeof(double) );
	bv = (double*) malloc( n * sizeof(double) );
	xv = (double*) malloc( n * sizeof(double) );
	vv = (double*) malloc( nv * n * sizeof(double) );
	alpha = (double*) malloc( ( p + 1 ) * sizeof(double) );
	beta = (double*) malloc( ( p + 1 ) * sizeof(double) );
	gamma = (double*) malloc( ( p + 1 ) * sizeof(double) );

	/* Build a fucking matrix to test */
	srand((unsigned)time(0));
	for(i=0;i<n*n;i++)
		amat[i] = ( rand() > RAND_MAX / 2 ? 1.0 : -1.0 ) * (double) rand() / (double) RAND_MAX;
	for(i=0;i<n;i++) /* Push it away from singularity */
		amat[i*n+i] = 2.0 + (double) rand() / (double) RAND_MAX;
	for(i=0;i<n;i++)
		bv[i] = ( rand() > RAND_MAX / 2 ? 1.0 : -1.0 ) * (double) rand() / (double) RAND_MAX;
	for(i=0;i<nv;i++)
		for(j=0;j<n;j++)
			vv[i*n+j] = ( rand() > RAND_MAX / 2 ? 1.0 : -1.0 ) * (double) rand() / (double) RAND_MAX;
	for(i=0;i<n;i++) 
                if( i == 4 ) 
                        vv[i] = 1.0;
                else 
                        vv[i] = 0.0;
	for(i=0;i<n;i++)
                if( i == 2 )
                        vv[n+i] = 1.0;
                else
                        vv[n+i] = 0.0;
	for(i=0;i<n;i++)
        	xv[i] = ( rand() > RAND_MAX / 2 ? 1.0 : -1.0 ) * (double) rand() / (double) RAND_MAX;
	fprintf( stdout, "# name: A\n# type: matrix\n# rows: %d\n# columns: %d\n", n, n );
	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			fprintf( stdout, "%15.7f", amat[i*n+j] );
	bicgstabp( n, 0, amat, bv, nv, vv, xv, 3 * n, 1e-9, 1, &ret );
	fprintf( stderr, "Returned %d\n", ret );

	/* Print solution */
	fprintf( stderr, "xv = " );
	for(i=0;i<n;i++)
		fprintf( stderr, "%15.7f", xv[i] );
	fprintf( stderr, "\n\n" );

	fprintf( stderr, "ax = " );
	for(i=0;i<n;i++)
	{
		val = 0.0;
		for(j=0;j<n;j++)
			val += amat[i*n+j] * xv[j];
		fprintf( stderr, "%15.7f", val );
	}
	fprintf( stderr, "\n\n" );

	fprintf( stderr, "bv = " );
	for(i=0;i<n;i++)
		fprintf( stderr, "%15.7f", bv[i] );
	fprintf( stderr, "\n" );

	free( amat );
	free( bv );
	free( xv );

	/* Build the lagrange basis here */
	long dim = 2;
	long deg = 3;
	long pdim = binomial( dim + deg, dim );
	double div = 1.0 / (double) deg;
	double nodes[pdim*dim];
	lagrange_t basis;

	/* Initialize */
	lagrange_init( &basis, dim, deg, NULL );
	lagrange_set_nodes_standard( &basis );
	lagrange_generate( &basis );

	/* Output the basis functions */
	long nd = 50;
	double dx = 1.0 / (double) ( nd - 1 );
	double pt[dim];
	//for(m=0;m<pdim;m++)
	//	for(i=0;i<nd;i++)
	//{
	//	pt[0] = (double) i * dx;
	//	for(j=0;j<nd-i;j++)
	//	{
	//		pt[1] = (double) j * dx;
	//		printf( "%15.7f%15.7f%15.7f\n", pt[0], pt[1], lagrange_evaluate( &basis, m, pt ) );
	//	}
	//}

	return 0;
}

