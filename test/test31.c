#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <unistd.h>
#include <time.h>
#include "lagrange.h"
#include "cgrad.h"

int main( int argc, char **argv )
{
	int i,j,k,m,ret,optc,ns,nc = 20,n = 400,nv = 20,p = 200,vb = 0;
	double val,*amat,*bmat,*bv,*tmpv,*ev;
	double *xv,*alpha,*beta,*gamma,*uu,*vv,a,b;
	long *ia,*ja,*ib,*jb;
	double *A,*B,*V,*W;

	while( ( optc = getopt( argc, argv, "v" ) ) != -1 )
	{
		switch( optc )
		{
			case 'v':
				vb = 1;
				break;
			case '?':
			default:
				fprintf( stderr, "I don't understand your silly request. Exiting.\n" );
				break;
		}
	}

	/* Sparse stuff */
	A = (double*) malloc( 3 * n * sizeof(double) );
	ja = (long*) malloc( 3 * n * sizeof(long) );
	ia = (long*) malloc( ( n + 1 ) * sizeof(long) );
	B = (double*) malloc( 3 * n * sizeof(double) );
	jb = (long*) malloc( 3 * n * sizeof(long) );
	ib = (long*) malloc( ( n + 1 ) * sizeof(long) );
	//srand((unsigned)time(0));
	k = 0, m = 0;
	ia[0] = 0;
	ja[k] = 0, A[k++] = (double) rand() / (double) RAND_MAX + 5.0;
	ja[k] = 1, A[k++] = (double) rand() / (double) RAND_MAX;
	ib[0] = 0;
	jb[m] = 0, B[m++] = (double) rand() / (double) RAND_MAX + 4.0;
	jb[m] = 1, B[m++] = (double) rand() / (double) RAND_MAX;
	for(i=1;i<n-1;i++)
	{
		ia[i] = k;
		ja[k] = i - 1, A[k++] = (double) rand() / (double) RAND_MAX;
		ja[k] = i, A[k++] = (double) rand() / (double) RAND_MAX + 2.0;
		ja[k] = i + 1, A[k++] = (double) rand() / (double) RAND_MAX;
		ib[i] = m;
		jb[m] = i - 1, B[m++] = (double) rand() / (double) RAND_MAX;
		jb[m] = i, B[m++] = (double) rand() / (double) RAND_MAX + 4.0;
		jb[m] = i + 1, B[m++] = (double) rand() / (double) RAND_MAX;
		
	}
	ia[n-1] = k;
	ja[k] = n - 2, A[k++] = (double) rand() / (double) RAND_MAX;
	ja[k] = n - 1, A[k++] = (double) rand() / (double) RAND_MAX + 3.0;
	ia[n] = k;
	ib[n-1] = m;
	jb[m] = n - 2, B[m++] = (double) rand() / (double) RAND_MAX;
	jb[m] = n - 1, B[m++] = (double) rand() / (double) RAND_MAX + 2.0;
	ib[n] = m;

	smsave( "amat.oct", n, n, ia, ja, A, 1 );
	smsave( "bmat.oct", n, n, ib, jb, B, 1 );

	/* Generate a set of constraints */
	V = (double*) malloc( nc * n * sizeof(double) );
	W = (double*) malloc( n * n * sizeof(double) );
	for(i=0;i<nc;i++)
	{
		nrandv( n, V + i * n );
		project( n, V + i * n, i, V );
		normalize( n, V + i * n );
	}
	fprintf( stdout, "# name: V\n# type: matrix\n# rows: %d\n# columns: %d\n", n, nc );
	for(i=0;i<n;i++)
	{
		for(j=0;j<nc;j++)
			fprintf( stdout, "%15.12f", V[j*n+i] );
		fprintf( stdout, "\n" );
	}

	/* Generate the projector */
	srand((unsigned)time(0)); /* Generate a new basis for same subspace and matrices A, B to test sensitivity to basis */
	for(i=0;i<n-nc;i++)
	{
		nrandv( n, W + i * n );
		project( n, W + i * n, nc, V );
		project( n, W + i * n, i, W );
		normalize( n, W + i * n );
	}
	fprintf( stdout, "# name: W\n# type: matrix\n# rows: %d\n# columns: %d\n", n, n - nc );
	for(i=0;i<n;i++)
	{
		for(j=0;j<n-nc;j++)
			fprintf( stdout, "%15.12f", W[j*n+i] );
		fprintf( stdout, "\n" );
	}

	/* Dense operators */
	amat = (double*) malloc( n * n * sizeof(double) );
	bmat = (double*) malloc( n * n * sizeof(double) );
	bv = (double*) malloc( n * sizeof(double) );
	xv = (double*) malloc( n * sizeof(double) );
	vv = (double*) malloc( ( nv > p ? nv : p ) * n * sizeof(double) );
	alpha = (double*) malloc( ( n + 2 ) * sizeof(double) );
	beta = (double*) malloc( ( n + 2 ) * sizeof(double) );
	gamma = (double*) malloc( ( n + 2 ) * sizeof(double) );
	uu = (double*) malloc( ( n + 2 ) * n * sizeof(double) );
	ev = (double*) malloc( 2 * p * sizeof(double) ); /* Complex eigenvalues as well */

	/* Build a fucking matrix to test */
	//srand((unsigned)time(0));
	for(i=0;i<n*n;i++)
		amat[i] = ( rand() > RAND_MAX / 2 ? 1.0 : -1.0 ) * (double) rand() / (double) RAND_MAX,
		bmat[i] = ( rand() > RAND_MAX / 2 ? 1.0 : -1.0 ) * (double) rand() / (double) RAND_MAX;
	for(i=0;i<n;i++) /* Push it away from singularity */
		amat[i*n+i] = 40.0 + (double) rand() / (double) RAND_MAX,
		bmat[i*n+i] = 40.0 + (double) rand() / (double) RAND_MAX,
		bv[i] = (double) rand() / (double) RAND_MAX;
	fprintf( stdout, "# name: A\n# type: matrix\n# rows: %d\n# columns: %d\n", n, n );
	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			fprintf( stdout, "%15.7f", amat[i*n+j] );
	fprintf( stdout, "\n# name: B\n# type: matrix\n# rows: %d\n# columns: %d\n", n, n );
	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			fprintf( stdout, "%15.7f", bmat[i*n+j] );
	fprintf( stdout, "\n# name: b\n# type: matrix\n# rows: %d\n# columns: 1\n", n );
	for(i=0;i<n;i++)
		fprintf( stdout, "%15.7f\n", bv[i] );
	fprintf( stderr, "p = %d\n", p );
	//gulanczos( n, p, amat, bmat, alpha, beta, gamma, NULL, NULL, 1e-9, 200, vb, &ret );
	//sgulanczos( n, p, ia, ja, A, ib, jb, B, alpha, beta, gamma, NULL, NULL, 1e-9, 200, vb, &ret );
	//ulanczos( n, p, amat, alpha, beta, gamma, NULL, NULL, &ret );
	csgulanczos( n, p, ia, ja, A, ib, jb, B, alpha, beta, gamma, NULL, NULL, nc, V, 1e-10, 200, vb, &ret );
	fprintf( stderr, "csgulanczos ret = %d\n", ret );
	fprintf( stdout, "\n# name: T\n# type: matrix\n# rows: %d\n# columns: %d\n", p, p );
	for(i=0;i<p;i++)
	{
		for(j=0;j<p;j++)
		{
			if( abs( i - j ) > 1 )
				fprintf( stdout, "%15.7f", 0.0 );
			else
			{
				if( j > i )
					fprintf( stdout, "%15.7f", gamma[i+1] );
				else if( j < i )
					fprintf( stdout, "%15.7f", beta[i] );
				else
					fprintf( stdout, "%15.7f", alpha[i+1] );
			}
		}
		fprintf( stdout, "\n" );
	}
	treigtridqds( p, alpha + 1, beta + 1, gamma + 1, 20000, ev, &ret, &ns );
	fprintf( stderr, "Reduction = %d Steps = %d\n", ret, ns );
	for(i=0;i<p;i++)
		fprintf( stderr, "%15.7f+%15.7fi\n", ev[2*i+0], ev[2*i+1] );

	free( amat );
	free( bv );
	free( xv );

	return 0;
}

