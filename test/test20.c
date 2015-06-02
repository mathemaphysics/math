#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include "mls.h"
#include "linalg.h"
#include "trees.h"

double cubic( int dim_in, double a_in, double *x_in )
{
	int i;
	double z = 0.0;
	for(i=0;i<dim_in;i++)
		z += x_in[i] * x_in[i];
	z = sqrt( z ) / a_in * 2.0;
	if( z > 2.0 )
		return 0.0;
	if( z > 1.0 )
		return ( 2.0 - z ) * ( 2.0 - z ) * ( 2.0 - z ) / 6.0 / a_in / a_in;
	if( z >= 0.0 )
		return ( 4.0 - 6.0 * z * z + 3.0 * z * z * z ) / 6.0 / a_in / a_in;
}

double density( double *x_in )
{
	return 0.01;
}

int main( int argc, char **argv )
{
	/* Basic variables */
	rkp_t rk;
	int i,j,k,p,m;
	int dim = 2;
	int deg = 2;
	int pdim;
	int radn = 10;
	int b_grd = 0; /* Load points from grid file */
	int b_iso = 0;
	char gfn[1024];
	double sum,x[2] = { 1.0, 1.0 };

	/* Take in the options */
	int optc;
	while( ( optc = getopt( argc, argv, "g:d:r:n" ) ) != -1 )
        {
                switch( optc )
                {
			case 'r':
				radn = atoi( optarg );
				break;
			case 'd':
				deg = atoi( optarg );
				break;
			case 'g':
				b_grd = 1;
				strcpy( gfn, optarg );
				break;
			case 'n':
				b_iso = 1;
				break;
			case '?':
			default:
				fprintf( stderr, "I don't understand the jibberish you are spouting. Goodbye.\n" );
				return 0;
				break;
		}
	}

	/* Allocate and generate grid points */
	int npts;
	double *pts,*dlt,*val,*gqw;

	/* Load points from file if switch is set */
	if( b_grd )
	{
		fprintf( stderr, "Loading points from file... " );
		load_points( gfn, &dim, &npts, &pts );
		fprintf( stderr, "dim = %d npts = %d. Done.\n", dim, npts );
	}
	pdim = binomial( deg + dim, dim );

	/* Generate supports */
	dlt = (double*) malloc( npts * sizeof(double) );
	val = (double*) malloc( npts * sizeof(double) );
	gqw = (double*) malloc( npts * sizeof(double) );
	for(i=0;i<npts;i++) 
                gqw[i] = 1.0; 
	generate_cloud_supports( 2, npts, pts, radn, 1.05, dlt );

	/* Initialize the RKP basis */
	rkp_init( &rk, npts, dim, deg, pts, dlt, gqw, val, cubic, 1.0 );
	fprintf( stderr, "dlt = %15.7f\n", dlt[93] );
	rkp_matrix_evaluate_scaled_const( &rk, pts + 93 * dim, dlt[93] );

	for(i=0;i<pdim;i++)
	{
		for(j=0;j<pdim;j++)
			fprintf( stderr, "%20.12f", rk.grm[i*pdim+j] );
		fprintf( stderr, "\n" );
	}

	free( dlt );
	free( val );
	free( gqw );

	return 0;
}

