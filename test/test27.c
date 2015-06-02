#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include "points.h"
#include "punity.h"

int main( int argc, char **argv )
{
	/* Basic variables */
	int i,j,k,p,m;
	int dim = 2;
	int deg = 3;
	int npts;
	int radn = 10;
	FILE *fp;

	/* Gaussian quadrature points */
	int quadn = 6;
	double qpts[6] = { -0.93246951, -0.66120939, -0.23861919, 0.23861919, 0.66120939, 0.93246951 };
	double qwts[6] = { 0.17132449, 0.360767157, 0.46791393, 0.46791393, 0.36076157, 0.17132449 };

	/* Particle configuration */
	char gfn[1024];
	double *pts,*dlt;
	punity_t pt;

	/* Configuration switches */
	int b_grd = 0;

	/* Take in the options */
	int optc;
	while( ( optc = getopt( argc, argv, "d:g:r:" ) ) != -1 )
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
			case '?':
			default:
				fprintf( stderr, "I don't understand the jibberish you are spouting. Goodbye.\n" );
				return 0;
				break;
		}
	}

	/* Load points from file if switch is set */
	if( b_grd )
	{
		fprintf( stderr, "Loading points from file... " );
		load_points( gfn, &dim, &npts, &pts, 0 );
		fprintf( stderr, "dim = %d npts = %d. Done.\n", dim, npts );
	}
	else
	{
		fprintf( stderr, "No point cloud given. Exiting.\n" );
		return 0;
	}

	/* Generate supports */
	dlt = (double*) malloc( npts * sizeof(double) );
	generate_cloud_supports_min( dim, npts, pts, radn, 1.05, dlt, 0.001 );
	punity_init( &pt, dim, npts, pts, dlt, &cubic_window, NULL );

	/* Variables for integration */
	double *ctr = (double*) malloc( dim * sizeof(double) );
	double *qbox = (double*) malloc( dim * dim * sizeof(double) );
	double sum,cr,qw,nqbox[dim*dim],qx[dim];
	
	/* Set up the index */
	long index[dim],size[dim];

	/* Intersect the spheres */
	for(i=0;i<npts;i++)
	{
		for(j=0;j<npts;j++)
		{
			if( j == i )
				continue;
			if( sphere_intersection( dim, pts + j * dim, dlt[j], pts + i * dim, dlt[i],
				ctr, &cr, qbox ) == 1 )
			{
				/* Create a normalized version of qbox */
				for(k=0;k<dim;k++)
				{
					sum = 0.0;
					for(m=0;m<dim;m++)
						sum += qbox[k*dim+m] * qbox[k*dim+m];
					sum = sqrt( sum );
					for(m=0;m<dim;m++)
						nqbox[k*dim+m] = qbox[k*dim+m] / sum;
				}
				for(k=0;k<dim;k++)
					index[k] = 0, size[k] = quadn - 1;
				do
				{
					lens_gauss_point( dim, pts + j * dim, dlt[j], pts + i * dim, dlt[i],
						cr, nqbox, index, qpts, qwts, qx, &qw );
					for(k=0;k<dim;k++)
						fprintf( stdout, "%15.7f", qx[k] );
					sum = punity_evaluate( &pt, i, qx ) * punity_evaluate( &pt, j, qx );
					fprintf( stdout, "%15.7f\n", sum );
				}
				while( arraynext( (long) dim, size, index ) != -1 );
				return 0;
			}
		}
	}

	free( ctr );
	free( qbox );
	free( pts );
	free( dlt );

	return 0;
}

