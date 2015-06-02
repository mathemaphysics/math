#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include "nuclei.h"
#include "mls.h"
#include "linalg.h"
#include "trees.h"
#include "sparse.h"
#include "points.h"

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

int on_boundary( double *pt_in )
{
        if( fabs( pt_in[0] ) < 0.04 || fabs( pt_in[1] ) < 0.04 || fabs( pt_in[0] - 5.0 ) < 0.04 || fabs( pt_in[1] - 5.0 ) < 0.04 )
                return 1;
        return 0;
}

double coulomb_potential( int dim_in, double *x_in, int nn_in, double *nuc_in )
{
	int i,k;
	double sum;
	double tot = 0.0;
	for(k=0;k<nn_in;k++)
	{
		for(sum=0.0,i=0;i<dim_in;i++)
			sum += pow( x_in[i] - nuc_in[k*dim_in+i], 2.0 );
		tot -= 1.0 / sqrt( sum );
	}
	return tot;
}

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
	double *pts,*dlt,*val,*gqw;

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
					fprintf( stdout, "%15.7f\n", qw );
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

