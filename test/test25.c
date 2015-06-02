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
	int radn = 10;

	/* Gaussian quadrature points */
	int quadn = 6;
	double qpts[6] = { -0.93246951, -0.66120939, -0.23861919, 0.23861919, 0.66120939, 0.93246951 };
	double qwts[6] = { 0.17132449, 0.360767157, 0.46791393, 0.46791393, 0.36076157, 0.17132449 };

	char gfn[1024];
	double sum,ssum;
	double x1,x2,res;
	int b_grd = 0; /* Load points from grid file */
	int npts,cgx = 50;
	double *pts,*dlt,*val,*gqw;
	FILE *fp;

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
	double *dr = (double*) malloc( dim * sizeof(double) );
	double qw,cr,ss;
	double tec[3],vec[3],wec[9],yec[3];
	
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
					ss = 0.0;
					for(m=0;m<dim;m++)
						ss += qbox[k*dim+m] * qbox[k*dim+m];
					ss = sqrt( ss );
					for(m=0;m<dim;m++)
						wec[k*dim+m] = qbox[k*dim+m] / ss;
				}

//				for(k=0;k<dim;k++)
//					index[k] = 0, size[k] = quadn - 1;
//				do
//				{
//					lens_gauss_point( dim, pts + j * dim, dlt[j], pts + i * dim, dlt[i],
//						cr, wec, index, qpts, qwts, tec, &qw );
//					for(k=0;k<dim;k++)
//						fprintf( stdout, "%15.7f", tec[k] );
//					fprintf( stdout, "\n" );
//				}
//				while( arraynext( (long) dim, size, index ) != -1 );

				/* Distance from center of j to center of i */
				ssum = 0.0;
				for(k=0;k<dim;k++)
					ssum += pow( pts[j*dim+k] - pts[i*dim+k], 2.0 );
				ssum = sqrt( ssum );

				/* Build the index array and iterate through all points in the plane of the lens */
				for(k=0;k<dim-1;k++)
					index[k] = 0;
				for(k=0;k<dim-1;k++)
					size[k] = quadn - 1;
				do
				{
					/* Calculate the limits for each dimension */
					for(k=0;k<dim-1;k++)
					{
						dr[k] = cr * cr;
						for(m=0;m<k;m++)
							dr[k] -= dr[m] * qpts[index[m]] * dr[m] * qpts[index[m]];
						dr[k] = sqrt( dr[k] );
						vec[k] = dr[k] * qpts[index[k]];
					}

					/* Now for the final dimension which is a function of sphere separation */
					x1 = dlt[j] * dlt[j];
					for(k=0;k<dim-1;k++)
						x1 -= vec[k] * vec[k];
					x1 = sqrt( x1 ); /* The right side of the left circle, j */
					x2 = dlt[i] * dlt[i];
					for(k=0;k<dim-1;k++)
						x2 -= vec[k] * vec[k];
					x2 = ssum - sqrt( x2 ); /* The left side of the right circle, i */

					/* Project the point onto the whole domain */
					for(m=0;m<dim;m++)
                                        	yec[m] = pts[j*dim+m]; /* Translate to the origin; the center of sphere j */
					for(k=1;k<dim;k++)
						for(m=0;m<dim;m++)
							yec[m] += vec[k-1] * wec[k*dim+m];
					for(k=0;k<quadn;k++)
					{
						for(m=0;m<dim;m++)
							tec[m] = yec[m] + ( 0.5 * ( x1 + x2 ) + 0.5 * ( x1 - x2 ) * qpts[k] ) * wec[m];
						qw = qwts[k];
						for(m=1;m<dim;m++)
							qw *= qwts[index[m-1]];
						/* The quadrature point is tec and the weight is in qw */
						for(m=0;m<dim;m++)
							fprintf( stdout, "%15.7f", tec[m] );
						fprintf( stdout, "\n" );
					}
				}
				while( arraynext( (long) dim - 1, size, index ) != -1 );
				return 0;
			}
		}
	}

	free( ctr );
	free( qbox );
	free( dr );
	free( pts );
	free( dlt );

	return 0;
}

