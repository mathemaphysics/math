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

int indicator_cloud( double *x_in )
{
	return 1;
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
	rkp_t rk;
	int i,j,k,p,m;
	int dim = 2;
	int deg = 3;
	int pdim = binomial( dim + deg, dim );
	int ul = 10;
	int nmin = 100;
	int radn = 10;
	int term;
	int nnz;
	int quadn = 6;
	//int quadn = 3;
	double qpts[6] = { -0.93246951, -0.66120939, -0.23861919, 0.23861919, 0.66120939, 0.93246951 };
	double qwts[6] = { 0.17132449, 0.360767157, 0.46791393, 0.46791393, 0.36076157, 0.17132449 };
	//double qpts[3] = { -0.774597, 0.0, 0.774597 };
	//double qwts[3] = { 0.555556, 0.888889, 0.555556 };
	double zrt = 1e-6;
	double ssum,sum,tum,scl = 0.3;
	double xx,yy,x1,x2,res;
	double ptmp[dim];
	double tol = zrt;
	int b_grd = 0; /* Load points from grid file */
	int b_lmat = 0; /* Load the matrices amat, bmat, dmat from files amat.sparse, etc. */
	char gfn[1024];
	int nb;
	int npts;
	double *pts,*dlt,*val,*gqw;
	FILE *fp;

	/* Take in the options */
	int optc;
	while( ( optc = getopt( argc, argv, "Lg:d:s:t:r:c" ) ) != -1 )
        {
                switch( optc )
                {
			case 'L':
				b_lmat = 1;
				break;
			case 'r':
				radn = atoi( optarg );
				break;
			case 'c':
				break;
			case 'd':
				deg = atoi( optarg );
				break;
			case 's':
				scl = atof( optarg );
				break;
			case 't':
				tol = atof( optarg );
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

	/* Set up the plot points */
        int pgx = 50;
        int pgy = 50;
        double x = 5.0;
        double y = 5.0;
        double pdx = x / (double) ( pgx - 1 );
        double pdy = y / (double) ( pgy - 1 );
        int nppts = pgx * pgy;
        double *ppts = (double*) malloc( nppts * dim * sizeof(double) );
        for(i=0;i<pgx;i++)
                for(j=0;j<pgy;j++)
                        ppts[i*pgy*dim+j*dim+0] = (double) i * pdx,
                        ppts[i*pgy*dim+j*dim+1] = (double) j * pdy;

	/* Generate supports */
	dlt = (double*) malloc( npts * sizeof(double) );
	val = (double*) malloc( npts * sizeof(double) );
	gqw = (double*) malloc( npts * sizeof(double) );
	generate_cloud_supports_min( 2, npts, pts, radn, 1.05, dlt, 0.001 );
	for(i=0;i<npts;i++) 
                gqw[i] = 1.0;

	/* Initialize the RKP basis */
	rkp_init( &rk, npts, dim, deg, pts, dlt, gqw, val, cubic, 1.0 );
	rkp_basis_generate_scaled_const( &rk );
	rkp_wavelet_basis_generate_scaled_const( &rk );

	/* Variables for integration */
	double *ctr = (double*) malloc( dim * sizeof(double) );
	double *qbox = (double*) malloc( dim * dim * sizeof(double) );
	double cr,dr;
	double qdx = 2.0 / (double) ( quadn - 1 );
	double vec[2],wec[2];
	int cgx = 50;
	char fname[512];
	
	/* Form the Laplace operator correctly with wavelets */
	int ord[2] = { 2, 4 }; /* Indexes which contain the xx and yy second derivatives */

	/* Intersect the spheres */
	for(i=0;i<npts;i++)
	{
		for(j=0;j<npts;j++)
		{
			if( sphere_intersection( dim, pts + j * dim, dlt[j] * rk.wrad, pts + i * dim, dlt[i] * rk.wrad,
				ctr, &cr, qbox ) == 1 )
			{
				ssum = 0.0;
				for(k=0;k<dim;k++)
					ssum += pow( pts[j*dim+k] - pts[i*dim+k], 2.0 );
				ssum = sqrt( ssum );

				/* Output the two circles */
				for(k=0;k<cgx;k++)
				{
					vec[0] = pts[j*dim+0] + rk.wrad * dlt[j] * cos( (double) k / (double) ( cgx - 1 ) * 2.0 * M_PI );
					vec[1] = pts[j*dim+1] + rk.wrad * dlt[j] * sin( (double) k / (double) ( cgx - 1 ) * 2.0 * M_PI );
					fprintf( stdout, "%15.7f%15.7f\n", vec[0], vec[1] );
				}
				for(k=0;k<cgx;k++)
                                {
                                        vec[0] = pts[i*dim+0] + rk.wrad * dlt[i] * cos( (double) k / (double) ( cgx - 1 ) * 2.0 * M_PI );
                                        vec[1] = pts[i*dim+1] + rk.wrad * dlt[i] * sin( (double) k / (double) ( cgx - 1 ) * 2.0 * M_PI );
					fprintf( stdout, "%15.7f%15.7f\n", vec[0], vec[1] );
                                }

				/* Iterate through all quadrature points to build stiffness matrix */
				for(k=0;k<quadn;k++)
				{
					xx = cr * qpts[k];

					/* Min and max calculated for x */
                                        x1 = sqrt( pow( rk.wrad * dlt[j], 2.0 ) - xx * xx );
                                        x2 = ssum - sqrt( pow( rk.wrad * dlt[i], 2.0 ) - xx * xx );

					/* Plot a sequence from minimum to maximum for this given value xx along the secondary axis */
					for(m=0;m<quadn;m++)
					{
						wec[0] = qbox[0] / sqrt( qbox[0] * qbox[0] + qbox[1] * qbox[1] );
						wec[1] = qbox[1] / sqrt( qbox[0] * qbox[0] + qbox[1] * qbox[1] );
						yy = ( x1 + x2 ) / 2.0 + ( x1 - x2 ) / 2.0 * qpts[m];
						vec[0] = pts[j*dim+0] + xx / cr * qbox[1*dim+0] + yy * wec[0];
						vec[1] = pts[j*dim+1] + xx / cr * qbox[1*dim+1] + yy * wec[1];
						fprintf( stdout, "%15.7f%15.7f\n", vec[0], vec[1] );
					}
				}
				if( j > 5 )
					return 0;
			}
			
		}
	}

	free( pts );
	free( dlt );
	free( val );
	free( gqw );
	free( ppts );

	return 0;
}

