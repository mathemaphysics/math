#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include "mls.h"
#include "linalg.h"
#include "trees.h"

double gaussian( int dim_in, double a_in, double *x_in )
{
	int i;
	double z = 0.0;
	for(i=0;i<dim_in;i++)
		z += x_in[i] * x_in[i];
	z /= ( a_in * a_in );
	return exp( -0.5 * z / a_in / a_in ) / a_in / sqrt( 2.0 * M_PI );
}

double cubic( int dim_in, double a_in, double *x_in )
{
	int i;
	double z = 0.0;
	for(i=0;i<dim_in;i++)
		z += x_in[i] * x_in[i];
	z = sqrt( z ) / a_in;
	if( z > 2.0 )
		return 0.0;
	if( z > 1.0 )
		return ( 2.0 - z ) * ( 2.0 - z ) * ( 2.0 - z ) / 6.0;
	if( z >= 0.0 )
		return ( 4.0 - 6.0 * z * z + 3.0 * z * z * z ) / 6.0;
}

void cubic_grad( int dim_in, double a_in, double *x_in, double *grad_out )
{
	int i;
	double r = 0.0;
	double sum = 0.0;
	for(i=0;i<dim_in;i++)
		r += x_in[i] * x_in[i];
	r = sqrt( r ) / a_in;
	for(i=0;i<dim_in;i++)
		grad_out[i] = x_in[i] / a_in;
	sum = 0.0;
	if( r > 2.0 )
		sum = 0.0;
	else if( r > 1.0 )
		sum = -0.5 * ( 2.0 - r ) * ( 2.0 - r ) / r;
	else if( r >= 0.0 )
		sum = 1.5 * r - 2.0;
	for(i=0;i<dim_in;i++)
		grad_out[i] *= sum;
}

double rhsf( int dim_in, double *x_in )
{
	return -1.0;
}

int main( int argc, char **argv )
{
	/* Basic variables */
	rkp_t rk;
	int i,j,k,p;
	int dim = 2;
	int deg = 3;
	int pdim = binomial( dim + deg, dim );
	int pfout = 0;
	int grd = 0;
	int ul = 10;
	int wav = 0;
	int wvo = deg - 1;
	double sum,scl = 0.3;
	double xx;
	double ptmp[dim];
	double *grad = (double*) malloc( dim * sizeof(double) );

	/* Take in the options */
	int optc;
	while( ( optc = getopt( argc, argv, "f:d:s:r:wg" ) ) != -1 )
        {
                switch( optc )
                {
			case 'g':
				grd = 1;
				break;
			case 'f':
				pfout = atoi( optarg );
				break;
			case 'd':
				deg = atoi( optarg );
				break;
			case 's':
				scl = atof( optarg );
				break;
			case 'r':
				wvo = atoi( optarg );
				break;
			case 'w':
				wav = 1;
				break;
			case '?':
			default:
				fprintf( stderr, "I don't understand the jibberish you are spouting. Goodbye.\n" );
				return 0;
				break;
		}
	}

	/* Allocate and generate grid points */
	double x = 1.0;
	double y = 1.0;
	int gx = 10; /* These are the number of points needed; NOT the spacing divisor */
	int gy = 10;
	double amat[(gx*gy)*(gx*gy)],bmat[(gx*gy)*(gx*gy)],smat[(gx*gy)*(gx*gy)],rhs[gx*gy];
	double dx = x / (double) ( gx - 1 );
	double dy = y / (double) ( gy - 1 );
	int npts = gx * gy;
	double *pts = (double*) malloc( npts * dim * sizeof(double) );
	for(i=0;i<gx;i++)
		for(j=0;j<gy;j++)
			pts[i*gy*dim+j*dim+0] = (double) i * dx,
			pts[i*gy*dim+j*dim+1] = (double) j * dy;
	double *dlt = (double*) malloc( npts * sizeof(double) );
	for(i=0;i<npts;i++)
		dlt[i] = scl;
	double *val = (double*) malloc( npts * sizeof(double) );
	double *gqw = (double*) malloc( npts * sizeof(double) );
	for(i=0;i<npts;i++)
		gqw[i] = 1.0;

	/* Fix the quadrature weights to start with; hopefully this will help */
	for(i=0;i<gx*gy;i++)
	{
		/* Point centered at pts + i * dim */
		for(j=0;j<gx*gy;j++)
		{
			xx = pow( pts[i*dim+0] - pts[j*dim+0], 2.0 ) + pow( pts[i*dim+1] - pts[j*dim+1], 2.0 );
			if( fabs( xx ) < scl )
				gqw[i] += 1.0;
		}
		gqw[i] /= (double) ( gx * gy );
	}

	/* Showing updated quadrature weights */
	fprintf( stderr, "Quadrature: " );
	for(i=0;i<gx*gy;i++)
		fprintf( stderr, "%15.7f", gqw[i] );
	fprintf( stderr, "\n" );

	/* Initialize the RKP basis */
	rkp_init( &rk, npts, dim, deg, pts, dlt, gqw, val, cubic, 2.0 ); /* The cubic has radius of 2.0 as defined above */
	rkp_basis_generate( &rk );
	rkp_wavelet_basis_generate( &rk );
	rk.wfsg = cubic_grad;

	/* Set up the plot points */
	int pgx = 50;
	int pgy = 50;
	double pdx = x / (double) ( pgx - 1 );
	double pdy = y / (double) ( pgy - 1 );
	int nppts = pgx * pgy;
	double *ppts = (double*) malloc( nppts * dim * sizeof(double) );
	for(i=0;i<pgx;i++)
                for(j=0;j<pgy;j++)
                        ppts[i*pgy*dim+j*dim+0] = (double) i * pdx,
                        ppts[i*pgy*dim+j*dim+1] = (double) j * pdy;

	int on_boundary( double *pt_in )
	{
		if( fabs( pt_in[0] ) < 1e-8 || fabs( pt_in[1] < 1e-8 || fabs( pt_in[0] - 1.0 ) < 1e-8 || fabs( pt_in[1] - 1.0 ) < 1e-8 ) )
			return 1;
		return 0;
	}

	/* Create the matrices */
	int ord[2] = { 2, 4 }; /* Indexes which contain the xx and yy second derivatives */
	for(i=0;i<gx*gy*gx*gy;i++)
		amat[i] = 0.0;
	for(i=0;i<gx*gy;i++)
	{
		for(j=0;j<gx*gy;j++)
		{
			if( on_boundary( pts + j * dim ) == 1 )
				amat[i*gx*gy+j] = gqw[i] * rkp_term_evaluate_node( &rk, i, j );
			else
				for(k=0;k<dim;k++)
					amat[i*gx*gy+j] += gqw[i] * rkp_wavelet_term_evaluate_node( &rk, i, ord[k], j );
		}
	}

	/* Checking */
	for(i=0;i<gx*gy;i++)
	{
		for(j=0;j<gx*gy;j++)
			fprintf( stderr, "%15.7f", amat[i*gx*gy+j] );
		fprintf( stderr, "\n" );
	}

	/* Build the righthand side vector from a combination of the boundary conditions and the forcing function of the system */
	for(i=0;i<gx*gy;i++)
	{
		if( on_boundary( pts + i * dim ) == 1 )
			rhs[i] = 0.0;
			//rhs[i] = 0.1 * pts[i*dim+0] * pts[i*dim+0];
		else
			rhs[i] = rhsf( 2, pts + i * dim );
	}

	/* Checking */
	for(i=0;i<gx*gy;i++)
		fprintf( stderr, "%15.7f\n", rhs[i] );

	/* Solve the system */
	i = gx * gy;
	j = 1;
	int ipiv[i],res;
	dgesv_( &i, &j, amat, &i, ipiv, rhs, &i, &res );
	fprintf( stderr, "res = %d\n", res );

	/* Copy the solution into the rk structure */
	for(i=0;i<gx*gy;i++)
		rk.vals[i] = rhs[i];

	/* Solution print */
	for(i=0;i<gx*gy;i++)
		fprintf( stderr, "val = %15.7f\n", rhs[i] );

	/* Output the function */
	FILE *fp = fopen( "plot", "w" );
	for(i=0;i<pgx*pgy;i++)
	{
		fprintf( fp, "%15.7f%15.7f%15.7f\n", ppts[i*dim+0], ppts[i*dim+1], rkp_evaluate( &rk, ppts + i * dim ) );
		if( ( i + 1 ) % pgy == 0 )
                        fprintf( fp, "\n" );
	}
	fclose( fp );

	free( grad );
	free( pts );
	free( dlt );
	free( val );
	free( gqw );
	free( ppts );

	return 0;
}

