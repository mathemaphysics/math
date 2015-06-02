#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include "points.h"
#include "punity.h"

#define SPARSE_MATRIX_SIZE_INC 1024

int on_boundary( double *x_in )
{
	
}

int indicator_square( double *x_in )
{
	if( x_in[0] > 0.0 && x_in[0] < 5.0 && x_in[1] > 0.0 && x_in[1] < 5.0 )
		return 1;
	else
		return 0;
}

int main( int argc, char **argv )
{
	/* Basic variables */
	int i,j,k,m,p,idx,dim,deg,npts,rdn;
	long *size,*index;
	double *pts,*dlt,*ctr,*qbox,*nqbox,*qx,qw,rad,x;
	punity_t pt;

	/* Switches */
	int b_grid = 0;

	/* Tensor product Gauss quadrature variables */
	int quadn = 6;
	double qpts[6] = { -0.93246951, -0.66120939, -0.23861919, 0.23861919, 0.66120939, 0.93246951 };
	double qwts[6] = {  0.17132449, 0.360767157, 0.46791393, 0.46791393, 0.36076157, 0.17132449 };

	/* Plot output variables */
	int gx = 100;
	int gy = 100;
	double lx = 3.0;
	double ly = 3.0;
	double dx = lx / (double) ( gx - 1 );
	double dy = ly / (double) ( gy - 1 );
	double sum,vec[dim];

	/* Take in the options */
        int optc;
        while( ( optc = getopt( argc, argv, "d:g:r:" ) ) != -1 )
        {
                switch( optc )
                {
                        case 'r':
                                rdn = atoi( optarg );
                                break;
                        case 'd':
                                deg = atoi( optarg );
                                break;
                        case 'g':
				load_points( optarg, &dim, &npts, &pts, 0 );
				dlt = (double*) malloc( npts * sizeof(double) );
				b_grid = 1;
				fprintf( stderr, "Loaded %d points in %d dimensions\n", npts, dim );
                                break;
                        case '?':
                        default:
                                fprintf( stderr, "There was an error or something. Have fun fixing it.\n" );
                                return 0;
                                break;
                }
        }

	/* Make sure a cloud was loaded */
	if( !b_grid )
	{
		fprintf( stderr, "Please select a cloud file. Exiting.\n" );
		return 0;
	}

	/* Test the partitioning functions */
	int drv[2] = { 1, 1 };
	double dr[2] = { 0.08, 0.15 };
	double ds[2] = { 1.55, 0.41 };

	/* Generate supports for each point in the cloud */
	generate_cloud_supports_min( dim, npts, pts, rdn, 1.05, dlt, 0.001 );
        punity_init( &pt, dim, npts, pts, dlt, &cubic_window, &cubic_window_deriv );
	generate_cloud_supports_min( dim, npts, pt.pts, rdn, 1.05, pt.dlt, 0.001 );
	fprintf( stderr, "rmax = %15.7f\n", pt.rmax );

	/* Allocate stuff for sphere intersection */
	ctr = (double*) malloc( dim * sizeof(double) );
	qbox = (double*) malloc( dim * dim * sizeof(double) );
	nqbox = (double*) malloc( dim * dim * sizeof(double) );
	size = (long*) malloc( dim * sizeof(long) );
	index = (long*) malloc( dim * sizeof(long) );
	qx = (double*) malloc( dim * sizeof(double) );

	/* Allocate space for sparse matrix entries */
	int salloc = SPARSE_MATRIX_SIZE_INC;
	double *amat = (double*) malloc( salloc * sizeof(double) );
	int *iamat = (int*) malloc( pt.npts * sizeof(int) );
	int *jamat = (int*) malloc( salloc * sizeof(int) );
	double *bmat = (double*) malloc( salloc * sizeof(double) );
	int *ibmat = (int*) malloc( pt.npts * sizeof(int) );
	int *jbmat = (int*) malloc( salloc * sizeof(int) );

	/* Calculate set boundary set */
	int nb = 0;
	for(i=0;i<pt.npts;i++)
		if( on_boundary( pt.pts + i * dim ) == 1 )
			++nb;
	double *v = (double*) malloc( nb * pt.npts * sizeof(double) );
	double *w = (double*) malloc( ( pt.npts - nb ) * sizeof(double) );

	/* Try a spherical integral to test */
        for(j=0;j<dim;j++)
                for(k=0;k<dim;k++)
                        if( j != k ) /* Alway have axis-aligned grid! */
                                nqbox[j*dim+k] = 0.0;
                        else
                                nqbox[j*dim+k] = 0.5;
	for(k=0;k<dim;k++)
        	index[k] = 0, size[k] = quadn - 1;
	ctr[0] = 0.0, ctr[1] = 0.0;
	sum = 0.0;
	do
        {
                sphere_gauss_point( dim, ctr, 0.5, nqbox, index, qpts, qwts, qx, &qw );
		sum += 1.0 * qw;
		fprintf( stderr, "qw = %15.7f\n", qw );
        }
        while( arraynext( (long) dim, size, index ) != -1 );
	fprintf( stderr, "sum = %15.7f\n", sum );

	/* Stop here for now for testing */
	//return 0;

	/* Initialize sparse structures */
	for(i=0;i<pt.npts;i++)
		iamat[i] = -1, ibmat[i] = -1;

	/* Iterate through pairs */
	p = 0; /* Start pointer to next matrix entry to zero in amat and bmat */
	for(i=0;i<npts;i++)
	{
		for(j=0;j<=i;j++)
		{
			if( sphere_intersection( dim, pt.pts + i * dim, pt.dlt[i], pt.pts + j * dim, pt.dlt[j], ctr, &rad, qbox ) == 1 )
			{
				fprintf( stderr, "%4d%4d", i, j );
				for(k=0;k<dim;k++)
				{
					sum = 0.0;
					for(m=0;m<dim;m++)
						sum += qbox[k*dim+m] * qbox[k*dim+m];
					sum = sqrt( sum );
					for(m=0;m<dim;m++)
						nqbox[k*dim+m] = qbox[k*dim+m] / sum;
				}

				/* Calculate the entries of the inner product matrix */
				for(k=0;k<dim;k++)
					index[k] = 0, size[k] = quadn - 1;
				sum = 0.0;
				do
                                {
                                        lens_gauss_point( dim, pt.pts + i * dim, pt.dlt[i], pt.pts + j * dim, pt.dlt[j],
                                                rad, nqbox, index, qpts, qwts, qx, &qw );
					if( indicator_square( qx ) == 1 )
					{
                                        	x = punity_evaluate( &pt, i, qx ) * punity_evaluate( &pt, j, qx );
						sum += qw * x;
					}
                                }
                                while( arraynext( (long) dim, size, index ) != -1 );
				fprintf( stderr, ": %15.12f", sum );

				/* Now put sum in the i,j entry in the sparse matrx of inner product entries */
				if( p > salloc )
				{
					salloc += SPARSE_MATRIX_SIZE_INC;
					amat = (double*) realloc( amat, salloc * sizeof(double) );
					jamat = (int*) realloc( jamat, salloc * sizeof(int) );
					bmat = (double*) realloc( bmat, salloc * sizeof(double) );
					jbmat = (int*) realloc( jbmat, salloc * sizeof(int) );
				}
				amat[p] = sum;
				jamat[p] = j;
				if( iamat[i] == -1 )
					iamat[i] = p;

				/* Calculate entries of the stiffness matrix */
				for(k=0;k<dim;k++)
					index[k] = 0;
				sum = 0.0;
				do
				{
					lens_gauss_point( dim, pt.pts + i * dim, pt.dlt[i], pt.pts + j * dim, pt.dlt[j],
						rad, nqbox, index, qpts, qwts, qx, &qw );
					for(k=0;k<dim;k++)
					{
						for(m=0;m<dim;m++)
							if( m == k )
								drv[m] = 1;
							else
								drv[m] = 0;
						if( indicator_square( qx ) == 1 )
						{
							x = punity_term_deriv_evaluate( &pt, i, drv, qx, NULL, 0 ) * punity_term_deriv_evaluate( &pt, j, drv, qx, NULL, 0 );
							sum += qw * x;
						}
					}
				}
				while( arraynext( (long) dim, size, index ) != -1 );
				fprintf( stderr, " %15.12f\n", sum );

				/* Put sum in the i,j entry of the sparse matrix of stiffness matrix */
				bmat[p] = sum;
				jbmat[p] = j;
				if( ibmat[i] == -1 )
					ibmat[i] = p;

				/* Increment sparse matrix pointer */
				++p;
			}
		}
	}
	iamat[pt.npts] = p;
	ibmat[pt.npts] = p;
	fprintf( stderr, "\n" );

	/* Generate the boundary projector */
	for(i=0,p=0;i<pt.npts;i++)
	{
		if( on_boundary( pt.pts + i * dim ) == 1 )
		{
			for(j=0;j<dim;j++)
				v[p*dim+j] = sparse_entry( i, j, pt.npts, pt.npts, iamat, jamat, amat );
			
		}
	}

	/* Start the solver routine */
	

	free( pts );
	free( dlt );

	return 0;
}

