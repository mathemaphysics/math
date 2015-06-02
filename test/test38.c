/**
 * This file is a basic template for use in writing any software which
 * loads all of its potentially complicated configuration information
 * from a single primary file. Other configuration routines may be
 * performed as a result of directives inside this primary configuration
 * file (which is loaded via a command line switch).
 */

/* Basic stuff we always use */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>

/* Includes for configuration */
#include "trees.h"
#include "parse.h"
#include "cfgread.h"

/* Includes for specific tasks for this template */
#include "points.h"
#include "punity.h"
#include "cgrad.h"
#include "nuclei.h"

/* Some variables defining data structure sizes */
#define TOKEN_BUFFER_LENGTH 1024
#define SPARSE_MAT_INCREMENT 1024
#define SYMTAB_INIT_SIZE 1024
#define NUM_EIG_DIMS_TO_PRINT 5
#define MAX_NUM_RANGES 512

int main( int argc, char **argv )
{
	/* Counting variables */
	int i,j,k,m,n,p,ret;
	double x,y,sum,tum;

	/* Some generic character buffers and file variables */
	FILE *fp;

	/* Variables for storing ranges */
	long nfr,rng[MAX_NUM_RANGES][2];

	/* Constants needed later */
	const long c_spmat_inc = SPARSE_MAT_INCREMENT;

	/* Point cloud variables */
	int dim,deg,rdn,npts,nbp;
	double *pts,*dlt;
	int *drv;
	punity_t pt;

	/* Quadrature variables */
	double *ctr,*qbox,*nqbox,*qx,qw,rad;
	long *size,*index,*tindex;
	int *ltg;
	shape_t dm;

	/* Quadrature itself */
	int quadn;
	double *qpts,*qwts;

	/* Sparse matrix variables themselves */
	long *ia,*ja,*ib,*jb,cc;
	double *A,*B,*alpha,*beta,*gamma,*omega,*ev,*eig;
	int mo,neig;

	/* Nuclei variables */
	nuclei_t nuc;

	/* Grid variables for plotting */
	int *dx;
	double *x0,*wx,*dd;

	/* Variables for the not so sparse stuff */
	double *Q;

	/* Variables for configuration */
	int optc;
	symtab_t sym;
	symbol_t *s;

	/* Load file containing variable declarations */
	#include "puksvar.h"

	/* Output who I am */
	fprintf( stderr, "\nPUKSHAM VERSION 0.1\n\n" );

	/* Do all the configuration operations here */
	#include "pukscfg.h"

	  //////////////////////////////////////////
         // Begin actual code here. Good luck... //
        //////////////////////////////////////////

	/* Allocate all space needed at least initially */
	if( !b_points_loaded )
	{
		fprintf( stderr, " * ERROR: No points loaded. Exiting.\n\n" );
		return 0;
	}
	if( !b_quad_loaded )
	{
		fprintf( stderr, " * ERROR: No quadrature loaded. Exiting.\n\n" );
		return 0;
	}
	if( !b_domain_set )
	{
		fprintf( stderr, " * ERROR: No domain indicator loaded. Exiting.\n\n" );
		return 0;
	}

	/* Check for discrepancies in Lanczos parameters */
	if( i_max_lanczos_steps > npts )
	{
		i_max_lanczos_steps = npts;
		fprintf( stderr, " ---> ERROR: Value of maxLanczosSteps > npts: Setting maxLanczosSteps to %d\n", npts );
	}

	/* Set up the partition of unity */
	dlt = (double*) malloc( npts * sizeof(double) );

	/* Setting up supports */
	fprintf( stderr, " * Generating supports\n\n" );
	generate_cloud_supports_min( dim, npts, pts, rdn, 1.05, dlt, 0.001 );

	/* Initializing Shepard partition of unity */
	fprintf( stderr, " * Initializing Shepard partition of unity " );
	punity_init( &pt, dim, npts, pts, dlt, &cubic_window, &cubic_window_deriv );
#ifdef PUNITY_USE_KDTREES
	generate_cloud_supports_min( dim, npts, pt.pts, rdn, 1.05, pt.dlt, 0.001 );
#endif
	fprintf( stderr, "with rmax = %15.7f\n\n", pt.rmax );

	/* Allocate stuff for integral evaluation */
	ctr = (double*) malloc( dim * sizeof(double) );
        qbox = (double*) malloc( dim * dim * sizeof(double) );
        nqbox = (double*) malloc( dim * dim * sizeof(double) );
        size = (long*) malloc( dim * sizeof(long) );
        index = (long*) malloc( dim * sizeof(long) );
	tindex = (long*) malloc( dim * sizeof(long) );
        qx = (double*) malloc( dim * sizeof(double) );

	/* Allocate space for the sparse matrices */
	ia = (long*) malloc( npts * sizeof(long) );
	ib = (long*) malloc( npts * sizeof(long) );
	ja = (long*) malloc( c_spmat_inc * sizeof(long) );
	jb = (long*) malloc( c_spmat_inc * sizeof(long) );
	A = (double*) malloc( c_spmat_inc * sizeof(double) );
	B = (double*) malloc( c_spmat_inc * sizeof(double) );
	cc = c_spmat_inc;

	/* Allocate space for Laplace operator */
	drv = (int*) malloc( dim * sizeof(int) );

	/* Deal with boundary functions */
	nbp = 0;
	for(i=0;i<npts;i++)
	{
		boundary_indicator( pt.pts + i * dim, &dm, f_bndry_tol, &ret );
		if( ret == 1 )
			++nbp, pt.bdry[i] = 1; /* Set this node to have a singularity to give delta property on boundary */
	}
	fprintf( stderr, " * Detected %d boundary points using given criteria\n\n", nbp );

	/* Now that boundary points are known, make final adjustments to Lanczos parameters */
	if( i_max_lanczos_steps > npts - nbp )
	{
		if( i_max_lanczos_steps > npts - nbp )
		{
			i_max_lanczos_steps = npts - nbp;
			fprintf( stderr, " * ERROR: Value of maxLanczosSteps > npts - nbp: Setting maxLanczosSteps to %d\n\n", npts - nbp );
		}
	}
	if( neig > i_max_lanczos_steps )
	{
		neig = i_max_lanczos_steps;
		fprintf( stderr, " * ERROR: Requested numberEigenvalues > maxLanczosSteps: Setting neig = %d\n\n", i_max_lanczos_steps );
	}

	/* Generate Q for saving Lanczos vectors later */
	Q = (double*) malloc( i_max_lanczos_steps * npts * sizeof(double) );

	/* Loading matrices */
	if( b_load_overl_mat )
	{
		fprintf( stderr, " * Loaded overlap matrix from \"%s\"\n\n", s_overl_fname );
	}
	if( b_load_stiff_mat )
	{
		fprintf( stderr, " * Loaded stiffness matrix from \"%s\"\n\n", s_stiff_fname );
	}

	/* Start building the matrices */
	if( !b_load_stiff_mat || !b_load_overl_mat )
	{
		fprintf( stderr, " * Building matrices " );
		/* Otherwise need to build system which does not contain the boundary nodes */
		ltg = (int*) malloc( ( npts - nbp ) * sizeof(int) );
		for(i=0,p=0;i<npts;i++)
			if( pt.bdry[i] == 0 )
				ltg[p++] = i;

		/* Initialize the matrices */
		for(i=0;i<pt.npts-nbp;i++)
			ia[i] = -1, ib[i] = -1;
		p = 0;
		for(i=0;i<npts-nbp;i++)
		{
			for(j=0;j<i/(npts/20);j++)
				fprintf( stderr, "=" );
			fprintf( stderr, ">%3d%%", i * 100 / ( npts - nbp ) );
			for(j=i;j<npts-nbp;j++)
			{
				if( sphere_intersection( dim, pt.pts + ltg[i] * dim, pt.dlt[ltg[i]], pt.pts + ltg[j] * dim, pt.dlt[ltg[j]], ctr, &rad, qbox ) == 1 )
				{
					/* Build the local basis for the intersection of spheres */
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
					sum = 0.0, tum = 0.0;
					do
                                	{
                                        	lens_gauss_point( dim, pt.pts + ltg[i] * dim, pt.dlt[ltg[i]], pt.pts + ltg[j] * dim, pt.dlt[ltg[j]],
                                                	rad, nqbox, index, qpts, qwts, qx, &qw );
						domain_indicator( qx, &dm, &ret );
						if( ret == 1 )
						{
							if( !b_load_overl_mat )
							{
                        	                		x = punity_evaluate_delta( &pt, ltg[i], qx, i_sing_order )
										* punity_evaluate_delta( &pt, ltg[j], qx, i_sing_order );
								sum += qw * x;
							}
							if( !b_load_stiff_mat )
							{
								y = 0.0;
								for(k=0;k<dim;k++)
								{
									for(m=0;m<dim;m++)
									{
										if( m == k )
											drv[m] = 1;
										else
											drv[m] = 0;
									}
									y += punity_term_deriv_evaluate_delta( &pt, ltg[i], drv,
													qx, i_sing_order, NULL, 0 )
									* punity_term_deriv_evaluate_delta( &pt, ltg[j], drv,
													qx, i_sing_order, NULL, 0 );
								}
								tum += qw * y;
							}
						}
					}
                                	while( arraynext( (long) dim, size, index ) != -1 );

					/* Now put sum in the i,j entry in the sparse matrx of inner product entries */
					if( p > cc )
					{
						cc += c_spmat_inc;
						if( !b_load_stiff_mat )
						{
							A = (double*) realloc( A, cc * sizeof(double) );
							ja = (long*) realloc( ja, cc * sizeof(long) );
						}
						if( !b_load_overl_mat )
						{
							B = (double*) realloc( B, cc * sizeof(double) );
							jb = (long*) realloc( jb, cc * sizeof(long) );
						}
					}
					if( !b_load_stiff_mat )
					{
						A[p] = tum;
						ja[p] = j;
						if( ia[i] == -1 )
                                                	ia[i] = p;
					}
					if( !b_load_overl_mat )
					{
						B[p] = sum;
						jb[p] = j;
						if( ib[i] == -1 )
                                                	ib[i] = p;
					}

					/* Step to the next entry */
					p++;
				}
			}
			parse_print_back( stderr, i / ( npts / 20 ) + 1 + 4 );
		}
			
		if( !b_load_stiff_mat )
			ia[npts-nbp] = p;
		if( !b_load_overl_mat )
			ib[npts-nbp] = p;
		fprintf( stderr, "\n\n" );

		/* Indicate that matrices are both built */
		if( !b_load_stiff_mat )
			b_have_stiff_mat = 1;
		if( !b_load_overl_mat )
			b_have_overl_mat = 1;
	}

	/* Don't go on if the matrices have not been loaded or built */
	if( !b_have_stiff_mat )
	{
		fprintf( stderr, " * ERROR: Stiffness matrix has not been built. This is really bad. Exiting angrily.\n\n" );
		return 0;
	}
	if( !b_have_overl_mat )
	{
		fprintf( stderr, " * ERROR: Overlap matrix has not been built. This is really bad. Exiting angrily.\n\n" );
		return 0;
	}

	/* Saving matrices to output files */
	if( b_save_stiff_mat || b_save_overl_mat )
		fprintf( stderr, " * Writing matrix output files\n\n" );
	if( b_save_overl_mat )
        {
		smsave( s_overl_fname, npts - nbp, npts - nbp, ib, jb, B, 1 );
                fprintf( stderr, " ---> Saved overlap matrix to output \"%s\"\n", s_overl_fname );
        }
	if( b_save_stiff_mat )
	{
		smsave( s_stiff_fname, npts - nbp, npts - nbp, ia, ja, A, 1 );
		fprintf( stderr, " ---> Saved stiffness matrix to output \"%s\"\n", s_stiff_fname );
	}
	if( b_save_overl_mat || b_save_stiff_mat )
		fprintf( stderr, "\n" );

	/* Solve the eigenvalue problem via Lanczos projection */
	if( !b_neig_set || ( b_neig_set && neig > npts ) )
	{
		fprintf( stderr, " * ERROR: Number of eigenvalues not specified or invalid. Exiting.\n\n " );
		return 0;
	}
	alpha = (double*) malloc( ( npts + 2 ) * sizeof(double) );
	beta = (double*) malloc( ( npts + 2 ) * sizeof(double) );
	gamma = (double*) malloc( ( npts + 2 ) * sizeof(double) );
	omega = (double*) malloc( ( npts + 2 ) * sizeof(double) );
	ev = (double*) malloc( 2 * npts * sizeof(double) );
	eig = (double*) malloc( neig * ( npts - nbp ) * sizeof(double) );

	/* Seed the random number generator */
	srand((unsigned)time(0));

	/* Lanczos procedure */
	fprintf( stderr, " * Beginning Lanczos procedure\n\n" );
	sgsilanczoscr( npts - nbp, i_max_lanczos_steps, neig, ib, jb, B, ia, ja, A, alpha, beta, omega, NULL, ev, Q, &mo, 1e-15, &x, 10000, 0, &ret );
	fprintf( stderr, " * Converged Lanczos procedure in %d steps with a collective residual of %5.9e over all EV's\n\n", mo, x );

	/* Now for sgsilanczos need to prepare alpha, beta by premultiplying by inverse of diag(omega) */
	copy( mo + 2, beta, gamma );
	for(i=0;i<mo-1;i++)
		alpha[i] /= omega[i], beta[i+1] /= omega[i], gamma[i+1] /= omega[i+1];
	alpha[mo-1] /= omega[mo-1];

	/* Solve the appropriate tridiagonal eigenvalue problem */
	if( strncmp( s_treig_routine, "lr", TOKEN_BUFFER_LENGTH ) == 0 )
		treiglr( mo, alpha, beta + 1, gamma + 1, 10000, ev, &ret, &n );
	else if( strncmp( s_treig_routine, "qds", TOKEN_BUFFER_LENGTH ) == 0 )
		treigqds( mo, alpha, beta + 1, gamma + 1, 10000, ev, &ret, &n );
	else if( strncmp( s_treig_routine, "dqds", TOKEN_BUFFER_LENGTH ) == 0 )
		treigdqds( mo, alpha, beta + 1, gamma + 1, 10000, ev, &ret, &n );
	else if( strncmp( s_treig_routine, "tridqds", TOKEN_BUFFER_LENGTH ) == 0 )
		treigtridqds( mo, alpha, beta + 1, gamma + 1, 10000, ev, &ret, &n );
	else
	{
		fprintf( stderr, " * ERROR: No valid tridiagonal eigenvalue solver selected: Given is \"%s\". Exiting.\n\n", s_treig_routine );
		return 0;
	}
	complex_bubble_sort( mo, ev );
	fprintf( stderr, " * Eigenvalues resulting from Lanczos process:\n\n" );
	for(i=0;i<neig;i++)
		fprintf( stderr, "%15.7f + %15.7fi\n", 1.0 / ev[2*i+0], ev[2*i+1] );
	fprintf( stderr, "\n" );

	/* Generate the eigenvectors */
	fprintf( stderr, " * Generating eigenvectors\n\n" );

	/* IMPORTANT: The shift mu used here should be an eigenvalue of the original system, not its inverse */
	fprintf( stderr, " ---> Vectors:\n" );
	for(i=0;i<neig;i++)
	{
		nrandv( npts - nbp, eig + i * ( npts - nbp ) );
		sgsinvi( npts - nbp, ia, ja, A, ib, jb, B, eig + i * ( npts - nbp ), 1.0 / ev[2*i+0], 1e-9, 10000, 0 );
		fprintf( stderr, " ----> %5d EV = %15.7f: ", i, 1.0 / ev[2*i+0] );
		for(j=0;j<NUM_EIG_DIMS_TO_PRINT;j++)
			fprintf( stderr, "%15.7f", eig[i*(npts-nbp)+j] );
		fprintf( stderr, "  ...\n" );
	}
	fprintf( stderr, "\n" );

	  /////////////////////////////////
	 // Normalize the wavefunctions //
	/////////////////////////////////

	//fprintf( stderr, " * Normalizing output wave functions\n\n" );
	//for(i=0;i<dim;i++)
	//	for(j=0;j<dim;j++)
	//		if( j == i )
	//			nqbox[i*dim+j] = 1.0;
	//		else
	//			nqbox[i*dim+j] = 0.0;
	//for(i=0;i<neig;i++)
	//{
	//	sum = 0.0;
	//	for(j=0;j<npts;j++)
	//	{
	//		for(k=0;k<dim;k++)
	//			index[k] = 0, size[k] = (long) quadn - 1;
	//		do
	//		{
	//			/* Calculate the spherical gauss point for this index */
	//			copy( dim, pt.pts + j * dim, ctr );
	//			rad = pt.dlt[j];
	//			sphere_gauss_point( dim, ctr, rad, nqbox, index, qpts, qwts, qx, &qw );
	//			domain_indicator( qx, &dm, &ret );
	//			if( ret == 1 )
	//				sum += qw * punity_evaluate( &pt, j, qx );
	//		}
	//		while( arraynext( (long) dim, size, index ) != -1 );
	//	}
	//}

	/* Do some plotting */
	if( b_out_eig_range_set )
	{
		/* Announce */
		fprintf( stderr, " * Plotting eigenvectors " );
		for(i=0;i<nfr;i++)
		{
			if( rng[i][0] == rng[i][1] )
				fprintf( stderr, "%d ", rng[i][0] );
			else if( rng[i][1] > rng[i][0] )
				fprintf( stderr, "%d-%d ", rng[i][0], rng[i][1] );
		}
		fprintf( stderr, "\n\n" );

		for(i=0;i<dim;i++)
			size[i] = (long) dx[i] - 1;
		for(i=0;i<dim;i++)
			index[i] = 0;
		do
		{
			/* Build the point at which to evaluate the eigenfunction */
			for(i=0;i<dim;i++)
				qx[i] = x0[i] + (double) index[i] * dd[i];
	
			/* Output the point for this line */
			for(i=0;i<dim;i++)
				fprintf( stdout, "%15.7f", qx[i] );
	
			/* Now evaluate the function at that point */
			for(n=0;n<nfr;n++)
			{
				for(m=rng[n][0];m<=rng[n][1];m++)
				{
					sum = 0.0;
					for(i=0;i<npts-nbp;i++)
					{
						y = 0.0;
						for(j=0;j<dim;j++)
							y += pow( qx[j] - pt.pts[ltg[i]*dim+j], 2.0 );
						y = sqrt( y );
						if( y < pt.dlt[ltg[i]] )
							sum += eig[m*(npts-nbp)+i] * punity_evaluate_delta( &pt, ltg[i], qx, i_sing_order );
					}

					/* Ouptut the values consecutively */
					fprintf( stdout, "%15.7f", sum );
				}
			}

			/* Close the line */
			fprintf( stdout, "\n" );

			/* Check to see if another newline is needed to maintain gnuplot's preferred format */
			for(i=0;i<dim;i++)
				tindex[i] = index[i];
			arraynext( (long) dim, size, tindex );
			for(i=0;i<dim;i++)
			{
				if( index[i] == dx[i] - 1 && tindex[i] != dx[i] - 1 )
				{
					fprintf( stdout, "\n" );
					break;
				}
			}
		}
		while( arraynext( (long) dim, size, index ) != -1 );
	}

	/* Announce exit */
	fprintf( stderr, " * Exiting happily. Have a nice day!\n\n" );

	/* Don't forget to clean up */
	free( pts ); free( dlt ); free( drv );
	free( qpts ); free( qwts );
	free( ia ); free( ib );
	free( ja ); free( jb );
	free( A ); free( B ); free( Q );
	free( alpha ); free( beta ); free( gamma ); free( omega );
	free( ctr ); free( qbox ); free( nqbox ); free( qx );
	free( size ); free( index );
	free( ev ); free( eig );
	free( ltg );

	return 0;
}

