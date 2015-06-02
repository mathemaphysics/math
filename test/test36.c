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

int main( int argc, char **argv )
{
	/* Counting variables */
	int i,j,k,m,n,p,ret;
	double x,y,sum,tum;

	/* Some generic character buffers and file variables */
	char buf[1024],bug[1024],*tok[1024],*tol[1024];
	FILE *fp;

	/* Constants needed later */
	const long c_spmat_inc = 1024;

	/* Point cloud variables */
	int dim,deg,rdn,npts,nbp;
	double *pts,*dlt;
	int *drv;
	punity_t pt;

	/* Quadrature variables */
	double *ctr,*qbox,*nqbox,*qx,qw,rad;
	long *size,*index;
	shape_t dm;

	/* Quadrature itself */
	int quadn;
	double *qpts,*qwts;

	/* Sparse matrix variables themselves */
	long *ia,*ja,*ib,*jb,cc;
	double *A,*B,*alpha,*beta,*gamma,*ev;
	int neig;

	/* Variables for the not so sparse stuff */
	double *V;

	/* Variables for configuration */
	int optc;
	symtab_t sym;
	symbol_t *s;

	/* Configuration switches */
	int b_points_loaded = 0;
	int b_config_loaded = 0;
	int b_quad_loaded = 0;
	int b_rdn_set = 0;
	int b_domain_set = 0;
	int b_use_periodic_boundary = 0;
	int b_use_basis_zero_on_boundary = 0;
	int b_neig_set = 0;
	int b_load_stiff_mat = 0;
	int b_load_overl_mat = 0;
	int b_save_stiff_mat = 0;
	int b_save_overl_mat = 0;

	/* String variables for configuration */
	char s_overl_fname[1024] = "amat.out";
	char s_stiff_fname[1024] = "bmat.out";

	/* Some floating point variables needed */
	double f_bndry_tol = 0.05;

	/* Output who I am */
	fprintf( stderr, "\nPUKSHAM VERSION 0.1\n\n" );

	/* Read command line options */
	while( ( optc = getopt( argc, argv, "C:p:q:r:e:" ) ) != -1 )
        {
                switch( optc )
                {
			case 'C':
				/* Initialize and load the symbol table */
				symtab_init( &sym, 1024 );
				cfgread_load_symbols_f( optarg, &sym );
				b_config_loaded = 1;

				/* Print symbol table quickly */
				fprintf( stderr, " * Loaded following symbols from \"%s\":\n\n", optarg );
				symtab_print( &sym );
				fprintf( stderr, "\n" );
				break;
			case 'p':
				/* Load points, number, and dimension */
				load_points( optarg, &dim, &npts, &pts, 0 );
				b_points_loaded = 1;

				/* Output notification */
				fprintf( stderr, " * Loaded %d points in %d dimensions from \"%s\"\n\n", npts, dim, optarg );
				break;
			case 'q':
				/* Load quadrature points from file */
				load_quadrature( optarg, &quadn, &qpts, &qwts, &ret );

				/* Make sure nothing went wrong and output what happens */
				if( ret == 0 )
				{
					b_quad_loaded = 1; /* Only mark loaded if there is no trouble */
					fprintf( stderr, " * Loaded quadrature with %d points from \"%s\"\n\n", quadn, optarg );
				}
				else
					fprintf( stderr, " * ERROR: Problem loading quadrature rule from \"%s\": Return value %d\n\n", optarg, ret );
				fprintf( stderr, "\n" );
				break;
			case 'r':
				/* Load radial neighbors number */
				rdn = atoi( optarg );
				b_rdn_set = 1;
				fprintf( stderr, " * Setting radial neighbor number to %d\n\n", rdn );
				break;
			case 'e':
				/* Load the number of eigenvalues to use */
				neig = atoi( optarg );
				if( neig > 0 )
				{
					b_neig_set = 1;
					fprintf( stderr, " * Setting the number of eigenvalues to %d\n\n", neig );
				}
				else
					fprintf( stderr, " * ERROR: Invalid number of eigenvalues: %d\n\n", neig );
				break;
			case '?':
			default:
				fprintf( stderr, "I don't understand the jibberish you are spouting.\n" );
				return 0;
		}
	}

	/* Process the configuration variables */
	if( b_config_loaded )
	{
		/* Print that searching for variables in configuration file */
		fprintf( stderr, " * Reading and executing configuration variables:\n\n" );

		/* Look for points file first if not loaded already */
		if( !b_points_loaded )
		{
			s = symtab_lookup( &sym, "pointCloud" );
			if( s != NULL )
			{
				load_points( (char*) s->data, &dim, &npts, &pts, 0 );
				fprintf( stderr, " ---> Loaded %d points in %d dimensions from \"%s\"\n", npts, dim, (char*) s->data );
				b_points_loaded = 1;
			}
		}
		if( !b_quad_loaded )
		{
			s = symtab_lookup( &sym, "quadrature" );
			if( s != NULL )
			{
				load_quadrature( (char*) s->data, &quadn, &qpts, &qwts, &ret );
				if( ret == 0 )
				{
					fprintf( stderr, " ---> Loaded quadrature rule with %d points from \"%s\"\n", quadn, (char*) s->data );
					b_quad_loaded = 1;
				}
				else
					fprintf( stderr, " ---> ERROR: Problem loading quadurature file in \"%s\": Return value %d\n", (char*) s->data, ret );
			}
		}
		if( !b_rdn_set )
		{
			s = symtab_lookup( &sym, "radialNeighbors" );
			if( s != NULL )
			{
				rdn = atoi( (char*) s->data );
				b_rdn_set = 1;
				fprintf( stderr, " ---> Loaded radial neighbors number %d\n", rdn );
			}
		}
		if( !b_domain_set )
		{
			s = symtab_lookup( &sym, "domain" );
			if( s != NULL )
			{
				load_domain( (char*) s->data, &dm, &ret );
				if( ret == 0 )
				{
					fprintf( stderr, " ---> Loaded domain information\n" );
					b_domain_set = 1;
				}
			}
		}
		s = symtab_lookup( &sym, "boundaryTolerance" );
		if( s != NULL )
		{
			f_bndry_tol = atof( (char*) s->data );
			fprintf( stderr, " ---> Loaded boundary tolerance of %15.7f\n", f_bndry_tol );
		}
		if( !b_neig_set )
		{
			s = symtab_lookup( &sym, "numberEigenvalues" );
			if( s != NULL )
			{
				neig = atoi( (char*) s->data );
				if( neig > 0 )
				{
					b_neig_set = 1;
					fprintf( stderr, " ---> Loaded number of eigenvalues to %d\n", neig );
				}
				else
					fprintf( stderr, " ---> ERROR: Invalid number of eigenvalues: %d\n", neig );
			}
		}
		s = symtab_lookup( &sym, "overlFileName" );
                if( s != NULL )
                {
                        strncpy( s_overl_fname, (char*) s->data, 1024 );
                        fprintf( stderr, " ---> File name of overlap matrix set to \"%s\"\n", s_overl_fname );
                }
                s = symtab_lookup( &sym, "stiffFileName" );
                if( s != NULL )
                {
                        strncpy( s_stiff_fname, (char*) s->data, 1024 );
                        fprintf( stderr, " ---> File name of stiffness matrix set to \"%s\"\n", s_stiff_fname );
                }
		s = symtab_lookup( &sym, "loadOverlMat" );
                if( s != NULL )
                {
                        if( strcmp( (char*) s->data, "yes" ) == 0 )
                        {
                                b_load_overl_mat = 1;
                                fprintf( stderr, " ---> Loading overlap matrix from file \"%s\"\n", s_overl_fname );
                        }
                }
		s = symtab_lookup( &sym, "loadStiffMat" );
		if( s != NULL )
		{
			if( strcmp( (char*) s->data, "yes" ) == 0 )
			{
				b_load_stiff_mat = 1;
				fprintf( stderr, " ---> Loading stiffness matrix from file \"%s\"\n", s_stiff_fname );
			}
		}
		s = symtab_lookup( &sym, "saveOverlMat" );
                if( s != NULL )
                {
                        if( strcmp( (char*) s->data, "yes" ) == 0 )
                        {
                                b_save_overl_mat = 1;
                                fprintf( stderr, " ---> Saving overlap matrix to file \"%s\"\n", s_overl_fname );
                        }
                }
                s = symtab_lookup( &sym, "saveStiffMat" );
                if( s != NULL )
                {
                        if( strcmp( (char*) s->data, "yes" ) == 0 )
                        {
                                b_save_stiff_mat = 1;
                                fprintf( stderr, " ---> Saving stiffness matrix to file \"%s\"\n", s_stiff_fname );
                        }
                }
		fprintf( stderr, "\n" );
	}

	/* Allocate all space needed at least initially */
	if( !b_points_loaded )
	{
		fprintf( stderr, "ERROR: No points loaded. Exiting.\n\n" );
		return 0;
	}
	if( !b_quad_loaded )
	{
		fprintf( stderr, "ERROR: No quadrature loaded. Exiting.\n\n" );
		return 0;
	}
	if( !b_domain_set )
	{
		fprintf( stderr, "ERROR: No domain indicator loaded. Exiting.\n\n" );
		return 0;
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
			++nbp;
	}
	V = (double*) malloc( nbp * npts * sizeof(double) );
	fprintf( stderr, " * Setting up constraint matrix for boundary functions with %d boundary points\n\n", nbp );
	for(p=0,i=0;i<npts;i++)
	{
		/* Make sure to call this thing with a tolerance, dumbass */
		boundary_indicator( pt.pts + i * dim, &dm, f_bndry_tol, &ret );
		if( ret == 0 )
			continue;

		/* Otherwise build the vector of all basis functions evaluated at each boundary point */
		for(j=0;j<npts;j++)
			V[p*npts+j] = punity_evaluate( &pt, j, pt.pts + i * dim );

		/* Project stuff */
		project( npts, V + p * npts, p, V );
		normalize( npts, V + p * npts );

		/* Advance to next vector */
		++p;
	}

	/* Temporarily output basis for checking */
	double *W = (double*) malloc( ( npts - nbp ) * npts * sizeof(double) );
	for(i=0;i<npts-nbp;i++)
	{
		nrandv( npts, W + i * npts );
		project( nbp, W + i * npts, nbp, V );
		project( npts, W + i * npts, i, W );
		normalize( npts, W + i * npts );
	}
	msave( "wmat.out", npts - nbp, npts, W, 1 );

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
		for(i=0;i<pt.npts;i++)
			ia[i] = -1, ib[i] = -1;
		p = 0;
		for(i=0;i<npts;i++)
		{
			for(j=0;j<i/(npts/20);j++)
				fprintf( stderr, "=" );
			fprintf( stderr, ">%3d%%", i * 100 / npts );
        	        for(j=i;j<npts;j++)
                	{
				if( sphere_intersection( dim, pt.pts + i * dim, pt.dlt[i], pt.pts + j * dim, pt.dlt[j], ctr, &rad, qbox ) == 1 )
	                        {
					/* If they intersect then output the two circles and then the quadrature points */
					  /////////////////////////////////
					 /* FIXME: Get rid of this ASAP */
					/////////////////////////////////
					if( i != j && i < 15 && j < 15 )
					{
						double xy[2];
						FILE *ff = fopen( "circles", "a" );
						for(k=0;k<100;k++)
						{
							xy[0] = pt.pts[i*dim+0] + cos( 2 * M_PI / (double) 100 * (double) k ) * pt.dlt[i];
							xy[1] = pt.pts[i*dim+1] + sin( 2 * M_PI / (double) 100 * (double) k ) * pt.dlt[i];
							fprintf( ff, "%15.7f%15.7f\n", xy[0], xy[1] );
							xy[0] = pt.pts[j*dim+0] + cos( 2 * M_PI / (double) 100 * (double) k ) * pt.dlt[j];
							xy[1] = pt.pts[j*dim+1] + sin( 2 * M_PI / (double) 100 * (double) k ) * pt.dlt[j];
							fprintf( ff, "%15.7f%15.7f\n", xy[0], xy[1] );
						}
						fclose( ff );
					}

					/* Build the normalized local basis */
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
                                        	lens_gauss_point( dim, pt.pts + i * dim, pt.dlt[i], pt.pts + j * dim, pt.dlt[j],
                                                	rad, nqbox, index, qpts, qwts, qx, &qw );
						domain_indicator( qx, &dm, &ret );
						if( ret == 1 )
						{
							  ////////////////////////////////////////////////////
							 /* FIXME: This is temporary; get rid of this ASAP */
							////////////////////////////////////////////////////
							if( i != j && i < 15 && j < 15 )
							{
								FILE *ff = fopen( "circles", "a" );
								fprintf( ff, "%15.7f%15.7f\n", qx[0], qx[1] );
								fclose( ff );
							}
							if( !b_load_overl_mat )
							{
                        	                		x = punity_evaluate( &pt, i, qx ) * punity_evaluate( &pt, j, qx );
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
									y += punity_term_deriv_evaluate( &pt, i, drv, qx, NULL, 0 ) * punity_term_deriv_evaluate( &pt, j, drv, qx, NULL, 0 );
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
			ia[npts] = p;
		if( !b_load_overl_mat )
			ib[npts] = p;
		fprintf( stderr, "\n\n" );
	}

	 ////////////////////////////////////////////////////////////////////
	/* FIXME: Temporarily exit for testing; get rid of this when done */
	return 0; /////////////////////////////////////////////////////////

	/* Saving matrices to output files */
	if( b_save_stiff_mat || b_save_overl_mat )
		fprintf( stderr, " * Writing matrix output files \"amat.oct\" and \"bmat.oct\"\n\n" );
	if( b_save_overl_mat )
        {
                smsave( s_overl_fname, npts, npts, ib, jb, B, 1 );
                fprintf( stderr, " ---> Saved overlap matrix to output \"%s\"\n", s_overl_fname );
        }
	if( b_save_stiff_mat )
	{
		smsave( s_stiff_fname, npts, npts, ia, ja, A, 1 );
		fprintf( stderr, " ---> Saved stiffness matrix to output \"%s\"\n", s_stiff_fname );
	}
	if( b_save_overl_mat || b_save_stiff_mat )
		fprintf( stderr, "\n" );

	/* Solve the eigenvalue problem via Lanczos projection */
	if( !b_neig_set || ( b_neig_set && neig > npts ) )
	{
		fprintf( stderr, "ERROR: Number of eigenvalues not specified or invalid. Exiting.\n\n " );
		return 0;
	}
	alpha = (double*) malloc( ( neig + 2 ) * sizeof(double) );
	beta = (double*) malloc( ( neig + 2 ) * sizeof(double) );
	gamma = (double*) malloc( ( neig + 2 ) * sizeof(double) );
	ev = (double*) malloc( 2 * neig * sizeof(double) );
	srand((unsigned)time(0));
	//sslanczos( npts, neig, ia, ja, A, alpha, beta, NULL, 1e-9, 0, &ret );
	sgsilanczos( npts, neig, ia, ja, A, ib, jb, B, alpha, beta, gamma, NULL, 1e-9, 10000, 0, &ret );
	treigtridqds( neig, alpha, beta + 1, beta + 1, 10000, ev, &ret, &n );
	complex_bubble_sort( neig, ev );
	fprintf( stderr, " * Eigenvalues resulting from Lanczos process:\n\n" );
	for(i=0;i<neig;i++)
		fprintf( stderr, "%15.7f + %15.7fi\n", ev[2*i+0], ev[2*i+1] );
	fprintf( stderr, "\n" );

	/* Don't forget to clean up */
	free( pts ); free( dlt ); free( drv );
	free( qpts ); free( qwts );
	free( ia ); free( ib );
	free( ja ); free( jb );
	free( A ); free( B );
	free( ctr ); free( qbox ); free( nqbox ); free( qx );
	free( size ); free( index );

	return 0;
}

