#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "points.h"
#include "punity.h"
#include "nuclei.h"

/**
 * Generate the overlap and Hamiltonian matrices using the given PU
 * and the polynomial basis for each element in lbase
 * @param pu PU basis
 * @param mdim Matrix dimension
 * @param nbp Number of boundary points
 * @param dm Domain shape
 * @param nlbase Number of basis functions in each local element
 * @param lbase Indexes of polynomials to use in each local basis
 * @param ltg Local to global map
 * @param nuc Nuclei list for the external potential
 * @param quadn Number of points in the one-dimensional quadrature rule
 * @param qpts Quadrature points
 * @param qwts Quadrature weights
 * @param iap Row start positions for A
 * @param jap Column index list for entries of A
 * @param Ap Entries of matrix A
 * @param ibp Row start positions for B
 * @param jbp Column index list for entries of B
 * @param Bp Entries of matrix B
 * @param cc Current number of doubles allocated in A/B and ja/jb
 * @param c_spmat_inc Increment size for growing A/B and ja/jb
 * @param b_use_external_potential Set 1 to use external potential, 0 to not
 * @param b_load_overl_mat Set 1 to not build overlap matrix, 0 to build
 * @param b_load_stiff_mat Set 1 to not build stiffness matrix, 0 to build
 * @param b_use_singular Set 1 to use singular particles on the boundary, 0 otherwise
 * @param i_sing_order Order of singularity to use if b_use_singular is 1
 * @param b_have_stiff_mat Returns 1 if successfully built stiffness matrix, 0 otherwise
 * @param b_have_overl_mat Returns 1 if successfully built overlap matrix, 0 otherwise
 */
void puksham_mbuild( punity_t *pu, int mdim, int nbp, shape_t *dm, int *nlbase, int **lbase, wfsd_t *rbfd, int *ltg, nuclei_t *nuc,
			int quadn, double *qpts, double *qwts,
			long **iap, long **jap, double **Ap, long **ibp, long **jbp, double **Bp, long cc, long c_spmat_inc,
				int b_use_external_potential, int b_load_overl_mat, int b_load_stiff_mat,
					int b_use_singular, int i_sing_order, int *b_have_stiff_mat, int *b_have_overl_mat )
{
	/* Counting variables */
	int i,ii,j,jj,k,m,p,pp,q,qq,ret;

	/* Important stuff */
	int dim = pu->dim;
	int npts = pu->npts;
	int *nb,nnb;
	int *expl,*expr,*drv;
	int idx;
	long *ia = *iap, *ja = *jap;
	long *ib = *ibp, *jb = *jbp;
	long *size,*index;
	double *A = *Ap, *B = *Bp;
	double *ctr;
	double x1,x2,y,rad;
	double *qbox,*nqbox;
	double qw,*qx,*tx;
	double rum,sum,tum;

	/* Allocate stuff */
	qbox = (double*) malloc( dim * dim * sizeof(double) );
	nqbox = (double*) malloc( dim * dim * sizeof(double) );
	ctr = (double*) malloc( dim * sizeof(double) );
	qx = (double*) malloc( dim * sizeof(double) );
	tx = (double*) malloc( dim * sizeof(double) );
	size = (long*) malloc( dim * sizeof(long) );
	index = (long*) malloc( dim * sizeof(long) );
	expl = (int*) malloc( dim * sizeof(int) );
	expr = (int*) malloc( dim * sizeof(int) );
	drv = (int*) malloc( dim * sizeof(int) );

	/* Initialize the matrices */
	for(i=0;i<mdim;i++)
		ia[i] = -1, ib[i] = -1;
	p = 0, pp = 0;
	for(i=0;i<npts-nbp;i++)
	{
		if( npts >= 20 )
		{
			for(j=0;j<i/(npts/20);j++)
				fprintf( stderr, "=" );
			fprintf( stderr, ">%3d%%", i * 100 / ( npts - nbp ) );
		}
		for(ii=0;ii<nlbase[ltg[i]];ii++)
		{
			/* Mark the start of the diagonal entry in the global scheme */
			qq = pp;

			/* Iterate through second function to calculate inner product */
			for(j=i;j<npts-nbp;j++)
			{
				for(jj=(j==i?ii:0);jj<nlbase[ltg[j]];jj++)
				{
					if( sphere_intersection( dim, pu->pts + ltg[i] * dim, pu->dlt[ltg[i]], pu->pts + ltg[j] * dim, pu->dlt[ltg[j]], ctr, &rad, qbox ) == 1 )
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
						rum = 0.0, sum = 0.0, tum = 0.0;
						do
						{
							lens_gauss_point( dim, pu->pts + ltg[i] * dim, pu->dlt[ltg[i]], pu->pts + ltg[j] * dim, pu->dlt[ltg[j]],
								rad, nqbox, index, qpts, qwts, qx, &qw );
							domain_indicator( qx, dm, &ret );
							if( ret == 1 )
							{
								/* Generate the neighbors list for point qx */
								nnb = punity_neighbors( pu, qx, &nb );

								/* Use nb to form these values */
								if( !b_load_overl_mat || !b_load_stiff_mat )
								{
									if( lbase[ltg[i]][ii] >= 0 )
									{
										global_polynomial_vector( lbase[ltg[i]][ii], dim, expl );
										if( b_use_singular )
											x1 = punity_eval_local_poly_delta( pu, ltg[i], expl, qx, i_sing_order );
										else
											x1 = punity_eval_local_poly( pu, ltg[i], expl, qx );
									}
									else
									{
										/* Deal with case in which lbase[ltg[i]][ii] < 0 */
										idx = -1 - lbase[ltg[i]][ii];

										/* Translate the vector to origin of particle ltg[i] */
										for(k=0;k<dim;k++)
											tx[k] = qx[k] - pu->pts[ltg[i]*dim+k];

										/* FIXME: Pick up here... */
										x1 = rbfd[idx]( dim, 0, pu->dlt[ltg[i]], tx ) * punity_evaluate( pu, ltg[i], qx );
									}
									if( lbase[ltg[j]][jj] >= 0 )
									{
										global_polynomial_vector( lbase[ltg[j]][jj], dim, expr );
										if( b_use_singular )
											x2 = punity_eval_local_poly_delta( pu, ltg[j], expr, qx, i_sing_order );
										else
											x2 = punity_eval_local_poly( pu, ltg[j], expr, qx );
									}
									else
									{
										/* Deal with case in which lbase[ltg[j]][jj] < 0 */
										idx = -1 - lbase[ltg[j]][jj];

										/* Translate the vector to origin of particle ltg[i] */
                                                                                for(k=0;k<dim;k++)
                                                                                        tx[k] = qx[k] - pu->pts[ltg[j]*dim+k];

										/* Call the right function pointer and localize it */
                                                                                x2 = rbfd[idx]( dim, 0, pu->dlt[ltg[j]], tx ) * punity_evaluate( pu, ltg[j], qx );
									}

									/* Do the evaluation */
									if( !b_load_overl_mat )
										sum += qw * x1 * x2;
									if( !b_load_stiff_mat && b_use_external_potential )
										rum += qw * x1 * x2 * nuclei_potential( nuc, qx );
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
										if( lbase[ltg[i]][ii] >= 0 )
										{
											global_polynomial_vector( lbase[ltg[i]][ii], dim, expl );
											if( b_use_singular )
												x1 = punity_term_deriv_local_poly_delta( pu, ltg[i], expl, drv, qx, i_sing_order, nb, nnb );
											else
												x1 = punity_term_deriv_local_poly( pu, ltg[i], expl, drv, qx, nb, nnb );
										}
										else
										{
											/* Deal with other case */
											idx = -1 - lbase[ltg[i]][ii];

											/* FIXME: Pick up here... */
											x1 = punity_term_deriv_local_radial( pu, ltg[i], rbfd[idx], drv, qx, nb, nnb );
										}
										if( lbase[ltg[j]][jj] >= 0 )
										{
											global_polynomial_vector( lbase[ltg[j]][jj], dim, expr );
											if( b_use_singular )
												x2 = punity_term_deriv_local_poly_delta( pu, ltg[j], expr, drv, qx, i_sing_order, nb, nnb );
											else
												x2 = punity_term_deriv_local_poly( pu, ltg[j], expr, drv, qx, nb, nnb );
										}
										else
										{
											/* Deal with other case */
											idx = -1 - lbase[ltg[j]][jj];

											/* Derivatives */
                                                                                        x2 = punity_term_deriv_local_radial( pu, ltg[j], rbfd[idx], drv, qx, nb, nnb );
										}
										y += 0.5 * x1 * x2;
									}
									tum += qw * y;
								}
								if( nnb > 0 )
									free( nb );
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
							if( b_use_external_potential )
								A[p] += rum; /* Only add this term in if reading in an external potential */
							ja[p] = qq;
							if( ia[pp] == -1 )
								ia[pp] = p;
						}
						if( !b_load_overl_mat )
						{
							B[p] = sum;
							jb[p] = qq;
							if( ib[pp] == -1 )
								ib[pp] = p;
						}

						/* Step to the next entry */
						p++;
					}
					qq++;
				}
			}
			pp++;
		}
		if( npts >= 20 )
			parse_print_back( stderr, i / ( npts / 20 ) + 1 + 4 );
	}
	if( !b_load_stiff_mat )
		ia[mdim] = p;
	if( !b_load_overl_mat )
		ib[mdim] = p;
	fprintf( stderr, "\n\n" );

	/* Indicate that matrices are both built */
	if( !b_load_stiff_mat )
		*b_have_stiff_mat = 1;
	if( !b_load_overl_mat )
		*b_have_overl_mat = 1;

	/* Pass the variables back */
	*iap = ia; *ibp = ib;
	*jap = ja; *jbp = jb;
	*Ap = A; *Bp = B;

	/* Don't forget to clean up memory */
	free( qx );
	free( tx );
}

