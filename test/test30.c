#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include "nuclei.h"
#include "mls.h"
#include "trees.h"
#include "sparse.h"
#include "cfgread.h"
#include "punity.h"

double density_integral();

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

int indicator( double *x_in )
{
	if( x_in[0] < 3.0 && x_in[1] < 3.0 && x_in[0] > 0.0 && x_in[1] > 0.0 )
		return 1;
	else
		return 0;
}

void evaluate_cloud( rkp_t *obj_in, int **tot_out )
{
	int i,j,k;
	double sum,d;

	*tot_out = (int*) malloc( obj_in->np * sizeof(int) );
	for(i=0;i<obj_in->np;i++)
		(*tot_out)[i] = 0;

	for(i=0;i<obj_in->np;i++)
	{
		for(j=0;j<obj_in->np;j++)
		{
			if( j == i )
				continue;
			sum = 0.0;
			for(k=0;k<obj_in->dim;k++)
				d = obj_in->pts[i*obj_in->dim+k] - obj_in->pts[j*obj_in->dim+k],
				sum += d * d;
			sum = sqrt( sum );
			if( sum < obj_in->dlts[i] * obj_in->wrad )
				(*tot_out)[i] = (*tot_out)[i] + 1;
		}
	}
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
	int pdim = binomial( dim + deg, dim );
	int plm = 0;
	int radn = 10;
	int nnz;
	double zrt = 1e-6;
	double sum;
	double xx,res;
	double tol = zrt;
	double sft = 0.0;

	/* Switches */
	int b_sft = 0; /* Use argument as shift */
	int b_prj = 0; /* Build projected system matrices */
	int b_pc = 0; /* Build a preconditioner */
	int b_grd = 0; /* Load points from grid file */
	int b_ext = 0; /* Load external potential */

	/* Quadrature */
	int quadn = 6;
	double qpts[6] = { -0.93246951, -0.66120939, -0.23861919, 0.23861919, 0.66120939, 0.93246951 };
	double qwts[6] = {  0.17132449, 0.360767157, 0.46791393, 0.46791393, 0.36076157, 0.17132449 };

	/* Other stuff */
	char gfn[1024];
	FILE *fp;
	rkp_t rk;
	punity_t pt;
	nuclei_t ext;
	symtab_t sym;
	symbol_t ss;

	/* Take in the options */
	int optc;
	while( ( optc = getopt( argc, argv, "C:g:d:s:t:S:r:x:pPc" ) ) != -1 )
        {
                switch( optc )
                {
			case 'C':
				symtab_init( &sym, 1024 );
				cfgread_load_symbols_f( optarg, &sym );
				symtab_print( &sym );
				break;
			case 'r':
				radn = atoi( optarg );
				break;
			case 'c':
				b_pc = 1;
				break;
			case 'd':
				deg = atoi( optarg );
				break;
			case 't':
				tol = atof( optarg );
				break;
			case 'S':
				sft = atof( optarg );
				b_sft = 1;
				break;
			case 'p':
				plm = 1;
				break;
			case 'P':
				b_prj = 1;
				break;
			case 'g':
				b_grd = 1;
				strcpy( gfn, optarg );
				break;
			case 'x':
				b_ext = 1;
				load_nuclei( optarg, &ext, 0 );
				nuclei_print( &ext );
				break;
			case '?':
			default:
				fprintf( stderr, "I don't understand the jibberish you are spouting. Goodbye.\n" );
				return 0;
				break;
		}
	}
	if( !b_ext )
	{
		fprintf( stderr, "Please specify an external potential. Exiting.\n" );
		return 0;
	}

	/* Allocate and generate grid points */
	int nb;
	int npts;
	double *pts,*dlt,*val,*gqw;

	/* Load points from file if switch is set */
	if( b_grd )
	{
		fprintf( stderr, "Loading points from file... " );
		load_points( gfn, &dim, &npts, &pts, 0 );
		fprintf( stderr, "dim = %d npts = %d. Done.\n", dim, npts );
		fflush( stderr );
	}

	/* Set up the plot points */
        int pgx = 70;
        int pgy = 70;
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

	/* Build variables */
	double f,g,alpha,beta,r,s,rr,a,b;
	double *vv,*ww,*pc,*xv,*xvh,*rv,*rvh,*rrv,*rrvh,*pv,*pvh; /* Variables for CG iteration */
	double *amat,*bmat,*cmat,*cmatt,*smat,*tmat,*temp,*dv,*rhs;
	long *iamat,*jamat,*ibmat,*jbmat,*icmat,*jcmat,*icmatt,*jcmatt,*ismat,*jsmat,*itmat,*jtmat,*list;

	/* Calculate the number of non-zeros */
	for(nnz=0,i=0;i<npts;i++)
	{
		for(j=0;j<npts;j++)
		{
			xx = 0.0;
			for(k=0;k<dim;k++)
				xx += pow( pts[i*dim+k] - pts[j*dim+k], 2.0 );
			if( xx < pow( dlt[j], 2.0 ) )
				++nnz;
		}
	}
	nnz = nnz + 1;

	/* Allocate this stuff static for now, but make dynamic ASAP */
	iamat = (long*) malloc( ( npts + 1 ) * sizeof(long) );
	jamat = (long*) malloc( nnz * sizeof(long) );
	ibmat = (long*) malloc( ( npts + 1 ) * sizeof(long) );
	jbmat = (long*) malloc( nnz * sizeof(long) );
	icmat = (long*) malloc( ( npts + 1 ) * sizeof(long) );
	jcmat = (long*) malloc( nnz * sizeof(long) );
	icmatt = (long*) malloc( ( npts + 1 ) * sizeof(long) );
	jcmatt = (long*) malloc( nnz * sizeof(long) );
	ismat = (long*) malloc( ( npts + 1 ) * sizeof(long) );
	jsmat = (long*) malloc( npts * npts * sizeof(long) );
	itmat = (long*) malloc( ( npts + 1 ) * sizeof(long) );
	jtmat = (long*) malloc( npts * npts * sizeof(long) );
	list = (long*) malloc( ( npts + 1 ) * sizeof(long) );
	temp = (double*) malloc( ( npts + 1 ) * sizeof(double) );

	/* Allocate stuff; can't use stack for this garbage!!! */
	vv = (double*) malloc( npts * npts * sizeof(double) );
        ww = (double*) malloc( npts * npts * sizeof(double) );
	pc = (double*) malloc( npts * sizeof(double) );
	xv = (double*) malloc( npts * sizeof(double) );
	xvh = (double*) malloc( npts * sizeof(double) );
	rv = (double*) malloc( npts * sizeof(double) );
	rvh = (double*) malloc( npts * sizeof(double) );
	rrv = (double*) malloc( npts * sizeof(double) );
	rrvh = (double*) malloc( npts * sizeof(double) );
	pv = (double*) malloc( npts * sizeof(double) );
	pvh = (double*) malloc( npts * sizeof(double) );

	/* Allocate matrix storage */
	amat = (double*) malloc( nnz * sizeof(double) );
	bmat = (double*) malloc( nnz * sizeof(double) );
	cmat = (double*) malloc( nnz * sizeof(double) );
	cmatt = (double*) malloc( nnz * sizeof(double) );
	smat = (double*) malloc( npts * npts * sizeof(double) );
	tmat = (double*) malloc( npts * npts * sizeof(double) );
	rhs = (double*) malloc( npts * sizeof(double) );
	dv = (double*) malloc( npts * sizeof(double) );

	/* Initialize the RKP basis */
	rkp_init( &rk, npts, dim, deg, pts, dlt, gqw, val, cubic, 1.0 );
	rkp_basis_generate_scaled_const( &rk );
	rkp_wavelet_basis_generate_scaled_const( &rk );

	/* Form the Laplace operator correctly with wavelets */
	int ord[2] = { 2, 4 }; /* Indexes which contain the xx and yy second derivatives */

	/* Build the matrices here to simplify things as much as possible */
	for(i=0;i<nnz;i++) /* Don't need this many entries anymore */
		amat[i] = 0.0, bmat[i] = 0.0;
	for(i=0;i<npts;i++)
		iamat[i] = -1, ibmat[i] = -1;
	for(p=0,i=0;i<npts;i++)
	{
		for(j=0;j<npts;j++)
		{
			/* Check to see if window function at j is large enough to be nonzero */
			xx = 0.0;
			for(k=0;k<dim;k++)
				xx += pow( pts[i*dim+k] - pts[j*dim+k], 2.0 );
			if( xx > pow( dlt[j] * rk.wrad, 2.0 ) )
				continue;
			if( iamat[i] == -1 )
				iamat[i] = p; /* The beginning of row i entries in amat */
			if( ibmat[i] == -1 )
				ibmat[i] = p;
			jamat[p] = j;
			jbmat[p] = j;
			for(k=0;k<dim;k++)
				amat[p] -= 2.0 * 0.5 * gqw[j] * rkp_wavelet_term_evaluate_node_scaled_const( &rk, j, ord[k], i ) / dlt[i] / dlt[i];
			amat[p] += gqw[j] * rkp_term_evaluate_node_scaled_const( &rk, j, i ) * coulomb_potential( dim, pts + i * dim, ext.np, ext.pts );
			bmat[p] = gqw[j] * rkp_term_evaluate_node_scaled_const( &rk, j, i );
			++p;
		}
	}
	iamat[npts] = p;
	ibmat[npts] = p;

	fprintf( stderr, "nnz = %d\n", p );
	fflush( stderr );

	/* Generate a random initial state vector */
	srand((unsigned)time(0));

	/* Count the number of boundary nodes */
	nb = 0;
	for(i=0;i<npts;i++)
		if( on_boundary( pts + i * dim ) == 1 )
			++nb;

	/* QR factorization first */
	fprintf( stderr, "Running QR factorization... " );
	fflush( stderr );
        eigenvalue_qrfactorize( dim, npts, pts, ibmat, jbmat, bmat, rhs, vv );
	fprintf( stderr, "Done.\n" );
	fflush( stderr );

	/* Now build the subspace projection matrix */
	fprintf( stderr, "Building subspace projection matrix... " );
	fflush( stderr );
	eigenvalue_subprojmat( npts, nb, vv, ww );
	fprintf( stderr, "Done.\n" );
	fflush( stderr );

	/* Build smat and tmat */
	fprintf( stderr, "Building smat and tmat... " );
	fflush( stderr );
	sparse_symbmm( (long) npts, (long) npts, (long) npts,
		ibmat, jbmat, iamat, jamat, ismat, jsmat, list );
	sparse_dgemm( (long) npts, (long) npts, (long) npts,
		ibmat, jbmat, bmat, iamat, jamat, amat, ismat, jsmat, smat, temp );
	sparse_symbmm( (long) npts, (long) npts, (long) npts,
		ibmat, jbmat, ibmat, jbmat, itmat, jtmat, list );
	sparse_dgemm( (long) npts, (long) npts, (long) npts,
		ibmat, jbmat, bmat, ibmat, jbmat, bmat, itmat, jtmat, tmat, temp );
	fprintf( stderr, "Done.\n" );
	fflush( stderr );

	/* Initialize to a normalized random vector */
	eigenvalue_initialize( npts, rhs );

	/* Now subtract out all vectors in vv from rhs */
        eigenvalue_project( npts, nb, vv, rhs );

	/* Normalize first using tmat */
	eigenvalue_normalize( npts, itmat, jtmat, tmat, rhs, temp );

	/* Initialize the preconditioner to the identity */
        for(i=0;i<npts;i++)
                pc[i] = 1.0;

	/* Solve the eigenvalue problem actually */
	eigenvalue_solver( npts, iamat, jamat, amat, ibmat, jbmat, bmat, icmat, jcmat, cmat,
				icmatt, jcmatt, cmatt, ismat, jsmat, smat, itmat, jtmat, tmat,
				rv, rvh, rrv, rrvh, pv, pvh, xv, xvh, dv, pc, rhs, temp, r, tol,
				nb, b_pc, sft, vv );
	
	/* Calculate actual eigenvalue */
	eigenvalue_calculate( npts, iamat, jamat, amat, ibmat, jbmat, bmat, rhs, temp, &f );
	fprintf( stderr, "RQ = %15.7f\n", f );

	/* Calculate the density and do the integral */
 	fprintf( stderr, "punity_init returns %d\n", punity_init( &pt, dim, npts, pts, rk.dlts, &cubic_window, &cubic_window_deriv ) );
	temp[0] = 1.0; temp[1] = 1.0;
	sum = density_integral( &rk, rhs, indicator, &pt, 1, temp, quadn, qpts, qwts );
	fprintf( stderr, "density integral = %15.7f\n", sum );

	/* Copy the solution into the rk structure */
	for(i=0;i<npts;i++)
		rk.vals[i] = rhs[i];

	/* Output the function */
	if( plm )
	{
		fp = fopen( "plot", "w" );
		for(i=0;i<pgx*pgy;i++)
                {
                        fprintf( fp, "%15.7f%15.7f%15.7f\n", ppts[i*dim+0], ppts[i*dim+1], rkp_evaluate_scaled_const( &rk, ppts + i * dim, 0.2 ), 0.1 );
                        if( ( i + 1 ) % pgy == 0 )
                                fprintf( fp, "\n" );
                }
		fclose( fp );
	}
	else
	{
		fp = fopen( "plot", "w" );
		for(i=0;i<npts;i++)
			fprintf( fp, "%15.7f%15.7f%15.7f\n", pts[i*dim+0], pts[i*dim+1], rkp_evaluate_node_scaled_const( &rk, i ) );
	}

	free( pts );
	free( dlt );
	free( val );
	free( gqw );
	free( ppts );

	return 0;
}

