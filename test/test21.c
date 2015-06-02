#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include "mls.h"
#include "linalg.h"
#include "trees.h"

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

double boundary( double *x_in )
{
	return 0.0;
}

int on_boundary( double *pt_in )
{
        if( fabs( pt_in[0] ) < 0.04 || fabs( pt_in[1] ) < 0.04 || fabs( pt_in[0] - 2.0 ) < 0.04 || fabs( pt_in[1] - 2.0 ) < 0.04 )
                return 1;
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

double standard_potential( int dim_in, double *x_in )
{
	int i;
	if( x_in[0] > 0.25 && x_in[0] < 0.75 && x_in[1] > 0.25 && x_in[1] < 0.75 )
		return -10.0;
	else
		return 0.0;
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
	int plm = 0;
	int radn = 10;
	double zrt = 1e-6;
	double sum,scl = 0.3;
	double xx,yy,res;
	double ptmp[dim];
	double tol = zrt;
	double sft = 0.0;
	int b_sft = 0; /* Use argument as shift */
	int b_prj = 0; /* Build projected system matrices */
	int b_pc = 0; /* Build a preconditioner */
	int b_grd = 0; /* Load points from grid file */
	char gfn[1024];
	FILE *fp;

	/* Take in the options */
	int optc;
	while( ( optc = getopt( argc, argv, "g:d:s:t:S:r:pPc" ) ) != -1 )
        {
                switch( optc )
                {
			case 'r':
				radn = atoi( optarg );
				break;
			case 'c':
				b_pc = 1; /* Apply automatic Jacobi preconditioner */
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
			case '?':
			default:
				fprintf( stderr, "I don't understand the jibberish you are spouting. Goodbye.\n" );
				return 0;
				break;
		}
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
	}

	/* Set up the plot points */
        int pgx = 70;
        int pgy = 70;
        double x = 2.0;
        double y = 2.0;
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
                gqw[i] = 1.0; //dlt[i];
	fp = fopen( "support.out", "w" );
	for(i=0;i<npts;i++)
	{
		sum = pow( pts[i*dim+0] - 1.0, 2.0 ) + pow( pts[i*dim+1] - 1.0, 2.0 );
		if( sum < 0.1 * 0.1 )
			fprintf( fp, "%15.7f%15.7f%15.7f\n", pts[i*dim+0], pts[i*dim+1], dlt[i] );
	}
	fclose( fp );

	/* Build variables */
	double f,g,alpha,beta,r,s,rr,a,b;
	double *vv,*ww,*pc,*xv,*xvh,*rv,*rvh,*rrv,*rrvh,*pv,*pvh; /* Variables for CG iteration */
	double *amat,*bmat,*amatc,*bmatc,*cmat,*smat,*tmat,*dv,*rhs;

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
	amat = (double*) malloc( npts * npts * sizeof(double) );
	bmat = (double*) malloc( npts * npts * sizeof(double) );
	amatc = (double*) malloc( npts * npts * sizeof(double) );
	bmatc = (double*) malloc( npts * npts * sizeof(double) );
	cmat = (double*) malloc( npts * npts * sizeof(double) );
	smat = (double*) malloc( npts * npts * sizeof(double) );
	tmat = (double*) malloc( npts * npts * sizeof(double) );
	rhs = (double*) malloc( npts * sizeof(double) );
	dv = (double*) malloc( npts * sizeof(double) );

	/* Initialize the RKP basis */
	rkp_init( &rk, npts, dim, deg, pts, dlt, gqw, val, cubic, 1.0 );
	rkp_basis_generate_scaled_const( &rk );
//	rkp_basis_generate_scaled( &rk );
	rkp_wavelet_basis_generate_scaled_const( &rk );
// 	rkp_wavelet_basis_generate( &rk );
	rk.wfsg = cubic_grad;

	/* Check for partition of unity */
	sum = 0.0;
	m = 23;
	ptmp[0] = 1.23;
	ptmp[1] = 0.58;
	for(i=0;i<npts;i++)
	{
		xx = rkp_term_evaluate_scaled_const( &rk, i, ptmp, 0.08 );
		sum += xx;
	}
	fprintf( stderr, "PU Check: %15.7f\n", sum );
	sum = 0.0;
	for(i=0;i<npts;i++)
	{
		xx = rkp_wavelet_term_evaluate_node_scaled_const( &rk, i, 4, m ) / dlt[m] / dlt[m];
		sum += xx;
	}
	fprintf( stderr, "PN Check: %15.7f\n", sum );

	/* Check the cloud */
	int *tot;
	evaluate_cloud( &rk, &tot );
	fprintf( stderr, "Occupations: " );
	for(i=0;i<npts;i++)
		fprintf( stderr, "%10d", tot[i] );
	fprintf( stderr, "\n" );
	free( tot );

	/* Print out */
	int term;
	srand((unsigned)time(0));
	if( plm )
	{
		m = 23;
		for(i=0;i<pgx*pgy;i++)
		{
			//xx = rkp_wavelet_term_evaluate_node_scaled_const( &rk, m, 4, i ) / dlt[i] / dlt[i];
			xx = rkp_term_evaluate_scaled_const( &rk, m, ppts + i * dim, 0.1 );
			printf( "%15.7f%15.7f%15.7f\n", ppts[i*dim+0], ppts[i*dim+1], xx );
			if( (i+1) % pgy == 0 )
				printf( "\n" );
		}
        	fflush( stdout );
	}
	else
	{
		term = rand() % npts;
		m = term;
		for(i=0;i<npts;i++)
		{
			xx = rkp_term_evaluate_node_scaled( &rk, m, i );
                        printf( "%15.7f%15.7f%15.7f\n", pts[i*dim+0], pts[i*dim+1], xx );
		}
		fflush( stdout );
	}

	/* Set up the nuclei */
	int nnuc = 1;
	double nuc[2] = { 1.0, 1.0 };

	/* Form the Laplace operator correctly with wavelets */
	int ord[2] = { 2, 4 }; /* Indexes which contain the xx and yy second derivatives */

	/* Build the matrices here to simplify things as much as possible */
	for(i=0;i<npts*npts;i++)
		amat[i] = 0.0, bmat[i] = 0.0;
	for(i=0;i<npts;i++)
	{
		for(j=0;j<npts;j++)
		{
			for(k=0;k<dim;k++)
				amat[i*npts+j] -= 0.5 * gqw[j] * rkp_wavelet_term_evaluate_node_scaled_const( &rk, j, ord[k], i ) / dlt[i] / dlt[i]; /* NOTE: This is not the Laplacian; it is half of it */
			amat[i*npts+j] += gqw[j] * rkp_term_evaluate_node_scaled_const( &rk, j, i ) * coulomb_potential( dim, pts + i * dim, nnuc, nuc );
			bmat[i*npts+j] = gqw[j] * rkp_term_evaluate_node_scaled_const( &rk, j, i );
		}
	}

	for(i=0;i<npts*npts;i++)
		smat[i] = 0.0, tmat[i] = 0.0;;
	for(i=0;i<npts;i++)
		for(j=0;j<npts;j++)
			for(k=0;k<npts;k++)
				smat[i*npts+j] += bmat[k*npts+i] * amat[k*npts+j],
				tmat[i*npts+j] += bmat[k*npts+i] * bmat[k*npts+j];

	/* Output amat and bmat for checking */
	fp = fopen( "amat", "w" );
	fprintf( fp , "# name: amat\n# type: matrix\n# rows: %d\n# columns: %d\n", npts, npts );
	for(i=0;i<npts;i++)
	{
		for(j=0;j<npts;j++)
			fprintf( fp, "%14.7f", amat[i*npts+j] );
		fprintf( fp, "\n" );
	}
	fclose( fp );
	fp = fopen( "bmat", "w" );
        fprintf( fp , "# name: bmat\n# type: matrix\n# rows: %d\n# columns: %d\n", npts, npts );
        for(i=0;i<npts;i++)
        {
                for(j=0;j<npts;j++)
                        fprintf( fp, "%14.7f", bmat[i*npts+j] );
                fprintf( fp, "\n" );
        }
        fclose( fp );

	/* Generate a random initial state vector */
	srand((unsigned)time(0));

	/* Count the number of boundary nodes */
	nb = 0;
	fp = fopen( "boundary.out", "w" );
	for(i=0;i<npts;i++)
		if( on_boundary( pts + i * dim ) == 1 )
			++nb, fprintf( fp, "%15.7f%15.7f\n", pts[i*dim+0], pts[i*dim+1] );
	fclose( fp );
	fprintf( stderr, "nb = %d\n", nb );

	/* QR factorization first */
	fprintf( stderr, "Running QR factorization... " );
	fflush( stderr );
	for(p=0,k=0;k<npts;k++)
        {
                /* Orthogonalize only the set of boundary points */
                if( on_boundary( pts + k * dim ) != 1 )
                        continue;

                /* Insert the next vector of bmat into vv */
                for(i=0;i<npts;i++)
                        rhs[i] = bmat[k*npts+i];

                /* Now project out all the previous orthogonal vectors stored in vv */
                for(i=0;i<p;i++)
                {
                        xx = 0.0;
                        for(j=0;j<npts;j++)
                                xx += vv[i*npts+j] * rhs[j];
                        for(j=0;j<npts;j++) 
                                rhs[j] -= xx * vv[i*npts+j];
                }

                /* Now normalize rhs and store it in vv[p*npts+...] */
                xx = 0.0;
                for(i=0;i<npts;i++) 
                        xx += rhs[i] * rhs[i];
                xx = sqrt( xx ); 
                for(i=0;i<npts;i++)
                        vv[p*npts+i] = rhs[i] / xx;

                /* Increment p since we got past the boundary test */
                ++p;
        }
	fprintf( stderr, "Done.\n" );
	fflush( stderr );

	/* Now build the subspace projection matrix */
	fprintf( stderr, "Building subspace projection matrix... " );
	fflush( stderr );
	for(i=0;i<npts-nb;i++)
	{
		for(j=0;j<npts;j++)
                	ww[i*npts+j] = ( rand() < RAND_MAX / 2 ? -1.0 : 1.0 ) * (double) rand() / (double) RAND_MAX;
		for(j=0;j<nb;j++)
		{
			xx = 0.0;
			for(k=0;k<npts;k++)
				xx += ww[i*npts+k] * vv[j*npts+k];
			for(k=0;k<npts;k++)
				ww[i*npts+k] -= xx * vv[j*npts+k];
		}
		for(j=0;j<i;j++)
		{
			xx = 0.0;
			for(k=0;k<npts;k++)
                                xx += ww[i*npts+k] * ww[j*npts+k];
                        for(k=0;k<npts;k++)
                                ww[i*npts+k] -= xx * ww[j*npts+k];
		}
		xx = 0.0;
		for(j=0;j<npts;j++)
			xx += ww[i*npts+j] * ww[i*npts+j];
		xx = sqrt( xx );
		for(j=0;j<npts;j++)
			ww[i*npts+j] /= xx;
	}
	fprintf( stderr, "Done.\n" );
	fflush( stderr );

	/* Build projected matrices for output */
	if( b_prj )
	{
		fprintf( stderr, "Calculating subspace projected operators amatc and bmatc... " );
		fflush( stderr );
		for(i=0;i<npts-nb;i++)
		{
			for(j=0;j<npts-nb;j++)
			{
				amatc[i*npts+j] = 0.0;
				bmatc[i*npts+j] = 0.0;
				for(k=0;k<npts;k++)
					for(m=0;m<npts;m++)
						amatc[i*npts+j] += ww[i*npts+k] * amat[k*npts+m] * ww[j*npts+m],
						bmatc[i*npts+j] += ww[i*npts+k] * bmat[k*npts+m] * ww[j*npts+m];
			}
		}
		fp = fopen( "amatc", "w" );
        	fprintf( fp , "# name: amatc\n# type: matrix\n# rows: %d\n# columns: %d\n", npts - nb, npts - nb );
        	for(i=0;i<npts-nb;i++)
        	{
                	for(j=0;j<npts-nb;j++)
                        	fprintf( fp, "%14.7f", amatc[i*npts+j] );
                	fprintf( fp, "\n" );
        	}
        	fclose( fp );
	        fp = fopen( "bmatc", "w" );
        	fprintf( fp , "# name: bmatc\n# type: matrix\n# rows: %d\n# columns: %d\n", npts - nb, npts - nb );
	        for(i=0;i<npts-nb;i++)
        	{
	                for(j=0;j<npts-nb;j++)
                	        fprintf( fp, "%14.7f", bmatc[i*npts+j] );
        	        fprintf( fp, "\n" );
	        }
        	fclose( fp );
		fprintf( stderr, "Done.\n" );
		fflush( stderr );
	}

	/* Build smat and tmat */
	fprintf( stderr, "Building smat and tmat... " );
	fflush( stderr );
	for(i=0;i<npts*npts;i++)
                smat[i] = 0.0, tmat[i] = 0.0;
        for(i=0;i<npts;i++)
                for(j=0;j<npts;j++)
                        for(k=0;k<npts;k++)
                                smat[i*npts+j] += bmat[k*npts+i] * amat[k*npts+j],
                                tmat[i*npts+j] += bmat[k*npts+i] * bmat[k*npts+j];
	fprintf( stderr, "Done.\n" );
	fflush( stderr );

	/* Initialize to a normalized random vector */
        for(i=0;i<npts;i++)
               	rhs[i] = ( rand() < RAND_MAX / 2 ? -1.0 : 1.0 ) * (double) rand() / (double) RAND_MAX;
	res = 0.0;
	for(i=0;i<npts;i++)
		res += rhs[i] * rhs[i];
	res = sqrt( res );
	for(i=0;i<npts;i++)
		rhs[i] /= res;

	/* Now subtract out all vectors in vv from rhs */
        for(k=0;k<nb;k++)
        {
                xx = 0.0;
                for(i=0;i<npts;i++)
                        xx += rhs[i] * vv[k*npts+i];
                for(i=0;i<npts;i++)
                        rhs[i] -= xx * vv[k*npts+i];
        }

	/* Normalize first; this is important */
	res = 0.0;
        for(i=0;i<npts;i++)
		for(j=0;j<npts;j++)
        		res += rhs[i] * tmat[i*npts+j] * rhs[j];
        res = sqrt( res );
        for(i=0;i<npts;i++)
                rhs[i] /= res;

	/* Project out again */
        for(k=0;k<nb;k++)
        {
                xx = 0.0;
                for(i=0;i<npts;i++)
                        xx += rhs[i] * vv[k*npts+i];
                for(i=0;i<npts;i++)
                        rhs[i] -= xx * vv[k*npts+i];
        }

	/* Initialize the preconditioner to the identity */
        for(i=0;i<npts;i++)
                pc[i] = 1.0;

        /* Normalize first; this is important */
        res = 0.0;
        for(i=0;i<npts;i++)
                for(j=0;j<npts;j++)
                        res += rhs[i] * tmat[i*npts+j] * rhs[j];
        res = sqrt( res );
        for(i=0;i<npts;i++)
                rhs[i] /= res;

	/* Now minimize the Rayleigh quotient over the subspace */
	for(m=0;m<100;m++)
	{
		/* Calculate actual eigenvalue */
		rr = r;
        	f = 0.0, g = 0.0;
	        for(i=0;i<npts;i++)
                	for(j=0;j<npts;j++)
        	                f += rhs[i] * smat[i*npts+j] * rhs[j],
	                        g += rhs[i] * tmat[i*npts+j] * rhs[j];
		r = f / g;
		if( fabs( rr - r ) < tol )
			break;

		/* Print step and rayleigh */
                fprintf( stderr, "%d: rq = %15.7f\n", m, r );

		/* Build the equation C y = d in cmat and dv */
		if( m < 5 )
			s = sft;
		else
			s = r;
		for(i=0;i<npts;i++)
			for(j=0;j<npts;j++)
				cmat[i*npts+j] = amat[i*npts+j] - s * bmat[i*npts+j];
		for(i=0;i<npts;i++)
		{
			dv[i] = 0.0;
			for(j=0;j<npts;j++)
				dv[i] += bmat[i*npts+j] * rhs[j];
		}

		/* Build a preconditioner */
		if( b_pc )
		{
			for(i=0;i<npts;i++)
				pc[i] = 0.0;
			for(i=0;i<npts;i++)
				for(j=0;j<npts;j++)
					if( j == 0 || fabs( cmat[i*npts+j] ) > pc[i] )
						pc[i] = fabs( cmat[i*npts+j] );

			/* Apply the preconditioner to cmat and dv */
			for(i=0;i<npts;i++)
				for(j=0;j<npts;j++)
					cmat[i*npts+j] /= pc[i];
			for(i=0;i<npts;i++)
				dv[i] /= pc[i];
		}

		/* Initialize the preconditioner to the identity */
                for(i=0;i<npts;i++)
                        pc[i] = 1.0;

		/* Initialize xv */
		for(i=0;i<npts;i++)
			xv[i] = rhs[i];

		/* Subtract constraints from xv */
                for(i=0;i<nb;i++)
                {
                        b = 0.0;
                        for(j=0;j<npts;j++)
                                b += xv[j] * vv[i*npts+j];
                        for(j=0;j<npts;j++)
                                xv[j] -= b * vv[i*npts+j];
                }

		/* Set conjugate xvh to xv to begin with */
		for(i=0;i<npts;i++)
			xvh[i] = xv[i];

		/* Build the residual and its conjugate */
		for(i=0;i<npts;i++)
			rv[i] = dv[i], rvh[i] = dv[i];
		for(i=0;i<npts;i++)
			for(j=0;j<npts;j++)
				rv[i] -= cmat[i*npts+j] * xv[j],
				rvh[i] -= cmat[j*npts+i] * xvh[j];

		/* Subtract constraints again */
		for(i=0;i<nb;i++)
                {
                        a = 0.0, b = 0.0;
                        for(j=0;j<npts;j++)
                                a += rvh[j] * vv[i*npts+j],
				b += rv[j] * vv[i*npts+j];
                        for(j=0;j<npts;j++)
				rvh[j] -= a * vv[i*npts+j],
                                rv[j] -= b * vv[i*npts+j];
                }

		/* Set the initial search direction */
		for(i=0;i<npts;i++)
                        pv[i] = rv[i] / pc[i],
			pvh[i] = rvh[i] / pc[i];

		/* Start bleeding internally and die... slowly... */
		for(p=0;p<40*npts;p++)
		{
			/* Calculate alpha */
			a = 0.0, b = 0.0;
			for(i=0;i<npts;i++)
				a += rvh[i] * rv[i] / pc[i];
			for(i=0;i<npts;i++)
				for(j=0;j<npts;j++)
					b += pvh[i] * cmat[i*npts+j] * pv[j];
			alpha = a / b;

			/* Take a step */
			for(i=0;i<npts;i++)
				xv[i] += alpha * pv[i],
				xvh[i] += alpha * pvh[i];

			/* Update the residuals */
			for(i=0;i<npts;i++)
				rrv[i] = rv[i],
				rrvh[i] = rvh[i];
			for(i=0;i<npts;i++)
				for(j=0;j<npts;j++)
					rrv[i] -= alpha * cmat[i*npts+j] * pv[j],
					rrvh[i] -=  alpha * cmat[j*npts+i] * pvh[j];

			/* Subtract constraints from xv */
                	for(i=0;i<nb;i++)
                	{
                        	a = 0.0, b = 0.0;
                        	for(j=0;j<npts;j++)
					a += rrvh[j] * vv[i*npts+j],
                         	     	b += rrv[j] * vv[i*npts+j];
                        	for(j=0;j<npts;j++)
					rrvh[j] -= a * vv[i*npts+j],
                        	        rrv[j] -= b * vv[i*npts+j];
                	}

			/* Calculate beta */
			a = 0.0, b = 0.0;
			for(i=0;i<npts;i++)
				a += rrvh[i] * rrv[i] / pc[i],
				b += rvh[i] * rv[i] / pc[i];
			beta = a / b;

			if( fabs( a ) < tol )
				break;

			/* Update search directions */
			for(i=0;i<npts;i++)
				pv[i] = rrv[i] / pc[i] + beta * pv[i],
				pvh[i] = rrvh[i] / pc[i] + beta * pvh[i];

			/* Update to most recent residual */
			for(i=0;i<npts;i++)
				rv[i] = rrv[i],
				rvh[i] = rrvh[i];
		}

		/* Print final iteration and residual */
		fprintf( stderr, "%d iterations: residual = %15.7f\n", p, a );

		/* Update rhs now with xv */
		for(i=0;i<npts;i++)
			rhs[i] = xv[i];

		/* Renormalize */
		res = 0.0;
		for(i=0;i<npts;i++)
			for(j=0;j<npts;j++)
				res += rhs[i] * tmat[i*npts+j] * rhs[j];
		res = sqrt( res );
		for(i=0;i<npts;i++)
			rhs[i] /= res;
	}

	/* Calculate actual eigenvalue */
	f = 0.0, g = 0.0;
	for(i=0;i<npts;i++)
		for(j=0;j<npts;j++)
			f += rhs[i] * amat[i*npts+j] * rhs[j],
			g += rhs[i] * bmat[i*npts+j] * rhs[j];
	fprintf( stderr, "RQ = %15.7f\n", f / g );

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

