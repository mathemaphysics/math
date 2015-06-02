#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
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
	z = sqrt( z ) / a_in * 2.0;
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

double boundary( double *x_in )
{
	return 0.0;
}

int on_boundary( double *pt_in )
{
        if( fabs( pt_in[0] ) < 1e-6 || fabs( pt_in[1] < 1e-6 || fabs( pt_in[0] - 1.0 ) < 1e-6 || fabs( pt_in[1] - 1.2 ) < 1e-6 ) )
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
	double zrt = 1e-6;
	double sum,scl = 0.3;
	double xx,res;
	double ptmp[dim];
	double tol = zrt;
	double sft = 0.0;
	int b_sft = 0;
	int b_prj = 0;
	int b_pc = 0;
	FILE *fp;

	/* Take in the options */
	int optc;
	while( ( optc = getopt( argc, argv, "d:s:t:S:pPc" ) ) != -1 )
        {
                switch( optc )
                {
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
			case '?':
			default:
				fprintf( stderr, "I don't understand the jibberish you are spouting. Goodbye.\n" );
				return 0;
				break;
		}
	}

	/* Allocate and generate grid points */
	double x = 1.0;
	double y = 1.2;
	int nb;
	int gx = 15; /* These are the number of points needed; NOT the spacing divisor */
	int gy = 15;
	double dx = x / (double) ( gx - 1 );
	double dy = y / (double) ( gy - 1 );
	int npts = gx * gy;
	double *pts = (double*) malloc( npts * dim * sizeof(double) );
	srand((unsigned)time(0));
	for(i=0;i<gx;i++)
		for(j=0;j<gy;j++)
		{
			pts[i*gy*dim+j*dim+0] = (double) i * dx;
			pts[i*gy*dim+j*dim+1] = (double) j * dy;
			//if( i != 0 && j != 0 && i != gx - 1 && j != gy - 1 )
			//{
			//	pts[i*gy*dim+j*dim+0] += 0.02 * ( rand() > RAND_MAX / 2 ? 1.0 : -1.0 ) * (double) rand() / (double) RAND_MAX;
			//	pts[i*gy*dim+j*dim+1] += 0.02 * ( rand() > RAND_MAX / 2 ? 1.0 : -1.0 ) * (double) rand() / (double) RAND_MAX;
			//}
		}
	double *dlt = (double*) malloc( npts * sizeof(double) );
	for(i=0;i<npts;i++)
		dlt[i] = scl;
	double *val = (double*) malloc( npts * sizeof(double) );
	double *gqw = (double*) malloc( npts * sizeof(double) );
	for(i=0;i<npts;i++)
		gqw[i] = 1.0;
	double f,g,alpha,beta,r,s,rr,a,b;
	double *vv,*ww,*pc,*xv,*xvh,*rv,*rvh,*rrv,*rrvh,*pv,*pvh; /* Variables for CG iteration */
	double *amat,*bmat,*amatc,*bmatc,*cmat,*smat,*tmat,*dv,*rhs;

	generate_cloud_supports( 2, npts, pts, 12, 1.2, dlt );
	//for(i=0;i<npts;i++)
	//	fprintf( stderr, "%15.7f", dlt[i] );
	

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

	/* Fix the quadrature weights to start with; hopefully this will help */
	for(i=0;i<gx*gy;i++)
	{
		/* Point centered at pts + i * dim */
		for(j=0;j<gx*gy;j++)
		{
			xx = pow( pts[i*dim+0] - pts[j*dim+0], 2.0 ) + pow( pts[i*dim+1] - pts[j*dim+1], 2.0 );
			if( sqrt( xx ) < scl )
				gqw[i] += 2.0;
		}
		gqw[i] /= (double) ( gx * gy );
	}

	/* Initialize the RKP basis */
	rkp_init( &rk, npts, dim, deg, pts, dlt, gqw, val, cubic, 1.0 );
	rkp_basis_generate( &rk );
	rkp_wavelet_basis_generate( &rk );
	rk.wfsg = cubic_grad;

	/* Do one round of correction by solving 1 / Vi = sum_j w( xi - xj ) */
	//double dl[npts];
	//for(i=0;i<npts;i++)
	//{
	//	sum = 0.0;
	//	for(j=0;j<npts;j++)
	//		sum += rkp_term_evaluate( &rk,i, pts + j * dim );
	//	dl[i] = 1.0 / sum;
	//}
	//fprintf( stderr, "dl = " );
	//for(i=0;i<npts;i++)
	//	fprintf( stderr, "%15.7f", dl[i] );
	//fprintf( stderr, "\n" );
	//for(i=0;i<npts;i++)
	//	rk.gqw[i] = dl[i];

	/* Reinitialize to correct */
	//rkp_basis_generate( &rk );
        //rkp_wavelet_basis_generate( &rk );

	/* Check the cloud */
	int *tot;
	evaluate_cloud( &rk, &tot );
	fprintf( stderr, "Occupations: " );
	for(i=0;i<npts;i++)
		fprintf( stderr, "%10d", tot[i] );
	fprintf( stderr, "\n" );
	free( tot );
	//for(i=0;i<npts;i++)
	//	if( on_boundary( pts + i * dim ) != 1 )
	//		rk.dlts[i] *= 0.6;
        //for(i=0;i<gx*gy;i++)
        //{
        //        /* Point centered at pts + i * dim */
        //        for(j=0;j<gx*gy;j++)
        //        {
        //                xx = pow( pts[i*dim+0] - pts[j*dim+0], 2.0 ) + pow( pts[i*dim+1] - pts[j*dim+1], 2.0 );
        //                if( xx < scl )
        //                        gqw[i] += 1.0;
        //        }
        //        gqw[i] /= (double) ( gx * gy );
        //}
	//rkp_basis_generate( &rk );
        //rkp_wavelet_basis_generate( &rk );
	//evaluate_cloud( &rk, &tot );
        //fprintf( stderr, "Occupations: " );
        //for(i=0;i<npts;i++)
        //        fprintf( stderr, "%10d", tot[i] );
        //fprintf( stderr, "\n" );
        //free( tot );

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

	/* Print out */
	int term = 143;
        for(i=0;i<pgx;i++)
        {
                for(j=0;j<pgy;j++)
                {
                        xx = rkp_term_evaluate( &rk, term, ppts + i * pgy * dim + j * dim );
                        printf( "%15.7f%15.7f%15.7f\n", ppts[i*pgy*dim+j*dim+0], ppts[i*pgy*dim+j*dim+1], xx );
                }
                printf( "\n" );
        }
        fflush( stdout );

	/* Form the Laplace operator correctly with wavelets */
	int ord[2] = { 2, 4 }; /* Indexes which contain the xx and yy second derivatives */

	/* Build the matrices here to simplify things as much as possible */
	for(i=0;i<gx*gy*gx*gy;i++)
		amat[i] = 0.0, bmat[i] = 0.0;
	for(i=0;i<gx*gy;i++)
	{
		for(j=0;j<gx*gy;j++)
		{
			for(k=0;k<dim;k++)
				amat[i*gx*gy+j] -= gqw[j] * rkp_wavelet_term_evaluate_node( &rk, j, ord[k], i );
			bmat[i*gx*gy+j] = gqw[j] * rkp_term_evaluate_node( &rk, j, i );
		}
	}

	for(i=0;i<gx*gy*gx*gy;i++)
		smat[i] = 0.0, tmat[i] = 0.0;;
	for(i=0;i<gx*gy;i++)
		for(j=0;j<gx*gy;j++)
			for(k=0;k<gx*gy;k++)
				smat[i*gx*gy+j] += bmat[k*gx*gy+i] * amat[k*gx*gy+j],
				tmat[i*gx*gy+j] += bmat[k*gx*gy+i] * bmat[k*gx*gy+j];

	/* Output amat and bmat for checking */
	fp = fopen( "amat", "w" );
	fprintf( fp , "# name: amat\n# type: matrix\n# rows: %d\n# columns: %d\n", gx * gy, gx * gy );
	for(i=0;i<gx*gy;i++)
	{
		for(j=0;j<gx*gy;j++)
			fprintf( fp, "%14.7f", amat[i*gx*gy+j] );
		fprintf( fp, "\n" );
	}
	fclose( fp );
	fp = fopen( "bmat", "w" );
        fprintf( fp , "# name: bmat\n# type: matrix\n# rows: %d\n# columns: %d\n", gx * gy, gx * gy );
        for(i=0;i<gx*gy;i++)
        {
                for(j=0;j<gx*gy;j++)
                        fprintf( fp, "%14.7f", bmat[i*gx*gy+j] );
                fprintf( fp, "\n" );
        }
        fclose( fp );

	/* Generate a random initial state vector */
	srand((unsigned)time(0));

	/* Count the number of boundary nodes */
	nb = 0;
	for(i=0;i<gx*gy;i++)
		if( on_boundary( pts + i * dim ) == 1 )
			++nb;

	/* QR factorization first */
	fprintf( stderr, "Running QR factorization... " );
	fflush( stderr );
	for(p=0,k=0;k<gx*gy;k++)
        {
                /* Orthogonalize only the set of boundary points */
                if( on_boundary( pts + k * dim ) != 1 )
                        continue;

                /* Insert the next vector of bmat into vv */
                for(i=0;i<gx*gy;i++)
                        rhs[i] = bmat[k*gx*gy+i];

                /* Now project out all the previous orthogonal vectors stored in vv */
                for(i=0;i<p;i++)
                {
                        xx = 0.0;
                        for(j=0;j<gx*gy;j++)
                                xx += vv[i*gx*gy+j] * rhs[j];
                        for(j=0;j<gx*gy;j++) 
                                rhs[j] -= xx * vv[i*gx*gy+j];
                }

                /* Now normalize rhs and store it in vv[p*gx*gy+...] */
                xx = 0.0;
                for(i=0;i<gx*gy;i++) 
                        xx += rhs[i] * rhs[i];
                xx = sqrt( xx ); 
                for(i=0;i<gx*gy;i++)
                        vv[p*gx*gy+i] = rhs[i] / xx;

                /* Increment p since we got past the boundary test */
                ++p;
        }
	fprintf( stderr, "Done.\n" );
	fflush( stderr );

	/* Now build the subspace projection matrix */
	fprintf( stderr, "Building subspace projection matrix... " );
	fflush( stderr );
	for(i=0;i<gx*gy-nb;i++)
	{
		for(j=0;j<gx*gy;j++)
                	ww[i*gx*gy+j] = ( rand() < RAND_MAX / 2 ? -1.0 : 1.0 ) * (double) rand() / (double) RAND_MAX;
		for(j=0;j<nb;j++)
		{
			xx = 0.0;
			for(k=0;k<gx*gy;k++)
				xx += ww[i*gx*gy+k] * vv[j*gx*gy+k];
			for(k=0;k<gx*gy;k++)
				ww[i*gx*gy+k] -= xx * vv[j*gx*gy+k];
		}
		for(j=0;j<i;j++)
		{
			xx = 0.0;
			for(k=0;k<gx*gy;k++)
                                xx += ww[i*gx*gy+k] * ww[j*gx*gy+k];
                        for(k=0;k<gx*gy;k++)
                                ww[i*gx*gy+k] -= xx * ww[j*gx*gy+k];
		}
		xx = 0.0;
		for(j=0;j<gx*gy;j++)
			xx += ww[i*gx*gy+j] * ww[i*gx*gy+j];
		xx = sqrt( xx );
		for(j=0;j<gx*gy;j++)
			ww[i*gx*gy+j] /= xx;
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
        	fprintf( fp , "# name: amatc\n# type: matrix\n# rows: %d\n# columns: %d\n", gx * gy - nb, gx * gy - nb );
        	for(i=0;i<gx*gy-nb;i++)
        	{
                	for(j=0;j<gx*gy-nb;j++)
                        	fprintf( fp, "%14.7f", amatc[i*gx*gy+j] );
                	fprintf( fp, "\n" );
        	}
        	fclose( fp );
	        fp = fopen( "bmatc", "w" );
        	fprintf( fp , "# name: bmatc\n# type: matrix\n# rows: %d\n# columns: %d\n", gx * gy - nb, gx * gy - nb );
	        for(i=0;i<gx*gy-nb;i++)
        	{
	                for(j=0;j<gx*gy-nb;j++)
                	        fprintf( fp, "%14.7f", bmatc[i*gx*gy+j] );
        	        fprintf( fp, "\n" );
	        }
        	fclose( fp );
		fprintf( stderr, "Done.\n" );
		fflush( stderr );
	}

	/* Build smat and tmat */
	fprintf( stderr, "Building smat and tmat... " );
	fflush( stderr );
	for(i=0;i<gx*gy*gx*gy;i++)
                smat[i] = 0.0, tmat[i] = 0.0;
        for(i=0;i<gx*gy;i++)
                for(j=0;j<gx*gy;j++)
                        for(k=0;k<gx*gy;k++)
                                smat[i*gx*gy+j] += bmat[k*gx*gy+i] * amat[k*gx*gy+j],
                                tmat[i*gx*gy+j] += bmat[k*gx*gy+i] * bmat[k*gx*gy+j];
	fprintf( stderr, "Done.\n" );
	fflush( stderr );

	/* Initialize to a normalized random vector */
        for(i=0;i<gx*gy;i++)
               	rhs[i] = ( rand() < RAND_MAX / 2 ? -1.0 : 1.0 ) * (double) rand() / (double) RAND_MAX;
	res = 0.0;
	for(i=0;i<gx*gy;i++)
		res += rhs[i] * rhs[i];
	res = sqrt( res );
	for(i=0;i<gx*gy;i++)
		rhs[i] /= res;

	/* Now subtract out all vectors in vv from rhs */
        for(k=0;k<nb;k++)
        {
                xx = 0.0;
                for(i=0;i<gx*gy;i++)
                        xx += rhs[i] * vv[k*gx*gy+i];
                for(i=0;i<gx*gy;i++)
                        rhs[i] -= xx * vv[k*gx*gy+i];
        }

	/* Normalize first; this is important */
	res = 0.0;
        for(i=0;i<gx*gy;i++)
		for(j=0;j<gx*gy;j++)
        		res += rhs[i] * tmat[i*gx*gy+j] * rhs[j];
        res = sqrt( res );
        for(i=0;i<gx*gy;i++)
                rhs[i] /= res;

	/* Project out again */
        for(k=0;k<nb;k++)
        {
                xx = 0.0;
                for(i=0;i<gx*gy;i++)
                        xx += rhs[i] * vv[k*gx*gy+i];
                for(i=0;i<gx*gy;i++)
                        rhs[i] -= xx * vv[k*gx*gy+i];
        }

	/* Initialize the preconditioner to the identity */
        for(i=0;i<gx*gy;i++)
                pc[i] = 1.0;

        /* Normalize first; this is important */
        res = 0.0;
        for(i=0;i<gx*gy;i++)
                for(j=0;j<gx*gy;j++)
                        res += rhs[i] * tmat[i*gx*gy+j] * rhs[j];
        res = sqrt( res );
        for(i=0;i<gx*gy;i++)
                rhs[i] /= res;

	/* Now minimize the Rayleigh quotient over the subspace */
	for(m=0;m<100;m++)
	{
		/* Calculate actual eigenvalue */
		rr = r;
        	f = 0.0, g = 0.0;
	        for(i=0;i<gx*gy;i++)
                	for(j=0;j<gx*gy;j++)
        	                f += rhs[i] * smat[i*gx*gy+j] * rhs[j],
	                        g += rhs[i] * tmat[i*gx*gy+j] * rhs[j];
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
		for(i=0;i<gx*gy;i++)
			for(j=0;j<gx*gy;j++)
				cmat[i*gx*gy+j] = amat[i*gx*gy+j] - s * bmat[i*gx*gy+j];
		for(i=0;i<gx*gy;i++)
		{
			dv[i] = 0.0;
			for(j=0;j<gx*gy;j++)
				dv[i] += bmat[i*gx*gy+j] * rhs[j];
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
			for(i=0;i<gx*gy;i++)
				for(j=0;j<gx*gy;j++)
					cmat[i*gx*gy+j] /= pc[i];
			for(i=0;i<gx*gy;i++)
				dv[i] /= pc[i];
		}

		/* Initialize the preconditioner to the identity */
                for(i=0;i<gx*gy;i++)
                        pc[i] = 1.0;

		/* Initialize xv */
		for(i=0;i<gx*gy;i++)
			xv[i] = rhs[i];

		/* Subtract constraints from xv */
                for(i=0;i<nb;i++)
                {
                        b = 0.0;
                        for(j=0;j<gx*gy;j++)
                                b += xv[j] * vv[i*gx*gy+j];
                        for(j=0;j<gx*gy;j++)
                                xv[j] -= b * vv[i*gx*gy+j];
                }

		/* Set conjugate xvh to xv to begin with */
		for(i=0;i<gx*gy;i++)
			xvh[i] = xv[i];

		/* Build the residual and its conjugate */
		for(i=0;i<gx*gy;i++)
			rv[i] = dv[i], rvh[i] = dv[i];
		for(i=0;i<gx*gy;i++)
			for(j=0;j<gx*gy;j++)
				rv[i] -= cmat[i*gx*gy+j] * xv[j],
				rvh[i] -= cmat[j*gx*gy+i] * xvh[j];

		/* Subtract constraints again */
		for(i=0;i<nb;i++)
                {
                        a = 0.0, b = 0.0;
                        for(j=0;j<gx*gy;j++)
                                a += rvh[j] * vv[i*gx*gy+j],
				b += rv[j] * vv[i*gx*gy+j];
                        for(j=0;j<gx*gy;j++)
				rvh[j] -= a * vv[i*gx*gy+j],
                                rv[j] -= b * vv[i*gx*gy+j];
                }

		/* Set the initial search direction */
		for(i=0;i<gx*gy;i++)
                        pv[i] = rv[i] / pc[i],
			pvh[i] = rvh[i] / pc[i];

		/* Start bleeding internally and die... slowly... */
		for(p=0;p<40*gx*gy;p++)
		{
			/* Calculate alpha */
			a = 0.0, b = 0.0;
			for(i=0;i<gx*gy;i++)
				a += rvh[i] * rv[i] / pc[i];
			for(i=0;i<gx*gy;i++)
				for(j=0;j<gx*gy;j++)
					b += pvh[i] * cmat[i*gx*gy+j] * pv[j];
			alpha = a / b;

			/* Take a step */
			for(i=0;i<gx*gy;i++)
				xv[i] += alpha * pv[i],
				xvh[i] += alpha * pvh[i];

			/* Update the residuals */
			for(i=0;i<gx*gy;i++)
				rrv[i] = rv[i],
				rrvh[i] = rvh[i];
			for(i=0;i<gx*gy;i++)
				for(j=0;j<gx*gy;j++)
					rrv[i] -= alpha * cmat[i*gx*gy+j] * pv[j],
					rrvh[i] -=  alpha * cmat[j*gx*gy+i] * pvh[j];

			/* Subtract constraints from xv */
                	for(i=0;i<nb;i++)
                	{
                        	a = 0.0, b = 0.0;
                        	for(j=0;j<gx*gy;j++)
					a += rrvh[j] * vv[i*gx*gy+j],
                         	     	b += rrv[j] * vv[i*gx*gy+j];
                        	for(j=0;j<gx*gy;j++)
					rrvh[j] -= a * vv[i*gx*gy+j],
                        	        rrv[j] -= b * vv[i*gx*gy+j];
                	}

			/* Calculate beta */
			a = 0.0, b = 0.0;
			for(i=0;i<gx*gy;i++)
				a += rrvh[i] * rrv[i] / pc[i],
				b += rvh[i] * rv[i] / pc[i];
			beta = a / b;

			if( fabs( a ) < tol )
				break;

			/* Update search directions */
			for(i=0;i<gx*gy;i++)
				pv[i] = rrv[i] / pc[i] + beta * pv[i],
				pvh[i] = rrvh[i] / pc[i] + beta * pvh[i];

			/* Update to most recent residual */
			for(i=0;i<gx*gy;i++)
				rv[i] = rrv[i],
				rvh[i] = rrvh[i];
		}

		/* Print final iteration and residual */
		fprintf( stderr, "%d iterations: residual = %15.7f\n", p, a );

		/* Update rhs now with xv */
		for(i=0;i<gx*gy;i++)
			rhs[i] = xv[i];

		/* Renormalize */
		res = 0.0;
		for(i=0;i<gx*gy;i++)
			for(j=0;j<gx*gy;j++)
				res += rhs[i] * tmat[i*gx*gy+j] * rhs[j];
		res = sqrt( res );
		for(i=0;i<gx*gy;i++)
			rhs[i] /= res;
	}

	/* Calculate actual eigenvalue */
	f = 0.0, g = 0.0;
	for(i=0;i<gx*gy;i++)
		for(j=0;j<gx*gy;j++)
			f += rhs[i] * amat[i*gx*gy+j] * rhs[j],
			g += rhs[i] * bmat[i*gx*gy+j] * rhs[j];
	fprintf( stderr, "RQ = %15.7f\n", f / g );

	/* Copy the solution into the rk structure */
	for(i=0;i<gx*gy;i++)
		rk.vals[i] = rhs[i];

	/* Output the function */
	fp = fopen( "plot", "w" );
	if( plm )
	{
		for(i=0;i<pgx*pgy;i++)
                {
                        fprintf( fp, "%15.7f%15.7f%15.7f\n", ppts[i*dim+0], ppts[i*dim+1], rkp_evaluate( &rk, ppts + i * dim ) );
                        if( ( i + 1 ) % pgy == 0 )
                                fprintf( fp, "\n" );
                }
	}
	else
	{
		for(i=0;i<gx*gy;i++)
		{
			fprintf( fp, "%15.7f%15.7f%15.7f\n", pts[i*dim+0], pts[i*dim+1], rkp_evaluate_node( &rk, i ) );
			if( ( i + 1 ) % gy == 0 )
        	                fprintf( fp, "\n" );
		}
	}
	fclose( fp );

	free( pts );
	free( dlt );
	free( val );
	free( gqw );
	free( ppts );

	return 0;
}

