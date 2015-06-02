#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "sparse.h"
#include "points.h"
#include "punity.h"
#include "mls.h"

void eigenvalue_qrfactorize( int dim, int npts, double *pts, long *ibmat, long *jbmat, double *bmat, double *rhs, double *vv )
{
	int i,j,k,p;
	double xx;

	for(p=0,k=0;k<npts;k++)
        {
                /* Orthogonalize only the set of boundary points */
                if( on_boundary( pts + k * dim ) != 1 )
                        continue;

                /* Insert the next vector of bmat into vv */
                for(i=0;i<npts;i++)
			rhs[i] = sparse_entry( k, i, npts, npts, ibmat, jbmat, bmat );

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
}

void eigenvalue_subprojmat( int npts, int nb, double *vv, double *ww )
{
	int i,j,k;
	double xx;

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
}

void eigenvalue_initialize( int npts, double *rhs )
{
	int i;
	double res;

	for(i=0;i<npts;i++)
               	rhs[i] = ( rand() < RAND_MAX / 2 ? -1.0 : 1.0 ) * (double) rand() / (double) RAND_MAX;
	res = 0.0;
	for(i=0;i<npts;i++)
		res += rhs[i] * rhs[i];
	res = sqrt( res );
	for(i=0;i<npts;i++)
		rhs[i] /= res;
}

void eigenvalue_normalize( int npts, long *itmat, long *jtmat, double *tmat, double *rhs, double *temp )
{
	int i;
	double res;

	sparse_dgemv( (long) npts, (long) npts, itmat, jtmat, tmat, 1, rhs, 1, temp );
	res = 0.0;
	for(i=0;i<npts;i++)
		res += rhs[i] * temp[i];
	res = sqrt( res );
	for(i=0;i<npts;i++)
		rhs[i] /= res;
}

void eigenvalue_project( int npts, int nb, double *vv, double *rhs )
{
	int i,k;
	double xx;

	for(k=0;k<nb;k++)
        {
                xx = 0.0;
                for(i=0;i<npts;i++)
                        xx += rhs[i] * vv[k*npts+i];
                for(i=0;i<npts;i++)
                        rhs[i] -= xx * vv[k*npts+i];
        }
}

void eigenvalue_solver( int npts, long *iamat, long *jamat, double *amat,
				  long *ibmat, long *jbmat, double *bmat,
				  long *icmat, long *jcmat, double *cmat,
				  long *icmatt, long *jcmatt, double *cmatt,
				  long *ismat, long *jsmat, double *smat,
				  long *itmat, long *jtmat, double *tmat,
				  double *rv, double *rvh, double *rrv, double *rrvh,
				  double *pv, double *pvh, double *xv, double *xvh, double *dv, double *pc,
				  double *rhs, double *temp, double rin, double tol,
				  int nb, int b_pc, double sft, double *vv )
{
	int i,j,k,p,m;
	double a,b,f,g,s,r,rr,res,alpha,beta;

	/* Now minimize the Rayleigh quotient over the subspace */
	r = rin;
	for(m=0;m<100;m++)
	{
		rr = r;
		f = 0.0, g = 0.0;
		sparse_dgemv( (long) npts, (long) npts, ismat, jsmat, smat, 1, rhs, 1, temp );
		for(i=0;i<npts;i++)
			f += rhs[i] * temp[i];
		sparse_dgemv( (long) npts, (long) npts, itmat, jtmat, tmat, 1, rhs, 1, temp );
		for(i=0;i<npts;i++)
			g += rhs[i] * temp[i];
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

		/* Build the shifted matrix */
		for(i=0;i<iamat[npts];i++)
		{
			cmat[i] = amat[i] - s * bmat[i];
			jcmat[i] = jamat[i];
		}
		for(i=0;i<=npts;i++)
			icmat[i] = iamat[i];
		sparse_dgemv( (long) npts, (long) npts, ibmat, jbmat, bmat, 1, rhs, 1, dv );

		/* Build a preconditioner */
		if( b_pc )
		{
			for(i=0;i<npts;i++)
				pc[i] = 0.0;
			for(i=0;i<npts;i++)
				for(j=icmat[i];j<icmat[i+1];j++)
					if( j == icmat[i] || fabs( cmat[j] ) > pc[i] )
						pc[i] = fabs( cmat[j] );

			/* Apply the preconditioner to cmat and dv */
			for(i=0;i<npts;i++)
				for(j=icmat[i];j<icmat[i+1];j++)
					cmat[j] /= pc[i];
			for(i=0;i<npts;i++)
				dv[i] /= pc[i];
		}

		/* Transpose cmat into cmatt */
		sparse_transp( 1, (long) npts, (long) npts, icmat, jcmat, cmat, icmatt, jcmatt, cmatt );

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
		sparse_dgemv( (long) npts, (long) npts, icmat, jcmat, cmat, 1, xv, 1, temp );
		for(i=0;i<npts;i++)
			rv[i] -= temp[i];
		sparse_dgemv( (long) npts, (long) npts, icmatt, jcmatt, cmatt, 1, xvh, 1, temp );
		for(i=0;i<npts;i++)
			rvh[i] -= temp[i];

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
			sparse_dgemv( (long) npts, (long) npts, icmat, jcmat, cmat, 1, pv, 1, temp );
			for(i=0;i<npts;i++)
				b += pvh[i] * temp[i];
			alpha = a / b;

			/* Take a step */
			for(i=0;i<npts;i++)
				xv[i] += alpha * pv[i],
				xvh[i] += alpha * pvh[i];

			/* Update the residuals */
			for(i=0;i<npts;i++)
				rrv[i] = rv[i],
				rrvh[i] = rvh[i];
			sparse_dgemv( (long) npts, (long) npts, icmat, jcmat, cmat, 1, pv, 1, temp );
			for(i=0;i<npts;i++)
				rrv[i] -= alpha * temp[i];
			sparse_dgemv( (long) npts, (long) npts, icmatt, jcmatt, cmatt, 1, pvh, 1, temp );
			for(i=0;i<npts;i++)
				rrvh[i] -= alpha * temp[i];

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
		sparse_dgemv( (long) npts, (long) npts, itmat, jtmat, tmat, 1, rhs, 1, temp );
		for(i=0;i<npts;i++)
			res += rhs[i] * temp[i];
		res = sqrt( res );
		for(i=0;i<npts;i++)
			rhs[i] /= res;
	}
}

void eigenvalue_calculate( int npts, long *iamat, long *jamat, double *amat,
				       long *ibmat, long *jbmat, double *bmat, double *rhs, double *temp, double *ev )
{
	int i;
	double f,g;

	f = 0.0, g = 0.0;
	sparse_dgemv( (long) npts, (long) npts, iamat, jamat, amat, 1, rhs, 1, temp );
	for(i=0;i<npts;i++)
		f += rhs[i] * temp[i];
	sparse_dgemv( (long) npts, (long) npts, ibmat, jbmat, bmat, 1, rhs, 1, temp );
	for(i=0;i<npts;i++)
		g += rhs[i] * temp[i];
	*ev = f / g;
}

double density_function( rkp_t *rk_in, int norb_in, double *coeffs_in, double *x_in, double rho_in, int (*ind_in)(double*) )
{
	int i,j;
	double tmp,sum = 0.0;

	if( ind_in( x_in ) == 0 )
		return 0.0;
	for(i=0;i<norb_in;i++)
	{
		for(j=0;j<rk_in->np;j++)
			rk_in->vals[j] = coeffs_in[i*rk_in->np+j];
		tmp = rkp_evaluate_scaled_const( rk_in, x_in, rho_in );
		sum += tmp * tmp;
	}
	return sum;
}

/**
 * Use a partition of unity to integrate density_function over
 * the entire domain; could use an RKP basis, but it is simpler
 * to use a more easily and quickly evaluated Shepard basis
 */
double density_integral( rkp_t *rk_in, double *coeffs_in, int (*ind_in)(double*), punity_t *pt_in, int norb_in, double *x_in, int qord_in, double *qpts_in, double *qwts_in )
{
	int i,j,k;
	long *index,*size;
	double nqbox[pt_in->dim*pt_in->dim],qp[pt_in->dim],qw,sum,davg,len,smp,ump,wmp,tmp,nrm;

	/* Make sure same dimension in each */
	assert( pt_in->dim == rk_in->dim );

	/* Calculate the average dilation factor; this might not work in general but just try for now */
	davg = 0.0;
	for(i=0;i<rk_in->np;i++)
		davg += rk_in->dlts[i];
	davg /= (double) rk_in->np;
	
	/* Allocate generalized enumeration index vectors */
	index = (long*) malloc( pt_in->dim * sizeof(long) );
	size = (long*) malloc( pt_in->dim * sizeof(long) );

	/* Iterate over all functions in the partition of unity */
	for(i=0;i<pt_in->dim;i++)
		index[i] = 0, size[i] = qord_in - 1;
	sum = 0.0, nrm = 0.0;
	for(i=0;i<pt_in->npts;i++)
	{
		/* Build nqbox local basis */
		for(j=0;j<pt_in->dim;j++)
			for(k=0;k<pt_in->dim;k++)
				if( j != k ) /* Alway have axis-aligned grid! */
					nqbox[j*pt_in->dim+k] = 0.0;
				else
					nqbox[j*pt_in->dim+k] = pt_in->dlt[i];
		do
		{
			sphere_gauss_point( pt_in->dim, pt_in->pts + i * pt_in->dim, pt_in->dlt[i],
					    nqbox, index, qpts_in, qwts_in, qp, &qw );
			len = 0.0;
			for(j=0;j<pt_in->dim;j++)
				len += pow( qp[j] - x_in[j], 2.0 );
			len = sqrt( len );
			ump = density_function( rk_in, norb_in, coeffs_in, qp, davg, ind_in );
			wmp = punity_evaluate( pt_in, i, qp );
			tmp = qw * ump * wmp / len;
			smp = qw * ump * wmp;
			if( !isnan( tmp ) )
				sum += tmp;
			if( !isnan( smp ) )
				nrm += smp;
			
			fprintf( stderr, "%15.7f\n", sum );
		}
		while( arraynext( (long) pt_in->dim, size, index ) != -1 );
	}

	/* Return the global integral */
	return sum / nrm;
}

/**
 * This function requires some sort of method for calculating the value of the
 * Hartree energy at a point internal to the cloud; this can be done using the
 * partition of unity formed by the reproducing kernel basis
 */
void poisson_solver()
{
	
}

