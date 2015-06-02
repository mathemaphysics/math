#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "trees.h"
#include "config.h"

#define ZERO_THRESHOLD 1e-6

typedef struct
{
	int dim;
	int type;
	double *orig;
	double *params;
} shape_t;

/**
 * Print a point to the given output
 */
void point_print( FILE *fp_in, int dim_in, double *x_in )
{
	int i;
	for(i=0;i<dim_in;i++)
		fprintf( fp_in, "%15.7f", x_in[i] );
	fprintf( fp_in, "\n" );
}

/**
 * Load simple domain information
 */
void load_domain( char *str_in, shape_t *dm_in, int *ret )
{
	int i,n;
	char buf[1024],bug[1024],*tok[1024],*tol[1024];

	strncpy( buf, str_in, 1024 );
	n = parse_stokenize( buf, tok, ":" );
	if( n != 2 )
		*ret = -1;
	else
	{
		if( strncmp( tok[0], "cube", 1024 ) == 0 )
		{
			dm_in->type = 0;
			strncpy( bug, tok[1], 1024 );
			n = parse_stokenize( bug, tol, "," );
			dm_in->dim = n - 1;
			dm_in->orig = (double*) malloc( ( n - 1 ) * sizeof(double) );
			dm_in->params = (double*) malloc( sizeof(double) );
			for(i=0;i<n-1;i++)
				dm_in->orig[i] = atof( tol[i] );
			dm_in->params[0] = atof( tol[n-1] );
			*ret = 0;
		}
		if( strncmp( tok[0], "box", 1024 ) == 0 )
		{
			dm_in->type = 1;
		}
		if( strncmp( tok[0], "sphere", 1024 ) == 0 )
        	{
			dm_in->type = 2;
		}
		if( strncmp( tok[0], "ellipsoid", 1024 ) == 0 )
		{
			dm_in->type = 3;
		}
	}
}

/**
 * Function which loads a tensor product grid
 * subdivision of variable dimension
 */
void load_grid( char *str_in, int **div, int *ret )
{
	int i,n;
	char buf[1024],*tok[1024];

	/* Initialize return to 0 */
	*ret = 0;

	/* Copy into new buffer so data not changed */
	strncpy( buf, str_in, 1024 );
	n = parse_stokenize( buf, tok, "," );
	if( n < 1 )
		*ret = -1;
	else
	{
		*div = (int*) malloc( n * sizeof(int) );
		if( *div == NULL )
			*ret = -2;
		else
		{
			for(i=0;i<n;i++)
				(*div)[i] = atoi( tok[i] );
			*ret = n;
		}
	}
}

/**
 * Generic boundary indicator function
 */
void boundary_indicator( double *x, shape_t *sh, double tol, int *ret )
{
	int i;
	double sum;

	switch( sh->type )
	{
		case 0:
			*ret = 0;
			for(i=0;i<sh->dim;i++)
			{
				if( fabs( x[i] - sh->orig[i] ) < tol || fabs( x[i] - ( sh->orig[i] + sh->params[0] ) ) < tol )
				{
					*ret = 1;
					break;
				}
			}
			break;
		case 1:
			*ret = 0;
			for(i=0;i<sh->dim;i++)
			{
				if( fabs( x[i] - sh->orig[i] ) < tol || fabs( x[i] - ( sh->orig[i] + sh->params[i] ) ) < tol )
				{
					*ret = 1;
					break;
				}
			}
			break;
		case 2:
			sum = 0.0;
			for(i=0;i<sh->dim;i++);
				sum += pow( x[i] - sh->orig[i], 2.0 );
			if( fabs( sum - sh->params[0] * sh->params[0] ) < tol * tol )
				*ret = 1;
			else
				*ret = 0;
			break;
		case 3:
			break;
		default:
			*ret = 0;
			break;
	}
}

/**
 * Generic domain indicator function
 */
void domain_indicator( double *x, shape_t *sh, int *ret )
{
	int i;
	double sum;

	switch( sh->type )
	{
		case 0:
			*ret = 1;
			for(i=0;i<sh->dim;i++)
			{
				if( x[i] < sh->orig[i] || x[i] > sh->params[0] + sh->orig[i] )
				{
					*ret = 0;
					break;
				}
			}
			break;
		case 1:
			*ret = 1;
			for(i=0;i<sh->dim;i++)
			{
				if( x[i] < sh->orig[i] || x[i] > sh->params[i] + sh->orig[i] )
				{
					*ret = 0;
					break;
				}
			}
			break;
		case 2:
			sum = 0.0;
			for(i=0;i<sh->dim;i++);
				sum += pow( x[i] - sh->orig[i], 2.0 );
			if( sum < sh->params[0] * sh->params[0] )
				*ret = 1;
			else
				*ret = 0;
			break;
		case 3:
			break;
		default:
			*ret = 0;
			break;
	}
}

/**
 * Generates a set of basis vectors which span the space
 * perpendicular to the input vector
 * @param dim_in The space dimension
 * @param vec_in Vector to which basis must be perpendicular
 * @param basis_out Output space for dim_in output vectors
 */
void local_axes( int dim_in, double *vec_in, double *basis_out )
{
	int i,j,k;
	double sum;

	sum = 0.0;
	for(i=0;i<dim_in;i++)
		basis_out[i] = vec_in[i], sum += basis_out[i] * basis_out[i];
	sum = sqrt( sum );
	for(i=0;i<dim_in;i++)
		basis_out[i] /= sum;
	for(i=1;i<dim_in;i++)
	{
		/* This is a horrible idea */
		//for(j=0;j<dim_in;j++)
		//	basis_out[i*dim_in+j] = ( j == i ? 1.0 : 0.0 );
		for(j=0;j<dim_in;j++)
			basis_out[i*dim_in+j] = 0.5 - ( (double) rand() / (double) RAND_MAX );
		sum = 0.0;
		for(j=0;j<dim_in;j++)
			sum += basis_out[i*dim_in+j] * basis_out[i*dim_in+j];
		sum = sqrt( sum );
		for(j=0;j<dim_in;j++)
			basis_out[i*dim_in+j] /= sum;
		for(j=0;j<i;j++)
		{
			sum = 0.0;
			for(k=0;k<dim_in;k++)
				sum += basis_out[i*dim_in+k] * basis_out[j*dim_in+k];
			for(k=0;k<dim_in;k++)
				basis_out[i*dim_in+k] -= sum * basis_out[j*dim_in+k];
		}
		sum = 0.0;
		for(j=0;j<dim_in;j++)
			sum += basis_out[i*dim_in+j] * basis_out[i*dim_in+j];
		sum = sqrt( sum );
		for(j=0;j<dim_in;j++)
			basis_out[i*dim_in+j] /= sum;
	}
}

/**
 * Calculate a rectangular box about the intersection of two spheres
 * of potentially different radii
 * @param dim_in Dimension in which the spheres live
 * @param ctr1_in Center of the first sphere
 * @param rad1_in Radius of the first sphere
 * @param ctr2_in Center of the second sphere
 * @param rad2_in Radius of the second sphere
 * @param ctr_out Center of the disc containing the intersection
 * @param rad_out Radius of the disc containing the intersection
 * @param qbox_out Contains a local basis covering the intersection of the spheres
 * @return Returns 1 if spheres intersect, 0 if not
 */
int sphere_intersection( int dim_in, double *ctr1_in, double rad1_in, double *ctr2_in, double rad2_in, double *ctr_out, double *rad_out, double *qbox_out )
{
	int i,j;
	double r,rr,sum;
	double *ax = (double*) malloc( dim_in * sizeof(double) );

	/* Calculate the origin and axis of the cylinder */
	sum = 0.0;
	for(i=0;i<dim_in;i++)
		ax[i] = ctr2_in[i] - ctr1_in[i], sum += ax[i] * ax[i];
	sum = sqrt( sum );
	if( sum > rad1_in + rad2_in )
		return 0;

	/* If circles are identical */
	if( sum < ZERO_THRESHOLD )
	{
		/* Then build a box around the smaller sphere */
		for(i=0;i<dim_in;i++)
			ctr_out[i] = ctr1_in[i];
		*rad_out = ( rad1_in < rad2_in ) ? rad1_in : rad2_in;
		for(i=0;i<dim_in;i++)
			for(j=0;j<dim_in;j++)
				qbox_out[i*dim_in+j] = ( i == j ) ? *rad_out : 0.0;
		return 1;
	}

	/* Otherwise */
	for(i=0;i<dim_in;i++)
		ax[i] /= sum;
	r = rad1_in + rad2_in - sum;
	for(i=0;i<dim_in;i++)
		ctr_out[i] = ctr1_in[i] + ( sum - rad2_in ) * ax[i];
	for(i=0;i<dim_in;i++)
		ax[i] *= r;

	/* Calculate the radius of the cylinder */
	rr = sqrt( fabs( rad1_in * rad1_in - rad2_in * rad2_in ) );
	if( sum < rr )
		*rad_out = ( rad1_in < rad2_in ) ? rad1_in : rad2_in;
	else
		*rad_out = sqrt( ( -sum + rad2_in - rad1_in ) * ( -sum - rad2_in + rad1_in )
			* ( -sum + rad2_in + rad1_in ) * ( sum + rad2_in + rad1_in ) ) / 2.0 / sum;

	/* Generate the local coordinates */
	local_axes( dim_in, ax, qbox_out );
	for(i=0;i<dim_in;i++)
		qbox_out[i] = ax[i];
	for(i=1;i<dim_in;i++)
		for(j=0;j<dim_in;j++)
			qbox_out[i*dim_in+j] *= (*rad_out);

	/* Move ctr_out to the middle of the box */
	for(i=0;i<dim_in;i++)
		qbox_out[i] *= 0.5, ctr_out[i] += qbox_out[i];

	free( ax );

	return 1;
}

void lens_gauss_point( int dim_in,
		double *ctr1_in, double rad1_in,
		double *ctr2_in, double rad2_in,
		double cr_in, double *nqbox_in, long *index_in,
		double *qpts_in, double *qwts_in,
		double *qp_out, double *qw_out )
{
	int i,j;
	double ssum,x1,x2;
	double dr[dim_in-1],vec[dim_in],wec[dim_in];

	/* Calculate distance between centers */
	ssum = 0.0;
	for(i=0;i<dim_in;i++)
		ssum += pow( ctr2_in[i] - ctr1_in[i], 2.0 );
	ssum = sqrt( ssum );
	
	/* Calculate the limits for each dimension */
	for(i=0;i<dim_in-1;i++)
	{
		dr[i] = cr_in * cr_in;
		for(j=0;j<i;j++)
			dr[i] -= dr[j] * qpts_in[index_in[j]] * dr[j] * qpts_in[index_in[j]];
		dr[i] = sqrt( dr[i] );
		vec[i] = dr[i] * qpts_in[index_in[i]];
	}

	/* Now for the final dimension which is a function of sphere separation */
	x1 = rad1_in * rad1_in;
	for(i=0;i<dim_in-1;i++)
		x1 -= vec[i] * vec[i];
	x1 = sqrt( x1 );
	x2 = rad2_in * rad2_in;
	for(i=0;i<dim_in-1;i++)
		x2 -= vec[i] * vec[i];
	x2 = ssum - sqrt( x2 );

	/* Project the point onto the whole domain */
	for(i=0;i<dim_in;i++)
      		wec[i] = ctr1_in[i]; /* Translate to the origin, ctr1_in */
	for(i=1;i<dim_in;i++)
		for(j=0;j<dim_in;j++)
			wec[j] += vec[i-1] * nqbox_in[i*dim_in+j];

	/* Index along the axis is given by index_in[dim_in-1] */
	for(i=0;i<dim_in;i++)
		qp_out[i] = wec[i] + ( 0.5 * ( x1 + x2 ) + 0.5 * ( x1 - x2 ) * qpts_in[index_in[dim_in-1]] ) * nqbox_in[i];
	*qw_out = qwts_in[index_in[dim_in-1]];
	for(i=1;i<dim_in;i++)
		*qw_out *= qwts_in[index_in[i-1]] * dr[i-1];
	*qw_out *= fabs( x2 - x1 ) / 2.0;
}

void sphere_gauss_point( int dim_in, double *ctr_in, double rad_in, double *nqbox_in, long *index_in, double *qpts_in, double *qwts_in, double *qp_out, double *qw_out )
{
	int i,j;
	double dr[dim_in],vec[dim_in];

	/* Calculate dimension limits */
	for(i=0;i<dim_in;i++)
	{
		dr[i] = rad_in * rad_in;
		for(j=0;j<i;j++)
			dr[i] -= dr[j] * qpts_in[index_in[j]] * dr[j] * qpts_in[index_in[j]];
		dr[i] = sqrt( dr[i] );
		vec[i] = dr[i] * qpts_in[index_in[i]];
	}

	/* Project the point onto the entire domain by affine transformation */
	for(i=0;i<dim_in;i++)
		qp_out[i] = ctr_in[i];
	for(i=0;i<dim_in;i++)
		for(j=0;j<dim_in;j++)
			qp_out[j] += vec[i] * nqbox_in[i*dim_in+j];
	*qw_out = 1.0;
	for(i=0;i<dim_in;i++)
		(*qw_out) *= dr[i] * qwts_in[index_in[i]];
}

void load_quadrature( char *fn, int *quadn, double **qpts, double **qwts, int *ret )
{
	int i,m,n;
	char buf[1024],*tok[1024];
	FILE *fp;

	fp = fopen( fn, "r" );
	if( fp == NULL )
		*ret = -1;
	else
	{
		n = parse_read_line( fp, buf );
		if( n <= 0 )
			*ret = -2;
		else
		{
			m = parse_stokenize( buf, tok, " \t" );
			if( m != 1 )
				*ret = -3;
			else
			{
				/* Now we know how many points and weights to expect */
				*quadn = atoi( tok[0] );

				/* Now allocate enough space inside qpts and qwts */
				*qpts = (double*) malloc( *quadn * sizeof(double) );
				*qwts = (double*) malloc( *quadn * sizeof(double) );

				/* Read the points first */
				n = parse_read_line( fp, buf );
				if( n <= 0 )
					*ret = -4;
				else
				{
					m = parse_stokenize( buf, tok, " \t" );
					if( m != *quadn )
						*ret = -5;
					else
					{
						for(i=0;i<*quadn;i++)
							(*qpts)[i] = atof( tok[i] );
						n = parse_read_line( fp, buf );
						if( n <= 0 )
							*ret = -6;
						else
						{
							m = parse_stokenize( buf, tok, " \t" );
							if( m != *quadn )
								*ret = -7;
							else
							{
								for(i=0;i<*quadn;i++)
									(*qwts)[i] = atof( tok[i] );
								*ret = 0;
							}
						}
					}
				}
			}
		}
	}
}

/**
 * Load set of points from ascii file.
 * @param fn_in File name
 * @param dim_out Dimension of the data loaded
 * @param np_out Number of points loaded
 * @param pts_out The points loaded
 * @param fmt_in Format (to be implemented)
 * @return Returns 0 if all went well
 */
int load_points( char *fn_in, int *dim_out, int *np_out, double **pts_out, int fmt_in )
{
        int i,j,n,m,dim,np;
        char buf[1024],*tok[1024];
        double *pts;
        FILE *fp;
 
        fp = fopen( fn_in, "r" );
	if( fp == NULL )
		return -1;
        n = parse_read_line( fp, buf );
        if( n <= 0 )
                return -1;
        m = parse_stokenize( buf, tok, " \t\n" );
        if( m != 2 )
                return -1;
        dim = atoi( tok[0] );
        np = atoi( tok[1] );
        pts = (double*) malloc( dim * np * sizeof(double) );
        if( pts == NULL )
                return -1;
        for(i=0;i<np;i++)
        {
                n = parse_read_line( fp, buf );
                if( n <= 0 )
                        continue;
                m = parse_stokenize( buf, tok, " \t\n" );
                if( m != dim )
                        continue;
                for(j=0;j<dim;j++)
                        pts[i*dim+j] = atof( tok[j] );
        }

	*dim_out = dim;
	*np_out = np;
	*pts_out = pts;

	return 0;
}

/**
 * Load set of points from ascii file.
 * @param fn_in File name
 * @param dim_in Dimension of the data to save
 * @param np_in Number of points to save
 * @param pts_out The points to save
 * @param fmt_in Format: 0 to save with header, 1 to save without
 * @return Returns 0 if all is well
 */ 
int write_points( char *fn_in, int dim_in, int np_in, double *pts_in, int fmt_in )
{
	int i,j;
	FILE *fp;

	switch( fmt_in )
	{
		case 0:
			fp = fopen( fn_in, "w" );
			if( fp == NULL )
				return -1;
			fprintf( fp, "%d %d\n", dim_in, np_in );
			for(i=0;i<np_in;i++)
			{
				for(j=0;j<dim_in;j++)
					fprintf( fp, "%15.7f", pts_in[i*dim_in+j] );
				fprintf( fp, "\n" );
			}
			fclose( fp );
			break;
		case 1:
			fp = fopen( fn_in, "w" );
                        if( fp == NULL )
                                return -1;
                        for(i=0;i<np_in;i++)
                        {
                                for(j=0;j<dim_in;j++)
                                        fprintf( fp, "%15.7f", pts_in[i*dim_in+j] );
                                fprintf( fp, "\n" );
                        }
                        fclose( fp );
			break;
	}

	return 0;
}

/**
 * Build a cloud of np_in points consistent with the indicator function
 * defining the domain and according to the density function.
 */
int generate_cloud_macqueen( int dim_in, int np_in, int nq_in, int (*ind_in)(double*,void*), double (*dns_in)(double*,void*), double *dom_in, double **pts_out, void *idat_in, void *ddat_in )
{
	int i,k,m,idx;
	int *j = (int*) malloc( np_in * sizeof(int) );
	double r,sum;
	double *pts = (double*) malloc( np_in * dim_in * sizeof(double) );
	double *x = (double*) malloc( dim_in * sizeof(double) );
	double *y = (double*) malloc( dim_in * sizeof(double) );

	/* Seed the generator */
	srand((unsigned)time(0));

	/* Generate the initial distribution of np_in points */
	for(i=0;i<np_in;)
	{
		for(k=0;k<dim_in;k++)
			x[k] = dom_in[2*k+0] + ( (double) rand() / (double) RAND_MAX ) * ( dom_in[2*k+1] - dom_in[2*k+0] );
		if( ind_in( x, idat_in ) != 1 )
			continue;
		else
		{
			for(k=0;k<dim_in;k++)
				pts[i*dim_in+k] = x[k];
			++i;
		}
	}
	for(i=0;i<np_in;i++)
		j[i] = 1;

	/* Now choose nq_in points according to the density function at random */
	for(i=0;i<nq_in;)
	{
		for(k=0;k<dim_in;k++)
			y[k] = dom_in[2*k+0] + ( (double) rand() / (double) RAND_MAX ) * ( dom_in[2*k+1] - dom_in[2*k+0] );
		if( ind_in( y, idat_in ) != 1 || dns_in( y, ddat_in ) < (double) rand() / (double) RAND_MAX )
			continue;
		++i;

		/* Find closest point in pts to y */
		for(k=0;k<np_in;k++)
		{
			sum = 0.0;
			for(m=0;m<dim_in;m++)
				sum += pow( pts[k*dim_in+m] - y[m], 2.0 );
			sum = sqrt( sum );
			if( k == 0 || sum < r )
				r = sum, idx = k;
		}

		/* Now know that point idx is closest to y */
		for(k=0;k<dim_in;k++)
			pts[idx*dim_in+k] = ( (double) j[idx] * pts[idx*dim_in+k] + y[k] ) / (double) ( j[idx] + 1 );
		++j[idx];
	}

	*pts_out = pts;

	free( x );
	free( y );

	return 0;
}

/**
 * Build a cloud of np_in points consistent with the indicator function
 * defining the domain and according to the density function.
 */
int smooth_cloud_macqueen( int dim_in, int np_in, int nq_in, int (*ind_in)(double*,void*), double (*dns_in)(double*,void*), double *dom_in, double *pts_in, void *idat_in, void *ddat_in )
{
	int i,k,m,idx;
	int *j = (int*) malloc( np_in * sizeof(int) );
	double r,sum;
	double *x = (double*) malloc( dim_in * sizeof(double) );
	double *y = (double*) malloc( dim_in * sizeof(double) );

	/* Seed the generator */
	srand((unsigned)time(0));

	/* Initialize */
	for(i=0;i<np_in;i++)
		j[i] = 1;

	/* Now choose nq_in points according to the density function at random */
	for(i=0;i<nq_in;)
	{
		for(k=0;k<dim_in;k++)
			y[k] = dom_in[2*k+0] + ( (double) rand() / (double) RAND_MAX ) * ( dom_in[2*k+1] - dom_in[2*k+0] );
		if( ind_in( y, idat_in ) != 1 || dns_in( y, ddat_in ) < (double) rand() / (double) RAND_MAX )
			continue;
		++i;

		/* Find closest point in pts to y */
		for(k=0;k<np_in;k++)
		{
			sum = 0.0;
			for(m=0;m<dim_in;m++)
				sum += pow( pts_in[k*dim_in+m] - y[m], 2.0 );
			sum = sqrt( sum );
			if( k == 0 || sum < r )
				r = sum, idx = k;
		}

		/* Now know that point idx is closest to y */
		for(k=0;k<dim_in;k++)
			pts_in[idx*dim_in+k] = ( (double) j[idx] * pts_in[idx*dim_in+k] + y[k] ) / (double) ( j[idx] + 1 );
		++j[idx];
	}

	free( x );
	free( y );

	return 0;
}

double lennard_jones( int dim_in, double *x_in, double *y_in, double eps_in, double sig_in, double rcut_in )
{
	int i;
	double rto,sum = 0.0;

	for(i=0;i<dim_in;i++)
		sum += pow( x_in[i] - y_in[i], 2.0 );
	if( sum > rcut_in * rcut_in )
		return 0.0;
	rto = sig_in * sig_in / sum;
	sum = 4.0 * eps_in * ( pow( rto, 6.0 ) - pow( rto, 3.0 ) );

	return sum;
}

double cloud_pair_energy( int dim_in, int np_in, double *pts_in, double *eps_in, double *sig_in, double rcut_in )
{
	int i,j,k;
	double sum;

	sum = 0.0;
	for(i=0;i<np_in;i++)
		for(j=0;j<i;j++)
			sum += lennard_jones( dim_in, pts_in + i * dim_in, pts_in + j * dim_in,
				0.5 * ( eps_in[i] + eps_in[j] ), 0.5 * ( sig_in[i] + sig_in[j] ), rcut_in );

	return sum;
}

double cloud_pair_energy_diff( int dim_in, int idx_in, int np_in, double *pts_in, double *dx_in, double *eps_in, double *sig_in, double rcut_in )
{
	int i,j;
	double sum;
	double x[dim_in];

	sum = 0.0;
	for(i=0;i<np_in;i++)
	{
		if( i == idx_in )
			continue;
		else
			sum -= lennard_jones( dim_in, pts_in + idx_in * dim_in, pts_in + i * dim_in,
				0.5 * ( eps_in[idx_in] + eps_in[i] ), 0.5 * ( sig_in[idx_in] + sig_in[i] ), rcut_in );
	}
	for(i=0;i<dim_in;i++)
                x[i] = pts_in[idx_in*dim_in+i] + dx_in[i];
	for(i=0;i<np_in;i++) 
        { 
                if( i == idx_in ) 
                        continue; 
                else 
                        sum += lennard_jones( dim_in, x, pts_in + i * dim_in, 
                                0.5 * ( eps_in[idx_in] + eps_in[i] ), 0.5 * ( sig_in[idx_in] + sig_in[i] ), rcut_in ); 
        }

	return sum;
}

int smooth_cloud_bubble_mc( int dim_in, int np_in, double *pts_in, int (*ind_in)(double*,void*), double (*dns_in)(double*,void*), double dx_in, double *eps_in, double *sig_in, double temp_in, double rcut_in, int max_in, int vmod_in, void *idat_in, void *ddat_in )
{
	int n,i,idx,jdx;
	int b_dns;
	long acc,tot;
	double sum,en;
	double *u = (double*) malloc( dim_in * sizeof(double) );
	double *v = (double*) malloc( dim_in * sizeof(double) );

	/* Determine whether or not to use density or external potential functions */
	if( dns_in == NULL )
                b_dns = 0;
        else
                b_dns = 1;

	srand((unsigned)time(0));
	tot = 0, acc = 0;
	jdx = -1;
	for(i=0;i<np_in;i++)
		eps_in[i] = 10.0;
	for(n=0;n<max_in;n++)
	{
		/* Out put the acceptance and rejection rates */
		if( n % vmod_in == 0 && tot > 0 )
		{
			sum = (double) acc / (double) tot;
			if( sum < 0.2 )
				dx_in = dx_in * 0.95;
			else if( sum > 0.4 )
				dx_in = dx_in * 1.05;
			fprintf( stderr, "%d: Accept: %4.1f%% dx = %15.7f en = %15.7f\n",
				n, 100.0 * (double) acc / (double) tot, dx_in,
				cloud_pair_energy( dim_in, np_in, pts_in, eps_in, sig_in, rcut_in ) );
			tot = 0, acc = 0;
		}

		/* Generate a trial move */
		idx = rand() % np_in;
		sum = 0.0;
		for(i=0;i<dim_in;i++)
			u[i] = (double) rand() / (double) RAND_MAX - 0.5, sum += u[i] * u[i];
		sum = sqrt( sum );
		for(i=0;i<dim_in;i++)
			u[i] = u[i] / sum * dx_in;
		for(i=0;i<dim_in;i++)
                        v[i] = pts_in[idx*dim_in+i] + u[i];
		if( ind_in( v, idat_in ) != 1 )
		{
			jdx = -1; /* No need to update any sigma parameters */
			continue;
		}
		if( b_dns == 1 ) /* Assume non-uniform values for sigma and epsilon */
		{
			if( jdx != -1 ) /* If jdx != -1, then jdx is the last updated point */
				sig_in[jdx] = dns_in( pts_in + jdx * dim_in, ddat_in );
		}
		en = cloud_pair_energy_diff( dim_in, idx, np_in, pts_in, u, eps_in, sig_in, rcut_in );
		sum = (double) rand() / (double) RAND_MAX;
		if( sum < exp( -1.0 / temp_in * en ) ) /* Reject if true, else keep */
		{
			for(i=0;i<dim_in;i++)
				pts_in[idx*dim_in+i] += u[i]; /* Undo the change */
			++acc;
			jdx = idx; /* Change was made to idx */
		}
		else
			jdx = -1; /* No change was made */
		++tot;
	}

	free( u );
	free( v );

	return 0;
}

int generate_cloud_supports( int dim_in, int np_in, double *pts_in, int deg_in, double gam_in, double *rad_in )
{
	int i,j,k,m;
	int max = 200;
	long tot;
	double sum;

	for(i=0;i<np_in;i++)
		rad_in[i] = 0.001;
	for(i=0;i<np_in;i++)
	{
		for(m=0;m<max;m++)
		{
			tot = 0;
			for(j=0;j<np_in;j++)
			{
				if( j == i )
					continue;
				sum = 0.0;
				for(k=0;k<dim_in;k++)
					sum += pow( pts_in[i*dim_in+k] - pts_in[j*dim_in+k], 2.0 );
				sum = sqrt( sum );
				if( sum < rad_in[i] )
					++tot;
			}
			if( tot < deg_in )
				rad_in[i] *= gam_in;
			else
				break;
		}
	}

	return 0;
}

int generate_cloud_supports_min( int dim_in, int np_in, double *pts_in, int deg_in, double gam_in, double *rad_in, double mrad_in )
{
        int i,j,k,m;
        int max = 200;
        long tot;
        double sum;

        for(i=0;i<np_in;i++)
                rad_in[i] = mrad_in;
        for(i=0;i<np_in;i++)
        {
                for(m=0;m<max;m++)
                {
                        tot = 0;
                        for(j=0;j<np_in;j++)
                        {
                                if( j == i )
                                        continue;
                                sum = 0.0;
                                for(k=0;k<dim_in;k++)
                                        sum += pow( pts_in[i*dim_in+k] - pts_in[j*dim_in+k], 2.0 );
                                sum = sqrt( sum );
                                if( sum < rad_in[i] )
                                        ++tot;
                        }
                        if( tot < deg_in )
                                rad_in[i] *= gam_in;
                        else
                                break;
                }
        }

        return 0;
}

/**
 * Data structure for kd trees for use in nearest neighbor
 * set finding algorithms
 */
struct kdnode_s
{
	int dim;
	int depth;
	int npts;
	int idx;
	double split;
	double *pos;
	struct kdnode_s *left;
	struct kdnode_s *right;
};

typedef struct kdnode_s kdnode_t;

void kdtree_heapify( int dim_in, double *pts_in, int d_in, int idx_in, int size_in )
{
	int i,l,r,lg;
	double tmp;

	l = 2 * ( idx_in + 1 ) - 1;
	r = 2 * ( idx_in + 1 );
	if( l < size_in && pts_in[idx_in*dim_in+d_in] < pts_in[l*dim_in+d_in] )
		lg = l;
	else
		lg = idx_in;
	if( r < size_in && pts_in[lg*dim_in+d_in] < pts_in[r*dim_in+d_in] )
		lg = r;
	if( lg != idx_in )
	{
		for(i=0;i<dim_in;i++)
		{
			tmp = pts_in[lg*dim_in+i];
			pts_in[lg*dim_in+i] = pts_in[idx_in*dim_in+i];
			pts_in[idx_in*dim_in+i] = tmp;
		}
		kdtree_heapify( dim_in, pts_in, d_in, lg, size_in );
	}
}

void kdtree_build_heap( int dim_in, double *pts_in, int d_in, int size_in )
{
	int i,n = size_in / 2 - 1;

	for(i=n;i>=0;i--)
		kdtree_heapify( dim_in, pts_in, d_in, i, size_in );
}

void kdtree_heapsort( int dim_in, double *pts_in, int d_in, int size_in )
{
	int i,j;
	double tmp;

	kdtree_build_heap( dim_in, pts_in, d_in, size_in );
	for(i=size_in-1;i>0;i--)
	{
		for(j=0;j<dim_in;j++)
		{
			tmp = pts_in[0*dim_in+j];
			pts_in[0*dim_in+j] = pts_in[i*dim_in+j];
			pts_in[i*dim_in+j] = tmp;
		}
		kdtree_heapify( dim_in, pts_in, d_in, 0, i );
	}
}

/**
 * @param root_in The k-d tree object
 * @param dim_in Dimension of the point space
 * @param npts_in Number of input points
 * @param pts_in List of points to subdivide by this root node
 * @param depth_in Depth of the tree; also indicates the dimension by which to split
 * @return Value is 0 if all went well
 */
int kdtree_build_median( kdnode_t *root_in, int dim_in, int npts_in, double *pts_in, int depth_in )
{
	int i,j,k,idx,nlp,nrp;
	double tmp,*lpts,*rpts;

	/* First find the median of the input parameter using a bubble sort (optimize later) */
	idx = depth_in % dim_in;
	root_in->dim = dim_in;
	root_in->depth = depth_in;
	root_in->npts = npts_in;
	kdtree_heapsort( dim_in, pts_in, idx, npts_in );

	/* NOTE: Actually only need to iterate i down to half the length of pts_in since we just want median */
	i = npts_in / 2;
	root_in->pos = (double*) malloc( dim_in * sizeof(double) );
	for(j=0;j<dim_in;j++)
		root_in->pos[j] = pts_in[i*dim_in+j];

	/* Build the two point sets to decompose for left and right */
	nlp = i;
	nrp = npts_in - ( i + 1 );
	lpts = pts_in;
	rpts = pts_in + ( i + 1 ) * dim_in;

	/* Now position is put into this root node */
	if( nlp > 0 )
	{
		root_in->left = (kdnode_t*) malloc( sizeof(kdnode_t) );
		kdtree_build_median( root_in->left, dim_in, nlp, lpts, depth_in + 1 );
	}
	if( nrp > 0 )
	{
		root_in->right = (kdnode_t*) malloc( sizeof(kdnode_t) );
		kdtree_build_median( root_in->right, dim_in, nrp, rpts, depth_in + 1 );
	}

	return 0;
}

int kdtree_build_average( kdnode_t *root_in, int dim_in, int npts_in, double *pts_in, int depth_in, int start_in )
{
        int i,j,k,idx,nlp,nrp,dir;
        double tmp,*lpts,*rpts;

	/* Just to be sure */
	if( npts_in <= 0 )
		return 0;

	/* Set basics */
	root_in->dim = dim_in;
	root_in->depth = depth_in;
	root_in->npts = npts_in;
	root_in->pos = NULL;

	/* Check if a leaf node */
	if( npts_in == 1 )
	{
		root_in->idx = start_in;
		root_in->pos = (double*) malloc( dim_in * sizeof(double) );
		for(i=0;i<dim_in;i++)
			root_in->pos[i] = pts_in[i];
		root_in->left = NULL;
		root_in->right = NULL;
		return 0;
	}

        /* First find the median of the input parameter using a bubble sort (optimize later) */
        idx = depth_in % dim_in;
	kdtree_heapsort( dim_in, pts_in, idx, npts_in );

	/* Calculate the average value for split */
	tmp = 0.0;
	for(i=0;i<npts_in;i++)
		tmp += pts_in[i*dim_in+idx];
	tmp /= (double) npts_in;
	root_in->split = tmp;

	/* Find the index to split the points */
	i = npts_in / 2;
	if( pts_in[i*dim_in+idx] >= tmp )
		dir = -1;
	else
		dir = 1;
	for(j=i;j>=0&&j<npts_in;)
	{
		if( dir == -1 )
		{
			if( pts_in[j*dim_in+idx] >= tmp && pts_in[(j-1)*dim_in+idx] < tmp )
				break;
			else
				--j;
		}
		else
		{
			if( pts_in[j*dim_in+idx] < tmp && pts_in[(j+1)*dim_in+idx] >= tmp )
			{
				++j; /* Want an exclusive upper limit so that the upper limit minus lower is the size */
				break;
			}
			else
				++j;
		}	
	}

        /* Build the two point sets to decompose for left and right */
        nlp = j;
        nrp = npts_in - j;
        lpts = pts_in;
        rpts = pts_in + j * dim_in;

        /* Now position is put into this root node */
        if( nlp > 0 ) /* Just for completeness; should never have nlp == 0 */
        {
                root_in->left = (kdnode_t*) malloc( sizeof(kdnode_t) );
                kdtree_build_average( root_in->left, dim_in, nlp, lpts, depth_in + 1, start_in );
        }
        if( nrp > 0 ) /* Same as above */
        {
                root_in->right = (kdnode_t*) malloc( sizeof(kdnode_t) );
                kdtree_build_average( root_in->right, dim_in, nrp, rpts, depth_in + 1, start_in + j );
        }

        return 0;
}

void kdtree_walk( kdnode_t *root_in )
{
	fprintf( stderr, "dim = %5d didx = %5d split = %15.7f npts = %5d depth = %5d idx = %5d\n", root_in->dim, root_in->depth % root_in->dim, root_in->split, root_in->npts, root_in->depth, root_in->idx );
	if( root_in->left != NULL )
		kdtree_walk( root_in->left );
	if( root_in->right != NULL )
		kdtree_walk( root_in->right );
}

int kdtree_free( kdnode_t *root_in )
{
	
}

#define KDTREE_RANGE_QUERY_INC 1024

/**
 * Return all points in the tree which are within r_in
 * units of the input point x_in
 */
void kdtree_range_query( kdnode_t *root_in, double *x_in, double r_in, int **idx_out, int *alloc_out, int *n_out )
{
	int i,p,idx,nidv,*idv;
	double *rn = (double*) malloc( 2 * root_in->dim * sizeof(double) );

	idx = root_in->depth % root_in->dim;
	for(i=0;i<root_in->dim;i++)
	{
		rn[2*i+0] = x_in[i] - r_in;
		rn[2*i+1] = x_in[i] + r_in;
	}
	if( root_in->pos == NULL )
	{
		if( rn[2*idx+0] < root_in->split && rn[2*idx+1] < root_in->split ) /* If range totally on left */
			kdtree_range_query( root_in->left, x_in, r_in, idx_out, alloc_out, n_out );
		else if( root_in->split < rn[2*idx+0] && root_in->split < rn[2*idx+1] ) /* If range totally on right */
			kdtree_range_query( root_in->right, x_in, r_in, idx_out, alloc_out, n_out );
		else /* Otherwise range is split */
		{
			kdtree_range_query( root_in->left, x_in, r_in, idx_out, alloc_out, n_out );
			kdtree_range_query( root_in->right, x_in, r_in, idx_out, alloc_out, n_out );
		}
	}
	else
	{
		p = 0;
		for(i=0;i<root_in->dim;i++)
		{
			if( root_in->pos[i] < rn[2*i+0] || rn[2*i+1] < root_in->pos[i] )
			{
				p = 1;
				break;
			}
		}
		if( p == 0 )
		{
			/* Then add to the list */
			if( (*n_out) + 1 > *alloc_out )
			{
				*idx_out = (int*) realloc( *idx_out, ( *alloc_out + KDTREE_RANGE_QUERY_INC ) * sizeof(int) );
				*alloc_out = *alloc_out + KDTREE_RANGE_QUERY_INC;
			}
			(*idx_out)[*n_out] = root_in->idx;
			*n_out = *n_out + 1;
		}
	}
	free( rn );
}

