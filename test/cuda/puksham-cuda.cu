/**
 * CUDA device code for combinadic functions
 */

typedef struct
{
	int dim;
	int np;
	double *pts;
	double *dlt;
	double *params;
	double *chg;
} nuclei_t;

__device__ void partition_init( int *s, int *m, int n )
{
	int i;
	for(i=0;i<n;i++)
		m[i] = 1, s[i] = 1;
}

__device__ int partition_next( int *s, int *m, int n )
{
	/* Update s: 1 1 1 1 -> 2 1 1 1 -> 1 2 1 1 -> 2 2 1 1 -> 3 2 1 1 -> 1 1 2 1 ... */
	int i = 0;
	++s[i];
	while ((i < n - 1) && (s[i] > m[i] + 1))
	{
		s[i] = 1;
		++i;
		++s[i];
	}

	/* If i is has reached n-1 th element, then the last unique partition has been found*/
	if (i == n - 1)
		return 0;

	/* Because all the first i elements are now 1, s[i] (i + 1 th element)
	is the largest. So we update max by copying it to all the first i
	positions in m.*/
	int max = s[i];
	for (i = i - 1; i >= 0; --i)
		m[i] = max;

	return 1;
}

__device__ int factorial( int num_in )
{
	int fct = 1;
	int i;
	if( num_in <= 0 )
		return 1;
	for(i=num_in;i>1;i--)
		fct *= i;
	return fct;
}

__device__ int truncfact( int n_in, int r_in )
{
	int i,prd = 1;
	for(i=0;i<r_in;i++)
		prd *= n_in - i;
	return prd;
}

__device__ int binomial( int n_in, int k_in )
{
	if( n_in < k_in )
		return 0;
	else
		return truncfact( n_in, k_in ) / factorial( k_in );
}

__device__ int combinadic_init( int lim_in, int dim_in, int *ptr_in )
{
	int i;
	if(lim_in<dim_in)
		return -1; // can't make dim_in-combinations of lim_in (< dim_in) objects
	for(i=0;i<dim_in;i++)
		ptr_in[i] = i;
	return 0;
}

__device__ int combinadic_next( int lim_in, int dim_in, int *vec_in )
{
	int i,j;
	for(i=0;i<dim_in;i++)
	{
		if( ( vec_in[i] < lim_in )
			&& ( ( i < dim_in - 1 && ( vec_in[i+1] - vec_in[i] ) > 1 ) || i == dim_in - 1 ) )
		{
			++vec_in[i];
			break;
		}
	}
	for(j=0;j<i;j++)
		vec_in[j] = j;
	return 0;
}

__device__ int combinadic_vector( int idx_in, int lim_in, int dim_in, int *vec_out )
{
	int i,j,idx,tmp;

	if( idx_in < 0 )
		return -1;

	idx = idx_in;
	for(i=dim_in-1;i>=0;i--)
	{
		for(j=i+1;j<=lim_in;j++)
		{
			tmp = binomial( j, i + 1 );
			if( tmp > idx )
				break;
		}
		idx -= binomial( j - 1, i + 1 );
		vec_out[i] = j - 1;
	}
	return 0;
}

__device__ int rcombinadic_vector( int idx_in, int lim_in, int dim_in, int *vec_out )
{
	return combinadic_vector( idx_in, lim_in + dim_in - 1, dim_in - 1, vec_out );
}

__device__ int rcombinadic_occupancy( int lim_in, int dim_in, int *vec_in, int *occ_out )
{
	int i;
	if( dim_in < 1 )
		return -1; /* Failure: Parameter out of acceptable domain */
	if( dim_in > 1 )
	{
		occ_out[0] = vec_in[0];
		for(i=1;i<dim_in-1;i++)
			occ_out[i] = vec_in[i] - vec_in[i-1] - 1;
		occ_out[dim_in-1] = lim_in + dim_in - 1 - vec_in[dim_in-2] - 1; /* minus addtl 'one' adjust difference of 1 to mean zero occupancy */
	}
	else
		occ_out[0] = lim_in;
	return 0;
}

__device__ int polynomial_exponents( int idx_in, int ord_in, int dim_in, int *occ_out )
{
	int *vec = (int*) malloc( ( dim_in-1 ) * sizeof(int) );
	rcombinadic_vector( idx_in, ord_in, dim_in, vec );
	rcombinadic_occupancy( ord_in, dim_in, vec, occ_out );
	free( vec );
	return 0;
}

__device__ int global_polynomial_vector( int idx_in, int dim_in, int *exp_out )
{
	int n,m,r;
	int *cmb = (int*) malloc( ( dim_in - 1 ) * sizeof(int) );

	for(n=0,r=0;;r++)
	{
		m = binomial( r + dim_in - 1, dim_in - 1 );
		if( n + m > idx_in )
			break;
		else
			n += m;
	}

	polynomial_exponents( idx_in - n, r, dim_in, exp_out );

	free( cmb );

	return 0;
}

/**
 * CUDA device code for stepping through a grid in d dimensions
 */

__device__ int arraynext( long dim_in, long *size_in, long *index_in )
{
	long i;
	for(i=dim_in-1;i>=0;i--)
	{
		if( index_in[i] < size_in[i] )
		{
			++index_in[i];
			return 0;
		}
		else /* Carry to the next spot */
			index_in[i] = 0;
	}
	return -1;
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define ZERO_THRESHOLD 1e-6

#define CUDA_RAND_MAX 2147483648

typedef struct
{
	int dim;
	int type;
	double *orig;
	double *params;
} shape_t;

__device__ void domain_indicator( double *x, shape_t *sh, int *ret )
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
 * A random number generator for use inside device functions
 */
__device__ unsigned int cuda_rand( unsigned int *m_z, unsigned int *m_w )
{
	*m_z = 36969 * ( (*m_z) & 65535 ) + ( (*m_z) >> 16 );
	*m_w = 18000 * ( (*m_w) & 65535 ) + ( (*m_w) >> 16 );

	return ( ( (*m_z) << 16 ) + (*m_w) ) % CUDA_RAND_MAX;
}

/**
 * Calculate a set of dim_in - 1 axes perpendicular to vec_in
 */
__device__ void local_axes( int dim_in, double *vec_in, double *basis_out )
{
	int i,j,k;
	unsigned int mz,mw;
	double sum;

	sum = 0.0;
	for(i=0;i<dim_in;i++)
		basis_out[i] = vec_in[i], sum += basis_out[i] * basis_out[i];
	sum = sqrt( sum );
	for(i=0;i<dim_in;i++)
		basis_out[i] /= sum;
	mz = 150; mw = 40;
	for(i=1;i<dim_in;i++)
	{
		/* This is a horrible idea */
		for(j=0;j<dim_in;j++)
			basis_out[i*dim_in+j] = 0.5 - ( (double) cuda_rand( &mz, &mw ) / (double) CUDA_RAND_MAX );
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
__device__ int sphere_intersection( int dim_in, double *ctr1_in, double rad1_in, double *ctr2_in, double rad2_in, double *ctr_out, double *rad_out, double *qbox_out )
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

/**
 * Generate a quadrature point in the lens of spherical intersection
 */
__device__ void lens_gauss_point( int dim_in,
		double *ctr1_in, double rad1_in,
		double *ctr2_in, double rad2_in,
		double cr_in, double *nqbox_in, long *index_in,
		double *qpts_in, double *qwts_in,
		double *qp_out, double *qw_out )
{
	int i,j;
	double ssum,x1,x2;
	double *dr,*vec,*wec;

	/* Must allocate memory via malloc in __device__ code */
	dr = (double*) malloc( ( dim_in - 1 ) * sizeof(double) );
	vec = (double*) malloc( dim_in * sizeof(double) );
	wec = (double*) malloc( dim_in * sizeof(double) );

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
		*qw_out *= qwts_in[index_in[i-1]] * 2.0 * dr[i-1];
	*qw_out *= fabs( x2 - x1 );
}

__device__ void sphere_gauss_point( int dim_in, double *ctr_in, double rad_in, double *nqbox_in, long *index_in, double *qpts_in, double *qwts_in, double *qp_out, double *qw_out )
{
	int i,j;
	double *dr,*vec;

	/* Allocate directly */
	dr = (double*) malloc( dim_in * sizeof(double) );
	vec = (double*) malloc( dim_in * sizeof(double) );

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
		(*qw_out) *= 2.0 * dr[i] * qwts_in[index_in[i]];
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define PUNITY_NEIGHBOR_INC 512
#define PUNITY_NB_INCREMENT 512

typedef struct
{
	/**
	 * Dimension of the point space
	 */
	int dim;

	/**
	 * The number of window functions
	 */
	int npts;

	/**
	 * List of the actual points
	 */
	double *pts;

	/**
	 * List of dilation factors for each point
	 */
	double *dlt;

	/**
	 * Window function pointer
	 */
	double (*wfs)(int,double,double*);

	/**
	 * Pointer to function returning gradient of
	 * the window function
	 */
	double (*wfsd)(int,int,double,double*);

	/**
	 * Maximum value of dilation factor
	 */
	double rmax;

	/**
	 * The kdtree storing the point list for
	 * rapid access
	 */
	void *kdt;

	/**
	 * List of ones and zeros; if bdry[i] = 1
	 * then function i centered at point i is
	 * on the boundary
	 */
	char *bdry;
} punity_t;

__device__ double cubic_window( int dim_in, double a_in, double *x_in )
{
        int i;
        double z = 0.0;
#ifdef PUNITY_NORMALIZE
	double vol = pow( a_in, (double) dim_in );
#else
	double vol = 1.0;
#endif
        for(i=0;i<dim_in;i++)
                z += x_in[i] * x_in[i];
        z = sqrt( z ) / a_in * 2.0;

        if( z > 2.0 )
                return 0.0 / vol;
        if( z > 1.0 )
                return ( 2.0 - z ) * ( 2.0 - z ) * ( 2.0 - z ) / 6.0 / vol;
        if( z >= 0.0 )
                return ( 4.0 - 6.0 * z * z + 3.0 * z * z * z ) / 6.0 / vol;
}

__device__ double cubic_window_deriv( int dim_in, int drv_in, double a_in, double *x_in )
{
	int i,j;
	double sum,prd,z = 0.0;

	/* Initialize polynomial coefficients correctly pre-scaled */
	double cf1[4] = { 4.0 / 3.0, 2.0 * -2.0 / a_in, 4.0 * 1.0 / a_in / a_in, 8.0 * -1.0 / 6.0 / a_in / a_in / a_in };
	double cf2[4] = { 2.0 / 3.0, 2.0 * 0.0 / a_in, 4.0 * -1.0 / a_in / a_in, 8.0 * 1.0 / 2.0 / a_in / a_in / a_in };

	/* Calculate distance from zero */
	for(i=0;i<dim_in;i++)
		z += x_in[i] * x_in[i];
        z = sqrt( z );

#ifdef PUNITY_NORMALIZE
	double vol = pow( a_in, (double) dim_in );
#else
	double vol = 1.0;
#endif

	/* Evaluate the drv_in derivative for each case */
	sum = 0.0;
	if( z / a_in > 1.0 )
		return 0.0 / vol;
	if( z / a_in > 0.5 )
	{
		/* Differentiate each term and evaluate at z */
		for(i=0;i<4;i++) /* Variable i doubles as the order of the current term */
		{
			prd = cf1[i];
			for(j=0;j<drv_in;j++)
				prd *= (double) ( i - j );
			prd *= pow( z, (double) ( i - drv_in ) );
			sum += prd;
		}
		return sum / vol;
	}
	if( z / a_in >= 0.0 )
	{
		/* Differentiate each term and evaluate at z */
                for(i=0;i<4;i++) /* Variable i doubles as the order of the current term */
                {
                        prd = cf2[i];
                        for(j=0;j<drv_in;j++)
                                prd *= (double) ( i - j );
                        prd *= pow( z, (double) ( i - drv_in ) );
                        sum += prd;
                }
		return sum / vol;
	}
}

__device__ double quartic_window( int dim_in, double a_in, double *x_in )
{
	int i;
	double z = 0.0;
	for(i=0;i<dim_in;i++)
		z += x_in[i] * x_in[i];
	z = sqrt( z ) / a_in;

#ifdef PUNITY_NORMALIZE
	double vol = pow( a_in, (double) dim_in );
#else
	double vol = 1.0;
#endif

	if( z > 1.0 )
		return 0.0 / vol;
	else
		return ( 1.0 - 6.0 * z * z + 8.0 * z * z * z - 3.0 * z * z * z * z ) / vol;
}

__device__ double quartic_window_deriv( int dim_in, int drv_in, double a_in, double *x_in )
{
	int i;

	
}

/**
 * Evaluate the window function; IMPORTANT: This function shifts the window
 * functions and places it at the points stored in obj_in->pts + idx_in points;
 * obj_in->wfs and obj_in->wfsd do not translation, but they do dilate
 * @param obj_in PU object
 * @param idx_in Window index
 * @param x_in Point at which to evaluate
 */
__device__ double punity_window_evaluate( punity_t *obj_in, int idx_in, double *x_in )
{
	int i;
	double y;
	double *x = (double*) malloc( obj_in->dim * sizeof(double) );

	for(i=0;i<obj_in->dim;i++)
		x[i] = x_in[i] - obj_in->pts[idx_in*obj_in->dim+i];
	y = obj_in->wfs( obj_in->dim, obj_in->dlt[idx_in], x );
	free( x );
	return y;
}

/**
 * Evaluate the window function placing a singularity at the
 * origin of the function for generating a partition of unity
 * with the Kronecker delta property; IMPORTANT: This function
 * also translates the origin to the point in obj_in->pts + idx_in
 * @param obj_in PU object
 * @param idx_in Window index
 * @param x_in Point at which to evaluate
 * @param exp_in Exponent to define the singularity
 */
__device__ double punity_window_evaluate_delta( punity_t *obj_in, int idx_in, double *x_in, int exp_in )
{
	int i;
	double r = 0.0;

	for(i=0;i<obj_in->dim;i++)
		r += pow( x_in[i] - obj_in->pts[idx_in*obj_in->dim+i], 2.0 );
	r = sqrt( r );

	return punity_window_evaluate( obj_in, idx_in, x_in ) / pow( r / obj_in->dlt[idx_in], (double) exp_in );
}

/**
 * Initialize the internal data structures for storage of points
 * and particle radii.
 * @param obj_in The partition of unity data structure
 * @param dim_in The dimension of the point space
 * @param npts_in Number of points to use to generate the partition
 * @param pts_in Set of points to use to generate the partition
 * @param dlt_in Set of dilation or scale factors or radii of particles
 * @param wfs_in Window function to use to generate the partition
 * @param wfsd_in Function returning specific derivatives of the window
 * @return Returns 0 if no error, -1 if memory issues
 */
__device__ int punity_init( punity_t *obj_in, int dim_in, int npts_in, double *pts_in, double *dlt_in, double (*wfs_in)(int,double,double*), double (*wfsd_in)(int,int,double,double*) )
{
	int i;

	/* Set all the basic dimensional information */
	obj_in->dim = dim_in;
	obj_in->npts = npts_in;
	obj_in->pts = (double*) malloc( dim_in * npts_in * sizeof(double) );
	obj_in->dlt = (double*) malloc( npts_in * sizeof(double) );
	if( obj_in->pts == NULL || obj_in->dlt == NULL )
		return -1;

	/* Copy the point data into the structure */
	for(i=0;i<dim_in*npts_in;i++)
		obj_in->pts[i] = pts_in[i];
	for(i=0;i<npts_in;i++)
		obj_in->dlt[i] = dlt_in[i];

	/* Set the function pointers */
	obj_in->wfs = wfs_in;
	obj_in->wfsd = wfsd_in;

	/* Calculate the maximum dilation factor in the system */
	for(i=0;i<npts_in;i++)
		if( i == 0 || dlt_in[i] > obj_in->rmax )
			obj_in->rmax = dlt_in[i];

	/* Initialize all points to internal points; set boundary points to 1 later */
	obj_in->bdry = (char*) malloc( npts_in * sizeof(char) );
	for(i=0;i<npts_in;i++)
		obj_in->bdry[i] = 0;

	return 0;
}

/**
 * Clean up memory allocated for punity_t
 * @param obj_in PU object to clean up
 */
__device__ int punity_free( punity_t *obj_in )
{
	free( obj_in->pts );
	free( obj_in->dlt );
	free( obj_in->bdry );
}

/**
 * Evaluate the particle function idx_in at point x_in.
 * @param obj_in Partition of unity object
 * @param idx_in Index of the function to evaluate
 * @param x_in Position at which to evaluate the function
 * @return Value of function idx_in at point x_in
 */
__device__ double punity_evaluate( punity_t *obj_in, int idx_in, double *x_in )
{
	int i,j;
	double res,sum,y;
	double *vec = (double*) malloc( obj_in->dim * sizeof(double) );

	/* If an index is given outside the range of number of points then return zero */
	if( idx_in > obj_in->npts - 1 )
		return 0.0;

	/* Don't bother with the denominator if the numerator is zero */
	sum = 0.0;
	for(i=0;i<obj_in->dim;i++)
		sum += pow( obj_in->pts[idx_in*obj_in->dim+i] - x_in[i], 2.0 );
	if( sum > obj_in->dlt[idx_in] * obj_in->dlt[idx_in] )
		return 0.0;

	/* NOTE: Change this to punity_neighbors_fast() */
	sum = 0.0;
	for(i=0;i<obj_in->npts;i++)
	{
		res = 0.0;
		for(j=0;j<obj_in->dim;j++)
			vec[j] = x_in[j] - obj_in->pts[i*obj_in->dim+j], res += vec[j] * vec[j];
		if( res < obj_in->dlt[i] * obj_in->dlt[i] )
			sum += obj_in->wfs( obj_in->dim, obj_in->dlt[i], vec );
	}
	for(i=0;i<obj_in->dim;i++)
		vec[i] = x_in[i] - obj_in->pts[idx_in*(obj_in->dim)+i];
	y = obj_in->wfs( obj_in->dim, obj_in->dlt[idx_in], vec ) / sum;
	free( vec );
	return y;
}

/**
 * Evaluate the partition of unity with singularities if needed; this function
 * reads the values of obj_in->bdry for all functions in the domain in order
 * to decide which windows to evaluate with singularities; this is mainly used
 * to make boundary conditions easier to implement, but it can be used to implement
 * interpolating partitions of unity in general (although there may be degradation
 * of interpolants); in the future, make the exponent a function of the function
 * index, i.e. store as obj_in->exp[i]
 * @param obj_in The PU object
 * @param idx_in Index of the function to evaluate
 * @param x_in Point at which to evaluate the functions
 * @param exp_in Exponent to use for singularity evaluation
 */
__device__ double punity_evaluate_delta( punity_t *obj_in, int idx_in, double *x_in, int exp_in )
{
	int i,j;
	double res,sum,y;
	double *vec = (double*) malloc( obj_in->dim * sizeof(double) );

	/* If an index is given outside the range of number of points then return zero */
	if( idx_in > obj_in->npts - 1 )
		return 0.0;

	/* Don't bother with the denominator if the numerator is zero */
	sum = 0.0;
	for(i=0;i<obj_in->dim;i++)
		sum += pow( obj_in->pts[idx_in*obj_in->dim+i] - x_in[i], 2.0 );
	if( sum > obj_in->dlt[idx_in] * obj_in->dlt[idx_in] )
		return 0.0;

	/* NOTE: Change this to punity_neighbors_fast() */
	sum = 0.0;
	for(i=0;i<obj_in->npts;i++)
	{
		res = 0.0;
		for(j=0;j<obj_in->dim;j++)
			vec[j] = x_in[j] - obj_in->pts[i*obj_in->dim+j], res += vec[j] * vec[j];
		res = sqrt( res );
		if( res < obj_in->dlt[i] )
		{
			if( obj_in->bdry[i] == 0 )
				sum += obj_in->wfs( obj_in->dim, obj_in->dlt[i], vec );
			else
				sum += obj_in->wfs( obj_in->dim, obj_in->dlt[i], vec ) / pow( res / obj_in->dlt[i], (double) exp_in );
		}
	}
	res = 0.0;
	for(i=0;i<obj_in->dim;i++)
		vec[i] = x_in[i] - obj_in->pts[idx_in*obj_in->dim+i], res += vec[i] * vec[i];
	res = sqrt( res );
	if( res > obj_in->dlt[idx_in] )
		return 0.0;
	else
	{
		if( obj_in->bdry[idx_in] == 0 )
			y = obj_in->wfs( obj_in->dim, obj_in->dlt[idx_in], vec ) / sum;
		else
			y = obj_in->wfs( obj_in->dim, obj_in->dlt[idx_in], vec ) / pow( res / obj_in->dlt[idx_in], (double) exp_in ) / sum;
	}
	free( vec );
	return y;
}

/**
 * Evaluate the derivative of r = | xi - xj | w.r.t. drv_in; FIXME: This function
 * returns values which become singular as r -> 0; figure out a way to get around this
 * @param dim_in Dimension of the vector space
 * @param drv_in Derivative integer vector
 * @param x_in Point at which to evaluate the derivative
 */
__device__ double radial_deriv_evaluate( int dim_in, int *drv_in, double *x_in )
{
	int b,c,i,j,p,n,*v,*w,*s,*m;
	double d,prd,sum = 0.0;

	/* Set up the vector to partition */
	n = 0;
	for(i=0;i<dim_in;i++)
		n += drv_in[i];
	v = (int*) malloc( n * sizeof(int) );
	w = (int*) malloc( n * sizeof(int) );
	s = (int*) malloc( n * sizeof(int) );
	m = (int*) malloc( n * sizeof(int) );
	for(i=0,p=0;i<dim_in;i++)
		for(j=0;j<drv_in[i];j++)
			v[p++] = i;

	/* Start iterating through the partitions */
	partition_init( s, m, n );
	do
	{
		/* Variable b contains the number of blocks in this partition */
		for(i=0;i<n;i++)
			if( i == 0 || s[i] > b )
				b = s[i];

		/* Generate each of the blocks from the membership vector s */
		prd = 1.0;
		for(i=0;i<b;i++) /* Index i is the current block */
		{
			/* The number of entries in block i is counted in c */
			for(j=0,c=0;j<n;j++)
				if( s[j] == i + 1 )
					w[c++] = v[j]; /* Take derivative of u = x1^2 + ... + xd^2 w.r.t. coordinate v[j] */

			/* Do the derivative */
			if( c > 2 )
				d = 0.0;
			else
			{
				if( c == 2 )
				{
					if( w[0] != w[1] )
						d = 0.0;
					else
						d = 2.0;
				}
				else if( c == 1 )
					d = 2.0 * x_in[w[0]];
				else /* The zero derivative should not occur in this sequence, but anyway... */
				{
					d = 0.0;
					for(j=0;j<dim_in;j++)
						d += x_in[j] * x_in[j];
				}
			}

			/* Multiply the derivative in d for this block into the total product */
			prd *= d; /* Variable d is derivative of u consistent with this block */
		}

		/* Calculate b derivative of radial function */
		for(i=0;i<b;i++)
			prd *= ( 0.5 - (double) i );
		d = 0.0;
		for(i=0;i<dim_in;i++)
			d += x_in[i] * x_in[i];
		prd *= pow( d, 0.5 - (double) b );
		sum += prd;
	}
	while( partition_next( s, m, n ) != 0 );

	free( v );
	free( w );
	free( s );
	free( m );

	return sum;
}

/**
 * Evaluates the derivatives of r**-p to any order in any dimension
 * using radial_deriv_evaluate
 * @param dim_in Dimension of vector space
 * @param drv_in The derivative integer vector
 * @param x_in Point at which to evaluate derivative
 * @param a_in Dilation factor to use; r -> (r/a)**-p
 * @param exp_in Exponent p to use
 */
__device__ double rational_radial_deriv_evaluate( int dim_in, int *drv_in, double *x_in, double a_in, int exp_in )
{
	int b,c,i,j,k,p,n,*v,*w,*s,*m,*q;
	double d,y,prd,sum = 0.0;

	/* Set up the vector to partition */
	n = 0;
	for(i=0;i<dim_in;i++)
		n += drv_in[i];
	v = (int*) malloc( n * sizeof(int) );
	w = (int*) malloc( n * sizeof(int) );
	s = (int*) malloc( n * sizeof(int) );
	m = (int*) malloc( n * sizeof(int) );
	q = (int*) malloc( dim_in * sizeof(int) );
	for(i=0,p=0;i<dim_in;i++)
		for(j=0;j<drv_in[i];j++)
			v[p++] = i;

	/* Form radial value */
	d = 0.0;
	for(i=0;i<dim_in;i++)
		d += x_in[i] * x_in[i];
	d = sqrt( d );

	/* Return now if derivative order is zero */
	if( n == 0 )
	{
		free( v );
		free( w );
		free( s );
		free( m );
		free( q );
		return pow( d / a_in, (double) ( -exp_in ) );
	}

	/* Start iterating through the partitions */
	partition_init( s, m, n );
	do
	{
		/* Variable b contains the number of blocks in this partition */
		for(i=0;i<n;i++)
			if( i == 0 || s[i] > b )
				b = s[i];

		/* Take the b-th derivative of r**-p w.r.t. r because b is the number of blocks in partition s */
		prd = 1.0;
		for(i=0;i<b;i++)
			prd *= (double) ( -exp_in - i );
		prd *= pow( d, (double) ( -exp_in - b ) );

		/* Now deal with all the other block derivatives of r w.r.t. x's */
		for(i=0;i<b;i++)
		{
			/* The number of entries in block i is counted in c */
			for(j=0,c=0;j<n;j++)
				if( s[j] == i + 1 )
					w[c++] = v[j]; /* Take derivative of u = x1^2 + ... + xd^2 w.r.t. coordinate v[j] */

			/* Now rebuild the derivative in drv_in format to pass to radial_deriv_evaluate */
			for(j=0;j<dim_in;j++)
				q[j] = 0;
			for(j=0;j<c;j++)
				q[w[j]]++; /* Everytime an index appears in w, increment its component once */
			y = radial_deriv_evaluate( dim_in, q, x_in );
			prd *= y;
		}
		sum += prd;
	}
	while( partition_next( s, m, n ) != 0 );

	free( v );
	free( w );
	free( s );
	free( m );
	free( q );

	return pow( a_in, (double) exp_in ) * sum;
}

/**
 * Evaluate the derivative of phi_i / ( sum_j phi_j ) w.r.t. the phi_j's themselves
 * where each phi_j is evaluated at the input point x_in
 */
__device__ double punity_comp_deriv_evaluate( punity_t *obj_in, int idx_in, int *drv_in, int q_in, int *nb_in, double *x_in )
{
	int i,j,k,n;
	double sum,prd1,prd2;

	/* Build the rank of the derivative w.r.t. all variables not equal to idx_in */
	n = 0;
	for(i=0;i<q_in;i++)
		n += drv_in[i];

	/* Which position in drv_in corresponds to idx_in */
	for(i=0,k=-1;i<q_in;i++)
		if( nb_in[i] == idx_in )
			k = i;
	if( k == -1 ) /* Zero/error because idx_in is not a neighbor to x_in */
		return 0.0;

	/* Calculate sum of all window functions */
	sum = 0.0;
        for(i=0;i<q_in;i++)
                sum += punity_window_evaluate( obj_in, nb_in[i], x_in );

	/* Simplified derivative only has two terms! */
	prd1 = punity_window_evaluate( obj_in, idx_in, x_in ) * ( n % 2 == 0 ? 1.0 : -1.0 )
		* (double) factorial( n ) * pow( sum, -1.0 * (double) ( n + 1 ) );

	/* If drv_in[k] == 0 then return prd1 as output */
	if( drv_in[k] == 0 )
		return prd1;

	/* Build the second term in the sum */
	prd2 = (double) drv_in[k] * ( ( n - 1 ) % 2 == 0 ? 1.0 : -1.0 ) * (double) factorial( n - 1 ) * pow( sum, -1.0 * (double) n );

	/* Form sum */
	return prd1 + prd2;
}

/**
 * This functions makes use of the Faa di Bruno formula for higher derivatives
 * of a composition of an r-dependent function with r as a function of
 * the individual x variables; r = sqrt( x1^2 + ... + xd^2 )
 */
__device__ double punity_window_deriv_evaluate( punity_t *obj_in, int idx_in, int *drv_in, double *x_in )
{
	int b,c,i,j,p,n,*v,*w,*s,*m,*t;
        double d,prd,sum = 0.0;

	/* Set up the vector to partition */
        n = 0;
        for(i=0;i<obj_in->dim;i++)
                n += drv_in[i];
        v = (int*) malloc( n * sizeof(int) );
        w = (int*) malloc( n * sizeof(int) );
        s = (int*) malloc( n * sizeof(int) );
        m = (int*) malloc( n * sizeof(int) );
	t = (int*) malloc( obj_in->dim * sizeof(int) );
        for(i=0,p=0;i<obj_in->dim;i++)
                for(j=0;j<drv_in[i];j++)
                        v[p++] = i;

	/* If derivative order is zero then return the function undifferentiated */
	if( n == 0 )
	{
		free( v );
		free( w );
		free( s );
		free( m );
		free( t );
		return obj_in->wfs( obj_in->dim, obj_in->dlt[idx_in], x_in );
	}

	/* Start iterating through the partitions */
        partition_init( s, m, n );
        do
        {
                /* Variable b contains the number of blocks in this partition */
                for(i=0;i<n;i++)
                        if( i == 0 || s[i] > b )
                                b = s[i];

                /* Generate each of the blocks from the membership vector s */
                prd = 1.0;
                for(i=0;i<b;i++) /* Index i is the current block */
                {
                        /* The number of entries in block i is counted in c */
                        for(j=0,c=0;j<n;j++)
                                if( s[j] == i + 1 )
                                        w[c++] = v[j]; /* Take derivative of u = x1^2 + ... + xd^2 w.r.t. coordinate v[j] */

			/* Build the derivative in terms of exponents */
			for(j=0;j<obj_in->dim;j++)
				t[j] = 0;
			for(j=0;j<c;j++)
				t[w[j]]++; /* Everytime an index appears in w, increment its component once */

                        /* Do the derivative */
			d = radial_deriv_evaluate( obj_in->dim, t, x_in ); /* No dilation factor here; all contained in wfsd */

                        /* Multiply the derivative in d for this block into the total product */
                        prd *= d; /* Variable d is derivative of u consistent with this block */
                }

		/* Calculate the b order derivative of the radial polynomial */
		prd *= obj_in->wfsd( obj_in->dim, b, obj_in->dlt[idx_in], x_in );

		/* Add the contribution from this partition to the total in sum */
                sum += prd;
	}
	while( partition_next( s, m, n ) != 0 );

	free( v );
	free( w );
	free( s );
	free( m );
	free( t );

	return sum;
}

/**
 * Evaluate arbitrary derivatives of partition of unity functions with
 * a singularity of order exp_in at the origin of the function
 */
__device__ double punity_window_deriv_evaluate_delta( punity_t *obj_in, int idx_in, int *drv_in, double *x_in, int exp_in )
{
	int i,j,k,m,n,p,q,r,*qc,*qd,*v,*cmb,*dmb;
	double d,prd,sum = 0.0;

	/* Do some setup */
	n = 0;
	for(i=0;i<obj_in->dim;i++)
		n += drv_in[i];
	v = (int*) malloc( n * sizeof(int) );
	qc = (int*) malloc( obj_in->dim * sizeof(int) );
	qd = (int*) malloc( obj_in->dim * sizeof(int) ); /* Complement of q w.r.t. the current cmb state */
	cmb = (int*) malloc( n * sizeof(int) );
	dmb = (int*) malloc( n * sizeof(int) );

	/* Calculate the derivative vector; a sequence of the numbers 0,...,dim-1 */
	for(i=0,p=0;i<obj_in->dim;i++)
		for(j=0;j<drv_in[i];j++)
			v[p++] = i;

	/* Calculate the radial vector */
	d = 0.0;
	for(i=0;i<obj_in->dim;i++)
		d += x_in[i] * x_in[i];
	d = sqrt( d );

	/* If derivative order is zero then return the function undifferentiated */
	if( n == 0 )
	{
		free( v );
		free( qc );
		free( qd );
		free( cmb );
		free( dmb );
		return obj_in->wfs( obj_in->dim, obj_in->dlt[idx_in], x_in ) / pow( d / obj_in->dlt[idx_in], (double) exp_in );
	}

	/* Iterate through all combinations sizes and for each size, iterate through all combinations */
	for(i=0;i<n;i++)
	{
		combinadic_init( n, i, cmb );
		k = binomial( n, i );
		for(j=0;j<k;j++)
		{
			/* Form the derivative vector for this case */
			for(m=0;m<obj_in->dim;m++)
				qc[m] = 0;
			for(m=0;m<i;m++)
				qc[v[cmb[m]]]++;

			/* Build the complement */
			if( i > 0 )
			{
				for(m=0;m<cmb[0];m++)
					dmb[m] = m;
				for(m=0,r=cmb[0];m<i-1;m++) /* Iterate consecutive integer pairs */
					for(q=cmb[m]+1;q<cmb[m+1];q++)
						dmb[r++] = q;
				for(m=cmb[i-1]+1;m<n;m++)
					dmb[r++] = m;
			}
			else
				for(m=0;m<n;m++)
					dmb[m] = m;

			/* Form the derivative of the window */
			for(m=0;m<obj_in->dim;m++)
				qd[m] = 0;
			for(m=0;m<n-i;m++)
				qd[v[dmb[m]]]++;

			/* Apply the derivative to the window function and evaluate it */
			prd = 1.0;
			if( i == 0 )
				prd *= pow( d, (double) ( -exp_in ) );
			else
				prd *= rational_radial_deriv_evaluate( obj_in->dim, qc, x_in, obj_in->dlt[idx_in], exp_in );
			if( n - i == 0 )
				prd *= obj_in->wfs( obj_in->dim, obj_in->dlt[idx_in], x_in );
			else
				prd *= punity_window_deriv_evaluate( obj_in, idx_in, qd, x_in );

			/* Add this shit up */
			sum += prd;

			/* Take a step to the next combinadic vector */
			combinadic_next( n, i, cmb );
		}
	}

	free( v );
	free( qc );
	free( qd );
	free( cmb );
	free( dmb );

	return sum;
}

__device__ int punity_neighbors( punity_t *obj_in, double *x_in, int **nb_out )
{
	int i,j,q,na,*nb,*nc;
	double prd;

	na = PUNITY_NB_INCREMENT;
        nb = (int*) malloc( na * sizeof(int) );
        for(i=0,q=0;i<obj_in->npts;i++)
        { 
                prd = 0.0;
                for(j=0;j<obj_in->dim;j++)
                        prd += pow( obj_in->pts[i*obj_in->dim+j] - x_in[j], 2.0 );
                prd = sqrt( prd );
                if( prd < obj_in->dlt[i] ) /* Then it is close enough to contribute */
                {
                        /* See if we need to allocate more space */
                        if( q + 1 > na )
                        {
                                na += PUNITY_NB_INCREMENT;
				nc = (int*) malloc( na * sizeof(int) );
				for(j=0;j<q;j++)
					nc[j] = nb[j];
				free( nb );
				nb = nc; /* Swap definitions to simulate realloc */
                        }
                        nb[q++] = i;
                }
        }
	*nb_out = nb;

	return q;
}

/**
 * Evaluate the derivative of a partition-of-unity function with given
 * index and given components
 * @param obj_in Partition of unity object
 * @param idx_in Index of the PU functions to evaluate (index to the point)
 * @param drv_in Derivative specified as the order w.r.t. each independent variable
 * @param x_in Point at which to evaluate the resulting derivative
 * @param nb_in Neighbors of point idx_in; if NULL it will be calculated
 * @return The value of the derivative at point x_in
 */
__device__ double punity_term_deriv_evaluate( punity_t *obj_in, int idx_in, int *drv_in, double *x_in, int *nb_in, int q_in )
{
	int b,c,i,j,k,p,q,n,bf,*v,*w,*s,*m,*t,*u,*nb;
	long *size,*index;
        double d,prd,tsum,*y,sum = 0.0;

	/* First check to see if x_in is too far from obj_in->pts + idx_in * dim */
	prd = 0.0;
	for(i=0;i<obj_in->dim;i++)
		prd += pow( obj_in->pts[idx_in*obj_in->dim+i] - x_in[i], 2.0 );
	prd = sqrt( prd );
	if( prd > obj_in->dlt[idx_in] )
		return 0.0;

	/* Temporary position vector for translation */
	y = (double*) malloc( obj_in->dim * sizeof(double) );

	/* Set up the vector to partition */
        n = 0;
        for(i=0;i<obj_in->dim;i++)
                n += drv_in[i];
        v = (int*) malloc( n * sizeof(int) );
        w = (int*) malloc( n * sizeof(int) );
        s = (int*) malloc( n * sizeof(int) );
        m = (int*) malloc( n * sizeof(int) );
	t = (int*) malloc( obj_in->dim * sizeof(int) );
        for(i=0,p=0;i<obj_in->dim;i++)
                for(j=0;j<drv_in[i];j++)
                        v[p++] = i;

	/* Let q count the number of functions contributing to the total particle function idx_in */
#ifdef PUNITY_USE_KDTREES
	if( nb_in == NULL )
		q = punity_neighbors_fast( obj_in, x_in, &nb ), bf = 1; /* Free nb at the end */
	else
		q = q_in, nb = nb_in, bf = 0; /* Don't free at the end */
#else
	if( nb_in == NULL )
		q = punity_neighbors( obj_in, x_in, &nb ), bf = 1; /* Free nb at the end */
	else
		q = q_in, nb = nb_in, bf = 0; /* Don't free nb at the end */
#endif
	u = (int*) malloc( q * sizeof(int) ); /* Serves as derivative vector */

	/* Maximum number of blocks possible is n */
	size = (long*) malloc( n * sizeof(long) );
	index = (long*) malloc( n * sizeof(long) );
	for(i=0;i<n;i++)
		size[i] = q - 1; /* Maximum value is q - 1; starting at 0 makes for q values */

	/* Start iterating through the partitions */
        partition_init( s, m, n );
        do
        {
                /* Variable b contains the number of blocks in this partition, s */
                for(i=0;i<n;i++)
                        if( i == 0 || s[i] > b )
                                b = s[i];

		/* Iterate through all length-b vectors with entries in {1,...,q} */
		for(i=0;i<b;i++)
			index[i] = 0;
		do
		{
			/* Build the derivative integer vector from the index vector in this loop */
			for(i=0;i<q;i++)
				u[i] = 0;
			for(i=0;i<b;i++)
				u[index[i]]++;

			/* Calculate the derivative w.r.t. the window functions */
			prd = punity_comp_deriv_evaluate( obj_in, idx_in, u, q, nb, x_in );

	                /* Generate each of the blocks from the membership vector */
	                for(i=0;i<b;i++) /* Index i is the current block */
	                {
	                        /* The number of entries in block i is counted in c */
	                        for(j=0,c=0;j<n;j++)
	                                if( s[j] == i + 1 )
	                                        w[c++] = v[j]; /* Take derivative of function w.r.t. coordinate v[j] */
				for(j=0;j<obj_in->dim;j++)
					t[j] = 0;
				for(j=0;j<c;j++)
					t[w[j]]++;

				/* Build the product of partition derivatives specified by the indexes in index of each block */
				for(j=0;j<obj_in->dim;j++)
					y[j] = x_in[j] - obj_in->pts[nb[index[i]]*obj_in->dim+j];
				tsum = punity_window_deriv_evaluate( obj_in, nb[index[i]], t, y );
				prd *= tsum;
	                }
			sum += prd;
		}
		while( arraynext( b, size, index ) != -1 );
	}
	while( partition_next( s, m, n ) != 0 );

	free( v );
	free( w );
	free( s );
	free( m );
	free( t );
	free( u );
	free( y );
	free( size ); free( index );
	if( bf )
		free( nb );

	return sum;
}

/**
 * Evaluate the derivative of a partition-of-unity function with given
 * index and given components
 * @param obj_in Partition of unity object
 * @param idx_in Index of the PU functions to evaluate (index to the point)
 * @param drv_in Derivative specified as the order w.r.t. each independent variable
 * @param x_in Point at which to evaluate the resulting derivative
 * @param exp_in Exponent to use for singularity, if it exists
 * @param nb_in Neighbors given as input; NULL means calculate it yourself
 * @return The value of the derivative at point x_in
 */
__device__ double punity_term_deriv_evaluate_delta( punity_t *obj_in, int idx_in, int *drv_in, double *x_in, int exp_in, int *nb_in, int q_in )
{
	int b,c,i,j,k,p,q,n,bf,*v,*w,*s,*m,*t,*u,*nb;
	long *size,*index;
        double d,prd,tsum,*y,sum = 0.0;

	/* First check to see if x_in is too far from obj_in->pts + idx_in * dim */
	prd = 0.0;
	for(i=0;i<obj_in->dim;i++)
		prd += pow( obj_in->pts[idx_in*obj_in->dim+i] - x_in[i], 2.0 );
	prd = sqrt( prd );
	if( prd > obj_in->dlt[idx_in] )
		return 0.0;

	/* Temporary position vector for translation */
	y = (double*) malloc( obj_in->dim * sizeof(double) );

	/* Set up the vector to partition */
        n = 0;
        for(i=0;i<obj_in->dim;i++)
                n += drv_in[i];
        v = (int*) malloc( n * sizeof(int) );
        w = (int*) malloc( n * sizeof(int) );
        s = (int*) malloc( n * sizeof(int) );
        m = (int*) malloc( n * sizeof(int) );
	t = (int*) malloc( obj_in->dim * sizeof(int) );
        for(i=0,p=0;i<obj_in->dim;i++)
                for(j=0;j<drv_in[i];j++)
                        v[p++] = i;

	/* Let q count the number of functions contributing to the total particle function idx_in */
#ifdef PUNITY_USE_KDTREES
	if( nb_in == NULL )
		q = punity_neighbors_fast( obj_in, x_in, &nb ), bf = 1;
	else
		q = q_in, nb = nb_in, bf = 0;
#else
	if( nb_in == NULL )
		q = punity_neighbors( obj_in, x_in, &nb ), bf = 1;
	else
		q = q_in, nb = nb_in, bf = 0;
#endif

	u = (int*) malloc( q * sizeof(int) ); /* Serves as derivative vector */

	/* Maximum number of blocks possible is n */
	size = (long*) malloc( n * sizeof(long) );
	index = (long*) malloc( n * sizeof(long) );
	for(i=0;i<n;i++)
		size[i] = q - 1; /* Maximum value is q - 1; starting at 0 makes for q values */

	/* Start iterating through the partitions */
        partition_init( s, m, n );
        do
        {
                /* Variable b contains the number of blocks in this partition, s */
                for(i=0;i<n;i++)
                        if( i == 0 || s[i] > b )
                                b = s[i];

		/* Iterate through all length-b vectors with entries in {1,...,q} */
		for(i=0;i<b;i++)
			index[i] = 0;
		do
		{
			/* Build the derivative integer vector from the index vector in this loop */
			for(i=0;i<q;i++)
				u[i] = 0;
			for(i=0;i<b;i++)
				u[index[i]]++;

			/* Calculate the derivative w.r.t. the window functions */
			prd = punity_comp_deriv_evaluate( obj_in, idx_in, u, q, nb, x_in );

	                /* Generate each of the blocks from the membership vector */
	                for(i=0;i<b;i++) /* Index i is the current block */
	                {
	                        /* The number of entries in block i is counted in c */
	                        for(j=0,c=0;j<n;j++)
	                                if( s[j] == i + 1 )
	                                        w[c++] = v[j]; /* Take derivative of function w.r.t. coordinate v[j] */
				for(j=0;j<obj_in->dim;j++)
					t[j] = 0;
				for(j=0;j<c;j++)
					t[w[j]]++;

				/* Build the product of partition derivatives specified by the indexes in index of each block */
				for(j=0;j<obj_in->dim;j++)
					y[j] = x_in[j] - obj_in->pts[nb[index[i]]*obj_in->dim+j];
				if( obj_in->bdry[nb[index[i]]] == 0 )
					tsum = punity_window_deriv_evaluate( obj_in, nb[index[i]], t, y );
				else
					tsum = punity_window_deriv_evaluate_delta( obj_in, nb[index[i]], t, y, exp_in );
				prd *= tsum;
	                }
			sum += prd;
		}
		while( arraynext( b, size, index ) != -1 );
	}
	while( partition_next( s, m, n ) != 0 );

	free( v );
	free( w );
	free( s );
	free( m );
	free( t );
	free( u );
	free( y );
	free( size ); free( index );
	if( bf )
		free( nb );

	return sum;
}

/**
 * This function is intended to be a faster version of punity_term_deriv_evaluate
 * for use calculating the first derivatives only
 */
__device__ double punity_term_first_deriv_evaluate( punity_t *obj_in, int idx_in, int drv_in, double *x_in )
{
	
}

/**
 * Evaluate particle-localized polynomial term
 */
__device__ double punity_eval_local_poly( punity_t *obj_in, int idx_in, int *pidx_in, double *x_in )
{
	int i;
	double prd = 1.0;

	/* Calculate the monomial term evaluated at x_in */
	for(i=0;i<obj_in->dim;i++) /* Evaluate everything on the local scale as well */
		prd *= pow( ( x_in[i] - obj_in->pts[idx_in*obj_in->dim+i] ) / obj_in->dlt[idx_in], (double) pidx_in[i] );

	return prd * punity_evaluate( obj_in, idx_in, x_in );
}

/**
 * Evaluate particle-localizaed polynomial term while observing
 * singular particle terms
 */
__device__ double punity_eval_local_poly_delta( punity_t *obj_in, int idx_in, int *pidx_in, double *x_in, int exp_in )
{
	int i;
	double prd = 1.0;

	/* Calculate the monomial term evaluated at x_in */
	for(i=0;i<obj_in->dim;i++)
		prd *= pow( ( x_in[i] - obj_in->pts[idx_in*obj_in->dim+i] ) / obj_in->dlt[idx_in], (double) pidx_in[i] );

	return prd * punity_evaluate_delta( obj_in, idx_in, x_in, exp_in );
}

/**
 * Evaluate derivatives of the particle-localized polynomial term
 * while ignoring singular particles
 */
__device__ double punity_term_deriv_local_poly( punity_t * obj_in, int idx_in, int *pidx_in, int *drv_in, double *x_in, int *nb_in, int q_in )
{
	int i,j,k,m,n,p,q,r,*qc,*qd,*v,*cmb,*dmb;
	double prd,coeff,sum = 0.0;

	/* Do some setup */
	n = 0;
	for(i=0;i<obj_in->dim;i++)
		n += drv_in[i];
	v = (int*) malloc( n * sizeof(int) );
	qc = (int*) malloc( obj_in->dim * sizeof(int) );
	qd = (int*) malloc( obj_in->dim * sizeof(int) ); /* Complement of q w.r.t. the current cmb state */
	cmb = (int*) malloc( n * sizeof(int) );
	dmb = (int*) malloc( n * sizeof(int) );

	/* Calculate the derivative vector; a sequence of the numbers 0,...,dim-1 */
	for(i=0,p=0;i<obj_in->dim;i++)
		for(j=0;j<drv_in[i];j++)
			v[p++] = i;

	/* If derivative order is zero then return the function undifferentiated */
	if( n == 0 )
	{
		free( v );
		free( qc );
		free( qd );
		free( cmb );
		free( dmb );
		return punity_eval_local_poly( obj_in, idx_in, pidx_in, x_in );
	}

	/* Iterate through all combinations sizes and for each size, iterate through all combinations */
	for(i=0;i<n;i++)
	{
		combinadic_init( n, i, cmb );
		k = binomial( n, i );
		for(j=0;j<k;j++)
		{
			/* Form the derivative vector for this case */
			for(m=0;m<obj_in->dim;m++)
				qc[m] = 0;
			for(m=0;m<i;m++)
				qc[v[cmb[m]]]++;

			/* Build the complement */
			if( i > 0 )
			{
				for(m=0;m<cmb[0];m++)
					dmb[m] = m;
				for(m=0,r=cmb[0];m<i-1;m++) /* Iterate consecutive integer pairs */
					for(q=cmb[m]+1;q<cmb[m+1];q++)
						dmb[r++] = q;
				for(m=cmb[i-1]+1;m<n;m++)
					dmb[r++] = m;
			}
			else
				for(m=0;m<n;m++)
					dmb[m] = m;

			/* Form the derivative of the window */
			for(m=0;m<obj_in->dim;m++)
				qd[m] = 0;
			for(m=0;m<n-i;m++)
				qd[v[dmb[m]]]++;

			/* Apply the derivative to the window function and evaluate it */
			prd = 1.0;
			if( i == 0 )
			{
				/* Put the value of the polynomial term in pidx_in in prd */
				for(m=0;m<obj_in->dim;m++)
					prd *= pow( ( x_in[m] - obj_in->pts[idx_in*obj_in->dim+m] ) / obj_in->dlt[idx_in], (double) pidx_in[m] );
			}
			else
			{
				/* Go through and do the derivative of the polynomial */
				for(m=0;m<obj_in->dim;m++)
				{
					p = pidx_in[m] - qc[m]; /* If qc[j] > pidx_in[j] then the whole thing is zero */
					if( p < 0 )
					{
						prd = 0.0;
						break; /* Break here because if one partial is zero then the whole term is; no need to do it */
					}
					else
					{
						prd *= pow( ( x_in[m] - obj_in->pts[idx_in*obj_in->dim+m] ) / obj_in->dlt[idx_in], (double) p ) / pow( obj_in->dlt[idx_in], (double) i );
						for(p=0;p<qc[m];p++)
							prd *= (double) ( pidx_in[m] - p );
					}
				}
			}
			if( n - i == 0 )
				prd *= punity_evaluate( obj_in, idx_in, x_in );
			else /* Evaluate the derivative of the partition of unity function */
			{
				prd *= punity_term_deriv_evaluate( obj_in, idx_in, qd, x_in, nb_in, q_in );
			}

			/* Add this shit up */
			sum += prd;

			/* Take a step to the next combinadic vector */
			combinadic_next( n, i, cmb );
		}
	}

	free( v );
	free( qc );
	free( qd );
	free( cmb );
	free( dmb );

	return sum;
}

/**
 * Evaluate derivatives of the particle-localized polynomial term
 * while ignoring singular particles
 */
__device__ double punity_term_deriv_local_poly_delta( punity_t * obj_in, int idx_in, int *pidx_in, int *drv_in, double *x_in, int exp_in, int *nb_in, int q_in )
{
	int i,j,k,m,n,p,q,r,*qc,*qd,*v,*cmb,*dmb;
	double prd,coeff,sum = 0.0;

	/* Do some setup */
	n = 0;
	for(i=0;i<obj_in->dim;i++)
		n += drv_in[i];
	v = (int*) malloc( n * sizeof(int) );
	qc = (int*) malloc( obj_in->dim * sizeof(int) );
	qd = (int*) malloc( obj_in->dim * sizeof(int) ); /* Complement of q w.r.t. the current cmb state */
	cmb = (int*) malloc( n * sizeof(int) );
	dmb = (int*) malloc( n * sizeof(int) );

	/* Calculate the derivative vector; a sequence of the numbers 0,...,dim-1 */
	for(i=0,p=0;i<obj_in->dim;i++)
		for(j=0;j<drv_in[i];j++)
			v[p++] = i;

	/* If derivative order is zero then return the function undifferentiated */
	if( n == 0 )
	{
		free( v );
		free( qc );
		free( qd );
		free( cmb );
		free( dmb );
		return punity_eval_local_poly_delta( obj_in, idx_in, pidx_in, x_in, exp_in );
	}

	/* Iterate through all combinations sizes and for each size, iterate through all combinations */
	for(i=0;i<n;i++)
	{
		combinadic_init( n, i, cmb );
		k = binomial( n, i );
		for(j=0;j<k;j++)
		{
			/* Form the derivative vector for this case */
			for(m=0;m<obj_in->dim;m++)
				qc[m] = 0;
			for(m=0;m<i;m++)
				qc[v[cmb[m]]]++;

			/* Build the complement */
			if( i > 0 )
			{
				for(m=0;m<cmb[0];m++)
					dmb[m] = m;
				for(m=0,r=cmb[0];m<i-1;m++) /* Iterate consecutive integer pairs */
					for(q=cmb[m]+1;q<cmb[m+1];q++)
						dmb[r++] = q;
				for(m=cmb[i-1]+1;m<n;m++)
					dmb[r++] = m;
			}
			else
				for(m=0;m<n;m++)
					dmb[m] = m;

			/* Form the derivative of the window */
			for(m=0;m<obj_in->dim;m++)
				qd[m] = 0;
			for(m=0;m<n-i;m++)
				qd[v[dmb[m]]]++;

			/* Apply the derivative to the window function and evaluate it */
			prd = 1.0;
			if( i == 0 )
			{
				/* Put the value of the polynomial term in pidx_in in prd */
				for(m=0;m<obj_in->dim;m++)
					prd *= pow( ( x_in[m] - obj_in->pts[idx_in*obj_in->dim+m] ) / obj_in->dlt[idx_in], (double) pidx_in[m] );

			}
			else
			{
				/* Go through and do the derivative of the polynomial */
				for(m=0;m<obj_in->dim;m++)
				{
					p = pidx_in[m] - qc[m]; /* If qc[j] > pidx_in[j] then the whole thing is zero */
					if( p < 0 )
					{
						prd = 0.0;
						break; /* Break here because if one partial is zero then the whole term is; no need to do it */
					}
					else
					{
						prd *= pow( ( x_in[m] - obj_in->pts[idx_in*obj_in->dim+m] ) / obj_in->dlt[idx_in], (double) p ) / pow( obj_in->dlt[idx_in], (double) i );
						for(p=0;p<qc[m];p++)
							prd *= (double) ( pidx_in[m] - p );
					}
				}
			}
			if( n - i == 0 )
				prd *= punity_evaluate_delta( obj_in, idx_in, x_in, exp_in );
			else /* Evaluate the derivative of the partition of unity function */
				prd *= punity_term_deriv_evaluate_delta( obj_in, idx_in, qd, x_in, exp_in, nb_in, q_in );

			/* Add this shit up */
			sum += prd;

			/* Take a step to the next combinadic vector */
			combinadic_next( n, i, cmb );
		}
	}

	free( v );
	free( qc );
	free( qd );
	free( cmb );
	free( dmb );

	return sum;
}

/**
 * Kernel function to generate the overlap and Hamiltonian matrices
 */
__global__ void puksham_mbuild_cuda( punity_t *pu, int mdim, int nbp, shape_t *dm, int *nlbase, int **lbase, int *ltg, nuclei_t *nuc, int quadn, double *qpts, double *qwts,
					 long **jap, double **Ap, long **jbp, double **Bp, long cc, long c_spmat_inc,
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
	double qw,*qx;
	double rum,sum,tum;

	/* Allocate stuff */
	qbox = (double*) malloc( dim * dim * sizeof(double) );
	nqbox = (double*) malloc( dim * dim * sizeof(double) );
	ctr = (double*) malloc( dim * sizeof(double) );
	qx = (double*) malloc( dim * sizeof(double) );
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
										idx = 1 - idx;
										
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
										idx = 1 - idx;
										
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
											idx = 1 - idx;
											
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
											idx = 1 - idx;
											
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
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
	
}

