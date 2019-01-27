#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "combinadic.h"
#include "trees.h"
#include "points.h"
#include "geom.h"
#include "config.h"

#define PUNITY_NEIGHBOR_INC 512
#define PUNITY_NB_INCREMENT 512
//#define PUNITY_NORMALIZE
//#define PUNITY_USE_KDTREES

/* Typedef this so it can be used to make an array */
typedef double (*wfsd_t)(int,int,double,double*);

/**
 * The partition of unity data structure
 */
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
	kdnode_t *kdt;

	/**
	 * List of ones and zeros; if bdry[i] = 1
	 * then function i centered at point i is
	 * on the boundary
	 */
	char *bdry;

	/**
	 * Adding in rectangular periodic boundary
	 * mode; when pbc == 1, assume the domain
	 * is periodic in each of its orthogonal
	 * dimensions; default is pbc = 0 (off)
	 */
	int pbc;

	/**
	 * Following variables define a rectangular
	 * set of boundaries in order to make interpolant
	 * periodic over the given rectangular domain;
	 * "rpb" = "rectangular periodic boundaries"
	 */
	double *rpb;
} punity_t;

typedef struct
{
	/**
	 * The partition of unity to use to generate the
	 * basis functions
	 */
	punity_t pu;

	/**
	 * Local-to-global map of particle basis function
	 * to a global function index; zero and positive
	 * entries in each ltg[i] refering to monomial terms
	 * whereas negative indexes refer to input bases
	 * stored in bfuncs
	 */
	int **lbase;

	/**
	 * The number of basis functions to use for each
	 * particular particle
	 */
	int *nlbase;

	/**
	 * Number of external basis functions loaded at init;
	 * this is the size of extf
	 */
	int nextf;

	/**
	 * A set of function pointers given as input for
	 * potential use in each particle; the first entry
	 * in extf is referenced as -1 in ltg; if, say
	 * ltg[2][0] = -1, then this indicates that the first
	 * basis functions to use is the one stored in extf[0];
	 * if ltg[2][1] = -3, then the second basis function
	 * on particle 2 is extf[2]
	 */
	double (**extf)(int*,double*);
} pubasis_t;

double cubic_window( int dim_in, double a_in, double *x_in )
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

/**
 * The analytical or symbolic calculation of derivatives in generic function
 * spaces would greatly benefti from algorithmic differentiation. One example
 * is symbolic python, the packages `sympy`.
 *
 * @param dim_in The dimension of the cubic window
 * @param drv_in Order of derivative in the radius
 * @param a_in Scaling factor
 * @param x_in The position at which to evaluate the derivative
 * @return The scalar value of the derivative
 */
double cubic_window_deriv( int dim_in, int drv_in, double a_in, double *x_in )
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

double quartic_window( int dim_in, double a_in, double *x_in )
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

double quartic_window_deriv( int dim_in, int drv_in, double a_in, double *x_in )
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
double punity_window_evaluate( punity_t *obj_in, int idx_in, double *x_in )
{
	int i;
	double x[obj_in->dim];

	for(i=0;i<obj_in->dim;i++)
		x[i] = x_in[i] - obj_in->pts[idx_in*obj_in->dim+i];
	return obj_in->wfs( obj_in->dim, obj_in->dlt[idx_in], x );
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
double punity_window_evaluate_delta( punity_t *obj_in, int idx_in, double *x_in, int exp_in )
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
int punity_init( punity_t *obj_in, int dim_in, int npts_in, double *pts_in, double *dlt_in, double (*wfs_in)(int,double,double*), double (*wfsd_in)(int,int,double,double*) )
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

    /* Set pbc options to off by default and put NULL in rpb */
    obj_in->pbc = 0;
    obj_in->rpb = NULL;

#ifdef PUNITY_USE_KDTREES
	/* Use kdtrees for storage of points for fast access if faster this way */
	obj_in->kdt = (kdnode_t*) malloc( sizeof(kdnode_t) );
	kdtree_build_average( obj_in->kdt, obj_in->dim, obj_in->npts, obj_in->pts, 0, 0 );
#endif

	return 0;
}

/**
 * Turn on periodic boundaries and read the bounds
 * as input arguments given in rpb_in
 */
int punity_use_pbc( punity_t *obj_in, double *rpb_in )
{
    int i;

    /* Set pbc to 1 (on) */
    obj_in->pbc = 1;

    /* Make sure obj_in->rpb not already allocated */
    if( obj_in->rpb != NULL )
        return -1; /* Something funny is going on */
    else
        obj_in->rpb = (double*) malloc( obj_in->dim * sizeof(double) );

    /* Set the actual boundaries for ea. dimension */
    for(i=0;i<obj_in->dim;i++)
        obj_in->rpb[i] = rpb_in[i];

    return 0;
}

/**
 * Clean up memory allocated for punity_t
 * @param obj_in PU object to clean up
 */
int punity_free( punity_t *obj_in )
{
    if( obj_in->pts != NULL )
	    free( obj_in->pts );
    if( obj_in->dlt != NULL )
	    free( obj_in->dlt );
    if( obj_in->dlt != NULL )
	    free( obj_in->bdry );
    if( obj_in->pbc != 0 && obj_in->rpb != NULL )
    {
        free( obj_in->rpb );
        obj_in->pbc = 0;
    }
    return 0;
}

/**
 * Evaluate the particle function idx_in at point x_in.
 * @param obj_in Partition of unity object
 * @param idx_in Index of the function to evaluate
 * @param x_in Position at which to evaluate the function
 * @return Value of function idx_in at point x_in
 */
double punity_evaluate( punity_t *obj_in, int idx_in, double *x_in )
{
	int i,j;
	double res,sum,vec[obj_in->dim];

	/* If an index is given outside the range of number of points then return zero */
	if( idx_in > obj_in->npts - 1 )
		return 0.0;

	/* Don't bother with the denominator if the numerator is zero */
    if( obj_in->pbc == 0 )  /* Case 1: No periodic boundaries   */
    {
        sum = 0.0;
        for(i=0;i<obj_in->dim;i++)
            sum += pow( obj_in->pts[idx_in*obj_in->dim+i] - x_in[i], 2.0 );
    }
    else                    /* Case 2: Periodic boundaries      */
    {
        assert( obj_in->dim == 3 );
        sum = pbc_dist_real3( x_in, &(obj_in->pts[idx_in*obj_in->dim]), obj_in->rpb );
        sum = sum * sum; /* FIXME: This is NOT ideal; calculate sqr version of dist */
    }
    if( sum > obj_in->dlt[idx_in] * obj_in->dlt[idx_in] )
        return 0.0;

	/* NOTE: Change this to punity_neighbors_fast() */
    if( obj_in->pbc == 0 )  /* Case 1: No periodic boundaries   */
    {
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
    }
    else
    {
        /**
         * Prototypes from pbc.c used here:
         * void pbc_vec_real3( t_real *x_in, t_real *y_in, t_real *bx_in, t_real *z_out )
         * t_real pbc_dist_real3( t_real *x_in, t_real *y_in, t_real *bx_in );
         */
        assert( obj_in->dim == 3 );
        sum = 0.0;
        for(i=0;i<obj_in->npts;i++)
        {
            res = pbc_dist_real3( x_in, &(obj_in->pts[i*obj_in->dim]), obj_in->rpb );
            res = res * res; /* FIXME: As above fix this to take no sqrt() at all!!! */
            pbc_vec_real3( x_in, &(obj_in->pts[i*obj_in->dim]), obj_in->rpb, vec );
            if( res < obj_in->dlt[i] * obj_in->dlt[i] )
                sum += obj_in->wfs( obj_in->dim, obj_in->dlt[i], vec );
        }
        pbc_vec_real3( x_in, &(obj_in->pts[idx_in*obj_in->dim]), obj_in->rpb, vec );
    }
	return obj_in->wfs( obj_in->dim, obj_in->dlt[idx_in], vec ) / sum;
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
double punity_evaluate_delta( punity_t *obj_in, int idx_in, double *x_in, int exp_in )
{
	int i,j;
	double res,sum,vec[obj_in->dim];

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
			return obj_in->wfs( obj_in->dim, obj_in->dlt[idx_in], vec ) / sum;
		else
			return obj_in->wfs( obj_in->dim, obj_in->dlt[idx_in], vec ) / pow( res / obj_in->dlt[idx_in], (double) exp_in ) / sum;
	}
}

/**
 * Evaluate the derivative of r = | xi - xj | w.r.t. drv_in; FIXME: This function
 * returns values which become singular as r -> 0; figure out a way to get around this
 * @param dim_in Dimension of the vector space
 * @param drv_in Derivative integer vector
 * @param x_in Point at which to evaluate the derivative
 */
double radial_deriv_evaluate( int dim_in, int *drv_in, double *x_in )
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
double rational_radial_deriv_evaluate( int dim_in, int *drv_in, double *x_in, double a_in, int exp_in )
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
double punity_comp_deriv_evaluate( punity_t *obj_in, int idx_in, int *drv_in, int q_in, int *nb_in, double *x_in )
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
double punity_window_deriv_evaluate( punity_t *obj_in, int idx_in, int *drv_in, double *x_in )
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
double punity_window_deriv_evaluate_delta( punity_t *obj_in, int idx_in, int *drv_in, double *x_in, int exp_in )
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

int punity_neighbors( punity_t *obj_in, double *x_in, int **nb_out )
{
	int i,j,q,na,*nb;
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
                                nb = (int*) realloc( nb, na * sizeof(int) );
                        }
                        nb[q++] = i;
                }
        }
	*nb_out = nb;

	return q;
}

/**
 * Find all particles, i, which are withing obj_in->dlt[i] units from x_in
 * @param obj_in PU object
 * @param x_in Point to query
 * @param nb_out Output neighbors; allocated on-the-fly; you need to clean this up
 */
int punity_neighbors_fast( punity_t *obj_in, double *x_in, int **nb_out )
{
	int i,j,n,a,nn,aa,*id;
	double sum;
	const double rsq = obj_in->rmax * obj_in->rmax;

	id = NULL, n = 0, a = 0;
	kdtree_range_query( obj_in->kdt, x_in, obj_in->rmax, &id, &a, &n );

	*nb_out = NULL, nn = 0, aa = 0;
	for(i=0;i<n;i++)
	{
		sum = 0.0;
		for(j=0;j<obj_in->dim;j++)
			sum += pow( x_in[j] - obj_in->pts[id[i]*obj_in->dim+j], 2.0 );
		if( sum < obj_in->dlt[id[i]] * obj_in->dlt[id[i]] )
		{
			if( nn + 1 > aa )
			{
				(*nb_out) = (int*) realloc( (*nb_out), ( aa + PUNITY_NEIGHBOR_INC ) * sizeof(int) );
				aa = aa + PUNITY_NEIGHBOR_INC;
			}
			(*nb_out)[nn++] = id[i];
		}
	}
	return nn;
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
double punity_term_deriv_evaluate( punity_t *obj_in, int idx_in, int *drv_in, double *x_in, int *nb_in, int q_in )
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

	/* If n is zero then nothing to do */
	if( n <= 0 )
		return punity_evaluate( obj_in, idx_in, x_in );

	/* Otherwise allocate stuff */
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
double punity_term_deriv_evaluate_delta( punity_t *obj_in, int idx_in, int *drv_in, double *x_in, int exp_in, int *nb_in, int q_in )
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

	/* Check to see if n = 0 */
	if( n <= 0 )
		return punity_evaluate_delta( obj_in, idx_in, x_in, exp_in );

	/* Allocate stuff if n > 0 */
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
double punity_term_first_deriv_evaluate( punity_t *obj_in, int idx_in, int drv_in, double *x_in )
{
	
}

/**
 * Evaluate particle-localized polynomial term
 */
double punity_eval_local_poly( punity_t *obj_in, int idx_in, int *pidx_in, double *x_in )
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
double punity_eval_local_poly_delta( punity_t *obj_in, int idx_in, int *pidx_in, double *x_in, int exp_in )
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
double punity_term_deriv_local_poly( punity_t * obj_in, int idx_in, int *pidx_in, int *drv_in, double *x_in, int *nb_in, int q_in )
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
	for(i=0;i<=n;i++)
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
						prd *= pow( x_in[m] - obj_in->pts[idx_in*obj_in->dim+m], (double) p ) / pow( obj_in->dlt[idx_in], (double) pidx_in[m] );
						for(p=0;p<qc[m];p++)
							prd *= (double) ( pidx_in[m] - p );
					}
				}
			}
			if( n - i == 0 )
				prd *= punity_evaluate( obj_in, idx_in, x_in );
			else /* Evaluate the derivative of the partition of unity function */
				prd *= punity_term_deriv_evaluate( obj_in, idx_in, qd, x_in, nb_in, q_in );

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
double punity_term_deriv_local_poly_delta( punity_t * obj_in, int idx_in, int *pidx_in, int *drv_in, double *x_in, int exp_in, int *nb_in, int q_in )
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
	for(i=0;i<=n;i++)
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
						prd *= pow( x_in[m] - obj_in->pts[idx_in*obj_in->dim+m], (double) p ) / pow( obj_in->dlt[idx_in], (double) pidx_in[m] );
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
 * A version of punity_window_deriv_evaluate allowing arbitrary specification
 * of the radial derivative function
 */
double radial_window_deriv_evaluate( int dim_in, double (*wfsd_in)(int,int,double,double*), double dlt_in, int *drv_in, double *x_in )
{
	int b,c,i,j,p,n,*v,*w,*s,*m,*t;
        double d,prd,sum = 0.0;

	/* Set up the vector to partition */
        n = 0;
        for(i=0;i<dim_in;i++)
                n += drv_in[i];
        v = (int*) malloc( n * sizeof(int) );
        w = (int*) malloc( n * sizeof(int) );
        s = (int*) malloc( n * sizeof(int) );
        m = (int*) malloc( n * sizeof(int) );
	t = (int*) malloc( dim_in * sizeof(int) );
        for(i=0,p=0;i<dim_in;i++)
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
		return wfsd_in( dim_in, 0, dlt_in, x_in );
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
			for(j=0;j<dim_in;j++)
				t[j] = 0;
			for(j=0;j<c;j++)
				t[w[j]]++; /* Everytime an index appears in w, increment its component once */

                        /* Do the derivative */
			d = radial_deriv_evaluate( dim_in, t, x_in ); /* No dilation factor here; all contained in wfsd */

                        /* Multiply the derivative in d for this block into the total product */
                        prd *= d; /* Variable d is derivative of u consistent with this block */
                }

		/* Calculate the b order derivative of the radial polynomial */
		prd *= wfsd_in( dim_in, b, dlt_in, x_in );

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
 * Generalized version of punity_window_deriv_evaluate_delta but with given
 * function pointer
 */
double radial_window_deriv_evaluate_delta( int dim_in, double (*wfsd_in)(int,int,double,double*), double dlt_in, int *drv_in, double *x_in, int exp_in )
{
	int i,j,k,m,n,p,q,r,*qc,*qd,*v,*cmb,*dmb;
	double d,prd,sum = 0.0;

	/* Do some setup */
	n = 0;
	for(i=0;i<dim_in;i++)
		n += drv_in[i];
	v = (int*) malloc( n * sizeof(int) );
	qc = (int*) malloc( dim_in * sizeof(int) );
	qd = (int*) malloc( dim_in * sizeof(int) ); /* Complement of q w.r.t. the current cmb state */
	cmb = (int*) malloc( n * sizeof(int) );
	dmb = (int*) malloc( n * sizeof(int) );

	/* Calculate the derivative vector; a sequence of the numbers 0,...,dim-1 */
	for(i=0,p=0;i<dim_in;i++)
		for(j=0;j<drv_in[i];j++)
			v[p++] = i;

	/* Calculate the radial vector */
	d = 0.0;
	for(i=0;i<dim_in;i++)
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
		return wfsd_in( dim_in, 0, dlt_in, x_in ) / pow( d / dlt_in, (double) exp_in );
	}

	/* Iterate through all combinations sizes and for each size, iterate through all combinations */
	for(i=0;i<=n;i++)
	{
		combinadic_init( n, i, cmb );
		k = binomial( n, i );
		for(j=0;j<k;j++)
		{
			/* Form the derivative vector for this case */
			for(m=0;m<dim_in;m++)
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
			for(m=0;m<dim_in;m++)
				qd[m] = 0;
			for(m=0;m<n-i;m++)
				qd[v[dmb[m]]]++;

			/* Apply the derivative to the window function and evaluate it */
			prd = 1.0;
			if( i == 0 )
				prd *= pow( d, (double) ( -exp_in ) );
			else
				prd *= rational_radial_deriv_evaluate( dim_in, qc, x_in, dlt_in, exp_in );
			if( n - i == 0 )
				prd *= wfsd_in( dim_in, 0, dlt_in, x_in );
			else
				prd *= radial_window_deriv_evaluate( dim_in, wfsd_in, dlt_in, qd, x_in );

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
 * Calculate arbitrary derivatives of localized radial function phi(x) * wfs(x);
 * radial because radial functions are simple, so can calculate any
 * derivative
 */
double punity_term_deriv_local_radial( punity_t * obj_in, int idx_in, double (*wfsd_in)(int,int,double,double*), int *drv_in, double *x_in, int *nb_in, int q_in )
{
	int i,j,k,m,n,p,q,r,*qc,*qd,*v,*cmb,*dmb;
	double d,prd,*y,sum = 0.0;

	/* Do some setup */
	n = 0;
	for(i=0;i<obj_in->dim;i++)
		n += drv_in[i];

	/* Allocate stuff */
	v = (int*) malloc( n * sizeof(int) );
	qc = (int*) malloc( obj_in->dim * sizeof(int) );
	qd = (int*) malloc( obj_in->dim * sizeof(int) ); /* Complement of q w.r.t. the current cmb state */
	cmb = (int*) malloc( n * sizeof(int) );
	dmb = (int*) malloc( n * sizeof(int) );
	y = (double*) malloc( obj_in->dim * sizeof(double) );

	/* Form y = x - x_i */
	for(i=0;i<obj_in->dim;i++)
		y[i] = x_in[i] - obj_in->pts[idx_in*obj_in->dim+i]; /* Scaling via obj_in->dlt[idx_in] is done via window function */

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
		d = wfsd_in( obj_in->dim, 0, obj_in->dlt[idx_in], y ) * punity_evaluate( obj_in, idx_in, x_in );
		free( y );
		return d;
	}

	/* Iterate through all combinations sizes and for each size, iterate through all combinations */
	for(i=0;i<=n;i++)
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
			prd *= radial_window_deriv_evaluate( obj_in->dim, wfsd_in, obj_in->dlt[idx_in], qc, y );
			prd *= punity_term_deriv_evaluate( obj_in, idx_in, qd, x_in, nb_in, q_in );

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
	free( y );

	return sum;
}

/**
 * Same as punity_term_deriv_local_radial() except use singular particles
 * where applicable
 */
double punity_term_deriv_local_radial_delta( punity_t * obj_in, int idx_in, double (*rfd_in)(int,double), int *drv_in, double *x_in, int exp_in, int *nb_in, int q_in )
{
	
}

/**
 * Partition of unity basis initializer
 */
int pubasis_init( pubasis_t *obj_in, int dim_in, int npts_in, double *pts_in, double *dlt_in, int *nlbase_in, int **lbase_in, double (*wfs_in)(int,double,double*), double (*wfsd_in)(int,int,double,double*), int nextf_in, double (**extf)(int,double,double*) )
{
	int i,j;

	/* Allocate the list for storage of local basis set indexes */
	obj_in->nlbase = (int*) malloc( npts_in * sizeof(int) );
	obj_in->lbase = (int**) malloc( npts_in * sizeof(int*) );
	for(i=0;i<npts_in;i++)
		obj_in->nlbase[i] = nlbase_in[i];
	for(i=0;i<npts_in;i++)
	{
		obj_in->lbase[i] = (int*) malloc( nlbase_in[i] * sizeof(int) );
		for(j=0;j<nlbase_in[i];j++)
			obj_in->lbase[i][j] = lbase_in[i][j];
	}
	
	return punity_init( &(obj_in->pu), dim_in, npts_in, pts_in, dlt_in, wfs_in, wfsd_in );
}

/**
 * Partition of unity basis cleanup
 */
int pubasis_free( pubasis_t *obj_in )
{
	
}

// vim: ts=4:sts=4:sw=4:et:sta
