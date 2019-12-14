#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "math_config.h"

/**
 * This file defines a data structure and supporting
 * functions which are used to generate and evaluate
 * a Lagrange basis on the n-simplex for any desired
 * dimension and polynomial space (order).
 */

typedef struct
{
    long dim; ///< Dimensionality of points
    long deg; ///< Degree of polynomial space
    long pdim; ///< Dimensionality of polynomial space
    double *nodes; ///< Nodes used to define lagrange basis
    double *coeffs; ///< Coefficients (pdim of them) defining each of the pdim basis functions
    double *vdm; ///< The Vandermonde matrix of the function/point system
    double *vdmlu; ///< The LU decomposition of the Vandermonde matrix in \c vdm
} lagrange_t;

/**
 * This function initializes all internal data
 * to the right data consistent with the dimension
 * and degree input.
 * \param lag_in The object
 * \param dim_in The dimensionality
 * \param deg_in The polynomial degree
 * \param nodes_in The nodes defining the Vandermonde matrix
 * \brief Initializer for lagrange_t
 */
int lagrange_init( lagrange_t *lag_in, long dim_in, long deg_in, double *nodes_in )
{
	/* Set parameters internally */
	lag_in->dim = dim_in;
	lag_in->deg = deg_in;
	lag_in->pdim = binomial_long( dim_in + deg_in, dim_in ); ///< Number of terms less than or equal to degree \c deg

	/* Allocate space for pdim vertices */
	lag_in->nodes = malloc( dim_in * lag_in->pdim * sizeof(double) );
	if( lag_in->nodes == NULL )
		return -1;

	/* Allocate space for coefficients defining the Lagrange shape functions */
	lag_in->coeffs = malloc( lag_in->pdim * lag_in->pdim * sizeof(double) );
	if( lag_in->coeffs == NULL )
		return -1;

	/* Allocate space for internal Vandermonde matrix */
	lag_in->vdm = malloc( lag_in->pdim * lag_in->pdim * sizeof(double) );
	if( lag_in->vdm == NULL )
		return -1;

	/* Allocate space for LU decomposition of vdm */
	lag_in->vdmlu = malloc( lag_in->pdim * lag_in->pdim * sizeof(double) );
	if( lag_in->vdmlu == NULL )
		return -1;

	/* Copy nodes if nodes_in not null */
	if( nodes_in != NULL )
		bcopy( nodes_in, lag_in->nodes, dim_in * lag_in->pdim * sizeof(double) );

	/* If no other problems, things are fine */
	return 0;
}

/**
 * This function deallocates all relevant data structures
 * \param lag_in Data structure
 * \brief Destructor for lagrange_t
 */
int lagrange_free( lagrange_t *lag_in )
{
	if( lag_in->nodes != NULL )
		free( lag_in->nodes );
	if( lag_in->coeffs != NULL )
		free( lag_in->coeffs );
	if( lag_in->vdm != NULL )
		free( lag_in->vdm );
	if( lag_in->vdmlu != NULL )
		free( lag_in->vdmlu );

	return 0;
}

/**
 * This function asserts that 
 */
int lagrange_set_nodes( lagrange_t *lag_in, double *nodes_in )
{
	if( nodes_in == NULL || lag_in->nodes == NULL  )
		return -1;
	bcopy( nodes_in, lag_in->nodes, lag_in->pdim );

	return 0;
}

/**
 * This function automatically loads the standard vertices
 * associated with the standard \c dim -simplex .
 */
int lagrange_set_nodes_standard( lagrange_t *lag_in )
{
	long i,j,k,m,index[lag_in->dim+lag_in->deg],exp[lag_in->dim];

	if( lag_in->nodes == NULL )
		return -1;

	for(m=0,k=0;m<=lag_in->deg;m++)
	{
		rcombinadic_init_long( m, lag_in->dim, index );
		for(i=0;i<binomial_long(lag_in->dim+m-1,lag_in->dim-1);i++)
		{
			rcombinadic_occupancy_long( m, lag_in->dim, index, exp );
			for(j=0;j<lag_in->dim;j++)
				lag_in->nodes[k*lag_in->dim+j] = (double) exp[j] / (double) lag_in->deg;
			rcombinadic_next_long( m, lag_in->dim, index );
			++k;
		}
	}

	return 0;
}

/**
 * Generates the coefficients defining the Lagrange basis
 * on the standard element. Calculates the Vandermonde matrix
 * and LU decomposes it, placing the result in \c vdmlu , and
 * fills in \c coeffs , each row of which contains the sequence
 * of \c pdim coefficients defining the shape functions.
 */
int lagrange_generate( lagrange_t *lag_in )
{
	/* Local variables */
	long i,k,m,n,p,index[lag_in->dim+lag_in->deg],exp[lag_in->dim];
	double x;

	/* FORTRAN variables for linear algebra */
	int q = (int) lag_in->pdim;
	int ipiv[lag_in->pdim];
	int info;

	/* Dealing with possible segmentation faults */
	if( lag_in->nodes == NULL || lag_in->coeffs == NULL || lag_in->vdm == NULL || lag_in->vdmlu == NULL )
		return -1;

	/* For each monomial order from 0 to lag_in->deg */
	for(m=0,k=0;m<=lag_in->deg;m++)
	{
		rcombinadic_init_long( m, lag_in->dim, index );
		for(i=0;i<binomial_long(lag_in->dim+m-1,lag_in->dim-1);i++)
		{
			rcombinadic_occupancy_long( m, lag_in->dim, index, exp );
			for(n=0;n<lag_in->pdim;n++)
			{
				x = 1.0;
				for(p=0;p<lag_in->dim;p++)
					x *= pow( lag_in->nodes[n*lag_in->dim+p], (double) exp[p] );
				lag_in->vdm[k*lag_in->pdim+n] = x;
			}
			rcombinadic_next_long( m, lag_in->dim, index );
			++k;
		}
	}

	/* Make the unitary matrix in coeffs */
	for(i=0;i<lag_in->pdim;i++)
		for(k=0;k<lag_in->pdim;k++)
			if( i == k )
				lag_in->coeffs[i*lag_in->pdim+k] = 1.0;
			else
				lag_in->coeffs[i*lag_in->pdim+k] = 0.0;

	/* Make a copy of vdm */
	bcopy( lag_in->vdm, lag_in->vdmlu, lag_in->pdim * lag_in->pdim * sizeof(double) );

	/* Do the solve */
	dgesv_( &q, &q, lag_in->vdmlu, &q, ipiv, lag_in->coeffs, &q, &info );

	/* Return value returned by dgesv_ */
	return info;
}

/**
 * Evaluate a Lagrange polynomial at a given point.
 * \param lag_in The Lagrange interpolant basis
 * \param index_in The index of the basis polynomial to evaluate
 * \param point_in The point at which to evaluate function \c index_in
 */
double lagrange_evaluate( lagrange_t *lag_in, long index_in, double *point_in )
{
	long i,j,k,m,index[lag_in->dim+lag_in->deg],exp[lag_in->dim];
	double x,sum=0.0;

	for(m=0,k=0;m<=lag_in->deg;m++)
	{
		rcombinadic_init_long( m, lag_in->dim, index );
		for(i=0;i<binomial_long(m+lag_in->dim-1,lag_in->dim-1);i++)
		{
			rcombinadic_occupancy_long( m, lag_in->dim, index, exp );
			x = 1.0;
			for(j=0;j<lag_in->dim;j++)
				x *= pow( point_in[j], (double) exp[j] );
			sum += lag_in->coeffs[index_in*lag_in->pdim+k] * x;
			rcombinadic_next_long( m, lag_in->dim, index );
			++k;
		}
	}

	return sum;
}

/**
 * Evaluate an interpolant
 */
double lagrange_interpolate( lagrange_t *lag_in, double *coeffs_in, double *point_in )
{
	
}

