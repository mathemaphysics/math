#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "math_config.h"

typedef struct
{
	long dim; ///< Dimensionality of points
	long deg; ///< Degree of polynomial to use
	long pdim; ///< Dimensionality of polynomial space
	double *nodes; ///< Nodes to use for integral evaluation
	double *coeffs; ///< Coefficients defining the quadrature rule
	double *vdm; ///< Stores Vandermonde matrix
	double *vdmlu; ///< Stores the LU decomposition of the Vandermonde matrix
} quadrature_t;

/**
 * I love this function. This function does dimension-generalized
 * (i.e. for any dimension) integration of arbitrary polynomial terms.
 * The dim_in term is the dimension, and the exp_in variable is a pointer
 * to dim_in long integers defining the exponents for each independent
 * variable.
 * \param dim_in The number of variables in the monomial term
 * \param exp_in The exponent vector defining the monomial term to integrate over the simplex
 * \brief Calculates the integral of \c x1^a1 \c + \c ... \c + \c xn^an for \c dim=n over the standard \c n -simplex
 */
double quadrature_term_integral( long dim_in, long *exp_in )
{
	long i,j;
	long num = 1;
	long dnm = 0;

	for(i=0;i<dim_in;i++)
	{
		num *= factorial( exp_in[i] );
		dnm += exp_in[i];
	}

	dnm += dim_in;
	return (double) num / (double) factorial( dnm );
}

/**
 * This function initializes the quadrature_t data structure,
 * which generates a quadrature rule over the n-simplex for
 * numerical integration
 * \param quad_in Pointer to the object
 * \param dim_in Dimensionality of the points
 * \param deg_in Degree of polynomial
 * \param nodes_in The nodes to use as quadrature
 * \brief Initializer for quadrature_t
 */
int quadrature_init( quadrature_t *quad_in, long dim_in, long deg_in, double *nodes_in )
{
	/* Set basics */
	quad_in->dim = dim_in;
	quad_in->deg = deg_in;
	quad_in->pdim = binomial_long( dim_in + deg_in, dim_in );

	/* Allocate space for node storage */
	quad_in->nodes = (double*) malloc( dim_in * quad_in->pdim * sizeof(double) );
	if( quad_in->nodes == NULL )
		return -1;

	/* Allocate space for coefficients */
	quad_in->coeffs = (double*) malloc( quad_in->pdim * sizeof(double) );
	if( quad_in->coeffs == NULL )
		return -1;

	/* Allocate space for Vandermonde matrix */
	quad_in->vdm = (double*) malloc( quad_in->pdim * quad_in->pdim * sizeof(double) );
	if( quad_in->vdm == NULL )
		return -1;

	/* Allocate space for LU decomposition of vdm */
	quad_in->vdmlu = (double*) malloc( quad_in->pdim * quad_in->pdim * sizeof(double) );
	if( quad_in->vdmlu == NULL )
		return -1;

	/* If possible, copy nodes from input */
	if( nodes_in != NULL )
		bcopy( nodes_in, quad_in->nodes, dim_in * quad_in->pdim * sizeof(double) );

	/* Otherwise nothing went wrong */
	return 0;
}

/**
 * This function frees all data internal to the quadrature_t
 * data structure
 * \param quad_in Pointer to the object
 * \brief Destructor for quadrature_t
 */
int quadrature_free( quadrature_t *quad_in )
{
	if( quad_in->nodes != NULL )
		free( quad_in->nodes );
	if( quad_in->coeffs != NULL )
		free( quad_in->coeffs );
	if( quad_in->vdm != NULL )
		free( quad_in->vdm );
	if( quad_in->vdmlu != NULL )
		free( quad_in->vdmlu );

	return 0;
}

int quadrature_set_nodes( quadrature_t *quad_in, double *nodes_in )
{
	if( nodes_in == NULL || quad_in->nodes == NULL )
		return -1;
	bcopy( nodes_in, quad_in->nodes, quad_in->pdim );

	return 0;
}

int quadrature_set_nodes_standard( quadrature_t *quad_in )
{
	long i,j,k,m,index[quad_in->dim+quad_in->deg],exp[quad_in->dim];

        if( quad_in->nodes == NULL )
                return -1;

        for(m=0,k=0;m<=quad_in->deg;m++)
        {
                rcombinadic_init_long( m, quad_in->dim, index );
                for(i=0;i<binomial_long(quad_in->dim+m-1,quad_in->dim-1);i++)
                {
                        rcombinadic_occupancy_long( m, quad_in->dim, index, exp );
                        for(j=0;j<quad_in->dim;j++)
                                quad_in->nodes[k*quad_in->dim+j] = (double) exp[j] / (double) quad_in->deg;
                        rcombinadic_next_long( m, quad_in->dim, index );
                        ++k;
                }
        }

        return 0;
}

/**
 * Sets up all the data inside by calculating the Vandermonde
 * matrix (transposed relative to the one in lagrange_t) and
 * inverting it to generate the set of coefficients for the domain.
 */
int quadrature_generate( quadrature_t *quad_in )
{
	/* Local variables */
        long i,k,m,n,p,index[quad_in->dim+quad_in->deg],exp[quad_in->dim];
        double x;

        /* FORTRAN variables for linear algebra */
        int q = (int) quad_in->pdim;
	int one = 1;
        int ipiv[quad_in->pdim];
        int info;

        /* Dealing with possible segmentation faults */
        if( quad_in->nodes == NULL || quad_in->coeffs == NULL || quad_in->vdm == NULL || quad_in->vdmlu == NULL )
                return -1;

        /* For each monomial order from 0 to lag_in->deg; this matrix is the transpose of lag->vdm */
        for(m=0,k=0;m<=quad_in->deg;m++)
        {
                rcombinadic_init_long( m, quad_in->dim, index );
                for(i=0;i<binomial_long(quad_in->dim+m-1,quad_in->dim-1);i++)
                {
                        rcombinadic_occupancy_long( m, quad_in->dim, index, exp );
                        for(n=0;n<quad_in->pdim;n++)
                        {
                                x = 1.0;
                                for(p=0;p<quad_in->dim;p++)
                                        x *= pow( quad_in->nodes[n*quad_in->dim+p], (double) exp[p] );
                                quad_in->vdm[n*quad_in->pdim+k] = x;
                        }

			/* Put the integral of the current term in y */
			quad_in->coeffs[k] = quadrature_term_integral( quad_in->dim, exp );

			/* Next exponent of order m */
                        rcombinadic_next_long( m, quad_in->dim, index );
                        ++k;
                }
        }

	/* Now invert the matrix in vdm */
	bcopy( quad_in->vdm, quad_in->vdmlu, quad_in->pdim * quad_in->pdim * sizeof(double) );

	/* Actually solve */
	dgesv_( &q, &one, quad_in->vdmlu, &q, ipiv, quad_in->coeffs, &q, &info );

	/* Return value from dgesv_ */
	return info;
}

/**
 * This function takes a pointer to a function whose integral over
 * the domain is desired and calculates the integral using the internal
 * quadrature rule.
 */
double quadrature_integrate( quadrature_t *quad_in, double (*func_in)(double*,void*) )
{

}

