#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include <math.h>
#include "config.h"

/**
 * This is an implementation of NURBS which which have
 * \c nparam input parameters which are each mutually
 * orthogonal, i.e. the case in which \c dim is greater
 * than 1 is achieved by tensor product NURBS. There is
 * another version of this code which makes use of
 * simplex-based interpolation in \c d dimensions, but
 * which suffers from the limitation that on the boundaries
 * of simplices all control points have the surface passing
 * exactly through them.
 */
typedef struct
{
    /**
     * The spatial dimension in which points are taken.
     */
    long dim;

    /**
     * This is the number of parameters defining
     * the objects as a line/surface/volume/etc.,
     */
    long nparam;

    /**
     * The degree of interpolation used to define
     * the spline curves in each of the \c nparam
     * orthogonal parameters. For example, the
     * the degree of interpolation in parameter \c i
     * is \c deg[i] where \c i ranges from 0 to
     * \c nparam - 1.
     */
    long *deg;

    /**
     * The start index of knots for parameter \c i .
     */
    long *iknots;

    /**
     * These are the actual knot locations for each
     * dimension in the parameter space.
     */
    double *knots;

    /**
     * The weights defining the interpolant. There
     * are \c n1 * \c ... \c ni , which are the
     * numbers of knots in each parameter, and \c i
     * is \c nparam . This array makes use of
     * multi-indexing in a one-dimensional array
     * for generalization of the number of parameters.
     */
    double *weights;

    /**
     * These are the control points. These points lie
     * in dimension \c dim . There are \c nparam of
     * indexes into \c cpoints . Space allocated to
     * \c cpoints must be \c n1 * ... * \c ni , where
     * \c ni is the number of knots in parameter \c i .
     */
    double *cpoints;
} nurbs_t;

/**
 * The input \c knots_in vector must contain the right number
 * of vertices for degree \c deg_in interpolation; need (at least) n + 1
 * vertices in the knot vector to specify such a basis function.
 * This is a basis outline of de Boor's algorithm for recursive
 * evaluation of B-splines.
 */
double spline( long deg_in, double u_in, double *knots_in )
{
    double f,g;

    if( deg_in > 0 )
    {
        if( knots_in[deg_in] - knots_in[0] > 0.0 )
            f = ( u_in - knots_in[0] ) / ( knots_in[deg_in] - knots_in[0] );
        else
            f = 0.0;
        if( knots_in[deg_in+1] - knots_in[1] > 0.0 )
            g = ( knots_in[deg_in+1] - u_in ) / ( knots_in[deg_in+1] - knots_in[1] );
        else
            g = 0.0;
        return f * spline( deg_in - 1, u_in, knots_in ) + g * spline( deg_in - 1, u_in, knots_in + 1 );
    }
    else
    {
        if( knots_in[0] <= u_in && u_in < knots_in[1] )
            return 1.0;
        else
            return 0.0;
    }
}

/**
 * This is the basic function which makes use of \c spline to
 * recursively compute values of non-uniform rational splines
 * given the particular knot index and the total list of knots
 * and weights.
 */
double nurbs( long deg_in, double u_in, long kdx_in, long nknots_in, double *knots_in, double *weights_in )
{
    long i;
    double num,dnm;

    assert( kdx_in < nknots_in - deg_in - 1 );
    num = weights_in[kdx_in] * spline( deg_in, u_in, knots_in + kdx_in );
    for(i=0,dnm=0.0;i<nknots_in-deg_in-1;i++)
        dnm += weights_in[i] * spline( deg_in, u_in, knots_in + i );
    return num / dnm;
}

/**
 * This function initializes all internal data in \c nurbs_t .
 * \brief Constructor for \c nurbs_t
 *
 */
int nurbs_init( nurbs_t *nurbs_in, long dim_in, long nparam_in, long *deg_in, long *iknots_in, double *knots_in, double *cpoints_in, double *weights_in )
{
    long i,idx;
    nurbs_in->dim = dim_in;
    nurbs_in->nparam = nparam_in;

    /* Allocate space for the degree integers for each parameter and copy the values from the input if available */
    nurbs_in->deg = malloc( nparam_in * sizeof(long) );
    bcopy( deg_in, nurbs_in->deg, nparam_in * sizeof(long) );

    /* Allocate space for \c iknots and copy the data over if input pointer is not null */
    nurbs_in->iknots = malloc( (nparam_in+1) * sizeof(long) );
    bcopy( iknots_in, nurbs_in->iknots, (nparam_in+1) * sizeof(long) );

    /* Copy the actual knots into the object */
    nurbs_in->knots = malloc( iknots_in[nparam_in] * sizeof(double) );
    bcopy( knots_in, nurbs_in->knots, iknots_in[nparam_in] * sizeof(double) );

    /* Calculate the product of all parameter knots to allocate space for \c cpoints and weights */
    for(i=0,idx=1;i<nparam_in;i++)
        idx = idx * ( iknots_in[i+1] - iknots_in[i] - deg_in[i] - 1 ); /* There are only n - deg - 1 basis functions */
    nurbs_in->cpoints = malloc( idx * dim_in * sizeof(double) );
    bcopy( cpoints_in, nurbs_in->cpoints, idx * dim_in * sizeof(double) ); // Need space for a length \c dim vector for each shape function
    nurbs_in->weights = malloc( idx * sizeof(double) ); // Need space for a single integer for each shape function
    bcopy( weights_in, nurbs_in->weights, idx * sizeof(double) );

    return 0;
}

/**
 * Frees all internal allocated data in \c nurbs_t .
 * \brief Destructor for \c nurbs_t
 *
 */
int nurbs_free( nurbs_t *nurbs_in )
{
    free( nurbs_in->deg );
    free( nurbs_in->iknots );
    free( nurbs_in->knots );
    free( nurbs_in->cpoints );
    free( nurbs_in->weights );
}

/**
 * \brief Loads the points and weights defining the interpolation from a file or array
 */
int nurbs_load( nurbs_t *nurbs_in, char *data_in )
{

}

/**
 * \brief Evaluate the function within \c nurbs_t
 * This method calculates the value of the parametrized
 * geometric object at a given parameter value. Important here
 * is to remember that if there are \c N knots in a parametric
 * direction and the basis is of order \c K , then there are
 * exactly \c N-K-1 total functions. This is due to the fact
 * that functions are indexed by the point at which they
 * begin, rather than the position at which they are centered;
 * order \c K spline basis functions cover the span of \c K+1
 * knot distances, i.e. it contacts \c K+2 different consecutive
 * knot nodes.
 */
long arrayindex( long, long *, long * );
void nurbs_evaluate( nurbs_t *nurbs_in, double *x_in, double *x_out )
{
    long i,j,index[nurbs_in->nparam],size[nurbs_in->nparam];
    double val,sum;

    /* Calculate number of knots in each parameter dimension; put it in size */
    for(i=0;i<nurbs_in->nparam;i++)
        size[i] = nurbs_in->iknots[i+1] - nurbs_in->iknots[i] - nurbs_in->deg[i] - 1; /* There are nknots - deg functions for nknots knots */

    /* Initialize the multi-index */
    for(i=0;i<nurbs_in->nparam;i++)
        index[i] = 0;

    /* Function needed to generalize the parameter dimension so number of loops can vary */
    inline int next()
    {
        for(j=nurbs_in->nparam-1;j>=0;j--)
        {
            if( index[j] < size[j] )
            {
                ++index[j];
                return 0;
            }
            else /* the carry */
                index[j] = 0; /* if j = 0, sets index[0] = 0 but does not carry, terminates returning -1 indicating reset */
        }
        return -1; /* Should return -1 iff counter has been reset to 0,0,...,0 */
    }

    /* Initialize x_out to zero */
    for(i=0;i<nurbs_in->dim;i++)
        x_out[i] = 0.0;
    sum = 0.0;
    do
    {
        /* Evaluate function indexed by index[] */
        val = 1.0;
        for(i=0;i<nurbs_in->nparam;i++)
            val *= spline( nurbs_in->deg[i], x_in[i], nurbs_in->knots + nurbs_in->iknots[i] + index[i] );

        /* Calculate the index into both weights and cpoints */
        j = arrayindex( nurbs_in->nparam, size, index );

        /* Keep running sum for denominator part of the NURBS */
        sum += val * nurbs_in->weights[j];

        /* Add contribution of control point indexed by index; control point coresponding to current weight and knot */
        for(i=0;i<nurbs_in->dim;i++)
            x_out[i] = x_out[i] + val * nurbs_in->weights[j] * nurbs_in->cpoints[j*nurbs_in->dim+i];
    }
    while( next() == 0 );

    /* Divide by the denominator */
    for(i=0;i<nurbs_in->nparam;i++)
        x_out[i] /= sum;
}

// vim: ts=4:sts=4:sw=4:et:sta
