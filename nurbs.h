#ifndef NURBS_H
#define NURBS_H

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

#endif
