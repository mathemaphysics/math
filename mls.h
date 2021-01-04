#ifndef MLS_H
#define MLS_H

#include "quadrature.h"

typedef struct
{
        /**
         * The number of points used in the interpolation;
         * this is also the number of interpolating basis
         * functions used in the reconstruction
         */
        int np;

        /**
         * The dimension of the points given as input
         */
        int dim;

        /**
         * The degree of interpolating polynomial to seek as
         * the weighted reconstruction
         */
        int deg;

        /**
         * The list of points used in the interpolation
         */
        double *pts;

		/**
		 * The dilation parameter for each point in pts
		 */
		double *dlts;

		/**
         * The list of function values at their corresponding
         * point in pts
         */
        double *vals;

        /**
         * Function pointer to the weight function; this function
         * is shifted to properly center in about each interpolating
         * point stored in pts; the first argument is the dimension;
         * the second is the actual point, which must be centered on
         * the appropriate point by the user, at which to evaluate
         */
        double (*wfs)(int,double,double*);
} mls_t;

typedef struct
{
	/**
	 * The total number of points to use as particle centers
	 */
	int np;

	/**
	 * The dimension of the points
	 */
	int dim;

	/**
	 * The polynomial degree to reproduce
	 */
	int deg;

	/**
	 * The points themselves
	 */
	double *pts;

	/**
	 * List of indexes into pts which lie on the
	 * boundary
	 */
	int *bpts;

	/**
	 * The dilation factors which define the radii
	 * of the window functions
	 */
	double *dlts;

	/**
	 * The values to interpolate; only necessary
	 * when calling global evaluate functions wuch
	 * as rkp_evaluate
	 */
	double *vals;

	/**
	 * The quadrature weights to use for producing
	 * the particle functions
	 */
	double *gqw;

	/**
	 * A pointer to the window function to use
	 * for producing correction polynomials with
	 * the moving least squares idea
	 */
	double (*wfs)(int,double,double*);

	/**
	 * A pointer to a function which gives the
	 * gradient of the window function above; this
	 * most likely will not be needed
	 */
	void (*wfsg)(int,double,double*,double*);

	/**
	 * Radius of the unscaled window function
	 * given as input; this is used to decide
	 * when function values will be identically
	 * zero since they are outside of the window
	 */
	double wrad;

	/**
	 * A pdim by pdim matrix which acts as internal
	 * storage for the Gram moment matrix when generating
	 * the correction coefficients for the basis functions
	 * and wavelets
	 */
	double *grm;

	/**
	 * Space allocated for LU factorization of the moment
	 * matrix
	 */
	double *grmlu;

	/**
	 * Coefficients for correction functions at each point
	 * listed in pts; this is used for fast evaluation of
	 * values for point collocation PDE solvers
	 */
	double *coeffs;

	/**
	 * These are the wavelet coefficients stored for
	 * each point for fast evaluation of derivatives in
	 * point collocation schemes; there are np * pdim * (pdim-1)
	 * entries
	 */
	double *wcoeffs;

	/**
	 * Mode which currently has no use; will likely
	 * be dealt with later to reflect, for example,
	 * which scaling mode is being used: none, scaled,
	 * or const
	 */
	int mode;
} rkp_t;

int mls_init( mls_t *, int, int, int, double *, double *, double *, double (*)(int,double,double*) );

int rkp_init( rkp_t *, int, int, int, double *, double *, double *, double *, double (*)(int,double,double*), double );

int mls_free( mls_t * );

int rkp_free( rkp_t * );

int mls_basis_evaluate( mls_t *, double *, double * );

int mls_matrix_evaluate( mls_t *, double, double *, double * );

int rkp_matrix_evaluate( rkp_t *, double * );

int rkp_matrix_evaluate_scaled( rkp_t *, double * );

int rkp_matrix_evaluate_scaled_const( rkp_t *, double *, double );

int rkp_matrix_generate( rkp_t *, double *, double * );

int rkp_matrix_generate_scaled( rkp_t *, double *, double * );

int rkp_matrix_generate_scaled_const( rkp_t *, double *, double *, double );

int rkp_wavelet_generate_order( rkp_t *, int, double *, double * );

int rkp_wavelet_generate_order_scaled( rkp_t *, int, double *, double * );

int rkp_wavelet_generate_order_scaled_const( rkp_t *, int, double *, double, double * );

int rkp_wavelet_generate( rkp_t *, double *, double * );

int rkp_wavelet_generate_scaled( rkp_t *, double *, double * );

int rkp_wavelet_generate_scaled_const( rkp_t *, double *, double, double * );

int rkp_basis_generate( rkp_t * );

int rkp_basis_generate_scaled( rkp_t * );

int rkp_wavelet_basis_generate( rkp_t * );

int rkp_wavelet_basis_generate_scaled( rkp_t * );

int rkp_wavelet_basis_generate_scaled_const( rkp_t * );

double mls_term_evaluate( mls_t *, int, double * );

double rkp_term_evaluate( rkp_t *, int, double * );

double rkp_term_evaluate_scaled( rkp_t *, int, double * );

double rkp_term_evaluate_scaled_const( rkp_t *, int, double *, double );

double rkp_term_evaluate_node( rkp_t *, int, int );

double rkp_term_evaluate_node_scaled( rkp_t *, int, int );

double rkp_term_evaluate_node_scaled_const( rkp_t *, int, int );

double rkp_wavelet_term_evaluate( rkp_t *, int, int, double *, double );

double rkp_wavelet_term_evaluate_scaled( rkp_t *, int, int, double * );

double rkp_wavelet_term_evaluate_scaled_const( rkp_t *, int, int, double *, double );

double rkp_wavelet_term_evaluate_node( rkp_t *, int, int, int );

double rkp_wavelet_term_evaluate_node_scaled( rkp_t *, int, int, int );

double rkp_wavelet_term_evaluate_node_scaled_const( rkp_t *, int, int, int );

int rkp_grad_term_evaluate( rkp_t *, int, double *, double * );

double mls_evaluate( mls_t *, double * );

double rkp_evaluate( rkp_t *, double * );

double rkp_evaluate_scaled( rkp_t *, double * );

double rkp_evaluate_scaled_const( rkp_t *, double *, double );

double rkp_evaluate_node( rkp_t *, int );

double rkp_evaluate_node_scaled( rkp_t *, int );

double rkp_evaluate_node_scaled_const( rkp_t *, int );

double rkp_wavelet_evaluate( rkp_t *, int, double *, double );

double rkp_wavelet_evaluate_scaled( rkp_t *, int, double *, double );

double rkp_wavelet_evaluate_node( rkp_t *, int, int );

double rkp_wavelet_evaluate_node_scaled( rkp_t *, int, int );

double rkp_evaluate_inner_prod( rkp_t *, int, int );

double rkp_evaluate_poisson_inner_prod( rkp_t *, int, int );

double rkp_evaluate_functional( rkp_t *, int, double (*)(int,double*) );

#endif

