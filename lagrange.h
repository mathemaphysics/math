#ifndef LAGRANGE_H
#define LAGRANGE_H

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

int lagrange_init( lagrange_t *, long, long, double * );

int lagrange_free( lagrange_t * );

int lagrange_set_nodes( lagrange_t *, double * );

int lagrange_set_nodes_standard( lagrange_t * );

int lagrange_generate( lagrange_t * );

double lagrange_evaluate( lagrange_t *, long, double * );

#endif

