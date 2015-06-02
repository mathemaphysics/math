#ifndef QUADRATURE_H
#define QUADRATURE_H

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

int quadrature_init( quadrature_t *, long, long, double * );

int quadrature_free( quadrature_t * );

int quadrature_set_nodes( quadrature_t *, double * );

int quadrature_set_nodes_standard( quadrature_t * );

int quadrature_generate( quadrature_t * );

#endif
