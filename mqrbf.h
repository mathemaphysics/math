#ifndef MQRBF_H
#define MQRBF_H

typedef struct
{
	int np;
	int dim;
	double *pts;
	double *prm;
	double *vals;
	double *coeffs;
	double *vdm;
} mqrbf_t;

double mqf( int, double, double * );

int mqrbf_init( mqrbf_t *, int, int, double *, double *, double * );

int mqrbf_free( mqrbf_t * );

int mqrbf_generate( mqrbf_t * );

#endif

