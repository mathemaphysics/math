#ifndef NUCLEI_H
#define NUCLEI_H

typedef struct
{
	int dim;
	int np;
	double *pts;
	double *dlt;
	double *params;
	double *chg;
} nuclei_t;

double gaussian( int, double, double *, double * );

int nuclei_init( nuclei_t *, int, int );

void nuclei_free( nuclei_t * );

double nuclei_potential( nuclei_t *, double * );

double nucleus_potential( nuclei_t *, double * );

double nucleus_length( nuclei_t *, double * );

int load_nuclei( char *, nuclei_t *, int );

int generate_nucleus_cloud( int, int, double *, int, double **, double, int );

double singular_wavefunction( int, double * );

#endif

