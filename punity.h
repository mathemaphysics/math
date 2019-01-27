#ifndef PUNITY_H
#define PUNITY_H

#include "points.h"

typedef double (*wfsd_t)(int,int,double,double*);

typedef struct
{
        int dim;
        int npts;
        double *pts;
        double *dlt;
        double (*wfs)(int,double,double*);
	double (*wfsd)(int,int,double,double*);
	double rmax;
	kdnode_t *kdt;
	char *bdry;
} punity_t;

double cubic_window( int, double, double * );

double cubic_window_deriv( int, int, double, double * );

double radial_deriv_evaluate( int, int *, double * );

double rational_radial_deriv_evaluate( int, int *, double *, double, int );

int punity_init( punity_t *, int, int, double *, double *, double (*)(int,double,double*), double (*)(int,int,double,double*) );

int punity_use_pbc( punity_t *, double * );

int punity_free( punity_t * );

int punity_neighbors( punity_t *, double *, int ** );

double punity_evaluate( punity_t *, int, double * );

double punity_evaluate_delta( punity_t *, int, double *, int );

double punity_window_deriv_evaluate( punity_t *, int, int *, double * );

double punity_window_deriv_evaluate_delta( punity_t *, int, int *, double *, int );

double punity_term_deriv_evaluate( punity_t *, int, int *, double *, int *, int );

double punity_term_deriv_evaluate_delta( punity_t *, int, int *, double *, int, int *, int );

double punity_eval_local_poly( punity_t *, int, int *, double * );

double punity_eval_local_poly_delta( punity_t *, int, int *, double *, int );

double punity_term_deriv_local_poly( punity_t *, int, int *, int *, double *, int *, int );

double punity_term_deriv_local_poly_delta( punity_t *, int, int *, int *, double *, int, int *, int );

double punity_term_deriv_local_radial( punity_t *, int, double (*)(int,int,double,double*), int *, double *, int *, int );

#endif

