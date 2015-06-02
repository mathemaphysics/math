#ifndef LINALG_H
#define LINALG_H

void dgemm_( char *, char *, int *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int * );

void dgemv_( char *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int * );

void daxpy_( int *, double *, int *, double *, int *, double * );

void dlacpy_( char *, int *, int *, double *, int *, double *, int * );

void dgepiv_( int *, double *, int *, int *, int * );

void dgelu_( char *, int *, int *, double *, int *, int *, int * );

void dgefb_( char *, int *, int *, double *, int *, int *, double *, int *, int * );

void dgedet_( char *, int *, double *, int *, int *, double *, int * );

void dgesv_( int *, int *, double *, int *, int *, double *, int *, int * );

void matrix_print( int *, int *, double * );

#endif

