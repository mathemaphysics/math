#ifndef SPARSE_H
#define SPARSE_H

double sparse_entry( long, long, long, long, long *, long *, double * );

int sparse_load( char *, long *, long *, long **, long **, double **, int );

int sparse_save( char *, long, long, long *, long *, double *, int );

void sparse_symbmm( long, long, long, long *, long *, long *, long *, long *, long *, long * );

void sparse_dgemm( long, long, long, long *, long *, double *, long *, long *, double *, long *, long *, double *, double * );

void sparse_dgemv( long, long, long *, long *, double *, long, double *, long, double * );

void sparse_dgemv_accum( long, long, long *, long *, double *, long, double *, long, double * );

void sparse_transp( char, long, long, long *, long *, double *, long *, long *, double * );

#endif
