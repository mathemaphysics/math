#ifndef	CGRAD_H
#define CGRAD_H

void copy( int, double *, double * );
void project( int, double *, int, double * );
void bicgstab( int, int, double *, double *, double *, int, double, int, int * );
void bicgstabp( int, int, double *, double *, int, double *, double *, int, double, int, int * );
void ulanczos( int, int, double *, double *, double *, double *, double *, double *, int * );

#endif
