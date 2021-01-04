#ifndef MPIFMTOUT_H
#define MPIFMTOUT_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <mpi.h>
#include "codes.h"

typedef struct
{
	int op;
	char pre[32];
	char name[1024];
	int rank;
	char msg[4096];
} t_mpi_msg;

void mffprintf( int, char *, char *, int, int, int, FILE *, const char *, ... );

void mpisetmsg( t_mpi_msg *, int, char *, char *, int, char * );

void mpiservercom( int, int, int, MPI_Comm * );

int mpiprintserver( int, int, int, FILE * );

int mpiprint( int, int, int, t_mpi_msg * );

#endif
