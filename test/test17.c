#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "trees.h"

int compare_double( void *dbl1_in, void *dbl2_in )
{
	double dbl1 = *((double*)dbl1_in);
	double dbl2 = *((double*)dbl2_in);
	if( dbl1 > dbl2 )
		return 1;
	else if( dbl1 < dbl2 )
		return -1;
	else
		return 0;
}

int main()
{
	int i,j;
	int n = 50;
	double dbl;
	heap_t hp;

	heap_init( &hp, sizeof(double), &compare_double );
	srand((unsigned)time(0));
	for(i=0;i<n;i++)
	{
		dbl = ( rand() > RAND_MAX / 2 ? 1.0 : -1.0 ) * (double) rand() / (double) RAND_MAX;
		heap_insert( &hp, &dbl );
	}

	heap_print( &hp );

	for(i=0;i<hp.size;i++)
	{
		dbl = *((double*)hp.start);
		printf( "%15.7f", dbl );
		heap_delete( &hp, 0 );
	}
	printf( "\n" );
}

