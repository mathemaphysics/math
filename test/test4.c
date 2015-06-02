#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "legendre.h"

int main( int argc, char **argv )
{
	int i,j;
	int div = 100;
	double ll = -1.0;
	double ul = 1.0;
	double step = ( ul - ll ) / (double) div;
	FILE *fp = fopen( "legendre.out", "w" );

	for(i=0;i<8;i++)
		for(j=0;j<=div;j++)
			fprintf( fp, "%10.4f %10.4f\n", (double) j * step + ll, legendre_evaluate( i, (double) j * step + ll ) );
	fclose( fp );

	return 0;
}

