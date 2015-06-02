#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "trees.h"

int main( int argc, char **argv )
{
	int i;
	long siz[3] = { 10, 10, 10 };
	long idx[3] = { 0, 0, 0 };

	do
	{
		printf( "%d %d %d\n", idx[0], idx[1], idx[2] );
	}
	while( arraynext( 3, siz, idx ) == 0 );

	return 0;
}

