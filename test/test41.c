#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "legendre.h"

int main( int argc, char **argv )
{
	int i,j;
	int div = 200;
	double ll = 0.0;
	double ul = 30.0;
	double x;
	double step = ( ul - ll ) / (double) div;
	FILE *fp = fopen( "laguerre.out", "w" );

	/* Set this shit */
	const int en = 2;
	const int el = 0;
	const double a0 = 0.529;

	/* Output some other shit */
	for(j=0;j<=div;j++)
	{
		x = radial_wave_evaluate( en, el, (double) j * step + ll, a0 );
		fprintf( fp, "%10.4f %10.4f\n", ((double) j * step + ll) / a0, /* pow( ((double) j * step + ll) / a0, 2.0 ) * */ x );
	}

	/* Close and clean up */
	fclose( fp );

	return 0;
}

