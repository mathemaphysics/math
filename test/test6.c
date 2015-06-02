#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* Domain is 0 to 1 */
double dist( double x_in )
{
	return sqrt( 2.0 * x_in );
}

double distp( double x_in )
{
	return x_in;
}

double func( double x_in )
{
	double x = ( 2.0 * x_in * x_in - 3.0 * x_in + 3.0 ) / 8.0;
        return x;
}

int main( int argc, char **argv )
{
	int i,k;
	double x,r,sum;
	FILE *fp = fopen( "nonuniform.out", "w" );

	srand((unsigned)time(0));
	sum = 0.0;
	for(i=0;i<100000;i++)
	{
		x = 0.5 * (double) rand() / (double) RAND_MAX; /* Integrate over 0 to 1 */
		fprintf( fp, "%d %10.5f\n", i, dist( x ) );
		sum = sum + func(dist(x)) / distp(dist(x)); //func( dist( x ) ) / distp( dist( x ) );
	}
	printf( "Sum = %14.10f\n", sum / (double) i * 0.5 );
	fclose( fp );

	return 0;
}

