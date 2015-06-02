#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "nurbs.h"
#include "points.h"

int main()
{
	int i,j,k,m,p;
	long dim = 2;
	long npm = 2;
	long deg[2] = { 3, 3 };
	double knots[38] = { 0.0, 0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.2, 1.2, 1.2,
                             0.0, 0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.2, 1.2, 1.2 };
	long iknots[3] = { 0, 19, 38 };
	double cpts[15*15*3];
	double wts[15*15];
	double x_in[2],x_out[2];
	for(i=0;i<15*15;i++)
		wts[i] = 1.0;
	srand((unsigned)time(0));
	for(p=0,i=0;i<15;i++)
	{
		for(j=0;j<15;j++)
		{
			cpts[p+0] = (double) i / (double) ( 15 - 1 ); // + 0.1 * ( (double) rand() / (double) RAND_MAX - 0.5 );
			cpts[p+1] = (double) j / (double) ( 15 - 1 ); // + 0.1 * ( (double) rand() / (double) RAND_MAX - 0.5 );
			p += 2;
		}
	}

	nurbs_t nb;
	nurbs_init( &nb, dim, npm, deg, iknots, knots, cpts, wts );

	for(i=0;i<70;i++)
	{
		x_in[0] = 1.195 * (double) i / (double) ( 70 - 1 );
		for(j=0;j<70;j++)
		{
			x_in[1] = 1.195 * (double) j / (double) ( 70 - 1 );
			nurbs_evaluate( &nb, x_in, x_out );
			for(m=0;m<2;m++)
				fprintf( stderr, "%15.7f", x_out[m] );
			fprintf( stderr, "\n" );
		}
		fprintf( stderr, "\n" );
	}

	return 0;
}
