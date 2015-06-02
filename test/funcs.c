#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double exp_func( int dim_in, int drv_in, double dlt_in, double *x_in )
{
	int i,n=1;
	double r,rho;
	const double a0 = 0.529; /* Angstroms */
	double tot;

	if( drv_in <= 0 )
		tot = 1.0;
	else
		tot = pow( -1.0 / (double) n / a0, (double) drv_in );

	r = 0.0;
	for(i=0;i<dim_in;i++)
		r += x_in[i] * x_in[i];
	r = sqrt( r );

	rho = 2.0 * r / (double) n / a0;

	return tot * exp( -0.5 * rho );
}
