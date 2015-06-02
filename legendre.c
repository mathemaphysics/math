#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "combinadic.h"
#include "config.h"

double legendre_evaluate( int deg_in, double x_in )
{
	int i,b;
	double val = 0.0;
	for(i=0;i<=deg_in;i++)
	{
		b = binomial( deg_in, i );
		val = val + (double) (b*b) * pow( x_in - 1.0, (double) (deg_in - i) ) * pow( x_in + 1.0, (double) i );
	}
	return val / pow( 2, (double) deg_in );
}

double laguerre_evaluate( int n_in, int alpha_in, double x_in )
{
	int i,k;
	double tmp,val = 0.0;

	k = 1;
	val = (double) ( binomial( n_in + alpha_in, n_in  ) );
	for(i=1;i<=n_in;i++)
	{
		tmp = (double) ( ( i % 2 == 0 ? 1 : -1 ) * binomial( n_in + alpha_in, n_in - i ) );
		k = k * i;
		val += tmp / (double) k * pow( x_in, (double) i );
	}
	return val;
}

double laguerre_deriv_evaluate( int n_in, int alpha_in, int k_in, double x_in )
{
	double tot = ( k_in % 2 == 0 ) ? 1.0 : -1.0;
	return tot * laguerre_evaluate( n_in - k_in, alpha_in + k_in, x_in );
}

double radial_exp_deriv_evaluate( int k_in, double x_in )
{
	double tot = pow( 0.5, (double) k_in );
	return tot * exp( -0.5 * x_in );
}

double radial_wave_evaluate( int n_in, int l_in, double r_in, double a0_in )
{
	double rho = 2.0 * r_in / a0_in / (double) n_in;
	return pow( rho, l_in ) * exp( -0.5 * rho ) * laguerre_evaluate( n_in - l_in - 1, 2 * l_in + 1, rho );
}

double radial_wave_deriv_evaluate( int n_in, int l_in, int k_in, double r_in, double a0_in )
{
	
}

