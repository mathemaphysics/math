#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mls.h"
#include "linalg.h"

double cubic( int dim_in, double a_in, double *x_in )
{
	double z = fabs( x_in[0] ) / a_in;
	if( z >= 0.0 && z < 1.0 )
		return ( 2.0 / 3.0 - z * z * ( 1.0 - z / 2.0 ) );
	else if( z >= 1.0 && z <= 2.0 )
		return 1.0 / 6.0 * pow( 2.0 - z, 3.0 );
	else
		return 0.0;
}

double gaussian( int dim_in, double a_in, double *x_in )
{
	int i;
	double val = 0.0;
	if( x_in[0] < -4.0 * a_in || x_in[0] > 4.0 * a_in )
		return 0.0;
	for(i=0;i<dim_in;i++)
		val += x_in[i] * x_in[i] / a_in / a_in;
	return exp( -0.5 * val ) / a_in / sqrt( 2.0 * M_PI );
}

double parabola( int dim_in, double a_in, double *x_in )
{
	if( x_in[0] < -1.0 * a_in || x_in[0] > 1.0 * a_in )
		return 0.0;
	else
		return 1 - x_in[0] * x_in[0] / a_in / a_in;
}

double step( int dim_in, double a_in, double *x_in )
{
	if( x_in[0] < 0 )
		return -1.0 * x_in[0];
	else
		return 1.0 * x_in[0];
}

int main( int argc, char **argv )
{
	rkp_t ml;
	int i,j,pdim=binomial(2+1,1);
	double pts[21*1] = { 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0 };
	double vals[21] = { 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0 };
	for(i=0;i<21;i++)
		vals[i] = vals[i] * vals[i] - 0.5 * vals[i];
	double xx,xt[1] = { 0.2 };
	double dlts[21];
	for(i=0;i<21;i++)
		dlts[i] = 0.09;
	//dlts[10] = 0.21;
	double val;
	double gqw[21];
	for(i=0;i<21;i++)
		gqw[i] = 1.0;
	double temp[300];
	for(i=0;i<300;i++)
		temp[i] = 0.0;

	rkp_init( &ml, 21, 1, 2, pts, dlts, gqw, vals, cubic, 2.0 );
	rkp_basis_generate( &ml );
	rkp_wavelet_basis_generate( &ml );

	/* Plot the basis */
	FILE *fp = fopen( "test5.out", "w" );
	double stp = 1.0 / 300.0;
	for(j=0;j<21;j++)
	for(i=0;i<21;i++)
		xt[0] = pts[i],
		val = rkp_term_evaluate_node( &ml, j, i ),
		fprintf( fp, "%f %f\n", xt[0], val ),
		temp[i] += val;
	fclose( fp );

	for(i=0;i<300;i++)
		temp[i] = 0.0;
	fp = fopen( "test5.out.2", "w" );
	for(j=0;j<21;j++)
	for(i=0;i<300;i++)
		xt[0] = (double) i * stp,
		val = rkp_term_evaluate( &ml, j, xt ),
		fprintf( fp, "%f %f\n", xt[0], val ),
		temp[i] += val;
	fclose( fp );

	fp = fopen( "test5.out.sum", "w" );
	for(i=0;i<300;i++)
		fprintf( fp, "%f %f\n", (double) i * stp, temp[i] );
	fclose( fp );

	/* Plot the wavelets */
	for(i=0;i<300;i++)
		temp[i] = 0.0;
	int order = 1;
	fp = fopen( "test6.out", "w" );
	order = 0;
	j = 1;
	for(j=0;j<21;j++)
	for(i=0;i<21/*300*/;i++)
// 		xt[0] = (double) i * stp,
		xt[0] = pts[i],
// 		val = rkp_wavelet_term_evaluate( &ml, j, order, xt ),
		val = rkp_wavelet_term_evaluate_node( &ml, j, order, i ),
		temp[i] += val,
		fprintf( fp, "%f %f\n", xt[0], val );
	fclose( fp );

	fp = fopen( "test6.out.2", "w" );
	for(i=0;i<21;i++)
		xt[0] = pts[i],
		val = rkp_wavelet_term_evaluate_node( &ml, 10, order, i ),
		fprintf( fp, "%f %f\n", xt[0], val );
	fclose( fp );

	fp = fopen( "test6.out.3", "w" );
	order = 2;
	for(j=0;j<21;j++)
	for(i=0;i<300;i++)
		xt[0] = (double) i * stp,
		val = rkp_wavelet_term_evaluate( &ml, j, order, xt, 0.1 ),
		fprintf( fp, "%f %f\n", xt[0], val ),
		temp[i] += val;
	fclose( fp );
	
	fp = fopen( "test6.out.sum", "w" );
	for(i=0;i<300;i++)
		fprintf( fp, "%f %f\n", (double) i * stp, temp[i] );
	fclose( fp );
	
	/* Output the plot itself */
	fp = fopen( "stuff-plot", "w" );
	for(i=0;i<300;i++)
		xt[0] = (double) i * stp, fprintf( fp, "%f %f\n", xt[0], rkp_evaluate( &ml, xt ) );
	fclose( fp );

	fp = fopen( "stuff-plot-2", "w" );
	for(i=0;i<300;i++)
	{
		xt[0] = (double) i * stp;
		xx = 0.0;
		xx += rkp_evaluate( &ml, xt );
		xx += rkp_wavelet_evaluate( &ml, 1, xt, 0.02 );
		fprintf( fp, "%f %f\n", xt[0], xx );
	}
	fclose( fp );

	return 0;
}

