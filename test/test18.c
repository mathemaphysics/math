#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <string.h>
#include "nurbs.h"
#include "points.h"
#include "nuclei.h"
#include "parse.h"

typedef struct
{
	int dim;
	double *corner;
	double side;
} square_t;

typedef struct
{
	int dim;
	double *center;
	double radius;
} circle_t;

typedef struct
{
	int dim;
	double *corner;
	double *sides;
} rectangle_t;

int square_init( square_t *obj_in, int dim_in )
{
	int i;
	obj_in->dim = dim_in;
	obj_in->corner = (double*) malloc( dim_in * sizeof(double) );

	return 0;
}

int square_free( square_t *obj_in )
{
	free( obj_in->corner );
}

int circle_init( circle_t *obj_in, int dim_in )
{
	int i;
	obj_in->dim = dim_in;
	obj_in->center = (double*) malloc( dim_in * sizeof(double) );

	return 0;
}

int circle_free( circle_t *obj_in )
{
	free( obj_in->center );
}

int rectangle_init( rectangle_t *obj_in, int dim_in )
{
	int i;
	obj_in->dim = dim_in;
	obj_in->corner = (double*) malloc( dim_in * sizeof(double) );
	obj_in->sides = (double*) malloc( dim_in * sizeof(double) );

	return 0;
}

int rectangle_free( rectangle_t *obj_in )
{
	free( obj_in->corner );
	free( obj_in->sides );
}

int indicator_square( double *x_in, void *args_in )
{
	int i;
	square_t *obj = (square_t*) args_in;

	for(i=0;i<obj->dim;i++)
		if( x_in[i] < obj->corner[i] || x_in[i] > obj->corner[i] + obj->side )
			return 0;
	return 1;
}

int indicator_circle( double *x_in, void *args_in )
{
	double x[2] = { 0.0, 0.0 };
        double sum = 0.0;
        sum = ( x_in[0] - x[0] ) * ( x_in[0] - x[0] ) + ( x_in[1] - x[1] ) * ( x_in[1] - x[1] );
        sum = sqrt( sum );
	if( sum < 1.0 )
		return 1;
	else
		return 0;
}

int indicator_rectangle( double *x_in, void *args_in )
{
	int i;
	rectangle_t *obj = (rectangle_t*) args_in;

	for(i=0;i<obj->dim;i++)
		if( x_in[i] < obj->corner[i] || x_in[i] > obj->corner[i] + obj->sides[i] )
			return 0;
	return 1;
}

double uniform_density( double *x_in, void *args_in )
{
	return 0.4;
}

double nucleus_density( double *x_in, void *args_in )
{
	nuclei_t *nuc = (nuclei_t*) args_in;
	return nucleus_length( nuc, x_in );
}

double distribution( double *x_in, void *args_in )
{
	double x[2] = { 2.5, 2.5 };
	double sum = 0.0;
	sum = ( x_in[0] - x[0] ) * ( x_in[0] - x[0] ) + ( x_in[1] - x[1] ) * ( x_in[1] - x[1] );
	sum = sqrt( sum );
	return exp( -0.5 * sum * sum );
}

double nucleus_distribution( double *x_in, void *args_in )
{
	nuclei_t *nuc = (nuclei_t*) args_in;
	return nucleus_potential( nuc, x_in );
}

int main( int argc, char **argv )
{
	int i,n,m;
	int dim = 2;
	int np = 600;
	int nq = 50000;
	int nstep = 10000;
	int vmod = 500;
	double *eps;
	double *sig;
	double *lim;
	double *pts;
	double rcut = 0.4;
	double temp = 1.5;
	double dx = 0.2;
	int (*ind)(double*,void*);
	double (*dns)(double*,void*);
	char buf[1024],bug[1024];
	char *tok[1024],*tol[1024];
	void *idata = malloc( 1024 );
	void *ddata = malloc( 1024 );
	char fname[1024];
	nuclei_t nuc;
	nuc.params = NULL;

	/* Switches */
	int b_post = 0; /* Do post processing */
	int b_load = 0;

	/* Take in the options */
	int optc;
	while( ( optc = getopt( argc, argv, "p:q:d:n:l:r:t:x:v:i:I:D:P" ) ) != -1 )
        {
                switch( optc )
                {
			case 'p':
				np = atoi( optarg );
				eps = (double*) malloc( np * sizeof(double) );
				sig = (double*) malloc( np * sizeof(double) );
				break;
			case 'q':
				nq = atoi( optarg );
				break;
			case 'd':
				dim = atoi( optarg );
				lim = (double*) malloc( 2 * dim * sizeof(double) );
				break;
			case 'n':
				nstep = atoi( optarg );
				break;
			case 'l':
				strncpy( buf, optarg, 1024 );
				n = parse_stokenize( buf, tok, "," );
				lim = (double*) malloc( n * sizeof(double) );
				for(i=0;i<n;i++)
					lim[i] = atof( tok[i] );
				break;
			case 'r':
				rcut = atof( optarg );
				break;
			case 't':
				temp = atof( optarg );
				break;
			case 'x':
				dx = atof( optarg );
				break;
			case 'v':
				vmod = atoi( optarg );
				break;
			case 'i':
				b_load = 1;
				load_points( optarg, &dim, &np, &pts, 1 );
				break;
			case 'f':
				
				break;
			case 'I':
				strncpy( buf, optarg, 1024 );
				n = parse_stokenize( buf, tok, ":" );
				if( n < 2 )
				{
					fprintf( stderr, "Bad domain arguments. Exiting.\n" );
					return 0;
				}
				if( strcmp( tok[0], "square" ) == 0 )
				{
					square_t *tmp = idata;
					strcpy( bug, tok[1] );
					m = parse_stokenize( bug, tol, "," );
					if( m < 2 )
					{
						fprintf( stderr, "Bad arguments to square. Exiting.\n" );
						return 0;
					}
					square_init( tmp, m - 1 );
					tmp->side = atof( tol[0] );
					for(i=0;i<tmp->dim;i++)
						tmp->corner[i] = atof( tol[i+1] );
					ind = &indicator_square;
				}
				if( strcmp( tok[0], "circle" ) == 0 )
				{
					circle_t *tmp = idata;
					
				}
				if( strcmp( tok[0], "rectangle" ) == 0 )
				{
					rectangle_t *tmp = idata;
					strcpy( bug, tok[1] );
					m = parse_stokenize( bug, tol, "," );
					if( m % 2 != 0 ) /* Should have an even number of arguments */
					{
						fprintf( stderr, "Bad arguments to rectangle. Exiting.\n" );
						return 0;
					}
					rectangle_init( tmp, m / 2 );
					for(i=0;i<m/2;i++)
					{
						tmp->sides[i] = atof( tol[i] );
						tmp->corner[i] = atof( tol[m/2+i] );
					}
					ind = &indicator_rectangle;
				}
				break;
			case 'D':
				load_nuclei( optarg, &nuc, 0 );
				nuclei_print( &nuc );
				dns = &nucleus_density;
				ddata = &nuc;
				break;
			case 'P':
				b_post = 1;
				break;
			case '?':
			default:
				fprintf( stderr, "I don't understand the jibberish you are spouting. Goodbye.\n" );
				return 0;
				break;
		}
	}

	/* Generate the initial distribution */
	if( !b_load )
	{
		generate_cloud_macqueen( dim, np, nq, ind, &nucleus_distribution, lim, &pts, idata, ddata );
		write_points( "cloud.out.pre", dim, np, pts, 0 );
		write_points( "cloud.out.pre.plot", dim, np, pts, 1 );
	}
	else
	{
		/* Add np points to the loaded distribution */
	}
	eps = (double*) malloc( np * sizeof(double) );
	sig = (double*) malloc( np * sizeof(double) );
	smooth_cloud_bubble_mc( dim, np, pts, ind, dns, dx, eps, sig, temp, rcut, nstep, vmod, idata, ddata );
	write_points( "cloud.out", dim, np, pts, 0 );
	write_points( "cloud.out.plot", dim, np, pts, 1 ); /* Output for view without header line */
	if( b_post )
	{
		smooth_cloud_macqueen( dim, np, nq, ind, &nucleus_distribution, lim, pts, idata, ddata );
		write_points( "cloud.out.post", dim, np, pts, 0 );
		write_points( "cloud.out.post.plot", dim, np, pts, 1 );
	}

	free( eps );
	free( sig );

	return 0;
}

