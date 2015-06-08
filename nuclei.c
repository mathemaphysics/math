#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "config.h"

typedef struct
{
	int dim;
	int np;
	double *pts;
	double *dlt;
	double *params;
	double *chg;
} nuclei_t;

int nuclei_init( nuclei_t *obj_in, int dim_in, int np_in )
{
	obj_in->dim = dim_in;
	obj_in->np = np_in;
	obj_in->pts = (double*) malloc( np_in * dim_in * sizeof(double) );
	obj_in->dlt = (double*) malloc( np_in * sizeof(double) );
	obj_in->params = (double*) malloc( 2 * np_in * sizeof(double) );
	obj_in->chg = (double*) malloc( np_in * sizeof(double) );

	return 0;
}

void nuclei_free( nuclei_t *obj_in )
{
	free( obj_in->pts );
	free( obj_in->dlt );
	free( obj_in->params );
	free( obj_in->chg );
}

void nuclei_print( nuclei_t *obj_in )
{
	int i,j;

	fprintf( stderr, " Nuclei: dim = %d np = %d\n", obj_in->dim, obj_in->np );
	for(i=0;i<obj_in->np;i++)
	{
		fprintf( stderr, " %d: ", i );
		for(j=0;j<obj_in->dim;j++)
			fprintf( stderr, "%15.7f", obj_in->pts[i*obj_in->dim+j] );
		fprintf( stderr, " dlt = %15.7f ", obj_in->dlt[i] );
		fprintf( stderr, "params = %15.7f%15.7f ", obj_in->params[2*i+0], obj_in->params[2*i+1] );
		fprintf( stderr, "charge = %15.7f\n", obj_in->chg[i] );
	}
}

double gaussian( int dim_in, double std_in, double *ctr_in, double *x_in )
{
	int i;
	double sum = 0.0;
	for(i=0;i<dim_in;i++)
		sum += pow( x_in[i] - ctr_in[i], 2.0 );
	return exp( -1.0 * sum / std_in / std_in );
}

int load_nuclei( char *fn_in, nuclei_t *nuc_in, int fmt_in )
{
    int i,j,n,m;
    char buf[1024],*tok[1024];
    FILE *fp;

    fp = fopen( fn_in, "r" );
    if( fp == NULL )
        return -1;
    n = parse_read_line( fp, buf );
    if( n <= 0 )
        return -1;
    m = parse_stokenize( buf, tok, " \t\n" );
    if( m != 2 )
        return -1;
    nuc_in->dim = atoi( tok[0] );
    nuc_in->np = atoi( tok[1] );
    nuc_in->pts = (double*) malloc( nuc_in->dim * nuc_in->np * sizeof(double) );
    nuc_in->dlt = (double*) malloc( nuc_in->np * sizeof(double) );
    nuc_in->params = (double*) malloc( 2 * nuc_in->np * sizeof(double) );
    nuc_in->chg = (double*) malloc( nuc_in->np * sizeof(double) );
    if( nuc_in->pts == NULL )
        return -1;
    for(i=0;i<nuc_in->np;i++)
    {
        n = parse_read_line( fp, buf );
        if( n <= 0 )
            continue;
        m = parse_stokenize( buf, tok, " \t\n" );
        if( m < nuc_in->dim + 4 )
            continue;
        for(j=0;j<nuc_in->dim;j++)
            nuc_in->pts[i*nuc_in->dim+j] = atof( tok[j] );
        nuc_in->chg[i] = atof( tok[nuc_in->dim] );
        nuc_in->dlt[i] = atof( tok[nuc_in->dim+1] );
        nuc_in->params[2*i+0] = atof( tok[nuc_in->dim+2] );
        nuc_in->params[2*i+1] = atof( tok[nuc_in->dim+3] );
    }

    return 0;
}

/**
 * Calculate the potential of an electron at x_in interacting with
 * all nuclei; FIXME: Need to implement charge parameters
 */
double nuclei_potential( nuclei_t *obj_in, double *x_in )
{
	int i,j;
	double r,sum;

	sum = 0.0; /* Don't forget this dumbass */
	for(i=0;i<obj_in->np;i++)
	{
		r = 0.0;
		for(j=0;j<obj_in->dim;j++)
			r += pow( x_in[j] - obj_in->pts[i*obj_in->dim+j], 2.0 );
		r = sqrt( r );
		sum -= obj_in->chg[i] / r; /* For now assume charge is one */
	}
	return sum;
}

/**
 * This builds a distribution function which is unity
 * at the given locations of nuclei
 */
double nucleus_potential( nuclei_t *obj_in, double *x_in )
{
	int i,j,idx;
	double sum = 0.0;
	double min;

	for(i=0;i<obj_in->np;i++)
	{
		sum = 0.0;
		for(j=0;j<obj_in->dim;j++)
			sum += pow( obj_in->pts[i*obj_in->dim+j] - x_in[j], 2.0 );
		if( i == 0 || sum < min )
			min = sum, idx = i;
	}
	sum = gaussian( obj_in->dim, obj_in->dlt[idx], obj_in->pts + idx * obj_in->dim, x_in );

	return sum;
}

/**
 * Returns a distances to be used for calculating mesh density;
 * specifically the distance is the van der Waals distance to use
 * at the point x_in
 */
double nucleus_length( nuclei_t *obj_in, double *x_in )
{ 
        int i,j,idx;
        double sum = 0.0;
        double min;
	double scl,ncl;
 
        for(i=0;i<obj_in->np;i++)
        {
                sum = 0.0;
                for(j=0;j<obj_in->dim;j++)
                        sum += pow( obj_in->pts[i*obj_in->dim+j] - x_in[j], 2.0 );
                if( i == 0 || sum < min )
                        min = sum, idx = i;
        }
	scl = obj_in->params[2*idx+0];
	ncl = obj_in->params[2*idx+1];
        sum = scl * ( 1.0 - gaussian( obj_in->dim, obj_in->dlt[idx], obj_in->pts + idx * obj_in->dim, x_in ) ) + ncl;

        return sum;
}

/**
 * Generate a point distribution for a set of nuclei.
 * @param dim_in Dimension of the space
 * @param nnuc_in The number of nuclei given as input
 * @param nuc_in The locations of the nuclei input
 * @param nppn_in The number of points per nucleus
 * @param pts_out The actual points generated as output
 * @return Returns 0 if success, otherwise an error
 */
int generate_nucleus_cloud( int dim_in, int nnuc_in, double *nuc_in, int nppn_in, double **pts_out, double step_in, int max_in )
{
	int i,j,k,m;
	double sum,vec[dim_in];

	/* Initialize all points to an assigned center */
	double *pts = (double*) malloc( nnuc_in * nppn_in * dim_in * sizeof(double) );
	for(i=0;i<nnuc_in;i++)
		for(j=0;j<nppn_in;j++)
			for(k=0;k<dim_in;k++)
				pts[i*nppn_in*dim_in+j*dim_in+k] = nuc_in[i*dim_in+k];

	/* Now start taking randomly oriented steps */
	srand((unsigned)time(0));
	for(i=0;i<max_in;i++)
	{
		for(j=0;j<nnuc_in;j++)
		{
			for(k=0;k<nppn_in;k++)
			{
				/* First generate a random direction */
				for(m=0;m<dim_in;m++)
					vec[m] = ( rand() > RAND_MAX / 2 ? 1.0 : -1.0 ) * (double) rand() / (double) RAND_MAX;
				sum = 0.0;
				for(m=0;m<dim_in;m++)
					sum += vec[m] * vec[m];
				sum = sqrt( sum );
				for(m=0;m<dim_in;m++)
					vec[m] /= sum;
				for(m=0;m<dim_in;m++)
					pts[j*nppn_in*dim_in+k*dim_in+m] += ( rand() > RAND_MAX / 2 ? 1.0 : -1.0 ) * step_in * vec[m];
			}
		}
	}

	/* Output the points */
	write_points( "cloud.out", dim_in, nnuc_in * nppn_in, pts, 1 );

	return 0;
}

int smooth_nucleus_cloud( int dim_in, int np_in, double *pts_in )
{
	
}

/**
 * This function is intended for use as a basis function for PU methods
 * to approximate electron density near a zero-angular-momentum nuclear
 * singularity
 */
double singular_wavefunction( int dim_in, double *x_in )
{
	int i;
	double sum = 0.0;

	for(i=0;i<dim_in;i++)
		sum += x_in[i] * x_in[i];
	sum = sqrt( sum );
	return exp( -sum );
}

