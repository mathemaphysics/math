/**
 * This file is a basic template for use in writing any software which
 * loads all of its potentially complicated configuration information
 * from a single primary file. Other configuration routines may be
 * performed as a result of directives inside this primary configuration
 * file (which is loaded via a command line switch).
 */

/* Basic stuff we always use */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>

/* Includes for configuration */
#include "trees.h"
#include "parse.h"
#include "cfgread.h"

/* Includes for specific tasks for this template */
#include "points.h"
#include "punity.h"

#define TOKEN_BUFFER_LENGTH 1024

int main( int argc, char **argv )
{
	int i,j,optc,ret;
	double x,y;
	symtab_t sym;
	symbol_t *s;

	/* General parsing variables */
	int n;
	char buf[TOKEN_BUFFER_LENGTH],*tok[TOKEN_BUFFER_LENGTH];

	/* Points data structures */
	int dim,npts,rdn,nbp;
	double *pts,*dlt;
	punity_t pt;
	shape_t dm;

	/* Grid variables for plotting */
	int *dx;
	long *size,*index,*tindex;
	double sum,*sx,*x0,*wx,*dd;

	/* Quadrature variables */
	int quadn;
	double *ctr,*qbox,*nqbox,*qpts,*qwts,*qx,qw,rad;

	/* Switches to indicate boolean things */
	int b_config_loaded = 0;
	int b_points_loaded = 0;
	int b_rdn_set = 0;
	int b_domain_set = 0;
	int b_grid_set_dx = 0;
	int b_grid_set_x0 = 0;
	int b_grid_set_wx = 0;
	int b_quad_loaded = 0;
	
	/* Floating points */
	double f_bndry_tol = 0.05;

	/* Announce */
	fprintf( stderr, "\nPU Template VERSION 0.1\n\n" );

	/* Load the options */
	while( ( optc = getopt( argc, argv, "C:p:r:" ) ) != -1 )
        {
                switch( optc )
                {
			case 'C':
				/* Load symbols */
				symtab_init( &sym, 1024 );
				ret = cfgread_load_symbols_f( optarg, &sym );
				if( ret == 0 )
				{
					b_config_loaded = 1;
					fprintf( stderr, " * Loaded following symbols from \"%s\":\n\n", optarg );
					symtab_print( &sym );
					fprintf( stderr, "\n" );
				}
				else
					fprintf( stderr, " * ERROR: Failed to load configuration variables from file \"%s\": Returned %d\n\n", optarg, ret );
				break;
			case 'p':
				if( load_points( optarg, &dim, &npts, &pts, 0 ) == 0 )
				{
					fprintf( stderr, " * Loaded %d points in %d dimensions from file \"%s\"\n\n", npts, dim, optarg );
					b_points_loaded = 1;
					dlt = (double*) malloc( npts * sizeof(double) );
					if( dlt == NULL )
					{
						fprintf( stderr, " * ERROR: Could not allocate space for dilation factors. Exiting.\n\n" );
						return 0;
					}
				}
				else
					fprintf( stderr, " * ERROR: Failed loading points from file \"%s\"\n\n", optarg );
				break;
			case 'r':
				rdn = atoi( optarg );
				b_rdn_set = 1;
				fprintf( stderr, " * Set radial neighbor number to %d\n\n", rdn );
				break;
			default:
				fprintf( stderr, "I don't understand the jibberish you are spouting. Goodbye.\n" );
				return 0;
		}
	}

	/* Process all the configuration variables by loading them one at a time */
	if( b_config_loaded )
	{
		fprintf( stderr, " * Reading and executing configuration variables:\n\n" );
		if( !b_points_loaded )
		{
			s = symtab_lookup( &sym, "pointCloud" );
			if( s != NULL )
			{
				if( load_points( (char*) s->data, &dim, &npts, &pts, 0 ) == 0 )
				{
					fprintf( stderr, " ---> Loaded %d points in %d dimensions from file \"%s\"\n", npts, dim, (char*) s->data );
					b_points_loaded = 1;
					dlt = (double*) malloc( npts * sizeof(double) );
					if( dlt == NULL )
					{
						fprintf( stderr, " * ERROR: Could not allocate space for dilation factors. Exiting.\n\n" );
						return 0;
					}
				}
				else
					fprintf( stderr, " ---> ERROR: Failed loading points from file \"%s\"\n", (char*) s->data );
			}
		}
		if( !b_rdn_set )
		{
			s = symtab_lookup( &sym, "radialNeighbors" );
			if( s != NULL )
			{
				rdn = atoi( (char*) s->data );
				b_rdn_set = 1;
				fprintf( stderr, " ---> Loaded radial neighbors number %d\n", rdn );
			}
		}
		if( !b_domain_set )
		{
			s = symtab_lookup( &sym, "domain" );
			if( s != NULL )
			{
				load_domain( (char*) s->data, &dm, &ret );
				if( ret == 0 )
				{
					fprintf( stderr, " ---> Loaded domain information\n" );
					b_domain_set = 1;
				}
			}
		}
		s = symtab_lookup( &sym, "boundaryTolerance" );
		if( s != NULL )
		{
			f_bndry_tol = atof( (char*) s->data );
			fprintf( stderr, " ---> Loaded boundary tolerance of %15.7f\n", f_bndry_tol );
		}
		if( !b_points_loaded )
			fprintf( stderr, " ---> ERROR: Cannot load grid without loading points; need dimension\n" );
		else
		{
			dx = (int*) malloc( dim * sizeof(int) );
			x0 = (double*) malloc( dim * sizeof(double) );
			wx = (double*) malloc( dim * sizeof(double) );
			dd = (double*) malloc( dim * sizeof(double) );
			sx = (double*) malloc( dim * sizeof(double) );
			size = (long*) malloc( dim * sizeof(long) );
			index = (long*) malloc( dim * sizeof(long) );
			tindex = (long*) malloc( dim * sizeof(long) );
			s = symtab_lookup( &sym, "gridDivision" );
			if( s != NULL )
			{
				strncpy( buf, (char*) s->data, TOKEN_BUFFER_LENGTH );
				n = parse_stokenize( buf, tok, "," );
				if( n < dim )
				{
					fprintf( stderr, " ---> ERROR: Too few entries for gridDivision\n" );
					for(i=0;i<n;i++)
						dx[i] = atoi( tok[i] );
					for(i=n;i<dim;i++)
						dx[i] = dx[n-1]; /* Fill it in with the last entry given */
				}
				else
					for(i=0;i<dim;i++)
						dx[i] = atoi( tok[i] );
				b_grid_set_dx = 1;
				fprintf( stderr, " ---> Set grid: dx = " );
				for(i=0;i<dim;i++)
					fprintf( stderr, "%5d", dx[i] );
				fprintf( stderr, "\n" );
			}
			s = symtab_lookup( &sym, "gridOrigin" );
			if( s != NULL )
			{
				strncpy( buf, (char*) s->data, TOKEN_BUFFER_LENGTH );
				n = parse_stokenize( buf, tok, "," );
				if( n < dim )
				{
					fprintf( stderr, " ---> ERROR: Too few entries for gridOrigin\n" );
					for(i=0;i<n;i++)
						x0[i] = atof( tok[i] );
					for(i=n;i<dim;i++)
						x0[i] = x0[n-1];
				}
				else
					for(i=0;i<dim;i++)
						x0[i] = atof( tok[i] );
				b_grid_set_x0 = 1;
				fprintf( stderr, " ---> Set grid: x0 = " );
				for(i=0;i<dim;i++)
					fprintf( stderr, "%12.5f", x0[i] );
				fprintf( stderr, "\n" );
			}
			s = symtab_lookup( &sym, "gridDimensions" );
			if( s != NULL )
			{
				strncpy( buf, (char*) s->data, TOKEN_BUFFER_LENGTH );
				n = parse_stokenize( buf, tok, "," );
				if( n < dim )
				{
					fprintf( stderr, " ---> ERROR: Too few entries for gridDimensions\n" );
					for(i=0;i<n;i++)
						wx[i] = atof( tok[i] );
					for(i=n;i<dim;i++)
						wx[i] = wx[n-1];
				}
				else
					for(i=0;i<dim;i++)
						wx[i] = atof( tok[i] );
				b_grid_set_wx = 1;
				fprintf( stderr, " ---> Set grid: wx = " );
				for(i=0;i<dim;i++)
					fprintf( stderr, "%12.5f", wx[i] );
				fprintf( stderr, "\n" );
			}
			for(i=0;i<dim;i++)
				dd[i] = wx[i] / (double) ( dx[i] - 1 );
		}
		if( !b_quad_loaded )
		{
			s = symtab_lookup( &sym, "quadrature" );
			if( s != NULL )
			{
				load_quadrature( (char*) s->data, &quadn, &qpts, &qwts, &ret );
				if( ret == 0 )
				{
					fprintf( stderr, " ---> Loaded quadrature rule with %d points from \"%s\"\n", quadn, (char*) s->data );
					b_quad_loaded = 1;
				}
				else
					fprintf( stderr, " ---> ERROR: Problem loading quadurature file in \"%s\": Return value %d\n", (char*) s->data, ret );
			}
		}
		fprintf( stderr, "\n" );
	}

	/* If no points loaded then say goodbye, come again */
	if( !b_points_loaded )
	{
		fprintf( stderr, " * ERROR: No point cloud loaded. Exiting.\n\n" );
		return 0;
	}
	if( !b_rdn_set )
	{
		fprintf( stderr, " * ERROR: No radial neighbor value set. Exiting.\n\n" );
		return 0;
	}
	if( !b_domain_set )
	{
		fprintf( stderr, " * ERROR: No domain indicator loaded. Exiting.\n\n" );
		return 0;
	}
	if( !b_quad_loaded )
	{
		fprintf( stderr, " * ERROR: No quadrature rule loaded. Exiting.\n\n" );
		return 0;
	}

	/* Allocate stuff with sizes given in data just loaded */
	ctr = (double*) malloc( dim * sizeof(double) );
	qbox = (double*) malloc( dim * dim * sizeof(double) );
	nqbox = (double*) malloc( dim * dim * sizeof(double) );
	qx = (double*) malloc( dim * sizeof(double) );

	/* Initialize everything */
	fprintf( stderr, " * Generating points supports\n\n" );
	generate_cloud_supports_min( dim, npts, pts, rdn, 1.05, dlt, 0.001 );
	punity_init( &pt, dim, npts, pts, dlt, &cubic_window, &cubic_window_deriv );

	/* Mark the boundary nodes */
	nbp = 0;
        for(i=0;i<pt.npts;i++)
        {
                boundary_indicator( pt.pts + i * pt.dim, &dm, f_bndry_tol, &ret );
                if( ret == 1 )
                        pt.bdry[i] = 1, ++nbp;
        }
        fprintf( stderr, " * Number of boundary points is %d\n\n", nbp );

	/* Quadrature loop prototype here */
	double lx1[2] = { 0.5, 0.0 };
	double lx2[2] = { 0.0, 0.5 };
	double lr1 = sqrt(2.0) * 0.5 * 1.5;
	double lr2 = sqrt(2.0) * 0.5 * 1.5;
	if( sphere_intersection( dim, lx1, lr1, lx2, lr2, ctr, &rad, qbox ) == 1 )
	{
		for(i=0;i<dim;i++)
		{
			sum = 0.0;
			for(j=0;j<dim;j++)
				sum += qbox[i*dim+j] * qbox[i*dim+j];
			sum = sqrt( sum );
			for(j=0;j<dim;j++)
				nqbox[i*dim+j] = qbox[i*dim+j] / sum;
		}
		for(i=0;i<dim;i++)
			index[i] = 0, size[i] = quadn - 1;
		sum = 0.0;
		do
		{
			lens_gauss_point( dim, lx1, lr1, lx2, lr2, rad, nqbox, index, qpts, qwts, qx, &qw );
			domain_indicator( qx, &dm, &ret );
			if( ret == 1 )
			{
				fprintf( stderr, "%15.7f%15.7f\n", qx[0], qx[1] );
				sum += qw * x;
			}
		}
		while( arraynext( (long) dim, size, index ) != -1 );
	}

	  //////////////////////////////////////
	 // Do your processing business here //
	//////////////////////////////////////

	/* Output grid iterated here */
	for(i=0;i<dim;i++)
		size[i] = (long) dx[i] - 1;
	for(i=0;i<dim;i++)
		index[i] = 0;
	do
	{
		/* Build the point at which to evaluate the eigenfunction */
		for(i=0;i<dim;i++)
			sx[i] = x0[i] + (double) index[i] * dd[i];

		/* Output the point for this line */
		for(i=0;i<dim;i++)
			fprintf( stdout, "%15.7f", sx[i] );

		/* Now evaluate the function at that point */
		sum = 0.0;

		  /////////////////////////////////////////////////////////////
		 // This is where you put your code to output your function //
		/////////////////////////////////////////////////////////////

		/* Ouptut the values consecutively */
		fprintf( stdout, "%15.7f", sum );

		/* Close the line */
		fprintf( stdout, "\n" );

		/* Check to see if another newline is needed to maintain gnuplot's preferred format */
		for(i=0;i<dim;i++)
			tindex[i] = index[i];
		arraynext( (long) dim, size, tindex );
		for(i=0;i<dim;i++)
		{
			if( index[i] == dx[i] - 1 && tindex[i] != dx[i] - 1 )
			{
				fprintf( stdout, "\n" );
				break;
			}
		}
	}
	while( arraynext( (long) dim, size, index ) != -1 );

	/* Don't forget to clean up */
	free( pts );
	free( dlt );
	free( dx ); free( x0 ); free( wx ); free( dd );
	free( size ); free( index ); free( tindex );
	free( dm.orig ); free( dm.params );
	free( sx );
	free( ctr ); free( qbox ); free( nqbox );

	/* Free the partition of unity data structures */
	punity_free( &pt );

	return 0;
}

