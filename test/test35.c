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
#include "combinadic.h"
#include "trees.h"

/* Define this stuff for parsing things */
#define TOKEN_BUFFER_LENGTH 1024

int main( int argc, char **argv )
{
	int optc,ret;
	symtab_t sym;
	symbol_t *s;

	/* Points data structures */
	int dim,npts,rdn,nbp,idx,jdx;
	double *pts,*dlt;
	punity_t pt;
	shape_t dm;

	/* Grid variables */
	int *dx;
	double y,*x0,*wx,*dd,*xx;
	long *size,*index,*tindex;
	FILE *ff;

	/* Variables for loading grid variables */
	int i,j,n;
	char buf[TOKEN_BUFFER_LENGTH],*tok[TOKEN_BUFFER_LENGTH];

	/* Switches to indicate boolean things */
	int b_config_loaded = 0;
	int b_points_loaded = 0;
	int b_rdn_set = 0;
	int b_domain_set = 0;
	int b_grid_set_dx = 0;
	int b_grid_set_x0 = 0;
	int b_grid_set_wx = 0;

	/* Float values needing to be loaded */
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
					fprintf( stderr, " * ERROR: Failed to load variables from file \"%s\": Returned %d\n\n", optarg, ret );
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
		if( !b_points_loaded )
			fprintf( stderr, " ---> ERROR: Cannot load grid without loading points; need dimension\n" );
		else
		{
			dx = (int*) malloc( dim * sizeof(int) );
			x0 = (double*) malloc( dim * sizeof(double) );
			wx = (double*) malloc( dim * sizeof(double) );
			dd = (double*) malloc( dim * sizeof(double) );
			xx = (double*) malloc( dim * sizeof(double) );
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
		s = symtab_lookup( &sym, "boundaryTolerance" );
		if( s != NULL )
		{
			f_bndry_tol = atof( (char*) s->data );
			fprintf( stderr, " ---> Loaded boundary tolerance of %15.7f\n", f_bndry_tol );
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
		fprintf( stderr, " * ERROR: No domain indicator set. Exiting.\n\n" );
		return 0;
	}
	if( !b_grid_set_dx || !b_grid_set_x0 || !b_grid_set_wx )
	{
		fprintf( stderr, " * ERROR: No grid information loaded. Exiting.\n\n" );
		return 0;
	}

	/* Initialize everything here */
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

	/* Do your business here */
	for(i=0;i<dim;i++)
		size[i] = (long) dx[i] - 1;
	for(i=0;i<dim;i++)
		index[i] = 0;

	/* Pick out a boundary node to use as reference */
	srand((unsigned)time(0));
	jdx = rand() % nbp;
	for(i=0,n=0;i<npts;i++)
	{
		if( pt.bdry[i] == 1 )
		{
			if( n == jdx )
			{
				jdx = i;
				break;
			}
			++n;
		}
	}

	/* Plot the center */
	ff = fopen( "bdry.out", "w" );
	for(i=0;i<dim;i++)
		size[i] = (long) dx[i] - 1;
	for(i=0;i<dim;i++)
		index[i] = 0;
	do
	{
		for(i=0;i<dim;i++)
			xx[i] = x0[i] + (double) index[i] * dd[i];
		y = punity_evaluate_delta( &pt, jdx, xx, 1 );
		for(i=0;i<dim;i++)
			fprintf( ff, "%15.7f", xx[i] );
		fprintf( ff, "%15.7f\n", y );
		for(i=0;i<dim;i++)
			tindex[i] = index[i];
		arraynext( (long) dim, size, tindex );
		for(i=0;i<dim;i++)
		{
			if( index[i] == dx[i] - 1 && tindex[i] != dx[i] - 1 )
			{
				fprintf( ff, "\n" );
				break;
			}
		}
	}
	while( arraynext( (long) dim, size, index ) != -1 );
	fclose( ff );

	/* Start plot iteration */
	int drv[2] = { 0, 1 };
	for(j=0,n=0;j<npts;j++)
	{
		if( pt.bdry[j] == 0 )
		{
			/* Calculate how particle i is from boundary particle jdx */
			y = 0.0;
			for(i=0;i<dim;i++)
				y += pow( 2.5 - pt.pts[j*dim+i], 2.0 );
			y = sqrt( y );

			/* Plot the internal particles closest to jdx */
			if( y < 0.6 )
			{
				/* Set particle number and output file name */
				idx = j;
				sprintf( buf, "plot-%d.out", n++ );
				ff = fopen( buf, "w" );

				/* Initialize the iteration */
				for(i=0;i<dim;i++)
					size[i] = (long) dx[i] - 1;
				for(i=0;i<dim;i++)
					index[i] = 0;

				/* Do the iteration */
				do
				{
					for(i=0;i<dim;i++)
						xx[i] = x0[i] + (double) index[i] * dd[i];
					y = punity_term_deriv_evaluate_delta( &pt, idx, drv, xx, 1, NULL, 0 );
					for(i=0;i<dim;i++)
						fprintf( ff, "%15.7f", xx[i] );
					fprintf( ff, "%15.7f\n", y );
					for(i=0;i<dim;i++)
						tindex[i] = index[i];
					arraynext( (long) dim, size, tindex );
					for(i=0;i<dim;i++)
					{
						if( index[i] == dx[i] - 1 && tindex[i] != dx[i] - 1 )
						{
							fprintf( ff, "\n" );
							break;
						}
					}
				}
				while( arraynext( (long) dim, size, index ) != -1 );
				fclose( ff );
			}
		}
	}

	/* Don't forget to clean up */
	free( pts );
	free( dlt );
	free( dx );
	free( x0 );
	free( wx );
	free( dd );
	free( xx );
	free( size );
	free( index );

	return 0;
}

