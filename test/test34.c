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

int main( int argc, char **argv )
{
	int i,p,dim,optc,ret,np,*gd;
	shape_t dm;
	symtab_t sym;
	symbol_t *s;

	/* Vector counting */
	long *size,*index;
	double *pts,*div;

	/* Some switches */
	int b_config_loaded = 0;
	int b_dim_set = 0;
	int b_domain_set = 0;
	int b_grid_loaded = 0;

	while( ( optc = getopt( argc, argv, "C:" ) ) != -1 )
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
					fprintf( stderr, "Loaded following symbols from \"%s\":\n\n", optarg );
					symtab_print( &sym );
					fprintf( stderr, "\n" );
				}
				else
					fprintf( stderr, " * ERROR: Failed to load variables from file \"%s\": Returned %d\n\n", optarg, ret );
				break;
			default:
				fprintf( stderr, "I don't understand the jibberish you are spouting. Goodbye.\n" );
				return 0;
		}
	}

	if( b_config_loaded )
	{
		/* Print that searching for variables in configuration file */
                fprintf( stderr, " * Reading and executing configuration variables:\n\n" );

		/* Reading dimension symbol */
		s = symtab_lookup( &sym, "dimension" );
		if( s != NULL )
		{
			dim = atoi( (char*) s->data );
			if( dim > 0 )
			{
				b_dim_set = 1;
				fprintf( stderr, " ---> Dimension explicitly set to %d\n", dim );
			}
			else
				fprintf( stderr, " ---> ERROR: Dimension out of bounds\n" );
		}

		/* Reading domain symbols */
		s = symtab_lookup( &sym, "domain" );
		if( s != NULL )
		{
			load_domain( (char*) s->data, &dm, &ret );
			if( ret == 0 )
			{
				if( b_dim_set )
				{
					if( dm.dim == dim )
						fprintf( stderr, " ---> Loaded domain information: Dimension = %d\n", dim );
					else
					{
						fprintf( stderr, " ---> Loaded domain information but dimension does not match configuration: " );
						fprintf( stderr, " Setting dimension = %d\n", dm.dim );
						dim = dm.dim; /* Variable "domain" value overrides "dimension" */
					}
				}
				else
				{
					fprintf( stderr, " ---> Loaded domain information: Dimension = %d\n", dm.dim );
					dim = dm.dim;
					b_dim_set = 1;
				}
				b_domain_set = 1;
			}
			else
				fprintf( stderr, " ---> ERROR: Failed to load domain indicator\n" );
		}

		/* Read grid information */
		s = symtab_lookup( &sym, "grid" );
		if( s != NULL )
		{
			load_grid( (char*) s->data, &gd, &ret );
			if( ret > 0 )
			{
				if( b_dim_set )
				{
					if( ret < dim ) /* This is okay, just too much grid info given */
						fprintf( stderr, " ---> ERROR: Not enough grid info specified\n" );
					else
					{
						fprintf( stderr, " ---> Grid data loaded: " );
						for(i=0;i<ret;i++)
							fprintf( stderr, "%5d", gd[i] );
						fprintf( stderr, "\n" );
						b_grid_loaded = 1;
					}
				}
				else
				{
					dim = ret;
					b_dim_set = 1;
				}
			}
			else
				fprintf( stderr, " ---> ERROR: Failed to load grid information\n" );
		}
	}

	if( !b_domain_set )
	{
		fprintf( stderr, " * ERROR: No domain indicator loaded. Exiting.\n\n" );
		return 0;
	}
	if( !b_grid_loaded )
	{
		fprintf( stderr, " * ERROR: No grid information loaded. Exiting.\n\n" );
	}
	fprintf( stderr, "\n" );

	/* Now make the grid and write the points to output */
	np = 1;
        for(i=0;i<dim;i++)
                np *= gd[i];
	div = (double*) malloc( dim * sizeof(double) );
	if( dm.type == 0 )
		for(i=0;i<dim;i++)
			div[i] = dm.params[0] / (double) ( gd[i] - 1 );
	if( dm.type == 1 )
		for(i=0;i<dim;i++)
			div[i] = dm.params[i] / (double) ( gd[i] - 1 );
	pts = (double*) malloc( dim * np * sizeof(double) );
	size = (long*) malloc( dim * sizeof(long) );
	index = (long*) malloc( dim * sizeof(long) );
	for(i=0;i<dim;i++)
		index[i] = 0, size[i] = gd[i] - 1;
	p = 0;
	do
	{
		for(i=0;i<dim;i++)
			pts[p*dim+i] = dm.orig[i] + (double) index[i] * div[i];
		++p;
	}
	while( arraynext( (long) dim, size, index ) != -1 );

	/* Write points to output file */
	write_points( "cloud.grid.out", dim, np, pts, 0 );

	/* Generate partition of unity */
	

	return 0;
}

