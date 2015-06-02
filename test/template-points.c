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
	int optc,ret;
	symtab_t sym;
	symbol_t *s;

	/* Points data structures */
	int dim,npts;
	double *pts;

	/* Switches to indicate boolean things */
	int b_config_loaded = 0;
	int b_points_loaded = 0;

	/* Announce */
	fprintf( stderr, "\nCloud Template VERSION 0.1\n\n" );

	/* Load the options */
	while( ( optc = getopt( argc, argv, "C:p:" ) ) != -1 )
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
				}
				else
					fprintf( stderr, " * ERROR: Failed loading points from file \"%s\"\n\n", optarg );
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
				}
				else
					fprintf( stderr, " ---> ERROR: Failed loading points from file \"%s\"\n", (char*) s->data );
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

	/* Do your business here */
	

	/* Don't forget to clean up */
	free( pts );

	return 0;
}

