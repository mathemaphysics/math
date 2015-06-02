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

int main( int argc, char **argv )
{
	int optc;
	int b_config_loaded = 0;
	symtab_t sym;
	symbol_t *s;

	/* Load the options */
	while( ( optc = getopt( argc, argv, "C:" ) ) != -1 )
        {
                switch( optc )
                {
			case 'C':
				/* Load symbols */
				symtab_init( &sym, 1024 );
				cfgread_load_symbols_f( optarg, &sym );
				b_config_loaded = 1;

				/* Let everyone know */
				fprintf( stderr, "Loaded following symbols from \"%s\":\n\n", optarg );
				symtab_print( &sym );
				fprintf( stderr, "\n" );
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
		s = symtab_lookup( &sym, "variable" );
		if( s != NULL )
		{
			
                }
	}
}

