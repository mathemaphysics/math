#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include "symtab.h"

int symbol_init( symbol_t *obj_in, char *name_in, int len_in )
{
    strncpy( obj_in->name, name_in, SYMBOL_MAX_NAME_LEN );
    obj_in->size = len_in;
    obj_in->data = malloc( len_in + 1 ); /* Allocate enough to null-terminate */
    bzero( obj_in->data, len_in + 1 ); /* This is important */

    return 0;
}

int symbol_free( symbol_t *obj_in )
{
    free( obj_in->data );
    obj_in->size = 0;

    return 0;
}

unsigned int symbol_hash( char *name_in )
{
    unsigned int hash = 0;
    unsigned c;

    while( c = *name_in++ )
        hash = hash * 9 ^ c;

    return hash;
}

int symbol_set_data( symbol_t *obj_in, void *data_in )
{
    bcopy( data_in, obj_in->data, obj_in->size );

    return 0;
}

int symtab_init( symtab_t *obj_in, long size_in )
{
    int i;
    obj_in->size = size_in;
    obj_in->nsyms = 0;
    obj_in->syms = (symbol_t*) malloc( size_in * sizeof(symbol_t) );
    obj_in->occs = (char*) malloc( size_in * sizeof(char) );
    for(i=0;i<size_in;i++)
        obj_in->occs[i] = 0;

    return 0;
}

int symtab_free( symtab_t *obj_in )
{
    int i;
    for(i=0;i<obj_in->nsyms;i++)
        symbol_free( obj_in->syms + i );
    free( obj_in->syms );
    free( obj_in->occs );
    obj_in->nsyms = 0;
    obj_in->size = 0;

    return 0;
}

int symtab_add( symtab_t *obj_in, char *name_in, long len_in, void *data_in )
{
    int i,offset,ulim;
    symbol_t *ptr;

    if( obj_in->nsyms >= obj_in->size )
        return -1; /* No space left in symbol table */

    offset = symbol_hash( name_in ) % (unsigned) obj_in->size;
    ulim = obj_in->size - offset;
    ptr = obj_in->syms + offset;
    for(i=0;i<ulim;i++)
    {
        if( obj_in->occs[offset+i] == 0 )
        {
            /* INSERT THE SYMBOL HERE */
            obj_in->occs[offset+i] = 1;
            ++obj_in->nsyms;
            symbol_init( obj_in->syms + offset + i, name_in, len_in );
            memcpy( ptr[i].data, data_in, len_in );
            return 0; /* Success */
        }
    }
    return -2; /* Didn't find a spot or something stupid happened */
}

int symtab_del( symtab_t *obj_in, char *var_in )
{

}

symbol_t *symtab_lookup( symtab_t *obj_in, char *var_in )
{
    int i,offset,ulim;
    symbol_t *sym;

    offset = symbol_hash( var_in ) % (unsigned) obj_in->size;
    ulim = obj_in->size - offset;
    sym = obj_in->syms + offset;

    for(i=0;i<ulim;i++)
        if( strncmp( sym[i].name, var_in, SYMBOL_MAX_NAME_LEN ) == 0 )
            return sym + i;
    return NULL;
}

void symtab_print( symtab_t *obj_in )
{
    long i,n,nx=0,m,mx=0,p=4;
    char fmt[64];
    for(i=0;i<obj_in->size;i++)
    {
        if( obj_in->occs[i] == 1 )
        {
            n = strlen( obj_in->syms[i].name );
            m = strlen( obj_in->syms[i].data );
            n = n / p + 1;
            m = m / p + 1;
            if( n > nx )
                nx = n;
            if( m > mx )
                mx = m;
        }
    }
    sprintf( fmt, " %%-%lds   %%-%lds\n", nx * p, mx * p );
    fprintf( stderr, fmt, "Variable", "Value" );
    fprintf( stderr, fmt, "========", "=====" );
    sprintf( fmt, " %%-%lds   %%-%lds\n", nx * p, mx * p );
    for(i=0;i<obj_in->size;i++)
        if( obj_in->occs[i] == 1 )
            fprintf( stderr, fmt, obj_in->syms[i].name, obj_in->syms[i].data );
}

// vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab:smarttab
