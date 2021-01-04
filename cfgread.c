#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "parse.h"
#include "symtab.h"

#define CFGREAD_BUFFER_SIZE 8192
#define CFGREAD_TOKEN_BUFFER_SIZE 4096

/**
 * Read all symbols set in a file named fn_in and
 * enter them into the symbol table smt_in, which
 * must already be initialized.
 * @param fn_in Name of file to read
 * @param smt_in Symbol table to add variables to
 * @return Returns 0 if all is well, -1 if file not open and -2 if failure to load symbols
 */
int cfgread_load_symbols_f( char *fn_in, symtab_t *smt_in )
{
    int i,n;
    char buf[CFGREAD_BUFFER_SIZE];
    char *tok[CFGREAD_TOKEN_BUFFER_SIZE];
    char *var,*val,*tp;
    symbol_t *tmp;

    FILE *fp = fopen( fn_in, "r" );
    if( fp == NULL )
        return -1; /* Failed to open file */
    while( parse_read_line( fp, buf ) != -1 )
    {
        n = parse_stokenize( buf, tok, "=" );
        if( n != 2 ) /* Possibly an error but... */
            continue;
        else
        {
            var = parse_strip_white( tok[0] );
            val = parse_strip_white( tok[1] );
            tmp = symtab_lookup( smt_in, var );
            if( tmp == NULL )
            {
                /* Symbol does not exist in the table so add it */
                if( symtab_add( smt_in, var, (long) strlen( val ), val ) < 0 )
                {
                    fclose( fp );
                    return -2;
                }
            }
            else
            {
                /* Symbol already exists so change its value */
                if( strlen( val ) > tmp->size )
                {
                    tp = (char*) realloc( tmp->data, strlen( val ) * sizeof(char) );
                    if( tp == NULL )
                        return -3;
                    tmp->data = tp;
                    tmp->size = strlen( val );
                }
                strncpy( tmp->data, val, tmp->size );
            }
        }
    }
    return 0;
}

/**
 * Does the same as cfgread_load_symbols_f but does so
 * by reading directly from a file point instead of opening
 * the file itself.
 * @param fp_in File pointer from which to read
 * @param smt_in The symbol table into which to deposit the variables found
 * @return Returns 0 if all is well and -1 if fail
 */
int cfgread_load_symbols_s( FILE *fp_in, symtab_t *smt_in )
{
    int i,n;
    char buf[CFGREAD_BUFFER_SIZE];
    char *tok[CFGREAD_TOKEN_BUFFER_SIZE];
    char *var,*val,*tp;
    symbol_t *tmp;

    if( fp_in == NULL )
        return -1; /* Failed to open file */
    while( parse_read_line( fp_in, buf ) != -1 )
    {
        n = parse_stokenize( buf, tok, "=" );
        if( n != 2 ) /* Possibly an error but... */
            continue;
        else
        {
            var = parse_strip_white( tok[0] );
            val = parse_strip_white( tok[1] );
            tmp = symtab_lookup( smt_in, var );
            if( tmp == NULL )
            {
                /* Symbol does not exist in the table so add it */
                if( symtab_add( smt_in, var, (long) strlen( val ), val ) < 0 )
                    return -2;
            }
            else
            {
                /* Symbol already exists so change its value */
                if( strlen( val ) > tmp->size )
                {
                    tp = (char*) realloc( tmp->data, strlen( val ) * sizeof(char) );
                    if( tp == NULL )
                        return -3;
                    tmp->data = tp;
                    tmp->size = strlen( val );
                }
                strncpy( tmp->data, val, tmp->size );
            }
        }
    }
    return 0;
}

/**
 * Saves symbols from a symbol table into a text file
 * in text format using "=" and smart spacing of variable
 * and value.
 * @param fn_in The name of the file to write to
 * @param smt_in The symbol table from which to write the symbols
 * @return Returns 0 if okay and -1 if fail
 */
int cfgread_save_symbols_f( char *fn_in, symtab_t *smt_in )
{
    int i,n,nx=0,m,mx=0,p=4;
    char fmt[64];

    FILE *fp = fopen( fn_in, "w" );
    if( fp == NULL )
        return -1;
    for(i=0;i<smt_in->size;i++)
    {
        if( smt_in->occs[i] == 1 )
        {
            n = strlen( smt_in->syms[i].name );
            m = strlen( smt_in->syms[i].data );
            n = n / p + 1;
            m = m / p + 1;
            if( n > nx )
                nx = n;
            if( m > mx )
                mx = m;
        }
    }
    sprintf( fmt, "%%-%ds = %%-%ds\n", nx * p, mx * p );
    for(i=0;i<smt_in->size;i++)
        if( smt_in->occs[i] == 1 )
            fprintf( fp, fmt, smt_in->syms[i].name, smt_in->syms[i].data );
    fclose( fp );
    return 0;
}

/**
 * Does same as above, but taken a stream pointer.
 * @param fp_in File pointer to write to
 * @param smt_in Symbol table containing variables
 * @return Returns 0 if okay and -1 if fail
 */
int cfgread_save_symbols_s( FILE *fp_in, symtab_t *smt_in )
{
    int i,n,nx=0,m,mx=0,p=4;
    char fmt[64];

    if( fp_in == NULL )
        return -1;
    for(i=0;i<smt_in->size;i++)
    {
        if( smt_in->occs[i] == 1 )
        {
            n = strlen( smt_in->syms[i].name );
            m = strlen( smt_in->syms[i].data );
            n = n / p + 1;
            m = m / p + 1;
            if( n > nx )
                nx = n;
            if( m > mx )
                mx = m;
        }
    }
    sprintf( fmt, "%%-%ds = %%-%ds\n", nx * p, mx * p );
    for(i=0;i<smt_in->size;i++)
        if( smt_in->occs[i] == 1 )
            fprintf( fp_in, fmt, smt_in->syms[i].name, smt_in->syms[i].data );
    return 0;	
}

// vim: ts=4:sts=4:sw=4:et:sta
