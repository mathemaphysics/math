#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "parse.h"
#include "utildefs.h"

#define MAX_WHITE_SPACE 1024
#define SNODE_CHILD_ALLOC_INCREMENT 128

/**
 * Find the position of the first non-whitespace character.
 * @param str_in The string of interest
 * @return Position of the first non-whitespace character
 */
int parse_get_start( char *str_in )
{
    int i;
    for(i=0;i<MAX_WHITE_SPACE;i++)
        if( str_in[i] != ' ' && str_in[i] != '\t' )
            break;
    return i;
}

/**
 * Use this function to see if chr_in is anywhere in the string
 * given in sep_in. This is needed to tokenizing strings.
 * @param chr_in Character to check
 * @param sep_in String of separators; see if chr_in is in here
 * @return Returns 1 if chr_in is in sep_in, 0 otherwise
 */
int parse_contains_char( char chr_in, char *sep_in )
{
    int i = 0;
    while( sep_in[i] != '\0' )
        if( chr_in == sep_in[i++] )
            return 1;
    return 0;
}

void parse_print_white( FILE *fp_in, int n_in )
{
	int i;
	for(i=0;i<n_in;i++)
		fprintf( fp_in, " " );
}

void parse_print_back( FILE *fp_in, int n_in )
{
	int i;
	for(i=0;i<n_in;i++)
		fprintf( fp_in, "\b" );
}

/**
 * Use this to reposition the start of the input string
 * str_in to cut off initial whitespace and trim it off
 * the end by inserting null characters.
 * @param str_in The input string
 * @return Returns a pointer to the whitespace-erased beginning of str_in
 */
char *parse_strip_white( char *str_in )
{
    int n,m;
    char *ptr;
    n = parse_get_start( str_in );
    m = strlen( str_in ) - 1;
    while( m >= 0 && parse_contains_char( str_in[m], " \t" ) == 1 )
        str_in[m--] = '\0';
    ptr = &str_in[n];
    return ptr;
}

/**
 * Read a line out of the file pointer into buf_out and return
 * the number of characters read or -1 if EOF.
 * @param fp_in The file pointer, already open
 * @param buf_out The buffer to use to store the line
 * @return Returns the number of characters read from the stream
 */
int parse_read_line( FILE *fp_in, char *buf_out )
{
    int i;
    char ch;

    if( !feof( fp_in ) )
    {
        buf_out[0] = '\0';
        i = 0;
        while( ( ch = fgetc( fp_in ) ) != -1 )
        {
            if( ch != '\n' )
                buf_out[i++] = (char) ch;
            else
            {
                buf_out[i] = '\0';
                break;
            }
        }
        return i;
    }
    else
        return -1;
}

/**
 * This function does the same thing as parse_ftokenize except
 * it does it by taking its data from a buffer, buf_in instead
 * of from a file stream pointer.
 * @param buf_in The originating buffer
 * @param tok_out Tokens given as output
 * @param sep_in The tokenizing characters in a null-terminated string
 * @return Returns the number of tokens found
 */
int parse_stokenize( char *buf_in, char **tok_out, char *sep_in )
{
    int i,j,k,dx;
    char cmt[] = "@#"; /* Lines beginning with these characters are ignored */

    i = strlen( buf_in );
    k = parse_get_start( buf_in );
    if( parse_contains_char( buf_in[k], cmt ) == 0 )
    {
        for(j=0;j<i-k;j++)
            if( parse_contains_char( buf_in[k+j], sep_in ) == 1 )
                buf_in[k+j] = '\0';
        dx = 0;
        if( buf_in[k] != '\0' )
            tok_out[dx++] = &buf_in[k];
        for(j=1;j<i-k-1;j++)
            if( buf_in[k+j] == '\0' && buf_in[k+j+1] != '\0' )
                tok_out[dx++] = &buf_in[k+j+1];
        return dx;
    }
    else
        return 0;
}

/**
 * Load a single line from the file into buf_out.
 * Function returns the integer number of tokens
 * found in the line.
 * @param fp_in The file pointer; should be open already
 * @param buf_out The pre-allocated buffer to use for tokenizing
 * @param tok_out The pointer to a series of char pointers pointing to each token inside buf_out
 * @return Returns the number of tokens found
 */
int parse_ftokenize( FILE *fp_in, char *buf_out, char **tok_out, char *sep_in )
{
    int i,j,k,ch,dx;
    char cmt[] = "@#"; /* Lines beginning with these characters are ignored */

    if( !feof( fp_in ) )
    {
        i = 0;
        buf_out[i] = '\0';
        while( ( ch = fgetc( fp_in ) ) != -1 )
        {
            if( ch != '\n' )
                buf_out[i++] = (char) ch;
            else
            {
                buf_out[i] = '\0';
                break;
            }
        }
	k = parse_get_start( buf_out );
	if( parse_contains_char( buf_out[k], cmt ) == 0 )
	{
            for(j=0;j<i-k;j++)
                if( parse_contains_char( buf_out[k+j], sep_in ) == 1 )
                    buf_out[k+j] = '\0';
            dx = 0;
            if( buf_out[k] != '\0' )
                tok_out[dx++] = &buf_out[k];
            for(j=1;j<i-k-1;j++)
                if( buf_out[k+j] == '\0' && buf_out[k+j+1] != '\0' )
                    tok_out[dx++] = &buf_out[k+j+1];
            return dx;
        }
	else
            return 0;
    }
    else
	return -1;
}

/**
 * Read a comma or space separated vector into memory
 */
int parse_read_vector( char *str_in, double **vec_out )
{
	int i;
}

/**
 * The 4-byte integer version of parse_read_ranges_long
 */
int parse_read_ranges( char *str_in, int ***range_out )
{
    int i,m;
    char buf[UTIL_BUFFER_LEN], bug[UTIL_BUFFER_LEN];
    char *tok[UTIL_BUFFER_LEN], *tol[UTIL_BUFFER_LEN];
    int n,cnfr,**range;

    /* Allocate range space */
    *range_out = (int**) malloc( UTIL_MAX_RANGES * sizeof(int*) );
    for(i=0;i<UTIL_MAX_RANGES;i++)
        (*range_out)[i] = (int*) malloc( 2 * sizeof(int) );

    /* Alias *range_out to range */
    range = *range_out;

    /* Copy input string into first buffer */
    strcpy( buf, str_in ); /* Cannot change str_in */

    /* Start counting */
    cnfr = 0;
    n = parse_stokenize( buf, tok, "," );
    for(i=0;i<n;i++)
    {
        strcpy( bug, tok[i] );
        m = parse_stokenize( bug, tol, "-" );
        if( m > 2 ) /* Only want definitive ranges like 127-439 */
            continue;
        if( m == 1 ) /* Then print single page */
        {
            range[cnfr][0] = atol( tol[0] );
            range[cnfr][1] = range[cnfr][0]; /* Same indicating a single frame */
            ++cnfr;
        }
        else
        {
            range[cnfr][0] = atol( tol[0] );
            range[cnfr][1] = atol( tol[1] );
            if( range[cnfr][1] < range[cnfr][0] )
                fprintf( stderr, "Warning: Bad central range given; ignoring.\n" );
            else
                ++cnfr;
        }
        if( cnfr > UTIL_MAX_RANGES )
        {
            fprintf( stderr, "ERROR: Too many ranges specified on command line. Sorry. Exiting." );
            exit(0);
        }
    }
    fprintf( stderr, "Central group: " );
    for(i=0;i<cnfr;i++)
    {
        if( range[i][1] == range[i][0] )
            fprintf( stderr, "%d ", range[i][0] );
        else
            fprintf( stderr, "%d-%d ", range[i][0], range[i][1] );
    }

    /* Return the number of ranges and values */
    return cnfr;
}
/**
 * Read integer range list from string
 */
long parse_read_ranges_long( char *str_in, long ***range_out )
{
    int i,m;
    char buf[UTIL_BUFFER_LEN], bug[UTIL_BUFFER_LEN];
    char *tok[UTIL_BUFFER_LEN], *tol[UTIL_BUFFER_LEN];
    long n,cnfr;
    long **range;

    /* Allocate range space */
    *range_out = (long**) malloc( UTIL_MAX_RANGES * sizeof(long*) );
    for(i=0;i<UTIL_MAX_RANGES;i++)
        (*range_out)[i] = (long*) malloc( 2 * sizeof(long) );

    /* Alias *range_out to range */
    range = *range_out;

    /* Copy input string into first buffer */
    strcpy( buf, str_in ); /* Cannot change str_in */

    /* Start counting */
    cnfr = 0;
    n = parse_stokenize( buf, tok, "," );
    for(i=0;i<n;i++)
    {
        strcpy( bug, tok[i] );
        m = parse_stokenize( bug, tol, "-" );
        if( m > 2 ) /* Only want definitive ranges like 127-439 */
            continue;
        if( m == 1 ) /* Then print single page */
        {
            range[cnfr][0] = atol( tol[0] );
            range[cnfr][1] = range[cnfr][0]; /* Same indicating a single frame */
            ++cnfr;
        }
        else
        {
            range[cnfr][0] = atol( tol[0] );
            range[cnfr][1] = atol( tol[1] );
            if( range[cnfr][1] < range[cnfr][0] )
                fprintf( stderr, "Warning: Bad central range given; ignoring.\n" );
            else
                ++cnfr;
        }
        if( cnfr > UTIL_MAX_RANGES )
        {
            fprintf( stderr, "ERROR: Too many ranges specified on command line. Sorry. Exiting." );
            exit(0);
        }
    }
    fprintf( stderr, "Central group: " );
    for(i=0;i<cnfr;i++)
    {
        if( range[i][1] == range[i][0] )
            fprintf( stderr, "%ld ", range[i][0] );
        else
            fprintf( stderr, "%ld-%ld ", range[i][0], range[i][1] );
    }

    /* Return the number of ranges and values */
    return cnfr;
}

/**
 * This function prints out a formatted line describing
 * the function of interest; it takes the flag name, the
 * value of that flag, if it has one, and the description
 * of the flag variable
 *
 * @param flag The string representing the flag, i.e. "-f"
 * @param val  The value of the flag's variable
 * @param desc Description of the variable's meaning
 */
void print_usage_line( char *flag, char *val, char *desc )
{
    int fw_flag = 5;
    int fw_val  = 12;
    int fw_desc = 50;
    int cnt;
    char tmp[128] = "";
    char ump[128] = "";
    snprintf( ump, (size_t) 127, "%+*s", fw_val, val );
    for(cnt=0;cnt<strlen(desc);cnt+=fw_desc)
    {
        snprintf( tmp, (size_t) 127, "%-*s", fw_desc, desc + cnt );
        printf( " %*s %+*s %*s\n", fw_flag, (cnt==0 ? flag : ""), fw_val, (cnt==0 ? ump : ""), fw_desc, tmp );
    }
}

//typedef struct snode
//{
//	int id;
//	int len;
//	char *str;
//	int deg;
//	int alloc;
//	struct snode *child;
//	char op[SNODE_OP_SIZE];
//	char optype;
//} snode_t;

int snode_init( snode_t *obj_in, int len_in )
{
	static int idx = 0;
	obj_in->id = idx++;
	if( len_in <= 0 )
	    return -1;
	obj_in->len = len_in;
	obj_in->str = (char*) malloc( len_in * sizeof(char) );
	if( obj_in->str == NULL )
	     return -2;
	obj_in->child = NULL;
	obj_in->deg = 0;
	obj_in->alloc = 0;
	return 0;
}

int snode_alloc( snode_t *obj_in, int len_in )
{
	char *tmp;
	if( len_in < 0 )
	    return -1;
	if( len_in == 0 )
	{
	    if( obj_in->str == NULL )
		return -2;
	    free( obj_in->str );
	    obj_in->str = NULL;
	    obj_in->len = 0;
            return 0;
	}
	tmp = (char*) realloc( obj_in->str, len_in * sizeof(char) );
	if( tmp == NULL )
	    return -3;
	obj_in->str = tmp;
	obj_in->len = len_in;

	return 0;
}

int snode_free( snode_t *obj_in )
{
	obj_in->deg = 0;
	obj_in->len = 0;
	return snode_alloc( obj_in, 0 );
}

int snode_set_string( snode_t *obj_in, char *str_in )
{
	if( strncpy( obj_in->str, str_in, obj_in->len ) == NULL )
		return -1;
	return 0;
}

int snode_add_child( snode_t *obj_in, char *str_in )
{
	snode_t *tmp;
	if( obj_in->deg + 1 > obj_in->alloc )
	{
		tmp = (snode_t*) realloc( obj_in->child, obj_in->alloc + SNODE_CHILD_ALLOC_INCREMENT * sizeof(snode_t) );
		if( tmp == NULL )
			return -1;
		obj_in->child = tmp;
		obj_in->alloc += SNODE_CHILD_ALLOC_INCREMENT;
	}
	snode_init( obj_in->child + obj_in->deg, strlen( str_in ) + 1 );
	snode_set_string( obj_in->child + obj_in->deg, str_in );
	obj_in->deg++;
	return 0;
}

/**
 * This function sets the operation to be performed on all
 * child nodes based on the parsing done in the tree above
 * this point. This is for use with syntax trees, which is
 * what this data type was intended to deal with.
 * @param obj_in The node object
 * @param op_in The string denoting the operation to perform on child nodes
 * @return Returns 1 and sets cop = NULL if op_in == NULL, returns -1 if fail to copy and 0 if all okay
 */
int snode_set_op( snode_t *obj_in, char *op_in, char optype_in )
{
	if( op_in == NULL )
	{
		op_in = NULL;
		return 1;
	}
	if( strncpy( obj_in->op, op_in, SNODE_OP_SIZE ) == NULL )
		return -1;
	obj_in->optype = optype_in;
	return 0;
}

/**
 * This function frees all nodes in the tree below and
 * including the obj_in node. This frees all children of
 * children, etc. Note that this function is recursive,
 * so watch for stack overflow.
 * @param obj_in The node object
 * @return Returns 0 if all is well or -1 if something fails
 */
int stree_free( snode_t *obj_in )
{
	int i;
	for(i=0;i<obj_in->deg;i++)
		if( stree_free( obj_in->child + i ) < 0 )
			return -1;
	return snode_free( obj_in );
}

int stree_print( snode_t *obj_in, int level_in )
{
	int i;
	for(i=0;i<level_in;i++)
		printf( "\t" );
	if( obj_in->str != NULL )
		printf( "%d: %s, op = %s, optype = %d\n", obj_in->id, obj_in->str, obj_in->op, obj_in->optype );
	else
		return -1;
	for(i=0;i<obj_in->deg;i++)
		stree_print( obj_in->child + i, level_in + 1 );
	return 0;
}

