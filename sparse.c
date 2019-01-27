#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parse.h"
#include "config.h"

int sparse_save( char *fn, long n, long k, long *ia, long *ja, double *a, int mode )
{
	long i,j;
	FILE *fp = fopen( fn, "w" );
	if( fp == NULL )
		return -1;

	switch( mode )
	{
		case 0:
			fprintf( fp, "%d %d\n", n, ia[n] );
			break;
		case 1:
			fprintf( fp, "# name: A\n# type: sparse matrix\n# nnz: %d\n# rows: %d\n# columns: %d\n", ia[n], n, k );
			break;
	}
	for(i=0;i<n;i++)
	{
		for(j=ia[i];j<ia[i+1];j++)
		{
			switch( mode )
			{
				case 0:
					fprintf( fp, "%d %d %17.9f\n", i, ja[j], a[j] );
					break;
				case 1:
					fprintf( fp, "%d %d %17.9f\n", ja[j] + 1, i + 1, a[j] );
					break;
			}
		}
	}

	fclose( fp );
	return 0;
}

int sparse_load( char *fn, long *n, long *k, long **ia, long **ja, double **a, int mode )
{
	int m;
	long i,j,nnz,at,bt;
	double ct;
	char buf[1024],*tok[1024];
	FILE *fp = fopen( fn, "r" );
	if( fp == NULL )
		return -1;

	m = parse_read_line( fp, buf );
	if( m <= 0 )
		return -1;
	m = parse_stokenize( buf, tok, " \t" );
	if( m != 2 )
		return -2;
	*n = atol( tok[0] );
	nnz = atol( tok[1] );
	*ia = (long*) malloc( ( nnz + 1 ) * sizeof(long) );
	*ja = (long*) malloc( nnz * nnz * sizeof(long) );
	*a = (double*) malloc( nnz * nnz * sizeof(double) );

	for(i=0;i<*n;i++)
		(*ia)[i] = -1;
	for(i=0;i<nnz;i++)
	{
		m = parse_read_line( fp, buf );
		if( m <= 0 )
		{
			fclose( fp );
			return -3;
		}
		m = parse_stokenize( buf, tok, " \t" );
		if( m != 3 )
		{
			fclose( fp );
			return -4;
		}
		at = atol( tok[0] );
		bt = atol( tok[1] );
		ct = atof( tok[2] );
		if( (*ia)[at] == -1 )
			(*ia)[at] = i;
		(*ja)[i] = bt;
		(*a)[i] = ct;
	}
	(*ia)[*n] = nnz;

	fclose( fp );
	return 0;
}

double sparse_entry( long i, long j, long n, long m, long *ia, long *ja, double *a )
{
	long k;
	for(k=ia[i];k<ia[i+1];k++)
	{
		if( ja[k] == j )
			return a[k];
	}
	return 0.0;
}

/**
 * A symbolic multiplication of two matrices to construct
 * the output structure.
 * @param n Number rows of A
 * @param k Number columns of A and rows of B
 * @param m Number columns of B
 * @param ia Indexes of start of each row in A
 * @param ja Columns of each entry in A
 * @param ib Indexes of start of each row in B
 * @param jb Columns of each entry in B
 * @param ic Indexes of start of each row in output C
 * @param jc Columns of each entry in output C
 * @param list Temporary integer space needed for computation 
 */
void sparse_symbmm( long n, long k, long m, long *ia, long *ja, long *ib, long *jb, long *ic, long *jc, long *list )
{
	long i,j,jj,kk,len,start;

	/* clear the index list */
	j = (n>k)?(n):(k); j = (j>m)?(j):(m);
	for(i=0;i<j;i++)
		list[i] = -1;

	/* start the first row off at index 0 of jc and c */
	ic[0] = 0;

	/* let's get this show on the road */
	for(i=0;i<n;i++)
	{
		/* start a new list */
		start = -2;
		len = 0;

		/* take union of all row lists include in row i's row list in a */
		for(jj=ia[i];jj<ia[i+1];jj++)
		{
			j = ja[jj]; /* the actual column number pointed to by jj */
			for(kk=ib[j];kk<ib[j+1];kk++)
			{
				if( list[jb[kk]] == -1 )
				{
					list[jb[kk]] = start;
					start = jb[kk];
					++len;
				}
			}
		}

		/* add the list to ic and jc */
		ic[i+1] = ic[i] + len;
		for(j=ic[i];j<ic[i+1];j++)
		{
			jc[j] = start;
			start = list[start];
			list[jc[j]] = -1; /* no other entry will point to jc[j] after this */
		}
	}
}

void sparse_dgemm( long n, long k, long m, long *ia, long *ja, double *a, long *ib, long *jb, double *b, long *ic, long *jc, double *c, double *temp )
{
	long i,j,jj,kk;
	double entry;

	j = (n>k)?(n):(k);j = (j>m)?(j):(m);
	for(i=0;i<j;i++)
		temp[i] = 0.0;

	for(i=0;i<n;i++)
	{
		for(jj=ia[i];jj<ia[i+1];jj++)
		{
			j = ja[jj];
			entry = a[jj]; /* saving entry A( i, ja[jj] ) */
			for(kk=ib[j];kk<ib[j+1];kk++) /* jb[kk] = column in B */
				temp[jb[kk]] += entry * b[kk];
		}
		for(j=ic[i];j<ic[i+1];j++)
		{
			c[j] = temp[jc[j]];
			temp[jc[j]] = 0.0;
		}
	}
}

/* A very nice and VERY simple algorithm! */
void sparse_dgemv( long n, long m, long *ia, long *ja, double *a, long rb, double *b, long rc, double *c )
{
	int i,j;
	for(i=0;i<n;i++)
	{
		c[i*rc] = 0.0;
		for(j=ia[i];j<ia[i+1];j++)
			c[i*rc] = c[i*rc] + a[j] * b[ja[j]*rb];
	}
}

void sparse_dgemv_accum( long n, long m, long *ia, long *ja, double *a, long rb, double *b, long rc, double *c )
{
	int i,j;
	for(i=0;i<n;i++)
	{
		for(j=ia[i];j<ia[i+1];j++)
			c[i*rc] = c[i*rc] + a[j] * b[ja[j]*rb];
	}
}

void sparse_transp( char move, long n, long m, long *ia, long *ja, double *a, long *ib, long *jb, double *b )
{
	long i,j,jj;
	for(i=0;i<m+1;i++)
		ib[i] = 0;
	if( move == 1 )
		for(i=0;i<m;i++)
			b[i] = 0.0;

	/* count number of new columns in each row */
	ib[0] = 0;
	for(i=0;i<n;i++)
		for(j=ia[i];j<ia[i+1];j++)
			ib[ja[j]+1] = ib[ja[j]+1] + 1;

	/* "Integrate" entries forward to get final end positions */
	for(i=0;i<m;i++)
		ib[i+1] = ib[i] + ib[i+1];

	/* Counting row sizes in b done; now construct jb */
	for(i=0;i<n;i++)
	{
		for(j=ia[i];j<ia[i+1];j++)
		{
			jj = ja[j];
			jb[ib[jj]] = i;
			if( move == 1 )
				b[ib[jj]] = a[j];
			ib[jj] = ib[jj] + 1;
		}
	}

	for(i=m;i>0;i--)
		ib[i] = ib[i-1];
	ib[0] = 0;
}

void sparse_cgrad(  )
{
	
}

void sparse_proj_cgrad(  )
{

}

// vim: ts=4:sts=4:sw=4:et:sta
