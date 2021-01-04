#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int polynomial_exponents( int, int, int, int * );
long int polynomial_exponents_long( long int, long int, long int, long int * );
int truncfact( int, int );
long int truncfact_long( long int, long int );

/**
 * This is the standard factorial function. Improvement can certainly be made
 * by rewriting it more intelligently.
 * @param num_in The number whose factorial you want
 * @return The factorial of num_in
 */
int factorial( int num_in )
{
	int fct = 1;
	int i;
	if( num_in <= 0 )
		return 1;
	for(i=num_in;i>1;i--)
		fct *= i;
	return fct;
}

long int factorial_long( long int num_in )
{
	long int fct = 1;
	long int i;
	if( num_in <= 0 )
		return 1;
	for(i=num_in;i>1;i--)
		fct *= i;
	return fct;
}

int binomial_slow( int n_in, int k_in )
{
	if( n_in < k_in )
		return 0;
	else
		return factorial( n_in ) / factorial( k_in ) / factorial( n_in - k_in );
}

long int binomial_slow_long( long int n_in, long int k_in )
{
	if( n_in < k_in )
		return 0;
	else
		return factorial( n_in ) / factorial( k_in ) / factorial( n_in - k_in );
}

int binomial( int n_in, int k_in )
{
	if( n_in < k_in )
		return 0;
	else
		return truncfact( n_in, k_in ) / factorial( k_in );
}

long int binomial_long( long int n_in, long int k_in )
{
	if( n_in < k_in )
		return 0;
	else
		return truncfact_long( n_in, k_in ) / factorial_long( k_in );
}

/**
 * This function returns the number of r_in-permutations of n_in objects.
 */
int truncfact( int n_in, int r_in )
{
	int i,prd = 1;
	for(i=0;i<r_in;i++)
		prd *= n_in - i;
	return prd;
}

long int truncfact_long( long int n_in, long int r_in )
{
	long int i,prd = 1;
	for(i=0;i<r_in;i++)
		prd *= n_in - i;
	return prd;
}

/**
 * Given a combinadic index, its corresponding vector is
 * the strictly increasing non-negative integer sequence,
 * {m_k}, of length dim_in and whose entries are less than
 * or equal to lim_in, and such that
 *   C(m_k,k+1) + C(m_{k-1},k) + ... + C(m_0,1) = idx_in
 */
int combinadic_vector( int idx_in, int lim_in, int dim_in, int *vec_out )
{
	int i,j,idx,tmp;

	if( idx_in < 0 )
		return -1;

	idx = idx_in;
	for(i=dim_in-1;i>=0;i--)
	{
		for(j=i+1;j<=lim_in;j++)
		{
			tmp = binomial( j, i + 1 );
			if( tmp > idx )
				break;
		}
		idx -= binomial( j - 1, i + 1 );
		vec_out[i] = j - 1;
	}
	return 0;
}

long int combinadic_vector_long( long int idx_in, long int lim_in, long int dim_in, long int *vec_out )
{
	long int i,j,idx,tmp;

	if( idx_in < 0 )
		return -1;

	idx = idx_in;
	for(i=dim_in-1;i>=0;i--)
	{
		for(j=i+1;j<=lim_in;j++)
		{
			tmp = binomial( j, i + 1 );
			if( tmp > idx )
				break;
		}
		idx -= binomial( j - 1, i + 1 );
		vec_out[i] = j - 1;
	}
	return 0;
}

int combinadic_vector_length( int idx_in, int dim_in, int lim_in )
{
	int i,sum;

	sum = 0;
	for(i=1;i<=dim_in;i++)
	{
		sum += binomial( lim_in, i );
		if( sum > idx_in )
			break;
	}
	return i;
}

long int combinadic_vector_length_long( long int idx_in, long int dim_in, long int lim_in )
{
	long int i,sum;

	sum = 0;
	for(i=1;i<=dim_in;i++)
	{
		sum += binomial( lim_in, i );
		if( sum > idx_in )
			break;
	}
	return i;
}

int combinadic_next( int lim_in, int dim_in, int *vec_in )
{
	int i,j;
	for(i=0;i<dim_in;i++)
	{
		if( ( vec_in[i] < lim_in )
			&& ( ( i < dim_in - 1 && ( vec_in[i+1] - vec_in[i] ) > 1 ) || i == dim_in - 1 ) )
		{
			++vec_in[i];
			break;
		}
	}
	for(j=0;j<i;j++)
		vec_in[j] = j;
}

long int combinadic_next_long( long int lim_in, long int dim_in, long int *vec_in )
{
	long int i,j;
	for(i=0;i<dim_in;i++)
	{
		if( ( vec_in[i] < lim_in )
			&& ( ( i < dim_in - 1 && ( vec_in[i+1] - vec_in[i] ) > 1 ) || i == dim_in - 1 ) )
		{
			++vec_in[i];
			break;
		}
	}
	for(j=0;j<i;j++)
		vec_in[j] = j;
}

int combinadic_init( int lim_in, int dim_in, int *ptr_in )
{
	int i;
	if(lim_in<dim_in)
		return -1; // can't make dim_in-combinations of lim_in (< dim_in) objects
	for(i=0;i<dim_in;i++)
		ptr_in[i] = i;
	return 0;
}

long int combinadic_init_long( long int lim_in, long int dim_in, long int *ptr_in )
{
	long int i;
	if(lim_in<dim_in)
		return -1; // can't make dim_in-combinations of lim_in (< dim_in) objects
	for(i=0;i<dim_in;i++)
		ptr_in[i] = i;
	return 0;
}

int combinadic_index( int dim_in, int *vec_in )
{
	int i;
	int idx = 0;
	for(i=0;i<dim_in;i++)
		idx += binomial( vec_in[i], i + 1 );
	return idx;
}

long int combinadic_index_long( long int dim_in, long int *vec_in )
{
	long int i;
	long int idx = 0;
	for(i=0;i<dim_in;i++)
		idx += binomial( vec_in[i], i + 1 );
	return idx;
}

int rcombinadic_init( int lim_in, int dim_in, int *ptr_in )
{
	return combinadic_init( lim_in + dim_in - 1, dim_in - 1, ptr_in );
}

long int rcombinadic_init_long( long int lim_in, long int dim_in, long int *ptr_in )
{
	return combinadic_init_long( lim_in + dim_in - 1, dim_in - 1, ptr_in );
}

int rcombinadic_index( int dim_in, int *vec_in )
{
	return combinadic_index( dim_in - 1, vec_in );
}

long int rcombinadic_index_long( long int dim_in, long int *vec_in )
{
	return combinadic_index_long( dim_in - 1, vec_in );
}

int rcombinadic_next( int lim_in, int dim_in, int *vec_in )
{
	return combinadic_next( lim_in + dim_in - 1, dim_in - 1, vec_in );
}

long int rcombinadic_next_long( long int lim_in, long int dim_in, long int *vec_in )
{
	return combinadic_next_long( lim_in + dim_in - 1, dim_in - 1, vec_in );
}

int rcombinadic_vector( int idx_in, int lim_in, int dim_in, int *vec_out )
{
	return combinadic_vector( idx_in, lim_in + dim_in - 1, dim_in - 1, vec_out );
}

long int rcombinadic_vector_long( long int idx_in, long int lim_in, long int dim_in, long int *vec_out )
{
	return combinadic_vector_long( idx_in, lim_in + dim_in - 1, dim_in - 1, vec_out );
}
/* takes a combinadic vector of length dim_in and outputs
   an occupancy vector of length dim_in + 1 */
int rcombinadic_occupancy( int lim_in, int dim_in, int *vec_in, int *occ_out )
{
	int i;
	if( dim_in < 1 )
		return -1; /* Failure: Parameter out of acceptable domain */
	if( dim_in > 1 )
	{
		occ_out[0] = vec_in[0];
		for(i=1;i<dim_in-1;i++)
			occ_out[i] = vec_in[i] - vec_in[i-1] - 1;
		occ_out[dim_in-1] = lim_in + dim_in - 1 - vec_in[dim_in-2] - 1; /* minus addtl 'one' adjust difference of 1 to mean zero occupancy */
	}
	else
		occ_out[0] = lim_in;
	return 0;
}

long int rcombinadic_occupancy_long( long int lim_in, long int dim_in, long int *vec_in, long int *occ_out )
{
	long int i;
	if( dim_in < 1 )
		return -1; /* Failure: Parameter out of acceptable domain */
	if( dim_in > 1 )
	{
		occ_out[0] = vec_in[0];
		for(i=1;i<dim_in;i++)
			occ_out[i] = vec_in[i] - vec_in[i-1] - 1;
		occ_out[dim_in-1] = lim_in + dim_in - 1 - vec_in[dim_in-2] - 1; /* minus addtl 'one' adjust difference of 1 to mean zero occupancy */
	}
	else
		occ_out[0] = lim_in;
	return 0;
}

/* takes an occupancy vector of length dim_in and outputs
   a vector of dim dim_in - 1 */
int rcombinadic_invoccup( int dim_in, int *occ_in, int *vec_out )
{
	int i;
	vec_out[0] = occ_in[0];
	for(i=1;i<dim_in;i++)
		vec_out[i] = vec_out[i-1] + occ_in[i] + 1;
	return 0;
}

long int rcombinadic_invoccup_long( long int dim_in, long int *occ_in, long int *vec_out )
{
	long int i;
	vec_out[0] = occ_in[0];
	for(i=1;i<dim_in;i++)
		vec_out[i] = vec_out[i-1] + occ_in[i] + 1;
	return 0;
}

int polynomial_index( int dim_in, int *exp_in )
{
	int cmb[dim_in-1];
	rcombinadic_invoccup( dim_in, exp_in, cmb );
	return rcombinadic_index( dim_in, cmb );
}

long int polynomial_index_long( long int dim_in, long int *exp_in )
{
	long int cmb[dim_in-1];
	rcombinadic_invoccup_long( dim_in, exp_in, cmb );
	return rcombinadic_index_long( dim_in, cmb );
}

int global_poly_index( int ord_in, int dim_in, int *exp_in )
{
	return polynomial_index( dim_in, exp_in ) + binomial( ord_in + dim_in - 1, dim_in );
}

long int global_poly_index_long( long int ord_in, long int dim_in, long int *exp_in )
{
	return polynomial_index_long( dim_in, exp_in ) + binomial_long( ord_in + dim_in - 1, dim_in );
}

/**
 * Generate the global polynomial vector idx_in of max order ord_in
 * and of dimension dim_in and put the exponent in exp_out
 */
int global_polynomial_vector( int idx_in, int dim_in, int *exp_out )
{
	int n,m,r,cmb[dim_in-1];

	for(n=0,r=0;;r++)
	{
		m = binomial( r + dim_in - 1, dim_in - 1 );
		if( n + m > idx_in )
			break;
		else
			n += m;
	}

	polynomial_exponents( idx_in - n, r, dim_in, exp_out );

	return 0;
}

long int global_polynomial_vector_long( long int idx_in, long int dim_in, long int *exp_out )
{
	long int n,m,r,cmb[dim_in-1];

        for(n=0,r=0;;r++)
        {
                m = binomial( r + dim_in - 1, dim_in - 1 );
                if( n + m > idx_in )
                        break;
                else
                        n += m;
        }

        polynomial_exponents_long( idx_in - n, r, dim_in, exp_out );

        return 0;
}

/* NOTE: it's probably best to stick to a two-index system: (term_order,local_index) */
/* takes an index of a specified order only */
int polynomial_exponents( int idx_in, int ord_in, int dim_in, int *occ_out )
{
	int vec[dim_in-1];
	rcombinadic_vector( idx_in, ord_in, dim_in, vec );
	rcombinadic_occupancy( ord_in, dim_in, vec, occ_out );
	return 0;
}

long int polynomial_exponents_long( long int idx_in, long int ord_in, long int dim_in, long int *occ_out )
{
	long int vec[dim_in-1];
	rcombinadic_vector_long( idx_in, ord_in, dim_in, vec );
	rcombinadic_occupancy_long( ord_in, dim_in, vec, occ_out );
	return 0;
}

/*
 * FACTORADIC OPERATIONS
 */

int factoradic_radix_index( int dim_in, int idx_in, int *rad_out )
{
	int i,fct,div,rem;

	rem = idx_in;
	for(i=dim_in-1;i>=0;i--)
	{
		fct = factorial( i );
		div = rem / fct;
		rem = rem % fct;
		rad_out[dim_in-1-i] = div;
	}
}

int factoradic_radix_vector( int dim_in, int *vec_in, int *rad_out )
{
	int c,i,j;
	for(i=0;i<dim_in;i++)
		rad_out[i] = 0;
	for(i=0;i<dim_in;i++)
	{
		for(j=i+1;j<dim_in;j++)
		{
			if( vec_in[j] < vec_in[i] )
				++rad_out[i];
		}
	}
	return 0l;
}

int factoradic_index( int dim_in, int *vec_in )
{
	int i,rad[dim_in],idx = 0;
	factoradic_radix_vector( dim_in, vec_in, rad );
	for(i=0;i<dim_in;i++)
		idx += rad[dim_in-1-i] * factorial( i );
	return idx;
}

int factoradic_init( int dim_in, int *vec_in )
{
	int i;

	for(i=0;i<dim_in;i++)
		vec_in[i] = i;

	return 0l;
}

int factoradic_next( int dim_in, int *vec_in )
{
	
}

int factoradic_vector( int idx_in, int dim_in, int *vec_out )
{
	int i,j,k,rad[dim_in];
	int fnd[dim_in];

	for(i=0;i<dim_in;i++)
		fnd[i] = 0;

	factoradic_radix_index( dim_in, idx_in, rad );

	for(i=0;i<dim_in;i++)
	{
		for(j=0,k=0;j<dim_in;j++)
		{
			if( fnd[j] == 0 )
				++k;
			if( k - 1 >= rad[i] )
				break;
		}
		fnd[j] = 1;
		vec_out[i] = j;
	}
	return 0;
}

/**
 * Given a factoradic index, this function computes the radix form of the
 * index. That is, it converts the base 10 input index into a non-constant
 * factorial base number.
 */
long int factoradic_radix_index_long( long int dim_in, long int idx_in, long int *rad_out )
{
	long int i,fct,div,rem;

	rem = idx_in;
	for(i=dim_in-1;i>=0;i--)
	{
		fct = factorial( i );
		div = rem / fct;
		rem = rem % fct;
		rad_out[dim_in-1-i] = div;
	}
}

long int factoradic_radix_vector_long( long int dim_in, long int *vec_in, long int *rad_out )
{
	long int c,i,j;
	for(i=0;i<dim_in;i++)
		rad_out[i] = 0;
	for(i=0;i<dim_in;i++)
	{
		for(j=i+1;j<dim_in;j++)
		{
			if( vec_in[j] < vec_in[i] )
				++rad_out[i];
		}
	}
	return 0l;
}

long int factoradic_index_long( long int dim_in, long int *vec_in )
{
	long int i,rad[dim_in],idx = 0;
	factoradic_radix_vector_long( dim_in, vec_in, rad );
	for(i=0;i<dim_in;i++)
		idx += rad[dim_in-1-i] * factorial( i );
	return idx;
}

long int factoradic_init_long( long int dim_in, long int *vec_in )
{
	long int i;

	for(i=0;i<dim_in;i++)
		vec_in[i] = i;

	return 0l;
}

long int factoradic_next_long( long int dim_in, long int *vec_in )
{
	
}

long int factoradic_vector_long( long int idx_in, long int dim_in, long int *vec_out )
{
	long int i,j,k,rad[dim_in];
	int fnd[dim_in];

	for(i=0;i<dim_in;i++)
		fnd[i] = 0;

	factoradic_radix_index_long( dim_in, idx_in, rad );

	for(i=0;i<dim_in;i++)
	{
		for(j=0,k=0;j<dim_in;j++)
		{
			if( fnd[j] == 0 )
				++k;
			if( k - 1 >= rad[i] )
				break;
		}
		fnd[j] = 1;
		vec_out[i] = j;
	}
	return 0;
}

void partition_init( int *s, int *m, int n )
{
	int i;
	for(i=0;i<n;i++)
		m[i] = 1, s[i] = 1;
}

int partition_next( int *s, int *m, int n )
{
	/* Update s: 1 1 1 1 -> 2 1 1 1 -> 1 2 1 1 -> 2 2 1 1 -> 3 2 1 1 -> 1 1 2 1 ... */
	int i = 0;
	++s[i];
	while ((i < n - 1) && (s[i] > m[i] + 1))
	{
		s[i] = 1;
		++i;
		++s[i];
	}

	/* If i is has reached n-1 th element, then the last unique partition has been found*/
	if (i == n - 1)
		return 0;

	/* Because all the first i elements are now 1, s[i] (i + 1 th element)
	is the largest. So we update max by copying it to all the first i
	positions in m.*/
	int max = s[i];
	for (i = i - 1; i >= 0; --i)
		m[i] = max;

	return 1;
}

