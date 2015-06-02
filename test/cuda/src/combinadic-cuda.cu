/**
 * CUDA device code for combinadic functions
 */

__device__ void partition_init( int *s, int *m, int n )
{
	int i;
	for(i=0;i<n;i++)
		m[i] = 1, s[i] = 1;
}

__device__ int partition_next( int *s, int *m, int n )
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

__device__ int factorial( int num_in )
{
	int fct = 1;
	int i;
	if( num_in <= 0 )
		return 1;
	for(i=num_in;i>1;i--)
		fct *= i;
	return fct;
}

__device__ int truncfact( int n_in, int r_in )
{
	int i,prd = 1;
	for(i=0;i<r_in;i++)
		prd *= n_in - i;
	return prd;
}

__device__ int binomial( int n_in, int k_in )
{
	if( n_in < k_in )
		return 0;
	else
		return truncfact( n_in, k_in ) / factorial( k_in );
}

__device__ int combinadic_init( int lim_in, int dim_in, int *ptr_in )
{
	int i;
	if(lim_in<dim_in)
		return -1; // can't make dim_in-combinations of lim_in (< dim_in) objects
	for(i=0;i<dim_in;i++)
		ptr_in[i] = i;
	return 0;
}

__device__ int combinadic_next( int lim_in, int dim_in, int *vec_in )
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
	return 0;
}

__device__ int combinadic_vector( int idx_in, int lim_in, int dim_in, int *vec_out )
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

__device__ int rcombinadic_vector( int idx_in, int lim_in, int dim_in, int *vec_out )
{
	return combinadic_vector( idx_in, lim_in + dim_in - 1, dim_in - 1, vec_out );
}

__device__ int rcombinadic_occupancy( int lim_in, int dim_in, int *vec_in, int *occ_out )
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

__device__ int polynomial_exponents( int idx_in, int ord_in, int dim_in, int *occ_out )
{
	int *vec = (int*) malloc( ( dim_in-1 ) * sizeof(int) );
	rcombinadic_vector( idx_in, ord_in, dim_in, vec );
	rcombinadic_occupancy( ord_in, dim_in, vec, occ_out );
	free( vec );
	return 0;
}

__device__ int global_polynomial_vector( int idx_in, int dim_in, int *exp_out )
{
	int n,m,r;
	int *cmb = (int*) malloc( ( dim_in - 1 ) * sizeof(int) );

	for(n=0,r=0;;r++)
	{
		m = binomial( r + dim_in - 1, dim_in - 1 );
		if( n + m > idx_in )
			break;
		else
			n += m;
	}

	polynomial_exponents( idx_in - n, r, dim_in, exp_out );

	free( cmb );

	return 0;
}

