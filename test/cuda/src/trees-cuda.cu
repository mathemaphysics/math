/**
 * CUDA device code for stepping through a grid in d dimensions
 */

__device__ int arraynext( long dim_in, long *size_in, long *index_in )
{
	long i;
	for(i=dim_in-1;i>=0;i--)
	{
		if( index_in[i] < size_in[i] )
		{
			++index_in[i];
			return 0;
		}
		else /* Carry to the next spot */
			index_in[i] = 0;
	}
	return -1;
}

