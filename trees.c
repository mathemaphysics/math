#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define LIST_SIZE_INC 20

typedef struct
{
	void *start;	///< Pointer to first element
	void *pointer;	///< The start of allocated memory
	int step;	    ///< The size (in bytes) of an element
	int size;	    ///< The number of elements stored
	int alloc;	    ///< The number of elements allocated
} list_t;

typedef char byte_t;

#define HEAP_SIZE_INC 512 ////////////

typedef struct
{
	void *start;
	int step;
	long size;
	long alloc;
	int (*cmp)( void *, void * );
} heap_t;

/**
 * This is the internal node data structure used within the 
 * linked list implementation.
 * \brief A node in a list
 */
typedef struct
{
	/**
	 * The customized data stored in the node
	 */
	void *data;

	/**
	 * Integer into the global list containing all nodes
	 * pointing to the location of the next node in the list
	 */
	long next;

	/**
	 * Integer pointing to the 
	 */
	long last;

	/**
	 * Integer indicating the type of node this is; this can
	 * be customized
	 */
	int type;

	/**
	 * This is a binary variable; 0 indicates available for
	 * use, 1 indicates in use
	 */
	char avail;
} node_t;

#define LLIST_SIZE_INC 20

typedef struct
{
	/**
	 * Start of the nodes data structures
	 */
	node_t *start;

	/**
	 * The size of the data stored in bytes
	 */
	long step;

	/**
	 * The index of the first element
	 */
	long begin;

	/**
	 * The index of the last element
	 */
	long end;

	/**
	 * The length of the list
	 */
	long size;

	/**
	 * The space currently allocated to start
	 */
	long alloc;

	/**
	 * Optionally stores the index of each node in the list from reference at start
	 */
	long *index;

	/**
	 * Function used to compare two data points
	 */
	int (*cmp)( void *, void * );
} llist_t;

#define llist_entry( list_in, i ) list_in.start[i].data
#define llist_print_entry( list_in, j, type, printfunc ) printfunc( (type*) llist_entry(list_in,j) )
#define BTREE_SIZE_INC 512

typedef struct
{
	/**
	 * Data which is particular to this node; can be anything
	 */
	void *data;

	/**
	 * The key associated with this node
	 */
	long key;

	/**
	 * Availability of the node; 0 for available
	 */
	char avail;

	/**
	 * Index of the parent of this node
	 */
	long parent;

	/**
	 * Index of the left child of this node
	 */
	long left;

	/**
	 * Index of the right child of this node
	 */
	long right;

	/**
	 * The currently color of this node
	 */
	char color;

	/**
	 * Global position of this data structure in memory; this
	 * may not be necessary but it is useful
	 */
	long index;
} bsnode_t;

typedef struct
{
	bsnode_t *start;
	long root;
	long step;
	long size;
	long alloc;
	int (*compf)(void*,void*);
	long min;
	long max;
} bstree_t;

typedef union
{
	void *pointer;
	long integer;
} target_t;

typedef struct
{
	void *data;
	long key;
	char avail;
	long parent;
	long left;
	long right;
	long priority;
} tpnode_t;

typedef struct
{
	tpnode_t *start;
	long root;
	long step;
	long size;
	long alloc;
} treap_t;

int binary_search_range( long int *, long int *, long int *, long int *, long int *, long int *, long int * );
int binary_closest_range( long int *, long int *, long int *, long int *, long int *, long int *, long int * );
int heap_grow( heap_t *, long );
int llist_grow( llist_t *, long );
int bstree_grow( bstree_t *, long );
int bstree_left_rotate( bstree_t *, long );
int bstree_right_rotate( bstree_t *, long );
int bstree_minimum( bstree_t *, long, long * );
int bstree_maximum( bstree_t *, long, long * );
int bstree_search( bstree_t *, long, target_t, long * );
int bstree_delete_index( bstree_t *, long, long * );
int bstree_delete_index_balanced( bstree_t *, long );

int arraynext( long dim_in, long *size_in, long *index_in )
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

long arrayindex( long dim_in, long *size_in, long *index_in )
{
        long i,j,k;

        /* The total index to output */
        long idx = 0;

        /* Sum contributions from each index */
        for(i=0;i<dim_in;i++)
        {
                for(j=i+1,k=1;j<dim_in;j++)
                        k *= size_in[j];
                idx += k * index_in[i];
        }

        /* Return index into one-dimensional array */
        return idx;
}

int binary_search( long int *i_in, long int *list_in, long int *llen_in, long int *idx_out, long int *its_out )
{
	long int left,right;
	left = 0;
	right = *llen_in - 1;
	return binary_search_range( i_in, &left, &right, list_in, llen_in, idx_out, its_out );
}

int binary_search_range( long int *i_in, long int *left_in, long int *right_in, long int *list_in, long int *llen_in, long int *idx_out, long int *its_out )
{
	long int left,right;
	*its_out = 0;
	if( *llen_in > 0 )
	{
		left = *left_in;
		right = *right_in;
		while(left<=right)
		{
			*idx_out = ( left + right ) / 2;
			if( list_in[*idx_out] == *i_in )
				return 0;
			else
			{
				if( list_in[*idx_out] < *i_in )
					left = *idx_out + 1;
				else
					right = *idx_out - 1;
			}
			*its_out = *its_out + 1;
		}
		return -1;
	}
	else
	{
		*idx_out = 0;
		return -1;
	}
}

int binary_closest( long int *i_in, long int *list_in, long int *llen_in, long int *idx_out, long int *its_out )
{
	long int left,right;
	left = 0;
	right = *llen_in - 1;
	return binary_closest_range( i_in, &left, &right, list_in, llen_in, idx_out, its_out );
}

int binary_closest_range( long int *i_in, long int *left_in, long int *right_in, long int *list_in, long int *llen_in, long int *idx_out, long int *its_out )
{
	long int idx;
	if( binary_search_range( i_in, left_in, right_in, list_in, llen_in, &idx, its_out ) == -1 )
	{
		if( *llen_in > 0 && idx > -1 )
			*idx_out = list_in[idx] > *i_in ? idx - 1 : idx;
		else
			*idx_out = -1;
		return -1;
	}
	else
	{
		*idx_out = idx;
		return 0;
	}
}

/**
 * Just a simple search-for-integer function. Returns -1 if not found.
 */
int simple_find( int target_in, int len_in, int *set_in )
{
	int i;
	for(i=0;i<len_in;i++)
	{
		if( set_in[i] == target_in )
			return i;
	}
	return -1;
}

int simple_is_subset( int len1_in, int *set1_in, int len2_in, int *set2_in )
{
	int i,j;
	if( len1_in > len2_in )
		return 0; /* If set1_in is bigger than set2_in then set1_in cannot be a subset of set2_in */
	for(i=0;i<len1_in;i++)
	{
		if( simple_find( set1_in[i], len2_in, set2_in ) == -1 )
			return 0; // set1_in is NOT a subset of set2_in
	}
	return 1; // set1_in is a subset of set2_in
}

int simple_union( int len1_in, int *set1_in, int len2_in, int *set2_in, int *set_out )
{
	int i,len;

	len = 0;
	for(i=0;i<len1_in;i++)
		if( simple_find( set1_in[i], len, set_out ) == -1 )
			set_out[len] = set1_in[i], ++len;
	for(i=0;i<len2_in;i++)
		if( simple_find( set2_in[i], len, set_out ) == -1 )
			set_out[len] = set2_in[i], ++len;

	return len;
}

int simple_set_intersect( int len1_in, int *set1_in, int len2_in, int *set2_in, int *set_out )
{
	int i,j;
	for(i=0,j=0;i<len1_in;i++)
	{
		if( simple_find( set1_in[i], len2_in, set2_in ) != -1 )
		{
			set_out[j] = set1_in[i];
			++j;
		}
	}
	return j;
}

void bubble_sort( int len_in, int *set_in )
{
	int i,j,tmp;
	if( len_in > 1 )
	{
		for(i=len_in-1;i>0;i--)
		{
			for(j=0;j<i;j++)
			{
				if( set_in[j] > set_in[j+1] )
				{
					tmp = set_in[j];
					set_in[j] = set_in[j+1];
					set_in[j+1] = tmp;
				}
			}
		}
	}
}

void complex_bubble_sort( int len_in, double *set_in )
{
        int i,j;
        double tmpr,tmpi;
        if( len_in > 1 )
        {
                for(i=len_in-1;i>0;i--)
                {
                        for(j=0;j<i;j++)
                        {
                                if( set_in[2*j] < set_in[2*(j+1)] )
                                {
                                        tmpr = set_in[2*j];
                                        set_in[2*j] = set_in[2*(j+1)];
                                        set_in[2*(j+1)] = tmpr;
                                        tmpi = set_in[2*j+1];
                                        set_in[2*j+1] = set_in[2*(j+1)+1];
                                        set_in[2*(j+1)+1] = tmpi;
                                }
                        }
                }
        }
}

void inverse_complex_bubble_sort( int len_in, double *set_in )
{
	int i,j;
        double tmpr,tmpi;
        if( len_in > 1 )
        {
                for(i=len_in-1;i>0;i--)
                {
                        for(j=0;j<i;j++)
                        {
                                if( 1.0 / set_in[2*j] > 1.0 / set_in[2*(j+1)] )
                                {
                                        tmpr = set_in[2*j];
                                        set_in[2*j] = set_in[2*(j+1)];
                                        set_in[2*(j+1)] = tmpr;
                                        tmpi = set_in[2*j+1];
                                        set_in[2*j+1] = set_in[2*(j+1)+1];
                                        set_in[2*(j+1)+1] = tmpi;
                                }
                        }
                }
        }
}

/**
 * Just a simple search-for-integer function. Returns -1 if not found.
 * This one is just the long version.
 */
long int simple_find_long( long int target_in, long int len_in, long int *set_in )
{
	long int i;
	for(i=0;i<len_in;i++)
	{
		if( set_in[i] == target_in )
			return i;
	}
	return -1;
}

long int simple_is_subset_long( long int len1_in, long int *set1_in, long int len2_in, long int *set2_in )
{
	long int i;
	if( len1_in > len2_in )
		return 0; /* If set1_in is bigger than set2_in then set1_in cannot be a subset of set2_in */
	for(i=0;i<len1_in;i++)
	{
		if( simple_find_long( set1_in[i], len2_in, set2_in ) == -1 )
			return 0; // set1_in is NOT a subset of set2_in
	}
	return 1; // set1_in is a subset of set2_in
}

long int simple_union_long( long len1_in, long *set1_in, long len2_in, long *set2_in, long *set_out )
{
	long i,len;

	len = 0;
	for(i=0;i<len1_in;i++)
		if( simple_find_long( set1_in[i], len, set_out ) == -1 )
			set_out[len] = set1_in[i], ++len;
	for(i=0;i<len2_in;i++)
		if( simple_find_long( set2_in[i], len, set_out ) == -1 )
			set_out[len] = set2_in[i], ++len;

	return len;
}

long int simple_set_intersect_long( long int len1_in, long int *set1_in, long int len2_in, long int *set2_in, long int *set_out )
{
	long int i,j;
	for(i=0,j=0;i<len1_in;i++)
	{
		if( simple_find_long( set1_in[i], len2_in, set2_in ) != -1 )
		{
			set_out[j] = set1_in[i];
			++j;
		}
	}
	return j;
}

void bubble_sort_long( long int len_in, long int *set_in )
{
	long int i,j,tmp;
	if( len_in > 1 )
	{
		for(i=len_in-1;i>0;i--)
		{
			for(j=0;j<i;j++)
			{
				if( set_in[j] > set_in[j+1] )
				{
					tmp = set_in[j];
					set_in[j] = set_in[j+1];
					set_in[j+1] = tmp;
				}
			}
		}
	}
}

void bubble_sort_double_long( long int len_in, double *set_in )
{
	long int i,j;
	double tmp;
	if( len_in > 1 )
	{
		for(i=len_in-1;i>0;i--)
		{
			for(j=0;j<i;j++)
			{
				if( set_in[j] > set_in[j+1] )
				{
					tmp = set_in[j];
					set_in[j] = set_in[j+1];
					set_in[j+1] = tmp;
				}
			}
		}
	}
}

/**
 * Create a new basic list.
 * @param list_in The list structure
 * @param step_in Size of the associated data type
 * @return Returns 0 if all is well or -1 if there was a problem
 */
int list_init( list_t *list_in, int step_in )
{
	void *tmp;

	list_in->start = NULL;
	list_in->pointer = NULL;
	list_in->step = step_in;
	list_in->size = 0;
	list_in->alloc = LIST_SIZE_INC;
	tmp = realloc( list_in->pointer, LIST_SIZE_INC * step_in );
	if( tmp == NULL )
		return -1;
	list_in->pointer = tmp;
	list_in->start = list_in->pointer;
	return 0;
}

int list_free( list_t *list_in )
{
	if( list_in->start != NULL && list_in->size > 0 )
		free( list_in->pointer );
	else
		return -1;
	list_in->pointer = NULL;
	list_in->start = NULL;
	return 0;
}

int list_reset( list_t *list_in )
{
	if( list_in->size <= 0 )
		return -1;
	memmove( list_in->pointer, list_in->start, list_in->size * list_in->step );
	list_in->start = list_in->pointer;
	return 0;
}

int list_shift( list_t *list_in, int shift_in )
{
	if( list_in->size <= 0
	 || ( ( (byte_t*) list_in->start + shift_in * list_in->step ) < (byte_t*) list_in->pointer )
	 || ( ( (byte_t*) list_in->start + ( list_in->size + shift_in ) * list_in->step ) > ( (byte_t*) list_in->pointer + list_in->alloc * list_in->step ) ) )
		return -1;
	if( shift_in == 0 )
		return 0;
	list_in->start = memmove( (byte_t*) list_in->start + shift_in * list_in->step, list_in->start, list_in->size * list_in->step );
	return 0;
}

int list_shift_set( list_t *list_in, int start_in, int len_in, int shift_in )
{
	if( (byte_t*) list_in->start + ( start_in + shift_in ) * list_in->step < (byte_t*) list_in->pointer )
		return -1;
	if( (byte_t*) list_in->start + ( start_in + len_in + shift_in ) * list_in->step > (byte_t*) list_in->pointer + list_in->alloc * list_in->step )
		return -1;
	if( shift_in == 0 )
		return 0;

	/* Just move the data; do not change start */
	memmove( (byte_t*) list_in->start + ( start_in + shift_in ) * list_in->step, (byte_t*) list_in->start + start_in * list_in->step, len_in * list_in->step );

	return 0;
}

int list_grow( list_t *list_in, int minc_in )
{
	void *tmp;
	list_reset( list_in );
	tmp = realloc( list_in->pointer, ( list_in->alloc + minc_in ) * list_in->step );
	if( tmp == NULL )
		return -1;
	list_in->pointer = tmp;
	list_in->start = list_in->pointer;
	list_in->alloc += minc_in;
	return list_in->alloc;
}

int list_shrink( list_t *list_in, int mdec_in )
{
	void *tmp;
	list_reset( list_in ); /* Force list_in->start to coincide with list_in->pointer */
	tmp = realloc( list_in->pointer, ( list_in->alloc - mdec_in ) * list_in->step );
	if( tmp == NULL )
		return -1;
	list_in->pointer = tmp;
	list_in->start = list_in->pointer;
	list_in->alloc -= mdec_in;
	return list_in->alloc;
}

int list_prepend( list_t *list_in, void *val_in )
{
	if( list_in->size + 1 > list_in->alloc )
	{
		if( list_grow( list_in, LIST_SIZE_INC ) == -1 )
			return -1;
		if( list_shift( list_in, 1 ) == -1 )
			return -1;
	}
	else
	{
		if( list_in->start <= list_in->pointer ) /* Meaning there is no memory allocated to the left, so make some */
			if( list_shift( list_in, 1 ) == -1 )
				return -1;
	}
	list_in->start = (byte_t*) list_in->start - list_in->step;
	bcopy( val_in, list_in->start, list_in->step );
	++list_in->size;
	return list_in->size;
}

int list_append( list_t *list_in, void *val_in )
{
	/* If there isn't enough memory allocated */
	if( list_in->size + 1 > list_in->alloc )
	{
		if( list_grow( list_in, LIST_SIZE_INC ) == -1 )
			return -1;
	}
	else /* If there is enough, then either there is some at the beginning or some at the end */
	{
		if( (byte_t*)list_in->start > (byte_t*)list_in->pointer )
			if( list_shift( list_in, -1 ) == -1 )
				return -1;
	}

	bcopy( val_in, (byte_t*) list_in->start + list_in->size * list_in->step, list_in->step );
	++list_in->size;
	return list_in->size;
}

int list_popf( list_t *list_in, void *ptr_out )
{
	if( list_in->size <= 0 )
		return -1;
	bcopy( list_in->start, ptr_out, list_in->step );
	if( list_in->size > 1 )
		list_in->start = (byte_t*) list_in->start + list_in->step;
	--list_in->size;
	if( list_in->alloc - list_in->size >= LIST_SIZE_INC )
		if( list_shrink( list_in, LIST_SIZE_INC ) == -1 )
			return -1;
	return 0;
}

int list_popb( list_t *list_in, void *ptr_out )
{
	if( list_in->size <= 0 )
		return -1;
	else
		bcopy( (byte_t*) list_in->start + ( list_in->size - 1 ) * list_in->step, ptr_out, list_in->step );
	--list_in->size;
	if( list_in->alloc - list_in->size >= LIST_SIZE_INC )
		if( list_shrink( list_in, LIST_SIZE_INC ) == -1 )
			return -1;
	return 0;
}

int list_insert( list_t *list_in, void *ptr_in, int idx_in )
{
	if( list_in->size + 1 > list_in->alloc )
	{
		if( list_grow( list_in, LIST_SIZE_INC ) == -1 )
			return -1;
	}
	if( (byte_t*) list_in->start - (byte_t*) list_in->pointer > 0 )
	{
		if( list_shift_set( list_in, 0, idx_in + 1, -1 ) == -1 )
			return -1;
		bcopy( ptr_in, (byte_t*) list_in->start + idx_in * list_in->step, list_in->step ); /* Remember actually inserting into idx_in + 1 position */
		list_in->start = (byte_t*) list_in->start - list_in->step;
	}
	else
	{
		if( list_shift_set( list_in, idx_in + 1, list_in->size - idx_in - 1, 1 ) == -1 )
			return -1;
		bcopy( ptr_in, (byte_t*) list_in->start + ( idx_in + 1 ) * list_in->step, list_in->step );
	}
	++list_in->size;
	return 0;
}

int list_insert_set( list_t *list_in, void *ptr_in, int idx_in, int len_in )
{
	
}

int list_delete( list_t *list_in, int idx_in )
{
	/* Decide which shift requires movement of smallest amount of data, from left or from right */
	if( idx_in > list_in->size - idx_in ) /* from right? */
	{
		if( list_shift_set( list_in, idx_in + 1, list_in->size - idx_in - 1, -1 ) == -1 )
			return -1;
	}
	else /* else from left */
	{
		if( list_shift_set( list_in, 0, idx_in, 1 ) == -1 )
			return -1;
		list_in->start = (byte_t*) list_in->start + list_in->step; /* Must manually adjust start */
	}
	--list_in->size;
	if( list_in->alloc - list_in->size >= LIST_SIZE_INC )
		if( list_shrink( list_in, LIST_SIZE_INC ) == -1 )
			return -1;
	return 0;
}

int list_delete_set( list_t *list_in, int idx_in, int len_in )
{
	
}

/**
 * Create a new heap structure.
 * @param heap_in The data structure
 * @param step_in The size of data in bytes
 * @param cmp_in Comparison function for the particular data type stored
 * @return 
 */
int heap_init( heap_t *heap_in, int step_in, int (*cmp_in)(void*,void*) )
{
	heap_in->cmp = cmp_in; /* Set the order relation function */
	heap_in->start = malloc( HEAP_SIZE_INC * step_in );
	heap_in->step = step_in;
	heap_in->size = 0;
	heap_in->alloc = 0;
	heap_grow( heap_in, HEAP_SIZE_INC );
	return 0;
}

int heap_free( heap_t *heap_in )
{
	free( heap_in->start );
	heap_in->size = 0;
	heap_in->alloc = 0;
	return 0;
}

int heap_append( heap_t *heap_in, void *data_in )
{
	if( heap_in->alloc < heap_in->size + 1 )
		heap_grow( heap_in, HEAP_SIZE_INC );
	bcopy( data_in, (byte_t*) heap_in->start + heap_in->step * heap_in->size, heap_in->step );
	++heap_in->size;
	return 0;
}

long heap_left( long idx_in ) /* Using zero-based indexing implicitly, converting to one-based and back */
{
	return ( ( idx_in + 1 ) << 1 ) - 1;
}

long heap_right( long idx_in ) /* Using zero-based indexing implicitly, converting to one-based and back */
{
	return ( ( idx_in + 1 ) << 1 ); /* Implicitly adding 1 by NOT subtracting 1 */
}

long heap_parent( long idx_in ) /* Using zero-based indexing implicitly, converting to one-based and back */
{
	return ( ( idx_in + 1 ) >> 1 ) - 1;
}

int heap_insert( heap_t *heap_in, void *data_in )
{
	int i,j,step = heap_in->step;
	byte_t *start = (byte_t*) heap_in->start;
	if( heap_in->alloc < heap_in->size + 1 )
		heap_grow( heap_in, HEAP_SIZE_INC );
	i = heap_in->size;
	while( i > 0 )
	{
		j = heap_parent( i );
		if( heap_in->cmp( start + j * step, data_in ) == 1 )
			break;
		bcopy( start + j * step, start + i * step, step );
		i = j; /* Take a step up the tree toward the root */
	}
	bcopy( data_in, start + i * step, step );
	++heap_in->size;
	return 0;
}

int heap_delete( heap_t *heap_in, long idx_in )
{
	return 0;
}

int heap_swap( heap_t *heap_in, long idx_in, long jdx_in )
{
	if( idx_in >= heap_in->size || jdx_in >= heap_in->size )
		return -1;
	byte_t *start = (byte_t*) heap_in->start;
	int step = heap_in->step;
	byte_t tmp[step];
	bcopy( start + idx_in * step, tmp, step );
	bcopy( start + jdx_in * step, start + idx_in * step, step );
	bcopy( tmp, start + jdx_in * step, step );
	return 0;
}

int heap_heapify( heap_t *heap_in, long idx_in )
{
	byte_t *start = (byte_t*) heap_in->start;
	int step = heap_in->step;
	long m = idx_in;
	const long l = heap_left( idx_in );
	const long r = heap_right( idx_in );
	if( ( l < heap_in->size ) && heap_in->cmp( start + l * step, start + idx_in * step ) == 1 )
		m = l;
	if( ( r < heap_in->size ) && heap_in->cmp( start + r * step, start + m * step ) == 1 )
		m = r;
	if( m != idx_in )
	{
		if( heap_swap( heap_in, idx_in, m ) != 0 )
			return -1;
		if( heap_heapify( heap_in, m ) != 0 )
			return -2;
	}
	return 0;
}

int heap_build( heap_t *heap_in )
{
	long j = ( heap_in->size >> 1 ); /* divide by 2... kind of */
	while( j >= 0 )
	{
		heap_heapify( heap_in, j );
		--j;
	}
	return 0;
}

int heap_find( heap_t *heap_in, void *data_in )
{
	return 0;
}

int heap_grow( heap_t *heap_in, long inc_in )
{
	void *ptr;
	ptr = realloc( heap_in->start, ( heap_in->alloc + inc_in ) * heap_in->step );
	if( ptr == NULL )
		return -1;
	heap_in->start = ptr;
	heap_in->alloc = heap_in->alloc + inc_in;
	return 0;
}

int heap_shrink( heap_t *heap_in, long dec_in )
{
	return 0;
}

void heap_print( heap_t *heap_in )
{
	double size = (double) heap_in->size;
	long i,j,lim = (long) ceil(log(size)/log(2.0));
	byte_t *start = (byte_t*) heap_in->start;
	int step = heap_in->step;
	for(i=0;i<lim;i++)
	{
		for(j=(1<<i)-1;j<(1<<(i+1))-1&&j<heap_in->size;j++)
			printf( "%15.7f ", ((double*)start)[j] );
		printf( "\n" );
	}
}

/**
 * Create a new linked list data type.
 * @param list_in The list object
 * @param step_in The size of the data in bytes
 * @param cmp_in Optional function which compares the data in each node to determine greater than or less than (use NULL for none)
 */
int llist_init( llist_t *list_in, int step_in, int (*cmp_in)(void*,void*) )
{
	list_in->start = NULL;
	list_in->step = step_in;
	list_in->cmp = cmp_in;
	list_in->size = 0;
	list_in->alloc = 0;
	list_in->begin = 0;
	list_in->end = 0;
	llist_grow( list_in, LLIST_SIZE_INC );
}

int llist_free( llist_t *list_in )
{
	int i;
	for(i=0;i<list_in->alloc;i++)
		free( list_in->start[i].data );
	free( list_in->start );
	list_in->size = 0;
	list_in->alloc = 0;
	list_in->begin = 0;
	list_in->end = 0;
}

int llist_avail( llist_t *list_in, long *idx_out )
{
	long i;
	for(i=0;i<list_in->alloc;i++)
	{
		if( list_in->start[i].avail == 0 )
		{
			*idx_out = i;
			break;
		}
	}
	return 0;
}

int llist_prepend( llist_t *list_in, void *data_in )
{
	long ins;
	if( list_in->alloc < list_in->size + 1 )
		llist_grow( list_in, LLIST_SIZE_INC );
	if( list_in->size > 0 )
	{
		llist_avail( list_in, &ins );
		list_in->start[list_in->begin].last = ins;
		bcopy( data_in, list_in->start[ins].data, list_in->step );
		list_in->start[ins].last = -1;
		list_in->start[ins].next = list_in->begin;
		list_in->begin = ins;
		list_in->start[ins].avail = 1;
	}
	else
	{
		list_in->start[0].last = -1;
		list_in->start[0].next = -1;
		bcopy( data_in, list_in->start[0].data, list_in->step );
		list_in->begin = 0;
		list_in->end = 0;
		list_in->start[0].avail = 1;
		ins = 0;
	}
	++list_in->size;
	return ins;
}

int llist_append( llist_t *list_in, void *data_in )
{
	long ins;
	if( list_in->alloc < list_in->size + 1 )
		llist_grow( list_in, LLIST_SIZE_INC );
	if( list_in->size > 0 )
	{
		llist_avail( list_in, &ins );
		list_in->start[list_in->end].next = ins;
		bcopy( data_in, list_in->start[ins].data, list_in->step );
		list_in->start[ins].last = list_in->end;
		list_in->start[ins].next = -1;
		list_in->end = ins;
		list_in->start[ins].avail = 1; /* set it to occupied */
	}
	else
	{
		list_in->start[0].last = -1;
		list_in->start[0].next = -1;
		bcopy( data_in, list_in->start[0].data, list_in->step );
		list_in->begin = 0;
		list_in->end = 0;
		list_in->start[0].avail = 1; /* set it to occupied */
		ins = 0;
	}
	++list_in->size;
	return ins;
}

int llist_popb( llist_t *list_in, void *data_out )
{
	if( list_in->size <= 0 )
		return -1;
	long last = list_in->start[list_in->end].last;
	bcopy( list_in->start[list_in->end].data, data_out, list_in->step );
	list_in->start[list_in->end].avail = 0; /* set to available */
	list_in->end = list_in->start[list_in->end].last; /* end is now set to the previous */
	list_in->start[last].next = -1; /* set past-the-end pointer to null */
	--list_in->size;
	return 0;
}

int llist_popf( llist_t *list_in, void *data_out )
{
	
}

int llist_insert_after( llist_t *list_in, long idx_in, void *data_in )
{
	if( idx_in >= list_in->size || idx_in < 0 )
		return -1;
	if( list_in->alloc < list_in->size + 1 )
		llist_grow( list_in, LLIST_SIZE_INC );
	long ins;
	long after = list_in->start[idx_in].next;
	llist_avail( list_in, &ins );
	if( after != -1 )
		list_in->start[after].last = ins;
	else
		list_in->end = ins;
	list_in->start[idx_in].next = ins;
	bcopy( data_in, list_in->start[ins].data, list_in->step );
	list_in->start[ins].next = after;
	list_in->start[ins].last = idx_in;
	list_in->start[ins].avail = 1; /* not available for use */
	return list_in->size++; /* Return first, then increment */
}

int llist_insert_before( llist_t *list_in, long idx_in, void *data_in )
{
	if( idx_in >= list_in->size || idx_in < 0 )
		return -1;
	if( list_in->alloc < list_in->size + 1 )
		llist_grow( list_in, LLIST_SIZE_INC );
	long ins;
	long before = list_in->start[idx_in].last;
	llist_avail( list_in, &ins );
	if( before != -1 )
		list_in->start[before].next = ins;
	else
		list_in->begin = ins;
	list_in->start[idx_in].last = ins;
	bcopy( data_in, list_in->start[ins].data, list_in->step );
	list_in->start[ins].next = idx_in;
	list_in->start[ins].last = before;
	list_in->start[ins].avail = 1; /* not available for use */
	return list_in->size++;
}

int llist_delete( llist_t *list_in, long idx_in )
{
	if( idx_in >= list_in->size || idx_in < 0 )
		return -1;
	long last = list_in->start[idx_in].last;
	long next = list_in->start[idx_in].next;
	if( last != -1 )
		list_in->start[last].next = next;
	else
		list_in->begin = next;
	if( next != -1 )
		list_in->start[next].last = last;
	else
		list_in->end = last;
	list_in->start[idx_in].avail = 0; /* this node is now available for use */
	return --list_in->size;
}

int llist_reorder( llist_t *list_in, llist_t *list_out )
{
	long i = list_in->begin;
	long j = 0;
	while( i != -1 && j < list_in->size )
	{
		if( llist_append( list_out, list_in->start[i].data ) < 0 )
			return -1;
		i = list_in->start[i].next;
		++j;
	}
	return 0;
}

int llist_grow( llist_t *list_in, long inc_in )
{
	int i;
	void *ptr;
	node_t *tmp;

	/* Allocate the nodes */
	tmp = (node_t*) realloc( list_in->start, ( list_in->alloc + inc_in + 1 ) * sizeof(node_t) );
	if( tmp == NULL )
		return -1;
	list_in->start = tmp;

	/* Now allocate data space for each node */
	for(i=list_in->alloc;i<list_in->alloc+inc_in+1;i++)
		list_in->start[i].data = malloc( list_in->step ), list_in->start[i].avail = 0; /* yes, available */

	/* Increment the shit */
	list_in->alloc += inc_in;

	return 0;
}

int llist_shrink( llist_t *list_in, long dec_in )
{
	int i;
	node_t *ptr;

	if( list_in->alloc < dec_in )
		dec_in = list_in->alloc; /* Free it all then */

	/* Deallocate the nodes */
	for(i=list_in->alloc-1;i>list_in->alloc-dec_in-1;i--)
		free( list_in->start[i].data );
	ptr = (node_t*) realloc( list_in->start, ( list_in->alloc - dec_in ) * list_in->step );
	if( ptr == NULL )
		return -1;
	list_in->start = ptr;
	list_in->alloc -= dec_in;
	if( list_in->size > list_in->alloc )
		list_in->size = list_in->alloc;
	return list_in->alloc;
}

int llist_find( llist_t *list_in, void *val_in, long *idx_out )
{
	long i = 0;
	long cur = list_in->begin;
	while( cur != -1 && i < list_in->size )
	{
		if( list_in->cmp( val_in, list_in->start[cur].data ) == 0 )
		{
			*idx_out = cur;
			return 0;
		}
		cur = list_in->start[cur].next;
		++i;
	}
	*idx_out = -1; /* Not found */
	return -1;
}

void llist_print( llist_t *list_in )
{
	int num = 0;
	int nxt = list_in->begin;
	while( num < list_in->size )
	{
		printf( "Entry %d, Index %d: %x\n", num, nxt, ((long*) (list_in->start[nxt].data)) );
		nxt = list_in->start[nxt].next;
		++num;
	}
}

#define BLACK 0
#define RED 1

/**
 * Initialize a node data structure for use in a binary
 * search tree.
 * @param node_in The node to set up and allocate
 * @param step_in The number of bytes to use to store data
 * @return Returns 0 if all is well
 */
int bsnode_init( bsnode_t *node_in, long step_in )
{
	node_in->data = malloc( step_in );
	node_in->key = -1;
	node_in->avail = 0; /* mark initially unoccupied */
	node_in->parent = -1;
	node_in->left = -1;
	node_in->right = -1;
	node_in->color = RED; /* Every node will be RED until changed */
	return 0;
}

int bsnode_free( bsnode_t *node_in )
{
	free( node_in->data );
	return 0;
}

int bsnode_reset( bsnode_t *node_in )
{
	node_in->key = -1;
	node_in->avail = 0; /* mark unoccupied */
	node_in->parent = -1;
	node_in->left = -1;
	node_in->right = -1;
	node_in->color = RED; /* Every node will be RED until changed */
	return 0;
}

void bsnode_print( bsnode_t *node_in, void *args_in )
{
	printf( "L: %d\tR: %d\tP: %d\tC: %d\n", node_in->left, node_in->right, node_in->parent, node_in->color );
}

/**
 * There are three cases: start[idx_in] has
 * no children, one child, or two children
 */
#define PARENT_LEFT(i) tree_in->start[tree_in->start[i].parent].left
#define PARENT_RIGHT(i) tree_in->start[tree_in->start[i].parent].right
#define LEFT_PARENT(i) tree_in->start[tree_in->start[i].left].parent
#define RIGHT_PARENT(i) tree_in->start[tree_in->start[i].right].parent
#define PARENT(i) tree_in->start[i].parent
#define LEFT(i) tree_in->start[i].left
#define RIGHT(i) tree_in->start[i].right
#define COLOR(i) tree_in->start[i].color
#define KEY(i) tree_in->start[i].key
#define DATA(i) tree_in->start[i].data

int bstree_init( bstree_t *tree_in, long step_in, int (*compf_in)(void*,void*) )
{
	tree_in->start = NULL;
	tree_in->root = -1;
	tree_in->step = step_in;
	tree_in->size = 0;
	tree_in->alloc = 0;
	tree_in->compf = compf_in;
	tree_in->min = -1;
	tree_in->max = -1;
	if( bstree_grow( tree_in, BTREE_SIZE_INC ) < 0 ) /* All non-success return codes will be < 0 */
		return -1;

	/* allocate the and initialize the NIL node */
	bsnode_init( tree_in->start - 1, tree_in->step );
	tree_in->start[-1].color = BLACK;

	return 0;
}

int bstree_free( bstree_t *tree_in )
{
	int i;
	for(i=-1;i<tree_in->alloc;i++)
		bsnode_free( tree_in->start + i );
	free( tree_in->start - 1 );
	return 0;
}

/**
 * FIXME: The addresses printed with %x specifier
 * here need to be generalized in terms of size, i.e.
 * depending on the system, the pointer size is
 * different
 */
void bstree_info( bstree_t *tree_in )
{
	printf( "%x: start = %x root = %ld step = %ld size = %ld alloc = %ld compf = %x min = %ld max = %ld",
			(unsigned int)tree_in, (unsigned int)tree_in->start, tree_in->root,
			tree_in->step, tree_in->size, tree_in->alloc,
			(unsigned int)tree_in->compf,
			tree_in->min, tree_in->max );
}

int bstree_avail( bstree_t *tree_in, long *node_out )
{
	long i;
	for(i=0;i<tree_in->alloc;i++)
	{
		if( tree_in->start[i].avail == 0 )
		{
			*node_out = i;
			return 0;
		}
	}
	return -1;
}

int bstree_insert( bstree_t *tree_in, void *data_in, long key_in, long *idx_out )
{
	/* Grow with bstree_grow() if necessary */
	if( tree_in->alloc < tree_in->size + 1 )
		if( bstree_grow( tree_in, BTREE_SIZE_INC ) < 0 )
			return -1;

	/* It's a sick sad world... I hate my fucking life */
	char usecomp;
	long cur = tree_in->root;
	long temp = -1;
	if( tree_in->compf == NULL )
	{
		usecomp = 0;
		while( cur != -1 )
		{
			temp = cur;
			if( key_in < tree_in->start[cur].key )
				cur = tree_in->start[cur].left;
			else
				cur = tree_in->start[cur].right;
		}
	}
	else
	{
		usecomp = 1;
		while( cur != -1 )
		{
			temp = cur;
			if( tree_in->compf( data_in, tree_in->start[cur].data ) < 0 )
				cur = tree_in->start[cur].left;
			else
				cur = tree_in->start[cur].right;
		}
	}

	/* Set the data */
	bstree_avail( tree_in, &cur ); /* Find a free spot and bsnode_reset() it */
	if( idx_out != NULL )
		*idx_out = cur;
	bsnode_reset( tree_in->start + cur ); /* makes sure left and right are NULL and key has no value until set */
	tree_in->start[cur].parent = temp;
	tree_in->start[cur].key = key_in; /* regardless of order relation, copy the key */
	if( data_in != NULL )
		bcopy( data_in, tree_in->start[cur].data, tree_in->step );
	else if( usecomp == 1 ) /* use compare function, compf */
		return -1; /* fail if data_in is NULL and usecomp indicates that order relation relies on data only, not key */
	tree_in->start[cur].avail = 1;

	/* now set cur's location relative to rest of tree */
	if( temp != -1 )
	{
		if( usecomp == 0 )
		{
			if( key_in < tree_in->start[temp].key )
				tree_in->start[temp].left = cur;
			else
				tree_in->start[temp].right = cur;
		}
		else
		{
			if( tree_in->compf( data_in, tree_in->start[temp].data ) < 0 )
				tree_in->start[temp].left = cur;
			else
				tree_in->start[temp].right = cur;
		}
	}
	else
		tree_in->root = cur;

	/* if added successfully increment size */
	++tree_in->size;
	
	/* now set min and max */
	if( tree_in->size > 1 )
	{
		if( tree_in->compf == NULL )
		{
			if( tree_in->start[cur].key < tree_in->start[tree_in->min].key )
				tree_in->min = cur;
			else if( tree_in->start[cur].key > tree_in->start[tree_in->max].key )
				tree_in->max = cur;
		}
		else
		{
			if( tree_in->compf( tree_in->start[cur].data, tree_in->start[tree_in->min].data ) < 0 )
				tree_in->min = cur;
			else if( tree_in->compf( tree_in->start[cur].data, tree_in->start[tree_in->max].data ) > 0 )
				tree_in->max = cur;
		}
	}
	else
	{
		/* initialize when no comparison can be made */
		tree_in->min = cur;
		tree_in->max = cur;
	}

	return 0;
}

/**
 * Function bstree_insert calls bsnode_reset which sets
 * node color to RED for the newly inserted node, so all
 * is well up there
 */
int bstree_balance_insert( bstree_t *tree_in, long idx_in )
{
	long y,z;
	z = idx_in;

	while( COLOR(PARENT(z)) == RED )
	{
		if( PARENT_LEFT(PARENT(z)) == PARENT(z) )
		{
			y = PARENT_RIGHT(PARENT(z));
			if( COLOR(y) == RED )
			{
				/* then just need to recolor some nodes */
				COLOR(PARENT(z)) = BLACK;
				COLOR(y) = BLACK;
				COLOR(PARENT(PARENT(z))) = RED;
				z = PARENT(PARENT(z));
			}
			else
			{
				if( z == RIGHT(PARENT(z)) ) /* if z is the child direction opposite that of its parent */
				{
					z = PARENT(z);
					bstree_left_rotate( tree_in, z );
				}
				COLOR(PARENT(z)) = BLACK;
				COLOR(PARENT(PARENT(z))) = RED;
				bstree_right_rotate( tree_in, PARENT(PARENT(z)) );
			}
		}
		else
		{
			y = PARENT_LEFT(PARENT(z));
			if( COLOR(y) == RED )
			{
				/* then just need to recolor some nodes */
				COLOR(PARENT(z)) = BLACK;
				COLOR(y) = BLACK;
				COLOR(PARENT(PARENT(z))) = RED;
				z = PARENT(PARENT(z));
			}
			else
			{
				if( z == LEFT(PARENT(z)) ) /* if z is the child direction opposite that of its parent */
				{
					z = PARENT(z);
					bstree_right_rotate( tree_in, z );
				}
				COLOR(PARENT(z)) = BLACK;
				COLOR(PARENT(PARENT(z))) = RED;
				bstree_left_rotate( tree_in, PARENT(PARENT(z)) );
			}
		}
	}
	COLOR(tree_in->root) = BLACK;
	return 0;
}

int bstree_insert_balanced( bstree_t *tree_in, void *data_in, long key_in, long *idx_out )
{
	if( bstree_insert( tree_in, data_in, key_in, idx_out ) < 0 )
		return -1;
	if( bstree_balance_insert( tree_in, *idx_out ) < 0 )
		return -1;
	return 0;
}

int bstree_successor( bstree_t *tree_in, long idx_in, long *idx_out )
{
	long p,q;
	if( tree_in->start[idx_in].right != -1 )
		bstree_minimum( tree_in, tree_in->start[idx_in].right, idx_out );
	else
	{
		q = idx_in;
		p = tree_in->start[q].parent;
		while( p != -1 && q == tree_in->start[p].right )
		{
			q = p;
			p = tree_in->start[p].parent;
		}
		*idx_out = p;
	}
	return 0;
}

int bstree_predecessor( bstree_t *tree_in, long idx_in, long *idx_out )
{
	long p,q;
	if( tree_in->start[idx_in].left != -1 )
		bstree_maximum( tree_in, tree_in->start[idx_in].left, idx_out );
	else
	{
		q = idx_in;
		p = tree_in->start[q].parent;
		while( p != -1 && q == tree_in->start[p].left )
		{
			q = p;
			p = tree_in->start[p].parent;
		}
		*idx_out = p;
	}
	return 0;
}

int bstree_minimum( bstree_t *tree_in, long start_in, long *idx_out )
{
	long cur = start_in;
	while( tree_in->start[cur].left != -1 )
		cur = tree_in->start[cur].left;
	*idx_out = cur;
	return 0;
}

int bstree_maximum( bstree_t *tree_in, long start_in, long *idx_out )
{
	long cur = start_in;
	while( tree_in->start[cur].right != -1 )
		cur = tree_in->start[cur].right;
	*idx_out = cur;
	return 0;
}

int bstree_delete_key( bstree_t *tree_in, target_t target_in, long *idx_out )
{
	if( bstree_search( tree_in, tree_in->root, target_in, idx_out ) >= 0 )
		return bstree_delete_index( tree_in, *idx_out, NULL );
	else
		return -1; /* did not exist so no need to delete it */
	return 0;
}

int bstree_delete_index( bstree_t *tree_in, long idx_in, long *idx_out )
{
	long x,y,z,idx;
	z = idx_in;
	if( tree_in->size > 1 )
	{
		if( idx_in == tree_in->min )
		{
			bstree_successor( tree_in, idx_in, &idx );
			tree_in->min = idx;
		}
		if( idx_in == tree_in->max )
		{
			bstree_predecessor( tree_in, idx_in, &idx );
			tree_in->max = idx;
		}
	}

	if( LEFT(z) == -1 || RIGHT(z) == -1 )
		y = z;
	else
		bstree_successor( tree_in, z, &y );
	if( LEFT(y) != -1 )
		x = LEFT(y);
	else
		x = RIGHT(y);
	PARENT(x) = PARENT(y);
	if( PARENT(y) == -1 )
		tree_in->root = x;
	else
	{
		if( y == LEFT(PARENT(y)) )
			LEFT(PARENT(y)) = x;
		else
			RIGHT(PARENT(y)) = x;
	}
	if( y != z )
	{
		KEY(z) = KEY(y);
		if( DATA(y) != NULL && DATA(z) != NULL )
			bcopy( DATA(y), DATA(z), tree_in->step );
	}
	if( idx_out != NULL )
		*idx_out = y;
	--tree_in->size;
	tree_in->start[y].avail = 0;
	return 0;
}

int bstree_balance_delete( bstree_t *tree_in, long idx_in )
{
	long x,w;
	x = idx_in;
	while( x != tree_in->root && COLOR(x) == BLACK )
	{
		if( PARENT_LEFT(x) == x )
		{
			w = PARENT_RIGHT(x);
			if( COLOR(w) == RED )
			{
				COLOR(w) = BLACK;
				COLOR(PARENT(x)) = RED;
				bstree_left_rotate( tree_in, PARENT(x) );
				w = PARENT_RIGHT(x);
			}
			if( COLOR(LEFT(w)) == BLACK && COLOR(RIGHT(w)) == BLACK )
			{
				/* Decrement the "blackness" of both x and w, making x black and w red */
				COLOR(w) = RED;
				x = PARENT(x);
			}
			else
			{
				if( COLOR(RIGHT(w)) == BLACK )
				{
					COLOR(LEFT(w)) = BLACK; /* Swap color pre-rotate */
					COLOR(w) = RED; /* Swap color pre-rotate */
					bstree_right_rotate( tree_in, w ); /* Rotate to maintain black count */
					w = PARENT_RIGHT(x); /* Set it back to where it was in the tree */
				}
				COLOR(w) = COLOR(PARENT(x));
				COLOR(PARENT(x)) = BLACK;
				COLOR(RIGHT(w)) = BLACK;
				bstree_left_rotate( tree_in, PARENT(x) );
				x = tree_in->root; /* Force termination; probably a better way */
			}
		}
		else
		{
			w = PARENT_LEFT(x);
			if( COLOR(w) == RED )
			{
				COLOR(w) = BLACK;
				COLOR(PARENT(x)) = RED;
				bstree_right_rotate( tree_in, PARENT(x) );
				w = PARENT_LEFT(x);
			}
			if( COLOR(LEFT(w)) == BLACK && COLOR(RIGHT(w)) == BLACK )
			{
				/* Decrement the "blackness" of both x and w, making x black and w red */
				COLOR(w) = RED;
				x = PARENT(x);
			}
			else
			{
				if( COLOR(LEFT(w)) == BLACK )
				{
					COLOR(RIGHT(w)) = BLACK; /* Swap color pre-rotate */
					COLOR(w) = RED; /* Swap color pre-rotate */
					bstree_left_rotate( tree_in, w ); /* Rotate to maintain black count */
					w = PARENT_LEFT(x); /* Set it back to where it was in the tree */
				}
				COLOR(w) = COLOR(PARENT(x));
				COLOR(PARENT(x)) = BLACK;
				COLOR(LEFT(w)) = BLACK;
				bstree_right_rotate( tree_in, PARENT(x) );
				x = tree_in->root; /* Force termination; probably a better way */
			}
		}
	}
	COLOR(x) = BLACK;
	return 0;
}

int bstree_delete_key_balanced( bstree_t *tree_in, target_t target_in, long *idx_out )
{
	if( bstree_search( tree_in, tree_in->root, target_in, idx_out ) >= 0 )
		return bstree_delete_index_balanced( tree_in, *idx_out );
	else
		return -1; /* Did not exist so no need to delete it */
	return 0;
}

int bstree_delete_index_balanced( bstree_t *tree_in, long idx_in )
{
	int ret;
	long temp;
	if( bstree_delete_index( tree_in, idx_in, &temp ) < 0 )
		return -1;
	if( COLOR(temp) == BLACK )
	{
		if( LEFT(temp) != -1 )
			temp = LEFT(temp);
		else
			temp = RIGHT(temp);
		if( bstree_balance_delete( tree_in, temp ) < 0 )
			return -1;
	}
	return 0;
}

int bstree_pop_front( bstree_t *tree_in, void *data_out, long *key_out )
{
	bcopy( tree_in->start[tree_in->min].data, data_out, tree_in->step );
	*key_out = tree_in->start[tree_in->min].key;
	return bstree_delete_index( tree_in, tree_in->min, NULL );
}

int bstree_pop_front_balanced( bstree_t *tree_in, void *data_out, long *key_out )
{
	bcopy( tree_in->start[tree_in->min].data, data_out, tree_in->step );
	*key_out = tree_in->start[tree_in->min].key;
	return bstree_delete_index_balanced( tree_in, tree_in->min );
}

int bstree_pop_back( bstree_t *tree_in, void *data_out, long *key_out )
{
	bcopy( tree_in->start[tree_in->max].data, data_out, tree_in->step );
	*key_out = tree_in->start[tree_in->max].key;
	return bstree_delete_index( tree_in, tree_in->max, NULL );
}

int bstree_pop_back_balanced( bstree_t *tree_in, void *data_out, long *key_out )
{
	bcopy( tree_in->start[tree_in->max].data, data_out, tree_in->step );
	*key_out = tree_in->start[tree_in->max].key;
	return bstree_delete_index_balanced( tree_in, tree_in->max );
}

int bstree_search( bstree_t *tree_in, long begin_in, target_t target_in, long *idx_out )
{
	if( begin_in != -1 )
	{
		if( tree_in->compf == NULL )
		{
			if( tree_in->start[begin_in].key == target_in.integer )
			{
				*idx_out = begin_in;
				return 0; /* zero indicates that it was found */
			}
			else
			{
				if( target_in.integer < tree_in->start[begin_in].key )
					return bstree_search( tree_in, tree_in->start[begin_in].left, target_in, idx_out );
				else
					return bstree_search( tree_in, tree_in->start[begin_in].right, target_in, idx_out );
				
			}
		}
		else
		{
			if( tree_in->compf( tree_in->start[begin_in].data, target_in.pointer ) == 0 )
			{
				*idx_out = begin_in;
				return 0; /* zero indicates that it was found */
			}
			else
			{
				if( tree_in->compf( target_in.pointer, tree_in->start[begin_in].data ) < 0 )
					return bstree_search( tree_in, tree_in->start[begin_in].left, target_in, idx_out );
				else
					return bstree_search( tree_in, tree_in->start[begin_in].right, target_in, idx_out );
				
			}
		}
	}
	return -1;
}

int bstree_search_key( bstree_t *tree_in, long begin_in, long key_in, long *idx_out )
{
	if( begin_in != -1 )
	{
		if( tree_in->start[begin_in].key == key_in )
		{
			*idx_out = begin_in;
			return 0; /* zero indicates that it was found */
		}
		else
		{
			if( key_in < tree_in->start[begin_in].key )
				return bstree_search_key( tree_in, tree_in->start[begin_in].left, key_in, idx_out );
			else
				return bstree_search_key( tree_in, tree_in->start[begin_in].right, key_in, idx_out );
			
		}
	}
	return -1;
}

int bstree_search_data( bstree_t *tree_in, long begin_in, void *data_in, long *idx_out )
{
	if( tree_in->compf == NULL )
		return -1; /* Fail if there is no partial order on the data field */
	if( begin_in != -1 )
	{
		if( tree_in->compf( tree_in->start[begin_in].data, data_in ) == 0 )
		{
			*idx_out = begin_in;
			return 0; /* Zero indicates that it was found */
		}
		else
		{
			if( tree_in->compf( data_in, tree_in->start[begin_in].data ) < 0 )
				return bstree_search_data( tree_in, tree_in->start[begin_in].left, data_in, idx_out );
			else
				return bstree_search_data( tree_in, tree_in->start[begin_in].right, data_in, idx_out );
		}
	}
	return -1;
}

int bstree_walk( bstree_t *tree_in, long node_in, long level_in )
{
	if( node_in != -1 )
	{
		bstree_walk( tree_in, tree_in->start[node_in].left, level_in + 1 );
		bstree_walk( tree_in, tree_in->start[node_in].right, level_in + 1 );
	}
	return 0;
}

int bstree_iterate( bstree_t *tree_in, long node_in, void (*func_in)(bsnode_t*,void*), void *args_in )
{
	if( node_in != -1 )
	{
		bstree_iterate( tree_in, tree_in->start[node_in].left, func_in, args_in );
		func_in( tree_in->start + node_in, args_in );
		bstree_iterate( tree_in, tree_in->start[node_in].right, func_in, args_in );
	}
	return 0;
}

int bstree_left_rotate( bstree_t *tree_in, long idx_in )
{
	long x,y;
	x = idx_in;
	y = tree_in->start[idx_in].right;
	if( PARENT(x) != -1 )
	{
		if( PARENT_LEFT(x) == x )
			PARENT_LEFT(x) = y;
		else
			PARENT_RIGHT(x) = y;
	}
	else
		tree_in->root = y;
	PARENT(y) = PARENT(x);
	PARENT(x) = y;
	RIGHT(x) = LEFT(y);
	if( LEFT(y) != -1 )
		LEFT_PARENT(y) = x;
	LEFT(y) = x;
}

int bstree_right_rotate( bstree_t *tree_in, long idx_in )
{
	long x,y;
	x = idx_in;
	y = tree_in->start[idx_in].left;
	if( PARENT(x) != -1 )
	{
		if( PARENT_LEFT(x) == x )
			PARENT_LEFT(x) = y;
		else
			PARENT_RIGHT(x) = y;
	}
	else
		tree_in->root = y;
	PARENT(y) = PARENT(x);
	PARENT(x) = y;
	LEFT(x) = RIGHT(y);
	if( RIGHT(y) != -1 )
		RIGHT_PARENT(y) = x;
	RIGHT(y) = x;
}

/**
 * Automatically allocate one more position at index -1
 * relative to tree_in->start so there is a place for NIL;
 * so actual pointer returned by realloc is start - 1
 * @param tree_in Tree in which to grow
 * @param inc_in Size to grow
 */
int bstree_grow( bstree_t *tree_in, long inc_in )
{
	long i;
	bsnode_t *ptr;

	if( tree_in->start != NULL )
		ptr = (bsnode_t*) realloc( tree_in->start - 1, ( tree_in->alloc + 1 + inc_in ) * sizeof(bsnode_t) );
	else
		ptr = (bsnode_t*) malloc( ( tree_in->alloc + 1 + inc_in ) * sizeof(bsnode_t) );
	if( ptr == NULL )
		return -1;
	tree_in->start = ptr + 1;
	for(i=tree_in->alloc;i<tree_in->alloc+inc_in;i++)
	{
		bsnode_init( tree_in->start + i, tree_in->step ); /* set to available and allocate data */
		tree_in->start[i].index = i;
	}
	tree_in->alloc += inc_in;
	return 0;
}

int bstree_shrink( bstree_t *tree_in, long dec_in )
{
	int i;
	int start = (tree_in->alloc-dec_in<0)?(0):(tree_in->alloc-dec_in);
	int end = tree_in->alloc;
	bsnode_t *ptr;
	for(i=start;i<end;i++)
		bsnode_free( tree_in->start + i );
	ptr = (bsnode_t*) realloc( tree_in->start - 1, ( tree_in->alloc - dec_in + 1 ) * sizeof(bsnode_t) );
	if( ptr == NULL )
		return -1;
	tree_in->start = ptr;
	tree_in->alloc -= dec_in;
	return 0;
}

void bstree_rb_check_iterate( bsnode_t *node_in, void *args_in )
{
	static long cc = 0;
	bstree_t *tree = (bstree_t*) args_in;
	if( node_in->color == RED && tree->start[node_in->parent].color == RED )
		printf( "TWO REDS IN A ROW!: p = %d : %d, color -1 = %d\n", node_in->parent, cc++, tree->start[-1].color );
}

int bstree_rb_check( bstree_t *tree_in )
{
	return bstree_iterate( tree_in, tree_in->root, &bstree_rb_check_iterate, tree_in );
}

void bstree_print( bstree_t *tree_in )
{
	bstree_iterate( tree_in, tree_in->root, &bsnode_print, NULL );
}

int tpnode_init(  )
{
	
}

int tpnode_free(  )
{
	
}

int treap_init(  )
{
	
}

int treap_free(  )
{
	
}

