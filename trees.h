#ifndef TREES_H
#define TREES_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int binary_search_int( int *, int *, int *, int *, int * ); ///< Search an integer array by bisection
int binary_search_range_int( int *, int *, int *, int *, int *, int *, int * ); ///< Search an integer array within a selected range of the array only
int binary_closest_int( int *, int *, int *, int *, int * ); ///< Search for the integer closest to the desired entry
int binary_closest_range_int( int *, int *, int *, int *, int *, int *, int * ); ///< Search for the integer closest to the desired entry within a given range only
int binary_set_intersect( int, int *, int, int *, int * );

int binary_search( long int *, long int *, long int *, long int *, long int * ); ///< The long version; see above
int binary_search_range( long int *, long int *, long int *, long int *, long int *, long int *, long int * ); ///< The long version; see above
int binary_closest( long int *, long int *, long int *, long int *, long int * ); ///< The long version; see above
int binary_closest_range( long int *, long int *, long int *, long int *, long int *, long int *, long int * ); ///< The long version; see above
int binary_set_intersect_long( long int, long int *, long int, long int *, long int * ); ///< The long version; see above

int simple_find( int, int, int * ); ///< A simple linear search function; quite slow
int simple_is_subset( int, int *, int, int * ); ///< Detect whether a set is a subset of another set
int simple_union( int, int *, int, int *, int * );
int simple_set_intersect( int, int *, int, int *, int * ); ///< Form the intersection of two sets; also slow
void bubble_sort( int, int * ); ///< A nice bubble sort algorithm for sorting integer arrays
void bubble_sort_double( int, double * ); ///< Same as above but for doubles, and using integer indexing
void complex_bubble_sort( int, double * );
void inverse_complex_bubble_sort( int, double * );
long int simple_find_long( long int, long int, long int * ); ///< The long version; see above
long int simple_is_subset_long( long int, long int *, long int, long int * ); ///< The long version; see above
long int simple_union_long( long int, long int *, long int, long int *, long int * );
long int simple_set_intersect_long( long int, long int *, long int, long int *, long int * ); ///< The long version; see above
void bubble_sort_long( long int, long int * ); ///< The long version; see above
void bubble_sort_double_long( long int, double * ); ///< The long version (indexing); see above

typedef struct
{
	void *start;	///< Pointer to first element
	void *pointer;	///< The start of allocated memory
	int step;	///< The size (in bytes) of an element
	int size;	///< The number of elements stored
	int alloc;	///< The number of elements allocated
} list_t;

typedef char byte_t;

int list_init( list_t *, int ); ///< Initialize the list for use with elements of a specified number of bytes
int list_free( list_t * ); ///< Free all data associated with the list
int list_reset( list_t * ); ///< Syncs actual data to the beginning of the allocated memory, i.e. list_t::pointer <= list_t::start
int list_shift( list_t *, int ); ///< Shifts the entire set of data the given number of element positions
int list_shift_set( list_t *, int, int, int ); ///< Shifts a certain subset of data
int list_grow( list_t *, int ); ///< Increase the allocated memory by a given amount
int list_shrink( list_t *, int ); ///< Decrease the allocated memory by a given amount
int list_prepend( list_t *, void * ); ///< Insert an element at the front of the list
int list_append( list_t *, void * ); ///< Insert an element at the end of the list
int list_popf( list_t *, void * ); ///< Get the value of the first entry in the list and remove it
int list_popb( list_t *, void * ); ///< The the value of the last entry in the list and remove it
int list_insert( list_t *, void *, int ); ///< Insert entry into an arbitrary position in a list
int list_insert_set( list_t *, void *, int, int ); ///< Insert a chunk of data
int list_delete( list_t *, int ); ///< Delete an arbitrary entry in a list
int list_delete_set( list_t *, int, int ); ///< Delete a chunk of data

typedef struct
{
	void *start;
	int step;
	long size;
	long alloc;
	int (*cmp)( void *, void * );
} heap_t;

int heap_init( heap_t *, int, int (*)(void*,void*) );
int heap_free( heap_t * );
int heap_append( heap_t *, void * );
int heap_insert( heap_t *, void * );
int heap_delete( heap_t *, long );
int heap_heapify( heap_t *, long );
int heap_build( heap_t * );
int heap_find( heap_t *, void * );
int heap_grow( heap_t *, long );
int heap_shrink( heap_t *, long );
void heap_print( heap_t * );

typedef struct
{
	void *data;
	long next;
	long last;
	int type;
	char avail; /* 0 if available for use, else not available */
} node_t;

typedef struct
{
	node_t *start; ///< Start of the nodes
	long step; ///< The size of the data structure stored
	long begin; ///< The index of the first element
	long end; ///< The index of the last element
	long size; ///< The length of the list
	long alloc; ///< The space currently allocated to start
	long *index; ///< Optionally stores the index of each node in the list from reference at llist_t::start
	int (*cmp)( void *, void * ); ///< Function used to compare data
} llist_t;

int llist_init( llist_t *, int, int (*)(void*,void*) );
int llist_free( llist_t * );
int llist_avail( llist_t *, long * ); ///< Returns an index to a free spot
int llist_prepend( llist_t *, void * );
int llist_append( llist_t *, void * );
int llist_popb( llist_t *, void * );
int llist_popf( llist_t *, void * );
int llist_insert_after( llist_t *, long, void * );
int llist_insert_before( llist_t *, long, void * );
int llist_delete( llist_t *, long );
int llist_reorder( llist_t *, llist_t * );
int llist_grow( llist_t *, long );
int llist_shrink( llist_t *, long );
int llist_find( llist_t *, void *, long * );
void llist_print( llist_t * );

typedef struct
{
	void *data;
	long key;
	char avail;
	long parent;
	long left;
	long right;
	char color;
	long index;
} bsnode_t;

int bsnode_init( bsnode_t *, long );
int bsnode_free( bsnode_t * );
int bsnode_reset( bsnode_t * );
void bsnode_print( bsnode_t *, void * );

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

int bstree_init( bstree_t *, long, int (*)(void*,void*) );
int bstree_free( bstree_t * );
int bstree_avail( bstree_t *, long * );
int bstree_insert( bstree_t *, void *, long, long * );
int bstree_balance_insert( bstree_t *, long );
int bstree_insert_balanced( bstree_t *, void *, long, long * );
int bstree_successor( bstree_t *, long, long * );
int bstree_predecessor( bstree_t *, long, long * );
int bstree_minimum( bstree_t *, long, long * );
int bstree_maximum( bstree_t *, long, long * );
int bstree_delete_key( bstree_t *, target_t, long * );
int bstree_delete_index( bstree_t *, long, long * );
int bstree_balance_delete( bstree_t *, long );
int bstree_delete_key_balanced( bstree_t *, target_t, long * );
int bstree_delete_index_balanced( bstree_t *, long );
int bstree_pop_front( bstree_t *, void *, long * );
int bstree_pop_front_balanced( bstree_t *, void *, long * );
int bstree_pop_back( bstree_t *, void *, long * );
int bstree_pop_back_balanced( bstree_t *, void *, long * );
int bstree_search( bstree_t *, long, target_t, long * );
int bstree_search_key( bstree_t *, long, long, long * );
int bstree_search_data( bstree_t *, long, void *, long * );
int bstree_walk( bstree_t *, long, long );
int bstree_iterate( bstree_t *, long, void (*)(bsnode_t*,void*), void * );
int bstree_left_rotate( bstree_t *, long );
int bstree_right_rotate( bstree_t *, long );
int bstree_grow( bstree_t *, long );
int bstree_shrink( bstree_t *, long );
int bstree_rb_check( bstree_t * );
void bstree_print( bstree_t * );
void bstree_info( bstree_t * );

typedef struct
{
	void *data;
	long key;
	byte_t avail;
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

int arraynext( long, long *, long * );

long arrayindex( long, long *, long * );

#endif
