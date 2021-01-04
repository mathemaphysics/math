#ifndef COMBINADIC_H
#define COMBINADIC_H

/* Typical versions */
int factorial( int ); ///< The standard factorial function
long int factorial_long( long int ); ///< The long (and more useful) version
int binomial_slow( int n_in, int k_in ); ///< A slower, simpler binomial coefficient function
long int binomial_slow_long( long int n_in, long int k_in ); ///< The long version; see above
int binomial( int, int ); ///< The binomial function which uses truncated factorials
long int binomial_long( long int, long int ); ///< The long version; later suffix with _long
int truncfact( int, int ); ///< The truncated factorial function
long int truncfact_long( long int, long int ); ///< The long version; later suffix with _long (pretty please?)

/* Combinadics */
int combinadic_index( int, int * ); ///< Calculate the combinadic index of an ordered set of integers
int combinadic_init( int, int, int * ); ///< Initialize the combinadic vector to the first combination
int combinadic_next( int, int, int * ); ///< Step to the next combination
int combinadic_vector( int, int, int, int * ); ///< Calculate the combinadic vector given the combinadic index
int combinadic_vector_length( int, int, int ); ///< Calculate the dimension of the combinadic with given index

/* Long combinadics */
long int combinadic_index_long( long int, long int * ); ///< The long version; later suffix with _long
long int combinadic_init_long( long int, long int, long int * ); ///< The long version; later suffix with _long
long int combinadic_next_long( long int, long int, long int * ); ///< The long version; later suffix with _long
long int combinadic_vector_long( long int, long int, long int, long int * ); ///< The long version; later suffix with _long
long int combinadic_vector_length_long( long int, long int, long int ); ///< The long version; later suffix with _long

/* Factoradics */
int factoradic_radix_index( int, int, int * ); ///< Calculate the index of a factoradic vector in radix form
int factoradic_radix_vector( int, int *, int * ); ///< Calculate the radix vector given its factoradic index
int factoradic_index( int, int * ); ///< Calculate the factoradic index from the actual permutation
int factoradic_init( int, int * ); ///< Initialize the factoradic vector
int factoradic_next( int, int * ); ///< Step one permutation forward; not yet written
int factoradic_vector( int, int, int * ); ///< Calculate the factoradic permutation vector from its index

/* Long factoradics */
long int factoradic_radix_index_long( long int, long int, long int * ); ///< The long version; see above
long int factoradic_radix_vector_long( long int, long int *, long int * ); ///< The long version; see above
long int factoradic_index_long( long int, long int * ); ///< The long version; see above
long int factoradic_init_long( long int, long int * ); ///< The long version; see above
long int factoradic_next_long( long int, long int * ); ///< The long version; see above
long int factoradic_vector_long( long int, long int, long int * ); ///< The long version; see above

/* Combinadic with repeat & distinguishable containers with indistinguishable objects */
int rcombinadic_index( int, int * );
int rcombinadic_init( int, int, int * );
int rcombinadic_next( int, int, int * );
int rcombinadic_vector( int, int, int, int * );
int rcombinadic_occupancy( int, int, int *, int * );
int rcombinadic_invoccup( int, int *, int * );

/* Long versions */
long int rcombinadic_index_long( long int, long int * );
long int rcombinadic_init_long( long int, long int, long int * );
long int rcombinadic_next_long( long int, long int, long int * );
long int rcombinadic_vector_long( long int, long int, long int, long int * );
long int rcombinadic_occupancy_long( long int, long int, long int *, long int * );
long int rcombinadic_invoccup_long( long int, long int *, long int * );

int polynomial_index( int, int * );
int polynomial_exponents( int, int, int, int * );
int global_poly_index( int, int, int * );

long int polynomial_index_long( long int, long int * );
long int polynomial_exponents_long( long int, long int, long int, long int * );
long int global_poly_index_long( long int, long int, long int * );

int global_polynomial_vector( int, int, int * );
long int global_polynomial_vector_long( long int, long int, long int * );

/* Stirling numbers of the second kind */
int stirling_index( int, int * ); ///< Calculate the stirling partition index from the partition vector
int stirling_vector( int, int, int * ); ///< Calculate the stirling partition vector from the index

/* Long Stirling numbers */
long int stirling_index_long( long int, long int * ); ///< The long version; see above
long int stirling_vector_long( long int, long int, long int * ); ///< The long version; see above

/* Partitioning functions */
void partition_init( int *, int *, int );
int partition_next( int *, int *, int );

#endif
