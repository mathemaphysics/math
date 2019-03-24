#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define CGRAD_ZERO 1e-15
#define CGRAD_TREIG_PREPROCESS /* Function trofrm will fail if no preprocessing is done */

/**
 * Sparse tridiagonal save in octave format
 */
int stmsave( char *fn, int n, double *alpha, double *beta, double *gamma )
{
    int i;
    FILE *fp = fopen( fn, "w" );
    fprintf( fp, "# name: T\n" );
    fprintf( fp, "# type: sparse matrix\n" );
    fprintf( fp, "# rows: %d\n", n );
    fprintf( fp, "# columns: %d\n", n );
    for(i=0;i<n;i++)
    {
        if( i > 0 )
            fprintf( fp, "%d %d %18.10f\n", i, i + 1, beta[i] );
        fprintf( fp, "%d %d %18.10f\n", i + 1, i + 1, alpha[i] );
        if( i < n - 1 )
            fprintf( fp, "%d %d %18.10f", i + 2, i + 1, gamma[i+1] );
        fprintf( fp, "\n" );
    }
    fclose( fp );
    return 0;
}

/**
 * Save sparse matrices in a readable format
 * @param fn File name to which output is directed
 * @param n Number of rows in matrix a
 * @param k Number of columns in matrix a
 * @param ia Sparse matrix parameters; list of beginning of each row
 * @param ja List of column numbers corresponding to entries in a
 * @param a Actual values corresponding to ja indexes
 * @param mode Mode of output; 0 means output coordinates, 1 means output in octave/matlab format
 * @return Returns -1 if error in opening the output file, 0 otherwise
 */
int smsave( char *fn, long n, long k, long *ia, long *ja, double *a, int mode )
{
    long i,j;
    static int idx = 1;
    char vn[512];
    FILE *fp = fopen( fn, "w" );
    if( fp == NULL )
        return -1;

    switch( mode )
    {
        case 0:
            fprintf( fp, "%d %d\n", n, ia[n] );
            break;
        case 1:
            snprintf( vn, 512, "A%d", idx );
            fprintf( fp, "# name: %s\n# type: sparse matrix\n# nnz: %d\n# rows: %d\n# columns: %d\n", vn, ia[n], n, k );
            ++idx;
            break;
    }
    for(i=0;i<n;i++)
    {
        for(j=ia[i];j<ia[i+1];j++)
        {
            switch( mode )
            {
                case 0:
                    fprintf( fp, "%d %d %17.12f\n", i, ja[j], a[j] );
                    break;
                case 1:
                    fprintf( fp, "%d %d %17.12f\n", ja[j] + 1, i + 1, a[j] );
                    break;
            }
        }
    }

    fclose( fp );
    return 0;
}

/**
 * Save a dense matrix into a file
 */
int msave( char *fn, long n, long k, double *a, int mode )
{
    long i,j;
    static int idx = 1;
    char vn[512];
    FILE *fp = fopen( fn, "w" );
    if( fp == NULL )
        return -1;

    switch( mode )
    {
        case 0:
            fprintf( fp, "%d %d\n", n, k );
            break;
        case 1:
            snprintf( vn, 512, "M%d", idx );
            fprintf( fp, "# name: %s\n# type: matrix\n# rows: %d\n# columns: %d\n", vn, n, k );
            ++idx;
            break;
    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<k;j++)
            fprintf( fp, "%17.12f", a[i*k+j] );
        fprintf( fp, "\n" );
    }

    fclose( fp );
    return 0;
}

/**
 * Print a vector to a stream
 * @param fp File stream pointer to which to direct output
 * @param n The number of entries in the vector x
 * @param x Pointer to the start of the entries of the vector
 */
void vprint( FILE *fp, int n, double *x )
{
    int i;
    for(i=0;i<n;i++)
        fprintf( fp, "%20.13f", x[i] );
    fprintf( fp, "\n" );
}

/**
 * Return 1 if x is within CGRAD_ZERO of zero
 * @param x Double precision floating point to check for zero
 * @return Returns 1 if zero, 0 otherwise
 */
int iszero( double x )
{
    if( fabs( x ) < CGRAD_ZERO )
        return 1;
    else
        return 0;
}

/**
 * Vector copy
 * @param n Length of the vector
 * @param u Origin vector
 * @param v Destination to which to copy u
 */
void copy( int n, double *u, double *v )
{
    int i;

    for(i=0;i<n;i++)
        v[i] = u[i];
}

/**
 * Vector copy and divide by constant
 * @param n Length of vector
 * @param u Origin vector to copy
 * @param c Constant by which to divide u
 * @param v Destination for u / c
 */
void copyvdiv( int n, double *u, double c, double *v )
{
    int i;

    for(i=0;i<n;i++)
        v[i] = u[i] / c;
}

/**
 * Divide every element of a vector for a constant
 * @param n Vector length
 * @param x Vector to divide
 * @param c Constant by which to divide x; result is in x as well
 */
void vdiv( int n, double *x, double c )
{
    int i;

    for(i=0;i<n;i++)
        x[i] /= c;
}

/**
 * Zero a double precision vector
 * @param n Length of vector
 * @param x Pointer to entries to zero
 */
void zerov( int n, double *x )
{
    int i;

    for(i=0;i<n;i++)
        x[i] = 0.0;
}

/**
 * Is norm of x less than or greater than input tol
 * @param n Length of vector
 * @param x Vector to check against tol
 * @param tol Tolerance to use
 * @param res Returns -1 if norm(x) < tol else 1
 */
void normchk( int n, double *x, double tol, int *res )
{
    int i;
    double sum = 0.0;
    const double tolsq = tol * tol;

    for(i=0;i<n;i++)
        sum += x[i] * x[i];
    if( sum < tolsq )
        *res = -1;
    else
        *res = 1;
}

/**
 * Is the dot product of x and y less than or greater than tol
 * @param n Dimension of vectors
 * @param x First vector
 * @param y Second vector
 * @param tol Tolerance with which to compare x dot y
 * @param res Returns -1 if x dot y < tol else 1
 */
void dotchk( int n, double *x, double *y, double tol, int *res )
{
    int i;
    double sum = 0.0;

    for(i=0;i<n;i++)
        sum += x[i] * y[i];
    if( fabs( sum ) < tol )
        *res = -1;
    else
        *res = 1;
}

/**
 * Calculate and return the norm of the vector
 * @param n Dimension of vector x
 * @param x Vector whose norm is calculated
 * @param res Value of the norm of x
 */
void norm( int n, double *x, double *res )
{
    int i;

    (*res) = 0.0;
    for(i=0;i<n;i++)
        (*res) += x[i] * x[i];
    (*res) = sqrt( *res );
}

/**
 * Calculate and return the dot product
 * @param n Dimension of vectors x and y
 * @param x First vector
 * @param y Second vector
 * @param res Returns x dot y
 */
void dotp( int n, double *x, double *y, double *res )
{
    int i;

    (*res) = 0.0;
    for(i=0;i<n;i++)
        (*res) += x[i] * y[i];
}

/**
 * Calculate and return the dot product divided by c
 * @param n Dimension of x and y
 * @param x First vector
 * @param y Second vector
 * @param c Scalar by which to divide x dot y
 * @param res Contains x dot y / c
 */
void dotpdiv( int n, double *x, double *y, double c, double *res )
{
    int i;

    (*res) = 0.0;
    for(i=0;i<n;i++)
        (*res) += x[i] * y[i];
    (*res) = *res / c;
}

/**
 * Normalize a double vector
 * @param n Dimension of x
 * @param x Vector to normalize
 */
void normalize( int n, double *x )
{
    int i;
    double sum = 0.0;

    for(i=0;i<n;i++)
        sum += x[i] * x[i];
    sum = sqrt( sum );
    for(i=0;i<n;i++)
        x[i] /= sum;
}

/**
 * Generate a random double precision n-vector normalized to unity;
 * make sure to call srand() before running this
 * @param n Vector length
 * @param x Vector
 */
void nrandv( int n, double *x )
{
    int i;

    for(i=0;i<n;i++)
        x[i] = 0.5 - ( (double) rand() / (double) RAND_MAX );
    normalize( n, x );
}

/**
 * Project a subspace out of the vector
 * @param n System dimension
 * @param x Vector from which to project vector space in V
 * @param nv Number of vectors in V to project out of x
 * @param V Vector space to project out of x
 */
void project( int n, double *x, int nv, double *V )
{
    int i,j;
    double f;

    for(i=0;i<nv;i++)
    {
        f = 0.0;
        for(j=0;j<n;j++)
            f += x[j] * V[i*n+j];
        for(j=0;j<n;j++)
            x[j] -= f * V[i*n+j];
    }
}

void bicgstab( int, int, double *, double *, double *, int, double, int, int * );

/**
 * Project a vector into a subspace via least-squares calculation
 */
void projectls( int n, double *x, int nv, double *V, int max )
{
    int i,j,k,res;
    double *C = (double*) malloc( n * n * sizeof(double) );
    double *b = (double*) malloc( nv * sizeof(double) );
    double *y = (double*) malloc( nv * sizeof(double) );

    /* Form the V**T V and put it in C */
    for(i=0;i<nv;i++)
    {
        for(j=0;j<nv;j++)
        {
            C[i*n+j] = 0.0;
            for(k=0;k<n;k++)
                C[i*n+j] += V[i*n+k] * V[j*n+k];
        }
    }
    for(i=0;i<nv;i++)
    {
        b[i] = 0.0;
        for(j=0;j<n;j++)
            b[i] += V[i*n+j] * x[j];
    }

    /* Solve the least squares problem in full rank form */
    nrandv( nv, y );
    bicgstab( nv, 0, C, b, y, max, 1e-9, 0, &res );

    /* Now calculate the actual projection as x - Cy */
    for(i=0;i<nv;i++)
        for(j=0;j<n;j++)
            x[j] -= y[i] * V[i*n+j];
}

/**
 * Check to make sure vector x contains no component in
 * the given vector subspace
 */
void spcheck( int n, double *x, int nv, double *V, double tol, int *res )
{
    int i,j;
    double f;

    *res = 0;
    for(i=0;i<nv;i++)
    {
        f = 0.0;
        for(j=0;j<n;j++)
            f += x[j] * V[i*n+j];
        if( fabs( f ) > tol )
            *res += 1;
    }
}

/**
 * Vector sum combination; calculate res = a*x + b*y
 * @param n Dimension of vectors x and y
 * @param a First constant
 * @param x First vector
 * @param b Second constant
 * @param y Second vector
 * @param res Returns a*x + b*y
 */
void vsum( int n, double a, double *x, double b, double *y, double *res )
{
    int i;

    for(i=0;i<n;i++)
        res[i] = a * x[i] + b * y[i];
}

/**
 * Matrix vector product
 */
void dgemv( int n, int tt, double *A, double *x, double *res )
{
    int i,j;

    for(i=0;i<n;i++)
        res[i] = 0.0;
    if( tt == 0 )
        for(i=0;i<n;i++)
            for(j=0;j<n;j++)
                res[i] += A[i*n+j] * x[j];
    else
        for(i=0;i<n;i++)
            for(j=0;j<n;j++)
                res[i] += A[j*n+i] * x[j];
}

/**
 *  A very nice and VERY simple algorithm!
 */
void sdgemv( long n, long m, long *ia, long *ja, double *a, long rb, double *b, long rc, double *c )
{
    int i,j;
    for(i=0;i<n;i++)
    {
        c[i*rc] = 0.0;
        for(j=ia[i];j<ia[i+1];j++)
            c[i*rc] = c[i*rc] + a[j] * b[ja[j]*rb];
    }
}

/**
 * Sparse symmetric dgemv; only lower triangular part is stored in
 * the sparse data structure
 */
void ssdgemv( long n, long *ia, long *ja, double *A, long rx, double *x, long ry, double *y )
{
    int i,j;
    double stemp,rtemp;

    for(i=0;i<n;i++)
        y[i*ry] = 0.0;
    for(i=0;i<n;i++)
    {
        j = ia[i];
        rtemp = x[i*rx];
        stemp = 0.0;
        if( ja[j] == i )
        {
            stemp = A[j] * x[i*rx];
            j++;
        }
        for(;j<ia[i+1];j++)
        {
            stemp += A[j] * x[ja[j]*rx];
            y[ja[j]*ry] += A[j] * rtemp;
        }
        y[i*ry] += stemp;
    }
}

/**
 * Sparse matrix transpose
 */
void stransp( char move, long n, long m, long *ia, long *ja, double *a, long *ib, long *jb, double *b )
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

/**
 * Compute an inner product w.r.t. a given operator
 * @param n Dimension of space
 * @param A Matrix defining the inner product
 * @param x Left vector
 * @param y Right vector
 * @param res Inner product value
 */
void mdotp( int n, double *A, double *x, double *y, double *res )
{
    int i,j;

    (*res) = 0.0;
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            (*res) += x[i] * A[i*n+j] * y[j];
}

/**
 * Calculate the norm of the difference between two vectors
 */
void vnormdiff( int n, double *x, double *y, double *res )
{
    int i;
    double sum = 0.0;

    for(i=0;i<n;i++)
        sum += pow( x[i] - y[i], 2.0 );
    *res = sqrt( sum );
}

/**
 * Calculate the residual of a linear system
 */
void residual( int n, double *A, double *x, double *b, double *r )
{
    int i,j;

    for(i=0;i<n;i++)
        r[i] = b[i];
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            r[i] -= A[i*n+j] * x[j];
}

/**
 * Calculate shifted-matrix-vector product with an addition;
 * res = (A-sI)**(t) x + c y
 * @param n Dimension of space
 * @param t No transpose if t = 0, otherwise transpose A
 * @param A Matrix input
 * @param s Value to shift the diagonal of A
 * @param x Vector input one
 * @param y Vector input two
 * @param c Multiplier parameter on y
 * @param res Output vector
 */
void smmadd( int n, int tt, double *A, double s, double *x, double *y, double c, double *res )
{
    int i,j;

    for(i=0;i<n;i++)
        res[i] = c * y[i];
    if( tt == 0 ) /* Not transposed */
    {
        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            {
                if( j == i )
                    res[i] += ( A[i*n+j] - s ) * x[j];
                else
                    res[i] += A[i*n+j] * x[j];
            }
        }
    }
    else /* If not zero assume transposed */
    {
        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            {
                if( j == i )
                    res[i] += ( A[j*n+i] - s ) * x[j];
                else
                    res[i] += A[j*n+i] * x[j];
            }
        }
    }
}

/**
 * Calculate shifted-matrix-vector product with an addition;
 * res = (A-sB)**(t) x + c y
 * @param n Dimension of space
 * @param t No transpose if t = 0, otherwise transpose A
 * @param A Matrix 1 input
 * @param B Matrix 2 input
 * @param s Value to shift the diagonal of A
 * @param x Vector input one
 * @param y Vector input two
 * @param c Multiplier parameter on y
 * @param res Output vector
 */
void gsmmadd( int n, int tt, double *A, double *B, double s, double *x, double *y, double c, double *res )
{
    int i,j;

    for(i=0;i<n;i++)
        res[i] = c * y[i];
    if( tt == 0 ) /* Not transposed */
    {
        for(i=0;i<n;i++)
            for(j=0;j<n;j++)
                res[i] += ( A[i*n+j] - s * B[i*n+j] ) * x[j];
    }
    else /* If not zero assume transposed */
    {
        for(i=0;i<n;i++)
            for(j=0;j<n;j++)
                res[i] += ( A[j*n+i] - s * B[j*n+i] ) * x[j];
    }
}

void spdiag( long n, long *ia, long *ja, double *A, long *d, int *ret )
{
    long i,j;

    /* Find indexes of the diagonal elements */
    for(i=0;i<n;i++)
    {
        d[i] = -1;
        for(j=ia[i];j<ia[i+1];j++)
        {
            if( ja[j] == i )
            {
                d[i] = j;
                break;
            }
        }
        if( d[i] == -1 )
        {
            *ret = -1; /* Missing a diagonal element */
            break;
        }
    }
}

/**
 * Apply a sparse lower triangular inverse preconditioner to a given input vector
 * and overwrite it with the result
 * @param n Dimension of the vector and sparse square operator
 * @param ia List of row start positions
 * @param ja List of column indexes
 * @param A Entries in matrix A
 * @param x Vector to which to apply lower triangular preconditioner
 * @param d Positions of the diagonal elements in each row
 */
void spclw( long n, long *ia, long *ja, double *A, double *x, long *d )
{
    long i,j;

    /* Do the forward substitution */
    for(i=0;i<n;i++)
    {
        for(j=ia[i];j<d[i];j++)
            x[i] = x[i] - A[j] * x[ja[j]];
        x[i] = x[i] / A[d[i]];
    }
}

/**
 * Apply a sparse upper triangular inverse preconditioner to a given input vector
 * and overwrite it with the result
 * @param n Dimension of the vector and sparse square operator
 * @param ia List of row start positions
 * @param ja List of column indexes
 * @param A Entries in matrix A
 * @param x Vector to which to apply lower triangular preconditioner
 * @param d Positions of the diagonal elements in each row
 */
void spcup( long n, long *ia, long *ja, double *A, double *x, long *d )
{
    long i,j;

    /* Do the back substitution */
    for(i=n-1;i>=0;i--)
    {
        for(j=ia[i+1]-1;j>d[i];j--)
            x[i] = x[i] - A[j] * x[ja[j]];
        x[i] = x[i] / A[d[i]];
    }
}

/**
 * Regular sparse symmetric conjugate gradient algorithm; A must be
 * in CRS format with only upper triangular and diagonal entries
 * @param n Dimension of sparse square matrix
 * @param ia Row beginning indexes
 * @param ja Column index list
 * @param A Entries of the matrix
 * @param b Righthand side
 * @param x Input/output solution
 * @param tol Tolerance to indicate when to converge
 * @param vb Verbose or not
 * @param ret Returns 0 if all went well, < 0 otherwise
 */
void scg( int n, long *ia, long *ja, double *A, double *b, double *x, int max, double tol, int vb, int *ret )
{
    int i,j;
    double alpha,tmp,rso,rsn,*r,*p,*t;

    /* Start out okay and change if necessary */
    *ret = 0;

    /* Allocate stuff */
    r = (double*) malloc( n * sizeof(double) );
    p = (double*) malloc( n * sizeof(double) );
    t = (double*) malloc( n * sizeof(double) );

    /* Set up initial vectors and coefficients */
    ssdgemv( (long) n, ia, ja, A, 1, x, 1, r );
    for(i=0;i<n;i++)
        r[i] = b[i] - r[i], p[i] = r[i];
    dotp( n, r, r, &rso );

    /* Start iterating */
    for(i=0;i<max;i++)
    {
        ssdgemv( (long) n, ia, ja, A, 1, p, 1, t );
        dotp( n, p, t, &tmp );
        alpha = rso / tmp;
        for(j=0;j<n;j++)
            x[j] = x[j] + alpha * p[j],
                r[j] = r[j] - alpha * t[j];
        dotp( n, r, r, &rsn );
        if( rsn < tol * tol )
            break;
        for(j=0;j<n;j++)
            p[j] = r[j] + rsn / rso * p[j];
        rso = rsn;
    }
    if( vb )
        fprintf( stderr, "%d: residual = %15.7f\n", i, sqrt( rsn ) );

    /* Clean up */
    free( r ); free( p ); free( t );
}

/**
 * Regular sparse symmetric conjugate gradient algorithm; A must be
 * in CRS format; Jacobi preconditioner applied
 * @param n Dimension of sparse square matrix
 * @param ia Row beginning indexes
 * @param ja Column index list
 * @param A Entries of the matrix
 * @param b Righthand side
 * @param x Input/output solution
 * @param tol Tolerance to indicate when to converge
 * @param vb Verbose or not
 * @param ret Returns 0 if all went well, < 0 otherwise
 */
void jpscg( int n, long *ia, long *ja, double *A, double *bb, double *x, int max, double tol, int vb, int *ret )
{
    int i,j,res;
    double alpha,tmp,rso,rsn,*r,*p,*t,*b,*pc;

    /* Start out okay and change if necessary */
    *ret = 0;

    /* Allocate stuff */
    r = (double*) malloc( n * sizeof(double) );
    p = (double*) malloc( n * sizeof(double) );
    t = (double*) malloc( n * sizeof(double) );
    b = (double*) malloc( n * sizeof(double) );
    pc = (double*) malloc( n * sizeof(double) );

    /* Create the Jacobi preconditioner */
    for(i=0;i<n;i++)
        for(j=ia[i];j<ia[i+1];j++)
            if( j == ia[i] || fabs( A[j] ) > fabs( pc[i] ) )
                pc[i] = fabs( A[j] );

    /* Set up the preconditioned system */
    copy( n, bb, b );
    for(i=0;i<n;i++)
        b[i] = b[i] / pc[i];

    /* Set up initial vectors and coefficients */
    ssdgemv( (long) n, ia, ja, A, 1, x, 1, r );
    for(i=0;i<n;i++)
        r[i] = r[i] / pc[i];
    for(i=0;i<n;i++)
        r[i] = b[i] - r[i], p[i] = r[i];
    dotp( n, r, r, &rso );

    /* Start iterating */
    for(i=0;i<max;i++)
    {
        ssdgemv( (long) n, ia, ja, A, 1, p, 1, t );
        for(j=0;j<n;j++)
            t[j] = t[j] / pc[j];
        dotp( n, p, t, &tmp );
        alpha = rso / tmp;
        for(j=0;j<n;j++)
            x[j] = x[j] + alpha * p[j],
                r[j] = r[j] - alpha * t[j];
        dotp( n, r, r, &rsn );
        if( rsn < tol * tol )
            break;
        for(j=0;j<n;j++)
            p[j] = r[j] + rsn / rso * p[j];
        rso = rsn;
        if( vb )
            fprintf( stderr, "%d: residual = %15.7f\n", i, sqrt( rsn ) );
    }

    /* Clean up */
    free( r ); free( p ); free( t ); free( b ); free( pc );
}

/**
 * Biconjugate gradient algorithm
 */
void bicg( int n, double *A, double *b, double *x, int max, double tol, int vb, int *ret )
{

}

/**
 * Preconditioned biconjugate gradient algorithm
 */
void pbicg( int n, double *A, double *b, double *M, double *x, int max, double tol, int vb, int *ret )
{

}

/**
 * Biconjugate gradient method projected
 */
void bicgp( int n, double *A, double *b, int nv, double *V, double *x, int max, double tol, int vb, int *ret )
{
    int i,j,m;
    double f,g,rho,rhon,alpha,beta,omega;
    double *xh,*r,*rh,*p,*ph;
}

/**
 * Preconditioned biconjugate gradient projected method
 */
void pbicgp( int n, double *A, double *b, double *M, int nv, double *V, double *x, int max, double tol, int vb, int *ret )
{

}

/**
 * Standard biconjugate gradient stabilized algorithm
 */
void bicgstab( int n, int tt, double *A, double *b, double *x, int max, double tol, int vb, int *ret )
{
    int i,j,m;
    double f,g,rho,rhon,alpha,beta,omega;
    double *r0,*r,*p,*v,*s,*t;

    r0 = (double*) malloc( n * sizeof(double) );
    r = (double*) malloc( n * sizeof(double) );
    p = (double*) malloc( n * sizeof(double) );
    v = (double*) malloc( n * sizeof(double) );
    s = (double*) malloc( n * sizeof(double) );
    t = (double*) malloc( n * sizeof(double) );

    /* Initialize */
    rho = 1.0;
    alpha = 1.0;
    omega = 1.0;
    for(i=0;i<n;i++)
    {
        r[i] = b[i];
        if( tt == 0 )
            for(j=0;j<n;j++)
                r[i] -= A[i*n+j] * x[j];
        else
            for(j=0;j<n;j++)
                r[i] -= A[j*n+i] * x[j];
        r0[i] = r[i];
        v[i] = 0.0;
        p[i] = 0.0;
    }

    /* Iterate */
    if( vb )
    {
        fprintf( stderr, "\n" );
        fprintf( stderr, "----------------------------------------------------------------------\n" );
        fprintf( stderr, "%10s%15s%15s%15s%15s\n", " Iteration", "Residual ", "Alpha  ", "Beta  ", "Omega  " );
        fprintf( stderr, "----------------------------------------------------------------------\n" );
    }
    *ret = 0;
    for(m=0;m<max;m++)
    {
        rhon = 0.0;
        for(i=0;i<n;i++)
            rhon += r0[i] * r[i];
        beta = ( rhon / rho ) * ( alpha / omega );
        rho = rhon;
        for(i=0;i<n;i++)
            p[i] = r[i] + beta * ( p[i] - omega * v[i] );
        for(i=0;i<n;i++)
        {
            v[i] = 0.0;
            if( tt == 0 )
                for(j=0;j<n;j++)
                    v[i] += A[i*n+j] * p[j];
            else
                for(j=0;j<n;j++)
                    v[i] += A[j*n+i] * p[j];
        }
        f = 0.0;
        for(i=0;i<n;i++)
            f += r0[i] * v[i];
        alpha = rho / f;
        if( isnan( alpha ) )
        {
            *ret = -1; /* Failed to converge */
            break;
        }
        for(i=0;i<n;i++)
            s[i] = r[i] - alpha * v[i];
        for(i=0;i<n;i++)
        {
            t[i] = 0.0;
            if( tt == 0 )
                for(j=0;j<n;j++)
                    t[i] += A[i*n+j] * s[j];
            else
                for(j=0;j<n;j++)
                    t[i] += A[j*n+i] * s[j];
        }
        f = 0.0, g = 0.0;
        for(i=0;i<n;i++)
            f += t[i] * s[i],
              g += t[i] * t[i];
        omega = f / g;
        for(i=0;i<n;i++)
            x[i] += alpha * p[i] + omega * s[i];

        /* Check for convergence here */
        f = 0.0;
        for(i=0;i<n;i++)
            f += r[i] * r[i];
        f = sqrt( f );
        if( vb )
            fprintf( stderr, "%10d%15.7f%15.7f%15.7f%15.7f\n", m, f, alpha, beta, omega );
        if( f < tol )
            break;

        for(i=0;i<n;i++)
            r[i] = s[i] - omega * t[i];
    }
    if( f > tol )
        *ret = -2;

    if( vb )
    {
        fprintf( stderr, "----------------------------------------------------------------------\n" );
        if( *ret == 0 )
            fprintf( stderr, " Converged                              System size:   %15d\n", n );
        else
            fprintf( stderr, " Failed                                 System size:   %15d\n", n );
        fprintf( stderr, "                                        Tolerance:     %15.3e\n", tol );
        fprintf( stderr, "                                        Iterations:    %15d\n", m );
        fprintf( stderr, "                                        Residual:      %15.3e\n", f );
        fprintf( stderr, "                                        ------------------------------\n" );
    }

    free( r0 );
    free( r );
    free( p );
    free( v );
    free( s );
    free( t );
}

/**
 * Standard biconjugate gradient stabilized algorithm using sparse matrices
 */
void sbicgstab( int n, long *ia, long *ja, double *A, double *b, double *x, int max, double tol, int vb, int *ret )
{
    int i,j,m;
    double f,g,rho,rhon,alpha,beta,omega;
    double *r0,*r,*p,*v,*s,*t;

    r0 = (double*) malloc( n * sizeof(double) );
    r = (double*) malloc( n * sizeof(double) );
    p = (double*) malloc( n * sizeof(double) );
    v = (double*) malloc( n * sizeof(double) );
    s = (double*) malloc( n * sizeof(double) );
    t = (double*) malloc( n * sizeof(double) );

    /* Initialize */
    rho = 1.0;
    alpha = 1.0;
    omega = 1.0;
    sdgemv( (long) n, (long) n, ia, ja, A, 1, x, 1, r );
    for(i=0;i<n;i++)
        r[i] = b[i] - r[i];
    for(i=0;i<n;i++)
    {
        r0[i] = r[i];
        v[i] = 0.0;
        p[i] = 0.0;
    }

    /* Iterate */
    if( vb )
    {
        fprintf( stderr, "\n" );
        fprintf( stderr, "----------------------------------------------------------------------\n" );
        fprintf( stderr, "%10s%15s%15s%15s%15s\n", " Iteration", "Residual ", "Alpha  ", "Beta  ", "Omega  " );
        fprintf( stderr, "----------------------------------------------------------------------\n" );
    }
    *ret = 0;
    for(m=0;m<max;m++)
    {
        rhon = 0.0;
        for(i=0;i<n;i++)
            rhon += r0[i] * r[i];
        beta = ( rhon / rho ) * ( alpha / omega );
        rho = rhon;
        for(i=0;i<n;i++)
            p[i] = r[i] + beta * ( p[i] - omega * v[i] );
        sdgemv( (long) n, (long) n, ia, ja, A, 1, p, 1, v );

        f = 0.0;
        for(i=0;i<n;i++)
            f += r0[i] * v[i];
        alpha = rho / f;
        if( isnan( alpha ) )
        {
            *ret = -1; /* Failed to converge */
            break;
        }
        for(i=0;i<n;i++)
            s[i] = r[i] - alpha * v[i];
        sdgemv( (long) n, (long) n, ia, ja, A, 1, s, 1, t );

        f = 0.0, g = 0.0;
        for(i=0;i<n;i++)
            f += t[i] * s[i],
              g += t[i] * t[i];
        omega = f / g;
        for(i=0;i<n;i++)
            x[i] += alpha * p[i] + omega * s[i];

        /* Check for convergence here */
        f = 0.0;
        for(i=0;i<n;i++)
            f += r[i] * r[i];
        f = sqrt( f );
        if( vb )
            fprintf( stderr, "%10d%15.7f%15.7f%15.7f%15.7f\n", m, f, alpha, beta, omega );
        if( f < tol )
            break;

        for(i=0;i<n;i++)
            r[i] = s[i] - omega * t[i];
    }
    if( f > tol )
        *ret = -2;

    if( vb )
    {
        fprintf( stderr, "----------------------------------------------------------------------\n" );
        if( *ret == 0 )
            fprintf( stderr, " Converged                              System size:   %15d\n", n );
        else
            fprintf( stderr, " Failed                                 System size:   %15d\n", n );
        fprintf( stderr, "                                        Tolerance:     %15.3e\n", tol );
        fprintf( stderr, "                                        Iterations:    %15d\n", m );
        fprintf( stderr, "                                        Residual:      %15.3e\n", f );
        fprintf( stderr, "                                        ------------------------------\n" );
    }

    free( r0 );
    free( r );
    free( p );
    free( v );
    free( s );
    free( t );
}

/**
 * Lower triangular preconditioned biconjugate gradient stabilized algorithm using sparse matrices
 */
void lpsbicgstab( int n, long *ia, long *ja, double *A, double *bb, double *x, int max, double tol, int vb, int *ret )
{
    int i,j,m,res;
    double f,g,rho,rhon,alpha,beta,omega;
    double *r0,*r,*p,*v,*s,*t,*b;
    long *d;

    r0 = (double*) malloc( n * sizeof(double) );
    r = (double*) malloc( n * sizeof(double) );
    p = (double*) malloc( n * sizeof(double) );
    v = (double*) malloc( n * sizeof(double) );
    s = (double*) malloc( n * sizeof(double) );
    t = (double*) malloc( n * sizeof(double) );
    b = (double*) malloc( n * sizeof(double) );
    d = (long*) malloc( n * sizeof(long) );

    /* Apply the preconditioner to the righthand side */
    copy( n, bb, b );
    spdiag( (long) n, ia, ja, A, d, &res );
    spclw( (long) n, ia, ja, A, b, d );

    /* Initialize */
    rho = 1.0;
    alpha = 1.0;
    omega = 1.0;
    sdgemv( (long) n, (long) n, ia, ja, A, 1, x, 1, r );
    spclw( n, ia, ja, A, r, d );
    for(i=0;i<n;i++)
        r[i] = b[i] - r[i];
    for(i=0;i<n;i++)
    {
        r0[i] = r[i];
        v[i] = 0.0;
        p[i] = 0.0;
    }

    /* Iterate */
    if( vb )
    {
        fprintf( stderr, "\n" );
        fprintf( stderr, "----------------------------------------------------------------------\n" );
        fprintf( stderr, "%10s%15s%15s%15s%15s\n", " Iteration", "Residual ", "Alpha  ", "Beta  ", "Omega  " );
        fprintf( stderr, "----------------------------------------------------------------------\n" );
    }
    *ret = 0;
    for(m=0;m<max;m++)
    {
        rhon = 0.0;
        for(i=0;i<n;i++)
            rhon += r0[i] * r[i];
        beta = ( rhon / rho ) * ( alpha / omega );
        rho = rhon;
        for(i=0;i<n;i++)
            p[i] = r[i] + beta * ( p[i] - omega * v[i] );
        sdgemv( (long) n, (long) n, ia, ja, A, 1, p, 1, v );
        spclw( n, ia, ja, A, v, d );

        f = 0.0;
        for(i=0;i<n;i++)
            f += r0[i] * v[i];
        alpha = rho / f;
        if( isnan( alpha ) )
        {
            *ret = -1; /* Failed to converge */
            break;
        }
        for(i=0;i<n;i++)
            s[i] = r[i] - alpha * v[i];
        sdgemv( (long) n, (long) n, ia, ja, A, 1, s, 1, t );
        spclw( n, ia, ja, A, t, d );

        f = 0.0, g = 0.0;
        for(i=0;i<n;i++)
            f += t[i] * s[i],
              g += t[i] * t[i];
        omega = f / g;
        for(i=0;i<n;i++)
            x[i] += alpha * p[i] + omega * s[i];

        /* Check for convergence here */
        f = 0.0;
        for(i=0;i<n;i++)
            f += r[i] * r[i];
        f = sqrt( f );
        if( vb )
            fprintf( stderr, "%10d%15.7f%15.7f%15.7f%15.7f\n", m, f, alpha, beta, omega );
        if( f < tol )
            break;

        for(i=0;i<n;i++)
            r[i] = s[i] - omega * t[i];
    }
    if( f > tol )
        *ret = -2;

    if( vb )
    {
        fprintf( stderr, "----------------------------------------------------------------------\n" );
        if( *ret == 0 )
            fprintf( stderr, " Converged                              System size:   %15d\n", n );
        else
            fprintf( stderr, " Failed                                 System size:   %15d\n", n );
        fprintf( stderr, "                                        Tolerance:     %15.3e\n", tol );
        fprintf( stderr, "                                        Iterations:    %15d\n", m );
        fprintf( stderr, "                                        Residual:      %15.3e\n", f );
        fprintf( stderr, "                                        ------------------------------\n" );
    }

    free( r0 );
    free( r );
    free( p );
    free( v );
    free( s );
    free( t );
    free( b );
    free( d );
}

/**
 * Projected biconjugate gradient stabilized algorithm
 */
void bicgstabp( int n, int tt, double *A, double *b, int nv, double *V, double *x, int max, double tol, int vb, int *ret )
{
    int i,j,m;
    double f,g,alpha,beta,omega;
    double *r0,*r,*u,*v,*d,*db,*s,*sb,*t;

    r0 = (double*) malloc( n * sizeof(double) );
    r = (double*) malloc( n * sizeof(double) );
    u = (double*) malloc( n * sizeof(double) );
    v = (double*) malloc( n * sizeof(double) );
    d = (double*) malloc( n * sizeof(double) );
    db = (double*) malloc( n * sizeof(double) );
    s = (double*) malloc( n * sizeof(double) );
    sb = (double*) malloc( n * sizeof(double) );
    t = (double*) malloc( n * sizeof(double) );

    /* Initialize */
    project( n, x, nv, V );
    for(i=0;i<n;i++)
    {
        r[i] = b[i];
        if( tt == 0 )
            for(j=0;j<n;j++)
                r[i] -= A[i*n+j] * x[j];
        else
            for(j=0;j<n;j++)
                r[i] -= A[j*n+i] * x[j];
        r0[i] = r[i];
        d[i] = r[i];
    }

    /* Project V out of r */
    project( n, r0, nv, V );

    /* Iterate */
    if( vb )
    {
        fprintf( stderr, "\n" );
        fprintf( stderr, "----------------------------------------------------------------------\n" );
        fprintf( stderr, "%10s%15s%15s%15s%15s\n", " Iteration", "Residual ", "Alpha  ", "Beta  ", "Omega  " );
        fprintf( stderr, "----------------------------------------------------------------------\n" );
    }
    *ret = 0;
    for(m=0;m<max;m++)
    {
        /* Project V out of d */
        copy( n, d, db );
        project( n, db, nv, V );

        /* Calculate alpha */
        f = 0.0, g = 0.0;
        for(i=0;i<n;i++)
            f += r0[i] * r[i];
        for(i=0;i<n;i++)
        {
            t[i] = 0.0;
            if( tt == 0 )
                for(j=0;j<n;j++)
                    t[i] += A[i*n+j] * db[j];
            else
                for(j=0;j<n;j++)
                    t[i] += A[j*n+i] * db[j];
        }
        for(i=0;i<n;i++)
            g += r0[i] * t[i];
        alpha = f / g;

        /* Update s */
        for(i=0;i<n;i++)
            s[i] = r[i] - alpha * t[i];

        /* Project V out of s */
        copy( n, s, sb );
        project( n, sb, nv, V );

        /* Calculate omega */
        f = 0.0, g = 0.0;
        for(i=0;i<n;i++)
        {
            v[i] = 0.0;
            if( tt == 0 )
                for(j=0;j<n;j++)
                    v[i] += A[i*n+j] * sb[j];
            else
                for(j=0;j<n;j++)
                    v[i] += A[j*n+i] * sb[j];
            f += sb[i] * v[i];
        }
        copy( n, v, u );
        project( n, u, nv, V );
        for(i=0;i<n;i++)
            g += v[i] * u[i];
        omega = f / g;

        /* Advance x and the residual and calculate beta */
        f = 0.0, g = 0.0;
        for(i=0;i<n;i++)
            g += r0[i] * r[i];
        for(i=0;i<n;i++)
        {
            x[i] += alpha * db[i] + omega * sb[i];
            r[i] = s[i] - omega * v[i];
            f += r0[i] * r[i];
        }
        beta = alpha / omega * f / g;

        /* Check for convergence */
        if( vb )
            fprintf( stderr, "%10d%15.7f%15.7f%15.7f%15.7f\n", m, f, alpha, beta, omega );
        if( fabs( f ) < tol )
            break;

        /* Advance d */
        for(i=0;i<n;i++)
            d[i] = r[i] + beta * ( d[i] - omega * t[i] );
    }
    if( f > tol )
        *ret = -2;

    if( vb )
    {
        fprintf( stderr, "----------------------------------------------------------------------\n" );
        if( *ret == 0 )
            fprintf( stderr, " Converged                              System size:   %15d\n", n );
        else
            fprintf( stderr, " Failed                                 System size:   %15d\n", n );
        fprintf( stderr, "                                        Tolerance:     %15.3e\n", tol );
        fprintf( stderr, "                                        Iterations:    %15d\n", m );
        fprintf( stderr, "                                        Residual:      %15.3e\n", f );
        fprintf( stderr, "                                        ------------------------------\n" );
    }

    free( r0 );
    free( r );
    free( u );
    free( v );
    free( d );
    free( db );
    free( s );
    free( sb );
    free( t );
}

/**
 * Do forward substitution for the lower triangular part of a given
 * sparse matrix by taking the result of multiplication by strictly
 * upper part; this is for use in Gauss-Seidel iteration
 * @param n Dimension of square system
 * @param ia List of row beginnings
 * @param ja Column indexing
 * @param A Entries of matrix A
 * @param x Current vector iterate
 * @param b Righthand side for forward substitution; b overwritten
 * @param d List of positions of the diagonal elements in each row
 */
void smfsub( long n, long *ia, long *ja, double *A, double *x, double *b, long *d )
{
    long i,j;

    /* Multiply by the strictly upper triangular part; x <- Fx + b */
    for(i=0;i<n;i++)
    {
        x[i] = b[i];
        for(j=d[i]+1;j<ia[i+1];j++)
            x[i] -= A[j] * x[ja[j]];
    }

    /* Do the forward substitution */
    for(i=0;i<n;i++)
    {
        for(j=ia[i];j<d[i];j++)
            x[i] = x[i] - A[j] * x[ja[j]];
        x[i] = x[i] / A[d[i]];
    }
}

/**
 * Gauss-Seidel fixed-point iteration; b is overwritten with the solution
 * @param n Dimension of the square matrix A
 * @param ia Start and end indicator for rows
 * @param ja Column indicators for each entry
 * @param A Elements of matrix A
 * @param b Righthand side of linear system
 * @param max Maximum number of GS iterations to take
 * @param tol Tolerance to use to detect convergence
 * @param vb Verbose if 1, not if 0
 * @param ret Return value; 0 if all went well
 */
void sfpgs( long n, long *ia, long *ja, double *A, double *x, double *b, int max, double tol, int vb, int *ret )
{
    int m,res;
    long i,j;
    long *d = (long*) malloc( n * sizeof(long) );

    /* Find diagonal entries only once */
    *ret = 0;
    for(i=0;i<n;i++)
    {
        d[i] = -1;
        for(j=ia[i];j<ia[i+1];j++)
        {
            if( ja[j] == i )
            {
                d[i] = j;
                break;
            }
        }
        if( d[i] == -1 )
        {
            *ret = -1; /* Missing a diagonal element */
            break;
        }
    }

    /* Iterate */
    if( *ret == 0 )
    {
        for(m=0;m<max;m++)
        {
            smfsub( n, ia, ja, A, x, b, d );
        }
    }
}

/**
 * Tridiagonal LU factorization routine; solves Tx = b; IMPORTANT: The size of
 * the tridiagonal system is actually n - 1 by n - 1
 * @param n Size of the square tridiagonal matrix
 * @param alpha Diagonal elements of tridiagonal matrix; length is n
 * @param beta Subdiagonal elements of tridiagonal matrix; length is n - 1
 * @param gamma Superdiagonal elements of tridiagonal matrix; length is n - 1
 * @param b Right-hand side of the linear system
 * @param x Output of the solution; don't need this; overwrite b instead; this is standard
 */
void trlusv( int n, double *alpha, double *beta, double *gamma, double *b, double *x )
{
    int i;
    double *d,*l,*u,*y;

    /* Allocate temp storage for LU factorzation vectors; try to eliminate this by overwriting alpha, beta, gamma */
    d = (double*) malloc( n * sizeof(double) );
    l = (double*) malloc( n * sizeof(double) );
    u = (double*) malloc( n * sizeof(double) );
    y = (double*) malloc( n * sizeof(double) );

    /* Calculate the LU decomposition of the tridiagonal system in alpha, beta, gamma */
    d[0] = alpha[0];
    u[0] = gamma[0];
    for(i=1;i<n;i++)
    {
        l[i-1] = beta[i-1] / d[i-1];
        d[i] = alpha[i] - l[i-1] * gamma[i-1];
        u[i] = gamma[i]; /* This is redundant, of course */
    }

    /* Do forward/back-substitution to solve the system */
    y[0] = b[0];
    for(i=1;i<n;i++)
        y[i] = b[i] - l[i-1] * y[i-1];
    x[n-1] = y[n-1] / d[n-1];
    for(i=n-2;i>=0;i--)
        x[i] = ( y[i] - u[i] * x[i+1] ) / d[i];

    /* Clean up */
    free( d ); free( l ) ; free( u ); free( y );
}

/**
 * Form a Givens rotation in cos/sin format
 */
void gvrot( double a11, double a21, double *c, double *s )
{
    double ss;
    double p;

    /* Calculate cos and sin components of rank two operation */
    p = a11 / a21;
    ss = 1.0 / ( 1.0 + p * p );
    *s = sqrt( ss );
    *c = sqrt( 1.0 - ss );
}

/**
 * Form Householder reflections into the vector v
 * @param n Vector dimension
 * @param idx Index below which to disappear all entries
 * @param ldx Leading dimension of x (how many entries to skip per index)
 * @param x Vector of inputs
 * @param m Number of vectors in y to which to apply I - 2 uuT
 * @param y Vectors to which to apply I - 2 uuT
 * @param Q Output matrix of QR transform
 */
void hhref( int n, int idx, int ldx, double *x, int m, double *y, int nq, double *Q )
{
    int i,j;
    double len,lem,*u,*w;

    /* Allocate */
    u = (double*) malloc( n * sizeof(double) );

    /* Compute length of original vector; len must be preserved by transformed vector */
    len = 0.0;
    for(i=0;i<n;i++)
        len += x[i*ldx] * x[i*ldx];
    len = sqrt( len );

    /* Make u */
    lem = 0.0;
    for(i=0;i<n;i++)
    {
        if( i == idx )
            u[i] = x[i*ldx] - len;
        else
            u[i] = x[i*ldx];
        lem += u[i] * u[i];
    }
    lem = sqrt( lem );
    for(i=0;i<n;i++)
        u[i] /= lem; /* This is the completed normalized reflection vector */

    /* Apply I - 2 uuT to x and put it in y */
    for(i=0;i<m;i++)
    {
        w = y + i;
        len = 0.0;
        for(j=0;j<n;j++)
            len += w[j*ldx] * u[j];
        for(j=0;j<n;j++)
            w[j*ldx] = w[j*ldx] - 2.0 * len * u[j];
    }
    for(i=0;i<nq;i++)
    {
        w = Q + i;
        len = 0.0;
        for(j=0;j<n;j++)
            len += w[j*nq] * u[j];
        for(j=0;j<n;j++)
            w[j*nq] = w[j*nq] - 2.0 * len * u[j];
    }
}

void hhqrls( int n, int m, double *A, int nq, double *Q )
{
    int i,j,k;
    double *q;

    /* Apply the sequence of reflections */
    for(k=0;k<m;k++)
        hhref( n - k, 0, m, A + k * ( m + 1 ), m - k, A + k * ( m + 1 ), nq, Q + k * nq );

    /* Now invert the smaller matrix using back substitution */
    for(k=0;k<nq;k++)
    {
        q = Q + k;
        for(i=m-1;i>=0;i--)
        {
            for(j=m-1;j>i;j--)
                q[i*nq] = q[i*nq] - A[i*m+j] * q[j*nq];
            q[i*nq] = q[i*nq] / A[i*m+i];
        }
    }
}

/**
 * Apply Givens rotation to a tridiagonal matrix; this
 * requires pentadiagonal storage in order to account
 * for new temporary non-zeros
 */
void tdgvrot( int n, double *A, int col, int row1, int row2 )
{
    int i;
    double a11,a21,c,s;

    /* Form matrix with right entries of A */
    a11 = A[row1*n+col];
    a21 = A[row2*n+col];
    gvrot( a11, a21, &c, &s );
}

/**
 * Basic Lanczos unsymmetric tridiagonalization method; store the
 * tridiagonal matrix in band format
 */
void ulanczos( int n, int m, double *A, double *alpha, double *beta, double *gamma, double *r0, double *s0, int *res )
{
    int i,j,res1,res2,res3;
    double *q,*qq,*p,*pp,*r,*s;

    /* Allocate storage for iteration vector */
    q = (double*) malloc( n * sizeof(double) );
    qq = (double*) malloc( n * sizeof(double) );
    p = (double*) malloc( n * sizeof(double) );
    pp = (double*) malloc( n * sizeof(double) );
    r = (double*) malloc( n * sizeof(double) );
    s = (double*) malloc( n * sizeof(double) );

    /* Generate an initially orthogonal pair of vectors in p, q */
    zerov( n, q );
    zerov( n, p );
    if( r0 != NULL )
        copy( n, r0, r );
    else
        nrandv( n, r );
    if( s0 != NULL )
        copy( n, s0, s );
    else
        nrandv( n, s );

    /* Loop until residual vector r, s vanish or maximum basis size, m, is reached */
    *res = 0;
    for(i=0;i<m;)
    {
        /* Make sure that r, s do not vanish and are not orthogonal */
        normchk( n, r, CGRAD_ZERO, &res1 );
        normchk( n, s, CGRAD_ZERO, &res2 );
        dotchk( n, r, s, CGRAD_ZERO, &res3 );

        /* Return code which tells which conditions failed */
        if( res1 == -1 || res2 == -1 || res3 == -1 )
        {
            if( res1 == -1 )
                *res |= ( 1 << 0 );
            if( res2 == -1 )
                *res |= ( 1 << 1 );
            if( res3 == -1 )
                *res |= ( 1 << 2 );
            break;
        }

        /* Build the next vectors and coefficients */
#ifdef CGRAD_LANCZOS_NORM
        norm( n, r, &beta[i] );
#else
        dotp( n, s, r, &beta[i] );
        beta[i] = sqrt( fabs( beta[i] ) );
#endif
        dotpdiv( n, s, r, beta[i], &gamma[i] );
        copyvdiv( n, r, beta[i], qq );
        copyvdiv( n, s, gamma[i], pp );

        /* Step forward one */
        ++i;

        /* Calcualte the next diagonal term */
        mdotp( n, A, pp, qq, &alpha[i] );
        smmadd( n, 0, A, alpha[i], qq, q, -1.0 * gamma[i-1], r );
        smmadd( n, 1, A, alpha[i], pp, p, -1.0 * beta[i-1], s );
        copy( n, qq, q );
        copy( n, pp, p );
    }

    free( q ); free( qq ); free( p ); free( pp ); free( r ); free( s );
}

/**
 * Basic Lanczos unsymmetric tridiagonalization method; store the
 * tridiagonal matrix in band format; this version saves the dual bases
 * along the way for later use
 */
void ulanczosb( int n, int m, double *A, double *alpha, double *beta, double *gamma, double *r0, double *s0, double *V, double *U, int *res )
{
    int i,j,res1,res2,res3;
    double *q,*qq,*p,*pp,*r,*s;

    /* Allocate storage for iteration vector */
    q = (double*) malloc( n * sizeof(double) );
    qq = (double*) malloc( n * sizeof(double) );
    p = (double*) malloc( n * sizeof(double) );
    pp = (double*) malloc( n * sizeof(double) );
    r = (double*) malloc( n * sizeof(double) );
    s = (double*) malloc( n * sizeof(double) );

    /* Generate an initially orthogonal pair of vectors in p, q */
    zerov( n, q );
    zerov( n, p );
    if( r0 != NULL )
        copy( n, r0, r );
    else
        nrandv( n, r );
    if( s0 != NULL )
        copy( n, s0, s );
    else
        nrandv( n, s );

    /* Loop until residual vector r, s vanish or maximum basis size, m, is reached */
    *res = 0;
    for(i=0;i<m;)
    {
        /* Make sure that r, s do not vanish and are not orthogonal */
        normchk( n, r, CGRAD_ZERO, &res1 );
        normchk( n, s, CGRAD_ZERO, &res2 );
        dotchk( n, r, s, CGRAD_ZERO, &res3 );

        /* Return code which tells which conditions failed */
        if( res1 == -1 || res2 == -1 || res3 == -1 )
        {
            if( res1 == -1 )
                *res |= ( 1 << 0 );
            if( res2 == -1 )
                *res |= ( 1 << 1 );
            if( res3 == -1 )
                *res |= ( 1 << 2 );
            break;
        }

        /* Build the next vectors and coefficients */
#ifdef CGRAD_LANCZOS_NORM
        norm( n, r, &beta[i] );
#else
        dotp( n, s, r, &beta[i] );
        beta[i] = sqrt( fabs( beta[i] ) );
#endif
        dotpdiv( n, s, r, beta[i], &gamma[i] );
        copyvdiv( n, r, beta[i], qq );
        copyvdiv( n, s, gamma[i], pp );

        /* Save in the basis set if no breakdown occurs */
        if( V != NULL )
            copy( n, qq, V + i * n );
        if( U != NULL )
            copy( n, pp, U + i * n );

        /* Step forward one */
        ++i;

        /* Calcualte the next diagonal term */
        mdotp( n, A, pp, qq, &alpha[i] );
        smmadd( n, 0, A, alpha[i], qq, q, -1.0 * gamma[i-1], r );
        smmadd( n, 1, A, alpha[i], pp, p, -1.0 * beta[i-1], s );
        copy( n, qq, q );
        copy( n, pp, p );
    }

    free( q ); free( qq ); free( p ); free( pp ); free( r ); free( s );
}

/**
 * Lanczos matrix inversion using tridiagonal operations and projection
 */
void gelzsvm( int n, int m, double *A, double *b, double *x, double *alpha, double *beta, double *gamma, double *V, double *U )
{
    int i,j,res;
    double beta0,tmp,*r0,*s0,*tb,*tx;

    r0 = (double*) malloc( n * sizeof(double) );
    s0 = (double*) malloc( n * sizeof(double) );
    tb = (double*) malloc( m * sizeof(double) );
    tx = (double*) malloc( m * sizeof(double) );

    /* Calculate the original residual for use in forming tb */
    residual( n, A, x, b, r0 );
    norm( n, r0, &beta0 );
    normalize( n, r0 );
    nrandv( n, s0 );
    dotp( n, s0, r0, &tmp );
    vdiv( n, s0, tmp );

    /* Build the Lanczos matrix and vectors */
    ulanczosb( n, m, A, alpha, beta, gamma, r0, s0, V, U, &res );

    /* Form the right-hand side of the projected tridiagonal problem */
    for(i=0;i<m;i++)
        if( i == 0 )
            tb[i] = beta0;
        else
            tb[i] = 0.0;

    /* Solve the tridiagonal system; put solution in tx */
    trlusv( m, alpha + 1, beta + 1, gamma + 1, tb, tx );

    /* Now put the solution back into the original vector space */
    for(i=0;i<m;i++)
        for(j=0;j<n;j++)
            x[j] += tx[i] * V[i*n+j];

    free( r0 ); free( s0 ); free( tb ); free( tx );
}

#define CGRAD_GELZSVT_INC 256

inline void vchkrealloc( int n, int size, int *alloc, double *x, double *b, double *V, double *U )
{
    if( size + 1 > *alloc )
    {
        x = (double*) realloc( x, ( *alloc + CGRAD_GELZSVT_INC ) * sizeof(double) );
        b = (double*) realloc( b, ( *alloc + CGRAD_GELZSVT_INC ) * sizeof(double) );
        V = (double*) realloc( V, ( *alloc + CGRAD_GELZSVT_INC ) * n * sizeof(double) );
        U = (double*) realloc( U, ( *alloc + CGRAD_GELZSVT_INC ) * n * sizeof(double) );
        *alloc += CGRAD_GELZSVT_INC;
    }
}

/**
 * Lanczos matrix inversion using tridiagonal operations and projection;
 * this version calculates a convergence estimate and terminates when the
 * tolerance criterion is met; this version does inline Lanczos basis calculation
 */
void gelzsv( int n, int *m, double *A, double *b, double *x, double *alpha, double *beta, double *gamma, double tol, int *res, int vb )
{
    int i,j,talloc,res1,res2,res3;
    double beta0,tmp,err,*r0,*s0,*tb,*tx,*q,*qq,*p,*pp,*r,*s,*V,*U;

    r0 = (double*) malloc( n * sizeof(double) );
    s0 = (double*) malloc( n * sizeof(double) );
    q = (double*) malloc( n * sizeof(double) );
    qq = (double*) malloc( n * sizeof(double) );
    p = (double*) malloc( n * sizeof(double) );
    pp = (double*) malloc( n * sizeof(double) );
    r = (double*) malloc( n * sizeof(double) );
    s = (double*) malloc( n * sizeof(double) );

    /* Initialize tridiagonal space */
    tx = (double*) malloc( CGRAD_GELZSVT_INC * sizeof(double) );
    tb = (double*) malloc( CGRAD_GELZSVT_INC * sizeof(double) );
    V = (double*) malloc( CGRAD_GELZSVT_INC * n * sizeof(double) );
    U = (double*) malloc( CGRAD_GELZSVT_INC * n * sizeof(double) );
    talloc = CGRAD_GELZSVT_INC;

    /* Write results if verbose on */
    if( vb )
    {
        fprintf( stderr, "\n" );
        fprintf( stderr, "-------------------------\n" );
        fprintf( stderr, "%10s%15s\n", " Iteration", "Residual " );
        fprintf( stderr, "-------------------------\n" );
    }

    /* Start setting up Lanczos steps */
    zerov( n, q );
    zerov( n, p );

    /* Calculate the original residual for use in forming tb */
    residual( n, A, x, b, r0 );
    norm( n, r0, &beta0 );
    normalize( n, r0 );
    nrandv( n, s0 );
    dotp( n, s0, r0, &tmp );
    vdiv( n, s0, tmp );
    copy( n, r0, r );
    copy( n, s0, s );

    /* Loop until residual vector r, s vanish or maximum basis size, m, is reached */
    *res = 0;
    tb[0] = beta0;
    for(i=0;i<n-1;)
    {
        /* Make sure that r, s do not vanish and are not orthogonal */
        normchk( n, r, CGRAD_ZERO, &res1 );
        normchk( n, s, CGRAD_ZERO, &res2 );
        dotchk( n, r, s, CGRAD_ZERO, &res3 );

        /* Return code which tells which conditions failed */
        if( res1 == -1 || res2 == -1 || res3 == -1 )
        {
            if( res1 == -1 )
                *res |= ( 1 << 0 );
            if( res2 == -1 )
                *res |= ( 1 << 1 );
            if( res3 == -1 )
                *res |= ( 1 << 2 );
            break;
        }

        /* Build the next vectors and coefficients */
        dotp( n, s, r, &beta[i] );
        beta[i] = sqrt( fabs( beta[i] ) );
        dotpdiv( n, s, r, beta[i], &gamma[i] );

        /* Copy new vectors into V and U, respectively */
        copyvdiv( n, r, beta[i], qq );
        copyvdiv( n, s, gamma[i], pp );

        /* Save in the basis set if no breakdown occurs */
        if( V != NULL )
            copy( n, qq, V + i * n );
        if( U != NULL )
            copy( n, pp, U + i * n );

        /* Step forward one */
        ++i;

        /* Calculate the next diagonal term */
        mdotp( n, A, pp, qq, &alpha[i] );

        /* Allocate more space for tx, tb if needed */
        vchkrealloc( n, i, &talloc, tx, tb, V, U );

        /* Form the right-hand side of the projected tridiagonal problem */
        for(j=1;j<i;j++)
            tb[j] = 0.0;

        /* Solve the tridiagonal system for an error estimate */
        copy( n, tx, r0 ); /* Save old solution */
        trlusv( i, alpha + 1, beta + 1, gamma + 1, tb, tx );
        vnormdiff( i, r0, tx, &err );
        if( err < tol )
            break;

        if( vb )
            fprintf( stderr, "%10d%15.7f\n", i, err );

        /* Prepare vectors qq, pp for next round */
        smmadd( n, 0, A, alpha[i], qq, q, -1.0 * gamma[i-1], r );
        smmadd( n, 1, A, alpha[i], pp, p, -1.0 * beta[i-1], s );
        copy( n, qq, q );
        copy( n, pp, p );
    }
    *m = i;

    if( vb )
    {
        fprintf( stderr, "-------------------------\n" );
        fprintf( stderr, "Tolerance: %14.3e\n", tol );
        fprintf( stderr, " Residual: %14.3e\n", err );
    }

    /* Now put the solution back into the original vector space */
    for(i=0;i<*m;i++)
        for(j=0;j<n;j++)
            x[j] += tx[i] * V[i*n+j];

    free( r0 ); free( s0 ); free( tb ); free( tx ); free( tx ); free( tb );
    free( V ); free( U ); free( q ); free( qq ); free( p ); free( pp );
}

/**
 * Basic generalized Lanczos unsymmetric tridiagonalization method; store the
 * tridiagonal matrix in band format; biconjugate gradient stabilized method is
 * used to apply the the inverse of B iteratively
 */
void gulanczos( int n, int m, double *A, double *B, double *alpha, double *beta, double *gamma, double *r0, double *s0, double tol, int max, int vb, int *res )
{
    /* Must apply the inverse of B at each step!!! */
    int i,j,res1,res2,res3;
    double *q,*qq,*p,*pp,*r,*s,*t,*u,*v,*w;

    /* Allocate storage for iteration vector */
    q = (double*) malloc( n * sizeof(double) );
    qq = (double*) malloc( n * sizeof(double) );
    p = (double*) malloc( n * sizeof(double) );
    pp = (double*) malloc( n * sizeof(double) );
    r = (double*) malloc( n * sizeof(double) );
    s = (double*) malloc( n * sizeof(double) );
    t = (double*) malloc( n * sizeof(double) );
    w = (double*) malloc( n * sizeof(double) );
    u = (double*) malloc( n * sizeof(double) );
    v = (double*) malloc( n * sizeof(double) );

    /* Generate an initially orthogonal pair of vectors in p, q */
    zerov( n, q );
    zerov( n, p );
    if( r0 != NULL )
        copy( n, r0, r );
    else
        nrandv( n, r );
    if( s0 != NULL )
        copy( n, s0, s );
    else
        nrandv( n, s );

    /* Loop until residual vector r, s vanish or maximum basis size, m, is reached */
    *res = 0;
    for(i=0;i<m;)
    {
        /* Make sure that r, s do not vanish and are not orthogonal */
        normchk( n, r, CGRAD_ZERO, &res1 );
        normchk( n, s, CGRAD_ZERO, &res2 );
        dotchk( n, r, s, CGRAD_ZERO, &res3 );

        /* Return code which tells which conditions failed */
        if( res1 == -1 || res2 == -1 || res3 == -1 )
        {
            if( res1 == -1 )
                *res |= ( 1 << 0 );
            if( res2 == -1 )
                *res |= ( 1 << 1 );
            if( res3 == -1 )
                *res |= ( 1 << 2 );
            break;
        }

        /* Build the next vectors and coefficients */
#ifdef CGRAD_LANCZOS_NORM
        norm( n, r, &beta[i] );
#else
        dotp( n, s, r, &beta[i] );
        beta[i] = sqrt( fabs( beta[i] ) );
#endif
        dotpdiv( n, s, r, beta[i], &gamma[i] );
        copyvdiv( n, r, beta[i], qq );
        copyvdiv( n, s, gamma[i], pp );

        /* Step forward one */
        ++i;

        /* Everything changes here because we need to multiply by the inverse of B */
        dgemv( n, 0, A, qq, t );
        bicgstab( n, 0, B, t, u, max, tol, vb, &res1 ); /* Vector uk is in r0 */
        bicgstab( n, 1, B, pp, w, max, tol, vb, &res2 ); /* Vector wk is in r0 */
        dgemv( n, 1, A, w, v );

        /* Now u, v are formed so calculate alpha[i], r[i], s[i] */
        dotp( n, pp, u, &alpha[i] );
        vsum( n, 1.0, u, -1.0 * alpha[i], qq, r );
        vsum( n, 1.0, r, -1.0 * gamma[i-1], q, r );
        vsum( n, 1.0, v, -1.0 * alpha[i], pp, s );
        vsum( n, 1.0, s, -1.0 * beta[i-1], p, s );
        copy( n, qq, q );
        copy( n, pp, p );
    }

    free( q ); free( qq ); free( p ); free( pp ); free( r ); free( s );
    free( t ); free( w ); free( u ); free( v );
}

void subdchk( int n, double *l, double *d, double *u, int *ret )
{
    int i;

    *ret = 0;
    for(i=0;i<n-1;i++)
        if( fabs( l[i] ) > CGRAD_ZERO && fabs( l[i+1] ) > CGRAD_ZERO )
            *ret = 1;
}

/**
 * Basic generalized Lanczos unsymmetric tridiagonalization method with sparse matrix format;
 * store the tridiagonal matrix in band format; biconjugate gradient stabilized method is
 * used to apply the the inverse of B iteratively
 * @param n Dimension of the matrix
 * @param m Dimension of the Krylov space to build
 * @param ia Row index vector of A
 * @param ja Column indicator of A
 * @param A Floating point values of A
 * @param ib Row index vector of B
 * @param jb Column indicator of B
 * @param B Floating point vector of B
 * @param alpha Output diagonal of the Krylov matrix
 * @param beta Output subdiagonal of the Krylov matrix
 * @param gamma Output superdiagonal of the Krylov matrix
 * @param r0 Initial r given as input; NULL means generate random
 * @param s0 Initial s given as input; NULL means generate random
 * @param tol Error tolerance to use
 * @param max Maximum number of iterations to take in the inversion of B via bicgstab
 * @param vb Verbose level; 0 means off, anything else means on
 * @param res Return value indicating success or failure
 */
void sgulanczos( int n, int m, long *ia, long *ja, double *A, long *ib, long* jb, double *B, double *alpha, double *beta, double *gamma, double *r0, double *s0, double tol, int max, int vb, int *res )
{
    /* Must apply the inverse of B at each step!!! */
    int i,j,res1,res2,res3;
    long *iat,*jat,*ibt,*jbt;
    double *q,*qq,*p,*pp,*r,*s,*t,*u,*v,*w,*At,*Bt;

    /* Allocate storage for iteration vector */
    q = (double*) malloc( n * sizeof(double) );
    qq = (double*) malloc( n * sizeof(double) );
    p = (double*) malloc( n * sizeof(double) );
    pp = (double*) malloc( n * sizeof(double) );
    r = (double*) malloc( n * sizeof(double) );
    s = (double*) malloc( n * sizeof(double) );
    t = (double*) malloc( n * sizeof(double) );
    w = (double*) malloc( n * sizeof(double) );
    u = (double*) malloc( n * sizeof(double) );
    v = (double*) malloc( n * sizeof(double) );

    /* Generate an initially orthogonal pair of vectors in p, q */
    zerov( n, q );
    zerov( n, p );
    if( r0 != NULL )
        copy( n, r0, r );
    else
        nrandv( n, r );
    if( s0 != NULL )
        copy( n, s0, s );
    else
        nrandv( n, s );

    /* Transposes of the A and B matrices */
    iat = (long*) malloc( ( n + 1 ) * sizeof(long) );
    jat = (long*) malloc( ia[n] * sizeof(long) );
    At = (double*) malloc( ia[n] * sizeof(double) );
    ibt = (long*) malloc( ( n + 1 ) * sizeof(long) );
    jbt = (long*) malloc( ib[n] * sizeof(long) );
    Bt = (double*) malloc( ib[n] * sizeof(double) );
    stransp( 1, (long) n, (long) n, ia, ja, A, iat, jat, At ); /* There must be a way to get rid of this */
    stransp( 1, (long) n, (long) n, ib, jb, B, ibt, jbt, Bt ); /* At least should not have to store the actual entries */

    /* Loop until residual vector r, s vanish or maximum basis size, m, is reached */
    *res = 0;
    for(i=0;i<m;)
    {
        /* Make sure that r, s do not vanish and are not orthogonal */
        normchk( n, r, CGRAD_ZERO, &res1 );
        normchk( n, s, CGRAD_ZERO, &res2 );
        dotchk( n, r, s, CGRAD_ZERO, &res3 );

        /* Return code which tells which conditions failed */
        if( res1 == -1 || res2 == -1 || res3 == -1 )
        {
            if( res1 == -1 )
                *res |= ( 1 << 0 );
            if( res2 == -1 )
                *res |= ( 1 << 1 );
            if( res3 == -1 )
                *res |= ( 1 << 2 );
            break;
        }

        /* Build the next vectors and coefficients */
#ifdef CGRAD_LANCZOS_NORM
        norm( n, r, &beta[i] );
#else
        dotp( n, s, r, &beta[i] );
        beta[i] = sqrt( fabs( beta[i] ) );
#endif
        dotpdiv( n, s, r, beta[i], &gamma[i] );
        copyvdiv( n, r, beta[i], qq );
        copyvdiv( n, s, gamma[i], pp );

        /* Step forward one */
        ++i;

        /* Everything changes here because we need to multiply by the inverse of B */
        sdgemv( n, n, ia, ja, A, 1, qq, 1, t );
        sbicgstab( n, ib, jb, B, t, u, max, tol, vb, &res1 ); /* Vector uk is in r0 */
        sbicgstab( n, ibt, jbt, Bt, pp, w, max, tol, vb, &res2 ); /* Vector wk is in r0 */
        sdgemv( n, n, iat, jat, At, 1, w, 1, v );

        /* Now u, v are formed so calculate alpha[i], r[i], s[i] */
        dotp( n, pp, u, &alpha[i] );
        vsum( n, 1.0, u, -1.0 * alpha[i], qq, r );
        vsum( n, 1.0, r, -1.0 * gamma[i-1], q, r );
        vsum( n, 1.0, v, -1.0 * alpha[i], pp, s );
        vsum( n, 1.0, s, -1.0 * beta[i-1], p, s );
        copy( n, qq, q );
        copy( n, pp, p );
    }

    free( q ); free( qq ); free( p ); free( pp ); free( r ); free( s );
    free( t ); free( w ); free( u ); free( v );
    free( iat ); free( jat ); free( At ); free( ibt ); free( jbt ); free( Bt );
}

/**
 * Symmetric version of sgulanczos
 */
void sslanczos( int n, int m, long *ia, long *ja, double *A, double *alpha, double *beta, double *r0, double tol, int vb, int *res )
{
    int i,j,ret;
    double *v,*vv,*w;

    /* Allocate everything */
    v = (double*) malloc( n * sizeof(double) );
    vv = (double*) malloc( n * sizeof(double) );
    w = (double*) malloc( n * sizeof(double) );

    /* Generate an initial vector in v */
    zerov( n, vv );
    if( r0 != NULL )
        copy( n, r0, v );
    else
        nrandv( n, v );
    beta[0] = 0.0;

    /* Start iteration */
    for(i=0;i<m;i++)
    {
        ssdgemv( (long) n, ia, ja, A, 1, v, 1, w );
        dotp( n, w, v, &alpha[i] );
        for(j=0;j<n;j++)
            w[j] = w[j] - alpha[i] * v[j] - beta[i] * vv[j];
        norm( n, w, &beta[i+1] );
        copy( n, v, vv );
        copyvdiv( n, w, beta[i+1], v );
    }

    free( v ); free( vv ); free( w );
}

/**
 * Sparse generalized symmetric indefinite Lanczos procedure; matrices A and B
 * are assumed to be symmetric but not necessarily definite; eigenvalues can be
 * extracted from the symmetric tridiagonal matrix (alpha,beta) but must be
 * scaled using omega for each eigenvalue
 * @param n Vector dimension
 * @param m Number of steps to take; also the final dimension of the tridiagonal matrix
 * @param ia Row index list of A
 * @param ja Column indexes of A
 * @param A Floating point entries of A
 * @param ib Row index list of B
 * @param jb Column indexes of B
 * @param B Floating point entries of B
 * @param alpha Diagonal entries of tridiagonal
 * @param beta Subdiagonal/superdiagonal entries of tridiagonal
 * @param omega Multiplier vector of A-norms of vectors
 * @param r0 Optional initial Krylov space vector
 * @param tol Tolerance determining convergence
 * @param max Maximum number of iterations to use in conjugate gradient
 * @param vb Verbose level
 * @param res Return value
 */
void sgsilanczos( int n, int m, long *ia, long *ja, double *A, long *ib, long* jb, double *B, double *alpha, double *beta, double *omega, double *r0, double tol, int max, int vb, int *res )
{
    int i,ret;
    double *w,*q,*s,*qq,tau;

    /* Allocate stuff */
    w = (double*) malloc( n * sizeof(double) );
    q = (double*) malloc( n * sizeof(double) );
    s = (double*) malloc( n * sizeof(double) );
    qq = (double*) malloc( n * sizeof(double) );

    /* If input start vector given, use it; otherwise generate random one */
    if( r0 != NULL )
        copy( n, r0, q );
    else
        nrandv( n, q );

    /* Multiply the initial vector */
    ssdgemv( (long) n, ia, ja, A, 1, q, 1, w );
    dotp( n, q, w, &omega[0] );

    /* Start iteration */
    for(i=0;i<m;i++)
    {
        /* Solve Bs = w to get s = inv(B) * w */
        scg( n, ib, jb, B, w, s, max, tol, vb, &ret );

        /* Only conditionally calculate the first sum */
        if( i > 0 )
            vsum( n, 1.0, s, -beta[i] / omega[i-1], qq, s ); /* qq is the previous step's vector */
        dotp( n, w, s, &alpha[i] );
        vsum( n, 1.0, s, -alpha[i] / omega[i], q, s );
        norm( n, s, &tau );

        /* Stop if tolerance is reached */
        if( fabs( tau ) < tol )
            break;

        /* Take the next step before jumping back to beginning of loop */
        copy( n, q, qq );
        copyvdiv( n, s, tau, q );

        /* Multiply again */
        ssdgemv( (long) n, ia, ja, A, 1, q, 1, w );
        dotp( n, w, q, &omega[i+1] );

        /* Stop if omega vanishes */
        if( fabs( omega[i+1] ) < tol )
            break;
        beta[i+1] = tau * omega[i+1];
    }

    free( w ); free( q ); free( s ); free( qq );
}

/**
 * Need this function for sgsilanczosc for convergence testing
 */
void treiglr( int, double *, double *, double *, int, double *, int *, int * );
void treigqds( int, double *, double *, double *, int, double *, int *, int * );
void treigdqds( int, double *, double *, double *, int, double *, int *, int * );
void treigtridqds( int, double *, double *, double *, int, double *, int *, int * );

/**
 * Same as sgsilanczos except also has argument neig; the old argument m
 * just indicates the maximum size of tridiagonal to build, whereas neig
 * is the number of eigenvalues to test for convergence; these will be
 * the eigenvalues of a potentially transformed system which have the
 * largest magnitude
 * @param n Vector dimension
 * @param m Maximum number of Lanczos iterations
 * @param neig The number of eigenvalues sought
 * @param ia Row index list of A
 * @param ja Column indexes of A
 * @param A Floating point entries of A
 * @param ib Row index list of B
 * @param jb Column indexes of B
 * @param B Floating point entries of B
 * @param alpha Diagonal entries of tridiagonal
 * @param beta Subdiagonal/superdiagonal entries of tridiagonal
 * @param omega Multiplier vector of A-norms of vectors
 * @param r0 Optional initial Krylov space vector
 * @param ev The output eigenvalues
 * @param ns Actual size of the output tridiagonal matrix generated
 * @param tol Tolerance determining convergence
 * @param max Maximum number of iterations to use in conjugate gradient
 * @param vb Verbose level
 * @param res Return value
 */
void sgsilanczosc( int n, int m, int neig, long *ia, long *ja, double *A, long *ib, long* jb, double *B, double *alpha, double *beta, double *omega, double *r0, double *ev, int *mo, double tol, int max, int vb, int *res )
{
    int i,j,ns,ret;
    double *w,*q,*s,*qq,*fv,tau,r,t;
    double *talpha,*tbeta,*tgamma;

    /* Allocate stuff */
    w = (double*) malloc( n * sizeof(double) );
    q = (double*) malloc( n * sizeof(double) );
    s = (double*) malloc( n * sizeof(double) );
    qq = (double*) malloc( n * sizeof(double) );
    fv = (double*) malloc( 2 * ( m + 2 ) * sizeof(double) );
    talpha = (double*) malloc( ( m + 2 ) * sizeof(double) );
    tbeta = (double*) malloc( ( m + 2 ) * sizeof(double) );
    tgamma = (double*) malloc( ( m + 2 ) * sizeof(double) );

    /* If input start vector given, use it; otherwise generate random one */
    if( r0 != NULL )
        copy( n, r0, q );
    else
        nrandv( n, q );

    /* Multiply the initial vector */
    ssdgemv( (long) n, ia, ja, A, 1, q, 1, w );
    dotp( n, q, w, &omega[0] );

    /* Start iteration */
    for(i=0;i<m;i++)
    {
        /* Solve Bs = w to get s = inv(B) * w */
        scg( n, ib, jb, B, w, s, max, tol, vb, &ret );

        /* Only conditionally calculate the first sum */
        if( i > 0 )
            vsum( n, 1.0, s, -beta[i] / omega[i-1], qq, s ); /* qq is the previous step's vector */
        dotp( n, w, s, &alpha[i] );
        vsum( n, 1.0, s, -alpha[i] / omega[i], q, s );
        norm( n, s, &tau );

        /* Stop if tolerance is reached */
        if( fabs( tau ) < tol )
            break;

        /* Take the next step before jumping back to beginning of loop */
        copy( n, q, qq );
        copyvdiv( n, s, tau, q );

        /* Multiply again */
        ssdgemv( (long) n, ia, ja, A, 1, q, 1, w );
        dotp( n, w, q, &omega[i+1] );

        /* Stop if omega vanishes */
        if( fabs( omega[i+1] ) < tol )
            break;
        beta[i+1] = tau * omega[i+1];

        /* Do a convergence test */
        if( ( i + 1 ) >= neig )
        {
            /* Generate the right matrix */
            copy( ( i + 1 ) + 2, alpha, talpha );
            copy( ( i + 1 ) + 2, beta, tbeta );
            copy( ( i + 1 ) + 2, tbeta, tgamma );
            for(j=0;j<(i+1)-1;j++)
                talpha[j] /= omega[j], tbeta[j+1] /= omega[j], tgamma[j+1] /= omega[j+1];
            talpha[(i+1)-1] /= omega[(i+1)-1];

            /* Save previous results and update */
            copy( 2 * ( i + 1 ), ev, fv );
            treiglr( i + 1, talpha, tbeta + 1, tgamma + 1, 10000, ev, &ret, &ns );
            inverse_complex_bubble_sort( i + 1, ev );

            /* Compare current to previous eigenvalues of interest */
            vnormdiff( 2 * neig, ev, fv, &r );
            t = 0.0; /* FIXME: Left off here; need to calculate change relative to magnitude of the eigenvalues, not absolute */
        }
    }
    *mo = i;

    free( w ); free( q ); free( s ); free( qq ); free( fv );
    free( talpha ); free( tbeta ); free( tgamma );
}

/**
 * Same as sgsilanczosc except also saves all vectors
 * @param n Vector dimension
 * @param m Maximum number of Lanczos iterations
 * @param neig The number of eigenvalues sought
 * @param ia Row index list of A
 * @param ja Column indexes of A
 * @param A Floating point entries of A
 * @param ib Row index list of B
 * @param jb Column indexes of B
 * @param B Floating point entries of B
 * @param alpha Diagonal entries of tridiagonal
 * @param beta Subdiagonal/superdiagonal entries of tridiagonal
 * @param omega Multiplier vector of A-norms of vectors
 * @param r0 Optional initial Krylov space vector
 * @param ev The output eigenvalues
 * @param Q Storage for vectors for explicit reorthogonalization
 * @param mo Number of iterations actually taken
 * @param tol Tolerance determining convergence
 * @param fr Final residual on termination
 * @param max Maximum number of iterations to use in conjugate gradient
 * @param vb Verbose level
 * @param res Return value
 */
void sgsilanczoscr( int n, int m, int neig, long *ia, long *ja, double *A, long *ib, long* jb, double *B, double *alpha, double *beta, double *omega, double *r0, double *ev, double *Q, int *mo, double tol, double *fr, int max, int vb, int *res )
{
    int i,j,ns,ret;
    double *w,*q,*s,*qq,*rr,*fv,tau,r,t;
    double *talpha,*tbeta,*tgamma;
    double a,b;

    /* Allocate stuff */
    w = (double*) malloc( n * sizeof(double) );
    q = (double*) malloc( n * sizeof(double) );
    s = (double*) malloc( n * sizeof(double) );
    qq = (double*) malloc( n * sizeof(double) );
    rr = (double*) malloc( n * sizeof(double) );
    fv = (double*) malloc( 2 * ( m + 2 ) * sizeof(double) );
    talpha = (double*) malloc( ( m + 2 ) * sizeof(double) );
    tbeta = (double*) malloc( ( m + 2 ) * sizeof(double) );
    tgamma = (double*) malloc( ( m + 2 ) * sizeof(double) );

    /* If input start vector given, use it; otherwise generate random one */
    if( r0 != NULL )
        copy( n, r0, q );
    else
        nrandv( n, q );

    /* Multiply the initial vector */
    copy( n, q, Q );
    ssdgemv( (long) n, ia, ja, A, 1, q, 1, w );
    dotp( n, q, w, &omega[0] );

    /* Start iteration */
    for(i=0;i<m;i++)
    {
        /* Solve Bs = w to get s = inv(B) * w */
        scg( n, ib, jb, B, w, s, max, tol, vb, &ret );

        /* Only conditionally calculate the first sum */
        if( i > 0 )
            vsum( n, 1.0, s, -beta[i] / omega[i-1], qq, s ); /* qq is the previous step's vector */
        dotp( n, w, s, &alpha[i] );
        vsum( n, 1.0, s, -alpha[i] / omega[i], q, s );
        norm( n, s, &tau );

        /* Stop if tolerance is reached */
        if( fabs( tau ) < tol )
            break;

        /* Take the next step before jumping back to beginning of loop */
        copy( n, q, qq );
        copyvdiv( n, s, tau, q );

        /* Do projection to restore orthogonality */
        for(j=0;j<i+1;j++)
        {
            ssdgemv( (long) n, ia, ja, A, 1, &Q[j*n], 1, rr );
            dotp( n, q, rr, &a );
            dotp( n, &Q[j*n], rr, &b );
            vsum( n, 1.0, q, -a / b, &Q[j*n], q );
        }
        copy( n, q, &Q[(i+1)*n] );

        /* Multiply again */
        ssdgemv( (long) n, ia, ja, A, 1, q, 1, w );
        dotp( n, w, q, &omega[i+1] );

        /* Stop if omega vanishes */
        if( fabs( omega[i+1] ) < tol )
            break;
        beta[i+1] = tau * omega[i+1];

        /* Do a convergence test */
        if( ( i + 1 ) >= neig )
        {
            /* Generate the right matrix */
            copy( ( i + 1 ) + 2, alpha, talpha );
            copy( ( i + 1 ) + 2, beta, tbeta );
            copy( ( i + 1 ) + 2, tbeta, tgamma );
            for(j=0;j<(i+1)-1;j++)
                talpha[j] /= omega[j], tbeta[j+1] /= omega[j], tgamma[j+1] /= omega[j+1];
            talpha[(i+1)-1] /= omega[(i+1)-1];

            /* Save previous results and update */
            copy( 2 * ( i + 1 ), ev, fv );
            treiglr( i + 1, talpha, tbeta + 1, tgamma + 1, 10000, ev, &ret, &ns );
            inverse_complex_bubble_sort( i + 1, ev );

            /* Compare current to previous eigenvalues of interest */
            vnormdiff( 2 * neig, ev, fv, &r );
            if( vb )
                fprintf( stderr, "%5d: EV Residual = %10.7e\n", i, r );
            if( r < tol )
            {
                *fr = r;
                break;
            }
        }
    }
    *mo = i;

    free( w ); free( q ); free( s ); free( qq ); free( fv );
    free( talpha ); free( tbeta ); free( tgamma );
}

/**
 * Constrained sparse generalized symmetric indefinite Lanczos procedure
 */
void csgsilanczos( int n, int m, long *ia, long *ja, double *A, long *ib, long* jb, double *B, double *alpha, double *beta, double *omega, double *r0, int nv, double *V, double tol, int max, int vb, int *res )
{
    int i,ret;
    double *w,*q,*s,*qq,tau;

    /* Allocate stuff */
    w = (double*) malloc( n * sizeof(double) );
    q = (double*) malloc( n * sizeof(double) );
    s = (double*) malloc( n * sizeof(double) );
    qq = (double*) malloc( n * sizeof(double) );

    /* If input start vector given, use it; otherwise generate random one */
    if( r0 != NULL )
        copy( n, r0, q );
    else
        nrandv( n, q );

    /* Multiply the initial vector */
    ssdgemv( (long) n, ia, ja, A, 1, q, 1, w );
    dotp( n, q, w, &omega[0] );

    /* Start iteration */
    for(i=0;i<m;i++)
    {
        /* Solve Bs = w to get s = inv(B) * w */
        scg( n, ib, jb, B, w, s, max, tol, vb, &ret );

        /* Only conditionally calculate the first sum */
        if( i > 0 )
            vsum( n, 1.0, s, -beta[i] / omega[i-1], qq, s ); /* qq is the previous step's vector */
        dotp( n, w, s, &alpha[i] );
        vsum( n, 1.0, s, -alpha[i] / omega[i], q, s );

        /* Calculate the norm to get omega */
        norm( n, s, &tau );

        /* Stop if tolerance is reached */
        if( fabs( tau ) < tol )
            break;

        /* Take the next step before jumping back to beginning of loop */
        copy( n, q, qq );
        copyvdiv( n, s, tau, q );

        /* Multiply again */
        ssdgemv( (long) n, ia, ja, A, 1, q, 1, w );
        dotp( n, w, q, &omega[i+1] );

        /* Stop if omega vanishes */
        if( fabs( omega[i+1] ) < tol )
            break;
        beta[i+1] = tau * omega[i+1];
    }

    free( w ); free( q ); free( s ); free( qq );
}

/**
 * Basic constrained generalized Lanczos unsymmetric tridiagonalization method with sparse matrix format;
 * store the tridiagonal matrix in band format; biconjugate gradient stabilized method is
 * used to apply the the inverse of B iteratively; this version does additional orthogonalization
 * at each iteration to produce a tridiagonal matrix on a projected subspace of the domain which
 * excludes the vectors given in V
 * @param n Dimension of the matrix
 * @param m Dimension of the Krylov space to build
 * @param ia Row index vector of A
 * @param ja Column indicator of A
 * @param A Floating point values of A
 * @param ib Row index vector of B
 * @param jb Column indicator of B
 * @param B Floating point vector of B
 * @param alpha Output diagonal of the Krylov matrix
 * @param beta Output subdiagonal of the Krylov matrix
 * @param gamma Output superdiagonal of the Krylov matrix
 * @param r0 Initial r given as input; NULL means generate random
 * @param s0 Initial s given as input; NULL means generate random
 * @param nv Number of constraint vectors
 * @param V List of entries of the constraint vectors sequentially
 * @param tol Error tolerance to use
 * @param max Maximum number of iterations to take in the inversion of B via bicgstab
 * @param vb Verbose level; 0 means off, anything else means on
 * @param res Return value indicating success or failure
 */
void csgulanczos( int n, int m, long *ia, long *ja, double *A, long *ib, long* jb, double *B, double *alpha, double *beta, double *gamma, double *r0, double *s0, int nv, double *V, double tol, int max, int vb, int *res )
{
    /* Must apply the inverse of B at each step!!! */
    int i,j,res1,res2,res3;
    long *iat,*jat,*ibt,*jbt;
    double *q,*qq,*p,*pp,*r,*s,*t,*u,*v,*w,*At,*Bt;

    /* Allocate storage for iteration vector */
    q = (double*) malloc( n * sizeof(double) );
    qq = (double*) malloc( n * sizeof(double) );
    p = (double*) malloc( n * sizeof(double) );
    pp = (double*) malloc( n * sizeof(double) );
    r = (double*) malloc( n * sizeof(double) );
    s = (double*) malloc( n * sizeof(double) );
    t = (double*) malloc( n * sizeof(double) );
    w = (double*) malloc( n * sizeof(double) );
    u = (double*) malloc( n * sizeof(double) );
    v = (double*) malloc( n * sizeof(double) );

    /* Generate an initially orthogonal pair of vectors in p, q */
    zerov( n, q );
    zerov( n, p );
    if( r0 != NULL )
        copy( n, r0, r );
    else
        nrandv( n, r );
    if( s0 != NULL )
        copy( n, s0, s );
    else
        nrandv( n, s );
    project( n, r, nv, V );
    project( n, s, nv, V );
    normalize( n, r );
    normalize( n, s );

    /* Transposes of the A and B matrices */
    iat = (long*) malloc( ( n + 1 ) * sizeof(long) );
    jat = (long*) malloc( ia[n] * sizeof(long) );
    At = (double*) malloc( ia[n] * sizeof(double) );
    ibt = (long*) malloc( ( n + 1 ) * sizeof(long) );
    jbt = (long*) malloc( ib[n] * sizeof(long) );
    Bt = (double*) malloc( ib[n] * sizeof(double) );
    stransp( 1, (long) n, (long) n, ia, ja, A, iat, jat, At ); /* There must be a way to get rid of this */
    stransp( 1, (long) n, (long) n, ib, jb, B, ibt, jbt, Bt ); /* At least should not have to store the actual entries */

    /* Loop until residual vector r, s vanish or maximum basis size, m, is reached */
    *res = 0;
    for(i=0;i<m;)
    {
        /* Make sure that r, s do not vanish and are not orthogonal */
        normchk( n, r, CGRAD_ZERO, &res1 );
        normchk( n, s, CGRAD_ZERO, &res2 );
        dotchk( n, r, s, CGRAD_ZERO, &res3 );

        /* Return code which tells which conditions failed */
        if( res1 == -1 || res2 == -1 || res3 == -1 )
        {
            if( res1 == -1 )
                *res |= ( 1 << 0 );
            if( res2 == -1 )
                *res |= ( 1 << 1 );
            if( res3 == -1 )
                *res |= ( 1 << 2 );
            break;
        }

        /* Build the next vectors and coefficients */
#ifdef CGRAD_LANCZOS_NORM
        norm( n, r, &beta[i] );
#else
        dotp( n, s, r, &beta[i] );
        beta[i] = sqrt( fabs( beta[i] ) );
#endif
        dotpdiv( n, s, r, beta[i], &gamma[i] );
        copyvdiv( n, r, beta[i], qq );
        copyvdiv( n, s, gamma[i], pp );

        /* Step forward one */
        ++i;

        /* Everything changes here because we need to multiply by the inverse of B */
        sdgemv( n, n, ia, ja, A, 1, qq, 1, t );
        sbicgstab( n, ib, jb, B, t, u, max, tol, vb, &res1 ); /* Vector uk is in r0 */
        sbicgstab( n, ibt, jbt, Bt, pp, w, max, tol, vb, &res2 ); /* Vector wk is in r0 */
        sdgemv( n, n, iat, jat, At, 1, w, 1, v );

        /* Now u, v are formed so calculate alpha[i], r[i], s[i] */
        dotp( n, pp, u, &alpha[i] );
        vsum( n, 1.0, u, -1.0 * alpha[i], qq, r );
        vsum( n, 1.0, r, -1.0 * gamma[i-1], q, r );
        vsum( n, 1.0, v, -1.0 * alpha[i], pp, s );
        vsum( n, 1.0, s, -1.0 * beta[i-1], p, s );

        /* Final projection */
        project( n, r, nv, V ); /* Find a better way to do this */
        project( n, s, nv, V );

        /* Final copy */
        copy( n, qq, q );
        copy( n, pp, p );
    }

    free( q ); free( qq ); free( p ); free( pp ); free( r ); free( s );
    free( t ); free( w ); free( u ); free( v );
    free( iat ); free( jat ); free( At ); free( ibt ); free( jbt ); free( Bt );
}

/**
 * Basic constrained generalized Lanczos unsymmetric tridiagonalization method with sparse matrix format;
 * store the tridiagonal matrix in band format; biconjugate gradient stabilized method is
 * used to apply the the inverse of B iteratively; this version does additional orthogonalization
 * at each iteration to produce a tridiagonal matrix on a projected subspace of the domain which
 * excludes the vectors given in V
 * @param n Dimension of the matrix
 * @param m Dimension of the Krylov space to build
 * @param ia Row index vector of A
 * @param ja Column indicator of A
 * @param A Floating point values of A
 * @param ib Row index vector of B
 * @param jb Column indicator of B
 * @param B Floating point vector of B
 * @param alpha Output diagonal of the Krylov matrix
 * @param beta Output subdiagonal of the Krylov matrix
 * @param gamma Output superdiagonal of the Krylov matrix
 * @param r0 Initial r given as input; NULL means generate random
 * @param s0 Initial s given as input; NULL means generate random
 * @param nv Number of constraint vectors
 * @param V List of entries of the constraint vectors sequentially
 * @param tol Error tolerance to use
 * @param max Maximum number of iterations to take in the inversion of B via bicgstab
 * @param vb Verbose level; 0 means off, anything else means on
 * @param res Return value indicating success or failure
 */
void csgulanczoslp( int n, int m, long *ia, long *ja, double *A, long *ib, long* jb, double *B, double *alpha, double *beta, double *gamma, double *r0, double *s0, int nv, double *V, double tol, int max, int vb, int *res )
{
    /* Must apply the inverse of B at each step!!! */
    int i,j,res1,res2,res3;
    long *iat,*jat,*ibt,*jbt;
    double *q,*qq,*p,*pp,*r,*s,*t,*u,*v,*w,*At,*Bt;

    /* Allocate storage for iteration vector */
    q = (double*) malloc( n * sizeof(double) );
    qq = (double*) malloc( n * sizeof(double) );
    p = (double*) malloc( n * sizeof(double) );
    pp = (double*) malloc( n * sizeof(double) );
    r = (double*) malloc( n * sizeof(double) );
    s = (double*) malloc( n * sizeof(double) );
    t = (double*) malloc( n * sizeof(double) );
    w = (double*) malloc( n * sizeof(double) );
    u = (double*) malloc( n * sizeof(double) );
    v = (double*) malloc( n * sizeof(double) );

    /* Generate an initially orthogonal pair of vectors in p, q */
    zerov( n, q );
    zerov( n, p );
    if( r0 != NULL )
        copy( n, r0, r );
    else
        nrandv( n, r );
    if( s0 != NULL )
        copy( n, s0, s );
    else
        nrandv( n, s );
    project( n, r, nv, V );
    project( n, s, nv, V );
    normalize( n, r );
    normalize( n, s );

    /* Transposes of the A and B matrices */
    iat = (long*) malloc( ( n + 1 ) * sizeof(long) );
    jat = (long*) malloc( ia[n] * sizeof(long) );
    At = (double*) malloc( ia[n] * sizeof(double) );
    ibt = (long*) malloc( ( n + 1 ) * sizeof(long) );
    jbt = (long*) malloc( ib[n] * sizeof(long) );
    Bt = (double*) malloc( ib[n] * sizeof(double) );
    stransp( 1, (long) n, (long) n, ia, ja, A, iat, jat, At ); /* There must be a way to get rid of this */
    stransp( 1, (long) n, (long) n, ib, jb, B, ibt, jbt, Bt ); /* At least should not have to store the actual entries */

    /* Loop until residual vector r, s vanish or maximum basis size, m, is reached */
    *res = 0;
    for(i=0;i<m;)
    {
        /* Make sure that r, s do not vanish and are not orthogonal */
        normchk( n, r, CGRAD_ZERO, &res1 );
        normchk( n, s, CGRAD_ZERO, &res2 );
        dotchk( n, r, s, CGRAD_ZERO, &res3 );

        /* Return code which tells which conditions failed */
        if( res1 == -1 || res2 == -1 || res3 == -1 )
        {
            if( res1 == -1 )
                *res |= ( 1 << 0 );
            if( res2 == -1 )
                *res |= ( 1 << 1 );
            if( res3 == -1 )
                *res |= ( 1 << 2 );
            break;
        }

        /* Build the next vectors and coefficients */
#ifdef CGRAD_LANCZOS_NORM
        norm( n, r, &beta[i] );
#else
        dotp( n, s, r, &beta[i] );
        beta[i] = sqrt( fabs( beta[i] ) );
#endif
        dotpdiv( n, s, r, beta[i], &gamma[i] );
        copyvdiv( n, r, beta[i], qq );
        copyvdiv( n, s, gamma[i], pp );

        /* Step forward one */
        ++i;

        /* Everything changes here because we need to multiply by the inverse of B */
        sdgemv( n, n, ia, ja, A, 1, qq, 1, t );
        lpsbicgstab( n, ib, jb, B, t, u, max, tol, vb, &res1 ); /* Vector uk is in r0 */
        lpsbicgstab( n, ibt, jbt, Bt, pp, w, max, tol, vb, &res2 ); /* Vector wk is in r0 */
        sdgemv( n, n, iat, jat, At, 1, w, 1, v );

        /* Now u, v are formed so calculate alpha[i], r[i], s[i] */
        dotp( n, pp, u, &alpha[i] );
        vsum( n, 1.0, u, -1.0 * alpha[i], qq, r );
        vsum( n, 1.0, r, -1.0 * gamma[i-1], q, r );
        vsum( n, 1.0, v, -1.0 * alpha[i], pp, s );
        vsum( n, 1.0, s, -1.0 * beta[i-1], p, s );

        /* Final projection */
        project( n, r, nv, V ); /* Find a better way to do this */
        project( n, s, nv, V );

        /* Final copy */
        copy( n, qq, q );
        copy( n, pp, p );
    }

    free( q ); free( qq ); free( p ); free( pp ); free( r ); free( s );
    free( t ); free( w ); free( u ); free( v );
    free( iat ); free( jat ); free( At ); free( ibt ); free( jbt ); free( Bt );
}

/**
 * This function transforms a tridiagonal matrix into its 1-form
 * as defined by thesis by Slemmons
 * @param n Dimension of the square tridiagonal matrix in alpha, beta, gamma
 * @param alpha Diagonal entries of the input matrix
 * @param beta Subdiagonal entries of the input matrix
 * @param gamma Superdiagonal entries of the input matrix
 * @param l Lower triangular part of factorization; the subdiagonal of L
 * @param u Upper triangular part of factorization; the diagonal of U
 */
void trofrm( int n, double *alpha, double *beta, double *gamma, double *l, double *u )
{
    int i;
    double *p;

    /* Allocate storage for diagonal similarity transformation */
    p = (double*) malloc( n * sizeof(double) );

    /* Form the entries of p */
    p[0] = 1.0;
    for(i=1;i<n;i++)
        p[i] = p[i-1] * gamma[i-1];

    /* Form the entries of DAD**-1 */
    for(i=0;i<n-1;i++)
    {
        l[i] = p[i+1] / p[i] * beta[i]; /* l[i] is entry DAD**-1(i+1,i) */
        u[i] = alpha[i];
    }
    u[n-1] = alpha[n-1];
}

/**
 * LU factorization for general nonsymmetric tridiagonal matrices
 * @param n Dimension of the square matrix
 * @param s Shift to use before factoring
 * @param dd Diagonal of the input matrix; length n
 * @param ll Subdiagonal of the input matrix; also make this length n (even though n-1 entries)
 * @param uu Superdiagonal of the input matrix; make this length n as well
 * @param d Output diagonal of U
 * @param l Output subdiagonal of L (diagonal of L is all ones)
 * @param u Output superdiagonal
 */
void trlu( int n, double s, double *dd, double *ll, double *uu, double *d, double *l, double *u )
{
    int i;

    /* Start LU factorization */
    d[0] = dd[0] - s;
    u[0] = uu[0];
    for(i=1;i<n;i++)
    {
        l[i-1] = ll[i-1] / d[i-1];
        d[i] = dd[i] - l[i-1] * uu[i-1] - s;
        u[i] = uu[i]; /* This is redundant, of course */
    }
}

/**
 * Same as trlu except it is assumed that the upper diagonal of the
 * tridiagonal matrix is all ones, i.e. it is assumed to be in its
 * 1-form (see Slemmons thesis "Toward solving the eigenproblem")
 * @param n Size of the square tridiagonal system
 * @param s Shift to use
 * @param dd Diagonal of the input tridiagonal matrix
 * @param ll Subdiagonal of the input tridiagonal matrix
 * @param d Diagonal of U in the LU factorization; superdiagonal is all ones
 * @param l Subdiagonal of L in the LU factorization; diagonal is all ones
 */
void trluof( int n, double s, double *dd, double *ll, double *d, double *l )
{
    int i;

    /* Start LU factorization */
    d[0] = dd[0] - s;
    for(i=1;i<n;i++)
    {
        l[i-1] = ll[i-1] / d[i-1];
        d[i] = dd[i] - l[i-1] - s;
    }
}

/**
 * Reverse multiply step of treig
 * @param n Size of the square tridiagonal matrix
 * @param s Shift to undo in the multiplication
 * @param d Diagonal of U
 * @param l Subdiagonal of L (diagonal of L is unit)
 * @param u Superdiagonal of U
 * @param dd Output diagonal of reverse multiplied matrix
 * @param ll Subdiagonal of reverse multiplied matrix
 * @param uu Superdiagonal of reverse multiplied matrix
 */
void trrm( int n, double s, double *d, double *l, double *u, double *dd, double *ll, double *uu )
{
    int i;

    /* Now have the LU decomposition in l, d, u; do reverse multiply */
    for(i=0;i<n-1;i++)
    {
        ll[i] = d[i+1] * l[i];
        dd[i] = d[i] + u[i] * l[i] + s;
        uu[i] = u[i]; /* This is redundant as well */
    }
    dd[n-1] = d[n-1] + s; /* Don't forget to do the last entry */
}

/**
 * Solves the two-by-two eigenvalue problem for use in shifting for treig
 * @param a Upper left entry in matrix
 * @param b Upper right entry in matrix
 * @param c Lower left entry in matrix
 * @param d Lower right entry in matrix
 * @param re1 Output of the real part of the first eigenvalue
 * @param im1 Output of the imaginary part of the first eigenvalue
 * @param re2 Output of the real part of the second eigenvalue
 * @param im2 Output of the imaginary part of the second eigenvalue
 */
void tteig( double a, double b, double c, double d, double *re1, double *im1, double *re2, double *im2 )
{
    double tr,dt,ds;

    /* If subdiagonal entry i is non-zero then submatrix i on diagonal is conjugate pair */
    tr = a + d;
    dt = a * d - b * c;
    ds = tr * tr - 4.0 * dt;
    *re1 = tr / 2.0;
    *re2 = *re1;
    *im1 = sqrt( fabs( ds ) ) / 2.0;
    *im2 = -(*im1);
    if( ds > 0.0 )
        *re1 = *re1 + *im1, *re2 = *re2 - *im1, *im1 = 0.0, *im2 = 0.0;
}

/**
 * QD algorithm of Rutishauer for tridiagonal matrices
 * @param n Size of the square tridiagonal system
 * @param alpha Diagonal entries
 * @param beta Subdiagonal entries
 * @param gamma Superdiagonal entries
 * @param max Maximum number of iterations to take
 * @param ev Output eigenvalues (real and imaginary parts)
 * @param res Return value; 0 if all eigenvalues found, if > 0, it is the number of eigenvalues not found
 * @param ns Return value; the number of steps taken before convergence
 */
void treiglr( int n, double *alpha, double *beta, double *gamma, int max, double *ev, int *res, int *ns )
{
    int k,m,nn,ret;
    double s,*d,*l,*u,*dd,*ll,*uu,tr,dt,ds,re1,re2,im1,im2;

    /* Allocate temp storage for LU factorzation vectors; try to eliminate this by overwriting alpha, beta, gamma */
    d = (double*) malloc( n * sizeof(double) );
    l = (double*) malloc( n * sizeof(double) );
    u = (double*) malloc( n * sizeof(double) );
    dd = (double*) malloc( n * sizeof(double) );
    ll = (double*) malloc( n * sizeof(double) );
    uu = (double*) malloc( n * sizeof(double) );

    copy( n, beta, ll );
    copy( n, alpha, dd );
    copy( n, gamma, uu );
    ll[n-1] = 0.0;
    uu[n-1] = 0.0;

    /* Calculate the LU decomposition of the tridiagonal system in alpha, beta, gamma and form UL */
    nn = n, k = 0;
    for(m=0;m<max;m++)
    {
        /* Shift strategy of some sort */
        s = dd[nn-1];

        /* Do LU factorization into l, d, u */
        trlu( nn, s, dd, ll, uu, d, l, u );

        /* Now have the LU decomposition in l, d, u; do reverse multiply into ll, dd, uu */
        trrm( nn, s, d, l, u, dd, ll, uu );

        /* Now check the lower righthand corner 2x2 block */
        if( ( nn > 1 && iszero( ll[nn-2] ) ) || nn == 1 ) /* Check in this order */
        {
            /* Then dd[nn-1] is an eigenvalue, so deflate by one */
            ev[2*k+0] = dd[nn-1], ev[2*k+1] = 0.0, ++k;
            nn -= 1;
        }
        else if( ( nn > 2 && iszero( ll[nn-3] ) ) || nn == 2 ) /* Check in this order */
        {
            tteig( dd[nn-2], uu[nn-2], ll[nn-2], dd[nn-1], &re1, &im1, &re2, &im2 );
            ev[2*k+0] = re1, ev[2*k+1] = im1, ++k;
            ev[2*k+0] = re2; ev[2*k+1] = im2, ++k;
            nn -= 2;
        }
        if( nn < 1 )
            break;
    }

    /* Output the reduction */
    *res = nn, *ns = m; /* Zero means all eigenvalues found */

    /* Clean up */
    free( d ); free( l ) ; free( u ); free( dd ); free( ll ); free( uu );
}

/**
 * Tridiagonal eigenvalue solver using quotient difference method with shifts (qds)
 */
void treigqds( int n, double *alpha, double *beta, double *gamma, int max, double *ev, int *res, int *ns )
{
    int i,k,m,nn,ret;
    double s,t,*uu,*ll,tr,dt,ds,re1,re2,im1,im2,avg;

    /* Allocate temp storage for LU factorzation vectors; try to eliminate this by overwriting alpha, beta, gamma */
    uu = (double*) malloc( n * sizeof(double) );
    ll = (double*) malloc( n * sizeof(double) );

#ifdef CGRAD_TREIG_PREPROCESS
    /* Preprocess here by normalizing to the average value of all matrix elements */
    avg = 0.0;
    for(i=0;i<n-1;i++)
        avg += alpha[i] + beta[i] + gamma[i];
    avg += alpha[n-1];
    avg /= (double) ( 3 * n - 2 ); /* The number of elements */
    for(i=0;i<n-1;i++)
        alpha[i] /= avg, beta[i] /= avg, gamma[i] /= avg;
    alpha[n-1] /= avg;
#endif

    /* Reduce alpha, beta, gamma to have a unit upper diagonal for the dqds algorithm */
    trofrm( n, alpha, beta, gamma, ll, uu );

    /* Do the initial LU factorization */
    trluof( n, 0.0, uu, ll, uu, ll );

    /* Start qds iteration */
    nn = n, t = 0.0, k = 0;
    for(m=0;m<max;m++)
    {
        /* Choose a shift */
        s = ll[nn-1] + uu[nn-1];
        t += s;

        /* Update to next qds vectors */
        ll[nn-1] = 0.0; /* This is important */
        uu[0] = uu[0] + ll[0] - s;
        for(i=0;i<nn-1;i++)
        {
            ll[i] = ll[i] * ( uu[i+1] / uu[i] );
            uu[i+1] = uu[i+1] + ll[i+1] - s - ll[i];
        }

        /* Now check the lower righthand corner 2x2 block */
        if( ( nn > 1 && iszero( ll[nn-2] * uu[n-2] ) ) || nn == 1 ) /* Check in this order */
        {
            /* Then dd[nn-1] is an eigenvalue, so deflate by one */
            ev[2*k+0] = uu[nn-1] + ll[nn-2] + t, ev[2*k+1] = 0.0, ++k;
            nn -= 1;
        }
        else if( ( nn > 2 && iszero( ll[nn-3] * uu[nn-2] ) ) || nn == 2 ) /* Check in this order */
        {
            if( nn - 2 == 0 )
                tteig( uu[nn-2], 1.0, ll[nn-2] * uu[nn-2], ll[nn-2] + uu[nn-1], &re1, &im1, &re2, &im2 );
            else
                tteig( ll[nn-3] + uu[nn-2], 1.0, ll[nn-2] * uu[nn-2], ll[nn-2] + uu[nn-1], &re1, &im1, &re2, &im2 );
            ev[2*k+0] = re1 + t, ev[2*k+1] = im1, ++k;
            ev[2*k+0] = re2 + t; ev[2*k+1] = im2, ++k;
            nn -= 2;
        }
        if( nn < 1 )
            break;
    }

#ifdef CGRAD_TREIG_PREPROCESS
    /* Postprocess eigenvalues */
    for(i=0;i<n;i++)
    {
        ev[2*i+0] *= avg;
        ev[2*i+1] *= avg;
    }
#endif

    /* Output the reduction */
    *res = nn, *ns = m; /* Zero means all eigenvalues found */

    /* Free everything */
    free( uu ); free( ll );
}

/**
 * Tridiagonal eigenvalue solver using differential quotient difference method with shifts (qds)
 */
void treigdqds( int n, double *alpha, double *beta, double *gamma, int max, double *ev, int *res, int *ns )
{
    int i,k,m,nn,ret;
    double s,t,*uu,*ll,*d,re1,re2,im1,im2,avg;

    /* Allocate temp storage for LU factorzation vectors; try to eliminate this by overwriting alpha, beta, gamma */
    uu = (double*) malloc( n * sizeof(double) );
    ll = (double*) malloc( n * sizeof(double) );
    d = (double*) malloc( n * sizeof(double) );

#ifdef CGRAD_TREIG_PREPROCESS
    /* Preprocess here by normalizing to the average value of all matrix elements */
    avg = 0.0;
    for(i=0;i<n-1;i++)
        avg += alpha[i] + beta[i] + gamma[i];
    avg += alpha[n-1];
    avg /= (double) ( 3 * n - 2 ); /* The number of elements */
    for(i=0;i<n-1;i++)
        alpha[i] /= avg, beta[i] /= avg, gamma[i] /= avg;
    alpha[n-1] /= avg;
#endif

    /* Reduce alpha, beta, gamma to have a unit upper diagonal for the dqds algorithm */
    trofrm( n, alpha, beta, gamma, ll, uu );

    /* Do the initial LU factorization */
    trluof( n, 0.0, uu, ll, uu, ll );

    /* Start qds iteration */
    nn = n, t = 0.0, k = 0;
    for(m=0;m<max;m++)
    {
        /* Choose a shift */
        s = ll[nn-2] + uu[nn-1];
        t += s;

        /* Update to next qds vectors */
        d[0] = uu[0] - s;
        for(i=0;i<nn-1;i++)
        {
            uu[i] = d[i] + ll[i];
            ll[i] = ll[i] * ( uu[i+1] / uu[i] );
            d[i+1] = d[i] * ( uu[i+1] / uu[i] ) - s;
        }
        uu[nn-1] = d[nn-1];

        /* Now check the lower righthand corner 2x2 block */
        if( ( nn > 1 && iszero( ll[nn-2] * uu[n-2] ) ) || nn == 1 ) /* Check in this order */
        {
            /* Then dd[nn-1] is an eigenvalue, so deflate by one */
            ev[2*k+0] = uu[nn-1] + ll[nn-2] + t, ev[2*k+1] = 0.0, ++k;
            nn -= 1;
        }
        else if( ( nn > 2 && iszero( ll[nn-3] * uu[n-3] ) ) || nn == 2 ) /* Check in this order */
        {
            if( nn - 2 == 0 )
                tteig( uu[nn-2], 1.0, ll[nn-2] * uu[nn-2], ll[nn-2] + uu[nn-1], &re1, &im1, &re2, &im2 );
            else
                tteig( ll[nn-3] + uu[nn-2], 1.0, ll[nn-2] * uu[nn-2], ll[nn-2] + uu[nn-1], &re1, &im1, &re2, &im2 );
            ev[2*k+0] = re1 + t, ev[2*k+1] = im1, ++k;
            ev[2*k+0] = re2 + t; ev[2*k+1] = im2, ++k;
            nn -= 2;
        }
        if( nn < 1 )
            break;
    }

#ifdef CGRAD_TREIG_PREPROCESS
    /* Postprocess eigenvalues */
    for(i=0;i<n;i++)
    {
        ev[2*i+0] *= avg;
        ev[2*i+1] *= avg;
    }
#endif

    /* Output the reduction */
    *res = nn, *ns = m; /* Zero means all eigenvalues found */

    /* Free everything */
    free( uu ); free( ll ); free( d );
}


/**
 * DQDS using implicit double complex conjugate shifts
 */
void treigtridqds( int n, double *alpha, double *beta, double *gamma, int max, double *ev, int *res, int *ns )
{
    int i,k,m,nn,ret;
    double sr,si,s,t,*uu,*ll,*d,re1,re2,im1,im2,avg;
    double xl,yl,xr,yr,zr; /* Bulge variables */

    /* Allocate temp storage for LU factorzation vectors; try to eliminate this by overwriting alpha, beta, gamma */
    uu = (double*) malloc( n * sizeof(double) );
    ll = (double*) malloc( n * sizeof(double) );
    d = (double*) malloc( n * sizeof(double) );

#ifdef CGRAD_TREIG_PREPROCESS
    /* Preprocess here by normalizing to the average value of all matrix elements */
    avg = 0.0;
    for(i=0;i<n-1;i++)
        avg += alpha[i] + beta[i] + gamma[i];
    avg += alpha[n-1];
    avg /= (double) ( 3 * n - 2 ); /* The number of elements */
    for(i=0;i<n-1;i++)
        alpha[i] /= avg, beta[i] /= avg, gamma[i] /= avg;
    alpha[n-1] /= avg;
#endif

    /* Reduce alpha, beta, gamma to have a unit upper diagonal for the dqds algorithm */
    trofrm( n, alpha, beta, gamma, ll, uu );

    /* Do the initial LU factorization */
    trluof( n, 0.0, uu, ll, uu, ll );

    /* Start qds iteration */
    nn = n, t = 0.0, k = 0;
    for(m=0;m<max;m++)
    {
        /* Calculate the shift and decide what to do */
        tteig( ll[nn-2] + uu[nn-2], 1.0, uu[nn-1] * ll[nn-2], uu[nn-1], &re1, &im1, &re2, &im2 );

        /* Set the shift */
        sr = re1; si = im1; /* Shift is not accumulated in tridqds */

        /* Choose which path to take */
        if( iszero( si ) || nn < 4 ) /* If nn < 4 then we can't do tridqds so default to dqds */
        {
            /* Choose a shift */
            s = re1;
            t += s;

            /* Update to next qds vectors */
            d[0] = uu[0] - s;
            for(i=0;i<nn-1;i++)
            {
                uu[i] = d[i] + ll[i];
                ll[i] = ll[i] * ( uu[i+1] / uu[i] );
                d[i+1] = d[i] * ( uu[i+1] / uu[i] ) - s;
            }
            uu[nn-1] = d[nn-1];
        }
        else
        {
            /* First step */
            xr = 1.0; yr = ll[0]; zr = 0.0;
            xr = xr * uu[0] + yr;
            xl = pow( uu[0] + ll[0], 2.0 ) + uu[1] * ll[0] - 2.0 * sr * ( uu[0] + ll[0] ) + sr * sr + si * si;
            yl = - ( uu[1] * ll[0] * uu[2] * ll[1] / xl );
            xl = - uu[1] * ll[0] * ( uu[0] + ll[0] + uu[1] + ll[1] - 2 * sr ) / xl;
            uu[0] = xr - xl;
            xr = yr - xl; yr = zr - yl - xl * ll[1]; zr = - yl * ll[2];
            xr = xr / uu[0]; yr = yr / uu[0]; zr = zr / uu[0];
            ll[0] = xl + yr + xr * uu[1];
            xl = yl + zr + yr * uu[2]; yl = zr * uu[3];
            xr = 1.0 - xr; yr = ll[1] - yr; zr = - zr;

            /* Do central loop */
            for(i=1;i<nn-3;i++)
            {
                xr = xr * uu[i] + yr;
                xl = - xl / ll[i-1]; yl = - yl / ll[i-1];
                uu[i] = xr - xl;
                xr = yr - xl; yr = zr - yl - xl * ll[i+1]; zr = - yl * ll[i+2];
                xr = xr / uu[i]; yr = yr / uu[i]; zr = zr / uu[i];
                ll[i] = xl + yr + xr * uu[i+1];
                xl = yl + zr + yr * uu[i+2]; yl = zr * uu[i+3];
                xr = 1.0 - xr; yr = ll[i+1] - yr; zr = - zr;
            }

            /* Step n - 3 */
            xr = xr * uu[nn-3] + yr;
            xl = - xl / ll[nn-4]; yl = - yl / ll[nn-4];
            uu[nn-3] = xr - xl;
            xr = yr - xl; yr = zr - yl - xl * ll[nn-2];
            xr = xr / uu[nn-3]; yr = yr / uu[nn-3];
            ll[nn-3] = xl + yr + xr * uu[nn-2];
            xl = yl + yr * uu[nn-1];
            xr = 1.0 - xr; yr = ll[nn-2] - yr;

            /* Step n - 2 */
            xr = xr * uu[nn-2] + yr;
            xl = - xl / ll[nn-3];
            uu[nn-2] = xr - xl;
            xr = yr - xl;
            xr = xr / uu[nn-2];
            ll[nn-2] = xl + xr * uu[nn-1];
            xr = 1.0 - xr;

            /* Step n - 1 */
            xr = xr * uu[nn-1];
            uu[nn-1] = xr;
        }

        /* Now check the lower righthand corner 2x2 block of L(i)*U(i) */
        if( ( nn > 1 && iszero( ll[nn-2] * uu[nn-2] ) ) || nn == 1 ) /* Check in this order */
        {
            /* Then dd[nn-1] is an eigenvalue, so deflate by one */
            ev[2*k+0] = uu[nn-1] + ll[nn-2] + t, ev[2*k+1] = 0.0, ++k;
            nn -= 1;
        }
        else if( ( nn > 2 && iszero( ll[nn-3] * uu[nn-3] ) ) || nn == 2 ) /* Check in this order */
        {
            if( nn - 2 == 0 )
                tteig( uu[nn-2], 1.0, ll[nn-2] * uu[nn-2], ll[nn-2] + uu[nn-1], &re1, &im1, &re2, &im2 );
            else
                tteig( ll[nn-3] + uu[nn-2], 1.0, ll[nn-2] * uu[nn-2], ll[nn-2] + uu[nn-1], &re1, &im1, &re2, &im2 );
            ev[2*k+0] = re1 + t, ev[2*k+1] = im1, ++k;
            ev[2*k+0] = re2 + t, ev[2*k+1] = im2, ++k;
            nn -= 2;
        }
        if( nn < 1 )
            break;
    }

#ifdef CGRAD_TREIG_PREPROCESS
    /* Postprocess eigenvalues */
    for(i=0;i<n;i++)
    {
        ev[2*i+0] *= avg;
        ev[2*i+1] *= avg;
    }
#endif

    /* Output the reduction */
    *res = nn, *ns = m; /* Zero means all eigenvalues found */

    /* Free everything */
    free( uu ); free( ll ); free( d );
}

/**
 * Do symmetric inverse iteration for a generalized eigenvalue problem
 * to recover eigenvectors of given eigenvalue
 */
void sgsinvi( int n, long *ia, long *ja, double *A, long *ib, long *jb, double *B, double *x, double mu, double tol, int max, int vb )
{
    int i,j,ret;
    long *ic,*jc;
    double r,*y,*xx,*C;

    /* Allocate stuff */
    y = (double*) malloc( n * sizeof(double) );
    xx = (double*) malloc( n * sizeof(double) );

    /* Build the sparse structure of the matrix to inverse, (A-mu*B) */
    ic = (long*) malloc( ( n + 1 ) * sizeof(long) );
    jc = (long*) malloc( ia[n] * sizeof(long) );
    C = (double*) malloc( ia[n] * sizeof(double) );

    /* Copy */
    for(i=0;i<=n;i++)
        ic[i] = ia[i];
    for(i=0;i<ia[n];i++)
        jc[i] = ja[i], C[i] = A[i] - mu * B[i];

    /* Repeatedly apply inverse transformed matrix to get eigenvector of interest */
    normalize( n, x );
    for(i=0;i<max;i++)
    {
        copy( n, x, xx );
        ssdgemv( n, ib, jb, B, 1, x, 1, y );
        scg( n, ic, jc, C, y, x, max, tol, vb, &ret );
        normalize( n, x );
        r = 0.0;
        for(j=0;j<n;j++)
            r += x[j] * xx[j];
        if( vb )
            fprintf( stderr, "%d: r = %15.7f\n", i, r );
        if( fabs( fabs( r ) - 1.0 ) < tol )
            break;
    }

    /* Clean up */
    free( y );
    free( xx );
    free( ic );
    free( jc );
    free( C );
}

/**
 * Generalized minimum residual method
 */
void gmres( int n, double *A, double *b, double *x, double *H, double *V, int max, double tol, int vb, int *ret )
{
    int i,j,k,m;
    double f,g,beta;
    double *v,*w;

    v = (double*) malloc( n * sizeof(double) );
    w = (double*) malloc( n * sizeof(double) );

    for(i=0;i<n;i++)
    {
        v[i] = b[i];
        for(j=0;j<n;j++)
            v[i] -= A[i*n+j] * x[j];
        f += v[i] * v[i];
    }
    f = sqrt( f );
    for(i=0;i<n;i++)
        v[i] /= f;

    for(i=0;i<max;i++)
    {
        for(j=0;j<n;j++)
        {
            w[j] = 0.0;
            for(k=0;k<n;k++)
                w[j] += A[j*n+k] * v[i*n+k];
        }
        for(j=0;j<=i;j++)
        {
            H[j*n+i] = 0.0;
            for(k=0;k<n;k++)
                H[j*n+i] += w[k] * v[j*n+k];
            for(k=0;k<n;k++)
                w[k] -= H[j*n+i] * v[j*n+k];
        }
        H[(i+1)*n+i] = 0.0;
        for(j=0;j<n;j++)
            H[(i+1)*n+i] += w[j] * w[j];
        H[(i+1)*n+i] = sqrt( H[(i+1)*n+i] );
        if( fabs( H[(i+1)*n+i] ) < tol )
        {
            m = i;
            break;
        }
        for(j=0;j<n;j++)
            v[(i+1)*n+j] = w[j] / H[(i+1)*n+i];
    }

    free( v );
    free( w );
}

// vim: ts=4:sts=4:sw=4:et:sta
