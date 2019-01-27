#ifndef POINTS_H
#define POINTS_H

#define DOMAIN_CUBE 0
#define DOMAIN_BOX 1
#define DOMAIN_SPHERE 2
#define DOMAIN_ELLIPSOID 3

typedef struct
{
	int dim;
	int type;
	double *orig;
	double *params;
} shape_t;

void point_print( FILE *, int, double * );

void domain_indicator( double *, shape_t *, int * );

int sphere_intersection( int, double *, double, double *, double, double *, double *, double * );

void lens_gauss_point( int, double *, double, double *, double, double , double *, long *, double *, double *, double *, double * );

double lens_volume( int, double, double, double );

void sphere_gauss_point( int, double *, double, double *, long *, double *, double *, double *, double * );

int load_points( char *, int *, int *, double **, int );

int write_points( char *, int, int, double *, int );

int generate_cloud_macqueen( int, int, int, int (*)(double*,void*), double (*)(double*,void*), double *, double **, void *, void * );

int smooth_cloud_macqueen( int, int, int, int (*)(double*,void*), double (*)(double*,void*), double *, double *, void *, void * );

int smooth_cloud_bubble_mc( int, int, double *, int (*)(double*,void*), double (*)(double*,void*), double, double *, double *, double, double, int, int, void *, void * );

int generate_cloud_supports( int, int, double *, int, double, double * );

int generate_cloud_supports_min( int, int, double *, int, double, double *, double );

/**
 * Data structure for kd trees for use in nearest neighbor
 * set finding algorithms
 */
struct kdnode_s
{
        int dim;
        int depth;
	int npts;
	int idx;
	double split;
        double *pos;
        struct kdnode_s *left;
        struct kdnode_s *right;
};

typedef struct kdnode_s kdnode_t;

int kdtree_build_median( kdnode_t *, int, int, double *, int );

int kdtree_build_average( kdnode_t *, int, int, double *, int, int );

int kdtree_range_query( kdnode_t *, double *, double, int **, int *, int * );

void kdtree_walk( kdnode_t * );

#endif

