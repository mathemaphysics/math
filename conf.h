// vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab:smarttab
#ifndef VC_CONF_H
#define VC_CONF_H

#include "config.h"

#ifndef MIN
#define MIN(A,B) ( ((A)<(B)) ? (A) : (B) )
#endif
#ifndef MAX
#define MAX(A,B) ( ((A)>(B)) ? (A) : (B) )
#endif

#define VC_FNAME_LEN 1024
#define VC_OUT_FMT_LEN 16
#define VC_MAX_NUM_GRO 32
#define VC_BUFFER_LEN 2048
#define VC_BUFFER_INCR 1024
#define VC_MAX_RANGES 512
#define VC_FRAME_CHUNK_SIZE 256
#define VC_UNWRAP_VEL_WARN
#define VC_GMX_GRO_TITLE_LEN 1024
#define VC_RESNAME_MAX_LEN 32
#define VC_GRPNAME_MAX_LEN 64
#define VC_NDX_LINE_MAX_LEN 64

#if HAVE_DOUBLE == 1  /* This is set in the configure script with --enable-double=yes */
	typedef double t_real;
#else
	typedef float t_real;
#endif

/**
 * Some output functions to simplify
 * printing info/warning/error messages
 * to the screen in an MPI environment
 */
#define mprintmsg(FP,PREFIX,...) mffprintf(0,PREFIX,PACKAGE_STRING,strlen(PACKAGE_STRING), \
                                         rank,size,FP,##__VA_ARGS__)
#define mprintdbg(...) mffprintf(0,"DEBUG",PACKAGE_STRING,strlen(PACKAGE_STRING),\
                                  rank,size,stderr,##__VA_ARGS__)
#define mprinterr(...) mffprintf(0,"ERROR",PACKAGE_STRING,strlen(PACKAGE_STRING),\
                                  rank,size,stderr,##__VA_ARGS__)
typedef struct
{
    int fs; /* Number of bytes per floating point number; only used to confirm format */
    int nt; /* Number of time steps stored */
    int np; /* Number of particles in (every) frame */
    int nf; /* Number of floats to store per time per particle */
} t_binhdr; /* Binary header type */

typedef struct
{
    int n;
} t_cfdata;

#endif

