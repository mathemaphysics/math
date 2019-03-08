#ifndef VC_GEOM_H
#define VC_GEOM_H

#include "config.h"
#include "conf.h"

#ifdef USE_GMXLIB_5
#include "gromacs/legacyheaders/tpxio.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/pbc.h"
#include "gromacs/legacyheaders/vec.h"
#else
#include "gromacs/tpxio.h"
#include "gromacs/txtdump.h"
#include "gromacs/pbc.h"
#include "gromacs/vec.h"
#endif

t_real vnorm_sqr( t_real * );

t_real vnorm( t_real * );

void normalize( t_real * );

t_real dist_real3( t_real *, t_real * );

t_real dist_real3_sqr( t_real *, t_real * );

void vavg_real3( t_real *, t_real *, t_real * );

void pbc_vec_real3( t_real *, t_real *, t_real *, t_real * );

#ifdef USE_GMX_PBC_DIST

void pbc_init_box_gmx( char *, int, int *, t_inputrec *, t_topology *, rvec *, t_pbc * );

t_real pbc_dist_real3_gmx( t_real *, t_real *, rvec *, t_pbc * );

t_real pbc_dist_real3_sqr_gmx( t_real *, t_real *, rvec *, t_pbc * );

void pbc_vec_real3_gmx( t_real *, t_real *, t_real *, t_pbc * );

#endif

t_real pbc_dist_real3( t_real *, t_real *, t_real * );

t_real pbc_dist_real3_sqr( t_real *, t_real *, t_real * );

t_real dotp_real3( t_real *, t_real * );

#endif

