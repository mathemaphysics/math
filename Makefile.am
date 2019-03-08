FC=gfortran
SUBDIRS=test
lib_LIBRARIES = libmath.a
include_HEADERS = lagrange.h linalg.h quadrature.h nurbs.h runge.h legendre.h mls.h mqrbf.h points.h nuclei.h sparse.h punity.h cgrad.h wavelet.h geom.h
libmath_a_SOURCES = lagrange.c quadrature.c linalg.f90 runge.f90 nurbs.c legendre.c mls.c cgrad.c mqrbf.c points.c nuclei.c sparse.c punity.c wavelet.c pbc.c anint.f90
libmath_a_CFLAGS = -I./. $(util_CFLAGS) $(libgmx_d_CFLAGS)
