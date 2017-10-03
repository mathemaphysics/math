FC=gfortran
SUBDIRS=test
lib_LIBRARIES = libmath.a
include_HEADERS = include/lagrange.h include/linalg.h include/quadrature.h include/nurbs.h include/runge.h include/legendre.h include/mls.h include/mqrbf.h include/points.h include/nuclei.h include/sparse.h include/punity.h include/cgrad.h include/wavelet.h
libmath_a_SOURCES = lagrange.c quadrature.c linalg.f90 runge.f90 nurbs.c legendre.c mls.c cgrad.c mqrbf.c points.c nuclei.c sparse.c punity.c wavelet.c pbc.c anint.f90
libmath_a_CFLAGS = -I./include $(util_CFLAGS)
