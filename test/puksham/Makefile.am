bin_PROGRAMS = puksham
include_HEADERS = puksvar.h pukscfg.h funcs.h
puksham_SOURCES = puksham.c mbuild.c funcs.c ../../points.c ../../punity.c ../../cgrad.c ../../nuclei.c ../../pbc.c ../../anint.f90
puksham_CFLAGS = -I../../include -I./ $(util_CFLAGS) $(libgmx_d_CFLAGS)
puksham_LDFLAGS = -lm
puksham_LDADD = $(util_LIBS)
