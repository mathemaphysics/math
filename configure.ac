#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(libmath, 0.1, ryan-p-daly@uiowa.edu)
AM_INIT_AUTOMAKE(libmath, 0.1)
AC_CONFIG_SRCDIR([lagrange.c])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_FC
AC_PROG_RANLIB
AM_PROG_CC_C_O

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h strings.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memmove pow sqrt])

PKG_CHECK_MODULES( util, util )
PKG_CHECK_MODULES( neural, neural )

# For pbc.c and anint.f90 need to tell it HAVE_DOUBLE always
AC_DEFINE( [HAVE_DOUBLE], [1], [Enable use of double precision] )

AC_CONFIG_FILES([Makefile] [test/Makefile] [test/puksham/Makefile] [math.pc])
AC_OUTPUT

