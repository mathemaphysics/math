cmake_minimum_required(VERSION 3.0.0)
project(math VERSION 0.1.0)

include(CTest)
enable_testing()

add_library(math lagrange.c quadrature.c linalg.f90 runge.f90 nurbs.c legendre.c mls.c cgrad.c mqrbf.c points.c nuclei.c sparse.c punity.c pbc.c anint.f90)

set(VARIABLE test2)
configure_file(${PROJECT_SOURCE_DIR}/math_config.h.in ${PROJECT_SOURCE_DIR}/math_config.h)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
