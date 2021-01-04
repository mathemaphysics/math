cmake_minimum_required(VERSION 3.0.0)
project(math VERSION 0.1.0)

include(CTest)
enable_testing()

enable_language(Fortran)
set(FC gfortran)

set(GROMACS_VERSION "2018")
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f90)

add_library(mathext linalg.f90 runge.f90 anint.f90)
add_library(math lagrange.c quadrature.c nurbs.c legendre.c mls.c cgrad.c mqrbf.c points.c nuclei.c sparse.c punity.c pbc.c anint.f90 combinadic.c parse.c trees.c symtab.c cfgread.c)
target_link_libraries(math gromacs mathext)
target_include_directories(math PUBLIC /usr/include)

configure_file(${PROJECT_SOURCE_DIR}/math_config.h.in ${PROJECT_SOURCE_DIR}/math_config.h)

add_subdirectory(test)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
