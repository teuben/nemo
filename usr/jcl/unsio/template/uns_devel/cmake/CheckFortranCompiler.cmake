# -*-cmake-*-
# check if fortran compiler is installed
if(DEFINED CMAKE_Fortran_COMPILER AND CMAKE_Fortran_COMPILER MATCHES "^$")
  set(CMAKE_Fortran_COMPILER CMAKE_Fortran_COMPILER-NOTFOUND)
endif(DEFINED CMAKE_Fortran_COMPILER AND CMAKE_Fortran_COMPILER MATCHES "^$")
enable_language(Fortran OPTIONAL)
message(STATUS "CMAKE_Fortran_COMPILER_WORKS = ${CMAKE_Fortran_COMPILER_WORKS}")

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME STREQUAL "gfortran")
  # gfortran
  set (EXTRA_Fortran_FLAGS " ${WARNF} -ffixed-line-length-none ")
endif()