
# check if fortran compiler is installed
if(DEFINED CMAKE_Fortran_COMPILER AND CMAKE_Fortran_COMPILER MATCHES "^$")
    set(CMAKE_Fortran_COMPILER CMAKE_Fortran_COMPILER-NOTFOUND)
endif(DEFINED CMAKE_Fortran_COMPILER AND CMAKE_Fortran_COMPILER MATCHES "^$")
enable_language(Fortran OPTIONAL)
message(STATUS "CMAKE_Fortran_COMPILER_WORKS = ${CMAKE_Fortran_COMPILER_WORKS}")
