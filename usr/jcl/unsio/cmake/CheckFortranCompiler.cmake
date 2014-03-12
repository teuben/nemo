# -*-cmake-*-
# ============================================================================
# Copyright Jean-Charles LAMBERT - 2008-2014
# e-mail:   Jean-Charles.Lambert@oamp.fr
# address:  Dynamique des galaxies
#           Centre de donneeS Astrophysique de Marseille (CeSAM)
#           Laboratoire d'Astrophysique de Marseille
#           Pole de l'Etoile, site de Chateau-Gombert
#           38, rue Frederic Joliot-Curie
#           13388 Marseille cedex 13 France
#           CNRS U.M.R 6110
# ============================================================================
# check if fortran compiler is installed
# ============================================================================

if ( NOT CMAKE_DISABLE_FORTRAN )
    MESSAGE("\n\nTrying do detect fortran support.....")
    MESSAGE("\nYou can change FORTRAN compiler name by setting FC variable before running cmake, like\nexport FC=ifort\n")
    MESSAGE("or run directly the command:\nFC=ifort cmake ..\n\n")
    MESSAGE("If you want to disable FORTRAN support, add following command to cmake:\n-DCMAKE_DISABLE_FORTRAN=1\n\n")
    enable_language(Fortran OPTIONAL)
else()
    MESSAGE("\n\nFortran support DISABLE\n\n")
endif()

if(DEFINED CMAKE_Fortran_COMPILER AND CMAKE_Fortran_COMPILER MATCHES "^$")
  set(CMAKE_Fortran_COMPILER CMAKE_Fortran_COMPILER-NOTFOUND)
endif(DEFINED CMAKE_Fortran_COMPILER AND CMAKE_Fortran_COMPILER MATCHES "^$")

if (CMAKE_Fortran_COMPILER)
  message(STATUS "CMAKE_Fortran_COMPILER_WORKS = ${CMAKE_Fortran_COMPILER_WORKS}")
  message(STATUS "Fortran compiler : " ${CMAKE_Fortran_COMPILER})

  # FFLAGS depend on the compiler
  get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
  
  if (Fortran_COMPILER_NAME STREQUAL "gfortran")
    # gfortran
    set (EXTRA_Fortran_FLAGS " ${WARNF} -ffixed-line-length-none ")
  endif()
endif()