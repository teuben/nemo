# ============================================================================
# Copyright Jean-Charles LAMBERT - 2009                                       
# e-mail:   Jean-Charles.Lambert@oamp.fr
# address:  Dynamique des galaxies                                            
#           Laboratoire d'Astrophysique de Marseille                          
#           Pole de l'Etoile, site de Chateau-Gombert                         
#           38, rue Frederic Joliot-Curie                                     
#           13388 Marseille cedex 13 France                                   
#           CNRS U.M.R 6110
# ============================================================================
# CMakeListst.txt file for UNS library                                        
# ============================================================================
cmake_minimum_required(VERSION 2.6)
#INCLUDE(CPack) 
# project name
project (UNSIO)
SET(CMAKE_VERBOSE_MAKEFILE OFF)
# extra path for NEMO
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_BINARY_DIR})

FIND_PACKAGE(NEMO)

IF (NOT NEMO_INSTALLED) 

  MESSAGE(FATAL_ERROR "NEMO environment not loaded !!!")

ELSE (NOT NEMO_INSTALLED) 

  set(NEMOLIB $ENV{NEMOLIB})
  set(NEMOINC $ENV{NEMOINC})

ENDIF (NOT NEMO_INSTALLED) 




MESSAGE( STATUS "NEMOLIB : " ${NEMOLIB} )

#execute_process(COMMAND  "grep 'FC  =' ${NEMOLIB}/makedefs | cut -d= -f2"
#	OUTPUT_VARIABLE FF OUTPUT_STRIP_TRAILING_WHITESPACE)
#MESSAGE( STATUS "FORTRAN: " ${FF} )
# contains the full path to the top level directory of your build tree 
MESSAGE( STATUS "PROJECT_BINARY_DIR: " ${PROJECT_BINARY_DIR} )

# contains the full path to the root of your project source directory,
# i.e. to the nearest directory where CMakeLists.txt contains the PROJECT() command 
MESSAGE( STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )

# Find all the sources for the LIB
FILE(GLOB LIBSOURCES ../src/*.cc)

# Destination path for the lib
SET(LIBRARY_OUTPUT_PATH ../lib)

# compilation options
add_definitions(-Wall -g -Wswitch -Wcast-align -Wreturn-type -Wno-unused-function  -fno-strict-aliasing -ggdb)

# create library "uns"
add_library (unsio SHARED ${LIBSOURCES})
#
# ------------------
# TEST PROGRAM      
# ------------------
# Enable Fortran support
ENABLE_LANGUAGE(Fortran)
SET(CMAKE_Fortran_COMPILER gfortran)


# Executbale output path
set(EXECUTABLE_OUTPUT_PATH ../bin)
MESSAGE( STATUS "EXECUTABLE_OUTPUT_PATH: " ${EXECUTABLE_OUTPUT_PATH} )

# Make sure the compiler can find include files from our Hello library.
include_directories (${UNSIO_SOURCE_DIR}/../src $ENV{NEMOINC} ${NEMOLIB})

# Make sure the linker can find the Hello library once it is built.
link_directories (${UNSIO_BINARY_DIR}/lib ${NEMOLIB} ${FC_LIB_PATH})

# Find all LIB's sources files
FILE(GLOB exe_sources ../test_src/*.F)


# ----------------------------------------------------------
# Make sure the compiler can find include files from our Hello library.
#include_directories (${UNS_SOURCE_DIR}/src $ENV{NEMOINC} ${NEMOLIB})

# Find all LIB's sources files
FILE(GLOB execpp_sources ../test_src/*.cc)

# build cpp executables according to the source
FOREACH(exe_cpp ${execpp_sources})

  get_filename_component(exe ${exe_cpp} NAME_WE)
  MESSAGE( STATUS "New executable ==> " ${exe})
  add_executable (${exe} ${exe_cpp})

  # Link the executable to the Hello library.
  target_link_libraries (${exe} unsio nemo sqlite3  X11 pgplot ${FC_LIB} stdc++)

ENDFOREACH(exe_cpp ${execpp_sources})

# build Fortran executables according to the source
# compilation options
FOREACH(exe_f ${exe_sources})

  get_filename_component(exe ${exe_f} NAME_WE)
  MESSAGE( STATUS "New executable ==> " ${exe})
  add_executable (${exe} ${exe_f})

  # Link the executable to the Hello library.
  #target_link_libraries (${exe} unsio nemomaing77 nemo g2c sqlite3 stdc++)
  target_link_libraries (${exe} unsio nemomaing77 nemo sqlite3 stdc++)

  SET_TARGET_PROPERTIES(${exe} PROPERTIES LINKER_LANGUAGE Fortran)
  # add "-lstdc++"
  #SET(CMAKE_Fortran_LINK_FLAGS "-L${NEMOLIB} ") #-lstdc++")
ENDFOREACH(exe_f ${exe_sources})

# ----------------------------------------------------------
# install target                                            

IF (NEMO_INSTALLED) 
  set (CMAKE_INSTALL_PREFIX $ENV{NEMO})

  # lib
  INSTALL(FILES  ${UNSIO_BINARY_DIR}/../lib/libunsio.so DESTINATION lib)
  # includes
  INSTALL(FILES  ${UNSIO_BINARY_DIR}/../src/uns.h DESTINATION inc/uns)
  INSTALL(FILES  ${UNSIO_BINARY_DIR}/../src/unsdebug.h DESTINATION inc/uns)
  INSTALL(FILES  ${UNSIO_BINARY_DIR}/../src/snapshotinterface.h DESTINATION inc/uns)
  INSTALL(FILES  ${UNSIO_BINARY_DIR}/../src/componentrange.h DESTINATION inc/uns)
  INSTALL(FILES  ${UNSIO_BINARY_DIR}/../src/ctools.h DESTINATION inc/uns)
  INSTALL(FILES  ${UNSIO_BINARY_DIR}/../src/userselection.h DESTINATION inc/uns)
ENDIF (NEMO_INSTALLED) 
# man
#INSTALL(FILES  ${UNSIO_BINARY_DIR}/../doc/man/man3/uns.3 DESTINATION man/man3)

#set (CMAKE_INSTALL_PREFIX $ENV{NEMOBIN})
#INSTALL(PROGRAMS  ${UNSIO_BINARY_DIR}/../bin/unscpp DESTINATION bin)
