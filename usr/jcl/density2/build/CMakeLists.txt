# ============================================================================
# Copyright Jean-Charles LAMBERT - 2008                                       
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
 
# project name
project (density)

set(NEMOLIB $ENV{NEMOLIB})

#SET(NEMO_DIR "${PROJECT_BINARY_DIR}" CACHE LOCATION "where was NEMO built?")
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_BINARY_DIR})
MESSAGE( STATUS "CMAKE_MODULE_PATH=" ${CMAKE_MODULE_PATH} )
FIND_PACKAGE(NEMO REQUIRED)

# Dehnen's dir
set (DEHNEN $ENV{NEMO}/usr/dehnen)

#FIND_PACKAGE(NEMO REQUIRED)
MESSAGE( STATUS "NEMOLIB : " ${NEMOLIB} )

# contains the full path to the top level directory of your build tree 
MESSAGE( STATUS "PROJECT_BINARY_DIR: " ${PROJECT_BINARY_DIR} )

# contains the full path to the root of your project source directory,
# i.e. to the nearest directory where CMakeLists.txt contains the PROJECT() command 
MESSAGE( STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )


# compilation options
add_definitions(-O2 -ggdb -DfalcON_NEMO -D_FILE_OFFSET_BITS=64 -mpreferred-stack-boundary=4 -funroll-loops -fforce-addr -rdynamic  -DfalcON_SINGLE)

# Executbale output path
set(EXECUTABLE_OUTPUT_PATH ../bin)
MESSAGE( STATUS "EXECUTABLE_OUTPUT_PATH: " ${EXECUTABLE_OUTPUT_PATH} )

# Make sure the compiler can find include files from our Hello library.
include_directories (${DEHNEN}/falcON/inc ${DEHNEN}/falcON/utils/inc $ENV{NEMOINC} ${NEMOLIB})

# Make sure the linker can find the Hello library once it is built.
link_directories ( ${NEMOLIB} ${DEHNEN}/falcON/lib ${DEHNEN}/utils/lib )


# ----------------------------------------------------------
# Make sure the compiler can find include files from our Hello library.
#include_directories (${UNS_SOURCE_DIR}/src $ENV{NEMOINC} ${NEMOLIB})

# Find all LIB's sources files
FILE(GLOB execpp_sources ../src/*.cc)

# build cpp executables according to the source
FOREACH(exe_cpp ${execpp_sources})

  get_filename_component(exe ${exe_cpp} NAME_WE)
  MESSAGE( STATUS "New executable ==> " ${exe})
  add_executable (${exe} ${exe_cpp})

  # Link the executable to the Hello library.
  target_link_libraries (${exe}  WDutils falcON nemo dl stdc++ )

ENDFOREACH(exe_cpp ${execpp_sources})

# ----------------------------------------------------------
# install target                                            

set (CMAKE_INSTALL_PREFIX $ENV{NEMO})

# 
# install bin targets
FOREACH(exe_cpp ${execpp_sources})

  get_filename_component(exe ${exe_cpp} NAME_WE)
  INSTALL(PROGRAMS  ${PROJECT_BINARY_DIR}/../bin/${exe} DESTINATION bin)

ENDFOREACH(exe_cpp ${execpp_sources})


