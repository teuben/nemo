# ============================================================================
# Copyright Jean-Charles LAMBERT - 2009-2012                                       
# e-mail:   Jean-Charles.Lambert@oamp.fr
# address:  Dynamique des galaxies
#           Centre de donneeS Astrophysique de Marseille (CeSAM)
#           Laboratoire d'Astrophysique de Marseille                          
#           Pole de l'Etoile, site de Chateau-Gombert                         
#           38, rue Frederic Joliot-Curie                                     
#           13388 Marseille cedex 13 France                                   
#           CNRS U.M.R 6110
# ============================================================================
# CMakeListst.txt file for nemo library                                        
# ============================================================================
cmake_minimum_required(VERSION 2.6)
 
# project name
project (nemo)

# set CMAKE Module PATH (NEMO)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_BINARY_DIR})
MESSAGE( STATUS "CMAKE_MODULE_PATH=" ${CMAKE_MODULE_PATH} )


# contains the full path to the top level directory of your build tree 
MESSAGE( STATUS "PROJECT_BINARY_DIR: " ${PROJECT_BINARY_DIR} )
MESSAGE( STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )

MESSAGE( STATUS "CMAKE_BUILD_TYPE =" ${CMAKE_BUILD_TYPE} )
MESSAGE( STATUS "CMAKE_CXX_FLAGS  =" ${CMAKE_CXX_FLAGS} )
#MESSAGE( STATUS "LIBRARY TYPE     =" ${LIBTYPE} )

MESSAGE( STATUS "nemolight OPT    =" ${OPT} )


# Lib : libnemo.so
# Find all the sources for the utils LIB
FILE(GLOB LIBUTILS src/*.c )
# create the library "libnemo.so"
#add_library (nemo ${LIBTYPE} ${LIBUTILS})
add_library (nemo STATIC ${LIBUTILS})

# Lib : libnemomaing77.so
# Find all the sources for old fortran compiler
FILE(GLOB LIBUTILS src/fortran_old_gcc/*.c )
# create the library "libnemomaing77.so"
#add_library (nemomaing77 ${LIBTYPE} ${LIBUTILS})
add_library (nemomaing77 STATIC ${LIBUTILS})

if(OSX)
   set_target_properties(nemo PROPERTIES LINK_FLAGS "-undefined suppress -flat_namespace")
   set_target_properties(nemomaing77 PROPERTIES LINK_FLAGS "-undefined suppress -flat_namespace")
endif(OSX)
# Destination path for the lib
SET(LIBRARY_OUTPUT_PATH lib)
# compilation options
#add_definitions(-Wall ${OPT} -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE)
#add_definitions(${OPT} ${WARNOPT}  -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE)
add_definitions(${OPT} -fPIC -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE)
SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${WARNCPP} )
SET(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${WARNC} )

if (CMAKE_Fortran_COMPILER)
  SET(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} ${WARNF})
endif()

# remove -DNDEBUG in RElease mode to keep activated assert calls
STRING(REPLACE "-DNDEBUG" ""  CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})

# SET NEMOPATH
#SET(NEMOLIB  ${PROJECT_SOURCE_DIR}/../lib)
#SET(NEMOINC  ${PROJECT_SOURCE_DIR}/src/inc)

# Make sure the compiler can find include files from our Hello library.
include_directories (${PROJECT_SOURCE_DIR}/src/inc /usr/include/malloc)

# 
# install bin targets
#
