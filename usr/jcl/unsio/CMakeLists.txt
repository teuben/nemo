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
#INCLUDE(CPack) 
# project name
project (UNSIO)


SET(CMAKE_VERBOSE_MAKEFILE OFF)
# extra path for NEMO
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# load setup flags
include(SetupFlags)

# contains the full path to the top level directory of your build tree
MESSAGE( STATUS "PROJECT_BINARY_DIR: " ${PROJECT_BINARY_DIR} )

# contains the full path to the root of your project source directory,
# i.e. to the nearest directory where CMakeLists.txt contains the PROJECT() command
MESSAGE( STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )

MESSAGE( STATUS "CMAKE_BUILD_TYPE =" ${CMAKE_BUILD_TYPE} )
MESSAGE( STATUS "CMAKE_CXX_FLAGS  =" ${CMAKE_CXX_FLAGS} )
MESSAGE( STATUS "LIBRARY TYPE     =" ${LIBTYPE} )

MESSAGE( STATUS "core OPT  =" ${OPT} )
MESSAGE( STATUS "OSX detected =" ${OSX} )
MESSAGE( STATUS "BUILD_TYPE =" ${RELEASE} )

FIND_PACKAGE(NEMO)
FIND_PACKAGE(SQLITE3)

#set (NEMO_INSTALLED FALSE) # fore NEMO not installed
IF (NOT NEMO_INSTALLED) 
  MESSAGE(STATUS "NEMO LIGHT will be build")

  add_subdirectory(${PROJECT_SOURCE_DIR}/3rdparty/nemolight/build  ${PROJECT_BINARY_DIR}/../3rdparty/nemolight/build )
  set(NEMOLIGHT ${PROJECT_SOURCE_DIR}/3rdparty/nemolight)
  set(NEMOLIB ${NEMOLIGHT}/lib)
  set(NEMOINC ${NEMOLIGHT}/src/inc)

ELSE (NOT NEMO_INSTALLED) 
  set(NEMOLIB $ENV{NEMOLIB})
  set(NEMOINC $ENV{NEMOINC})
ENDIF (NOT NEMO_INSTALLED) 

MESSAGE( STATUS "NEMOLIB : " ${NEMOLIB} )
MESSAGE( STATUS "NEMOINC : " ${NEMOINC} )

# Find all the sources for the LIB
FILE(GLOB LIBSOURCES src/*.cc)

# Destination path for the lib
SET(LIBRARY_OUTPUT_PATH lib)

# compilation options
#dd_definitions(-Wall ${OPT} -Wswitch -Wcast-align -Wreturn-type -Wno-unused-function  -fno-strict-aliasing -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE ${DNOSQLITE3})
add_definitions(${OPT} -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE ${DNOSQLITE3})
SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${WARNCPP} )
SET(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${WARNC} )
SET(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} ${WARNF})

# remove -DNDEBUG in RElease mode to keep activated assert calls
STRING(REPLACE "-DNDEBUG" ""  CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})

# create library "uns"
add_library (unsio ${LIBTYPE} ${LIBSOURCES} )

if(OSX) 
  set_target_properties(unsio PROPERTIES LINK_FLAGS "-L${NEMOLIB} -lnemo -lsqlite3")
endif(OSX)
#
# ------------------
# TEST PROGRAM      
# ------------------
if(DEFINED CMAKE_Fortran_COMPILER AND CMAKE_Fortran_COMPILER MATCHES "^$")
    set(CMAKE_Fortran_COMPILER CMAKE_Fortran_COMPILER-NOTFOUND)
endif(DEFINED CMAKE_Fortran_COMPILER AND CMAKE_Fortran_COMPILER MATCHES "^$")
enable_language(Fortran OPTIONAL)
message(STATUS "CMAKE_Fortran_COMPILER_WORKS = ${CMAKE_Fortran_COMPILER_WORKS}")

# Enable Fortran support
#ENABLE_LANGUAGE(Fortran)
#if (${CMAKE_Fortran_COMPILER} STREQUAL CMAKE_Fortran_COMPILER-NOTFOUND)
#  MESSAGE(STATUS "NO FORTRAN")
#else ()
#  SET(CMAKE_Fortran_COMPILER gfortran)
#endif()

# Executbale output path
set(EXECUTABLE_OUTPUT_PATH bin)
MESSAGE( STATUS "EXECUTABLE_OUTPUT_PATH: " ${EXECUTABLE_OUTPUT_PATH} )

# Make sure the compiler can find include files from our Hello library.
include_directories (${UNSIO_SOURCE_DIR}/src ${NEMOINC} ${NEMOLIB} ${SQLITE3_H_PATH} /usr/include/malloc)

# Make sure the linker can find the Hello library once it is built.
link_directories (${UNSIO_BINARY_DIR}/lib /usr/X11/lib ${SQLITE3_LIB_PATH} ${NEMOLIB} ${FC_GFORT_PATH} ${FC_G77_PATH})

# Find all LIB's sources files
FILE(GLOB exe_sources test_src/*.F)


# ----------------------------------------------------------
# Make sure the compiler can find include files from our Hello library.
#include_directories (${UNS_SOURCE_DIR}/src $ENV{NEMOINC} ${NEMOLIB})

# Find all LIB's sources files
FILE(GLOB execpp_sources test_src/*.cc)

# build cpp executables according to the source
FOREACH(exe_cpp ${execpp_sources})

  get_filename_component(exe ${exe_cpp} NAME_WE)
  MESSAGE( STATUS "New executable ==> " ${exe})
  add_executable (${exe} ${exe_cpp})

  # Link the executable to the Hello library.
  target_link_libraries (${exe} unsio nemo ${SQLITE3_LIB} ${FC_GFORT_LIB} ${FC_G77_LIB} stdc++)

ENDFOREACH(exe_cpp ${execpp_sources})

# build Fortran executables according to the source
# compilation options
if (CMAKE_Fortran_COMPILER)
FOREACH(exe_f ${exe_sources})

  get_filename_component(exe ${exe_f} NAME_WE)
  MESSAGE( STATUS "New executable ==> " ${exe})
  add_executable (${exe} ${exe_f})

  # Link the executable to the Hello library.
  #target_link_libraries (${exe} unsio nemomaing77 nemo g2c sqlite3 stdc++)
  target_link_libraries (${exe} unsio nemomaing77 nemo ${SQLITE3_LIB}  ${FC_LIB} stdc++)

  SET_TARGET_PROPERTIES(${exe} PROPERTIES LINKER_LANGUAGE Fortran)
  # add "-lstdc++"
  #SET(CMAKE_Fortran_LINK_FLAGS "-L${NEMOLIB} ") #-lstdc++")
ENDFOREACH(exe_f ${exe_sources})
endif()


# ----------------------------------------------------------
# Install procedure
# ----------------------------------------------------------
mark_as_advanced(CMAKE_INSTALL_PREFIX)

IF (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  IF (NEMO_INSTALLED)
     set   (CMAKE_INSTALL_PREFIX "$ENV{NEMO}"  CACHE PATH "" FORCE)        #if NEMO exist  => default install into NEMO path
  ELSE()
     set   (CMAKE_INSTALL_PREFIX "$ENV{HOME}/local/unsio" CACHE PATH "" FORCE)  #if NEMO NOT exist => default install into ${HOME}/local/unsio
  ENDIF()
ENDIF ()#NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

MESSAGE( STATUS "--------> CMAKE_INSTALL_PREFIX =" ${CMAKE_INSTALL_PREFIX})

# unsio lib
INSTALL(FILES  ${UNSIO_BINARY_DIR}/lib/libunsio.${LIBEXT} DESTINATION lib)

  
# UNS specific includes
INSTALL(FILES  ${UNSIO_BINARY_DIR}/../src/uns.h DESTINATION inc/uns)
INSTALL(FILES  ${UNSIO_BINARY_DIR}/../src/unsdebug.h DESTINATION inc/uns)
INSTALL(FILES  ${UNSIO_BINARY_DIR}/../src/snapshotinterface.h DESTINATION inc/uns)
INSTALL(FILES  ${UNSIO_BINARY_DIR}/../src/componentrange.h DESTINATION inc/uns)
INSTALL(FILES  ${UNSIO_BINARY_DIR}/../src/ctools.h DESTINATION inc/uns)
INSTALL(FILES  ${UNSIO_BINARY_DIR}/../src/userselection.h DESTINATION inc/uns)

# NEMO specific includes
IF (NOT NEMO_INSTALLED)

  # install nemo lib
  INSTALL(FILES  ${NEMOLIGHT}/lib/libnemo.${LIBEXT} DESTINATION lib)
  INSTALL(FILES  ${NEMOLIGHT}/lib/libnemomaing77.${LIBEXT} DESTINATION lib)

  # Find all NEMOLIGHT header files
  FILE(GLOB nlight_files ${NEMOLIGHT}/src/*.h ${NEMOLIGHT}/src/inc/*.h)
  FOREACH(source_header ${nlight_files})
     INSTALL(FILES  ${source_header} DESTINATION inc)
  ENDFOREACH()

  FILE(GLOB nlight_files ${NEMOLIGHT}/src/inc/snapshot/*.{c,h})
  FOREACH(source_header ${nlight_files})
     INSTALL(FILES  ${source_header} DESTINATION inc/snapshot)
  ENDFOREACH()

ENDIF()
# man
INSTALL(FILES  ${UNSIO_BINARY_DIR}/../man/man3/unsiof.3 DESTINATION man/man3)
# some binaries

#set (CMAKE_INSTALL_PREFIX $ENV{NEMOBIN})
INSTALL(PROGRAMS  ${UNSIO_BINARY_DIR}/bin/uns_info DESTINATION bin)
