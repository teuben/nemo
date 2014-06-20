# ============================================================================
# Copyright Jean-Charles LAMBERT - 2010-2012
# e-mail:   Jean-Charles.Lambert@oamp.fr                                      
# address:  Dynamique des galaxies
#           Centre de donneeS Astrophysique de Marseille (CeSAM)
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
project (uns_projects)

# extra path for NEMO
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH} )
MESSAGE( STATUS "CMAKE_MODULE_PATH=" ${CMAKE_MODULE_PATH} )

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

# require NEMO
FIND_PACKAGE(NEMO REQUIRED)
# require UNSIO
FIND_PACKAGE(UNSIO REQUIRED)
# Sqlite3
FIND_PACKAGE(SQLITE3)
# Boost (deprecated)
FIND_PACKAGE(BOOST)

#FIND_PACKAGE(NEMO REQUIRED)
MESSAGE( STATUS "NEMOLIB : " ${NEMOLIB} )

# contains the full path to the top level directory of your build tree 
MESSAGE( STATUS "PROJECT_BINARY_DIR: " ${PROJECT_BINARY_DIR} )
#MESSAGE ( STATUS "Boost_LIBRARY_DIRS : " ${Boost_LIBRARY_DIRS}) 
# contains the full path to the root of your project source directory,
# i.e. to the nearest directory where CMakeLists.txt contains the PROJECT() command 
MESSAGE( STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )


# Find all the sources for the utils LIB

if ( NEMO_INSTALLED )
   FILE(GLOB LIBUTILS ${PROJECT_SOURCE_DIR}/lib/utils/*.cc ${PROJECT_SOURCE_DIR}/lib/utils/*.c ${PROJECT_SOURCE_DIR}/lib/utils/nemodep/*.cc)
else ()
   FILE(GLOB LIBUTILS ${PROJECT_SOURCE_DIR}/lib/utils/*.cc ${PROJECT_SOURCE_DIR}/lib/utils/*.c)
endif ()
# create the library "JCLutils"
add_library (JCLutils SHARED ${LIBUTILS})
if(OSX)
set_target_properties(JCLutils PROPERTIES LINK_FLAGS "-undefined suppress -flat_namespace")
endif(OSX)
# Find all the sources for the projects LIB
if ( NEMO_INSTALLED )
   FILE(GLOB LIBPROJECTS ${PROJECT_SOURCE_DIR}/lib/projects/*.cc ${PROJECT_SOURCE_DIR}/lib/projects/nemodep/*.cc)
else ()
   FILE(GLOB LIBPROJECTS ${PROJECT_SOURCE_DIR}/lib/projects/*.cc)
endif ()
# create the library "JCLutils"
add_library (JCLprojects SHARED ${LIBPROJECTS})
if(OSX)
set_target_properties(JCLprojects PROPERTIES LINK_FLAGS "-undefined suppress -flat_namespace")
endif(OSX)
# Destination path for the lib
SET(LIBRARY_OUTPUT_PATH lib)


# compilation options
#add_definitions(${OPT} -DNO_CUDA -DfalcON_NEMO -DfalcON_SINGLE  -funroll-loops -Wall ${DNOBOOST})
add_definitions(${OPT} -DNO_CUDA -DfalcON_NEMO -DfalcON_SINGLE  -Wall ${DNOBOOST})
SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${WARNCPP} )
SET(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${WARNC} )
SET(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} ${WARNF})

# remove -DNDEBUG in RElease mode to keep activated assert calls
STRING(REPLACE "-DNDEBUG" ""  CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
message (STATUS "CMAKE_CXX_FLAGS_RELEASE = " ${CMAKE_CXX_FLAGS_RELEASE}  )

# Executbale output path
set(EXECUTABLE_OUTPUT_PATH bin)
MESSAGE( STATUS "EXECUTABLE_OUTPUT_PATH: " ${EXECUTABLE_OUTPUT_PATH} )

# Make sure the compiler can find include files from our Hello library.
include_directories (${PROJECT_SOURCE_DIR}/src )
include_directories (${PROJECT_SOURCE_DIR}/lib/utils ${PROJECT_SOURCE_DIR}/lib/utils/nemodep )
include_directories (${PROJECT_SOURCE_DIR}/lib/projects ${PROJECT_SOURCE_DIR}/lib/projects/nemodep)
if ( FALCON_INSTALLED )
  include_directories (${DEHNEN}/falcON/utils/inc ${DEHNEN}/falcON/inc)
endif()
include_directories (${UNSIOPATH}/lib ${UNSIOPATH}/inc ${UNSIOPATH}/inc/uns)
include_directories (/usr/include/malloc)


# Make sure the linker can find the Hello library once it is built.
link_directories ( ${UNSIOPATH}/lib ${DEHNEN}/falcON/lib ${DEHNEN}/utils/lib $ENV{PGPLOT_DIR}  ${G2C_DIR} ${PROJECT_BINARY_DIR}/lib /usr/lib64 /usr/X11/lib ${FC_GFORT_PATH} ${FC_G77_PATH})

# ----------------------------------------------------------
# Install SETUP
# ----------------------------------------------------------

mark_as_advanced(CMAKE_INSTALL_PREFIX)

set (CMAKE_INSTALL_PREFIX ${UNSIOPATH})
MESSAGE( STATUS "--------> CMAKE_INSTALL_PREFIX =" ${CMAKE_INSTALL_PREFIX})

if (OSX)
  MESSAGE(STATUS "\nAfter running \"make install\", do not forget to set your environment variable DYLD_LIBRARY_PATH using command:\n\n export DYLD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib\n")
endif(OSX)

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

# ----------------------------------------------------------
# Make sure the compiler can find include files from our Hello library.
#include_directories (${UNS_SOURCE_DIR}/src $ENV{NEMOINC} ${NEMOLIB})

# Find all exes sources files

if ( NEMO_INSTALLED )
   FILE(GLOB execpp_sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/src/nemodep/*.cc)
else ()
   FILE(GLOB execpp_sources ${PROJECT_SOURCE_DIR}/src/*.cc)
endif ()
# build cpp executables according to the source
FOREACH(exe_cpp ${execpp_sources})

  get_filename_component(exe ${exe_cpp} NAME_WE)
  MESSAGE( STATUS "New executable ==> " ${exe})
  add_executable (${exe} ${exe_cpp})

 #
  if ( NEMO_INSTALLED )
     target_link_libraries (${exe}  JCLutils JCLprojects cpgplot pgplot.a ${FC_GFORT_LIB} ${FC_G77_LIB}  X11 unsio WDutils.gcc falcON nemo ${SQLITE3_LIB} gomp ${LIBBOOST} pthread dl)
  else ()
     target_link_libraries (${exe}  JCLutils JCLprojects  ${FC_GFORT_LIB} ${FC_G77_LIB} unsio nemo ${SQLITE3_LIB} gomp ${LIBBOOST} dl)
  endif ()
  INSTALL(TARGETS ${exe} RUNTIME  DESTINATION bin)

ENDFOREACH(exe_cpp ${execpp_sources})

# ----------------------------------------------------------
# install target                                            



IF   (OSX) # Mac OSX
  SET(SOEXT "dylib")
ELSE (OSX) # Linux
  SET(SOEXT "so")
ENDIF(OSX)
# 
# install bin targets
FOREACH(exe_cpp ${execpp_sources})

  get_filename_component(exe ${exe_cpp} NAME_WE)
  #INSTALL(PROGRAMS  ${PROJECT_BINARY_DIR}/bin/${exe}     DESTINATION bin)
ENDFOREACH(exe_cpp ${execpp_sources})
INSTALL(FILES  ${PROJECT_BINARY_DIR}/lib/libJCLutils.${SOEXT} DESTINATION lib)
INSTALL(FILES  ${PROJECT_BINARY_DIR}/lib/libJCLprojects.${SOEXT} DESTINATION lib)
INSTALL(FILES  ${PROJECT_SOURCE_DIR}/man/uns_2dplot.1 DESTINATION man/man1)
INSTALL(FILES  ${PROJECT_SOURCE_DIR}/man/uns_stack.1 DESTINATION man/man1)
INSTALL(FILES  ${PROJECT_SOURCE_DIR}/man/uns2uns.1 DESTINATION man/man1)
INSTALL(FILES  ${PROJECT_SOURCE_DIR}/man/uns_projects.1 DESTINATION man/man1)
INSTALL(FILES  ${PROJECT_SOURCE_DIR}/lib/utils/nemodep/cfalcon.h DESTINATION inc/uns)
INSTALL(FILES  ${PROJECT_SOURCE_DIR}/lib/utils/csnaptools.h DESTINATION inc/uns)
INSTALL(FILES  ${PROJECT_SOURCE_DIR}/lib/utils/ctimer.h DESTINATION inc/uns)
#
