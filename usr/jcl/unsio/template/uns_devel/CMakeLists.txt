# -*-cmake-*-
# ============================================================================
# Copyright Jean-Charles LAMBERT - 2008-2014
#           Centre de donneeS Astrophysiques de Marseille (CeSAM)       
# e-mail:   Jean-Charles.Lambert@lam.fr                                      
# address:  Aix Marseille Universite, CNRS, LAM 
#           Laboratoire d'Astrophysique de Marseille                          
#           Pole de l'Etoile, site de Chateau-Gombert                         
#           38, rue Frederic Joliot-Curie                                     
#           13388 Marseille cedex 13 France                                   
#           CNRS UMR 7326    
# ============================================================================
# CMakeListst.txt file for a developpement project based in UNSIO api
# ============================================================================

cmake_minimum_required(VERSION 2.6)

# project name
project (uns_devel)

# extra path for NEMO
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH} )
MESSAGE( STATUS "CMAKE_MODULE_PATH=" ${CMAKE_MODULE_PATH} )

# load setup flags
include(SetupFlags)

# require NEMO
FIND_PACKAGE(NEMO REQUIRED)
# require UNSIO
FIND_PACKAGE(UNSIO REQUIRED)
# Sqlite3
FIND_PACKAGE(SQLITE3)

# Check fortran compiler
include(CheckFortranCompiler)

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
MESSAGE( STATUS "UNSIOPATH =" ${UNSIOPATH} )

# Boost (deprecated)
#FIND_PACKAGE(BOOST)

#FIND_PACKAGE(NEMO REQUIRED)
MESSAGE( STATUS "NEMOLIB : " ${NEMOLIB} )

# contains the full path to the top level directory of your build tree 
MESSAGE( STATUS "PROJECT_BINARY_DIR: " ${PROJECT_BINARY_DIR} )
#MESSAGE ( STATUS "Boost_LIBRARY_DIRS : " ${Boost_LIBRARY_DIRS}) 
# contains the full path to the root of your project source directory,
# i.e. to the nearest directory where CMakeLists.txt contains the PROJECT() command 
MESSAGE( STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )


# Find all the sources for the utils LIB
if (CMAKE_Fortran_COMPILER)
  FILE(GLOB LIBUTILS ${PROJECT_SOURCE_DIR}/lib/src/*.cc ${PROJECT_SOURCE_DIR}/lib/src/*.c ${PROJECT_SOURCE_DIR}/lib/src/*.F )
else ()
  FILE(GLOB LIBUTILS ${PROJECT_SOURCE_DIR}/lib/src/*.cc ${PROJECT_SOURCE_DIR}/lib/src/*.c)
endif ()
# create the library "MYutils"
add_library (MYutils SHARED ${LIBUTILS})
if(OSX)
  set_target_properties(JCLutils PROPERTIES LINK_FLAGS "-undefined suppress -flat_namespace")
endif(OSX)

# Destination path for the lib
SET(LIBRARY_OUTPUT_PATH lib)


# compilation options
#add_definitions(${OPT} -DNO_CUDA -DfalcON_NEMO -DfalcON_SINGLE  -funroll-loops -Wall ${DNOBOOST})
add_definitions(${OPT} -DNO_CUDA -DfalcON_NEMO -DfalcON_SINGLE  -Wall ${DNOBOOST})
SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${WARNCPP} )
SET(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${WARNC} )
#SET(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} ${WARNF} )
SET(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} ${EXTRA_Fortran_FLAGS})
message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
message ("Fortran flags   : " ${CMAKE_Fortran_FLAGS})

# remove -DNDEBUG in RElease mode to keep activated assert calls
STRING(REPLACE "-DNDEBUG" ""  CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
#message (STATUS "CMAKE_CXX_FLAGS_RELEASE = " ${CMAKE_CXX_FLAGS_RELEASE}  )

# Executbale output path
set(EXECUTABLE_OUTPUT_PATH bin)
MESSAGE( STATUS "EXECUTABLE_OUTPUT_PATH: " ${EXECUTABLE_OUTPUT_PATH} )

# Make sure the compiler can find include files from our Hello library.
include_directories (${PROJECT_SOURCE_DIR}/src )
include_directories (${PROJECT_SOURCE_DIR}/lib/src)
include_directories (${UNSIOPATH}/lib ${UNSIOPATH}/inc ${UNSIOPATH}/inc/uns)

# Make sure the linker can find libraries once it is built.
link_directories ( ${UNSIOPATH}/lib $ENV{PGPLOT_DIR}  ${G2C_DIR} ${PROJECT_BINARY_DIR}/lib /usr/lib64 /usr/X11/lib ${FC_GFORT_PATH} ${FC_G77_PATH})

# ----------------------------------------------------------
# Install SETUP
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
# Find all c/c++ exes sources files

FILE(GLOB execpp_sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/src/*.c)

# build cpp executables according to the source
FOREACH(exe_cpp ${execpp_sources})

  get_filename_component(exe ${exe_cpp} NAME_WE)
  MESSAGE( STATUS "New executable ==> " ${exe})
  add_executable (${exe} ${exe_cpp})

  #
  target_link_libraries (${exe}  MYutils ${FC_GFORT_LIB} ${FC_G77_LIB} unsio nemo ${SQLITE3_LIB} dl)
  INSTALL(TARGETS ${exe} RUNTIME  DESTINATION bin)
ENDFOREACH(exe_cpp ${execpp_sources})


#
# find all fortran source files
#
if (CMAKE_Fortran_COMPILER)
  # Find all LIB's sources files
  FILE(GLOB exe_sources ${PROJECT_SOURCE_DIR}/src/*.F)
  FOREACH(exe_f ${exe_sources})

    get_filename_component(exe ${exe_f} NAME_WE)
    MESSAGE( STATUS "New executable ==> " ${exe})
    add_executable (${exe} ${exe_f})

    # Link the executable to the Hello library.
    #target_link_libraries (${exe} unsio nemomaing77 nemo g2c sqlite3 stdc++)
    target_link_libraries (${exe} MYutils unsio nemomaing77 nemo ${SQLITE3_LIB}  ${FC_LIB} stdc++)

    SET_TARGET_PROPERTIES(${exe} PROPERTIES LINKER_LANGUAGE Fortran)
    INSTALL(TARGETS ${exe} RUNTIME  DESTINATION bin)
    # add "-lstdc++"
    #SET(CMAKE_Fortran_LINK_FLAGS "-L${NEMOLIB} ") #-lstdc++")
  ENDFOREACH(exe_f ${exe_sources})
endif()

# ----------------------------------------------------------
# install target                                            

IF   (OSX) # Mac OSX
  SET(SOEXT "dylib")
ELSE (OSX) # Linux
  SET(SOEXT "so")
ENDIF(OSX)

INSTALL(FILES  ${PROJECT_BINARY_DIR}/lib/libMYutils.${SOEXT} DESTINATION lib)
#
