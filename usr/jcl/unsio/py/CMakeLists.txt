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
# CMakeListst.txt file for UNS python library                                     
# ============================================================================

# This is a CMake example for Python
cmake_minimum_required(VERSION 2.6) 

# project name
project (py_unsio)

# extra path for NEMO
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake )
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/../cmake )

# load setup flags
include(SetupFlags)

# required package
FIND_PACKAGE(SWIG REQUIRED)
FIND_PACKAGE(NumPy REQUIRED)
FIND_PACKAGE(UNSIO REQUIRED)
FIND_PACKAGE(SQLITE3 REQUIRED)

INCLUDE(${SWIG_USE_FILE})

FIND_PACKAGE(PythonLibs)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

# some status
MESSAGE( STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )
MESSAGE( STATUS "PROJECT_BINARY_DIR: " ${PROJECT_BINARY_DIR} )
MESSAGE( STATUS "CMAKE_BUILD_TYPE =" ${CMAKE_BUILD_TYPE} )
MESSAGE( STATUS "CMAKE_CXX_FLAGS  =" ${CMAKE_CXX_FLAGS} )
MESSAGE( STATUS "LIBRARY TYPE     =" ${LIBTYPE} )
MESSAGE( STATUS "core OPT  =" ${OPT} )
MESSAGE( STATUS "OSX detected =" ${OSX} )
MESSAGE( STATUS "BUILD_TYPE =" ${RELEASE} )

MESSAGE(STATUS "PYTHON_INCLUDE_PATH:" ${PYTHON_INCLUDE_PATH})
MESSAGE(STATUS "PYTHON_INCLUDE_DIRS:" ${PYTHON_INCLUDE_DIRS})
MESSAGE(STATUS "PYTHON_LIBRARY :" ${PYTHON_LIBRARY})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/swig ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/../src ${NUMPY_INCLUDE_DIRS})


# Set definitions
set (INTERFACE swig/py_unsio.i)
set (MODULENAME py_unsio)

# Find all LIB's sources files
FILE(GLOB execpp_sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h )

SET(CMAKE_SWIG_FLAGS "-c++")

SET_SOURCE_FILES_PROPERTIES(${INTERFACE} PROPERTIES CPLUSPLUS ON)

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/swig)

set(CMAKE_SWIG_OUTDIR ${PROJECT_BINARY_DIR}/swig) 

IF   (OSX) # Mac OSX
  SET(SOEXT "dylib")
ELSE (OSX) # Linux
  SET(SOEXT "so")
ENDIF(OSX)

SWIG_ADD_MODULE(${MODULENAME} python  ${INTERFACE} ${execpp_sources} ${UNSIOPATH}/lib/libunsio.${SOEXT} ${UNSIOPATH}/lib/libnemo.a )

SWIG_LINK_LIBRARIES(${MODULENAME} ${PYTHON_LIBRARIES} ${UNSIOPATH}/lib/libunsio.${SOEXT} ${UNSIOPATH}/lib/libnemo.a ${SQLITE3_LIB_PATH}/libsqlite3.${SOEXT} )
# ----------------------------------------------------------
# Install SETUP
# ----------------------------------------------------------

mark_as_advanced(CMAKE_INSTALL_PREFIX)



set (CMAKE_INSTALL_PREFIX ${UNSIOPATH})
MESSAGE( STATUS "--------> CMAKE_INSTALL_PREFIX =" ${CMAKE_INSTALL_PREFIX})

INSTALL(FILES  ${PROJECT_BINARY_DIR}/_${MODULENAME}.so DESTINATION py/modules)
INSTALL(FILES  ${PROJECT_BINARY_DIR}/swig/${MODULENAME}.py DESTINATION py/modules)

MESSAGE( STATUS "\nAfter running \"make install\", do no forget to set PYTHONPATH variable with command:\n  export PYTHONPATH=${CMAKE_INSTALL_PREFIX}/py/modules\n\n")
#
