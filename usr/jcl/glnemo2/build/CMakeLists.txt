# ============================================================================
# Copyright Jean-Charles LAMBERT - 2007-2011                                  
# e-mail:   Jean-Charles.Lambert@oamp.fr                                      
# address:  Dynamique des galaxies                                            
#           Laboratoire d'Astrophysique de Marseille                          
#           Pôle de l'Etoile, site de Château-Gombert                         
#           38, rue Frédéric Joliot-Curie                                     
#           13388 Marseille cedex 13 France                                   
#           CNRS U.M.R 6110                                                   
# ============================================================================
# See the complete license in LICENSE and/or "http://www.cecill.info".        
# ============================================================================

cmake_minimum_required(VERSION 2.6)
 
# project name
project (glnemo2)

# check for nemo
SET(NEMO_INSTALLED FALSE)
FILE(GLOB GLOB_TEMP_VAR $ENV{NEMO})
IF(GLOB_TEMP_VAR)
  SET(NEMO_INSTALLED TRUE)
ENDIF(GLOB_TEMP_VAR)

# find QT4 files
FIND_PACKAGE(Qt4 4.7.4 REQUIRED)
SET(QT_USE_QTOPENGL 1)
SET(QT_USE_QTUITOOLS 1)
SET(QT_USE_QTNETWORK 1)
INCLUDE(${QT_USE_FILE})

# find OpenGL package
FIND_PACKAGE(OpenGL REQUIRED)

# contains the full path to the top level directory of your build tree 
MESSAGE( STATUS "PROJECT_BINARY_DIR: " ${PROJECT_BINARY_DIR} )

# contains the full path to the root of your project source directory,
# i.e. to the nearest directory where CMakeLists.txt contains the PROJECT() command 
MESSAGE( STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )

# Find nemo files
FILE(GLOB lib_files 
     ../plugins/nemolight/*.c 
     ../plugins/*.cc
     ../plugins/ftm/*.cc 
     ../plugins/gadget/*.cc
     ../plugins/network/*.cc
     ../plugins/network/*.cpp
     ../plugins/ramses/*.cc
     ../plugins/zlib/*.c
     ../utils/*.cc
     ../3rdparty/pfntlib/*.cxx
    )

# Find plugins qt headers
FILE(GLOB hdrs1_files ../plugins/*.h ../plugins/ftm/*.h ..plugins/gadget/*.h 
     ../plugins/network/*.h ../plugins/ramses/*.h
     )

# Destination path for the lib
SET(LIBRARY_OUTPUT_PATH LIB)

# compilation options
add_definitions(-Wall -O2 -ggdb)
ADD_DEFINITIONS(-DQT_PLUGIN)
ADD_DEFINITIONS(-DQT_STATICPLUGIN)
ADD_DEFINITIONS(-DQT_OPENGL_LIB)
ADD_DEFINITIONS(-DQT_GUI_LIB)
ADD_DEFINITIONS(-DQT_CORE_LIB)
ADD_DEFINITIONS(${QT_DEFINITIONS})

# Make sure the compiler can find include files from our Hello library.
include_directories (${PROJECT_SOURCE_DIR}/../plugins ${PROJECT_SOURCE_DIR}/../plugins/nemolight/inc ${PROJECT_SOURCE_DIR}/../src ${PROJECT_SOURCE_DIR}/../plugins/ftm ${PROJECT_SOURCE_DIR}/../plugins/gadget ${PROJECT_SOURCE_DIR}/../src/glew ${PROJECT_SOURCE_DIR}/../plugins/network ${PROJECT_SOURCE_DIR}/../plugins/ramses ${PROJECT_SOURCE_DIR}/../plugins/zlib ${PROJECT_SOURCE_DIR}/../3rdparty/pfntlib ${PROJECT_SOURCE_DIR}/../utils)

# run moc files
QT4_WRAP_CPP(moc1_src ${hdrs1_files})

# create library "uns"
add_library (uns STATIC ${moc1_src} ${lib_files})

#

# -----------------------
# Glnemo2 binary PROGRAM 
# -----------------------

# Make sure the linker can find the Hello library once it is built.
link_directories (${PROJECT_SOURCE_DIR}/../LIB)

# Glew file
set (glew ../src/glew/glew.c) 

# Find all cpp sources files
FILE(GLOB src_files ../src/*.cc)

# find all header's file
FILE(GLOB hdrs_files ../src/*.h)

# Find all UI's files
FILE(GLOB ui_files ../src/*.ui)

# Find Resource file
SET(rcs_file ../src/glnemo.qrc)

# compute ressources
QT4_ADD_RESOURCES( rc_src ${rcs_file} )

# run uic files
QT4_WRAP_UI(ui_headers ${ui_files})

# run moc files
QT4_WRAP_CPP(moc_src ${hdrs_files})

# we need this to be able to include headers produced by uic in our code
# (CMAKE_BINARY_DIR holds a path to the build directory, while INCLUDE_DIRECTORIES() 
# works just like INCLUDEPATH from qmake)
INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR} )

# Executbale output path
set(EXECUTABLE_OUTPUT_PATH ../bin)

add_executable (glnemo2 ${src_files} ${rc_src} ${ui_headers} ${moc_src} ${glew})

# Link the executable to the Hello library.
target_link_libraries (glnemo2 uns ${OPENGL_LIBRARY}  ${QT_LIBRARIES} )

# INSTALLATION

IF (NEMO_INSTALLED) 
  set (CMAKE_INSTALL_PREFIX $ENV{NEMO})

  INSTALL(PROGRAMS  ${EXECUTABLE_OUTPUT_PATH}/glnemo2     DESTINATION bin)
  INSTALL(FILES  ${PROJECT_BINARY_DIR}/../man/man1/glnemo2.1 DESTINATION man/man1)

ENDIF (NEMO_INSTALLED) 
#
