# ============================================================================
# Copyright Jean-Charles LAMBERT - 2010                                       
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
project (uns_projects)

# set CMAKE Module PATH (NEMO)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_BINARY_DIR})
MESSAGE( STATUS "CMAKE_MODULE_PATH=" ${CMAKE_MODULE_PATH} )

# require NEMO
FIND_PACKAGE(NEMO REQUIRED)
# require boost
find_package( Boost 1.36.0 REQUIRED)

# Some necessary variable
set(DEHNEN  $ENV{NEMO}/usr/dehnen)
set(NEMOLIB $ENV{NEMOLIB})

#FIND_PACKAGE(NEMO REQUIRED)
MESSAGE( STATUS "NEMOLIB : " ${NEMOLIB} )

# contains the full path to the top level directory of your build tree 
MESSAGE( STATUS "PROJECT_BINARY_DIR: " ${PROJECT_BINARY_DIR} )
#MESSAGE ( STATUS "Boost_LIBRARY_DIRS : " ${Boost_LIBRARY_DIRS}) 
# contains the full path to the root of your project source directory,
# i.e. to the nearest directory where CMakeLists.txt contains the PROJECT() command 
MESSAGE( STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )

# find lib g2c
execute_process(COMMAND g77 -print-file-name=libg2c.a 
	OUTPUT_VARIABLE G2CFPATH OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND dirname ${G2CFPATH} 
	OUTPUT_VARIABLE G2C_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)

MESSAGE( STATUS "G2C dir : " ${G2C_DIR})

# Find all the sources for the utils LIB
FILE(GLOB LIBUTILS ../lib/utils/*.cc)
# create the library "JCLutils"
add_library (JCLutils SHARED ${LIBUTILS})

# Find all the sources for the projects LIB
FILE(GLOB LIBPROJECTS ../lib/projects/*.cc)
# create the library "JCLutils"
add_library (JCLprojects SHARED ${LIBPROJECTS})

# Destination path for the lib
SET(LIBRARY_OUTPUT_PATH ../lib)


# compilation options
add_definitions( -DfalcON_NEMO -DfalcON_SINGLE -g -ggdb -funroll-loops -Wall)
#add_definitions(-O2 -ggdb -DfalcON_NEMO -D_FILE_OFFSET_BITS=64 -mpreferred-stack-boundary=4 -funroll-loops -fforce-addr -rdynamic  -DfalcON_SINGLE)


# Executbale output path
set(EXECUTABLE_OUTPUT_PATH ../bin)
MESSAGE( STATUS "EXECUTABLE_OUTPUT_PATH: " ${EXECUTABLE_OUTPUT_PATH} )

# Make sure the compiler can find include files from our Hello library.
include_directories (../lib/utils ../lib/projects ${DEHNEN}/falcON/inc ${DEHNEN}/falcON/utils/inc $ENV{NEMOINC} ${NEMOLIB} $ENV{NEMOINC}/uns)


# Make sure the linker can find the Hello library once it is built.
link_directories (${NEMOLIB} ${DEHNEN}/falcON/lib ${DEHNEN}/utils/lib $ENV{PGPLOT_DIR}  ${G2C_DIR} ${PROJECT_BINARY_DIR}/../lib /usr/lib64)


# ----------------------------------------------------------
# Make sure the compiler can find include files from our Hello library.
#include_directories (${UNS_SOURCE_DIR}/src $ENV{NEMOINC} ${NEMOLIB})

# Find all LIB's sources files
FILE(GLOB execpp_sources ../src/*.cc ../src/*.c)

# build cpp executables according to the source
FOREACH(exe_cpp ${execpp_sources})

  get_filename_component(exe ${exe_cpp} NAME_WE)
  MESSAGE( STATUS "New executable ==> " ${exe})
  add_executable (${exe} ${exe_cpp})

  # Link the executable to the Hello library.
#  target_link_libraries (${exe}  JCLutils JCLprojects cpgplot pgplot X11 g2c unsio nemo sqlite3 ${Boost_LIBRARIES} boost_thread-mt.a boost_system-mt.a)
 target_link_libraries (${exe}  JCLutils JCLprojects cpgplot pgplot X11 g2c unsio WDutils falcON nemo sqlite3  boost_thread-mt boost_system-mt dl)


ENDFOREACH(exe_cpp ${execpp_sources})

# ----------------------------------------------------------
# install target                                            

set (CMAKE_INSTALL_PREFIX $ENV{NEMO})

# 
# install bin targets
FOREACH(exe_cpp ${execpp_sources})

  get_filename_component(exe ${exe_cpp} NAME_WE)
  INSTALL(PROGRAMS  ${PROJECT_BINARY_DIR}/../bin/${exe}     DESTINATION bin)
ENDFOREACH(exe_cpp ${execpp_sources})
INSTALL(FILES  ${PROJECT_BINARY_DIR}/../lib/libJCLutils.so    DESTINATION lib)
INSTALL(FILES  ${PROJECT_BINARY_DIR}/../lib/libJCLprojects.so DESTINATION lib)
INSTALL(FILES  ${PROJECT_BINARY_DIR}/../man/uns_2dplot.1 DESTINATION man/man1)
INSTALL(FILES  ${PROJECT_BINARY_DIR}/../lib/utils/cfalcon.h DESTINATION inc/uns)
INSTALL(FILES  ${PROJECT_BINARY_DIR}/../lib/utils/csnaptools.h DESTINATION inc/uns)
INSTALL(FILES  ${PROJECT_BINARY_DIR}/../lib/utils/ctimer.h DESTINATION inc/uns)
#
