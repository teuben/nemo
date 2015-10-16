#-*-cmake-*- 
# ============================================================================
# Copyright Jean-Charles LAMBERT - 2007-2015                                  
# e-mail:   Jean-Charles.Lambert@lam.fr                                      
# address:  Centre de donneeS Astrophysique de Marseille (CeSAM)                                          
#           Laboratoire d'Astrophysique de Marseille                          
#           Pôle de l'Etoile, site de Château-Gombert                         
#           38, rue Frédéric Joliot-Curie                                     
#           13388 Marseille cedex 13 France                                   
#           CNRS U.M.R 7326                                                   
# ============================================================================
# See the complete license in LICENSE and/or "http://www.cecill.info".        
# ============================================================================

cmake_minimum_required(VERSION 2.8)

# project name
project (glnemo2)

# extra path for .cmake .files
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# Cpack Generator detection
include(DetectCpackGenerator)

if    (WIN32)
  MESSAGE(STATUS "\n\n\n- - - WINDOWS  - - -\n\n\nWARNING : DO NOT FORGET TO RUN CMAKE LIKE THIS\n\n ")
  MESSAGE(STATUS "cmake -G \"MinGW Makefiles\" .. -DCMAKE_PREFIX_PATH=c:\\QtSDK-x86_64\n\n\n")
endif ()
if    (APPLE)
  MESSAGE(STATUS "\n\n\n- - - Mac OSX  - - -\n\n\nWARNING : DO NOT FORGET TO RUN CMAKE LIKE THIS\n\n ")
  MESSAGE(STATUS "cmake .. -DCMAKE_PREFIX_PATH=/Users/jcl/Qt/5.4/clang_64/lib/cmake\n\n\n")
endif ()
# check QTDIR5
FILE(GLOB GLOB_TEMP_VAR1 $ENV{QTDIR5})
IF(NOT GLOB_TEMP_VAR1)
  MESSAGE(STATUS "\n\nWARNING,  environment variable >>> QTDIR5 <<< is not defined")
  MESSAGE(STATUS "If it does not compile, try to set it to your QT5 location directory\n\n")
  SET(QTDIR5 "/usr")
ELSE (NOT GLOB_TEMP_VAR1)
  SET(QTDIR5 $ENV{QTDIR5})
  MESSAGE(STATUS "QTDIR5 location = $ENV{QTDIR5}")
ENDIF()

# check for nemo
SET(NEMO_INSTALLED FALSE)
FILE(GLOB GLOB_TEMP_VAR $ENV{NEMO})
IF(GLOB_TEMP_VAR)
  SET(NEMO_INSTALLED TRUE)
ENDIF()

# MODULEs PATH
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake ${QTDIR5}) 

# Instruct CMake to run moc automatically when needed.
#set(CMAKE_AUTOMOC ON)

# find QT5 files
find_package(Qt5Widgets)
find_package(Qt5Core)
find_package(Qt5OpenGL)
find_package(Qt5Network)
find_package(Qt5PrintSupport)
find_package(Qt5Gui)

# find OpenGL package
find_package(OpenGL REQUIRED)

# contains the full path to the top level directory of your build tree 
MESSAGE( STATUS "PROJECT_BINARY_DIR: " ${PROJECT_BINARY_DIR} )

# contains the full path to the root of your project source directory,
# i.e. to the nearest directory where CMakeLists.txt contains the PROJECT() command 
MESSAGE( STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )

# Library's files
FILE(GLOB lib_files 
  "plugins/nemolight/*.c"
  "plugins/*.cc"
  "plugins/ftm/*.cc"
  "plugins/gadget/*.cc"
  "plugins/tipsy/*.cc"
  "plugins/network/*.cc"
  "plugins/network/*.cpp"
  "plugins/ramses/*.cc"
  "utils/*.cc"
  "3rdparty/pfntlib/*.cxx"
  )

# Find plugins qt headers
FILE(GLOB hdrs1_files plugins/*.h plugins/ftm/*.h ..plugins/gadget/*.h 
  plugins/network/*.h plugins/ramses/*.h  plugins/tipsy/*.h 
  )

# Tipsy is not suported by  WIN32
if   (WIN32) 
  list (REMOVE_ITEM lib_files   ${CMAKE_CURRENT_SOURCE_DIR}/plugins/snapshottipsy.cc ${CMAKE_CURRENT_SOURCE_DIR}/plugins/tipsy/tipsyio.cc)   
  list (REMOVE_ITEM hdrs1_files ${CMAKE_CURRENT_SOURCE_DIR}/plugins/snapshottipsy.h  ${CMAKE_CURRENT_SOURCE_DIR}/plugins/tipsy/tipsyio.h) 
endif()

# Destination path for the lib
SET(LIBRARY_OUTPUT_PATH LIB)

# compilation options
add_definitions(-Wall -O2 -D_GLIBCXX_PARALLEL -ggdb)
ADD_DEFINITIONS(-DQT_PLUGIN)
ADD_DEFINITIONS(-DQT_STATICPLUGIN)
ADD_DEFINITIONS(-DQT_OPENGL_LIB)
ADD_DEFINITIONS(-DQT_GUI_LIB)
ADD_DEFINITIONS(-DQT_CORE_LIB)
ADD_DEFINITIONS(-DQT_PRINTSUPPORT_LIB)
ADD_DEFINITIONS(${QT_DEFINITIONS})
if   (WIN32)
  ADD_DEFINITIONS(-DUNICODE -DGLEW_STATIC)
endif()

if (NOT APPLE)
  ADD_DEFINITIONS(-fopenmp)
endif()

# Add compiler flags for building executables (-fPIE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")

# Make sure the compiler can find include files from our Hello library.
include_directories (${PROJECT_SOURCE_DIR}/plugins ${PROJECT_SOURCE_DIR}/plugins/nemolight/inc ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/plugins/ftm ${PROJECT_SOURCE_DIR}/plugins/gadget ${PROJECT_SOURCE_DIR}/plugins/tipsy ${PROJECT_SOURCE_DIR}/src/glew ${PROJECT_SOURCE_DIR}/plugins/network ${PROJECT_SOURCE_DIR}/plugins/ramses ${PROJECT_SOURCE_DIR}/plugins/zlib ${PROJECT_SOURCE_DIR}/3rdparty/pfntlib ${PROJECT_SOURCE_DIR}/utils)
include_directories (${Qt5Widgets_INCLUDE_DIRS} ${Qt5Core_INCLUDE_DIRS} ${Qt5OpenGL_INCLUDE_DIRS} ${Qt5Network_INCLUDE_DIRS} {Qt5PrintSupport_INCLUDE_DIRS})
include_directories(/usr/include/CCfits /usr/include/cfitsio /usr/local/include/CCfits)

if (APPLE)
  include_directories(/usr/include/malloc /usr/X11/include)
endif()
# run moc files
QT5_WRAP_CPP(moc1_src ${hdrs1_files})

# create library "uns"
add_library (uns STATIC ${moc1_src} ${lib_files})

#

# -----------------------
# Glnemo2 binary PROGRAM 
# -----------------------

# Make sure the linker can find the Hello library once it is built.
link_directories (${CMAKE_BINARY_DIR}/LIB /usr/local/lib )

# Glew file
set (glew src/glew/glew.c) 

# Find all cpp sources files
FILE(GLOB src_files src/*.cc)

# find all header's file
FILE(GLOB hdrs_files src/*.h)

# Find all UI's files
FILE(GLOB ui_files src/*.ui)

# Find Resource file
SET(rcs_file src/glnemo.qrc)

# compute ressources
QT5_ADD_RESOURCES( rc_src ${rcs_file} )

# run uic files
QT5_WRAP_UI(ui_headers ${ui_files})

# run moc files
QT5_WRAP_CPP(moc_src ${hdrs_files})

# we need the following to be able to include headers produced by uic in our code
# (CMAKE_BINARY_DIR holds a path to the build directory, while INCLUDE_DIRECTORIES() 
# works just like INCLUDEPATH from qmake)
INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR} )

# Executbale output path
if (UNIX AND NOT APPLE ) 
  set(EXECUTABLE_OUTPUT_PATH bin)
endif ()

if (APPLE)
  set(EXECUTABLE_OUTPUT_PATH "")
endif ()

add_executable (glnemo2 WIN32 MACOSX_BUNDLE ${src_files} ${rc_src} ${ui_headers} ${moc_src} ${glew})

# Link the executable to the Hello library.
set (GOMP "gomp")
if (APPLE)
  set (GOMP "")
endif()
target_link_libraries (glnemo2 uns z ${OPENGL_LIBRARY}  ${QT_LIBRARIES} ${GOMP} CCfits cfitsio)

qt5_use_modules(glnemo2 Gui Widgets Core OpenGL Network PrintSupport)

# INSTALLATION
IF (NEMO_INSTALLED) 
  set (CMAKE_INSTALL_PREFIX $ENV{NEMO})
  INSTALL(FILES  ${PROJECT_SOURCE_DIR}/man/man1/glnemo2.1 DESTINATION man/man1)
ELSE()
  set(CMAKE_INSTALL_PREFIX  "/usr")
ENDIF () 

MESSAGE(STATUS "CMAKE_INSTALL_PREFIX :" ${CMAKE_INSTALL_PREFIX} )

# VERSION
if (EXISTS "${PROJECT_SOURCE_DIR}/src/version.h") # detect version from "version.h" file
  set(CURRENT_VERSION "")
  file(READ "${PROJECT_SOURCE_DIR}/src/version.h" CURRENT_VERSION_H)

  # GLNEMO2_MAJOR
  string(REGEX MATCH "GLNEMO2_MAJOR \"([a-zA-Z0-9.-]+)\"" XXXX "${CURRENT_VERSION_H}")
  set(PROJECT_VERSION_MAJOR ${CMAKE_MATCH_1})

  # GLNEMO2_MINOR
  string(REGEX MATCH "GLNEMO2_MINOR \"([a-zA-Z0-9.-]+)\"" XXXX "${CURRENT_VERSION_H}")
  set(PROJECT_VERSION_MINOR ${CMAKE_MATCH_1})

  #  GLNEMO2_PATCH
  string(REGEX MATCH "GLNEMO2_PATCH \"([a-zA-Z0-9.-]+)\"" XXXX "${CURRENT_VERSION_H}")
  set(PROJECT_VERSION_PATCH ${CMAKE_MATCH_1})

  # GLNEMO2_EXTRA
  string(REGEX MATCH "GLNEMO2_EXTRA \"([a-zA-Z0-9.-]+)\"" XXXX "${CURRENT_VERSION_H}")
  set(PROJECT_VERSION_EXTRA ${CMAKE_MATCH_1})
  
  # PROJECT_VERSION
  set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}${PROJECT_VERSION_EXTRA}")

  message(STATUS "Project Version = " ${PROJECT_VERSION})
endif()

#set(PROJECT_VERSION_MAJOR "1")
#set(PROJECT_VERSION_MINOR "9")
#set(PROJECT_VERSION_PATCH "0")
#set(PROJECT_VERSION_EXTRA "-dev")
#set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}${PROJECT_VERSION_EXTRA}")


if (APPLE)
  #set(CMAKE_INSTALL_PREFIX "/Applications")
  set(CMAKE_INSTALL_PREFIX "")
  set(BIN_INSTALL_DIR ".")
  set(DOC_INSTALL_DIR ".")
endif()

message("${PROJECT_NAME} will be installed to ${CMAKE_INSTALL_PREFIX}")
#install(TARGETS ${PROJECT_NAME} DESTINATION ${BIN_INSTALL_DIR})

if (UNIX AND NOT APPLE)
  INSTALL(PROGRAMS   ${PROJECT_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/glnemo2   DESTINATION bin)
  #INSTALL(CODE "execute_process(COMMAND ${CMAKE_STRIP} bin/glnemo2)") 
  #INSTALL(CODE "execute_process(COMMAND ${CMAKE_STRIP} _CPack_Packages/${CPACK_TOPLEVEL_TAG}/${CPACK_GENERATOR}/${PACK_PACKAGE_FILE_NAME}/${CMAKE_INSTALL_PREFIX}/bin/glnemo2)")
  #include(PostInstall)
  if (NOT NEMO_INSTALLED)
  	INSTALL(FILES  ${PROJECT_SOURCE_DIR}/man/man1/glnemo2.1 DESTINATION share/man/man1)
  endif()
endif ()

if (APPLE)
  INSTALL(PROGRAMS  ${PROJECT_BINARY_DIR}/glnemo2.app   DESTINATION   . )
endif()

if(APPLE)
  set(MACOSX_BUNDLE_INFO_STRING "${PROJECT_NAME} ${PROJECT_VERSION}")
  set(MACOSX_BUNDLE_BUNDLE_VERSION "${PROJECT_NAME} ${PROJECT_VERSION}")
  set(MACOSX_BUNDLE_LONG_VERSION_STRING "${PROJECT_NAME} ${PROJECT_VERSION}")
  set(MACOSX_BUNDLE_SHORT_VERSION_STRING "${PROJECT_VERSION}")
  #set(MACOSX_BUNDLE_COPYRIGHT "${PROJECT_COPYRIGHT_YEAR} ${PROJECT_VENDOR}")
  set(MACOSX_BUNDLE_ICON_FILE "glnemo2.icns")
  #set(MACOSX_BUNDLE_GUI_IDENTIFIER "${PROJECT_DOMAIN_SECOND}.${PROJECT_DOMAIN_FIRST}")
  set(MACOSX_BUNDLE_BUNDLE_NAME "${PROJECT_NAME}")
  set(MACOSX_BUNDLE_RESOURCES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.app/Contents/Resources")
  set(MACOSX_BUNDLE_ICON "${CMAKE_CURRENT_SOURCE_DIR}/res/images/glnemo2.icns")
  execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${MACOSX_BUNDLE_RESOURCES})
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${MACOSX_BUNDLE_ICON} ${MACOSX_BUNDLE_RESOURCES})
endif()

# CPACK infos
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_NAME "glnemo2")
set(CPACK_PACKAGE_RELEASE 1)
set(CPACK_PACKAGE_CONTACT "Jean-Charles LAMBERT")
set(CPACK_PACKAGE_VENDOR "LAM - CeSAM")
set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}.${LINUX_NAME}.${CMAKE_SYSTEM_PROCESSOR}")
#set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_PACKAGE_RELEASE}.${LINUX_NAME}.${CMAKE_SYSTEM_PROCESSOR}")

# CPACK RPM specific
set(CPACK_RPM_CHANGELOG_FILE "${PROJECT_SOURCE_DIR}/ChangeLog")
set(CPACK_RPM_PACKAGE_LICENSE "CeCILL2")
set(CPACK_RPM_PACKAGE_DESCRIPTION "An interactive 3D visualization program for N-body snaphots")

# CPACK debian specific
set(CPACK_DEBIAN_PACKAGE_NAME ${CPACK_PACKAGE_NAME})
set(CPACK_DEBIAN_PACKAGE_VERSION ${CPACK_PACKAGE_VERSION}-${CPACK_PACKAGE_RELEASE}.${LINUX_NAME})
IF(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
	SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE amd64)
ELSEIF(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "i686")
	SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE i386)
ENDIF(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
set(CPACK_DEBIAN_PACKAGE_DESCRIPTION ${CPACK_RPM_PACKAGE_DESCRIPTION})


if (APPLE)
  set(CMAKE_INSTALL_PREFIX "/Applications")
  set(CPACK_GENERATOR "DragNDrop")
  set(CPACK_DMG_FORMAT "UDBZ")
  set(CPACK_DMG_VOLUME_NAME "${PROJECT_NAME}")
  set(CPACK_SYSTEM_NAME "OSX")
  set(CPACK_PACKAGE_FILE_NAME "${PROJECT_NAME}-${PROJECT_VERSION}")
  #set(CPACK_PACKAGE_ICON "${ICONS_DIR}/DMG.icns")
  #set(CPACK_DMG_DS_STORE "${ICONS_DIR}/DMGDSStore")
  #set(CPACK_DMG_BACKGROUND_IMAGE "${ICONS_DIR}/DMGBackground.png")
endif()

# stripping
if (NOT WIN32)
	#SET(CPACK_STRIP_FILES {CMAKE_INSTALL_PREFIX}/bin/glnemo2)
	#SET(CPACK_STRIP_FILES "True")
	#SET(CPACK_SET_DESTDIR On)
        set(CMAKE_INSTALL_DO_STRIP TRUE)  
	#SET(CPACK_INSTALL_PREFIX "/usr")	
	#SET(CMAKE_INSTALL_PREFIX ${CPACK_INSTALL_PREFIX})
endif()

# This must always be last!
include(CPack)

set(CMAKE_INSTALL_SYSTEM_RUNTIME_DESTINATION "${BIN_INSTALL_DIR}")
include(InstallRequiredSystemLibraries)


if(APPLE)
  set(EXECUTABLE "${PROJECT_NAME}.app")
endif()

if(APPLE) # OR WIN32)
  include(DeployQt5)

  #set(LIST_PLUGIN "Qt5::AccessibleFactory;Qt5::QJpegPlugin;Qt5::QICNSPlugin;Qt5::QCocoaPrinterSupportPlugin")
  # The following part give me a big headeach, seriously....
  # it's intended to add plugin in DMG package in order to have all dependencies resolved.
  # So far, I am not abble to add PrintingSupport though... :( ...which should come from
  # Qt5::QCocoaPrinterSupportPlugin.... but fails when I run cmake command...
  #
 
  #set(LIST_PLUGIN "Qt5::AccessibleFactory;Qt5::QJpegPlugin;Qt5::QICNSPlugin")
  set(LIST_PLUGIN "Qt5::QJpegPlugin;Qt5::QICNSPlugin")
  set (QT_ADD_PLUGIN "")
  foreach (plugin ${LIST_PLUGIN})
    get_property(QT_GET_PLUGIN TARGET ${plugin} PROPERTY LOCATION_RELEASE) 
    message( ${plugin} " : " ${QT_GET_PLUGIN})
    get_filename_component(xx ${QT_GET_PLUGIN} DIRECTORY)
    list (APPEND QT_ADD_PLUGIN ${QT_GET_PLUGIN})
  endforeach()
  
  # See macdeployqt source code to find a hack :
  # https://qt.gitorious.org/qt/qttools/source/4b711fec48c395ab4e486dd1a789309d4d488e9d:src/macdeployqt/shared/shared.cpp
  # Ugly hack to add printerplugin support
  
  get_filename_component(xx ${QT_GET_PLUGIN} DIRECTORY)
  get_filename_component(xx ${xx} DIRECTORY)
  set (PrinterPlugin "${xx}/printsupport/libcocoaprintersupport.dylib")
  message ("PrinterPlugin = " ${PrinterPlugin})
  list (APPEND QT_ADD_PLUGIN ${PrinterPlugin})
  
  install_qt5_executable("${EXECUTABLE}" "${QT_ADD_PLUGIN}")
endif()
