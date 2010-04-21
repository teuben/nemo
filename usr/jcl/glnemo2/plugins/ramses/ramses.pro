# include global configuration
include(../../config.arch)

OBJECTS_DIR = .obj/$${ARCH}/$$COMPILEMODE
DESTDIR = lib/$${ARCH}/$$COMPILEMODE
TARGET = ramses

CONFIG += $$GLOBAL warn_on \
          lib  \
	  static \
          staticlib
          
TEMPLATE = lib

INCLUDEPATH += \
../../src \
../../utils

# Input
HEADERS += camr.h cpart.h
SOURCES += camr.cc cpart.cc
