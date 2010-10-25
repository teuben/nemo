######################################################################

INCLUDEPATH += .
# include global configuration
include(../../config.arch)

OBJECTS_DIR = .obj/$${ARCH}/$$COMPILEMODE
TARGET  = pfntlib
DESTDIR = lib/$${ARCH}/$$COMPILEMODE
CONFIG += $$GLOBAL warn_on \
          lib  \
	  static \
          staticlib
          
TEMPLATE = lib 


# Input
HEADERS += fnt.h fntLocal.h sg.h ul.h ulLocal.h ulRTTI.h
SOURCES += fnt.cxx \
           fntBitmap.cxx \
           fntTXF.cxx \
           sg.cxx \
           sgd.cxx \
           sgdIsect.cxx \
           sgIsect.cxx \
           sgPerlinNoise.cxx \
           ul.cxx \
           ulClock.cxx \
           ulError.cxx \
           ulLinkedList.cxx \
           ulList.cxx \
           ulRTTI.cxx
