# GyrfalcON
GYRFALCON = $(NEMO)/usr/dehnen/falcON
include(../../config.arch)
MOC_DIR = .moc/$${ARCH}/$$COMPILEMODE
UI_DIR = .ui/$${ARCH}/$$COMPILEMODE
OBJECTS_DIR = .obj/$${ARCH}/$$COMPILEMODE
DESTDIR = ../lib/$${ARCH}/$$COMPILEMODE
TEMPLATE = lib
TARGET = glnemo2
DEPENDPATH += .
INCLUDEPATH += ./ \
    ./inc \
    $$GYRFALCON/inc \
    $$GYRFALCON/utils/inc \
    $$GYRFALCON/../utils/inc \
    ../../plugins/network \
    $(NEMOINC) \
    $(NEMOLIB)
QMAKE_CFLAGS += -DfalcON_INDI \
    -DfalcON_NEMO -DfalcON_SINGLE

QMAKE_CXXFLAGS += -DfalcON_INDI \
    -DfalcON_NEMO -DfalcON_SINGLE

unix {
   QMAKE_CFLAGS   += -rdynamic
   QMAKE_CXXFLAGS += -rdynamic
}
QMAKE_LFLAGS = -shared
macx {
  QMAKE_LIBDIR =  $(NEMO)/lib $(NEMO)/usr/dehnen/utils
  #QMAKE_LFLAGS += -Wl,-rpath,$(NEMO)/usr/dehnen/utils
  #QMAKE_LFLAGS += -Wl,-rpath,$(NEMO)/lib
  LIBS +=  -lWDutils -lfalcon
}

CONFIG += console
QT += network
# INSTALLS
target.path  = $$(NEMO)/obj/manip
unix:target.extra = cp ../lib/$${ARCH}/$$COMPILEMODE/libglnemo2.so $$(NEMO)/obj/manip/glnemo2.so
macx:target.extra = cp ../lib/$${ARCH}/$$COMPILEMODE/libglnemo2.1.0.0.dylib $$(NEMO)/obj/manip/glnemo2.so
INSTALLS += target

SOURCES += glnemo2.cc \
    master_server_thread.cc \
    client.cc \
    ../../plugins/network/sendreseau.cpp \
    ../../plugins/network/readreseau.cpp \
    ../../plugins/network/global.cpp \
    masterserver.cc
HEADERS += master_server_thread.h \
    client.h \
    ../../plugins/network/sendreseau.h \
    ../../plugins/network/readreseau.h \
    ../../plugins/network/global.h \
    masterserver.h
