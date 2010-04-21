######################################################################

######################################################################

# include global configuration
include(../config.arch)
MOC_DIR = .moc/$${ARCH}/$$COMPILEMODE
UI_DIR = .ui/$${ARCH}/$$COMPILEMODE
OBJECTS_DIR = .obj/$${ARCH}/$$COMPILEMODE
TEMPLATE = lib
TARGET = utils
DESTDIR = lib/$${ARCH}/$$COMPILEMODE

CONFIG += $$GLOBAL \
    warn_on \
    static \
    staticlib

DEPENDPATH += .
INCLUDEPATH += . \
    $$NEMOINC \
    $$NEMOLIB

# Input
HEADERS += cfortio.h
SOURCES += cfortio.cc
