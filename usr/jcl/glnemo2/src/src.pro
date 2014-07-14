# -------------------------------------------
# Subdir relative project main directory: ./src
# Target is an application: ../bin/architecture/release|debug/glnemo2
# include global configuration
include(../config.arch)
DEFINES += GLEW_STATIC
FORMS += formobjectcontrol.ui \
    formabout.ui \
    formscreenshot.ui \
    formselectpart.ui \
    formoptions.ui \
    formconnect.ui \
    formhelp.ui
HEADERS += mainwindow.h \
    glwindow.h \
    globject.h \
    glgridobject.h \
    globaloptions.h \
    snapshotinterface.h \
    pluginsmanage.h \
    particlesobject.h \
    particlesselectrange.h \
    particlesdata.h \
    globjectparticles.h \
    loadingthread.h \
    formobjectcontrol.h \
    formabout.h \
    gltexture.h \
    componentrange.h \
    userselection.h \
    glselection.h \
    frustumculling.h \
    tools3d.h \
    vec3d.h \
    gltextobject.h \
    globjectosd.h \
    formscreenshot.h \
    formselectpart.h \
    formoptions.h \
    orbits.h \
    gloctree.h \
    densityhisto.h \
    colormap.h \
    densitycolorbar.h \
    camera.h \
    catmull_rom_spline.h \
    formconnect.h \
    glcubeobject.h \
    formhelp.h \
    cshader.h \
    glcolorbar.h \
    glaxesobject.h
SOURCES += glnemo.cc \
    mainwindow.cc \
    glwindow.cc \
    globject.cc \
    glgridobject.cc \
    globaloptions.cc \
    pluginsmanage.cc \
    particlesobject.cc \
    particlesselectrange.cc \
    particlesdata.cc \
    globjectparticles.cc \
    loadingthread.cc \
    formobjectcontrol.cc \
    formabout.cc \
    gltexture.cc \
    componentrange.cc \
    userselection.cc \
    glselection.cc \
    frustumculling.cc \
    tools3d.cc \
    vec3d.cc \
    gltextobject.cc \
    globjectosd.cc \
    formscreenshot.cc \
    formselectpart.cc \
    formoptions.cc \
    orbits.cc \
    glew/glew.c \
    gloctree.cc \
    densityhisto.cc \
    colormap.cc \
    densitycolorbar.cc \
    camera.cc \
    catmull_rom_spline.cc \
    snapshotinterface.cc \
    formconnect.cc \
    glcubeobject.cc \
    formhelp.cc \
    cshader.cc \
    glcolorbar.cc \
    glaxesobject.cc
RESOURCES = glnemo.qrc
CONFIG += $$GLOBAL \
    warn_on \
    opengl \
    thread
CONFIG(debug, debug|release) { 
    TARGET = ../bin/$$ARCH/$$COMPILEMODE/glnemo2
    win32 { 
        DESTDIR = ../bin/$$COMPILEMODE/$$ARCH
        TARGET = glnemo2
    }
    unix:
}

else { 
    TARGET = ../bin/$$ARCH/$$COMPILEMODE/glnemo2
    win32 { 
        DESTDIR = ../bin/$$COMPILEMODE/$$ARCH
        TARGET = glnemo2
    }
}
TEMPLATE = app
QT      += opengl
QT      += network
QT      += printsupport
#QT      += widgets
MOC_DIR     = .moc/$$ARCH/$$COMPILEMODE
UI_DIR      = ._ui/$$ARCH/$$COMPILEMODE
OBJECTS_DIR = .obj/$$ARCH/$$COMPILEMODE
RCC_DIR     = .res/$$ARCH/$$COMPILEMODE
QMAKE_LIBDIR = \
    ../utils/lib/$$ARCH/$$COMPILEMODE \
    ../3rdparty/pfntlib/lib/$$ARCH/$$COMPILEMODE \
    ../plugins/lib/$$ARCH/$$COMPILEMODE \
    $$NEMOLIB \
    ../plugins/ftm/lib/$$ARCH/$$COMPILEMODE \
    ../plugins/ramses/lib/$$ARCH/$$COMPILEMODE \
    ../plugins/gadget/lib/$$ARCH/$$COMPILEMODE \
    ../plugins/zlib/lib/$$ARCH/$$COMPILEMODE \
    ../plugins/tipsy/lib/$$ARCH/$$COMPILEMODE \
    ../plugins/network/lib/$$ARCH/$$COMPILEMODE

# Icons
macx {
  ICON = ../res/images/glnemo2.icns
}
# INSTALLS for Linux and Mac OS X
MYNEMO = $$(NEMO)
!win32 {
   !isEmpty( MYNEMO ) {  # NEMO repository exist and loaded
      target.path = $$(NEMO)/bin
      man.path    = $$(NEMO)/man/man1
      man.files   = ../man/man1/glnemo2.1
      INSTALLS   += man
   } else            {  # NEMO repository does not exist           
      target.path = $$(HOME)/bin
   }
   INSTALLS      += target
}
INCLUDEPATH += \
    ../3rdparty/pfntlib/ \
    ../plugins \
    ../plugins/ftm \
    ../plugins/nemolight \
    ../plugins/gadget \
    ../plugins/zlib \
    ../plugins/network \
    ../plugins/tipsy \
    ../src \
    $$NEMOLIB \
    $$NEMOINC \
    glew
LIBS += \
    -lpfntlib \
    -lsnapshot \
    -lftm \
    -lnemo \
    -lgadget \
    -lramses \
    -lnetwork \
    -ltipsy \
    -lutils
win32 {
    LIBS += -lzlib
    LIBS += -lopengl32
}
unix {
    LIBS += -lz
    LIBS += -lGLU
}
macx {
    LIBS += -lz
    LIBS -= -lGLU
}
POST_TARGETDEPS += \
    ../3rdparty/pfntlib/lib/$$ARCH/$$COMPILEMODE/libpfntlib.a \
    ../plugins/nemolight/lib/$$ARCH/$$COMPILEMODE/libnemo.a \
    ../plugins/lib/$$ARCH/$$COMPILEMODE/libsnapshot.a \
    ../plugins/ftm/lib/$$ARCH/$$COMPILEMODE/libftm.a \
    ../plugins/gadget/lib/$$ARCH/$$COMPILEMODE/libgadget.a \
    ../plugins/ramses/lib/$$ARCH/$$COMPILEMODE/libramses.a \
    ../plugins/network/lib/$$ARCH/$$COMPILEMODE/libnetwork.a \
    ../plugins/zlib/lib/$$ARCH/$$COMPILEMODE/libzlib.a \
    ../plugins/tipsy/lib/$$ARCH/$$COMPILEMODE/libtipsy.a \
    ../utils/lib/$$ARCH/$$COMPILEMODE/libutils.a
DISTFILES += ../ChangeLog

# No TPSY support for windows
win32 {
QMAKE_LIBDIR -= ../plugins/tipsy/lib/$$ARCH/$$COMPILEMODE
INCLUDEPATH  -= ../plugins/tipsy
LIBS         -= -ltipsy
POST_TARGETDEPS -= ../plugins/tipsy/lib/$$ARCH/$$COMPILEMODE/libtipsy.a

}


