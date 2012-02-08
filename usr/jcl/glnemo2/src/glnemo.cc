// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2012                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <QtGui>
#include <QApplication>
#include <GL/glew.h>
#include <QtOpenGL>
#include <QGLFormat>
#include <QDesktopWidget>
#include <iostream>
#include <QSplashScreen>
// Nemo stuffs
#define _vectmath_h // put this statement to avoid conflict with C++ vector class
#include <nemo.h>

#include "mainwindow.h"
using namespace std;

#define RELEASE_VERSION "1.40"

// ============================================================================
// NEMO parameters                                                             
  const char * defv[] = {  
    "in=\n             Input snapshot (Nemo,Gadget 2 & 1, Ramses, phiGrape, ftm, list of files)",
    "select=\n         Select particles using:                                         \n"
    "                   1) component name ex: gas,halo,disk,stars,bulge,bndry           \n"
    "                   2) range operator ex: 0:999,1000:1999 would select two sets     \n"
    "                   of 1000 particles                                              ",
    "server=\n         Running simulation server hostname                              ",
    "keep_all=f\n      keep all the particles despite the selection     ",
    "times=all\n       Select time                                      ",
    "xmin=0.0\n        xmin box (for ramses input)                      ",
    "xmax=1.0\n        xmax box (for ramses input)                      ",
    "ymin=0.0\n        ymin box (for ramses input)                      ",
    "ymax=1.0\n        ymax box (for ramses input)                      ",
    "zmin=0.0\n        zmin box (for ramses input)                      ",
    "zmax=1.0\n        zmax box (for ramses input)                      ",
    "lmin=0.0\n        level min (for ramses amr input)                 ",
    "lmax=0.0\n        level max (for ramses amr input)                 ",
    "scale=1000.\n        ramses rescaling factor                          ",
    "vel=f\n           load velocity coordinates                        ",
    "disp_vel=f\n      display velocity vectors                         ",
    "blending=t\n      Activate blending colors                         ",
    "dbuffer=f\n       Activate OpenGL depth buffer                     ",
    "perspective=t\n   false means orthographic                         ",
    "bestzoom=t\n      automatic zoom                                   ",
    "auto_render=t\n   automatic rendering mode                         ",
    "play=f\n          automatically load and display next snapshot     ",
    "glsl=t\n          try to initialyze GLSL engine                    ",
    "ortho_range=6.0\n xy range if orthographic projection              ",
    "zoom=-14\n        zoom value                                       ",
    "xrot=0.0\n        rotation angle on X axis                         ",
    "yrot=0.0\n        rotation angle on Y axis                         ",
    "zrot=0.0\n        rotation angle on Z axis                         ",
    "xtrans=0.0\n      translation on X                                 ",
    "ytrans=0.0\n      translation on Y                                 ",
    "ztrans=0.0\n      translation on Z                                 ",
    "grid=t\n          Show grid                                        ",
    "nb_meshs=28\n     #meshs for the grid                              ",
    "mesh_size=1.0\n   grid's size of one mesh                          ",
    "xyg=t\n           display a grid in XY plan                        ",
    "yzg=f\n           display a grid in YZ plan                        ",
    "xzg=f\n           display a grid in XZ plan                        ",
    "cube=f\n          display a cube                                   ",
    "osd=t\n           Show On Screen Display                           ",
    "osdtime=t\n       Show time on OSD                                 ",
    "osdnbody=t\n      Show nbody on OSD                                ",
    "osdzoom=t\n       Show zoom on OSD                                 ",
    "osdrot=t\n        Show rotattion on OSD                            ",
    "osdtrans=t\n      Show transformation on OSD                       ",
    "osdtitle=t\n      Show title on OSD                                ",
    "osddata=t\n       Show data type on OSD                            ",
    "osd_set_title=\n  Set an explicit title on OSD                     ",
    "osdfs=13.0\n Size of OSD's font                                    ",
    "axis=t\n          display axis                                     ",
    "cb=t\n            display Color Bar (CB) on the screen             ",
    "cblog=f\n         display real or log of the physical value on CB  ",
    "cbloc=3\n         CB location, 0:top 1:right 2:bottom 3:left       ",
    "cbdigits=1\n      CB #digits                                       ",
    "cboffset=35\n     CB #offset pixels from the border location       ",
    "cbpw=0.03\n       CB size in percentage of the OpenGL windows width",
    "cbph=0.65\n       CB size in percentage of the OpenGL windows height",
    "cbfs=13\n   size of the fonts used to display CB             ",
    "com=t\n           centering according Center Of Mass               ", 
    "point=f\n         show particles as points                         ",
    "selphys=1\n        select physical quantity to display\n           "
    "             (1:density, 2:temperature, 3:pressure)                ",
    "minphys=\n        set minimal physical value                       ",
    "maxphys=\n        set maximal physical value                       ",
    "auto_ts=t\n       automatic texture size                           ",
    "texture=t\n       show particles as textures                       ",
    "texture_s=1.0\n   texture size of gaz particle                     ",
    "texture_a=1.0\n   texture alpha of gaz particle                    ",
    "cmapindex=0\n     Color map index                                  ",
    "psize=1.0\n       Set particles size                               ",
    "port=4000\n       Server's communication port                      ",
    "wsize=931\n       Windows's width size                             ",
    "hsize=750\n       Windows's height size                            ",
    "screenshot=\n     Screenshot name                                  ",
    "shot_ext=jpg\n    Screenshot's extension jpg|png                    ",
    "smooth_gui=t\n    if true it allows a smoother interactivity with  ",
    "                   the GUI, but it **double** the memory usage. \n",
    "VERSION="RELEASE_VERSION"\n    "__DATE__"  - JCL  compiled at <"__TIME__">      ",
    NULL
  };
  const char * usage="Interactive 3D OpenGL NBody simulation Snapshots rendering program";
  Q_IMPORT_PLUGIN(nemoplugin);
  Q_IMPORT_PLUGIN(ftmplugin);
  Q_IMPORT_PLUGIN(gadgetplugin);
  Q_IMPORT_PLUGIN(phigrapeplugin); 
  Q_IMPORT_PLUGIN(ramsesplugin); 
  Q_IMPORT_PLUGIN(listplugin);
  Q_IMPORT_PLUGIN(networkplugin);

// ============================================================================
//  The main program is here                                                   
int main(int argc, char *argv[])
{
  QApplication::setDesktopSettingsAware(true);
  QApplication app(argc, argv);
  if ( !QGLFormat::hasOpenGL() ) {
    qWarning( "This system has no OpenGL support. Exiting." );
    return -1;
  }
//   QGLFormat f;
//   f.setDirectRendering(false);
//   QGLFormat::setDefaultFormat(f);
  // Get screen coordinates
  QDesktopWidget  * desktop1 = QApplication::desktop();
  std::cerr << "# screens = " << desktop1->numScreens()<< "\n";
  std::cerr << "Is virtual desktop : " << desktop1->isVirtualDesktop() << "\n";
  QWidget  * desktop = desktop1->screen(desktop1->primaryScreen());

  // initialyze NEMO engine
  initparam(argv,const_cast<char**>(defv));
  // CAUTION !!! do not call getparam function after **MainWindow** object    
  // instantiation bc nemo engine get corrupted by io_nemo afterwhile the     
  // snapshot has been loaded.                                                
  const int wsize=getiparam((char *)"wsize");
  const int hsize=getiparam((char* )"hsize");
  std::string shot;
  bool interact=true;
  bool play=getbparam((char* )"play");
  if ( hasvalue((char *) "screenshot") )  {
    shot           = getparam((char *) "screenshot");
    interact       = false;
  }  else
    shot="";
  Q_INIT_RESOURCE(glnemo); // load resources
//  QPixmap pixmap(glnemo::GlobalOptions::RESPATH+"/images/glnemo2.png");
//  QSplashScreen splash(pixmap);
//  if (interact) {    
//    splash.show();
//    app.processEvents();
//  }
  glnemo::MainWindow main_win(RELEASE_VERSION); // main window object

  // compute window size
  const int x=((desktop->width()  - wsize)/2);
  const int y=((desktop->height() - hsize)/2);
  main_win.setGeometry(x,y,wsize,hsize);
  // move to the center of the screen
  main_win.move(x,y);

  
  if (interact) {
    main_win.show();
    
  }
  main_win.start(shot);
//  splash.finish(&main_win);
  finiparam();  // garbage collecting for nemo

  //if (interact) return app.exec();
  if (interact || play) 
    return app.exec();
}
