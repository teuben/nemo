// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2008                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           P�le de l'Etoile, site de Ch�teau-Gombert                         
//           38, rue Fr�d�ric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
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
// Nemo stuffs
#define _vectmath_h // put this statement to avoid conflict with C++ vector class
#include <nemo.h>

#include "mainwindow.h"
using namespace std;

#define RELEASE_VERSION "1.preview.2009-Jul-02"

// ============================================================================
// NEMO parameters                                                             
  const char * defv[] = {  
    "in=\n             Nemo input snapshot                                             ",
    "select=\n      Select particles using range operator separated               \n"
    "                   by a comma. E.g 0:999,1000:1999 would select two sets        \n"
    "                   of 1000 particles and give them a different color              ",
    "server=\n         Running simulation server hostname                              ",
    "range_visib=t\n   toggle visibility  for the particles selected via \'select=\' \n"
    "                   options. Can be usefull if you only want to display particles\n"
    "                   selected via \'select_list\' options(see below). In that case\n"
    "                   you should set \'f\'                                           ",
    "select_list=\n    Select particles from list of indexes            ",
    "keep_all=f\n      keep all the particles despite the selection     ",
    "comp=\n           select particles from component name, separated\n"
    "                   by a comma. E.g \"halo,disk,gas\" would select\n"
    "                   Halo, Disk and Gas components.\n                         ",
    "times=all\n       Select time                                      ",
    "vel=f\n           load velocity coordinates                        ",
    "disp_vel=f\n      display velocity vectors                         ",
    "blending=t\n      Activate blending colors                         ",
    "dbuffer=f\n       Activate OpenGL depth buffer                     ",
    "perspective=t\n   false means orthographic                         ",
    "bestzoom=t\n      automatic zoom                                   ",
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
    "osd=t\n           Show On Screen Display                           ",
    "part=f\n          show particles as points                         ",
    "texture=t\n       show particles as textures                       ",
    "texture_s=0.15\n  texture size of gaz particle                     ",
    "texture_ac=125\n  texture alpha color of gaz particle              ",
    "psize=1.0\n       Set particles size                               ",
    "port=4444\n       Server's communication port                      ",
    "wsize=931\n       Windows's width size                             ",
    "hsize=750\n       Windows's height size                            ",
    "screenshot=\n     Screenshot name                                  ",
    "anim_file=\n      Animation filename                               ",
    "anim_play=t\n     play animation ?                                 ",
    "anim_bench=t\n    play animation in benchmark mode ?               ",
    "smooth_gui=t\n    if true it allows a smoother interactivity with  ",
    "                   the GUI, but it **doubles** the memory usage. \n",
    "VERSION="RELEASE_VERSION"\n    "__DATE__"  - JCL  compiled at <"__TIME__">      ",
    NULL
  };
  const char * usage="Interactive 3D OpenGL NBody simulation Snapshots rendering program";
  Q_IMPORT_PLUGIN(listplugin);
  Q_IMPORT_PLUGIN(nemoplugin);
  Q_IMPORT_PLUGIN(ftmplugin);
  Q_IMPORT_PLUGIN(gadgetplugin);
  Q_IMPORT_PLUGIN(phigrapeplugin); 


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
  if ( hasvalue((char *) "screenshot") )  {
    shot           = getparam((char *) "screenshot");
    interact       = false;
  }  else
    shot="";
  Q_INIT_RESOURCE(glnemo); // load resources
  glnemo::MainWindow main_win(RELEASE_VERSION); // main window object

  // compute window size
  const int x=((desktop->width()  - wsize)/2);
  const int y=((desktop->height() - hsize)/2);
  main_win.setGeometry(x,y,wsize,hsize);
  // move to the center of the screen
  main_win.move(x,y);

  if (interact) main_win.show();
  main_win.start(shot);
  finiparam();  // garbage collecting for nemo

  if (interact) return app.exec();
}
