// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2006                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
//                                                                             
// Glnemo is  an interactive 3D plotting program  , using OpenGL and trolltech 
// QT API. It allows to plot in 3D the positions x,y,z of a NEMO snapshot.     
//                                                                             
// This software is governed by the CeCILL  license under French law and       
// abiding by the rules of distribution of free software.  You can  use,       
// modify and/ or redistribute the software under the terms of the CeCILL      
// license as circulated by CEA, CNRS and INRIA at the following URL           
// "http://www.cecill.info".                                                   
//                                                                             
// As a counterpart to the access to the source code and  rights to copy,      
// modify and redistribute granted by the license, users are provided only     
// with a limited warranty  and the software's author,  the holder of the      
// economic rights,  and the successive licensors  have only  limited          
// liability.                                                                  
//                                                                             
// In this respect, the user's attention is drawn to the risks associated      
// with loading,  using,  modifying and/or developing or reproducing the       
// software by the user in light of its specific status of free software,      
// that may mean  that it is complicated to manipulate,  and  that  also       
// therefore means  that it is reserved for developers  and  experienced       
// professionals having in-depth computer knowledge. Users are therefore       
// encouraged to load and test the software's suitability as regards their     
// requirements in conditions enabling the security of their systems and/or    
// data to be ensured and,  more generally, to use and operate it in the       
// same conditions as regards security.                                        
//                                                                             
// The fact that you are presently reading this means that you have had        
// knowledge of the CeCILL license and that you accept its terms.              
// ============================================================================
//                                                                             
// main program                                                                
//                                                                             
// ============================================================================

#include <iostream>
#include <qapplication.h>
#include <qgl.h>

// Nemo stuffs
#define _vectmath_h // put this statement to avoid conflict with C++ vector class
#include <nemo.h>

#include "globjwin.h"

using namespace std;
// ============================================================================
// NEMO parameters                                                             
::string defv[] = {  // use `::'string because of 'using namespace std'
  "in=\n             Nemo input snapshot                                             ",
  "server=\n         Running simulation server hostname                              ",
  "select=all\n      Select particles from the : range operator , separeted\n"
  "                   by a comma. E.g 0:999,1000:1999 would select two sets\n"
  "                   of 1000 particles and give them a different color              ",
  "range_visib=t\n   toggle visibility  for the particles selected via \'select=\'\n"
  "                   options. Can be usefull if you only want to display particles\n"
  "                   selected via \'select_list\' options(see below). In that case\n"
  "                   you should set \'f\'                                           ",
  "select_list=\n    Select particles from list of indexes          ",
  "times=all\n       Select time                                    ",
  "vel=t\n           load velocity coordinates                      ",
  "disp_vel=f\n      display velocity vectors                       ",
  "blending=t\n      Activate blending colors                       ",
  "dbuffer=f\n       Activate OpenGL depth buffer                   ",
  "perspective=t\n   false means orthographic                       ",
  "bestzoom=t\n      automatic zoom                                 ",
  "play=f\n          automatically load and display next snapshot   ",
  "ortho_range=6.0\n xy range if orthographic projection            ",
  "zoom=-14\n        zoom value                                     ",  
  "xrot=0.0\n        rotation angle on X axis                       ",
  "yrot=0.0\n        rotation angle on Y axis                       ",
  "zrot=0.0\n        rotation angle on Z axis                       ",
  "xtrans=0.0\n      translation on X                               ",
  "ytrans=0.0\n      translation on Y                               ",
  "ztrans=0.0\n      translation on Z                               ",
  "grid=f\n          Show grid                                      ",
  "gas=f\n           gas like particle effect                       ",
  "texture_s=0.15\n  texture size of gaz particle                   ",
  "texture_ac=125\n  texture alpha color of gaz particle            ",
  "psize=1.0\n       Set particles size                             ",
  "port=4444\n       Server's communication port                    ",
  "wsize=925\n       Windows's width size                           ",
  "hsize=685\n       Windows's height size                          ",
  "screenshot=\n     Screenshot name                                ",
  "VERSION=0.93\n    "__DATE__"  - JCL  compiled at <"__TIME__">   ",
  NULL
};
::string usage="3D OpenGL Nemo Snapshot rendering";

// ============================================================================
//  The main program is here                                                   
int main( int argc, char **argv )
{
  QApplication::setColorSpec( QApplication::CustomColor );
  QApplication a(argc,argv);			

  if ( !QGLFormat::hasOpenGL() ) {
    qWarning( "This system has no OpenGL support. Exiting." );
    return -1;
  }
  // Get screen coordinates
  QDesktopWidget  * desktop = QApplication::desktop();
  
  // initialyze NEMO engine
  initparam(argv,defv);

  // CAUTION !!! do not call getparam function after **GLObjectWindow** object
  // instantiation bc nemo engine get corrupted by io_nemo afterwhile the     
  // snapshot has been loaded.                                                
  const int wsize=getiparam("wsize");
  const int hsize=getiparam("hsize");
  ::string screenshot=NULL;
  bool has_screenshot=false;
  if ( hasvalue("screenshot")) {
    screenshot=getparam("screenshot");
    has_screenshot=true;
  }
  // Create widget window
  GLObjectWindow w(0,"glnemo");

  // move window to the center of the screen
  w.setGeometry(desktop->width() /2 - (wsize/2),
                desktop->height()/2 - (hsize/2),
                wsize,hsize);

  a.setMainWidget( &w );
  w.show();
  
  w.resize(wsize-1,hsize-1); // trick to update GL buffer with right W & H
  w.glbox->enable_s=true;
  if (has_screenshot) {
    w.showMinimized();
    w.takeScreenshot(screenshot);
  }
  finiparam();  // garbage collecting for nemo
  if (! has_screenshot) {
    return a.exec();
  }
}
// ============================================================================
