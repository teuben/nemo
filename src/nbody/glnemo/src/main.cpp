// ============================================================================
//  $Id $
//                                                                             
// main program                                                                
//                                                                             
//                                                                             
// ============================================================================

#include "globjwin.h"
#include <qapplication.h>
#include <qgl.h>

/*
  The main program is here. 
*/

int main( int argc, char **argv )
{
  QApplication::setColorSpec( QApplication::CustomColor );
  QApplication a(argc,argv);			

  if ( !QGLFormat::hasOpenGL() ) {
    qWarning( "This system has no OpenGL support. Exiting." );
    return -1;
  }

  GLObjectWindow w(0,"glnemo",argc,argv);
  w.resize( 750, 750 );
  a.setMainWidget( &w );
  w.show();
  return a.exec();
}
