// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
//                                                                             
// MovieThread class implementation                                            
//                                                                             
// Grab GL frame buffer to make movie                                          
// ============================================================================
#include "movie_thread.h"
#include <qimage.h>
#include <qmessagebox.h>
// ============================================================================
// Constructor
MovieThread::MovieThread(int fps, GLBox * _glbox)
 : QThread()
{
  my_timer = 1000/fps;
  glbox = _glbox;
}
// ============================================================================
// Destructor
MovieThread::~MovieThread()
{
  
}

// ============================================================================
// run method
void MovieThread::run()
{
  play_timer = startTimer( my_timer ); 
}
// ============================================================================
// stop method
void MovieThread::stop()
{
 killTimer( play_timer );
 terminate();
}
// ============================================================================
//
void MovieThread::timerEvent( QTimerEvent *e )
{
  static int no;
  if (e->timerId() == play_timer) {    
    QImage img=glbox->grabFrameBuffer(); // grab OpenGL buffer
    char shot[100];
    sprintf(shot,"frame.%05d.png",no);
    std::cerr <<  ">> shot = " << shot << "\n";
#if 1    
    if ( !img.save( shot, "PNG" ) ) { // save the image
      //QMessageBox::warning( this, "Save failed", "Error saving file" );
    }
#endif    
    no++;
  }
}
// ============================================================================
