// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2007                                  
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
#include <unistd.h>

// ============================================================================
// Constructor                                                                 
MovieThread::MovieThread(int fps, GLBox * _glbox)
 : QThread()
{
  my_timer = 1000/fps;   // 1 second = 1 000 000 microseconds
  glbox = _glbox;
  enable=false;

}
// ============================================================================
// Destructor                                                                  
MovieThread::~MovieThread()
{
  //store_image.clear();
}
// ============================================================================
// MovieThread::timerEvent()                                                   
// grab open gl frame buffer at each timer event                               
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
// MovieThread::run()                                                          
// method launch in parallel in a thread                                       
void MovieThread::run()
{
#if 0
  enable=true;
  while (enable) {
      this->grabFrameBuffer();
      //callMe();
      usleep(my_timer);
  }
#else
  play_timer = startTimer( my_timer );
  
#endif  
}
// ============================================================================
// MovieThread::stop()                                                         
// stop the running thread                                                     
void MovieThread::stop()
{
#if 0
 enable=false;     // stop running loop       
 usleep(my_timer); // wait a timer delay      
 terminate();      // terminate running thread
 for (unsigned int i=0;i<store_image.size(); i++) {
  char shot[100];
  sprintf(shot,"frame.%05d.png",i);
  store_image[i].save(shot,"PNG");  
 }
 store_image.erase(store_image.begin(),store_image.end()); 
#else
  killTimer( play_timer );
  terminate();
#endif 
}
// ============================================================================
// MovieThread::grabFrameBuffer()                                              
// grab opengl frame buffer                                                    
void MovieThread::grabFrameBuffer()
{
  static int no;
  QImage img=glbox->grabFrameBuffer(TRUE); // grab OpenGL buffer
  char shot[100];
  sprintf(shot,"frame.%05d.png",no);
  std::cerr <<  ">> shot = " << shot << "\n";
  store_image.push_back(img);

  no++;
}
// ============================================================================
