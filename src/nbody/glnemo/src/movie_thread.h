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
// MovieThread class definition                                                
//                                                                             
// Grab GL frame buffer to make movie                                          
// ============================================================================
#ifndef MOVIE_THREAD_H
#define MOVIE_THREAD_H

#include <qthread.h>
#include <qimage.h>
#include <qobject.h>
#include "glbox.h"
class MovieThread : public QThread, QObject {
  //Q_OBJECT
public:
  MovieThread(int fps, GLBox * _glbox);
  ~MovieThread();
  
  void run();
  void stop();

private:
  int my_timer, play_timer;
  GLBox * glbox;
  vector <QImage> store_image;
  bool enable;
  
  void grabFrameBuffer();
  void callMe();    
  void timerEvent( QTimerEvent *e );
};
#endif
// ============================================================================
