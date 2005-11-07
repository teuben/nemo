// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// generic_thread.h                                                            |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Jean-Charles LAMBERT - 2005                                       |
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      |
// address:  Dynamique des galaxies                                            |
//           Laboratoire d'Astrophysique de Marseille                          |
//           2, place Le Verrier                                               |
//           13248 Marseille Cedex 4, France                                   |
//           CNRS U.M.R 6110                                                   |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                              
// definition of the GenericThread class                                        
//                                                                              
// This class implement a generic way to use simple program                     
// using thread                                                                 
//------------------------------------------------------------------------------
#ifndef GENERIC_THREAD_H
#define GENERIC_THREAD_H

#include <pthread.h>

extern "C" void *threadStartup(void *); // helper method which    
                                        // will start the thread  
                                        // via run(). It's used   
                                        // to build the interface 
                                        // between C an C++       

class GenericThread  {
 public:
  GenericThread();
  virtual ~GenericThread();
  
  int  init();         // create the thread                 
  void start();        // start the virtual "run"ning thread
  void exit();         // exit                              
  int join();          // wait the end of the running thread
  bool isFinished() { return finished;} ;   
  bool isRunning()  { return running ;} ;
  bool isStarted()  { return started ;} ;
  void setFinished(const bool b) { finished = b; if (b) started=false;} ;
  void setRunning(const bool b) { running = b;} ;
  int mutexLock();     // lock the mutex  
  int mutexUnLock();   // unlock the mutex
  int getMyid() { return id;};
  // "run()" is a virtual method which must be implemented
  // in the class which derive GenericThread class        
  virtual void  run() { std::cerr << "In GenericThread::run()\n";};

 private:
  pthread_t id;            // pthread id                       
  pthread_mutex_t mutex;   // mutex for thread control         
  bool started,            // true if thread has been started   
       finished,           // true is thread has been terminated
       running;            // true is thread is running         

};
#endif
//------------------------------------------------------------------------------
