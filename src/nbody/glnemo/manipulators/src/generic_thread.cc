// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// generic_thread.cc                                                           |
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
// implementation of the GenericThread class                                    
//                                                                              
// This class implement a generic way to use simple program                     
// using thread                                                                 
//------------------------------------------------------------------------------

#include <iostream>
#include <assert.h>
#include "generic_thread.h"

#define LOCAL_DEBUG 0
#include "print_debug.h"
using namespace std;

//------------------------------------------------------------------------------
// Constructor                                                                  
// create a mutex to control the thread                                         
// initialyze control variables                                                 
GenericThread::GenericThread()
{
  int status;

  // initialyse thread mutex control
  status = pthread_mutex_init(&mutex,NULL );

  if (status != 0 ) {
    std::cerr << "GenericThread::GenericThread() - Error in [pthread_mutex_init], aborted...\n";
    std::exit(1);
  }  

  finished = false;
  running  = false;
  started  = false;
}
//------------------------------------------------------------------------------
// Destructor                                                                   
// free mutex ressources                                                        
GenericThread::~GenericThread()
{
  pthread_mutex_destroy(&mutex);
}
//------------------------------------------------------------------------------
// init():                                                                      
// lock the mutex                                                               
// create the thread                                                            
int GenericThread::init()
{
  int status,status_create;  
  assert(started  == false); 
  // lock the mutex to prevent the thread to start (lock the thread)
  status = mutexLock();
  if (status != 0 ) {
    std::cerr << "GenericThread::init() - ERROR in mutexLock, aborted..\n";
    std::exit(1);
  }
  // spwan thread
  status_create = pthread_create(&id, NULL,threadStartup, this);
  if (status_create != 0 ) {
    std::cerr << "GenericThread::init() - ERROR in [pthread_create]...\n";
    std::cerr << "Status = [" << status_create << "]\n";
    //    std::exit(1);
  }  
  if (status_create==0 ){
    started = true; // thread has been started
  } 
  else { // problem occur during thread creation
    status = mutexUnLock(); // release mutex
    if (status != 0 ) {
      std::cerr << "GenericThread::init() - ERROR  mutexUnLock, aborted...\n";
      std::exit(1);
    }    
  }
  return status_create;
}

//------------------------------------------------------------------------------
// start():                                                                     
//  release mutex to allow "run()" method to start                              
void GenericThread::start()
{
  int status;
  // free the mutex
  status = pthread_mutex_unlock(&mutex);
  if (status != 0 ) {
    cerr << "[GenericThread::start()] Ploblem occured in mutex unlock, aborted...\n";
    std::exit(1);
  }
  running = true;
}
//------------------------------------------------------------------------------
// mutexLock():                                                                 
// Lock the mutex                                                               
int GenericThread::mutexLock()
{
   int status = pthread_mutex_lock(&mutex);
   return status;
}
//------------------------------------------------------------------------------
// mutexUnLock():                                                               
// unlock the mutex                                                             
int GenericThread::mutexUnLock()
{
   int status = pthread_mutex_unlock(&mutex);
   return status;
}
//------------------------------------------------------------------------------
// join():                                                                      
// wait for the end of the thread                                               
//------------------------------------------------------------------------------
int GenericThread::join()
{
  int status=pthread_join(id,NULL);
  return status;
}
//------------------------------------------------------------------------------
// threadStartup(void * arg):                                                   
// I: arg =>  pointer to the derived class                                      
// (helper function)                                                            
// this is the running // thread                                                
// wait for the mutex                                                           
// call the "run" method                                                        
//------------------------------------------------------------------------------
void * threadStartup(void * arg)
{
  int status;
  
  GenericThread * t = (GenericThread *) arg;

  // get the mutex to start the "run" method
  PRINT_D std::cerr << "[threadStartup()] Waiting for mutex.....\n";
  status = t->mutexLock();
  if (status != 0 ) {
    std::cerr << "[threadStartup()] Ploblem occured in mutex lock, aborted...\n";
    std::exit(1);
  }
  PRINT_D std::cerr << "[threadStartup()] got the mutex !!!! \n";

  // start "run" method
  t->run();

  // release the mutex
  status = t->mutexUnLock();
  if (status != 0 ) {
    std::cerr << "[threadStartup()] Ploblem occured in mutex unlock, aborted...\n";
    std::exit(1);
  }
  PRINT_D std::cerr << "[threadStartup()] released the mutex !!!! \n";

  // reset variables
  t->setFinished(true);
  t->setRunning(false);

  // I am finished
  int ret;
  pthread_exit(&ret);
}
//------------------------------------------------------------------------------
