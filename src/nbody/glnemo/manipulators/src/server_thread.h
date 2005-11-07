// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// server_thread.h                                                             |
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

#ifndef SERVER_THREAD_H
#define SERVER_THREAD_H
#include <iostream>

#include <body.h>  // from gyrfalcon distribution

#include "generic_thread.h"
#include <string>

//using namespace std;

class ServerThread : public GenericThread
{
 public:
  //
  ServerThread( const std::string _sim_name, const falcON::snapshot * S,
		pthread_mutex_t  * _mut, pthread_cond_t  * _cond);
  ~ServerThread();
  //
  void start(int);
  void run();
  
  void kill();

 private:
  std::string sim_name;
  const falcON::snapshot * my_snapshot;
  float * timu;
  float * selected_pos;   // store x,y,z particles    
  int     selected_nbody; // #bodies selected         
  int   * selected_index; // bodies's indexes selected
  pthread_mutex_t  * mut;
  pthread_cond_t   * cond;
  pthread_mutex_t  condition_mutex;
  int    newSd;

  int parseSelectedString(char * select_string);
};

#endif

//------------------------------------------------------------------------------
