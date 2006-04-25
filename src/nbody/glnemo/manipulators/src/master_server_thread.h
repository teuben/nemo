// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// master_server_thread.h                                                      |
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
// main thread which manage the others concurrent server thread                |
// ----------------------------------------------------------------------------+

#ifndef MASTER_SERVER_THREAD_H
#define MASTER_SERVER_THREAD_H
#include <iostream>
#include <pthread.h>
#include <string>

#include "generic_thread.h"
#include "server_thread.h"

#include <body.h>  // from gyrfalcon distribution
//#include <public/nbody.h>


using namespace std;
#define NB_SERVER 2
class MasterServerThread : public GenericThread {

 public:
  //
  MasterServerThread( const std::string _sim_name,
		      const int _port, const int _max_port, const falcON::snapshot * S);

  ~MasterServerThread();
  
  void updateData(const falcON::snapshot * S);

 private:
  std::string sim_name;
  const falcON::snapshot * my_snapshot;
  pthread_mutex_t   condition_mutex;
  pthread_cond_t    condition_cond;
  int port, max_port;
  ServerThread * server_t[NB_SERVER];
  bool stoploop;
  //
  void start();
  void run();

};

#endif
// -----------------------------------------------------------------------------
