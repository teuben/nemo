// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// master_server_thread.cc                                                     |
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

#include <iostream>
#include <sys/times.h>
#include <assert.h>
#include "master_server_thread.h"
#include "MessageBuffer.h"
#include "ServerClientDef.h"
#include "SocketServer.h"


using namespace falcON;
#define LOCAL_DEBUG 0
#include "print_debug.h"


// -----------------------------------------------------------------------------
// Constructor                                                                  
// Initialyze mutex, launch master server thread then return to gyrfalcON       
MasterServerThread::MasterServerThread(const std::string _sim_name,
				       const int _port, const int _max_port, 
				       const snapshot * S):GenericThread()
{
  // get data pointers
  port       =  _port;
  max_port   =  _max_port;
  sim_name   =  _sim_name;
  my_snapshot = S;

  // Initialyze mutex
  int status = pthread_mutex_init(&condition_mutex,NULL);
  if (status != 0 ) {
    std::cerr << "MasterServerThread::MasterServerThread() - ERROR " 
	      << "during [pthread_mutex_init], aborted...\n";
    std::exit(1);
  }
  // Initialyze condiiton mutex
  status = pthread_cond_init(&condition_cond,NULL );
  if (status != 0 ) {
    std::cerr << "MasterServerThread::MasterServerThread() - ERROR "
	      << "during [pthread_cond_init], aborted...\n";
    std::exit(1);
  }
  
  start(); // Start parallel task
}
// -----------------------------------------------------------------------------
// Destructor                                                                   
MasterServerThread::~MasterServerThread()
{
  stoploop = true;          // terminate the Master thread  
  for (int i=0; i<NB_SERVER; i++) {
    if ( server_t[i]->isRunning()) {
      std::cerr << 
	"MasterServerThread::~MasterServerThread() -"
	"killing server ["<<i<<"] thread...\n";
      server_t[i]->kill(); // kill the running server thread
    }    
    delete server_t[i];
  }  
}
// -----------------------------------------------------------------------------
// updateData:                                                                  
// at each time step, the manipulator call this method and the master thread    
// can inform the slave server thread that a new time step is ready to transmit 
void MasterServerThread::updateData(const snapshot * S)
{
  //my_snapshot = S; // update data from snapshot
  //std::cerr << "Update time [" << my_snapshot->time() << "]\n";
  pthread_cond_broadcast( &condition_cond ); // grant autorization to Thread servers
}
// -----------------------------------------------------------------------------
// start:                                                                       
// create/launch the thread and return                                          
void MasterServerThread::start()
{
  init();                 // create the running thread
  GenericThread::start(); // start the run method
}

// -----------------------------------------------------------------------------
// run:                                                                         
// this is the running Master Server Thread                                     
// - listen the communication port for a new glnemo client                      
// - launch a slave server thread which will communicate with the new glnemo    
//    client                                                                    
void MasterServerThread::run()
{
  /* share the job */
  //  SocketServer * my_server = new SocketServer(SOCK_STREAM,SERV_PORT,BACKLOG);
  SocketServer * my_server = new SocketServer(SOCK_STREAM,port,max_port,BACKLOG);

  //ServerThread * server_t[NB_SERVER];

  // initialyse servers
  for (int i=0; i<NB_SERVER; i++) {
    server_t[i] = new ServerThread(sim_name,my_snapshot,&condition_mutex,&condition_cond);
  }
  stoploop = false;

  while (!stoploop) {
    PRINT_D std::cerr << "\n*** MasterServerThread::run() - Waiting for a new client...***\n";
    int newSd=my_server->acceptClient();
    PRINT_D std::cerr << "\n*** MasterServerThread::run() => new client accepted\n";

    if ( ! stoploop) {
      // look for a free server
      bool find=false;
      for (int i=0; (i<NB_SERVER && !find) ; i++) {
	if (! server_t[i]->isStarted() ) {  // server not started
	  find=true;
	  server_t[i]->start(newSd);        // launch new server 
	  if (! server_t[i]->isStarted()) { // failed to start   
	    std::cerr << "MasterServerThread::run() - WARNING Failed to start server_thread\n";
	    close(newSd);
	  }
	}
      }
      if (!find) {
	std::cerr << "\n\n*** MasterServerThread::run() - TOO MANY SERVERS RUNNING ! ***\n\n";
	close(newSd);
      }
    }
  }    
}
// -----------------------------------------------------------------------------

