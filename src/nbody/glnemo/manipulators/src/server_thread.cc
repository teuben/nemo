// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// server_thread.cc                                                            |
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
#include "server_thread.h"
#include <string.h>
#define _vectmath_h // put this statement to avoid conflict with C++ vector class
#include <nemo.h>


#include <iostream>
#include <sys/times.h>
#include <assert.h>
#include <unistd.h>  // seems to be already included from <iostream> under gcc 3.3.1

#include "MessageBuffer.h"

#include <string>
#define LOCAL_DEBUG 0
#include "print_debug.h"

using namespace falcON;
//------------------------------------------------------------------------------
// Constructor                                                                  
// Initialyze mutex and arrays to store particles data                          
ServerThread::ServerThread(const std::string _sim_name,
			   const snapshot * S,
			   pthread_mutex_t * _mut,
			   pthread_cond_t  * _cond):GenericThread()
{
  // get data pointers
  my_snapshot = S;
  mut   = _mut;             // global mutex
  cond  = _cond;            // condition variable
  selected_3d = NULL;
  selected_index = NULL;
  sim_name = _sim_name;

  int nbody = my_snapshot->N_bodies();
  selected_index= new int[nbody];     // keep memory
  selected_3d   = new float[nbody*3]; // keep memory


}
//------------------------------------------------------------------------------
// Destructor                                                                   
// deallocate memory                                                            
ServerThread::~ServerThread()
{
  if (selected_3d) {
    delete[] selected_3d;
  }
  if (selected_index) {
    delete[] selected_index;
  }
}
//------------------------------------------------------------------------------
// kill:                                                                        
// if running kill the thread by closing its socket descriptor which launch an  
// exception.                                                                   
void ServerThread::kill()
{
  if (isRunning()) {
    close (newSd);
  }
}
//------------------------------------------------------------------------------
// start:                                                                       
// launch in parallel the "run()" method                                        
void ServerThread::start(int sock_fd)
{
  newSd = sock_fd;
  if (GenericThread::init()==0) {  // create the running thread
    GenericThread::start(); // start the run method
  }
}
//------------------------------------------------------------------------------
// run()::                                                                      
// this is the parallel running method                                          
void ServerThread::run()
{ 
  int len;
  clock_t t_start,t_stop;
  float tps;
  struct tms qq;
  float cps=1./(float)  (sysconf(_SC_CLK_TCK)); // Clock per second  
  bool stop     = false;
  bool first_connect = true;

  // Instantiate a new MessageBuffer object dealing with 
  // the new socket "newSd"                              
  MessageBuffer * serverMB=new MessageBuffer(newSd,8192*2);

  std::string hello("[gyrfalcon] "+sim_name);
  try {
    PRINT_D std::cerr <<"["<<getMyid()<<"]" << "Sending hello\n";
    serverMB->sendData(MessageBuffer::Hello,hello.length()+1,hello.c_str());
    PRINT_D std::cerr <<"["<<getMyid()<<"]"<< "Sending nbody ["<<my_snapshot->N_bodies()<<"] \n";
    serverMB->sendData(MessageBuffer::Nbody,1,(char *) &my_snapshot->N_bodies());
    PRINT_D std::cerr <<"["<<getMyid()<<"]"<< "end sending nbodies...\n";
  }
  catch (int n) {
    switch(n) {
    case -1 : PRINT_D std::cerr <<"["<<getMyid()<<"]"
				<< "ServerThread::run() - WARNING: Catch error during SEND\n";
      delete serverMB;
      stop    =true; // not necessary to keep running
      break;
    default :
      assert(1);
    }
  }      
  // Control variable for select_pos array
  char * selected_string;
  bool lock_mut=false;;
  float mytime;
  while (!stop) {
    try {
      serverMB->recvData();
      switch (int get_buf=*(serverMB->getTagBuffer())) {  // Get Tag buffer
      case MessageBuffer::Select :     // Select string
      case MessageBuffer::SelectV :    // Select string
	// do something with the selected string
	selected_string = serverMB->getDatBuffer() ;
	PRINT_D std::cerr <<"["<<getMyid()<<"]"<< "Got Select from client [" << selected_string << "]\n";
	parseSelectedString(selected_string);

	// ********** >> Protect this area with conditionnal mutex
	PRINT_D std::cerr <<"["<<getMyid()<<"]"<< "before LOCK A\n\n";
	pthread_mutex_lock(mut);
	PRINT_D std::cerr <<"["<<getMyid()<<"]"<< "After LOCK A\n\n";
	lock_mut = true;
	if (first_connect) {      // this is the first connexion, we 
	  first_connect = false;  // can get the data without waiting
	} 
	else { // here I am waiting authorization from main thread (updateData())
	  PRINT_D std::cerr <<"["<<getMyid()<<"]"<< "Waiting authorization...\n";
	  pthread_cond_wait( ServerThread::cond, mut);
	  PRINT_D std::cerr <<"["<<getMyid()<<"]"<< "Got authorization\n";
	}
	PRINT_D std::cerr <<"["<<getMyid()<<"]"<< "Sending time => [" << my_snapshot->time() << "]\n";
	mytime=my_snapshot->time();
	serverMB->sendData(MessageBuffer::Time,1,(char *) &(mytime));                   // send time     
	t_start=times(&qq);
	PRINT_D std::cerr <<"start = " << t_start <<"\n";
	PRINT_D std::cerr <<"["<<getMyid()<<"]"<< "Sending selected_nbody [" << selected_nbody << "] and pos\n";
	fill3DArray(); // put positions in 3D array
	serverMB->sendData(MessageBuffer::Pos,selected_nbody*3,(char *) selected_3d);    // send positions 
	if ( get_buf == MessageBuffer::SelectV) {
	  PRINT_D std::cerr <<"["<<getMyid()<<"]"<< "Sending selected_nbody [" << selected_nbody << "] and vel\n";
	  fill3DArray(true); // put velocities in 3D array
	  serverMB->sendData(MessageBuffer::Vel,selected_nbody*3,(char *) selected_3d);  // send velocities
	}
	t_stop=times(&qq);
	PRINT_D std::cerr <<"stop = " << t_stop <<"\n";
	pthread_mutex_unlock(mut);
	lock_mut = false;  
	// ********** << Protect this area with conditionnal mutex

	len = sizeof(int) * selected_nbody * 3;
	tps=((float)(t_stop-t_start));
	tps=( tps<=0 ? 1 : tps);
	PRINT_D std::cerr <<"["<<getMyid()<<"]"<<"Transfer Time ["<< (tps )
			  << "] seconds =  ["<< len/(tps*cps)/1024/1024 << "] MB/s\n"; 
	
	break;
      case MessageBuffer::Stop :
	stop = 1;
	PRINT_D std::cerr <<"["<<getMyid()<<"]"<< "ServerThread::run() - Client ask me to stop the connection....";
	delete serverMB;
	close(newSd);
	break;
      default :
	stop = 1;
	std::cerr <<"["<<getMyid()<<"]"<< "ServerThread::run() - Unknown received data [" 
		  << *(serverMB->getTagBuffer()) << "]\n";
	std::cerr <<"["<<getMyid()<<"]"<< "ServerThread::run() - Closing connection with client..\n";
	delete serverMB;
	close(newSd);
	break;
      }
      //std::exit(1);
    }
    // catch exceptions
    catch (int n) {  // catch errors thrown by SEND and RECV
      switch(n) {
      case -1 : PRINT_D std::cerr <<"["<<getMyid()<<"]"
				  << "ServerThread::run() - WARNING: Catch error during SEND\n";
	break;
      case -2 : PRINT_D std::cerr <<"["<<getMyid()<<"]"
				  << "ServerThread::run() - WARNING: Catch error during RECV\n";
	break;
      case -3 : PRINT_D std::cerr <<"["<<getMyid()<<"]"
				  << "ServerThread::run() - WARNING: Catch error during RECV, nbuffer=0!\n";
	break;
      case -40: PRINT_D std::cerr <<"["<<getMyid()<<"]"
				  << "ServerThread::run() - WARNING: "
	  "Catch error during ServerThread::parseSelectedString\n";
	break;
      default :
	assert(1);
      }
      delete serverMB;  // destroy object
      close(newSd);
      stop = true;      // get out of loop
      if (lock_mut) {   // must release the mut, bc it was caught
	PRINT_D std::cerr <<"["<<getMyid()<<"]"
			  << "ServerThread::run() - WARNING: release mutex trapped in a Try and Catch !!!!\n";
	pthread_mutex_unlock(mut );
      }
    }
  }

}
//------------------------------------------------------------------------------
// fillPosVelArray:                                                             
// fill Positins and velocities array                                           
int ServerThread::fill3DArray(bool vel)
{
  int nbody=my_snapshot->N_bodies();

  // fill selected_pos array according to selected_index
  PRINT_D std::cerr << "ServerThread::parseSelectedString, filling select_pos [" 
		    <<  selected_nbody << "]\n";
  int nbody_out=0;
  for (int i=0; i<selected_nbody; i++) {
    int p_index=selected_index[i];
    if (p_index < nbody) {
      body b_current = my_snapshot->bodyNo(p_index);
      vect x;
      if ( vel ) {
	x=b_current.vel();
      }
      else {
	x=b_current.pos();
      }
      selected_3d[nbody_out*3  ] = x[0];//my_snapshot->pos(p_index)[0]; // x
      selected_3d[nbody_out*3+1] = x[1];//my_snapshot->pos(p_index)[1]; // y
      selected_3d[nbody_out*3+2] = x[2];//my_snapshot->pos(p_index)[2]; // z
      nbody_out++;
    }
  }
  selected_nbody = nbody_out;
  //  std::cerr << "In parseSelected npart =["<< npart << "]\n";
  PRINT_D std::cerr << "End of ServerThread::parseSelectedString\n";
  return selected_nbody;
  
}
//------------------------------------------------------------------------------
// parseSelectedString:                                                         
// parse the selected_string to fill selected_pos array                         
int ServerThread::parseSelectedString(char * selected_string)
					
{
  int nbody=my_snapshot->N_bodies();

  PRINT_D std::cerr << "ServerThread::parseSelectedString ["<< selected_string << "]\n";
  if (  strcmp(selected_string,"all")) {
    selected_nbody = nemoinpi(selected_string, selected_index, nbody);
    if (selected_nbody <=0 ) {
      std::cerr << "nemoinpi = [" << selected_string << "] selected_nbody = "
		<< selected_nbody <<"\n";
      throw(-40);
    }
  }
  else {  // select all the particles
    selected_nbody=nbody;
    PRINT_D std::cerr << "ServerThread::parseSelectedString , select all[" 
		      << selected_nbody << "]\n";
    PRINT_D std::cerr << "Pointer selected_index = [" << selected_index<< "]\n";
    for (int i=0; i<nbody; i++) {
      selected_index[i] = i;
    }
  }

}
//------------------------------------------------------------------------------
