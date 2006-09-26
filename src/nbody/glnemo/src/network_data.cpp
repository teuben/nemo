// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2006                                  
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
// implementation of the NetworkData Class                                     
//                                                                             
//                                                                             
// ============================================================================

extern "C" { // Must include NEMO's header before calling [using namespace std]
             // to avoid confict on 'string' in <stdinc.h>
#include <stdinc.h>
}
#include "network_data.h"
#include "virtual_particles_select.h"
#include "particles_range.h"
#include <iostream>
#include <assert.h>

#define LOCAL_DEBUG 0
#include "print_debug.h"
using namespace std;

// ============================================================================
// NetworkData constructor                                                     
// Connect a socket to a server                                                
// create a message buffer for communication                                   
NetworkData::NetworkData(const char * hostname,const int _port)
{
  int port = _port; // the communication port  
  // get a new client instance
  client = new SocketClient(hostname,SOCK_STREAM,port);
  // connect to the server
  sock_fd = client->connectServer();
  
  clientMB = NULL;
  if (sock_fd == -1) {   // unable to connect
    is_connected = FALSE;
  } else {
    is_connected = TRUE;
    // create a new message buffer interface
    clientMB= new MessageBuffer(sock_fd ,8192*2);
    PRINT_D cerr << "Adress clienMB=" << clientMB << "\n";
  }
  
  // allocate memory for pointers
  //nbody = new int;
  //timu  = new float;
  part_data->allocVar(); // allocate nbody, timu
  // Initialize global class data
  *part_data->nbody = -1;
  *part_data->timu  = -1.0;
  is_parsed = FALSE;
  is_loading_thread = FALSE;
  is_new_data_loaded = FALSE;
}
// ============================================================================
// Destructor                                                                  
NetworkData::~NetworkData()
{
  if (client) {
    PRINT_D std::cerr << "NetworkData killing client\n";
    delete client;
    client = NULL;
  } 
  if (clientMB) {
    PRINT_D std::cerr << "NetworkData killing clientMB\n";
    delete clientMB;
    clientMB = NULL;
  }
  delete part_data;
  //delete nbody;
  //delete timu;
}
// ============================================================================
// NetworkData::loadPos()                                                      
// get from the server new positions data                                      
int NetworkData::loadPos(ParticlesSelectVector * psv, const bool vel)
{
  int n;   // happy red hat
  if (n) ; // remove compiler warning
  
  if (clientMB) {
    try {
      if (vel) {
        // send se lected string to the server
        clientMB->sendData(MessageBuffer::SelectV,qselect_part.length()+1,qselect_part);
      } 
      else {
        clientMB->sendData(MessageBuffer::Select,qselect_part.length()+1,qselect_part);  
      }
      
      // receive time from server
      clientMB->recvData();
      if (*(clientMB->getTagBuffer()) !=  MessageBuffer::Time) {
        throw (-4); // Bad protocol error
      }
      // convert buffer to TIME value
      *part_data->timu =  *((float *) (clientMB->getDatBuffer()));
      
      // receive position (correct nbody according selected range)
      clientMB->recvData();
      if (*(clientMB->getTagBuffer()) !=  MessageBuffer::Pos) {
        throw (-4); // Bad protocol error
      } 
      int nbytes=*((int *) (clientMB->getSizBuffer()));  // #bytes transferred
      int new_body = (int) ((float) (nbytes) / 3.0 / sizeof(float));
      assert (new_body <= *part_data->nbody);
      
      // allocate memory to store positions           
      if ( ! part_data->pos ||               // not yet allocated
            new_body > *part_data->nbody) {   // new_body bigger 
        part_data->pos = (float *) part_data->mallocate((char *) part_data->pos, sizeof(int) * new_body * 3, true);
      }
      *part_data->nbody = new_body;
      memcpy((float *) part_data->pos, (char *) clientMB->getDatBuffer(), nbytes);  
      PRINT_D cerr << "New nbody = " << *part_data->nbody << "\n";
        
      if (vel) {
        // receive velocities (correct nbody according selected range)
        clientMB->recvData();
        if (*(clientMB->getTagBuffer()) !=  MessageBuffer::Vel) {
          throw (-4); // Bad protocol error
        } 
        int nbytes=*((int *) (clientMB->getSizBuffer()));  // #bytes transferred
        int new_body = (int) ((float) (nbytes) / 3.0 / sizeof(float));
        assert (new_body <= *part_data->nbody);
        
        // allocate memory to store velocities           
        if ( ! part_data->vel ||               // not yet allocated
              new_body > *part_data->nbody) {   // new_body bigger 
          part_data->vel = (float *) part_data->mallocate((char *) part_data->vel, sizeof(int) * new_body * 3, true);
        }
        *part_data->nbody = new_body;
        memcpy((float *) part_data->vel, (char *) clientMB->getDatBuffer(), nbytes);
        
        // compute velocity vector norm
        part_data->computeVelNorm();
      }
      //computeCooMax();
      if (! is_parsed) {  // selected range string not yet parsed
        is_parsed = TRUE; // parse only one time
        int nobject=VirtualData::fillParticleRange(psv,full_nbody,
                          qselect_part);
        //VirtualParticlesSelect * vps = new VirtualParticlesSelect();
        //nobject=vps->storeParticles<ParticlesRange>(psv,full_nbody,qselect_part);
        //int nobject=vps->A<int>();
        if (nobject); // do nothing (remove compiler warning)
      }
      if (! is_loading_thread) {
        // a running Thread MUST not emit data to the GLBox      
        // because it's not possible the share the OpenGl Display
        // it crashs the aplication                              
        computeCooMax();
        emit loadedData(part_data,psv);
        is_new_data_loaded = FALSE;
      } 
      else {
        is_new_data_loaded = TRUE;
      }
      return 1;
    } // end of try
    catch (int n) {
       switch(n) {
       case -1 : cerr << "Catch error during SEND\n";
        break;
       case -2 : cerr << "Catch error during RECV\n";
        break;
       case -3 : cerr << "Catch error during RECV, nbuffer=0!\n";
        break;
       default :
	assert(1);     
       }
       is_end_of_data = TRUE;
       is_connected = FALSE;
       return 0;
    } // end of catch
  } // end if clientMB
  else {
    return 0;
  }
}
// ============================================================================
// NetworkData::uploadGlData()                                                 
// send data to the opengl object                                              
void NetworkData::uploadGlData(ParticlesSelectVector * psv)
{
   if (is_new_data_loaded) {
    emit loadedData(part_data,psv);
    is_new_data_loaded = FALSE;
  } 
}
// ============================================================================
// NetworkData::getNbody()                                                     
// return #bodies from the server                                              
int NetworkData::getNbody()
{
  PRINT_D cerr << "get Nbody Adress clienMB=" << clientMB << "\n";
  if ( (*part_data->nbody)<0 && is_connected) { // Nbody not got
  
    // Wait hello Data
    clientMB->recvData();   
    if (*(clientMB->getTagBuffer()) !=  MessageBuffer::Hello) {    
      throw (-4); // Bad protocol error
    }
    data_name = clientMB->getDatBuffer();
    PRINT_D cerr << "Got hello from server : [" << clientMB->getDatBuffer() << "]\n";
    
    // Wait Nbody data
    clientMB->recvData();   
    if (*(clientMB->getTagBuffer()) !=  MessageBuffer::Nbody) {
      throw (-4); // Bad protocol error
    }
    *part_data->nbody = *((int *) (clientMB->getDatBuffer()));
    full_nbody = *part_data->nbody;
    PRINT_D cerr << "Got Nbody from server : [" << *part_data->nbody << "]\n";
  } 
  return *part_data->nbody;
}
// ============================================================================
// NetworkData::getTime()                                                      
// return time                                                                 
float NetworkData::getTime()
{
  return *part_data->timu;
}
// ============================================================================
// NetworkData::endOfDataMessage()                                             
// return message in case of lost communication                                
QString NetworkData::endOfDataMessage()
{
  QString message="Lost connexion with server";
  return message;
}
// ============================================================================
