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
#include <iostream>
#include <errno.h>
#include <netdb.h>
#include "Socket.h"

#define LOCAL_DEBUG 0
#include "print_debug.h"

using namespace std; // prevent writing statment like 'std::cerr'

// ============================================================================
// Constructor                                                                 
// Create a socket                                                             
Socket::Socket(int type)
{
  PRINT_D cerr << "I am in Socket constructor...\n";
  if ( type != SOCK_DGRAM && type != SOCK_STREAM) {
    std::cerr << "Unknown socket type  aborted..\n";
    exit(1);
  }
  /* get a socket descriptor */
  sd = socket(AF_INET, type,0);
  if ( sd == -1 ) {
    perror("Socket function failed in [get_socket]");
    //exit(1);
  }
}
// ============================================================================
// Socket::sockOpt()                                                           
// set Socket Options                                                          
int Socket::sockOpt()
{
  int yes=1, status = 1;
  if (setsockopt(sd,SOL_SOCKET, SO_REUSEADDR,&yes,
		 sizeof(int))== -1 )  {
    perror("sock_opt");
    status = -1;
    //exit(1);
  }
  return status;
}
// ============================================================================
// Socket::getSocketAddr()                                                     
// get Socket Adress desctiption according to the host                         
int Socket::getSocketAddr(const char * host_name,             // hostname      
			    struct sockaddr_in *remoteAddr ,  // remote Addr   
			    int port)                         // listening port
{
  struct hostent * h;
  int status=1;

  if ((h =gethostbyname(host_name)) == NULL) {
    perror("gethostbyname");
    status = -1;
  }
  else {
    /* set host Information */
    remoteAddr->sin_family = AF_INET;
    remoteAddr->sin_addr = *((struct in_addr *) h->h_addr);
    remoteAddr->sin_port = htons(port);
    memset(&(remoteAddr->sin_zero),'\0',8);
  }
  return status;
}
// ============================================================================
