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
#include <cstdio>
#include <netinet/in.h>
#include "SocketClient.h"


using namespace std; // prevent writing statment like 'std::cerr'

// ============================================================================
// Constructor                                                                 
// Create a socket                                                             
SocketClient::SocketClient(const char * server, 
		 	   int type,
			   int listen_port):Socket(type)
{
  getSocketAddr(server,&serAddr,listen_port);
}
// ============================================================================
// Destructor                                                                  
// close socket                                                                
SocketClient::~SocketClient()
{
  Socket::close();
}
// ============================================================================
// SocketClient::connectServer()                                               
// Connect socket to the server                                                
int SocketClient::connectServer()
{
  if (connect(sd,(struct sockaddr *) &serAddr, 
	      sizeof(serAddr)) == -1 ) {
    perror("connect");
    sd = -1;
  }
  return sd;
}
// ============================================================================
