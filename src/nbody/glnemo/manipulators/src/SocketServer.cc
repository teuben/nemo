// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// SocketServer.cc                                                             |
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
// This class provide all the needed method to run a Server                     
// using message socket protocol                                                
//------------------------------------------------------------------------------
#include <iostream>
#include <errno.h>

#include "SocketServer.h"


using namespace std; // prevent writing statment like 'std::cerr'

//------------------------------------------------------------------------------
// Constructor                                                                  
// - Create a socket                                                            
// - Bind socket to the listen port                                             
// - Listen socket                                                              
SocketServer::SocketServer(int type, 
			   int listen_port,
			   int backlog):Socket(type)
{
  port    = listen_port;    
  sockOpt();             // specify sockets options
  bindSocket();          // bind socket to the listen port
  listenSocket(backlog); // listen socket
}
//------------------------------------------------------------------------------
// bindSocket:                                                                  
// Bind the socket to the listening port                                        
int SocketServer::bindSocket()
{
  struct sockaddr_in my_addr;
  int status;

  // set host Information 
  my_addr.sin_family = AF_INET;                 // socket family
  my_addr.sin_addr.s_addr = htonl(INADDR_ANY);
  my_addr.sin_port = htons(port);
  memset(&(my_addr.sin_zero),'\0',8);

  status = bind (sd, (struct sockaddr *) &my_addr,sizeof(my_addr));
  if(status<0) {
    std::cerr << "SocketServer::bindSocket() - WARNING cannot bind port number ["
	      << port << "]\n";
    perror("SocketServer::bindSocket() - PERROR bind:");
  }

  return status;
}
//------------------------------------------------------------------------------
// listenSocket:                                                                
// listen socket                                                                
int SocketServer::listenSocket(int backlog)
{
  int status;
  status = listen (sd, backlog);
  if(status<0) {
    perror("SocketServer::listenSocket() - PERROR listen:");
  }
  return status;
}
//------------------------------------------------------------------------------
// acceptClient:                                                                
// establish a connection with a new Client                                     
// return socket descritor                                                      
int SocketServer::acceptClient()
{
  int cliLen;
  cliLen = sizeof(cliAddr);

  sock_accept_client = accept(sd, (struct sockaddr *) &cliAddr, 
			      //#ifdef linux
			      (socklen_t *) &cliLen);
  //#else
  //  &cliLen);
  //#endif

  if(sock_accept_client<0) {
    perror("SocketServer::acceptClient() - PERROR accept:");
  }
  return sock_accept_client;
}
//------------------------------------------------------------------------------

