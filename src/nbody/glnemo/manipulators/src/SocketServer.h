// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// SocketServer.h                                                              |
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
#ifndef MY_SOCKET_SERVER_H
#define MY_SOCKET_SERVER_H

#include "Socket.h"

class SocketServer : public Socket {
public:
  SocketServer(int , int, int);
  ~SocketServer();
  int acceptClient();

private:
  int port;               // listening port                 
  int sock_accept_client; // socket id returned after accept
  struct sockaddr_in cliAddr;

  int bindSocket();
  int listenSocket(int);
};
#endif
//------------------------------------------------------------------------------
