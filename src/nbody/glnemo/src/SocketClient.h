// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2007                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique de galaxies                                             
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================

#ifndef MY_SOCKET_CLIENT_H
#define MY_SOCKET_CLIENT_H

#include "Socket.h"

class SocketClient : public Socket {
public:
  SocketClient(const char * , int, int); // server's name, sock's type, port
  ~SocketClient();
  int connectServer();
private:
  struct sockaddr_in serAddr;   // server's specification
};
#endif

// ============================================================================
