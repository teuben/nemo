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

#ifndef MY_SOCKET_H
#define MY_SOCKET_H

// Sockets's Headers
#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>

class Socket {
public:
  Socket(int);
  //~Socket();
      
  int sockOpt();             // specify socket option        
  int sockFd() { return sd;} // return socket file descriptor


protected:
  int sd;                     // socket file descriptor      
  int getSocketAddr(const char *,           // hostname      
		    struct sockaddr_in *,   // remote Addr   
		    int);                   // listening port

  void close() { ::close(sd); };
};
#endif
// ============================================================================
