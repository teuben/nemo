/* -*-c++-*- */
// -------------------------------------------------------------
// -------------------------------------------------------------

#ifndef MY_SERVER_CLIENT_DEF_H
#define MY_SERVER_CLIENT_DEF_H

#define SERV_PORT 5431
#define BACKLOG 25

extern "C" {
  int io_nemo(char * , char *, ...);
}
#endif
// -------------------------------------------------------------
