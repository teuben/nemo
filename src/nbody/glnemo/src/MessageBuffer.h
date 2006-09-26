// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2005                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique de galaxies                                             
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
//                                                                             
// Definition of the class MessageBuffer                                       
//                                                                             
// This class allow to send and receive data of any size over                  
// the network and through socket mechanism                                    
// ============================================================================
#ifndef MESSAGE_BUFFER_H
#define MESSAGE_BUFFER_H

#include "ChunkBuffer.h"

#define FLOAT_BYTE 4
#define INT_BYTE   4
#define CHAR_BYTE  1

class MessageBuffer  {
  
public :

  enum TAG {
    Hello    = -0x0001,
    Nbody    = -0x0002,
    Time     = -0x0003,
    Pos      = -0x0004,
    Vel      = -0x0014,
    Select   = -0x0005,
    SelectV  = -0x0015,
    Stop     = -0x0006,
    ACK      = -0x0007,
    EndOfMes = -0x0666
  };

  enum TAG_TYPE {
    HelloB   = CHAR_BYTE  ,
    NbodyB   = INT_BYTE   ,
    TimeB    = FLOAT_BYTE ,
    DataB    = FLOAT_BYTE ,
    SelectB  = CHAR_BYTE  ,
    StopB    = INT_BYTE
  };
 
  MessageBuffer(int sd, int length=8192);
  ~MessageBuffer();
  int recvData(); 
  char * newBuffer(int length);

  int   * getTagBuffer();
  int   * getSizBuffer();
  int     getSizDatBuffer(const TAG tag);
  char  * getDatBuffer();

  int sendAll(int, char *, int *);  // sd, buffer, len
  int sendAll(char *, int *);       // buffer, len    
  int sendData(const TAG ,int ,const char * );
protected :
  char * buffer;                   // buffer to store received data 

  ChunkBuffer * tag_chunk;         // Store tag buffer 
  ChunkBuffer * siz_chunk;         // Store size buffer
  ChunkBuffer * dat_chunk;         // Store Data buffer

private :
  //variables related to buffer
  int    buffer_length;  // max buffer size                          
  char * p_buffer;       // pointer to the current byte in the buffer
  int    n_buffer;       // #bytes received in the buffer            
  int    send_flags;     // flags for send(2)

  bool   is_empty;       // true if nothing in the buffer            
  
  // variables related to socket
  int sock_fd;          // socket file descriptor                    

  // Method
  bool isEndOfMessage();
  int sizeOfData(const TAG);
};
#endif // MESSAGE_BUFFER_H
// ============================================================================

