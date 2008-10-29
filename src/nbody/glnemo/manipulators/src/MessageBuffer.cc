// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// MessageBuffer.cc                                                            |
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
//                                                                             |
// Implementation of the class MessageBuffer                                   |
//                                                                             |
// This class allow to send and receive data of any size over                  |
// the network and through socket mechanism                                    |
//-----------------------------------------------------------------------------+
#include <iostream>
#include <sys/socket.h>
#include <netinet/in.h>
#include <errno.h>
#include <unistd.h>
#include <assert.h>
#include <cstdlib>
#include "MessageBuffer.h"

using namespace std;

//------------------------------------------------------------------------------

//#define BUFFER_SIZE 8192*2
#define MIN(A,B) ((A)<(B)?(A):(B))

#define LOCAL_DEBUG 0
#include "print_debug.h"

//------------------------------------------------------------------------------
// Constructor                                                                  
// store socket file descriptors and allocate a working buffer of               
// size length                                                                  
//------------------------------------------------------------------------------
MessageBuffer::MessageBuffer(int sd, int length)
{
  //
  buffer = NULL;
  // Initialyze buffer
  newBuffer(length);

  // init value
  sock_fd   = sd;

  // init Chunk
  tag_chunk = new ChunkBuffer(4);
  assert(tag_chunk != NULL);
  siz_chunk = new ChunkBuffer(4);
  assert(siz_chunk != NULL);
  dat_chunk = NULL;
}
//------------------------------------------------------------------------------
// destructor                                                                   
// deallocate buffers                                                           
//------------------------------------------------------------------------------
MessageBuffer::~MessageBuffer()
{
  PRINT_D std::cerr << "MessageBuffer::~MessageBuffer(): Kill p_buffer\n";
  if (p_buffer) {
    delete[] p_buffer;
  }
  PRINT_D std::cerr << "MessageBuffer::~MessageBuffer(): Kill Tag\n";
  if (tag_chunk) {
    delete tag_chunk;
  }
  PRINT_D std::cerr << "MessageBuffer::~MessageBuffer(): Kill Siz\n";
  if (siz_chunk) {
    delete  siz_chunk;
  }
  PRINT_D std::cerr << "MessageBuffer::~MessageBuffer(): Kill Dat\n";
  if (dat_chunk) {
    delete dat_chunk;
  }
}
//------------------------------------------------------------------------------
// MessageBuffer::recvData                                                      
// allocate the receive buffer with a new size                                  
//------------------------------------------------------------------------------
char * MessageBuffer::newBuffer(int length)
{
  // allocate memory for the buffer
  buffer_length = length;
  buffer = new char[buffer_length];
  assert(buffer != NULL);
  if (!buffer) {
    std::cerr << "MessageBuffer::newBuffer() - ERROR allocation failed for the buffer\n";
    //std::exit(1);
  }
  p_buffer = buffer;

  // init value
  is_empty  = true;

  return buffer;
}

//------------------------------------------------------------------------------
// MessageBuffer::sendAll                                                       
// Send data stored in buffer, through socket sd                                
int  MessageBuffer::sendAll( char * s_buffer, int * len)
{
  return sendAll(sock_fd,s_buffer,len);
}
//------------------------------------------------------------------------------
// MessageBuffer::sendAll                                                       
// Send data stored in buffer, through socket sd                                
int MessageBuffer::sendAll(int sd, char * s_buffer, int * len)
{
  int total=0;      // how many bytes we have sent
  int b_left= *len; // how many bytes we have left to send
  int n_send;
  int n;

  while (total < *len) {
    //n_send = MIN(b_left,BUFFER_SIZE);
    n_send = MIN(b_left,buffer_length);
    //fprintf(stderr,"N-send expected ==> [%d]\n",n_send);
    n = send(sd, s_buffer+total, n_send,MSG_NOSIGNAL);
    if ( n == -1) {
      perror("MessageBuffer::sendAll - PERROR send:");
      throw(-1); // throw execption 
      return -1;
      break; // error during send
    }
    //fprintf(stderr,"N-send real ==> [%d]\n",n);
    total += n;

    b_left -= n;
  }
  //cerr << "End of sendall [" << total << "]\n";
  return (n==-1?-1:0); // -1 on failure, 0 on success 
  
}
//------------------------------------------------------------------------------
// MessageBuffer::sendData                                                      
// send a message structured as following :                                     
// TAG              : 4 bytes                                                   
// Size of buffer   : 4 bytes                                                   
// buffer           : array of bytes                                            
int MessageBuffer::sendData(const TAG tag,int n_element,const char * s_buffer )
{
  int status;

  // Send [tag] 
  int len = sizeof(int);                          // #bytes
  status = sendAll((char *) &tag, &len);

  if (status >= 0 ) {
    // Send [n_element]
    int len = sizeof(int);
    int lendata = sizeOfData(tag) * n_element;    // #bytes      
    status = sendAll((char *) &lendata, &len);    // !!! important, we send the #bytes of
                                                  // the buffer size !!!
    
    if (status >=0 ) {
      // Send [s_buffer]
      //cerr << "Send s_buffer in INT = " << * ( (int*) s_buffer)<< "\n";
      int lendata = sizeOfData(tag) * n_element;  // #bytes      
      status = sendAll((char *) s_buffer, &lendata);

#if 0
      if (status >= 0 ) {
	// Send end of message
	int lentag=sizeof(int);
	int tag=EndOfMes;
	status = sendAll((char *) &tag,&lentag );
      }
#endif
    }
  }
  return status;
}
//------------------------------------------------------------------------------
// MessageBuffer::recvData                                                      
// Listen to the socket                                                         
// fill up the receiv buffer                                                    
// parse the buffer :                                                           
//       TAG --> SIZEofBUFFER --> BUFFER                                        
//------------------------------------------------------------------------------
int MessageBuffer::recvData()
{
  bool end_of_message=false;

  if (dat_chunk != NULL) {                                    // exist but completed
    delete dat_chunk;                                         // delete
    dat_chunk = NULL;
  }

  while (!end_of_message) {
    //cerr << "Recv : is_empty = " << is_empty << "\n";
    if (is_empty) {                                       // empty buffer         
      buffer = p_buffer;                                  // set to the beginning 
      if ((n_buffer=                                      // store #bytes received
	   recv(sock_fd,buffer,buffer_length,MSG_NOSIGNAL))==-1) {   // recv data 
	perror("MessageBuffer::recvData - PERROR recv:");            // failed               
	close(sock_fd);                                   // close sd             
	throw(-2);
	return -1;
      }
      if (n_buffer) {   // we have read something              
	//std::cerr << "recv: n_buffer [" << n_buffer << "]\n";
	is_empty=false; // buffer is no more empty             
      } else {
	throw(-3);
      }
    } // if .. is_empty

    
    //parseBuffer(); // Find out TAG
    if ( ! tag_chunk->complete) {
      PRINT_D std::cerr << "< Processing TAG....\n";
      buffer = tag_chunk->parseBuffer(buffer,n_buffer,is_empty);
      PRINT_D std::cerr << "> Processing TAG...."<< n_buffer <<"]\n";
    } 
    else {
      //std::cerr << "Tag = " << tag_chunk->getIntValue() << "\n";
      if (  ! siz_chunk->complete) {
	PRINT_D std::cerr << "<< Processing SIZ....\n";
	buffer = siz_chunk->parseBuffer(buffer,n_buffer,is_empty);
	PRINT_D std::cerr << ">> Processing SIZ....["<< n_buffer <<"]\n";
      } 
      else {
	//std::cerr << "Siz = " << siz_chunk->getIntValue() << "\n";
	if (siz_chunk->getIntValue() > 0) { 
	  //std::cerr << "siz_chunk->getIntValue() =[ " << siz_chunk->getIntValue() << "]\n";

	  if (dat_chunk == NULL ||                                      // does not exist 
	      (dat_chunk->complete )) {                                 // already completed

	    if (dat_chunk != NULL) {                                    // exist but copleted
	      delete dat_chunk;                                         // delete
	    }
	    //int size_dat = siz_chunk->getIntValue() * sizeof(float) * 3;

	    //int size_dat = getSizDatBuffer((MessageBuffer::TAG) (* getTagBuffer()));
	    int size_dat = (* getSizBuffer());
	    PRINT_D cerr << "TAG = " << (* getTagBuffer()) << " SIZE = " << (* getSizBuffer()) << "\n";
	    PRINT_D cerr << "Size Dat Buffer = ["<<size_dat<<"]";
	    dat_chunk = new ChunkBuffer( size_dat );
	    assert(dat_chunk);
	    PRINT_D std::cerr << "Processing DAT....["<< n_buffer <<"]\n";
	  } 
	  if ( ! dat_chunk->complete ) {
	    buffer = dat_chunk->parseBuffer(buffer,n_buffer,is_empty);
	    //std::cerr << "after parsed Data Buffer ["<< n_buffer <<"] is_empty ?" << is_empty << "\n";
	  } else {
	    PRINT_D std::cerr << "Data COMPLETE !!! \n";
	    tag_chunk->razVar();
	    siz_chunk->razVar();
	  }
	}
      }
    }
    end_of_message = isEndOfMessage();
  }
  return 1;
}

//---------------------------------------------------------------------------------
// MessageBuffer::isEndOfMessage                                                   
// test if it's the end of message                                                 
bool MessageBuffer::isEndOfMessage()
{
  //if (( tag_chunk->complete)  && tag_chunk->getIntValue()==EndOfMes) {
  if ((dat_chunk != NULL && dat_chunk->complete) ) {
    PRINT_D std::cerr << "isEndOfMessage=TRUE, Tag_chunk buffer = " 
		      << tag_chunk->getIntValue() << "\n";
    //cerr << "TAG buffer equal in INT " << *((int *) (getDatBuffer())) << "\n";
    tag_chunk->razVar();
    siz_chunk->razVar();
    return true;
  }
  else
    return false;
}
//---------------------------------------------------------------------------------
// MessageBuffer::getTagBuffer                                                     
// return the TAG value                                                            
int * MessageBuffer::getTagBuffer(){
  return ((int *) tag_chunk->getBuffer());
}
//---------------------------------------------------------------------------------
// MessageBuffer::getSizBuffer                                                     
// return the #bytes in the data buffer                                            
int * MessageBuffer::getSizBuffer(){
  return ((int *) siz_chunk->getBuffer());
}
//---------------------------------------------------------------------------------
// MessageBuffer::getSizDatBuffer                                                  
// return the data size store in the buffer                                        
int   MessageBuffer::getSizDatBuffer(const TAG tag){
  return (  *((int *) siz_chunk->getBuffer()) * sizeOfData(tag)) ;
}

//---------------------------------------------------------------------------------
// MessageBuffer::getDatBuffer()                                                   
// return the Data buffer (array of bytes)
char * MessageBuffer::getDatBuffer(){
  //return ((float *) dat_chunk->getBuffer());
  return ((char *) dat_chunk->getBuffer());
}

//---------------------------------------------------------------------------------
// MessageBuffer::sizeOfData                                                       
// return the Data size according to the TAG                                       
int MessageBuffer::sizeOfData(const TAG tag)
{
  int size=0;
  switch (tag) {
  case Hello    :
  case Select   :
  case SelectV  :
    size = sizeof(char);
    break;
  case Nbody    :
  case EndOfMes :
  case Stop     :
    size = sizeof(int);
    break;
  case Time     :
  case Pos      :
  case Vel      :
    size = sizeof(float);
    break;
  default       :
    cerr << "Unknown TAG, exiting....\n";
    exit(1);
    break;
  }
  return size;
}
//---------------------------------------------------------------------------------
