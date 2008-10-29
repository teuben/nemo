// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// ChunkBuffer.cc                                                              |
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
#include <iostream>
#include <assert.h>
#include <cstring>
#include "ChunkBuffer.h"
using namespace std;
//------------------------------------------------------------------------------
// constructor                                                                  
// Allocate a chunk buffer with size_buffer bytes                               
ChunkBuffer::ChunkBuffer(int size_buffer)
{
  //
  razVar();
  // allocate memory for buffer
  l_chunk = size_buffer;
  chunk = new char[l_chunk];
}
// -----------------------------------------------------------------------------
// destructor                                                                   
// deallocate chunk buffer.                                                     
ChunkBuffer::~ChunkBuffer()
{
  delete[] chunk;
}
//------------------------------------------------------------------------------
// char * ChunkBuffer::parseBuffer                                              
// parse the received buffer and return a pointer on the position of            
// the reveived buffer not yet parsed                                           
// fill up chunk buffer according to #bytes received                            
char *  ChunkBuffer::parseBuffer(char * recv_buffer, int & nrecv, bool & is_empty)
{
  int pop_bytes;

  if (nrecv>0) {
    if (nrecv >= (l_chunk-n_chunk) ) {    // full tag buffer               
      if ((n_chunk > l_chunk)) {
	//std::cerr << n_chunk << " " << l_chunk << "\n";
	assert((n_chunk > l_chunk));
      }
      pop_bytes=(l_chunk-n_chunk);        // whole length of tag buffer    
      assert(pop_bytes >= 0);
      complete = true;                    // current chunk buffer completed
    } 
    else {
      pop_bytes=nrecv;                    // got only nrecv bytes          
    }

    memcpy(chunk+n_chunk,recv_buffer,pop_bytes); // fill up chunk                     
    recv_buffer += pop_bytes;                    // pop 'pop_bytes' bytes from buffer,
                                                 // incr buffer pointer               

    nrecv -= pop_bytes;  // #bytes remaining in the buffer  
    assert(nrecv >= 0);  // nrecv < 0 : ERROR               
    if (nrecv == 0 ) {   // no more bytes in reveived buffer
      is_empty = true;   // buffer is empty                 
    }

    n_chunk  += pop_bytes;      // size of chunk buffer increase 
    assert(n_chunk <= l_chunk); // n_chunk > l_chunk : ERROR     
    if (n_chunk == l_chunk) {   // received == chunk's size      
      complete = true;          // current chunk buffer completed
    } 
  }
  else {  // should not arrive here, means that nrecv=0, impossible !!
    assert(1);
  }
  return recv_buffer;
}
//------------------------------------------------------------------------------
// char * ChunkBuffer::getIntValue()                                            
// return int value stored in a buffer                                          
int  ChunkBuffer::getIntValue()
{
  int * p = (int *) chunk;
  return (*p);
}
//------------------------------------------------------------------------------
// char * ChunkBuffer::getBuffer()                                              
char * ChunkBuffer::getBuffer()
{
  return chunk;
}
//------------------------------------------------------------------------------
// ChunkBuffer::razVar()                                                        
// reset buffer                                                                 
int  ChunkBuffer::razVar()
{
  complete = false;
  n_chunk  = 0;

  return 1;
}
//-----------------------------------------------------------------------------
