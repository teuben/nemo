// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// ChunkBuffer.h                                                               |
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
#ifndef CHUNK_BUFFER_H
#define CHUNK_BUFFER_H


class ChunkBuffer {
public:
  ChunkBuffer(int length_buffer);
  ~ChunkBuffer();

  bool   complete;           // true if buffer is complete

  char * parseBuffer(char *, int &, bool &);
  int getIntValue();         // return the int value of the buffer's 4 first bytes
  char * getBuffer();        // return the content of the chunk buffer            
  int razVar();              // raz variables                                     
  int getLength() { return l_chunk; };

protected:
  char * chunk;              // buffer to store data
  int    l_chunk;            // buffer length       
  int    n_chunk;            // #bytes in the buffer

};

#endif  // CHUNK_BUFFER_H
//------------------------------------------------------------------------------
