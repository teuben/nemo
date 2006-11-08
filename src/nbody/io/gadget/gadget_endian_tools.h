// -----------------------------------------------------------------------------
// gadget_endian_tools.h                                                        
// Definitition of the class GadgetEndianTools                                  
//                                                                              
// This class perform data IO operation according to the endianness             
//                                                                              
// Intel/Amd/Alpha are little-endian compliant                                  
//                                                                              
// PowerPC IBM/Sparc/SGI are big-endian compliant                               
// -----------------------------------------------------------------------------
#ifndef GADGET_ENDIAN_TOOLS_H
#define GADGET_ENDIAN_TOOLS_H

#include <stdio.h>

#include "gadget_data_structure.h"

class GadgetEndianTools {
 public:
  enum ioop { READ,WRITE };

  //constructor
  GadgetEndianTools(const FILE * _fd,const  bool _swap);
  //destructor
  ~GadgetEndianTools();

  int ioHeader(t_io_header_1 * header,const ioop op);
  int ioData(char * ptr,int size_bytes, int items,const ioop op);

 private:
  bool swap;                // control swapping variable
  const FILE * fd;          // file descriptor

  void swapBytes(void * x,const int size);

};
#endif
// -----------------------------------------------------------------------------
