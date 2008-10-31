// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/public/lib/io.cc                                                
///                                                                             
/// \brief  contains definitions of methods declared in inc/public/io.h         
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2000-2008                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2000-2008 Walter Dehnen                                        
//                                                                              
// This program is free software; you can redistribute it and/or modify         
// it under the terms of the GNU General Public License as published by         
// the Free Software Foundation; either version 2 of the License, or (at        
// your option) any later version.                                              
//                                                                              
// This program is distributed in the hope that it will be useful, but          
// WITHOUT ANY WARRANTY; without even the implied warranty of                   
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
// General Public License for more details.                                     
//                                                                              
// You should have received a copy of the GNU General Public License            
// along with this program; if not, write to the Free Software                  
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                    
//                                                                              
////////////////////////////////////////////////////////////////////////////////
#include <public/io.h>
#include <fstream>
#include <cstdlib>
#include <cstring>

#include <public/nemo++.h>

using namespace falcON;
////////////////////////////////////////////////////////////////////////////////
namespace {
  //----------------------------------------------------------------------------
  int openstdout = 0;
  int openstdin  = 0;
}
void output::open_std() falcON_THROWING {
  if( ++openstdout > 1 )
    falcON_THROW("trying to open more than one output to stdout");
}
void output::close_std() {
  if(openstdout) --openstdout;
}
void input::open_std() falcON_THROWING {
  if( ++openstdin > 1 )
    falcON_THROW("trying to open more than one input from stdin");
}
void input::close_std() {
  if(openstdin) --openstdin;
}
// /////////////////////////////////////////////////////////////////////////////
//                                                                              
// falcON::FileSize()                                                           
//                                                                              
// /////////////////////////////////////////////////////////////////////////////
size_t falcON::FileSize(const char*sFileName)
{
  std::ifstream f;
  f.open(sFileName, std::ios_base::binary | std::ios_base::in);
  if (!f.good() || f.eof() || !f.is_open()) { return 0; }
  f.seekg(0, std::ios_base::beg);
  std::ifstream::pos_type begin_pos = f.tellg();
  f.seekg(0, std::ios_base::end);
  return f.tellg() - begin_pos;
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::output                                                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
void output::__open(bool append)
{
  APPENDING = false;
  if     (0 == FILE    ||
	  0 == FILE[0] ||
	  0 == std::strcmp(FILE,".") ) {
    OUT = 0;
    DebugInfo(2,"output: open sink\n");
  } else if(0 == std::strcmp(FILE,"-") ) {
    open_std();
    OUT = &std::cout;
    DebugInfo(2,"output: open stdout\n");
  } else {
    std::ofstream *FOUT = new std::ofstream();
    if(append) {
      FOUT->open(FILE,std::ios::out | std::ios::app);
      if(FOUT->is_open()) {
	APPENDING = true;
	DebugInfo(2,"output: append to file \"%s\"\n",FILE);
      }
    }
    if(!FOUT->is_open() )
      FOUT->open(FILE,std::ios::out);
    if( FOUT->is_open() ) {
      OUT = FOUT;
      DebugInfo(2,"output: open file \"%s\"\n",FILE);
    } else {
      DebugInfo(2,"output: could not open file \"%s\"\n",FILE);
      OUT = 0;
      falcON_DEL_O(FOUT);
    }
  }
}
//------------------------------------------------------------------------------
void output::close() {
  if(FREC) {
    if(FILE)
      falcON_Warning("closing FortranORec before output from file \"%s\"\n",
		     FILE);
    else
      falcON_Warning("closing FortranORec before output\n");
    FREC->close();
  }
  if(OUT) {
    DebugInfo(2,"output: closing\n");
    if(OUT == &std::cout) close_std();
    else falcON_DEL_O(OUT);
  }
  APPENDING = false;
  OUT = 0;
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::input                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
void input::__open() {
  if     (0 == FILE || FILE[0] == 0) {
    IN = 0;
    DebugInfo(2,"input: empty file\n");
  } else if(0 == std::strcmp(FILE,"-") ) {
    open_std();
    IN= &std::cin;
    DebugInfo(2,"input: stdin\n");
  } else {
    std::ifstream *FIN = new std::ifstream(FILE);
    if( FIN->is_open() ) {
      IN = FIN;
      DebugInfo(2,"input: open file \"%s\"\n",FILE);
    } else {
      DebugInfo(2,"input: could not open file \"%s\"\n",FILE);
      IN = 0;
      falcON_DEL_O(FIN);
    }
  }
}
//------------------------------------------------------------------------------
void input::close() {
  if(FREC) {
    if(FILE)
      falcON_Warning("closing FortranIRec before input from file \"%s\"\n",
		     FILE);
    else
      falcON_Warning("closing FortranIRec before input\n");
    FREC->close();
  }
  DebugInfo(2,"input: closing\n");
  if(IN == &std::cin) close_std();
  else if(IN) falcON_DEL_O(IN);
  IN = 0;
}
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class falcON::FortranIRec                                                    
//                                                                              
////////////////////////////////////////////////////////////////////////////////
unsigned falcON::FortranIRec::read_size() throw(falcON::exception)
{
  if(HSZE == 4) {
    uint32 S;
    IN.read(static_cast<char*>(static_cast<void*>(&S)),sizeof(uint32));
    if(SWAP) swap_bytes(S);
    return S;
  } else if(HSZE == 8) {
    uint64 S;
    IN.read(static_cast<char*>(static_cast<void*>(&S)),sizeof(uint64));
    if(SWAP) swap_bytes(S);
    return S;
  } else 
    throw exception("FortranIRec: header size must be 4 or 8\n");
  return 0;
}
//------------------------------------------------------------------------------
falcON::FortranIRec::FortranIRec(input& in, unsigned rec, bool swap)
  throw(falcON::exception) : IN(in), HSZE(rec), SWAP(swap), READ(0)
{
  DebugInfo(8,"FortranIRec: opening ... \n");
  if(!IN) throw exception("FortranIRec::FortranIRec(): input corrupted");
  if(IN.FREC)
    throw exception("trying to open 2nd FortranIRec to same input\n");
  IN.FREC = this;
  SIZE = read_size();
  DebugInfo(6,"FortranIRec: opened with %u bytes\n",SIZE);
}
//------------------------------------------------------------------------------
unsigned falcON::FortranIRec::read_bytes(char*buf, unsigned n)
  throw(falcON::exception)
{
  if(!IN) throw exception("FortranIRec::read_bytes(): input corrupted");
  if(READ + n > SIZE) {
    falcON_Warning("FortranIRec::read(): cannot read %u, but only %u bytes\n",
		   n, SIZE-READ);
    n = SIZE - READ;
  }
  IN.read(buf,n);
  if(!IN) throw exception("FortranIRec: input corrupted");
  READ += n;
  DebugInfo(6,"FortranIRec: read %u bytes\n",n);
  return n;
}
//------------------------------------------------------------------------------
void falcON::FortranIRec::skip_bytes(unsigned n) {
  if(READ + n > SIZE) n = SIZE - READ;
  if(n && !IN) throw exception("FortranIRec::skip_bytes(): input corrupted");
  for(char C; n; --n,++READ) IN.read(&C,1);
}
//------------------------------------------------------------------------------
void falcON::FortranIRec::close() throw(falcON::exception)
{
  if(!IN) throw exception("FortranIRec::close(): input corrupted");
  if(READ != SIZE) {
    falcON_Warning("FortranIRec: only %u of %u bytes read on closing record\n",
		   READ, SIZE);
    for(char C; READ!=SIZE; ++READ) IN.read(&C,1);
  }
  unsigned S = read_size();
  IN.FREC = 0;
  if(S != SIZE) throw exception("FortranIRec: record size mismatch");
  DebugInfo(6,"FortranIRec: closed with %u bytes\n",SIZE);
}
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class falcON::FortranORec                                                    
//                                                                              
////////////////////////////////////////////////////////////////////////////////
void falcON::FortranORec::write_size() throw(falcON::exception)
{
  if(HSZE == 4) {
    uint32 S = SIZE;
    OUT.write(static_cast<char*>(static_cast<void*>(&S)),sizeof(uint32));
  } else if(HSZE == 8) {
    uint64 S;
    OUT.write(static_cast<char*>(static_cast<void*>(&S)),sizeof(uint64));
  } else 
    throw exception("FortranORec: header size must be 4 or 8\n");
}
//------------------------------------------------------------------------------
falcON::FortranORec::FortranORec(output& out, unsigned size, unsigned rec)
  throw(falcON::exception) : OUT(out), HSZE(rec), SIZE(size), WRITTEN(0)
{
  if(!OUT) throw exception("FortranORec: output corrupted");
  if(OUT.FREC)
    throw exception("trying to open 2nd FortranORec to same output\n");
  OUT.FREC = this;
  write_size();
  DebugInfo(6,"FortranORec: opened for %u bytes\n",SIZE);
}
//------------------------------------------------------------------------------
unsigned falcON::FortranORec::write_bytes(const char*buf, unsigned n)
  throw(falcON::exception)
{
  if(!OUT) throw exception("FortranORec: output corrupted");
  if(WRITTEN + n > SIZE) {
    falcON_Warning("FortranORec::write(): cannot read %u, but only %u bytes\n",
		   n, SIZE-WRITTEN);
    n = SIZE - WRITTEN;
  }
  OUT.write(buf,n);
  if(!OUT) throw exception("FortranORec: ostream corrupted");
  WRITTEN += n;
  DebugInfo(6,"FortranORec: written %u bytes\n",n);
  return n;
}
//------------------------------------------------------------------------------
void falcON::FortranORec::fill_bytes(unsigned n, char C) {
  if(WRITTEN + n > SIZE) n = SIZE - WRITTEN;
  for(; n; --n,++WRITTEN) OUT.write(&C,1);
}
//------------------------------------------------------------------------------
void falcON::FortranORec::close() throw(falcON::exception)
{
  if(!OUT) throw exception("FortranORec: output corrupted");
  if(WRITTEN!=SIZE) {
    falcON_Warning("FortranORec: only %u of %u bytes written on closing record"
		   " ... padding with 0\n", WRITTEN, SIZE);
    for(char C=0; WRITTEN!=SIZE; ++WRITTEN) OUT.write(&C,1);
  }
  write_size();
  OUT.FREC = 0;
  DebugInfo(6,"FortranORec: closed with %u bytes\n",SIZE);
}
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// struct falcON::GadgetHeader                                                  
//                                                                              
////////////////////////////////////////////////////////////////////////////////
falcON::GadgetHeader::GadgetHeader() :
  time(0.), redshift(0.), flag_sfr(0), flag_feedback(0), flag_cooling(0),
  num_files(0), BoxSize(0.), Omega0(0.), OmegaLambda(0.), HubbleParam(0.),
  flag_stellarage(0), flag_metals(0), flag_entropy_instead_u(0)
{
  for(int k=0; k!=6; ++k) {
    npart[k] = 0;
    npartTotal[k] = 0;
    npartTotalHighWord[k] = 0;
    masstab[k] = 0.;
  }
}
//------------------------------------------------------------------------------
bool falcON::GadgetHeader::Read(input& in, unsigned rec, bool& swap)
  throw(falcON::exception)
{
  swap = 0;
  // read record header and determine swapping necessity
  if(rec == 4) {
    uint32 S;
    in.read(static_cast<char*>(static_cast<void*>(&S)), sizeof(uint32));
    if(S != sizeof(GadgetHeader)) {
      swap_bytes(S);
      if(S == sizeof(GadgetHeader)) swap = 1;
      else return false;
    }
  } else if(rec == 8) {
    uint64 S;
    in.read(static_cast<char*>(static_cast<void*>(&S)), sizeof(uint64));
    if(S != sizeof(GadgetHeader)) {
      swap_bytes(S);
      if(S == sizeof(GadgetHeader)) swap = 1;
      else return false;
    }
  } else
    throw falcON::exception("Fortran header size must be 4 or 8\n");
  // read full GadgetHeader
  in.read(static_cast<char*>
	  (static_cast<void*>(this)), sizeof(GadgetHeader));
  // if required swap bytes
  if(swap) {
    swap_bytes(npart,6);
    swap_bytes(masstab,6);
    swap_bytes(time);
    swap_bytes(redshift);
    swap_bytes(flag_sfr);
    swap_bytes(flag_feedback);
    swap_bytes(npartTotal,6);
    swap_bytes(flag_cooling);
    swap_bytes(num_files);
    swap_bytes(BoxSize);
    swap_bytes(Omega0);
    swap_bytes(OmegaLambda);
    swap_bytes(HubbleParam);
    swap_bytes(flag_stellarage);
    swap_bytes(flag_metals);
    swap_bytes(npartTotalHighWord,6);
    swap_bytes(flag_entropy_instead_u);
  }
  // read record trailer and check for consistency
  if(rec == 4) {
    uint32 S;
    in.read(static_cast<char*>(static_cast<void*>(&S)), sizeof(uint32));
    if(swap) swap_bytes(S);
    if(S != sizeof(GadgetHeader)) {
      falcON_Warning("GadgetHeader::Read(): record size mismatch\n");
      return false;
    }
  } else if(rec == 8) {
    uint64 S;
    in.read(static_cast<char*>(static_cast<void*>(&S)), sizeof(uint64));
    if(swap) swap_bytes(S);
    if(S != sizeof(GadgetHeader)) {
      falcON_Warning("GadgetHeader::Read(): record size mismatch\n");
      return false;
    }
  }
  return true;
}
//------------------------------------------------------------------------------
bool falcON::GadgetHeader::mismatch(GadgetHeader const&H) const {
  bool okay = true;
#define CHECK_I(FIELD,FIELDNAME)				\
  if(FIELD != H.FIELD) {					\
    okay = false;						\
    falcON_Warning("GadgetHeader \"%s\" mismatch (%u vs %u)\n",	\
		   FIELDNAME, FIELD, H.FIELD);			\
  }
#define CHECK_D(FIELD,FIELDNAME)				\
  if(FIELD != H.FIELD) {					\
    okay = false;						\
    falcON_Warning("GadgetHeader\"%s\" mismatch (%f vs %f)\n",	\
		   FIELDNAME, FIELD, H.FIELD);			\
  }
  CHECK_D(time,"time");
  CHECK_D(redshift,"redshift");
  CHECK_I(flag_sfr,"flag_sfr");
  CHECK_I(flag_feedback,"flag_feedback");
  CHECK_I(npartTotal[0],"npartTotal[0]");
  CHECK_I(npartTotal[1],"npartTotal[1]");
  CHECK_I(npartTotal[2],"npartTotal[2]");
  CHECK_I(npartTotal[3],"npartTotal[3]");
  CHECK_I(npartTotal[4],"npartTotal[4]");
  CHECK_I(npartTotal[5],"npartTotal[5]");
  CHECK_I(flag_cooling,"flag_cooling");
  CHECK_I(num_files,"num_files");
  CHECK_D(BoxSize,"BoxSize");
  CHECK_D(Omega0,"Omega0");
  CHECK_D(OmegaLambda,"OmegaLambda");
  CHECK_D(HubbleParam,"HubbleParam");
  CHECK_I(flag_stellarage,"flag_stellarage");
  CHECK_I(flag_metals,"flag_metals");
  CHECK_I(flag_entropy_instead_u,"flag_entropy_instead_u");
  CHECK_I(npartTotalHighWord[0],"npartTotalHighWord[0]");
  CHECK_I(npartTotalHighWord[1],"npartTotalHighWord[1]");
  CHECK_I(npartTotalHighWord[2],"npartTotalHighWord[2]");
  CHECK_I(npartTotalHighWord[3],"npartTotalHighWord[3]");
  CHECK_I(npartTotalHighWord[4],"npartTotalHighWord[4]");
  CHECK_I(npartTotalHighWord[5],"npartTotalHighWord[5]");
  return !okay;
#undef CHECK_I
#undef CHECK_D
}
//------------------------------------------------------------------------------
void falcON::GadgetHeader::check_simple_npart_error() {
  for(int k=0; k!=6; ++k)
    if(npart[k] > npartTotal[k]) {
      falcON_Warning("GadgetHeader: npart[%u]=%u > npartTotal[%u]=%u: "
		     "we will try to fix by setting npartTotal[%u]=%u\n",
		     k,npart[k],k,npartTotal[k],k,npart[k]);
      npartTotal[k] = npart[k];
    }
}
//------------------------------------------------------------------------------
void falcON::GadgetHeader::dump(std::ostream&out) const {
  out<<" gadget header dump:";
  for(int k=0; k!=6; ++k)
    out<<"\n type "<<k
       <<": npart="<<std::setw(8)<<npart[k]
       <<" npartTotal="<<std::setw(8)<<npartTotal[k]
       <<" masstab="<<masstab[k];
  out<<"\n redshift               = "<<redshift
     <<"\n flag_sfr               = "<<flag_sfr
     <<"\n flag_feedback          = "<<flag_feedback
     <<"\n flag_cooling           = "<<flag_cooling
     <<"\n num_files              = "<<num_files
     <<"\n BoxSize                = "<<BoxSize
     <<"\n Omega0                 = "<<Omega0
     <<"\n OmegaLambda            = "<<OmegaLambda
     <<"\n HubbleParam            = "<<HubbleParam
     <<"\n flag_stellarage        = "<<flag_stellarage
     <<"\n flag_metals            = "<<flag_metals
     <<"\n flag_entropy_instead_u = "<<flag_entropy_instead_u
     <<std::endl;
}
////////////////////////////////////////////////////////////////////////////////
