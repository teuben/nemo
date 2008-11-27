// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/src/io.cc                                                
///
/// \brief  contains definitions of methods declared in utils/inc/io.h         
///
/// \author Walter Dehnen
/// \date   2000-2008
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2008 Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but          
// WITHOUT ANY WARRANTY; without even the implied warranty of                   
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
// General Public License for more details.                                     
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 675
// Mass Ave, Cambridge, MA 02139, USA.
//
////////////////////////////////////////////////////////////////////////////////
#include <io.h>
#include <memory.h>
#include <fstream>
#include <cstdlib>
#include <cstring>

using namespace WDutils;
////////////////////////////////////////////////////////////////////////////////
namespace {
  //----------------------------------------------------------------------------
  int openstdout = 0;
  int openstdin  = 0;
}
void output::open_std() WDutils_THROWING {
  if( ++openstdout > 1 )
    WDutils_THROW("trying to open more than one output to stdout");
}
void output::close_std() {
  if(openstdout) --openstdout;
}
void input::open_std() WDutils_THROWING {
  if( ++openstdin > 1 )
    WDutils_THROW("trying to open more than one input from stdin");
}
void input::close_std() {
  if(openstdin) --openstdin;
}
// /////////////////////////////////////////////////////////////////////////////
//
// WDutils::FileSize()
//
// /////////////////////////////////////////////////////////////////////////////
size_t WDutils::FileSize(const char*sFileName)
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
//
// class WDutils::output
//
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
      WDutils_DEL_O(FOUT);
    }
  }
}
//------------------------------------------------------------------------------
void output::close() {
  if(FREC) {
    if(FILE)
      WDutils_Warning("closing FortranORec before output from file \"%s\"\n",
		     FILE);
    else
      WDutils_Warning("closing FortranORec before output\n");
    FREC->close();
  }
  if(OUT) {
    DebugInfo(2,"output: closing\n");
    if(OUT == &std::cout) close_std();
    else WDutils_DEL_O(OUT);
  }
  APPENDING = false;
  OUT = 0;
}
////////////////////////////////////////////////////////////////////////////////
//
// class WDutils::input
//
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
      WDutils_DEL_O(FIN);
    }
  }
}
//------------------------------------------------------------------------------
void input::close() {
  if(FREC) {
    if(FILE)
      WDutils_Warning("closing FortranIRec before input from file \"%s\"\n",
		     FILE);
    else
      WDutils_Warning("closing FortranIRec before input\n");
    FREC->close();
  }
  DebugInfo(2,"input: closing\n");
  if(IN == &std::cin) close_std();
  else if(IN) WDutils_DEL_O(IN);
  IN = 0;
}
////////////////////////////////////////////////////////////////////////////////
//
// class WDutils::FortranIRec
//
////////////////////////////////////////////////////////////////////////////////
unsigned WDutils::FortranIRec::read_size() throw(WDutils::exception)
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
WDutils::FortranIRec::FortranIRec(input& in, unsigned rec, bool swap)
  throw(WDutils::exception) : IN(in), HSZE(rec), SWAP(swap), READ(0)
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
unsigned WDutils::FortranIRec::read_bytes(char*buf, unsigned n)
  throw(WDutils::exception)
{
  if(!IN) throw exception("FortranIRec::read_bytes(): input corrupted");
  if(READ + n > SIZE) {
    WDutils_Warning("FortranIRec::read(): cannot read %u, but only %u bytes\n",
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
void WDutils::FortranIRec::skip_bytes(unsigned n) {
  if(READ + n > SIZE) n = SIZE - READ;
  if(n && !IN) throw exception("FortranIRec::skip_bytes(): input corrupted");
  for(char C; n; --n,++READ) IN.read(&C,1);
}
//------------------------------------------------------------------------------
void WDutils::FortranIRec::close() throw(WDutils::exception)
{
  if(!IN) throw exception("FortranIRec::close(): input corrupted");
  if(READ != SIZE) {
    WDutils_Warning("FortranIRec: only %u of %u bytes read on closing record\n",
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
// class WDutils::FortranORec
//
////////////////////////////////////////////////////////////////////////////////
void WDutils::FortranORec::write_size() throw(WDutils::exception)
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
WDutils::FortranORec::FortranORec(output& out, unsigned size, unsigned rec)
  throw(WDutils::exception) : OUT(out), HSZE(rec), SIZE(size), WRITTEN(0)
{
  if(!OUT) throw exception("FortranORec: output corrupted");
  if(OUT.FREC)
    throw exception("trying to open 2nd FortranORec to same output\n");
  OUT.FREC = this;
  write_size();
  DebugInfo(6,"FortranORec: opened for %u bytes\n",SIZE);
}
//------------------------------------------------------------------------------
unsigned WDutils::FortranORec::write_bytes(const char*buf, unsigned n)
  throw(WDutils::exception)
{
  if(!OUT) throw exception("FortranORec: output corrupted");
  if(WRITTEN + n > SIZE) {
    WDutils_Warning("FortranORec::write(): cannot read %u, but only %u bytes\n",
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
void WDutils::FortranORec::fill_bytes(unsigned n, char C) {
  if(WRITTEN + n > SIZE) n = SIZE - WRITTEN;
  for(; n; --n,++WRITTEN) OUT.write(&C,1);
}
//------------------------------------------------------------------------------
void WDutils::FortranORec::close() throw(WDutils::exception)
{
  if(!OUT) throw exception("FortranORec: output corrupted");
  if(WRITTEN!=SIZE) {
    WDutils_Warning("FortranORec: only %u of %u bytes written on closing record"
		   " ... padding with 0\n", WRITTEN, SIZE);
    for(char C=0; WRITTEN!=SIZE; ++WRITTEN) OUT.write(&C,1);
  }
  write_size();
  OUT.FREC = 0;
  DebugInfo(6,"FortranORec: closed with %u bytes\n",SIZE);
}
//------------------------------------------------------------------------------
