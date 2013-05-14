// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/src/io.cc                                                
///
/// \brief  contains definitions of methods declared in utils/inc/io.h         
///
/// \author Walter Dehnen
/// \date   2000-2009,2013
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2009,2013 Walter Dehnen
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
void output::open_std()
{
  if( ++openstdout > 1 )
    WDutils_THROW("trying to open more than one output to stdout");
}
void output::close_std()
{
  if(openstdout) --openstdout;
}
void input::open_std()
{
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
  std::ifstream::pos_type end_pos = f.tellg();
  return size_t(end_pos > begin_pos? end_pos-begin_pos : 0);
}
////////////////////////////////////////////////////////////////////////////////
//
// class WDutils::output
//
////////////////////////////////////////////////////////////////////////////////
void output::__open(bool append)
{
  DebugInfo(8,"output::__open(%d): FILE=%s\n",append,FILE);
  APPENDING = false;
  if     (0 == FILE    ||
	  0 == FILE[0] ||
	  0 == std::strcmp(FILE,".") ) {
    OUT = 0;
    DebugInfo(5,"output: open sink\n");
  } else if(0 == std::strcmp(FILE,"-") ) {
    open_std();
    OUT = &std::cout;
    DebugInfo(5,"output: open stdout\n");
  } else {
    DebugInfo(10,"output::__open(%d): FILE=%s\n",append,FILE);
    std::ofstream *FOUT = new std::ofstream();
    if(append) {
      FOUT->open(FILE,std::ios::out | std::ios::app);
      if(FOUT->is_open()) {
	APPENDING = true;
	DebugInfo(4,"output: append to file \"%s\"\n",FILE);
      }
    }
    if(!FOUT->is_open() )
      FOUT->open(FILE,std::ios::out);
    if( FOUT->is_open() ) {
      OUT = FOUT;
      DebugInfo(5,"output: open file \"%s\"\n",FILE);
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
    DebugInfo(6,"output: closing\n");
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
size_t WDutils::FortranIRec::read_size() throw(WDutils::exception)
{
  if(HSZE == 4) {
    uint32_t S;
    IN.read(reinterpret_cast<char*>(&S),sizeof(uint32_t));
    if(SWAP) swap_bytes(S);
    return size_t(S);
  } else if(HSZE == 8) {
    uint64_t S;
    IN.read(reinterpret_cast<char*>(&S),sizeof(uint64_t));
    if(SWAP) swap_bytes(S);
    return size_t(S);
  } else 
    throw exception("FortranIRec: header size must be 4 or 8\n");
#ifndef __PGI
  return 0; // unreachable; avoids warnings about not returning a value
#endif
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
  DebugInfo(6,"FortranIRec: opened with %lu bytes\n",SIZE);
}
//------------------------------------------------------------------------------
size_t WDutils::FortranIRec::read_bytes(char*buf, size_t n)
  throw(WDutils::exception)
{
  if(!IN) throw exception("FortranIRec::read_bytes(): input corrupted");
  if(READ + n > SIZE) {
    WDutils_Warning("FortranIRec::read(): "
		    "can only read %lu bytes, not %lu\n",SIZE-READ,n);
    n = SIZE - READ;
  }
  IN.read(buf,n);
  if(!IN) throw exception("FortranIRec: input corrupted");
  READ += n;
  DebugInfo(6,"FortranIRec: read %lu bytes\n",n);
  return n;
}
//------------------------------------------------------------------------------
void WDutils::FortranIRec::skip_bytes(size_t n) {
  if(READ + n > SIZE) n = SIZE - READ;
  if(n && !IN) throw exception("FortranIRec::skip_bytes(): input corrupted");
  for(char C; n; --n,++READ) IN.read(&C,1);
}
//------------------------------------------------------------------------------
void WDutils::FortranIRec::close() throw(WDutils::exception)
{
  if(!IN) throw exception("FortranIRec::close(): input corrupted");
  if(READ != SIZE) {
    WDutils_Warning("FortranIRec: only %lu of %lu bytes read "
		    "on closing record\n", READ, SIZE);
    for(char C; READ!=SIZE; ++READ) IN.read(&C,1);
  }
  size_t S = read_size();
  IN.FREC = 0;
  if(S != SIZE) throw exception("FortranIRec: record size mismatch");
  DebugInfo(6,"FortranIRec: closed with %lu bytes\n",SIZE);
}
////////////////////////////////////////////////////////////////////////////////
//
// class WDutils::FortranORec
//
////////////////////////////////////////////////////////////////////////////////
void WDutils::FortranORec::write_size() throw(WDutils::exception)
{
  if(HSZE == 4) {
    uint32_t S = uint32_t(SIZE);
    OUT.write(reinterpret_cast<const char*>(&S),sizeof(uint32_t));
  } else if(HSZE == 8) {
    uint64_t S = uint64_t(SIZE);
    OUT.write(reinterpret_cast<const char*>(&S),sizeof(uint64_t));
  } else 
    throw exception("FortranORec: header size must be 4 or 8\n");
}
//------------------------------------------------------------------------------
WDutils::FortranORec::FortranORec(output& out, size_t rsize, unsigned rec)
  throw(WDutils::exception) : OUT(out), HSZE(rec), SIZE(rsize), WRITTEN(0)
{
  if(!OUT) throw exception("FortranORec: output corrupted");
  if(OUT.FREC)
    throw exception("trying to open 2nd FortranORec to same output\n");
  OUT.FREC = this;
  write_size();
  DebugInfo(6,"FortranORec: opened for %lu bytes\n",SIZE);
}
//------------------------------------------------------------------------------
size_t WDutils::FortranORec::write_bytes(const char*buf, size_t n)
  throw(WDutils::exception)
{
  if(!OUT) throw exception("FortranORec: output corrupted");
  if(WRITTEN + n > SIZE) {
    WDutils_Warning("FortranORec::write(): can only write %lu bytes, not %lu\n",
		    SIZE-WRITTEN,n);
    n = SIZE - WRITTEN;
  }
  OUT.write(buf,n);
  if(!OUT) throw exception("FortranORec: ostream corrupted");
  WRITTEN += n;
  DebugInfo(6,"FortranORec: written %lu bytes\n",n);
  return n;
}
//------------------------------------------------------------------------------
void WDutils::FortranORec::fill_bytes(size_t n, char C) {
  if(WRITTEN + n > SIZE) n = SIZE - WRITTEN;
  for(; n; --n,++WRITTEN) OUT.write(&C,1);
}
//------------------------------------------------------------------------------
void WDutils::FortranORec::close() throw(WDutils::exception)
{
  if(!OUT) throw exception("FortranORec: output corrupted");
  if(WRITTEN!=SIZE) {
    WDutils_Warning("FortranORec: only %lu of %lu bytes written "
		    "on closing record ... padding with 0\n", WRITTEN, SIZE);
    for(char C=0; WRITTEN!=SIZE; ++WRITTEN) OUT.write(&C,1);
  }
  write_size();
  OUT.FREC = 0;
  DebugInfo(6,"FortranORec: closed with %lu bytes\n",SIZE);
}
//------------------------------------------------------------------------------
