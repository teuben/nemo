// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   src/public/exe/manipulate.cc
///
/// \author Walter Dehnen
///
/// \date   2005-2008,2011
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005-2008,2011 Walter Dehnen
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
//
// history: 
//
// v 1.0   28/04/2005  WD created
// v 1.1   20/05/2005  WD new manip.cc; added manippath option
// v 2.0   14/06/2005  WD new falcON
// v 2.1   27/07/2006  WD write all data out if write not given
// v 2.1.1 10/09/2008  WD happy gcc 4.3.1
// v 2.2   22/06/2011  WD allow for times=first and times=last
////////////////////////////////////////////////////////////////////////////////
#define falcON_VERSION   "2.2"
#define falcON_VERSION_D "22-jun-2011 Walter Dehnen                          "
//
#ifndef falcON_NEMO                                // this is a NEMO program
#  error You need NEMO to compile "manipulate"
#endif
#define falcON_RepAction 0                         // no action reporting
#define falcON_NOT_USING_PROPER                    // not using PROPRIETARY code
//
#include <public/manip.h>                          // N-body manipulators
#include <main.h>                                  // main & NEMO stuff
//
const char*defv[] = {
  "in=???\n           input file                                         ",
  "out=\n             output file (default: none)                        ",
  "manipname=\n       name(s) of run-time manipulator(s)                 ",
  "manippars=\n       parameter set(s) for manipulator(s)                ",
  "manipfile=\n       data file(s) required by manipulator(s)            ",
  "manippath=\n       path to search for manipulators first              ",
  "times=all\n        time range(s) to copy (\"first\", \"last\" are allowed)",
  "write=\n           what to write (if output desired; default: all)    ",
  "stop=f\n           stop if manipulator(s) returns 0                   ",
  falcON_DEFV, NULL };
//
const char*usage = "manipulate -- use manipulators on nemo snapshots";
//
namespace {
  using namespace falcON;
  struct ManipOut {
    const Manipulator MANIP;
    nemo_out          OUT;
    const fieldset    WRITE;
    const fieldset    WANT;
    const bool        STOP;
    //
    ManipOut()
      : MANIP( getparam  ("manipname"), getparam_z("manippars"),
	       getparam_z("manipfile"), getparam_z("manippath")),
	WRITE( getioparam_a("write") ),
	WANT ( MANIP.need() | WRITE ),
	STOP ( getbparam("stop") )
    { if(!MANIP) falcON_THROW("empty manipulator"); }
    //
    bool operator()(snapshot&shot, fieldset read)
    {
      check_sufficient(read, MANIP.need());
      shot.add_fields(MANIP.provide());
      bool stop = MANIP(shot);
      if(hasvalue("out") && strcmp(getparam("out"),".") ) {
	if(!OUT && !OUT.open(getparam("out")))
	  falcON_THROW("cannot open file \"%s\" for output",getparam("out"));
	shot.write_nemo(OUT,WRITE);
      }
      return stop && STOP;
    }
  };
}
//
void falcON::main() falcON_THROWING
{
  nemo_in  IN(getparam("in"));
  ManipOut MANO;
  snapshot SHOT;
  fieldset READ;
  if(!IN.has_snapshot())
    falcON_THROW("no snapshots found in input file\n");
  if(0==strcmp(getparam("times"),"first")) {
    // special case times=first
    SHOT.read_nemo(IN,READ,MANO.WANT,0,0);
    MANO(SHOT,READ);
  } else if(0==strcmp(getparam("times"),"last")) {
    // special case times=last
    while(IN.has_snapshot())
      SHOT.read_nemo(IN,READ,MANO.WANT,0,0);
    SHOT.del_fields(~READ);
    MANO(SHOT,READ);
  } else {
    // general case for times
    while(IN.has_snapshot()) {
      if(!SHOT.read_nemo(IN,READ,MANO.WANT,getparam("times"),0)) continue;
      if(MANO(SHOT,READ)) break;
    }
  }
}
//
  
  
