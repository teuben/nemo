// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// manipulate.cc                                                               |
//                                                                             |
// Copyright (C) 2005 Walter Dehnen                                            |
//                                                                             |
// This program is free software; you can redistribute it and/or modify        |
// it under the terms of the GNU General Public License as published by        |
// the Free Software Foundation; either version 2 of the License, or (at       |
// your option) any later version.                                             |
//                                                                             |
// This program is distributed in the hope that it will be useful, but         |
// WITHOUT ANY WARRANTY; without even the implied warranty of                  |
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU           |
// General Public License for more details.                                    |
//                                                                             |
// You should have received a copy of the GNU General Public License           |
// along with this program; if not, write to the Free Software                 |
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                   |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 1.0   28/04/2005  WD created                                              |
// v 1.1   20/05/2005  WD new manip.cc; added manippath option                 |
// v 2.0   14/06/2005  WD new falcON                                           |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "2.0"
#define falcON_VERSION_D "13-jun-2005 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile "manipulate"
#endif
#define falcON_RepAction 0                         // no action reporting       
#define falcON_NOT_USING_PROPER                    // not using PROPRIETARY code
//-----------------------------------------------------------------------------+
#include <public/manip.h>                          // N-body manipulators       
#include <main.h>                                  // main & NEMO stuff         
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n           input file                                         ",
  "out=\n             output file (default: none)                        ",
  "manipname=\n       name(s) of run-time manipulator(s)                 ",
  "manippars=\n       parameter set(s) for manipulator(s)                ",
  "manipfile=\n       data file(s) required by manipulator(s)            ",
  "manippath=\n       path to search for manipulators first              ",
  "times=all\n        times to process                                   ",
  "write=\n           what to write (if output desired; default: all)    ",
  "stop=f\n           stop if manipulator(s) returns 0                   ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string usage = "manipulate -- use manipulators on nemo snapshots";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  const bool stop(getbparam("stop"));
  nemo_in in(getparam("in"));
  nemo_out out;
  const Manipulator manip(getparam  ("manipname"),
			  getparam_z("manippars"),
			  getparam_z("manipfile"),
			  getparam_z("manippath"));
  if(!manip)
    falcON_THROW("empty manipulator");
  const fieldset need ( manip.need() );
  const fieldset write( hasvalue("write")?
			fieldset(getparam("write")) :
			need|manip.provide());
  snapshot shot;
  bool goon (true);
  while(goon && in.has_snapshot()) {
    fieldset read;
    if(!shot.read_nemo(in,read,need,getparam("times"),0)) continue;
    check_sufficient(read, need);
    shot.add_fields(manip.provide());
    goon = manip(shot) || !stop;
    if(hasvalue("out")) {
      if(!out && !out.open(getparam("out")))
	falcON_THROW("cannot open file \"%s\" for output",getparam("out"));
      shot.write_nemo(out,write);
    }
  }
}
//------------------------------------------------------------------------------
  
  
