// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// scale_eps.cc                                                                |
//                                                                             |
// Copyright (C) 2003-2008 Walter Dehnen                                       |
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
// v 1.0    23/10/2003  WD created                                             |
// v 1.1    07/22/2003  WD marginal changes due to changes in body.h           |
// v 1.2    30/04/2004  WD more marginal changes due to changes in body.h      |
// v 2.0    19/05/2004  WD option give (previously write), no default.         |
// v 2.1    20/05/2005  WD several minor updates, option give -> write         |
// v 3.0    14/06/2005  WD new falcON                                          |
// v 3.0.1  10/09/2008  WD happy gcc 4.3.1                                     |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "3.0.1"
#define falcON_VERSION_D "10-sep-2008 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#error You need NEMO to compile "scale_eps"
#endif
#define falcON_RepAction 0                         // no action reporting       
#define falcON_NOT_USING_PROPER                    // not using PROPRIETARY code
//-----------------------------------------------------------------------------+
#include <iostream>                                // C++ I/O                   
#include <fstream>                                 // C++ file I/O              
#include <cmath>                                   // standard math             
#include <main.h>                                  // main & NEMO stuff         
//------------------------------------------------------------------------------
const char*defv[] = {
  "in=???\n           input file                                         ",
  "out=???\n          output file                                        ",
  "fac=???\n          scale factor for eps_i                             ",
  "write=\n           output only these; if not given, output all we got ",
  "times=all\n        times of snapshots to scale                        ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
const char*usage = "scale_eps -- scales individual softening lengths";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  nemo_in         in (getparam("in"));
  nemo_out        out;
  fieldset        read, write;
  const fieldset  give (getioparam_a("write"));
  snapshot        shot;
  const real      fac(getrparam("fac"));
  while(in.has_snapshot()) {
    if(! shot.read_nemo(in,read,give,getparam("times"),0)) continue;
    if(! read.contain(give))
    falcON_Warning("cannot write '%s' data (only read '%s')",
		   word(give & ~read), word(read));
    write = read & give;
    if(write.contain(fieldbit::e))
      LoopAllBodies(&shot,Bi) Bi.eps() *= fac;
    if(!out.is_open()) out.open(getparam("out"));
    shot.write_nemo(out,write);
  }
}
