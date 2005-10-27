// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// set_time.cc                                                                 |
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
// sets time of snapshots                                                      |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 1.0.0  16/05/2005  WD created                                             |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "1.6"
#define falcON_VERSION_D "16-may-2005 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#error You need NEMO to compile addgravity
#endif
#define falcON_RepAction 0                         // no action reporting       
//-----------------------------------------------------------------------------+
#include <body.h>                                  // bodies etc..              
#include <public/io.h>                             // WDs C++ NEMO I/O          
#include <main.h>                                  // main & NEMO stuff         
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n           input file                                         ",
  "out=???\n          output file                                        ",
  "dtime=???\n        add this to every snapshots' time                  ",
  "times=all\n        time range                                         ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string usage = "addgravity -- adds gravity to a snapshot; using falcON";
//------------------------------------------------------------------------------
void falcON::main() throw(falcON::exception)
{
  nemo_in  IN (getparam("in"));
  nemo_out OUT(getparam("out"));
  snapshot SHOT;
  io       READ;
  double   DT (getdparam("dtime"));
  while(IN.is_present(nemo_io::snap)) {
    if(! SHOT.read_nemo(IN, READ, io::all, getparam("times"), 0))
      continue;
    SHOT.advance_time_by(DT);
    SHOT.write_nemo(OUT,READ);
  }
}
