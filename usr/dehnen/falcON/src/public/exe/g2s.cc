// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// g2s.cc                                                                      |
//                                                                             |
// Copyright (C) 2007,2008 Walter Dehnen                                       |
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
// Gadget to NEMO snapshot file converter                                      |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 0.0    14/02/2007  WD created & seriously debugged                        |
// v 0.1    14/02/2007  WD added param header                                  |
// v 0.2    15/02/2007  WD can handle different endianess; added param first   |
// v 0.3    16/02/2007  WD checks on file size, only writes what it got        |
// v 1.0    16/02/2007  WD moved to public part of falcON                      |
// v 1.0.1  23/02/2007  WD renamed to g2s (previously gadget2nemo)             |
// v 1.0.2  10/09/2008  WD happy gcc 4.3.1                                     |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "1.0.2"
#define falcON_VERSION_D "10-sep-2008 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#include <main.h>                                  // main & NEMO stuff         
#include <fstream>                                 // C++ file I/O              
#include <iomanip>     
#include <string>
////////////////////////////////////////////////////////////////////////////////
const char*defv[] = {	
  "in=???\n           input base file (GADGET format)                    ",
  "out=???\n          output file (nemo snapshot format)                 ",
  "nshot=1\n          number of snapshots to read                        ",
  "copy=mxvU\n        data to copy (maximum: mxvkURHpa)                  ",
  "first=0\n          first snapshot to read                             ",
  "header=4\n         header size or unfio (4 or 8)                      ",
  falcON_DEFV, NULL};
const char*usage="GADGET to NEMO converter (better than \"gadget2nemo\")";
////////////////////////////////////////////////////////////////////////////////
void falcON::main() falcON_THROWING
{
  const char* fbase=getparam("in");
  unsigned nshot=getuparam("nshot");
  unsigned header=getuparam("header");
  nemo_out out(getparam("out"));
  fieldset copy=getioparam("copy") & fieldset("mxvkURHpa"), got;
  snapshot shot;
  if(nshot > 1) {
    for(unsigned i=0; i!=nshot; ++i) {
      char file[256];
      sprintf(file, "%s%03d", fbase, i+getuparam("first"));
      got = shot.read_gadget(file, copy, header);
      shot.write_nemo(out, copy & got);
    }
  } else {
    got = shot.read_gadget(fbase, copy, header);
    shot.write_nemo(out, copy & got);
  }
}
////////////////////////////////////////////////////////////////////////////////
