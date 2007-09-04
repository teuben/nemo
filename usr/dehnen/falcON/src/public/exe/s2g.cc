// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// s2g.cc                                                                      |
//                                                                             |
// Copyright (C) 2007 Walter Dehnen                                            |
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
// NEMO to Gadget snapshot file converter                                      |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 0.0    14/02/2007  WD created & seriously debugged                        |
// v 0.1    14/02/2007  WD added param header                                  |
// v 1.0    16/02/2007  WD moved to public part of falcON                      |
// v 1.0.1  23/02/2007  WD renamed to s2g (previously nemo2gadget)             |
// v 1.0.2  05/03/2007  WD avoid warning from bodies::read_snapshot()          |
// v 1.0.3  30/08/2007  WD debugged error message                              |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "1.0.3"
#define falcON_VERSION_D "30-aug-2007 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#include <body.h>                                  // bodies                    
#include <public/io.h>                             // NEMO file I/O             
#include <main.h>                                  // main & NEMO stuff         
#include <fstream>                                 // C++ file I/O              
#include <iomanip>     
#include <string>
////////////////////////////////////////////////////////////////////////////////
string defv[] = {	
  "in=???\n           input file (in nemo snapshot format)               ",
  "out=???\n          base name for output files (out000 out001 ...)     ", 
  "copy=mxvkU\n       data to copy (minimum mxvkU, maximum mxvkURHpa)    ",
  "times=first\n      time range(s) to select snapshots from;\n"
  "                    if equals \"first\", write only one file \"out\",\n"
  "                    otherwise allow for multiple time steps           ",
  "first=0\n          index for first snapshot to write (if more than 1) ",
  "warn=f\n           warn about missing data (which will be zeroed)     ",
  "header=4\n         header size for unfio (4 or 8)                     ",
  falcON_DEFV, NULL};
string usage="NEMO to GADGET converter (better than \"nemo2gadget\")";
////////////////////////////////////////////////////////////////////////////////
void falcON::main() falcON_THROWING
{
  const char* fbase=getparam("out");
  nemo_in in(getparam("in"));
  if(!in) falcON_THROW("cannot open file \"%s\" for input\n",getparam("in"));
  fieldset copy = getioparam("copy") & fieldset("mxvkURHpa"), read;
  snapshot shot;
  if(0 == strcmp(getparam("times"), "first")) {
    if(!in.has_snapshot())
      falcON_THROW("file \"%s\" contains no snapshot\n",getparam("in"));
    shot.read_nemo(in,read,copy,0,0);
    output out(fbase);
    if(!out) falcON_THROW("cannot open file \"%s\" for output\n",fbase);
    shot.write_gadget(out,copy,getbparam("warn"),getuparam("header"));
  } else {
    unsigned ifile = getuparam("first");
    char file[256];
    strcpy(file,fbase);
    strcat(file,"%03d");
    output out;
    while(in.has_snapshot()) {
      if(! shot.read_nemo(in,read,copy,getparam("times"),0) ) continue;
      out.reopen(file,ifile++);
      if(!out) falcON_THROW("cannot open file \"%s\" for output\n",out.file());
      shot.write_gadget(out,copy,getbparam("warn"),getuparam("header"));
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
