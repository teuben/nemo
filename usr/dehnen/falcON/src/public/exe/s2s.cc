// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/public/exe/s2s.cc                                               
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2007                                                                
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2007 Walter Dehnen                                             
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
// v 0.0   27/02/2007  WD created.                                              
// v 0.1   28/02/2007  WD renamed keyword 'write' -> 'copy'                     
////////////////////////////////////////////////////////////////////////////////
#define falcON_VERSION   "0.1"
#define falcON_VERSION_D "28-feb-2007 Walter Dehnen                          "
//------------------------------------------------------------------------------
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile "s2s"
#endif
#define falcON_RepAction 0                         // no action reporting       
//------------------------------------------------------------------------------
#include <body.h>                                  // bodies etc..              
#include <public/io.h>                             // WDs C++ NEMO I/O          
#ifdef falcON_PROPER
#  include <public/bodyfunc.h>                     // body functions            
#endif
#include <main.h>                                  // NEMO basics & main        
#include <cstdio>                                  // C std I/O                 
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n         snapshot input file                                ",
  "out=???\n        snapshot output file                               ",
  "times=all\n      times to copy                                      ",
  "copy=\n          select data to write out (default: all)            ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string usage = "s2s -- Walter's alternative to snapcopy";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING {
  nemo_in        IN   (getparam("in"));
  nemo_out       OUT  (getparam("out"));
  const fieldset KEEP (getioparam_a("copy"));
  fieldset       READ;
  snapshot       SHOT;
  while(IN.has_snapshot())
    if(SHOT.read_nemo(IN,READ,KEEP,getparam("times"),0))
      SHOT.write_nemo(OUT,READ);
}
