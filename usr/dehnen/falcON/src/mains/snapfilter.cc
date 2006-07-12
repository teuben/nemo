// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/mains/snapfilter.cc                                             
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2005-2006                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2005-2006 Walter Dehnen                                        
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
// v 0.0    21/07/2005 WD created as testbed for bodies::remove()               
// v 0.1    28/04/2006 WD bodyfunc replaced by BodyFunc<>                       
// v 1.0    04/07/2006 WD made public (along with bodyfunc)                     
// v 1.0.1  06/07/2006 WD slight changes in bodyfunc error handling             
////////////////////////////////////////////////////////////////////////////////
#define falcON_VERSION   "1.0.1"
#define falcON_VERSION_D "06-jul-2006 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile "snapfilter"
#endif
#define falcON_RepAction 0                         // no action reporting       
//-----------------------------------------------------------------------------+
#include <public/bodyfunc.h>                       // body functions            
#include <main.h>                                  // main & NEMO stuff         
using namespace falcON;
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n           input file                                         ",
  "out=???\n          output file for filtered bodies                    ",
  "filter=???\n       boolean bodyfunc expression (man page): filter     ",
  "params=\n          parameters, must match requirements from filter    ",
//   "remain=\n          output file for remaining bodies                   ",
  "write=\n           which data to write [default: all read]            ",
  "times=all\n        times to process                                   ",
  "zeromissing=f\n    zero missing body properties (or error out)?       ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string
usage = "snapfilter -- filters bodies from snapshot";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  // set up parameters                                                          
  const bool     Z    (getbparam("zeromissing"));
  const nemo_in  IN   (getparam("in"));
  BodyFunc<bool> BF   (getparam_z("filter"), getparam_z("params"));
  const fieldset WRITE(getioparam_a("write"));
  const fieldset WANT (WRITE | (BF? BF.need() : fieldset(fieldset::o)) );
  snapshot       SHOT (0., 0u, WANT | fieldset::f);
  nemo_out       OUT;
  // loop snapshots and process them                                            
  while(IN.has_snapshot()) {
    fieldset GOT;
    if(!SHOT.read_nemo(IN,GOT,WANT,getparam("times"),0)) continue;
    if(!GOT.contain(BF.need())) {
      fieldset MISS = GOT.missing(BF.need());
      if(Z) {
	warning("data '%s' missing; will zero them",word(MISS));
	SHOT.reset_data(MISS);
      } else
	falcON_THROW("data on '%s' missing",word(MISS));
    }
    SHOT.add_field(fieldbit::f);
    SHOT.reset_flags();
    unsigned COUNT = 0;
    if(BF) LoopAllBodies(&SHOT,B)
      if(! BF(B,SHOT.time())) {
	B.flag_for_removal();
	COUNT++;
      }
    if(COUNT) {
      SHOT.remove();
      debug_info(1,"%d bodies removed at time %f\n",COUNT,SHOT.time());
    }
    if(SHOT.N_bodies() > 0) {
      if(OUT || OUT.open(getparam("out")))
	SHOT.write_nemo(OUT,GOT&WRITE);
    } else
      warning("no body filtered at time %f",SHOT.time());
  }
}
//------------------------------------------------------------------------------
