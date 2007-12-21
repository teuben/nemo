// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/proper/exe/addprop.cc                                           
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
// v 0.0   20/12/2007  WD created.                                              
// v 0.1   21/12/2007  WD allow for replacement, changes in bodyfunc            
////////////////////////////////////////////////////////////////////////////////
#define falcON_VERSION   "0.1"
#define falcON_VERSION_D "21-dec-2007 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile "addprop"
#endif
#define falcON_RepAction 0                         // no action reporting       
//-----------------------------------------------------------------------------+
#include <body.h>                                  // bodies etc..              
#include <public/io.h>                             // WDs C++ NEMO I/O          
#include <public/bodyfunc.h>                       // body functions            
#include <main.h>                                  // NEMO basics & main        
#include <cstdio>                                  // C std I/O                 
extern "C" {
#  include  <stdinc.h>                             // for nemoinpd, nemoinpf    
}
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n         snapshot input file                                ",
  "out=???\n        snapshot output file                               ",
  "times=all\n      times to process                                   ",
  "write=\n         select data to write out (default: all)            ",
  "add=???\n        single character: body field to add                ",
  "value=???\n      bodyfunc: value to assign to property added        ",
  "pars=\n          parameters, if any, for value                      ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string usage = "addprop -- add/replace a single body property to/in snapshots";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING {
  nemo_in        IN(getparam("in"));
  nemo_out       OUT;
  Bodyfunc       FUNC(getparam("value"),getparam("pars"));
  real           P[bodyfunc::MAXPAR], *PARS = 0;
  fieldset       NEED = FUNC.need();
  fieldbit       ADD(getparam("add")[0]);
  // sanity checks for added property and the value to be assigned to it
  if(ADD == fieldbit::invalid)
    error("property '%c' not recognised \n",getparam("add")[0]);
  // loop snapshots ...
  fieldset WRITE = getioparam_a("write") | fieldset(ADD);
  fieldset GET   = NEED | getioparam_a("write"), GOT;
  snapshot SHOT;
  while(IN.has_snapshot()) {
    if(!SHOT.read_nemo(IN,GOT,GET,getparam("times"),0)) continue;
    SHOT.add_field(ADD);
    if(!GOT.contain(NEED)) {
      fieldset miss = GOT.missing(NEED);
      error("data '%s' required for value are missing in snapshot\n",
	    word(miss));
    }
    switch(value(ADD)) {
#ifdef ASSIGN_PROP
#  error C-macro ASSIGN_PROP already #defined
#endif
#define ASSIGN_PROP(BIT,NAME)			\
    case BIT: {					\
      BodyProp<BIT> PROP(&FUNC,SHOT.time());	\
      LoopAllBodies(&SHOT,B)			\
	B.NAME() = PROP(B);			\
    } break;
      DEF_NAMED(ASSIGN_PROP);
#undef ASSIGN_PROP
    }
    if(!OUT) OUT.open(getparam("out"));
    SHOT.write_nemo(OUT,WRITE);
  }
}
