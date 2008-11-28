// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   src/public/exe/addprop.c
///
/// \author Walter Dehnen
/// \date   2007-2008
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007-2008 Walter Dehnen
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
// v 0.0   20/12/2007  WD created 
// v 0.1   21/12/2007  WD allow for replacement, changes in bodyfunc
// v 0.1.1 10/09/2008  WD happy gcc 4.3.1
// v 0.1.2 27/11/2008  WD WDutils::io.h
////////////////////////////////////////////////////////////////////////////////
#define falcON_VERSION   "0.1.2"
#define falcON_VERSION_D "27-nov-2008 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile "addprop"
#endif
#define falcON_RepAction 0                         // no action reporting       
//-----------------------------------------------------------------------------+
#include <main.h>                                  // NEMO basics & main        
#include <public/bodyfunc.h>                       // body functions            
#include <cstdio>                                  // C std I/O                 
//------------------------------------------------------------------------------
const char*defv[] = {
  "in=???\n         snapshot input file                                ",
  "out=???\n        snapshot output file                               ",
  "times=all\n      times to process                                   ",
  "write=\n         select data to write out (default: all)            ",
  "add=???\n        single character: body field to add                ",
  "value=???\n      bodyfunc: value to assign to property added        ",
  "pars=\n          parameters, if any, for value                      ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
const char*usage =
    "addprop -- add/replace a single body property to/in snapshots";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING {
  nemo_in        IN(getparam("in"));
  nemo_out       OUT;
  Bodyfunc       FUNC(getparam("value"),getparam("pars"));
  fieldset       NEED = FUNC.need();
  fieldbit       ADD(getparam("add")[0]);
  // sanity checks for added property and the value to be assigned to it
  if(ADD == fieldbit::invalid)
    falcON_Error("property '%c' not recognised \n",getparam("add")[0]);
  // loop snapshots ...
  fieldset WRITE = getioparam_a("write") | fieldset(ADD);
  fieldset GET   = NEED | getioparam_a("write"), GOT;
  snapshot SHOT;
  while(IN.has_snapshot()) {
    if(!SHOT.read_nemo(IN,GOT,GET,getparam("times"),0)) continue;
    SHOT.add_field(ADD);
    if(!GOT.contain(NEED)) {
      fieldset miss = GOT.missing(NEED);
      falcON_Error("data '%s' required for value are missing in snapshot\n",
		   word(miss));
    }
    switch(value(ADD)) {
#undef ASSIGN_PROP
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
