// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   mksingle.cc
///
/// \brief  creates 1-body snapshot
///
/// \author Walter Dehnen
///
/// \date   2011
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 675
// Mass Ave, Cambridge, MA 02139, USA.
//
////////////////////////////////////////////////////////////////////////////////
//
// history:
//
// v 0.0    28/03/2011  WD  created
////////////////////////////////////////////////////////////////////////////////
#define falcON_VERSION   "0.0"
#define falcON_VERSION_D "28-mar-2011 Walter Dehnen                          "
//
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile mksingle
#endif
#define falcON_RepAction 0                         // no action reporting
//
#include <main.h>                                  // main & NEMO stuff
////////////////////////////////////////////////////////////////////////////////
const char*defv[] = {
  "out=???\n          output file                                        ",
  "m=???\n            mass of particle                                   ",
  "pos=0,0,0\n        position of particle                               ",
  "vel=0,0,0\n        velocity of particle                               ",
  "sink=t\n           sink (or std) particle?                            ",
  falcON_DEFV, NULL };
//
const char*usage = "mksingle -- create single-particle snapshot";
////////////////////////////////////////////////////////////////////////////////
void falcON::main() falcON_THROWING
{
  unsigned N[bodytype::NUM] = {0u};
  N[getbparam("sink")? bodytype::sink : bodytype::std] = 1u;
  snapshot shot(0.0,N);
  body b = shot.begin_all_bodies();
  b.mass() = getrparam("m");
  b.pos () = getvrparam("pos");
  b.vel () = getvrparam("vel");
  nemo_out out(getparam("out"));
  shot.write_nemo(out);
}
