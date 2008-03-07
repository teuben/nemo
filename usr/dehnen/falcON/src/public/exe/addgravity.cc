// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// addgravity.cc                                                               |
//                                                                             |
// Copyright (C) 2002-2008 Walter Dehnen                                       |
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
// adds gravity (pot & acc) to a snapshot                                      |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 1.0.0  19/10/2001  WD created. initial version reads all snapshots.       |
// v 1.0.1  23/10/2001  WD added option times= to select time range(s)         |
// v 1.0.2  28/08/2002  WD adapted to various changes in falcON                |
// v 1.1    30/08/2002  WD adapted this file for usage of MPI otherwise        |
// v 1.1.1  04/02/2003  WD default falcON parameters automized                 |
// v 1.2    20/03/2003  WD gravity, action reporting                           |
// v 1.3    23/05/2003  WD automated NEMO history                              |
// v 1.3.1  23/05/2003  WD automated version & compile information             |
// v 1.4    06/05/2004  WD new body.h; write=read + gravity                    |
// v 1.5    25/08/2004  WD allowing for individual softening lengths           |
// v 1.6    16/05/2005  WD added external potential                            |
// v 2.0    13/06/2005  WD new falcON, new body.h, new nemo I/O                |
// v 2.1    13/06/2005  WD changes in fieldset                                 |
// v 2.2    06/03/2008  WD debugged (problem when using external potential)    |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "2.2"
#define falcON_VERSION_D "06-mar-2008 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#error You need NEMO to compile addgravity
#endif
//-----------------------------------------------------------------------------+
#include <body.h>                                  // bodies etc..              
#include <forces.h>                                // falcON                    
#include <public/io.h>                             // WDs C++ NEMO I/O          
#include <externacc.h>                             // external potential        
#include <main.h>                                  // main & NEMO stuff         
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n           input file                                         ",
  "out=???\n          output file                                        ",
  "times=all\n        time range                                         ",
  "eps=0.05\n         softening length [ <0: use eps_i from snapshot]    ",
  "kernel="falcON_KERNEL_TEXT
  "\n                 softening kernel                                   ",
  "theta="falcON_THETA_TEXT
  "\n                 tolerance parameter at M=M_tot                     ",
  "Ncrit="falcON_NCRIT_TEXT
  "\n                 max # bodies in un-split cells                     ",
  "Grav=1\n           Newton's constant of gravity (0-> no self-gravity) ",
  "root_center=\n     if given (3 numbers), forces tree-root centering   ",
  "accname=\n         name of external acceleration field                ",
  "accpars=\n         parameters of external acceleration field          ",
  "accfile=\n         file required by external acceleration field       ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string usage = "addgravity -- adds gravity to a snapshot; using falcON";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  nemo_in  IN (getparam("in"));
  nemo_out OUT;
  unsigned NCRIT(getiparam("Ncrit"));
  fieldset       READ, NEED(fieldset::m|fieldset::x);
#ifdef falcON_INDI
  bool     SOFT(getrparam("eps") < 0);
  if(SOFT) NEED |= fieldset::e;
#endif
  vect     X0, *RC(getvparam_z("root_center",X0));
  snapshot SHOT;
  forces   FALCON(&SHOT,
		  getrparam("eps"),
		  getrparam("theta"),
		  kern_type(getiparam("kernel")),
#ifdef falcON_INDI
		  SOFT,
#endif
		  getrparam("Grav") );
  acceleration *ACCEXT = hasvalue("accname") ?
    new nemo_acc(getparam("accname"), getparam("accpars"), getparam("accfile"))
    : 0;
  if(ACCEXT) SHOT.add_fields(fieldset::q);
  while(IN.has_snapshot()) {
    if(! SHOT.read_nemo(IN, READ, fieldset::all, getparam("times"), 0))
      continue;
    if(READ.contain(NEED)) {
      if(!OUT.is_open()) OUT.open(getparam("out"));
      SHOT.add_fields(fieldset::gravity);
      if(ACCEXT) SHOT.add_fields(fieldset::q);
      SHOT.reset_flags();
      if(FALCON.NewtonsG() != zero) {
	FALCON.grow(NCRIT, RC);
	FALCON.approximate_gravity(1,1);
      }
      if(ACCEXT)
	ACCEXT->set(&SHOT,true, FALCON.NewtonsG()? 2 : 0);
      SHOT.write_nemo(OUT,READ|fieldset(fieldset::a|fieldset::p|fieldset::q));
    } else
      warning("data '%s' missing at time %f: cannot compute gravity",
	      word(READ.missing(NEED)), SHOT.time());
  }
}
