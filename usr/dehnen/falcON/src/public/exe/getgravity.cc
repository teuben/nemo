// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// getgravity.cc                                                               |
//                                                                             |
// Copyright (C) 2002-2010 Walter Dehnen                                       |
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
// v 0.0    23/11/2002  WD created.                                            |
// v 0.1    04/02/2003  WD default falcON parameters automized                 |
// v 0.2    20/03/2003  WD gravity, action reporting                           |
// v 0.3    23/05/2003  WD automated NEMO history                              |
// v 1.0    20/05/2005  WD several minor updates                               |
// v 2.0    14/06/2005  WD new falcON, new body.h, new nemo I/O                |
// v 2.1    22/06/2005  WD changes in nemo I/O support                         |
// v 2.2    13/06/2005  WD changes in fieldset                                 |
// v 2.3    13/06/2005  WD changes in fieldset and body.h                      |
// v 2.3.1  19/09/2007  WD ??                                                  |
// v 2.3.2  20/02/2008  WD more changes in body.h                              |
// v 2.3.3  10/09/2008  WD happy gcc 4.3.1                                     |
// v 2.4    19/03/2010  WD sink particle gravity extra tree, epssink, no fsink |
// v 2.4.1  25/03/2010  WD debugged sink particle gravity, fsink back
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "2.4.1"
#define falcON_VERSION_D "25-mar-2010 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO
#error You need NEMO to compile "src/mains/getgravity.cc"
#endif
#include <main.h>                                  // main & NEMO stuff         
#include <forces.h>                                // falcON forces             
//------------------------------------------------------------------------------
const char*defv[] = {
  "srce=???\n         input file: sources [m,x]       ",
  "sink=???\n         input file: sinks   [x]         ",
  "out=???\n          output file         [x,a,p]     ",
  "times=all\n        time range (for srce only)      ",
  "eps=0.05\n         softening length                ",
#ifdef falcON_PROPER
  "epssink=\n         softening length for sink particles (default: eps) ",
  "fsink=0.2\n        theta_sink/theta <= 1                              ",
#endif
  "kernel=" falcON_KERNEL_TEXT
  "\n                 softening kernel                ",
  "theta=" falcON_THETA_TEXT
  "\n                 tolerance parameter at M=M_tot  ",
  "Ncrit=" falcON_NCRIT_TEXT
  "\n                 max # bodies in un-split cells  ",
  "Grav=1\n           Newton's constant of gravity (0-> no self-gravity) ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
const char*usage =
    "getgravity -- computes gravity at sink positions; using falcON";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  nemo_in        srce(getparam("srce"));
  nemo_out       out(getparam("out"));
  const fieldset write(fieldset::x | fieldset::p | fieldset::a);
  snapshot       shot(fieldset::gravity);
  bool           SOFT(getrparam("eps") < 0);
  real           TH(getrparam("theta"));
  forces         falcon(&shot,
			getrparam("eps"),
			TH,
			kern_type(getiparam("kernel")),
			SOFT,
			getrparam("Grav"),
			TH< zero? const_theta : theta_of_M,
#ifdef falcON_PROPER
			getrparam_z("epssink"),
			getrparam  ("fsink")
#else
			zero, one
#endif
      );
  const fieldset srcedata(fieldset::m|fieldset::x);
  while(srce.has_snapshot()) {
    // open snapshot with sources and check for time in range, if both given
    snap_in srce_in(srce);
    if(srce_in.has_time() && !time_in_range(srce_in.time(),getparam("times")))
      continue;
    // open snapshot with sinks and ensure we have enough bodies
    nemo_in sink   (getparam("sink"));
    snap_in sink_in(sink);
    unsigned nbod[bodytype::NUM] = {0};
    for(bodytype t; t; ++t)
      nbod[t] = srce_in.Nbod(t) + sink_in.Nbod(t);
    shot.resetN(nbod);
    // read sources
    const body sources(shot.begin_all_bodies());
    fieldset read = shot.read_part(srce_in, srcedata, sources);
    if(!read.contain(srcedata))
      falcON_THROW("sources must have mx data");
    // read sinks
    const body sinks(sources, srce_in.Ntot());
    read = shot.read_part(sink_in, fieldset::x, sinks);
    if(!read.contain(fieldbit::x))
      falcON_THROW("sinks must have x data");
    shot.set_time(srce_in.time());
    // loop sources, get their mass and flag them to be inactive
    real M(zero);
    for(body b(sources); b!=sinks; ++b) {
      b.unflag_active();
      M += mass(b);
    }
    // loop sinks, set their mass and flag them to be active
    M *= 1.e-10/real(srce_in.Ntot());
    for(body b(sinks); b; ++b) {
      b.flag_as_active();
      b.mass() = M;
    }
    // compute gravity
    falcon.grow(getiparam("Ncrit"));
    falcon.approximate_gravity();
    // write sink data to output
    if(out)
      shot.write_nemo(out,fieldset(fieldset::x|fieldset::a|fieldset::p),sinks);
  }
}
