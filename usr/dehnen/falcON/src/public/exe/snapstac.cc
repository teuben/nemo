// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// snapstac.cc                                                                 |
//                                                                             |
// Copyright (C) 2005  Walter Dehnen                                           |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// This is a non-public part of the code.                                      |
// It is property of its author and not to be made public without his written  |
// consent.                                                                    |
//                                                                             |
//-----------------------------------------------------------------------------+
// history:                                                                    |
//                                                                             |
// v 0.0    21/09/2005 WD created                                              |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "0.0"
#define falcON_VERSION_D "22-sep-2005 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile "snapfilter"
#endif
#define falcON_RepAction 0                         // no action reporting       
//-----------------------------------------------------------------------------+
#include <body.h>                                  // bodies & I/O              
#include <main.h>                                  // main & NEMO stuff         
using namespace falcON;
//------------------------------------------------------------------------------
string defv[] = {
  "in1=???\n          input file name                                    ",
  "in2=???\n          input file name                                    ",
  "out=???\n          output file name                                   ",
  "deltax=0,0,0\n     position of in1 wrt in2                            ",
  "deltav=0,0,0\n     velocity of in1 wrt in2                            ",
  "zerocm=f\n         zero center of mass                                ",
  "write=\n           which data to write [default: all read]            ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string
usage = "snapstack -- stack two N-body systems on top of each other\n";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  fieldset       got(getioparam_a("write"));
  const fieldset write(getioparam_z("write"));
  nemo_out       out(getparam("out"));
  nemo_in        in1(getparam("in1")), in2(getparam("in2"));
  if(!in1.has_snapshot())
    falcON_THROW("cannot load snapshot from file %s\n",getparam("in1"));
  if(!in2.has_snapshot())
    falcON_THROW("cannot load snapshot from file %s\n",getparam("in2"));
  snap_in        snap1(in1), snap2(in2);
  const double   time(snap1.has_time()? snap1.time() :
		      snap2.has_time()? snap2.time() : 0.);
  snapshot       shot(time,
		      snap1.Nbod()+snap2.Nbod(),
		      fieldset::o,
		      snap1.Nsph()+snap2.Nsph());
  // read snapshots
  const body     b1(shot.begin_all_bodies());
  got &= shot.read_nemo(snap1, got, b1, 0, 0);
  const body     b2(b1, snap1.Nbod());
  got &= shot.read_nemo(snap2, got, b2, 0, 0);
  if(write && !got.contain(write))
    warning("couldn't read %s from both %s and %s\n",
	    word(got.missing(write)), getparam("in1"), getparam("in2"));
  // shift centre of snapshot 1
  const vect dx(getvparam("deltax")), dv(getvparam("deltav"));
  if(dx != zero || dv != zero)
    for(body b(b1); b!=b2; ++b) {
      b.pos() += dx;
      b.vel() += dv;
    }
  // center to centre of mass
  if(getbparam("zerocm")) {
    real m0(zero);
    vect x0(zero), v0(zero);
    LoopAllBodies(&shot,b) {
      m0 += mass(b);
      x0 += pos(b) * mass(b);
      v0 += vel(b) * mass(b);
    }
    if(m0) {
      x0 /= m0;
      v0 /= m0;
      LoopAllBodies(&shot,b) {
	b.pos() -= x0;
	b.vel() -= v0;
      }
    }
  }
  // output
  if(out)
    shot.write_nemo(out,got);
}
//------------------------------------------------------------------------------

