// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/public/exe/snapstac.cc                                          
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2005-2011                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2005-2011 Walter Dehnen                                        
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
// v 0.0    21/09/2005 WD created                                               
// v 0.1    10/09/2007 WD fixed bug with history buffer                         
// v 0.2    19/09/2007 WD changes in fields.h, body.h, io.h; keyword time       
// v 0.2.1  10/09/2008 WD happy gcc 4.3.1
// v 0.2.2  15/06/2011 WD keyword mirror
// v 0.3    30/06/2011 WD keywords dx1,dx2,dv1,dv2 replace deltax and deltav
////////////////////////////////////////////////////////////////////////////////
#define falcON_VERSION   "0.3"
#define falcON_VERSION_D "30-jun-2011 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile "snapstac"
#endif
#define falcON_RepAction 0                         // no action reporting       
//-----------------------------------------------------------------------------+
#include <main.h>                                  // main & NEMO stuff         
using namespace falcON;
//------------------------------------------------------------------------------
const char*defv[] = {
  "in1=???\n          input file name                                    ",
  "in2=???\n          input file name                                    ",
  "out=???\n          output file name                                   ",
  "dx1=0,0,0\n        shift positions  of snapshot 1 by dx1              ",
  "dv1=0,0,0\n        shift velocities of snapshot 1 by dv1              ",
  "dx2=0,0,0\n        shift positions  of snapshot 2 by dx1              ",
  "dv2=0,0,0\n        shift velocities of snapshot 2 by dv1              ",
  "mirror=f\n         reflect snapshot 2 w.r.t. origin (before shifting) ",
  "time=\n            set simulation time                                ",
  "zerocm=f\n         zero center of mass (after all shifting etc.)      ",
  "write=\n           which data to write [default: all read]            ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
const char*usage = "snapstac -- stack two N-body systems on top of each other";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  fieldset       got(getioparam_a("write"));
  const fieldset write(getioparam_z("write"));
  nemo_in        in1(getparam("in1")), in2(getparam("in2"));
  if(!in1.has_snapshot())
    falcON_THROW("cannot load snapshot from file %s\n",getparam("in1"));
  if(!in2.has_snapshot())
    falcON_THROW("cannot load snapshot from file %s\n",getparam("in2"));
  snap_in        snap1(in1), snap2(in2);
  const double   time(snap1.has_time()? snap1.time() :
		      snap2.has_time()? snap2.time() : 0.);
  unsigned nbod[bodytype::NUM] = {0};
  for(bodytype t; t; ++t)
    nbod[t] = snap1.Nbod(t) + snap2.Nbod(t);
  snapshot shot(time, nbod, fieldset::empty);
  // read snapshots
  const body b1(shot.begin_all_bodies());
  got &= shot.read_part(snap1, got, b1, 0);
  const body b2(b1, snap1.Ntot());
  got &= shot.read_part(snap2, got, b2, 0);
  if(write && !got.contain(write))
    falcON_Warning("couldn't read %s from both %s and %s\n",
		   word(got.missing(write)), getparam("in1"), getparam("in2"));
  // shift centre of snapshot 1
  vect dx, dv;
  dx = getvparam("dx1");
  dv = getvparam("dv1");
  if(dx != zero || dv != zero)
    for(body b(b1); b!=b2; ++b) {
      b.pos() += dx;
      b.vel() += dv;
    }
  // mirror snapshot 2
  if(getbparam("mirror"))
    for(body b(b2); b!=shot.end_all_bodies(); ++b) {
      b.pos() = -pos(b);
      b.vel() = -vel(b);
    }
  // shift centre of snapshot 2
  dx = getvparam("dx2");
  dv = getvparam("dv2");
  if(dx != zero || dv != zero)
    for(body b(b2); b!=shot.end_all_bodies(); ++b) {
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
  if(hasvalue("time"))
    shot.set_time(getdparam("time"));
  // output
  nemo_out out(getparam("out"));
  if(out) shot.write_nemo(out,got);
}
//------------------------------------------------------------------------------

