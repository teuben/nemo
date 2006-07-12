// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// mover.cc                                                                    |
//                                                                             |
// Copyright (C) 2002-2005 Walter Dehnen                                       |
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
// v 1.0   18/11/2002  PJM created from density_centre input/output            |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "3.1"
#define falcON_VERSION_D "13-jul-2005 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile "mover"
#endif
#define falcON_RepAction 0                         // no action reporting       
#define falcON_NOT_USING_PROPER                    // not using PROPRIETARY code
//-----------------------------------------------------------------------------+
#include <iostream>                                // C++ I/O                   
#include <fstream>                                 // C++ file I/O              
#include <iomanip>                                 // C++ I/O formatting        
#include <body.h>                                  // the bodies                
#include <public/io.h>                             // my NEMO I/O               
#include <public/tools.h>                          // my tools for bodies       
#include <utils/inline_io.h>                      // my I/O utilities          
#include <main.h>                                  // main & NEMO stuff         
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n           input file                                         ",
  "out=\n             output file for snapshots                          ",
  "give=\n            output only these; if not given, output all we got ",
  "xrot=\n            rotation about x axis, done first  (degrees)       ",
  "yrot=\n            rotation about y axis, done second (degrees)       ",
  "zrot=\n            rotation about z axis, done third  (degrees)       ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string
usage = "mover - rotates snapshots\n";
//------------------------------------------------------------------------------

void rotatex (const double _cxr,
	      const double _sxr,
	      double &_y,  double &_z,
	      double &_vy, double &_vz){

  register double newy=_cxr*_y + _sxr*_z, newvy=_cxr*_vy + _sxr*_vz;

  _z  = _cxr*_z  - _sxr*_y;
  _vz = _cxr*_vz - _sxr*_vy;
  _y  = newy;
  _vy = newvy;
}

void rotatey (const double _cyr,
	      const double _syr,
	      double &_x,  double &_z,
	      double &_vx, double &_vz){

  register double newx=_cyr*_x - _syr*_z, newvx=_cyr*_vx - _syr*_vz;

  _z  = _cyr*_z  + _syr*_x;
  _vz = _cyr*_vz + _syr*_vx;
  _x  = newx;
  _vx = newvx;
}

void rotatez (const double _czr,
	      const double _szr,
	      double &_x,  double &_y,
	      double &_vx, double &_vy){

  register double newx=_czr*_x + _szr*_y, newvx=_czr*_vx + _szr*_vy;

  _y  = _czr*_y  - _szr*_x;
  _vy = _czr*_vy - _szr*_vx;
  _x  = newx;
  _vx = newvx;
}

void falcON::main() falcON_THROWING
{
  // open I/O                                                                   
  const bool     DO_OUT   (file_for_output("out"));
  const nemo_in  IN   (getparam("in"));
  const fieldset GIVE (DO_OUT   ? getioparam_a("give") : fieldset::o);
  nemo_out       OUT;
  // loop snapshots & process them                                              
  snapshot       SHOT;
  fieldset       READ;
  double xr=0.017453292*(getdparam("xrot")),
    yr=0.017453292*(getdparam("yrot")),
    zr=0.017453292*(getdparam("zrot"));

  double cxr=cos(xr),sxr=sin(xr),
    cyr=cos(yr),syr=sin(yr),
    czr=cos(zr),szr=sin(zr);

  while(IN.has_snapshot()) {
    // read time, read snapshot if in times                                     
    const bool IN_TIMES = SHOT.read_nemo(IN,READ,GIVE);

    // output of snapshot                                                       
    if(DO_OUT) {
      if(!OUT.is_open()) {
	OUT.open(getparam_z("out"));
	if(!OUT) ::error("cannot open out\n");
      }
      LoopAllBodies(&SHOT,Bi) {
	double tmpx=Bi.pos()[0], tmpy=Bi.pos()[1], tmpz=Bi.pos()[2],
	  tmpvx=Bi.vel()[0], tmpvy=Bi.vel()[1], tmpvz=Bi.vel()[2];

	if(xr) rotatex(cxr,sxr,tmpy,tmpz,tmpvy,tmpvz);
	if(yr) rotatey(cyr,syr,tmpx,tmpz,tmpvx,tmpvz);
	if(zr) rotatey(czr,szr,tmpx,tmpy,tmpvx,tmpvy);

	Bi.pos()[0]  = tmpx;
	Bi.pos()[1]  = tmpy;
	Bi.pos()[2]  = tmpz;
	Bi.vel()[0]  = tmpvx;
	Bi.vel()[1]  = tmpvy;
	Bi.vel()[2]  = tmpvz;
      }
      SHOT.write_nemo(OUT,READ&GIVE);
    }
  }
}
