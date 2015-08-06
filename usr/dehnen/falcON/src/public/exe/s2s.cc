// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   src/public/exe/s2s.cc
///
/// \author Walter Dehnen
///
/// \date   2007-12,15
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007-2012 Walter Dehnen 
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
// v 0.0   27/02/2007  WD created
// v 0.1   28/02/2007  WD renamed keyword 'write' -> 'copy'
// v 0.2   09/10/2007  WD added keyword 'time' 
// v 0.3   30/10/2007  WD allowed times=first and times=last
// v 1.0   31/10/2007  WD added filter, etc, makes snapfilter redundant
// v 1.0.1 25/03/2008  WD warn if no snapshot matched times
// v 1.0.2 10/09/2008  WD happy gc 4.3.1
// v 2.0   09/03/2010  WD sorting, new filter (in body.h)
// v 2.0.1 22/06/2011  WD check for snapshot first
// v 2.1   04/12/2012  WD added rotaxis and rotangle
// v 2.2   28/01/2015  WD added dx,dv
////////////////////////////////////////////////////////////////////////////////
#define falcON_VERSION   "2.2"
#define falcON_VERSION_D "28-jan-2015 Walter Dehnen                          "
//
#ifndef falcON_NEMO                                // this is a NEMO program
#  error You need NEMO to compile "s2s"
#endif
#define falcON_RepAction 0                         // no action reporting
//
#include <public/nemo++.h>                         // WDs C++ NEMO I/O
#include <public/bodyfunc.h>                       // body functions
#include <utils/matr33.h>                          // matrix
#include <main.h>                                  // NEMO basics & main
#include <cstdio>                                  // C std I/O
//
const char*defv[] = {
  "in=???\n         snapshot input file                                ",
  "out=???\n        snapshot output file                               ",
  "times=all\n      time range(s) to copy (\"first\", \"last\" are allowed)",
  "dx=0,0,0\n       shift positions  of snapshot by dx                 ",
  "dv=0,0,0\n       shift velocities of snapshot by dv                 ",
  "filter=\n        boolean bodyfunc expression (man page): filter     ",
  "params=\n        parameters, must match requirements for filter     ",
  "sorting=\n       scalar bodyfunc expression: property to sort       ",
  "sortpars=\n      parameters, must match requirements for sorting    ",
  "zeromissing=f\n  set body properties missing for filter to zero?    ",
  "rotaxis=\n       rotate snapshot around this axis                   ",
  "rotangle=0\n     angle [degrees] to rotate snapshot                 ",
  "time=\n          set time of snapshots to this value (at output)    ",
  "copy=\n          select data to write out (default: all read)       ",
  falcON_DEFV, NULL };
const char*usage = "s2s -- Walter's alternative to snapcopy";
//
typedef WDutils::Matrix33<falcON::real> Matrix;
template<typename T> struct rotate_helper {
  static void job(T&, Matrix const&) {}
};
template<> struct rotate_helper<falcON::vect> {
  static void job(falcON::vect&x, Matrix const&R) { x = R*x; }
};
template<typename T> inline
void rotate_it(T&x, Matrix const&R)
{
  rotate_helper<T>::job(x,R);
}
//
template<int BIT> struct rotate_body {
  static void job(falcON::fieldset r, falcON::body&B, Matrix const&R)
  {
    if(r.contain(falcON::fieldbit(BIT)))
      rotate_it(B.template datum<BIT>(), R);
    rotate_body<BIT-1>::job(r,B,R);
  }
};
template<> struct rotate_body<0> {
  static void job(falcON::fieldset, falcON::body&, Matrix const&) {}
};
//
inline void rotate(falcON::fieldset r, falcON::body&B, Matrix const&R)
{
  rotate_body<falcON::fieldbit::NQUANT-1>::job(r,B,R);
}
//
namespace {
  using namespace falcON;
  struct FilterSortRotateShiftWrite
  {
    BodyFilter     Filter;
    BodyFunc<real> SortFunc;
    fieldset       Copy,Need;
    nemo_out       Out;
    bool           keys,zm,rot;
    vect           dx,dv;
    Matrix33<real> RotMat;

    FilterSortRotateShiftWrite()
      : Filter(getparam_z("filter"),getparam_z("params"))
      , SortFunc(getparam_z("sorting"),getparam_z("sortpars"))
      , Copy(getioparam_a("copy"))
      , Need(Copy | Filter.need() | SortFunc.need())
      , keys(getioparam_z("copy").contain(fieldbit::k))
      , zm(getbparam("zeromissing"))
      , rot(hasvalue("rotaxis"))
      , dx(getvparam("dx"))
      , dv(getvparam("dv"))
    { 
      if(keys) Copy |= fieldset::k;
      if(getrparam("rotangle")==0) {
	if(rot)
	  falcON_Warning("rotaxis given but rotangle=%g\n",
			 getrparam("rotangle"));
	rot = 0;
      } else if(!rot)
	falcON_Warning("rotangle =%g but no rotaxis given\n",
		       getrparam("rotangle"));
      if(rot) {
	vect u = getvparam("rotaxis");
	if(norm(u)==0) {
	  rot = 0;
	  falcON_Warning("rotaxis=0,0,0: cannot rotate\n");
	  return;
	}
	u /= abs(u);
	real the = getrparam("rotangle")*Pi/180;
	real cth = std::cos(the);
	real sth = std::sin(the);
	RotMat.dyadic(u,u);
	RotMat    *= 1-cth;
	RotMat(0) += cth;
	RotMat(1) -= u[2]*sth;
	RotMat(2) += u[1]*sth;
	RotMat(3) += u[2]*sth;
	RotMat(4) += cth;
	RotMat(5) -= u[0]*sth;
	RotMat(6) -= u[1]*sth;
	RotMat(7) += u[0]*sth;
	RotMat(8) += cth;
      }
    }
    void operator()(snapshot&shot)
    {
      if(keys) shot.add_field(fieldbit::k);
      shot.apply_filter(Filter,zm);
      if(shot.N_bodies()) {
	shot.apply_sort(SortFunc,Copy,zm);
	if(hasvalue("time")) shot.set_time(getdparam("time"));
	if(rot) {
	  fieldset VecCopy = Copy & fieldset::vectors;
	  VecCopy &= shot.all_data();
	  LoopAllBodies(&shot,B)
	    rotate(VecCopy,B,RotMat);
	}
	if(norm(dx))
	  LoopAllBodies(&shot,B)
	    B.pos() += dx;
	if(norm(dv))
	  LoopAllBodies(&shot,B)
	    B.vel() += dv;
	if(!Out) Out.open(getparam("out"));
	shot.write_nemo(Out,Copy);
      }
    }
  };
}
//
void falcON::main() falcON_THROWING
{
  nemo_in  In(getparam("in"));
  fieldset Read;
  snapshot Shot;
  FilterSortRotateShiftWrite FSRSW;
  if(!In.has_snapshot())
    falcON_THROW("no snapshots found in input file\n");
  if(0==strcmp(getparam("times"),"first")) {
    // special case times=first
    Shot.read_nemo(In,Read,FSRSW.Need,0,0);
    FSRSW(Shot);
  } else if(0==strcmp(getparam("times"),"last")) {
    // special case times=last
    while(In.has_snapshot())
      Shot.read_nemo(In,Read,FSRSW.Need,0,0);
    Shot.del_fields(~Read);
    FSRSW(Shot);
  } else {
    // general case for times
    while(In.has_snapshot())
      if(Shot.read_nemo(In,Read,FSRSW.Need,getparam("times"),0)) {
	Shot.del_fields(~Read);
	FSRSW(Shot);
      }
  }
  if(!FSRSW.Out)
    falcON_Warning("no snapshot matching \"times=%s\" found in input\n",
		   getparam("times"));
}
//
