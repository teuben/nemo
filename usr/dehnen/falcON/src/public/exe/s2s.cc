// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   src/public/exe/s2s.cc
///
/// \author Walter Dehnen
/// \date   2007-2010
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
// v 0.0   27/02/2007  WD created
// v 0.1   28/02/2007  WD renamed keyword 'write' -> 'copy'
// v 0.2   09/10/2007  WD added keyword 'time' 
// v 0.3   30/10/2007  WD allowed times=first and times=last
// v 1.0   31/10/2007  WD added filter, etc, makes snapfilter redundant
// v 1.0.1 25/03/2008  WD warn if no snapshot matched times
// v 1.0.2 10/09/2008  WD happy gc 4.3.1
// v 2.0   09/03/2010  WD sorting, new filter (in body.h)
////////////////////////////////////////////////////////////////////////////////
#define falcON_VERSION   "2.0"
#define falcON_VERSION_D "09-mar-2010 Walter Dehnen                          "
//
#ifndef falcON_NEMO                                // this is a NEMO program
#  error You need NEMO to compile "s2s"
#endif
#define falcON_RepAction 0                         // no action reporting
//
#include <public/nemo++.h>                         // WDs C++ NEMO I/O
#include <public/bodyfunc.h>                       // body functions
#include <main.h>                                  // NEMO basics & main
#include <cstdio>                                  // C std I/O
//
const char*defv[] = {
  "in=???\n         snapshot input file                                ",
  "out=???\n        snapshot output file                               ",
  "times=all\n      time range(s) to copy (\"first\", \"last\" are allowed)",
  "filter=\n        boolean bodyfunc expression (man page): filter     ",
  "params=\n        parameters, must match requirements for filter     ",
  "sorting=\n       scalar bodyfunc expression: property to sort       ",
  "sortpars=\n      parameters, must match requirements for sorting    ",
  "zeromissing=f\n  set body properties missing for filter to zero?    ",
  "time=\n          set time of snapshots to this value (at output)    ",
  "copy=\n          select data to write out (default: all read)       ",
  falcON_DEFV, NULL };
const char*usage = "s2s -- Walter's alternative to snapcopy";
//
namespace {
  using namespace falcON;
  struct FilterSortWrite
  {
    BodyFilter     Filter;
    BodyFunc<real> SortFunc;
    fieldset       Copy,Need;
    nemo_out       Out;
    bool           zm;

    FilterSortWrite()
      : Filter(getparam_z("filter"),getparam_z("params")),
	SortFunc(getparam_z("sorting"),getparam_z("sortpars")),
	Copy(getioparam_a("copy")),
	Need(Copy | Filter.need() | SortFunc.need()),
	zm(getbparam("zeromissing"))
    {}
    void operator()(snapshot&shot)
    {
      if(Need.contain(fieldbit::k)) shot.add_field(fieldbit::k);
      shot.apply_filter(Filter,zm);
      if(shot.N_bodies()) {
	shot.apply_sort(SortFunc,Copy,zm);
	if(hasvalue("time")) shot.set_time(getdparam("time"));
	if(!Out) Out.open(getparam("out"));
      }
    }
  };
}
//
void falcON::main() falcON_THROWING {
  nemo_in         In(getparam("in"));
  fieldset        Read;
  snapshot        Shot;
  FilterSortWrite FSW;
  if(0==strcmp(getparam("times"),"first")) {
    // special case times=first
    if(In.has_snapshot()) {
      Shot.read_nemo(In,Read,FSW.Need,0,0);
      FSW(Shot);
    }
  } else if(0==strcmp(getparam("times"),"last")) {
    // special case times=last
    if(In.has_snapshot()) {
      while(In.has_snapshot())
	Shot.read_nemo(In,Read,FSW.Need,0,0);
      Shot.del_fields(~Read);
      FSW(Shot);
    }
  } else {
    // general case for times
    while(In.has_snapshot())
      if(Shot.read_nemo(In,Read,FSW.Need,getparam("times"),0)) {
	Shot.del_fields(~Read);
	FSW(Shot);
      }
  }
  if(!FSW.Out)
    falcON_Warning("no snapshot matching \"times=%s\" found in input\n",
		   getparam("times"));
}
