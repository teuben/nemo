// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/public/exe/s2s.cc                                               
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
// v 0.0   27/02/2007  WD created.                                              
// v 0.1   28/02/2007  WD renamed keyword 'write' -> 'copy'                     
// v 0.2   09/10/2007  WD added keyword 'time'                                  
// v 0.3   30/10/2007  WD allowed times=first and times=last                    
// v 1.0   31/10/2007  WD added filter, etc, makes snapfilter redundant         
// v 1.0.1 25/03/2008  WD warn if no snapshot matched times                     
////////////////////////////////////////////////////////////////////////////////
#define falcON_VERSION   "1.0.1"
#define falcON_VERSION_D "25-mar-2008 Walter Dehnen                          "
//------------------------------------------------------------------------------
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile "s2s"
#endif
#define falcON_RepAction 0                         // no action reporting       
//------------------------------------------------------------------------------
#include <body.h>                                  // bodies etc..              
#include <public/io.h>                             // WDs C++ NEMO I/O          
#include <public/bodyfunc.h>                       // body functions            
#include <main.h>                                  // NEMO basics & main        
#include <cstdio>                                  // C std I/O                 
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n         snapshot input file                                ",
  "out=???\n        snapshot output file                               ",
  "times=all\n      time range(s) to copy (\"first\", \"last\" are allowed)",
  "filter=\n        boolean bodyfunc expression (man page): filter     ",
  "params=\n        parameters, must match requirements for filter     ",
  "zeromissing=f\n  set body properties missing for filter to zero?    ",
  "time=\n          set time of snapshots to this value (after filter) ",
  "copy=\n          select data to write out (default: all read)       ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string usage = "s2s -- Walter's alternative to snapcopy";
//------------------------------------------------------------------------------
namespace {
  using namespace falcON;
  using falcON::real;
  BodyFilter    *BF=0;    // filter to apply (if any)
  BodyFunc<real>*BS=0;    // expression to sort (if any)
  bool           ZM=0;    // set missing body data to zero?
  /// apply filter (if any) to snapshot, return true if any body remains
  /// \return did any body remain after filtering?
  /// \param ss (pter to) snapshot to be filtered
  /// \param got body data obtained from data file.
  bool apply_filter(snapshot*ss, fieldset got) {
    if(BF==0) return true;
    BF->set_time(ss->time());
    ss->add_field(fieldbit::f);
    ss->reset_flags();
    got |= fieldset::f;
    if(!got.contain(BF->need())) {
      fieldset miss = got.missing(BF->need());
      if(ZM) {
	warning("data '%s' required for filter but missing; "
		"will zero them\n",word(miss));
	ss->add_fields(miss);
	ss->reset_data(miss);
      } else
	error("data '%s' required for filter but missing "
	      "(use 'zeromissing=t' if you want me to reset them to 0)\n",
	      word(miss));
    }
    LoopAllBodies(ss,B)
      if(! (*BF)(B) ) B.flag_for_removal();
    ss->remove();
    debug_info(1,"%d bodies removed at time %f\n",ss->N_del(),ss->time());
    return ss->N_bodies() > 0;
  }
  void apply_sort(snapshot*ss, fieldset copy) {
    if(BS==0) return;
  }
}
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING {
  nemo_in  IN(getparam("in"));
  nemo_out OUT;
  fieldset COPY(getioparam_a("copy")), NEED(COPY), READ;
  snapshot SHOT;
  // initialize bodyfiltering
  if(hasvalue("filter")) {
    BF    = new BodyFilter(getparam("filter"),getparam_z("params"));
    ZM    = getbparam("zeromissing");
    NEED |= BF->need();
  }
  if(0==strcmp(getparam("times"),"first")) {
    // special case times=first
    if(IN.has_snapshot()) {
      SHOT.read_nemo(IN,READ,NEED,0,0);
      if(apply_filter(&SHOT,READ)) {
	if(hasvalue("time")) SHOT.set_time(getdparam("time"));
	if(!OUT) OUT.open(getparam("out"));
	SHOT.write_nemo(OUT,COPY);
      }
    } else
      warning("no snapshot found in input\n");
  } else if(0==strcmp(getparam("times"),"last")) {
    // special case times=last
    if(IN.has_snapshot()) {
      while(IN.has_snapshot())
	SHOT.read_nemo(IN,READ,NEED,0,0);
      if(apply_filter(&SHOT,READ)) {
	if(hasvalue("time")) SHOT.set_time(getdparam("time"));
	if(!OUT) OUT.open(getparam("out"));
	SHOT.write_nemo(OUT,COPY);
      }
    } else
      warning("no snapshot found in input\n");
  } else {
    // general case for times
    while(IN.has_snapshot())
      if(SHOT.read_nemo(IN,READ,NEED,getparam("times"),0)) {
	if(apply_filter(&SHOT,READ)) {
	  if(hasvalue("time")) SHOT.set_time(getdparam("time"));
	  if(!OUT) OUT.open(getparam("out"));
	  SHOT.write_nemo(OUT,COPY);
	}
      }
    if(!OUT)
      warning("no snapshot matching \"times=%s\" found in input\n",
	      getparam("times"));
  }
  if(BF) falcON_DEL_O(BF);
}
