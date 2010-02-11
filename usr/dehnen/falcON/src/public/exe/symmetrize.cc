// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// symmetrize.cc                                                               |
//                                                                             |
// Copyright (C) 2002-2006, 2008 Walter Dehnen                                 |
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
// v 0.0   18/11/2002  WD created                                              |
// v 0.1   19/11/2002  WD first bug ...                                        |
// v 0.2   19/11/2002  WD action reporting                                     |
// v 0.3   23/05/2003  WD automated NEMO history                               |
// v 0.4   23/10/2003  WD automated version, compiler in version etc           |
// v 0.5   01/05/2004  WD happy icc 8.0; new body.h & changes in this file     |
// v 0.6   01/05/2004  WD removed bug (Seg fault) made at v 0.5; acceleration  |
// v 1.0   19/05/2004  WD removed bug (v0.6); allowed for missing x,v, or a    |
// v 1.1   20/05/2005  WD several minor updates                                |
// v 2.0   14/06/2005  WD new falcON                                           |
// v 2.1   28/06/2005  WD deBUGged                                             |
// v 2.1.1 31/03/2006  WD BD_NQUANT -> BodyData::NQUANT                        |
// v 2.1.2 20/02/2008  WD change in body.h (removed old-style constructors)    |
// v 2.1.3 10/09/2008  WD happy gcc 4.3.1                                      |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "2.1.3"
#define falcON_VERSION_D "10-sep-2008 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile "symmetrize"
#endif
#define falcON_RepAction 0                         // no action reporting       
#define falcON_NOT_USING_PROPER                    // not using PROPRIETARY code
//-----------------------------------------------------------------------------+
#include <iostream>                                // C++ I/O                   
#include <fstream>                                 // C++ file I/O              
#include <iomanip>                                 // C++ I/O formatting        
#include <main.h>                                  // main & NEMO stuff         
//------------------------------------------------------------------------------
const char*defv[] = {
  "in=???\n        input file                                         ",
  "out=???\n       output file                                        ",
  "times=all\n     times to process                                   ",
  "use=2\n         use every ith particle                             ",
  "copy=1\n        make 2^h symmetized copies, h in [0,1,2]           ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
const char*
usage = "symmetrize -- symmetrizes snapshots\n"
        "every ith particle is taken and replaced by its symmetric\n"
        "wrt origin and z=0 copies. can also be used to reduce N\n";
//------------------------------------------------------------------------------
namespace { using namespace falcON;
  // class that determines how a datum's copy No C looks like
  // default: identical
  template<int C, typename T> struct copy_type {
    static void c(T&to, T const&from) { to = from; }
  };
  // vectors: 2nd copy (C=1): reflected w.r.t. origin
  template<> struct copy_type<1,vect> {
    static void c(vect&to, vect const&from) { 
      to[0]=-from[0];
      to[1]=-from[1];
      to[2]=-from[2];
    }
  };
  // vectors: 3rd copy (C=2): reflected w.r.t. equatorial plane
  template<> struct copy_type<2,vect> {
    static void c(vect&to, vect const&from) { 
      to[0]= from[0];
      to[1]= from[1];
      to[2]=-from[2];
    }
  };
  // vectors: 4th copy (C=3): reflected w.r.t. z-axis
  template<> struct copy_type<3,vect> {
    static void c(vect&to, vect const&from) {
      to[0]=-from[0];
      to[1]=-from[1];
      to[2]= from[2];
    }
  };
  //----------------------------------------------------------------------------
  template<int C> struct copy {
    template<typename T> static void c(T&to, T const&from) {
      copy_type<C,T>::c(to,from);
    }
  };
  //----------------------------------------------------------------------------
  template<int BIT, int COPIES, int C=0> struct copy_datum {
    static void c(body const&from, body&to) {
      copy<C>::c(to. template datum<BIT>(), const_datum<BIT>(from));
      copy_datum<BIT, COPIES, C+1>::c(from,++to);
    }
  };
  template<int BIT, int C> struct copy_datum<BIT,C,C> {
    static void c(body const&, body&) {}
  };
  //----------------------------------------------------------------------------
  template<int COPIES, int BIT=0> struct copy_data {
    static void c(body const&from, body const&to, fieldset copy) {
      if(copy.contain(fieldbit(BIT))) {
	body f(from), t(to);
	copy_datum<BIT,COPIES>::c(f,t);
      }
      copy_data<COPIES,BIT+1>::c(from,to,copy);
    }
  };
  template<int COPIES> struct copy_data<COPIES,BodyData::NQUANT> {
    static void c(body const&, body const&, fieldset) {}
  };
} // namespace {
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  nemo_in  in(getparam("in"));
  nemo_out out;
  unsigned use (max(1u,getuparam("use")));
  unsigned copy(min(2u,getuparam("copy")));
  snapshot shin;
  fieldset read,want(fieldset::nemoin);
  real     Mfac(use==(1u<<copy) ? one : double(use)/double(1<<copy));
  while(in.has_snapshot()) {
    if(! shin.read_nemo(in,read,want,getparam("times"),0)) continue;
    if(use==1u && copy == 0u) {                    // IF output==input          
      if(!out.is_open()) out.open(getparam("out"));
      shin.write_nemo(out,read);
    } else {                                       // ELSE (output != input)    
      unsigned Nout[bodytype::NUM];
      for(bodytype t; t; ++t) {
	Nout[t] = shin.N_bodies(t)/use;
	if(shin.N_bodies(t) % use) ++(Nout[t]);
	Nout[t] *= 1<<copy;
      }
      snapshot shou(shin.time(),Nout,read);
      if       (copy == 0u)
	for(body
	      from=shin.begin_all_bodies(),
	      to  =shou.begin_all_bodies(); from; from+=use, ++to)
	  copy_data<1>::c(from,to,read);
      else if(copy == 1u)
	for(body
	      from=shin.begin_all_bodies(),
	      to  =shou.begin_all_bodies(); from; from+=use, to+=2)
	  copy_data<2>::c(from,to,read);
      else
	for(body
	      from=shin.begin_all_bodies(),
	      to  =shou.begin_all_bodies(); from; from+=use, to+=4)
	  copy_data<4>::c(from,to,read);
      if(!out.is_open()) out.open(getparam("out"));
      if(read.contain(fieldbit::m) && Mfac != one)
	LoopAllBodies(&shou,b) b.mass() *= Mfac;
      shou.write_nemo(out,read);
    }
  }
}
