// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// symmetrize.cc                                                               |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
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
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "0.6"
#define falcON_VERSION_D "10-may-2004 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile "src/mains/symmetrize.cc"
#endif
#define falcON_RepAction 0                         // no action reporting       
#include <iostream>                                // C++ I/O                   
#include <fstream>                                 // C++ file I/O              
#include <iomanip>                                 // C++ I/O formatting        
#include <body.h>                                  // the bodies                
#include <public/nmio.h>                           // my NEMO I/O               
#include <public/ionl.h>                           // my I/O utilities          
#include <main.h>                                  // main & NEMO stuff         
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n        input file                                         ",
  "out=???\n       output file                                        ",
  "times=all\n     times to process                                   ",
  "use=2\n         use every ith particle                             ",
  "copy=1\n        make 2^h symmetized copies, h in [0,1,2]           ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string
usage = "symmetrize -- symmetrizes snapshots\n"
        "every ith particle is taken and replaced by its symmetric\n"
        "wrt origin and z=0 copies. can also be used to reduce N\n";
//------------------------------------------------------------------------------
namespace { using namespace nbdy;
  template<int IO, int COPIES, int C=0> struct copy_datum {
    static void c(bodies const&Bf, bodies const&Bt, int f, int t) {
      Bt. template datum_io<IO>(t+C) = Bf. template const_datum_io<IO>(f);
      copy_datum<IO,COPIES,C+1>::c(Bf,Bt,f,t);
    } };
  template<int IO, int C> struct copy_datum<IO,C,C> {
    static void c(bodies const&, bodies const&, int, int) {} };
  //----------------------------------------------------------------------------
  template<int COPIES, int BIT=0> struct copy_data {
    static void c(bodies const&Bf, bodies const&Bt,
		  io copy, int f, int t) {
      if(copy & 1<<BIT) copy_datum<1<<BIT,COPIES>::c(Bf,Bt,f,t);
      copy_data<COPIES,BIT+1>::c(Bf,Bt,copy,f,t);
    } };
  template<int COPIES> struct copy_data<COPIES,IO_NQUANT> {
    static void c(bodies const&, bodies const&, io, int, int) {} };
}
//------------------------------------------------------------------------------
void nbdy::main()
{
  nemo_in      IN  (getparam("in"));
  nemo_out     OUT;
  double       TIME;
  int          USE (getiparam("use"));
  if(USE  < 1) USE  = 1;
  int          COPY(getiparam("copy"));
  if(COPY < 0) COPY = 0;
  if(COPY > 2) COPY = 2;
  bodies       BIN;
  io           READ, WANT = bodies::NIOBits;
  real         MFAC = (USE == 1<<COPY)? 1. : double(USE) / double(1<<COPY);
  while(IN.is_present(nemo_io::snap)) {
    if(! BIN.read_nemo_snapshot(IN,READ,&TIME,WANT,getparam("times"),false))
      continue;
    if(USE == 1 && COPY == 0)                      // IF output==input          
      BIN.write_nemo_snapshot(OUT,&TIME,READ);
    else {                                         // ELSE (output != input)    
      int NOUT = BIN.N_bodies()/USE;
      if(BIN.N_bodies() % USE) ++NOUT;
      NOUT *= 1<<COPY;
      bodies BOUT(NOUT,READ);
      if(COPY == 0) {
	for(register int i=0,j=0; i<BIN.N_bodies(); i+=USE, ++j) {
	  copy_data<1>::c(BIN,BOUT,READ,i,j);
	}
      } else if(COPY==1) {
	if(READ & io::a) {
	  for(register int i=0,j=0; i<BIN.N_bodies(); i+=USE, j+=2) {
	    copy_data<2>::c(BIN,BOUT,READ,i,j);
	    BOUT.pos(j+1) *=-one;                  // -x,-y,-z                  
	    BOUT.vel(j+1) *=-one;                  // -u,-v,-w                  
	    BOUT.acc(j+1) *=-one;                  // -ax,-ay,-az               
	  }
	} else {
	  for(register int i=0,j=0; i<BIN.N_bodies(); i+=USE, j+=2) {
	    copy_data<2>::c(BIN,BOUT,READ,i,j);
	    BOUT.pos(j+1) *=-one;                  // -x,-y,-z                  
	    BOUT.vel(j+1) *=-one;                  // -u,-v,-w                  
	  }
	}
      } else {
	if(READ & io::a) {
	  for(register int i=0,j=0; i<BIN.N_bodies(); i+=USE, j+=4) {
	    copy_data<4>::c(BIN,BOUT,READ,i,j);
	    BOUT.pos(j+1)    *=-one;               // -x,-y,-z                  
	    BOUT.vel(j+1)    *=-one;               // -u,-v,-w                  
	    BOUT.pos(j+2)[2] *=-one;               //  x, y,-z                  
	    BOUT.vel(j+2)[2] *=-one;               //  u, v,-w                  
	    BOUT.pos(j+3)[0] *=-one;               // -x                        
	    BOUT.vel(j+3)[0] *=-one;               // -u                        
	    BOUT.pos(j+3)[1] *=-one;               // -y                        
	    BOUT.vel(j+3)[1] *=-one;               // -v                        
	  }
	} else {
	  for(register int i=0,j=0; i<BIN.N_bodies(); i+=USE, j+=4) {
	    copy_data<4>::c(BIN,BOUT,READ,i,j);
	    BOUT.pos(j+1)    *=-one;               // -x,-y,-z                  
	    BOUT.vel(j+1)    *=-one;               // -u,-v,-w                  
	    BOUT.acc(j+1)    *=-one;               // -ax,-ay,-az               
	    BOUT.pos(j+2)[2] *=-one;               //  x, y,-z                  
	    BOUT.vel(j+2)[2] *=-one;               //  u, v,-w                  
	    BOUT.acc(j+2)[2] *=-one;               // ax,ay,-az                 
	    BOUT.pos(j+3)[0] *=-one;               // -x                        
	    BOUT.vel(j+3)[0] *=-one;               // -u                        
	    BOUT.acc(j+3)[0] *=-one;               // -ax                       
	    BOUT.pos(j+3)[1] *=-one;               // -y                        
	    BOUT.vel(j+3)[1] *=-one;               // -v                        
	    BOUT.acc(j+3)[1] *=-one;               // -ay                       
	  }
	}
      }
      if(!OUT.is_open()) OUT.open(getparam("out"));
      if(READ & io::m && USE != 1<<COPY)
	LoopBodies(bodies,&BOUT,B) B.mass() *= MFAC;
      BOUT.write_nemo_snapshot(OUT,&TIME,READ);
    }
  }
}
