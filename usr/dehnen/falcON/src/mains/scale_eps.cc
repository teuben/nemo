// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// scale_eps.cc                                                                |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2003-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 1.0    23/10/2003  WD created                                             |
// v 1.1    07/22/2003  WD marginal changes due to changes in body.h           |
// v 1.2    30/04/2004  WD marginal changes due to changes in body.h           |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "1.2"
#define falcON_VERSION_D "30-apr-2004 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#error You need NEMO to compile "src/mains/scale_eps.cc"
#endif
#define falcON_RepAction 0                         // no action reporting       
#include <iostream>                                // C++ I/O                   
#include <fstream>                                 // C++ file I/O              
#include <cmath>                                   // standard math             
#include <body.h>                                  // the bodies                
#include <public/nmio.h>                           // my NEMO I/O               
#include <main.h>                                  // main & NEMO stuff         
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n      input file                            ",
  "out=???\n     output file                           ",
  "fac=???\n     scale factor for eps_i                ",
  "write=mxve\n  what to write to out                  ",
  "times=all\n   times of snapshots to scale           ",
  falcON_DEFV, NULL};
//------------------------------------------------------------------------------
string
usage = "scale_eps -- scales individual softening lengths";
//------------------------------------------------------------------------------
void nbdy::main()
{
  nemo_in   IN (getparam("in"));
  nemo_out  OUT;
  io        READ, NEED = io(getparam("write")) | io::e;
  double    TIME;
  bodies    BODIES;
  const real FAC(getdparam("fac"));
  while(IN.is_present(nbdy::nemo_io::snap)) {
    if(! BODIES.read_nemo_snapshot(IN,READ,&TIME,NEED,getparam("times")))
      continue;
    if(! READ.contains(NEED)) error("insufficient input data");
    LoopBodies(bodies,(&BODIES),Bi) Bi.eps() *= FAC;
    if(!OUT.is_open()) OUT.open(getparam("out"));
    BODIES.write_nemo_snapshot(OUT,&TIME,READ);
  }
}
