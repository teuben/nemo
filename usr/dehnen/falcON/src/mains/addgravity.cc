// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// addgravity.cc                                                               |
//                                                                             |
// copyright Walter Dehnen, 2002-2003                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// adds gravity (pot & acc) to a snapshot                                      |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 1.0.0  19/10/2001  WD created. initial version reads all snapshots.       |
// v 1.0.1  23/10/2001  WD added option times= to select time range(s)         |
// v 1.0.2  28/08/2002  WD adapted to various changes in nbdy                  |
// v 1.1    30/08/2002  WD adapted this file for usage of MPI otherwise        |
// v 1.1.1  04/02/2003  WD default falcON parameters automized                 |
// v 1.2    20/03/2003  WD gravity, action reporting                           |
// v 1.3    23/05/2003  WD automated NEMO history                              |
// v 1.3.1  23/05/2003  WD automated version & compile information             |
// v 1.4    06/05/2004  WD new body.h; write=read + gravity                    |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "1.4"
#define falcON_VERSION_D "06-may-2004 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO
#error You need NEMO to compile "addgravity.cc"
#endif
#include <body.h>                                     // bodies etc..           
#include <falcON.h>                                   // falcON                 
#include <public/nmio.h>                              // WDs C++ NEMO I/O       
#include <main.h>                                     // main & NEMO stuff      
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n           input file                                         ",
  "out=???\n          output file                                        ",
  "times=all\n        time range                                         ",
  "eps=0.05\n         softening length                                   ",
  "kernel="falcON_KERNEL_TEXT
  "\n                 softening kernel                                   ",
  "theta="falcON_THETA_TEXT
  "\n                 tolerance parameter at M=M_tot                     ",
  "Ncrit="falcON_NCRIT_TEXT
  "\n                 max # bodies in un-split cells                     ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string usage = "addgravity -- adds pot & acc to a snapshot using\n"
               "falcON = Force ALgorithm with Complexity O(N)\n";
//------------------------------------------------------------------------------
void nbdy::main()
{
  nemo_in  IN (getparam("in"));
  nemo_out OUT;
  unsigned NCRIT(getiparam("Ncrit"));
  io       READ;
  double   TIME;
  bodies   BODIES(0,io::mx|io::ap|io::f);
  falcON   FALCON(&BODIES,
		  getdparam("eps"),
		  getdparam("theta"),
		  kern_type(getiparam("kernel")));
  while(IN.is_present(nemo_io::snap)) {
    if(! BODIES.read_nemo_snapshot(IN,READ,&TIME,io::all,getparam("times"),0))
      continue;
    if(READ & io::mx) {
      if(!OUT.is_open()) OUT.open(getparam("out"));
      BODIES.flag_all_as_active();
      FALCON.grow(NCRIT);
      FALCON.approximate_gravity();
      BODIES.write_nemo_snapshot(OUT,&TIME,READ|io::ap);
    } else
      warning("bodies' mx not read at time %f: cannot compute gravity",TIME);
  }
}
