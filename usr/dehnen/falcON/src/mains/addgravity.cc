// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// addgravity.cc                                                               |
//                                                                             |
// copyright Walter Dehnen, 2002-2003                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// adds gravity (pot & acc) to a snapshot                                      |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 1.0.0  19/10/2001  WD  created. initial version reads all snapshots.      |
// v 1.0.1  23/10/2001  WD  added option times= to select time range(s)        |
// v 1.0.2  28/08/2002  WD  adapted to various changes in nbdy                 |
// v 1.1    30/08/2002  WD  adapted this file for usage of MPI otherwise       |
// v 1.1.1  04/02/2003  WD  default falcON parameters automized                |
// v 1.2    20/03/2003  WD  gravity, action reporting                          |
// v 1.3    23/05/2003  WD automated NEMO history                              |
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
  "in=???\n            input file                      ",
  "out=???\n           output file                     ",
  "times=all\n         time range                      ",
  "eps=0.05\n          softening length                ",
  "kernel="falcON_KERNEL_TEXT
  "\n                  softening kernel                ",
  "theta="falcON_THETA_TEXT
  "\n                  tolerance parameter at M=M_tot  ",
  "Ncrit="falcON_NCRIT_TEXT
  "\n                  max # bodies in un-split cells  ",
  "VERSION=1.3"
#ifdef falcON_PROPER
  "P"
#endif
#ifdef falcON_SSE
  "S"
#endif
#ifdef falcON_INDI
  "I"
#endif
             "\n       23-may-2003 WD\n"
  "                   compiled " __DATE__ ", " __TIME__ "  ",
  NULL};
string usage = "addgravity -- adds pot & acc to a snapshot using\n"
               "falcON = Force ALgorithm with Complexity O(N)\n";
//------------------------------------------------------------------------------
void nbdy::main()
{
  nbdy::nemo_in  IN (getparam("in"));
  nbdy::nemo_out OUT;
  unsigned       NCRIT(getiparam("Ncrit"));
  nbdy::io       READ;
  const nbdy::io WRITE = nbdy::io::mxv | nbdy::io::p | nbdy::io::a;
  nbdy::real     TIME;
  nbdy::bodies   BODIES;
  nbdy::falcON   FALCON(&BODIES,
			nbdy::real(getdparam("eps")),
			nbdy::real(getdparam("theta")),
			nbdy::kern_type(getiparam("kernel")));
  while(IN.is_present(nbdy::nemo_io::snap)) {
    if(! BODIES.read_nemo_snapshot(IN,READ,&TIME,
				   nbdy::io::mxv,getparam("times")))
      continue;
    if(READ & nbdy::io::mxv) {
      if(!OUT.is_open()) OUT.open(getparam("out"));
      BODIES.flag_all_as_active();
      FALCON.grow(NCRIT);
      FALCON.approximate_gravity();
      BODIES.write_nemo_snapshot(OUT,&TIME,WRITE);
    }
  }
}
