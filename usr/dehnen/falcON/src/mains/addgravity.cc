// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// addgravity.cc                                                               |
//                                                                             |
// copyright Walter Dehnen, 2002                                               |
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
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef ALLOW_NEMO
#error You need "ALLOW_NEMO" to compile __FILE__
#endif
#include <body.h>                                     // bodies etc..           
#include <falcON.h>                                   // fakcIB                 
#include <public/nmio.h>                              // WDs C++ NEMO I/O       
#include <nemo.h>                                     // NEMO basics & main     
#include <nemomain.h>                                 // NEMO main              
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n            input file                      ",
  "out=???\n           output file                     ",
  "times=all\n         time range                      ",
  "eps=0.05\n          softening length                ",
  "kernel=1\n          softening kernel                ",
  "theta=0.6\n         tolerance parameter at M=M_tot  ",
  "Ncrit=6\n           max # bodies in un-split cells  ",
  "VERSION=1.1\n       30/August/2002 WD\n"
  "                   compiled  " __DATE__ ", " __TIME__ "                    ",
  NULL};
string usage = "addgravity -- adds pot & acc to a snapshot using\n"
               "falcON = Force ALgorithm with Complexity O(N)\n";
//------------------------------------------------------------------------------
void nemo::main()
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
  IN.read_history();
  while(IN.is_present(nbdy::nemo_io::snap)) {
    if(! BODIES.read_nemo_snapshot(IN,READ,&TIME,
				   nbdy::io::mxv,getparam("times")))
      continue;
    if(READ & nbdy::io::mxv) {
      if(!OUT.is_open()) {
	OUT.open(getparam("out"));
	OUT.write_history();
      }
      BODIES.flag_all_as_sink();
      FALCON.grow(NCRIT);
      FALCON.approximate_gravity();
      BODIES.write_nemo_snapshot(OUT,&TIME,WRITE);
    }
  }
}
