// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// getgravity.cc                                                               |
//                                                                             |
// copyright Walter Dehnen, 2002-2003                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// computes gravity (pot & acc) at sink position and due to other sources.     |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 0.0    23/11/2002  WD  created.                                           |
// v 0.1    04/02/2003  WD  default falcON parameters automized                |
// v 0.2    20/03/2003  WD  gravity, action reporting                          |
// v 0.3    23/05/2003  WD automated NEMO history                              |
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO
#error You need NEMO to compile "src/mains/getgravity.cc"
#endif
#include <body.h>                                  // bodies etc..              
#include <falcON.h>                                // fakcIB                    
#include <public/nmio.h>                           // WDs C++ NEMO I/O          
#include <main.h>                                  // main & NEMO stuff         
//------------------------------------------------------------------------------
string defv[] = {
  "srce=???\n          input file: sources [m,x]       ",
  "sink=???\n          input file: sinks   [x]         ",
  "out=???\n           output file         [x,a,p]     ",
  "times=all\n         time range (for srce only)      ",
  "eps=0.05\n          softening length                ",
  "kernel="falcON_KERNEL_TEXT
  "\n                  softening kernel                ",
  "theta="falcON_THETA_TEXT
  "\n                  tolerance parameter at M=M_tot  ",
  "Ncrit="falcON_NCRIT_TEXT
  "\n                  max # bodies in un-split cells  ",
  "VERSION=0.3"
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
string usage = "getgravity -- computes gravity at test positions using\n"
               "falcON = Force ALgorithm with Complexity O(N)\n";
//------------------------------------------------------------------------------
void nbdy::main()
{
  nbdy::nemo_in  IN[2];
  IN[0].open(getparam("srce"));
  nbdy::nemo_out OUT;
  unsigned       NCRIT(getiparam("Ncrit"));
  nbdy::io       READ[2],TOREAD[2] = {nbdy::io::mx, nbdy::io::x};
  const nbdy::io WRITE = nbdy::io::x | nbdy::io::p | nbdy::io::a;
  nbdy::real     TIME;
  nbdy::uint     NN[2];
  nbdy::sbodies  BODIES;
  nbdy::falcON   FALCON(&BODIES,
			nbdy::real(getdparam("eps")),
			nbdy::real(getdparam("theta")),
			nbdy::kern_type(getiparam("kernel")));
  while(IN[0].is_present(nbdy::nemo_io::snap)) {
    IN[1].open(getparam("sink"));
    if(BODIES.read_nemo_snapshots(IN,2,READ,NN,&TIME,TOREAD,getparam("times"))
       && READ[0] & TOREAD[0]
       && READ[1] & TOREAD[1]) {
      register nbdy::real mass = BODIES.mass(0);
      for(register nbdy::uint i=0; i!=NN[0]; ++i) {
	BODIES.flg(i).un_set(nbdy::flag::ACTIVE);
	if(BODIES.mass(i) < mass) mass = BODIES.mass(i);
      }
      mass *= 1.e-10/nbdy::real(NN[1]);
      for(register nbdy::uint i=NN[0]; i!=BODIES.N_bodies(); ++i) {
	BODIES.flg(i).add(nbdy::flag::ACTIVE);
	BODIES.mass(i) = mass;
      }
      FALCON.grow(NCRIT);
      FALCON.approximate_gravity();
      if(!OUT.is_open()) OUT.open(getparam("out"));
      BODIES.write_nemo_snapshot(OUT,&TIME,WRITE,NN[1],NN[0]);
    }
  }
}
