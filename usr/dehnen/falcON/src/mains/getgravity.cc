// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// getgravity.cc                                                               |
//                                                                             |
// copyright Walter Dehnen, 2002                                               |
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
  "srce=???\n          input file: sources [m,x]       ",
  "sink=???\n          input file: sinks   [x]         ",
  "out=???\n           output file         [x,a,p]     ",
  "times=all\n         time range (for srce only)      ",
  "eps=0.05\n          softening length                ",
  "kernel=1\n          softening kernel                ",
  "theta=0.6\n         tolerance parameter at M=M_tot  ",
  "Ncrit=6\n           max # bodies in un-split cells  ",
  "VERSION=0.0\n       23/November/2002 WD             ",
  NULL};
string usage = "getgravity -- computes gravity at test positions using\n"
               "falcON = Force ALgorithm with Complexity O(N)\n";
//------------------------------------------------------------------------------
void nemo::main()
{
  nbdy::nemo_in  IN[2];
  IN[0].open(getparam("srce"));
  IN[0].read_history();

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
  IN[0].read_history();
  while(IN[0].is_present(nbdy::nemo_io::snap)) {
    IN[1].open(getparam("sink"));
    IN[1].read_history();
    if(BODIES.read_nemo_snapshots(IN,2,READ,NN,&TIME,TOREAD,getparam("times"))
       && READ[0] & TOREAD[0]
       && READ[1] & TOREAD[1]) {
      register nbdy::real mass = BODIES.mas(0);
      for(register nbdy::uint i=0; i!=NN[0]; ++i) {
	BODIES.flg(i).un_set(nbdy::flag::SINK);
	if(BODIES.mas(i) < mass) mass = BODIES.mas(i);
      }
      mass *= 1.e-10/nbdy::real(NN[1]);
      for(register nbdy::uint i=NN[0]; i!=BODIES.N_bodies(); ++i) {
	BODIES.flg(i).add(nbdy::flag::SINK);
	BODIES.mas(i) = mass;
      }
      FALCON.grow(NCRIT);
      FALCON.approximate_gravity();
      if(!OUT.is_open()) {
	OUT.open(getparam("out"));
	OUT.write_history();
      }
      BODIES.write_nemo_snapshot(OUT,&TIME,WRITE,NN[1],NN[0]);
    }
  }
}
