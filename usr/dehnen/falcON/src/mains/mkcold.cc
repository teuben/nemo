// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// mkcold.cc                                                                   |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2002                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 0.0     06/12/2002  WD created                                            |
// v 0.1     09/12/2002  WD added option for figure rotation                   |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef ALLOW_NEMO                                 // this is a NEMO program    
#  error You need "ALLOW_NEMO" to compile __FILE__
#endif
#include <body.h>                                  // N-body bodies             
#if NDIM != 3                                      // this is a NEMO program    
#  error You need NDIM == 3 to compile __FILE__
#endif
#include <public/nmio.h>                           // my NEMO I/O               
#include <proper/rand.h>                           // random deviates           
#include <nemomain.h>                              // NEMO main                 
extern "C" {
#include <nemo.h>                                  // NEMO basics & main        
#include <mathfns.h>                               // NEMO random numbers       
}

using namespace nbdy;
//------------------------------------------------------------------------------
string defv[] = {
  "out=???\n           output file                              ",
  "nbody=???\n         number of bodies                         ",
  "a=1\n               semi major axis                          ",
  "b=1\n               semi intermediate axis                   ",
  "c=1\n               semi minor axis                          ",
  "Mtot=1\n            total mass                               ",
  "seed=0\n            seed for the randum number generator     ",
  "time=0\n            simulation time of snapshot              ",
  "sym=t\n             enforce 8-fold symmetry?                 ",
  "sigma=0\n           velocity dispersion (Maxwellian)         ",
  "omega=0\n           rigid pattern speed about z-axis         ",
  "VERSION=0.0\n       06/December/2002 WD\n"
  "                   compiled  " __DATE__ ", " __TIME__ "          ",
  NULL};
string usage = "mkcold -- cold triaxial Gaussian initial conditions\n";
//------------------------------------------------------------------------------
namespace nbdy {
  class NemoRNG : public RandomNumberGenerator {
  public:
    explicit NemoRNG(const char* seed)
    { init_xrandom(const_cast<char*>(seed)); }
    double RandomDouble() { return xrandom(0.,1.); }
  };
}
//------------------------------------------------------------------------------
void nemo::main()
{
  bool symmetric = getbparam("sym");
  unsigned N     = getiparam("nbody");
  if(symmetric & N%8) {
    warning(" %d not dividable by 8; using nbody=%d instead",N,8*(N/8+1));
    N = 8 * (N/8+1);
  }
  bodies BB(N);
  double mass  = getdparam("Mtot") / double(N);
  NemoRNG RNG (getparam("seed"));
  double sigma = getdparam("sigma");
  double omega = getdparam("omega");
  Gaussian Gx (&RNG,&RNG,getdparam("a"));
  Gaussian Gy (&RNG,&RNG,getdparam("b"));
  Gaussian Gz (&RNG,&RNG,getdparam("c"));
  Gaussian Gv (&RNG,&RNG,sig);
  register nbdy::real x[3], v[3] = {0.,0.,0.};
  for(register int n=0; n!=N; ) {
    x[0] = Gx();
    x[1] = Gy();
    x[2] = Gz();
    if(sigma) {
      v[0] = Gv();
      v[1] = Gv();
      v[2] = Gv();
    }
    if(symmetric)
      for(register int i=0; i!=8; ++i, ++n) {
	BB.mas(n)    = mass;
	BB.pos(n)[0] = i&1 ? x[0] : -x[0];
	BB.pos(n)[1] = i&2 ? x[1] : -x[1];
	BB.pos(n)[2] = i&4 ? x[2] : -x[2];
 	BB.vel(n)[0] = (i&1 ? v[0] : -v[0]) - omega * (i&2 ? x[1] : -x[1]);
	BB.vel(n)[1] = (i&2 ? v[1] : -v[1]) + omega * (i&1 ? x[0] : -x[0]);
	BB.vel(n)[2] = i&4 ? v[2] : -v[2];
	if(omega) {
	  BB.vel(n)[0] += rot[0];
	  BB.vel(n)[1] += rot[1];
	}
     }
    else {
      BB.mas(n)    = mass;
      BB.pos(n)    = x;
      BB.vel(n)[0] = v[0] - omega * x[1];
      BB.vel(n)[1] = v[1] + omega * x[0];
      BB.vel(n)[2] = v[2];
      n++;
    }
  }
  nbdy::nemo_out out(getparam("out"));
  out.write_history();
  nbdy::real time = getdparam("time");
  BB.write_nemo_snapshot(out,&time,nbdy::io::mxv);
}
