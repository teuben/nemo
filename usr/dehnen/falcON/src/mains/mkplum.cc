// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// mkplum.cc                                                                   |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 0.0   28/10/2003  WD created                                              |
// v 0.1   27/01/2004  WD velocity sampling always by pseudo RNG               |
// v 1.0   08/03/2004  WD changed to using nsam.h; added Ossipkov-Merritt DF   |
// v 1.0.1 02/04/2004  WD moved stuff into plum.h                              |
// v 1.1   12/05/2004  WD made PUBLIC; changed 'a' -> 'r_s'                    |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "1.0.1"
#define falcON_VERSION_D "02-apr-2004 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile mkplum
#endif
#define falcON_RepAction 0                         // no action reporting       
#include <public/plum.h>                           // Plummer model & sampler   
#if falcON_NDIM != 3                               // this is a NEMO program    
#  error You need falcON_NDIM == 3 to compile mkplum
#endif
#include <public/nmio.h>                           // my NEMO I/O               
#include <public/rand.h>                           // random deviates           
#include <main.h>                                  // main & NEMO stuff         
////////////////////////////////////////////////////////////////////////////////
string defv[] = {
  "out=???\n          output file                                        ",
  "nbody=???\n        number of bodies                                   ",
  "r_s=1\n            scale radius                                       ",
  "mass=1\n           total mass of Plummer model                        ",
#ifdef falcON_PROPER
  "r_a=\n             Ossipkov-Merritt anisotropy radius                 ",
#endif
  "seed=0\n           seed for the randum number generator               ",
  "q-ran=f\n          use quasi- instead of pseudo-random numbers        ",
  "time=0\n           simulation time of snapshot                        ",
  "f_pos=0.5\n        fraction of bodies with positive sense of rotation ",
  "rmax=\n            if given, only emit bodies with r <= rmax          ",
  "WD_units=f\n       input:  kpc, M_Sun\n"
  "                   output: kpc, kpc/Gyr, G=1 (-> mass unit)           ",
  falcON_DEFV, NULL };
string usage = "mkplum -- initial conditions from a Plummer model";
////////////////////////////////////////////////////////////////////////////////
void nbdy::main()
{
  const double mf= 2.2228847e5;                    // mass unit for WD_units    
  const bool   WD  (getbparam("WD_units"));        // using WD_units?           
  const Random Ran (getparam("seed"),6);           // random-number-generators  
  if(hasvalue("rmax") && getdparam("rmax")<=0.) ::error("rmax <= 0");
#ifdef falcON_PROPER
  if(hasvalue("r_a")  && getdparam("r_a") <=0.) ::error("r_a <= 0");
#endif
  bodies BB(getiparam("nbody"));
  PlummerModelSampler PS(getdparam("r_s"),
			 WD? getdparam("mass")/mf : getdparam("mass"),
#ifdef falcON_PROPER
			 getdparam_z("r_a"),
#endif
			 getdparam_z("rmax"));
  PS.sample(BB,getbparam("q-ran"),Ran,getdparam("f_pos"));
  nemo_out out(getparam("out"));
  double time = getdparam("time");
  BB.write_nemo_snapshot(out,&time,io::mxv);
}
