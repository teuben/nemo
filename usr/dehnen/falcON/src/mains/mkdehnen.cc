// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// mkdehnen.cc                                                                 |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2004                                               |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// creates N-body initial conditions from a Dehnen model                       |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 1.0   26/02/2004  WD created                                              |
// v 1.1   27/02/2004  WD added some consistency checks & error messages       |
// v 1.2   01/03/2004  WD corrected asymptotic of f at Q=0 with OM anisotropy  |
// v 1.3   01/03/2004  WD added option f_pos; re-arranged order of options     |
// v 2.0   03/03/2004  WD using nsam.h; fixed bug with scaling & Rad[]         |
// v 2.1   05/03/2004  WD re-written gama; improved scaling support            |
// v 2.2   26/03/2004  WD bug fixed in this file (f(E=0)=nan -> sampling error)|
// v 2.2.1 02/04/2004  WD put class DehnenModelSampler into gama.h             |
// v 2.3   04/05/2004  WD happy icc 8.0; new body.h; new make.proper           |
// v 2.4   12/05/2004  WD made partly PUBLIC                                   |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "2.4"
#define falcON_VERSION_D "12-may-2004 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#error You need NEMO to compile mkdehnen
#endif
#define falcON_RepAction 0                         // no action reporting       
//-----------------------------------------------------------------------------+
#include <body.h>                                  // N-body bodies             
#if falcON_NDIM != 3                               // this is a NEMO program    
#  error You need falcON_NDIM == 3 to compile "src/mains/mkdehnen.cc"
#endif
#include <public/nmio.h>                           // my NEMO file I/O          
#include <public/gama.h>                           // my gamma models           
#include <main.h>                                  // main & NEMO stuff         
////////////////////////////////////////////////////////////////////////////////
string defv[] = {
  "out=???\n          output file                                        ",
  "nbody=???\n        number of bodies                                   ",
  "gamma=???\n        inner power-law slope of density, gamma in[0,2.5[  ",
  "mass=1\n           total mass                                         ",
  "r_s=1\n            scale radius                                       ",
#ifdef falcON_PROPER
  "r_a=\n             Ossipkov-Merritt anisotropy radius                 ",
#endif
  "seed=0\n           seed for the randum number generator               ",
  "q-ran=f\n          use quasi- instead of pseudo-random numbers        ",
  "time=0\n           simulation time of snapshot                        ",
  "f_pos=0.5\n        fraction of bodies with positive sense of rotation ",
  "rmax=\n            if given, only emit bodies with r <= rmax          ",
#ifdef falcON_PROPER
  "Rp=\n              for mass adaption: list of R in increasing order   ",
  "fac=1.2\n          for mass adaption: factor between mass bins        ",
  "peri=f\n           for mass adaption: R_peri(E,L) rather than R_c(E)  ",
  "epar=\n            if given, set eps_i = epar * sqrt(m_i/M_tot)       ",
#endif
  "WD_units=f\n       input:  kpc, solar masses\n"
  "                   output: kpc, kpc/Gyr, G=1 (-> mass unit)           ",
  falcON_DEFV, NULL };
////////////////////////////////////////////////////////////////////////////////
string usage = "mkdehnen -- construct a Dehnen (1993, MNRAS, 265, 250) model"
#ifdef falcON_PROPER
               "\n            possibly with Osipkov-Merritt anisotropy"
               "\n            and mass adaption (proprietary version only)"
#endif
;
////////////////////////////////////////////////////////////////////////////////
void nbdy::main()
{
  const double mf = 2.2228847e5;                   // M_sun in WD_units         
  //----------------------------------------------------------------------------
  // 1. set some parameters                                                     
  //----------------------------------------------------------------------------
  if(hasvalue("rmax") && getdparam("rmax")<=0.) ::error("rmax <= 0");
#ifdef falcON_PROPER
  if(hasvalue("r_a")  && getdparam("r_a") <=0.) ::error("r_a <= 0");
#endif
  const bool   WD  (getbparam("WD_units"));        // using WD_units?           
  const Random Ran (getparam("seed"),6);           // random-number-generators  
  const io data = 
#ifdef falcON_PROPER
    hasvalue("epar")? io::e  | io::mxv :
#endif
    io::mxv;
#ifdef falcON_PROPER
  const int    nbmax(100);
  double       Rad[nbmax];
  int          nb=0;
  if(hasvalue("Rp")) {
    nb=nemoinpd(getparam("Rp"),Rad,nbmax);
    if(nb+1>nbmax)
      ::error("exceeding expected number of radii in mass adaption");
    Rad[nb] = Rad[nb-1] * 1.e20;
    nb++;
  }
#endif
  //----------------------------------------------------------------------------
  // 2. create initial conditions from a Dehnen model using mass adaption       
  //----------------------------------------------------------------------------
  DehnenModelSampler GS(getdparam("gamma"),
			getdparam("r_s"),
			WD? getdparam("mass")/mf : getdparam("mass"),
#ifdef falcON_PROPER
			getdparam_z("r_a"),
#endif
			getdparam_z("rmax"),
			9999,1.e-8
#ifdef falcON_PROPER
		       ,Rad,nb,
			getdparam("fac"),
			getbparam("peri")
#endif
			);
  bodies BB(getiparam("nbody"), data);
#ifdef falcON_PROPER
#endif
  GS.sample(BB,getbparam("q-ran"),Ran,getdparam("f_pos")
#ifdef falcON_PROPER
	    ,getdparam("epar")
#endif
	    );
  //----------------------------------------------------------------------------
  // 3. output of snapshot                                                      
  //----------------------------------------------------------------------------
  nemo_out out(getparam("out"));
  double time = getdparam("time");
  BB.write_nemo_snapshot(out,&time,data);
}
