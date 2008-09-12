// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// mkking.cc                                                                   |
//                                                                             |
// Copyright (C) 2000-2005, 2008 Walter Dehnen                                 |
//                                                                             |
// This program is free software; you can redistribute it and/or modify        |
// it under the terms of the GNU General Public License as published by        |
// the Free Software Foundation; either version 2 of the License, or (at       |
// your option) any later version.                                             |
//                                                                             |
// This program is distributed in the hope that it will be useful, but         |
// WITHOUT ANY WARRANTY; without even the implied warranty of                  |
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU           |
// General Public License for more details.                                    |
//                                                                             |
// You should have received a copy of the GNU General Public License           |
// along with this program; if not, write to the Free Software                 |
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                   |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// creates N-body initial conditions from a King model                         |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 1.0   11/09/2001  WD created for nemo 3.0.4                               |
// v 1.0.1 17/09/2001  WD changed to init_xrandom (nemo 3.0.8)                 |
// v 1.0.2 19/09/2001  WD bug fixed; centrate centre of mass (option zerocm)   |
// v 1.0.3 21/09/2001  WD changes in body & io                                 |
// v 1.0.4 26/09/2001  WD option WD_units; more bugs fixed                     |
// v 1.0.5 27/02/2002  WD changed core_radius -> r_c, tidal_radius -> r_t      |
// v 1.0.6 02/04/2002  WD improved outputs to log                              |
// v 1.0.7 07/06/2002  WD added output of history                              |
// v 1.1   30/08/2002  WD adapted this file for usage of MPI otherwise         |
// v 1.1.1 05/12/2002  WD nemo::main, new file layout ...                      |
// v 1.2   20/03/2003  WD action reporting                                     |
// v 1.2.1 20/03/2003  WD report core density with outputs                     |
// v 1.2.3 09/04/2003  WD half-mass radius, table output                       |
// v 1.3   23/05/2003  WD automated NEMO history                               |
// v 1.3.1 17/06/2003  WD added option eps                                     |
// v 1.4   28/10/2003  WD automatic version; quasi-random numbers              |
// v 1.4.1 30/04/2004  WD abandoned zerocm; new body.h; happy icc 8.0          |
// v 1.4.2 18/05/2004  WD rearranged options                                   |
// v 1.4.3 20/05/2005  WD several minor updates                                |
// v 2.0   14/06/2005  WD new falcON                                           |
// v 2.1   13/06/2005  WD changes in fieldset                                  |
// v 2.1.1 20/02/2008  WD change in body.h (removed old-style constructors)    |
// v 2.1.2 10/09/2008  WD happy gcc 4.3.1                                      |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "2.1.2"
#define falcON_VERSION_D "10-sep-2008 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#error You need NEMO to compile "mkking"
#endif
#define falcON_RepAction 0                         // no action reporting       
#define falcON_NOT_USING_PROPER                    // not using PROPRIETARY code
//-----------------------------------------------------------------------------+
#include <public/king.h>                           // my king model             
#include <body.h>                                  // my bodies                 
#include <public/io.h>                             // my NEMO file I/O          
#include <public/random.h>
#include <iomanip>                                 // C++ I/O manipulators      
#include <cmath>                                   // C++ math                  
#include <main.h>                                  // main & NEMO stuff         
using namespace falcON;
////////////////////////////////////////////////////////////////////////////////
const char*defv[] = {
  "out=???\n          output file                                        ",
  "nbody=???\n        number of bodies                                   ",
  "W0=???\n           W0 = Psi0/sigma                                    ",
  "seed=0\n           seed for the randum number generator               ",
  "time=0\n           simulation time of snapshot                        ",
  "mass=1\n           total mass                                         ",
  "r_c=1\n            core radius                                        ",
  "r_t=\n             tidal radius (don't use r_c)                       ",
  "eps=\n             output softening length eps with every body        ",
  "q-ran=f\n          use quasi- instead of pseudo-random numbers        ",
  "WD_units=f\n       input:  kpc, solar masses\n"
  "                   output: kpc, kpc/Gyr, G=1 (-> mass unit)           ",
  "outputs=f\n        give some global quantities to stderr              ",
  "table=\n           file for table of r,Phi,rho,M(<r)                  ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
const char*usage = "mkking -- initial conditions from a King model";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  //----------------------------------------------------------------------------
  // 1. create king model, i.e. integrate to get Phi(r) etc.                    
  //----------------------------------------------------------------------------
  king_model KM(getdparam("W0"));
  //----------------------------------------------------------------------------
  // 2. rescale model to mass and core/tidal radius wanted                      
  //----------------------------------------------------------------------------
  const    double vf = 0.977775320024919, mf = 2.2228847e5;
  register double mass = getdparam("mass");
  if(getbparam("WD_units")) mass /= mf;
  if(hasvalue("r_t"))
    KM.reset_scales_tidal(mass, getdparam("r_t"));
  else
    KM.reset_scales_core (mass, getdparam("r_c"));
  if(getbparam("outputs"))
    if(getbparam("WD_units"))
      std::clog<<  "# ------------------------------------------------"
	       <<"\n# King model: "
	       <<"\n# W_0   = "<<std::setw(12)<<getdparam("W0")
	       <<"\n# M_tot = "<<std::setw(12)<<KM.total_mass()
	       <<" = "<<std::setw(12)<<mf * KM.total_mass()<<" M_sun"
	       <<"\n# r_0   = "<<std::setw(12)<<KM.core_radius()    
	       <<"                kpc"
	       <<"\n# r_t   = "<<std::setw(12)<<KM.tidal_radius()
	       <<"                kpc"
	       <<"\n# r_rms = "<<std::setw(12)<<KM.rms_radius()
	       <<"                kpc"
	       <<"\n# r_h   = "<<std::setw(12)<<KM.half_mass_radius()
	       <<"                kpc"
	       <<"\n# sigma = "<<std::setw(12)<<KM.sigma()
	       <<" = "<<std::setw(12)<<vf * KM.sigma()     <<" km/s"
	       <<"\n# E_tot =" <<std::setw(13)<<KM.Etot()
	       <<" =" <<std::setw(13)<<vf*vf*mf * KM.Etot()<<" M_sun(km/s)^2"
	       <<"\n# rho_0 = "<<std::setw(12)<<KM.core_density()
	       <<" = "<<std::setw(12)<<mf*KM.core_density()<<" M_sun/kpc^3"
	       <<"\n# ------------------------------------------------"
	       <<std::endl;
    else
      std::clog<<  "# ---------------------------------"
	       <<"\n# King model: "
	       <<"\n#             W_0   = "<<getdparam("W0")
	       <<"\n#             M_tot = "<<KM.total_mass()
	       <<"\n#             r_0   = "<<KM.core_radius()
	       <<"\n#             r_t   = "<<KM.tidal_radius()
	       <<"\n#             r_rms = "<<KM.rms_radius()
	       <<"\n#             r_h   = "<<KM.half_mass_radius()
	       <<"\n#             sigma = "<<KM.sigma()
	       <<"\n#             E_tot =" <<KM.Etot()
	       <<"\n#             rho_0 = " <<KM.core_density()
	       <<"\n# ---------------------------------"<<std::endl;
  if(hasvalue("table")) KM.write_table(getparam("table"));
  //----------------------------------------------------------------------------
  // 3. create initial conditions from King model                               
  //----------------------------------------------------------------------------
  const int N = getuparam("nbody");
  const fieldset data(hasvalue("eps")? fieldset::basic | fieldset::e :
		      fieldset::basic);
  unsigned nbod[BT_NUM]={0}; nbod[bodytype::std] = N;
  snapshot        shot(getdparam("time"),nbod,data);
  const    bool   q(getbparam("q-ran"));
  const    Random Ran(getparam("seed"),6);
  const    double m = mass/double(N);
  double r,v,cth,R,phi;
  bool   again;
  int    errors=0;
  LoopAllBodies(&shot,Bi) {
    //                                                                          
    // 3.1 set mass and get r and v                                             
    //                                                                          
    Bi.mass() = m;
    do {
      again = false;
      try {
	if(q) KM.random(Ran(0), Ran(1), r, v);
	else  KM.random(Ran( ), Ran( ), r, v);
      } catch(WDutils::exception E) {
	if(++errors > 1000)
	    falcON_Error("exceeding 1000 errors \"%s\" in sampling\n", text(E));
	again = true;
      }
    } while(again);
    //                                                                          
    // 3.2 get x,y,z                                                            
    //                                                                          
    cth         = q? Ran(2,-1.,1.) : Ran(-1.,1.);
    R           = r * std::sqrt(1.-cth*cth);
    phi         = q? Ran(3,0.,TPi) : Ran(0.,TPi);
    Bi.pos()[0] = R * std::cos(phi);
    Bi.pos()[1] = R * std::sin(phi);
    Bi.pos()[2] = r * cth;
    //                                                                          
    // 3.3 get vx,vy,vz                                                         
    //                                                                          
    cth         = q? Ran(4,-1.,1.) : Ran(-1.,1.);
    R           = v * std::sqrt(1.-cth*cth);
    phi         = q? Ran(5,0.,TPi) : Ran(0.,TPi);
    Bi.vel()[0] = R * std::cos(phi);
    Bi.vel()[1] = R * std::sin(phi);
    Bi.vel()[2] = v * cth;
  }
  //                                                                            
  // 3.4 set eps                                                                
  //                                                                            
  if(hasvalue("eps")) {
    const real eps = getdparam("eps");
    LoopAllBodies(&shot,Bi) Bi.eps() = eps;
  }
  //----------------------------------------------------------------------------
  // 4. output of snapshot                                                      
  //----------------------------------------------------------------------------
  nemo_out out(getparam("out"));
  shot.write_nemo(out,data);
};
