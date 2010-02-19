// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// mkdehnen.cc                                                                 |
//                                                                             |
// Copyright (C) 2004-2008 Walter Dehnen                                       |
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
// creates N-body initial conditions from a Dehnen (1993) model                |
//                                                                             |
//-----------------------------------------------------------------------------+
//
// history:
//
// v 1.0   26/02/2004  WD created
// v 1.1   27/02/2004  WD added some consistency checks & error messages
// v 1.2   01/03/2004  WD corrected asymptotic of f at Q=0 with OM anisotropy
// v 1.3   01/03/2004  WD added option f_pos; re-arranged order of options
// v 2.0   03/03/2004  WD using sample.h; fixed bug with scaling & Rad[]
// v 2.1   05/03/2004  WD re-written gama; improved scaling support
// v 2.2   26/03/2004  WD bug fixed in this file (f(E=0)=nan -> sampling error)
// v 2.2.1 02/04/2004  WD put class DehnenModelSampler into gamma.h
// v 2.3   04/05/2004  WD happy icc 8.0; new body.h; new make.proper
// v 2.4   12/05/2004  WD made partly PUBLIC
// v 2.5   27/10/2004  WD option giveF
// v 2.5.1 20/05/2005  WD several minor updates
// v 3.0   10/06/2005  WD new falcON: new bodies, new nemo_io 
// v 3.1   06/07/2005  WD default rmax=1000 r_s
// v 3.2   13/06/2005  WD changes in fieldset
// v 3.3   02/05/2007  WD made Ossipkov-Merritt anisotropic models public
// v 3.3.1 23/01/2008  WD DF in phden (previous: aux) if giveF=true
// v 3.3.2 20/02/2008  WD change in body.h (removed old-style constructors)
// v 3.3.3 10/09/2008  WD happy gcc 4.3.1
// v 3.3.4 13/11/2008  WD new mass adaption (proprietary only)
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "3.3.4"
#define falcON_VERSION_D "13-nov-2008 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#error You need NEMO to compile mkdehnen
#endif
//-----------------------------------------------------------------------------+
#include <public/gamma.h>                          // WD's gamma models         
#include <main.h>                                  // main & NEMO stuff         
////////////////////////////////////////////////////////////////////////////////
const char*defv[] = {
  "out=???\n          output file                                        ",
  "nbody=???\n        number of bodies                                   ",
  "gamma=???\n        inner power-law slope of density, gamma in[0,2.5[  ",
  "mass=1\n           total mass                                         ",
  "r_s=1\n            scale radius                                       ",
  "r_a=\n             Ossipkov-Merritt anisotropy radius                 ",
  "seed=0\n           seed for the randum number generator               ",
  "q-ran=f\n          use quasi- instead of pseudo-random numbers        ",
  "time=0\n           simulation time of snapshot                        ",
  "f_pos=0.5\n        fraction of bodies with positive sense of rotation ",
  "rmax=1000\n        if != 0, only emit bodies with r <= rmax * r_s     ",
  "giveF=f\n          give distribution function in phden data?          ",
#ifdef falcON_PROPER
  "MA_rs=\n           mass adaption: scale radius                        ",
  "MA_eta=\n          mass adaption: shape parameter                     ",
  "MA_mmm=\n          mass adaption: ration m_max/m_min                  ",
  "MA_nmax=1\n        mass adaption: max n per (E,L)                     ",
  "MA_peri=f\n        mass adaption: use R_peri(E,L) rather than R_c(E)  ",
  "epar=\n            if given, set eps_i = epar * sqrt(m_i/<m>)         ",
#endif
  "WD_units=f\n       input:  kpc, solar masses\n"
  "                   output: kpc, kpc/Gyr, G=1 (-> mass unit)           ",
  falcON_DEFV, NULL };
////////////////////////////////////////////////////////////////////////////////
const char*usage = "mkdehnen -- initial conditions from a Dehnen (1993) model"
                 "\n            possibly with Osipkov-Merritt anisotropy"
#ifdef falcON_PROPER
                 "\n            and mass adaption (proprietary version only):"
                 "\n\n                     m_min + (r/rs)^eta m_max"
                 "\n            m propto ------------------------"
                 "\n                        1  + (r/rs)^eta\n";
#endif
;
////////////////////////////////////////////////////////////////////////////////
void falcON::main() falcON_THROWING
{
  const double mf = 2.2228847e5;                   // M_sun in WD_units         
  //----------------------------------------------------------------------------
  // 1. set some parameters                                                     
  //----------------------------------------------------------------------------
  if(hasvalue("r_a")  && getdparam("r_a") <=0.)
    falcON_THROW("r_a <= 0");
  const bool     WD  (getbparam("WD_units"));      // using WD_units?           
  const Random   Ran (getparam("seed"),6);         // random-number-generators  
  const fieldset data( 
    (getbparam("giveF")? fieldset::d : fieldset::empty) |
#ifdef falcON_PROPER
    (hasvalue("epar")  ? fieldset::e : fieldset::empty) |
#endif
    fieldset::basic);
  //----------------------------------------------------------------------------
  // 2. create initial conditions from a Dehnen model using mass adaption       
  //----------------------------------------------------------------------------
  double rs(getdparam("r_s"));
  double rmax(getdparam("rmax"));
  if(rmax > 0.) rmax *= rs;
  else          rmax  = 0.;
  DehnenModelSampler GS(getdparam("gamma"), rs,
			WD? getdparam("mass")/mf : getdparam("mass"),
			getdparam_z("r_a"),
			rmax,9999,1.e-8
#ifdef falcON_PROPER
		       ,getdparam_z("MA_rs"),
			getdparam_z("MA_eta"),
			getdparam_z("MA_mmm"),
			getdparam  ("MA_nmax"),
			getbparam  ("MA_peri")
#endif
			);
  unsigned nbod[bodytype::NUM]={0}; nbod[bodytype::std] = getuparam("nbody");
  snapshot shot(getdparam("time"), nbod, data);
  GS.sample(shot,getbparam("q-ran"),Ran,getdparam("f_pos"),
#ifdef falcON_PROPER
	    getdparam("epar"),
#endif
	    getbparam("giveF"));
  //----------------------------------------------------------------------------
  // 3. output of snapshot                                                      
  //----------------------------------------------------------------------------
  nemo_out out(getparam("out"));
  shot.write_nemo(out,data);
}
