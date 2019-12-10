// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   mkhalo.cc
///
/// \brief  creates N-body initial conditions for a quite general spherical
///         halo with an optional external potential
///
/// \author Paul McMillan
/// \author Walter Dehnen
/// \date   2000-2011
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005      Walter Dehnen, Paul McMillan
// Copyright (C) 2005-2011 Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 675
// Mass Ave, Cambridge, MA 02139, USA.
//
////////////////////////////////////////////////////////////////////////////////
//
// history:
//
// v 1.0    21/03/2005  PJM bastardized mkcdm to create mkanyhalo
// v 1.1    29/09/2005  PJM added Ossipkov-Merritt type anisotropy radius
// v 1.2    09/02/2007  WD  bastardized mkanyhalo to create this
// v 1.3    02/05/2007  WD  somewhat stable version, improved documentation 
// v 2.0    02/05/2007  WD  tabfile, core radius
// v 2.0.1  03/05/2007  WD  debugged (qbulir)
// v 2.0.2  08/05/2007  WD  r_2
// v 2.0.3  09/05/2007  WD  Mhalo -> M, beta->b 
// v 2.1    13/09/2007  WD  added parameter eps; made public
// v 2.2    19/09/2007  WD  ensured total mass
// v 2.3    01/11/2007  WD  WD_units
// v 2.3.1  02/11/2007  WD  debugged r_2 -> r_s conversion
// v 2.4    23/11/2007  WD  steeper truncation for r_t<0
// v 2.4.1  23/01/2008  WD  DF into phden, not aux data field
// v 2.4.2  08/02/2008  WD  removed minor bug with RNG seed
// v 2.4.3  20/02/2008  WD  change in body.h (removed old-style constructors)
// v 2.4.4  11/06/2008  WD  changes in error/exception treatment (falcON)
// v 2.4.5  10/09/2008  WD  happy gcc 4.3.1
// v 2.4.6  13/11/2008  WD  new mass adaption (proprietary only)
// v 2.4.7  15/12/2008  WD  debugged r_2 -> r_s conversion
// v 2.4.8  12/03/2009  WD  warning if accpars or accfile given but no accname
// v 2.4.9  03/07/2009  WD  keyword r_max
// v 2.5    08/07/2009  WD  automatically careful if DF non-monotonic
// v 2.5.1  15/03/2010  WD  inserted explanation of max_r
// v 2.6    04/08/2010  WD  keyword model
// v 2.6.1  09/12/2011  WD  improved explanation
////////////////////////////////////////////////////////////////////////////////
#define falcON_VERSION   "2.6.2"
#define falcON_VERSION_D "10-dec-2019 PJT/Walter Dehnen                        "
//------------------------------------------------------------------------------
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile mkhalo
#endif
#define falcON_RepAction 0                         // no action reporting
//------------------------------------------------------------------------------
#include <main.h>                                  // main & NEMO stuff
#include <public/halo.h>                           // halo model
#include <fstream>                                 // C++ file I/O 
#include <iomanip>                                 // C++ I/O formatting
////////////////////////////////////////////////////////////////////////////////
const char*defv[] = {
  "out=???\n          output file                                        ",
  "nbody=???\n        number of bodies                                   ",
  "model=Zhao\n       (untruncated) density model; recognised values:\n"
  "                    Plummer:   1/Model= (1+x^2)^(5/2)\n"
  "                    Jaffe:     1/Model= x^2 (1+x)^2\n"
  "                    Hernquist: 1/Model= x (1+x)^3\n"
  "                    Dehnen:    1/Model= x^inner (1+x)^(4-inner)\n"
  "                    Zhao:      1/Model= x^inner (1+x^eta)^((outer-inner)/eta)\n"
  "                    NFW:       1/Model= x (1+x)^3\n"
  "                    Moore:     1/Model= x^(3/2) (1+x^(3/2))\n"
  "                    DM:        1/Model= x^(7/9) (1+x^(4/9))^6         ",
  "inner=\n           inner density exponent                             ",
  "outer=\n           outer density exponent                             ",
  "eta=\n             transition strength between inner and outer        ",
  "r_s=1\n            scale radius                                       ",
  "r_2=\n             use r_s = r_2*((2-inner)/(outer-2))^(-1/eta)       ",
  "M=1\n              total mass of halo                                 ",
  "r_c=0\n            core radius                                        ",
  "r_t=0\n            truncation radius (0 maps to infinity)             ",
  "b=0\n              anisotropy parameter; -1.5 <= b <= min(1,g/2)      ",
  "r_a=0\n            anisotropy radius (0 maps to infinity)             ",
  "max_r=\n           radius beyond which total Pot ~= GM/r\n"
  "                   (NOT a truncation radius)\n"
  "                   use when external potential's mass extends much\n"
  "                   further than that of halo sampled                  ",
  "seed=0\n           seed for the randum number generator               ",
  "q-ran=f\n          use quasi- instead of pseudo-random numbers        ",
  "time=0\n           simulation time of snapshot                        ",
  "f_pos=0.5\n        fraction of bodies with positive sense of rotation ",
  "eps=\n             if given, set eps_i = eps                          ",
#ifdef falcON_PROPER
  "MA_rs=\n           mass adaption: scale radius                        ",
  "MA_eta=\n          mass adaption: shape parameter                     ",
  "MA_mmm=\n          mass adaption: ratio m_max/m_min                   ",
  "MA_nmax=0\n        mass adaption: max n per (E,L) (default: MA_mmm)   ",
  "MA_peri=f\n        mass adaption: use R_peri(E,L) rather than R_c(E)  ",
  "epar=\n            if given, set eps_i = epar * sqrt(m_i/<m>)         ",
#endif
  "giveF=f\n          give distribution function in phden data?          ",
  "giveP=f\n          give Phi(r) (halo+external) in pot data?           ",
  "giveA=f\n          give -dPhi/dr in acc data?                         ",
  "accname=\n         name of external monopole acceleration field       ",
  "accpars=\n         parameters of external acceleration field          ",
  "accfile=\n         file required by external acceleration field       ",
  "tabfile=\n         write table with r, M(r), rho, Psi, g(Q) to file   ",
  "prec=8\n           # significant digits in table                      ",
  "WD_units=f\n       input  units: kpc, Msun\n"
  "                   output units: kpc, 222288*Msun, Gyr (=> G=1)\n"
  "                   NOTE: external potential MUST also use these units ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
const char*usage =
  "mkhalo -- initial conditions from an equilibrium distribution function\n"
  "          with spherical density proportional to\n"
  "             Model(x) * Trunc(r/|r_t|),  x=sqrt(r^2+r_c^2)/r_s\n"
  "          with\n"
  "                         -inner   eta    (inner-outer)/eta\n"
  "             Model(x) = x       (x    + 1),\n"
  "             Trunc(z) = 1                      if r_t = 0,\n"
  "                      = sech(z)                if r_t > 0,\n"
  "                      = 2/(sech(z)+1/sech(z))  if r_t < 0.\n"
  "          The distribution function is of the form (Cuddeford 1991)\n"
  "             f(E,L) = g(Q)/L^2b\n"
  "          with  Q=-E - L^2 / 2 r^2_a  (Ossipkov 1979, Merritt 1985).\n"
  "          These models have velocity anisotropy\n"
  "             beta(r) = (r^2 + b r_a^2)/(r^2 + r_a^2).\n"
  "          If an external potential is given, the initial conditions will\n"
  "          be in equilibrium with the total potential."
#ifdef falcON_PROPER
  "\n"
  "          Individual masses are supported with :\n"
  "                                     eta                   eta\n"
  "             m propto (m_min + [r/rs]    m_max)/(1 + [r/rs]    )\n"
  "          with rs and eta independent of those for the density."
#endif
  ;
//------------------------------------------------------------------------------
namespace {
  inline double getshapeparam(const char*keyword)
  { return falcON::hasvalue(keyword)? falcON::getdparam(keyword) : -1.; }
}
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  const double mf = 2.2228847e5;                      // M_sun in WD_units
  // 1. set some parameters
  const bool WD   (getbparam("WD_units"));
  const Random Ran(getiparam("seed"),6);
  const fieldset data((
#ifdef falcON_PROPER
     hasvalue ("epar") |
#endif
     hasvalue ("eps")   ? fieldset::e : fieldset::empty) |
    (getbparam("giveF") ? fieldset::d : fieldset::empty) |
    (getbparam("giveP") ? fieldset::p : fieldset::empty) |
    (getbparam("giveA") ? fieldset::a : fieldset::empty) |
    fieldset::basic );
  const nemo_acc*mono=0;
  if(hasvalue("accname"))                             // IF(accname given) THEN
    mono = new nemo_acc(getparam  ("accname"),        //   initialize external
			getparam_z("accpars"),        //   accelerations
			getparam_z("accfile"));
  else if(hasvalue("accpars"))
    falcON_Warning("'accpars' given but no 'accname': "
		   "will have no external potential\n");
  else if(hasvalue("accfile"))
    falcON_Warning("'accfile' given but no 'accname': "
		   "will have no external potential\n");
  // 2. create initial conditions from a halo model using mass adaption
  DoublePowerLawHalo::Model model=DoublePowerLawHalo::model(getparam("model"));
  double gi = DoublePowerLawHalo::inner_value(model,getshapeparam("inner"));
  double go = DoublePowerLawHalo::outer_value(model,getshapeparam("outer"));
  double et = DoublePowerLawHalo::trans_value(model,getshapeparam("eta"));
  double rs = getdparam("r_s");
  double rc = getdparam("r_c");
  double rt = getdparam("r_t");
  double ra = getdparam("r_a");
  double be = getdparam("b");
  double Mt = getdparam("M");
  if(WD) Mt/= mf;
  if(hasvalue("r_2")) {
    if(rs!=1.)
      falcON_Warning("'r_2' given: will ignore non-default 'r_s' value\n");
    if(go<=2.)
      falcON_Error("outer<=2 ==> gamma(r)<2 at finite r; cannot use r_2\n");
    rs = getdparam("r_2") * pow((2-gi)/(go-2),-1./et);
  }
  ModifiedDoublePowerLawHalo Halo(rs,rc,rt,Mt,gi,go,et);
  HaloSampler HaloSample(Halo,be,ra,mono,
#ifdef falcON_PROPER
			 getdparam_z("MA_rs"),
			 getdparam_z("MA_mmm"),
			 getdparam_z("MA_eta"),
			 getdparam  ("MA_nmax"),
			 getbparam  ("MA_peri"),
#endif
			 getdparam_z("max_r"));
  // 3. sample snapshot and write to output
  unsigned N = getuparam("nbody");
  if(N) {
    unsigned nbod[bodytype::NUM]={0}; nbod[bodytype::std] = N;
    snapshot shot(getdparam("time"), nbod, data);
    HaloSample.sample(shot, getbparam("q-ran"), Ran,
		      getdparam("f_pos"), 
#ifdef falcON_PROPER
		      getdparam("epar"),
#endif
		      getbparam("giveF"),getbparam("giveP"),getbparam("giveA"));
    if(hasvalue("eps")
#ifdef falcON_PROPER
       && !hasvalue("epar")
#endif
       ) {
      real epsi=getrparam("eps");
      LoopAllBodies(&shot,B) B.eps()=epsi;
    }
    if(Mt != HaloSample.Mt()) {
      const real mfac = Mt/HaloSample.Mt(), vfac=sqrt(mfac);
      LoopAllBodies(&shot,B) {
	B.mass() *= mfac;
	B.vel () *= vfac;
      }
    }
    nemo_out out(getparam("out"));
    shot.write_nemo(out,data);
  }
  // 4. optional output of table
  if(hasvalue("tabfile")) {
    HaloModel const&HM (HaloSample.Model());
    output   Tab(getparam("tabfile"));
    if(!Tab) return;
    const int p(getiparam("prec"));
    const int w(p+5);
    char space[20] = "               ";
    space[w>10? w-10 : 0] = 0;
    Tab  <<"#---------------------------------------------------------------\n";
    if(RunInfo::cmd_known())
      Tab<< "# file generated by command:\n"
	 << "# \""<<RunInfo::cmd() <<"\"\n";
    Tab  << "# run at "<<RunInfo::time()<<'\n';
    if(RunInfo::user_known())
      Tab<< "#     by \""<<RunInfo::user()<<"\"\n";
    if(RunInfo::host_known())
      Tab<< "#     on \""<<RunInfo::host()<<"\"\n";
    if(RunInfo::pid_known())
      Tab<< "#     pid "<<RunInfo::pid()<<'\n';
    Tab  << "#\n"
	 << "# table for halo model with density\n"
	 << "#\n"
	 << (Halo.trunc_radius()? 
	    "#                        C sech(r/r_t)\n" :
            "#                              C\n")
	 << "#     rho(r) = ----------------------------------\n"
	 << "#               inner   eta     [outer-inner]/eta\n"
	 << "#              x      (x    + 1)\n#\n"
	 << (Halo.core_radius()?
	    "# with x = sqrt(r^2 + r_c^2) / r_s\n" :
	    "# with x = r/r_s\n")
	 << "# and\n";
    if(Halo.core_radius()) 
      Tab<< "#     r_c   = "<<Halo.core_radius()<<(WD? " kpc\n" : "\n");
    Tab  << "#     r_s   = "<<Halo.scale_radius()<<(WD? " kpc\n" : "\n");
    if(Halo.trunc_radius())
      Tab<< "#     r_t   = "<<Halo.trunc_radius()<<(WD? " kpc\n" : "\n");
    Tab  << "#     inner = "<<gi<<"\n"
	 << "#     outer = "<<go<<"\n"
	 << "#     eta   = "<<et<<"\n"
	 << "#     C     = "<<Halo.rho0();
    if(WD)
      Tab<< " = "<<(Halo.rho0()*mf)<<" M_sun/kpc^3";
    Tab  << '\n'
	 << "#     M_tot = "<<Halo.total_mass();
    if(WD)
      Tab<< " = "<<(Halo.total_mass()*mf)<<" M_sun";
    Tab  << '\n';
    if(hasvalue("accname")) {
      Tab<<"# and embedded in the external potential \"" << getparam("accname")
	 <<"\"\n";
      if(hasvalue("accpars"))
	Tab<<"# with parameters \"" << getparam("accpars") <<"\"\n";
      if(hasvalue("accfile"))
	Tab<<"# with auxiliary file \"" << getparam("accfile") <<"\"\n";
    }
    Tab  << "#\n"
         << "# Psi(0) = "<< setprecision(p) << HM.Ps(0.);
    if(WD)
      Tab<< "(kpc/Gyr)^2";
    Tab  << '\n'
	 << "#\n#"
	 << space << "      r  "
	 << space << "   M_h(<r) "
	 << space << "   M_t(<r) "
	 << space << "    rho(r) "
	 << space << "    Psi(r) "
	 << space << "ln[g(Q=Psi)]\n"
	 << "#\n";
    for(int n=0; n!=HM.Ntab(); ++n)
      Tab<< setw(w) << setprecision(p) << HM.rad(n) << ' '
	 << setw(w) << setprecision(p) << HM.mhr(n) << ' '
	 << setw(w) << setprecision(p) << HM.mtr(n) << ' '
	 << setw(w) << setprecision(p) << Halo(HM.rad(n)) << ' '
	 << setw(w) << setprecision(p) << HM.psi(n) << ' '
	 << setw(w) << setprecision(p) << HM.lng(n) << '\n';
  }
}
////////////////////////////////////////////////////////////////////////////////
