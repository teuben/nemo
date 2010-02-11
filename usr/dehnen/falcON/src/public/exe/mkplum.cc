// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// mkplum.cc                                                                   |
//                                                                             |
// Copyright (C) 2000-2008 Walter Dehnen                                       |
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
// history:                                                                    |
//                                                                             |
// v 0.0   28/10/2003  WD created                                              |
// v 0.1   27/01/2004  WD velocity sampling always by pseudo RNG               |
// v 1.0   08/03/2004  WD changed to using sample.h; added Ossipkov-Merritt DF |
// v 1.0.1 02/04/2004  WD moved stuff into plum.h                              |
// v 1.1   12/05/2004  WD made PUBLIC; changed 'a' -> 'r_s'                    |
// v 1.1.1 20/05/2005  WD several minor updates                                |
// v 2.0   14/06/2005  WD new falcON                                           |
// v 2.1   13/06/2005  WD changes in fieldset                                  |
// v 2.2   02/05/2007  WD made Ossipkov-Merritt anisotropic model public       |
// v 2.2.1 20/02/2008  WD change in body.h (removed old-style constructors)    |
// v 2.2.2 10/09/2008  WD happy gcc 4.3.1                                      |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "2.2.2"
#define falcON_VERSION_D "10-sep-2008 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile mkplum
#endif
#define falcON_RepAction 0                         // no action reporting       
//-----------------------------------------------------------------------------+
#include <cmath>                                   // math standard library     
#include <public/sample.h>                         // N-body sampling           
#include <main.h>                                  // main & NEMO stuff         
////////////////////////////////////////////////////////////////////////////////
const char*defv[] = {
  "out=???\n          output file                                        ",
  "nbody=???\n        number of bodies                                   ",
  "r_s=1\n            scale radius                                       ",
  "mass=1\n           total mass of Plummer model                        ",
  "r_a=\n             Ossipkov-Merritt anisotropy radius                 ",
  "seed=0\n           seed for the randum number generator               ",
  "q-ran=f\n          use quasi- instead of pseudo-random numbers        ",
  "time=0\n           simulation time of snapshot                        ",
  "f_pos=0.5\n        fraction of bodies with positive sense of rotation ",
  "rmax=\n            if given, only emit bodies with r <= rmax          ",
  "WD_units=f\n       input:  kpc, M_Sun\n"
  "                   output: kpc, kpc/Gyr, G=1 (-> mass unit)           ",
  falcON_DEFV, NULL };
////////////////////////////////////////////////////////////////////////////////
const char*usage = "mkplum -- initial conditions from a Plummer model";
////////////////////////////////////////////////////////////////////////////////
namespace {
  using namespace falcON;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class PlummerModel                                                       //
  //                                                                          //
  // a scale-free Plummer model                                               //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class PlummerModel {
    const double ra, iraq, umiraq, fciraq;         // anisotropy radius, etc    
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    PlummerModel(double a= 0.) falcON_THROWING :
      ra(a), iraq(ra==0.? 0. : 1./(ra*ra)), umiraq(1-iraq), fciraq(0.4375*iraq)
    {
      if(a!=0. && a<0.75)
	falcON_THROW("PlummerModel: f(Q=1) negative for r_a/r_s < 3/4");
      if(a!=0. && a*a < 0.8125)
	falcON_Warning("PlummerModel: "
		       "f(Q) not monotonic for (r_a/r_s)^2 < 13/16");
    }
    //--------------------------------------------------------------------------
    double Ps(double r) const { return 1/sqrt(1+r*r); }
    double Rh(double r) const { return iFPit/power<5>(sqrt(1+r*r)); }
    double Mr(double r) const { return power<3>(r/sqrt(1+r*r)); }
    //--------------------------------------------------------------------------
    double rM(double m) const {
      const double x=pow(m,0.666666666666666666667);
      return sqrt(x/(1-x));
    }
    //--------------------------------------------------------------------------
    double DF(double q) const {
      const double fac = 0.15637905395236658903;
      return fac * pow(q,1.5) * (umiraq * q*q + fciraq);
    }
#ifdef falcON_PROPER
    //--------------------------------------------------------------------------
    double Re(double) const falcON_THROWING {
      falcON_THROW("PlummerModel::Re() not implemented"); return 0.; }
    //--------------------------------------------------------------------------
    double Rp(double, double) const falcON_THROWING { 
      falcON_THROW("PlummerModel::Rp() not implemented"); return 0.; }
    //--------------------------------------------------------------------------
#endif
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class ScaledPlummerModel                                                 //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class ScaledPlummerModel: 
    private PlummerModel
  {
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
    const double sR, iR;                           // scale radius, inverse     
    const double sM, iM;                           // scale mass,   inverse     
    const double sE, iE;                           // scale energy, inverse     
    const double sV, iV;                           // scale velocity, inverse   
    const double sL, iL;                           // scale ang.mom., inverse   
    const double sF;                               // scale DF                  
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    ScaledPlummerModel(double r,                   // I: scale radius           
		       double m,                   // I: total mass             
		       double a= 0.) :             //[I: anisotropy radius]     
      PlummerModel ( a/r ),
      sR  ( r               ), iR ( 1/sR ),
      sM  ( m               ), iM ( 1/sM ),
      sE  ( sM/sR           ), iE ( 1/sE ),
      sV  ( sqrt(sE)        ), iV ( 1/sV ),
      sL  ( sV * sR         ), iL ( 1/sL ),
      sF  ( sM * power<3>(iV*iR) ) 
    {}
    //--------------------------------------------------------------------------
    double Ps(double r) const { return sE*PlummerModel::Ps(r*iR); }
    double Mr(double r) const { return sM*PlummerModel::Mr(r*iR); }
    double rM(double m) const { return sR*PlummerModel::rM(m*iM); }
    double DF(double q) const { return sF*PlummerModel::DF(q*iE); }
#ifdef falcON_PROPER
    double Re(double e) const { return sR*PlummerModel::Re(e*iE); }
    double Rp(double e,
	      double l) const { return sR*PlummerModel::Rp(e*iE,l*iL); }
#endif
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class PlummerModelSampler                                                //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class PlummerModelSampler : 
    private ScaledPlummerModel,
    public  SphericalSampler
  {
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    PlummerModelSampler(double const&radius,       // I: scale radius           
			double const&Mtot,         // I: GM (untruncated)       
			double const&r_a,          // I: anisotropy radius      
			double const&rmax) :       // I: maximum radius         
      ScaledPlummerModel ( radius, Mtot , r_a ),
      SphericalSampler   ( rmax>0? Mr(rmax) : Mtot , r_a )
    {}
    //--------------------------------------------------------------------------
    // provide the abstract functions                                           
    //--------------------------------------------------------------------------
    double Ps(double r) const { return ScaledPlummerModel::Ps(r); }
    double Mr(double r) const { return ScaledPlummerModel::Mr(r); }
    double rM(double m) const { return ScaledPlummerModel::rM(m); }
    double DF(double q) const { return ScaledPlummerModel::DF(q); }
#ifdef falcON_PROPER
    double Re(double e) const { return ScaledPlummerModel::Re(e); }
    double Rp(double e,
	      double l) const { return ScaledPlummerModel::Rp(e,l); }
#endif
  };
  //////////////////////////////////////////////////////////////////////////////
} // namespace {
////////////////////////////////////////////////////////////////////////////////
void falcON::main() falcON_THROWING
{
  const double mf= 2.2228847e5;                    // mass unit for WD_units    
  const bool   WD  (getbparam("WD_units"));        // using WD_units?           
  const Random Ran (getparam("seed"),6);           // random-number-generators  
  if(hasvalue("rmax") && getdparam("rmax")<=0.) falcON_THROW("rmax <= 0");
  if(hasvalue("r_a")  && getdparam("r_a") <=0.) falcON_THROW("r_a <= 0");
  unsigned nbod[bodytype::NUM]={0}; nbod[bodytype::std] = getuparam("nbody");
  snapshot shot(getdparam("time"), nbod, fieldset(fieldset::basic));
  PlummerModelSampler PS(getdparam("r_s"),
			 WD? getdparam("mass")/mf : getdparam("mass"),
			 getdparam_z("r_a"),
			 getdparam_z("rmax"));
  PS.sample(shot,getbparam("q-ran"),Ran,getdparam("f_pos"));
  nemo_out out(getparam("out"));
  shot.write_nemo(out,fieldset::basic);
}
////////////////////////////////////////////////////////////////////////////////
