// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// plum.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2004                                               |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
// defines:                                                                    |
// class nbdy::PlummerModel;                                                   |
// class nbdy::ScaledPlummerModel;                                             |
// class nbdy::PlummerModelSampler;                                            |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_plum_h
#define falcON_included_plum_h

#ifndef falcON_included_iostream
#  include <iostream>
#  define falcON_included_iostream
#endif
#ifndef falcON_included_Pi_h
#  include <public/Pi.h>
#endif
#ifndef falcON_included_inln_h
#  include <public/inln.h>
#endif
#ifndef falcON_included_exit_h
#  include <public/exit.h>
#endif
#ifndef falcON_included_nsam_h
#  include <public/nsam.h>
#endif
#ifndef falcON_included_cmath
#  include <cmath>
#  define falcON_included_cmath
#endif
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
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
    PlummerModel(double a= 0.) :
      ra(a), iraq(ra==0.? 0. : 1./(ra*ra)), umiraq(1-iraq), fciraq(0.4375*iraq)
    {
      if(a!=0. && a<0.75)
	error("PlummerModel: f(Q=1) negative for r_a/r_s < 3/4");
      if(a!=0. && a*a < 0.8125)
	warning("PlummerModel: f(Q) not monotonic for (r_a/r_s)^2 < 13/16");
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
    double Re(double) const {
      error("PlummerModel::Re() not implemented"); return 0.; }
    //--------------------------------------------------------------------------
    double Rp(double, double) const { 
      error("PlummerModel::Rp() not implemented"); return 0.; }
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
#ifdef falcON_PROPER
			double const&r_a,          // I: anisotropy radius      
#endif
			double const&rmax) :       // I: maximum radius         
      ScaledPlummerModel ( radius, Mtot
#ifdef falcON_PROPER
			   , r_a
#endif
			   ),
      SphericalSampler   ( rmax>0? Mr(rmax) : Mtot
#ifdef falcON_PROPER
			   , r_a
#endif
			   )
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
}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_included_plum_h    
