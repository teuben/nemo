// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// timer.h                                                                     |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2003-2004,2013                                     |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// timer functions for external potentials                                     |
//                                                                             |
// various functional forms are supported:                                     |
//                                                                             |
// index   name         functional form                                        |
// -----+--------------+-----------------------------------------------------  |
//  0   | adiabatic    |  a smooth transition from 0 at t<t0 to 1 at t>t0+tau  |
//  1   | saturate     |  1-exp([t-t0]/tau)                                    |
//  2   | quasi-linear |  sqrt(x^2+1) - 1, x=[t-t0]/tau for t>t0, 0 for t<t0   |
//  3   | linear       |  [t-t0]/tau for t>t0, 0 for t<t0                      |
//                                                                             |
// initialization: void   init_timer(int index, double t0, double tau);        |
// call:           double timer_double(double t);                              |
// call:           float  timer_float (float  t);                              |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// Versions                                                                    |
// 0.0   20-jan-2004 WD   created from $NEMO/src/orbits/potential/data/timer.h |
// 0.1   30-jun-2004 WD   changed to non-template, anonymous namespace         |
// 0.2   22-nov-2013 WD   added exponential growth                             |
//-----------------------------------------------------------------------------+
#ifndef included_timer_h
#define included_timer_h

#include <cmath>
////////////////////////////////////////////////////////////////////////////////
namespace {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class timer                                                              //
  //                                                                          //
  // implements four different ways to grow an amplitude A(t).                //
  //                                                                          //
  // 1. ADIABATIC                                                             //
  // a slow change from A=0 at t<=t0 to A=1 at t>=t1:=t0+tau:                 //
  //                                                                          //
  //      3   5    5   3    15       1           t - t0                       //
  // A = --- x  - --- x  + ---- x + ---;  x = 2 -------- - 1;                 //
  //     16        8        16       2            tau                         //
  //                                                                          //
  // for t in [t0,t0+tau]                                                     //
  //                                                                          //
  // 2. SATURATE                                                              //
  // an initially linear increase that saturates at t=oo to A=1               //
  //                                                                          //
  //              t0 - t                                                      //
  // A = 1 - exp(--------);                                                   //
  //               tau                                                        //
  //                                                                          //
  // 3. QUASI LINEAR                                                          //
  // an eventually linear increase which starts of with zero dA/dt            //
  //                                                                          //
  //            2                  t - t0                                     //
  // A = sqrt( x  + 1 ) - 1;  x = --------;                                   //
  //                                tau                                       //
  // 4. LINEAR                                                                //
  // a simple linear grow for t>t0                                            //
  //                                                                          //
  //          t - t0                                                          //
  // A = x = --------;                                                        //
  //            tau                                                           //
  //                                                                          //
  // 5. EXPONENTIAL                                                           //
  // an exponential grow minus initial                                        //
  //                                                                          //
  //   = 0                              for t0 < t                            //
  // A = fac * exp([t-t0]/tau)-1        for t0 < t < t1                       //
  //   = 1                              for      t > t1                       //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class timer {
  public:
    enum index {
      adiabatic    = 0,
      saturate     = 1,
      quasi_linear = 2,
      linear       = 3,
      exponential  = 4,
      constant     = 9
    };
    //--------------------------------------------------------------------------
    // some useful inline functions                                             
    //--------------------------------------------------------------------------
    static double twice (double const&x) { return x+x; }
    static double square(double const&x) { return x*x; }
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
  private:
    index  in;
    double t0, tau, itau, t1, fac;                 // t_0, tau, 1/tau, t_1      
    //--------------------------------------------------------------------------
    // functions for A(t)                                                       
    //--------------------------------------------------------------------------
    double _m_adiabatic(double t) const {
      if(t <= t0) return 0.;
      if(t >= t1) return 1.;
      const double a0=0.1875, a1=0.625, a2=0.9375;
      register double
	xi  = twice(t-t0)*itau-1,
	xq  = xi*xi;
      xi *= ((a0*xq - a1)*xq + a2);
      return 0.5 + xi;
    }
    //--------------------------------------------------------------------------
    double _m_saturate(double t) const {
      if(t <= t0) return 0.;
      return 1 - std::exp((t0-t)*itau);
    }
    //--------------------------------------------------------------------------
    double _m_quasi_linear(double t) const {
      if(t <= t0) return 0.;
      return std::sqrt(square((t-t0)*itau)+1)-1;
    }
    //--------------------------------------------------------------------------
    double _m_linear(double t) const {
      return t<=t0? 0. : (t-t0)*itau;
    }
    //--------------------------------------------------------------------------
    double _m_exponential(double t) const {
      return 
	t <= t0?  0 :
	t >= t1?  1 :
	fac * (std::exp(itau*(t-t0))-1);
    }
    //--------------------------------------------------------------------------
    // construction and initialization                                          
    //--------------------------------------------------------------------------
  public:
    void init(index  _in,
	      double _t0,
	      double _ta,
	      double _t1=0.0)
    {
      in   = _in;
      t0   = _t0;
      tau  = _ta;
      itau = tau==0.? 0. : 1./tau;
      t1   = in==exponential? _t1 : t0+tau;
      fac  = in==exponential? 1/(std::exp(itau*(t1-t0))-1) : 1;
    }
    timer() {}
    timer(index  _in,
	  double _t0,
	  double _ta,
	  double _t1=0.0)
    { init(_in,_t0,_ta,_t1); }
    //--------------------------------------------------------------------------
    // returning the amplitude function                                         
    //--------------------------------------------------------------------------
    double operator()(double t) const {
      if (tau== 0) return 1.;
      switch(in) {
      case saturate:     return _m_saturate(t);
      case quasi_linear: return _m_quasi_linear(t);
      case linear:       return _m_linear(t);
      case adiabatic:    return _m_adiabatic(t);
      case exponential:  return _m_exponential(t);
      default:           return 1.;
      }
    }
    //--------------------------------------------------------------------------
    // a description of the amplitude function                                  
    //--------------------------------------------------------------------------
    static const char* describe(index i) {
      switch(i) {
      case saturate:     
	return "saturate: A(t) = 1-exp([t0-t]/tau)";
      case quasi_linear: 
	return "quasi-linear: A(t) = sqrt(x^2+1)-1, x=(t-t0)/tau";
      case linear:       
	return "linear: A(t) = (t-t0)/tau";
      case adiabatic:           
	return "adiabatic: a smooth transition from A(t<t0)=0 to A(t>t0+tau)=1";
      case exponential:
	return "exponential A(t) = f*{exp([t-t0]/tau)-1}, for t0<t<t1; A(t1)=1";
      default:
	return "none: A == 1 at all times";
      }
    }
    //--------------------------------------------------------------------------
    const char* describe() const {
      return describe(in);
    }
    //--------------------------------------------------------------------------
    double const& T0 () const { return t0; }
    double const& TAU() const { return tau; }
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
}                                                  // anonymous namespace       
////////////////////////////////////////////////////////////////////////////////
#endif                                             // included_timer_h          

