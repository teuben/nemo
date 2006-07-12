// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// truncatedCDM.cc                                                             |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2004                                          |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// This is a non-public part of the code.                                      |
// It is property of its author and not to be made public without his written  |
// consent.                                                                    |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <pjm/pjmCDM.h>
#include <public/basic.h>
#include <utils/Pi.h>
#include <utils/inline_io.h>
#include <utils/tupel.h>
#include <utils/numerics.h>
#include <utils/spline.h>
#include <utils/WDMath.h>
#include <fstream>
#include <iomanip>
#include <ctime>

using std::pow;
using std::exp;
using std::log;

using namespace falcON;

#define SECH_TRUNC 1                               // sech(r/b) truncation      
// #undef SECH_TRUNC                               // exp (r/b) truncation      

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// auxiliary methods                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace {
  typedef tupel<4,double> quad_d;                  // quad of doubles           
  //----------------------------------------------------------------------------
  using namespace falcON;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // struct TruncFac<>                                                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template <TruncCDMModel::truncfac> struct TruncFac;
  //////////////////////////////////////////////////////////////////////////////
  // struct TruncFac<TruncCDMModel::expn>                                     //
  template<> struct TruncFac<TruncCDMModel::expn> {
    static double fac (double x) {
      return exp(-x);
    }
    //--------------------------------------------------------------------------
    static double fac (double x,
		       double&d1) {
      register double ex = exp(-x);
      d1=-ex;
      return ex;
    }
    //--------------------------------------------------------------------------
    static double fac (double x,
		       double&d1,
		       double&d2) {
      register double ex = exp(-x);
      d1=-ex;
      d2= ex;
      return ex;
    }
    //--------------------------------------------------------------------------
    static double fac (double x,
		       double&d1,
		       double&d2,
		       double&d3) {
      register double ex = exp(-x);
      d1=-ex;
      d2= ex;
      d3= d1;
      return ex;
    }
    //--------------------------------------------------------------------------
    static double fac (double x,
		       double&d1,
		       double&d2,
		       double&d3,
		       double&d4) {
      register double ex = exp(-x);
      d1=-ex;
      d2= ex;
      d3=-ex;
      d4= ex;
      return ex;
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  // struct TruncFac<TruncCDMModel::sech>                                     //
  template<> struct TruncFac<TruncCDMModel::sech> {
    static double fac (double x) {
      register double y=exp(-x);
      return twice(y)/(1+y*y);
    }
    //--------------------------------------------------------------------------
    static double fac (double x,
		       double&d1) {
      register double y=exp(-x), y2=y+y, z=y*y, f=1/(1+z);
      d1= y2*(z-1)*square(f);
      return y2*f;
    }
    //--------------------------------------------------------------------------
    static double fac (double x,
		       double&d1,
		       double&d2) {
      register double y=exp(-x), y2=y+y, z=y*y, f=1/(1+z), d=f*f;
      d1= y2*( z-1)     * d;
      d2= y2*((z-6)*z+1)*(d*=f);
      return y2*f;
    }
    //--------------------------------------------------------------------------
    static double fac (double x,
		       double&d1,
		       double&d2,
		       double&d3) {
      register double y=exp(-x), y2=y+y, z=y*y, f=1/(1+z), d=f*f;
      d1= y2*(  z- 1)           * d;
      d2= y2*(( z- 6)*z+ 1)     *(d*=f);
      d3= y2*(((z-23)*z+23)*z-1)*(d*=f);
      return y2*f;
    }
    //--------------------------------------------------------------------------
    static double fac (double x,
		       double&d1,
		       double&d2,
		       double&d3,
		       double&d4) {
      register double y = exp(-x), y2=y+y, z=y*y, f=1/(1+z), d=f*f;
      d1= y2*(   z- 1)                  * d;
      d2= y2*((  z- 6)*z+  1)           *(d*=f);
      d3= y2*((( z-23)*z+ 23)*z- 1)     *(d*=f);
      d4= y2*((((z-76)*z+230)*z-76)*z+1)*(d*=f);
      return y2*f;
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  // struct TruncFac<TruncCDMModel::suph>                                     //
  template<> struct TruncFac<TruncCDMModel::suph> {
    static double fac (double x) {
      register double y=TruncFac<TruncCDMModel::sech>::fac(x);
      return twice(y)/(1+y*y);
    }
    //--------------------------------------------------------------------------
    static double fac (double x,
		       double&d1) {
      register double
	y1,
	y = TruncFac<TruncCDMModel::sech>::fac(x,y1),
	z = y*y, f=1/(1+z);
      d1  = twice(y1)*(1-z)*f*f;
      return twice(y)*f;
    }
    //--------------------------------------------------------------------------
    static double fac (double x,
		       double&d1,
		       double&d2) {
      register double
	y1,y2,
	y = TruncFac<TruncCDMModel::sech>::fac(x,y1,y2),
	z = y*y, f=1/(1+z), d=f*f;
      d1  = twice(y1)*(1-z)*d;
      d2  = twice(y2*(1-z*z)-(y+y)*y1*y1*(3-z))*(d*=f);
      return twice(y)*f;
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class HaloModel                                                          //
  //                                                                          //
  // model for a density                                                      //
  //                                                                          //
  //                      1                                                   //
  //      rho(r) = ----------------                                           //
  //               r^g (r+1)^(b-g)                                            //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class HaloModel {
    const double g, b, bg;
    //--------------------------------------------------------------------------
  public:
    HaloModel(const double g_,
	      const double b_) : g(g_), b(b_), bg(b-g) {}
    //--------------------------------------------------------------------------
    double operator() (const double r) const
    {
      if(g == 0) return pow(r+1,-bg);
      if(g == 1) return pow(r+1,-bg) / r;
      if(g == 2) return pow(r+1,-bg) / (r*r);
      return pow(r,-g) * pow(r+1,-bg);
    }
    //--------------------------------------------------------------------------
    double operator() (const double r, double &rh1) const
    {
      register double r1=r+1, rh=pow(r,-g)*pow(r1,-bg);
      rh1 =-rh * (g/r + bg/r1);
      return rh;
    }
    //--------------------------------------------------------------------------
    double operator() (const double r, double &rh1, double &rh2) const
    {
      register double r1=r+1, rh=pow(r,-g)*pow(r1,-bg);
      rh1 =-rh * (g/r + bg/r1);
      rh2 = rh * (g*(g+1)/(r*r) + 2*g*bg/(r*r1) + bg*(bg+1)/(r1*r1));
      return rh;
    }
    //--------------------------------------------------------------------------
    double operator() (const double r,
		       double &rh1, double &rh2, double &rh3) const
    {
      register double r1=r+1, rh=pow(r,-g)*pow(r1,-bg);
      rh1 =-rh * (g/r + bg/r1);
      register double t1=g*(g+1)/(r*r), t2=bg*(bg+1)/(r1*r1);
      rh2 = rh * (t1 + 2*g*bg/(r*r1) + t2);
      rh3 = rh * (t1*((g+2)/2+3*bg/r1) + t2*(3*g/r+(bg+2)/r1));
      return rh;
    }
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class ReductionFactor                                                    //
  //                                                                          //
  // not yet used;                                                            //
  // to be used with complete DF : L^-2b g(Q=Eps-u^2*L^2/2)                   //
  //                                                                          //
  //                          1-b   2b                                        //
  //      f(r) = (1 + r^2 u^2)     r                                          //
  //                                            b                             //
  //           = (1 + r^2 u^2) (r^2/[1+r^2 u^2])                              //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class ReductionFactor {
    const double uq, b, tb;
    //--------------------------------------------------------------------------
  public:
    ReductionFactor(double u_, double b_) : uq(u_*u_), b(b_), tb(b+b) {}
    //--------------------------------------------------------------------------
    double operator() (double r) const
    {
      if(b == 0.) return 1+r*r*uq;
      if(uq== 0.) return pow(r,tb);
      const double rq=r*r, ft=1+rq*uq;
      return ft * pow(rq/ft,b);
    }
    //--------------------------------------------------------------------------
    double operator() (double r, double &d1) const
    {
      if(b == 0.) {
	d1 = twice(r*uq);
	return 1+r*r*uq;
      }
      if(uq== 0.) {
	const double p=(tb==1.)? 1. : pow(r,tb-1);
	d1 = tb*p;
	return r*p;
      }
      const double rq=r*r, ft=1+rq*uq;
      const double fc=ft * pow(rq/ft,b);
      d1 = 2*r*fc*((1-b)*uq/ft+b/rq);
      return fc;
    }
//     //--------------------------------------------------------------------------
//     double operator() (double r, double &d1, double &d2) const;
//     //--------------------------------------------------------------------------
//     double operator() (double r, double &d1, double &d2, double &d3) const;
//     //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // auxiliary variables                                                      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  double                  AB_;                     // a/b                       
  TruncCDMModel::truncfac TF_;                     // type of truncation factor 
  HaloModel               *RHO_, *RED_;            // full & reduced density    
  //----------------------------------------------------------------------------
  // truncation factor                                                          
  //----------------------------------------------------------------------------
  inline double TruncF(double x) {
    switch(TF_) {
    case TruncCDMModel::expn:
      return TruncFac<TruncCDMModel::expn>::fac(x);
    case TruncCDMModel::sech: default:
      return TruncFac<TruncCDMModel::sech>::fac(x);
    }
  };
  //----------------------------------------------------------------------------
  inline double TruncF(double x, double&d1) {
    switch(TF_) {
    case TruncCDMModel::expn:
      return TruncFac<TruncCDMModel::expn>::fac(x,d1);
    case TruncCDMModel::sech: default:
      return TruncFac<TruncCDMModel::sech>::fac(x,d1);
    }
  };
  //----------------------------------------------------------------------------
  inline double TruncF(double x, double&d1, double&d2) {
    switch(TF_) {
    case TruncCDMModel::expn:
      return TruncFac<TruncCDMModel::expn>::fac(x,d1,d2);
    case TruncCDMModel::sech: default:
      return TruncFac<TruncCDMModel::sech>::fac(x,d1,d2);
    }
  };
  //----------------------------------------------------------------------------
  inline double TruncF(double x, double&d1, double&d2, double&d3) {
    switch(TF_) {
    case TruncCDMModel::expn:
      return TruncFac<TruncCDMModel::expn>::fac(x,d1,d2,d3);
    case TruncCDMModel::sech: default:
      return TruncFac<TruncCDMModel::sech>::fac(x,d1,d2,d3);
    }
  };
  //----------------------------------------------------------------------------
  inline double TruncF(double x, double&d1, double&d2, double&d3, double&d4) {
    switch(TF_) {
    case TruncCDMModel::expn:
      return TruncFac<TruncCDMModel::expn>::fac(x,d1,d2,d3,d4);
    case TruncCDMModel::sech: default:
      return TruncFac<TruncCDMModel::sech>::fac(x,d1,d2,d3,d4);
    }
  };
  //----------------------------------------------------------------------------
  // truncated CDM model density, using static a/b & RHO                        
  //----------------------------------------------------------------------------
  inline double rho_tcdm(                          // trunced density           
			 double x)                 // I: x = r/a                
  { 
    return TruncF (x*AB_) * (*RHO_)(x);            // using global a/b          
  }
  //----------------------------------------------------------------------------
  inline double rho_tcdm(                          // trunced density           
			 double x,                 // I: x = r/a                
			 double&d1)                // O: drho/dx                
  {
    register double 
      rf1,rf=TruncF (x*AB_,rf1),
      rh1,rh=(*RHO_)(x,    rh1);
    rf1*= AB_;
    d1 = rf1*rh + rf*rh1;
    return rf*rh;
  }
  //----------------------------------------------------------------------------
  inline double rho_tcdm(                          // trunced density           
			 double x,                 // I: x = r/a                
			 double&d1,                // O: drho/dx                
			 double&d2)                // O: d^2rho/dx^2            
  {
    register double 
      rf1,rf2,rf=TruncF (x*AB_,rf1,rf2),
      rh1,rh2,rh=(*RHO_)(x,    rh1,rh2);
    rf1*= AB_;
    rf2*= square(AB_);
    d1 = rf1*rh + rf*rh1;
    d2 = rf2*rh + 2*rf1*rh1 + rf*rh2;
    return rf*rh;
  }
  //----------------------------------------------------------------------------
  inline double rho_tcdm(                          // trunced density           
			 double x,                 // I: x = r/a                
			 double&d1,                // O: drho/dx                
			 double&d2,                // O: d^2rho/dx^2            
			 double&d3)                // O: d^3rho/dx^3            
  {
    register double 
      rf1,rf2,rf3,rf=TruncF (x*AB_,rf1,rf2,rf3),
      rh1,rh2,rh3,rh=(*RHO_)(x,    rh1,rh2,rh3);
    rf1*= AB_;
    rf2*= square(AB_);
    rf3*= cube(AB_);
    d1 = rf1*rh + rf*rh1;
    d2 = rf2*rh + 2*rf1*rh1 + rf*rh2;
    d3 = rf3*rh + 3*(rf2*rh1+rf1*rh2) + rf*rh3;
    return rf*rh;
  }
  //----------------------------------------------------------------------------
  // reduced truncated CDM model density, using static a/b & RED                
  //----------------------------------------------------------------------------
  inline double red_tcdm(                          // trunced density           
			 double x)                 // I: x = r/a                
  { 
    return TruncF (x*AB_) * (*RED_)(x);            // using global a/b          
  }
  //----------------------------------------------------------------------------
  inline double red_tcdm(                          // trunced density           
			 double x,                 // I: x = r/a                
 			 double&d1)                // O: drho/dx                
  {
    register double 
      rf1,rf=TruncF (x*AB_,rf1),
      rh1,rh=(*RED_)(x,    rh1);
    rf1*= AB_;
    d1 = rf1*rh + rf*rh1;
    return rf*rh;
  }
  //----------------------------------------------------------------------------
  inline double red_tcdm(                          // trunced density           
			 double x,                 // I: x = r/a                
			 double&d1,                // O: drho/dx                
			 double&d2)                // O: d^2rho/dx^2            
  {
    register double 
      rf1,rf2,rf=TruncF (x*AB_,rf1,rf2),
      rh1,rh2,rh=(*RED_)(x,    rh1,rh2);
    rf1*= AB_;
    rf2*= square(AB_);
    d1 = rf1*rh + rf*rh1;
    d2 = rf2*rh + 2*rf1*rh1 + rf*rh2;
    return rf*rh;
  }
  //----------------------------------------------------------------------------
  inline double red_tcdm(                          // trunced density           
			 double x,                 // I: x = r/a                
			 double&d1,                // O: drho/dx                
			 double&d2,                // O: d^2rho/dx^2            
			 double&d3)                // O: d^3rho/dx^3            
  {
    register double 
      rf1,rf2,rf3,rf=TruncF (x*AB_,rf1,rf2,rf3),
      rh1,rh2,rh3,rh=(*RED_)(x,    rh1,rh2,rh3);
    rf1*= AB_;
    rf2*= square(AB_);
    rf3*= cube(AB_);
    d1 = rf1*rh + rf*rh1;
    d2 = rf2*rh + 2*rf1*rh1 + rf*rh2;
    d3 = rf3*rh + 3*(rf2*rh1+rf1*rh2) + rf*rh3;
    return rf*rh;
  }
  //============================================================================
  inline double dM(double const&lr,                // dM/dln r                  
		   double const&M)
  {
    register double r = exp(lr);                   // ln r as independent var   
    return cube(r)*rho_tcdm(r);
  }
  //----------------------------------------------------------------------------
  // static spline and function returning it.                                   
  //----------------------------------------------------------------------------
  spline<double> *SPLINE;                          // global spline             
  inline double glob_spline(const double M)
  {
    return (*SPLINE)(M);
  }
  //----------------------------------------------------------------------------
  // integrand for dPsi/dlnr and d(r^2beta rho sigma_r^2) / dlnr                
  //----------------------------------------------------------------------------
  double TBET;
  //----------------------------------------------------------------------------
  inline quad_d dP(double const&l, quad_d const&F)
  {
    const double x = exp(l), rd=red_tcdm(x), rh=rho_tcdm(x);
    register quad_d D;
    D[0] =-(*SPLINE)(l) / x;                       // dP       = M/r * dln r    
    D[1] = rd * D[0];                              // d(rd*sq) = rd * M/r dln r 
    D[2] = rh * F[0] * cube(x);                    // dW       = rh * dP * r^3  
    D[3] = F[1] * pow(x,TBET) * TBET;              // dK       = Pr * r^3 (3-2b)
    return D;
  }
}                                                  // END: unnamed namespace    
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class TruncCDMModel                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
TruncCDMModel::TruncCDMModel(double al,            // density exponent          
			     double a_,            // break radius              
			     double b_,            // sech^2 truncation radius  
			     double m_,            // maximal v_circ            
			     double be,            // anisotropy parameter      
			     truncfac t_) :        // type of truncation factor 
  A    ( al ),
  B    ( be ),
  a    ( a_ ),
  b    ( b_ ),
  tf   ( t_ ),
  RHO  ( new HaloModel(A,    3    ) ),
  RED  ( new HaloModel(A-2*B,3-2*B) ),
  ia   ( 1./a ),
  atb  ( A - B - B ),
  ab   ( a/b ),
  ln_a ( log(a) )
{
  if(B) {
    if( atb < 0 ) error("TruncCDMModel: beta > 2*alpha: unphysical model\n");
    if( B > 1.0 ) error("TruncCDMModel: beta > 1: unphysical model\n");
    if( B <-1.5 ) error("TruncCDMModel: beta < -3/2 not supported\n");
  }
  AB_  = ab;
  TF_  = tf;
  if(tf == suph)
    warning("TruncCDMModel: truncation factor suph not yet implemented;"
	    "will use sech instead\n");
  RHO_ = RHO;
  RED_ = RED;
  const double A3=3-A;
  // 1. initialize grids: lr,r,rh,m;                                            
  const double 
    ba   = b/a,
    rmin = 0.001*min(ba,1.),
    rmax = 200*max(ba,1.);
  n1 = int(200*log10(rmax/rmin));
  n  = 1+n1;
  const double dlr= log(rmax/rmin)/double(n1);
  double M  = pow(rmin,A3)/A3*(1-square(A3)/(4-A)*rmin);
  double *rd= new double[n];
  lr    = new double[n];
  r     = new double[n];
  m     = new double[n];
  r [0] = rmin;
  lr[0] = log(rmin);
  rd[0] = red_tcdm(rmin);
  m [0] = FPi * M;
  vqm  = m[0] / r[0];
  for(int i=1; i!=n; ++i) {
    lr[i] = lr[0] + i * dlr;
    r [i] = exp(lr[i]);
    rd[i] = red_tcdm(r[i]);
    M     = rk4(M,lr[i-1],dlr,dM);
    m[i]  = FPi * M;
    if(m[i] > vqm * r[i]) vqm = m[i] / r[i];
  }

  // 2. set scale factors                                          
  c   = m_ / m[n1];                                // c    : mass               
  ca  = c / a;                                     // c/a  : pot, vc^2, sigma^2 
  sca = sqrt(ca);                                  //      : velocities         
  ac  = a * c;                                     // c*a  : ang mom^2          
  cac = ca / square(a);                            // c/a^3: density            

  // 3. potential & radial velocity dispersion, Ec                              
  double
  s0     = FPi * power<3>(r[0] ) * rho_tcdm(r[ 0]),
  sN     = FPi * power<3>(r[n1]) * rho_tcdm(r[n1]);
  SPLINE = new spline<double>(n,lr,m,&s0,&sN);     // static spline: m(ln r)    
  ps     = new double[n];
  sq     = new double[n];
  ec     = new double[n];
  quad_d y;
  y[0]   = m[n1] / r[n1];
  y[1]   = m[n1]*pow(ab,2*B-4)*exp(LogGamma(2*B-4,r[n1]*ab));
  y[2]   = y[3] = 0.;
  ps[n1] = y[0];
  sq[n1] = rd[n1]? y[1] / rd[n1] : 0;
  TBET = 3-2*B;
  for(int i=n-2; i!=-1; --i) {
    y     = rk4(y,lr[i+1],-dlr,dP);
    ps[i] = y[0];
    sq[i] = rd[i]? y[1] / rd[i] : 0.;
    ec[i] = ps[i] - 0.5*m[i]/r[i];
  }
  ps0  = A<2? ps[0] + m[0]/(r[0]*(2-A)) : 0.;
  sq0  = sq[0]*rd[0] + FPi/(3-A)*Ipow(r[0],1-A-atb);
  wtot = TPi * y[2];
  ktot =-TPi * y[3];
  delete   SPLINE;
  delete[] rd;
  for(nm=0; nm!=n1; ++nm) if(m[nm+1] == m[nm]) break;
  nm++;
  // 4. compute distibution function                                            
  setgE();
}
////////////////////////////////////////////////////////////////////////////////
namespace {
  double E_,M_,aM_,B_;                             // Eps, Mtot, Mtot*a/b,beta  
  double nu;                                       // 1/(1-p)                   
  inline double intge_of_q(double q)
  {
    //   d Psi       Q^(1-p)                                                    
    // ----------- = ------- dq    with  Psi = Q (1-q^[1/(1-p)])                
    // (Q - Psi)^p    1-p                                                       
    register double psi = E_*(1-pow(q,nu));        // Psi = Eps * (1 - q^nu)    
    if     (psi <= 0.)
      return 0.;
    else if(psi < SPLINE->last_X()) {
      const double aMps=aM_/psi;
      const double Mops=M_/psi;
      const double tbet=3-B_-B_;
      const double rred=pow(Mops,-tbet)*exp(-aMps);
      const double tmp1=(tbet+aMps);;
      if(B_ > 0.5) return tmp1 * rred / psi;
      const double tmp2=-(tbet+2*aMps);
      if(B_ >-0.5) return (square(tmp1) + tmp2) * rred / square(psi);
      const double tmp3= (2*tbet+6*aMps);
      return (cube(tmp1) + 3*tmp1*tmp2 + tmp3) * rred / cube(psi);
    } else
      return (*SPLINE)(psi);
  }
  //----------------------------------------------------------------------------
  inline double intgE(double z)                    // q = z^4; dq = 4 z^3 dz    
  {
    register double f=cube(z);
    return times<4>(f)*intge_of_q(z*f);
  }
}                                                  // END: unnamced namespace   
//------------------------------------------------------------------------------
inline void TruncCDMModel::setgE()
{
  // 1. consider possible cases for beta, tabulate integrand on grid            
  double *in = new double[n];
  if       (B >= 0.5) {        //  0.5 <= beta <  1  :  tabulate d red/ dpsi    
    double rd,rd1,ps1;
    for(int i=0; i!=n; ++i) {
      rd    = red_tcdm(r[i],rd1);                 // red, red'                  
      ps1   =-m[i]/square(r[i]);                  // psi'                       
      in[i] = rd1/ps1;                            // dred/dpsi                  
    }
  } else if(B >=-0.5) {        // -0.5 <= beta <  0.5:  tabulate d^2 red/ dpsi^2
    double rh,rd,rd1,rd2,ps1,ps2;
    for(int i=0; i!=n; ++i) {
      rh    = rho_tcdm(r[i]);                     // rho                        
      rd    = red_tcdm(r[i],rd1,rd2);             // red, red', red"            
      ps1   =-m[i]/square(r[i]);                  // psi'                       
      ps2   =-2*ps1/r[i]-FPi*rh;                  // psi"                       
      in[i] = (rd2*ps1 - rd1*ps2)/cube(ps1);      // d^2red/dpsi^2              
    }
  } else if(B >=-1.5) {        // -1.5 <= beta < -0.5:  tabulate d^3 red/ dpsi^3
    double rh,rh1,rd,rd1,rd2,rd3,ps1,ps2,ps3;
    for(int i=0; i!=n; ++i) {
      rh    = rho_tcdm(r[i],rh1);                 // rho, rho'                  
      rd    = red_tcdm(r[i],rd1,rd2,rd3);         // red, red', red", red"'     
      ps1   =-m[i]/square(r[i]);                  // psi'                       
      ps2   =-2*ps1/r[i]-FPi*rh;                  // psi"                       
      ps3   = 2*(ps1/r[i]-ps2)/r[i]-FPi*rh1;      // psi"'                      
      in[i] = (ps1*(rd3*ps1-rd1*ps3)-             // d^3red/dpsi^3              
	       3*ps2*(rd2*ps1-rd1*ps2))/pow(ps1,4);
    }
  } else
    error("TruncCDMModel: beta < -1.5 not supported\n");
  // 2. compute g(E) on grid                                                    
  lg = new double[n];
  const double
    Logalfa = (1.5-B) * LogofTwo + LogBeta(0.5,1-B) - LogofPi +
              log( (B > 0.5 ? 1 : ( 0.5-B)) *
	           (B >-0.5 ? 1 : (-0.5-B)) );
  if(B == 0.5  ||  B ==-0.5 ||  B ==-1.5 )
    for(int i=0; i!=n; ++i) lg[i] = log(in[i]) - Logalfa;
  else {
    const int mm = int(1.5-B);
    const double
    p      = 1.5-B-mm,
    p1     = 1-p,
    lfc    = log(sin(p1*Pi)) - LogofPi - Logalfa,
    s0     = (atb==0 && A <2)? (mm+1-1/(2-A))*in[0]*r[0]/m[0]   :
             (atb==0 && A==2)? -in[0]/FPi                       :
             (A <2)          ? (mm+1-atb/(2-A))*in[0]*r[0]/m[0] :
             (A==2)          ? atb*in[0]/FPi                    :
                               (mm+1-atb/(2-A))*in[0]/ps[0]     ;
    nu     = 1./p1;
    SPLINE = new spline<double>(n,ps,in,&s0);
    M_     = m[n1];
    B_     = B;
    aM_    = ab*m[n1];
    for(register int i=0; i!=n; ++i) {
      E_       = ps[i];
      double err, g = qbulir(intgE,0.,1.,1.e-7,&err,0,50);
      if(g  < 0.)
	error("TruncCDMModel: very g(E) < 0");
      if(err>1.e-3)
	warning("TruncCDMModel: very inaccurate integration for g(E)\n");
      lg[i] = lfc + p1*log(E_) + log(nu*g);
    }
    delete SPLINE;
  }
  delete[] in;
  fac    = c / pow(c*a,1.5);
  ln_fac = log(fac);
}
//------------------------------------------------------------------------------
TruncCDMModel::~TruncCDMModel()
{
  if(r)   delete[] r;
  if(lr)  delete[] lr;
  if(m)   delete[] m;
  if(ps)  delete[] ps;
  if(sq)  delete[] sq;
  if(lg)  delete[] lg;
  if(ec)  delete[] ec;
  if(RHO) delete   RHO;
  if(RED) delete   RED;
}
////////////////////////////////////////////////////////////////////////////////
inline double TruncCDMModel::redx(double x) const {
  switch(tf) {
  case expn:
    return TruncFac<expn>::fac(x*ab) * (*RED)(x);
  case sech: default:
    return TruncFac<sech>::fac(x*ab) * (*RED)(x);
  }
}
//------------------------------------------------------------------------------
inline double TruncCDMModel::rhox(double x) const {
  switch(tf) {
  case expn:
    return TruncFac<expn>::fac(x*ab) * (*RHO)(x);
  case sech: default:
    return TruncFac<sech>::fac(x*ab) * (*RHO)(x);
  }
}
//------------------------------------------------------------------------------
double TruncCDMModel::rho(double rad) const {
  return cac * rhox(rad*ia);
}
//------------------------------------------------------------------------------
inline double TruncCDMModel::psix(double x) const {
  if(x<= 0.)   return ps0;
  if(x< r[0] ) {
    if(A< 2.) return  ps[0]-m[0]/(r[0]*(2-A))*pow(x/r[0],2-A);
    if(A==2.) return  FPi * log(x);
    else      return  m[0]/(r[0]*(A-2))*pow(x/r[0],2-A);
  }
  if(x> r[n1]) return m[n1]/x;
  else         return polev(log(x),lr,ps,n);
}
//------------------------------------------------------------------------------
double TruncCDMModel::psi(double rad) const {
  return ca*psix(rad*ia);
}
//------------------------------------------------------------------------------
inline double TruncCDMModel::cumx(double x) const {
  if(x<= 0.)  return 0.;
  if(x<r[0] ) return m[0]*pow(x/r[0],3-A);
  if(x>r[n1]) return m[n1];
  else        return polev(log(x),lr,m,n);
}
//------------------------------------------------------------------------------
double TruncCDMModel::cum(double rad) const {
  return c*cumx(rad*ia);
}
//------------------------------------------------------------------------------
inline double TruncCDMModel::xcum(double mr) const {
  if(mr<= 0.)   return 0.;
  if(mr<= m[0]) return r[0]*pow(mr/m[0],1/(3-A));
  if(mr>  m[n1]) error("TruncCDMModel: M>Mtot\n");
  return exp(polev(mr,m,lr,n));
}
//------------------------------------------------------------------------------
double TruncCDMModel::rM(double M) const {
  return a*xcum(M/c);
}
//------------------------------------------------------------------------------
double TruncCDMModel::vcq(double rad) const {
  if(rad <= 0.) return 0.;
  return c*cumx(rad*ia)/rad;
}
//------------------------------------------------------------------------------
double TruncCDMModel::omq(double rad) const {
  if(rad <= 0.) return 0.;
  return cum(rad)/cube(rad);
}
//------------------------------------------------------------------------------
double TruncCDMModel::kpq(double rad) const {
  if(rad <= 0.) return 0.;
  return omq(rad) + FPi * rho(rad);
}
//------------------------------------------------------------------------------
inline double TruncCDMModel::gamx(double x) const {
  if(x < r[0] ) return 2./sqrt(4-A);
  if(x > r[n1]) return 1.;
  double oq = cumx(x)/cube(x);
  return 2*sqrt(oq/(oq + FPi*rhox(x)));
}
//------------------------------------------------------------------------------
double TruncCDMModel::gam(double rad) const {
  return gamx(rad/a);
}
//------------------------------------------------------------------------------
inline double TruncCDMModel::epcx(double x) const {
  if(x<=0.)    return ps0;
  if(x< r[0])  return psix(x) - 0.5*vcq(x);
  if(x> r[n1]) return 0.5*m[n1]/x;
  else         return polev(log(x),lr,ec,n);
}
//------------------------------------------------------------------------------
double TruncCDMModel::Epc(double rad) const {
  return ca*epcx(rad/a);
}
//------------------------------------------------------------------------------
inline double TruncCDMModel::xepc(double e) const {
  if(e>=ps0)    return 0.;
  if(e> ec[0]) {
    if(A< 2.)   return r[0]*pow((ps[0]-e)/(ps[0]-ec[0]),1/(2-A));
    if(A==2.)   return exp(0.5+e/FPi);
    else        return r[0]*pow(e/ec[0],1/(2-A));
  }
  if(e< ec[n1]) return 0.5*m[n1]/e;
  else          return exp(polev(e,ec,lr,n));
}
//------------------------------------------------------------------------------
double TruncCDMModel::RcE(double e) const {
  return a*xepc(e/ca);
}
//------------------------------------------------------------------------------
inline double TruncCDMModel::xap(double e, double lq, double ce) const {
  register double
    xc  = xepc(e),
    ecc = sqrt(1.-lq/(cumx(xc)*xc));
  return xc * pow(1+ecc*ce, 0.5*gamx(xc));
}
//------------------------------------------------------------------------------
double TruncCDMModel::Rp(double e, double q) const {
  return a * xap(e/ca, q/ac, -1.);
}
//------------------------------------------------------------------------------
double TruncCDMModel::Ra(double e, double q) const {
  return a * xap(e/ca, q/ac, 1.);
}
//------------------------------------------------------------------------------
inline double TruncCDMModel::sqrx(double x) const {
  if(x <= 0.) return 0.;
  if(x<r[0] ) return (sq0-FPi/(3-A)*Ipow(x,1-A-atb))/redx(x);
  if(x>r[n1]) return m[n1]*pow(x,2*B-4)*exp(LogGamma(2*B-4,x*ab))/redx(x);
              return polev(log(x),lr,sq,n);
}
//------------------------------------------------------------------------------
double TruncCDMModel::sqr(double rad) const {
  return ca*sqrx(rad*ia);
}
//------------------------------------------------------------------------------
double TruncCDMModel::sqt(double rad) const {
  return ca*sqrx(rad*ia)*(1-B);
}
//------------------------------------------------------------------------------
inline double TruncCDMModel::lngE(double e) const {
  if(e <= ps[n] || A>2 && e>ps0) return -1.e100;
  if(e >ps[0]) {
    if(atb==0)
      return lg[0] -0.5*((A-5)*A+4)/(2-A)*log((ps0-e)/(ps0-ps[0]));
    if(A  < 2)
      return lg[0] + (B-0.5*(6-A-4*B)/(2-A))*log((ps0-e)/(ps0-ps[0]));
    if(A  ==2)
      return lg[0] + 2*(1-B)*(e-ps[0]);
    if(A  > 2)
      return lg[0] + (B+0.5*(6-A-4*B)/(A-2))*log(e/ps[0]);
  }
  return polev(e,ps,lg,n);
}
//------------------------------------------------------------------------------
inline double TruncCDMModel::gofE(double e) const {
  if(e <= ps[n] || A>2 && e>ps0) return 0.;
  return exp(lngE(e));
}
//------------------------------------------------------------------------------
double TruncCDMModel::g_E(double Eps) const {
  return fac * gofE(Eps/ca);
}
//------------------------------------------------------------------------------
double TruncCDMModel::fEL(double Eps, double Lq) const {
  if(B) return fac * gofE(Eps/ca) * pow(Lq/ac,-B);
  else  return fac * gofE(Eps/ca);
}
////////////////////////////////////////////////////////////////////////////////
#if 0
int main(int argc, char* argv[])
{
  register double alfa = (argc>1)? atof(argv[1]) : 1.;
  register double b    = (argc>2)? atof(argv[2]) : 10.;
  register double beta = (argc>3)? atof(argv[3]) : 0.;
  register double r;
  TruncCDMModel CDM(alfa,1.,b,1.,beta);
  cout << "\n" 
       << " E_tot =" << CDM.total_energy() << "\n\n";
  r = 0.5 * CDM.rad_grid(0);
  cout << "     "
       << std::setw(10) << r  << " "
       << std::setw(13) << CDM.vcq(r) * r << " "
       << std::setw(13) <<-CDM.pot(r) << " "
       << std::setw(13) << CDM.vcq(r) << " "
       << std::setw(13) << CDM.srq(r) << " "
       << std::setw(10) << 1-CDM.stq(r)/CDM.srq(r)  << " "
       << std::setw(13) << log(CDM.g_E(CDM.pot(r))) << "\n";
  for(register int i=0; i!=CDM.N_grid(); ++i)
    cout << std::setw(4)  << i <<" "
	 << std::setw(10) << CDM.rad_grid(i) << " "
	 << std::setw(13) << CDM.cum_grid(i) << " "
	 << std::setw(13) << CDM.psi_grid(i) << " "
	 << std::setw(13) << CDM.vcq_grid(i) << " "
	 << std::setw(13) << CDM.srq_grid(i) << " "
	 << std::setw(10) << 1-CDM.stq(r)/CDM.srq(r)  << " "
	 << std::setw(13) << CDM.lng_grid(i) << "\n";
  
  r = 2 * CDM.rad_grid(CDM.N_grid()-1);
  cout << "     "
       << std::setw(10) << r  << " "
       << std::setw(13) << CDM.vcq(r) * r << " "
       << std::setw(13) <<-CDM.pot(r) << " "
       << std::setw(13) << CDM.vcq(r) << " "
       << std::setw(13) << CDM.srq(r) << " "
       << std::setw(10) << 1-CDM.stq(r)/CDM.srq(r)  << " "
       << std::setw(13) << log(CDM.g_E(CDM.pot(r))) << "\n";
}
#endif
