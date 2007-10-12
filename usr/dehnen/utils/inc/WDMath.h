// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    utils/inc/WDMath.h                                                 
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    1994-2006                                                          
///                                                                             
/// \todo    complete doxygen documentation                                     
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 1994-2006  Walter Dehnen                                       
//                                                                              
// This program is free software; you can redistribute it and/or modify         
// it under the terms of the GNU General Public License as published by         
// the Free Software Foundation; either version 2 of the License, or (at        
// your option) any later version.                                              
//                                                                              
// This program is distributed in the hope that it will be useful, but          
// WITHOUT ANY WARRANTY; without even the implied warranty of                   
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
// General Public License for more details.                                     
//                                                                              
// You should have received a copy of the GNU General Public License            
// along with this program; if not, write to the Free Software                  
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                    
//                                                                              
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Contents:                                                                    
//                                                                              
//   mathematical constants                                                     
//   volume of unit sphere in N dimensions                                      
//   integral of power                                                          
//   logarithms and exponentiations: natural and to bases 2 and 10              
//   sech(x)                                                                    
//   sincos: sin(x), cos(x)                                                     
//   log of sin,cos,sinh,cosh for complex arguments                             
//   complete and incomplete Gamma functions                                    
//   complete and incomplete Beta functions                                     
//   exponential integrals Ei and En                                            
//   Bessel functions J, Y, I, K                                                
//   Hermite polynomials and their normalisization                              
//                                                                              
////////////////////////////////////////////////////////////////////////////////
#ifndef WDutils_included_WDMath_h
#define WDutils_included_WDMath_h

#ifndef WDutils_included_cmath
# include <cmath>
# define WDutils_included_cmath
#endif
#ifndef WDutils_included_cstdlib
# include <cstdlib> 
# define WDutils_included_cstdlib
#endif

#if  defined(__COMPLEX__)					\
  || defined(_CPP_COMPLEX)					\
  || defined(__STD_COMPLEX)					\
  || defined(__PGCC__) && defined(_STLP_template_complex)	\
  || defined(__GNUC__) && defined(_GLIBCXX_COMPLEX)		\
  || defined(WDutils_included_complex)
#  define WDutils_COMPLEX
#endif

#ifndef WDutils_included_Pi_h
# include <Pi.h>
#endif
#ifndef WDutils_included_basic_h
# include <exception.h>
#endif
#ifndef WDutils_included_traits_h
# include <traits.h>
#endif
#ifndef WDutils_included_inline_h
# include <inline.h>
#endif

#ifdef __GNUC__
extern "C" {
  // from math.h
  double hypot(double,double);
  float  hypotf(float,float);
}
#endif

namespace WDutils {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// \name some constants                                                      
  //@{                                                                          
  // ///////////////////////////////////////////////////////////////////////////
  const double EulerGamma  = 0.577215664901532860606512090082; ///< Euler's Gam
  const double LogofTwo    = 0.693147180559945309417232121458; ///< ln(2)
  const double LogofPi     = 1.144729885849400174143427351353; ///< ln(Pi)
  const double LogofTwoInv = 1.442695040888963407359924681002; ///< 1/ln(2)
  const double LogofTen    = 2.302585092994045684017991454684; ///< ln(10)
  const double LogofTenInv = 0.434294481903251827651128918917; ///< 1/ln(10)
  //@}
  // ///////////////////////////////////////////////////////////////////////////
  /// volume of unit sphere in n dimensions                                     
  double SphVol(int);
  // ///////////////////////////////////////////////////////////////////////////
  /// integral of power (inline)                                                
  inline double Ipow(double x, double p)
  {
    register double p1=p+1;
    if(p1) return std::pow(x,p1)/p1;
    else   return std::log(x);
  }
  // ///////////////////////////////////////////////////////////////////////////
  /// hypotenus of x,y                                                          
#ifdef __GNUC__
  inline double hypot(double x, double y) { return ::hypot (x,y); }
  inline float  hypot(float  x, float  y) { return ::hypotf(x,y); }
#else
  template<typename X> inline
  X hypot(X x, X y) {
    X ax=abs(x), ay=abs(y);
    return ax>ay? ax*sqrt(1+square(ay/ax)) : ay*sqrt(1+square(ax/ay));
  }
#endif
  // ///////////////////////////////////////////////////////////////////////////
  /// fast inverse square root.
  ///
  /// implementation of a fast inverse sqrt function.
  /// 
  /// described in Chris Lomont, "Fast Inverse Square Root", see
  /// www.lomont/Math/Papers/2003/InvSqrt.pdf
  ///
  /// \version       2003 CL  initial C routine
  /// \version 26-09-2007 WD  adopted to C++
  /// \version 26-09-2007 WD  using union to give good code under gcc -O2
  /// \version 26-09-2007 WD  two Newton-Raphson iterations
  ///
  /// \author Chris Lomont, Walter Dehnen
  ///
  /// \param x argument
  /// \return approximation to 1/sqrt(x), 5-6 digits accurate
  inline float invsqrt(float x)
  {
    register union { float y; int i; } R; // union to manipulate bits via int i
    R.y = x;
    register float xhalf = 0.5f*R.y;      // take x/2
    R.i = 0x5f375a86 - (R.i>>1);          // manipulate bits: get initial guess
    R.y*= 1.5f-xhalf*R.y*R.y;             // 1st Newton Raphson step
    R.y*= 1.5f-xhalf*R.y*R.y;             // 2nd Newton Raphson step
    return R.y;                           // more steps not sensible for floats
  }
  /// inverse square root in double precision: 1/std::sqrt() is fastest
  inline double invsqrt(double x) {
    return 1./std::sqrt(x);
  }
//   inline double invsqrt(double x)
//   {
//     register union { double y; long long int i; } R;
//     R.y = x;
//     register double xhalf = 0.5*R.y;
//     R.i = 0x5fe6ec85e7de30daLL - (R.i>>1);
//     R.y*= 1.5-xhalf*R.y*R.y;
//     R.y*= 1.5-xhalf*R.y*R.y;
//     R.y*= 1.5-xhalf*R.y*R.y;
//     return R.y;
//   }
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// \name Log's and Exp's (inlines)                                           
  //@{                                                                          
  // ///////////////////////////////////////////////////////////////////////////
  /// natural logarithm
  inline double ln(double x) {
    if(x<=0) WDutils_ErrorF("argument <= 0","ln()");
    return std::log(x);
  }
  /// logarithm to base 2
  inline double ld(double x) {
    if(x<=0) WDutils_ErrorF("argument <= 0","ld()");
    return LogofTwoInv*std::log(x);
  }
  /// logarithm to base 10
  inline double lg(double x) {
    if(x<=0) WDutils_ErrorF("argument <= 0","lg()");
    return std::log10(x);
  }
  /// ten to the power \a x
  inline double Tento(double x) {
    return std::exp(LogofTen*x);
  }
  /// two to the power \a x
  inline double Twoto(double x) {
    return std::exp(LogofTwo*x);
  }
  /// secans hyperbolicus
  inline double sech(double x) {
    register double ex = x<0 ? exp(x) : exp(-x);
    return 2.*ex/(1+ex*ex);
  }
  template<typename T> struct __sincos {
    static void sc(T x, T&, T&) {
      WDutils_THROW("sincos() of \"%s\" called\n",traits<T>::name());
    }
  };
  template<> struct __sincos<float> {
    static void sc(float x, float&s, float&c) {
#if defined(__GNUC__) || defined (__INTEL_COMPILER) || defined (__PGCC__)
    __asm __volatile__ ("fsincos" : "=t" (s), "=u" (c) : "0" (x) );
#else
    s = std::sin(x);
    c = std::cos(x);
#endif
    }
  };
  template<> struct __sincos<double> {
    static void sc(double x, double&s, double&c) {
#if defined(__GNUC__) || defined (__INTEL_COMPILER) || defined (__PGCC__)
    __asm __volatile__ ("fsincos" : "=t" (s), "=u" (c) : "0" (x) );
#else
    s = std::sin(x);
    c = std::cos(x);
#endif
    }
  };
  /// sinus and cosinus of real-valued scalar simultaneously
  template<typename REAL> inline void sincos(REAL x, REAL&s, REAL&c) {
    __sincos<REAL>::sc(x,s,c);
  }
  //@}
#ifdef WDutils_COMPLEX
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// \name logarithms of complex trigonometric and hyperbolic functions        
  //@{                                                                          
  // ///////////////////////////////////////////////////////////////////////////
  std::complex<double> lnsin (std::complex<double> const&); ///< ln(sin(x))
  std::complex<double> lncos (std::complex<double> const&); ///< ln(cos(x))
  std::complex<double> lnsinh(std::complex<double> const&); ///< ln(sinh(x))
  std::complex<double> lncosh(std::complex<double> const&); ///< ln(cosh(x))
  //@}
#endif
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// \name Gamma functions                                                     
  ///                                                                           
  /// \note Gamma(z)   = int( exp(-t) t^(z-1), t=0..oo )                        
  /// \note Gamma(a,x) = int( exp(-t)*t^(a-1), t=x..oo )                        
  /// \note gamma(a,x) = int( exp(-t)*t^(a-1), t=0..x  )                        
  //@{                                                                          
  // ///////////////////////////////////////////////////////////////////////////
  double LogGamma(double);   ///< ln(Gamma(x))
#ifdef WDutils_COMPLEX
  ///< ln(Gamma(z))    
  std::complex<double> LogGamma(std::complex<double> const&);
#endif
  double LogGamma(double, double);  ///< ln(Gamma(a,x))
  double Loggamma(double, double);  ///< ln(gamma(a,x))  a>0
  double GammaP  (double, double);  ///< P(a,x):=gamma(a,x)/Gamma(a)
  double GammaQ  (double, double);  ///< Q(a,x):=Gamma(a,x)/Gamma(a)
  //@}
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// \name Beta functions                                                      
  ///                                                                           
  /// \note Beta  (a,b) = int( t^(a-1) (1-t)^(b-1), t=0..1 )                    
  /// \note Beta_x(a,b) = int( t^(a-1) (1-t)^(b-1), t=0..x )                    
  //@{                                                                          
  // ///////////////////////////////////////////////////////////////////////////
  double LogBeta (double a, double b);           ///< ln(Beta(a,b))
  double Beta    (double a, double b);           ///< Beta(a,b)
  double Beta    (double a, double b, double x); ///< Beta_x(a,b)
  /// class for incomplete Beta function at fixed (a,b)
  class BetaFunc {
    double a,b,B,x0;
  public:
    /// construction: get a,b
    BetaFunc(double a, double b);
    double const& operator() () const { return B;} ///< Beta(a,b)
    double operator() (double) const;              ///< Beta_x(a,b)
  };
  //@}
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// \name  Exponential integrals                                              
  //@{                                                                          
  // ///////////////////////////////////////////////////////////////////////////
  double En(int, double);                  ///< E_n(x)
  double Ei(double);                       ///< Ei(x)
  //@}
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// \name Bessel function                                                     
  //@{                                                                          
  // ///////////////////////////////////////////////////////////////////////////
  double J0(double);                       ///< J_0(x)
  double J1(double);                       ///< J_1(x)
  double Jn(unsigned, double);             ///< J_n(x)
  double Y0(double);                       ///< Y_0(x) [GR: N_0(x)]
  double Y1(double);                       ///< Y_1(x) [GR: N_1(x)]
  double Yn(unsigned, double);             ///< Y_n(x) [GR: N_n(x)]
  double I0(double);                       ///< I_0(x)
  double I1(double);                       ///< I_1(x)
  double In(unsigned, double);             ///< I_n(x)
  double K0(double);                       ///< K_0(x)
  double K1(double);                       ///< K_1(x)
  double Kn(unsigned, double);             ///< K_n(x)
  //@}
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// \name Hermite polynomials                                                 
  //@{                                                                          
  // ///////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  double HermiteH(unsigned, double);         ///< H_n(x)              
  void   HermiteH(unsigned, double, double*);///< H_i(x), i=0...n     
  /// gives the inverse squared normalization constants for the H_n(x):
  /// Int dx Exp[-x^2] H_n(x) H_m(x) = N_n delta_{nm}
  void   NormSqHermite(unsigned, double*);
  /// returns the nth Hermite polynomial
  /// normalized to be orthonormal w.r.t. the weight function exp(-x^2)
  double HermiteH_normalized(unsigned, double);
  /// evaluates the Hermite polynomials 0 to n
  /// normalized to be orthonormal w.r.t. the weight function exp(-x^2)
  void HermiteH_normalized(unsigned, double, double*);
  //@}
  // ///////////////////////////////////////////////////////////////////////////
} // namespace WDutils {
////////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_WDMath_h
