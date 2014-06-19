// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    src/WDMath.cc                                                      
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    1994-2005                                                          
///                                                                             
/// \todo    add doxygen documentation                                          
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 1994-2005  Walter Dehnen                                       
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
#include <complex>                                  // order important!         
#include <WDMath.h>                                 // have __COMPLEX__ defined 
#include <inline.h>
#include <iostream>

#ifdef __INTEL_COMPILER
#pragma warning (disable:981) /* operands are evaluated in unspecified order */
#endif

using std::exp;
using std::log;
using std::cos;
using std::sin;
using std::complex;

#define MathError(A,B)   WDutils_Error("in %s: %s",B,A)
#define MathWarning(A,B) WDutils_Warning("in %s: %s",B,A)
#define maxit            100
#define fpmin            1.e-40
#define eps              1.e-10

using namespace WDutils;
////////////////////////////////////////////////////////////////////////////////
// volume of unit sphere                                                        
////////////////////////////////////////////////////////////////////////////////
double WDutils::SphVol(int d)
{
  if(d==1) return 2;
  if(d==2) return Pi;
  if(d==3) return FPit;
  int k,n;
  double 
    cn = Pi,
    ae = 2.,
    ao = Pih;
  for(k=1,n=2;;k++) {
    ae *= n/double(n+1);
    cn *= ae;
    if(++n==d) break;
    ao *= n/double(n+1);
    cn *= ao;
    if(++n==d) break;
  }
  return cn;
}
////////////////////////////////////////////////////////////////////////////////
// logarithms of complex trigonometric and hyperbolic functions                 
////////////////////////////////////////////////////////////////////////////////
complex<double> WDutils::lnsin(complex<double> const&x)
{
#ifdef __GNUC__
  double s,c;
  sincos(std::real(x),s,c);
  double 
  ep = exp(-2*abs(std::imag(x))),
  em = s*(1.+ep);
  ep = c*(1.-ep);
#else
  double 
  ep = exp(-2*abs(std::imag(x))),
  em = sin(std::real(x))*(1.+ep);
  ep = cos(std::real(x))*(1.-ep);
#endif
  return complex<double>
    (abs(std::imag(x))+0.5*log(0.25*(em*em+ep*ep)),atan2(sign(std::imag(x))*ep,em));
}
//------------------------------------------------------------------------------
complex<double> WDutils::lncos(complex<double> const&x)
{
#ifdef __GNUC__
  double s,c;
  sincos(std::real(x),s,c);
  double 
  ep = exp(-2*abs(std::imag(x))),
  em = c*(1.+ep);
  ep = s*(1.-ep);
#else
  double 
  ep = exp(-2*abs(std::imag(x))),
  em = cos(std::real(x))*(1.+ep);
  ep = sin(std::real(x))*(1.-ep);
#endif
  return complex<double>
    (abs(std::imag(x))+0.5*log(0.25*(em*em+ep*ep)),atan2(-sign(std::imag(x))*ep,em));
}
//------------------------------------------------------------------------------
complex<double> WDutils::lnsinh(complex<double> const&x)
{
#ifdef __GNUC__
  double s,c;
  sincos(std::imag(x),s,c);
  double 
  ep = exp(-2*abs(std::real(x))),
  em = s*(1.+ep);
  ep = c*(1.-ep);
#else
  double 
  ep = exp(-2*abs(std::real(x))),
  em = sin(std::imag(x))*(1.+ep);
  ep = cos(std::imag(x))*(1.-ep);
#endif
  return complex<double>
    (abs(std::real(x))+0.5*log(0.25*(em*em+ep*ep)),atan2(em,sign(std::real(x))*ep));
}
//------------------------------------------------------------------------------
complex<double> WDutils::lncosh(complex<double> const&x)
{
#ifdef __GNUC__
  double s,c;
  sincos(std::imag(x),s,c);
  double 
  ep = exp(-2*abs(std::real(x))),
  em = c*(1.+ep);
  ep = s*(1.-ep);
#else
  double 
  ep = exp(-2*abs(std::real(x))),
  em = cos(std::imag(x))*(1.+ep);
  ep = sin(std::imag(x))*(1.-ep);
#endif
  return complex<double>(
			 abs(std::real(x))+0.5*log(0.25*(em*em+ep*ep)) ,
			 atan2(sign(std::real(x))*ep,em) );
}
////////////////////////////////////////////////////////////////////////////////
// Gamma functions                                                              
////////////////////////////////////////////////////////////////////////////////
namespace {
  inline double lnGam_pos(double x)  // x >= 0
  {
    const double cof[6]={ 76.18009172947146, -86.50532032941677,
			  24.01409824083091, -1.231739572450155,
			  1.208650973866179e-3, -5.395239384953e-6 };
    double ser=1.000000000190015, y=x, tmp=y+5.5;
    tmp-= (y+0.5) * log(tmp);
    for(int j=0; j<6; j++) ser+= cof[j]/++y;
    return -tmp + log(STPi*ser/x);
  }
  //----------------------------------------------------------------------------
  inline double lnGam(                           // R: Gamma(x)                 
		      double x,                  // I: x                        
		      const char* func=0)        //[I: name of caller]          
  {
    if(x<= 0.) {
      if( is_integral(-x) ) MathError("negative integer argument",func);
      return log(Pi/sin(Pi*x)) - lnGam_pos(1.-x);
    }
    return lnGam_pos(x);
  }
  //----------------------------------------------------------------------------
  inline double lngam_ser(                        // R: ln gamma(a,x) via series
			  double a,               // I: a >= 0                  
			  double x,               // I: x >= 0                  
			  const char  *func = 0)  //[I: name of caller]         
  {
    double
      ap  = a,
      del = 1.0/a,
      sum = del;
    for(int n=1; n<=maxit; n++) {
      ap  += 1.;
      del *= x/ap;
      sum += del;
      if(abs(del) < abs(sum)*eps)
	return log(sum)-x+a*log(x);
    }
    MathError("a too large or maxit too small in lngam_ser()",func);
    return 0.;
  }
  //----------------------------------------------------------------------------
  inline double lnGam_cfr(                        // R: ln Gamma(a,x) via c'd fr
			  double a,         // I: a (can be negative)     
			  double x,         // I: x >= 0                  
			  const char  *func = 0)  //[I: name of caller]         
  {
    double an,del,
      b = x+1.-a,
      c = 1./fpmin,
      d = 1./b,
      h = d;
    for(int i=1; i<=maxit; i++) {
      an =-i*(i-a);
      b += 2.;
      d  = an*d+b; if(abs(d)<fpmin) d=fpmin;
      c  = b+an/c; if(abs(c)<fpmin) c=fpmin;
      d  = 1./d;
      del= d*c;
      h *= del;
      if(abs(del-1.) < eps)
	return log(h)-x+a*log(x);
    }
    MathError("a too large or maxit too small in lnGam_cfr()",func);
    return 0.;
  }
}
//==============================================================================
double WDutils::LogGamma(double x)
{
  return lnGam(x,"LogGamma(x)");
}
//------------------------------------------------------------------------------
complex<double> WDutils::LogGamma(complex<double> const&z)
{
  const double c[6]={ 76.18009172947146, -86.50532032941677,
		      24.01409824083091, -1.231739572450155,
		      1.208650973866179e-3, -5.395239384953e-6 };
  if(iszero(std::imag(z)) && std::real(z) <= 0. && is_integral(-std::real(z))) 
    MathError("z=-n","LogGamma(z)");
  bool turn = std::real(z)<1;
  complex<double> 
    ser = complex<double>(1.000000000190015),
    y   = turn? complex<double>(2.0)-z : z,
    tmp = y+complex<double>(4.5);
  tmp-= (y-complex<double>(0.5))*log(tmp);
  for(int j=0; j<6; j++) {
    ser+= c[j]/y;
    y  += 1.;
  }
  if(turn) {
    y    = Pi*z-Pi;
    tmp -= log(STPi*ser/y) + lnsin(y);
  } else 
    tmp = log(STPi*ser) - tmp;
  while(std::imag(tmp)> Pi) tmp-= complex<double>(0,TPi);
  while(std::imag(tmp)<-Pi) tmp+= complex<double>(0,TPi);
  return tmp;
}
//------------------------------------------------------------------------------
double WDutils::GammaP(double a, double x)
{
  if(x<=0. || a<=0.) MathError("invalid arguments","GammaP(a,x)");
  if(x<(a+1.)) return      exp(lngam_ser(a,x,"GammaP(a,x)") - lnGam_pos(a));
  else         return 1. - exp(lnGam_cfr(a,x,"GammaP(a,x)") - lnGam_pos(a));
}
//------------------------------------------------------------------------------
double WDutils::GammaQ(double a, double x)
{
  if(x<0. || a<=0.) MathError("invalid arguments","GammaQ(a,x)");
  if(x<(a+1.)) return 1. - exp(lngam_ser(a,x,"GammaQ(a,x)") - lnGam_pos(a));
  else         return      exp(lnGam_cfr(a,x,"GammaQ(a,x)") - lnGam_pos(a));
}
//------------------------------------------------------------------------------
double WDutils::Loggamma(double a, double x)
{
  if(x<=0) MathError("x <= 0","Loggamma(a,x)");
  if(a<=0) MathError("a <= 0","Loggamma(a,x)");
  if(x < a+1 )
    return lngam_ser(a,x,"Loggamma(a,x)");
  else
    return log(exp(lnGam_pos(a)) - exp(lnGam_cfr(a,x,"Loggamma(a,x)")));
}
//------------------------------------------------------------------------------
double WDutils::LogGamma(double a, double x)
{
  if(iszero(x)) return lnGam(a,"LogGamma(a,x)");
  if(x < 0 ) MathError("x < 0","LogGamma(a,x)");
  if(x < a+1 && a > 0)
    return log(exp(lnGam_pos(a)) - exp(lngam_ser(a,x,"LogGamma(a,x)")));
  else
    return lnGam_cfr(a,x,"LogGamma(a,x)");
}
////////////////////////////////////////////////////////////////////////////////
// Beta functions                                                               
////////////////////////////////////////////////////////////////////////////////
namespace {
  inline double betacf(double a, double b, double x)
  {
    int    m,m2;
    double aa,del,h,
      qab=a+b,
      qap=a+1,
      qam=a-1,
      c  =1,
      d  =1-qab*x/qap;
    if(abs(d) < fpmin) d=fpmin;
    d = 1/d;
    h = d;
    for(m=1; m<=maxit; ++m) {
      m2 = m+m;
      aa = m*(b-m)*x/((qam+m2)*(a+m2));
      d  = 1+aa*d; if(abs(d)<fpmin) d=fpmin;
      c  = 1+aa/c; if(abs(c)<fpmin) c=fpmin;
      d  = 1/d;
      h *= d*c;
      aa =-(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
      d  = 1+aa*d; if(abs(d)<fpmin) d=fpmin;
      c  = 1+aa/c; if(abs(c)<fpmin) c=fpmin;
      d  = 1/d;
      del= d*c;
      h *= del;
      if(abs(del-1) < eps) return h;
    }
    MathError("a or b too big, or maxit too small","Beta(a,b,x)");
    return 0;
  }
}
//------------------------------------------------------------------------------
double WDutils::LogBeta(double a, double b)
{
  if(a<=0.) MathError("a<=0","LogBeta(a,b)");
  if(b<=0.) MathError("b<=0","LogBeta(a,b)");
  return lnGam_pos(a) +lnGam_pos(b) -lnGam_pos(a+b);
}
//------------------------------------------------------------------------------
double WDutils::Beta(double a, double b)
{
  if(a<=0.) MathError("a<=0","Beta(a,b)");
  if(b<=0.) MathError("b<=0","Beta(a,b)");
  return exp(lnGam_pos(a) +lnGam_pos(b) -lnGam_pos(a+b));
}
//------------------------------------------------------------------------------
double WDutils::Beta(double a, double b, double x)
{
  if(a <= 0) MathError("a <=0","Beta(a,b,x)");
  if(b <= 0) MathError("b <=0","Beta(a,b,x)");
  if(x <  0) MathError("x < 0","Beta(a,b,x)");
  if(x >  1) MathError("x > 1","Beta(a,b,x)");
  if(iszero(x) ) return 0.;
  if(equal(x,1.)) return exp(lnGam_pos(a)+lnGam_pos(b)-lnGam_pos(a+b));
  if(x < (a+1)/(a+b+2))
    return
      exp(a*log(x)+b*log(1-x))*betacf(a,b,x)/a;
  else 
    return
      exp(lnGam_pos(a)+lnGam_pos(b)-lnGam_pos(a+b)) -
      exp(a*log(x)+b*log(1-x))*betacf(b,a,1-x)/b;
}
//------------------------------------------------------------------------------
WDutils::BetaFunc::BetaFunc(double __a, double __b)
  : a ( __a ), 
    b ( __b ), 
    B ( exp(lnGam_pos(a)+lnGam_pos(b)-lnGam_pos(a+b)) ),
    x0( (a+1)/(a+b+2) ) {}
double WDutils::BetaFunc::operator() (double x) const {
  if(x <  0) MathError("x < 0","BetaFunc(x)");
  if(x >  1) MathError("x > 1","BetaFunc(x)");
  if(iszero(x) ) return 0;
  if(equal(x,1.)) return B;
  return x < x0 ?
        exp(a*log(x)+b*log(1-x))*betacf(a,b,  x)/a : 
    B - exp(a*log(x)+b*log(1-x))*betacf(b,a,1-x)/b ;
}
////////////////////////////////////////////////////////////////////////////////
// Exponential integrals                                                        
////////////////////////////////////////////////////////////////////////////////
double WDutils::En(int n, double x)
{
  if(n<0 || x<0. || (iszero(x) && n<=1)) MathError("bad argumends","En()");
  if(n==0)  return exp(-x)/x;
  if(iszero(x)) return 1./double(n-1);
  double ans;
  if(x>1.) {
    int    i,nm1=n-1;
    double a,b,c,d,del,h;
    b = x+n;
    c = 1./fpmin;
    d = 1./b;
    h = d;
    for(i=1; i<=maxit; i++) {
      a   =-i*(nm1+i);
      b  += 2.;
      d   = 1./(a*d+b);
      c   = b+a/c;
      del = c*d;
      h  *= del;
      if(abs(del-1.) < eps) return h*exp(-x);
    }
    ans = h*exp(-x);
    MathWarning("continued fraction failed","En()");
  } else {
    int    i,ii,nm1=n-1;
    double del,fac,psi;
    ans = nm1? 1./double(nm1) : -log(x)-EulerGamma();
    fac = 1.;
    for(i=1; i<=maxit; i++) {
      fac *=-x/double(i);
      if(i!=nm1)
	del =-fac/double(i-nm1);
      else {
	psi =-EulerGamma();
	for(ii=1; ii<=nm1; ii++)
	  psi+= 1./double(ii);
	del = fac*(psi-log(x));
      }
      ans += del;
      if(abs(del) < abs(ans)*eps) return ans;
    }
    MathWarning("series failed","En()");
  }
  return ans;
}
//------------------------------------------------------------------------------
double WDutils::Ei(double x)
{
  if(x<=0.)   return -En(1,-x);
  if(x<fpmin) return log(x)+EulerGamma();
  int    k;
  const double logeps =-20.72326583694641115616192309216;
  double fact=1.,sum=0.,term=1.;
  if(x<=-logeps) {
    for(k=1; k<=maxit; k++) {
      fact*= x/k;
      term = fact/k;
      sum += term;
      if(term<eps*sum) break;
    }
    if(k>maxit) MathError("series failed","Ei()");
    return sum+log(x)+EulerGamma();
  }
  for(k=1; k<=maxit; k++) {
    fact = term;
    term*= k/x;
    if(term<eps) break;
    if(term<fact) sum+= term;
    else {
      sum -= fact;
      break;
    }
  }
  if(k>maxit) MathError("series failed","Ei()");
  return exp(x)*(1.0+sum)/x;
}
////////////////////////////////////////////////////////////////////////////////
// Bessel functions                                                             
////////////////////////////////////////////////////////////////////////////////
double WDutils::J0(double x)
{
  double ax=abs(x),y,ans1,ans2;
  if(ax < 8.) {
    y    = x*x;
    ans1 = 57568490574.0+y*(-13362590354.0+y*(651619640.7
	   +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
    ans2 = 57568490411.0+y*(1029532985.0+y*(9494680.718
	   +y*(59272.64853+y*(267.8532712+y*1.0))));
    return ans1/ans2;
  } else {
    double z=8./ax, xx=ax-0.785398164;
    y   =z*z;
    ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
         +y*(-0.2073370639e-5+y*0.2093887211e-6)));
    ans2=       -0.1562499995e-1+y*(0.1430488765e-3
         +y*(-0.6911147651e-5+y*(0.7621095161e-6-y*0.934935152e-7)));
    return sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
  }
}
//------------------------------------------------------------------------------
double WDutils::J1(double x)
{
  double ax=abs(x),y,ans1,ans2;
  if(ax < 8.) {
    y    = x*x;
    ans1 = x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
	   +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
    ans2 = 144725228442.0+y*(2300535178.0+y*(18583304.74
	   +y*(99447.43394+y*(376.9991397+y*1.0))));
    return ans1/ans2;
  } else {
    double z=8./ax, xx=ax-2.356194491;
    y    = z*z;
    ans1 = 1.0+y*(0.183105e-2+y*(-0.3516396496e-4
	   +y*(0.2457520174e-5+y*(-0.240337019e-6))));
    ans2 =        0.04687499995+y*(-0.2002690873e-3
	   +y*(0.8449199096e-5+y*(-0.88228987e-6 +y*0.105787412e-6)));
    return sign(x) * sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
  }
}
//------------------------------------------------------------------------------
double WDutils::Jn(unsigned n, double x)
{
  if(n==0)  return J0(x);
  if(n==1)  return J1(x);

  const double acc=60., bigno=1.e10, bigni=1.e-10;
  unsigned m;
  bool jsum;
  double ax,bj,bjm,bjp,sum,tox,ans;

  ax=abs(x);
  if(iszero(ax)) return 0.;
  if(ax>double(n)) {
    tox = 2./ax;
    bjm = J0(ax);
    bj  = J1(ax);
    for(unsigned j=1; j!=n; ++j) {
      bjp=j*tox*bj-bjm;
      bjm=bj;
      bj=bjp;
    }
    ans=bj;
  } else {
    tox  = 2./ax;
    m    = 2*((n+unsigned(sqrt(acc*n))/2));
    jsum = false;
    bjp  = ans = sum = 0.;
    bj   = 1.;
    for(unsigned j=m; j; --j) {
      bjm = j*tox*bj-bjp;
      bjp = bj;
      bj  = bjm;
      if(abs(bj) > bigno) {
	bj  *= bigni;
	bjp *= bigni;
	ans *= bigni;
	sum *= bigni;
      }	
      if(jsum) sum += bj;
      jsum=!jsum;
      if(j==n) ans=bjp;
    }
    sum=2.0*sum-bj;
    ans /= sum;
  }
  return x < 0.0 && (n & 1) ? -ans : ans;
}
//------------------------------------------------------------------------------
double WDutils::Y0(double x)
{
  if(x<0.) MathError("negative argument","Y0(x)");
  double y,ans1,ans2;
  if(x < 8.0) {
    y    = x*x;
    ans1 =-2957821389.0+y*(7062834065.0+y*(-512359803.6
	   +y*(10879881.29+y*(-86327.92757+y*228.4622733))));
    ans2 = 40076544269.0+y*(745249964.8+y*(7189466.438
	   +y*(47447.26470+y*(226.1030244+y*1.0))));
    return (ans1/ans2)+0.636619772*J0(x)*log(x);
  } else {
    double z=8./x, xx=x-0.785398164;
    y    = z*z;
    ans1 = 1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
	   +y*(-0.2073370639e-5+y*0.2093887211e-6)));
    ans2 =        -0.1562499995e-1+y*(0.1430488765e-3
	   +y*(-0.6911147651e-5+y*(0.7621095161e-6
	   +y*(-0.934945152e-7))));
    return sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
  }
}
//------------------------------------------------------------------------------
double WDutils::Y1(double x)
{
  if(x<0.) MathError("negative argument","Y1(x)");
  double y,ans1,ans2;
  if(x < 8.) {
    y    = x*x;
    ans1 = x*(-0.4900604943e13+y*(0.1275274390e13
	   +y*(-0.5153438139e11+y*(0.7349264551e9
	   +y*(-0.4237922726e7+y*0.8511937935e4)))));
    ans2 = 0.2499580570e14+y*(0.4244419664e12
	   +y*(0.3733650367e10+y*(0.2245904002e8
	   +y*(0.1020426050e6+y*(0.3549632885e3+y)))));
    return (ans1/ans2)+0.636619772*(J1(x)*log(x)-1.0/x);
  } else {
    double z=8./x, xx=x-2.356194491;
    y    = z*z;
    ans1 = 1.0+y*(0.183105e-2+y*(-0.3516396496e-4
	   +y*(0.2457520174e-5+y*(-0.240337019e-6))));
    ans2 =    0.04687499995+y*(-0.2002690873e-3
	   +y*(0.8449199096e-5+y*(-0.88228987e-6+y*0.105787412e-6)));
    return sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
  }
}
//------------------------------------------------------------------------------
double WDutils::Yn(unsigned n, double x)
{
  if(x<0.) MathError("negative argument","Yn(x)");
  if(n==0) return Y0(x);
  if(n==1) return Y1(x);
  double by=Y1(x),bym=Y0(x),byp,tox=2./x;
  for(unsigned j=1; j!=n; ++j) {
    byp = j*tox*by-bym;
    bym = by;
    by  = byp;
  }
  return by;
}
//------------------------------------------------------------------------------
double WDutils::I0(double x)
{
  double ax=abs(x),y;
  if(ax < 3.75) {
    y = x/3.75;
    y*= y;
    return 1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
		  +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
  } else {
    y = 3.75/ax;
    return (exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
	    +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
	    +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
	    +y*0.392377e-2))))))));
  }
}
//------------------------------------------------------------------------------
double WDutils::I1(double x)
{
  double ans,ax=abs(x),y;
  if(ax < 3.75) {
    y = x/3.75;
    y*= y;
    ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
	+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
  } else {
    y=3.75/ax;
    ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1-y*0.420059e-2));
    ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
		      +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
    ans *= (exp(ax)/sqrt(ax));
  }
  return x < 0.0 ? -ans : ans;
}
//------------------------------------------------------------------------------
double WDutils::In(unsigned n, double x)
{
  const double acc=60., bigno=1.e10, bigni=1.e-10;
  if(n==0)  return I0(x);
  if(n==1)  return I1(x);
  if(iszero(x)) return 0.;
  double bi,bim,bip,tox,ans;
  tox=2.0/fabs(x);
  bip=ans=0.0;
  bi=1.0;
  for(unsigned j=2*(n+unsigned(sqrt(acc*n))); j; --j) {
    bim = bip+j*tox*bi;
    bip = bi;
    bi  = bim;
    if(abs(bi) > bigno) {
      ans *= bigni;
      bi  *= bigni;
      bip *= bigni;
    }
    if(j==n) ans=bip;
  }
  ans *= I0(x)/bi;
  return x < 0.0 && (n & 1) ? -ans : ans;
}
//------------------------------------------------------------------------------
double WDutils::K0(double x)
{
  if(x<0.) MathError("negative argument","K0(x)");
  double y;
  if(x <= 2.) {
    y = x*x/4.;
    return (-log(x/2.0)*I0(x))+(-0.57721566+y*(0.42278420
	    +y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
	    +y*(0.10750e-3+y*0.74e-5))))));
  } else {
    y = 2./x;
    return (exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
	    +y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
	    +y*(-0.251540e-2+y*0.53208e-3))))));
  }
}
//------------------------------------------------------------------------------
double WDutils::K1(double x)
{
  if(x<0.) MathError("negative argument","K1(x)");
  double y;
  if(x <= 2.) {
    y=x*x/4.0;
    return (log(x/2.0)*I1(x))+(1.0/x)*(1.0+y*(0.15443144
	    +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
	    +y*(-0.110404e-2+y*(-0.4686e-4)))))));
  } else {
    y = 2./x;
    return (exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
	    +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
	    +y*(0.325614e-2+y*(-0.68245e-3)))))));
  }
}
//------------------------------------------------------------------------------
double WDutils::Kn(unsigned n, double x)
{
  if(x<0.) MathError("negative argument","Kn(x)");
  if(n==0) return K0(x);
  if(n==1) return K1(x);
  double bk,bkm,bkp,tox;
  tox = 2./x;
  bkm = K0(x);
  bk  = K1(x);
  for(unsigned j=1; j!=n; ++j) {
    bkp = bkm+j*tox*bk;
    bkm = bk;
    bk  = bkp;
  }
  return bk;
}
////////////////////////////////////////////////////////////////////////////////
// Hermite Polynomials                                                          
////////////////////////////////////////////////////////////////////////////////
double WDutils::HermiteH(unsigned n, double x)
{
  if(n==0) return 1;
  if(n==1) return 2.*x;
  double h0=1., h1=2*x, hi=h1;
  for(unsigned i=1; i!=n; ++i) {
    hi = 2. * (x*h1 - i*h0);
    h0 = h1;
    h1 = hi;
  }
  return hi;
}
//------------------------------------------------------------------------------
void WDutils::HermiteH(unsigned n, double x, double *H)
{
  H[0] = 1.; 	if(n==0) return;
  H[1] = 2*x;	if(n==1) return;
  for(unsigned i=1; i!=n; ++i) H[i+1] = 2*(x*H[i]-2*H[i-1]);
}
//------------------------------------------------------------------------------
void WDutils::NormSqHermite(unsigned n, double *N)
{
  N[0] = SPi;	if(n==0) return;   	                // Sqrt[Pi]             
  N[1] = 2*SPi;	if(n==1) return;	                // 2*Sqrt[Pi]           
  for(unsigned i=2; i<=n; ++i) N[i] = 2*i*N[i-1];
}
//------------------------------------------------------------------------------
double WDutils::HermiteH_normalized(unsigned n, double x)
{
  if(n==0) return 1.   / sqrt(Pi);
  if(n==1) return 2.*x / sqrt(TPi);
  double   h0=1., h1=2*x, hi=h1;
  unsigned N=2;
  for(unsigned i=1; i<n; ++i, N*=2*i) {
    hi = 2. * (x*h1 - i*h0);
    h0 = h1;
    h1 = hi;
  }
  return hi / sqrt(N*Pi);
}
//------------------------------------------------------------------------------
void WDutils::HermiteH_normalized(unsigned n, double x, double *H)
{
  H[0] = 1.;
  if(n>0)  H[1] = 2*x;
  for(unsigned i=1; i<n; ++i)
    H[i+1] = 2 * (x*H[i] - 2*H[i-1]);
  for(unsigned i=0,N=1; i<=n; ++i,N*=2*i) 
    H[i] /= sqrt(N*Pi);
}
////////////////////////////////////////////////////////////////////////////////
