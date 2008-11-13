// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// gamma.cc                                                                    |
//                                                                             |
// Copyright (C) 1994, 1995, 2004, 2005,2008  Walter Dehnen                    |
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
#include <public/gamma.h>
using namespace falcON;
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// implementing DehnenModel::YcofE();                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace {
  const double tiny=1.e-16;
  double       g,ga,g12,g1,g2,ig2,g3,g4,e,eg2,lq;
  //----------------------------------------------------------------------------
  inline void subyce2(double y, double&f, double&df)
  {
    f  = 0.5*y-std::log(y)-0.5-e;
    df = 0.5-1./y;
  }
  //----------------------------------------------------------------------------
  inline void subyceg(double y, double&f, double&df)
  { 
    register double yg1 = pow(y,g1);
    f  = 0.5*y*yg1*(g2*y-g4) + 1. - eg2;
    df = 0.5*g2*(g3*y-g4)*yg1;
  }
}
//------------------------------------------------------------------------------
double DehnenModel::YcofE(double E) const
{
  if(E==0.) return 1.;
  ::e   = E;
  ::eg2 = e*g2; 
  if((g2>0 && eg2 >1) || E<0.) falcON_Error("DehnenModel: Eps out of range");
  if( g2>0. && eg2==1.) return 0.;
  if(g==1.) return 0.5*(3.-sqrt(8*E+1.));
  if(g==2.) return rtsafe(&subyce2,tiny,1.,eps);
  ::g1 = g1;
  ::g2 = g2;
  ::g3 = g3;
  ::g4 = g4;
  return rtsafe(&subyceg,tiny,1.,eps);
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// implementing DehnenModel::YcofL();                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace {
  inline void subyclg(double y, double&f, double&df)
  { 
    register double yg3 = pow(y,g3);
    f  = (yg3+lq)*y-lq;
    df = g4*yg3+lq;
  }
}
//------------------------------------------------------------------------------
double DehnenModel::YcofLq(double Lq) const
{
  if(Lq==0.) return 0.;
  ::lq = Lq;
  if(g==2.) return 0.5*(sqrt(Lq*(Lq+4))-lq);
  ::g3 = g3;
  ::g4 = g4;
  return rtsafe(&subyclg,0.,1.,eps);
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// implementing DehnenModel::SigIsotropicSquared();                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
double DehnenModel::SigIsoQ(double x) const
{
  if(g>2. && x==0.)
    falcON_Error("DehnenModel: SigIsotropicSquared() diverges at x=0");
  register double x1=x+1;
  if(g==0.)  return (6.*x+1.)/(30.*x1*x1);
  register double y=x/x1;
  if(g==0.5) return 0.2*sqrt(y)/x1;
  if(x>=10.) { // use formula (9) of Tremaine et al. (1994)
    register double xi=1./x, al=2.*g-6., fac=1., xi5=power<5>(xi), sum=0.2*xi5;
    for(register int k=1; k<10; k++)
      {
	al  -= 1.;
	fac *= al/double(k);
	xi5 *= xi;
	sum += fac / (5.+k) * xi5;
      }
    return pow(y,g) * power<4>(x1) * sum;
  }
  if(g==1.) 
    return y * power<4>(x1) * ((((-0.25*y+4./3.)*y-3.)*y+4.)*y-25./12.-log(y));
  if(g==1.5)
    return pow(y,1.5) * power<4>(x1)
      * (((-1./3.*y+2.)*y-6.)*y+1./y+10./3.+4.*log(y));
  if(g==2.) {
    if(x==0.) return 0.5;
    register double y2 = y * y;
    return y2 * power<4>(x1) * (0.5/y2*(1.-y2)*((y-8.)*y+1)-6.*log(y));
  }
  if(g==2.5)
    return power<4>(x1) / sqrt(y)
      * ((((-y+(4.*log(y)-10./3.))*y+6.)*y-2.)*y+1./3.);
  register double g2n=g2-g, ygn=pow(y,g2n), sum;
  sum = (1.-ygn) / g2n;
  sum-= 4.*(1.-(ygn*=y)) / (g2n+=1.);
  sum+= 6.*(1.-(ygn*=y)) / (g2n+=1.);
  sum-= 4.*(1.-(ygn*=y)) / (g2n+=1.);
  sum+= (1.-(ygn*=y)) / (g2n+=1.);
  return pow(y,g) * power<4>(x1) * sum;
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// implementing                                                               //
//   DehnenModel::SurfaceDensity();                                           //
//   DehnenModel::CumSurfaceDensity();                                        //
//   DehnenModel::DSurfaceDensityDR();                                        //
//   DehnenModel::EffectiveRadius();                                          //
//   DehnenModel::SigIsotropicProjected();                                    //
//   DehnenModel::SigCircProjected();                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace {
  double s,sa,ts,fsg;
  //----------------------------------------------------------------------------
  inline double subsur(double t)
  { 
    register double tq = t*t;
    return pow(s+tq,g1) * power<2>(1-tq) / sqrt(ts+sa*tq);
  }
}
//------------------------------------------------------------------------------
double DehnenModel::SurfaceDensity(double R) const
{
  if(g>=1. && R==0.) falcON_Error("DehnenModel: SurfaceDensity() diverges at R=0");
  if(R==0.) return 1 / (Pi * g1 * g2);
  ::s  = R;
  ::sa = 1-s;
  ::ts = s+s;
  ::g1 = g1;
  return g3 / Pi / pow(s+1.,3.5-g) * qbulir(&subsur,0.,1.,eps);
}
//==============================================================================
namespace {
  inline double subdsu(double t)
  { 
    register double tq = t*t;
    return (g4*tq + fsg) * power<4>(1.-tq) / (pow(s+tq,ga) * sqrt(ts+sa*tq));
  }
}
//------------------------------------------------------------------------------
double DehnenModel::DSurfaceDensityDR(double R) const
{
  if(R==0.) falcON_Error("DehnenModel: DSurfaceDensityDR() diverges at R=0");
  ::s   = R;
  ::sa  = 1-s;
  ::ts  = s+s;
  ::fsg = 4*s+g;
  ::ga  = 1+g;
  ::g4  = g4;
  return -g3 * s / Pi * pow(s+1.,g-4.5) * qbulir(&subdsu,0.,1.,eps);
}
//==============================================================================
namespace {
  inline double subcsd(double t)
  {
    register double tq = t*t;
    return pow(s+tq,g1) * tq * sqrt(ts+sa*tq);
  }
}
//------------------------------------------------------------------------------
double DehnenModel::CumSurfaceDensity(double R) const
{
  if(R==0.) return 0.;
  ::s  = R;
  ::sa = 1-s;
  ::ts = s+s;
  ::g1 = g1;
  return 1 - 2 * g3 * pow(s+1,g-2.5) * qbulir(&subcsd,0.,1.,eps);
}
//==============================================================================
namespace {
  double __ep;
  void subxef(double R, double&f, double&df)
  {
    s  = R;
    sa = 1-s;
    ts = s+s;
    f  = 0.5 - 2 * g3 / pow(s+1,2.5-g) * qbulir(&subcsd,0.,1.,__ep);
    df =       2 * g3 / pow(s+1,3.5-g) * qbulir(&subsur,0.,1.,__ep);
  }
}
//------------------------------------------------------------------------------
double DehnenModel::EffectiveRadius() const
{
  ::g = g;
  ::g1= g1;
  ::g3= g3;
  __ep= eps;
  const double a = X<mm>(0.5) * (((-0.00182*g+0.00322)*g-0.00439)*g+0.7549);
  return rtsafe(&subxef,0.9*a,1.1*a,eps);
}
//==============================================================================
namespace {
  inline double subsip(double t)
  {
    register double tq = t*t;
    return pow(s+tq,g12) * power<3>(1-tq) * tq * sqrt(ts+sa*tq);
  }
}
//------------------------------------------------------------------------------
double DehnenModel::SigIsotropicProj(double R) const
{
  if(g>2.  && R==0.)
    falcON_Error("DehnenModel: SigIsotropicProj() diverges at R=0");
  if(g==2. && R==0.) return sqrt(0.5);
  if(g>=1. && R==0.) return 0.;
  if(R==0.) return 0.75 * g1 / (3-2*g) / (5-2*g);
  ::s   = R;
  ::sa  = 1-s;
  ::ts  = s+s;
  ::g12 = 1-g-g;
  return sqrt( qbulir(&subsip,0.,1.,eps) /
	       qbulir(&subsur,0.,1.,eps) / pow(s+1,g2) );
}
//==============================================================================
namespace {
  inline double subscp(double t)
  {
    register double tq = t*t;
    return pow(s+tq,g12) * power<5>(1-tq) / sqrt(ts+sa*tq);
  }
}
//------------------------------------------------------------------------------
double DehnenModel::SigCircProj(double R) const 
{
  if(g>2.  && R==0.) falcON_Error("DehnenModel: SigCircProj() diverges at R=0");
  if(g==2. && R==0.) return 0.5;
  if(R==0.)          return 0.;
  ::s   = R;
  ::sa  = 1-s;
  ::ts  = s+s;
  ::g12 = 1-g-g;
  return s * sqrt( 0.5 * qbulir(&subscp,0.,1.,eps) /
		         qbulir(&subsur,0.,1.,eps) / pow(s+1,g3) );
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// implementing                                                               //
//   DehnenModel::F();                                                        //
//   DehnenModel::G();                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace {
  const double sqrt2 = sqrt(2.);
  const double Pi2   = square(Pi);
  const double Pi3   = cube(Pi);
  double Ag,Bg,Cg,Pg;
  //----------------------------------------------------------------------------
  inline double subyofp(double p)
  {
    return p==0.? 1. : g2==0.? exp(-p) : pow(1-p*g2, ig2);
  }
  //----------------------------------------------------------------------------
  inline double subfei(double t)
  {
    register double y = subyofp(e*(1-t*t));
    return power<2>(1-y) * (Ag+y*(Bg+Cg*y)) / pow(y,Pg);
  }
}
//------------------------------------------------------------------------------
double DehnenModel::Fsub(double E, double G) const
{
  if(E==0.) return 0.;
  ::e   = E;
  ::eg2 = e*g2;
  if(g2>0. && eg2>1.)
    falcON_Error("DehnenModel: Eps=%f > %f=1/(2-g) in DfIsotropic()",E,ig2);
  if(e<0.)
    falcON_Error("DehnenModel: Eps=%f < 0 in DfIsotropic()",E);
  if(g2>0. && eg2==1.)
    falcON_Error("DehnenModel: DfIsotropic() diverging at E=Psi(0)");
  if(G<=0.) G=g;
  double dg = g-G;
  if(dg<0.)
    falcON_Error("DehnenModel: G > gamma in F(E,G)");
  if(g==0. && dg==0.) {
    register double e1=sqrt(e+e), e2=1-e-e, e3=sqrt(e2);
    return 1.5/Pi3*(e1*(2+1./e2)-3.*std::log((1+e1)/e3));
  }
  ::g2 = g2;
  ::ig2= ig2;
  ::Ag = G*(2-dg);
  ::Bg = twice((2-G)*(2-dg)+2*(G-1));
  ::Cg = (4-G)*(2+dg);
  ::Pg = 4-g-dg;
  return (3-G) * 0.25 /(sqrt2*Pi3) *sqrt(e) *qbulir(&subfei,0.,1.,eps,0,0,50);
}
//==============================================================================
namespace {
  double uq;
  //----------------------------------------------------------------------------
  inline double subfom(double t)
  {
    register double y = subyofp(e*(1-t*t));
    return (power<2>(1-y)*(g+y*(2+g4*y))/power<3>(y) + uq*(g4*y-g3))/pow(y,g1);
  }
}
//------------------------------------------------------------------------------
double DehnenModel::F(double Q, double ra) const
{
  if(Q==0.) return 0.;
  ::e   = Q;
  ::eg2 = e*g2;
  if((g2>0. &&eg2>1.) ||e<0.)
    falcON_Error("DehnenModel::F(): Q out of range");
  if(g2>0. && eg2==1.)
    falcON_Error("DehnenModel::F(): diverging at Q=Psi(0)");
  if(ra==0.)
    falcON_Error("DehnenModel::F(): zero anisotropy radius");
  ::uq = 1./(ra*ra);
  ::g  = g;
  ::g1 = g1;
  ::g2 = g2;
  ::ig2= ig2;
  ::g3 = g3;
  ::g4 = g4;
  return g3 * 0.5 / (sqrt2*Pi3) * sqrt(e) * qbulir(&subfom,0.,1.,eps,0,0,50);
}
//==============================================================================
namespace {
  inline double sub__geA(double t)
  {
    // for gamma < 2
    register double tq = t*t;
    register double z  = eg2-tq;
    if(z<=0. || z>=1.) return 0;
    return pow(z,g1) * tq/power<4>(1-pow(z,ig2));
  }
  inline double sub__geB(double t)
  {
    // for gamma > 2
    register double tq = t*t;
    register double y  = pow(eg2*(1-tq),-ig2);
    if(y<=0. || y>=1.) return 0;
    return pow(y,g1) * tq/power<4>(1-y);
  }
  inline double sub__geC(double z)
  {
    // for gamma = 2
    if(z<=0. || z>=1.) return 0.;
    register double z1 = 1-z;
    register double y  = exp(-e-square(z/z1));
    return y*square(y*z/square(z1*(1-y)));
  }
}
//------------------------------------------------------------------------------
double DehnenModel::G(double E) const
{
  ::e   = E;
  ::eg2 = e*g2;
  if(g2>0. && eg2==1.) return 0.;
  if((g2>0. && eg2>1.) || E<0.)
    falcON_Error("DehnenModel::G() E out of range");
  if(E==0.) falcON_Error("DehnenModel::G() diverging at E=0");
  if(g==2.) {
    return 32 * sqrt2 * Pi2 * qbulir(&sub__geC,0.,1.,eps,0,0,50);
  } else if(g<2.) {
    ::g1  = (g+1)*ig2;
    ::eg2 = 1-e*g2;
    ::ig2 = ig2;
    return 32 * sqrt2 * Pi2 * pow(ig2,1.5) * 
      qbulir(&sub__geA,0.,sqrt(eg2),eps,0,0,50);
  } else {
    ::g1  = 6-3*g/2;
    ::eg2 = 1./(1-e*g2);
    ::ig2 = ig2;
    return 32 * sqrt2 * Pi2 * pow(-ig2,1.5) * eg2 *
      qbulir(&sub__geB,0.,1.,eps,0,0,50);
  }
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// implementing class falcON::DehnenModelSampler                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
DehnenModelSampler::
DehnenModelSampler(double gamma,             // I: gamma                  
		   double radius,            // I: scale radius           
		   double Mtot,              // I: GM (untruncated)       
		   double r_a,               //[I: anisotropy radius]     
		   double rmax,              //[I: maximum radius]        
		   int    N,                 //[I: # points on table f(y)]
		   double eps                //[I: numerical precision]   
#ifdef falcON_PROPER
		  ,double __rs,              //[I: mass adaption: scale radius]
		   double __mm,              //[I: mass adaption: mass ratio]
		   double __et,              //[I: mass adaption: shape param]
		   double __nm,              //[I: mass adaption: n_max]
		   bool   __pr               //[I: mass adaption: R_-/Re]
#endif
		   ) :
  ScaledDehnenModel ( gamma,radius,Mtot,eps ),
  SphericalSampler  ( rmax>0? Mr(rmax):Mtot, r_a
#ifdef falcON_PROPER
		      ,__rs,__mm,__et,__nm,__pr
#endif
                    ),
  n                 ( N ),
  y                 ( falcON_NEW(double,n) ),
  f                 ( falcON_NEW(double,n) ),
  fi                ( gamma>0.? 0.5*(gamma-6) : -2. ),
  fo                ( r_a<=0? 2.5 : 0.5 )
{
  if(radius <= 0.0) falcON_Error("DehnenModel: scale radius <= 0\n");
  if(Mtot   <= 0.0) falcON_Error("DehnenModel: total mass <= 0\n");
  if(gamma  <  0.0) falcON_Error("DehnenModel: gamma < 0\n");
  if(gamma  >= 2.5) falcON_Error("DehnenModel: gamma >= 5/2: self-energy=oo\n");
  const double dy=1./double(n+1);
  for(int i=0; i!=n; ++i) { 
    y[i] = (i+1)*dy;
    f[i] = log(r_a>0? F(Psy(y[i]), r_a) : F(Psy(y[i])) );
  }
}
////////////////////////////////////////////////////////////////////////////////
