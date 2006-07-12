// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// king.cc                                                                     |
//                                                                             |
// Copyright (C) 2000, 2001, 2002, 2003  Walter Dehnen                         |
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
#include <public/king.h>
#include <public/basic.h>
#include <numerics.h>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace falcON;
//==============================================================================
// falcON::king_model::set_up() & auxiliary data & routines                     
//==============================================================================
namespace {
  double c00;                                           // used by set_up() etc 
  typedef tupel<2,double> vec2D;
  //----------------------------------------------------------------------------
  inline double rhok(const double Psi)
    // B&T eq. 4.131 with rho_1 = 1 and sigma = 1
  {
    const double f1  =1.1283791670955125739, f2=2./3.;  // sqrt(4/Pi)           
    if(Psi<=0.) return 0.;
    register double sPsi=sqrt(Psi);
    return exp(Psi)*erf(sPsi) - f1*sPsi*(1+f2*Psi);
  }
  //--------------------------------------------------------------------------
  inline vec2D integrand1(double const&r, vec2D const&Y)
    // takes x=r as independent, Y[0]=Psi, Y[1]=Psi'
  {
    register vec2D F;
    F[0] = Y[1];
    F[1] = (r>0.)? -FPi*rhok(Y[0])-2*Y[1]/r : -c00;
    return F;
  }
  //----------------------------------------------------------------------------
  inline vec2D integrand2(double const&Psi, vec2D const&Y)
    // takes Psi as independent, Y[0]=r, Y[1]=Psi'
  {
    register vec2D F;
    F[0] = 1./Y[1];
    F[1] = -FPi*rhok(Psi)/Y[1] -2/Y[0];
    return F;
  }
}
//------------------------------------------------------------------------------
void king_model::setup(const unsigned n)
{
  if(N != n) {
    if(r)  falcON_DEL_A(r);
    if(ps) falcON_DEL_A(ps);
    if(m)  falcON_DEL_A(m);
    if(rh) falcON_DEL_A(rh);
    N  = n;
    r  = falcON_NEW(double,N);
    ps = falcON_NEW(double,N);
    m  = falcON_NEW(double,N);
    rh = falcON_NEW(double,N);
  }
  r [0] = 0.;
  ps[0] = Psi0;
  m [0] = 0.;
  rh[0] = rhok(Psi0);
  rh0   = rh[0];
  r0    = 1.5/sqrt(Pi*rh0);            // r_0 = sqrt(9/4PiG rho_0)
  c00   = FPi*rh0/3.;                  // 4PiG rho_0/3

  register unsigned        i,Ni=N/2;
  register double          dx=2.*r0/double(Ni-1);
  register vec2D y;
  y[0] = ps[0];
  y[1] = 0.;
  for(i=1; i<Ni; i++) {
    y     = rk4(y,r[i-1],dx,&integrand1);
    r [i] = r[i-1] + dx;
    ps[i] = y[0];
    m [i] =-y[1]*r[i]*r[i];
    rh[i] = rhok(ps[i]);
  }
  y[0] = r [i-1];
  dx   =-ps[i-1] / double(N-i);
  for(; i<N; i++) {
    y     = rk4(y,ps[i-1],dx,&integrand2);
    r [i] = y[0];
    ps[i] = ps[i-1] + dx;
    m [i] =-y[1]*r[i]*r[i];
    rh[i] = rhok(ps[i]);
  }
  c  = log10(r[N-1]/r0);
  rh[N-1] = 0.0;
  P0 =-m[N-1]/r[N-1];
}
//==============================================================================
// falcON::king_model: reset                                                    
//==============================================================================
void king_model::reset()
{
  N = 0;
  if(r)  { falcON_DEL_A(r);  r =0; }
  if(ps) { falcON_DEL_A(ps); ps=0; }
  if(m)  { falcON_DEL_A(m);  m =0; }
  if(rh) { falcON_DEL_A(rh); rh=0; }
}
//==============================================================================
// falcON::king_model: re-scaling                                               
//==============================================================================
void king_model::reset_scales_tidal(const double M, const double Rt)
{
  rscal = Rt/r[N-1];
  Pscal = M/(m[N-1]*rscal);
  vscal = sqrt(Pscal);
  rhscl = Pscal/(rscal*rscal);
  drscl = rhscl/rscal;
  mscal = Pscal*rscal;
  sdscl = Pscal/rscal;
}
//------------------------------------------------------------------------------
void king_model::reset_scales_core(const double M, const double R0)
{
  rscal = R0/r0;
  Pscal = M/(m[N-1]*rscal);
  vscal = sqrt(Pscal);
  rhscl = Pscal/(rscal*rscal);
  drscl = rhscl/rscal;
  mscal = Pscal*rscal;
  sdscl = Pscal/rscal;
}
//==============================================================================
// falcON::king_model::rms_radius()                                             
//==============================================================================
double king_model::rms_radius() const
{
  register double dM,y=0.;
  for(register unsigned i=0,j=1; j!=N; ++i,++j) {
    dM  = m[j]-m[i];
    y  += (square(r[j])+square(r[i]))*dM;
  }
  return rscal*sqrt(0.5*y/m[N-1]);
}
//==============================================================================
// falcON::king_model::half_mass_radius()                                       
//==============================================================================
double king_model::half_mass_radius() const
{
  return rscal*polev(0.5*m[N-1],m,r,N);
}
//==============================================================================
// falcON::king_model::Etot()                                                   
//==============================================================================
double king_model::Etot() const
{
  register double dM, psi,e=0.;
  for(register unsigned i=0,j=1; j!=N; ++i,++j) {
    dM  = m[j]-m[i];
    psi = ps[j]+ps[i];
    e  += dM*psi;
  }
  return -0.25*Pscal*mscal*(square(m[N-1])/r[N-1] + e*0.5);
}
//==============================================================================
// falcON::king_model::write_table()                                            
//==============================================================================
using std::ios;
using std::setw;
using std::setprecision;
using std::ofstream;
//------------------------------------------------------------------------------
void king_model::write_table(const char* file) const
{
  ofstream table;
  if(! open(table,file) )  return;
  table.setf(ios::left, ios::adjustfield);
  table<<"#        r     Phi(r)     rho(r)       M(r)\n";
  for(register unsigned i=0; i<N; i++)
    table<<setw(10)<<setprecision(6)<<rscal*r[i]      <<" "  // r     
	 <<setw(10)<<setprecision(6)<<Pscal*(P0-ps[i])<<" "  // Phi(r)
	 <<setw(10)<<setprecision(6)<<rhscl*rh[i]     <<" "  // rho(r)
	 <<setw(10)<<setprecision(6)<<mscal*m[i]      <<"\n";// M(r)  
  table.close();
}
//==============================================================================
// falcON::king_model::random() & auxiliary data & methods                      
//==============================================================================
namespace {
  double psi,rhi;
  //----------------------------------------------------------------------------
  inline double feps(const double Eps) {
    // B&T eq. 4.130 with rho_1 = 1 and sigma = 1                               
    const double fac = 0.063493635934240969785;    // (2*Pi)^(-3/2)             
    if(Eps <= 0.) return 0.;
    return fac * (exp(-Eps)-1.);
  }
  //----------------------------------------------------------------------------
  inline void cumdens(const double v, double& P, double& p) {
    // given v and r (via psi(r)), we return mass(<v|r)-rho and                 
    // its derivative p(v|r) = 4Pi f(v|r) v^2                                   
    const double sqh = 0.70710678118654752440,
                 sq2p= 0.79788456080286535588;     // sqrt(2/Pi)                
    register double vq=v*v;
    P = exp(psi)*erf(v*sqh)-sq2p*(vq*v/3.+v*exp(psi-0.5*vq)) - rhi;
    p = FPi * feps(psi-0.5*vq) * vq;
  }
}
//------------------------------------------------------------------------------
double king_model::random(const double R1, const double R2,
			  double &rad, double &vel) const {
  // 1. find radius from cumulative mass                                        
  register double mi = R1*m[N-1];             // mass(<rad)                     
  rad = polev(mi,m,r, N);                     // corresponding radius           
  // 2. find velocity from cumulative mass at given radius                      
  psi = polev(mi,m,ps,N);                     // psi(rad)                       
  rhi = R2*rhok(psi);                         // mass(<vel|rad)                 
  vel = rtsafe(&cumdens,0.,sqrt(2*psi),1.e-8);
  // 3. scale them to appropriate units                                         
  rad*= rscal;
  vel*= vscal;
  // 4. return potential                                                        
  return (P0-psi)*Pscal;
}
//==============================================================================
// falcON::king_model::SurfDens() & auxiliary data & methods                    
//==============================================================================
namespace {
  double *rad, *rho, Rq, rt;
  int     K;
  //----------------------------------------------------------------------------
  inline double sd_integrand(const double z)
  {
    register double radius = sqrt(Rq+z*z);
    return (radius >= rt)? 0.0 : polev(radius,rad,rho,K);
  }
}
//------------------------------------------------------------------------------
double king_model::SurfDens(double R) const
{
  R    /= rscal;
  ::rt  = r[N-1];
  ::K   = N;
  if(R >= rt) return 0.;
  ::Rq  = R*R;
  ::rad = r;
  ::rho = rh;
  double S = qbulir(sd_integrand,0.0,sqrt(rt*rt-Rq),1.e-6);
  return 2*Pscal/rscal*S;
}
////////////////////////////////////////////////////////////////////////////////
