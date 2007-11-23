// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
/// \file src/public/lib/halo.cc                                               |
//                                                                             |
// Copyright (C) 2000-2007  Walter Dehnen, Paul McMillan                       |
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
#include <public/halo.h>
#include <public/basic.h>
#include <utils/spline.h>
#include <utils/WDMath.h>
using namespace falcON;
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::HaloModifier                                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace {
  /// z(y) = 2/(y+1/y) = 2*y/(1+y^2)
  /// \return z(y)
  /// \param y y
  inline double z(double y) {
    return (y+y)/(1+y*y);
  }
  /// z(y) = 2/(y+1/y) = 2*y/(1+y^2)
  /// \param z1 dz(y)/dy
  inline double zd(double y, double&z1) {
    double yq=y*y, yq1=1/(1+yq);
    z1 = twice(1-yq)*yq1*yq1;
    return (y+y)*yq1;
  }
  /// z(y) = 2/(y+1/y) = 2*y/(1+y^2)
  /// \param z1 dz(y)/dy
  /// \param z2 d^2z(y)/dy^2
  inline double zdd(double y, double&z1, double&z2) {
    double yq=y*y, yq1=1/(1+yq), yq2=yq1*yq1;
    z1 = twice(1-yq)*yq2;
    z2 = 4*y*(yq-3)*yq2*yq1;
    return (y+y)*yq1;
  }
  /// z(y) = 2/(y+1/y) = 2*y/(1+y^2)
  /// \param y1 (input) dy/dx
  /// \param z1 dz/dx
  inline double z(double y, double y1, double&z1) {
    double r=zd(y,z1);
    z1 *= y1;
    return r;
  }
  /// z(y) = 2/(y+1/y) = 2*y/(1+y^2)
  /// \param y1 (input) dy/dx
  /// \param y2 (input) d^2y/dx^2
  /// \param z1 dz/dx
  /// \param z2 d^2z/dx^2
  inline double z(double y,
		  double y1, double y2,
		  double&z1, double&z2) {
    double r=zdd(y,z1,z2);
    z2 *= y1*y1;
    z2 += z1*y2;
    z1 *= y1;
    return r;
  }
} // namespace {
//------------------------------------------------------------------------------
inline double HaloModifier::trunc(double r) const {
  if(sechtr)
    return ::z(exp(-r*irt));
  else
    return ::z(::z(exp(-r*irt)));
}
//------------------------------------------------------------------------------
inline double HaloModifier::trunc(double r, double&t1) const {
  double y=exp(-r*irt), y1=-irt*y;
  if(sechtr)
    return ::z(y,y1,t1);
  else {
    double z1,z=::z(y,y1,z1);
    return ::z(z,z1,t1);
  }
}
//------------------------------------------------------------------------------
inline double HaloModifier::trunc(double r, double&t1, double&t2) const {
  double y=exp(-r*irt), y1=-irt*y, y2=-irt*y1;
  if(sechtr)
    return ::z(y,y1,y2,t1,t2);
  else {
    double z1,z2,z=::z(y,y1,y2,z1,z2);
    return ::z(z,z1,z2,t1,t2);
  }
}
//------------------------------------------------------------------------------
inline double HaloModifier::core(double r) const {
  return sqrt(r*r+rcq);
}
//------------------------------------------------------------------------------
inline double HaloModifier::core(double r, double&t1) const {
  double c=sqrt(r*r+rcq);
  t1 = r/c;
  return c;
}
//------------------------------------------------------------------------------
inline double HaloModifier::core(double r, double&t1, double&t2) const {
  double c=sqrt(r*r+rcq),ic=1./c;
  t1 = r*ic;
  t2 = rcq*ic*ic*ic;
  return c;
}
//------------------------------------------------------------------------------
HaloModifier::HaloModifier(double c, double t) falcON_THROWING
: rc(abs(c)), rcq(c*c),
  rt(abs(t)), irt(rt? 1/rt : 0.), sechtr(t>=0)
{
  if(isinf(t)) falcON_THROW("HaloModifier: truncation radius == inf\n");
  if(isnan(t)) falcON_THROW("HaloModifier: truncation radius == nan\n");
  if(c    <0.) warning("HaloModifier: core radius = %g<0; will use %g\n",c,rc);
//   std::cerr<<" Testing HaloModifier::trunc(): (irt="
// 	   <<irt<<", sechtr="<<sechtr<<")\n";
//   for(;;) {
//     double r;
//     std::cout<<" r="; std::cin>>r;
//     if(r<0) break;
//     double dr=0.0001*r;
//     double rl=r-dr, rh=r+dr;
//     double t =trunc(r),tl=trunc(rl),th=trunc(rh);
//     double t1,t1l,t1h,t1r,t2;
//     double tr=trunc(r,t1),trl=trunc(rl,t1l),trh=trunc(rh,t1h);
//     double trr=trunc(r,t1r,t2);
//     std::cerr<<" tr="<<t<<" ="<<tr<<" ="<<trr<<'\n'
// 	     <<" t1="<<t1<<" ="<<t1r
// 	     <<" ="<<((th-tl)/(dr+dr))<<'\n'
// 	     <<" t2="<<t2
// 	     <<" ="<<((t+t-tl-th)/(dr*dr))
// 	     <<" ="<<((t1h-t1l)/(dr+dr))<<'\n';
//   }
}
//------------------------------------------------------------------------------
inline double HaloModifier::cored(HaloDensity const&Model,
				  double r) const
{
  return rcq? Model(core(r)) : Model(r);
}
//------------------------------------------------------------------------------
inline double HaloModifier::cored(HaloDensity const&Model,
				  double r, double&rh1) const
{
  if(rcq==0.) return Model(r,rh1);
  double dx1,x =core (r,dx1);
  double dr1,rh=Model(x,dr1);
  rh1 = dr1*dx1;
  return rh;
}
//------------------------------------------------------------------------------
inline double HaloModifier::cored(HaloDensity const&Model,
				  double r, double&rh1, double&rh2) const
{
  if(rcq==0.) return Model(r,rh1,rh2);
  double dx1,dx2,x =core (r,dx1,dx2);
  double dr1,dr2,rh=Model(x,dr1,dr2);
  rh1 = dr1*dx1;
  rh2 = dr2*dx1*dx1 + dr1*dx2;
  return rh;
}
//------------------------------------------------------------------------------
double HaloModifier::operator()(HaloDensity const&Model,
				double r) const
{
  return irt? trunc(r)*cored(Model,r) : cored(Model,r);
}
//------------------------------------------------------------------------------
double HaloModifier::operator()(HaloDensity const&Model,
				double r, double&rh1) const
{
  if(!irt) return cored(Model,r,rh1);
  double dt1,tr=trunc(r,dt1);
  double dr1,rh=cored(Model,r,dr1);
  rh1 = dr1*tr + rh*dt1;
  return rh*tr;
}
//------------------------------------------------------------------------------
double HaloModifier::operator()(HaloDensity const&Model,
				double r, double&rh1, double&rh2) const
{
  if(!irt) return cored(Model,r,rh1,rh2);
  double dt1,dt2,tr=trunc(r,dt1,dt2);
  double dr1,dr2,rh=cored(Model,r,dr1,dr2);
  rh1 = dr1*tr + rh*dt1;
  rh2 = dr2*tr + twice(dr1*dt1) + rh*dt2;
  return rh*tr;
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::DoublePowerLawHalo                                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
DoublePowerLawHalo::DoublePowerLawHalo(double inner, double outer, double trans)
  : gi(inner), go(outer), et(trans), gg(go-gi), al(gg/et)
{
  if(gi < 0.)
    falcON_THROW("DoublePowerHalo: inner power-law slope < 0\n");
  if(gi > go)
    falcON_THROW("DoublePowerHalo: inner power-law slope > outer\n");
  if(et <= 0.)
    falcON_THROW("DoublePowerHalo: transition steepness <= 0\n");
}
//------------------------------------------------------------------------------
double DoublePowerLawHalo::operator()(double x) const {
  return pow(x,-gi) * pow(1.+pow(x,et),-al);
}
//------------------------------------------------------------------------------
double DoublePowerLawHalo::operator()(double x, double&rh1) const {
  double q = pow(x,et), q1=1./(1+q);
  double f = pow(x,-gi) * pow(q1,al);
  double g = gi + gg*q*q1;
  rh1 = -f*g/x;
  return f;
}
//------------------------------------------------------------------------------
double DoublePowerLawHalo::operator()(double x, double&rh1, double&rh2) const {
  double q = pow(x,et), q1=1./(1+q);
  double f = pow(x,-gi) * pow(q1,al);
  double g = gg*q*q1;
  double gx= et*q1*g;
  g  += gi;
  rh1 = f/x;
  rh2 = rh1*(g*(1+g)-gx)/x;
  rh1*=-g;
  return f;
}
//------------------------------------------------------------------------------
namespace {
  const HaloModifier* __HM;
  double __z0,__iE,__Ai,__Ao,__Rc,__Yo;
  inline double __dM(double z)
  {
    if(z<=__z0) return 0.;
    double z1 = 1-z;
    if(z1<=0) return __Yo;
    double u = std::pow(z/z1,__iE);
    if(u<=__Rc) return 0.;
    double r = __Rc? sqrt(u*u-__Rc*__Rc) : u;
    double y = (r/u) * std::pow(z,__Ai) * std::pow(z1,__Ao);
    return __HM->truncated()? __HM->trunc(r) * y : y;
  }
}
double DoublePowerLawHalo::Mtot(const HaloModifier&hm) const
  falcON_THROWING
{
  if(hm.cored() || hm.truncated()) {
    if(go < 3+et && !hm.truncated())
      falcON_THROW("DoublePowerLawHalo::Mtot(): cannot compute total mass"
		   "for outer<3+eta and no truncation\n");
    __HM = &hm;
    __iE = 1./et;
    __Ai = __iE*(3-gi)-1;
    __Ao = __iE*(go-3)-1;
    __Yo = __Ao > 1.e-7 || __HM->truncated()? 0. : 1.;
    __Rc = hm.r_c();
    double tm = std::pow(__Rc,et);
    __z0 = tm/(1.+tm);
    return FPi*qbulir(&__dM,__z0,1.,1.e-9)/et;
  } else {
    if(gi >= 3. || go <= 3.)
      falcON_THROW("DoublePowerLawHalo::Mtot(): total mass diverges\n");
    return FPi*Beta((3-gi)/et,(go-3)/et)/et;
  }
}
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// auxiliary stuff for class HaloPotential                                      
//                                                                              
////////////////////////////////////////////////////////////////////////////////
namespace {
  /// for the maximum table radius in case of trunction: gamma(Rmax)=gam_trunc
  const double gam_trunc = 100.;
  /// for the table extrema in case of no trunction: 
  ///    |gamma(r)-gamma_asymptotic| < eps_gamma * gamma_asymptotic
  const double eps_gamma = 0.002;
  //////////////////////////////////////////////////////////////////////////////
  /// find an appropriate maximum radius for the tables in HaloModel
  double Rmax(HaloDensity const&halo)
  {
    double r = halo.scale_radius();
    double gam,rh=halo(r,gam);
    gam *= -r/rh;
    if(halo.trunc_radius() > 0.) {
      while(gam < gam_trunc) {
	r  += r;
	rh  = halo(r,gam);
	gam*=-r/rh;
      }
    } else {
      const double gamo = halo.outer_gamma();
      const double egmo = eps_gamma;
      while(abs(gam-gamo) > egmo) {
	r  += r;
	rh  = halo(r,gam);
	gam*=-r/rh;
      }
    }
    return r;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// find an appropriate minimum radius for the tables in HaloModel
  double Rmin(HaloDensity const&halo)
  {
    double r = halo.scale_radius();
    double gam,rh=halo(r,gam);
    gam *= -r/rh;
    const double gam0 = halo.inner_gamma();
    const double egm0 = eps_gamma;
    while(abs(gam-gam0) > egm0) {
      r  *= 0.5;
      rh  = halo(r,gam);
      gam*=-r/rh;
    }
    return r;
  }
  //////////////////////////////////////////////////////////////////////////////
  const HaloDensity *RHO;
  spline<double>    *SPLINE;
  //----------------------------------------------------------------------------
  /// give dM/dln r for HaloDensity pointed to by RHO
  inline double dM(double lr, double const&M) {
    const double r = exp(lr);
    return cube(r)*(*RHO)(r);
  }
  //----------------------------------------------------------------------------
  /// give dPsi_h / dln r, assuming SPLINE holds M(ln r)
  inline double dP(double lr, double const&P) {
    return -(*SPLINE)(lr) * exp(-lr);              // dPsi = -M/r * dlnr        
  }
} // namespace {
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class falcON::HaloPotential                                                  
//                                                                              
////////////////////////////////////////////////////////////////////////////////
HaloPotential::HaloPotential(HaloDensity const&model,
			     const acceleration*mono)
  : DEN(model), MON(mono),
    Ah(model.inner_gamma()), At(Ah), Ch(model.outer_gamma()), PS(0)
{
  // 0. compute min/max radius for tables
  const double rmin = Rmin(model);
  const double rmax = Rmax(model);
  debug_info(2,"HaloPotential: minimum & maximum tabulated radii = %g, %g\n",
	     rmin,rmax);

  // 1. set r,lr,mh
  n1 = int(200*log10(rmax/rmin));
  n  = 1+n1;
  const double dlr= log(rmax/rmin)/double(n1);
  r .reset(n);
  lr.reset(n);
  mh.reset(n);
  {
    RHO =&DEN;
    double M = DEN(rmin)*cube(rmin)/(3-Ah);     // M_h(<rmin)/4Pi
    r [0] = rmin;
    lr[0] = log(rmin);
    mh[0] = FPi*M;
    for(int i=1; i!=n; ++i) {
      lr[i] = lr[0] + i * dlr;                  // grid is linear in ln r
      r [i] = exp(lr[i]);                       // set r_i
      M     = rk4(M,lr[i-1],dlr,dM);            // integrate mass outwards
      mh[i] = FPi*M;                            // set cumulative halo mass
    }
    debug_info(4,"HaloPotential: M_halo(<%g)=%g\n",r[n1],mh[n1]);
  }
  // 2. set mt,rh,ps,ec
  // 2.1 set effects of external potential
  if(!MON) {
    ps.reset(n,0.);
    mt.reset(n,0.);
    rh.reset(n,0.);
  } else {
    ps.reset(n);
    mt.reset(n);
    rh.reset(n);
    //    find mass and potential of external monopole
    Array<double> pot_e(n);
    Array<vect_d> pos_e(n), acc_e(n);
    for(int i=0; i!=n; ++i) 
      pos_e[i] = vect_d(r[i],0.,0.);
    MON->set(0.,n,0,pos_e.array(),0,0,
	     pot_e.array(),acc_e.array(),0);    // get external monopole
    for(int i=0; i!=n; ++i) {
      ps[i] =-pot_e[i];                         // get Psi_e(r)
      mt[i] =-square(r[i])*acc_e[i][0];         // get M_e(<r)
      if(i && ps[i]>ps[i-1])
	error("HaloPotential: external Ps(%g)=%g > Ps(%g)=%g\n",
	      r[i],ps[i],r[i-1],ps[i-1]);
    }
    //    find density generating external monopole
    double extA  = (log(mt[1])-log(mt[0]))/dlr; // gamma_0 for monopole
    if(extA > At) At = extA;                    // gamma_0 for potential
    double dmdlr = mt[0] * extA;                // dM_e/dlnr assuming power law
    spline<double> Sext(lr,mt,&dmdlr);
    for(int i=0; i!=n; ++i) {
      Sext(lr[i], &(rh[i]));                    // G*rho = d(G*M)/dlnr
      rh[i] /= FPi * cube(r[i]);                //       / (4 Pi r^3)
    }
    debug_info(4,"HaloPotential: M_mono(<%g)=%g\n",r[n1],mt[n1]);
  }
  // 2.2 make mt hold the total cumulative mass and rh the total density
  for(int i=0; i!=n; ++i) {
    mt[i] += mh[i];
    rh[i] += DEN(r[i]);
  }
  // 2.3 find total potential & Ec: add psi_halo to ps[]; get ec
  ec.reset(n);
  {
    double P = mh[n1] / r[n1],                  // Psi_h(r_max)
      s0h = FPi * cube(r[ 0]) *  DEN(r[ 0]),
      sNh = FPi * cube(r[n1]) *  DEN(r[n1]);
    SPLINE = new spline<double>(lr,mh,&s0h,&sNh);
    ps[n1] += P;
    for(int i=n-2; i!=-1; --i) {
      P     = rk4(P,lr[i+1],-dlr,dP);
      ps[i]+= P;
      ec[i] = ps[i] - 0.5*mt[i]/r[i];
    }
    falcON_DEL_O(SPLINE);
  }
  ps0 = At<2? ps[0] + mt[0]/(r[0]*(2-At)) : 0.;
  // 2.4 find smallest index beyond which the total mass does not seem to change
  for(nm=0; nm!=n1; ++nm) if(mt[nm+1] == mt[nm]) break;
  nm++;
  // 3   make penta spline for potential and force
  dp.reset(n);
  for(int i=0; i!=n; ++i)
    dp[i] = -mt[i]/r[i];
  PS = new Pspline<double>(lr,ps,dp);
}
//------------------------------------------------------------------------------
// destructor
HaloPotential::~HaloPotential() {
  if(PS) falcON_DEL_O(PS);
}
//------------------------------------------------------------------------------
// potential Phi(r) and (-dPhi/dr)/r using penta spline
double HaloPotential::PotAcc(double Rq, double&A) const {
  double P, lR=0.5*log(Rq);
  if     (Rq <= 0.) { A=0.; P=ps0; }
  else if(lR<lr[0]) {
    if     (At< 2.) { P = (ps[0]-ps0)*exp((2-At)*(lR-lr[0]));
                      A = (2-At)*P/Rq;
		      P+= ps0; }
    else if(At==2.) { P = ps[0]*(lR-lr[0]);
                      A = ps[0]/Rq; }
    else            { P = ps[0]*exp((2-At)*(lR-lr[0]));
                      A = (2-At)*P/Rq; }
  }
  else if(lR>lr[n1]){ P = mt[n1]/sqrt(Rq);
                      A =-P/Rq; }
  else              { P = (*PS)(lR,&A);
                      A/= Rq; }
  return -P;
}
//------------------------------------------------------------------------------
// potential Psi(r)
double HaloPotential::Ps(double R) const {
  if(R<= 0.)   return ps0;
  if(R< r[0] ) {
    if(At< 2.) return ps0 - (ps0-ps[0]) * pow(R/r[0],2-At);
    if(At==2.) return ps[0] * (log(R)-lr[0]);
    else       return ps[0] * pow(R/r[0],2-At);
  }
  if(R> r[n1]) return mt[n1]/R;
  else         return polev(log(R),lr,ps);
}
//------------------------------------------------------------------------------
// log R_psi(E)
double HaloPotential::lnRPsi(double P) const {
  if(P> ps[0]) {
    if(At< 2.) return lr[0]+log((ps0-P)/(ps0-ps[0]))/(2-At);
    if(At==2.) return lr[0]-P/ps[0];
    else       return lr[0]+log(P/ps[0])/(2-At);
  }
  if(P<ps[n1]) return log(mt[n1]/P);
  else         return polev(P,ps,lr);
}
//------------------------------------------------------------------------------
// R_psi(E)
double HaloPotential::RPsi(double P) const {
  if(P> ps[0]) {
    if(At< 2.) return r[0]*pow((ps0-P)/(ps0-ps[0]),1/(2-At));
    if(At==2.) return r[0]*exp(-P/ps[0]);
    else       return r[0]*pow(P/ps[0],1/(2-At));
  }
  if(P<ps[n1]) return mt[n1]/P;
  else         return exp(polev(P,ps,lr));
}
//------------------------------------------------------------------------------
// total cumulative mass
double HaloPotential::Mt(double R) const {
  if(R<= 0.)  return 0.;
  if(R<r[0] ) return mt[0]*pow(R/r[0],3-At);
  if(R>r[n1]) return mt[n1];
  else        return polev(log(R),lr,mt);
}
//------------------------------------------------------------------------------
// cumulative halo mass
double HaloPotential::Mh(double R) const {
  if(R<= 0.)  return 0.;
  if(R<r[ 0]) return mh[0]*pow(R/r[0],3-Ah);
  if(R>r[n1]) return mh[n1];
  else        return polev(log(R),lr,mh);
}
//------------------------------------------------------------------------------
// total mass density
double HaloPotential::rhot(double R) const {
  if(R<= 0.)  return 0.;
  if(R<r[0] ) return rh[0]*pow(R/r[0],-At);
  if(R>r[n1]) return 0.;
  else        return polev(log(R),lr,rh);
}
//------------------------------------------------------------------------------
// v_circ^2(r)
double HaloPotential::vcq(double R) const {
  if(R<= 0.) return 0.;
  return Mt(R)/R;
}
//------------------------------------------------------------------------------
// Omega^2(r)
double HaloPotential::omq(double R) const {
  if(R <= 0.) return 0.;
  return Mt(R)/cube(R);
}
//------------------------------------------------------------------------------
// kappa^2(r)
double HaloPotential::kpq(double R) const {
  if(R <= 0.) return 0.;
  return omq(R) + FPi * rhot(R);
}
//------------------------------------------------------------------------------
// gamma := 2*Omega/kappa 
double HaloPotential::gam(double R) const {
  if(R < r[ 0]) return 2./sqrt(4.-At);
  if(R > r[n1]) return 1.;
  double oq = Mt(R)/cube(R);
  return 2.*sqrt(oq/(oq + FPi*rhot(R)));
}
//------------------------------------------------------------------------------
// Eps_c(r)
double HaloPotential::Epc(double R) const {
  if(R<=0.)    return ps0;
  if(R< r[ 0]) return Ps(R) - 0.5*vcq(R);
  if(R> r[n1]) return 0.5*mt[n1]/R;
  else         return polev(log(R),lr,ec);
}
//------------------------------------------------------------------------------
// R_circ(E)
double HaloPotential::RcE(double E) const {
  if(At>=2 && E>=ps0) return 0.;
  if(E> ec[0]) {
    if(At< 2.)  return r[0]*pow((ps0-E)/(ps0-ec[0]),1/(2-At));
    if(At==2.)  return r[0]*exp(0.5*(E/ec[0]-1));
    else        return r[0]*pow(E/ec[0],1/(2-At));
  }
  if(E< ec[n1]) return 0.5*mt[n1]/E;
  else          return exp(polev(E,ec,lr));
}
//------------------------------------------------------------------------------
// estimate for R(E, L^2, cos[eta])
double HaloPotential::Rap(double E, double Lq, double ce) const {
  const double
    Rc  = RcE(E),
    ecc = sqrt(1-Lq/(Mt(Rc)*Rc));
  return Rc * pow(1+ecc*ce, 0.5*gam(Rc));
}
//------------------------------------------------------------------------------
// estimate for R_peri(E, L^2)
double HaloPotential::Rp(double E, double Lq) const {
  return Rap(E,Lq,-1.);
}
//------------------------------------------------------------------------------
// estimate for R_apo(E, L^2)
double HaloPotential::Ra(double E, double Lq) const {
  return Rap(E,Lq,1.);
}
//------------------------------------------------------------------------------
// radius given halo mass
double HaloPotential::RMh(double M) const {
  if(M<= 0.)   return 0.;
  if(M<= mh[ 0]) return r[0]*pow(M/mh[0],1/(3-Ah));
  if(M>  mh[n1]) error("HaloPotential::rMh(): M>M_halo(oo)\n");
  return exp(polev(M,mh,lr));
}
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// auxiliary stuff for class HaloModel                                          
//                                                                              
////////////////////////////////////////////////////////////////////////////////
namespace {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class ReducedDensity                                                       
  //                                                                            
  // to be used with complete DF : L^-2b g(Q=Eps-u^2*L^2/2)                     
  //                                                                            
  // returns reduced density  f(r) * rho(r)  with reduction factor              
  //                                                                            
  //                          1-b   2b                                          
  //      f(r) = (1 + r^2 u^2)     r                                            
  //                                            b                               
  //           = (1 + r^2 u^2) (r^2/[1+r^2 u^2])                                
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class ReducedDensity : public HaloDensity {
    const HaloDensity&Model;
    const double uq, b, tb;
    //--------------------------------------------------------------------------
    double reduc(double r) const {
      if(b == 0.) return 1+r*r*uq;
      if(uq== 0.) return pow(r,tb);
      const double rq=r*r, ft=1+rq*uq;
      return ft * pow(rq/ft,b);
    }
    //--------------------------------------------------------------------------
    double reduc(double r, double &d1) const {
      if(b == 0.) {
	double t=r*uq;
	d1 = t+t;
	return 1+r*t;
      }
      if(uq == 0.) {
	const double p=(tb==1.)? 1. : pow(r,tb-1);
	d1 = tb*p;
	return r*p;
      }
      const double rq=r*r, qq=rq*uq, ft=1+qq;
      const double fc=ft * pow(rq/ft,b);
      d1 = twice(b+qq)*fc/(r*ft);
      return fc;
    }
    //--------------------------------------------------------------------------
    double reduc(double r, double&d1, double&d2) const {
      if(b == 0.) {
	double t=r*uq;
	d1 = t+t;
	d2 = uq+uq;
	return 1+r*t;
      }
      if(uq== 0.) {
	const double p=(tb==1.)? 1. : pow(r,tb-1);
	d1 = tb*p;
	d2 = d1*(tb-1)/r;
	return r*p;
      }
      const double rq=r*r, qq=rq*uq, ft=1+qq, ift=1/ft, al=twice(b+qq)*ift;
      const double fc=ft * pow(rq*ift,b);
      d1 = fc/r;
      d2 = d1/r*(al*(al-1)+4*qq*(1-b)*ift*ift);
      d1*= al;
      return fc;
    }
    //--------------------------------------------------------------------------
  public:
    ReducedDensity(const HaloDensity&m, double r_a, double beta)
      : Model(m), uq(r_a>0? 1/square(r_a):0), b(beta), tb(b+b) {}
    //--------------------------------------------------------------------------
    double inner_gamma() const {
      return Model.inner_gamma() - tb;
    }
    //--------------------------------------------------------------------------
    double scale_radius() const {
      return Model.scale_radius();
    }
    //--------------------------------------------------------------------------
    double trunc_radius() const {
      return Model.trunc_radius();
    }
    //--------------------------------------------------------------------------
    double outer_gamma() const {
      return Model.outer_gamma() - (uq? 2 : tb);
    }
    //--------------------------------------------------------------------------
    double operator() (double r) const {
      return reduc(r)*Model(r);
    }
    //--------------------------------------------------------------------------
    double operator() (double r, double &d1) const {
      double r1,re=reduc(r,r1);
      double m1,mo=Model(r,m1);
      d1 = re*m1 + r1*mo;
      return re*mo;
    }
    //--------------------------------------------------------------------------
    double operator()(double r, double&d1, double&d2) const {
      double r1,r2,re=reduc(r,r1,r2);
      double m1,m2,mo=Model(r,m1,m2);
      d1 = re*m1 + r1*mo;
      d2 = re*m2 + twice(r1*m1) + r2*mo;
      return re*mo;
    }
  };
  //----------------------------------------------------------------------------
  //  Clever way of integrating the equation for the distribution function      
  //  N.B. the part for when psi < last point on spline IS used, and works      
  //  for all density profiles provided ReducedDensity is done right.           
  //----------------------------------------------------------------------------
  const ReducedDensity  *RED;
  double nu,MT,Be,Q;                               // 1/(1-p), Mtot, beta, Q    
  inline double intgQ_of_q(double q) {
    //                                                                          
    //   d Psi       Q^(1-p)                                                    
    // ----------- = ------- dq    with  Psi = Q (1-q^[1/(1-p)])                
    // (Q - Psi)^p    1-p                                                       
    //                                                                          
    // Q and p are both constants in this integration, making life very easy    
    const double psi = Q * (1-pow(q,nu));          // Psi = Q * (1 - q^nu)      
    if     (psi <= 0.)
      return 0.;
    else if(psi < SPLINE->last_X()) {
      if(Be > 0.5) {
	double r=MT/psi,rd1;
	(*RED)(r,rd1);
	return -r*rd1/psi;
      }
      if(Be >-0.5){
	double r=MT/psi,rd1,rd2;
	(*RED)(r,rd1,rd2);
 	return r/(psi*psi)*(twice(rd1) + r*rd2);
      }
    } else
      return (*SPLINE)(psi);
  }
  //----------------------------------------------------------------------------
  inline double intgQ(double z) {                  // q = z^4; dq = 4 z^3 dz    
    const double f=cube(z);
    return 4*f*intgQ_of_q(z*f);
  }
}
////////////////////////////////////////////////////////////////////////////////
falcON_TRAITS(ReducedDensity,"ReducedDensity");
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class falcON::HaloModel                                                      
//                                                                              
////////////////////////////////////////////////////////////////////////////////
HaloModel::HaloModel(HaloDensity const&model,
		     double beta, double r_a,
		     const acceleration*mono)
  : HaloPotential(model,mono), B(beta), RA(r_a)
{
  // 0. perform some sanity checks
  if(Ah < B+B) falcON_THROW("HaloModel: 2*beta > gamma_0: unphysical model\n");
  if(B >  1.0) falcON_THROW("HaloModel: beta > 1: unphysical model\n");
  if(B < -0.5) falcON_THROW("HaloModel: beta < -1/2 not supported\n");
  // 1 compute distribution function                                           
  // 1.1 tabulate integrand on grid considering possible cases for beta
  Array<double,1> in(n);
  RED = new ReducedDensity(DEN,B,RA);
  if       (B >= 0.5) {        //  0.5 <= beta <  1  :  tabulate d red/ dpsi    
    bool integrand_negative = false;
    for(int i=0; i!=n; ++i) {
      double rd1,
	ps1 = -mt[i]/square(r[i]);                // psi'                       
      (*RED)(r[i],rd1);                           // red'                       
      in[i] = rd1/ps1;                            // dred/dpsi                  
      if(in[i]<0.) integrand_negative = true;
    }
    if(integrand_negative)
      warning("HaloModel: integrand for DF is negative;\n");
  } else if(B >=-0.5) {        // -0.5 <= beta <  0.5:  tabulate d^2 red/ dpsi^2
    bool integrand_negative = false;
    for(int i=0; i!=n; ++i) {
      double rd1,rd2,
	ps1 =-mt[i]/square(r[i]),                 // psi'                       
	ps2 =-twice(ps1)/r[i]-FPi*(rh[i]);        // psi"                       
      (*RED)(r[i],rd1,rd2),                       // red', red"                 
      in[i] = (rd2*ps1 - rd1*ps2)/cube(ps1);      // d^2red/dpsi^2              
      if(in[i]<0.) integrand_negative = true;
    }
    if(integrand_negative)
      warning("HaloModel: integrand for DF is negative;\n");
  } else
    error("HaloModel: beta < -0.5 not supported\n");
  // 1.2 compute g(Q) of f=L^-2B * g(Q)                                         
  lg.reset(n);
  const double
    Logalfa = (1.5-B) * LogofTwo + LogBeta(0.5,1-B) - LogofPi +
              log( (B >= 0.5 ? 1 : ( 0.5-B)) *
	           (B >=-0.5 ? 1 : (-0.5-B)) );
  if(B == 0.5  ||  B ==-0.5 ||  B ==-1.5)     // trivial cases: no integration  
    for(int i=0; i!=n; ++i)
      lg[i] = log(in[i]) - Logalfa;
  else {
    const int mm = int(1.5-B);
    const double
      p    = 1.5-B-mm,
      p1   = 1-p,
      lfc  = log(sin(p1*Pi)) - LogofPi - Logalfa;
    nu     = 1/p1;
    Be     = B;
    MT     = mt[n1];
    SPLINE = new spline<double>(ps,in);
    for(int i=0; i!=n; ++i) {
      Q = ps[i];
      double err, g = qbulir(intgQ,0.,1.,1.e-8,&err,0,50);
      if(g<0.) error("HaloModel: g(Q=%g)=%g < 0, err=%g\n",Q,g,err);
      if(err>1.e-3) 
	warning("HaloModel: inaccurate integration for g(Q) at Q=%g\n",Q);
      lg[i] = lfc + p1*log(Q) + log(nu*g);
      if(i && lg[i]>lg[i-1])
	warning("HaloModel: non-monotinic DF at E=%g\n",ps[i]);
    }
    delete SPLINE;
  }
  falcON_DEL_O(RED); RED=0;
}
//------------------------------------------------------------------------------
// g(Q)
double HaloModel::lnG(double Q) const {
  if(Q < ps[n1])
    // 1 at small Q (large radii)
    if(RA) // 1.1 finite anisotropy radius
      return lg[n1] + (3.5-Ch-B) * (lnRPsi(Q)-lr[n1]);
    else   // 1.2 infinite anisotropy radius (recognised as RA=0)
      return lg[n1] + (1.5-Ch+B) * (lnRPsi(Q)-lr[n1]);
  else if(Q > ps[0])
    // 2 at large Q (small radii)
    return lg[0] + (B+B-Ah-(1.5-B)*(2-At)) * (lnRPsi(Q)-lr[0]);
  else
    // 3 withing tabulated interval
    return polev(Q,ps,lg);
}
//------------------------------------------------------------------------------
// distribution function f(Eps,L^2)
double HaloModel::fEL (double E, double Lq) const {
  return pow(Lq,-B) * exp(lnG(RA? E-0.5*Lq/(RA*RA) : E));
}
//------------------------------------------------------------------------------
