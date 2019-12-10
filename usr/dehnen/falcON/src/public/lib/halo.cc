// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
//
/// \file src/public/lib/halo.cc
//
// Copyright (C) 2000-2004 Walter Dehnen
// Copyright (C) 2005-2007 Walter Dehnen, Paul McMillan
// Copyright (C) 2008-2011 Walter Dehnen
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
#include <public/halo.h>
#include <public/basic.h>
#include <utils/spline.h>
#include <utils/WDMath.h>
#include <utils/timer.h>
#include <cassert>
#include <strings.h>

using namespace falcON;
//
// class falcON::HaloModifier
//
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
//
inline double HaloModifier::trunc(double r, double&t1) const {
  double y=exp(-r*irt), y1=-irt*y;
  if(sechtr)
    return ::z(y,y1,t1);
  else {
    double z1,z=::z(y,y1,z1);
    return ::z(z,z1,t1);
  }
}
//
inline double HaloModifier::trunc(double r, double&t1, double&t2) const {
  double y=exp(-r*irt), y1=-irt*y, y2=-irt*y1;
  if(sechtr)
    return ::z(y,y1,y2,t1,t2);
  else {
    double z1,z2,z=::z(y,y1,y2,z1,z2);
    return ::z(z,z1,z2,t1,t2);
  }
}
//
inline double HaloModifier::core(double r) const {
  return sqrt(r*r+rcq);
}
//
inline double HaloModifier::core(double r, double&t1) const {
  double c=sqrt(r*r+rcq);
  t1 = r/c;
  return c;
}
//
inline double HaloModifier::core(double r, double&t1, double&t2) const {
  double c=sqrt(r*r+rcq),ic=1./c;
  t1 = r*ic;
  t2 = rcq*ic*ic*ic;
  return c;
}
//
HaloModifier::HaloModifier(double c, double t) falcON_THROWING
: rc(abs(c)), rcq(c*c),
  rt(abs(t)), irt(rt? 1/rt : 0.), sechtr(t>=0)
{
  if(std::isinf(t)) falcON_THROW("HaloModifier: truncation radius == inf");
  if(std::isnan(t)) falcON_THROW("HaloModifier: truncation radius == nan");
  if(c<0.) falcON_Warning("HaloModifier: core radius = %g<0; will use %g\n",
			  c,rc);
}
//
inline double HaloModifier::cored(HaloDensity const&Model,
				  double r) const
{
  return rcq? Model(core(r)) : Model(r);
}
//
inline double HaloModifier::cored(HaloDensity const&Model,
				  double r, double&rh1) const
{
  if(rcq==0.) return Model(r,rh1);
  double dx1,x =core (r,dx1);
  double dr1,rh=Model(x,dr1);
  rh1 = dr1*dx1;
  return rh;
}
//
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
//
double HaloModifier::operator()(HaloDensity const&Model,
				double r) const
{
  return irt? trunc(r)*cored(Model,r) : cored(Model,r);
}
//
double HaloModifier::operator()(HaloDensity const&Model,
				double r, double&rh1) const
{
  if(!irt) return cored(Model,r,rh1);
  double dt1,tr=trunc(r,dt1);
  double dr1,rh=cored(Model,r,dr1);
  rh1 = dr1*tr + rh*dt1;
  return rh*tr;
}
//
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
//
// class falcON::DoublePowerLawHalo
//
DoublePowerLawHalo::Model DoublePowerLawHalo::model(const char*mod, bool warn)
{
  if(mod==0) return Default;
  if(0==strcasecmp(mod,"Plummer"))   return Plummer;
  if(0==strcasecmp(mod,"Jaffe"))     return Jaffe;
  if(0==strcasecmp(mod,"Hernquist")) return Hernquist;
  if(0==strcasecmp(mod,"Dehnen"))    return Dehnen;
  if(0==strcasecmp(mod,"Zhao"))      return Default;
  if(0==strcasecmp(mod,"NFW"))       return NFW;
  if(0==strcasecmp(mod,"Moore"))     return Moore;
  if(0==strcasecmp(mod,"DM"))        return DM;
  if(warn) falcON_WarningN("DoublePowerLawHalo: "
			   "unknown model '%s'; assume default\n",mod);
  return Default;
}
//
const char* DoublePowerLawHalo::name(DoublePowerLawHalo::Model model)
{
  switch(model){
  case Default:   return "Zhao";
  case Plummer:   return "Plummer";
  case Jaffe:     return "Jaffe";
  case Hernquist: return "Hernquist";
  case Dehnen:    return "Dehnen";
  case NFW:       return "NFW";
  case Moore:     return "Moore";
  case DM:        return "DM";
  default:        return "unknown model";
  }
}
//
namespace {
  using namespace falcON;
  const double req = DoublePowerLawHalo::null_value();
  //
  const double Value[8][3] = {
    { req, req, req },         // Default:       gi=???  go=???   et=???
    { 0.0, 5.0, 2.0 },         // Plummer:       gi=0    go=5     et=2
    { 2.0, 4.0, 1.0 },         // Jaffe:         gi=2    go=4     et=1
    { 1.0, 4.0, 1.0 },         // Hernquist:     gi=1    go=4     et=1
    { req, 4.0, 1.0 },         // Dehnen:        gi=???  go=4     et=1
    { 1.0, 3.0, 1.0 },         // NFW:           gi=1    go=3     et=1
    { 1.5, 3.0, 1.5 },         // Moore:         gi=3/2  go=3     et=3/2
    { 0.777777777777777778,    // DM:            gi=7/9  go=31/9  et=4/9
      3.444444444444444444,
      0.444444444444444444 } };
  const char* Name[3] = {"inner", "outer", "eta"};
  //
  inline double value(DoublePowerLawHalo::Model m, double x, int v)
    falcON_THROWING
  {
    if(Value[m][v] == req) {
      if(x == req)
	falcON_THROWN("parameter '%s' required for model '%s'\n",
		      Name[v],DoublePowerLawHalo::name(m));
      return x;
    } else {
      if(x != req && x != Value[m][v])
	falcON_WarningN("%s=%g ignored, using %s=%g for model '%s'\n",
			Name[v],x,Name[v],Value[m][v],
			DoublePowerLawHalo::name(m));
      return Value[m][v];
    }
  }
}
//
double DoublePowerLawHalo::inner_value(Model m, double x) falcON_THROWING
{ return ::value(m,x,0); }
double DoublePowerLawHalo::outer_value(Model m, double x) falcON_THROWING
{ return ::value(m,x,1); }
double DoublePowerLawHalo::trans_value(Model m, double x) falcON_THROWING
{ return ::value(m,x,2); }
//
double DoublePowerLawHalo::inner_default(Model m)
{ return ::Value[m][0]; }
double DoublePowerLawHalo::outer_default(Model m)
{ return ::Value[m][1]; }
double DoublePowerLawHalo::trans_default(Model m)
{ return ::Value[m][2]; }
//
DoublePowerLawHalo::DoublePowerLawHalo(Model mod,
				       double inner, double outer, double trans)
  falcON_THROWING :
  go(outer_value(mod,outer)),
  gi(inner_value(mod,inner)),
  et(trans_value(mod,trans)),
  gg(go-gi), al(gg/et)
{
  if(gi < 0.)
    falcON_THROW("DoublePowerHalo: inner power-law slope = %g < 0",gi);
  if(gi > go)
    falcON_THROW("DoublePowerHalo: inner power-law slope = %g > outer = %g",
		 gi,go);
  if(et <= 0.)
    falcON_THROW("DoublePowerHalo: transition steepness = %g <= 0",et);
}
//
DoublePowerLawHalo::DoublePowerLawHalo(double inner, double outer, double trans)
  : go(outer), gi(inner), et(trans), gg(go-gi), al(gg/et)
{
  if(gi < 0.)
    falcON_THROW("DoublePowerHalo: inner power-law slope = %g < 0",gi);
  if(gi > go)
    falcON_THROW("DoublePowerHalo: inner power-law slope = %g > outer = %g",
		 gi,go);
  if(et <= 0.)
    falcON_THROW("DoublePowerHalo: transition steepness = %g <= 0",et);
}
//
double DoublePowerLawHalo::operator()(double x) const {
  return pow(x,-gi) * pow(1.+pow(x,et),-al);
}
//
double DoublePowerLawHalo::operator()(double x, double&rh1) const {
  double q = pow(x,et), q1=1./(1+q);
  double f = pow(x,-gi) * pow(q1,al);
  double g = gi + gg*q*q1;
  rh1 = -f*g/x;
  return f;
}
//
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
//
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
		   "for outer<3+eta and no truncation");
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
      falcON_THROW("DoublePowerLawHalo::Mtot(): total mass diverges");
    return FPi*Beta((3-gi)/et,(go-3)/et)/et;
  }
}
//
// class HaloPotential
//
namespace {
  /// for the maximum table radius in case of truncation: gamma(Rmax)=gam_trunc
  const double gam_trunc = 100.;
  /// for the table extrema in case of no truncation: 
  ///    |gamma(r)-gamma_asymptotic| < eps_gamma * gamma_asymptotic
  const double eps_gamma = 0.01;
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
    return min(r,0.01*halo.scale_radius());
  }
  //
  const HaloDensity *RHO;
  spline<double>    *SPLINE;
  /// give dM/dln r for HaloDensity pointed to by RHO
  inline double dM(double lr, double const&) {
    const double r = exp(lr);
    return cube(r)*(*RHO)(r);
  }
  /// give dPsi_h / dln r, assuming SPLINE holds M(ln r)
  inline double dP(double lr, double const&) {
    return -(*SPLINE)(lr) * exp(-lr);              // dPsi = -M/r * dlnr        
  }
} // namespace {
//
HaloPotential::HaloPotential(HaloDensity const&model,
			     const acceleration*mono,
			     double r_max) falcON_THROWING
  : DEN(model), MON(mono),
    Ah(model.inner_gamma()), At(Ah), Ch(model.outer_gamma()), PS(0)
{
  // 0. compute min/max radius for tables
  const double rmin = Rmin(model);
  const double rmax = r_max? max(r_max,Rmax(model)) : Rmax(model);
  DebugInfo(2,"HaloPotential: minimum & maximum tabulated radii = %g, %g\n",
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
    DebugInfo(4,"HaloPotential: M_halo(<%g)=%g\n",r[n1],mh[n1]);
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
    WDutils::Timer Tim;
    MON->set(0.,n,0,pos_e.array(),0,0,
	     pot_e.array(),acc_e.array(),0);    // get external monopole
    DebugInfo(4,"HaloPotential: call to external potential for %d particles "
	      " took %f seconds",n,Tim.stop());
    for(int i=0; i!=n; ++i) {
      ps[i] =-pot_e[i];                         // get Psi_e(r)
      mt[i] =-square(r[i])*acc_e[i][0];         // get M_e(<r)
      if(!(ps[i]>0))
	falcON_THROW("HaloPotential: external Ps(%g)=%20.16g",
	      r[i],ps[i]);
      if(!(mt[i]>0))
	falcON_THROW("HaloPotential: external M(%g)=%20.16g",
		     r[i],mt[i]);
      if(i && !(ps[i]<=ps[i-1]))
	falcON_THROW("HaloPotential: external Ps(%g)=%20.16g > Ps(%g)=%20.16g",
		     r[i],ps[i],r[i-1],ps[i-1]);	
    }
    //    find density generating external monopole
    double extA  = (log(mt[1])-log(mt[0]))/dlr; // dlnM/dlnr for monopole
    if(3-extA > At) At = 3-extA;                // gamma_0 for potential
    double dmdlr = mt[0] * extA;                // dM_e/dlnr assuming power law
    spline<double> Sext(lr,mt,&dmdlr);
    int pjt = 0;                                // warn if pjt becomes non-0
    for(int i=0; i!=n; ++i) {
      Sext(lr[i], &(rh[i]));                    // G*rho = d(G*M)/dlnr
      rh[i] /= FPi * cube(r[i]);                //       / (4 Pi r^3)
      DebugInfo(9,"PJT: %d %g %g\n", i, r[i], rh[i]);
#if 0      
      if(!(rh[i]>0))
	falcON_THROW("HaloPotential: external Rh(%g)=%20.16g i=%d n=%d",
		     r[i],rh[i],i,n);
#else
      if (rh[i] < 0) {                          // dunno why this is now
	if (pjt==0) pjt=i;                      // needed; float/double?
	rh[i] = 0.0;
      }
#endif      
    }
    if (pjt) falcON_Warning("HaloModel: PJT patched at %d/%d\n",pjt,n);
    DebugInfo(4,"HaloPotential: M_mono(<%g)=%g\n",r[n1],mt[n1]);
  }
  // 2.2 make mt hold the total cumulative mass and rh the total density
  for(int i=0; i!=n; ++i) {
    mt[i] += mh[i];
    assert(mt[i]>0);
    assert(i==0 || mt[i]>=mt[i-1]);
    rh[i] += DEN(r[i]);
    assert(rh[i]>0);
  }
  // 2.3 find the smallest index beyond which mh does not seem to change
  for(nm=0; nm!=n1; ++nm)
    if(mh[nm+1]==mh[nm]) break;
  // 2.4 set parameters for extrapolation at r>rmax
  tr  = DEN.trunc_radius() > 0;
  go  = DEN.outer_gamma();
  g3  = 3-go;
  g4  = 4-go;
  fmt = tr? 0. : -FPi*  rh[n1]  *r[n1]*r[n1]*r[n1]/g3;
  fmh = tr? 0. : -FPi*DEN(r[nm])*r[nm]*r[nm]*r[nm]/g3;
  Mtt = mt[n1]+fmt;
  Mht = mh[nm]+fmh;
  // 2.5 find total potential & Ec: add psi_halo to ps[]; get ec
  ec.reset(n);
  {
    double P = Mht/r[n1],                      // Psi_h(r_max)
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
  // 2.6 find largest index below which ps and ec do change
  for(n0=1; n0!=n1; ++n0) 
    if(ps[n0-1]!=ps[n0] && ps[n0]!=ps[n0+1] &&
       ec[n0-1]!=ec[n0] && ec[n0]!=ec[n0+1]) break;
  n0--;
  // 3   make penta spline for potential and force
  dp.reset(n);
  for(int i=0; i!=n; ++i)
    dp[i] = -mt[i]/r[i];
  PS = new Pspline<double>(lr,ps,dp);
}
// destructor
HaloPotential::~HaloPotential() {
  if(PS) falcON_DEL_O(PS);
}
// potential Psi(r) and (dPsi/dr)/r using penta spline
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
  else if(lR>lr[n1]){ P = Mtt/sqrt(Rq);
                      A =-P/Rq; }
  else              { P = (*PS)(lR,&A);
                      A/= Rq; }
  return -P;
}
// potential Psi(r)
double HaloPotential::Ps(double R) const {
  if(R<= 0.)   return ps0;
  if(R< r[0] ) {
    if(At< 2.) return ps0 - (ps0-ps[0]) * pow(R/r[0],2-At);
    if(At==2.) return ps[0] * (log(R)-lr[0]);
    else       return ps[0] * pow(R/r[0],2-At);
  }
  if(R> r[n1]) return Mtt/R;
  else         return Polev(log(R),lr,ps);
}
// log R_psi(E)
double HaloPotential::lnRPsi(double P) const {
  if(P> ps[n0]) {
    if(At< 2.) return lr[n0]+log((ps0-P)/(ps0-ps[n0]))/(2-At);
    if(At==2.) return lr[n0]-P/ps[n0];
    else       return lr[n0]+log(P/ps[n0])/(2-At);
  }
  if(P<ps[n1]) return log(Mtt/P);
  else         return Polev(P,ps.array()+n0,lr.array()+n0,n-n0);
}
// R_psi(E)
double HaloPotential::RPsi(double P) const {
  if(P> ps[n0]) {
    if(At< 2.) return r[n0]*pow((ps0-P)/(ps0-ps[n0]),1/(2-At));
    if(At==2.) return r[n0]*exp(-P/ps[n0]);
    else       return r[n0]*pow(P/ps[n0],1/(2-At));
  }
  if(P<ps[n1]) return Mtt/P;
  else         return exp(Polev(P,ps.array()+n0,lr.array()+n0,n-n0));
}
// total cumulative mass
double HaloPotential::Mt(double R) const {
  if(R<= 0.)  return 0.;
  if(R<r[0] ) return mt[0]*pow(R/r[0],3-At);
  if(R>r[n1]) return Mtt - fmt*pow(R/r[n1],g3);
  else        return Polev(log(R),lr,mt);
}
// total density
double HaloPotential::rht(double R) const {
  if(R<= 0.)  return 0.;
  if(R<r[0] ) return rh[ 0]*pow(r[ 0]/R,At);
  if(R>r[n1]) return rh[n1]*pow(r[n1]/R,go);
  else        return Polev(log(R),lr,rh);
}
// cumulative halo mass
double HaloPotential::Mh(double R) const {
  if(R<= 0.)  return 0.;
  if(R<r[ 0]) return mh[0]*pow(R/r[0],3-Ah);
  if(R>r[nm]) return Mht - fmh*pow(R/r[nm],g3);
  else        return Polev(log(R),lr,mh);
}
// total mass density
double HaloPotential::rhot(double R) const {
  if(R<=  0.) return 0.;
  if(R<r[ 0]) return rh[ 0]*pow(R/r[ 0],-At);
  if(R>r[n1]) return rh[n1]*pow(R/r[n1],-go);
  else        return Polev(log(R),lr,rh);
}
// v_circ^2(r)
double HaloPotential::vcq(double R) const {
  if(R<= 0.) return 0.;
  return Mt(R)/R;
}
// Omega^2(r)
double HaloPotential::omq(double R) const {
  if(R <= 0.) return 0.;
  return Mt(R)/cube(R);
}
// kappa^2(r)
double HaloPotential::kpq(double R) const {
  if(R <= 0.) return 0.;
  return omq(R) + FPi * rhot(R);
}
// gamma := 2*Omega/kappa 
double HaloPotential::gam(double R) const {
  if(R < r[ 0]) return 2./sqrt(4.-At);
  if(R > r[n1]) return 1.;
  double oq = Mt(R)/cube(R);
  return 2.*sqrt(oq/(oq + FPi*rhot(R)));
}
// Eps_c(r)
double HaloPotential::Epc(double R) const {
  if(R<=0.)    return ps0;
  if(R< r[ 0]) return Ps(R) - 0.5*vcq(R);
  if(R> r[n1]) return 0.5*Mtt/R;
  else         return Polev(log(R),lr,ec);
}
// R_circ(E)
double HaloPotential::RcE(double E) const {
  if(At>=2 && E>=ps0) return 0.;
  if(E >ec[n0]) {
    if(At< 2.)  return r[n0]*pow((ps0-E)/(ps0-ec[n0]),1/(2-At));
    if(At==2.)  return r[n0]*exp(0.5*(E/ec[n0]-1));
    else        return r[n0]*pow(E/ec[n0],1/(2-At));
  }
  if(E< ec[n1]) return Mtt/(E+E);
  else          return exp(Polev(E,ec.array()+n0,lr.array()+n0,n-n0));
}
// estimate for R(E, L^2, cos[eta])
double HaloPotential::Rap(double E, double Lq, double ce) const {
  const double
    Rc  = RcE(E),
    ecc = sqrt(1-Lq/(Mt(Rc)*Rc));
  return Rc * pow(1+ecc*ce, 0.5*gam(Rc));
}
// estimate for R_peri(E, L^2)
double HaloPotential::Rp(double E, double Lq) const {
  return Rap(E,Lq,-1.);
}
// estimate for R_apo(E, L^2)
double HaloPotential::Ra(double E, double Lq) const {
  return Rap(E,Lq,1.);
}
// radius given halo mass
double HaloPotential::RMh(double M) const falcON_THROWING {
  if(M<= 0.)     return 0.;
  if(M<= mh[ 0]) return r[0]*pow(M/mh[0],1/(3-Ah));
  if(M>  Mht   ) falcON_THROW("HaloPotential::RMh(): M>M_halo(oo)");
  if(M>  mh[nm]) return r[nm]*pow((Mht-M)/fmh,1/g3);
  return exp(Polev(M,mh.array(),lr.array(),nm+1));
}
//
// class HaloModel
//
namespace {
  //
  // class ReducedDensity                                                       
  //
  // to be used with complete DF : L^-2b g(Q=Eps-u^2*L^2/2)
  //                                                                            
  // returns reduced density f(r) * rho(r) with reduction factor
  //
  //                          1-b   2b
  //      f(r) = (1 + r^2 u^2)     r 
  //                                            b
  //           = (1 + r^2 u^2) (r^2/[1+r^2 u^2])
  //
  class ReducedDensity : public HaloDensity {
    const HaloDensity&Model;
    const double uq, b, tb, tb1;
    //
    double reduc(double r) const {
      if(b == 0.) return 1+r*r*uq;
      if(uq== 0.) return pow(r,tb);
      const double rq=r*r, ft=1+rq*uq;
      return ft * pow(rq/ft,b);
    }
    //
    double reduc(double r, double &d1) const {
      if(b == 0.) {
	double t=r*uq;
	d1 = t+t;
	return 1+r*t;
      }
      if(uq == 0.) {
	double p=tb1? pow(r,tb1) : 1.;
	d1 = tb*p;
	return r*p;
      }
      double rq=r*r, qq=rq*uq, ft=1+qq;
      double fc=pow(rq/ft,b);
      d1 = twice(b+qq)*fc/r;
      return ft*fc;
    }
    //
    double reduc(double r, double&d1, double&d2) const {
      if(b == 0.) {
	double t=r*uq;
	d1 = t+t;
	d2 = uq+uq;
	return 1+r*t;
      }
      if(uq== 0.) {
	const double p=tb1? pow(r,tb1) : 1.;
	d1 = tb*p;
	d2 = d1*tb1/r;
	return r*p;
      }
      const double rq=r*r, qq=rq*uq, ft=1+qq, ift=1/ft, al=twice(b+qq)*ift;
      const double fc=ft * pow(rq*ift,b);
      d1 = fc/r;
      d2 = d1/r*(al*(al-1)+4*qq*(1-b)*ift*ift);
      d1*= al;
      return fc;
    }
    //
  public:
    ReducedDensity(const HaloDensity&m, double r_a, double beta)
      : Model(m), uq(r_a>0? 1/square(r_a):0), b(beta), tb(b+b), tb1(tb-1)
    {}
    //
    double inner_gamma() const {
      return Model.inner_gamma() - tb;
    }
    //
    double scale_radius() const {
      return Model.scale_radius();
    }
    //
    double trunc_radius() const {
      return Model.trunc_radius();
    }
    //
    double outer_gamma() const {
      return Model.outer_gamma() - (uq? 2 : tb);
    }
    //
    double operator() (double r) const {
      return reduc(r)*Model(r);
    }
    //
    double operator() (double r, double &d1) const {
      double r1,re=reduc(r,r1);
      double m1,mo=Model(r,m1);
      d1 = re*m1 + r1*mo;
      return re*mo;
    }
    //
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
  //  26-05-2011
  //  NOTE We re-normalise the integral by dividing the integrand (usually
  //       d^2rho/dPsi^2) by its maximum value. Strangely, this improves the
  //       performance of qbulir (hinting at a problem with qbulir).
  //----------------------------------------------------------------------------
  const ReducedDensity*RED=0;
  double nu,MT,Be,Q,iI;         // 1/(1-p), Mtot, beta, Q, 1/in[i]
  inline double intgQ_of_q(double q) {
    //
    //   d Psi       Q^(1-p)                                      1
    // ----------- = ------- dq    with  Psi = Q [1-q^nu]   nu = ---
    // (Q - Psi)^p     1-p                                       1-p
    //                                                                          
    // Q and nu=1/(1-p) are both constants
    //                                                                          
    const double psi = Q * (1-pow(q,nu));
    if     (psi <= 0.)
      // Psi < 0 should never happen
      return 0.;
    if(psi < SPLINE->last_X()) {
      // Psi < Psi(r_max)
      // assuming  Psi = MT/r  to compute integrand
      if(Be > 0.5) {
	double r=MT/psi,rd1;
	RED->operator()(r,rd1);
	return -r*iI*rd1/psi;
      }
      if(Be >-0.5){
	double r=MT/psi,rd1,rd2;
	RED->operator()(r,rd1,rd2);
 	return r*iI/(psi*psi)*(twice(rd1) + r*rd2);
      }
    } 
    // Psi(r_max) < Psi < Psi(r_min) 
    return iI*(*SPLINE)(psi);
  }
  //
  inline double intgQ(double z) {                  // q = z^4; dq = 4 z^3 dz    
    const double f=cube(z);
    return 4*f*intgQ_of_q(z*f);
  }
}
//
falcON_TRAITS(ReducedDensity,"ReducedDensity");
//
// class falcON::HaloModel
//
HaloModel::HaloModel(HaloDensity const&model,
		     double beta, double r_a,
		     const acceleration*mono, double r_max) falcON_THROWING
  : HaloPotential(model,mono,r_max), B(beta), RA(r_a),
    g_nonmon(false)
{
  // 0. perform some sanity checks
  if(Ah < B+B) falcON_THROW("HaloModel: 2*beta > gamma_0: unphysical model");
  if(B >  1.0) falcON_THROW("HaloModel: beta > 1: unphysical model");
  if(B < -0.5) falcON_THROW("HaloModel: beta < -1/2 not supported");
  // 1 compute distribution function                                           
  // 1.1 tabulate integrand on grid considering possible cases for beta
  Array<double,1> in(n);
  RED = new ReducedDensity(DEN,RA,B);
  double negative_v=0;
  int    negative_i=0;
  if       (B >= 0.5) {        //  0.5 <= beta <  1  :  tabulate d red/ dpsi    
    bool integrand_negative = false;
    for(int i=0; i!=n; ++i) {
      double rd1,
	ps1 = -mt[i]/square(r[i]);                // psi'
      (*RED)(r[i],rd1);                           // red'
      in[i] = rd1/ps1;                            // dred/dpsi
      assert(!std::isinf(in[i]) && !std::isnan(in[i]));
      if(in[i]<0) {
	integrand_negative = true;
	if(in[i]<negative_v) {
	  negative_v = in[i];
	  negative_i = i;
	}
      }
    }
    if(integrand_negative)
      falcON_Warning("HaloModel: drho/dPsi<0 min=%g at Ps(%g) r=%g",
		     negative_v,ps[negative_i],r[negative_i]);
  } else if(B >=-0.5) {        // -0.5 <= beta <  0.5:  tabulate d^2 red/ dpsi^2
    bool integrand_negative = false;
    for(int i=0; i!=n; ++i) {
      double rd1,rd2,
	ps1 =-mt[i]/square(r[i]),                 // psi'                       
	ps2 =-twice(ps1)/r[i]-FPi*(rh[i]);        // psi"                       
      (*RED)(r[i],rd1,rd2);                       // red', red"                 
      in[i] = (rd2*ps1 - rd1*ps2)/cube(ps1);      // d^2red/dpsi^2              
      assert(!std::isinf(in[i]) && !std::isnan(in[i]));
      if(in[i]<0.) {
	integrand_negative = true;
	if(in[i]<negative_v) {
	  negative_v = in[i];
	  negative_i = i;
	}
      }
    }
    if(integrand_negative)
      falcON_Warning("HaloModel: d^2rho/dPsi^2<0 min=%g at Ps(%g) r=%g",
		     negative_v,ps[negative_i],r[negative_i]);
  } else
    falcON_THROW("HaloModel: beta < -0.5 not supported");
  // 1.2 compute g(Q) of f=L^-2B * g(Q)                                         
  lg.reset(n);
  const double
    Logalfa = (1.5-B) * LogofTwo() + LogBeta(0.5,1-B) - LogofPi() +
              log( (B >= 0.5 ? 1 : ( 0.5-B)) *
	           (B >=-0.5 ? 1 : (-0.5-B)) );
  bool   g_negative = false;
  negative_i=0;
  negative_v=0;
  if(B == 0.5  ||  B ==-0.5 ||  B ==-1.5)     // trivial cases: no integration  
    for(int i=0; i!=n; ++i) {
      assert(!std::isinf(in[i]) && !std::isnan(in[i]));
      if(in[i] < 0) {
	g_negative = true;
	if(in[i] < negative_v) { 
	  negative_v = in[i];
	  negative_i = i;
	}
	falcON_Warning("HaloModel: g(Q) < 0 at Q=%20.14g =Psi(%g)\n",Q,r[i]);
      }
      lg[i] = log(in[i]) - Logalfa;
      if(i && lg[i]>lg[i-1] && !g_nonmon) {
	g_nonmon=true;
	falcON_Warning("HaloModel: g(Q) non-monotinic at Q=%20.14g\n",ps[i]);
      }
    }
  else {
    //
    // edited 8-Sep-2011 WD
    //
    // on examination of the compute F(E), it appeared it was exactly Pi^2 too
    // large (the density computed at any radius for a DM model was too large
    // by this exact factor for all radii). Therefore, I have edited below to
    // reduce lfc by 2*LogofPi.
    //
    const int mm = int(1.5-B);
    const double
      p    = 1.5-B-mm,
      p1   = 1-p,
      lfc  = log(sin(p1*Pi)) - 3*LogofPi() - Logalfa;
    nu     = 1/p1;
    Be     = B;
    MT     = mt[n1];
    SPLINE = new spline<double>(ps,in);
    for(int i=0; i!=n; ++i) {
      Q = ps[i];
      iI= 1/in[i];
      double err, g = qbulir(intgQ,0.,1.,1.e-8,&err,0,50);
      assert(!std::isinf(g) && !std::isnan(g));
      if(g<0.) {
	g_negative = true;
	g *= pow(Q,p1);
	if(g < negative_v) { 
	  negative_v = g;
	  negative_i = i;
	}
	falcON_Warning("HaloModel: g(Q) < 0 at Q=%20.14g =Psi(%g)\n",Q,r[i]);
      }
      else if(err>1.e-3) 
	falcON_Warning("HaloModel: inaccurate integration for g(Q=%20.14g)\n",
		       Q);
      lg[i] = lfc + p1*log(Q) + log(nu*g*in[i]);
      if(i && lg[i]>lg[i-1] && !g_nonmon) {
	g_nonmon=true;
	falcON_Warning("HaloModel: g(Q) non-monotinic at Q=%20.14g\n",Q);
      }
    }
    if(g_negative) negative_v *= exp(lfc);
    delete SPLINE;
  }
  if(g_negative)
    falcON_THROW("HaloModel: g(Q)<0 g_min=%g at Q=%g = Ps(%g)",
		 negative_v,ps[negative_i],r[negative_i]);
  falcON_DEL_O(RED); RED=0;
}
// g(Q)
double HaloModel::lnG(double Q) const falcON_THROWING {
  if(Q < ps[n1]) {
    // 1 at small Q (large radii)
    if(Q < 0.) falcON_THROW("HaloModel::lnG(): Q=%g < 0",Q);
    if(RA) // 1.1 finite anisotropy radius
      return lg[n1] + (3.5-Ch-B) * (lnRPsi(Q)-lr[n1]);
    else   // 1.2 infinite anisotropy radius (recognised as RA=0)
      return lg[n1] + (1.5-Ch+B) * (lnRPsi(Q)-lr[n1]);
  } else if(Q > ps[0])
    // 2 at large Q (small radii)
    return lg[0] + (3+Ah-1.5*At-B*(4-At)) * (lr[0]-lnRPsi(Q));
  else
    // 3 within tabulated interval
    return Polev(Q,ps,lg);
}
// distribution function f(Eps,L^2)
double HaloModel::fEL (double E, double Lq) const {
  return pow(Lq,-B) * exp(lnG(RA? E-0.5*Lq/(RA*RA) : E));
}
////////////////////////////////////////////////////////////////////////////////
