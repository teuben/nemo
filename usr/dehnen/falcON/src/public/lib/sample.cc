// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file    src/public/lib/sample.cc
///
/// \author  Walter Dehnen
///
/// \date    2004-2009
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004-2009  Walter Dehnen
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
#include <public/sample.h>
#include <public/basic.h>
#include <numerics.h>
#include <cmath>

using namespace falcON;

//
// class falcON::SphericalSampler
//
namespace {
  const int Ne1=1000, Ne=Ne1+1;
  double __p;
  inline double dI(double eta, double const&I)
  {
    return pow(sin(eta),__p);
  }
}
//
inline void SphericalSampler::setis()
  // set table with int(sin(eta)^(1-2b), eta=0..x)
{
  if(!beta) return;
  Is[0]    = 0.;
  Xe[0][0] = 0.;
  Xe[0][1] = 1.;
  const double de =Pi/Ne1;
  double __eta = 0.;
  __p = 1-b0-b0;
  for(register int i=1; i!=Ne; ++i) {
    Is[i]    = rk4(Is[i-1],__eta,de,dI);
    __eta   += de;
    Xe[i][0] = sin(__eta);
    Xe[i][1] = cos(__eta);
  }
}
//
#ifdef falcON_PROPER
#  include <proper/sample.cc>
#endif
SphericalSampler::SphericalSampler(double __mt,
				   double __ra,
				   double __b0,
				   bool   __c) :
  careful ( __c ),
  OM      ( __ra>0. ),
  beta    ( __b0 != 0. ),
  Mt      ( __mt ),
  ra      ( __ra ),
  iraq    ( __ra>0? 1./(__ra*__ra) : 0. ),
  b0      ( __b0 ),
  ibt     ( 1./(3. - __b0 - __b0) ),
  Xe      ( beta? Ne : 0 ),
  Is      ( beta? Ne : 0 )
#ifdef falcON_PROPER
  ,
  adapt_masses ( false ),
  Peri         ( false ),
  irs          ( 0 ),
  eta          ( 1 ),
  mmm          ( 1 ),
  nmax         ( 1 )
#endif
{
  if(beta) setis();
}
//
inline double SphericalSampler::F0(double Psi) const
{
  if(careful) {
    // be careful in finding a value f0 so that f(Eps) < f0 for Eps in [0,Psi]
    // but not that f0 is MUCH larger than the maximum of f in that range
    const int    n = 11;
    const double s[n] = {0.975,0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1};
    // 1st find maximum of f several points
    double f0=DF(Psi);
    for(int i=0; i!=n; ++i) {
      double fi = DF(s[i]*Psi);
      if(fi > f0) f0 = fi;
    }
    // 2nd allow for an extra 200% higher maximum (extra safe)
    return f0+f0+f0;
  } else
    return DF(Psi);
}
//
void SphericalSampler::sample(body   const&B0,
			      unsigned     N,
			      bool         q,
			      Random const&R,
			      double       f_,
#ifdef falcON_PROPER
			      double       epar,
#endif
			      bool         givef,
			      bool         giveP,
			      bool         giveA) const falcON_THROWING
  {
  if(givef && !has_phden(B0))
    falcON_Warning("SphericalSampler: "
		   "bodies not supporting phden: cannot give DF");
  if(giveP && !has_pot(B0))
    falcON_Warning("SphericalSampler: "
		   "bodies not supporting pot: cannot give Phi");
  if(giveA && !has_acc(B0))
    falcON_Warning("SphericalSampler: "
		   "bodies not supporting acc: cannot give dPhi/dr");
  if(!(B0+(N-1)).is_valid())
    falcON_THROW("SphericalSampler::sample(): not enough bodies free");
  if(q && R.Nsob() < 6)
    falcON_THROW("SphericalSampler::sample(): "
		 "too few quasi-random number generators\n");
  const body BN(B0,N);
  givef = givef && has_aux(B0);
  giveP = giveP && has_pot(B0);
  giveA = giveA && has_acc(B0);
  const double m  (Mt/double(N));
  const double fp (2*f_-1);
  int    n, k=0;
  double Ncum=0.,Mcum=0.,mu=m;
  for(body Bi(B0); Bi!=BN; ) {
    //
    // 1. get Mr,r,Psi,vr,vt,f
    //
    double Mr = (q? R(0):R())*Mt;
    double r  = rM(Mr);
    double Psi= Ps(r),Q;
    double w, we=sqrt(2*Psi);
    double f0=F0(Psi),f;
    do {
      do {
	w = we * pow(R(),ibt);
	Q = Psi-0.5*w*w;
      } while(Q<=0.);
      f = DF(Q);
      if(WDutils::isnan(f))
	falcON_THROW("SphericalSampler::sample(): "
		     "%s is NaN; Eps=%g [M=%g, r=%g, we=%g, w=%g]\n",
		     (beta? (iraq==0? "g(E)" : "g(Q)"):
		      (iraq==0? "f(E)" : "f(Q)")), Q,Mr,r,we,w);
      if(f>f0)
	falcON_THROW("SphericalSampler::sample(): DF non-monotonic"
		     ": f(Psi=%g)=%g > f(Eps=%g)=%g [r=%g]\n", Q,f,Psi,f0,r);
    } while(f0 * R() > f);
    double vr,vt;
    if(beta) {
      pair_d SC=polev((q? R(1):R())*Is[Ne1],Is,Xe);
      vt = w * SC[0]/sqrt(1+r*r*iraq);
      vr = w * SC[1];
      if(b0<zero || vt>zero) f *= pow(vt*r,-b0-b0);
      if(givef && WDutils::isnan(f))
	falcON_THROW("SphericalSampler::sample(): f(E,L^2) is NaN"
		     ": Eps=%g, [vt=%g, r=%g]\n",Psi-0.5*w*w,vt,r);
    } else {
      double ce = q? R(1,-1.,1.) : R(-1.,1.);
      vr = w * ce;
      vt = w * sqrt((1-ce*ce)/(1+r*r*iraq));
    }
    //
    // 2. establish number of bodies at (r,vr,vt), depending on mass adaption
    //#
#ifdef falcON_PROPER
    if(adapt_masses)
      falcON_SAMPLE_ADAPT
    else
#endif
    {
      ++Ncum;
      n  = 1;
    }
    if(n<1) continue;
    //
    // 3. set mass, position, and velocity of n bodies with this (r,vr,vt)
    //
    for(int i=0; i!=n && Bi!=BN; ++i,++k,++Bi) {
      Bi.mass() = mu;
      Mcum     += mu;
      register double
	cth = q? R(2,-1.,1.):R(-1.,1.),
	sth = std::sqrt(1.-cth*cth),
	phi = q? R(3,0.,TPi):R(0.,TPi),
	cph = std::cos(phi),
	sph = std::sin(phi);
      Bi.pos()[0] = r * sth * cph;
      Bi.pos()[1] = r * sth * sph;
      Bi.pos()[2] = r * cth;
      register double
	psi = q? R(4,0.,TPi):R(0.,TPi),
        vth = vt * std::cos(psi),
        vph = vt * std::sin(psi),
        vm  = vr * sth + vth * cth;
      if       (fp>0.) {
	if(fp>  (q? R(5):R())) vph= abs(vph);
      } else if(fp<0.) {
	if(fp< -(q? R(5):R())) vph=-abs(vph);
      }
      Bi.vel()[0] = vm * cph - vph * sph;
      Bi.vel()[1] = vm * sph + vph * cph;
      Bi.vel()[2] = vr * cth - vth * sth;
      if(givef) Bi.phden() = f;
      if(giveP) Bi.pot() =-Psi;
      if(giveA) Bi.acc() = pos(Bi)*(-Mr/(r*r));
    }
  }
  //
  // 4. adjust total mass
  //
  if(Mcum != Mt) {
    const double mfac = Mt/Mcum;
    for(body Bi(B0); Bi!=BN; ++Bi)
      Bi.mass()*= mfac;
  }
  //
  // 5. set eps_i                                                               
  //
#ifdef falcON_PROPER
  if(epar>0. && has_eps(B0)) {
    const double iMt  = double(N)/Mt;
    for(body Bi(B0); Bi!=BN; ++Bi)
      Bi.eps () = epar * sqrt(mass(Bi)*iMt);
  }
#endif
}
//
void SphericalSampler::sample_pos(body const  &B0,
				  unsigned     N,
				  bool         q,
				  Random const&R) const falcON_THROWING
  {
  if(!(B0+(N-1)).is_valid())
    falcON_THROW("SphericalSampler::sample_pos(): not enough bodies free");
  if(q && R.Nsob() < 6)
    falcON_THROW("SphericalSampler::sample_pos(): "
		 "too few quasi-random number generators\n");
  const body BN(B0,N);
  const double m  (Mt/double(N));
  for(body Bi(B0); Bi!=BN; ++Bi) {
    double Mr = (q? R(0):R())*Mt;
    double r  = rM(Mr);
    Bi.mass() = m;
    double
      cth = q? R(1,-1.,1.):R(-1.,1.),
      sth = std::sqrt(1.-cth*cth),
      phi = q? R(2,0.,TPi):R(0.,TPi);
    Bi.pos()[0] = r * sth * std::cos(phi);
    Bi.pos()[1] = r * sth * std::sin(phi);
    Bi.pos()[2] = r * cth;
  }
}
//
