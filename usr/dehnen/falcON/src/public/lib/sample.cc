// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// sample.cc                                                                   |
//                                                                             |
// Copyright (C) 2004-2007  Walter Dehnen                                      |
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
#include <public/sample.h>
#include <public/basic.h>
#include <numerics.h>
#include <cmath>

using namespace falcON;

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::SphericalSampler                                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace {
  const int Ne1=1000, Ne=Ne1+1;
  double __p;
  inline double dI(double eta, double const&I)
  {
    return pow(sin(eta),__p);
  }
}
////////////////////////////////////////////////////////////////////////////////
inline void SphericalSampler::setis()
  // set table with int(sin(eta)^(1-2b), eta=0..x)                              
{
  if(!beta) return;
  Is[0]    = 0.;
  Xe[0][0] = 0.;
  Xe[0][1] = 1.;
  const double    de =Pi/Ne1;
  double eta = 0.;
  __p = 1-b0-b0;
  for(register int i=1; i!=Ne; ++i) {
    Is[i]    = rk4(Is[i-1],eta,de,dI);
    eta     += de;
    Xe[i][0] = sin(eta);
    Xe[i][1] = cos(eta);
  }
}
////////////////////////////////////////////////////////////////////////////////
#ifdef falcON_PROPER
#  include <proper/sample.cc>
#else
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
{
  if(beta) setis();
}
#endif
////////////////////////////////////////////////////////////////////////////////
inline double SphericalSampler::F0(double Psi) const
{
  if(careful) {
    // be careful in finding a value f0 so that f(Eps) < f0 for Eps in [0,Psi]
    // but not that f0 is MUCH larger than the maximum of f in that range
    const int    n = 5;
    const double s[n] = {0.9,0.8,0.7,0.6,0.5};
    // 1st find maximum of f at 0.5, 0.6, 0.7, 0.8, 0.9, and 1.0 times Psi
    double f0=DF(Psi);
    for(int i=0; i!=n; ++i) {
      double fi = DF(s[i]*Psi);
      if(fi > f0) f0 = fi;
    }
    // 2nd allow for an extra 50% higher maximum (extra safe)
    return 1.5*f0;
  } else
    return DF(Psi);
}
////////////////////////////////////////////////////////////////////////////////
double SphericalSampler::
set_radvel(                                        // R: Psi(r)                 
	   bool         q,                         // I: quasi (or random) RNG? 
	   Random const&R,                         // I: pseudo & quasi RNG     
	   double      &r,                         // O: radius                 
	   double      &vr,                        // O: radial velocity        
	   double      &vt,                        // O: tangential velocity    
	   double      &f) const                   // O: f(E,L)                 
{
  //                                                                            
  // 1. get r & Psi(r)                                                          
  //                                                                            
  r = rM(q? R(0)*Mt : R()*Mt);                     // get Lagrange radius from M
  double Psi = Ps(r);                              // get potential             
  //                                                                            
  // 2. get w                                                                   
  //                                                                            
  double w, we=sqrt(2*Psi);                        // escape velocity           
  double f0=F0(Psi);                               // DF(w=0)                   
  do {                                             // DO                        
    w = we * pow(R(),ibt);                         //   sample velocity 0<w<we  
    f = DF(Psi-0.5*w*w);                           //   get g(Q)                
    if(f>f0)                                       //   IF f>f(w=0)             
      error("sampling error: DF non-monotonic"     //     ERROR                 
	    ": f(Psi=%g)=%g < f(Eps=%g)=%g\n",
	    f,Psi-0.5*w*w,f0,Psi);                 //     ERROR                 
  } while(f0 * R() > f);                           // WHILE ( rejected )        
  //                                                                            
  // 3. get angle in velocity space -> vr, vt  and set f=L^-2b*g(Q)             
  //                                                                            
  if(beta) {                                       // IF b0 != 0                
    pair_d SC=polev((q? R(1):R())*Is[Ne1],Is,Xe);  //   sample eta              
    vt = w * SC[0]/sqrt(1+r*r*iraq);               //   tangential velocity     
    vr = w * SC[1];                                //   radial velocity         
    if(b0<zero || vt>zero) f *= pow(vt*r,-b0-b0);  //   distribution function   
  } else {                                         // ELSE                      
    double ce = q? R(1,-1.,1.) : R(-1.,1.);        //   sample cos(eta)         
    vr = w * ce;                                   //   radial velocity         
    vt = w * sqrt((1-ce*ce)/(1+r*r*iraq));         //   tangential velocity     
  }                                                // ENDIF                     
  return Psi;                                      // return potential          
}
////////////////////////////////////////////////////////////////////////////////
void SphericalSampler::sample(body   const&B0,     // I: first body to sample   
			      unsigned     N,      // I: # bodies to sample     
			      bool         q,      // I: quasi random?          
			      Random const&R,      // I: pseudo & quasi RNG     
			      double       f_,     //[I: fraction with vphi>0]  
#ifdef falcON_PROPER
			      double       epar,   //[I: factor: setting eps_i] 
#endif
			      bool         givef)  //[I: write DF into aux?]    
  const {
  if(givef && !has_aux(B0))
    warning("SphericalSampler: bodies not supporting aux: cannot give DF");
  if(!(B0+(N-1)).is_valid())
    error("SphericalSampler: not enough bodies free");
  const body BN(B0,N);
  givef = givef && has_aux(B0);
  const double m  (Mt/double(N));                  // Mt/N: body mass           
  const double fp (2*f_-1);                        // fraction: swap sign(vphi) 
  int    n, k=0;                                   // # bodies / (r,vr,vt)      
  double Ncum=0.,Mcum=0.,mu=m;                     // cumulated num, mu         
  for(body Bi(B0); Bi!=BN; ) {                     // LOOP bodies               
    //                                                                          
    // 1. get r,vr,vt,f                                                         
    //                                                                          
    double r,vr,vt,f;                              //   radius, v_tan, v_rad, DF
    double Psi = set_radvel(q,R,r,vr,vt,f);        //   get Psi,r,vr,vt,f(E,L)  
    //                                                                          
    // 2. establish number of bodies at (r,vr,vt), depending on mass adaption   
    //                                                                          
#ifdef falcON_PROPER
    if(adapt_masses) falcON_SAMPLE_ADAPT else      //   IF mass adaption        
#endif
      n  = 1;                                      //   ELSE set n  = 1         
    if(n<1) continue;                              //   IF n<1: continue        
    //                                                                          
    // 3. set mass, position, and velocity of n bodies with this (r,vr,vt)      
    //                                                                          
    for(int i=0; i!=n && Bi!=BN; ++i,++k,++Bi) {   //   LOOP n bodies           
      if(givef) Bi.aux() = f;                      //     set DF                
      Bi.mass() = mu;                              //     set mass              
      Mcum     += mu;                              //     cumulate total mass   
      register double                              //     some auxiliary vars:  
	cth = q? R(2,-1.,1.):R(-1.,1.),            //     sample cos(theta)     
	sth = std::sqrt(1.-cth*cth),               //     sin(theta)            
	phi = q? R(3,0.,TPi):R(0.,TPi),            //     sample azimuth phi    
	cph = std::cos(phi),                       //     cos(phi)              
	sph = std::sin(phi);                       //     sin(phi)              
      Bi.pos()[0] = r * sth * cph;                 //     x=r*sin(th)*cos(ph)   
      Bi.pos()[1] = r * sth * sph;                 //     y=r*sin(th)*sin(ph)   
      Bi.pos()[2] = r * cth;                       //     z=r*cos(th)           
      register double                              //     some auxiliary vars:  
	psi = q? R(4,0.,TPi):R(0.,TPi),            //     sample angle psi      
        vth = vt * std::cos(psi),                  //     v_theta               
        vph = vt * std::sin(psi),                  //     v_phi                 
        vm  = vr * sth + vth * cth;                //     v_meridional          
      if       (fp>0.) {                           //     IF flip sign to pos?  
	if(fp>  (q? R(5):R())) vph= abs(vph);      //       draw -> flip        
      } else if(fp<0.) {                           //     ELIF flip sign to neg?
	if(fp< -(q? R(5):R())) vph=-abs(vph);      //       draw -> flip        
      }                                            //     ENDIF                 
      Bi.vel()[0] = vm * cph - vph * sph;          //     v_x= ...              
      Bi.vel()[1] = vm * sph + vph * cph;          //     v_y= ...              
      Bi.vel()[2] = vr * cth - vth * sth;          //     v_z= ...              
    }                                              //   END LOOP                
  }                                                // END LOOP                  
  //                                                                            
  // 4. adjust total mass [& set eps]                                           
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
    const double iMt  = 1./Mt;
    for(body Bi(B0); Bi!=BN; ++Bi)
      Bi.eps () = epar * sqrt(mass(Bi)*iMt);
  }
#endif
}
