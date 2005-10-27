// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// sample.cc                                                                   |
//                                                                             |
// Copyright (C) 2004-2005  Walter Dehnen                                      |
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
#include <public/numerics.h>
#include <public/Pi.h>
#include <cmath>

using namespace falcON;

#ifdef falcON_PROPER
#  include <proper/sample.cc>
#endif

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::SphericalSampler                                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
double SphericalSampler::
pseudo_random(                                     // R: Psi(r)                 
	      PseudoRandom const&RNG,              // I: RNG for all purposes   
	      double            &r,                // O: radius                 
	      double            &vr,               // O: radial velocity        
	      double            &vt,               // O: tangential velocity    
	      double            &f) const          // O: f(E) or g(Q)           
{
  //                                                                            
  // 1. get r & Psi(r)                                                          
  //                                                                            
  r = rM(RNG()*Mt);                                // get Lagrange radius from M
  double Psi = Ps(r);                              // get potential             
  //                                                                            
  // 2. get w                                                                   
  //                                                                            
  double w, we=sqrt(2*Psi);                        // escape velocity           
  double f0=DF(Psi);                               // DF(w=0)                   
  do {                                             // DO                        
    w = we * pow(RNG(),ibt);                       //   sample velocity 0<w<we  
    f = DF(Psi-0.5*w*w);                           //   get DF                  
    if(f>f0)                                       //   IF f>f(w=0)             
      error("sampling error: DF non-monotonic"     //     ERROR                 
	    ": f(%f)=%f < f(%f)=%f\n",
	    f,Psi-0.5*w*w,f0,Psi);                 //     ERROR                 
  } while(f0 * RNG() > f);                         // WHILE ( rejected )        
  //                                                                            
  // 3. get angle in velocity space -> vr, vt                                   
  //                                                                            
#ifdef falcON_PROPER
  if(beta) falcON_SAMPLE_SETV else                 // IF b0 != 0                
#endif
  {                                                // ELSE                      
    double ce = RNG(-1.,1.);                       //   sample cos(eta)         
    vr = w * ce;                                   //   radial velocity         
    vt = w * sqrt((1-ce*ce)
#ifdef falcON_PROPER
		  /(1+r*r*iraq)
#endif
		  );                               //   tangential velocity     
  }                                                // ENDIF                     
  return Psi;                                      // return potential          
}
////////////////////////////////////////////////////////////////////////////////
double SphericalSampler::
quasi_random(                                      // R: Psi(r)                 
	     Random const&RNG,                     // I: RNG for all purposes   
	     double      &r,                       // O: radius                 
	     double      &vr,                      // O: radial velocity        
	     double      &vt,                      // O: tangential velocity    
	     double      &f) const                 // O: f(E) or g(Q)           
{
  //                                                                            
  // 1. get r & Psi(r)                                                          
  //                                                                            
  r = rM(RNG(0)*Mt);                               // get Lagrange radius from M
  double Psi = Ps(r);                              // get potential             
  //                                                                            
  // 2. get w                                                                   
  //                                                                            
  double w, we=sqrt(2*Psi);                        // escape velocity           
  double f0=DF(Psi);                               // DF(w=0)                   
  do {                                             // DO                        
    double __ran = RNG();
    w = we * pow(__ran,ibt);                       //   sample velocity 0<w<we  
    double Eps = Psi-0.5*w*w;
    f = DF(Eps);                                   //   get DF                  
    if(f>f0)                                       //   IF f>f(w=0)             
      error("sampling error: DF non-monotonic"     //     ERROR                 
	    ": f(Eps) > f0:=f(Psi) where\n"
	    "\t\t Psi=%f, f0=%f, Eps=%f, f=%f, ran=%f, we=%f, w=%f\n",
	    Psi,f0,Eps,f,__ran,we,w);              //     ERROR                 
  } while(f0 * RNG() > f);                         // WHILE ( rejected )        
  //                                                                            
  // 3. get angle in velocity space -> vr, vt                                   
  //                                                                            
#ifdef falcON_PROPER
  if(beta) falcON_SAMPLE_SETV else                 // IF b0 != 0                
#endif
  {                                                // ELSE                      
    double ce = RNG(1,-1.,1.);                     //   sample cos(eta)         
    vr = w * ce;                                   //   radial velocity         
    vt = w * sqrt((1-ce*ce)
#ifdef falcON_PROPER
		  /(1+r*r*iraq)
#endif
		  );                               //   tangential velocity     
  }                                                // ENDIF                     
  return Psi;                                      // return potential          
}
////////////////////////////////////////////////////////////////////////////////
void SphericalSampler::sample(bodies const&B,      // I/O: bodies to sample     
			      bool   const&q,      // I: quasi random?          
			      Random const&RNG,    // I: pseudo & quasi RNG     
			      double const&f_,     //[I: fraction with vphi>0]  
#ifdef falcON_PROPER
			      double const&epar,   //[I: factor: setting eps_i] 
#endif
			      bool         givef)  //[I: write DF into aux?]    
  const {
  if(givef && !B.have(fieldbit::y))
    warning("SphericalSampler: bodies not supporting aux: cannot give DF");
  givef = givef && B.have(fieldbit::y);
  const double m  (Mt/double(B.N_bodies()));       // Mt/N: body mass           
  const double fp (2*f_-1);                        // fraction: swap sign(vphi) 
  int    n, k=0;                                   // # bodies / (r,vr,vt)      
  double Ncum=0.,Mcum=0.,mu=m;                     // cumulated num, mu         
  for(body Bi=B.begin_all_bodies(); Bi; ) {        // LOOP bodies               
    //                                                                          
    // 1. get r,vr,vt,f                                                         
    //                                                                          
    double r,vr,vt,f;                              //   radius, v_tan, v_rad, DF
    double Psi = q?                                //   quasi-randomly?         
      quasi_random (RNG,r,vr,vt,f) :               //     YES: get r,vr,vt,f    
      pseudo_random(RNG,r,vr,vt,f);                //     NO:  get r,vr,vt,f    
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
    for(int i=0; i!=n && Bi; ++i,++k,++Bi) {       //   LOOP n bodies           
      if(givef) Bi.aux() = f;                      //     set DF                
      Bi.mass() = mu;                              //     set mass              
      Mcum    += mu;                               //     cumulate total mass   
      register double                              //     some auxiliary vars:  
	cth = q? RNG(2,-1.,1.):RNG(-1.,1.),        //     sample cos(theta)     
	sth = std::sqrt(1.-cth*cth),               //     sin(theta)            
	phi = q? RNG(3,0.,TPi):RNG(0.,TPi),        //     sample azimuth phi    
	cph = std::cos(phi),                       //     cos(phi)              
	sph = std::sin(phi);                       //     sin(phi)              
      Bi.pos()[0] = r * sth * cph;                 //     x=r*sin(th)*cos(ph)   
      Bi.pos()[1] = r * sth * sph;                 //     y=r*sin(th)*sin(ph)   
      Bi.pos()[2] = r * cth;                       //     z=r*cos(th)           
      register double                              //     some auxiliary vars:  
	psi = q? RNG(4,0.,TPi):RNG(0.,TPi),        //     sample angle psi      
        vth = vt * std::cos(psi),                  //     v_theta               
        vph = vt * std::sin(psi),                  //     v_phi                 
        vm  = vr * sth + vth * cth;                //     v_meridional          
      if       (fp>0.) {                           //     IF flip sign to pos?  
	if(fp>  (q? RNG(5):RNG())) vph= abs(vph);  //       draw -> flip        
      } else if(fp<0.) {                           //     ELIF flip sign to neg?
	if(fp< -(q? RNG(5):RNG())) vph=-abs(vph);  //       draw -> flip        
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
    LoopAllBodies(&B,Bi)
      Bi.mass()*= mfac;
  }
  //                                                                            
  // 5. set eps_i                                                               
  //                                                                            
#ifdef falcON_PROPER
  if(epar>0. && B.have(fieldbit::e)) {
    const double iMt  = 1./Mt;
    LoopAllBodies(&B,Bi)
      Bi.eps () = epar * sqrt(mass(Bi)*iMt);
  }
#endif
}
