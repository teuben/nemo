// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// nsam.cc                                                                     |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2004                                               |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <public/nsam.h>
#include <public/exit.h>
#include <public/nums.h>
#include <public/Pi.h>
#include <cmath>

using namespace nbdy;

#ifdef falcON_PROPER
#  include <proper/nsam.cc>
#endif

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class nbdy::SphericalSampler                                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
double SphericalSampler::
pseudo_random(                                     // R: Psi(r)                 
	      PseudoRandom const&RNG,              // I: RNG for all purposes   
	      double            &r,                // O: radius                 
	      double            &vr,               // O: radial velocity        
	      double            &vt) const         // O: tangential velocity    
{
  double w,f;
  //                                                                            
  // 1. get r & Psi(r)                                                          
  //                                                                            
  r = rM(RNG()*Mt);                                // get Lagrange radius from M
  double Psi = Ps(r);                              // get potential             
  //                                                                            
  // 2. get w                                                                   
  //                                                                            
  double we  = sqrt(2*Psi);                        // escape velocity           
  double f0  = DF(Psi);                            // DF(w=0)                   
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
	     double      &vt) const                // O: tangential velocity    
{
  double w,f;
  //                                                                            
  // 1. get r & Psi(r)                                                          
  //                                                                            
  r = rM(RNG(0)*Mt);                               // get Lagrange radius from M
  double Psi = Ps(r);                              // get potential             
  //                                                                            
  // 2. get w                                                                   
  //                                                                            
  double we  = sqrt(2*Psi);                        // escape velocity           
  double f0  = DF(Psi);                            // DF(w=0)                   
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
			      double const&f_      //[I: fraction with vphi>0]  
#ifdef falcON_PROPER
			     ,double const&epar    //[I: factor: setting eps_i] 
#endif
			      ) const
{
  const double m  (Mt/double(B.N_bodies()));       // Mt/N: body mass           
  const double fp (2*f_-1);                        // fraction: swap sign(vphi) 
  int    n;                                        // # bodies / (r,vr,vt)      
  double Ncum=0.,Mcum=0.,mu=m;                     // cumulated num, mu         
  for(int k=0; k!=B.N_bodies();) {                 // LOOP bodies               
    //                                                                          
    // 1. get r,vr,vt                                                           
    //                                                                          
    double r,vr,vt;                                //   radius, rad. & tan. vel.
    double Psi = q?                                //   quasi-randomly?         
      quasi_random (RNG,r,vr,vt) :                 //     YES: get r,vr,vt      
      pseudo_random(RNG,r,vr,vt);                  //     NO:  get r,vr,vt      
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
    for(int i=0; i!=n && k!=B.N_bodies();++i,++k) {//   LOOP n bodies           
      B.mass(k) = mu;                              //     set mass              
      Mcum     += mu;                              //     cumulate total mass   
      register double                              //     some auxiliary vars:  
	cth = q? RNG(2,-1.,1.):RNG(-1.,1.),        //     sample cos(theta)     
	sth = std::sqrt(1.-cth*cth),               //     sin(theta)            
	phi = q? RNG(3,0.,TPi):RNG(0.,TPi),        //     sample azimuth phi    
	cph = std::cos(phi),                       //     cos(phi)              
	sph = std::sin(phi);                       //     sin(phi)              
      B.pos(k)[0] = r * sth * cph;                 //     x=r*sin(th)*cos(ph)   
      B.pos(k)[1] = r * sth * sph;                 //     y=r*sin(th)*sin(ph)   
      B.pos(k)[2] = r * cth;                       //     z=r*cos(th)           
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
      B.vel(k)[0] = vm * cph - vph * sph;          //     v_x= ...              
      B.vel(k)[1] = vm * sph + vph * cph;          //     v_y= ...              
      B.vel(k)[2] = vr * cth - vth * sth;          //     v_z= ...              
    }                                              //   END LOOP                
  }                                                // END LOOP                  
  //                                                                            
  // 4. adjust total mass [& set eps]                                           
  //                                                                            
  if(Mcum != Mt) {
    const double mfac = Mt/Mcum;
    LoopBodies(bodies,&B,Bi)
      Bi.mass()*= mfac;
  }
  //                                                                            
  // 5. set eps_i                                                               
  //                                                                            
#ifdef falcON_PROPER
  if(epar>0. && B.has(io::e)) {
    const double iMt  = 1./Mt;
    LoopBodies(bodies,&B,Bi)
      Bi.eps () = epar * sqrt(mass(Bi)*iMt);
  }
#endif
}
