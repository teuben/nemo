// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// nsam.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2004                                               |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_nsam_h
#define falcON_included_nsam_h

#ifndef falcON_included_body_h
#  include <body.h>
#endif
#ifndef falcON_included_rand_h
#  include <public/rand.h>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::SphericalSampler                                             //
  //                                                                          //
  // PUBLIC version:                                                          //
  //                                                                          //
  // We consider isotropic DFs, i.e. f=f(E)                                   //
  //                                                                          //
  //                                                                          //
  // PROPRIETARY version:                                                     //
  //                                                                          //
  // We consider DFs of the general form (Eps := -E = Psi-v^2/2)              //
  //                                                                          //
  //     f = L^(-2b) g(Q)                                                     //
  //                                                                          //
  // with Q := Eps - L^2/(2a^2)                                               //
  //                                                                          //
  // The Binney anisotropy parameter beta for these models is                 //
  //                                                                          //
  //     beta = (r^2 + b*a^2) / (r^2 + a^2),                                  //
  //                                                                          //
  // which is beta=b at r=0 and beta=1 at r=oo.                               //
  // For b=0, we obtain the Ossipkov-Merritt model. For a=oo, we obtain a     //
  // model with constant anisotropy. For b=0 and a=oo, we get the isotropic   //
  // DF.                                                                      //
  //                                                                          //
  // NOTE  1. Since a=0 makes no sense, we interprete a=0 as a=oo.            //
  //       2. If the model has a different type of DF, one may not use these  //
  //          methods but superseed the sampling routines below.              //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class SphericalSampler {
  private:
    const double Mt;
    const double ibt;
    typedef tupel<2,double> pair_d;
    //--------------------------------------------------------------------------
    const bool   adapt_masses, Peri, OM, beta;
    const double ra,iraq;
    const double b0;
    const int    nR;
    const double*Rad, fac;
    double      *num;
    pair_d      *Xe;
    double      *Is;
    int          Ne;
    //--------------------------------------------------------------------------
    // private methods                                                          
    //--------------------------------------------------------------------------
#ifdef falcON_PROPER
    inline void setis();
#endif
  public:
    //--------------------------------------------------------------------------
    // construction & destruction                                               
    //--------------------------------------------------------------------------
    SphericalSampler(double const& mt              // I: total mass             
#ifdef falcON_PROPER
		    ,double const& =0.,            //[I: a: anisotropy radius]  
		     double const& =0.,            //[I: b0]                    
		     const double* =0,             //[I: mass adaption: radii]  
		     int    const& =0,             //[I: mass adaption: # -- ]  
		     double const& =1.2,           //[I: mass adaption: factor] 
		     bool   const& =0              //[I: mass adaption: R_-/Re] 
#endif
		     )
#ifdef falcON_PROPER
      ;
#else
      : Mt(mt), ibt(1./3.) {}
#endif
    //--------------------------------------------------------------------------
    ~SphericalSampler() { 
#ifdef falcON_PROPER
      if(num) delete[] num;
      if(Xe)  delete[] Xe;
      if(Is)  delete[] Is;
#endif
    }
    //--------------------------------------------------------------------------
    // 1. Methods required for sampling radius and velocity: abstract           
    //--------------------------------------------------------------------------
  protected:
    virtual double DF(double)         const=0;     // f(E) or g(Q)              
#ifdef falcON_PROPER
    virtual double Re(double)         const=0;     // Rcirc(Eps)                
    virtual double Rp(double, double) const=0;     // R_peri(Eps, L)            
#endif
    virtual double Ps(double)         const=0;     // Psi(R)                    
    virtual double rM(double)         const=0;     // r(M), M in [0,M_total]    
    //--------------------------------------------------------------------------
    // 2. Sampling of phase space: virtual                                      
    //                                                                          
    // We implement two routines for sampling radius, radial and tangential     
    // velocity from the model. They differ in the treatment of the random      
    // numbers. The first accepts any type nbdy::PseudoRandom random number     
    // generator and the second a nbdy::Random, of which it uses the first two  
    // quasi-random number generators.                                          
    //                                                                          
    // NOTE  These methods are virtual but not abstract. That is you may super- 
    //       seed them if you want to provide your own.                         
    //--------------------------------------------------------------------------
    double pseudo_random(                          // R: Psi(r)                 
			 PseudoRandom const&,      // I: pseudo RNG             
			 double&,                  // O: radius                 
			 double&,                  // O: radial velocity        
			 double&) const;           // O: tangential velocity    
    //--------------------------------------------------------------------------
    double quasi_random (                          // R: Psi(r)                 
			 Random const&,            // I: pseudo & quasi RNG     
			 double&,                  // O: radius                 
			 double&,                  // O: radial velocity        
			 double&) const;           // O: tangential velocity    
    //--------------------------------------------------------------------------
    // 3. Sampling of full phase-space: non-virtual                             
    //                                                                          
    // NOTE on scales                                                           
    // We assume that the randomly generated (r,vr,vt) are already properly     
    // scaled as well as the Psi returnd. Also the routines for pericentre and  
    // Rc(E) are assummed to take arguments and return values that are both     
    // scaled.                                                                  
    //--------------------------------------------------------------------------
  public:
    void sample(bodies const&,                     // I/O: bodies to sample     
		bool   const&,                     // I: quasi random?          
		Random const&,                     // I: pseudo & quasi RNG     
		double const& =0.5                 //[I: fraction with vphi>0]  
#ifdef falcON_PROPER
	       ,double const& =0.0
#endif
		) const;                           //[I: factor: setting eps_i] 
    //--------------------------------------------------------------------------
  };
}
////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_included_nsam_h    
////////////////////////////////////////////////////////////////////////////////


