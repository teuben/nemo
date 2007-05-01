// -*- C++ -*-                                                                 |
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   inc/public/sample.h                                                 
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2000-2007                                                           
///                                                                             
/// \brief  code for sampling spherical stellar dynamical equilibria            
///                                                                             
/// \todo   finish doxygen documentation                                        
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2004-2007  Walter Dehnen                                       
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
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// some history                                                                 
// 27/10/2004 WD  added support for DF evaluation                               
// 16/02/2006 WD  added support for non-monotonic DF                            
//                                                                              
////////////////////////////////////////////////////////////////////////////////
#ifndef falcON_included_sample_h
#define falcON_included_sample_h

#ifndef falcON_included_body_h
#  include <body.h>
#endif
#ifndef falcON_included_random_h
#  include <public/random.h>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::SphericalSampler                                             
  //                                                                            
  /// provides routines for sampling phase-space densities                      
#ifndef falcON_PROPER
  /// PUBLIC version:  we consider isotropic DFs, i.e. f=f(E)                   
#else
  /// PROPRIETARY version: we consider DFs of the general form\n                
  ///    f = L^(-2b) g(Q)\n                                                     
  /// with\n                                                                    
  ///    Q := Eps - L^2/(2a^2)\n                                                
  /// where\n                                                                   
  ///    Eps := -E = Psi-v^2/2.\n                                               
  /// The Binney anisotropy parameter beta for these models is\n                
  ///    beta = (r^2 + b*a^2) / (r^2 + a^2),\n                                  
  /// which is beta=b at r=0 and beta=1 at r=oo.\n                              
  /// For b=0, we obtain the Ossipkov-Merritt model. For a=oo, we obtain a      
  /// model with constant anisotropy. For b=0 and a=oo, the DF is isotropic.    
  /// \note Since a=0 makes no sense, we interprete a=0 as a=oo.                
  /// \note If the model has a different type of DF, one may not use these      
  ///       methods but superseed the sampling routines below.                  
#endif
  //////////////////////////////////////////////////////////////////////////////
  class SphericalSampler {
  private:
    const double Mt;
    const double ibt;
    const bool   careful;
    typedef tupel<2,double> pair_d;
    //--------------------------------------------------------------------------
#ifdef falcON_PROPER
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
    inline void setis();
#endif
    inline double F0(double) const;
  public:
    //--------------------------------------------------------------------------
    // construction & destruction                                               
    //--------------------------------------------------------------------------
    explicit 
    SphericalSampler(double mt,                    // I: total mass             
#ifdef falcON_PROPER
		     double   =0.,                 //[I: a: anisotropy radius]  
		     double   =0.,                 //[I: b0]                    
		     const double* =0,             //[I: mass adaption: radii]  
		     int      =0,                  //[I: mass adaption: # -- ]  
		     double   =1.2,                //[I: mass adaption: factor] 
		     bool     =0,                  //[I: mass adaption: R_-/Re] 
#endif
		     bool  c  =0)                  //[I: allow non-monotonic f] 
#ifdef falcON_PROPER
      ;
#else
      : 
    Mt(mt), ibt(1./3.), careful(c) {}
#endif
    //--------------------------------------------------------------------------
    ~SphericalSampler() { 
#ifdef falcON_PROPER
      if(num) falcON_DEL_A(num);
      if(Xe)  falcON_DEL_A(Xe);
      if(Is)  falcON_DEL_A(Is);
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
    // 2. Sampling of phase space                                               
    //                                                                          
    // Sampling radius, radial and tangential velocity from the model. The      
    // boolean template specifies whether quasi-random number are used (true)   
    // for getting the radius as well as the angles in velocity space.          
    //                                                                          
    //--------------------------------------------------------------------------
    template<bool QUASI>
    double set_radvel(                             // R: Psi(r)                 
		      Random const&,               // I: pseudo & quasi RNG     
		      double&,                     // O: radius                 
		      double&,                     // O: radial velocity        
		      double&,                     // O: tangential velocity    
		      double&) const;              // O: f(E) or g(Q)           
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
    void sample(body   const&,                     // I: first body to sample   
		unsigned     ,                     // I: # bodies to sample     
		bool         ,                     // I: quasi random?          
		Random const&,                     // I: pseudo & quasi RNG     
		double       =0.5,                 //[I: fraction with vphi>0]  
#ifdef falcON_PROPER
	        double       =0.0,                 //[I: factor: setting eps_i] 
#endif
		bool         =false                //[I: write DF into aux?]    
		) const;          
    //--------------------------------------------------------------------------
    void sample(bodies const&B,                    // I/O: bodies to sample     
		bool         q,                    // I: quasi random?          
		Random const&R,                    // I: pseudo & quasi RNG     
		double       f =0.5,               //[I: fraction with vphi>0]  
#ifdef falcON_PROPER
	        double       e =0.0,               //[I: factor: setting eps_i] 
#endif
		bool         g =false              //[I: write DF into aux?]    
		) const
    {
      sample(B.begin_all_bodies(), B.N_bodies(), q,R,f,
#ifdef falcON_PROPER
	     e,
#endif
	     g);
    }
    //--------------------------------------------------------------------------
  };
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#endif// falcON_included_sample_h
////////////////////////////////////////////////////////////////////////////////


