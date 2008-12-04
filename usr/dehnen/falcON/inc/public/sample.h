// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   inc/public/sample.h                                                 
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2000-2007                                                           
///                                                                             
/// \brief  code for sampling spherical stellar dynamical equilibria            
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
// 02/05/2007 WD  made support for Cuddeford (1991) DF public                   
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
  /// provides routines for sampling phase-space densities.                     
  /// We consider DFs of the general form\n                                     
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
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class SphericalSampler {
    typedef tupel<2,double> pair_d;
  private:
    const bool    careful,OM,beta;
    const double  Mt,ra,iraq,b0,ibt;
    Array<pair_d> Xe;
    Array<double> Is;
    //--------------------------------------------------------------------------
#ifdef falcON_PROPER
    const bool    adapt_masses, Peri;
    const double  irs,eta,mmm,nmax;
#endif
    //--------------------------------------------------------------------------
    // private methods                                                          
    //--------------------------------------------------------------------------
    inline void setis();
    inline double F0(double) const;
  public:
    //--------------------------------------------------------------------------
#ifdef falcON_PROPER
    /// ctor
    /// \param mt total mass to be sampled
    /// \param ra Ossipkov (1979) - Merritt (1985) anisotropy radius
    /// \param b0 Cuddeford (1991) central anisotropy
    /// \param rs mass adaption: scale radiu
    /// \param mm mass adaption: m_max/m_min
    /// \param et mass adaption: shape parameter
    /// \param nm mass adaption: n_max per (E,L) (default: mm)
    /// \param pr mass adaption: using R_peri (or R_circ) for mass adaption ?
    /// \param cf be extra careful when g(Q) may be non-monotonic.
    explicit 
    SphericalSampler(double mt, double ra, double b0,
		     double rs, double mm=1, double et=1, double nm=0,
		     bool pr=0, bool cf=0) falcON_THROWING;
#endif
    /// ctor
    /// \param mt total mass to be sampled
    /// \param ra Ossipkov (1979) - Merritt (1985) anisotropy radius
    /// \param b0 Cuddeford (1991) central anisotropy
    /// \param cf be extra careful when g(Q) may be non-monotonic.
    explicit 
    SphericalSampler(double mt, double ra=0., double b0=0., bool cf=0);
    //--------------------------------------------------------------------------
    /// \name abstract methods required for sampling radius and velocity
    //@{
  protected:
    /// g(Q) (reduces tof(E) for isotropuc case)
    virtual double DF(double) const=0;
    /// Psi(R)
    virtual double Ps(double) const=0;
    /// r(M), M in [0,M_total]
    virtual double rM(double) const=0;
#ifdef falcON_PROPER
    /// R_circ(Eps)
    virtual double Re(double) const=0;
    /// R_peri(Eps,L)
    virtual double Rp(double, double) const=0;
#endif
    //@}
  public:
    //--------------------------------------------------------------------------
    /// sampling of full phase-space: non-virtual
    /// \note assuming the randomly generated (r,vr,vt) are already properly
    ///       scaled as well as the Psi returnd. Also the routines for
    ///       R_peri and R_circ are assummed to take arguments and return
    ///       values that are both scaled.
    /// \param[in] B0  first body to sample
    /// \param[in] N   number of bodies to sample
    /// \param[in] q   use quasi (or pseudo) random numbers?
    /// \param[in] R   quasi and pseudo random number generator
    /// \param[in] f   fraction with vphi>0
#ifdef falcON_PROPER
    /// \param[in] e   factor: setting eps_i
#endif
    /// \param[in] gF  write DF into aux?
    /// \param[in] gP  write Phi into pot?
    /// \param[in] gA  write -dPhi/dr into acc?
    void sample(body const&B0, unsigned N, bool q, Random const&R, double f=0.5,
#ifdef falcON_PROPER
	        double e=0.0,
#endif
		bool gF=false, bool gP=false, bool gA=false)
      const falcON_THROWING;
    /// sampling of full phase-space: non-virtual
    /// \note assuming the randomly generated (r,vr,vt) are already properly
    ///       scaled as well as the Psi returnd. Also the routines for
    ///       R_peri and R_circ are assummed to take arguments and return
    ///       values that are both scaled.
    /// \param[in] B   bodies to sample
    /// \param[in] q   use quasi (or pseudo) random numbers?
    /// \param[in] R   quasi and pseudo random number generator
    /// \param[in] f   fraction with vphi>0
#ifdef falcON_PROPER
    /// \param[in] e   factor: setting eps_i
#endif
    /// \param[in] gF  write DF into aux?
    /// \param[in] gP  write Phi into pot?
    /// \param[in] gA  write -dPhi/dr into acc?
    void sample(bodies const&B, bool q, Random const&R, double f=0.5,
#ifdef falcON_PROPER
	        double e=0.0,
#endif
		bool gF=false, bool gP=false, bool gA=false)
      const falcON_THROWING
    {
      sample(B.begin_all_bodies(), B.N_bodies(), q,R,f,
#ifdef falcON_PROPER
	     e,
#endif
	     gF,gP,gA);
    }
    //--------------------------------------------------------------------------
    /// sampling positions only: non-virtual
    /// \param[in] B0  first body to sample
    /// \param[in] N   number of bodies to sample
    /// \param[in] q   use quasi (or pseudo) random numbers?
    /// \param[in] R   quasi and pseudo random number generator
    void sample_pos(body const&B0, unsigned N, bool q, Random const&R)
      const falcON_THROWING;
    /// sampling positions only: non-virtual
    /// \param[in] B   bodies to sample
    /// \param[in] q   use quasi (or pseudo) random numbers?
    /// \param[in] R   quasi and pseudo random number generator
    void sample_pos(bodies const&B, bool q, Random const&R)
      const falcON_THROWING
    {
      sample_pos(B.begin_all_bodies(), B.N_bodies(), q,R);
    }
    /// noop dtor
    virtual~SphericalSampler() {}
    //--------------------------------------------------------------------------
  };// class SphericalSampler
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#endif// falcON_included_sample_h
////////////////////////////////////////////////////////////////////////////////


