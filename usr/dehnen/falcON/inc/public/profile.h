// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   inc/public/profile.h                                                
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2000-2006                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2000-2006 Walter Dehnen                                        
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
#ifndef falcON_included_profile_h
#define falcON_included_profile_h

#ifndef falcON_included_body_h
#  include <body.h>
#endif

namespace falcON {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::spherical_profile                                            
  //                                                                            
  /// All bodies in_subset() are binned in radial shells to estimate radial     
  /// profiles                                                                  
  ///                                                                           
  /// Radial bins are of minimal width in ln(r) but contain at least a minimum  
  /// number of bodies. For each bin, the following quantities are estimated:\n 
  /// - the median radius within the spherical shell                         \n 
  /// - the cumulative mass, gravitational potential, and circular speed     \n 
  /// - the mean density                                                     \n 
  /// - the mean angular momentum                                            \n 
  /// - the mean rotational velocity and direction of rotation               \n 
  /// - the mean radial, azimuthal (in direction of rotation) and meridional    
  ///   velocity dispersion                                                  \n 
  /// - the axis ratios and principal axes (derived from the moment of inertia) 
  ///                                                                           
  /// \note                                                                     
  ///     In fact, we make the bins smaller by a factor of two and then consider
  ///     the merger of two adjacent smaller bins. This guarantees a somewhat   
  ///     higher resolution, but implies that adjacent larger bins share about  
  ///     half of their bodies.                                                 
  ///                                                                           
  /// \note                                                                     
  ///     used by manipulator falcON::Manipulate::spherprof                     
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class spherical_profile {
  private:
    typedef bodies::index index;                   // type used to address body 
    typedef double *pdouble;
    typedef vect_d *pvect_d;
    const int    kmin;
    const double dmax;
    int          nb,n;
    double       mt;
    pdouble      rr,mr,rh,ps,vr,sr,st,sp,ca,ba;
    pvect_d      am,vp,xa,xi;
    //--------------------------------------------------------------------------
    // construction                                                             
    //--------------------------------------------------------------------------
  public:
    /// constructor
    /// \param B bodies 
    /// \param W minimum bodies per bin
    /// \param D minimum bin width in ln(r)
    /// \param X (optional) centre position (default: origin)
    /// \param V (optional) centre velocity (default: origin)
    spherical_profile(const bodies*B,
		      unsigned     W=20,
		      double       D=0.1,
		      const vect  *X=0,
		      const vect  *V=0) falcON_THROWING;
    //--------------------------------------------------------------------------
    ~spherical_profile() {
#define DELETE_IT(P) if(P) { falcON_DEL_A(P); P=0; }
      DELETE_IT(rr);
      DELETE_IT(mr);
      DELETE_IT(rh);
      DELETE_IT(vr);
      DELETE_IT(vp);
      DELETE_IT(sr);
      DELETE_IT(st);
      DELETE_IT(sp);
      DELETE_IT(am);
      DELETE_IT(ca);
      DELETE_IT(ba);
      DELETE_IT(xa);
      DELETE_IT(xi);
#undef DELETE_IT
    }
    //--------------------------------------------------------------------------
    // data access                                                              
    //--------------------------------------------------------------------------
    /// number of bodies used
    const int   &Nb  ()      const { return nb; }
    /// total mass
    const double&Mtot()      const { return mt; }
    /// number of radial bins
    const int   &N   ()      const { return n; }
    /// do we have velocity information?
    bool         has_vels()  const { return vr!=0; }
    /// radius of bin
    const double&rad (int i) const { return rr[i]; }
    /// cumulative mass at radius
    const double&Mr  (int i) const { return mr[i]; }
    /// circular speed squared at radius
    double vcq (int i) const { return mr[i]/rr[i]; }
    /// density at radius
    const double&rho (int i) const { return rh[i]; }
    /// (negative) gravitational potential, assuming sphericity, at radius
    const double&psi (int i) const { return ps[i]; }
    /// mean radial motion at radius
    const double&vrad(int i) const { return vr[i]; }
    /// mean rotation vector (mean of er^v) at radius
    const vect_d&vrot(int i) const { return vp[i]; }
    /// mean rotational velocity at radius
    double vphi(int i) const { return abs(vp[i]); }
    /// radial velocity dispersion at radius
    const double&sigr(int i) const { return sr[i]; }        // sigma_r          
    /// velocity dispersion in rotational direction
    const double&sigp(int i) const { return sp[i]; }
    /// tangential velocity dispersion away from rotational direction
    const double&sigt(int i) const { return st[i]; }
    /// mean angular momentum in bin
    const vect_d&angm(int i) const { return am[i]; }
    /// Binney's anisotropy parameter beta
    double beta(int i) const {
      return 1.-twice(square(sr[i])/(square(st[i])+square(sp[i])));
    }
    /// minor to major axis ratio
    const double&cova(int i) const { return ca[i]; }
    /// intermediate to major axis ratio
    const double&bova(int i) const { return ba[i]; }
    /// direction of rotation axis
    vect_d drot(int i) const { 
      vect_d t(vp[i]);
      return t.normalize();
    }
    /// direction of major axis
    const vect_d&major(int i) const { return xa[i]; }
    /// direction of minor axis
    const vect_d&minor(int i) const { return xi[i]; }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::projected_profile                                            
  //                                                                            
  /// All bodies in_subset() are binned in projected radius (assuming an        
  /// infinite distance) to estimate projected profiles.                        
  ///                                                                           
  /// Bins in projected radius are of minimal width in ln(r) but contain at     
  /// least a minimum number of bodies. For each bin, the following quantities  
  /// are estimated:                                                         \n 
  /// - the median projected radius within the circular shell                \n 
  /// - the mean surface (projected) density                                 \n 
  /// - the mean line-of-sight velocity                                      \n 
  /// - the apparent rotation (los velocity feature consistent with rotation)\n 
  /// - the mean line-of-sight velocity dispersion                           \n 
  /// - the axis ratio                                                       \n 
  /// - the position angle and the rotation angle                            \n 
  ///                                                                           
  /// \note                                                                     
  ///     In fact, we make the bins smaller by a factor of two and then consider
  ///     the merger of two adjacent smaller bins. This guarantees a somewhat   
  ///     higher resolution, but implies that adjacent larger bins share about  
  ///     half of their bodies.                                                 
  ///                                                                           
  /// \note                                                                     
  ///     used by manipulator falcON::Manipulate::projprof                      
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class projected_profile {
  private:
    typedef bodies::index index;                   // type used to address body 
    typedef double *pdouble;
    const int    kmin;
    const double dmax;
    const vect_d elos;
    int          nb,n;
    double       mt;
    pdouble      mr,rr,sd,vl,vr,sl,ba,ph,al;
    //--------------------------------------------------------------------------
    // construction                                                             
    //--------------------------------------------------------------------------
  public:
    /// constructor
    /// \param B bodies
    /// \param P direction of the line-of-sight
    /// \param W minimum bodies per bin
    /// \param D minimum bin width in ln(r)
    /// \param X (optional) centre position (default: origin)
    /// \param V (optional) centre velocity (default: origin)
    projected_profile(const bodies*B,
		      vect const  &P,
		      unsigned     W=20,
		      double       D=0.1,
		      const vect  *X=0,
		      const vect  *V=0) falcON_THROWING;
    //--------------------------------------------------------------------------
    ~projected_profile() {
#define DELETE_IT(P) if(P) { falcON_DEL_A(P); P=0; }
      DELETE_IT(mr);
      DELETE_IT(rr);
      DELETE_IT(sd);
      DELETE_IT(vl);
      DELETE_IT(vr);
      DELETE_IT(sl);
      DELETE_IT(ba);
      DELETE_IT(ph);
      DELETE_IT(al);
#undef DELETE_IT
    }
    //--------------------------------------------------------------------------
    // data access                                                              
    //--------------------------------------------------------------------------
    /// number of bodies used
    const int   &Nb  ()      const { return nb; }
    /// total mass
    const double&Mtot()      const { return mt; }
    /// number of radial bins
    const int   &N   ()      const { return n; }
    /// do we have velocity information?
    bool         has_vels()  const { return vr!=0; }
    /// unit vector along line-of-sight
    vect_d const&los() const { return elos; }
    /// radius of bin
    const double&rad (int i) const { return rr[i]; }
    /// cumulative projected mass at radius
    const double&Mr  (int i) const { return mr[i]; }
    /// surface density at radius
    const double&Sig (int i) const { return sd[i]; }
    /// mean los motion at radius
    const double&vlos(int i) const { return vl[i]; }
    /// mean rotation velocity at radius
    const double&vrot(int i) const { return vr[i]; }
    /// radial los velocity dispersion at radius
    const double&sigl(int i) const { return sl[i]; }
    /// projected axis ratio
    const double&bova(int i) const { return ba[i]; }
    /// orientation angle
    const double&posa(int i) const { return ph[i]; }
    /// orientation of rotation
    const double&rota(int i) const { return al[i]; }
  };
  //////////////////////////////////////////////////////////////////////////////
} // namespace falcON {
#ifdef falcON_PROPER
#  include <proper/profile.h>
#endif
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_profile_h
