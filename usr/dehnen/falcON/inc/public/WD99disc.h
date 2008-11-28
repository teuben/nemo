// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   inc/public/WD99disc.h                                               
///                                                                             
/// \brief  generating initial conditions from Dehnen (1999) disc model         
///                                                                             
/// \author Paul McMillan                                                       
/// \author Walter Dehnen                                                       
/// \date   2000-2007                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2000-2007  Walter Dehnen, Paul McMillan                        
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
#ifndef falcON_included_WD99disc_h
#define falcON_included_WD99disc_h

#ifndef falcON_included_body_h
#  include <body.h>
#endif
#ifndef falcON_included_random_h
#  include <public/random.h>
#endif
#ifndef falcON_included_externacc_h
#  include <externacc.h>
#endif

namespace falcON {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class nbdy::WD99disc                                                       
  //                                                                            
  /// class for sampling a Dehnen (1999) disc.                                  
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class WD99disc {
  private:
    int          n1,n;                          ///< Table sizes, other
    double      *lr,*pot,*dpdr;                 ///< tables: log(r),Psi,dPsi/dr
    double       sig0;                          ///< parameter
    const int    No;                            ///< parameter
    const double Rd,iRd,Dens0,Rsig,Qmin,Z0,Eps,Mt,rmin,rmax; ///< paramters
    const        acceleration *Pe;              ///< External combi-potential
    const        ExpDisk Disc;                  ///< disc for particles

    //--------------------------------------------------------------------------
    /// sub-class used for all calculations for an individual orbit in the disc
    ///
    /// It is given two inputs specific for any particular orbit, the radius A 
    /// and Xi (used to determine the eccentricity), and integrates position 
    /// and velocity in R over a full period of radial motion.\n
    ///
    /// Once this has been done, the function WD99disc::PlanarOrbit::sample()
    /// provides a randomly sampled point on the orbit (with randomised phi).
    //--------------------------------------------------------------------------
    class PlanarOrbit {
      friend class WD99disc;
      bool   range;
      int    i,Npoints;
      double re,xi,rs,s0,ire,POT,CENACC,Q,sdens;
      double Ec,Lc,Omc,kap,gam,sigre,Lorb,e,g2,omr,Tr,dT;
      double*ttable,*Rtable,*vRtable;
      vect_d*W;
      vect_d Dt(vect_d const&);
      vect_d DR(vect_d const&);
      void step_back_to_R_equal(double, vect_d& );
      double CashKarp(vect_d&, double, double );
      void CashKarpStep(vect_d&, double);
    public:
      /// ctor: integrates orbit and sets tables
      PlanarOrbit(double R, 
		  double Xi, 
		  double Rs, 
		  double S0, 
		  double Qmin, 
		  double Sdens,
		  double Sigcorr);
      /// samples a point on the orbit
      /// \param RNG random number generator
      /// \param q quasi random numbers?
      /// \param rad (output) radius
      /// \param vrad (output) radial velocity
      /// \param phi (output) azimuth
      /// \param vphi (output) azimuthal velocity
      void sample(Random const&RNG, bool q,
		  double&rad, double&vrad, double&phi, double&vphi) const;
      /// dtor: de-allocates memory
      ~PlanarOrbit();
    };
    //--------------------------------------------------------------------------
  public:
    /// ctor.
    /// Creates table of Potential and acceleration in the z=0 plane of the
    /// given external field as a function of ln(R)
    /// \param no Number of orbits
    /// \param rd Disc scale radius
    /// \param dens0 Disc scale surface density
    /// \param rsig Vdisp scale radius
    /// \param q Toomre's Q / sigma_0
    /// \param z0 scale height
    /// \param eps particle smoothing length
    /// \param pe External potential
    WD99disc(int no, double rd, double dens0, double rsig,
	     double q, double z0, double eps, const acceleration*pe);
    /// sampler: populates the disc (except in the case where Q=0).
    /// 
    /// Pick a radius R (in such a way that the final disc has the surface
    /// density desired), define an orbit with 1) energy equal to that of a
    /// circular orbit at radius R, z=0, 2) Angular momentum picked from a
    /// distribution defined so that sig_R is as defined. Place 'Nbpero' bodies
    /// on that orbit. Repeat; see Dehnen (1999) for more details.\n
    /// Since this is effectively sampling particles in energy, rather than
    /// radius, the N-body distribution found will not have exactly the desired
    /// density distribution or velocity dispersion.  An iterative scheme is
    /// used to ensure that the density distribution and velocity dispersion of
    /// the resulting N-body distribution is as close to the required as
    /// possible.
    /// \param B bodies to sample
    /// \param Ni # iterations
    /// \param q use quasi-random number?
    /// \param RNG random number generator
    /// \param g write DF to bodies?
    void sample(bodies const&B, int Ni, bool q, Random const&RNG, bool g) const;
    /// sampling a cold disc
    /// \param B bodies to sample
    /// \param q use quasi-random number?
    /// \param RNG random number generator
    void coldsample(bodies const&B, bool q, Random const&RNG) const;
    /// iterates the density distribution closer to the target
    void iterate(int tsize,
		 int Np,
		 bool q,
		 Random const&RNG,
		 double*rtar=0,
		 double*minp=0,
	         double*star=0,
	         double*sinp=0) const;
    /// dtor: de-allocate memory
    ~WD99disc();
  };
}

////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_WD99disc_h    

