// *-* C++ *-*                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// GalPot_in.h                                                                 |
//                                                                             |
// C++ code written by Walter Dehnen, 1996-2002,                               |
//                                                                             |
// present address:                                                            |
// Astrophysikalisches Institut Potsdam                                        |
// An der Sternwarte 16, D-14482 Potsdam, Germany                              |
// e@mail: wdehnen@aip.de                                                      |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// file containing the inline functions of class GalaxyPotential               |
//                                                                             |
// TO BE INCLUDED BY GalPot.h WHICH IN TURN IS INCLUDED BY THE ENDUSER         |
//                                                                             |
// Version 0.0    15. July      1997                                           |
// Version 0.1    24. March     1998                                           |
// Version 0.2    22. September 1998                                           |
// Version 0.3    07. June      2001                                           |
// Version 0.4    22. April     2002                                           |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef GalPot_in
#define GalPot_in 1

#include "GalPot.h"	                          // make sure GalPot.h is known
namespace GalPot {                                // v0.4                       
  //////////////////////////////////////////////////////////////////////////////
  // Constructors and related (all inline)                                      
  //////////////////////////////////////////////////////////////////////////////
  inline double GalaxyPotential::Residual(const double a,
					  const double b,
					  const double c) const {
    return Disks::Residual(a,b,c) + Spheroids::Residual(a,b,c);
  }
  //----------------------------------------------------------------------------
  inline GalaxyPotential::GalaxyPotential(std::istream &from)
    : Disks(from),
      Spheroids(from),
      M(NRAD,RMIN,RMAX,Spheroids::gamma(),Spheroids::beta(),this)
  {}
  //----------------------------------------------------------------------------
  inline GalaxyPotential::GalaxyPotential(const int Nd, const DiskPar* pd,
					  const int Ns, const SphrPar* ps,
					  const double rmin, const double rmax,
					  const int k)
    : Disks(Nd,pd),
      Spheroids(Ns,ps),
      M(k,rmin,rmax,Spheroids::gamma(),Spheroids::beta(),this)
  {}
  //----------------------------------------------------------------------------
  inline void GalaxyPotential::reset(const int Nd, const DiskPar* pd,
				     const int Ns, const SphrPar* ps,
				     const double rmin, const double rmax) {
    Disks::reset(Nd,pd);
    Spheroids::reset(Ns,ps);
    M.reset(rmin,rmax,Spheroids::gamma(),Spheroids::beta(),this);
  }
  //////////////////////////////////////////////////////////////////////////////
  // information on the disks alone (all inline)                                
  //////////////////////////////////////////////////////////////////////////////
  inline double GalaxyPotential::DisksDensity(const double R, const double z)
    const {
    return Disks::Density(R,z);
  }
  //----------------------------------------------------------------------------
  inline double GalaxyPotential::DisksSurfaceDensity(const double R) const {
    return Disks::SurfaceDensity(R);
  }
  //----------------------------------------------------------------------------
  inline double GalaxyPotential::DisksMass(const double R) const {
    return Disks::Mass(R);
  }
  //----------------------------------------------------------------------------
  inline int GalaxyPotential::NumberofDisks() const {
    return Disks::NumberofDisks();
  }
  //----------------------------------------------------------------------------
  inline DiskPar GalaxyPotential::DiskParameter(const int i) const {
    return Disks::Parameter(i);
  }
  //////////////////////////////////////////////////////////////////////////////
  // information on the spheroids alone (all inline)                            
  //////////////////////////////////////////////////////////////////////////////
  inline double GalaxyPotential::SpheroidsDensity(const double R, const double z)
    const {
    return Spheroids::Density(R,z);
  }
  //----------------------------------------------------------------------------
  inline double GalaxyPotential::SpheroidsMass(const double R) const {
    return Spheroids::Mass(R);
  }
  //----------------------------------------------------------------------------
  inline int GalaxyPotential::NumberofSpheroids() const {
    return Spheroids::NumberofSpheroids();
  }
  //----------------------------------------------------------------------------
  inline SphrPar GalaxyPotential::SpheroidParameter(const int i) const {
    return Spheroids::Parameter(i);
  }
  //////////////////////////////////////////////////////////////////////////////
  // information on the total (most are non-inline functions)                   
  //////////////////////////////////////////////////////////////////////////////
  inline double GalaxyPotential::Density(const double R, const double z) const {
    return Disks::Density(R,z) + Spheroids::Density(R,z);
  }
  //----------------------------------------------------------------------------
  inline double GalaxyPotential::Mass(const double R) const {
    return Disks::Mass(R) + Spheroids::Mass(R);
  }
  //----------------------------------------------------------------------------
  inline double GalaxyPotential::vcsquare(const double R) const {
    return M.vcsquare(R);
  }
}                                                     // namespace GalPot       
#endif                                                // GalPot_in_h            
// end of file GalPot_in.h /////////////////////////////////////////////////////
