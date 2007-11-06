// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// GalPot.h                                                                    |
//                                                                             |
// Copyright (C) 1996-2007 Walter Dehnen                                       |
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
//                                                                             |
// Header file for class GalaxyPotential                                       |
//                                                                             |
// TO BE INCLUDED BY THE ENDUSER                                               |
//                                                                             |
// Version 0.0    15. July      1997                                           |
// Version 0.1    24. March     1998                                           |
// Version 0.2    22. September 1998                                           |
// Version 0.3    07. June      2001                                           |
// Version 0.4    22. April     2002                                           |
// Version 0.5    05. December  2002                                           |
// Version 0.6    05. February  2003                                           |
// Version 0.7    23. September 2004  fixed "find(): x out of range" error     |
// Version 0.8    04. February  2005  debugged error in GaussLegendre()        |
//                                    used WDutils::tupel                      |
// Version 0.9    06. November  2007  brought in line with GalPot package      |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// About the unit system                                                       |
// The unit system used throughout the classes and functions defined here and  |
// in included files is based on the following basic units                     |
// 	unit of length:		1 kilo parsec (kpc)                            |
//	unit of time:		1 mega year   (Myr)                            |
// 	unit of mass:		1 solar mass  (Msun)                           |
// This implies the following dimensions                                       |
//                                                                             |
// quantity        dimension / seize                using other units          |
//-----------------------------------------------------------------------------+
// angular vel.  1 Myr^-1	                   = 977.775320024919 km/s/kpc |
// velocity      1 kpc/Myr                         = 977.775320024919 km/s     |
// action/mass   1 kpc^2/Myr                       = 977.775320024919 kpc*km/s |
// potential     1 kpc^2/Myr^2                     = 956044.576449833 (km/s)^2 |
// acceleration  1 kpc/Myr^2                                                   |
// G             4.49865897 E-12 kpc^3/Myr^2/Msun                              |
// 4 Pi G        5.65318158 E-11 kpc^3/Myr^2/Msun                              |
//-----------------------------------------------------------------------------+
#ifndef GalPot_h
#define GalPot_h

#include <iostream>

namespace GalPot {              // added in V0.4                                
  //----------------------------------------------------------------------------
  // define the maximum l used in the multipole expansion                       
  // After adjustment re-compile GalPot.cc for the change to take effect        
  //----------------------------------------------------------------------------

#ifdef GalPot_cc
  const int LMAX=80;		// maximum l for the multipole expansion        
#endif

  //----------------------------------------------------------------------------
  // Set some values used as default in the constructors below                  
  //----------------------------------------------------------------------------

  const int    NRAD=201;	// DEFAULT number of radial points in Multipole 
  const double RMIN=1.e-4,	// DEFAULT min radius of logarithmic radial grid
               RMAX=1.e3;	// DEFAULT max radius of logarithmic radial grid
  
  //----------------------------------------------------------------------------
  // Include the definitions for all auxiliary functions. These contain a long  
  // tail of dependencies needed. It seems unavoidable to have them all visible 
  // to the enduser, if GalaxyPotential below is defined as a proper class, e.g.
  // the construction of not just one object is possible etc .                  
  //                                                                            
  // Most important for the enduser is the meaning of DiskPar and SphrPar.      
  // These are Vectors of 5 and 6 doubles, respectively, holding the parameters 
  // for one disk or spheroid component. The meaning of them is as follows      
  // DiskPar[0]   is the surface density normalisation Sigma_0 [Msun/kpc^2]     
  // DiskPar[1]   is the scale length R_d [kpc]                                 
  // DiskPar[2]   is the scale height h [kpc]. For h<0 an isothermal (sech^2)   
  //		  profile is used, for h>0 an exponential one, and for h=0 the  
  //		  disk is infinitesimal thin.                                   
  // DiskPar[3]   is the inner cut-off radius R_m [kpc]                         
  // DiskPar[4]   is eps. A term eps*cos(R/R_d) is added to the exponent.       
  //                                                                            
  // SphrPar[0]   is the density normalization rho_0 [Msun/kpc^3]               
  // SphrPar[1]   is the axis ration q                                          
  // SphrPar[2]   is the inner power slope gamma                                
  // SphrPar[3]   is the outer power slope beta                                 
  // SphrPar[4]   is the transition radius r_0 [kpc]                            
  // SphrPar[5]   is the outer cut-off radius r_t [kpc]                         
  //----------------------------------------------------------------------------

  const double Grav = 4.498658966346282e-12,            // Newtons G            
               fPiG = 5.653181583871732e-11;            // 4 Pi G               
}                                                       // v0.4                 

#include "GalPot_pre.h"

namespace GalPot {                                      // v0.4                 
//------------------------------------------------------------------------------
// finally define class GalaxyPotential for the Galactic potential              
//------------------------------------------------------------------------------

  class GalaxyPotential : public  PotResidual,
			  private Disks,
			  private Spheroids {
  private:

    Multipole M;

    double Residual      (double, double, double) const;
    // this function is not intended for the enduser, it is needed in the       
    // construction of GalaxyPotential itself                                   

  public:

    virtual ~GalaxyPotential() {}

    //--------------------------------------------------------------------------
    // constructors and related functions                                       
    //--------------------------------------------------------------------------
    GalaxyPotential(std::istream&);
    // constructor from input stream. Use this constructor to establish the     
    // potential of one of the models as shown in the following code fragment.  
    //                                                                          
    // ifstream from("Model.pot"); 	// file `Model.pot' contains the data   
    // GalaxyPotenial Phi(from);	// read from file and construct object  
    // from.close();			// close file                           

    GalaxyPotential(int, const DiskPar*,	// No & parameters of disks     
		    int, const SphrPar*,	// No & parameters of spheroids 
		    double = RMIN,		// min radius of radial grid    
		    double = RMAX, 		// max radius of radial grid    
		    int    = NRAD);		// No of radial grid points     
    // constructor from parameters. See above for the meaning of the parameters 
    
    void   reset   (int, const DiskPar*,	// No & parameters of disks     
		    int, const SphrPar*, 	// No & parameters of spheroids 
		    double = RMIN, 		// min radius of radial grid    
		    double = RMAX); 		// max radius of radial grid    
    // resets to new parameters (arguments as the above constructor), the       
    // number of radial grid points remains fixed at the old value              

    //--------------------------------------------------------------------------
    // applications (are all const member functions)                            
    //--------------------------------------------------------------------------
    // information on the disks alone                                           
    //--------------------------------------------------------------------------

    double DisksDensity(double, double) const;
    // returns the disks' volume density at some (R,z) (1st & 2nd argument)     

    double DisksSurfaceDensity(double) const;
    // returns the disks' surface density at some radius (1st argument)         
    // = density integrated over z from -oo to oo.                              

    double DisksMass(double=0.) const;
    // returns the disks' mass inside some radius (1st argument)                
    // = 2*Pi*R*density integrated over z from -oo to oo, and R from 0 to R.    
    // if called with no or zero argument, the total mass is returned           

    int NumberofDisks() const;
    // returns the number of disk components                                    

    DiskPar DiskParameter(int) const;
    // returns the parameters of the ith (1st argument) disk                    

    //--------------------------------------------------------------------------
    // information on the spheroids alone                                       
    //--------------------------------------------------------------------------

    double SpheroidsDensity(double, double) const;
    // returns the spheroids' volume density at some (R,z) (1st & 2nd argument) 

    double SpheroidsMass(double) const;
    // returns the spheroids' total mass inside some radius (1st argument)      
    // = 4*Pi*q*m^2*density integrated over m from 0 to R                       

    int NumberofSpheroids() const;
    // returns the number of spheroid components                                

    SphrPar SpheroidParameter(int) const;
    // returns the parameters of the ith (1st argument) spheroid                

    //--------------------------------------------------------------------------
    // information on the total                                                 
    //--------------------------------------------------------------------------

    double Density(double, double) const;
    // returns the (input) density in Msun/kpc^3 given (R,z) in kpc.            

    double vcsquare(double) const;
    // return the circular speed squared in (kpc/Myr)^2 at R given in kpc and   
    // z=0.  vc^2 is defined by R*dPhi/dR, which can become negative, e.g. in   
    // the central parts of a hollow disk                                       

    double Mass(double) const;
    // returns the total mass inside some radius (1st argument)                 
    // this is simply the sum of DisksMass() and SpheroidsMass() above          

    void OortConstants(double, double&, double&) const;
    // given R in kpc (1st argument) returns Oort's constants A, B as 2nd and   
    // 3rd argument. Units of the latter are 1/Myr                              

    double operator() (double, double) const;
    // returns Phi in (kpc/Myr)^2 at (R,z) given in kpc                         

    double operator() (double, double, double&, double&) const;
    // returns Phi in (kpc/Myr)^2 at (R,z) given in kpc.                        
    // additionally, on return the 3rd and 4th argument contain the derivatives 
    // dPhi/dR and dPhi/dz in units of kpc/Myr^2                                

    double Laplace(double, double) const;
    // returns Laplace(Phi) in 1/Myr^2 given (R,z) in kpc. This is not          
    // necessarily identical to 4 Pi G times GalaxyPot::Density() above, as the 
    // latter returns the input density, while this routine gives the Laplace   
    // of the potential as evaluated. However, both numbers should agree in the 
    // limit of infinite many radial grid points, multipoles, and infinite      
    // numerical accuracy                                                       

    Frequs KapNuOm  (double) const;
    // given R, the epicycle frequencies kappa (radial), nu (vertical), and     
    // omega (azimuthal) for the circular orbit through (R,z=0) are returned.   
    // units for the frequencies are 1/Myr                                      

  };// class GalaxyPotential
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // inline members methods of class GalaxyPotential                            
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  // Constructors and related (all inline)                                      
  //////////////////////////////////////////////////////////////////////////////
  inline double GalaxyPotential::Residual(double a,
					  double b,
					  double c) const {
    return Disks::Residual(a,b,c) + Spheroids::Residual(a,b,c);
  }
  //----------------------------------------------------------------------------
  inline GalaxyPotential::GalaxyPotential(std::istream &from)
    : Disks(from),
      Spheroids(from),
      M(NRAD,RMIN,RMAX,Spheroids::gamma(),Spheroids::beta(),this)
  {}
  //----------------------------------------------------------------------------
  inline GalaxyPotential::GalaxyPotential(int Nd, const DiskPar* pd,
					  int Ns, const SphrPar* ps,
					  double rmin, double rmax,
					  int k)
    : Disks(Nd,pd),
      Spheroids(Ns,ps),
      M(k,rmin,rmax,Spheroids::gamma(),Spheroids::beta(),this)
  {}
  //----------------------------------------------------------------------------
  inline void GalaxyPotential::reset(int Nd, const DiskPar* pd,
				     int Ns, const SphrPar* ps,
				     double rmin, double rmax) {
    Disks::reset(Nd,pd);
    Spheroids::reset(Ns,ps);
    M.reset(rmin,rmax,Spheroids::gamma(),Spheroids::beta(),this);
  }
  //////////////////////////////////////////////////////////////////////////////
  // information on the disks alone (all inline)                                
  //////////////////////////////////////////////////////////////////////////////
  inline double GalaxyPotential::DisksDensity(double R, double z)
    const {
    return Disks::Density(R,z);
  }
  //----------------------------------------------------------------------------
  inline double GalaxyPotential::DisksSurfaceDensity(double R) const {
    return Disks::SurfaceDensity(R);
  }
  //----------------------------------------------------------------------------
  inline double GalaxyPotential::DisksMass(double R) const {
    return Disks::Mass(R);
  }
  //----------------------------------------------------------------------------
  inline int GalaxyPotential::NumberofDisks() const {
    return Disks::NumberofDisks();
  }
  //----------------------------------------------------------------------------
  inline DiskPar GalaxyPotential::DiskParameter(int i) const {
    return Disks::Parameter(i);
  }
  //////////////////////////////////////////////////////////////////////////////
  // information on the spheroids alone (all inline)                            
  //////////////////////////////////////////////////////////////////////////////
  inline double GalaxyPotential::SpheroidsDensity(double R, double z)
    const {
    return Spheroids::Density(R,z);
  }
  //----------------------------------------------------------------------------
  inline double GalaxyPotential::SpheroidsMass(double R) const {
    return Spheroids::Mass(R);
  }
  //----------------------------------------------------------------------------
  inline int GalaxyPotential::NumberofSpheroids() const {
    return Spheroids::NumberofSpheroids();
  }
  //----------------------------------------------------------------------------
  inline SphrPar GalaxyPotential::SpheroidParameter(int i) const {
    return Spheroids::Parameter(i);
  }
  //////////////////////////////////////////////////////////////////////////////
  // information on the total (most are non-inline functions)                   
  //////////////////////////////////////////////////////////////////////////////
  inline double GalaxyPotential::Density(double R, double z) const {
    return Disks::Density(R,z) + Spheroids::Density(R,z);
  }
  //----------------------------------------------------------------------------
  inline double GalaxyPotential::Mass(double R) const {
    return Disks::Mass(R) + Spheroids::Mass(R);
  }
  //----------------------------------------------------------------------------
  inline double GalaxyPotential::vcsquare(double R) const {
    return M.vcsquare(R);
  }
  //////////////////////////////////////////////////////////////////////////////
} // namespace GalPot {
#endif                                             // GalPot_h                  
// end of file GalPot.h ////////////////////////////////////////////////////////
