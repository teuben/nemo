// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   inc/public/halo.h                                                   
///                                                                             
/// \brief  classes for generating spherical halo models                        
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
#ifndef falcON_included_halo_h
#define falcON_included_halo_h

#ifndef falcON_included_externacc_h
#  include <externacc.h>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace WDutils { template<class C, class B> class Pspline; }
namespace falcON {
  using namespace WDutils;
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::HaloDensity                                                  
  //                                                                            
  /// abstract base class for a spherical halo model                            
  /// used in the construction of halo equilibrium and initial conditions       
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class HaloDensity {
  public:
    /// pure virtual function: density at given radius
    /// \return density
    /// \param r (input) radius
    virtual double operator()(double r) const = 0;
    /// pure virtual function: negative logarithmic density slope at r->0
    /// \return gamma(r->0)
    virtual double inner_gamma() const = 0;
    /// pure virtual function: scale radius
    /// \return scale radius
    virtual double scale_radius() const = 0;
    /// pure virtual function: truncation radius (if any)
    /// \return truncation radius, zero if density decays like power law
    virtual double trunc_radius() const = 0;
    /// pure virtual function: negative logarithmic density slope at r->oo
    /// \return gamma(r->oo)
    virtual double outer_gamma() const = 0;
    /// pure virtual function: density and its first derivative
    /// \return density
    /// \param r (input) radius
    /// \param rh1 (output) drho/dr
    virtual double operator()(double r, double&rh1) const = 0;
    /// pure virtual function: density and its first two derivatives
    /// \return density
    /// \param r (input) radius
    /// \param rh1 (output) drho/dr
    /// \param rh2 (output) d^2rho/dr^2
    virtual double operator()(double r, double&rh1, double&rh2) const = 0;
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::HaloPotential                                                
  //                                                                            
  /// given a HaloDensity, we find the potential generated by it and an optional
  /// external monopole potential.                                              
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class HaloPotential {
  protected:
    const HaloDensity &DEN;                        ///< density of halo         
    const acceleration*MON;                        ///< external monopole       
    double             Ah,At;                      ///< gamma_0 for rho, psi    
    double             Ch;                         ///< gamma_00 for rho        
    double             ps0;                        ///< Psi(0) (if finite)      
    int                n,n1,nm;                    ///< size of tables          
    Array<double,1>    r,lr;                       ///< r_i, ln(r_i)            
    Array<double,1>    mh;                         ///< M_halo(<r_i)            
    Array<double,1>    mt;                         ///< M_tot(<r_i)             
    Array<double,1>    ps;                         ///< Psi_tot(r_i)            
    Array<double,1>    rh;                         ///< rho_tot(r_i)            
    Array<double,1>    ec;                         ///< E_c(r_i)                
    Array<double,1>    dp;                         ///< -v_c^2(r_i) = dPsi/dlnr 
    Pspline<double,double>   *PS;                  ///< penta spline: Psi(r)    
  public:
    /// constructor
    /// \param halo density model for halo
    /// \param mono external monopole potential
    HaloPotential(HaloDensity const&halo, const acceleration*mono);
    /// destructor
    ~HaloPotential();
    /// potential Psi(r) using polynomial interpolation
    double Ps(double R) const;
    /// potential Phi(r), and (-dPhi/dr)/r using penta spline
    /// \return Phi(r)
    /// \param  rq  (input) radius squared
    /// \param  acx (output) (-dPhi/dr)/r, so that acc = pos * acx
    double PotAcc(double rq, double&acx) const;
    /// total cumulative mass
    double Mt(double R) const;
    /// halo cumulative mass
    double Mh(double R) const;
    /// total halo mass
    double Mh_tot() const { return mh[n1]; }
    /// total mass density
    double rhot(double R) const;
    /// v_circ^2(r)
    double vcq(double R) const;
    /// Omega^2(r)
    double omq(double R) const;
    /// kappa^2(r)
    double kpq(double R) const;
    /// gamma := 2*Omega/kappa 
    double gam(double R) const;
    // ln R_psi(E)
    double lnRPsi(double P) const;
    // R_psi(E)
    double RPsi(double P) const;
    /// Eps_c(r)
    double Epc(double R) const;
    /// R_circ(Eps)
    double RcE(double E) const;
    /// estimate for R(Eps, L^2, cos[eta])
    double Rap(double E, double Lq, double ce) const;
    /// estimate for R_peri(Eps, L^2)
    double Rp(double E, double Lq) const;
    /// estimate for R_apo(Eps, L^2)
    double Ra(double E, double Lq) const;
    /// radius given cumulative halo mass
    double RMh(double M) const;
    //--------------------------------------------------------------------------
    /// \name constant access to some data
    //{@
    /// size of tables
    int Ntab() const { return n; }
    /// r_i
    double const&rad(int i) const { return r[i]; }
    /// ln(r_i)
    double const&lnr(int i) const { return lr[i]; }
    /// M_halo(r<r_i)
    double const&mhr(int i) const { return mh[i]; }
    /// Psi_tot(r_i)
    double const&psi(int i) const { return ps[i]; }
    /// M_tot(r<r_i)
    double const&mtr(int i) const { return mt[i]; }
    /// Eps_circ(r_i)
    double const&epc(int i) const { return ec[i]; }
    //@}
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::HaloModel                                                    
  //                                                                            
  /// given a HaloDensity, we construct a full spherical halo model with        
  /// Cuddeford (1991) distribution function (default: isotropic) in the        
  /// presence of an external spherical potential.                              
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class HaloModel : public HaloPotential {
  private:
    const double       B;                          ///< beta parameter of DF    
    const double       RA;                         ///< anisotropy radius of DF 
    Array<double,1>    lg;                         ///< table: log(g(Q))        
  public:
    /// constructor
    /// \param halo density model for halo
    /// \param beta anisotropy parameter beta 
    /// \param r_a anisotropy radius (0 maps to infinity)
    /// \param mono external monopole potential
    HaloModel(HaloDensity const&halo,
	      double beta, double r_a,
	      const acceleration*mono);
    /// ln g(Q)
    double lnG(double Q) const;
    /// distribution function f(Eps,L^2)
    double fEL(double E, double Lq) const;
    /// ln g(Q=Psi_tot(r_i))
    double const&lng(int i) const { return lg[i]; }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::HaloModifier                                                 
  //                                                                            
  /// given a HaloDensity rho(r), we modify it to rho(sqrt[r^2+c^2]) sech(r/rt) 
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class HaloModifier {
    const double rc,rcq;              ///< core radius, core radius squared
    const double rt,irt;              ///< truncation radius & its inverse
  public:
    /// sqrt(r^2+rc^2)
    double core(double r) const;
    /// sqrt(r^2+rc^2) and its first derivative
    double core(double r, double&d1) const;
    /// sqrt(r^2+rc^2) and its first two derivatives
    double core(double r, double&d1, double&d2) const;
    /// cored but untrancated density at given radius
    double cored(HaloDensity const&m, double r) const;
    /// cored but untrancated density and its first derivative
    double cored(HaloDensity const&m, double r, double&rh1) const;
    /// cored but untrancated density and its first two derivatives
    double cored(HaloDensity const&m, double r, double&rh1, double&rh2) const;
    /// truncation factor
    double trunc(double r) const;
    /// truncation factor and its first derivative
    double trunc(double r, double&t1) const;
    /// truncation factor and its first two derivatives
    double trunc(double r, double&t1, double&t2) const;
    /// core radius
    double const&r_c() const { return rc; }
    /// truncation radius
    double const&r_t() const { return rt; }
    /// constructor
    /// \param c core radius
    /// \param t truncation radius
    HaloModifier(double c, double t) falcON_THROWING;
    /// modified density at given radius
    double operator()(HaloDensity const&m, double r) const;
    /// modified density and its first derivative
    double operator()(HaloDensity const&m, double r, double&rh1) const;
    /// modified density and its first two derivatives
    double operator()(HaloDensity const&m, double r, double&rh1, double&rh2)
      const;
    /// is the core radius actually non-zero?
    bool cored() const { return rcq; }
    /// is the truncation radius finite?
    bool truncated() const { return irt; }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::DoublePowerLawHalo                                           
  //                                                                            
  /// implements a halo density model with density\n                            
  /// rho(r) = r^(-gamma_i) * (1 + r^eta)^((gamma_o-gamma_i)/eta)               
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class DoublePowerLawHalo : public HaloDensity {
    const double go,gi,et;
    const double gg,al;
  public:
    /// constructor
    /// \param inner inner power law slope gamma_i of density
    /// \param outer outer power law slope gamma_o of density
    /// \param trans transition steepness eta
    DoublePowerLawHalo(double inner, double outer, double trans);
    /// negative logarithmic density slope at r->0
    double inner_gamma() const { return gi; }
    /// scale radius: one
    double scale_radius() const { return 1.; }
    /// transition steepness
    double transition() const { return et; }
    /// no truncation radius
    double trunc_radius() const { return 0.; }
    /// negative logarithmic density slope at r->oo
    double outer_gamma() const { return go; }
    /// density at given (scaled) radius
    double operator()(double x) const;
    /// density and its first derivative
    double operator()(double x, double&rh1) const;
    /// density and its first two derivatives
    double operator()(double x, double&rh1, double&rh2) const;
    /// total mass if multiplied with a truncation factor and used with core
    /// \param M modifier, specifying truncation and core
    double Mtot(const HaloModifier&M) const;
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::ModifiedDoublePowerLawHalo                                   
  //                                                                            
  /// implements a halo with truncated double power-law density                 
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class ModifiedDoublePowerLawHalo : public HaloDensity {
    const DoublePowerLawHalo Model;
    const HaloModifier       Modif;
    const double             rsc,irs,mt,rh0,fc1,fc2;
  public:
    /// constructor
    /// It is required that the total mass is finite; otherwise, we error out.
    /// The total mass is finite if outer > 3 or trunc > 0.
    /// \param scale scale radius
    /// \param core  core radius
    /// \param trunc truncation radius (zero -> infinity, ie. no truncation)
    /// \param mtot  total mass
    /// \param inner inner power law slope of density
    /// \param outer outer power law slope of density
    /// \param trans transition steepness
    ModifiedDoublePowerLawHalo(double scale, double core, double trunc,
			       double mtot,
			       double inner, double outer, double trans)
      falcON_THROWING :
      Model(inner,outer,trans),
      Modif(core/scale,trunc/scale),
      rsc  (scale),
      irs  (1/scale),
      mt   (mtot),
      rh0  (cube(irs)*mtot/Model.Mtot(Modif)),
      fc1  (rh0*irs),
      fc2  (fc1*irs) 
    {
      if(isinf(rsc)) falcON_THROW("ModifiedDoublePowerLawHalo: r_s=inf\n");
      if(rsc   ==0.) falcON_THROW("ModifiedDoublePowerLawHalo: r_s==0\n");
      if(rsc   < 0.) falcON_THROW("ModifiedDoublePowerLawHalo: r_s=%g<0\n",rsc);
      if(isinf(mt) ) falcON_THROW("ModifiedDoublePowerLawHalo: M=inf\n");
      if(mt    ==0.) falcON_THROW("ModifiedDoublePowerLawHalo: M==0\n");
      if(mt     <0.) falcON_THROW("ModifiedDoublePowerLawHalo: M=%g<0\n",mt);
      debug_info(2,"ModifiedDoublePowerLawHalo: rh0=%f\n",rh0);
     }
    /// total mass
    double const&total_mass() const {
      return mt;
    }
    /// density normalisation
    double const&rho0() const {
      return rh0;
    }
    /// negative logarithmic density slope at r->0
    double inner_gamma() const {
      return Modif.cored()? 0. : Model.inner_gamma();
    }
    /// scale radius
    double scale_radius() const {
      return rsc;
    }
    /// core radius
    double core_radius() const {
      return rsc*Modif.r_c();
    }
    /// truncation radius
    double trunc_radius() const {
      return rsc*Modif.r_t();
    }
    /// negative logarithmic density slope at r->oo
    double outer_gamma() const {
      return Modif.truncated()? 100. : Model.outer_gamma();
    }
    /// transition steepness
    double transition() const {
      return Model.transition();
    }
    /// density at given (scaled) radius
    double operator()(double r) const {
      return rh0*Modif(Model,r*irs);
    }
    /// density and its first derivative
    double operator()(double r, double&rh1) const {
      double rho=rh0*Modif(Model,r*irs,rh1);
      rh1 *= fc1;
      return rho;
    }
    /// density and its first two derivatives
    double operator()(double r, double&rh1, double&rh2) const {
      double rho=rh0*Modif(Model,r*irs,rh1,rh2);
      rh1 *= fc1;
      rh2 *= fc2;
      return rho;
    }
  };
////////////////////////////////////////////////////////////////////////////////
} // namespace falcON
falcON_TRAITS(falcON::HaloDensity,"falcON::HaloDensity");
falcON_TRAITS(falcON::HaloModel,"falcON::HaloModel");
falcON_TRAITS(falcON::HaloModifier,"falcON::HaloModifier");
falcON_TRAITS(falcON::DoublePowerLawHalo,"falcON::DoublePowerLawHalo");
falcON_TRAITS(falcON::ModifiedDoublePowerLawHalo,"falcON::ModifiedDoublePowerLawHalo");
////////////////////////////////////////////////////////////////////////////////
#endif
