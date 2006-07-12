// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// anyhalo.h                                                                   |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2005                                          |
//           Paul McMillan, 2004-2005                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
//           paul.mcmillan@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// class for a self-consistent spherical stellar dynamical model with radial   |
// density profile defined by user (within the file, for the time being)       |
//                                                                             |
//                                                                             |
// and distribution function                                                   |
//                                                                             |
//                      2 -beta                                                |
//                  / (L )       for beta < 1,                                 |
// f(E,L) = g(E) x |         2                                                 |
//                  \ delta(L )  for beta = 1,                                 |
//                                                                             |
// where:                                                                      |
//                                                                             |
// beta = anisotropy parameter,                                                |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_anyhalo_h
#define falcON_included_anyhalo_h
 
#ifndef falcON_included_externacc_h
#  include<externacc.h>
#endif

////////////////////////////////////////////////////////////////////////////////
namespace {
  class HaloModel;                                // forward declaration        
} 
namespace falcON {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::AnyHalo                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class AnyHalo {
    //--------------------------------------------------------------------------
    // private data                                                             
    //--------------------------------------------------------------------------
  private:
    double        A,B,RA,c,MM;                        // defining parameters       
    HaloModel    *RHO;                             // auxiliary class           
    int           n,n1,nm,k,k1,tester;             // size of tables            
    double       *r,*lr,*rhe,*mg,*mh,*ps,*sq,*lg,*ec;  // tables                
    double        ia,atb,ps0,sq0,vqm,wtot,ktot,brokeps;// derived parameters    
    double        ic,fac,ln_a,ln_fac;              // scaling factors
    const acceleration  *external;                 // External monopole      
    //--------------------------------------------------------------------------
    // private methods                                                          
    //--------------------------------------------------------------------------
    void  setgE();                                 // tabulate g(Eps)           
    double lngE(double) const;                     // log(g(Eps))               
    double gofE(double) const;                     // g(Eps)                    
    double psix(double) const;                     // Psi(x=r/a)                
    double epcx(double) const;                     // Eps_c(x=r/a)              
    double xepc(double) const;                     // x_c(Eps/ca)               
    double rhox_g(double) const;                   // rho(x=r/a) grav           
    double rhox_h(double) const;                   // rho(x=r/a) halo
    double rhox_e(double) const;                   // rho(x=r/a) monopole
    double redx(double) const;                     // x^2beta * rho(x=r/a)      
    double cumx_g(double) const;                   // M(<x=r/a) grav            
    double cumx_h(double) const;                   // M(<x=r/a) halo            
    double xcum_h(double) const;                   // x(M/c)                    
    double sqrx(double) const;                     // sigma_r^2(x)              
    double gamx(double) const;                     // gamma(x=r/a)              
    double xap (double, double, double) const;     // x(e,lq,cos(eta))          
    typedef tupel<3,double> vectd;
    //--------------------------------------------------------------------------
    // construction & destruction                                               
    //--------------------------------------------------------------------------
  public:
    AnyHalo      (double,                          // I: Mdisc/Mhalo         
		  double,                          // I: inner density exponent
		  double,                          // I: outer density exponent
		  double,                          // I: scale radius        
		  double,                          // I: trunc. radius       
		  double,                          // I: anisotropy param beta  
		  double,                          // I: anisotropy radius
		  const acceleration * = 0);       //[I: external monopole]    
    //--------------------------------------------------------------------------
    ~AnyHalo();                              // destruction               
    //--------------------------------------------------------------------------
    // const methods: global properties                                         
    //--------------------------------------------------------------------------
    double const&inner_slope () const { return A; }
    double const&anisotropy  () const { return B; }
    double const&mass_normal () const { return c; }
    double       central_pot () const { return -c*ps0; }
    double       total_mass  () const { return c*mh[n1]; }
    double       vcirc_square() const { return c*vqm; }
    double       pot_energy  () const { return c*c*wtot; }
    double       kin_energy  () const { return c*c*ktot; }
    double       total_energy() const { return c*c*(wtot+ktot); }
    //--------------------------------------------------------------------------
    // const methods: local properties                                          
    //--------------------------------------------------------------------------
    double rho_g (double) const;                   // rho(r) of gravitating mass
    double rho_h (double) const;                   // rho(r) of halo only       
    double cum_g (double) const;                   // gravitating M(<r)         
    double cum_h (double) const;                   // halo's M(<r)              
    double rM_h  (double) const;                   // r(M_halo)                
    double psi (double) const;                     // Psi(r)                    
    double Epc (double) const;                     // Eps_c(r)                  
    double RcE (double) const;                     // Rc(Eps)                   
    double vcq (double) const;                     // vc^2(r)                   
    double omq (double) const;                     // omega^2(r)                
    double kpq (double) const;                     // kappa^2(r)                
    double gam (double) const;                     // gamma(r) = 2*omage/kappa  
    double sqr (double) const;                     // sigma_r^2(r)              
    double sqt (double) const;                     // sigma_theta^2(r)          
    double g_E (double) const;                     // g(Eps)                    
    double fEL (double, double) const;             // f(Eps,L^2)                
    //--------------------------------------------------------------------------
    // estimate peri and apocentre of an orbit                                  
    //--------------------------------------------------------------------------
    double Rp  (double, double) const;             // R_peri(Eps,L^2)           
    double Ra  (double, double) const;             // R_apo (Eps,L^2)           
    //--------------------------------------------------------------------------
    // direct const access to grids                                             
    //--------------------------------------------------------------------------
    int const&N_grid   ()      const { return n; }
    double    rad_grid (int i) const { return r[i]; }
    double    lnr_grid (int i) const { return lr[i]; }
    double    cum_grid (int i) const { return c  * mh[i]; }
    double    psi_grid (int i) const { return c * ps[i]; }
    double    vcq_grid (int i) const { return c * mg[i] / r[i]; }
    double    sqr_grid (int i) const { return c * sq[i]; }
    double    lng_grid (int i) const { return ln_fac + lg[i]; }
    double    epc_grid (int i) const { return c * ec[i]; }
  };
}
////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_included_anyhalo_h    

