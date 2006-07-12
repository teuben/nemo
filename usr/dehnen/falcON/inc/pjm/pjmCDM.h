// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// truncatedCDM.h                                                              |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2004                                          |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// This is a non-public part of the code.                                      |
// It is property of its author and not to be made public without his written  |
// consent.                                                                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// class for a self-consistent spherical stellar dynamical model with radial   |
// density profile                                                             |
//                                                                             |
//                  sech(r/b)                                                  |
// rho(r) = C ---------------------                                            |
//            r^alfa (r+a)^(3-alfa)                                            |
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
// alfa = inner density exponent,                                              |
// beta = anisotropy parameter,                                                |
// a    = scale- or break-radius,                                              |
// b    = truncation radius, and                                               |
// C    = normalization constant                                               |
//                                                                             |
// C is set such that the maximum circular velocity equals a given value.      |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_truncatedCDM_h
#define falcON_included_truncatedCDM_h

////////////////////////////////////////////////////////////////////////////////
namespace {
  class trunc_factor;                             // forward declaration        
  class HaloModel;                                // forward declaration        
} 
namespace falcON {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::TruncCDMModel                                              //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class TruncCDMModel {
  public:
    enum truncfac {
      expn = 0,                                    // exp(-r/b)                 
      sech = 1,                                    // sech(r/b)                 
      suph = 2                                     // 2/[cosh(r/b) + sech(r/b)] 
    };
    //--------------------------------------------------------------------------
    // private data                                                             
    //--------------------------------------------------------------------------
  private:
    double        A,B,c,a,b;                       // defining parameters       
    truncfac      tf;
    HaloModel    *RHO,*RED;                        // auxiliary class           
    int           n,n1,nm,k,k1;                    // size of tables            
    double       *r,*lr,*m,*ps,*sq,*lg,*ec;        // tables                    
    double        ia,atb,ps0,sq0,vqm,wtot,ktot;    // derived parameters        
    double        ab,ac,ca,sca,cac,fac,ln_a,ln_fac;// scaling factors           
    //--------------------------------------------------------------------------
    // private methods                                                          
    //--------------------------------------------------------------------------
    void  setgE();                                 // tabulate g(Eps)           
    double lngE(double) const;                     // log(g(Eps))               
    double gofE(double) const;                     // g(Eps)                    
    double psix(double) const;                     // Psi(x=r/a)                
    double epcx(double) const;                     // Eps_c(x=r/a)              
    double xepc(double) const;                     // x_c(Eps/ca)               
    double rhox(double) const;                     // rho(x=r/a)                
    double redx(double) const;                     // x^2beta * rho(x=r/a)      
    double cumx(double) const;                     // M(<x=r/a)                 
    double xcum(double) const;                     // x(M/c)                    
    double sqrx(double) const;                     // sigma_r^2(x)              
    double gamx(double) const;                     // gamma(x=r/a)              
    double xap (double, double, double) const;     // x(e,lq,cos(eta))          
    //--------------------------------------------------------------------------
    // construction & destruction                                               
    //--------------------------------------------------------------------------
  public:
    TruncCDMModel(double,                          // I: inner density exponent 
		  double,                          // I: break radius a         
		  double,                          // I: truncation radius b    
		  double,                          // I: mass                   
		  double,                          // I: anisotropy param beta  
		  truncfac = sech);                //[I: type of truncation fac]
    //--------------------------------------------------------------------------
    ~TruncCDMModel();                              // destruction               
    //--------------------------------------------------------------------------
    // const methods: global properties                                         
    //--------------------------------------------------------------------------
    double const&inner_slope () const { return A; }
    double const&anisotropy  () const { return B; }
    double const&break_radius() const { return a; }
    double const&trunc_radius() const { return b; }
    double const&mass_normal () const { return c; }
    double       central_pot () const { return -ca*ps0; }
    double       total_mass  () const { return c*m[n1]; }
    double       vcirc_square() const { return ca*vqm; }
    double       pot_energy  () const { return c*ca*wtot; }
    double       kin_energy  () const { return c*ca*ktot; }
    double       total_energy() const { return c*ca*(wtot+ktot); }
    //--------------------------------------------------------------------------
    // const methods: local properties                                          
    //--------------------------------------------------------------------------
    double rho (double) const;                     // rho(r)                    
    double psi (double) const;                     // Psi(r)                    
    double Epc (double) const;                     // Eps_c(r)                  
    double RcE (double) const;                     // Rc(Eps)                   
    double cum (double) const;                     // M(<r)                     
    double rM  (double) const;                     // r(M)                      
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
    double    rad_grid (int i) const { return a  * r[i]; }
    double    lnr_grid (int i) const { return ln_a + lr[i]; }
    double    cum_grid (int i) const { return c  * m[i]; }
    double    psi_grid (int i) const { return ca * ps[i]; }
    double    vcq_grid (int i) const { return ca * m[i] / r[i]; }
    double    sqr_grid (int i) const { return ca * sq[i]; }
    double    lng_grid (int i) const { return ln_fac + lg[i]; }
    double    epc_grid (int i) const { return ca * ec[i]; }
  };
}
////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_included_tcdm_h    

