// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// profile.h                                                                   |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// This is a non-public part of the code.                                      |
// It is property of its author and not to be made public without his written  |
// consent.                                                                    |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_profile_h
#define falcON_included_profile_h

#ifndef falcON_included_body_h
#  include <body.h>
#endif

namespace falcON {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::spherical_model                                              
  //                                                                            
  // given either a set of bodies or a set of radii,                            
  // we compute/estimate assuming spherical symmetry:                           
  // - the radial cumulative mass profile (smoothing radii)                     
  // - the circular speed and frequency                                         
  // - the density and epicyclic freqency                                       
  // - the potential and the energy of the circular orbit                       
  // We then provide:                                                           
  // - direct access to the tabled values                                       
  // - interpolations from the tables                                           
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class spherical_model {
  private:
    typedef bodies::index      index;              // type used to address body 
    const int k;                                   // # neighb's on either side 
    const int n;                                   // size of tables            
    double    mt;                                  // total mass                
    double   *r,*lr,*mr,*lq,*oq,*kq,*rh,*ps,*ec;   // tables made by setup_mass 
    void setup_mass(const bodies      *,
		    Array<real>  const&,
		    Array<index> const&) falcON_THROWING;
    double   *vr,*vt,*vp,*sr,*st,*sp;              // tables made by setup_vels 
    double   *skr,*skt,*skp,*kur,*kut,*kup;        // tables made by setup_vels
    vect_d   *am;                                  // angular momentum in window
    void setup_vels(const bodies      *,
		    Array<real>  const&,
		    Array<index> const&) falcON_THROWING;
    //--------------------------------------------------------------------------
    // construction                                                             
    //--------------------------------------------------------------------------
  public:
    spherical_model(const bodies*,                 // I: bodies                 
		    int          =20,              // I: size of window         
		    bool         =1)               // I: do velocities as well? 
      falcON_THROWING;
    ~spherical_model() {
      falcON_DEL_A(r);
      falcON_DEL_A(lr);
      falcON_DEL_A(mr);
      falcON_DEL_A(lq);
      falcON_DEL_A(oq);
      falcON_DEL_A(kq);
      falcON_DEL_A(rh);
      falcON_DEL_A(ps);
      falcON_DEL_A(ec);
      if(vr) {
	falcON_DEL_A(vr);
	falcON_DEL_A(vt);
	falcON_DEL_A(vp);
	falcON_DEL_A(sr);
	falcON_DEL_A(st);
	falcON_DEL_A(sp);
	falcON_DEL_A(am);
      }
    }
    //--------------------------------------------------------------------------
    // data access                                                              
    // 1 tables setup if body positions and masses are provided                 
    //--------------------------------------------------------------------------
    const int   &N   ()      const { return n; }
    const double&rad (int i) const { return r [i]; }        // radius           
    const double&lrd (int i) const { return lr[i]; }        // ln(radius)       
    const double&Mr  (int i) const { return mr[i]; }        // M(<r)            
    const double&Omq (int i) const { return oq[i]; }        // Omega^2(r)       
    const double&kaq (int i) const { return kq[i]; }        // kappa^2(r)       
    const double&rho (int i) const { return rh[i]; }        // rho(r)           
    const double&psi (int i) const { return ps[i]; }        // Psi(r)           
    const double&Epc (int i) const { return ec[i]; }        // Eps_circ(r)      
    const double&Lcq (int i) const { return lq[i]; }        // L_circ^2(r)      
    const double&Mtot()      const { return mt; }           // M_total          
    //--------------------------------------------------------------------------
    // 2 tables setup if body velocities are provided, too                      
    //--------------------------------------------------------------------------
    bool     has_vels()      const { return vr != 0; }
    const double&vrad(int i) const { return vr[i]; }        // mean v_r         
    const double&vthe(int i) const { return vt[i]; }        // mean v_theta     
    const double&vphi(int i) const { return vp[i]; }        // mean v_phi       
    const double&sigr(int i) const { return sr[i]; }        // sigma_r          
    const double&sigt(int i) const { return st[i]; }        // sigma_theta      
    const double&sigp(int i) const { return sp[i]; }        // sigma_phi        
    const double&sker(int i) const { return skr[i]; }       // skew_vr          
    const double&sket(int i) const { return skt[i]; }       // skew_vtheta      
    const double&skep(int i) const { return skp[i]; }       // skew_vphi        
    const double&kurr(int i) const { return kur[i]; }       // kurtosis_r       
    const double&kurt(int i) const { return kut[i]; }       // kurtosis_theta   
    const double&kurp(int i) const { return kup[i]; }       // kurtosis_phi     
    const vect_d&angm(int i) const { return am[i]; }        // mean L           
    const double beta(int i) const {                        // Binney's beta    
      return 1.-twice(square(sr[i])/(square(st[i])+square(sp[i])));
    }
    //--------------------------------------------------------------------------
    // sum over several bins                                                    
    // args: I: central bin, width, O: mean and dispersion velocity             
    //--------------------------------------------------------------------------
    void binup_vr(int,int, double&,double&) const;
    void binup_vt(int,int, double&,double&) const;
    void binup_vp(int,int, double&,double&) const;
  private:
    void binup_v(int,int, const double*,const double*, double&,double&) const;
    //--------------------------------------------------------------------------
    // table interpolations                                                     
    //                                                                          
    // Note: we will not check for out-of-range. In such a case, table          
    // extrapolation will be used, but may give non-sense results.              
    // You have been warned!                                                    
    //--------------------------------------------------------------------------
  public:
    enum arg {
      radial,                                               // f = f(r)         
      lograd,                                               // f = f(log[r])    
      psient,                                               // f = f(Psi)       
      epcirc,                                               // f = f(Eps_circ)  
      lcircq,                                               // f = f(L_circ^2)  
      omegaq,                                               // f = f(Omega^2)   
    };
    double R_c      (const double&, const arg=epcirc) const;// f = r = R_circ   
    double density  (const double&, const arg=radial) const;// f = rho          
    double Mofr     (const double&, const arg=radial) const;// f = M(<r)        
    double v_c_sq   (const double&, const arg=radial) const;// f = v_circ^2     
    double L_c_sq   (const double&, const arg=radial) const;// f = L_circ^2     
    double Omega_sq (const double&, const arg=radial) const;// f = Omega^2      
    double kappa_sq (const double&, const arg=radial) const;// f = kappa^2      
    double gamma_sq (const double&, const arg=radial) const;// f = gamma^2      
    double Psi      (const double&, const arg=radial) const;// f = Psi==-Phi    
    double Eps_c    (const double&, const arg=radial) const;// f = Eps_circ     
    double mean_vr  (const double&, const arg=radial) const;// f = <v_r>        
    double mean_vth (const double&, const arg=radial) const;// f = <v_theta>    
    double mean_vphi(const double&, const arg=radial) const;// f = <v_phi>      
    double sigma_r  (const double&, const arg=radial) const;// f = sigma_r      
    double sigma_th (const double&, const arg=radial) const;// f = sigma_theta  
    double sigma_phi(const double&, const arg=radial) const;// f = sigma_phi    
    vect_d angmom   (const double&, const arg=radial) const;// f = L            
  private:
    const double* table(const arg) const;
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::spherical_profile                                            
  //                                                                            
  // given a set of bodies, we bin the bodies in radial shells of width Dlr in  
  // ln(r) but containing at at least Nmin bodies. For each bin, we estimate    
  // - the radial cumulative mass profile and gravitational potential           
  // - the density                                                              
  // - the mean radial, azimuthal and meridional velocity and dispersions       
  // - the mean angular momentum                                                
  // We then provide:                                                           
  // - direct access to the tabled values                                       
  //                                                                            
  // NOTE in fact, we make the bins smaller by a factor of two and then consider
  //      the merger of two adjacent smaller bins. This guarantees a somewhat   
  //      higher resolution, but implies that adjacent larger bins share about  
  //      half of their bodies.                                                 
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class spherical_profile {
  private:
    typedef bodies::index index;                   // type used to address body 
    typedef double *pdouble;
    typedef vect_d *pvect_d;
    const int    kmin;
    const double dmax;
    int          n;
    double       mt;
    pdouble      rr,mr,rh,ps,vr,vt,vp,sr,st,sp,skr,skt,skp,kur,kut,kup;
    pvect_d      am;
    //--------------------------------------------------------------------------
    // construction                                                             
    //--------------------------------------------------------------------------
  public:
    spherical_profile(const bodies*,               // I: bodies                 
		      int         =20,             // I: mininum #bodies/bin    
		      double      =0.1,            // I: maximum ln(bin size)   
		      bool        =false)          // I: do velocities?         
      falcON_THROWING;
    ~spherical_profile() {
      falcON_DEL_A(rr);
      falcON_DEL_A(mr);
      falcON_DEL_A(rh);
      if(vr) {
	falcON_DEL_A(vr);
	falcON_DEL_A(vt);
	falcON_DEL_A(vp);
	falcON_DEL_A(sr);
	falcON_DEL_A(st);
	falcON_DEL_A(sp);
	falcON_DEL_A(am);
      }
    }
    //--------------------------------------------------------------------------
    // data access                                                              
    //--------------------------------------------------------------------------
    const int   &N   ()      const { return n; }
    bool         has_vels()  const { return vr!=0; }
    const double&rad (int i) const { return rr[i]; }        // radius           
    const double&Mr  (int i) const { return mr[i]; }        // M(<r)            
          double vcq (int i) const { return mr[i]/rr[i]; }  // v_circ^2(r)      
    const double&rho (int i) const { return rh[i]; }        // rho(r)           
    const double&psi (int i) const { return ps[i]; }        // Psi(r)           
    const double&Mtot()      const { return mt; }           // M_total          
    const double&vrad(int i) const { return vr[i]; }        // mean v_r         
    const double&vthe(int i) const { return vt[i]; }        // mean v_theta     
    const double&vphi(int i) const { return vp[i]; }        // mean v_phi       
    const double&sigr(int i) const { return sr[i]; }        // sigma_r          
    const double&sigt(int i) const { return st[i]; }        // sigma_theta      
    const double&sigp(int i) const { return sp[i]; }        // sigma_phi        
    const double&sker(int i) const { return skr[i]; }       // skew_vr          
    const double&sket(int i) const { return skt[i]; }       // skew_vtheta      
    const double&skep(int i) const { return skp[i]; }       // skew_vphi        
    const double&kurr(int i) const { return kur[i]; }       // kurtosis_r       
    const double&kurt(int i) const { return kut[i]; }       // kurtosis_theta   
    const double&kurp(int i) const { return kup[i]; }       // kurtosis_phi     
    const vect_d&angm(int i) const { return am[i]; }        // mean L           
          double beta(int i) const {                        // Binney's beta    
	    return 1.-twice(square(sr[i])/(square(st[i])+square(sp[i]))); }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::DistributionProfile                                          
  //                                                                            
  // given a table of F (phase-space densities) and masses (optional), we bin   
  // the bodies in F with maximum width dlF in ln(F) but containing at least    
  // Nmin bodies. For each bin, we estimate                                     
  // - the cumulative mass M(F>f)                                               
  // - the volume distribution function v(f)                                    
  // - the cumulative volume V(F>f)                                             
  // We then provide                                                            
  // - access to the tabulated values                                           
  //                                                                            
  // NOTE in fact, we make the bins smaller by a factor of two and then consider
  //      the merger of two adjacent smaller bins. This guarantees a somewhat   
  //      higher resolution, but implies that adjacent larger bins share about  
  //      half of their bodies.                                                 
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class  DistributionProfile {
  private:
    typedef double *pdouble;
    const int    kmin;
    const double dmax;
    double       Mt;
    int          n;
    pdouble      fm,vf,Vf,Mf;
    //--------------------------------------------------------------------------
    // construction                                                             
    //--------------------------------------------------------------------------
  public:
    DistributionProfile(int          ,             // I: size of tables         
			const double*,             // I: phase-space densities  
			const double*,             // I: optional: masses       
			int          =50,          // I: mininum #bodies/bin    
			double       =0.1);        // I: maximum ln(bin size)   
    //--------------------------------------------------------------------------
    ~DistributionProfile() {
      delete[] fm;
      delete[] vf;
      delete[] Vf;
      delete[] Mf;
    }
    //--------------------------------------------------------------------------
    int    const& bins()   const { return n; }
    double const& Mtot()   const { return Mt; }
    double const& f(int i) const { return fm[i]; }
    double const& v(int i) const { return vf[i]; }
    double const& V(int i) const { return Vf[i]; }
    double const& M(int i) const { return Mf[i]; }
    double        D(int i) const { return Mf[i] - fm[i]*Vf[i]; }
  };
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_profile_h
