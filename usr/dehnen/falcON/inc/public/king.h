// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// king.h                                                                      |
//                                                                             |
// Copyright (C) 2000-2003, 2005, 2007  Walter Dehnen                          |
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
#ifndef falcON_included_king_h
#define falcON_included_king_h

#include <public/basic.h>

////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::king_model                                                 //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class king_model {
    //--------------------------------------------------------------------------
    // private data                                                             
    //--------------------------------------------------------------------------
  private:
    double   Psi0,P0,r0,                          // Psi(0), Phi_0, r_0,        
             c,rh0;                               // log10(rt/r0), rho(0)       
    unsigned N;                                   // N: # grid points           
    double   *r,*ps,*m,*rh;                       // grids in r, Psi, M(r),rho  
    double   vscal,rscal,Pscal,rhscl,mscal,sdscl, // scaling factors            
             drscl;
    //--------------------------------------------------------------------------
    // private methods                                                          
    //--------------------------------------------------------------------------
    void setup(const unsigned);
    //--------------------------------------------------------------------------
    // protected methods                                                        
    //--------------------------------------------------------------------------
  protected:
    king_model() :
      Psi0(0.0), N(0), r(0), ps(0), m(0), rh(0),
      vscal(1), rscal(1), Pscal(1), rhscl(1), mscal(1), sdscl(1), drscl(1)
      {}
    //--------------------------------------------------------------------------
    void set(const double w, const unsigned n) {
      Psi0 = w;
      setup(n);
    }
    //--------------------------------------------------------------------------
    void reset();
    //--------------------------------------------------------------------------
    // construction & destruction                                               
    //--------------------------------------------------------------------------
  public:
    king_model(const double W,                    // W_0 := Psi_0 / sigma^2     
	       const unsigned n=1000) :           // # grid points              
      Psi0(W), N(0), r(0), ps(0), m(0), rh(0),
      vscal(1), rscal(1), Pscal(1), rhscl(1), mscal(1), sdscl(1), drscl(1)
    {
      if(W <= 0.0) error("king_model: W0 [%f] must be positive\n",W);
      setup(n);
    }
    //--------------------------------------------------------------------------
    ~king_model() { reset(); }
    //--------------------------------------------------------------------------
    // re-scaling                                                               
    //--------------------------------------------------------------------------
    void reset_scales_tidal(                      // set scales (with G=1):     
			    const double,         // I: total mass   is this    
			    const double);        // I: tidal radius is this    
    //--------------------------------------------------------------------------
    void reset_scales_core (                      // set scales (with G=1)      
			    const double,         // I: total mass   is this,   
			    const double);        // I: core radius  is this    
    //--------------------------------------------------------------------------
    // routines returning const global informations                             
    //--------------------------------------------------------------------------
    double tidal_radius    () const { return rscal*r[N-1]; }
    double total_mass      () const { return mscal*m[N-1]; }
    double core_radius     () const { return rscal*r0; }
    double sigma           () const { return vscal; }
    double concentration   () const { return c; }
    double core_density    () const { return rhscl*rh0; }
    double Phi0            () const { return Pscal*(P0-Psi0); }
    double W0              () const { return Psi0; }
    double half_mass_radius() const;              // half-mass radius           
    double rms_radius      () const;              // sqrt<r^2>                  
    double Etot            () const;              // self energy of equilibrium 
    //--------------------------------------------------------------------------
    // direct const access to grids                                             
    //--------------------------------------------------------------------------
    double rad_grid     (const int i) const { return rscal*r[i]; }
    double pot_grid     (const int i) const { return Pscal*(P0-ps[i]); }
    double cum_grid     (const int i) const { return mscal*m[i]; }
    double rho_grid     (const int i) const { return rhscl*rh[i]; }
    double drh_grid     (const int i) const {
      if(i            ==0) return 0.0;
      if(unsigned(i+1)==N) return drscl * rh[N-1]/(r[N-1]-r[N-2]);
      return drscl*(rh[i+1]-rh[i-1])/(r[i+1]-r[i-1]);
    }
    double sur_grid     (const int i) const {
      return sdscl*SurfDens(rscal*r[i]);
    }
    //--------------------------------------------------------------------------
    // miscellaneous                                                            
    //--------------------------------------------------------------------------
    double SurfDens     (double) const;           // surface density profile    
    void   write_table  (const char*) const;      // output r,Phi,rho,M(<r)     
    //--------------------------------------------------------------------------
    // create random radius & velocity                                          
    //--------------------------------------------------------------------------
    double random(                                // R: Phi(r)                  
		  const double, const double,     // I: RNs in (0,1)            
		  double&, double&) const;        // O: r & v                   
  };
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_king_h
