// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// pext.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_pext_h
#define falcON_included_pext_h

#ifndef falcON_included_auxx_h
#  include <public/auxx.h>
#endif
#ifndef falcON_included_flag_h
#  include <public/flag.h>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::extpot                                                       //
  //                                                                          //
  // is an abstract base class for a potential that can be used as external   //
  // potential in a N-body code, or as stand-alone potential in orbit         //
  // integration                                                              //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class extpot {
  protected:
    typedef double tensor[Ndim][Ndim];
    //--------------------------------------------------------------------------
  public:
    virtual vect force      (                      // R: acceleration           
			     vect const&,          // I: position               
			     real const&) const=0; // I: time                   
    //--------------------------------------------------------------------------
    virtual real pot_f      (                      // R: potential              
			     vect      &,          // O: acceleration           
			     vect const&,          // I: position               
			     real const&) const=0; // I: time                   
    //--------------------------------------------------------------------------
    virtual real rho        (                      // R: density                
			     vect const&,          // I: position               
			     real const&) const=0; // I: time                   
    //--------------------------------------------------------------------------
    virtual void energies   (                      // contribution to totals in:
			     double&u,             // scalar pot energy         
			     tensor&k,             // 2*kin energy tensor       
			     tensor&w,             //   pot energy tensor       
			     vect_d&x,             // dipole                    
			     vect_d&v,             // momentum                  
			     amom_d&l,             // angular momentum          
			     double&m) const  {}   // mass                      
    //--------------------------------------------------------------------------
    virtual bool is_empty() const { return false; }// no potential?             
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::gravity                                                      //
  //                                                                          //
  // is an abstract base class for a class that computes gravity              //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class gravity {
  public:
    //--------------------------------------------------------------------------
    virtual bool is_empty() const { return false; }// no potential?             
    //--------------------------------------------------------------------------
    // set potential & acceleration                                             
    //--------------------------------------------------------------------------
    virtual void set(real       const&,            // I: time                   
		     uint       const&,            // I: # bodies               
		     const real*const&,            // I: body masses (can be 0) 
		     const vect*const&,            // I: body positions         
		     const flag*const&,            // I: body flags             
		     real      *const&,            // O: potential to set       
		     vect      *const&,            // O: accelerations to set   
		     bool       const& =0) const=0;//[I: all or active only?]   
    //--------------------------------------------------------------------------
    // add potential & acceleration                                             
    //--------------------------------------------------------------------------
    virtual void add(real       const&,            // I: time                   
		     uint       const&,            // I: # bodies               
		     const real*const&,            // I: body masses (can be 0) 
		     const vect*const&,            // I: body positions         
		     const flag*const&,            // I: body flags             
		     real      *const&,            // O: potential to add to    
		     vect      *const&,            // O: accelerations to add to
		     bool       const& =0) const=0;//[I: all or active only?]   
    //--------------------------------------------------------------------------
    // set potential & add acceleration                                         
    //--------------------------------------------------------------------------
    virtual void sad(real       const&,            // I: time                   
		     uint       const&,            // I: # bodies               
		     const real*const&,            // I: body masses (can be 0) 
		     const vect*const&,            // I: body positions         
		     const flag*const&,            // I: body flags             
		     real      *const&,            // O: potential to set       
		     vect      *const&,            // O: accelerations to add to
		     bool       const& =0) const=0;//[I: all or active only?]   
    //--------------------------------------------------------------------------
    // set mass density                                                         
    //--------------------------------------------------------------------------
    virtual 
    void set_rho    (real       const&,            // I: time                   
		     uint       const&,            // I: # bodies               
		     const real*const&,            // I: body masses (can be 0) 
		     const vect*const&,            // I: body positions         
		     const flag*const&,            // I: body flags             
		     real      *const&,            // O: density to set         
		     bool       const& =0) const=0;//[I: all or active only?]   
    //--------------------------------------------------------------------------
    // add mass density                                                         
    //--------------------------------------------------------------------------
    virtual 
    void add_rho    (real       const&,            // I: time                   
		     uint       const&,            // I: # bodies               
		     const real*const&,            // I: body masses (can be 0) 
		     const vect*const&,            // I: body positions         
		     const flag*const&,            // I: body flags             
		     real      *const&,            // O: density to add to      
		     bool       const& =0) const=0;//[I: all or active only?]   
  };
}                                                  // END namespace nbdy        
////////////////////////////////////////////////////////////////////////////////
#ifdef falcON_NEMO
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// forward declaration of some nemo stuff needed in construction of nemo_pot  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
  typedef void (*pp_d)(const int*,const double*,double*,double*,const double*);
  typedef void (*pp_f)(const int*,const float *,float *,float *,const float *);
  pp_d get_potential_double(const char*, const char*, const char*);
  pp_f get_potential_float (const char*, const char*, const char*);
  void warning(char*, ...);
}
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::nemo_pot                                                     //
  //                                                                          //
  // is a specialization of class potential using the dynamically linkable    //
  // NEMO potential. That is we wrap the NEMO-routine potential_float or      //
  // potential_double into a nbdy::extpot.                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class nemo_pot : public extpot {
  private:
    //--------------------------------------------------------------------------
    // type                                                                     
    //--------------------------------------------------------------------------
    typedef void (*potproc_real)(                  // type of NEMO routine:     
				 const int*,       // I: # dimensions           
				 const real*,      // I: position               
				       real*,      // O: acceleration           
				       real*,      // O: potential              
				 const real*);     // I: time                   
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
    const potproc_real MYPOT;                      // pter to NEMO routine      
    //--------------------------------------------------------------------------
  public:
    //--------------------------------------------------------------------------
    // construction                                                             
    //                                                                          
    // A sensible usage of the constructor in a nemo_main which has the options 
    // potname, potpars, and potfile would be:                                  
    // ...                                                                      
    // nemo_pot *pot = hasvalue("potname")?                                     
    //   new nemo_pot(getparam("potname"),                                      
    //                hasvalue("potpars")? getparam("potpars") : 0,             
    //                hasvalue("potfile")? getparam("potfile") : 0)             
    //   : 0;                                                                   
    // ...                                                                      
    //--------------------------------------------------------------------------
    nemo_pot(const char*potname,
	     const char*potpars=0,
	     const char*potfile=0) :
#ifdef falcON_REAL_IS_FLOAT
      MYPOT ( get_potential_float (potname,potpars,potfile) )
#else
      MYPOT ( get_potential_double(potname,potpars,potfile) )
#endif
    {
      if(MYPOT==0) ::warning("nemo_pot(): no external potential");
    }
    //--------------------------------------------------------------------------
    // potential and acceleration                                               
    //--------------------------------------------------------------------------
    vect force(vect const&X, real const&T) const {
      register vect F;
      register real P;
      MYPOT(&Ndim, (const real*)X, (real*)F, &P, &T);
      return F;
    }
    //--------------------------------------------------------------------------
    real pot_f(vect&F, const vect&X, real const&T) const {
      register real P;
      MYPOT(&Ndim, (const real*)X, (real*)F, &P, &T);
      return P;
    }
    //--------------------------------------------------------------------------
    bool is_empty() const { return MYPOT == 0; }
    //--------------------------------------------------------------------------
    real rho  (vect const&X, real const&T) const { return zero; }
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::nemo_grav                                                    //
  //                                                                          //
  // is a simple implementation of gravity using class nbdy::nemo_pot         //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class nemo_grav : 
    private nemo_pot,
    public  gravity
  {
  public:
    //--------------------------------------------------------------------------
    nemo_grav(const char*potname,                  // I: potname                
	      const char*potpars =0,               //[I: potpars]               
	      const char*potfile =0) :             //[I: potfile]               
      nemo_pot(potname,potpars,potfile) {}
    //--------------------------------------------------------------------------
    bool is_empty() { return nemo_pot::is_empty(); }
    //--------------------------------------------------------------------------
#define LoopAll    for(int n=0; n!=N; ++n)
#define LoopActive for(int n=0; n!=N; ++n) if(is_active(F[n]))
    //--------------------------------------------------------------------------
    void set(real       const&T,                   // I: time                   
	     uint       const&N,                   // I: # bodies               
	     const real*const&M,                   // I: body masses (can be 0) 
	     const vect*const&X,                   // I: body positions         
	     const flag*const&F,                   // I: body flags  (can be 0) 
	     real      *const&P,                   // O: potential to set       
	     vect      *const&A,                   // O: accelerations to set   
	     bool       const&all=0) const         //[I: all or active only?]   
    {
      if(all || F==0) {                            // IF all are active         
	LoopAll                                    //   LOOP all bodies         
	  P[n] = pot_f(A[n],X[n],T);               //     set pot & acc         
      } else {                                     // ELSE (only active)        
	LoopActive                                 //   LOOP active bodies      
	  P[n] = pot_f(A[n],X[n],T);               //     set pot & acc         
      }                                            // ENDIF                     
    }
    //--------------------------------------------------------------------------
    void add(real       const&T,                   // I: time                   
	     uint       const&N,                   // I: # bodies               
	     const real*const&M,                   // I: body masses (can be 0) 
	     const vect*const&X,                   // I: body positions         
	     const flag*const&F,                   // I: body flags  (can be 0) 
	     real      *const&P,                   // O: potential to add to    
	     vect      *const&A,                   // O: accelerations to add to
	     bool       const&all=0) const         //[I: all or active only?]   
    {
      register vect a;                             // to hold acceleration      
      if(all || F==0) {                            // IF all are active         
	LoopAll {                                  //   LOOP all bodies         
	  P[n]+= pot_f(a,X[n],T);                  //     add potential         
	  A[n]+= a;                                //     add acceleration      
	}                                          //   END LOOP                
      } else {                                     // ELSE (only active)        
	LoopActive {                               //   LOOP active bodies      
	  P[n]+= pot_f(a,X[n],T);                  //     add potential         
	  A[n]+= a;                                //     add acceleration      
	}                                          //   END LOOP                
      }                                            // ENDIF                     
    }
    //--------------------------------------------------------------------------
    void sad(real       const&T,                   // I: time                   
	     uint       const&N,                   // I: # bodies               
	     const real*const&M,                   // I: body masses (can be 0) 
	     const vect*const&X,                   // I: body positions         
	     const flag*const&F,                   // I: body flags  (can be 0) 
	     real      *const&P,                   // O: potential to set       
	     vect      *const&A,                   // O: accelerations to add to
	     bool       const&all=0) const         //[I: all or active only?]   
    {
      register vect a;                             // to hold acceleration      
      if(all || F==0) {                            // IF all are active         
	LoopAll {                                  //   LOOP all bodies         
	  P[n] = pot_f(a,X[n],T);                  //     set potential         
	  A[n]+= a;                                //     add acceleration      
	}                                          //   END LOOP                
      } else {                                     // ELSE (only active)        
	LoopActive {                               //   LOOP active bodies      
	  P[n] = pot_f(a,X[n],T);                  //     set potential         
	  A[n]+= a;                                //     add acceleration      
	}                                          //   END LOOP                
      }                                            // ENDIF                     
    }
    //--------------------------------------------------------------------------
    void set_rho(real       const&T,               // I: time                   
		 uint       const&N,               // I: # bodies               
		 const real*const&M,               // I: body masses (can be 0) 
		 const vect*const&X,               // I: body positions         
		 const flag*const&F,               // I: body flags  (can be 0) 
		 real      *const&R,               // O: density to set         
		 bool       const&all=0) const     //[I: all or active only?]   
    {
      if(all || F==0) {                            // IF all are active         
	LoopAll                                    //   LOOP all bodies         
	  R[n] = rho(X[n],T);                      //     set density           
      } else {                                     // ELSE (only active)        
	LoopActive                                 //   LOOP active bodies      
	  R[n] = rho(X[n],T);                      //     set density           
      }                                            // ENDIF                     
    }
    //--------------------------------------------------------------------------
    void add_rho(real       const&T,               // I: time                   
		 uint       const&N,               // I: # bodies               
		 const real*const&M,               // I: body masses (can be 0) 
		 const vect*const&X,               // I: body positions         
		 const flag*const&F,               // I: body flags  (can be 0) 
		 real      *const&R,               // O: density to add to      
		 bool       const&all=0) const     //[I: all or active only?]   
    {
      if(all || F==0) {                            // IF all are active         
	LoopAll                                    //   LOOP all bodies         
	  R[n]+= rho(X[n],T);                      //     set density           
      } else {                                     // ELSE (only active)        
	LoopActive                                 //   LOOP active bodies      
	  R[n]+= rho(X[n],T);                      //     set density           
      }                                            // ENDIF                     
    }
#undef LoopAll
#undef LoopActive
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
}                                                  // END namespace nbdy        
////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_NEMO               
////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_included_pext_h    
