// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// pext.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2002                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef included_pext_h
#define included_pext_h

#ifndef included_auxx_h
#  include <public/auxx.h>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  // class nbdy::potential                                                      
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // this is an abstract base class for a potential that can be used as external
  // potential in a N-body code, or as stand-alone potential in orbit           
  // integration                                                                
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class extpot {
  protected:
    typedef real tensor[NDIM][NDIM];
    //--------------------------------------------------------------------------
  public:
    virtual vect force      (                      // R: acceleration           
			     const vect&,          // I: position               
			     const real) const=0;  // I: time                   
    //--------------------------------------------------------------------------
    virtual real pot_f      (                      // R: potential              
			     vect&,                // O: acceleration           
			     const vect&,          // I: position               
			     const real)const=0;   // I: time                   
    //--------------------------------------------------------------------------
    virtual real rho        (                      // R: density                
			     const vect&,          // I: position               
			     const real) const=0;  // I: time                   
    //--------------------------------------------------------------------------
    virtual void energies   (                      // contribution to totals in:
			     real  &u,             // scalar pot energy         
			     tensor&w,             //   pot energy tensor       
			     tensor&t,             // 2*kin energy tensor       
			     vect  &x,             // dipole                    
			     vect  &v,             // momentum                  
			     amom  &l,             // angular momentum          
			     real  &m) const  {}   // mass                      
    //--------------------------------------------------------------------------
    virtual bool is_empty() const { return false; }// no potential?             
    //--------------------------------------------------------------------------
  };
}                                                  // END namespace nbdy        
#ifdef ALLOW_NEMO
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  // class nbdy::nemo_pot                                                       
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // this is a specialization of class potential using the dynamically linkable 
  // NEMO potential. That is we wrap the NEMO-routine potential_float or        
  // potential_double into a nbdy::potential.                                   
  //                                                                            
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
    const int          ndim;                       // NDIM (need pointer to)    
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
    nemo_pot(char*,                                // I: potname                
	     char* =0,                             //[I: potpars]               
	     char* =0);                            //[I: potfile]               
    //--------------------------------------------------------------------------
    // potential and acceleration                                               
    //--------------------------------------------------------------------------
    vect force(const vect& X, const real T) const {
      register vect F;
      register real P;
      MYPOT(&ndim, X.const_pointer(), F.pointer(), &P, &T);
      return F;
    }
    //--------------------------------------------------------------------------
    real pot_f(vect&F, const vect&X, const real T) const {
      register real P;
      MYPOT(&ndim, X.const_pointer(), F.pointer(), &P, &T);
      return P;
    }
    //--------------------------------------------------------------------------
    bool is_empty() const { return MYPOT == 0; }
    //--------------------------------------------------------------------------
    real rho  (const vect& X, const real T) const { return zero; }
    //--------------------------------------------------------------------------
  };
}                                                  // END namespace nbdy        
//------------------------------------------------------------------------------
// forward declaration of some nemo stuff needed in construction of nemo_pot    
//------------------------------------------------------------------------------
extern "C" {
  typedef void (*pp_d)(const int*,const double*,double*,double*,const double*);
  typedef void (*pp_f)(const int*,const float *,float *,float *,const float *);
  pp_d get_potential_double(char*, char*, char*);
  pp_f get_potential_float (char*, char*, char*);
  void warning(char*, ...);
}
//------------------------------------------------------------------------------
// inline nemo_pot::nemo_pot()                                                  
//------------------------------------------------------------------------------
inline nbdy::nemo_pot::nemo_pot(char* potname, char* potpars, char* potfile)
  : ndim  ( NDIM ),
#if defined(DOUBLE_SINGLE) || defined(DOUBLE_DOUBLE)
    MYPOT ( get_potential_double(potname,potpars,potfile) )
#else
    MYPOT ( get_potential_float (potname,potpars,potfile) )
#endif
{
  if(MYPOT==0) warning("nemo_pot(): no external potential");
}
////////////////////////////////////////////////////////////////////////////////
#endif                                             // ALLOW_NEMO                
////////////////////////////////////////////////////////////////////////////////
#endif                                             // included_pext_h           
