// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// kern.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2002                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// versions/changes                                                            |
//                                                                             |
// 09/11/2002 - created by drawing from grav.h                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// defines                                                                     |
//                                                                             |
// class grav_kern                                                             |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef included_kern_h
#define included_kern_h

#ifndef included_grat_h
#  include <public/grat.h>
#endif

////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class nbdy::grav_kern                                                      
  //                                                                            
  // This class implements the direct summation and approximate computation     
  // of gravity between tree nodes.                                             
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class grav_kern {
    //--------------------------------------------------------------------------
    // static data                                                              
    //--------------------------------------------------------------------------
  protected:
    static const int
      N_C1    = 1,
      N_C2    = N_C1 + ten1::NDAT,
      N_C3    = N_C2 + ten2::NDAT,
#if   P_ORDER > 3
      N_C4    = N_C3 + ten3::NDAT,
# if  P_ORDER > 4
      N_C5    = N_C4 + ten4::NDAT,
#  if P_ORDER > 5
      N_C6    = N_C5 + ten5::NDAT,
      N_COEFF = N_C6 + ten6::NDAT;
#  else
      N_COEFF = N_C5 + ten5::NDAT;
#  endif
# else
      N_COEFF = N_C4 + ten4::NDAT;
# endif
#else
      N_COEFF = N_C3 + ten3::NDAT;
#endif
    //--------------------------------------------------------------------------
    // types                                                                    
    //--------------------------------------------------------------------------
  private:
    typedef grav_tree::cell_iterator  cell_iter;   // cell iterator             
    typedef grav_tree::soul_iterator  soul_iter;   // soul iterator             
    //--------------------------------------------------------------------------
    // protected data                                                           
    //--------------------------------------------------------------------------
  protected:
    kern_type             KERN;                    // softening kernel          
#ifdef ALLOW_INDI
    soft_type             SOFT;                    // global/individual         
#endif
    bool                  FULL_POT;                // Pth pole in pot           
    mutable real          EPS, EPQ, HEQ, QEQ;      // eps & powers of eps       
    //--------------------------------------------------------------------------
    // private methods                                                          
    //--------------------------------------------------------------------------
  private:
#define ARGS__ soul_iter const&, uint const&, soul_iter const&, uint const&
    void many_AA(ARGS__) const;
    void many_AS(ARGS__) const;
    void many_AN(ARGS__) const;
    void many_SA(ARGS__) const;
    void many_SS(ARGS__) const;
    void many_SN(ARGS__) const;
    void many_NA(ARGS__) const;
    void many_NS(ARGS__) const;
#undef ARGS__
    //--------------------------------------------------------------------------
    // protected methods                                                        
    //--------------------------------------------------------------------------
  protected:
    real&      c0(real*         const&T) const { return *(T); }
    real&      c0(cell_iter     const&C) const { return *(C->coeffs()); }
    //--------------------------------------------------------------------------
    ten1       c1(real*         const&T) const { return ten1(T+N_C1); }
    ten1       c1(cell_iter     const&C) const { return c1(C->coeffs()); }
    //--------------------------------------------------------------------------
    ten2       c2(real*         const&T) const { return ten2(T+N_C2); }
    ten2       c2(cell_iter     const&C) const { return c2(C->coeffs()); }
    //--------------------------------------------------------------------------
    ten3       c3(real*         const&T) const { return ten3(T+N_C3); }
    ten3       c3(cell_iter     const&C) const { return c3(C->coeffs()); }
    //--------------------------------------------------------------------------
#if P_ORDER > 3
    ten4       c4(real*         const&T) const { return ten4(T+N_C4); }
    ten4       c4(cell_iter     const&C) const { return c4(C->coeffs()); }
    //--------------------------------------------------------------------------
#if P_ORDER > 4
    ten5       c5(real*         const&T) const { return ten5(T+N_C5); }
    ten5       c5(cell_iter     const&C) const { return c5(C->coeffs()); }
    //--------------------------------------------------------------------------
#if P_ORDER > 5
    ten6       c6(real*         const&T) const { return ten6(T+N_C6); }
    ten6       c6(cell_iter     const&C) const { return c6(C->coeffs()); }
#endif
#endif
#endif
    //--------------------------------------------------------------------------
    // main purpose methods                                                     
    //--------------------------------------------------------------------------
    grav_kern(
	      kern_type       k,                   // I: type of kernel         
	      real            e,                   //[I: softening length]      
#ifdef ALLOW_INDI
	      soft_type       s,                   //[I: type of softening]     
#endif
	      bool            fp    = false) :     //[I: use Pth pole in pot]   
      KERN     ( k ),                              // set softening kernel      
#ifdef ALLOW_INDI
      SOFT     ( s ),                              // set softening type        
#endif
      FULL_POT ( fp ),                             // Pth pole in pot?          
      EPS      ( e ),                              // set softening length      
      EPQ      ( EPS*EPS ),                        // set eps^2                 
      HEQ      ( half*EPQ ),                       // set 0.50 * eps^2          
      QEQ      ( half*HEQ ) {}                     // set 0.25 * eps^2          
    //--------------------------------------------------------------------------
    // single soul-soul interaction                                             
    //--------------------------------------------------------------------------
    void single(soul_iter const&, soul_iter const&) const;
    //--------------------------------------------------------------------------
    // cell-soul interaction via direct summation                               
    //--------------------------------------------------------------------------
    void many(cell_iter const&, soul_iter const&) const;
    //--------------------------------------------------------------------------
    // cell-cell interaction via direct summation                               
    //--------------------------------------------------------------------------
    void many(cell_iter const&, cell_iter const&) const;
    //--------------------------------------------------------------------------
    // cell-self interaction via direct summation                               
    //--------------------------------------------------------------------------
    void many(cell_iter const&) const;
    //--------------------------------------------------------------------------
    // cell-soul interaction via approximation                                  
    //--------------------------------------------------------------------------
    void grav(cell_iter const&, soul_iter const&, vect&, real const&) const;
    //--------------------------------------------------------------------------
    // cell-cell interaction via approximation                                  
    //--------------------------------------------------------------------------
    void grav(cell_iter const&, cell_iter const&, vect&, real const&) const;
    //--------------------------------------------------------------------------
    void flush_buffers() const {}
    //--------------------------------------------------------------------------
  public:
    const real&current_eps  ()     const { return EPS; }
    const real&current_epsq ()     const { return EPQ; }
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  inline void grav_kern::many(cell_iter const&CA, cell_iter const&CB) const {

    const uint      NA=number(CA),       NB=number(CB);
    const soul_iter A0=CA.begin_souls(), B0=CB.begin_souls();
    if(NA%4 > NB%4) {
      if       (al_sink(CA))
	if     (al_sink(CB)) many_AA(A0,NA,B0,NB); // sinks: all  A, all  B    
	else if(is_sink(CB)) many_AS(A0,NA,B0,NB); // sinks: all  A, some B    
	else                 many_AN(A0,NA,B0,NB); // sinks: all  A, no   B    
      else if  (is_sink(CA))
	if     (al_sink(CB)) many_SA(A0,NA,B0,NB); // sinks: some A, all  B    
	else if(is_sink(CB)) many_SS(A0,NA,B0,NB); // sinks: some A, some B    
	else                 many_SN(A0,NA,B0,NB); // sinks: some A, no   B    
      else
	if     (al_sink(CB)) many_NA(A0,NA,B0,NB); // sinks: no   A, all  B    
	else if(is_sink(CB)) many_NS(A0,NA,B0,NB); // sinks: no   A, some B    
    } else {
      if       (al_sink(CB))
	if     (al_sink(CA)) many_AA(B0,NB,A0,NA); // sinks: all  B, all  A    
	else if(is_sink(CA)) many_AS(B0,NB,A0,NA); // sinks: all  B, some A    
	else                 many_AN(B0,NB,A0,NA); // sinks: all  B, no   A    
      else if(is_sink(CB))
	if     (al_sink(CA)) many_SA(B0,NB,A0,NA); // sinks: some B, all  A    
	else if(is_sink(CA)) many_SS(B0,NB,A0,NA); // sinks: some B, some A    
	else                 many_SN(B0,NB,A0,NA); // sinks: some B, no   A    
      else
	if     (al_sink(CA)) many_NA(B0,NB,A0,NA); // sinks: no   B, all  A    
	else if(is_sink(CA)) many_NS(B0,NB,A0,NA); // sinks: no   B, some A    
    }
  }
  //////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////
#endif                                             // included_kern_h           
