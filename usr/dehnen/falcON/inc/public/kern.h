// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// kern.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2003                                          |
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
// class grav_kern_all                                                         |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_kern_h
#define falcON_included_kern_h

#ifndef falcON_included_grat_h
#  include <public/grat.h>
#endif

////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::grav_kern_base                                               //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class grav_kern_base {
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
  public:
    kern_type             KERN;                    // softening kernel          
#ifdef falcON_INDI
    soft_type             SOFT;                    // global/individual         
#endif
    bool                  FULL_POT;                // Pth pole in pot           
    mutable real          EPS, EPQ, HEQ, QEQ;      // eps & powers of eps       
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
    real&      c0(real*           const&T) const { return *(T); }
    real&      c0(grav::cell_iter const&C) const { return *(C->coeffs()); }
    //--------------------------------------------------------------------------
    ten1       c1(real*           const&T) const { return ten1(T+grav::N_C1); }
    ten1       c1(grav::cell_iter const&C) const { return c1(C->coeffs()); }
    //--------------------------------------------------------------------------
    ten2       c2(real*           const&T) const { return ten2(T+grav::N_C2); }
    ten2       c2(grav::cell_iter const&C) const { return c2(C->coeffs()); }
    //--------------------------------------------------------------------------
    ten3       c3(real*           const&T) const { return ten3(T+grav::N_C3); }
    ten3       c3(grav::cell_iter const&C) const { return c3(C->coeffs()); }
    //--------------------------------------------------------------------------
#if falcON_ORDER > 3
    ten4       c4(real*           const&T) const { return ten4(T+grav::N_C4); }
    ten4       c4(grav::cell_iter const&C) const { return c4(C->coeffs()); }
    //--------------------------------------------------------------------------
#if falcON_ORDER > 4
    ten5       c5(real*           const&T) const { return ten5(T+grav::N_C5); }
    ten5       c5(grav::cell_iter const&C) const { return c5(C->coeffs()); }
    //--------------------------------------------------------------------------
#if falcON_ORDER > 5
    ten6       c6(real*           const&T) const { return ten6(T+grav::N_C6); }
    ten6       c6(grav::cell_iter const&C) const { return c6(C->coeffs()); }
#endif
#endif
#endif
    //--------------------------------------------------------------------------
    // protected methods                                                        
    //--------------------------------------------------------------------------
    grav_kern_base(
		   kern_type const&k,              // I: type of kernel         
		   real      const&e,              // I: softening length       
#ifdef falcON_INDI
		   soft_type const&s,              // I: type of softening      
#endif
		   bool      const&fp) :           // I: use Pth pole in pot    
      KERN     ( k ),                              // set softening kernel      
#ifdef falcON_INDI
      SOFT     ( s ),                              // set softening type        
#endif
      FULL_POT ( fp ),                             // Pth pole in pot?          
      EPS      ( e ),                              // set softening length      
      EPQ      ( EPS*EPS ),                        // set eps^2                 
      HEQ      ( half*EPQ ),                       // set 0.50 * eps^2          
      QEQ      ( half*HEQ ) {}                     // set 0.25 * eps^2          
    //--------------------------------------------------------------------------
  protected:
#define ARGS__ grav::soul_iter const&, uint const&, 	\
               grav::soul_iter const&, uint const&
    void many_AA(ARGS__) const;
    void many_AS(ARGS__) const;
    void many_AN(ARGS__) const;
    void many_SA(ARGS__) const;
    void many_SS(ARGS__) const;
    void many_SN(ARGS__) const;
    void many_NA(ARGS__) const;
    void many_NS(ARGS__) const;
#undef ARGS__
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::grav_kern                                                    //
  //                                                                          //
  // This class implements the direct summation and approximate computation   //
  // of gravity between tree nodes.                                           //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class grav_kern : protected grav_kern_base {
    //--------------------------------------------------------------------------
    // main purpose methods                                                     
    //--------------------------------------------------------------------------
  protected:
    grav_kern(
	      kern_type const&k,                   // I: type of kernel         
	      real      const&e,                   // I: softening length       
#ifdef falcON_INDI
	      soft_type const&s,                   // I: type of softening      
#endif
	      bool      const&fp=false) :          // I:[use Pth pole in pot]   
      grav_kern_base(k,e,
#ifdef falcON_INDI
		     s,
#endif
		     fp) {}
    //--------------------------------------------------------------------------
    // single soul-soul interaction                                             
    //--------------------------------------------------------------------------
    void single(grav::soul_iter const&, grav::soul_iter const&) const;
    //--------------------------------------------------------------------------
    // cell-soul interaction via direct summation                               
    //--------------------------------------------------------------------------
    void many(grav::cell_iter const&, grav::soul_iter const&) const;
    //--------------------------------------------------------------------------
    // cell-cell interaction via direct summation                               
    //--------------------------------------------------------------------------
    void many(grav::cell_iter const&, grav::cell_iter const&) const;
    //--------------------------------------------------------------------------
    // cell-self interaction via direct summation                               
    //--------------------------------------------------------------------------
    void many(grav::cell_iter const&) const;
    //--------------------------------------------------------------------------
    // cell-soul interaction via approximation                                  
    //--------------------------------------------------------------------------
    void grav(grav::cell_iter const&, grav::soul_iter const&,
	      vect&, real const&) const;
    //--------------------------------------------------------------------------
    // cell-cell interaction via approximation                                  
    //--------------------------------------------------------------------------
    void grav(grav::cell_iter const&, grav::cell_iter const&,
	      vect&, real const&) const;
    //--------------------------------------------------------------------------
    void flush_buffers() const {}
    //--------------------------------------------------------------------------
  public:
    const real&current_eps  ()     const { return EPS; }
    const real&current_epsq ()     const { return EPQ; }
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  inline
  void grav_kern::many(grav::cell_iter const&CA, grav::cell_iter const&CB) const
  {
    const uint            NA=number(CA),       NB=number(CB);
    const grav::soul_iter A0=CA.begin_souls(), B0=CB.begin_souls();
    if(NA%4 > NB%4) {
      if       (al_active(CA))
	if     (al_active(CB)) many_AA(A0,NA,B0,NB); // actives: all  A, all  B 
	else if(is_active(CB)) many_AS(A0,NA,B0,NB); // actives: all  A, some B 
	else                   many_AN(A0,NA,B0,NB); // actives: all  A, no   B 
      else if  (is_active(CA))
	if     (al_active(CB)) many_SA(A0,NA,B0,NB); // actives: some A, all  B 
	else if(is_active(CB)) many_SS(A0,NA,B0,NB); // actives: some A, some B 
	else                   many_SN(A0,NA,B0,NB); // actives: some A, no   B 
      else
	if     (al_active(CB)) many_NA(A0,NA,B0,NB); // actives: no   A, all  B 
	else if(is_active(CB)) many_NS(A0,NA,B0,NB); // actives: no   A, some B 
    } else {
      if       (al_active(CB))
	if     (al_active(CA)) many_AA(B0,NB,A0,NA); // actives: all  B, all  A 
	else if(is_active(CA)) many_AS(B0,NB,A0,NA); // actives: all  B, some A 
	else                   many_AN(B0,NB,A0,NA); // actives: all  B, no   A 
      else if(is_active(CB))
	if     (al_active(CA)) many_SA(B0,NB,A0,NA); // actives: some B, all  A 
	else if(is_active(CA)) many_SS(B0,NB,A0,NA); // actives: some B, some A 
	else                   many_SN(B0,NB,A0,NA); // actives: some B, no   A 
      else
	if     (al_active(CA)) many_NA(B0,NB,A0,NA); // actives: no   B, all  A 
	else if(is_active(CA)) many_NS(B0,NB,A0,NA); // actives: no   B, some A 
    }
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::grav_kern_all                                                //
  //                                                                          //
  // Like grav_kern, except that all cells and souls are assumed active.      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class grav_kern_all : protected grav_kern_base {
    //--------------------------------------------------------------------------
    // main purpose methods                                                     
    //--------------------------------------------------------------------------
  protected:
    grav_kern_all(
		  kern_type const&k,               // I: type of kernel         
		  real      const&e,               // I: softening length       
#ifdef falcON_INDI
		  soft_type const&s,               // I: type of softening      
#endif
		  bool      const&fp=false) :      // I:[use Pth pole in pot]   
      grav_kern_base(k,e,
#ifdef falcON_INDI
		     s,
#endif
		     fp) {}
    //--------------------------------------------------------------------------
    // single soul-soul interaction                                             
    //--------------------------------------------------------------------------
    void single(grav::soul_iter const&, grav::soul_iter const&) const;
    //--------------------------------------------------------------------------
    // cell-soul interaction via direct summation                               
    //--------------------------------------------------------------------------
    void many(grav::cell_iter const&, grav::soul_iter const&) const;
    //--------------------------------------------------------------------------
    // cell-cell interaction via direct summation                               
    //--------------------------------------------------------------------------
    void many(grav::cell_iter const&, grav::cell_iter const&) const;
    //--------------------------------------------------------------------------
    // cell-self interaction via direct summation                               
    //--------------------------------------------------------------------------
    void many(grav::cell_iter const&) const;
    //--------------------------------------------------------------------------
    // cell-soul interaction via approximation                                  
    //--------------------------------------------------------------------------
    void grav(grav::cell_iter const&, grav::soul_iter const&,
	      vect&, real const&) const;
    //--------------------------------------------------------------------------
    // cell-cell interaction via approximation                                  
    //--------------------------------------------------------------------------
    void grav(grav::cell_iter const&, grav::cell_iter const&,
	      vect&, real const&) const;
    //--------------------------------------------------------------------------
    void flush_buffers() const {}
    //--------------------------------------------------------------------------
  public:
    const real&current_eps  ()     const { return EPS; }
    const real&current_epsq ()     const { return EPQ; }
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  inline  void grav_kern_all::many(grav::cell_iter const&CA,
				   grav::cell_iter const&CB) const
  {
    const uint NA=number(CA), NB=number(CB);
    if(NA%4 > NB%4) many_AA(CA.begin_souls(),NA,CB.begin_souls(),NB);
    else            many_AA(CB.begin_souls(),NB,CA.begin_souls(),NA);
  }
  //////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_included_kern_h    
