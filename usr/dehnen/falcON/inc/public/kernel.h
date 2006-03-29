// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
/// /file inc/public/kernel.h                                                  |
//                                                                             |
// Copyright (C) 2000-2005  Walter Dehnen                                      |
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
// design note: experiments have shown that, strange enough, a code with       |
//              GravKern and GravKernAll being templates with template         |
//              parameter being the kern_type, is somewhat slower (with gcc    |
//              3.2.2: 4% for SSE code). This affects only the approximate     |
//              part of gravity.                                               |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_kernel_h
#define falcON_included_kernel_h

#ifndef falcON_included_gravity_h
#  include <public/gravity.h>
#endif

#if defined(falcON_SSE_CODE) && !defined(falcON_included_simd_h)
#  include <proper/simd.h>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::TaylorSeries                                               //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class TaylorSeries {
    vect       X;                                  // expansion center          
    grav::Cset C;                                  // expansion coefficients    
  public:
    explicit
    TaylorSeries(vect         const&x) : X(x), C(zero) {} 
    TaylorSeries(TaylorSeries const&T) : X(T.X) { C = T.C; }
    //--------------------------------------------------------------------------
    friend bool is_empty(TaylorSeries const&T) {
      return T.C == zero; }
    //--------------------------------------------------------------------------
    inline void shift_and_add(const grav::cell*const&);
    inline void extract_grav (grav::leaf*const&) const;
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::GravKernBase                                                 
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class GravKernBase {
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
  protected:
    const kern_type       KERN;                    // softening kernel          
#ifdef falcON_INDI
    const bool            INDI_SOFT;               // use individual eps?       
#endif
    mutable real          EPS, EQ;                 // eps, eps^2                
#ifdef  falcON_SSE_CODE
    mutable fvec4         fHQ, fQQ;                // eps^2/2, eps^2/4          
#else
    mutable real          HQ, QQ;                  // eps^2/2, eps^2/4          
#endif
    //--------------------------------------------------------------------------
  private:
    falcON::pool         *COEFF_POOL;              // pool for TaylorCoeffs     
    mutable int           NC, MAXNC;               // # TaylorCoeffs ever used  
    //--------------------------------------------------------------------------
    // blocking of 4 cell-leaf interactions                                     
    //--------------------------------------------------------------------------
#ifdef  falcON_SSE_CODE
  public:
    struct acl_block {
      int              NR;
      grav::cell_pter  A[4];
      grav::leaf_pter  B[4];
      vect             dX[4];
      fvec4            D0,D1;
      void load_g(grav::cell*a, grav::leaf*b, vect&dR, real Rq)
      {
	A [NR] = a;
	B [NR] = b;
	dX[NR] = dR;
	D0[NR] = mass(a)*mass(b);
	D1[NR] = Rq;
	++NR;
      }
#ifdef falcON_INDI
      fvec4            EQ;
      void load_i(grav::cell*a, grav::leaf*b, vect&dR, real Rq)
      {
	A [NR] = a;
	B [NR] = b;
	dX[NR] = dR;
	EQ[NR] = square(eph(a)+eph(b));
	D0[NR] = mass(a)*mass(b);
	D1[NR] = Rq+EQ[NR];
	++NR;
      }
#endif
      acl_block() : NR(0) {}
      bool is_full () const { return NR==4; }
      bool is_empty() const { return NR==0; }
      void reset   ()       { NR=0; }
    };
    mutable acl_block ACL;
    //--------------------------------------------------------------------------
    // blocking of 4 cell-cell interactions                                     
    //--------------------------------------------------------------------------
    struct acc_block {
      int              NR;
      grav::cell_pter  A[4], B[4];
      vect             dX[4];
      fvec4            D0,D1;
      void load_g(grav::cell*a, grav::cell*b, vect&dR, real Rq)
      {
	A [NR] = a;
	B [NR] = b;
	dX[NR] = dR;
	D0[NR] = mass(a)*mass(b);
	D1[NR] = Rq;
	++NR;
      }
#ifdef falcON_INDI
      fvec4            EQ;
      void load_i(grav::cell*a, grav::cell*b, vect&dR, real Rq)
      {
	A [NR] = a;
	B [NR] = b;
	dX[NR] = dR;
	EQ[NR] = square(eph(a)+eph(b));
	D0[NR] = mass(a)*mass(b);
	D1[NR] = Rq+EQ[NR];
	++NR;
      }
#endif
      acc_block() : NR(0) {}
      bool is_full () const { return NR==4; }
      bool is_empty() const { return NR==0; }
      void reset   ()       { NR=0; }
    };
    mutable acc_block ACC;
#endif
    //--------------------------------------------------------------------------
    // protected methods                                                        
    //--------------------------------------------------------------------------
  protected:
    GravKernBase(
		 kern_type const&k,                // I: type of kernel         
		 real      const&e,                // I: softening length       
#ifdef falcON_INDI
		 bool      const&s,                // I: type of softening      
#endif
		 unsigned  const&np) :             // I: initial pool size      
      KERN       ( k ),                            // set softening kernel      
#ifdef falcON_INDI
      INDI_SOFT  ( s ),                            // set softening type        
#endif
      EPS        ( e ),                            // set softening length      
      EQ         ( e*e ),
#ifdef falcON_SSE_CODE
      fHQ        ( half * EQ ),
      fQQ        ( quarter * EQ ),
#else
      HQ         ( half * EQ ),
      QQ         ( quarter * EQ ),
#endif
      COEFF_POOL ( new falcON::pool(max(4u,np), grav::NCOEF*sizeof(real)) ),
      NC         ( 0 ),
      MAXNC      ( 0 ) {}
    //--------------------------------------------------------------------------
    ~GravKernBase() {
      if(COEFF_POOL) delete COEFF_POOL;
    }
    //--------------------------------------------------------------------------
    void give_coeffs(grav::cell_pter const&C) const {
      if(COEFF_POOL && !hasCoeffs(C)) {
	register grav::Cset*X = static_cast<grav::Cset*>(COEFF_POOL->alloc());
	X->set_zero();
	C->setCoeffs(X);
	++NC;
      }
    }
    //--------------------------------------------------------------------------
    void take_coeffs(grav::cell_pter const&C) const {
      if(COEFF_POOL && hasCoeffs(C)) {
	COEFF_POOL->free(C->returnCoeffs());
	C->resetCoeffs();
	update_max(MAXNC,NC--);
      }
    }
    //--------------------------------------------------------------------------
#define ARGS__ grav::leaf_iter const&, unsigned const&, 	\
               grav::leaf_iter const&, unsigned const&
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
    void eval_grav    (grav::cell_iter const&, TaylorSeries const&) const;
    void eval_grav_all(grav::cell_iter const&, TaylorSeries const&) const;
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
  public:
    int      const&coeffs_used() const { return MAXNC; }
    unsigned       chunks_used() const {
      return COEFF_POOL?  COEFF_POOL->N_chunks() : 0u; }
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::GravKern                                                     
  //                                                                            
  // This class implements the direct summation and approximate computation     
  // of gravity between tree nodes.                                             
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class GravKern : public GravKernBase
  {
#ifdef falcON_SSE_CODE
    //--------------------------------------------------------------------------
    // blocking of 4 cell-node interactions                                     
    //--------------------------------------------------------------------------
    void flush_acl() const;
    void flush_acc() const;
#endif
    //--------------------------------------------------------------------------
    // main purpose methods                                                     
    //--------------------------------------------------------------------------
  protected:
    GravKern(kern_type const&k,                    // I: type of kernel         
	     real      const&e,                    // I: softening length       
#ifdef falcON_INDI
	     bool      const&s,                    // I: type of softening      
#endif
	     unsigned  const&np) :                 // I: initial pool size      
#ifdef falcON_INDI
      GravKernBase(k,e,s,np) {}
#else
      GravKernBase(k,e,np) {}
#endif
    //--------------------------------------------------------------------------
    // single leaf-leaf interaction                                             
    //--------------------------------------------------------------------------
    void single(grav::leaf_iter const&, grav::leaf_iter const&) const;
    //--------------------------------------------------------------------------
    // cell-leaf interaction via direct summation                               
    //--------------------------------------------------------------------------
    void direct(grav::cell_iter const&, grav::leaf_iter const&) const;
    //--------------------------------------------------------------------------
    // cell-cell interaction via direct summation                               
    //--------------------------------------------------------------------------
    void direct(grav::cell_iter const&, grav::cell_iter const&) const;
    //--------------------------------------------------------------------------
    // cell-self interaction via direct summation                               
    //--------------------------------------------------------------------------
    void direct(grav::cell_iter const&) const;
    //--------------------------------------------------------------------------
    // cell-leaf interaction via approximation                                  
    //--------------------------------------------------------------------------
    void approx(grav::cell_iter const&, grav::leaf_iter const&,
		vect&, real const&) const;
    //--------------------------------------------------------------------------
    // cell-cell interaction via approximation                                  
    //--------------------------------------------------------------------------
    void approx(grav::cell_iter const&, grav::cell_iter const&,
		vect&, real const&) const;
    //--------------------------------------------------------------------------
    void flush_buffers() const {
#ifdef falcON_SSE_CODE
      if(!ACL.is_empty()) flush_acl();
      if(!ACC.is_empty()) flush_acc();
#endif
    }
    //--------------------------------------------------------------------------
    // destruction: flush buffers                                               
    //--------------------------------------------------------------------------
    ~GravKern() { flush_buffers(); }
    //--------------------------------------------------------------------------
  public:
    const real&current_eps  ()     const { return EPS; }
    const real&current_epsq ()     const { return EQ; }
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::GravKernAll                                                  
  //                                                                            
  // Like GravKern, except that all cells and leafs are assumed active.         
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class GravKernAll : public GravKernBase
  {
#ifdef falcON_SSE_CODE
    //--------------------------------------------------------------------------
    // blocking of 4 cell-node interactions                                     
    //--------------------------------------------------------------------------
    void flush_acl() const;
    void flush_acc() const;
#endif
    //--------------------------------------------------------------------------
    // main purpose methods                                                     
    //--------------------------------------------------------------------------
  protected:
    GravKernAll(kern_type const&k,                 // I: type of kernel         
		real      const&e,                 // I: softening length       
#ifdef falcON_INDI
		bool      const&s,                 // I: type of softening      
#endif
		unsigned  const&np) :              // I: initial pool size      
#ifdef falcON_INDI
      GravKernBase(k,e,s,np) {}
#else
      GravKernBase(k,e,np) {}
#endif
    //--------------------------------------------------------------------------
    // single leaf-leaf interaction                                             
    //--------------------------------------------------------------------------
    void single(grav::leaf_iter const&, grav::leaf_iter const&) const;
    //--------------------------------------------------------------------------
    // cell-leaf interaction via direct summation                               
    //--------------------------------------------------------------------------
    void direct(grav::cell_iter const&, grav::leaf_iter const&) const;
    //--------------------------------------------------------------------------
    // cell-cell interaction via direct summation                               
    //--------------------------------------------------------------------------
    void direct(grav::cell_iter const&, grav::cell_iter const&) const;
    //--------------------------------------------------------------------------
    // cell-self interaction via direct summation                               
    //--------------------------------------------------------------------------
    void direct(grav::cell_iter const&) const;
    //--------------------------------------------------------------------------
    // cell-leaf interaction via approximation                                  
    //--------------------------------------------------------------------------
    void approx(grav::cell_iter const&, grav::leaf_iter const&,
		vect&, real const&) const;
    //--------------------------------------------------------------------------
    // cell-cell interaction via approximation                                  
    //--------------------------------------------------------------------------
    void approx(grav::cell_iter const&, grav::cell_iter const&,
		vect&, real const&) const;
    //--------------------------------------------------------------------------
    void flush_buffers() const {
#ifdef falcON_SSE_CODE
      if(!ACL.is_empty()) flush_acl();
      if(!ACC.is_empty()) flush_acc();
#endif
    }
    //--------------------------------------------------------------------------
    // destruction: flush buffers                                               
    //--------------------------------------------------------------------------
    ~GravKernAll() { flush_buffers(); }
    //--------------------------------------------------------------------------
  public:
    const real&current_eps  ()     const { return EPS; }
    const real&current_epsq ()     const { return EQ; }
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // inline functions                                                         //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  inline void GravKern::direct(grav::cell_iter const&CA,
			       grav::cell_iter const&CB) const
  {
#define ARGS_A CA.begin_leafs(),NA,CB.begin_leafs(),NB
#define ARGS_B CB.begin_leafs(),NB,CA.begin_leafs(),NA
    const unsigned NA=number(CA), NB=number(CB);
    if(NA%4 > NB%4) {
      if       (al_active(CA))
	if     (al_active(CB)) many_AA(ARGS_A);    // active: all  A, all  B    
	else if(is_active(CB)) many_AS(ARGS_A);    // active: all  A, some B    
	else                   many_AN(ARGS_A);    // active: all  A, no   B    
      else if  (is_active(CA))
	if     (al_active(CB)) many_SA(ARGS_A);    // active: some A, all  B    
	else if(is_active(CB)) many_SS(ARGS_A);    // active: some A, some B    
	else                   many_SN(ARGS_A);    // active: some A, no   B    
      else
	if     (al_active(CB)) many_NA(ARGS_A);    // active: no   A, all  B    
	else if(is_active(CB)) many_NS(ARGS_A);    // active: no   A, some B    
    } else {
      if       (al_active(CB))
	if     (al_active(CA)) many_AA(ARGS_B);    // active: all  B, all  A    
	else if(is_active(CA)) many_AS(ARGS_B);    // active: all  B, some A    
	else                   many_AN(ARGS_B);    // active: all  B, no   A    
      else if(is_active(CB))
	if     (al_active(CA)) many_SA(ARGS_B);    // active: some B, all  A    
	else if(is_active(CA)) many_SS(ARGS_B);    // active: some B, some A    
	else                   many_SN(ARGS_B);    // active: some B, no   A    
      else
	if     (al_active(CA)) many_NA(ARGS_B);    // active: no   B, all  A    
	else if(is_active(CA)) many_NS(ARGS_B);    // active: no   B, some A    
    }
  }
  //----------------------------------------------------------------------------
  inline void GravKernAll::direct(grav::cell_iter const&CA,
				  grav::cell_iter const&CB) const
  {
    const unsigned  NA=number(CA), NB=number(CB);
    if(NA%4 > NB%4) many_AA(ARGS_A);
    else            many_AA(ARGS_B);
  }
#undef ARGS_A
#undef ARGS_B
  //////////////////////////////////////////////////////////////////////////////
#ifdef falcON_SSE_CODE
  //----------------------------------------------------------------------------
  // blocking of approximate interactions: just fill the blocks                 
  //----------------------------------------------------------------------------
  inline void GravKern::approx(grav::cell_iter const&A,
			       grav::leaf_iter const&B,
			       vect                 &dR,
			       real            const&Rq) const
  {
#if falcON_INDI
    if(INDI_SOFT) ACL.load_i(A,B,dR,Rq); else
#endif
      ACL.load_g(A,B,dR,Rq+EQ);
    if(ACL.is_full()) flush_acl();
  }
  //----------------------------------------------------------------------------
  inline void GravKern::approx(grav::cell_iter const&A,
			       grav::cell_iter const&B,
			       vect                 &dR,
			       real            const&Rq) const
  {
#if falcON_INDI
    if(INDI_SOFT) ACC.load_i(A,B,dR,Rq); else
#endif
      ACC.load_g(A,B,dR,Rq+EQ);
    if(ACC.is_full()) flush_acc();
  }
  //----------------------------------------------------------------------------
  inline void GravKernAll::approx(grav::cell_iter const&A,
				  grav::leaf_iter const&B,
				  vect                 &dR,
				  real            const&Rq) const
  {
#if falcON_INDI
    if(INDI_SOFT) ACL.load_i(A,B,dR,Rq); else
#endif
      ACL.load_g(A,B,dR,Rq+EQ);
    if(ACL.is_full()) flush_acl();
  }
  //----------------------------------------------------------------------------
  inline void GravKernAll::approx(grav::cell_iter const&A,
				  grav::cell_iter const&B,
				  vect                 &dR,
				  real            const&Rq) const
  {
#if falcON_INDI
    if(INDI_SOFT) ACC.load_i(A,B,dR,Rq); else
#endif
      ACC.load_g(A,B,dR,Rq+EQ);
    if(ACC.is_full()) flush_acc();
  }
#endif
  //////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_included_kernel.h  
