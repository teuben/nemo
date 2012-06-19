// -*- C++ -*-
// /////////////////////////////////////////////////////////////////////////////
//
/// \file    src/public/kernel.cc
//
/// \brief   implements inc/public/kernel.h
/// \author  Walter Dehnen
/// \date    2000-2010,2012
//
// /////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008-2010  Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 675
// Mass Ave, Cambridge, MA 02139, USA.
//
// /////////////////////////////////////////////////////////////////////////////
#include <public/types.h>
#include <public/kernel.h>
#include <public/tensor_set.h>
#include <utils/WDMath.h>

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Macros facilitating the Cell-Leaf and Cell-Cell interactions               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifdef falcON_SSE_CODE
#  define __ARG_D D,J
#else
#  define __ARG_D D
#endif

#define CellLeaf(A,B,D,J,R) {						       \
  grav::Cset F;                                    /* to hold F^(n)        */  \
  if(is_active(A)) {                               /* IF A is active       */  \
    set_dPhi(F,R,__ARG_D);                         /*   F^(n) = d^nPhi/dR^n*/  \
    add_C_B2C(A->Coeffs(),F);                      /*   C_A   = ...        */  \
    if(is_active(B)) {                             /*   IF B is active, too*/  \
      F.flip_sign_odd();                           /*     flip sign:F^(odd)*/  \
      add_C_C2B(B->Coeffs(),F,A->poles());         /*     C_B   = ...      */  \
    }                                              /*   ENDIF              */  \
  } else if(is_active(B)) {                        /* ELIF B is active     */  \
    R.negate();                                    /*   flip sign: R       */  \
    set_dPhi(F,R,__ARG_D);                         /*   F^(n) = d^nPhi/dR^n*/  \
    add_C_C2B(B->Coeffs(),F,A->poles());           /*   C_B   = ...        */  \
  }                                                /* ENDIF                */  \
}

#define CellLeafAll(A,B,D,J,R) {					       \
  grav::Cset F;                                    /* to hold F^(n)        */  \
  set_dPhi(F,R,__ARG_D);                           /* F^(n) = d^nPhi/dR^n  */  \
  add_C_B2C(A->Coeffs(),F);                        /* C_A   = ...          */  \
  F.flip_sign_odd();                               /* F^(n) = d^nPhi/dR^n  */  \
  add_C_C2B(B->Coeffs(),F,A->poles());             /* C_B   = ...          */  \
}

#define CellCell(A,B,D,J,R) {						       \
  grav::Cset F;                                    /* to hold F^(n)        */  \
  if(is_active(A)) {                               /* IF A is active       */  \
    set_dPhi(F,R,__ARG_D);                         /*   F^(n) = d^nPhi/dR^n*/  \
    add_C_C2C(A->Coeffs(),F,B->poles());           /*   C_A   = ...        */  \
    if(is_active(B)) {                             /*   IF B is active, too*/  \
      F.flip_sign_odd();                           /*     flip sign:F^(odd)*/  \
      add_C_C2C(B->Coeffs(),F,A->poles());         /*     C_B   = ...      */  \
    }                                              /*   ENDIF              */  \
  } else if(is_active(B)) {                        /* ELIF B is active     */  \
    R.negate();                                    /*   flip sign: R       */  \
    set_dPhi(F,R,__ARG_D);                         /*   F^(n) = d^nPhi/dR^n*/  \
    add_C_C2C(B->Coeffs(),F,A->poles());           /*   C_B   = ...        */  \
  }                                                /* ENDIF                */  \
}

#define CellCellAll(A,B,D,J,R) {					       \
  grav::Cset F;                                    /* to hold F^(n)        */  \
  set_dPhi(F,R,__ARG_D);                           /* F^(n) = d^nPhi/dR^n  */  \
  add_C_C2C(A->Coeffs(),F,B->poles());             /* C_A   = ...          */  \
  F.flip_sign_odd();                               /* F^(n) = d^nPhi/dR^n  */  \
  add_C_C2C(B->Coeffs(),F,A->poles());             /* C_B   = ...          */  \
}

////////////////////////////////////////////////////////////////////////////////

using namespace falcON;
using namespace falcON::grav;
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class falcON::TaylorSeries                                                   
//                                                                              
////////////////////////////////////////////////////////////////////////////////
inline void TaylorSeries::
shift_and_add(const grav::cell*const&c) {          // I: cell & its coeffs      
  if(hasCoeffs(c)) {                               // IF(cell has had iaction)  
    vect dX = cofm(c) - X;                         //   vector to shift by      
    if(dX != zero && C != zero) {                  //   IF(dX != 0 AND C != 0)  
      shift_by(C,dX);                              //     shift expansion       
    }                                              //   ENDIF                   
    X = cofm(c);                                   //   set X to new position   
    C.add_times(Coeffs(c), one/mass(c));           //   add cell's coeffs in    
  }                                                // ENDIF                     
}
//------------------------------------------------------------------------------
inline void TaylorSeries::
extract_grav(leaf_iter const&L) const {            // I: leaf to get grav to    
  eval_expn(L->Coeffs(),C,cofm(L)-X);              // evaluate expansion        
}
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class falcON::GravKernBase                                                   
//                                                                              
////////////////////////////////////////////////////////////////////////////////
void GravKernBase::eval_grav(cell_iter    const&C,
			     TaylorSeries const&T) const
{
  TaylorSeries G(T);                               // G = copy of T             
  G.shift_and_add(C);                              // shift G; G+=T_C           
  take_coeffs(C);                                  // free memory: C's coeffs   
  LoopLeafKids(cell_iter,C,l) if(is_active(l)) {   // LOOP C's active leaf kids 
    l->normalize_grav();                           //   pot,acc/=mass           
    if(!is_empty(G)) G.extract_grav(l);            //   add pot,acc due to G    
  }                                                // END LOOP                  
  LoopCellKids(cell_iter,C,c) if(is_active(c))     // LOOP C's active cell kids 
    eval_grav(c,G);                                //   recursive call          
}
//------------------------------------------------------------------------------
void GravKernBase::eval_grav_all(cell_iter    const&C,
				 TaylorSeries const&T) const
{
  TaylorSeries G(T);                               // G = copy of T             
  G.shift_and_add(C);                              // shift G; G+=T_C           
  take_coeffs(C);                                  // free memory: C's coeffs   
  LoopLeafKids(cell_iter,C,l) {                    // LOOP C's leaf kids        
    l->normalize_grav();                           //   pot,acc/=mass           
    if(!is_empty(G)) G.extract_grav(l);            //   add pot,acc due to G    
  }                                                // END LOOP                  
  LoopCellKids(cell_iter,C,c)                      // LOOP C's cell kids        
    eval_grav_all(c,G);                            //   recursive call          
}
//------------------------------------------------------------------------------
real GravKernBase::Psi(kern_type k, real Xq, real Eq)
{
  switch(k) {
  case p1: {
    real   x = one/(Xq+Eq), d0=sqrt(x), d1=d0*x, hq=half*Eq;
    return d0 + hq*d1;
  }
  case p2: {
    real   x = one/(Xq+Eq), d0=sqrt(x), d1=d0*x, d2=3*d1*x, hq=half*Eq;
    return d0 + hq*(d1+hq*d2);
  }
  case p3: {
    real   x = one/(Xq+Eq), d0=sqrt(x), d1=d0*x, d2=3*d1*x, d3=5*d2*x, 
      hq=half*Eq;
    return d0 + hq*(d1+half*hq*(d2+hq*d3));
  }
  default:
      return WDutils::invsqrt(Xq+Eq);
  }
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// 1. Single body-body interaction (these are extremely rare)                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#define P0(MUM)					\
  x  = one/(Rq+EQ);				\
  D0 = MUM*sqrt(x);				\
  R *= D0*x;
//------------------------------------------------------------------------------
#define P1(MUM)					\
  x  = one/(Rq+EQ);				\
  D0 = MUM*sqrt(x);				\
  register real D1 = D0*x;			\
  register real hq = half*EQ;			\
  D0+= hq*D1;					\
  D1+= hq*3*D1*x;				\
  R *= D1;
//------------------------------------------------------------------------------
#define P2(MUM)						\
  x  = one/(Rq+EQ);					\
  D0 = MUM*sqrt(x);					\
  register real D1 = D0*x, D2= 3*D1*x, D3= 5*D2*x;	\
  register real hq = half*EQ;				\
  D0+= hq*(D1+hq*D2);					\
  D1+= hq*(D2+hq*D3);					\
  R *= D1;
//------------------------------------------------------------------------------
#define P3(MUM)							\
  x  = one/(Rq+EQ);						\
  D0 = MUM*sqrt(x);						\
  register real D1 = D0*x, D2= 3*D1*x, D3= 5*D2*x, D4= 7*D3*x;	\
  register real hq = half*EQ;					\
  register real qq = half*hq;					\
  D0+= ((hq*D3+D2)*qq+D1)*hq;					\
  D1+= ((hq*D4+D3)*qq+D2)*hq;					\
  R *= D1;
//------------------------------------------------------------------------------
#define P0_I(MUM)				\
  EQ = square(eph(A)+eph(B));			\
  P0(MUM)	
//------------------------------------------------------------------------------
#define P1_I(MUM)				\
  EQ = square(eph(A)+eph(B));			\
  P1(MUM)	
//------------------------------------------------------------------------------
#define P2_I(MUM)				\
  EQ = square(eph(A)+eph(B));			\
  P2(MUM)	
//------------------------------------------------------------------------------
#define P3_I(MUM)				\
  EQ = square(eph(A)+eph(B));			\
  P3(MUM)	
//==============================================================================
void GravKern::single(leaf_iter const &A, leaf_iter const&B) const
{
  vect R  = cofm(A)-cofm(B);
  real Rq = norm(R),x,D0;
  if(INDI_SOFT)
    switch(KERN) {
    case p1: { P1_I(mass(A)*mass(B)) } break;
    case p2: { P2_I(mass(A)*mass(B)) } break;
    case p3: { P3_I(mass(A)*mass(B)) } break;
    default: { P0_I(mass(A)*mass(B)) } break;
    }
  else
    switch(KERN) {
    case p1: { P1(mass(A)*mass(B)) } break;
    case p2: { P2(mass(A)*mass(B)) } break;
    case p3: { P3(mass(A)*mass(B)) } break;
    default: { P0(mass(A)*mass(B)) } break;
    }
  if(is_active(A)) { A->pot()-=D0; A->acc()-=R; }
  if(is_active(B)) { B->pot()-=D0; B->acc()+=R; }
}
//------------------------------------------------------------------------------
void GravKernAll::single(leaf_iter const &A, leaf_iter const&B) const
{
  vect R  = cofm(A)-cofm(B);
  real Rq = norm(R),x,D0;
  if(INDI_SOFT)
    switch(KERN) {
    case p1: { P1_I(mass(A)*mass(B)) } break;
    case p2: { P2_I(mass(A)*mass(B)) } break;
    case p3: { P3_I(mass(A)*mass(B)) } break;
    default: { P0_I(mass(A)*mass(B)) } break;
    }
  else
    switch(KERN) {
    case p1: { P1(mass(A)*mass(B)) } break;
    case p2: { P2(mass(A)*mass(B)) } break;
    case p3: { P3(mass(A)*mass(B)) } break;
    default: { P0(mass(A)*mass(B)) } break;
    }
  A->pot()-=D0; A->acc()-=R;
  B->pot()-=D0; B->acc()+=R;
}
//------------------------------------------------------------------------------
#undef P0
#undef P1
#undef P2
#undef P3
#undef P0_I
#undef P1_I
#undef P2_I
#undef P3_I
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// 2. Cell-Node interactions with SSE instructions                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#ifdef falcON_SSE_CODE
#  include <proper/kernel_SSE.h>
#else
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// 2. Cell-Node interactions without SSE instructions                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//==============================================================================
//                                                                              
// macros for computing the gravity                                             
//                                                                              
// DSINGL       called after loading a leaf-leaf interaction                    
//                                                                              
//==============================================================================
//                                                                              
// for direct summation code, we assume:                                        
// - real  D0          contains Mi*Mj   on input                                
// - real  D1          contains R^2+e^2 on input                                
// - real  EP[3]       are set to e^2, e^2/2 and e^2/4 for global softening     
// - real  EQ          is set to e^2                   for individual softening 
//                                                                              
//==============================================================================
//                                                                              
// NOTE that we changed the definition of the D_n by a sign:                    
//                                                                              
// D_n = (-1/r d/dr)^n g(r) at r=|R|                                            
//                                                                              
//==============================================================================
# define DSINGL_P0_G				\
  register real					\
  XX  = one/D1;					\
  D0 *= sqrt(XX);				\
  D1  = XX * D0;
# define DSINGL_P0_I				\
  DSINGL_P0_G
//------------------------------------------------------------------------------
# define DSINGL_P1_G				\
  register real					\
  XX  = one/D1;					\
  D0 *= sqrt(XX);				\
  D1  = XX * D0;				\
  XX *= 3  * D1;          /* XX == T2 */	\
  D0 += HQ * D1;				\
  D1 += HQ * XX;
# define DSINGL_P1_I				\
  HQ  = half*EQ;				\
  DSINGL_P1_G
//------------------------------------------------------------------------------
# define DSINGL_P2_G				\
  register real					\
  XX  = one/D1;					\
  D0 *= sqrt(XX);				\
  D1  = XX * D0;				\
  register real					\
  D2  = 3 * XX * D1;				\
  XX *= 5 * D2;           /* XX == T3 */	\
  D0 += HQ*(D1+HQ*D2);				\
  D1 += HQ*(D2+HQ*XX);
# define DSINGL_P2_I				\
  HQ = half*EQ;					\
  DSINGL_P2_G
//------------------------------------------------------------------------------
# define DSINGL_P3_G				\
  register real					\
  XX  = one/D1;					\
  D0 *= sqrt(XX);				\
  D1  =     XX * D0;				\
  register real					\
  D2  = 3 * XX * D1;				\
  register real					\
  D3  = 5 * XX * D2;				\
  XX *= 7 * D3;           /* XX == T4 */	\
  D0 += HQ*(D1+QQ*(D2+HQ*D3));			\
  D1 += HQ*(D2+QQ*(D3+HQ*XX));
# define DSINGL_P3_I				\
  HQ = half*EQ;					\
  QQ = half*HQ;					\
  DSINGL_P3_G
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// 2.1 direct summation of many-body interactions                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// macros for organizing the gravity computation via direct summation.          
//                                                                              
// we have to take care for active and non-active                               
//                                                                              
//==============================================================================
//                                                                              
// These macros assume:                                                         
//                                                                              
// - vect    X0     position of left leaf                                       
// - real    M0     mass of left leaf                                           
// - vect    F0     force for left leaf (if used)                               
// - real    P0     potential for left leaf (if used)                           
// - vect    dR     to be filled with  X0-X_j                                   
// - real    D0     to be filled with  M0*M_j                                   
// - real    D1     to be filled with  norm(dR[J])+eps^2  on loading            
//                                                                              
//==============================================================================
#define LOAD_G					\
  dR = X0 - cofm(B);				\
  D1 = norm(dR) + EQ;				\
  D0 = M0 * mass(B);
//------------------------------------------------------------------------------
#define LOAD_I					\
  dR = X0 - cofm(B);				\
  EQ = square(E0+eph(B));			\
  D1 = norm(dR) + EQ;				\
  D0 = M0 * mass(B);
//------------------------------------------------------------------------------
#define PUT_LEFT				\
  dR *= D1;					\
  P0 -= D0;					\
  F0 -= dR;
//------------------------------------------------------------------------------
#define PUT_RGHT				\
  dR       *= D1;				\
  B->pot() -= D0;				\
  B->acc() += dR;
//------------------------------------------------------------------------------
#define PUT_BOTH				\
  dR       *= D1;				\
  P0       -= D0;				\
  F0       -= dR;				\
  B->pot() -= D0;				\
  B->acc() += dR;
//------------------------------------------------------------------------------
#define PUT_SOME				\
  dR *= D1;					\
  P0 -= D0;					\
  F0 -= dR;					\
  if(is_active(B)) {				\
    B->pot() -= D0;				\
    B->acc() += dR;				\
  }
//------------------------------------------------------------------------------
#define GRAV_ALL(LOAD,DSINGL,PUT)		\
for(register leaf_iter B=B0; B!=BN; ++B) {	\
  LOAD						\
  DSINGL					\
  PUT						\
} 
//------------------------------------------------------------------------------
#define GRAV_FEW(LOAD,DSINGL)			\
for(register leaf_iter B=B0; B!=BN; ++B)	\
  if(is_active(B)) {				\
    LOAD					\
    DSINGL					\
    PUT_RGHT					\
  } 
//------------------------------------------------------------------------------
#define START_G					\
  const    real      M0=mass(A);		\
  const    vect      X0=cofm(A);		\
           vect      dR;			\
  register real      D0,D1;
//------------------------------------------------------------------------------
#define START_I					\
  const    real      E0=eph(A);			\
  const    real      M0=mass(A);		\
  const    vect      X0=cofm(A);		\
           vect      dR;			\
  register real      D0,D1;
//==============================================================================
// now defining auxiliary inline functions for the computation of  N            
// interactions. There are the following 10 cases:                              
// - each for the cases YA, YS, YN, NA, NS                                      
// - each for global and individual softening                                   
//==============================================================================
namespace {
  using namespace falcON; using namespace falcON::grav;
  //////////////////////////////////////////////////////////////////////////////
#define DIRECT(START,LOAD,DSINGL)				\
    static void many_YA(ARGS) {					\
      START; register real P0(zero); vect F0(zero);		\
      GRAV_ALL(LOAD,DSINGL,PUT_BOTH)				\
      A->pot()+=P0;  A->acc()+=F0;				\
    }								\
    static void many_YS(ARGS) {					\
      START; register real P0(zero); vect F0(zero);		\
      GRAV_ALL(LOAD,DSINGL,PUT_SOME)				\
      A->pot()+=P0; A->acc()+=F0;				\
    }								\
    static void many_YN(ARGS) {					\
      START; register real P0(zero); vect F0(zero);		\
      GRAV_ALL(LOAD,DSINGL,PUT_LEFT)				\
      A->pot()+=P0; A->acc()+=F0;				\
    }								\
    static void many_NA(ARGS) {					\
      START;							\
      GRAV_ALL(LOAD,DSINGL,PUT_RGHT)				\
    }								\
    static void many_NS(ARGS) {					\
      START;							\
      GRAV_FEW(LOAD,DSINGL)					\
    }
  //////////////////////////////////////////////////////////////////////////////
  template<kern_type, bool> struct __direct;
  //----------------------------------------------------------------------------
#define ARGS					\
  leaf_iter const&A,				\
  leaf_iter const&B0,				\
  leaf_iter const&BN,				\
  real&EQ, real&, real& 

  template<> struct __direct<p0,0> {
    DIRECT(START_G,LOAD_G,DSINGL_P0_G);
  };
  template<> struct __direct<p0,1> {
    DIRECT(START_I,LOAD_I,DSINGL_P0_I);
  };
  //----------------------------------------------------------------------------
#undef  ARGS
#define ARGS					\
  leaf_iter const&A,				\
  leaf_iter const&B0,				\
  leaf_iter const&BN,				\
  real&EQ, real&HQ, real& 

  template<> struct __direct<p1,0> {
    DIRECT(START_G,LOAD_G,DSINGL_P1_G);
  };
  template<> struct __direct<p1,1> {
    DIRECT(START_I,LOAD_I,DSINGL_P1_I);
  };
  //----------------------------------------------------------------------------
  template<> struct __direct<p2,0> {
    DIRECT(START_G,LOAD_G,DSINGL_P2_G);
  };
  template<> struct __direct<p2,1> {
    DIRECT(START_I,LOAD_I,DSINGL_P2_I);
  };
  //----------------------------------------------------------------------------
#undef  ARGS
#define ARGS					\
  leaf_iter const&A,				\
  leaf_iter const&B0,				\
  leaf_iter const&BN,				\
  real&EQ, real&HQ, real&QQ 

  template<> struct __direct<p3,0> {
    DIRECT(START_G,LOAD_G,DSINGL_P3_G);
  };
  template<> struct __direct<p3,1> {
    DIRECT(START_I,LOAD_I,DSINGL_P3_I);
  };
#undef LOAD_G
#undef LOAD_I
#undef START_G
#undef START_I
#undef DIRECT
  //////////////////////////////////////////////////////////////////////////////
  template<bool I> struct Direct {
    static void many_YA(kern_type KERN, ARGS) {
      switch(KERN) {
      case p1: __direct<p1,I>::many_YA(A,B0,BN,EQ,HQ,QQ); break;
      case p2: __direct<p2,I>::many_YA(A,B0,BN,EQ,HQ,QQ); break;
      case p3: __direct<p3,I>::many_YA(A,B0,BN,EQ,HQ,QQ); break;
      default: __direct<p0,I>::many_YA(A,B0,BN,EQ,HQ,QQ); break;
      } }
    static void many_YS(kern_type KERN, ARGS) {
      switch(KERN) {
      case p1: __direct<p1,I>::many_YS(A,B0,BN,EQ,HQ,QQ); break;
      case p2: __direct<p2,I>::many_YS(A,B0,BN,EQ,HQ,QQ); break;
      case p3: __direct<p3,I>::many_YS(A,B0,BN,EQ,HQ,QQ); break;
      default: __direct<p0,I>::many_YS(A,B0,BN,EQ,HQ,QQ); break;
      } }
    static void many_YN(kern_type KERN, ARGS) {
      switch(KERN) {
      case p1: __direct<p1,I>::many_YN(A,B0,BN,EQ,HQ,QQ); break;
      case p2: __direct<p2,I>::many_YN(A,B0,BN,EQ,HQ,QQ); break;
      case p3: __direct<p3,I>::many_YN(A,B0,BN,EQ,HQ,QQ); break;
      default: __direct<p0,I>::many_YN(A,B0,BN,EQ,HQ,QQ); break;
      } }
    static void many_NA(kern_type KERN, ARGS) {
      switch(KERN) {
      case p1: __direct<p1,I>::many_NA(A,B0,BN,EQ,HQ,QQ); break;
      case p2: __direct<p2,I>::many_NA(A,B0,BN,EQ,HQ,QQ); break;
      case p3: __direct<p3,I>::many_NA(A,B0,BN,EQ,HQ,QQ); break;
      default: __direct<p0,I>::many_NA(A,B0,BN,EQ,HQ,QQ); break;
      } }
    static void many_NS(kern_type KERN, ARGS) {
      switch(KERN) {
      case p1: __direct<p1,I>::many_NS(A,B0,BN,EQ,HQ,QQ); break;
      case p2: __direct<p2,I>::many_NS(A,B0,BN,EQ,HQ,QQ); break;
      case p3: __direct<p3,I>::many_NS(A,B0,BN,EQ,HQ,QQ); break;
      default: __direct<p0,I>::many_NS(A,B0,BN,EQ,HQ,QQ); break;
      } }
  };
  //////////////////////////////////////////////////////////////////////////////
}                                                  // END: unnamed namespace    
#undef ARGS
//=============================================================================
// we now can define the cell-leaf and cell-self interaction via direct sums    
//==============================================================================
#define ARGS KERN,B,A.begin_leafs(),A.end_leaf_desc(),EQ,HQ,QQ
void GravKern::direct(cell_iter const&A, leaf_iter const&B) const
{
  if(INDI_SOFT)
    if(is_active(B)) {
      if     (al_active(A)) Direct<1>::many_YA(ARGS);
      else if(is_active(A)) Direct<1>::many_YS(ARGS);
      else                  Direct<1>::many_YN(ARGS);
    } else {
      if     (al_active(A)) Direct<1>::many_NA(ARGS);
      else if(is_active(A)) Direct<1>::many_NS(ARGS);
    }
  else
  {
    if(is_active(B)) {
      if     (al_active(A)) Direct<0>::many_YA(ARGS);
      else if(is_active(A)) Direct<0>::many_YS(ARGS);
      else                  Direct<0>::many_YN(ARGS);
    } else {
      if     (al_active(A)) Direct<0>::many_NA(ARGS);
      else if(is_active(A)) Direct<0>::many_NS(ARGS);
    }
  }
}
//------------------------------------------------------------------------------
void GravKernAll::direct(cell_iter const&A, leaf_iter const&B) const
{
  if(INDI_SOFT)
    Direct<1>::many_YA(ARGS);
  else
    Direct<0>::many_YA(ARGS);
}
#undef ARGS
//------------------------------------------------------------------------------
#define ARGS KERN,A,A+1,A+1+Nk,EQ,HQ,QQ
void GravKern::direct(cell_iter const&C) const
{
  const unsigned  N1 = number(C)-1;
  leaf_iter       A  = C.begin_leafs();
  if(INDI_SOFT)
    if(al_active(C))
      for(unsigned Nk=N1; Nk; --Nk,++A) Direct<1>::many_YA(ARGS);
    else
      for(unsigned Nk=N1; Nk; --Nk,++A)
	if(is_active(A))                Direct<1>::many_YS(ARGS);
	else                            Direct<1>::many_NS(ARGS);
  else
  {
    if(al_active(C))
      for(unsigned Nk=N1; Nk; --Nk,++A) Direct<0>::many_YA(ARGS);
    else
      for(unsigned Nk=N1; Nk; --Nk,++A)
	if(is_active(A))                Direct<0>::many_YS(ARGS);
	else                            Direct<0>::many_NS(ARGS);
  }
}
//------------------------------------------------------------------------------
void GravKernAll::direct(cell_iter const&C) const
{
  const unsigned  N1 = number(C)-1;
  leaf_iter       A  = C.begin_leafs();
  if(INDI_SOFT)
    for(unsigned Nk=N1; Nk; --Nk,++A) Direct<1>::many_YA(ARGS);
  else
    for(unsigned Nk=N1; Nk; --Nk,++A) Direct<0>::many_YA(ARGS);
}
#undef ARGS
//==============================================================================
// we now define non-inline functions for the computation of many direct        
// interactions between NA and NB leafs                                         
// there are 8 cases, depending on whether all, some, or none of either A or    
// B are active (case none,none is trivial).                                    
//                                                                              
// these functions are called by GravKern::direct(cell,cell), which is inline   
// in kernel.h, or by GravKern::flush_scc() below.                              
//==============================================================================
void GravKernBase::many_AA(leaf_iter const&A0, unsigned NA,
			   leaf_iter const&B0, unsigned NB) const
{
  const leaf_iter AN=A0+NA, BN=B0+NB;
  if(INDI_SOFT)
    for(leaf_iter A=A0; A!=AN; ++A) Direct<1>::many_YA(KERN,A,B0,BN,EQ,HQ,QQ);
  else
    for(leaf_iter A=A0; A!=AN; ++A) Direct<0>::many_YA(KERN,A,B0,BN,EQ,HQ,QQ);
}
//------------------------------------------------------------------------------
void GravKernBase::many_AS(leaf_iter const&A0, unsigned NA,
			   leaf_iter const&B0, unsigned NB) const
{
  const    leaf_iter AN=A0+NA, BN=B0+NB;
  if(INDI_SOFT)
    for(leaf_iter A=A0; A!=AN; ++A) Direct<1>::many_YS(KERN,A,B0,BN,EQ,HQ,QQ);
  else
    for(leaf_iter A=A0; A!=AN; ++A) Direct<0>::many_YS(KERN,A,B0,BN,EQ,HQ,QQ);
}
//------------------------------------------------------------------------------
void GravKernBase::many_AN(leaf_iter const&A0, unsigned NA,
			   leaf_iter const&B0, unsigned NB) const
{
  const leaf_iter AN=A0+NA, BN=B0+NB;
  if(INDI_SOFT)
    for(leaf_iter A=A0; A!=AN; ++A) Direct<1>::many_YN(KERN,A,B0,BN,EQ,HQ,QQ);
  else
    for(leaf_iter A=A0; A!=AN; ++A) Direct<0>::many_YN(KERN,A,B0,BN,EQ,HQ,QQ);
}
//------------------------------------------------------------------------------
void GravKernBase::many_SA(leaf_iter const&A0, unsigned NA,
			   leaf_iter const&B0, unsigned NB) const
{
  const leaf_iter AN=A0+NA, BN=B0+NB;
  if(INDI_SOFT) {
    for(leaf_iter A=A0; A!=AN; ++A)
      if(is_active(A)) Direct<1>::many_YA(KERN,A,B0,BN,EQ,HQ,QQ);
      else             Direct<1>::many_NA(KERN,A,B0,BN,EQ,HQ,QQ);
  } else
    for(leaf_iter A=A0; A!=AN; ++A)
      if(is_active(A)) Direct<0>::many_YA(KERN,A,B0,BN,EQ,HQ,QQ);
      else             Direct<0>::many_NA(KERN,A,B0,BN,EQ,HQ,QQ);
}
//------------------------------------------------------------------------------
void GravKernBase::many_SS(leaf_iter const&A0, unsigned NA,
			   leaf_iter const&B0, unsigned NB) const
{
  const leaf_iter AN=A0+NA, BN=B0+NB;
  if(INDI_SOFT) {
    for(leaf_iter A=A0; A!=AN; ++A)
      if(is_active(A)) Direct<1>::many_YS(KERN,A,B0,BN,EQ,HQ,QQ);
      else             Direct<1>::many_NS(KERN,A,B0,BN,EQ,HQ,QQ);
  } else
    for(leaf_iter A=A0; A!=AN; ++A) 
      if(is_active(A)) Direct<0>::many_YS(KERN,A,B0,BN,EQ,HQ,QQ);
      else             Direct<0>::many_NS(KERN,A,B0,BN,EQ,HQ,QQ);
}
//------------------------------------------------------------------------------
void GravKernBase::many_SN(leaf_iter const&A0, unsigned NA,
			   leaf_iter const&B0, unsigned NB) const
{
  const leaf_iter AN=A0+NA, BN=B0+NB;
  if(INDI_SOFT) {
    for(leaf_iter A=A0; A!=AN; ++A)
      if(is_active(A)) Direct<1>::many_YN(KERN,A,B0,BN,EQ,HQ,QQ);
  } else
    for(leaf_iter A=A0; A!=AN; ++A) 
      if(is_active(A)) Direct<0>::many_YN(KERN,A,B0,BN,EQ,HQ,QQ);
}
//------------------------------------------------------------------------------
void GravKernBase::many_NA(leaf_iter const&A0, unsigned NA,
			   leaf_iter const&B0, unsigned NB) const
{
  const leaf_iter AN=A0+NA, BN=B0+NB;
  if(INDI_SOFT)
    for(leaf_iter A=A0; A!=AN; ++A) Direct<1>::many_NA(KERN,A,B0,BN,EQ,HQ,QQ);
  else
    for(leaf_iter A=A0; A!=AN; ++A) Direct<0>::many_NA(KERN,A,B0,BN,EQ,HQ,QQ);
}
//------------------------------------------------------------------------------
void GravKernBase::many_NS(leaf_iter const&A0, unsigned NA,
			   leaf_iter const&B0, unsigned NB) const
{
  const leaf_iter AN=A0+NA, BN=B0+NB;
  if(INDI_SOFT)
    for(leaf_iter A=A0; A!=AN; ++A) Direct<1>::many_NS(KERN,A,B0,BN,EQ,HQ,QQ);
  else
    for(leaf_iter A=A0; A!=AN; ++A) Direct<0>::many_NS(KERN,A,B0,BN,EQ,HQ,QQ);
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// 2.2 approximate gravity                                                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace {
  using namespace falcON; using namespace falcON::grav;
  //////////////////////////////////////////////////////////////////////////////
#define LOAD_G					\
  real D[ND];					\
  real XX=one/(Rq+EQ);				\
  D[0] = mass(A)*mass(B);
  //////////////////////////////////////////////////////////////////////////////
#define LOAD_I					\
  real D[ND];					\
  EQ   = square(eph(A)+eph(B));			\
  real XX=one/(Rq+EQ);				\
  D[0] = mass(A)*mass(B);			\
  __setE<P>::s(EQ,HQ,QQ); 
  //////////////////////////////////////////////////////////////////////////////
#define ARGS_B					\
  cell_iter const&A,				\
  leaf_iter const&B,				\
  vect           &R,				\
  real const     &Rq,				\
  real&EQ, real&HQ, real&QQ
  //////////////////////////////////////////////////////////////////////////////
#define ARGS_C					\
  cell_iter const&A,				\
  cell_iter const&B,				\
  vect           &R,				\
  real const     &Rq,				\
  real&EQ, real&HQ, real&QQ
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class __setE<kern_type>                                                  //
  // class __block<kern_type, order>                                          //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<kern_type>          struct __setE;
  template<kern_type,int>      struct __block;
#define sv   static void
  //////////////////////////////////////////////////////////////////////////////
  // kern_type = p0                                                             
  //////////////////////////////////////////////////////////////////////////////
  template<> struct __setE<p0> {
    sv s(real,real&,real&) {
    } };
  //----------------------------------------------------------------------------
  template<> struct __block<p0,1> : public __setE<p0> {
    enum { ND=2 };
    sv b(real&X, real D[ND], real, real, real) {
      D[0] *= sqrt(X);
      D[1]  = X * D[0];
    } };
  template<int K> struct __block<p0,K> : public __setE<p0> {
    enum { ND=K+1, F=K+K-1 };
    sv b(real&X, real D[ND], real EQ, real HQ, real QQ) {
      __block<p0,K-1>::b(X,D,EQ,HQ,QQ);
      D[K] = int(F) * X * D[K-1];
    } };
  //////////////////////////////////////////////////////////////////////////////
  // kern_type = p1                                                             
  //////////////////////////////////////////////////////////////////////////////
  template<> struct __setE<p1> {
    sv s(real EQ, real&HQ, real&) {
      HQ = half * EQ;
    } };
  //----------------------------------------------------------------------------
  template<> struct __block<p1,1> : public __setE<p1> {
    enum { ND=3 };
    sv b(real&X, real D[ND], real, real HQ, real) {
      D[0] *= sqrt(X);
      D[1]  =     X * D[0];
      D[2]  = 3 * X * D[1];
      D[0] += HQ*D[1];
      D[1] += HQ*D[2];
    } };
  template<int K> struct __block<p1,K> : public __setE<p1> {
    enum { ND=K+2, F=K+K+1 };
    sv b(real&X, real D[ND], real EQ, real HQ, real QQ) {
      __block<p1,K-1>::b(X,D,EQ,HQ,QQ);
      D[K+1] = int(F) * X * D[K];
      D[K]  += HQ  * D[K+1];
    } };
  //////////////////////////////////////////////////////////////////////////////
  // kern_type = p2                                                             
  //////////////////////////////////////////////////////////////////////////////
  template<> struct __setE<p2> {
    sv s(real EQ, real&HQ, real&) {
      HQ = half * EQ;
    } };
  //----------------------------------------------------------------------------
  template<> struct __block<p2,1> : public __setE<p2> {
    enum { ND=4 };
    sv b(real&X, real D[ND], real, real HQ, real) {
      D[0] *= sqrt(X);
      D[1]  =     X * D[0];
      D[2]  = 3 * X * D[1];
      D[3]  = 5 * X * D[2];
      D[0] += HQ*(D[1]+HQ*D[2]);
      D[1] += HQ*(D[2]+HQ*D[3]);
    } };
  template<int K> struct __block<p2,K> : public __setE<p2> {
    enum { ND=K+3, F=K+K+3 };
    sv b(real&X, real D[ND], real EQ, real HQ, real QQ) {
      __block<p2,K-1>::b(X,D,EQ,HQ,QQ);
      D[K+2] = int(F) * X * D[K+1];
      D[K]  += HQ*(D[K+1]+HQ*D[K+2]);
    } };
  //////////////////////////////////////////////////////////////////////////////
  // kern_type = p3                                                             
  //////////////////////////////////////////////////////////////////////////////
  template<> struct __setE<p3> {
    sv s(real EQ, real&HQ, real&QQ) {
      HQ = half * EQ;
      QQ = half * QQ;
    } };
  //----------------------------------------------------------------------------
  template<> struct __block<p3,1> : public __setE<p3> {
    enum { ND=5 };
    sv b(real&X, real D[ND], real, real HQ, real QQ) {
      D[0] *= sqrt(X);
      D[1]  =     X * D[0];
      D[2]  = 3 * X * D[1];
      D[3]  = 5 * X * D[2];
      D[4]  = 7 * X * D[3];
      D[0] += HQ*(D[1]+QQ*(D[2]+HQ*D[3]));
      D[1] += HQ*(D[2]+QQ*(D[3]+HQ*D[4]));
    } };
  template<int K> struct __block<p3,K> : public __setE<p3> {
    enum { ND=K+4, F=K+K+5 };
    sv b(real&X, real D[ND], real EQ, real HQ, real QQ) {
      __block<p3,K-1>::b(X,D,EQ,HQ,QQ);
      D[K+3] = int(F) * X * D[K+2];
      D[K]  += HQ*(D[K+1]+QQ*(D[K+2]+HQ*D[K+3]));
    } };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class kernel<kern_type, order, all, indi_soft>                           //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<kern_type,int,bool,bool=0> struct kernel;
  //////////////////////////////////////////////////////////////////////////////
  template<kern_type P, int K> struct kernel<P,K,0,0> : private __block<P,K> {
    enum { ND = __block<P,K>::ND };
    sv a(ARGS_B) { LOAD_G kernel::b(XX,D,EQ,HQ,QQ); CellLeaf(A,B,D,0,R); }
    sv a(ARGS_C) { LOAD_G kernel::b(XX,D,EQ,HQ,QQ); CellCell(A,B,D,0,R); }
  };
  template<kern_type P, int K> struct kernel<P,K,0,1> : private __block<P,K> {
    enum { ND = __block<P,K>::ND };
    sv a(ARGS_B) { LOAD_I kernel::b(XX,D,EQ,HQ,QQ); CellLeaf(A,B,D,0,R); }
    sv a(ARGS_C) { LOAD_I kernel::b(XX,D,EQ,HQ,QQ); CellCell(A,B,D,0,R); }
  };
  template<kern_type P, int K> struct kernel<P,K,1,0> : private __block<P,K> {
    enum { ND = __block<P,K>::ND };
    sv a(ARGS_B) { LOAD_G kernel::b(XX,D,EQ,HQ,QQ); CellLeafAll(A,B,D,0,R); }
    sv a(ARGS_C) { LOAD_G kernel::b(XX,D,EQ,HQ,QQ); CellCellAll(A,B,D,0,R); }
  };
#if defined(__GNUC__) && (__GNUC__ == 4) && (__GNUC_MINOR__ == 1) 
  // gcc 4.1.2 gives crashing code, if this is inlined.
  template<kern_type P, int K> struct kernel<P,K,1,1> : private __block<P,K> {
    enum { ND = __block<P,K>::ND };
    sv a(ARGS_B);
    sv a(ARGS_C);
  };
  template<kern_type P, int K> void kernel<P,K,1,1>::a(ARGS_B)
    { LOAD_I kernel::b(XX,D,EQ,HQ,QQ); CellLeafAll(A,B,D,0,R); }
  template<kern_type P, int K> void kernel<P,K,1,1>::a(ARGS_C)
    { LOAD_I kernel::b(XX,D,EQ,HQ,QQ); CellCellAll(A,B,D,0,R); }
#else
  template<kern_type P, int K> struct kernel<P,K,1,1> : private __block<P,K> {
    enum { ND = __block<P,K>::ND };
    sv a(ARGS_B) { LOAD_I kernel::b(XX,D,EQ,HQ,QQ); CellLeafAll(A,B,D,0,R); }
    sv a(ARGS_C) { LOAD_I kernel::b(XX,D,EQ,HQ,QQ); CellCellAll(A,B,D,0,R); }
  };
#endif
#undef sv
  //////////////////////////////////////////////////////////////////////////////
}                                                  // END: unnamed namespace    
////////////////////////////////////////////////////////////////////////////////
#define ARGS A,B,R,Rq,EQ,HQ,QQ
void GravKern::approx(cell_iter const&A, leaf_iter const&B,
		      vect           &R, real      Rq) const
{
  if(is_active(A)) give_coeffs(A);
  if(INDI_SOFT) 
    switch(KERN) {
    case p1: kernel<p1,ORDER,0,1>::a(ARGS); break;
    case p2: kernel<p2,ORDER,0,1>::a(ARGS); break;
    case p3: kernel<p3,ORDER,0,1>::a(ARGS); break;
    default: kernel<p0,ORDER,0,1>::a(ARGS); break;
    }
  else
    switch(KERN) {
    case p1: kernel<p1,ORDER,0,0>::a(ARGS); break;
    case p2: kernel<p2,ORDER,0,0>::a(ARGS); break;
    case p3: kernel<p3,ORDER,0,0>::a(ARGS); break;
    default: kernel<p0,ORDER,0,0>::a(ARGS); break;
    }
}
//------------------------------------------------------------------------------
void GravKernAll::approx(cell_iter const&A, leaf_iter const&B,
			 vect           &R, real      Rq) const
{
  give_coeffs(A);
  if(INDI_SOFT) 
    switch(KERN) {
    case p1: kernel<p1,ORDER,1,1>::a(ARGS); break;
    case p2: kernel<p2,ORDER,1,1>::a(ARGS); break;
    case p3: kernel<p3,ORDER,1,1>::a(ARGS); break;
    default: kernel<p0,ORDER,1,1>::a(ARGS); break;
    }
  else
    switch(KERN) {
    case p1: kernel<p1,ORDER,1,0>::a(ARGS); break;
    case p2: kernel<p2,ORDER,1,0>::a(ARGS); break;
    case p3: kernel<p3,ORDER,1,0>::a(ARGS); break;
    default: kernel<p0,ORDER,1,0>::a(ARGS); break;
    }
}
//------------------------------------------------------------------------------
void GravKern::approx(cell_iter const&A, cell_iter const&B,
		      vect           &R, real      Rq) const
{
  if(is_active(A)) give_coeffs(A);
  if(is_active(B)) give_coeffs(B);
  if(INDI_SOFT) 
    switch(KERN) {
    case p1: kernel<p1,ORDER,0,1>::a(ARGS); break;
    case p2: kernel<p2,ORDER,0,1>::a(ARGS); break;
    case p3: kernel<p3,ORDER,0,1>::a(ARGS); break;
    default: kernel<p0,ORDER,0,1>::a(ARGS); break;
    }
  else
    switch(KERN) {
    case p1: kernel<p1,ORDER,0,0>::a(ARGS); break;
    case p2: kernel<p2,ORDER,0,0>::a(ARGS); break;
    case p3: kernel<p3,ORDER,0,0>::a(ARGS); break;
    default: kernel<p0,ORDER,0,0>::a(ARGS); break;
    }
}
//------------------------------------------------------------------------------
void GravKernAll::approx(cell_iter const&A, cell_iter const&B,
			 vect           &R, real      Rq) const
{
  give_coeffs(A);
  give_coeffs(B);
  if(INDI_SOFT) 
    switch(KERN) {
    case p1: kernel<p1,ORDER,1,1>::a(ARGS); break;
    case p2: kernel<p2,ORDER,1,1>::a(ARGS); break;
    case p3: kernel<p3,ORDER,1,1>::a(ARGS); break;
    default: kernel<p0,ORDER,1,1>::a(ARGS); break;
    }
  else
    switch(KERN) {
    case p1: kernel<p1,ORDER,1,0>::a(ARGS); break;
    case p2: kernel<p2,ORDER,1,0>::a(ARGS); break;
    case p3: kernel<p3,ORDER,1,0>::a(ARGS); break;
    default: kernel<p0,ORDER,1,0>::a(ARGS); break;
    }
}
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_SSE_CODE                                                       
