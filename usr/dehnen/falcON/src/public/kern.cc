// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// kern.cc                                                                     |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2003                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <public/auxx.h>
#include <public/kern.h>
#include <public/coef.h>

using namespace nbdy;
using namespace nbdy::grav;
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class nbdy::TaylorSeries                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
inline void TaylorSeries::
shift_and_add(const grav_cell*const&c) {           // I: cell & its coeffs      
  if(hasCoeffs(c)) {                               // IF(cell has had iaction)  
    register vect dX = cofm(c) - X;                //   vector to shift by      
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
//                                                                            //
// class nbdy::grav_kern_base                                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
void grav_kern_base::eval_grav(cell_iter    const&C,
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
void grav_kern_base::eval_grav_all(cell_iter    const&C,
				   TaylorSeries const&T) const
{
  TaylorSeries G(T);                               // G = copy of T             
  G.shift_and_add(C);                              // shift G; G+=T_C           
  take_coeffs(C);                                  // free memory: C's coeffs   
  LoopLeafKids(cell_iter,C,l) {                    // LOOP C's leaf kids        
    l->normalize_grav();                           //   pot,acc/=mass           
//     // TEST
// #define __BODY 5
//     if(mybody(l)==__BODY) {
//       std::cerr<<" body #"<<__BODY<<":\n"
// 	       <<"    direct interations only: p="
// 	       <<pot(l)<<" a="<<acc(l)<<'\n';
//     }
// #undef __BODY
//     // TSET
    if(!is_empty(G)) G.extract_grav(l);            //   add pot,acc due to G    
  }                                                // END LOOP                  
  LoopCellKids(cell_iter,C,c)                      // LOOP C's cell kids        
    eval_grav_all(c,G);                            //   recursive call          
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
void grav_kern::single(leaf_iter const &A, leaf_iter const&B) const
{
  register vect R  = cofm(A)-cofm(B);
  register real Rq = norm(R),x,D0;
#ifdef falcON_INDI
  if(INDI_SOFT)
    switch(KERN) {
    case p1: { P1_I(mass(A)*mass(B)) } break;
    case p2: { P2_I(mass(A)*mass(B)) } break;
    case p3: { P3_I(mass(A)*mass(B)) } break;
    default: { P0_I(mass(A)*mass(B)) } break;
    }
  else
#endif
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
void grav_kern_all::single(leaf_iter const &A, leaf_iter const&B) const
{
  register vect R  = cofm(A)-cofm(B);
  register real Rq = norm(R),x,D0;
#ifdef falcON_INDI
  if(INDI_SOFT)
    switch(KERN) {
    case p1: { P1_I(mass(A)*mass(B)) } break;
    case p2: { P2_I(mass(A)*mass(B)) } break;
    case p3: { P3_I(mass(A)*mass(B)) } break;
    default: { P0_I(mass(A)*mass(B)) } break;
    }
  else
#endif
    switch(KERN) {
    case p1: { P1(mass(A)*mass(B)) } break;
    case p2: { P2(mass(A)*mass(B)) } break;
    case p3: { P3(mass(A)*mass(B)) } break;
    default: { P0(mass(A)*mass(B)) } break;
    }
  A->pot()-=D0; A->acc()-=R;
  B->pot()-=D0; B->acc()+=R;
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// 2. Cell-Node interactions                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#ifdef falcON_SSE_CODE
#  include <proper/kern_sse.cc>
#else
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
  XX *= 3  * D1;     /* XX == T2 */		\
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
#define LOAD_G(RGHT)				\
  dR = X0 - cofm(RGHT);				\
  D1 = norm(dR) + EQ;				\
  D0 = M0 * mass(RGHT)
//------------------------------------------------------------------------------
#define LOAD_I(RGHT)				\
  dR = X0 - cofm(RGHT);				\
  EQ = square(E0+eph(RGHT));			\
  D1 = norm(dR) + EQ;				\
  D0 = M0 * mass(RGHT)
//------------------------------------------------------------------------------
#define PUT_LEFT(RGHT)				\
  dR *= D1;					\
  P0 -= D0;					\
  F0 -= dR;
//------------------------------------------------------------------------------
#define PUT_RGHT(RGHT)				\
  dR          *= D1;				\
  RGHT->pot() -= D0;				\
  RGHT->acc() += dR;
//------------------------------------------------------------------------------
#define PUT_BOTH(RGHT)				\
  dR          *= D1;				\
  P0          -= D0;				\
  F0          -= dR;				\
  RGHT->pot() -= D0;				\
  RGHT->acc() += dR;
//------------------------------------------------------------------------------
#define PUT_SOME(RGHT)				\
  dR   *= D1;					\
  P0   -= D0;					\
  F0   -= dR;					\
  if(is_active(RGHT)) {				\
    RGHT->pot() -= D0;				\
    RGHT->acc() += dR;				\
  }
//------------------------------------------------------------------------------
#define GRAV_ALL(LOAD,DSINGL,PUT)		\
for(register leaf_iter B=B0; B!=BN; ++B) {	\
  LOAD(B);					\
  DSINGL					\
  PUT(B)					\
} 
//------------------------------------------------------------------------------
#define GRAV_FEW(LOAD,DSINGL)			\
for(register leaf_iter B=B0; B!=BN; ++B)	\
  if(is_active(B)) {				\
    LOAD(B);					\
    DSINGL					\
    PUT_RGHT(B)					\
  } 
//------------------------------------------------------------------------------
#define ARGS					\
  leaf_iter const&A,				\
  leaf_iter const&B0,				\
  leaf_iter const&BN,				\
  real&EQ, real&HQ, real&QQ 
//------------------------------------------------------------------------------
#define START_G					\
  const    real      M0=mass(A);		\
  const    vect      X0=cofm(A);		\
  register vect      dR;			\
  register real      D0,D1;
//------------------------------------------------------------------------------
#define START_I					\
  const    real      E0=eph(A);			\
  const    real      M0=mass(A);		\
  const    vect      X0=cofm(A);		\
  register vect      dR;			\
  register real      D0,D1;
//==============================================================================
// now defining auxiliary inline functions for the computation of  N            
// interactions. There are the following 10 cases:                              
// - each for the cases YA, YS, YN, NA, NS                                      
// - each for global and individual softening                                   
//==============================================================================
namespace {
  using namespace nbdy; using nbdy::uint; using namespace nbdy::grav;
  //////////////////////////////////////////////////////////////////////////////
#define DIRECT(START,LOAD,DSINGL)				\
    static void many_YA(ARGS) {					\
      START; register real P0=zero; register vect F0=zero;	\
      GRAV_ALL(LOAD,DSINGL,PUT_BOTH)				\
      A->pot()+=P0;  A->acc()+=F0;				\
    }								\
    static void many_YS(ARGS) {					\
      START; register real P0=zero; register vect F0=zero;	\
      GRAV_ALL(LOAD,DSINGL,PUT_SOME)				\
      A->pot()+=P0; A->acc()+=F0;				\
    }								\
    static void many_YN(ARGS) {					\
      START; register real P0=zero; register vect F0=zero;	\
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
  template<> struct __direct<p0,0> {
    DIRECT(START_G,LOAD_G,DSINGL_P0_G);
  };
#ifdef falcON_INDI
  template<> struct __direct<p0,1> {
    DIRECT(START_I,LOAD_I,DSINGL_P0_I);
  };
#endif
  //----------------------------------------------------------------------------
  template<> struct __direct<p1,0> {
    DIRECT(START_G,LOAD_G,DSINGL_P1_G);
  };
#ifdef falcON_INDI
  template<> struct __direct<p1,1> {
    DIRECT(START_I,LOAD_I,DSINGL_P1_I);
  };
#endif
  //----------------------------------------------------------------------------
  template<> struct __direct<p2,0> {
    DIRECT(START_G,LOAD_G,DSINGL_P2_G);
  };
#ifdef falcON_INDI
  template<> struct __direct<p2,1> {
    DIRECT(START_I,LOAD_I,DSINGL_P2_I);
  };
#endif
  //----------------------------------------------------------------------------
  template<> struct __direct<p3,0> {
    DIRECT(START_G,LOAD_G,DSINGL_P3_G);
  };
#ifdef falcON_INDI
  template<> struct __direct<p3,1> {
    DIRECT(START_I,LOAD_I,DSINGL_P3_I);
  };
#endif
#undef LOAD_G
#undef LOAD_I
#undef START_G
#undef START_I
#undef DIRECT
  //////////////////////////////////////////////////////////////////////////////
  template<bool I> struct Direct {
    static void many_YA(kern_type const&KERN, ARGS) {
      switch(KERN) {
      case p1: __direct<p1,I>::many_YA(A,B0,BN,EQ,HQ,QQ); break;
      case p2: __direct<p2,I>::many_YA(A,B0,BN,EQ,HQ,QQ); break;
      case p3: __direct<p3,I>::many_YA(A,B0,BN,EQ,HQ,QQ); break;
      default: __direct<p0,I>::many_YA(A,B0,BN,EQ,HQ,QQ); break;
      } }
    static void many_YS(kern_type const&KERN, ARGS) {
      switch(KERN) {
      case p1: __direct<p1,I>::many_YS(A,B0,BN,EQ,HQ,QQ); break;
      case p2: __direct<p2,I>::many_YS(A,B0,BN,EQ,HQ,QQ); break;
      case p3: __direct<p3,I>::many_YS(A,B0,BN,EQ,HQ,QQ); break;
      default: __direct<p0,I>::many_YS(A,B0,BN,EQ,HQ,QQ); break;
      } }
    static void many_YN(kern_type const&KERN, ARGS) {
      switch(KERN) {
      case p1: __direct<p1,I>::many_YN(A,B0,BN,EQ,HQ,QQ); break;
      case p2: __direct<p2,I>::many_YN(A,B0,BN,EQ,HQ,QQ); break;
      case p3: __direct<p3,I>::many_YN(A,B0,BN,EQ,HQ,QQ); break;
      default: __direct<p0,I>::many_YN(A,B0,BN,EQ,HQ,QQ); break;
      } }
    static void many_NA(kern_type const&KERN, ARGS) {
      switch(KERN) {
      case p1: __direct<p1,I>::many_NA(A,B0,BN,EQ,HQ,QQ); break;
      case p2: __direct<p2,I>::many_NA(A,B0,BN,EQ,HQ,QQ); break;
      case p3: __direct<p3,I>::many_NA(A,B0,BN,EQ,HQ,QQ); break;
      default: __direct<p0,I>::many_NA(A,B0,BN,EQ,HQ,QQ); break;
      } }
    static void many_NS(kern_type const&KERN, ARGS) {
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
//==============================================================================
// we now can define the cell-leaf and cell-self interaction via direct sums    
//==============================================================================
#define ARGS KERN,B,A.begin_leafs(),A.end_leaf_desc(),EQ,HQ,QQ
void grav_kern::direct(cell_iter const&A, leaf_iter const&B) const
{
#ifdef falcON_INDI
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
#endif
    if(is_active(B)) {
      if     (al_active(A)) Direct<0>::many_YA(ARGS);
      else if(is_active(A)) Direct<0>::many_YS(ARGS);
      else                  Direct<0>::many_YN(ARGS);
    } else {
      if     (al_active(A)) Direct<0>::many_NA(ARGS);
      else if(is_active(A)) Direct<0>::many_NS(ARGS);
    }
}
//------------------------------------------------------------------------------
void grav_kern_all::direct(cell_iter const&A, leaf_iter const&B) const
{
#ifdef falcON_INDI
  if(INDI_SOFT)
    Direct<1>::many_YA(ARGS);
  else
#endif
    Direct<0>::many_YA(ARGS);
}
#undef ARGS
//------------------------------------------------------------------------------
#define ARGS KERN,A,A+1,A+1+Nk,EQ,HQ,QQ
void grav_kern::direct(cell_iter const&C) const
{
  const    uint      N1 = number(C)-1;
  register uint      k, Nk;
  register leaf_iter A = C.begin_leafs();
#ifdef falcON_INDI
  if(INDI_SOFT)
    if(al_active(C))
      for(k=0, Nk=N1; Nk!=0; --Nk, ++k, ++A) Direct<1>::many_YA(ARGS);
    else
      for(k=0, Nk=N1; Nk!=0; --Nk, ++k, ++A)
	if(is_active(A))                     Direct<1>::many_YS(ARGS);
	else                                 Direct<1>::many_NS(ARGS);
  else
#endif
    if(al_active(C))
      for(k=0, Nk=N1; Nk!=0; --Nk, ++k, ++A) Direct<0>::many_YA(ARGS);
    else
      for(k=0, Nk=N1; Nk!=0; --Nk, ++k, ++A)
	if(is_active(A))                     Direct<0>::many_YS(ARGS);
	else                                 Direct<0>::many_NS(ARGS);
}
//------------------------------------------------------------------------------
void grav_kern_all::direct(cell_iter const&C) const
{
  const    uint      N1 = number(C)-1;
  register uint      k, Nk;
  register leaf_iter A = C.begin_leafs();
#ifdef falcON_INDI
  if(INDI_SOFT)
    for(k=0, Nk=N1; Nk!=0; --Nk, ++k, ++A) Direct<1>::many_YA(ARGS);
  else
#endif
    for(k=0, Nk=N1; Nk!=0; --Nk, ++k, ++A) Direct<0>::many_YA(ARGS);
}
#undef ARGS
//==============================================================================
// we now define non-inline functions for the computation of many direct        
// interactions between NA and NB leafs                                         
// there are 8 cases, depending on whether all, some, or none of either A or    
// B are active (case none,none is trivial).                                    
//                                                                              
// these functions are called by grav_kern::direct(cell,cell), which is inline  
// in kern.h, or by grav_kern::flush_scc() below.                               
//==============================================================================
void grav_kern_base::many_AA(leaf_iter const&A0, uint const&NA,
			     leaf_iter const&B0, uint const&NB) const
{
  const leaf_iter AN=A0+NA, BN=B0+NB;
#ifdef falcON_INDI
  if(INDI_SOFT)
    for(leaf_iter A=A0; A!=AN; ++A) Direct<1>::many_YA(KERN,A,B0,BN,EQ,HQ,QQ);
  else
#endif
    for(leaf_iter A=A0; A!=AN; ++A) Direct<0>::many_YA(KERN,A,B0,BN,EQ,HQ,QQ);
}
//------------------------------------------------------------------------------
void grav_kern_base::many_AS(leaf_iter const&A0, uint const&NA,
			     leaf_iter const&B0, uint const&NB) const
{
  const    leaf_iter AN=A0+NA, BN=B0+NB;
#ifdef falcON_INDI
  if(INDI_SOFT)
    for(leaf_iter A=A0; A!=AN; ++A) Direct<1>::many_YS(KERN,A,B0,BN,EQ,HQ,QQ);
  else
#endif
    for(leaf_iter A=A0; A!=AN; ++A) Direct<0>::many_YS(KERN,A,B0,BN,EQ,HQ,QQ);
}
//------------------------------------------------------------------------------
void grav_kern_base::many_AN(leaf_iter const&A0, uint const&NA,
			     leaf_iter const&B0, uint const&NB) const
{
  const leaf_iter AN=A0+NA, BN=B0+NB;
#ifdef falcON_INDI
  if(INDI_SOFT)
    for(leaf_iter A=A0; A!=AN; ++A) Direct<1>::many_YN(KERN,A,B0,BN,EQ,HQ,QQ);
  else
#endif
    for(leaf_iter A=A0; A!=AN; ++A) Direct<0>::many_YN(KERN,A,B0,BN,EQ,HQ,QQ);
}
//------------------------------------------------------------------------------
void grav_kern_base::many_SA(leaf_iter const&A0, uint const&NA,
			     leaf_iter const&B0, uint const&NB) const
{
  const leaf_iter AN=A0+NA, BN=B0+NB;
#ifdef falcON_INDI
  if(INDI_SOFT) {
    for(leaf_iter A=A0; A!=AN; ++A)
      if(is_active(A)) Direct<1>::many_YA(KERN,A,B0,BN,EQ,HQ,QQ);
      else             Direct<1>::many_NA(KERN,A,B0,BN,EQ,HQ,QQ);
  } else
#endif
    for(leaf_iter A=A0; A!=AN; ++A)
      if(is_active(A)) Direct<0>::many_YA(KERN,A,B0,BN,EQ,HQ,QQ);
      else             Direct<0>::many_NA(KERN,A,B0,BN,EQ,HQ,QQ);
}
//------------------------------------------------------------------------------
void grav_kern_base::many_SS(leaf_iter const&A0, uint const&NA,
			     leaf_iter const&B0, uint const&NB) const
{
  const leaf_iter AN=A0+NA, BN=B0+NB;
#ifdef falcON_INDI
  if(INDI_SOFT) {
    for(leaf_iter A=A0; A!=AN; ++A)
      if(is_active(A)) Direct<1>::many_YS(KERN,A,B0,BN,EQ,HQ,QQ);
      else             Direct<1>::many_NS(KERN,A,B0,BN,EQ,HQ,QQ);
  } else
#endif
    for(leaf_iter A=A0; A!=AN; ++A) 
      if(is_active(A)) Direct<0>::many_YS(KERN,A,B0,BN,EQ,HQ,QQ);
      else             Direct<0>::many_NS(KERN,A,B0,BN,EQ,HQ,QQ);
}
//------------------------------------------------------------------------------
void grav_kern_base::many_SN(leaf_iter const&A0, uint const&NA,
			     leaf_iter const&B0, uint const&NB) const
{
  const leaf_iter AN=A0+NA, BN=B0+NB;
#ifdef falcON_INDI
  if(INDI_SOFT) {
    for(leaf_iter A=A0; A!=AN; ++A)
      if(is_active(A)) Direct<1>::many_YN(KERN,A,B0,BN,EQ,HQ,QQ);
  } else
#endif
    for(leaf_iter A=A0; A!=AN; ++A) 
      if(is_active(A)) Direct<0>::many_YN(KERN,A,B0,BN,EQ,HQ,QQ);
}
//------------------------------------------------------------------------------
void grav_kern_base::many_NA(leaf_iter const&A0, uint const&NA,
			     leaf_iter const&B0, uint const&NB) const
{
  const leaf_iter AN=A0+NA, BN=B0+NB;
#ifdef falcON_INDI
  if(INDI_SOFT)
    for(leaf_iter A=A0; A!=AN; ++A) Direct<1>::many_NA(KERN,A,B0,BN,EQ,HQ,QQ);
  else
#endif
    for(leaf_iter A=A0; A!=AN; ++A) Direct<0>::many_NA(KERN,A,B0,BN,EQ,HQ,QQ);
}
//------------------------------------------------------------------------------
void grav_kern_base::many_NS(leaf_iter const&A0, uint const&NA,
			     leaf_iter const&B0, uint const&NB) const
{
  const leaf_iter AN=A0+NA, BN=B0+NB;
#ifdef falcON_INDI
  if(INDI_SOFT)
    for(leaf_iter A=A0; A!=AN; ++A) Direct<1>::many_NS(KERN,A,B0,BN,EQ,HQ,QQ);
  else
#endif
    for(leaf_iter A=A0; A!=AN; ++A) Direct<0>::many_NS(KERN,A,B0,BN,EQ,HQ,QQ);
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// 2.2 approximate gravity                                                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace {
  using namespace nbdy; using nbdy::uint; using namespace nbdy::grav;
  //////////////////////////////////////////////////////////////////////////////
#define LOAD_G					\
  register real D[ND];				\
  register real XX=one/(Rq+EQ);			\
  D[0] = mass(A)*mass(B);
  //////////////////////////////////////////////////////////////////////////////
#define LOAD_I					\
  register real D[ND];				\
  EQ   = square(eph(A)+eph(B));			\
  register real XX=one/(Rq+EQ);			\
  D[0] = mass(A)*mass(B);			\
  s(EQ,HQ,QQ); 
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
  // class kernel<kern_type, order, all, indi_soft>                           //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<kern_type>          struct __setE;
  template<kern_type,int>      struct __block;
#define sv   static void
  //////////////////////////////////////////////////////////////////////////////
  // kern_type = p0                                                             
  //////////////////////////////////////////////////////////////////////////////
  template<> struct __setE<p0> {
    sv s(real const&,real&,real&) {
    } };
  //----------------------------------------------------------------------------
  template<> struct __block<p0,1> : public __setE<p0> {
    enum { ND=2 };
    sv b(real&X, real*D, real const&EQ, real const&HQ, real const&QQ) {
      D[0] *= sqrt(X);
      D[1]  = X * D[0];
    } };
  template<int K> struct __block<p0,K> : public __setE<p0> {
    enum { ND=K+1, F=K+K-1 };
    sv b(real&X, real*D, real const&EQ, real const&HQ, real const&QQ) {
      __block<p0,K-1>::b(X,D,EQ,HQ,QQ);
      D[K] = F * X * D[K-1];
    } };
  //////////////////////////////////////////////////////////////////////////////
  // kern_type = p1                                                             
  //////////////////////////////////////////////////////////////////////////////
  template<> struct __setE<p1> {
    sv s(real const&EQ, real&HQ, real&QQ) {
      HQ = half * EQ;
    } };
  //----------------------------------------------------------------------------
  template<> struct __block<p1,1> : public __setE<p1> {
    enum { ND=3 };
    sv b(real&X, real*D, real const&EQ, real const&HQ, real const&QQ) {
      D[0] *= sqrt(X);
      D[1]  =     X * D[0];
      D[2]  = 3 * X * D[1];
      D[0] += HQ*D[1];
      D[1] += HQ*D[2];
    } };
  template<int K> struct __block<p1,K> : public __setE<p1> {
    enum { ND=K+2, F=K+K+1 };
    sv b(real&X, real*D, real const&EQ, real const&HQ, real const&QQ) {
      __block<p1,K-1>::b(X,D,EQ,HQ,QQ);
      D[K+1] = F * X * D[K];
      D[K]  += HQ  * D[K+1];
    } };
  //////////////////////////////////////////////////////////////////////////////
  // kern_type = p2                                                             
  //////////////////////////////////////////////////////////////////////////////
  template<> struct __setE<p2> {
    sv s(real const&EQ, real&HQ, real&QQ) {
      HQ = half * EQ;
    } };
  //----------------------------------------------------------------------------
  template<> struct __block<p2,1> : public __setE<p2> {
    enum { ND=4 };
    sv b(real&X, real*D, real const&EQ, real const&HQ, real const&QQ) {
      D[0] *= sqrt(X);
      D[1]  =     X * D[0];
      D[2]  = 3 * X * D[1];
      D[3]  = 5 * X * D[2];
      D[0] += HQ*(D[1]+HQ*D[2]);
      D[1] += HQ*(D[2]+HQ*D[3]);
    } };
  template<int K> struct __block<p2,K> : public __setE<p2> {
    enum { ND=K+3, F=K+K+3 };
    sv b(real&X, real*D, real const&EQ, real const&HQ, real const&QQ) {
      __block<p2,K-1>::b(X,D,EQ,HQ,QQ);
      D[K+2] = F * X * D[K+1];
      D[K]  += HQ*(D[K+1]+HQ*D[K+2]);
    } };
  //////////////////////////////////////////////////////////////////////////////
  // kern_type = p3                                                             
  //////////////////////////////////////////////////////////////////////////////
  template<> struct __setE<p3> {
    sv s(real const&EQ, real&HQ, real&QQ) {
      HQ = half * EQ;
      QQ = half * QQ;
    } };
  //----------------------------------------------------------------------------
  template<> struct __block<p3,1> : public __setE<p3> {
    enum { ND=5 };
    sv b(real&X, real*D, real const&EQ, real const&HQ, real const&QQ) {
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
    sv b(real&X, real*D, real const&EQ, real const&HQ, real const&QQ) {
      __block<p3,K-1>::b(X,D,EQ,HQ,QQ);
      D[K+3] = F * X * D[K+2];
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
    enum { ND=K+1 };
    sv a(ARGS_B) { LOAD_G b(XX,D,EQ,HQ,QQ); CellLeaf(A,B,D,0,R); }
    sv a(ARGS_C) { LOAD_G b(XX,D,EQ,HQ,QQ); CellCell(A,B,D,0,R); }
  };
  template<kern_type P, int K> struct kernel<P,K,0,1> : private __block<P,K> {
    enum { ND=K+1 };
    sv a(ARGS_B) { LOAD_I b(XX,D,EQ,HQ,QQ); CellLeaf(A,B,D,0,R); }
    sv a(ARGS_C) { LOAD_I b(XX,D,EQ,HQ,QQ); CellCell(A,B,D,0,R); }
  };
  template<kern_type P, int K> struct kernel<P,K,1,0> : private __block<P,K> {
   enum { ND=K+1 };
    sv a(ARGS_B) { LOAD_G b(XX,D,EQ,HQ,QQ); CellLeafAll(A,B,D,0,R); }
    sv a(ARGS_C) { LOAD_G b(XX,D,EQ,HQ,QQ); CellCellAll(A,B,D,0,R); }
  };
  template<kern_type P, int K> struct kernel<P,K,1,1> : private __block<P,K> {
    enum { ND=K+1 };
    sv a(ARGS_B) { LOAD_I b(XX,D,EQ,HQ,QQ); CellLeafAll(A,B,D,0,R); }
    sv a(ARGS_C) { LOAD_I b(XX,D,EQ,HQ,QQ); CellCellAll(A,B,D,0,R); }
  };
#undef sv
  //////////////////////////////////////////////////////////////////////////////
}                                                  // END: unnamed namespace    
////////////////////////////////////////////////////////////////////////////////
#define ARGS A,B,R,Rq,EQ,HQ,QQ
void grav_kern::approx(cell_iter const&A, leaf_iter const&B,
		       vect           &R, real      const&Rq) const
{
  if(is_active(A)) give_coeffs(A);
#ifdef falcON_INDI
  if(INDI_SOFT) 
    switch(KERN) {
    case p1: kernel<p1,ORDER,0,1>::a(ARGS); break;
    case p2: kernel<p2,ORDER,0,1>::a(ARGS); break;
    case p3: kernel<p3,ORDER,0,1>::a(ARGS); break;
    default: kernel<p0,ORDER,0,1>::a(ARGS); break;
    }
  else
#endif
    switch(KERN) {
    case p1: kernel<p1,ORDER,0,0>::a(ARGS); break;
    case p2: kernel<p2,ORDER,0,0>::a(ARGS); break;
    case p3: kernel<p3,ORDER,0,0>::a(ARGS); break;
    default: kernel<p0,ORDER,0,0>::a(ARGS); break;
    }
}
//------------------------------------------------------------------------------
void grav_kern_all::approx(cell_iter const&A, leaf_iter const&B,
			   vect           &R, real      const&Rq) const
{
  give_coeffs(A);
#ifdef falcON_INDI
  if(INDI_SOFT) 
    switch(KERN) {
    case p1: kernel<p1,ORDER,1,1>::a(ARGS); break;
    case p2: kernel<p2,ORDER,1,1>::a(ARGS); break;
    case p3: kernel<p3,ORDER,1,1>::a(ARGS); break;
    default: kernel<p0,ORDER,1,1>::a(ARGS); break;
    }
  else
#endif
    switch(KERN) {
    case p1: kernel<p1,ORDER,1,0>::a(ARGS); break;
    case p2: kernel<p2,ORDER,1,0>::a(ARGS); break;
    case p3: kernel<p3,ORDER,1,0>::a(ARGS); break;
    default: kernel<p0,ORDER,1,0>::a(ARGS); break;
    }
}
//------------------------------------------------------------------------------
void grav_kern::approx(cell_iter const&A, cell_iter const&B,
		       vect           &R, real      const&Rq) const
{
  if(is_active(A)) give_coeffs(A);
  if(is_active(B)) give_coeffs(B);
#ifdef falcON_INDI
  if(INDI_SOFT) 
    switch(KERN) {
    case p1: kernel<p1,ORDER,0,1>::a(ARGS); break;
    case p2: kernel<p2,ORDER,0,1>::a(ARGS); break;
    case p3: kernel<p3,ORDER,0,1>::a(ARGS); break;
    default: kernel<p0,ORDER,0,1>::a(ARGS); break;
    }
  else
#endif
    switch(KERN) {
    case p1: kernel<p1,ORDER,0,0>::a(ARGS); break;
    case p2: kernel<p2,ORDER,0,0>::a(ARGS); break;
    case p3: kernel<p3,ORDER,0,0>::a(ARGS); break;
    default: kernel<p0,ORDER,0,0>::a(ARGS); break;
    }
}
//------------------------------------------------------------------------------
void grav_kern_all::approx(cell_iter const&A, cell_iter const&B,
			   vect           &R, real      const&Rq) const
{
  give_coeffs(A);
  give_coeffs(B);
#ifdef falcON_INDI
  if(INDI_SOFT) 
    switch(KERN) {
    case p1: kernel<p1,ORDER,1,1>::a(ARGS); break;
    case p2: kernel<p2,ORDER,1,1>::a(ARGS); break;
    case p3: kernel<p3,ORDER,1,1>::a(ARGS); break;
    default: kernel<p0,ORDER,1,1>::a(ARGS); break;
    }
  else
#endif
    switch(KERN) {
    case p1: kernel<p1,ORDER,1,0>::a(ARGS); break;
    case p2: kernel<p2,ORDER,1,0>::a(ARGS); break;
    case p3: kernel<p3,ORDER,1,0>::a(ARGS); break;
    default: kernel<p0,ORDER,1,0>::a(ARGS); break;
    }
}
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_SSE_CODE                                                       
