// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// tn2D.cc                                                                     |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2003-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_tn2D_cc
#define falcON_included_tn2D_cc

#ifndef falcON_included_tn2D_h
#  include <public/tn2D.h>
#endif

namespace meta {
  //////////////////////////////////////////////////////////////////////////////
#define tm    template
#define tKX   template<int K,typename X>
#define tSY   tKX inline symt2D<K,X>
#define tS2   tX inline symt2D<2,X>
#define V     tupel<2,X>
#define tQ    template<int Q>
#define sv    static void
#define si    static int
#define sci   static const int
#define scb   static const bool
#define sX    static X
#define cX    const  X
#define tX    template<typename X>
#define cA    const  A
#define cB    const  B
#define cC    const  C
#define cD    const  D
#define tA    template<typename A>
#define tAB   template<typename A, typename B>
#define tABC  template<typename A, typename B, typename C>
#define tABCD template<typename A, typename B, typename C, typename D>
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // macros to be used for minimum and maximum of const integer expressions   //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
#undef  __min
#undef  __max
#ifdef __GNUC__
#  define __min(A,B) (A)<?(B)
#  define __max(A,B) (A)>?(B)
#else
#  define __min(A,B) (A)<(B)? (A):(B)
#  define __max(A,B) (A)>(B)? (A):(B)
#endif
////////////////////////////////////////////////////////////////////////////////
namespace meta2D {
  using namespace meta;
  using namespace nbdy;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // struct ONE2D<N,X>                                                        //
  //                                                                          //
  // generic constants, types and methods                                     //
  //                                                                          //
  // ND:  # independent elements in symmetric tensor of order N               //
  // CD:  # independent elements in symmetric tensors of orders 0 to N        //
  // I0:  index of first element of Nth tensor of set with orders 0 to N      //
  // PD:  # independent elements in symmetric tenrors of orders 2 to N        //
  // J0:  index of first element of Nth tensor of set with orders 2 to N      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  tm<int N, typename X=real> struct ONE2D : ONE<N> {
    enum {
      ND = N+1,                                    // N+1                       
      CD = ((N+1)*(N+2))/2,                        // (N+1)(N+2)/2              
      I0 = (N*(N+1))/2,                            // N(N+1)/2                  
      PD = CD-3,                                   // (N+1)(N+2)/2-3            
      J0 = I0-3                                    // N(N+1)/2-3                
    };
    typedef symt2D<N,X> Tensor;                    // sym 2D tensor of order N  
    static Tensor      &tens(      X*a) { return *static_cast<Tensor*>
					    (static_cast<void*>(a+I0)); }
    static Tensor const&tens(const X*a) { return *static_cast<const Tensor*>
					    (static_cast<const void*>(a+I0)); }
    static Tensor      &pole(      X*a) { return *static_cast<Tensor*>
					    (static_cast<void*>(a+J0)); }
    static Tensor const&pole(const X*a) { return *static_cast<const Tensor*>
					    (static_cast<const void*>(a+J0)); }
  };
  tm<typename X> struct ONE2D<1,X> : ONE<1> {
    enum { ND=2, CD=3, I0=1 };
    typedef tupel<2,X> Tensor;
    static Tensor      &tens(      X*a) { return *static_cast<Tensor*>
					    (static_cast<void*>(a+I0)); }
    static Tensor const&tens(const X*a) { return *static_cast<const Tensor*>
					    (static_cast<const void*>(a+I0)); }
  };
  tm<typename X> struct ONE2D<0,X> : ONE<0> {
    enum { ND=1, CD=1, I0=0 };
    typedef X Tensor;
    static Tensor      &tens(      X*a) { return *a; }
    static Tensor const&tens(const X*a) { return *a; }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // struct LoopL<Perform,Z,K>                                                //
  //                                                                          //
  // loop l  and  perform some templated action                               //
  //                                                                          //
  // template arguments of P: Z(anything),K,L                                 //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
#define T3i template<int,int,int> class
  // struct l for looping L=0...K, M=0...L                                      
  tm<T3i P, int Z, int K, int L> struct l {                     // K,L          
    tA    sv loop(A&a) { 
      P<  Z,K,L  >::job (a);
      l<P,Z,K,L+1>::loop(a); }
    tAB   sv loop(A&a,cB&b) { 
      P<  Z,K,L  >::job (a,b);
      l<P,Z,K,L+1>::loop(a,b); }
    tABC  sv loop(A&a,cB&b,cC&c) { 
      P<  Z,K,L  >::job (a,b,c);
      l<P,Z,K,L+1>::loop(a,b,c); }
    tABCD sv loop(A&a,cB&b,cC&c,cD&d) { 
      P<  Z,K,L  >::job (a,b,c,d);
      l<P,Z,K,L+1>::loop(a,b,c,d); }
  };
  //----------------------------------------------------------------------------
  tm<T3i P, int Z, int K> struct l<P,Z,K,K> {                   // L=K          
    tA    sv loop(A&a)                { P<Z,K,K>::job(a); }
    tAB   sv loop(A&a,cB&b)           { P<Z,K,K>::job(a,b); }
    tABC  sv loop(A&a,cB&b,cC&c)      { P<Z,K,K>::job(a,b,c); }
    tABCD sv loop(A&a,cB&b,cC&c,cD&d) { P<Z,K,K>::job(a,b,c,d); }
  };
  //----------------------------------------------------------------------------
  tm<T3i P, int K, int Z> struct LoopL {
    tA    sv loop(A&a)                { l<P,Z,K,0>::loop(a); }
    tAB   sv loop(A&a,cB&b)           { l<P,Z,K,0>::loop(a,b); }
    tABC  sv loop(A&a,cB&b,cC&c)      { l<P,Z,K,0>::loop(a,b,c); }
    tABCD sv loop(A&a,cB&b,cC&c,cD&d) { l<P,Z,K,0>::loop(a,b,c,d); }
  };
#undef T3i
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // K fold outer product of vector                                           //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  tm<int K> struct vop {
    tX sv as(X*a,cX*u,cX*v) {
      taux<X,K-1>::v_ast(a,u,v[0]);
      a[K] = u[K-1] * v[1];
    }
    tX sv as(X*a,cX*u,cX*v,cX&x) {
      taux<X,K-1>::v_ast(a,u,x*v[0]);
      a[K] = u[K-1] * x*v[1];
    }
    tX sv ad(X*a,cX*u,cX*v) {
      taux<X,K-1>::v_adt(a,u,v[0]);
      a[K]+= u[K-1] * v[1];
    }
    tX sv ad(X*a,cX*u,cX*v,cX&x) {
      taux<X,K-1>::v_adt(a,u,x*v[0]);
      a[K]+= u[K-1] * x*v[1];
    }
    tX sv su(X*a,cX*u,cX*v) {
      taux<X,K-1>::v_sut(a,u,v[0]);
      a[K]-= u[K-1] * v[1];
      taux<X,K-1>::v_sut(a,u,x*v[0]);
      a[K]-= u[K-1] * x*v[1];
    }
    tX sv su(X*a,cX*u,cX*v,cX&x) {
      taux<X,K-1>::v_sut(a,u,x*v[0]);
      a[K]-= u[K-1] * x*v[1];
    }
  };
  //----------------------------------------------------------------------------
  tm<> struct vop<2> {
    tX sv as(X*a,cX*v) {
      taux<X,1>::v_ast(a,v,v[0]);
      a[2] = v[1] * v[1];
    }
    tX sv as(X*a,cX*v,cX&x) {
      taux<X,1>::v_ast(a,v,x*v[0]);
      a[2] = x*v[1] * v[1];
    }
    tX sv ad(X*a,cX*v) {
      taux<X,1>::v_adt(a,v,v[0]);
      a[2]+= v[1] * v[1];
    }
    tX sv ad(X*a,cX*v,cX&x) {
      taux<X,1>::v_adt(a,v,x*v[0]);
      a[2]+= x*v[1] * v[1];
    }
    tX sv su(X*a,cX*v) {
      taux<X,1>::v_sut(a,v,v[0]);
      a[2]-= v[1] * v[1];
    }
    tX sv su(X*a,cX*v,cX&x) {
      taux<X,1>::v_sut(a,v,x*v[0]);
      a[2]-= x*v[1] * v[1];
    }
  };
}                                                  // END namespace meta2D      
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // implementing                                                             //
  //                                                                          //
  // nbdy::symt2D<K,S>::ass_out_prd();                                        //
  // nbdy::symt2D<K,S>::add_out_prd();                                        //
  // nbdy::symt2D<K,S>::sub_out_prd();                                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  tSY& symt2D<K,X>::ass_out_prd(symt2D<K-1,X>const&u, V const&v) {
    meta2D::vop<K>::as(a,pX(u),pX(v)); return*this; }
  tSY& symt2D<K,X>::add_out_prd(symt2D<K-1,X>const&u, V const&v) {
    meta2D::vop<K>::ad(a,pX(u),pX(v)); return*this; }
  tSY& symt2D<K,X>::sub_out_prd(symt2D<K-1,X>const&u, V const&v) {
    meta2D::vop<K>::su(a,pX(u),pX(v)); return*this; }
  tSY& symt2D<K,X>::ass_out_prd(symt2D<K-1,X>const&u, V const&v, cX&s) {
    meta2D::vop<K>::as(a,pX(u),pX(v),s); return*this; }
  tSY& symt2D<K,X>::add_out_prd(symt2D<K-1,X>const&u, V const&v, cX&s) {
    meta2D::vop<K>::ad(a,pX(u),pX(v),s); return*this; }
  tSY& symt2D<K,X>::sub_out_prd(symt2D<K-1,X>const&u, V const&v, cX&s) {
    meta2D::vop<K>::su(a,pX(u),pX(v),s); return*this; }
  //////////////////////////////////////////////////////////////////////////////
  tS2& symt2D<2,X>::ass_out_prd(V const&v) {
    meta2D::vop<2>::as(a,pX(v)); return*this; }
  tS2& symt2D<2,X>::add_out_prd(V const&v) {
    meta2D::vop<2>::ad(a,pX(v)); return*this; }
  tS2& symt2D<2,X>::sub_out_prd(V const&v) {
    meta2D::vop<2>::su(a,pX(v)); return*this; }
  tS2& symt2D<2,X>::ass_out_prd(V const&v, cX&s) {
    meta2D::vop<2>::as(a,pX(v),s); return*this; }
  tS2& symt2D<2,X>::add_out_prd(V const&v, cX&s) {
    meta2D::vop<2>::ad(a,pX(v),s); return*this; }
  tS2& symt2D<2,X>::sub_out_prd(V const&v, cX&s) {
    meta2D::vop<2>::su(a,pX(v),s); return*this; }
}                                                  // END namespace nbdy        
namespace meta2D {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // symmetric outer product                                                  //
  //                                                                          //
  // coded as follows:                                                        //
  //                                                                          //
  // 1. loop L of product using LoopL<>                                       //
  // 2. at each L    loop L1 = max(0,L-K2) ... min(L,K1)                      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  //----------------------------------------------------------------------------
  //                                                                            
  // struct syp_p<> computes single contribution to symmetric outer product     
  //                                                                            
  // static const int P:  factor in symmetric outer product                     
  // p(a,b): return contribution to symmetric outer product                     
  //                                                                            
  //----------------------------------------------------------------------------
  tm<int K, int L, int K1, int L1>
  struct syp_p {
    enum { L2 = L-L1,  P = TWO<K-L,K1-L1>::B * TWO<L,L1>::B };
    tX sX p(cX*b,cX*c) { return times<P>( b[L1] * c[L2] ); }
  };
  //----------------------------------------------------------------------------
  //                                                                            
  // struct syp_l<>::l<>() for the loop over L1                                 
  //                                                                            
  //----------------------------------------------------------------------------
  tm<int K,int L,int K1,int L1,int L1M>
  struct syp_l {
    tX sX l(cX*b,cX*c) {
      return
	syp_p<K,L,K1,L1      >::p(b,c) + 
	syp_l<K,L,K1,L1+1,L1M>::l(b,c);
    }
  };
  //----------------------------------------------------------------------------
  tm<int K,int L,int K1,int L1>
  struct syp_l<K,L,K1,L1,L1> {
    tX sX l(cX*b,cX*c) { return syp_p<K,L,K1,L1>::p(b,c); }
  };
  //----------------------------------------------------------------------------
  //                                                                            
  // struct syp_as<>: trigger loop over L1 via syp_l<>::l()                     
  // struct syp_ad<>: trigger loop over L1 via syp_l<>::l()                     
  // struct syp_su<>: trigger loop over L1 via syp_l<>::l()                     
  //                                                                            
  //----------------------------------------------------------------------------
  tm<int K1,int K,int L> struct syp_as {
    tX sv job(X*a,cX*b,cX*c) {
      a[L] =    syp_l<K,L,K1, __max(0,L-K+K1), __min(L,K1)>::l(b,c); }
    tX sv job(X*a,cX*b,cX*c,cX&d) {
      a[L] =  d*syp_l<K,L,K1, __max(0,L-K+K1), __min(L,K1)>::l(b,c); }
  };
  //----------------------------------------------------------------------------
  tm<int K1,int K,int L> struct syp_ad {
    tX sv job(X*a,cX*b,cX*c) {
      a[L]+=    syp_l<K,L,K1, __max(0,L-K+K1), __min(L,K1)>::l(b,c); }
    tX sv job(X*a,cX*b,cX*c,cX&d) {
      a[L]+=  d*syp_l<K,L,K1, __max(0,L-K+K1), __min(L,K1)>::l(b,c); }
  };
  //----------------------------------------------------------------------------
  tm<int K1,int K,int L> struct syp_su {
    tX sv job(X*a,cX*b,cX*c) {
      a[L]-=    syp_l<K,L,K1, __max(0,L-K+K1), __min(L,K1)>::l(b,c); }
    tX sv job(X*a,cX*b,cX*c,cX&d) {
      a[L]-=  d*syp_l<K,L,K1, __max(0,L-K+K1), __min(L,K1)>::l(b,c); }
  };
  //----------------------------------------------------------------------------
  //                                                                            
  // struct syp<>: perform symmetric outer product                              
  //                                                                            
  //----------------------------------------------------------------------------
  tm<int K1,int K2> struct syp {
    tX sv as (X*a,cX*b,cX*c)      { LoopL<syp_as,K1+K2,K1>::loop(a,b,c); }
    tX sv ad (X*a,cX*b,cX*c)      { LoopL<syp_ad,K1+K2,K1>::loop(a,b,c); }
    tX sv su (X*a,cX*b,cX*c)      { LoopL<syp_su,K1+K2,K1>::loop(a,b,c); }
    tX sv ast(X*a,cX*b,cX*c,cX&d) { LoopL<syp_as,K1+K2,K1>::loop(a,b,c,d); }
    tX sv adt(X*a,cX*b,cX*c,cX&d) { LoopL<syp_ad,K1+K2,K1>::loop(a,b,c,d); }
    tX sv sut(X*a,cX*b,cX*c,cX&d) { LoopL<syp_su,K1+K2,K1>::loop(a,b,c,d); }
  };
}                                                  // END namespace meta2D      
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // implementing                                                             //
  //                                                                          //
  // nbdy::symt2D<K,X>::ass_sym_prd();                                        //
  // nbdy::symt2D<K,X>::add_sym_prd();                                        //
  // nbdy::symt2D<K,X>::sub_sym_prd();                                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  tKX tQ inline symt2D<K,X>& symt2D<K,X>::ass_sym_prd(symt2D<K-Q,X>const&p,
						      symt2D<Q  ,X>const&q) {
    meta2D::syp<K-Q,Q>::as(a,pX(p),pX(q)); return*this; }
  tKX tQ inline symt2D<K,X>& symt2D<K,X>::add_sym_prd(symt2D<K-Q,X>const&p,
						      symt2D<Q  ,X>const&q) {
    meta2D::syp<K-Q,Q>::ad(a,pX(p),pX(q)); return*this; }
  tKX tQ inline symt2D<K,X>& symt2D<K,X>::sub_sym_prd(symt2D<K-Q,X>const&p,
						      symt2D<Q  ,X>const&q) {
    meta2D::syp<K-Q,Q>::su(a,pX(p),pX(q)); return*this; }
  //----------------------------------------------------------------------------
  tKX tQ inline symt2D<K,X>& symt2D<K,X>::ass_sym_prd(symt2D<K-Q,X>const&p,
						      symt2D<Q  ,X>const&q,
						      X            const&s) {
    meta2D::syp<K-Q,Q>::as(a,pX(p),pX(q),s); return*this; }
  tKX tQ inline symt2D<K,X>& symt2D<K,X>::add_sym_prd(symt2D<K-Q,X>const&p,
						      symt2D<Q  ,X>const&q,
						      X            const&s) {
    meta2D::syp<K-Q,Q>::ad(a,pX(p),pX(q),s); return*this; }
  tKX tQ inline symt2D<K,X>& symt2D<K,X>::sub_sym_prd(symt2D<K-Q,X>const&p,
						      symt2D<Q  ,X>const&q,
						      X            const&s) {
    meta2D::syp<K-Q,Q>::su(a,pX(p),pX(q),s); return*this; }
  //----------------------------------------------------------------------------
  tSY& symt2D<K,X>::ass_sym_prd(symt2D<K-1,X>const&p,V const&q) {
    meta2D::syp<K-1,1>::as(a,pX(p),pX(q)); return*this; }
  tSY& symt2D<K,X>::add_sym_prd(symt2D<K-1,X>const&p,V const&q) {
    meta2D::syp<K-1,1>::ad(a,pX(p),pX(q)); return*this; }
  tSY& symt2D<K,X>::sub_sym_prd(symt2D<K-1,X>const&p,V const&q) {
    meta2D::syp<K-1,1>::su(a,pX(p),pX(q)); return*this; }
  tSY& symt2D<K,X>::ass_sym_prd(V const&q,symt2D<K-1,X>const&p) {
    meta2D::syp<K-1,1>::as(a,pX(p),pX(q)); return*this; }
  tSY& symt2D<K,X>::add_sym_prd(V const&q,symt2D<K-1,X>const&p) {
    meta2D::syp<K-1,1>::ad(a,pX(p),pX(q)); return*this; }
  tSY& symt2D<K,X>::sub_sym_prd(V const&q,symt2D<K-1,X>const&p) {
    meta2D::syp<K-1,1>::su(a,pX(p),pX(q)); return*this; }
  //----------------------------------------------------------------------------
  tSY& symt2D<K,X>::ass_sym_prd(symt2D<K-1,X>const&p,V const&q,cX&s) {
    meta2D::syp<K-1,1>::as(a,pX(p),pX(q),s); return*this; }
  tSY& symt2D<K,X>::add_sym_prd(symt2D<K-1,X>const&p,V const&q,cX&s) {
    meta2D::syp<K-1,1>::ad(a,pX(p),pX(q),s); return*this; }
  tSY& symt2D<K,X>::sub_sym_prd(symt2D<K-1,X>const&p,V const&q,cX&s) {
    meta2D::syp<K-1,1>::su(a,pX(p),pX(q),s); return*this; }
  tSY& symt2D<K,X>::ass_sym_prd(V const&q,symt2D<K-1,X>const&p,cX&s) {
    meta2D::syp<K-1,1>::as(a,pX(p),pX(q),s); return*this; }
  tSY& symt2D<K,X>::add_sym_prd(V const&q,symt2D<K-1,X>const&p,cX&s) {
    meta2D::syp<K-1,1>::ad(a,pX(p),pX(q),s); return*this; }
  tSY& symt2D<K,X>::sub_sym_prd(V const&q,symt2D<K-1,X>const&p,cX&s) {
    meta2D::syp<K-1,1>::su(a,pX(p),pX(q),s); return*this; }
  //////////////////////////////////////////////////////////////////////////////
  tS2& symt2D<2,X>::ass_sym_prd(V const&p,V const&q) {
    meta2D::syp<1,1>::as(a,pX(p),pX(q)); return*this; }
  tS2& symt2D<2,X>::add_sym_prd(V const&p,V const&q) {
    meta2D::syp<1,1>::ad(a,pX(p),pX(q)); return*this; }
  tS2& symt2D<2,X>::sub_sym_prd(V const&p,V const&q) {
    meta2D::syp<1,1>::su(a,pX(p),pX(q)); return*this; }
  //----------------------------------------------------------------------------
  tS2& symt2D<2,X>::ass_sym_prd(V const&p,V const&q,cX&s) {
    meta2D::syp<1,1>::as(a,pX(p),pX(q),s); return*this; }
  tS2& symt2D<2,X>::add_sym_prd(V const&p,V const&q,cX&s) {
    meta2D::syp<1,1>::ad(a,pX(p),pX(q),s); return*this; }
  tS2& symt2D<2,X>::sub_sym_prd(V const&p,V const&q,cX&s) {
    meta2D::syp<1,1>::su(a,pX(p),pX(q),s); return*this; }
}                                                  // END namespace nbdy        
namespace meta2D {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // symmetric outer product with delta_ij^(N)                                //
  //                                                                          //
  // coded as follows:                                                        //
  //                                                                          //
  // 1. loop L of product using LoopL<>                                       //
  // 2. at each L    loop l1 = max(0,(L-K2+1)/2) ... min(L,K1)/2              //
  //                                                                          //
  // where l1 = L1/2.                                                         //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  //----------------------------------------------------------------------------
  //                                                                            
  // N fold outer product of delta_ij                                           
  //                                                                            
  //----------------------------------------------------------------------------
  tm<int k, int l> struct D {
    enum {
      K = k+k,
      L = l+l,
      Z = ONE<k-l>::H * ONE<l>::H
    }; };
  //----------------------------------------------------------------------------
  //                                                                            
  // struct syd_p<K,L,k1,l1>                                                    
  //                                                                            
  // static const int P:  factor in symmetric outer product                     
  // p(a,b): return contribution to symmetric outer product                     
  //                                                                            
  //----------------------------------------------------------------------------
  tm<int K,int L, int k1, int l1> struct syd_p {
    enum {
      K1 = k1+k1,
      L1 = l1+l1,
      P  = D<k1,l1>::Z * TWO<K-L,K1-L1>::B * TWO<L,L1>::B,
      L2 = L-L1
    };
    tX sX p(cX*b) { return times<P>( b[L2] ); }
  };
  //----------------------------------------------------------------------------
  //                                                                            
  // struct syd_l<>::l<>() for the loop over l1                                 
  //                                                                            
  //----------------------------------------------------------------------------
  tm<int K,int L,int k1,int l1,int l1M> struct syd_l {
    tX sX l(cX*b) {
      return
	syd_p<K,L,k1,l1      >::p(b) +
	syd_m<K,L,k1,l1+1,l1M>::l(b);
    }
  };
  //----------------------------------------------------------------------------
  tm<int K,int L,int k1,int l1>
  struct syd_m<K,L,k1,l1,l1> {
    tX sX l(cX*b) { return syd_p<K,L,k1,l1>::p(b); }
  };
  //----------------------------------------------------------------------------
  //                                                                            
  // struct syd_as<>: trigger loop over L1 via syd_l<>::l()                     
  // struct syd_ad<>: trigger loop over L1 via syd_l<>::l()                     
  // struct syd_su<>: trigger loop over L1 via syd_l<>::l()                     
  //                                                                            
  //----------------------------------------------------------------------------
  tm<int k1,int K,int L> struct syd__ {
    enum {
      K2     = K-k1-k1,
      l1_lo  = __max(0,(L-K2+1)/2),
      l1_hi  = __min(L/2,k1) };
    scb loop = l1_lo <= l1_hi;
  };
  //----------------------------------------------------------------------------
  tm<int k1,int K,int L> struct syd_as : public syd__<k1,K,L> {
    tX sv job(X*a,cX*b) {
      a[L]= loop?   syd_l<K,L,k1,l1_lo,l1_hi>::l(b) : X(0) ;
    }
    tX sv job(X*a,cX*b,cX&d) {
      a[L]= loop? d*syd_l<K,L,k1,l1_lo,l1_hi>::l(b) : X(0) ;
    }
  };
  //----------------------------------------------------------------------------
  tm<int k1,int K,int L> struct syd_ad : public syd__<k1,K,L> {
    tX sv job(X*a,cX*b) {
      if(loop) a[L]+=   syd_l<K,L,k1,l1_lo,l1_hi>::l(b);
    }
    tX sv job(X*a,cX*b,cX&d) {
      if(loop) a[L]+= d*syd_l<K,L,k1,l1_lo,l1_hi>::l(b);
    }
  };
  //----------------------------------------------------------------------------
  tm<int k1,int K,int L> struct syd_su : public syd__<k1,K,L> {
    tX sv job(X*a,cX*b) {
      if(loop) a[L]-=   syd_l<K,L,k1,l1_lo,l1_hi>::l(b);
    }
    tX sv job(X*a,cX*b,cX&d) {
      if(loop) a[L]-= d*syd_l<K,L,k1,l1_lo,l1_hi>::l(b);
    }
  };
  //----------------------------------------------------------------------------
  //                                                                            
  // struct syd_as0: for equating to delta^(k)                                  
  // struct syd_ad0: for equating to delta^(k)                                  
  // struct syd_su0: for equating to delta^(k)                                  
  //                                                                            
  //----------------------------------------------------------------------------
  tm<int Q,int K,int L> struct syd_as0 : public D<K/2,L/2> {
    enum { k = K/2, l=L/2, m=M/2 };
    tX sv job(X*a)      { a[L] = (!(L&1))?    Z : X(0) ; }
    tX sv job(X*a,cX&d) { a[L] = (!(L&1))? d* Z : X(0) ; }
  };
  //----------------------------------------------------------------------------
  tm<int Q,int k,int l> struct syd_ad0 : public D<k,l> {
    tX sv job(X*a)      { a[L]+=    Z; }
    tX sv job(X*a,cX&d) { a[L]+= d* Z; }
  };
  //----------------------------------------------------------------------------
  tm<int Q,int k,int l> struct syd_su0 : public D<k,l> {
    tX sv job(X*a)      { a[L]-=    Z; }
    tX sv job(X*a,cX&d) { a[L]-= d* Z; }
  };
  //----------------------------------------------------------------------------
  //                                                                            
  // struct syd<>: perform symmetric outer product with delta^k                 
  //                                                                            
  //----------------------------------------------------------------------------
  tm<int k, int K2> struct syd {
    sci K=k+k+K2;
    tX sv asp(X*a,cX*b)      { LoopL<syd_as,K,k>::loop(a,b); }
    tX sv adp(X*a,cX*b)      { LoopL<syd_ad,K,k>::loop(a,b); }
    tX sv sup(X*a,cX*b)      { LoopL<syd_su,K,k>::loop(a,b); }
    tX sv asp(X*a,cX*b,cX&s) { LoopL<syd_as,K,k>::loop(a,b,s); }
    tX sv adp(X*a,cX*b,cX&s) { LoopL<syd_ad,K,k>::loop(a,b,s); }
    tX sv sup(X*a,cX*b,cX&s) { LoopL<syd_su,K,k>::loop(a,b,s); }
    tX sv as(symt2D<K,X>&A,symt2D<K2,X> const&B)      { asp((X*)A,(cX*)B); }
    tX sv ad(symt2D<K,X>&A,symt2D<K2,X> const&B)      { adp((X*)A,(cX*)B); }
    tX sv su(symt2D<K,X>&A,symt2D<K2,X> const&B)      { sup((X*)A,(cX*)B); }
    tX sv as(symt2D<K,X>&A,symt2D<K2,X> const&B,cX&s) { asp((X*)A,(cX*)B,s); }
    tX sv ad(symt2D<K,X>&A,symt2D<K2,X> const&B,cX&s) { adp((X*)A,(cX*)B,s); }
    tX sv su(symt2D<K,X>&A,symt2D<K2,X> const&B,cX&s) { sup((X*)A,(cX*)B,s); }
  };
  //----------------------------------------------------------------------------
  tm<int k> struct syd<k,1> {
    sci K=k+k+1;
    tX sv asp(X*a,cX*b)      { LoopL<syd_as,K,k>::loop(a,b); }
    tX sv adp(X*a,cX*b)      { LoopL<syd_ad,K,k>::loop(a,b); }
    tX sv sup(X*a,cX*b)      { LoopL<syd_su,K,k>::loop(a,b); }
    tX sv asp(X*a,cX*b,cX&s) { LoopL<syd_as,K,k>::loop(a,b,s); }
    tX sv adp(X*a,cX*b,cX&s) { LoopL<syd_ad,K,k>::loop(a,b,s); }
    tX sv sup(X*a,cX*b,cX&s) { LoopL<syd_su,K,k>::loop(a,b,s); }
    tX sv as(symt2D<K,X>&A,V const&B)      { asp((X*)A,(cX*)B); }
    tX sv ad(symt2D<K,X>&A,V const&B)      { adp((X*)A,(cX*)B); }
    tX sv su(symt2D<K,X>&A,V const&B)      { sup((X*)A,(cX*)B); }
    tX sv as(symt2D<K,X>&A,V const&B,cX&s) { asp((X*)A,(cX*)B,s); }
    tX sv ad(symt2D<K,X>&A,V const&B,cX&s) { adp((X*)A,(cX*)B,s); }
    tX sv su(symt2D<K,X>&A,V const&B,cX&s) { sup((X*)A,(cX*)B,s); }
  };
  //----------------------------------------------------------------------------
  tm<int k> struct syd<k,0> {
    sci K=k+k;
    tX sv asp(X*a)      { LoopL<syd_as0,K,0>::loop(a); }
    tX sv asp(X*a,cX&s) { LoopL<syd_as0,K,0>::loop(a,s); }
    tX sv adp(X*a)      { LoopL<syd_ad0,k,0>::loop(a); }
    tX sv adp(X*a,cX&s) { LoopL<syd_ad0,k,0>::loop(a,s); }
    tX sv sup(X*a)      { LoopL<syd_su0,k,0>::loop(a); }
    tX sv sup(X*a,cX&s) { LoopL<syd_su0,k,0>::loop(a,s); }
    tX sv as(symt2D<K,X>&A)      { asp((X*)A); }
    tX sv as(symt2D<K,X>&A,cX&s) { asp((X*)A,s); }
    tX sv ad(symt2D<K,X>&A)      { adp((X*)A); }
    tX sv ad(symt2D<K,X>&A,cX&s) { adp((X*)A,s); }
    tX sv su(symt2D<K,X>&A)      { sup((X*)A); }
    tX sv su(symt2D<K,X>&A,cX&s) { sup((X*)A,s); }
  };
}                                                  // END namespace meta2D      
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // implementing                                                             //
  //                                                                          //
  // nbdy::symt2D<K,X>::ass_del_prd();                                        //
  // nbdy::symt2D<K,X>::add_del_prd();                                        //
  // nbdy::symt2D<K,X>::sub_del_prd();                                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  tKX tQ inline symt2D<K,X>& symt2D<K,X>::ass_del_prd(symt2D<Q,X> const&p) {
    meta2D::syd<(K-Q)/2,Q>::as(*this,p); return*this; }
  tKX tQ inline symt2D<K,X>& symt2D<K,X>::add_del_prd(symt2D<Q,X> const&p) {
    meta2D::syd<(K-Q)/2,Q>::ad(*this,p); return*this; }
  tKX tQ inline symt2D<K,X>& symt2D<K,X>::sub_del_prd(symt2D<Q,X> const&p) {
    meta2D::syd<(K-Q)/2,Q>::su(*this,p); return*this; }
  tKX tQ inline symt2D<K,X>& symt2D<K,X>::ass_del_prd(symt2D<Q,X> const&p,
						      X           const&s) {
    meta2D::syd<(K-Q)/2,Q>::as(*this,p,s); return*this; }
  tKX tQ inline symt2D<K,X>& symt2D<K,X>::add_del_prd(symt2D<Q,X> const&p,
						      X           const&s) {
    meta2D::syd<(K-Q)/2,Q>::ad(*this,p,s); return*this; }
  tKX tQ inline symt2D<K,X>& symt2D<K,X>::sub_del_prd(symt2D<Q,X> const&p,
						      X           const&s) {
    meta2D::syd<(K-Q)/2,Q>::su(*this,p,s); return*this; }
  tSY& symt2D<K,X>::ass_del_prd(V const&p) {
    meta2D::syd<(K-1)/2,1>::as(*this,p); return*this; }
  tSY& symt2D<K,X>::add_del_prd(V const&p) {
    meta2D::syd<(K-1)/2,1>::ad(*this,p); return*this; }
  tSY& symt2D<K,X>::sub_del_prd(V const&p) {
    meta2D::syd<(K-1)/2,1>::su(*this,p); return*this; }
  tSY& symt2D<K,X>::ass_del_prd(V const&p, X const&s) {
    meta2D::syd<(K-1)/2,1>::as(*this,p,s); return*this; }
  tSY& symt2D<K,X>::add_del_prd(V const&p, X const&s) {
    meta2D::syd<(K-1)/2,1>::ad(*this,p,s); return*this; }
  tSY& symt2D<K,X>::sub_del_prd(V const&p, X const&s) {
    meta2D::syd<(K-1)/2,1>::su(*this,p,s); return*this; }
  tSY& symt2D<K,X>::ass_del_prd() {
    meta2D::syd<K/2,0>::as(*this); return*this; }
  tSY& symt2D<K,X>::add_del_prd() {
    meta2D::syd<K/2,0>::ad(*this); return*this; }
  tSY& symt2D<K,X>::sub_del_prd() {
    meta2D::syd<K/2,0>::su(*this); return*this; }
  tSY& symt2D<K,X>::ass_del_prd(X const&s) {
    meta2D::syd<K/2,0>::as(*this,s); return*this; }
  tSY& symt2D<K,X>::add_del_prd(X const&s) {
    meta2D::syd<K/2,0>::ad(*this,s); return*this; }
  tSY& symt2D<K,X>::sub_del_prd(X const&s) {
    meta2D::syd<K/2,0>::su(*this,s); return*this; }
  //////////////////////////////////////////////////////////////////////////////
  tS2& symt2D<2,X>::ass_del_prd() {
    meta2D::syd<1,0>::as(*this); return*this; }
  tS2& symt2D<2,X>::add_del_prd() {
    meta2D::syd<1,0>::ad(*this); return*this; }
  tS2& symt2D<2,X>::sub_del_prd() {
    meta2D::syd<1,0>::su(*this); return*this; }
  tS2& symt2D<2,X>::ass_del_prd(X const&s) {
    meta2D::syd<1,0>::as(*this,s); return*this; }
  tS2& symt2D<2,X>::add_del_prd(X const&s) {
    meta2D::syd<1,0>::ad(*this,s); return*this; }
  tS2& symt2D<2,X>::sub_del_prd(X const&s) {
    meta2D::syd<1,0>::su(*this,s); return*this; }
}                                                  // END namespace nbdy        
namespace meta2D {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // complete inner product                                                   //
  //                                                                          //
  // coded as follows:                                                        //
  //                                                                          //
  // 1. loop L1=0...K1 of product using LoopL<>                               //
  // 2. at each L1  loop (kx,ky): kx>=ky, kx+ky=K2                            //
  // 3. at each (kx,ky) loop over all permutations of (kx,ky)                 //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // struct inp_p<>::p() for the product b[i]*c[i2]                             
  //                                                                            
  //----------------------------------------------------------------------------
  tm<int L1,int L2> struct inp_p {
    tX sX p(cX*b,cX*c) { return b[L1+L2] * c[L2]; }
  };
  //----------------------------------------------------------------------------
  //                                                                            
  // struct inp_a<> give some constants                                         
  //                                                                            
  //----------------------------------------------------------------------------
  tm<int Kx,int Ky> struct inp_a {
    enum {
      K = Kx+Ky,
      L = Ky,
      P = TWO<K,L>::B };
  };
  //----------------------------------------------------------------------------
  //                                                                            
  // struct inp_k<>::k() loop over all permutations in (Kx,Ky) space            
  //                                                                            
  //----------------------------------------------------------------------------
  tm<int L1,int Kx,int Ky>
  struct inp_k : public inp_a<Kx,Ky> {                         // Kx,Ky         
    tX sX k(cX*b,cX*c) {
      return times<P>( inp_p<L1,Ky>::p(b,c) + inp_p<L1,Kx>::p(b,c) );
    }
  };
  //----------------------------------------------------------------------------
  tm<int L1,int Kx>
  struct inp_k<L1,Kx,Kx> : public inp_a<Kx,Kx> {               // Ky=Kx         
    tX sX k(cX*b,cX*c) {
      return times<P>( inp_p<L1,Kx>::p(b,c) );
    }
  };
  //----------------------------------------------------------------------------
  //                                                                            
  // struct inpk<>::k() loop over all Kx>=Ky                                    
  //                                                                            
  //----------------------------------------------------------------------------
  tm<int L1,int Kx,int Ky,int KyM>
  struct inpk {
    tX sX k(cX*b,cX*c) {
      return
	inp_k<L1,Kx,  Ky      >::k(b,c) +
	inpk <L1,Kx-1,Ky+1,KyM>::k(b,c);
    }
  };
  //----------------------------------------------------------------------------
  tm<int L1,int Kx,int Ky>                                    // Ky=KyM         
  struct inpk<L1,Kx,Ky,Ky> {
    tX sX k(cX*b,cX*c) {
      return inp_k<L1,Kx,Ky>::k(b,c);
    }
  };
  //----------------------------------------------------------------------------
  //                                                                            
  // struct inp_as<>::job<>(): trigger loop over Kx,Ky via inpk<>::k()          
  // struct inp_ad<>::job<>(): trigger loop over Kx,Ky via inpk<>::k()          
  // struct inp_au<>::job<>(): trigger loop over Kx,Ky via inpk<>::k()          
  //                                                                            
  //----------------------------------------------------------------------------
  tm<int K2,int K1,int L1> struct inp_as {
    tX sv job(X*a,cX*b,cX*c) {
      a[L1] =   inpk<L1,K2,0,K2/2>::k(b,c); }
    tX sv job(X*a,cX*b,cX*c,cX&d) {
      a[L1] = d*inpk<L1,K2,0,K2/2>::k(b,c); }
  };
  //----------------------------------------------------------------------------
  tm<int K2,int K1,int L1> struct inp_ad {
    tX sv job(X*a,cX*b,cX*c) {
      a[L1]+=   inpk<L1,K2,0,K2/2>::k(b,c); }
    tX sv job(X*a,cX*b,cX*c,cX&d) {
      a[L1]+= d*inpk<L1,K2,0,K2/2>::k(b,c); }
  };
  //----------------------------------------------------------------------------
  tm<int K2,int K1,int L1> struct inp_su {
    tX sv job(X*a,cX*b,cX*c) {
      a[L1]-=   inpk<L1,K2,0,K2/2>::k(b,c); }
    tX sv job(X*a,cX*b,cX*c,cX&d) {
      a[L1]-= d*inpk<L1,K2,0,K2/2>::k(b,c); }
  };
  //----------------------------------------------------------------------------
  //                                                                            
  // struct inp<>: perform complete inner product                               
  //                                                                            
  //----------------------------------------------------------------------------
  tm<int K1,int K2> struct inp {
    tX sv as( X*a,cX*b,cX*c)      { LoopL<inp_as,K1,K2>::loop(a,b,c); }
    tX sv ad( X*a,cX*b,cX*c)      { LoopL<inp_ad,K1,K2>::loop(a,b,c); }
    tX sv su( X*a,cX*b,cX*c)      { LoopL<inp_su,K1,K2>::loop(a,b,c); }
    tX sv as( X*a,cX*b,cX*c,cX&d) { LoopL<inp_as,K1,K2>::loop(a,b,c,d); }
    tX sv ad( X*a,cX*b,cX*c,cX&d) { LoopL<inp_ad,K1,K2>::loop(a,b,c,d); }
    tX sv su( X*a,cX*b,cX*c,cX&d) { LoopL<inp_su,K1,K2>::loop(a,b,c,d); }
  };
  tm<int K2> struct inp<0,K2> {
    tX sX ras (cX*b,cX*c) {
      return inpk<0,K2,0,K2/2>::k(b,c);
    }
  };
}                                                  // END namespace meta3D      
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // implementing                                                             //
  //                                                                          //
  // nbdy::symt2D<K,X>::ass_inn_prd();                                        //
  // nbdy::symt2D<K,X>::add_inn_prd();                                        //
  // nbdy::symt2D<K,X>::sub_inn_prd();                                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  tKX inline void symt2D<K,X>::ass_inn_prd(V const&p, symt2D<K-1,X>&q) const {
    meta2D::inp<K-1,1>::as(pX(q),a,pX(p)); }
  tKX inline void symt2D<K,X>::add_inn_prd(V const&p, symt2D<K-1,X>&q) const {
    meta2D::inp<K-1,1>::ad(pX(q),a,pX(p)); }
  tKX inline void symt2D<K,X>::sub_inn_prd(V const&p, symt2D<K-1,X>&q) const {
    meta2D::inp<K-1,1>::su(pX(q),a,pX(p)); }

  tKX tQ inline void symt2D<K,X>::ass_inn_prd(symt2D<Q,X>  const&p,
					      symt2D<K-Q,X>     &q) const {
    meta2D::inp<K-Q,Q >::as(pX(q),a,pX(p)); }
  tKX tQ inline void symt2D<K,X>::add_inn_prd(symt2D<Q,X>  const&p,
					      symt2D<K-Q,X>     &q) const {
    meta2D::inp<K-Q,Q >::ad(pX(q),a,pX(p)); }
  tKX tQ inline void symt2D<K,X>::sub_inn_prd(symt2D<Q,X>  const&p,
					      symt2D<K-Q,X>     &q) const {
    meta2D::inp<K-Q,Q >::su(pX(q),a,pX(p)); }

  tKX inline void symt2D<K,X>::ass_inn_prd(symt2D<K-1,X>const&p, V&q) const {
    meta2D::inp<1,K-1>::as(pX(q),a,pX(p)); }
  tKX inline void symt2D<K,X>::add_inn_prd(symt2D<K-1,X>const&p, V&q) const {
    meta2D::inp<1,K-1>::ad(pX(q),a,pX(p)); }
  tKX inline void symt2D<K,X>::sub_inn_prd(symt2D<K-1,X>const&p, V&q) const {
    meta2D::inp<1,K-1>::su(pX(q),a,pX(p)); }

  tKX inline void symt2D<K,X>::ass_inn_prd(symt2D<K,X>  const&p, X&q) const {
    q = meta2D::inp<0,K>::ras(a,pX(p)); }
  tKX inline void symt2D<K,X>::add_inn_prd(symt2D<K,X>  const&p, X&q) const {
    q+= meta2D::inp<0,K>::ras(a,pX(p)); }
  tKX inline void symt2D<K,X>::sub_inn_prd(symt2D<K,X>  const&p, X&q) const {
    q-= meta2D::inp<0,K>::ras(a,pX(p)); }

  tKX inline X    symt2D<K,X>::ass_inn_prd(symt2D<K,X>  const&p) const {
    meta2D::inp<0,K>::ras(a,pX(p)); }
  //////////////////////////////////////////////////////////////////////////////
  tX inline void symt2D<2,X>::ass_inn_prd(V const&p, V&q) const {
    meta2D::inp<1,1>::as(pX(q),a,pX(p)); }
  tX inline void symt2D<2,X>::add_inn_prd(V const&p, V&q) const {
    meta2D::inp<1,1>::ad(pX(q),a,pX(p)); }
  tX inline void symt2D<2,X>::sub_inn_prd(V const&p, V&q) const {
    meta2D::inp<1,1>::su(pX(q),a,pX(p)); }

  tX inline void symt2D<2,X>::ass_inn_prd(symt2D<2,X> const&p, X&q) const {
    q = meta2D::inp<0,2>::ras(a,pX(p)); }
  tX inline void symt2D<2,X>::add_inn_prd(symt2D<2,X> const&p, X&q) const {
    q+= meta2D::inp<0,2>::ras(a,pX(p)); }
  tX inline void symt2D<2,X>::sub_inn_prd(symt2D<2,X> const&p, X&q) const {
    q-= meta2D::inp<0,2>::ras(a,pX(p)); }

  tX inline X    symt2D<2,X>::ass_inn_prd(symt2D<2,X> const&p) const {
    meta2D::inp<0,2>::ras(a,pX(p)); }
}                                                  // END namespace nbdy        
////////////////////////////////////////////////////////////////////////////////
#undef __min
#undef __max
#undef sci
#undef scb
#undef si
#undef sv
#undef sX
#undef cX
#undef tKX
#undef tSY
#undef tS2
#undef V
#undef tQ
#undef tX
#undef tm
#undef cA
#undef cB
#undef cC
#undef cD
#undef tA
#undef tAB
#undef tABC
#undef tABCD
////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_included_tn2D_cc   
