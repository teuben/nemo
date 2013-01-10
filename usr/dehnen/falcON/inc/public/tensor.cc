// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// tensor.cc                                                                   |
//                                                                             |
// Copyright (C) 2003-2005   Walter Dehnen                                     |
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
#ifndef falcON_included_tensor_cc
#define falcON_included_tensor_cc

#ifndef falcON_included_tensor_h
#  include <public/tensor.h>
#endif
#ifndef Wdutils_included_meta_h
#  include <utils/meta.h>
#endif

////////////////////////////////////////////////////////////////////////////////
#define tm    template
#define tKX   template<int K,typename X>
#define tSY   tKX inline symt3D<K,X>
#define tS2   tX inline symt3D<2,X>
#define V     falcONVec<3,X>
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
//
// macros to be used for minimum and maximum of const integer expressions
//
#undef  __min
#undef  __max
#define __min(A,B) ((A)<(B)? (A):(B))
#define __max(A,B) ((A)>(B)? (A):(B))
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  namespace meta3D {
    using meta::taux;
    using meta::ONE;
    using meta::TWO;
    //
    // struct ONE3D<N,X>
    //
    // generic constants, types and methods
    //
    // ND:  # independent elements in symmetric tensor of order N
    // CD:  # independent elements in symmetric tensors of orders 0 to N
    // I0:  index of first element of Nth tensor of set with orders 0 to N
    // PD:  # independent elements in symmetric tenrors of orders 2 to N
    // J0:  index of first element of Nth tensor of set with orders 2 to N
    //
    tm<int N, typename X=real> struct ONE3D : ONE<N> {
      enum {
	ND = ((N+1)*(N+2))/2,                      // (N+1)(N+2)/2
	CD = ((N+1)*(N+2)*(N+3))/6,                // (N+1)(N+2)(N+3)/6
	I0 = (N*(N+1)*(N+2))/6,                    // N(N+1)(N+2)/6
	PD = CD-4,                                 // (N+1)(N+2)(N+3)/6-4
	J0 = I0-4                                  // N(N+1)(N+2)/6-4
      };
      typedef symt3D<N,X> Tensor;                  // sym 3D tensor of order N
      static Tensor      &tens(      X*a)
      { return*static_cast<      Tensor*>(static_cast<void*>(a+I0)); }
      static Tensor const&tens(const X*a)
      { return*static_cast<const Tensor*>(static_cast<const void*>(a+I0)); }
      static Tensor      &pole(      X*a)
      { return*static_cast<      Tensor*>(static_cast<void*>(a+J0)); }
      static Tensor const&pole(const X*a)
      { return*static_cast<const Tensor*>(static_cast<const void*>(a+J0)); }
    };
    tm<typename X> struct ONE3D<1,X> : ONE<1> {
      enum { ND=3, CD=4, I0=1 };
      typedef falcONVec<3,X> Tensor;
      static Tensor      &tens(      X*a)
      { return reinterpret_cast<Tensor&>(*(a+I0)); }
      static Tensor const&tens(const X*a)
      { return reinterpret_cast<const Tensor&>(*(a+I0)); }
    };
    tm<typename X> struct ONE3D<0,X> : ONE<0> {
      enum { ND=1, CD=1, I0=0 };
      typedef X Tensor;
      static Tensor      &tens(      X*a) { return *a; }
      static Tensor const&tens(const X*a) { return *a; }
    };
    //
    // struct LoopLM<Perform,Z,K>
    //
    // loop l,m  and  perform some templated action
    //
    // template arguments of P: Z(anything),K,L,M
    //
#define T4i template<int,int,int,int> class
    // struct l for looping L=0...K, M=0...L
    tm<T4i P, int Z, int K, int L, int M> struct l {            // K,L,M
      tA    sv loop(A&a) { 
	P<  Z,K,L,M  >::job (a);
	l<P,Z,K,L,M+1>::loop(a); }
      tAB   sv loop(A&a,cB&b) { 
	P<  Z,K,L,M  >::job (a,b);
	l<P,Z,K,L,M+1>::loop(a,b); }
      tABC  sv loop(A&a,cB&b,cC&c) { 
	P<  Z,K,L,M  >::job (a,b,c);
	l<P,Z,K,L,M+1>::loop(a,b,c); }
      tABCD sv loop(A&a,cB&b,cC&c,cD&d) { 
	P<  Z,K,L,M  >::job (a,b,c,d);
	l<P,Z,K,L,M+1>::loop(a,b,c,d); }
    };
    //
    tm<T4i P, int Z, int K, int L> struct l<P,Z,K,L,L> {        // M=L
      tA    sv loop(A&a) { 
	P<  Z,K,L,  L>::job (a);
	l<P,Z,K,L+1,0>::loop(a); }
      tAB   sv loop(A&a,cB&b) { 
	P<  Z,K,L,  L>::job (a,b);
	l<P,Z,K,L+1,0>::loop(a,b); }
      tABC  sv loop(A&a,cB&b,cC&c) { 
	P<  Z,K,L,  L>::job (a,b,c);
	l<P,Z,K,L+1,0>::loop(a,b,c); }
      tABCD sv loop(A&a,cB&b,cC&c,cD&d) { 
	P<  Z,K,L,  L>::job (a,b,c,d);
	l<P,Z,K,L+1,0>::loop(a,b,c,d); }
    };
    //
    tm<T4i P, int Z, int K> struct l<P,Z,K,K,K> {               // M=L=K
      tA    sv loop(A&a)                { P<Z,K,K,K>::job(a); }
      tAB   sv loop(A&a,cB&b)           { P<Z,K,K,K>::job(a,b); }
      tABC  sv loop(A&a,cB&b,cC&c)      { P<Z,K,K,K>::job(a,b,c); }
      tABCD sv loop(A&a,cB&b,cC&c,cD&d) { P<Z,K,K,K>::job(a,b,c,d); }
    };
    //
    tm<T4i P, int K, int Z> struct LoopLM {
      tA    sv loop(A&a)                { l<P,Z,K,0,0>::loop(a); }
      tAB   sv loop(A&a,cB&b)           { l<P,Z,K,0,0>::loop(a,b); }
      tABC  sv loop(A&a,cB&b,cC&c)      { l<P,Z,K,0,0>::loop(a,b,c); }
      tABCD sv loop(A&a,cB&b,cC&c,cD&d) { l<P,Z,K,0,0>::loop(a,b,c,d); }
    };
#undef T4i
    //
    // K fold outer product of vector
    //
    tm<int K> struct vop {
      enum { 
	KK = ((K+1)*(K+2))/2,
	K0 = ( K   *(K+1))/2,
	K1 = ((K-1)* K   )/2 };
      tX sv as(X*a,cX*u,cX*v) {
	taux<X,K0-1>::v_ast(a,   u,   v[0]);
	taux<X,K -1>::v_ast(a+K0,u+K1,v[1]);
	a[KK-1] = u[K0-1]*v[2];
      }
      tX sv as(X*a,cX*u,cX*v,cX&x) {
	taux<X,K0-1>::v_ast(a,   u,   x*v[0]);
	taux<X,K -1>::v_ast(a+K0,u+K1,x*v[1]);
	a[KK-1] = u[K0-1]*x*v[2];
      }
      tX sv ad(X*a,cX*u,cX*v) {
	taux<X,K0-1>::v_adt(a,   u,   v[0]);
	taux<X,K -1>::v_adt(a+K0,u+K1,v[1]);
	a[KK-1]+= u[K0-1]*v[2];
      }
      tX sv ad(X*a,cX*u,cX*v,cX&x) {
	taux<X,K0-1>::v_adt(a,   u,   x*v[0]);
	taux<X,K -1>::v_adt(a+K0,u+K1,x*v[1]);
	a[KK-1]+= u[K0-1]*x*v[2];
      }
      tX sv su(X*a,cX*u,cX*v) {
	taux<X,K0-1>::v_sut(a,   u,   v[0]);
	taux<X,K -1>::v_sut(a+K0,u+K1,v[1]);
	a[KK-1]-= u[K0-1]*v[2];
      }
      tX sv su(X*a,cX*u,cX*v,cX&x) {
	taux<X,K0-1>::v_sut(a,   u,   x*v[0]);
	taux<X,K -1>::v_sut(a+K0,u+K1,x*v[1]);
	a[KK-1]-= u[K0-1]*x*v[2];
      }
    };
    //
    tm<> struct vop<2> {
      tX sv as(X*a,cX*v) {
	taux<X,2>::v_ast(a,  v,  v[0]);
	taux<X,1>::v_ast(a+3,v+1,v[1]);
	a[5] = v[2]*v[2];
      }
      tX sv as(X*a,cX*v,cX&x) {
	taux<X,2>::v_ast(a,  v,  x*v[0]);
	taux<X,1>::v_ast(a+3,v+1,x*v[1]);
	a[5] = x*v[2]*v[2];
      }
      tX sv ad(X*a,cX*v) {
	taux<X,2>::v_adt(a,  v,  v[0]);
	taux<X,1>::v_adt(a+3,v+1,v[1]);
	a[5]+= v[2]*v[2];
      }
      tX sv ad(X*a,cX*v,cX&x) {
	taux<X,2>::v_adt(a,  v,  x*v[0]);
	taux<X,1>::v_adt(a+3,v+1,x*v[1]);
	a[5]+= x*v[2]*v[2];
      }
      tX sv su(X*a,cX*v) {
	taux<X,2>::v_sut(a,  v,  v[0]);
	taux<X,1>::v_sut(a+3,v+1,v[1]);
	a[5]-= v[2]*v[2];
      }
      tX sv su(X*a,cX*v,cX&x) {
	taux<X,2>::v_sut(a,  v,  x*v[0]);
	taux<X,1>::v_sut(a+3,v+1,x*v[1]);
	a[5]-= x*v[2]*v[2];
      }
    };
  } // namespace meta3D
  //
  // implementing
  //
  // falcON::symt3D<K,S>::ass_out_prd();
  // falcON::symt3D<K,S>::add_out_prd();
  // falcON::symt3D<K,S>::sub_out_prd();
  //
  tSY& symt3D<K,X>::ass_out_prd(symt3D<K-1,X>const&u, V const&v) {
    meta3D::vop<K>::as(a,pX(u),pX(v)); return*this; }
  tSY& symt3D<K,X>::add_out_prd(symt3D<K-1,X>const&u, V const&v) {
    meta3D::vop<K>::ad(a,pX(u),pX(v)); return*this; }
  tSY& symt3D<K,X>::sub_out_prd(symt3D<K-1,X>const&u, V const&v) {
    meta3D::vop<K>::su(a,pX(u),pX(v)); return*this; }
  tSY& symt3D<K,X>::ass_out_prd(symt3D<K-1,X>const&u, V const&v, cX&s) {
    meta3D::vop<K>::as(a,pX(u),pX(v),s); return*this; }
  tSY& symt3D<K,X>::add_out_prd(symt3D<K-1,X>const&u, V const&v, cX&s) {
    meta3D::vop<K>::ad(a,pX(u),pX(v),s); return*this; }
  tSY& symt3D<K,X>::sub_out_prd(symt3D<K-1,X>const&u, V const&v, cX&s) {
    meta3D::vop<K>::su(a,pX(u),pX(v),s); return*this; }
  //
  tS2& symt3D<2,X>::ass_out_prd(V const&v) {
    meta3D::vop<2>::as(a,pX(v)); return*this; }
  tS2& symt3D<2,X>::add_out_prd(V const&v) {
    meta3D::vop<2>::ad(a,pX(v)); return*this; }
  tS2& symt3D<2,X>::sub_out_prd(V const&v) {
    meta3D::vop<2>::su(a,pX(v)); return*this; }
  tS2& symt3D<2,X>::ass_out_prd(V const&v, cX&s) {
    meta3D::vop<2>::as(a,pX(v),s); return*this; }
  tS2& symt3D<2,X>::add_out_prd(V const&v, cX&s) {
    meta3D::vop<2>::ad(a,pX(v),s); return*this; }
  tS2& symt3D<2,X>::sub_out_prd(V const&v, cX&s) {
    meta3D::vop<2>::su(a,pX(v),s); return*this; }
  //
  namespace meta3D {
    //
    // symmetric outer product
    //
    // coded as follows:
    //
    // 1. loop (L,M) of product using LoopLM<>
    // 2. at each (L,M)    loop L1 = max(0,L-K2) ... min(L,K1)
    // 3. at each (L,M,L1) loop M1 = max(0,M-L2) ... min(M,L1)
    //

    //
    // struct syp_p<> computes single contribution to symmetric outer product
    //
    // static const int P:  factor in symmetric outer product
    // p(a,b): return contribution to symmetric outer product
    //
    tm<int K, int L, int M, int K1, int L1, int M1>
    struct syp_p {
      sci P  = TWO<K-L,K1-L1>::B * TWO<L-M,L1-M1>::B * TWO<M,M1>::B;
      sci I1 = TWO<L1,  M1  >::I;
      sci I2 = TWO<L-L1,M-M1>::I;
      tX sX p(cX*b,cX*c) { return times<P>( b[I1] * c[I2] ); }
    };
    //
    // struct syp_m<>::m<>() for the loop over M1
    //
    tm<int K,int L,int M,int K1,int L1,int M1,int M1M>
    struct syp_m {
      tX sX m(cX*b,cX*c) {
	return
	  syp_p<K,L,M,K1,L1,M1      >::p(b,c) + 
	  syp_m<K,L,M,K1,L1,M1+1,M1M>::m(b,c);
      }
    };
    //
    tm<int K,int L,int M,int K1,int L1,int M1>
    struct syp_m<K,L,M,K1,L1,M1,M1> {
      tX sX m(cX*b,cX*c) { return syp_p<K,L,M,K1,L1,M1>::p(b,c); }
    };
    //
    // struct syp_l<>::l<>() for the loop over L1
    //
    tm<int K, int L, int M, int K1, int L1, int L1M> struct syp_l {
      sci M1_LO=__max(0,M-L+L1), M1_HI=__min(M,L1);
      tX sX l(cX*b,cX*c) {
	return
	  syp_m<K,L,M,K1,L1,M1_LO,M1_HI>::m(b,c)	+
	  syp_l<K,L,M,K1,L1+1,L1M      >::l(b,c);
      }
    };
    //
    tm<int K, int L, int M, int K1, int L1> struct syp_l<K,L,M,K1,L1,L1> {
      sci M1_LO=__max(0,M-L+L1), M1_HI=__min(M,L1);
      tX sX l(cX*b,cX*c) {
	return syp_m<K,L,M, K1,L1,M1_LO,M1_HI>::m(b,c);
      }
    };
    //
    // struct syp_as<>: trigger loop over L1 via syp_l<>::l()
    // struct syp_ad<>: trigger loop over L1 via syp_l<>::l()
    // struct syp_su<>: trigger loop over L1 via syp_l<>::l()
    //
    tm<int K1,int K,int L,int M> struct syp_as {
      sci I = TWO<L,M>::I;
      tX sv job(X*a,cX*b,cX*c) {
	a[I] =    syp_l<K,L,M,K1, __max(0,L-K+K1), __min(L,K1)>::l(b,c); }
      tX sv job(X*a,cX*b,cX*c,cX&d) {
	a[I] =  d*syp_l<K,L,M,K1, __max(0,L-K+K1), __min(L,K1)>::l(b,c); }
    };
    //
    tm<int K1,int K,int L,int M> struct syp_ad {
      sci I = TWO<L,M>::I;
      tX sv job(X*a,cX*b,cX*c) {
	a[I]+=    syp_l<K,L,M,K1, __max(0,L-K+K1), __min(L,K1)>::l(b,c); }
      tX sv job(X*a,cX*b,cX*c,cX&d) {
	a[I]+=  d*syp_l<K,L,M,K1, __max(0,L-K+K1), __min(L,K1)>::l(b,c); }
    };
    //
    tm<int K1,int K,int L,int M> struct syp_su {
      sci I = TWO<L,M>::I;
      tX sv job(X*a,cX*b,cX*c) {
	a[I]-=    syp_l<K,L,M,K1, __max(0,L-K+K1), __min(L,K1)>::l(b,c); }
      tX sv job(X*a,cX*b,cX*c,cX&d) {
	a[I]-=  d*syp_l<K,L,M,K1, __max(0,L-K+K1), __min(L,K1)>::l(b,c); }
    };
    //
    // struct syp<>: perform symmetric outer product
    //
    tm<int K1,int K2> struct syp {
      tX sv as (X*a,cX*b,cX*c)      { LoopLM<syp_as,K1+K2,K1>::loop(a,b,c); }
      tX sv ad (X*a,cX*b,cX*c)      { LoopLM<syp_ad,K1+K2,K1>::loop(a,b,c); }
      tX sv su (X*a,cX*b,cX*c)      { LoopLM<syp_su,K1+K2,K1>::loop(a,b,c); }
      tX sv ast(X*a,cX*b,cX*c,cX&d) { LoopLM<syp_as,K1+K2,K1>::loop(a,b,c,d); }
      tX sv adt(X*a,cX*b,cX*c,cX&d) { LoopLM<syp_ad,K1+K2,K1>::loop(a,b,c,d); }
      tX sv sut(X*a,cX*b,cX*c,cX&d) { LoopLM<syp_su,K1+K2,K1>::loop(a,b,c,d); }
    };
  } // namespace meta3D
  //
  // implementing
  //
  // falcON::symt3D<K,X>::ass_sym_prd();
  // falcON::symt3D<K,X>::add_sym_prd();
  // falcON::symt3D<K,X>::sub_sym_prd();
  //
  tKX tQ inline symt3D<K,X>& symt3D<K,X>::ass_sym_prd(symt3D<K-Q,X>const&p,
						      symt3D<Q  ,X>const&q) {
    meta3D::syp<K-Q,Q>::as(a,pX(p),pX(q)); return*this; }
  tKX tQ inline symt3D<K,X>& symt3D<K,X>::add_sym_prd(symt3D<K-Q,X>const&p,
						      symt3D<Q  ,X>const&q) {
    meta3D::syp<K-Q,Q>::ad(a,pX(p),pX(q)); return*this; }
  tKX tQ inline symt3D<K,X>& symt3D<K,X>::sub_sym_prd(symt3D<K-Q,X>const&p,
						      symt3D<Q  ,X>const&q) {
    meta3D::syp<K-Q,Q>::su(a,pX(p),pX(q)); return*this; }
  //
  tKX tQ inline symt3D<K,X>& symt3D<K,X>::ass_sym_prd(symt3D<K-Q,X>const&p,
						      symt3D<Q  ,X>const&q,
						      X            const&s) {
    meta3D::syp<K-Q,Q>::as(a,pX(p),pX(q),s); return*this; }
  tKX tQ inline symt3D<K,X>& symt3D<K,X>::add_sym_prd(symt3D<K-Q,X>const&p,
						      symt3D<Q  ,X>const&q,
						      X            const&s) {
    meta3D::syp<K-Q,Q>::ad(a,pX(p),pX(q),s); return*this; }
  tKX tQ inline symt3D<K,X>& symt3D<K,X>::sub_sym_prd(symt3D<K-Q,X>const&p,
						      symt3D<Q  ,X>const&q,
						      X            const&s) {
    meta3D::syp<K-Q,Q>::su(a,pX(p),pX(q),s); return*this; }
  //
  tSY& symt3D<K,X>::ass_sym_prd(symt3D<K-1,X>const&p,V const&q) {
    meta3D::syp<K-1,1>::as(a,pX(p),pX(q)); return*this; }
  tSY& symt3D<K,X>::add_sym_prd(symt3D<K-1,X>const&p,V const&q) {
    meta3D::syp<K-1,1>::ad(a,pX(p),pX(q)); return*this; }
  tSY& symt3D<K,X>::sub_sym_prd(symt3D<K-1,X>const&p,V const&q) {
    meta3D::syp<K-1,1>::su(a,pX(p),pX(q)); return*this; }
  tSY& symt3D<K,X>::ass_sym_prd(V const&q,symt3D<K-1,X>const&p) {
    meta3D::syp<K-1,1>::as(a,pX(p),pX(q)); return*this; }
  tSY& symt3D<K,X>::add_sym_prd(V const&q,symt3D<K-1,X>const&p) {
    meta3D::syp<K-1,1>::ad(a,pX(p),pX(q)); return*this; }
  tSY& symt3D<K,X>::sub_sym_prd(V const&q,symt3D<K-1,X>const&p) {
    meta3D::syp<K-1,1>::su(a,pX(p),pX(q)); return*this; }
  //
  tSY& symt3D<K,X>::ass_sym_prd(symt3D<K-1,X>const&p,V const&q,cX&s) {
    meta3D::syp<K-1,1>::as(a,pX(p),pX(q),s); return*this; }
  tSY& symt3D<K,X>::add_sym_prd(symt3D<K-1,X>const&p,V const&q,cX&s) {
    meta3D::syp<K-1,1>::ad(a,pX(p),pX(q),s); return*this; }
  tSY& symt3D<K,X>::sub_sym_prd(symt3D<K-1,X>const&p,V const&q,cX&s) {
    meta3D::syp<K-1,1>::su(a,pX(p),pX(q),s); return*this; }
  tSY& symt3D<K,X>::ass_sym_prd(V const&q,symt3D<K-1,X>const&p,cX&s) {
    meta3D::syp<K-1,1>::as(a,pX(p),pX(q),s); return*this; }
  tSY& symt3D<K,X>::add_sym_prd(V const&q,symt3D<K-1,X>const&p,cX&s) {
    meta3D::syp<K-1,1>::ad(a,pX(p),pX(q),s); return*this; }
  tSY& symt3D<K,X>::sub_sym_prd(V const&q,symt3D<K-1,X>const&p,cX&s) {
    meta3D::syp<K-1,1>::su(a,pX(p),pX(q),s); return*this; }
  //
  tS2& symt3D<2,X>::ass_sym_prd(V const&p,V const&q) {
    meta3D::syp<1,1>::as(a,pX(p),pX(q)); return*this; }
  tS2& symt3D<2,X>::add_sym_prd(V const&p,V const&q) {
    meta3D::syp<1,1>::ad(a,pX(p),pX(q)); return*this; }
  tS2& symt3D<2,X>::sub_sym_prd(V const&p,V const&q) {
    meta3D::syp<1,1>::su(a,pX(p),pX(q)); return*this; }
  //
  tS2& symt3D<2,X>::ass_sym_prd(V const&p,V const&q,cX&s) {
    meta3D::syp<1,1>::as(a,pX(p),pX(q),s); return*this; }
  tS2& symt3D<2,X>::add_sym_prd(V const&p,V const&q,cX&s) {
    meta3D::syp<1,1>::ad(a,pX(p),pX(q),s); return*this; }
  tS2& symt3D<2,X>::sub_sym_prd(V const&p,V const&q,cX&s) {
    meta3D::syp<1,1>::su(a,pX(p),pX(q),s); return*this; }
  //////////////////////////////////////////////////////////////////////////////
  namespace meta3D {
    //
    // symmetric outer product with delta_ij^(N)
    //
    // coded as follows:
    //
    // 1. loop (L,M) of product using LoopLM<>
    // 2. at each (L,M)    loop l1 = max(0,(L-K2+1)/2) ... min(L,K1)/2
    // 3. at each (L,M,l1) loop m1 = max(0,(m-L2+1)/2) ... min(M,L1)/2
    //
    // where l1,m1 = (L1,M1)/2.
    //

    //
    // N fold outer product of delta_ij
    //
    tm<int k, int l,int m> struct D {
      enum {
	K = k+k,
	L = l+l,
	M = m+m,
	Z = ONE<k-l>::H * ONE<l-m>::H * ONE<m>::H
      }; };
    //
    // struct syd_p<K,L,M,k1,l1,m1>
    //
    // static const int P:  factor in symmetric outer product
    // p(a,b): return contribution to symmetric outer product
    //
    tm<int K,int L, int M, int k1, int l1, int m1> struct syd_p {
      sci K1=k1+k1, L1=l1+l1, M1=m1+m1;
      sci P =D<k1,l1,m1>::Z *TWO<K-L,K1-L1>::B *TWO<L-M,L1-M1>::B *TWO<M,M1>::B;
      sci I2=TWO<L-L1,M-M1>::I;
      tX sX p(cX*b) { return times<P>( b[I2] ); }
    };
    //
    // struct syd_m<>::m<>() for the loop over m1
    //
    tm<int K,int L,int M,int k1,int l1,int m1,int m1M> struct syd_m {
      tX sX m(cX*b) {
	return
	  syd_p<K,L,M,k1,l1,m1      >::p(b) +
	  syd_m<K,L,M,k1,l1,m1+1,m1M>::m(b);
      }
    };
    //
    tm<int K,int L,int M,int k1,int l1,int m1>
    struct syd_m<K,L,M,k1,l1,m1,m1> {
      tX sX m(cX*b) { return syd_p<K,L,M,k1,l1,m1>::p(b); }
    };
    //
    // struct syd_l<>::l<>() for the loop over l1
    //
    tm<int,int,int,int,int,int> struct syd_l;
    tm<int,int,int,int,int,int,int,int,bool> struct syd__l;
    //
    tm<int K,int L,int M,int k1,int l1,int l1M,int m1l,int m1h>
    struct syd__l<K,L,M,k1,l1,l1M,m1l,m1h,true> {
      tX sX l(cX*b) { return
	  syd_m<K,L,M,k1,l1,m1l,m1h>::m(b) +
	  syd_l<K,L,M,k1,l1+1,l1M  >::l(b); }
    };
    //
    tm<int K,int L,int M,int k1,int l1,int l1M,int m1l,int m1h>
    struct syd__l<K,L,M,k1,l1,l1M,m1l,m1h,false> {
      tX sX l(cX*b) { return syd_l<K,L,M,k1,l1+1,l1M>::l(b); }
    };
    //
    tm<int K,int L,int M,int k1,int l1,int m1l,int m1h>
    struct syd__l<K,L,M,k1,l1,l1,m1l,m1h,true> {
      tX sX l(cX*b) { return syd_m<K,L,M,k1,l1,m1l,m1h>::m(b); }
    };
    //
    tm<int K,int L,int M,int k1,int l1,int m1l,int m1h>
    struct syd__l<K,L,M,k1,l1,l1,m1l,m1h,false> {
      tX sX l(cX*b) { return X(0); }
    };
    //
    tm<int K, int L, int M, int k1, int l1, int l1M> struct syd_l {
      sci L2  =L-l1-l1;
      sci m1l =__max(0,(M-L2+1)/2), m1h=__min(M/2,l1);
      scb loop=m1l <= m1h;
      tX sX l(cX*b) { return syd__l<K,L,M,k1,l1,l1M,m1l,m1h,loop>::l(b); }
    };
    //
    // struct syd_as<>: trigger loop over L1 via syd_l<>::l()
    // struct syd_ad<>: trigger loop over L1 via syd_l<>::l()
    // struct syd_su<>: trigger loop over L1 via syd_l<>::l()
    //
    tm<int k1,int K,int L,int M> struct syd_as {
      sci I     = TWO<L,M>::I;
      sci K2    = K-k1-k1;
      sci l1_lo = __max(0,(L-K2+1)/2), l1_hi=__min(L/2,k1);
      scb loop  = l1_lo<=l1_hi;
      tX sv job(X*a,cX*b) {
	a[I]= loop?   syd_l<K,L,M,k1,l1_lo,l1_hi>::l(b) : X(0) ;
      }
      tX sv job(X*a,cX*b,cX&d) {
	a[I]= loop? d*syd_l<K,L,M,k1,l1_lo,l1_hi>::l(b) : X(0) ;
      }
    };
    //
    tm<int k1,int K,int L,int M> struct syd_ad {
      sci I     = TWO<L,M>::I;
      sci K2    = K-k1-k1;
      sci l1_lo = __max(0,(L-K2+1)/2), l1_hi=__min(L/2,k1);
      scb loop  = l1_lo<=l1_hi;
      tX sv job(X*a,cX*b) {
	if(loop) a[I]+=   syd_l<K,L,M,k1,l1_lo,l1_hi>::l(b);
      }
      tX sv job(X*a,cX*b,cX&d) {
	if(loop) a[I]+= d*syd_l<K,L,M,k1,l1_lo,l1_hi>::l(b);
      }
    };
    //
    tm<int k1,int K,int L,int M> struct syd_su {
      sci I     = TWO<L,M>::I;
      sci K2    = K-k1-k1;
      sci l1_lo = __max(0,(L-K2+1)/2), l1_hi=__min(L/2,k1);
      scb loop  = l1_lo<=l1_hi;
      tX sv job(X*a,cX*b) {
	if(loop) a[I]-=   syd_l<K,L,M,k1,l1_lo,l1_hi>::l(b);
      }
      tX sv job(X*a,cX*b,cX&d) {
	if(loop) a[I]-= d*syd_l<K,L,M,k1,l1_lo,l1_hi>::l(b);
      }
    };
    //
    // struct syd_as0: for equating to delta^(k)
    // struct syd_ad0: for equating to delta^(k)
    // struct syd_su0: for equating to delta^(k)
    //
    tm<int Q,int K,int L,int M> struct syd_as0 : public D<K/2,L/2,M/2>{
      sci k = K/2, l=L/2, m=M/2;
      sci I = TWO<L,M>::I;
      tX sv job(X*a)      {
	a[I] = (!(L&1 || M&1))?    D<K/2,L/2,M/2>::Z : X(0) ; }
      tX sv job(X*a,cX&d) {
	a[I] = (!(L&1 || M&1))? d* D<K/2,L/2,M/2>::Z : X(0) ; }
    };
    //
    tm<int Q,int k,int l,int m> struct syd_ad0 : public D<k,l,m> {
      sci I = TWO<D<k,l,m>::L, D<k,l,m>::M>::I;
      tX sv job(X*a)      { a[I]+=    D<k,l,m>::Z; }
      tX sv job(X*a,cX&d) { a[I]+= d* D<k,l,m>::Z; }
    };
    //
    tm<int Q,int k,int l,int m> struct syd_su0 : public D<k,l,m> {
      sci I = TWO<D<k,l,m>::L,D<k,l,m>::M>::I;
      tX sv job(X*a)      { a[I]-=    D<k,l,m>::Z; }
      tX sv job(X*a,cX&d) { a[I]-= d* D<k,l,m>::Z; }
    };
    //
    // struct syd<>: perform symmetric outer product with delta^k
    //
    tm<int k, int K2> struct syd {
      sci K=k+k+K2;
      tX sv asp(X*a,cX*b)      { LoopLM<syd_as,K,k>::loop(a,b); }
      tX sv adp(X*a,cX*b)      { LoopLM<syd_ad,K,k>::loop(a,b); }
      tX sv sup(X*a,cX*b)      { LoopLM<syd_su,K,k>::loop(a,b); }
      tX sv asp(X*a,cX*b,cX&s) { LoopLM<syd_as,K,k>::loop(a,b,s); }
      tX sv adp(X*a,cX*b,cX&s) { LoopLM<syd_ad,K,k>::loop(a,b,s); }
      tX sv sup(X*a,cX*b,cX&s) { LoopLM<syd_su,K,k>::loop(a,b,s); }
      tX sv as(symt3D<K,X>&A,symt3D<K2,X> const&B)      { asp((X*)A,(cX*)B); }
      tX sv ad(symt3D<K,X>&A,symt3D<K2,X> const&B)      { adp((X*)A,(cX*)B); }
      tX sv su(symt3D<K,X>&A,symt3D<K2,X> const&B)      { sup((X*)A,(cX*)B); }
      tX sv as(symt3D<K,X>&A,symt3D<K2,X> const&B,cX&s) { asp((X*)A,(cX*)B,s); }
      tX sv ad(symt3D<K,X>&A,symt3D<K2,X> const&B,cX&s) { adp((X*)A,(cX*)B,s); }
      tX sv su(symt3D<K,X>&A,symt3D<K2,X> const&B,cX&s) { sup((X*)A,(cX*)B,s); }
    };
    //
    tm<int k> struct syd<k,1> {
      sci K=k+k+1;
      tX sv asp(X*a,cX*b)      { LoopLM<syd_as,K,k>::loop(a,b); }
      tX sv adp(X*a,cX*b)      { LoopLM<syd_ad,K,k>::loop(a,b); }
      tX sv sup(X*a,cX*b)      { LoopLM<syd_su,K,k>::loop(a,b); }
      tX sv asp(X*a,cX*b,cX&s) { LoopLM<syd_as,K,k>::loop(a,b,s); }
      tX sv adp(X*a,cX*b,cX&s) { LoopLM<syd_ad,K,k>::loop(a,b,s); }
      tX sv sup(X*a,cX*b,cX&s) { LoopLM<syd_su,K,k>::loop(a,b,s); }
      tX sv as(symt3D<K,X>&A,V const&B)      { asp((X*)A,(cX*)B); }
      tX sv ad(symt3D<K,X>&A,V const&B)      { adp((X*)A,(cX*)B); }
      tX sv su(symt3D<K,X>&A,V const&B)      { sup((X*)A,(cX*)B); }
      tX sv as(symt3D<K,X>&A,V const&B,cX&s) { asp((X*)A,(cX*)B,s); }
      tX sv ad(symt3D<K,X>&A,V const&B,cX&s) { adp((X*)A,(cX*)B,s); }
      tX sv su(symt3D<K,X>&A,V const&B,cX&s) { sup((X*)A,(cX*)B,s); }
    };
    //
    tm<int k> struct syd<k,0> {
      sci K=k+k;
      tX sv asp(X*a)      { LoopLM<syd_as0,K,0>::loop(a); }
      tX sv asp(X*a,cX&s) { LoopLM<syd_as0,K,0>::loop(a,s); }
      tX sv adp(X*a)      { LoopLM<syd_ad0,k,0>::loop(a); }
      tX sv adp(X*a,cX&s) { LoopLM<syd_ad0,k,0>::loop(a,s); }
      tX sv sup(X*a)      { LoopLM<syd_su0,k,0>::loop(a); }
      tX sv sup(X*a,cX&s) { LoopLM<syd_su0,k,0>::loop(a,s); }
      tX sv as(symt3D<K,X>&A)      { asp((X*)A); }
      tX sv as(symt3D<K,X>&A,cX&s) { asp((X*)A,s); }
      tX sv ad(symt3D<K,X>&A)      { adp((X*)A); }
      tX sv ad(symt3D<K,X>&A,cX&s) { adp((X*)A,s); }
      tX sv su(symt3D<K,X>&A)      { sup((X*)A); }
      tX sv su(symt3D<K,X>&A,cX&s) { sup((X*)A,s); }
    };
  } // namespace meta3D
  //
  // implementing
  //
  // falcON::symt3D<K,X>::ass_del_prd();
  // falcON::symt3D<K,X>::add_del_prd();
  // falcON::symt3D<K,X>::sub_del_prd();
  //
  tKX tQ inline symt3D<K,X>& symt3D<K,X>::ass_del_prd(symt3D<Q,X> const&p) {
    meta3D::syd<(K-Q)/2,Q>::as(*this,p); return*this; }
  tKX tQ inline symt3D<K,X>& symt3D<K,X>::add_del_prd(symt3D<Q,X> const&p) {
    meta3D::syd<(K-Q)/2,Q>::ad(*this,p); return*this; }
  tKX tQ inline symt3D<K,X>& symt3D<K,X>::sub_del_prd(symt3D<Q,X> const&p) {
    meta3D::syd<(K-Q)/2,Q>::su(*this,p); return*this; }
  tKX tQ inline symt3D<K,X>& symt3D<K,X>::ass_del_prd(symt3D<Q,X> const&p,
						      X           const&s) {
    meta3D::syd<(K-Q)/2,Q>::as(*this,p,s); return*this; }
  tKX tQ inline symt3D<K,X>& symt3D<K,X>::add_del_prd(symt3D<Q,X> const&p,
						      X           const&s) {
    meta3D::syd<(K-Q)/2,Q>::ad(*this,p,s); return*this; }
  tKX tQ inline symt3D<K,X>& symt3D<K,X>::sub_del_prd(symt3D<Q,X> const&p,
						      X           const&s) {
    meta3D::syd<(K-Q)/2,Q>::su(*this,p,s); return*this; }
  tSY& symt3D<K,X>::ass_del_prd(V const&p) {
    meta3D::syd<(K-1)/2,1>::as(*this,p); return*this; }
  tSY& symt3D<K,X>::add_del_prd(V const&p) {
    meta3D::syd<(K-1)/2,1>::ad(*this,p); return*this; }
  tSY& symt3D<K,X>::sub_del_prd(V const&p) {
    meta3D::syd<(K-1)/2,1>::su(*this,p); return*this; }
  tSY& symt3D<K,X>::ass_del_prd(V const&p, X const&s) {
    meta3D::syd<(K-1)/2,1>::as(*this,p,s); return*this; }
  tSY& symt3D<K,X>::add_del_prd(V const&p, X const&s) {
    meta3D::syd<(K-1)/2,1>::ad(*this,p,s); return*this; }
  tSY& symt3D<K,X>::sub_del_prd(V const&p, X const&s) {
    meta3D::syd<(K-1)/2,1>::su(*this,p,s); return*this; }
  tSY& symt3D<K,X>::ass_del_prd() {
    meta3D::syd<K/2,0>::as(*this); return*this; }
  tSY& symt3D<K,X>::add_del_prd() {
    meta3D::syd<K/2,0>::ad(*this); return*this; }
  tSY& symt3D<K,X>::sub_del_prd() {
    meta3D::syd<K/2,0>::su(*this); return*this; }
  tSY& symt3D<K,X>::ass_del_prd(X const&s) {
    meta3D::syd<K/2,0>::as(*this,s); return*this; }
  tSY& symt3D<K,X>::add_del_prd(X const&s) {
    meta3D::syd<K/2,0>::ad(*this,s); return*this; }
  tSY& symt3D<K,X>::sub_del_prd(X const&s) {
    meta3D::syd<K/2,0>::su(*this,s); return*this; }
  //
  tS2& symt3D<2,X>::ass_del_prd() {
    meta3D::syd<1,0>::as(*this); return*this; }
  tS2& symt3D<2,X>::add_del_prd() {
    meta3D::syd<1,0>::ad(*this); return*this; }
  tS2& symt3D<2,X>::sub_del_prd() {
    meta3D::syd<1,0>::su(*this); return*this; }
  tS2& symt3D<2,X>::ass_del_prd(X const&s) {
    meta3D::syd<1,0>::as(*this,s); return*this; }
  tS2& symt3D<2,X>::add_del_prd(X const&s) {
    meta3D::syd<1,0>::ad(*this,s); return*this; }
  tS2& symt3D<2,X>::sub_del_prd(X const&s) {
    meta3D::syd<1,0>::su(*this,s); return*this; }
  //
  namespace meta3D {
    //
    // complete inner product
    //
    // coded as follows:
    //
    // 1. loop (l1,m1) of product using LoopLM<>
    // 2. at each (l1,m1)  loop (kx,ky,kz): kx>=ky>=kz, kx+ky+kz=K2
    // 3. at each (kx,ky,kz) loop over all permutations of (kx,ky,kz)
    //

    //
    // struct inp_p<>::p() for the product b[i]*c[i2]
    //
    tm<int L1,int M1,int Ky,int M2> struct inp_p {
      enum {
	L2 = Ky+M2,
	I  = TWO<L1+L2,M1+M2>::I,
	I2 = TWO<L2,M2>::I
      };
      tX sX p(cX*b,cX*c) { return b[I] * c[I2]; }
    };
    //
    // struct inp_a<> give some constants
    //
    tm<int Kx,int Ky,int Kz> struct inp_a {
      enum {
	K  = Kx+Ky+Kz,
	L  = Ky+Kz,
	M  = Kz,
	kl = TWO<K,L>::B,
	lm = TWO<L,M>::B,
	P  = kl*lm
      };
    };
    //
    // struct inp_k<>::k() loop over all permutations in (Kx,Ky,Kz) space
    //
    tm<int L1,int M1,int Kx,int Ky,int Kz>
    struct inp_k {                                             // Kx,Ky,Kz
      tX sX k(cX*b,cX*c) {
	return times<inp_a<Kx,Ky,Kz>::P>( inp_p<L1,M1,Ky,Kz>::p(b,c) +
					  inp_p<L1,M1,Kz,Kx>::p(b,c) +
					  inp_p<L1,M1,Kx,Ky>::p(b,c) +
					  inp_p<L1,M1,Ky,Kx>::p(b,c) +
					  inp_p<L1,M1,Kx,Kz>::p(b,c) +
					  inp_p<L1,M1,Kz,Ky>::p(b,c) );
      }
    };
    //
    tm<int L1,int M1,int Kx,int Ky>
    struct inp_k<L1,M1,Kx,Ky,Kx> {                             // Kz=Kx != Ky
      tX sX k(cX*b,cX*c) {
	return times<inp_a<Kx,Ky,Kx>::P>( inp_p<L1,M1,Ky,Kx>::p(b,c) +
					  inp_p<L1,M1,Kx,Kx>::p(b,c) +
					  inp_p<L1,M1,Kx,Ky>::p(b,c) );
      }
    };
    //
    tm<int L1,int M1,int Kx,int Ky>
    struct inp_k<L1,M1,Kx,Ky,Ky> {                             // Kz=Ky != Kx
      tX sX k(cX*b,cX*c) {
	return times<inp_a<Kx,Ky,Ky>::P>( inp_p<L1,M1,Ky,Ky>::p(b,c) +
					  inp_p<L1,M1,Ky,Kx>::p(b,c) +
					  inp_p<L1,M1,Kx,Ky>::p(b,c) );
      }
    };
    //
    tm<int L1,int M1,int Kx,int Kz>
    struct inp_k<L1,M1,Kx,Kx,Kz> {                             // Ky=Kx != Kz
      tX sX k(cX*b,cX*c) {
	return times<inp_a<Kx,Kx,Kz>::P>( inp_p<L1,M1,Kx,Kz>::p(b,c) +
					  inp_p<L1,M1,Kz,Kx>::p(b,c) +
					  inp_p<L1,M1,Kx,Kx>::p(b,c) );
      }
    };
    //
    tm<int L1,int M1,int Kx>
    struct inp_k<L1,M1,Kx,Kx,Kx> {                             // Kx=Ky=Kz
      tX sX k(cX*b,cX*c) {
	return times<inp_a<Kx,Kx,Kx>::P>( inp_p<L1,M1,Kx,Kx>::p(b,c) );
      }
    };
    //
    // struct inpk<>::k() loop over all Kx>=Ky>=Kz
    //
    tm<int L1,int M1,int Kx,int Ky,int Kz,int KyM,int KzM>
    struct inpk {
      tX sX k(cX*b,cX*c) {
	return
	  inp_k<L1,M1,Kx,  Ky,  Kz        >::k(b,c) +
	  inpk <L1,M1,Kx-1,Ky+1,Kz,KyM,KzM>::k(b,c);
      }
    };
    //
    tm<int L1,int M1,int Kx,int Ky,int Kz,int KzM>            // Ky=KyM
    struct inpk<L1,M1,Kx,Ky,Kz,Ky,KzM> {
      tX sX k(cX*b,cX*c) {
	return
	  inp_k<L1,M1,Kx,        Ky,  Kz                  >::k(b,c) +
	  inpk <L1,M1,Kx+Ky-Kz-2,Kz+1,Kz+1,(Kx+Ky-1)/2,KzM>::k(b,c);
      }
    };
    //
    tm<int L1,int M1,int Kx,int Ky,int Kz>                    // Kz=KzM, Ky=KyM
    struct inpk<L1,M1,Kx,Ky,Kz,Ky,Kz> {
      tX sX k(cX*b,cX*c) {
	return inp_k<L1,M1,Kx,Ky,Kz>::k(b,c);
      }
    };
    //
    // struct inp_as<>::job<>(): trigger loop over Kx,Ky,Kz via inpk<>::k()
    // struct inp_ad<>::job<>(): trigger loop over Kx,Ky,Kz via inpk<>::k()
    // struct inp_au<>::job<>(): trigger loop over Kx,Ky,Kz via inpk<>::k()
    //
    tm<int K2,int K1,int L1,int M1> struct inp_as {
      sci I1 = TWO<L1,M1>::I;
      tX sv job(X*a,cX*b,cX*c) {
	a[I1] =   inpk<L1,M1,K2,0,0,K2/2,K2/3>::k(b,c); }
      tX sv job(X*a,cX*b,cX*c,cX&d) {
	a[I1] = d*inpk<L1,M1,K2,0,0,K2/2,K2/3>::k(b,c); }
    };
    //
    tm<int K2,int K1,int L1,int M1> struct inp_ad {
      sci I1 = TWO<L1,M1>::I;
      tX sv job(X*a,cX*b,cX*c) {
	a[I1]+=   inpk<L1,M1,K2,0,0,K2/2,K2/3>::k(b,c); }
      tX sv job(X*a,cX*b,cX*c,cX&d) {
	a[I1]+= d*inpk<L1,M1,K2,0,0,K2/2,K2/3>::k(b,c); }
    };
    //
    tm<int K2,int K1,int L1,int M1> struct inp_su {
      sci I1 = TWO<L1,M1>::I;
      tX sv job(X*a,cX*b,cX*c) {
	a[I1]-=   inpk<L1,M1,K2,0,0,K2/2,K2/3>::k(b,c); }
      tX sv job(X*a,cX*b,cX*c,cX&d) {
	a[I1]-= d*inpk<L1,M1,K2,0,0,K2/2,K2/3>::k(b,c); }
    };
    //
    // struct inp<>: perform complete inner product
    //
    tm<int K1,int K2> struct inp {
      tX sv as( X*a,cX*b,cX*c)      { LoopLM<inp_as,K1,K2>::loop(a,b,c); }
      tX sv ad( X*a,cX*b,cX*c)      { LoopLM<inp_ad,K1,K2>::loop(a,b,c); }
      tX sv su( X*a,cX*b,cX*c)      { LoopLM<inp_su,K1,K2>::loop(a,b,c); }
      tX sv as( X*a,cX*b,cX*c,cX&d) { LoopLM<inp_as,K1,K2>::loop(a,b,c,d); }
      tX sv ad( X*a,cX*b,cX*c,cX&d) { LoopLM<inp_ad,K1,K2>::loop(a,b,c,d); }
      tX sv su( X*a,cX*b,cX*c,cX&d) { LoopLM<inp_su,K1,K2>::loop(a,b,c,d); }
    };
    //
    tm<int K2> struct inp<0,K2> {
      tX sX ras (cX*b,cX*c) {
	return inpk<0,0,K2,0,0,K2/2,K2/3>::k(b,c);
      }
    };
  } // namespace meta3D
  //
  // implementing
  //
  // falcON::symt3D<K,X>::ass_inn_prd();
  // falcON::symt3D<K,X>::add_inn_prd();
  // falcON::symt3D<K,X>::sub_inn_prd();
  //
  tKX inline void symt3D<K,X>::ass_inn_prd(V const&p, symt3D<K-1,X>&q) const {
    meta3D::inp<K-1,1>::as(pX(q),a,pX(p)); }
  tKX inline void symt3D<K,X>::add_inn_prd(V const&p, symt3D<K-1,X>&q) const {
    meta3D::inp<K-1,1>::ad(pX(q),a,pX(p)); }
  tKX inline void symt3D<K,X>::sub_inn_prd(V const&p, symt3D<K-1,X>&q) const {
    meta3D::inp<K-1,1>::su(pX(q),a,pX(p)); }

  tKX tQ inline void symt3D<K,X>::ass_inn_prd(symt3D<Q,X>  const&p,
					      symt3D<K-Q,X>     &q) const {
    meta3D::inp<K-Q,Q >::as(pX(q),a,pX(p)); }
  tKX tQ inline void symt3D<K,X>::add_inn_prd(symt3D<Q,X>  const&p,
					      symt3D<K-Q,X>     &q) const {
    meta3D::inp<K-Q,Q >::ad(pX(q),a,pX(p)); }
  tKX tQ inline void symt3D<K,X>::sub_inn_prd(symt3D<Q,X>  const&p,
					      symt3D<K-Q,X>     &q) const {
    meta3D::inp<K-Q,Q >::su(pX(q),a,pX(p)); }

  tKX inline void symt3D<K,X>::ass_inn_prd(symt3D<K-1,X>const&p, V&q) const {
    meta3D::inp<1,K-1>::as(pX(q),a,pX(p)); }
  tKX inline void symt3D<K,X>::add_inn_prd(symt3D<K-1,X>const&p, V&q) const {
    meta3D::inp<1,K-1>::ad(pX(q),a,pX(p)); }
  tKX inline void symt3D<K,X>::sub_inn_prd(symt3D<K-1,X>const&p, V&q) const {
    meta3D::inp<1,K-1>::su(pX(q),a,pX(p)); }

  tKX inline void symt3D<K,X>::ass_inn_prd(symt3D<K,X>  const&p, X&q) const {
    q = meta3D::inp<0,K>::ras(a,pX(p)); }
  tKX inline void symt3D<K,X>::add_inn_prd(symt3D<K,X>  const&p, X&q) const {
    q+= meta3D::inp<0,K>::ras(a,pX(p)); }
  tKX inline void symt3D<K,X>::sub_inn_prd(symt3D<K,X>  const&p, X&q) const {
    q-= meta3D::inp<0,K>::ras(a,pX(p)); }

  tKX inline X    symt3D<K,X>::ass_inn_prd(symt3D<K,X>  const&p) const {
    meta3D::inp<0,K>::ras(a,pX(p)); }
  //
  tX inline void symt3D<2,X>::ass_inn_prd(V const&p, V&q) const {
    meta3D::inp<1,1>::as(pX(q),a,pX(p)); }
  tX inline void symt3D<2,X>::add_inn_prd(V const&p, V&q) const {
    meta3D::inp<1,1>::ad(pX(q),a,pX(p)); }
  tX inline void symt3D<2,X>::sub_inn_prd(V const&p, V&q) const {
    meta3D::inp<1,1>::su(pX(q),a,pX(p)); }

  tX inline void symt3D<2,X>::ass_inn_prd(symt3D<2,X> const&p, X&q) const {
    q = meta3D::inp<0,2>::ras(a,pX(p)); }
  tX inline void symt3D<2,X>::add_inn_prd(symt3D<2,X> const&p, X&q) const {
    q+= meta3D::inp<0,2>::ras(a,pX(p)); }
  tX inline void symt3D<2,X>::sub_inn_prd(symt3D<2,X> const&p, X&q) const {
    q-= meta3D::inp<0,2>::ras(a,pX(p)); }

  tX inline X    symt3D<2,X>::ass_inn_prd(symt3D<2,X> const&p) const {
    meta3D::inp<0,2>::ras(a,pX(p)); }
} // namespace falcON
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
#endif // falcON_included_tensor_cc
