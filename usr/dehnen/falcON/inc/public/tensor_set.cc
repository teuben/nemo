// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// tensor_set.cc                                                               |
//                                                                             |
// Copyright (C) 2003-2005  Walter Dehnen                                      |
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
// implements tensor_set.h via (inlined) meta template programming methods     |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_tensor_set_cc
#define falcON_included_tensor_set_cc

#ifndef falcON_included_tensor_set_h
#  include <public/tensor_set.h>
#endif
////////////////////////////////////////////////////////////////////////////////
#define tm    template
#define tNX   template<int N,typename X>
#define sv    static void
#define si    static int
#define sci   static const int
#define scb   static const bool
#define sX    static X
#define cX    const  X
#define cA    const  A
#define cB    const  B
#define cC    const  C
#define cD    const  D
#define tX    template<typename X>
#define tA    template<typename A>
#define tAB   template<typename A, typename B>
#define tABC  template<typename A, typename B, typename C>
#define tABCD template<typename A, typename B, typename C, typename D>
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// macros to be used for minimum and maximum of const integer expressions     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#undef  __min
#undef  __max
#define __min(A,B) ((A)<(B)? (A):(B))
#define __max(A,B) ((A)>(B)? (A):(B))
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  namespace meta3D {
    using meta::taux;
    ////////////////////////////////////////////////////////////////////////////
    const double iN[16] = { 1.,                    // 1/N                       
			    1.,
			    0.5,
			    0.3333333333333333333333333,
			    0.25,
			    0.2,
			    0.1666666666666666666666667,
			    0.1428571428571428571428571,
			    0.125,
			    0.1111111111111111111111111,
			    0.1,
			    0.09090909090909090909090909,
			    0.08333333333333333333333333,
			    0.076923076923076923076923076,
			    0.071428571428571428571428571,
			    0.066666666666666666666666667 };
    ////////////////////////////////////////////////////////////////////////////
    const double iF[16] = { 1.,                    // 1/N!                      
			    1.,
			    0.5,
			    0.166666666666666666666667,
			    0.416666666666666666666667e-1,
			    0.833333333333333333333333e-2,
			    0.138888888888888888888889e-2,
			    0.198412698412698412698413e-3,
			    0.248015873015873015873016e-4,
			    0.275573192239858906525573e-5,
			    0.275573192239858906525573e-6,
			    0.250521083854417187750521e-7,
			    0.208767569878680989792101e-8,
			    0.160590438368216145993924e-9,
			    0.114707455977297247138517e-10,
			    0.764716373181981647590113e-12 };
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // struct LoopK<Perform,N,Z,K0=0>                                         //
    //                                                                        //
    // loop K=K0..N  and  perform some templated action                       //
    //                                                                        //
    // template arguments of P: Z(anything),K                                 //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
#define T2i template<int,int> class
    tm<T2i P, int Z, int K, int N> struct _k {
  tA    sv loop(A&a) {
    P<Z,K>::job(a);       _k<P,Z,K+1,N>::loop(a); }
  tAB   sv loop(A&a,cB&b) {
    P<Z,K>::job(a,b);     _k<P,Z,K+1,N>::loop(a,b); }
  tABC  sv loop(A&a,cB&b,cC&c) {
    P<Z,K>::job(a,b,c);   _k<P,Z,K+1,N>::loop(a,b,c); }
  tABCD sv loop(A&a,cB&b,cC&c,cD&d) { 
    P<Z,K>::job(a,b,c,d); _k<P,Z,K+1,N>::loop(a,b,c,d); }
};
    //--------------------------------------------------------------------------
    tm<T2i P, int Z, int K> struct _k<P,Z,K,K> {
      tA    sv loop(A&a)                { P<Z,K>::job(a); }
      tAB   sv loop(A&a,cB&b)           { P<Z,K>::job(a,b); }
      tABC  sv loop(A&a,cB&b,cC&c)      { P<Z,K>::job(a,b,c); }
      tABCD sv loop(A&a,cB&b,cC&c,cD&d) { P<Z,K>::job(a,b,c,d); }
    };
    //--------------------------------------------------------------------------
    tm<T2i P, int N, int Z, int K0=0> struct LoopK {
      tA    sv loop(A&a)                { _k<P,Z,K0,N>::loop(a); }
      tAB   sv loop(A&a,cB&b)           { _k<P,Z,K0,N>::loop(a,b); }
      tABC  sv loop(A&a,cB&b,cC&c)      { _k<P,Z,K0,N>::loop(a,b,c); }
      tABCD sv loop(A&a,cB&b,cC&c,cD&d) { _k<P,Z,K0,N>::loop(a,b,c,d); }
    };
#undef T2i
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // struct LoopKLM<Perform,N,Z,K0=0>                                       //
    //                                                                        //
    // loop K=K0..N, L=0..K, M=0..L  and  perform some templated action       //
    //                                                                        //
    // template arguments of P: Z(integer),K,L,M                              //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
#define T4i template<int,int,int,int> class
    tm<T4i P, int Z, int K, int N> struct k {
  tA    sv loop(A&a) { 
    l<P,Z,K,0,0>::loop(a);
    k<P,Z,K+1,N>::loop(a);
  }
  tAB   sv loop(A&a,cB&b) { 
    l<P,Z,K,0,0>::loop(a,b);
    k<P,Z,K+1,N>::loop(a,b);
  }
  tABC  sv loop(A&a,cB&b,cC&c) { 
    l<P,Z,K,0,0>::loop(a,b,c);
    k<P,Z,K+1,N>::loop(a,b,c);
  }
  tABCD sv loop(A&a,cB&b,cC&c,cD&d) { 
    l<P,Z,K,0,0>::loop(a,b,c,d);
    k<P,Z,K+1,N>::loop(a,b,c,d);
  }
};
    //--------------------------------------------------------------------------
    tm<T4i P, int Z, int K> struct k<P,Z,K,K> {
      tA    sv loop(A&a)                { l<P,Z,K,0,0>::loop(a); }
      tAB   sv loop(A&a,cB&b)           { l<P,Z,K,0,0>::loop(a,b); }
      tABC  sv loop(A&a,cB&b,cC&c)      { l<P,Z,K,0,0>::loop(a,b,c); }
      tABCD sv loop(A&a,cB&b,cC&c,cD&d) { l<P,Z,K,0,0>::loop(a,b,c,d); }
    };
    //--------------------------------------------------------------------------
    tm<T4i P, int N, int Z, int K0=0> struct LoopKLM {
      tA    sv loop(A&a)                { k<P,Z,K0,N>::loop(a); }
      tAB   sv loop(A&a,cB&b)           { k<P,Z,K0,N>::loop(a,b); }
      tABC  sv loop(A&a,cB&b,cC&c)      { k<P,Z,K0,N>::loop(a,b,c); }
      tABCD sv loop(A&a,cB&b,cC&c,cD&d) { k<P,Z,K0,N>::loop(a,b,c,d); }
    };
#undef T4i
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // implementing                                                           //
    //                                                                        //
    // falcON::poles3D<N,X>::ass_out_prd();                                   //
    // falcON::symset3D<N,X>::ass_out_prd();                                  //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    tm<int Z,int K> struct o {
      tX sv job(X*a, cX*v) {
	vop<K>::as(a+ONE3D<K>::I0,a+ONE3D<K-1>::I0,v);
      } };
    tm<int Z> struct o<Z,0> {
      tX sv job(X*a, cX*v) {
	a[0]=X(1);
      } };
    tm<int Z> struct o<Z,1> {
      tX sv job(X*a, cX*v) {
	a[1]=v[0];
	a[2]=v[1];
	a[3]=v[2];
      } };
    tm<int Z> struct o<Z,2> {
      tX sv job(X*a, cX*v) {
	a[4]=v[0]*v[0];
	a[5]=v[1]*v[0];
	a[6]=v[2]*v[0];
	a[7]=v[1]*v[1];
	a[8]=v[2]*v[1];
	a[9]=v[2]*v[2];
      } };
    ////////////////////////////////////////////////////////////////////////////
    tm<int N, int K0=0> struct op {
      tX sv as(X*a, cX*v) { LoopK<o,N,0,K0>::loop(a,v); }
    };
  } // namespace meta3D {
  //////////////////////////////////////////////////////////////////////////////
  tNX inline symset3D<N,X>& symset3D<N,X>::ass_out_prd(tupel<3,X> const&v) {
    meta3D::op<N>::as(a,pX(v)); return*this; }
  tNX inline poles3D<N,X>& poles3D<N,X>::ass_out_prd(tupel<3,X> const&v) {
    meta3D::op<N,2>::as(a-4,pX(v)); return*this; }
  //////////////////////////////////////////////////////////////////////////////
  namespace meta3D {
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // implementing                                                           //
    //                                                                        //
    // falcON::poles3D<N,X>::add_body();                                      //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // in place m times outer products of vector for K>=2 only                //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    tm<int K> struct ipop {                        // advances from K-1 to K    
      sci KK = ((K+1)*(K+2))/2;
      sci K0 = ( K   *(K+1))/2;
      sci K1 = ((K-1)* K   )/2;
      tX sv as(X*a, cX*v) {
	a[KK-1] = a[K0-1] * v[2];
	taux<X,K -1>::v_ast(a+K0,a+K1,v[1]);
	taux<X,K0-1>::s_ml (a,v[0]);
      } };
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // add m times outer products of vector for K>=2 only                     //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    tm<int K, int N> struct po {
      sci KK = ((K+1)*(K+2))/2;
      tX sv ad(X*a, X*b, cX*v) {                   // a += (v op b)   & cont    
	ipop<K>::as(b,v);
	taux<X,KK-1>::v_ad(a+ONE3D<K>::J0, b);
	po<K+1,N>::ad(a,b,v);
      } };
    tm<int K> struct po<K,K> {
      sci KK = ((K+1)*(K+2))/2;
      tX sv ad(X*a, X*b, cX*v) {                   // a += (v op b)             
	ipop<K>::as(b,v);
	taux<X,KK-1>::v_ad(a+ONE3D<K>::J0, b);
      } };
    //--------------------------------------------------------------------------
    tm<int N> struct po<2,N> {                     // a+=b=m*(v op v) & cont    
      tX sv ad(X*a, X*b, cX*v, cX&m) {
	register X tmp  = m * v[0];
	a[0] += b[0] = v[0] * tmp;
	a[1] += b[1] = v[1] * tmp;
	a[2] += b[2] = v[2] * tmp;
	tmp   = m * v[1];
	a[3] += b[3] = v[1] * tmp;
	a[4] += b[4] = v[2] * tmp;
	a[5] += b[5] = m * v[2] * v[2];
	po<3,N>::ad(a,b,v);
      } };
    tm<> struct po<2,2> {                          // a+=m*(v op v)  DONE       
      tX sv ad(X*a, cX*v, cX&m) {
	register X tmp  = m * v[0];
	a[0] += v[0] * tmp;
	a[1] += v[1] * tmp;
	a[2] += v[2] * tmp;
	tmp   = m * v[1];
	a[3] += v[1] * tmp;
	a[4] += v[2] * tmp;
	a[5] += m * v[2] * v[2];
      } };
    ////////////////////////////////////////////////////////////////////////////
    tm<int N> struct pol {
      tX sv ad(X*a, cX*v, cX&m) {
	X b[ONE3D<N>::ND];
	po<2,N>::ad(a,b,v,m);
      } };
    tm<> struct pol<2> {
      tX sv ad(X*a, cX*v, cX&m) {
	po<2,2>::ad(a,v,m);
      } };
  } // namespace meta3D {
  //////////////////////////////////////////////////////////////////////////////
  tNX inline poles3D<N,X>&
  poles3D<N,X>::add_body(tupel<3,X> const&v, X const&m) {
    meta3D::pol<N>::ad(a,pX(v),m); return*this; }
  //////////////////////////////////////////////////////////////////////////////
  namespace meta3D {
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // implementing                                                           //
    //                                                                        //
    // falcON::poles3D<N,X>::add_cell();                                      //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // we have to compute M_k += x^(l) o P_(k-l)  for l=0,...,N               //
    // note that P_1 vanishes, P_0 is the mass m                              //
    //                                                                        //
    // we do this in 3 steps:                                                 //
    //                                                                        //
    // 1.     M_k += x^(k) * P_0   (P_0 = mass) via add_body();               //
    // 2.     M_k += P_k           for all k  is simple.                      //
    // 3.     loop l=1...N-2                                                  //
    //            compute x^(l) in place in memory pointed to by 'b'          //
    //            loop k=l+2...N                                              //
    //                M_k += x^(l) o P_(k-l)                                  //
    //                                                                        //
    // steps 2 and 3 are done by __cell_l<N,L>                                //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    tm<int L, int K, int N> struct __cell_k {      // general case              
      tX sv ad(X*a, cX*b, cX*c) {                  //   M_K += x^(L) o P_(K-L)  
	syp<L,K-L>::ad(a+ONE3D<K>::J0,b,c+ONE3D<K-L>::J0);
	__cell_k<L,K+1,N>::ad(a,b,c);
      } };
    tm<int L, int K> struct __cell_k<L,K,K> {      // terminator: K=N           
      tX sv ad(X*a, cX*b, cX*c) {                  //   M_K += x^(L) o P_(K-L)  
	syp<L,K-L>::ad(a+ONE3D<K>::J0,b,c+ONE3D<K-L>::J0);
      } };
    ////////////////////////////////////////////////////////////////////////////
    tm<int N, int L, int HL> struct __cell_l {     // general case HL==N-2      
      tX sv ad(X*a, X*b, cX*v, cX*c) {
	ipop<L>::as(b,v);                          //   compute x^(L) in place  
	__cell_k<L,L+2,N>::ad(a,b,c);              //   loop K=L+2..N { ... }   
	__cell_l<N,L+1,HL>::ad(a,b,v,c);           //   go ahead with L+1       
      } };
    tm<int N, int L> struct __cell_l<N,L,L> {      // terminator L=HL==N-2      
      tX sv ad(X*a, X*b, cX*v, cX*c) {
	ipop<L>::as(b,v);                          //   compute x^(L) in place  
	__cell_k<L,N,N>::ad(a,b,c);                //   K=L+2=N (no loop)       
      } };
    tm<int N, int HL> struct __cell_l<N,0,HL> {    // special case L=0; N       
      tX sv ad(X*a, cX*v, cX*c) {
	taux<X,ONE3D<N>::PD-1>::v_ad(a,c);         //   add all poles           
	__cell_l<N,1,HL>::ad(a,v,c);               //   go ahead with L=1       
      } };
    tm<int N, int HL> struct __cell_l<N,1,HL> {    // special case L=1; N       
      tX sv ad(X*a, cX*v, cX*c) {
	X b[ONE3D<N>::ND];                         //   make memory for b       
	b[0] = v[0]; b[1]=v[1]; b[2]=v[2];         //   set b to v              
	__cell_k<1,3,N>::ad(a,v,c);                //   loop K=3..N { ... }     
	__cell_l<N,2,HL>::ad(a,b,v,c);             //   go ahead with l=2       
      } };
    tm<> struct __cell_l<2,0,0> {                  // special case L=0; N=2     
      tX sv ad(X*a, cX*, cX*c) {
	taux<X,ONE3D<2>::PD-1>::v_ad(a,c);         //   add all poles; DONE     
      } };
    tm<> struct __cell_l<3,1,1> {                  // special case L=1, N=3     
      tX sv ad(X*a, cX*v, cX*c) {
	__cell_k<1,3,3>::ad(a,v,c);                //   K=3 (no loop)           
      } };
    ////////////////////////////////////////////////////////////////////////////
    tm<int N> struct __cell {
      tX sv ad(X*a, cX*v, cX*c) {
	__cell_l<N,0,N-2>::ad(a,v,c);
      } };
  } // namespace meta3D {
  //////////////////////////////////////////////////////////////////////////////
  tNX inline poles3D<N,X>&
  poles3D<N,X>::add_cell(tupel<3,X> const&v, X const&m, poles3D<N,X>const&p) {
    meta3D::__cell<N>::ad(a,pX(v),pX(p));
    return add_body(v,m);
  }
  //////////////////////////////////////////////////////////////////////////////
  namespace meta3D {
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // implementing                                                           //
    //                                                                        //
    // falcON::poles3D<N,X>::normalize()                                      //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    tm<int K, int N> struct __norm {
      tX sv job(poles3D<N,X>&M, cX&ix) {
	M.template pole<K>() *= (ix*iF[K]);
	__norm<K+1,N>::job(M,ix);
      } };
    tm<int K> struct __norm<K,K> {
      tX sv job(poles3D<K,X>&M, cX&ix) {
	M.template pole<K>() *= (ix*iF[K]);
      } };
    tm<int N> struct __norm<2,N> {
      tX sv job(poles3D<N,X>&M, cX&x) {
	register X ix = X(1)/x;
	M.template pole<2>() *= (ix*iF[2]);
	__norm<3,N>::job(M,ix);
      } };
    tm<> struct __norm<2,2> {
      tX sv job(poles3D<2,X>&M, cX&x) {
	M.template pole<2>() *= iF[2]/x;
      } };
  } // namespace meta3D {
  //////////////////////////////////////////////////////////////////////////////
  tNX inline poles3D<N,X>& 
  poles3D<N,X>::normalize(X const&x) {
    meta3D::__norm<2,N>::job(*this,x);
    return*this;
  }
  //////////////////////////////////////////////////////////////////////////////
  namespace meta3D {
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // implementing                                                           //
    //                                                                        //
    // falcON::symset3D<N,X>::flip_sign_odd()                                 //
    // falcON::symset3D<N,X>::flip_sign_even()                                //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    tm<int k, int n, int N> struct __flip_o {
      sci K=k+k+1;
      tX sv job(symset3D<N,X>&s) {
	s.template tensor<K>().negate();
	__flip_o<k+1,n,N>::job(s);
      } };
    tm<int k, int N> struct __flip_o<k,k,N> {
      sci K=k+k+1;
      tX sv job(symset3D<N,X>&s) {
	s.template tensor<K>().negate();
      } };
    ////////////////////////////////////////////////////////////////////////////
    tm<int k, int n, int N> struct __flip_e {
      sci K=k+k;
      tX sv job(symset3D<N,X>&s) {
	s.template tensor<K>().negate();
	__flip_e<k+1,n,N>::job(s);
      } };
    tm<int k, int N> struct __flip_e<k,k,N> {
      sci K=k+k;
      tX sv job(symset3D<N,X>&s) {
	s.template tensor<K>().negate();
      } };
  } // namespace meta3D {
  //////////////////////////////////////////////////////////////////////////////
  tNX inline symset3D<N,X>&
  symset3D<N,X>::flip_sign_odd () { 
    meta3D::__flip_o<0,(N-1)/2,N>::job(*this);
    return*this;
  }
  tNX inline symset3D<N,X>&
  symset3D<N,X>::flip_sign_even() {
    meta3D::__flip_e<0,N/2,N>::job(*this);
    return*this;
  }
  //////////////////////////////////////////////////////////////////////////////
  namespace meta3D {
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // implementing                                                           //
    //                                                                        //
    // falcON::symset3D<>::ass_inn_prd();                                     //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    tm<int K, int U> struct __inn_prd_set {
      tX sv as(cX*r, cX*v, X*l) {
	inp<K,1>::as(l+ONE3D<K>::I0, r+ONE3D<K+1>::I0 ,v);
	__inn_prd_set<K+1,U>::as(r,v,l);
      } };
    tm<int K> struct __inn_prd_set<K,K> {
      tX sv as(cX*r, cX*v, X*l) {
	inp<K,1>::as(l+ONE3D<K>::I0, r+ONE3D<K+1>::I0 ,v);
      } };
    tm<int U> struct __inn_prd_set<0,U> {
      tX sv as(cX*r, cX*v, X*l) {
	l[0] = v[0]*r[1] + v[1]*r[2] + v[2]*r[3];
	__inn_prd_set<1,U>::as(r,v,l);
      } };
    tm<> struct __inn_prd_set<0,0> {
      tX sv as(cX*r, cX*v, X*l) {
	l[0] = v[0]*r[1] + v[1]*r[2] + v[2]*r[3];
      } };
    //--------------------------------------------------------------------------
    tm<int N> struct __inn_prd_Set {
      tX sv as(cX*r, cX*v, X*l) {
	__inn_prd_set<0,N-1>::as(r,v,l);
      } };
  } // namespace meta3D {
  //////////////////////////////////////////////////////////////////////////////
  tNX inline void
  symset3D<N,X>::ass_inn_prd(tupel<3,X> const&v, symset3D<N-1,X>&l) const {
    meta3D::__inn_prd_Set<N>::as(a,pX(v),pX(l));
  }
  //////////////////////////////////////////////////////////////////////////////
  namespace meta3D {
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // implementing                                                           //
    //                                                                        //
    // falcON::shift_by();                                                    //
    // falcON::eval_expn();                                                   //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    tm<int N, int K, int I> struct __shift {
      tX sv at(X*n, X*o, cX*v) {
	n[0]-=o[0];
	n[1]+=o[1];  n[2]+=o[2];  n[3]+=o[3];
	X __v[3];
	taux<X,2>::v_ast(__v,v,iN[I+1]);
	__inn_prd_Set<K>::as(o,__v,o);
	__shift<N,K-1,I+1>::at(n,o,v);
      }
      tX sv by(X*n, X*o, cX*v) {
	taux<X,ONE3D<K>::CD-1>::v_ad(n,o);
	X __v[3];
	taux<X,2>::v_ast(__v,v,iN[I+1]);
	__inn_prd_Set<K>::as(o,__v,o);
	__shift<N,K-1,I+1>::by(n,o,v);
      } };
    tm<int N, int I> struct __shift<N,1,I> {
      tX sv at(X*n, cX*o, cX*v) {
	n[0]-=o[0] + iN[N] * (v[0]*o[1]+v[1]*o[2]+v[2]*o[3]);
	n[1]+=o[1];  n[2]+=o[2];  n[3]+=o[3];
      }
      tX sv by(X*n, cX*o, cX*v) {
	n[0]+=o[0] + iN[N] * (v[0]*o[1]+v[1]*o[2]+v[2]*o[3]);
	n[1]+=o[1];  n[2]+=o[2];  n[3]+=o[3];
      } };
    tm<int N> struct __shift<N,N,0> {
      tX sv at(X*n, cX*c, cX*v) {
	n[0]-=c[0];
	n[1]+=c[1];  n[2]+=c[2];  n[3]+=c[3];
	X o[ONE3D<N-1>::CD];
	__inn_prd_Set<N>::as(c,v,o);
	__shift<N,N-1,1>::at(n,o,v);
      }
      tX sv by(X*n, cX*v) {
	X o[ONE3D<N-1>::CD];
	__inn_prd_Set<N>::as(n,v,o);
	__shift<N,N-1,1>::by(n,o,v);
      } };
  } // namespace meta3D {
  //////////////////////////////////////////////////////////////////////////////
  tNX inline void
  shift_by(symset3D<N,X>&C, tupel<3,X>&x) {
    meta3D::__shift<N,N,0>::by(static_cast<X*>(C), static_cast<X*>(x));
  }
  tNX inline void
  eval_expn(symset3D<1,X>&G, symset3D<N,X> const&C, tupel<3,X> const&x) {
    meta3D::__shift<N,N,0>::at(static_cast<X*>(G),
			       static_cast<const X*>(C),
			       static_cast<const X*>(x));
  }
  //////////////////////////////////////////////////////////////////////////////
  namespace meta3D {
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // implementing                                                           //
    //                                                                        //
    // falcON::set_dPhi();                                                    //
    // falcON::add_dPhi();                                                    //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // C_n = Sum(R^(n-2*k) o D^(k) * (-1)^(n-k) * d_[n-k], k=0...[n/2])       //
    //                                                                        //
    // we do this as follows:                                                 //
    //                                                                        //
    // - set C_n = R^(n)                                                      //
    // - loop k=N...0                                                         //
    //      - loop l=1...(N-k)/2                                              //
    //           - n = N+2*L                                                  //
    //           - C_n += C_k o D^l * (-1)^(k+l) * d[k+l]                     //
    //      - C_k *= (-1)^k * d[k]                                            //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    tm<int N> struct __sign_times {
      tX sX s(X const&x) { return -x; } };
    tm<> struct __sign_times<0> {
      tX sX s(X const&x) { return  x; } };
    ////////////////////////////////////////////////////////////////////////////
    // sign_times<K>(x); return (-1)^K * x                                      
    template<int K, typename X>
    X sign_times(X const&x) { return __sign_times<K%2>::s(x); }
    ////////////////////////////////////////////////////////////////////////////
#ifdef falcON_SSE_CODE
#  define D_ARG      const fvec4*const&D, int const&J
#  define ARG_D      D,J
#  define __D(I,J)   D[I][J]
#  define __SD(I,J)  sign_times<I>(D[I][J])
#else
#  define D_ARG      const X*const&D
#  define ARG_D      D
#  define __D(I,J)   D[I]
#  define __SD(I,J)  sign_times<I>(D[I])
#endif
    ////////////////////////////////////////////////////////////////////////////
    tm<int ORDER,int K,int L,int U> struct __grav_l {
      sci N=K+L+L;
      tX sv ad(symset3D<ORDER,X>&C, D_ARG) {
	C.template tensor<N>().add_del_prd(C.template tensor<K>(),__SD(K+L,J));
	__grav_l<ORDER,K,L+1,U>::ad(C,ARG_D);
      } };
    tm<int ORDER,int K,int L> struct __grav_l<ORDER,K,L,L> {
      sci N=K+L+L;
      tX sv ad(symset3D<ORDER,X>&C, D_ARG) {
	C.template tensor<N>().add_del_prd(C.template tensor<K>(),__SD(K+L,J));
      } };
    tm<int ORDER,int L,int U> struct __grav_l<ORDER,0,L,U> {
      sci N=L+L;
      tX sv ad(symset3D<ORDER,X>&C, D_ARG) {
	C.template tensor<N>().add_del_prd(__SD(L,J));
	__grav_l<ORDER,0,L+1,U>::ad(C,ARG_D);
      } };
    tm<int ORDER,int L> struct __grav_l<ORDER,0,L,L> {
      sci N=L+L;
      tX sv ad(symset3D<ORDER,X>&C, D_ARG) {
	C.template tensor<N>().add_del_prd(__SD(L,J));
      } };
    ////////////////////////////////////////////////////////////////////////////
    tm<int,int,int,bool> struct __grav_dol;
    tm<int ORDER,int K,int U> struct __grav_dol<ORDER,K,U,1> {
      tX sv ad(symset3D<ORDER,X>&C, D_ARG) {
	__grav_l<ORDER,K,1,U>::ad(C,ARG_D);
      } };
    tm<int ORDER,int K,int U> struct __grav_dol<ORDER,K,U,0> {
      tX sv ad(symset3D<ORDER,X>&, D_ARG) {}
    };
    ////////////////////////////////////////////////////////////////////////////
    tm<int ORDER,int K> struct __grav_k {
      sci U = (ORDER-K)/2;
      scb Q = 1 <= U;
      tX sv ad(symset3D<ORDER,X>&C, D_ARG) {
	__grav_dol<ORDER,K,U,Q>::ad(C,ARG_D);
	C.template tensor<K>() *= __SD(K,J);
	__grav_k<ORDER,K-1>::ad(C,ARG_D);
      } };
    tm<int ORDER> struct __grav_k<ORDER,0> {
      sci U = ORDER/2;
      scb Q = 1 <= U;
      tX sv ad(symset3D<ORDER,X>&C, D_ARG) {
	__grav_dol<ORDER,0,U,Q>::ad(C,ARG_D);
	C.template tensor<0>() = __D(0,J);
      }
    };
    ////////////////////////////////////////////////////////////////////////////
    tm<int ORDER> struct __grav {
      tX sv ass(symset3D<ORDER,X>&C, const X*v, D_ARG) {
	C.ass_out_prd(v);
	__grav_k<ORDER,ORDER>::ad(C,ARG_D);
      }
      tX sv add(symset3D<ORDER,X>&C, const X*v, D_ARG) {
	symset3D<ORDER,X>&F;
	ass(F,v,ARG_D);
	C += F;
      }
    };
    ////////////////////////////////////////////////////////////////////////////
    // special treatment for small orders                                       
    ////////////////////////////////////////////////////////////////////////////
    tm<> struct __grav<1> {
      tX sv __ass(X*a, const X*v, D_ARG) {
	a[ 0] = __D(0,J);
	register X t =-__D(1,J);
	a[ 1] = v[0]*t;
	a[ 2] = v[1]*t;
	a[ 3] = v[2]*t;
      }
      tX sv ass(symset3D<1,X>&C, const X*v, D_ARG) { __ass((X*)(C),v,ARG_D); }
      tX sv __add(X*a, const X*v, D_ARG) {
	a[ 0]+= __D(0,J);
	register X t =-__D(1,J);
	a[ 1]+= v[0]*t;
	a[ 2]+= v[1]*t;
	a[ 3]+= v[2]*t;
      }
      tX sv add(symset3D<1,X>&C, const X*v, D_ARG) { __add((X*)(C),v,ARG_D); }
    };
    tm<> struct __grav<2> {
      tX sv __ass(X*a, const X*v, D_ARG) {
	a[ 0] = __D(0,J);
	register X t =-__D(1,J);
	a[ 1] = v[0]*t;
	a[ 2] = v[1]*t;
	a[ 3] = v[2]*t;
	t     = __D(2,J)*v[0];
	a[ 4] = t*v[0] - __D(1,J);
	a[ 5] = t*v[1];
	a[ 6] = t*v[2];
	t     = __D(2,J)*v[1];
	a[ 7] = t*v[1] - __D(1,J);
	a[ 8] = t*v[2];
	t     = __D(2,J)*v[2];
	a[ 9] = t*v[2] - __D(1,J);
      }
      tX sv ass(symset3D<2,X>&C, const X*v, D_ARG) { __ass((X*)(C),v,ARG_D); }
      tX sv __add(X*a, const X*v, D_ARG) {
	a[ 0]+= __D(0,J);
	register X t =-__D(1,J);
	a[ 1]+= v[0]*t;
	a[ 2]+= v[1]*t;
	a[ 3]+= v[2]*t;
	t     = __D(2,J)*v[0];
	a[ 4]+= t*v[0] - __D(1,J);
	a[ 5]+= t*v[1];
	a[ 6]+= t*v[2];
	t     = __D(2,J)*v[1];
	a[ 7]+= t*v[1] - __D(1,J);
	a[ 8]+= t*v[2];
	t     = __D(2,J)*v[2];
	a[ 9]+= t*v[2] - __D(1,J);
      }
      tX sv add(symset3D<2,X>&C, const X*v, D_ARG) { __add((X*)(C),v,ARG_D); }
    };
    tm<> struct __grav<3> {
      tX sv __ass(X*a, const X*v, D_ARG) {
	a[ 0] = __D(0,J);
	register X t =-__D(1,J);
	a[ 1] = v[0]*t;
	a[ 2] = v[1]*t;
	a[ 3] = v[2]*t;
	t     = __D(2,J)*v[0];
	a[ 4] = t*v[0] - __D(1,J);
	a[ 5] = t*v[1];
	a[ 6] = t*v[2];
	t     = __D(2,J)*v[1];
	a[ 7] = t*v[1] - __D(1,J);
	a[ 8] = t*v[2];
	t     = __D(2,J)*v[2];
	a[ 9] = t*v[2] - __D(1,J);
	register X D2_3 = times<3>(__D(2,J));
	t     = __D(3,J)*v[0]*v[0];
	a[10] = (D2_3-t)*v[0];
	a[11] = (__D(2,J)-t)*v[1];
	a[12] = (__D(2,J)-t)*v[2];
	t     =  __D(3,J)*v[1]*v[1];
	a[13] = (__D(2,J)-t)*v[0];
	a[16] = (D2_3-t)*v[1];
	a[17] = (__D(2,J)-t)*v[2];
	t     =  __D(3,J)*v[2]*v[2];
	a[15] = (__D(2,J)-t)*v[0];
	a[18] = (__D(2,J)-t)*v[1];
	a[19] = (D2_3-t)*v[2];
	a[14] =-__D(3,J)*v[0]*v[1]*v[2];
      }
      tX sv ass(symset3D<3,X>&C, const X*v, D_ARG) { __ass((X*)(C),v,ARG_D); }
      tX sv __add(X*a, const X*v, D_ARG) {
	a[ 0]+= __D(0,J);
	register X t = __D(1,J);
	a[ 1]-= v[0]*t;
	a[ 2]-= v[1]*t;
	a[ 3]-= v[2]*t;
	t     = __D(2,J)*v[0];
	a[ 4]+= t*v[0] - __D(1,J);
	a[ 5]+= t*v[1];
	a[ 6]+= t*v[2];
	t     = __D(2,J)*v[1];
	a[ 7]+= t*v[1] - __D(1,J);
	a[ 8]+= t*v[2];
	t     = __D(2,J)*v[2];
	a[ 9]+= t*v[2] - __D(1,J);
	register X D2_3 = times<3>(__D(2,J));
	t     = __D(3,J)*v[0]*v[0];
	a[10]+= (D2_3-t)*v[0];
	a[11]+= (__D(2,J)-t)*v[1];
	a[12]+= (__D(2,J)-t)*v[2];
	t     =  __D(3,J)*v[1]*v[1];
	a[13]+= (__D(2,J)-t)*v[0];
	a[16]+= (D2_3-t)*v[1];
	a[17]+= (__D(2,J)-t)*v[2];
	t     =  __D(3,J)*v[2]*v[2];
	a[15]+= (__D(2,J)-t)*v[0];
	a[18]+= (__D(2,J)-t)*v[1];
	a[19]+= (D2_3-t)*v[2];
	a[14]-= __D(3,J)*v[0]*v[1]*v[2];
      }
      tX sv add(symset3D<3,X>&C, const X*v, D_ARG) { __add((X*)(C),v,ARG_D); }
    };
#undef D_ARG
#undef ARG_D
#undef __D
#undef __SD
  } // namespace meta3D {
  //////////////////////////////////////////////////////////////////////////////
#ifdef falcON_SSE_CODE
  tNX inline void
  set_dPhi(symset3D<N,X>&C, tupel<3,X> const&v, const fvec4 D[N+1], int const&J)
  { meta3D::__grav<N>::ass(C,static_cast<const X*>(v),D,J); }
  tNX inline void
  add_dPhi(symset3D<N,X>&C, tupel<3,X> const&v, const fvec4 D[N+1], int const&J)
  { meta3D::__grav<N>::add(C,static_cast<const X*>(v),D,J); }
#else
  tNX inline void
  set_dPhi(symset3D<N,X>&C, tupel<3,X> const&v, X const D[N+1])
  { meta3D::__grav<N>::ass(C,static_cast<const X*>(v),D); }
  tNX inline void
  add_dPhi(symset3D<N,X>&C, tupel<3,X> const&v, X const D[N+1])
  { meta3D::__grav<N>::add(C,static_cast<const X*>(v),D); }
#endif
  //////////////////////////////////////////////////////////////////////////////
  namespace meta3D {
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // implementing                                                           //
    //                                                                        //
    // falcON::add_C_C2C();                                                   //
    // falcON::add_C_C2B();                                                   //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // C_k += F_k + M^2 . F^(k+2) - M^3 . F^(k+3) + ....                      //
    //                                                                        //
    // - loop k=0...N                                                         //
    //      C_k += F_k                                                        //
    // - loop m=2...ORDER-1                                                   //
    //      - loop k=0...min(N, ORDER-m)                                      //
    //           C_k += (-1)^m M^m . F^(k+m)                                  //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    tm<int,int,int,int,int,bool> struct __coef_k;
    //      loop k=0...KU : C_k += M^m . F^(k+m)                                
    tm<int ORDR, int N, int M, int K, int KU> struct __coef_k<ORDR,N,M,K,KU,1> {
      tX sv job(symset3D<N,X>&c, symset3D<ORDR,X> const&f, symt3D<M,X> const&p)
      {
	f.template tensor<K+M>().add_inn_prd(p, c.template tensor<K>());
	__coef_k<ORDR,N,M,K+1,KU,1>::job(c,f,p);
      } };
    tm<int ORDR, int N, int M, int K> struct __coef_k<ORDR,N,M,K,K,1> {
      tX sv job(symset3D<N,X>&c, symset3D<ORDR,X> const&f, symt3D<M,X> const&p)
      {
	f.template tensor<K+M>().add_inn_prd(p, c.template tensor<K>());
      } };
    //      loop k=0...KU : C_k -= M^m . F^(k+m)                                
    tm<int ORDR, int N, int M, int K, int KU> struct __coef_k<ORDR,N,M,K,KU,0> {
      tX sv job(symset3D<N,X>&c, symset3D<ORDR,X> const&f, symt3D<M,X> const&p)
      {
	f.template tensor<K+M>().sub_inn_prd(p, c.template tensor<K>());
	__coef_k<ORDR,N,M,K+1,KU,0>::job(c,f,p);
      } };
    tm<int ORDR, int N, int M, int K> struct __coef_k<ORDR,N,M,K,K,0> {
      tX sv job(symset3D<N,X>&c, symset3D<ORDR,X> const&f, symt3D<M,X> const&p)
      {
	f.template tensor<K+M>().sub_inn_prd(p, c.template tensor<K>());
      } };
    ////////////////////////////////////////////////////////////////////////////
    tm<int ORDR, int N, int M, int MU, bool S> struct __coef_m {
      sci KU = __min(N,ORDR-M), P=ORDR-1;
      tX sv job(symset3D<N,X>&c, symset3D<ORDR,X> const&f, poles3D<P,X> const&m)
      {
	__coef_k<ORDR,N,M,0,KU, S>::job(c,f,m.template pole<M>());
	__coef_m<ORDR,N,M+1,MU,!S>::job(c,f,m);
      } };
    tm<int ORDR, int N, int M, bool S> struct __coef_m<ORDR,N,M,M,S> {
      sci KU = __min(N,ORDR-M), P=ORDR-1;
      tX sv job(symset3D<N,X>&c, symset3D<ORDR,X> const&f, poles3D<P,X> const&m)
      {
	__coef_k<ORDR,N,M,0,KU, S>::job(c,f,m.template pole<M>());
      } };
    ////////////////////////////////////////////////////////////////////////////
    tm<int ORD> struct __c_cell {
      sci P=ORD-1;
      tX sv ad(symset3D<ORD,X>&c, symset3D<ORD,X> const&f, poles3D<P,X> const&m)
      {
	c += f;
	__coef_m<ORD,ORD,2,P,1>::job(c,f,m);
      } };
    tm<> struct __c_cell<2> {
      tX sv ad(symset3D<2,X>&c, symset3D<2,X> const&f) {
	c += f;
      } };
    ////////////////////////////////////////////////////////////////////////////
    tm<int ORD> struct __c_body {
      sci P=ORD-1;
      tX sv ad(symset3D<1,X>&c, symset3D<ORD,X> const&f, poles3D<P,X> const&m) {
	symset3D<1,X> cc=f.template subset<1>();
	__coef_m<ORD,1,2,P,1>::job(cc,f,m);
	c.template tensor<0>() -= cc.template tensor<0>();
	c.template tensor<1>() += cc.template tensor<1>();
      } };
    tm<> struct __c_body<2> {
      tX sv ad(symset3D<1,X>&c, symset3D<2,X> const&f) {
	c.template tensor<0>() -= f.template tensor<0>();
	c.template tensor<1>() += f.template tensor<1>();
      } };
  } // namespace meta3D {
  //////////////////////////////////////////////////////////////////////////////
  tNX inline void 
  add_C_C2C(symset3D<N,X>&c, symset3D<N,X> const&f, poles3D<N-1,X> const&m) {
    meta3D::__c_cell<N>::ad(c,f,m);
  }
  tNX inline void 
  add_C_C2B(symset3D<1,X>&c, symset3D<N,X> const&f, poles3D<N-1,X> const&m) {
    meta3D::__c_body<N>::ad(c,f,m);
  }
  //////////////////////////////////////////////////////////////////////////////
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#undef __min
#undef __max
#undef tm
#undef tNX
#undef sv
#undef si
#undef sci
#undef scb
#undef sX
#undef cX
#undef cA
#undef cB
#undef cC
#undef cD
#undef tX
#undef tA
#undef tAB
#undef tABC
#undef tABCD
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_tensor_set_cc
