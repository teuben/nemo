// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    inc/tupel.cc                                                       
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    1996-2005                                                          
///                                                                             
/// \brief   definition of some member of template class tupel<>                
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 1996-2005  Walter Dehnen                                       
//                                                                              
// This program is free software; you can redistribute it and/or modify         
// it under the terms of the GNU General Public License as published by         
// the Free Software Foundation; either version 2 of the License, or (at        
// your option) any later version.                                              
//                                                                              
// This program is distributed in the hope that it will be useful, but          
// WITHOUT ANY WARRANTY; without even the implied warranty of                   
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
// General Public License for more details.                                     
//                                                                              
// You should have received a copy of the GNU General Public License            
// along with this program; if not, write to the Free Software                  
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                    
//                                                                              
////////////////////////////////////////////////////////////////////////////////
#ifndef WDutils_included_tupel_cc
#define WDutils_included_tupel_cc

#ifndef WDutils_included_algorithm
#  include <algorithm>
#  define WDutils_included_algorithm
#endif
#ifndef WDutils_included_cmath
#  include <cmath>
#  define WDutils_included_cmath
#endif
////////////////////////////////////////////////////////////////////////////////
#define c_      const
#define sv      static void
#define sb      static bool
#define sX      static X
#define cS      const S
#define cT      const T
#define tX      template<typename X>
#define tS      template<typename S>
#define tST     template<typename S, typename T>
#define tP      template<typename P>
#define tZS     template<int Z,typename S>
////////////////////////////////////////////////////////////////////////////////
namespace WDutils { namespace meta {
  //////////////////////////////////////////////////////////////////////////////
  using std::min;
  using std::max;
  using std::abs;
  using std::isinf;
  using std::isnan;
  //----------------------------------------------------------------------------
  template<int N> struct times__ {
    template<typename X> static X is(X const&x) { return N * x; } };
  template<> struct times__<2> {
    template<typename X> static X is(X const&x) { return x + x; } };
  template<> struct times__<1> {
    template<typename X> static X is(X const&x) { return x; } };
  template<> struct times__<0> {
    template<typename X> static X is(X const&x) { return X(0); } };
  template<> struct times__<-1> {
    template<typename X> static X is(X const&x) { return -x; } };
  template<> struct times__<-2> {
    template<typename X> static X is(X const&x) { return -x-x; } };
  //----------------------------------------------------------------------------
  template<int N, typename X> inline
  X times (const X&x) { return times__<N>::is(x); }
  //----------------------------------------------------------------------------
  template<typename X> inline
  static X square(X const&x) { return x*x; }
  //----------------------------------------------------------------------------
  template<typename X> inline
  void update_max(X&x, const X&y) { if(y>x) x=y; }
  //----------------------------------------------------------------------------
  template<typename X> inline
  void update_min(X&x, const X&y) { if(y<x) x=y; }
  //----------------------------------------------------------------------------
  template<typename X> inline
  void update_max(X&x, const X&y, const X&a) { if((y+a)>x) x=y+a; }
  //----------------------------------------------------------------------------
  template<typename X> inline
  void update_min(X&x, const X&y, const X&a) { if((y-a)<x) x=y-a; }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class meta::taux<scalar,N,I>                                             //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename X, int N, int I=0> class taux {
    typedef const X       cX;
    typedef taux<X,N,I+1> M;
  public:
    tS  sv s_as  ( X*a,cS&x)      { a[I] =x; M::s_as  (a,x); }
        sv s_ze  ( X*a)           { a[I] =X(0); M::s_ze(a); }
    tS  sv s_ml  ( X*a,cS&x)      { a[I]*=x; M::s_ml(a,x); }
    tS  sv v_ml  ( X*a,cS*b)      { a[I]*=b[I]; M::v_ml(a,b); }
    tS  sv v_dv  ( X*a,cS*b)      { a[I]/=b[I]; M::v_dv(a,b); }
    tS  sv v_as  ( X*a,cS*b)      { a[I] =b[I]; M::v_as(a,b); }
    tS  sv v_ad  ( X*a,cS*b)      { a[I]+=b[I]; M::v_ad(a,b); }
    tS  sv v_su  ( X*a,cS*b)      { a[I]-=b[I]; M::v_su(a,b); }
    tST sv v_ast ( X*a,cS*b,cT&x) { a[I] =x*b[I]; M::v_ast(a,b,x); }
    tST sv v_adt ( X*a,cS*b,cT&x) { a[I]+=x*b[I]; M::v_adt(a,b,x); }
    tST sv v_sut ( X*a,cS*b,cT&x) { a[I]-=x*b[I]; M::v_sut(a,b,x); }
    tZS sv v_asi ( X*a,cX*b)      { a[I] =times<Z>(b[I]); M::v_asi<Z>(a,b); }
    tZS sv v_adi ( X*a,cX*b)      { a[I]+=times<Z>(b[I]); M::v_adi<Z>(a,b); }
    tZS sv v_sui ( X*a,cX*b)      { a[I]-=times<Z>(b[I]); M::v_sui<Z>(a,b); }
    tZS sv v_asti( X*a,cX*b,cS&x) { a[I] =times<Z>(x*b[I]);
                                    M::v_asti<Z>(a,b,x); }
    tZS sv v_adti( X*a,cX*b,cS&x) { a[I]+=times<Z>(x*b[I]);
                                    M::v_adti<Z>(a,b,x); }
    tZS sv v_suti( X*a,cX*b,cS&x) { a[I]-=times<Z>(x*b[I]);
                                    M::v_suti<Z>(a,b,x); }
        sv v_neg ( X*a)           { a[I]=-a[I]; M::v_neg(a); }
    tS  sv v_nega( X*a,cS*b)      { a[I]=-b[I]; M::v_nega(a,b); }
    tS  sv v_sum ( X*a,cX*b,cS*c) { a[I] =b[I]+c[I]; M::v_sum(a,b,c); }
    tS  sv v_dif ( X*a,cX*b,cS*c) { a[I] =b[I]-c[I]; M::v_dif(a,b,c); }
        sb s_eq  (cX*a,cX&b)      { return a[I]==b && M::s_eq(a,b); }
        sb s_neq (cX*a,cX&b)      { return a[I]!=b || M::s_neq(a,b); }
        sb v_eq  (cX*a,cX*b)      { return a[I]==b[I] && M::v_eq(a,b); }
        sb v_neq (cX*a,cX*b)      { return a[I]!=b[I] || M::v_neq(a,b); }
        sX v_norm(cX*a)           { return a[I] *a[I] + M::v_norm(a); }
    tS  sX v_dot (cX*a,cS*b)      { return a[I] *b[I] + M::v_dot(a,b); }
    tS  sX v_diq (cX*a,cS*b)      { return square(a[I]-b[I])+M::v_diq(a,b); }
    tS  sX v_suq (cX*a,cS*b)      { return square(a[I]+b[I])+M::v_suq(a,b); }
        sX v_vol (cX*a)           { return a[I] * M::v_vol(a); }
        sX v_min (cX*a)           { return min(a[I], M::v_min(a)); }
        sX v_max (cX*a)           { return max(a[I], M::v_max(a)); }
        sX v_amax(cX*a)           { return max(abs(a[I]), M::v_amax(a)); }
        sX v_amin(cX*a)           { return min(abs(a[I]), M::v_amin(a)); }
        sb v_nan (cX*a)           { return isnan(a[I]) || M::v_nan(a); }
        sb v_inf (cX*a)           { return isinf(a[I]) || M::v_inf(a); }
        sv v_uma ( X*a,cX*b)      { update_max(a[I],b[I]); M::v_uma(a,b); } 
        sv v_umi ( X*a,cX*b)      { update_min(a[I],b[I]); M::v_umi(a,b); } 
        sv v_umax( X*a,cX*b,cX&x) { update_max(a[I],b[I],x); M::v_umax(a,b,x);}
        sv v_umix( X*a,cX*b,cX&x) { update_min(a[I],b[I],x); M::v_umix(a,b,x);}
    tS  sv v_umia( S*a, S*b,cX*x) { update_min_max(a[I],b[I],x[I]);
                                    M::v_umia(a,b,x); }
        sv v_appl( X*a,cX*b, X(*f)(X)) { a[I]=f(b[I]); M::v_appl(a,b,f); } 
        sv v_outw(std::ostream&o, cX*a, unsigned w)
                                  { o.width(w); o<<a[I]<<' '; M::v_outw(o,a,w);}
        sv v_out (std::ostream&o, cX*a)
                                  { v_outw(o,a,o.width()); }
        sv v_in  (std::istream&i,  X*a) { i>>a[I]; M::v_in(i,a); }
    tP  sv s_app(P &p, cX&x)      { p[I] = x; M::s_app(p,x); }
    tP  sv s_apa(P &p, cX&x)      { p[I]+= x; M::s_apa(p,x); }
    tP  sv s_aps(P &p, cX&x)      { p[I]-= x; M::s_aps(p,x); }
    tP  sv v_ass(X*a,P c_&p)      { a[I] =p[I]; M::v_ass(a,p); }
    tP  sv v_asa(X*a,P c_&p)      { a[I]+=p[I]; M::v_asa(a,p); }
    tP  sv v_asu(X*a,P c_&p)      { a[I]-=p[I]; M::v_asu(a,p); }
    tP  sv v_app(P &p, cX*a)      { p[I] =a[I]; M::v_app(p,a); }
    tP  sv v_apa(P &p, cX*a)      { p[I]+=a[I]; M::v_apa(p,a); }
    tP  sv v_aps(P &p, cX*a)      { p[I]-=a[I]; M::v_aps(p,a); }
    tP  sv v_asst(X*a,P c_&p,cX&x){ a[I] = x*p[I]; M::v_asst(a,p,x); }
    tP  sv v_asat(X*a,P c_&p,cX&x){ a[I]+= x*p[I]; M::v_asat(a,p,x); }
    tP  sv v_asut(X*a,P c_&p,cX&x){ a[I]-= x*p[I]; M::v_asut(a,p,x); }
    tP  sv v_appt(P &p, cX*a,cX&x){ p[I] = x*a[I]; M::v_appt(p,a,x); }
    tP  sv v_apat(P &p, cX*a,cX&x){ p[I]+= x*a[I]; M::v_apat(p,a,x); }
    tP  sv v_apst(P &p, cX*a,cX&x){ p[I]-= x*a[I]; M::v_apst(p,a,x); }
  };
  //----------------------------------------------------------------------------
  template<typename X, int I> class taux<X,I,I> {
    typedef const X      cX;
  public:
    tS  sv s_as  ( X*a,cS&x)      { a[I] =x; }
        sv s_ze  ( X*a)           { a[I] =X(0); }
    tS  sv s_ml  ( X*a,cS&x)      { a[I]*=x; }
    tS  sv v_ml  ( X*a,cS*b)      { a[I]*=b[I]; }
    tS  sv v_dv  ( X*a,cS*b)      { a[I]/=b[I]; }
    tS  sv v_as  ( X*a,cS*b)      { a[I] =b[I]; }
    tS  sv v_ad  ( X*a,cS*b)      { a[I]+=b[I]; }
    tS  sv v_su  ( X*a,cS*b)      { a[I]-=b[I]; }
    tST sv v_ast ( X*a,cS*b,cT&x) { a[I] =x*b[I]; }
    tST sv v_adt ( X*a,cS*b,cT&x) { a[I]+=x*b[I]; }
    tST sv v_sut ( X*a,cS*b,cT&x) { a[I]-=x*b[I]; }
    tZS sv v_asi ( X*a,cX*b)      { a[I] =times<Z>(b[I]); }
    tZS sv v_adi ( X*a,cX*b)      { a[I]+=times<Z>(b[I]); }
    tZS sv v_sui ( X*a,cX*b)      { a[I]-=times<Z>(b[I]); }
    tZS sv v_asti( X*a,cX*b,cS&x) { a[I] =times<Z>(x*b[I]); }
    tZS sv v_adti( X*a,cX*b,cS&x) { a[I]+=times<Z>(x*b[I]); }
    tZS sv v_suti( X*a,cX*b,cS&x) { a[I]-=times<Z>(x*b[I]); }
        sv v_neg ( X*a)           { a[I]=-a[I]; }
    tS  sv v_nega( X*a,cS*b)      { a[I]=-b[I]; }
    tS  sv v_sum ( X*a,cX*b,cS*c) { a[I] =b[I]+c[I]; }
    tS  sv v_dif ( X*a,cX*b,cS*c) { a[I] =b[I]-c[I]; }
        sb s_eq  (cX*a,cX&b)      { return a[I]==b; }
        sb s_neq (cX*a,cX&b)      { return a[I]!=b; }
        sb v_eq  (cX*a,cX*b)      { return a[I]==b[I]; }
        sb v_neq (cX*a,cX*b)      { return a[I]!=b[I]; }
        sX v_norm(cX*a)           { return a[I] *a[I]; }
    tS  sX v_dot (cX*a,cS*b)      { return a[I] *b[I]; }
    tS  sX v_diq (cX*a,cS*b)      { return square(a[I]-b[I]); }
    tS  sX v_suq (cX*a,cS*b)      { return square(a[I]+b[I]); }
        sX v_vol (cX*a)           { return a[I]; }
        sX v_min (cX*a)           { return a[I]; }
        sX v_max (cX*a)           { return a[I]; }
        sX v_amax(cX*a)           { return abs(a[I]); }
        sX v_amin(cX*a)           { return abs(a[I]); }
        sb v_nan (cX*a)           { return isnan(a[I]); }
        sb v_inf (cX*a)           { return isinf(a[I]); }
        sv v_uma ( X*a,cX*b)      { update_max(a[I],b[I]); } 
        sv v_umi ( X*a,cX*b)      { update_min(a[I],b[I]); } 
        sv v_umax( X*a,cX*b,cX&x) { update_max(a[I],b[I],x); }
        sv v_umix( X*a,cX*b,cX&x) { update_min(a[I],b[I],x); }
    tS  sv v_umia( S*a, S*b,cX*x) { update_min_max(a[I],b[I],x[I]); }
        sv v_appl( X*a,cX*b, X(*f)(X)) { a[I]=f(b[I]); } 
        sv v_out (std::ostream&o, cX*a) { o<<a[I]; }
        sv v_outw(std::ostream&o, cX*a, unsigned w)
                                  { o.width(w); o<<a[I]<<' '; }
        sv v_in  (std::istream&i,  X*a) { i>>a[I]; }
    tP  sv s_app (P &p, cX&x)     { p[I] = x; }
    tP  sv s_apa (P &p, cX&x)     { p[I]+= x; }
    tP  sv s_aps (P &p, cX&x)     { p[I]-= x; }
    tP  sv v_ass (X*a,P c_&p)     { a[I] = p[I]; }
    tP  sv v_asa (X*a,P c_&p)     { a[I]+= p[I]; }
    tP  sv v_asu (X*a,P c_&p)     { a[I]-= p[I]; }
    tP  sv v_app (P &p, cX*a)     { p[I] = a[I]; }
    tP  sv v_apa (P &p, cX*a)     { p[I]+= a[I]; }
    tP  sv v_aps (P &p, cX*a)     { p[I]-= a[I]; }
    tP  sv v_asst(X*a,P c_&p,cX&x){ a[I] = x*p[I]; }
    tP  sv v_asat(X*a,P c_&p,cX&x){ a[I]+= x*p[I]; }
    tP  sv v_asut(X*a,P c_&p,cX&x){ a[I]-= x*p[I]; }
    tP  sv v_appt(P &p, cX*a,cX&x){ p[I] = x*a[I]; }
    tP  sv v_apat(P &p, cX*a,cX&x){ p[I]+= x*a[I]; }
    tP  sv v_apst(P &p, cX*a,cX&x){ p[I]-= x*a[I]; }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // struct ONE<N>                                                            //
  //                                                                          //
  // F: factorial N!!                                                         //
  // G: N!!                                                                   //
  // H: (2*N-1)!!                                                             //
  // P: 2^N                                                                   //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<int N> struct ONE {
    enum {
      F  = N*ONE<N-1>::F,
      G  = N*ONE<N-2>::G,
      H  = (N+N-1)*ONE<N-1>::H,
      P  = 4*ONE<N-2>::P
    }; };
  template<> struct ONE<1> { enum { F=1, G=1, H=1, P=2 }; };
  template<> struct ONE<0> { enum { F=1, G=1, H=1, P=1 }; };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // struct TWO<L,M>                                                          //
  //                                                                          //
  // B: binomial (L,M)                                                        //
  // I: index of (L,M)                                                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<int L, int M> struct TWO {
    enum { B = TWO<L-1,M-1>::B + TWO<L-1,M>::B,
	   I = (L*(L+1))/2+M
    }; };
  template<int M> struct TWO<0,M> { enum { B=1, I=0 }; };
  template<int L> struct TWO<L,0> { enum { B=1, I=(L*(L+1))/2 }; };
  template<int L> struct TWO<L,L> { enum { B=1, I=(L*(L+3))/2 }; };
  template<>      struct TWO<0,0> { enum { B=1, I=0 }; };
  //////////////////////////////////////////////////////////////////////////////
} // namespace meta {
} // namespace WDutils {
////////////////////////////////////////////////////////////////////////////////
#undef c_
#undef sv
#undef sb
#undef sX
#undef cS
#undef cT
#undef tX
#undef tS
#undef tST
#undef tP
#undef tZS
////////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_tupel_cc
