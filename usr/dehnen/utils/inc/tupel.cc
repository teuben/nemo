// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    utils/inc/tupel.cc                                                 
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    2003-2006                                                          
///                                                                             
/// \brief   definition of auxiliary methods for template class tupel<>         
///                                                                             
/// \version aug-2003: created                                                  
/// \version sep-2006: made human readable; unused code commented out           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 1996-2006  Walter Dehnen                                       
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
  //                                                                            
  /// template metaprogramming: unroll loops over dimensions for class tupel    
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  template<typename X, int N, int I=0> class taux {
    typedef const X       cX;
    typedef taux<X,N,I+1> M;
  public:
    /// used in tupel::tupel, in tupel::operator=(scalar), and in falcON
    template<typename S>
    static void s_as(X*a, S const&x) {
      a[I] = x;
      M::s_as(a,x);
    }
    /// used in tupel::reset(), and in falcON
    static void s_ze(X*a) {
      a[I] = X(0);
      M::s_ze(a);
    }
    /// used in tupel::operator*=(), and in falcON
    template<typename S>
    static void s_ml(X*a, S const &x) {
      a[I] *= x;
      M::s_ml(a,x);
    }
    /// used in tupel::ass_mul(tupel)
    template<typename S>
    static void v_ml(X*a, const S*b) {
      a[I] *= b[I];
      M::v_ml(a,b);
    }
    /// used in tupel::ass_div(tupel)
    template<typename S>
    static void v_dv(X*a, const S*b) {
      a[I] /= b[I];
      M::v_dv(a,b);
    }
    /// used in tupel::tupel, in tupel::operator=(tupel), and in falcON
    template<typename S>
    static void v_as(X*a, const S*b) {
      a[I] = b[I];
      M::v_as(a,b);
    }
    /// used in tupel::operator+=(tupel), and in falcON
    template<typename S>
    static void v_ad(X*a, const S*b) {
      a[I] += b[I];
      M::v_ad(a,b);
    }
    /// used in tupel::operator-=(tupel), and in falcON
    template<typename S>
    static void v_su(X*a, const S*b) {
      a[I] -= b[I];
      M::v_su(a,b);
    }
    /// used in tupel::operator*(scalar), and in falcON
    template<typename S, typename T>
    static void v_ast(X*a, const S*b, T const&x) {
      a[I] = x*b[I];
      M::v_ast(a,b,x);
    }
    /// used in falcON
    template<typename S, typename T>
    static void v_adt(X*a, const S*b, T const&x) {
      a[I] += x*b[I];
      M::v_adt(a,b,x);
    }
    /// used in falcON
    template<typename S, typename T>
    static void v_sut(X*a, const S*b, T const&x) {
      a[I] -= x*b[I];
      M::v_sut(a,b,x);
    }
#if(0)
    /// NOT USED
    template<int Z, typename S>
    static void v_asi(X*a, cX*b) {
      a[I] = times<Z>(b[I]);
      M::v_asi<Z>(a,b);
    }
    /// NOT USED
    template<int Z, typename S>
    static void v_adi(X*a, cX*b) {
      a[I] += times<Z>(b[I]);
      M::v_adi<Z>(a,b);
    }
    /// NOT USED
    template<int Z, typename S>
    static void v_sui(X*a, cX*b) {
      a[I] -= times<Z>(b[I]);
      M::v_sui<Z>(a,b);
    }
    /// NOT USED
    template<int Z, typename S>
    static void v_asti(X*a, cX*b, S const&x) {
      a[I] = times<Z>(x*b[I]);
      M::v_asti<Z>(a,b,x);
    }
    /// NOT USED
    template<int Z, typename S>
    static void v_adti(X*a, cX*b, S const&x) {
      a[I] += times<Z>(x*b[I]);
      M::v_adti<Z>(a,b,x);
    }
    /// NOT USED
    template<int Z, typename S>
    static void v_suti(X*a, cX*b, S const&x) {
      a[I] -= times<Z>(x*b[I]);
      M::v_suti<Z>(a,b,x);
    }
#endif
    /// used in tupel::negate(), and in falcON
    static void v_neg(X*a) {
      a[I]=-a[I];
      M::v_neg(a);
    }
    /// used in tupel::operator-() const
    template<typename S>
    static void v_nega(X*a, const S*b) {
      a[I] = -b[I];
      M::v_nega(a,b);
    }
    /// used in tupel::operator+(tupel) const
    template<typename S>
    static void v_sum(X*a, cX*b, const S*c) {
      a[I] = b[I]+c[I];
      M::v_sum(a,b,c);
    }
    /// used in tupel::operator-(tupel) const
    template<typename S>
    static void v_dif(X*a, cX*b, const S*c) {
      a[I] = b[I]-c[I];
      M::v_dif(a,b,c);
    }
    /// used in tupel::operator==(scalar) const, and in falcON
    static bool s_eq(cX*a, cX&b) {
      return a[I]==b && M::s_eq(a,b);
    }
    /// used in tupel::operator!=(scalar) const, and in falcON
    static bool s_neq(cX*a, cX&b) {
      return a[I]!=b || M::s_neq(a,b);
    }
    /// used in tupel::operator==(tupel) const
    static bool v_eq(cX*a, cX*b) {
      return a[I]==b[I] && M::v_eq(a,b);
    }
    /// used in tupel::operator!=(tupel) const
    static bool v_neq (cX*a,cX*b) {
      return a[I]!=b[I] || M::v_neq(a,b);
    }
    /// used in tupel::norm() const
    static X v_norm(cX*a) {
      return a[I] *a[I] + M::v_norm(a);
    }
    /// used in tupel::operator*(tupel) const
    template<typename S>
    static X v_dot(cX*a, const S*b) {
      return a[I] *b[I] + M::v_dot(a,b);
    }
    /// used in tupel::dist_sq(tupel) const
    template<typename S>
    static X v_diq(cX*a, const S*b) {
      return square(a[I]-b[I])+M::v_diq(a,b);
    }
    /// used in tupel::sum_sq(tupel) const
    template<typename S>
    static X v_suq(cX*a, const S*b) {
      return square(a[I]+b[I])+M::v_suq(a,b);
    }
    /// used in tupel::volume() const
    static X v_vol(cX*a) {
      return a[I] * M::v_vol(a);
    }
    /// used in tupel::min() const
    static X v_min(cX*a) {
      return min(a[I], M::v_min(a));
    }
    /// used in tupel::max() const
    static X v_max(cX*a) {
      return max(a[I], M::v_max(a));
    }
    /// used in tupel::maxnorm() const
    static X v_amax(cX*a) {
      return max(abs(a[I]), M::v_amax(a));
    }
    /// used in tupel::minnorm() const
    static X v_amin(cX*a) {
      return min(abs(a[I]), M::v_amin(a));
    }
    /// used in tupel::isnan() const
    static bool v_nan(cX*a) {
      return isnan(a[I]) || M::v_nan(a);
    }
    /// used in tupel::isinf() const
    static bool v_inf(cX*a) {
      return isinf(a[I]) || M::v_inf(a);
    }
    /// used in tupel::up_max(tupel)
    static void v_uma(X*a, cX*b) {
      update_max(a[I],b[I]);
      M::v_uma(a,b);
    } 
    /// used in tupel::up_min(tupel)
    static void v_umi(X*a, cX*b) {
      update_min(a[I],b[I]);
      M::v_umi(a,b);
    }
    /// used in tupel::up_max(tupel, scalar)
    static void v_umax(X*a, cX*b, cX&x) {
      update_max(a[I],b[I],x);
      M::v_umax(a,b,x);
    }
    /// used in tupel::up_min(tupel, scalar)
    static void v_umix(X*a, cX*b, cX&x) {
      update_min(a[I],b[I],x);
      M::v_umix(a,b,x);
    }
    /// used in tupel::up_min_max(tupel,tupel) const
    template<typename S>
    static void v_umia(S*a, S*b, cX*x) {
      update_min_max(a[I],b[I],x[I]);
      M::v_umia(a,b,x);
    }
    /// used in tupel::apply(func)
    static void v_appl(X*a, cX*b, X(*f)(X)) {
      a[I]=f(b[I]);
      M::v_appl(a,b,f);
    }
    /// used in v_out
    static void v_outw(std::ostream&o, cX*a, unsigned w) {
      o.width(w);
      o<<a[I]<<' ';
      M::v_outw(o,a,w);
    }
    /// used in operator<< (std::ostream, tupel)
    static void v_out(std::ostream&o, cX*a) {
      v_outw(o,a,o.width());
    }
    /// used in operator>> (std::istream, tupel)
    static void v_in(std::istream&i, X*a) {
      i>>a[I];
      M::v_in(i,a);
    }
#if(0)
    /// NOT USED
    template<typename P>
    static void s_app(P&p, cX&x) {
      p[I] = x;
      M::s_app(p,x);
    }
    /// NOT USED
    template<typename P>
    static void s_apa(P&p, cX&x) {
      p[I] += x;
      M::s_apa(p,x);
    }
    /// NOT USED
    template<typename P>
    static void s_aps(P&p, cX&x) {
      p[I] -= x;
      M::s_aps(p,x);
    }
    /// NOT USED
    template<typename P>
    static void v_ass(X*a, P const&p) {
      a[I] = p[I];
      M::v_ass(a,p);
    }
    /// NOT USED
    template<typename P>
    static void v_asa(X*a, P const&p) {
      a[I] += p[I];
      M::v_asa(a,p);
    }
    /// NOT USED
    template<typename P>
    static void v_asu(X*a, P const&p) {
      a[I] -= p[I];
      M::v_asu(a,p);
    }
    /// NOT USED
    template<typename P>
    static void v_app(P&p, cX*a) {
      p[I] = a[I];
      M::v_app(p,a);
    }
    /// NOT USED
    template<typename P>
    static void v_apa(P&p, cX*a) {
      p[I] += a[I];
      M::v_apa(p,a);
    }
    /// NOT USED
    template<typename P>
    static void v_aps(P&p, cX*a) {
      p[I] -= a[I];
      M::v_aps(p,a);
    }
    /// NOT USED
    template<typename P>
    static void v_asst(X*a,P const&p,cX&x) {
      a[I] = x*p[I];
      M::v_asst(a,p,x);
    }
    /// NOT USED
    template<typename P>
    static void v_asat(X*a, P const&p, cX&x) {
      a[I] += x*p[I];
      M::v_asat(a,p,x);
    }
    /// NOT USED
    template<typename P>
    static void v_asut(X*a, P const&p, cX&x) {
      a[I] -= x*p[I];
      M::v_asut(a,p,x);
    }
    /// NOT USED
    template<typename P>
    static void v_appt(P&p, cX*a, cX&x) {
      p[I] = x*a[I];
      M::v_appt(p,a,x);
    }
    /// NOT USED
    template<typename P>
    static void v_apat(P&p, cX*a, cX&x) {
      p[I] += x*a[I];
      M::v_apat(p,a,x);
    }
    /// NOT USED
    template<typename P>
    static void v_apst(P&p, cX*a, cX&x) {
      p[I] -= x*a[I];
      M::v_apst(p,a,x);
    }
#endif
  };
  //----------------------------------------------------------------------------
  // case N==I: truncation of loop                                              
  //----------------------------------------------------------------------------
  template<typename X, int I> class taux<X,I,I> {
    typedef const X      cX;
  public:
    template<typename S>
    static void s_as(X*a, S const&x) {
      a[I] = x;
    }
    static void s_ze(X*a) {
      a[I] = X(0);
    }
    template<typename S>
    static void s_ml(X*a, S const &x) {
      a[I] *= x;
    }
    template<typename S>
    static void v_ml(X*a, const S*b) {
      a[I] *= b[I];
    }
    template<typename S>
    static void v_dv(X*a, const S*b) {
      a[I] /= b[I];
    }
    template<typename S>
    static void v_as(X*a, const S*b) {
      a[I] = b[I];
    }
    template<typename S>
    static void v_ad(X*a, const S*b) {
      a[I] += b[I];
    }
    template<typename S>
    static void v_su(X*a, const S*b) {
      a[I] -= b[I];
    }
    template<typename S, typename T>
    static void v_ast(X*a, const S*b, T const&x) {
      a[I] = x*b[I];
    }
    template<typename S, typename T>
    static void v_adt(X*a, const S*b, T const&x) {
      a[I] += x*b[I];
    }
    template<typename S, typename T>
    static void v_sut(X*a, const S*b, T const&x) {
      a[I] -= x*b[I];
    }
#if(0)
    template<int Z, typename S>
    static void v_asi(X*a, cX*b) {
      a[I] = times<Z>(b[I]);
    }
    template<int Z, typename S>
    static void v_adi(X*a, cX*b) {
      a[I] += times<Z>(b[I]);
    }
    template<int Z, typename S>
    static void v_sui(X*a, cX*b) {
      a[I] -= times<Z>(b[I]);
    }
    template<int Z, typename S>
    static void v_asti(X*a, cX*b, S const&x) {
      a[I] = times<Z>(x*b[I]);
    }
    template<int Z, typename S>
    static void v_adti(X*a, cX*b, S const&x) {
      a[I] += times<Z>(x*b[I]);
    }
    template<int Z, typename S>
    static void v_suti(X*a, cX*b, S const&x) {
      a[I] -= times<Z>(x*b[I]);
    }
#endif
    static void v_neg(X*a) {
      a[I]=-a[I];
    }
    template<typename S>
    static void v_nega(X*a, const S*b) {
      a[I] = -b[I];
    }
    template<typename S>
    static void v_sum(X*a, cX*b, const S*c) {
      a[I] = b[I]+c[I];
    }
    template<typename S>
    static void v_dif(X*a, cX*b, const S*c) {
      a[I] = b[I]-c[I];
    }
    static bool s_eq(cX*a, cX&b) {
      return a[I]==b;
    }
    static bool s_neq(cX*a, cX&b) {
      return a[I]!=b;
    }
    static bool v_eq(cX*a, cX*b) {
      return a[I]==b[I];
    }
    static bool v_neq (cX*a,cX*b) {
      return a[I]!=b[I];
    }
    static X v_norm(cX*a) {
      return a[I] *a[I];
    }
    template<typename S>
    static X v_dot(cX*a, const S*b) {
      return a[I] *b[I];
    }
    template<typename S>
    static X v_diq(cX*a, const S*b) {
      return square(a[I]-b[I]);
    }
    template<typename S>
    static X v_suq(cX*a, const S*b) {
      return square(a[I]+b[I]);
    }
    static X v_vol(cX*a) {
      return a[I];
    }
    static X v_min(cX*a) {
      return a[I];
    }
    static X v_max(cX*a) {
      return a[I];
    }
    static X v_amax(cX*a) {
      return abs(a[I]);
    }
    static X v_amin(cX*a) {
      return abs(a[I]);
    }
    static bool v_nan(cX*a) {
      return isnan(a[I]);
    }
    static bool v_inf(cX*a) {
      return isinf(a[I]);
    }
    static void v_uma(X*a, cX*b) {
      update_max(a[I],b[I]);
    } 
    static void v_umi(X*a, cX*b) {
      update_min(a[I],b[I]);
    } 
    static void v_umax(X*a, cX*b, cX&x) {
      update_max(a[I],b[I],x);
    }
    static void v_umix(X*a, cX*b, cX&x) {
      update_min(a[I],b[I],x);
    }
    template<typename S>
    static void v_umia(S*a, S*b, cX*x) {
      update_min_max(a[I],b[I],x[I]);
    }
    static void v_appl(X*a, cX*b, X(*f)(X)) {
      a[I]=f(b[I]);
    } 
    static void v_outw(std::ostream&o, cX*a, unsigned w) {
      o.width(w);
      o<<a[I];
    }
    static void v_out(std::ostream&o, cX*a) {
      v_outw(o,a,o.width());
    }
    static void v_in(std::istream&i, X*a) {
      i>>a[I];
    }
#if(0)
    template<typename P>
    static void s_app(P&p, cX&x) {
      p[I] = x;
    }
    template<typename P>
    static void s_apa(P&p, cX&x) {
      p[I] += x;
    }
    template<typename P>
    static void s_aps(P&p, cX&x) {
      p[I] -= x;
    }
    template<typename P>
    static void v_ass(X*a, P const&p) {
      a[I] = p[I];
    }
    template<typename P>
    static void v_asa(X*a, P const&p) {
      a[I] += p[I];
    }
    template<typename P>
    static void v_asu(X*a, P const&p) {
      a[I] -= p[I];
    }
    template<typename P>
    static void v_app(P&p, cX*a) {
      p[I] = a[I];
    }
    template<typename P>
    static void v_apa(P&p, cX*a) {
      p[I] += a[I];
    }
    template<typename P>
    static void v_aps(P&p, cX*a) {
      p[I] -= a[I];
    }
    template<typename P>
    static void v_asst(X*a,P const&p,cX&x) {
      a[I] = x*p[I];
    }
    template<typename P>
    static void v_asat(X*a, P const&p, cX&x) {
      a[I] += x*p[I];
    }
    template<typename P>
    static void v_asut(X*a, P const&p, cX&x) {
      a[I] -= x*p[I];
    }
    template<typename P>
    static void v_appt(P&p, cX*a, cX&x) {
      p[I] = x*a[I];
    }
    template<typename P>
    static void v_apat(P&p, cX*a, cX&x) {
      p[I] += x*a[I];
    }
    template<typename P>
    static void v_apst(P&p, cX*a, cX&x) {
      p[I] -= x*a[I];
    }
#endif
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
#endif // WDutils_included_tupel_cc
