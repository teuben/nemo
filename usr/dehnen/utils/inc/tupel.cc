// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    utils/inc/tupel.cc                                                 
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    2003-2011                                                          
///                                                                             
/// \brief   definition of auxiliary methods for template class tupel<>         
///                                                                             
/// \version aug-2003: created                                                  
/// \version sep-2006: made human readable; unused code commented out           
/// \version aug-2011: taux<X,N,I> replaced by taux<X,I>
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 1996-2011  Walter Dehnen                                       
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
#ifndef WDutils_included_limits
#  include <limits>
#  define WDutils_included_limits
#endif
#ifndef WDutils_included_cmath
#  include <cmath>
#  define WDutils_included_cmath
#endif
////////////////////////////////////////////////////////////////////////////////
namespace WDutils {
namespace meta {
  //////////////////////////////////////////////////////////////////////////////
#ifndef WDutils_included_inline_h
  using std::min;
  using std::max;
#ifndef __PGI
  using std::abs;
#else
  // the pgCC compiler is faulty: std::abs(float) returns double
  template<typename T> T abs(T x) { return x<T(0)? -x:x; }
#endif
#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
  using std::isnan;
  using std::isinf;
#else
  using ::isnan;   // with some compilers (pgCC, CC) these are not in std
  using ::isinf;   // with icc, the std:: version causes compiler error
#endif
#endif
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
  template<typename X, int I> class taux
  {
    typedef const X     cX;
    typedef taux<X,I-1> M;
  public:
    /// used in tupel::tupel, in tupel::operator=(scalar), and in falcON
    template<typename S>
    static void s_as(X*a, S const&x)
    { 
      M::s_as(a,x);
      a[I] = x;
    }
    /// used in tupel::reset(), and in falcON
    static void s_ze(X*a)
    {
      M::s_ze(a);
      a[I] = X(0);
    }
    /// used in tupel::operator*=(), and in falcON
    template<typename S> static void s_ml(X*a, S const &x)
    {
      M::s_ml(a,x);
      a[I] *= x;
    }
    /// used in tupel::ass_mul(tupel)
    template<typename S>
    static void v_ml(X*a, const S*b)
    {
      M::v_ml(a,b);
      a[I] *= b[I];
    }
    /// used in tupel::ass_div(tupel)
    template<typename S>
    static void v_dv(X*a, const S*b) {
      M::v_dv(a,b);
      a[I] /= b[I];
    }
    /// used in tupel::tupel, in tupel::operator=(tupel), and in falcON
    template<typename S>
    static void v_as(X*a, const S*b)
    {
      M::v_as(a,b);
      a[I] = b[I];
    }
    /// used in tupel::operator+=(tupel), and in falcON
    template<typename S>
    static void v_ad(X*a, const S*b)
    {
      M::v_ad(a,b);
      a[I] += b[I];
    }
    /// used in tupel::operator-=(tupel), and in falcON
    template<typename S>
    static void v_su(X*a, const S*b)
    {
      M::v_su(a,b);
      a[I] -= b[I];
    }
    /// used in tupel::operator*(scalar), and in falcON
    template<typename S, typename T>
    static void v_ast(X*a, const S*b, T const&x)
    {
      M::v_ast(a,b,x);
      a[I] = x*b[I];
    }
    /// used in falcON
    template<typename S, typename T>
    static void v_adt(X*a, const S*b, T const&x)
    {
      M::v_adt(a,b,x);
      a[I] += x*b[I];
    }
    /// used in falcON
    template<typename S, typename T>
    static void v_sut(X*a, const S*b, T const&x)
    {
      M::v_sut(a,b,x);
      a[I] -= x*b[I];
    }
#if(0)
    /// NOT USED
    template<int Z, typename S>
    static void v_asi(X*a, cX*b)
    {
      M::v_asi<Z>(a,b);
      a[I] = times<Z>(b[I]);
    }
    /// NOT USED
    template<int Z, typename S>
    static void v_adi(X*a, cX*b)
    {
      M::v_adi<Z>(a,b);
      a[I] += times<Z>(b[I]);
    }
    /// NOT USED
    template<int Z, typename S>
    static void v_sui(X*a, cX*b)
    {
      M::v_sui<Z>(a,b);
      a[I] -= times<Z>(b[I]);
    }
    /// NOT USED
    template<int Z, typename S>
    static void v_asti(X*a, cX*b, S const&x)
    {
      M::v_asti<Z>(a,b,x);
      a[I] = times<Z>(x*b[I]);
    }
    /// NOT USED
    template<int Z, typename S>
    static void v_adti(X*a, cX*b, S const&x)
    {
      M::v_adti<Z>(a,b,x);
      a[I] += times<Z>(x*b[I]);
    }
    /// NOT USED
    template<int Z, typename S>
    static void v_suti(X*a, cX*b, S const&x)
    {
      M::v_suti<Z>(a,b,x);
      a[I] -= times<Z>(x*b[I]);
    }
#endif
    /// used in tupel::negate(), and in falcON
    static void v_neg(X*a)
    {
      M::v_neg(a);
      a[I]=-a[I];
    }
    /// used in tupel::operator-() const
    template<typename S>
    static void v_nega(X*a, const S*b)
    {
      M::v_nega(a,b);
      a[I] = -b[I];
    }
    /// used in tupel::operator+(tupel) const
    template<typename S>
    static void v_sum(X*a, cX*b, const S*c)
    {
      M::v_sum(a,b,c);
      a[I] = b[I]+c[I];
    }
    /// used in tupel::operator-(tupel) const
    template<typename S>
    static void v_dif(X*a, cX*b, const S*c)
    {
      M::v_dif(a,b,c);
      a[I] = b[I]-c[I];
    }
    /// used in tupel::operator==(scalar) const, and in falcON
    static bool s_eq(cX*a, cX&b)
    { return M::s_eq(a,b) && a[I]==b; }
    /// used in tupel::operator!=(scalar) const, and in falcON
    static bool s_neq(cX*a, cX&b)
    { return M::s_neq(a,b) || a[I]!=b; }
    /// used in tupel::operator==(tupel) const
    static bool v_eq(cX*a, cX*b)
    { return a[I]==b[I] && M::v_eq(a,b); }
    /// used in tupel::operator!=(tupel) const
    static bool v_neq (cX*a,cX*b)
    { return M::v_neq(a,b) || a[I]!=b[I]; }
    /// used in tupel::norm() const
    static X v_norm(cX*a)
    { return M::v_norm(a) + a[I]*a[I]; }
    /// used in tupel::operator*(tupel) const
    template<typename S>
    static X v_dot(cX*a, const S*b)
    { return M::v_dot(a,b) + a[I]*b[I]; }
    /// used in tupel::dist_sq(tupel) const
    template<typename S>
    static X v_diq(cX*a, const S*b)
    { return M::v_diq(a,b) + square(a[I]-b[I]); }
    /// used in tupel::sum_sq(tupel) const
    template<typename S>
    static X v_suq(cX*a, const S*b)
    { return M::v_suq(a,b) + square(a[I]+b[I]); }
    /// used in tupel::volume() const
    static X v_vol(cX*a)
    { return M::v_vol(a) * a[I]; }
    /// used in tupel::min() const
    static X v_min(cX*a)
    { return min(M::v_min(a),a[I]); }
    /// used in tupel::max() const
    static X v_max(cX*a)
    { return max(M::v_max(a),a[I]); }
    /// used in tupel::maxnorm() const
    static X v_amax(cX*a)
    { return max(M::v_amax(a),std::abs(a[I])); }
    /// used in tupel::minnorm() const
    static X v_amin(cX*a)
    { return min(M::v_amin(a),std::abs(a[I])); }
    /// used in tupel::isnan() const
    static bool v_nan(cX*a)
    { return M::v_nan(a) || isnan(a[I]); }
    /// used in tupel::isinf() const
    static bool v_inf(cX*a)
    { return M::v_inf(a) || isinf(a[I]); }
    /// used in tupel::up_max(tupel)
    static void v_uma(X*a, cX*b)
    {
      M::v_uma(a,b);
      update_max(a[I],b[I]);
    } 
    /// used in tupel::up_min(tupel)
    static void v_umi(X*a, cX*b)
    {
      M::v_umi(a,b);
      update_min(a[I],b[I]);
    }
    /// used in tupel::up_max(tupel, scalar)
    static void v_umax(X*a, cX*b, cX&x)
    {
      M::v_umax(a,b,x);
      update_max(a[I],b[I],x);
    }
    /// used in tupel::up_min(tupel, scalar)
    static void v_umix(X*a, cX*b, cX&x)
    {
      M::v_umix(a,b,x);
      update_min(a[I],b[I],x);
    }
    /// used in tupel::up_min_max(tupel,tupel) const
    template<typename S>
    static void v_umia(S*a, S*b, cX*x)
    {
      M::v_umia(a,b,x);
      update_min_max(a[I],b[I],x[I]);
    }
    /// used in tupel::apply(func)
    static void v_appl(X*a, cX*b, X(*f)(X))
    {
      M::v_appl(a,b,f);
      a[I]=f(b[I]);
    }
    /// used in v_out
    static void v_outwp(std::ostream&o, cX*a, unsigned w, unsigned p)
    {
      M::v_outwp(o,a,w,p);
      o<<' ';
      o.width(w);
      o.precision(p);
      o<<a[I];
    }
    /// used in operator<< (std::ostream, tupel)
    static void v_out(std::ostream&o, cX*a)
    { v_outwp(o,a,o.width(),o.precision()); }
    /// used in operator>> (std::istream, tupel)
    static void v_in(std::istream&i, X*a)
    {
      M::v_in(i,a);
      i>>a[I];
    }
#if(0)
    /// NOT USED
    template<typename P>
    static void s_app(P&p, cX&x)
    {
      M::s_app(p,x);
      p[I] = x;
    }
    /// NOT USED
    template<typename P>
    static void s_apa(P&p, cX&x)
    {
      M::s_apa(p,x);
      p[I] += x;
    }
    /// NOT USED
    template<typename P>
    static void s_aps(P&p, cX&x)
    {
      M::s_aps(p,x);
      p[I] -= x;
    }
    /// NOT USED
    template<typename P>
    static void v_ass(X*a, P const&p)
    {
      M::v_ass(a,p);
      a[I] = p[I];
    }
    /// NOT USED
    template<typename P>
    static void v_asa(X*a, P const&p)
    {
      M::v_asa(a,p);
      a[I] += p[I];
    }
    /// NOT USED
    template<typename P>
    static void v_asu(X*a, P const&p)
    {
      M::v_asu(a,p);
      a[I] -= p[I];
    }
    /// NOT USED
    template<typename P>
    static void v_app(P&p, cX*a)
    {
      M::v_app(p,a);
      p[I] = a[I];
    }
    /// NOT USED
    template<typename P>
    static void v_apa(P&p, cX*a)
    {
      M::v_apa(p,a);
      p[I] += a[I];
    }
    /// NOT USED
    template<typename P>
    static void v_aps(P&p, cX*a)
    {
      M::v_aps(p,a);
      p[I] -= a[I];
    }
    /// NOT USED
    template<typename P>
    static void v_asst(X*a,P const&p,cX&x)
    {
      M::v_asst(a,p,x);
      a[I] = x*p[I];
    }
    /// NOT USED
    template<typename P>
    static void v_asat(X*a, P const&p, cX&x)
    {
      M::v_asat(a,p,x);
      a[I] += x*p[I];
    }
    /// NOT USED
    template<typename P>
    static void v_asut(X*a, P const&p, cX&x)
    {
      M::v_asut(a,p,x);
      a[I] -= x*p[I];
    }
    /// NOT USED
    template<typename P>
    static void v_appt(P&p, cX*a, cX&x)
    {
      M::v_appt(p,a,x);
      p[I] = x*a[I];
    }
    /// NOT USED
    template<typename P>
    static void v_apat(P&p, cX*a, cX&x)
    {
      M::v_apat(p,a,x);
      p[I] += x*a[I];
    }
    /// NOT USED
    template<typename P>
    static void v_apst(P&p, cX*a, cX&x)
    {
      M::v_apst(p,a,x);
      p[I] -= x*a[I];
    }
#endif
  };
  //----------------------------------------------------------------------------
  // case N==0: truncation of loop                                              
  //----------------------------------------------------------------------------
  template<typename X> class taux<X,0>
  {
    typedef const X cX;
  public:
    //
    template<typename S>
    static void s_as(X*a, S const&x)
    { a[0] = x; }
    //
    static void s_ze(X*a)
    { a[0] = X(0); }
    //
    template<typename S>
    static void s_ml(X*a, S const &x)
    { a[0] *= x; }
    //
    template<typename S>
    static void v_ml(X*a, const S*b)
    { a[0] *= b[0]; }
    //
    template<typename S>
    static void v_dv(X*a, const S*b)
    { a[0] /= b[0]; }
    //
    template<typename S>
    static void v_as(X*a, const S*b)
    { a[0] = b[0]; }
    //
    template<typename S>
    static void v_ad(X*a, const S*b)
    { a[0] += b[0]; }
    //
    template<typename S>
    static void v_su(X*a, const S*b)
    { a[0] -= b[0]; }
    //
    template<typename S, typename T>
    static void v_ast(X*a, const S*b, T const&x)
    { a[0] = x*b[0]; }
    //
    template<typename S, typename T>
    static void v_adt(X*a, const S*b, T const&x)
    { a[0] += x*b[0]; }
    //
    template<typename S, typename T>
    static void v_sut(X*a, const S*b, T const&x)
    { a[0] -= x*b[0]; }
#if(0)
    //
    template<int Z, typename S>
    static void v_asi(X*a, cX*b)
    { a[0] = times<Z>(b[0]); }
    //
    template<int Z, typename S>
    static void v_adi(X*a, cX*b)
    { a[0] += times<Z>(b[0]); }
    //
    template<int Z, typename S>
    static void v_sui(X*a, cX*b)
    { a[0] -= times<Z>(b[0]); }
    //
    template<int Z, typename S>
    static void v_asti(X*a, cX*b, S const&x)
    { a[0] = times<Z>(x*b[0]); }
    //
    template<int Z, typename S>
    static void v_adti(X*a, cX*b, S const&x)
    { a[0] += times<Z>(x*b[0]); }
    //
    template<int Z, typename S>
    static void v_suti(X*a, cX*b, S const&x)
    { a[0] -= times<Z>(x*b[0]); }
#endif
    //
    static void v_neg(X*a)
    { a[0] = -a[0]; }
    //
    template<typename S>
    static void v_nega(X*a, const S*b)
    { a[0] = -b[0]; }
    //
    template<typename S>
    static void v_sum(X*a, cX*b, const S*c)
    { a[0] = b[0]+c[0]; }
    //
    template<typename S>
    static void v_dif(X*a, cX*b, const S*c)
    { a[0] = b[0]-c[0]; }
    //
    static bool s_eq(cX*a, cX&b)
    { return a[0]==b; }
    //
    static bool s_neq(cX*a, cX&b)
    { return a[0]!=b; }
    //
    static bool v_eq(cX*a, cX*b)
    { return a[0]==b[0]; }
    //
    static bool v_neq (cX*a,cX*b)
    { return a[0]!=b[0]; }
    //
    static X v_norm(cX*a)
    { return a[0]*a[0]; }
    //
    template<typename S>
    static X v_dot(cX*a, const S*b)
    { return a[0]*b[0]; }
    //
    template<typename S>
    static X v_diq(cX*a, const S*b)
    { return square(a[0]-b[0]); }
    //
    template<typename S>
    static X v_suq(cX*a, const S*b)
    { return square(a[0]+b[0]); }
    //
    static X v_vol(cX*a)
    { return a[0]; }
    //
    static X v_min(cX*a)
    { return a[0]; }
    //
    static X v_max(cX*a)
    { return a[0]; }
    //
    static X v_amax(cX*a)
    { return std::abs(a[0]); }
    //
    static X v_amin(cX*a)
    { return std::abs(a[0]); }
    //
    static bool v_nan(cX*a)
    { return isnan(a[0]); }
    //
    static bool v_inf(cX*a)
    { return isinf(a[0]); }
    //
    static void v_uma(X*a, cX*b) 
    { update_max(a[0],b[0]); } 
    //
    static void v_umi(X*a, cX*b)
    { update_min(a[0],b[0]); } 
    //
    static void v_umax(X*a, cX*b, cX&x)
    { update_max(a[0],b[0],x); }
    //
    static void v_umix(X*a, cX*b, cX&x)
    { update_min(a[0],b[0],x); }
    //
    template<typename S>
    static void v_umia(S*a, S*b, cX*x)
    { update_min_max(a[0],b[0],x[0]); }
    //
    static void v_appl(X*a, cX*b, X(*f)(X))
    { a[0]=f(b[0]); } 
    //
    static void v_outwp(std::ostream&o, cX*a, unsigned w, unsigned p)
    {
      o.width(w);
      o.precision(p);
      o<<a[0];
    }
    //
    static void v_out(std::ostream&o, cX*a)
    { o<<a[0]; }
    //
    static void v_in(std::istream&i, X*a)
    { i>>a[0]; }
#if(0)
    //
    template<typename P>
    static void s_app(P&p, cX&x)
    { p[0] = x; }
    //
    template<typename P>
    static void s_apa(P&p, cX&x)
    { p[0] += x; }
    //
    template<typename P>
    static void s_aps(P&p, cX&x)
    { p[0] -= x; }
    //
    template<typename P>
    static void v_ass(X*a, P const&p)
    { a[0] = p[0]; }
    //
    template<typename P>
    static void v_asa(X*a, P const&p)
    { a[0] += p[0]; }
    //
    template<typename P>
    static void v_asu(X*a, P const&p)
    { a[0] -= p[0]; }
    //
    template<typename P>
    static void v_app(P&p, cX*a)
    { p[0] = a[0]; }
    //
    template<typename P>
    static void v_apa(P&p, cX*a)
    { p[0] += a[0]; }
    //
    template<typename P>
    static void v_aps(P&p, cX*a)
    { p[0] -= a[0]; }
    //
    template<typename P>
    static void v_asst(X*a,P const&p,cX&x)
    { a[0] = x*p[0]; }
    //
    template<typename P>
    static void v_asat(X*a, P const&p, cX&x)
    { a[0] += x*p[0]; }
    //
    template<typename P>
    static void v_asut(X*a, P const&p, cX&x)
    { a[0] -= x*p[0]; }
    //
    template<typename P>
    static void v_appt(P&p, cX*a, cX&x)
    { p[0] = x*a[0]; }
    //
    template<typename P>
    static void v_apat(P&p, cX*a, cX&x)
    { p[0] += x*a[0]; }
    //
    template<typename P>
    static void v_apst(P&p, cX*a, cX&x)
    { p[0] -= x*a[0]; }
#endif
  };
} // namespace meta {
} // namespace WDutils {
////////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_tupel_cc
