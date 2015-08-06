// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file    utils/inc/inline.h
///
/// \author  Walter Dehnen
///
/// \date    1994-2011
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 1994-2011  Walter Dehnen
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
////////////////////////////////////////////////////////////////////////////////
#ifndef WDutils_included_inline_h
#define WDutils_included_inline_h

#ifndef WDutils_included_cstdlib
#  include <cstdlib>
#  define WDutils_included_cstdlib
#endif
#ifndef WDutils_included_cmath
#  include <cmath>
#  define WDutils_included_cmath
#endif
#ifndef WDutils_included_limits
#  include <limits>
#  define WDutils_included_limits
#endif
#ifndef WDutils_included_algorithm
#  include <algorithm>
#  define WDutils_included_algorithm
#endif
//
namespace WDutils {
  /// provide min(a,b) for built-in types
  using std::min;
  /// provide max(a,b) for built-in types
  using std::max;
#ifndef __PGI
  /// provide abs(a) for built-in types
  using std::abs;
#else
  /// the pgCC compiler is faulty: std::abs(float) returns double
  template<typename T> T abs(T x) { return x<T(0)? -x:x; }
#endif
#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
  /// provide isnan(a) for built-in types
  using std::isnan;
  /// provide isinf(a) for built-in types
  using std::isinf;
#else
  /// provide isnan(a) for built-in types
  using ::isnan;   // with some compilers (pgCC, CC) these are not in std
  /// provide isinf(a) for built-in types
  using ::isinf;
#endif
  /// sign(a) for scalar tyupes
  template<typename scalar_type> inline
  int sign(scalar_type x)
  { return (x<0)? -1:((x>0)? 1:0 ); }
  /// x*sign(s)
  template<typename scalar_type> inline
  scalar_type sign (scalar_type x, scalar_type s)
  { return ( s>=0 )?  abs(x) : -abs(x); }
  /// min(a,b,c)
  template<typename scalar_type> inline
  scalar_type min(scalar_type x, scalar_type y, scalar_type z)
  { return min(x,min(y,z)); }
  /// max(a,b,c)
  template<typename scalar_type> inline
  scalar_type max(scalar_type x, scalar_type y, scalar_type z)
  { return max(x,max(y,z)); }
  /// division remainder: mod(x,y) = x-y*int(x/y)
  template<typename scalar_type> inline
  scalar_type mod(scalar_type x, scalar_type y)
  { return x-y*int(x/y); }
  /// 2*x
  template<typename scalar_type> inline
  scalar_type twice(scalar_type x)
  { return x+x; }
  /// 3*x
  template<typename scalar_type> inline
  scalar_type trice(scalar_type x)
  { return 3*x; }
  /// 4*x
  template<typename scalar_type> inline
  scalar_type times4(scalar_type x)
  { return 4*x; }
  //
  template<int N> struct times_ {
    template<typename S> static S is(S const&x) { return N * x; } };
  template<> struct times_<2> {
    template<typename S> static S is(S const&x) { return x + x; } };
  template<> struct times_<1> {
    template<typename S> static S is(S const&x) { return x; } };
  template<> struct times_<0> {
    template<typename S> static S is(S const& ) { return S(0); } };
  template<> struct times_<-1> {
    template<typename S> static S is(S const&x) { return -x; } };
  template<> struct times_<-2> {
    template<typename S> static S is(S const&x) { return -x-x; } };
  /// N*x: uses faster code for |N| < 3
  template<int N, typename scalar_type> inline
  scalar_type times (scalar_type x) { return times_<N>::is(x); }
  //
  template<unsigned, bool ISODD> struct powU;
  template<unsigned N> struct powU<N,0> {
    static const unsigned K = N>>1;
    template<typename S> static S is(S x) { return powU<K,K%2>::is(x*x); } };
  template<unsigned N> struct powU<N,1> {
    static const unsigned K = N>>1;
    template<typename S> static S is(S x) { return x*powU<K,K%2>::is(x*x); } };
  template<> struct powU<2,0> {
    template<typename S> static S is(S x) { return x*x; } };
  template<> struct powU<1,1> {
    template<typename S> static S is(S x) { return x; } };
  template<> struct powU<0,0> {
    template<typename S> static S is(S  ) { return S(1); } };

  template<unsigned N> struct powerU {
    template<typename S> static S is(S x) { return powU<N,N%2>::is(x); } };

  template<int, bool ISNEG> struct powI;
  template<int N> struct powI<N,0> {
    template<typename S> static S is(S x) { return powerU<N>::is(x); } };
  template<int N> struct powI<N,1> {
    template<typename S> static S is(S x) { return powerU<-N>::is(S(1)/x); } };

  template<int N> struct powerI {
    template<typename S> static S is(S x) { return powI<N,N<0>::is(x); } };
  /// integer power of scalar (integer known at compiler time)
  template<int N, typename scalar_type>
  scalar_type power (scalar_type x) { return powerI<N>::is(x); }
  /// x*x
  template<typename scalar_type> inline
  scalar_type square(scalar_type x)
  { return x*x; }
  /// x^3
  template<typename scalar_type> inline
  scalar_type cube (scalar_type x)
  { return x*x*x; }
  /// sqrt(max(x,0))
  template<typename scalar_type> inline
  scalar_type sqrt0(scalar_type x)
  { return x <= scalar_type(0)? scalar_type(0) : std::sqrt(x); }
  /// provide pow(x,z)
  using std::pow;
#ifdef WDutils_non_standard_math
  // integer power of floating point number                                     
  template<typename _Tp> inline
  _Tp pow(const _Tp &x, unsigned int n) {
    if(n==0) return _Tp(1);
    register _Tp z=x, y=(n%2)? x : _Tp(1);
    for(register unsigned int i=n>>1; i; i>>=1) { z*=z; if(i%2) y*=z; }
    return y;
  }
  template<typename _Tp> inline
  _Tp pow(const _Tp &x, int n) {
    if(n==0) return _Tp(1);
    if(n <0) return pow(_Tp(1)/x, unsigned (-n));
    return pow(x, unsigned (n));
  }
#else
  /// single-precision power with unsigned integer exponent
  inline float  pow(float  x, unsigned i) { return float(std::pow(x, int(i))); }
  /// double-precision power with unsigned integer exponent
  inline double pow(double x, unsigned i) { return std::pow(x, int(i)); }
#endif
  /// provide sqrt(x) for built-in floating-point types
  using std::sqrt;
  /// provide exp(x) for built-in floating-point types
  using std::exp;
  /// provide log(x) for built-in floating-point types
  using std::log;
#if defined(__linux) || defined(__DARWIN_UNIX03)
  /// provide cube-root(x) for built-in floating-point types
  using ::cbrt;
#else
  /// single-precision provide cube-root(x)
  inline float cbrt(float x)
  { return float( std::pow( double(x), 0.333333333333333333333 ) ); }
  /// double-precision provide cube-root(x)
  inline double cbrt(double x)
  { return std::pow( x, 0.333333333333333333333 ); }
#endif
  /// update maximum: x = max(x,y)
  template<typename scalar_type> inline
  void update_max(scalar_type&x, scalar_type y)
  { if(y>x) x=y; }
  /// update minimum: x = min(x,y)
  template<typename scalar_type> inline
  void update_min(scalar_type&x, scalar_type y)
  { if(y<x) x=y; }
//   /// update maximum with add: x = max(x,y+a)
//   template<typename scalar_type> inline
//   void update_max(scalar_type&x, scalar_type y, scalar_type a)
//   { update_max(x,y+a); }
//   /// update minimum with sub: x = min(x,y-a)
//   template<typename scalar_type> inline
//   void update_min(scalar_type&x, scalar_type y, scalar_type a)
//   { update_min(x,y-a); }
  /// update minimum and maximum:  min=min(min,x), max=max(max,x)
  template<typename S, typename T> inline
  void update_min_max(S&min, S&max, T x)
  {
    if     (x < min) min = x;
    else if(x > max) max = x;
  }
  /// is floating-point number integral?
  inline bool is_integral(float x)
  { return floor(abs(x))==abs(x);  }
  /// is floating-point number integral?
  inline bool is_integral(double x)
  { return floor(abs(x))==abs(x);  }
} // namespace WDutils {
//
#endif // WDutils_included_inline_h
