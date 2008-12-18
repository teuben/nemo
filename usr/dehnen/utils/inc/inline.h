// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    utils/inc/inline.h                                                 
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    1994-2005, 2007                                                    
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 1994-2005, 2007  Walter Dehnen                                 
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
#ifndef WDutils_included_algorithm
#  include <algorithm>
#  define WDutils_included_algorithm
#endif
//------------------------------------------------------------------------------
namespace WDutils {
  using std::abs;
  using std::min;
  using std::max;
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  int sign    (const scalar_type&x)
  { return (x<0)? -1:((x>0)? 1:0 ); }
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  scalar_type sign (const scalar_type&x, const scalar_type&s)
  { return ( s>=0 )?  abs(x) : -abs(x); }
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  const scalar_type&min(const scalar_type&x,
			const scalar_type&y,
			const scalar_type&z)
  { return min(x,min(y,z)); }
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  const scalar_type&max(const scalar_type&x,
			const scalar_type&y,
			const scalar_type&z)
  { return max(x,max(y,z)); }
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  scalar_type mod(const scalar_type&x,
		  const scalar_type&y)
  { return x-y*int(x/y); }
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  scalar_type twice(const scalar_type&x)
  { return x+x; }
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  scalar_type trice(const scalar_type&x)
  { return 3*x; }
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  scalar_type times4(const scalar_type&x)
  { return 4*x; }
  //----------------------------------------------------------------------------
  template<int N> struct times__ {
    template<typename S> static S is(S const&x) { return N * x; } };
  template<> struct times__<2> {
    template<typename S> static S is(S const&x) { return x + x; } };
  template<> struct times__<1> {
    template<typename S> static S is(S const&x) { return x; } };
  template<> struct times__<0> {
    template<typename S> static S is(S const&x) { return S(0); } };
  template<> struct times__<-1> {
    template<typename S> static S is(S const&x) { return -x; } };
  template<> struct times__<-2> {
    template<typename S> static S is(S const&x) { return -x-x; } };

  template<int N, typename scalar_type> inline
  scalar_type times (const scalar_type&x) { return times__<N>::is(x); }
  //----------------------------------------------------------------------------
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
    template<typename S> static S is(S x) { return S(1); } };

  template<unsigned N> struct powerU {
    template<typename S> static S is(S x) { return powU<N,N%2>::is(x); } };

  template<int, bool ISNEG> struct powI;
  template<int N> struct powI<N,0> {
    template<typename S> static S is(S x) { return powerU<N>::is(x); } };
  template<int N> struct powI<N,1> {
    template<typename S> static S is(S x) { return powerU<-N>::is(S(1)/x); } };

  template<int N> struct powerI {
    template<typename S> static S is(S x) { return powI<N,N<0>::is(x); } };

  template<int N, typename scalar_type>
  scalar_type power (const scalar_type&x) { return powerI<N>::is(x); }
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  scalar_type square(const scalar_type&x)
  { return x*x; }
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  scalar_type cube (const scalar_type&x)
  { return x*x*x; }
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  scalar_type sqrt0(scalar_type const&x)
  { return x <= scalar_type(0)? scalar_type(0) : std::sqrt(x); }
  //----------------------------------------------------------------------------
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
  inline float  pow(float  x, unsigned i) { return std::pow(x, int(i)); }
  inline double pow(double x, unsigned i) { return std::pow(x, int(i)); }
#endif
  using std::sqrt;
  using std::exp;
  using std::log;
#ifdef linux
  using ::cbrt;
#else
  inline float cbrt(float x)          { 
    return float( std::pow( double(x), 0.333333333333333333333 ) );
  }
  inline double cbrt(double x)        { 
    return std::pow( x, 0.333333333333333333333 );
  }
#endif
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  void update_max(scalar_type&x, const scalar_type&y)
  {
    if(y>x) x=y;
  }
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  void update_min(scalar_type&x, const scalar_type&y)
  { 
    if(y<x) x=y;
  }
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  void update_max(scalar_type&x, const scalar_type&y, const scalar_type&a)
  {
    if((y+a)>x) x=y+a;
  }
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  void update_min(scalar_type&x, const scalar_type&y, const scalar_type&a)
  { 
    if((y-a)<x) x=y-a;
  }
  //----------------------------------------------------------------------------
  template<typename S, typename T> inline
  void update_min_max(S&min, S&max, T const&x)
  {
    if     (x < min) min = x;
    else if(x > max) max = x;
  }
  //----------------------------------------------------------------------------
  inline bool is_integral(const float&x)
  { return floor(abs(x))==abs(x);  }
  //----------------------------------------------------------------------------
  inline bool is_integral(const double&x)
  { return floor(abs(x))==abs(x);  }
} // namespace WDutils {
//------------------------------------------------------------------------------
#endif // WDutils_included_inline_h
