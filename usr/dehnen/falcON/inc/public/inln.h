// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// inln.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 1994-2003                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_inln_h
#define falcON_included_inln_h

#ifndef falcON_included_cstdlib
#  include <cstdlib>
#  define falcON_included_cstdlib
#endif
#ifndef falcON_included_cmath
#  include <cmath>
#  define falcON_included_cmath
#endif
#ifndef falcON_included_algorithm
#  include <algorithm>
#  define falcON_included_algorithm
#endif
#ifndef falcON_included_frst_h
#  include <public/frst.h>
#endif
#ifndef falcON_included_exit_h
#  include <public/exit.h>
#endif
//------------------------------------------------------------------------------
namespace nbdy {
  using std::abs;
#ifndef __GNUC__
  template<typename scalar_type> inline
  scalar_type   abs     (const scalar_type&x)
  { return (x<0)? -x : x; }
  //----------------------------------------------------------------------------
  template<> inline
  double abs<double> (const double&x)
  { return fabs(x); }
  //----------------------------------------------------------------------------
  template<> inline
  float abs<float>   (const float &x) 
  { return fabs(x); }
#endif
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  int sign    (const scalar_type&x)
  { return (x<0)? -1:((x>0)? 1:0 ); }
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  scalar_type sign (const scalar_type&x, const scalar_type&s)
  { return ( s>0 )?  abs(x) : -abs(x); }
  //----------------------------------------------------------------------------
#if(0)
  template<typename scalar_type> inline
  const scalar_type&min(const scalar_type&x,
			const scalar_type&y)
  { 
#ifdef __GNUC__
    return x <? y;
#else
    return (x<y)? x : y;
#endif
  }
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  const scalar_type&max(const scalar_type&x,
			const scalar_type&y)
  {
#ifdef __GNUC__
    return x >? y;
#else
    return (x>y)? x : y;
#endif
  }
  //----------------------------------------------------------------------------
#else // !0: use min & max from <algorithm>
  using std::min;
  using std::max;
#endif
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

  template<unsigned N, typename scalar_type>
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
  void      swap    (scalar_type&a, scalar_type&b)
  { register scalar_type t=a; a=b; b=t; }
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  scalar_type sqrt0(scalar_type const&x)
  { return x <= scalar_type(0)? scalar_type(0) : std::sqrt(x); }
  //----------------------------------------------------------------------------
#ifdef falcON_non_standard_math
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
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  void update_max(scalar_type&x, const scalar_type&y)
  {
#ifdef __GNUC__
    x = x >? y;
#else
    if(y>x) x=y;
#endif
  }
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  void update_min(scalar_type&x, const scalar_type&y)
  { 
#ifdef __GNUC__
    x = x <? y;
#else
    if(y<x) x=y;
#endif
  }
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  void update_max(scalar_type&x, const scalar_type&y, const scalar_type&a)
  {
#ifdef __GNUC__
    x = x >? (y+a);
#else
    if((y+a)>x) x=y+a;
#endif
  }
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  void update_min(scalar_type&x, const scalar_type&y, const scalar_type&a)
  { 
#ifdef __GNUC__
    x = x <? (y-a);
#else
    if((y-a)<x) x=y-a;
#endif
  }
  //----------------------------------------------------------------------------
  template<typename scalar_type> inline
  char* negspace(const scalar_type&x)
  { return (x<scalar_type(0))? " " : "  "; }
  //----------------------------------------------------------------------------
  inline bool is_integral(const float&x)
  { return floor(abs(x))==abs(x);  }
  //----------------------------------------------------------------------------
  inline bool is_integral(const double&x)
  { return floor(abs(x))==abs(x);  }
}
//------------------------------------------------------------------------------
#endif // falcON_included_inln_h
