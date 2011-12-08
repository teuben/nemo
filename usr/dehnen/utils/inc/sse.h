// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/sse.h
///
/// \brief  support for SSE coding and simple SSE supported code
///
/// \author Walter Dehnen
///
/// \date   2009-2011
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009-2011 Walter Dehnen
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
#ifndef WDutils_included_sse_h
#define WDutils_included_sse_h

#ifdef __SSE__
#  ifndef WDutils_included_xmmintrin_h
#    define WDutils_included_xmmintrin_h

#    if defined(__INTEL_COMPILER) && \
        defined(__GNUC__) && \
        defined(_MM_MALLOC_H_INCLUDED)
#      warning
#      warning The intel compiler has seen GNU's _mm_malloc.h which declares _mm_malloc() and _mm_free() to have different linking than those declared in INTEL's xmmintrin.h header file, which we are going to include now. This may cause a compiler error, which can be prevented by ensuring that _mm_malloc.h is not explicitly included when using the intel compiler.
#      warning
#    endif // __INTEL_COMPILER etc

extern "C" {
#    include <xmmintrin.h>
}
#  endif // WDutils_included_xmmintrin_h

# ifdef __SSE2__
#  ifdef __SSE4_1__
#   ifndef WDutils_included_smmintrin_h
#    define WDutils_included_smmintrin_h
extern "C" {
#    include <smmintrin.h>
}
#   endif // WDutils_included_smmintrin_h
#  else
#   ifndef WDutils_included_emmintrin_h
#    define WDutils_included_emmintrin_h
extern "C" {
#    include <emmintrin.h>
}
#   endif // WDutils_included_emmintrin_h
#  endif  // __SSE4_1__
//
// macros for |x|, -x, -|x|, |x-y|, and sign(x)*|y| of packed double
//
#  define _mm_abs_pd(__x)						\
    _mm_and_pd(__x,(__m128d)_mm_set_epi32(0x7fffffff,0xffffffff,	\
					  0x7fffffff,0xffffffff))
#  define _mm_neg_pd(__x)						\
    _mm_xor_pd(__x,(__m128d)_mm_set_epi32(0x80000000,0x0,0x80000000,0x0))
#  define _mm_nabs_pd(__x)						\
    _mm_or_pd (__x,(__m128d)_mm_set_epi32(0x80000000,0x0,0x80000000,0x0))
#  define _mm_diff_pd(__x,__y) _mm_abs_pd(_mm_sub_pd(__x,__y))
#  define _mm_signmask_pd(__x)						\
  _mm_and_pd (__x,(__m128d)_mm_set_epi32(0x80000000,0x0,0x80000000,0x0))
#  define _mm_signmove_pd(__x,__y)					\
  _mm_or_pd(_mm_signmask_pd(__x),_mm_abs_pd(__y))

//
// conversion:  two __m128d  [D,D],[D,D]  -->  one __m128  [S,S,S,S]
// 
# define _mm_cvt2pd_ps(__A,__B) \
  _mm_movelh_ps(_mm_cvtpd_ps(__A),_mm_cvtpd_ps(__B))

//
// macros for |x|, -x, -|x|, |x-y|, and sign(x)*|y| of packed single 
//
#  define _mm_abs_ps(__x) _mm_and_ps(__x,(__m128)_mm_set1_epi32(0x7fffffff))
#  define _mm_neg_ps(__x) _mm_xor_ps(__x,(__m128)_mm_set1_epi32(0x80000000))
#  define _mm_nabs_ps(__x) _mm_or_ps(__x,(__m128)_mm_set1_epi32(0x80000000))
#  define _mm_diff_ps(__x,__y) _mm_abs_ps(_mm_sub_ps(__x,__y))
#  define _mm_signmask_ps(__x)				\
  _mm_and_ps(__x,(__m128)_mm_set1_epi32(0x80000000))
#  define _mm_signmove_ps(__x,__y)			\
  _mm_or_ps(_mm_signmask_ps(__x),_mm_abs_ps(__y))
# else // __SSE2__
// in case we have no __SSE2__, we cannot use _mm_set1_epi32
namespace WDutils {
  namespace meta {
    union FandI {
      float __F;
      int   __I;
      FandI(int i) : __I(i) {}
    };
    const FandI __val_mask(0x7fffffff);
    const FandI __sgn_mask(0x80000000);
  }
}
#  define _mm_abs_ps(__x)					\
  _mm_and_ps(__x,_mm_set1_ps(WDutils::meta::__val_mask.__F))
#  define _mm_neg_ps(__x)					\
  _mm_xor_ps(__x,_mm_set1_ps(WDutils::meta::__sgn_mask.__F))
#  define _mm_nabs_ps(__x)					\
  _mm_or_ps(__x,_mm_set1_ps(WDutils::meta::__sgn_mask.__F))
#  define _mm_diff_ps(__x,__y) _mm_abs_ps(_mm_sub_ps(__x,__y))
#  define _mm_signmask_ps(__x)					\
  _mm_and_ps(__x,_mm_set1_ps(WDutils::meta::__neg_mask.__F))
#  define _mm_signmove_ps(__x,__y)				\
  _mm_or_ps(_mm_signmask_ps(__x),_mm_abs_ps(__y))
# endif // __SSE2__

#ifdef __INTEL_COMPILER
  inline __m128&operator*=(__m128&x, __m128 const&y)
  { return x=_mm_mul_ps(x,y); }
  inline __m128&operator/=(__m128&x, __m128 const&y)
  { return x=_mm_div_ps(x,y); }
  inline __m128&operator+=(__m128&x, __m128 const&y)
  { return x=_mm_add_ps(x,y); }
  inline __m128&operator-=(__m128&x, __m128 const&y)
  { return x=_mm_sub_ps(x,y); }
# ifdef __SSE2__
  inline __m128i&operator+=(__m128i&x, __m128i const&y)
  { return x=_mm_add_epi32(x,y); }
#endif // __SSE2__
  inline float xmm0(__m128 __A)
  {
    union { float f; int i; } tmp;
    tmp.i = _mm_extract_ps(__A,0);
    return tmp.f;
  }
  inline float xmm1(__m128 __A)
  {
    union { float f; int i; } tmp;
    tmp.i = _mm_extract_ps(__A,1);
    return tmp.f;
  }
  inline float xmm2(__m128 __A)
  {
    union { float f; int i; } tmp;
    tmp.i = _mm_extract_ps(__A,2);
    return tmp.f;
  }
  inline float xmm3(__m128 __A)
  {
    union { float f; int i; } tmp;
    tmp.i = _mm_extract_ps(__A,3);
    return tmp.f;
  }
#elif defined(__GNUC__)  // __INTEL_COMPILER / __GNUC__
  inline float xmm0(__m128 __A)
  { return __builtin_ia32_vec_ext_v4sf(__A,0); }
  inline float xmm1(__m128 __A)
  { return __builtin_ia32_vec_ext_v4sf(__A,1); }
  inline float xmm2(__m128 __A)
  { return __builtin_ia32_vec_ext_v4sf(__A,2); }
  inline float xmm3(__m128 __A)
  { return __builtin_ia32_vec_ext_v4sf(__A,3); }
#endif // __INTEL_COMPILER / __GNUC__
#endif // __SSE__

//

#ifndef WDutils_included_meta_h
#  include <meta.h>
#endif
#ifndef WDutils_included_memory_h
#  include <memory.h>
#endif
#ifndef WDutils_included_exception_h
#  include <exception.h>
#endif
#ifndef WDutils_included_cstring
#  define WDutils_included_cstring
#  include <cstring>
#endif

namespace WDutils {

#ifdef __SSE__
  /// horizontal sum: A0+A1+A2+A3.
  /// compute in each SPFP value the sum of the four SPFP values from argument:
  /// result: [A0+A1+A2+A3, A0+A1+A2+A3, A0+A1+A2+A3, A0+A1+A2+A3]
  inline __m128 _mm_sum_ps(__m128 __A)
  {
#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
    __v4sf __a = (__v4sf)__A, __t,__s;
    __t  = __builtin_ia32_shufps(__a,__a,_MM_SHUFFLE (2,3,0,1));
    __s  = __builtin_ia32_addps(__a,__t);
    __t  = __builtin_ia32_shufps(__s,__s,_MM_SHUFFLE (1,0,3,2));
    return (__m128) __builtin_ia32_addps(__s,__t);
#else // gcc
    __m128 __S;
    __S =  _mm_add_ps(__A,_mm_shuffle_ps(__A,__A,_MM_SHUFFLE (2,3,0,1)));
    return _mm_add_ps(__S,_mm_shuffle_ps(__S,__S,_MM_SHUFFLE (1,0,3,2)));
#endif// gcc
  }

  /// horizontal sum: A0+A1+A2+A3.
  /// extract the sum of the four SPFP values from argument
  inline float _mm_getsum_ps(__m128 __A)
  {
#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
  __v4sf __a = (__v4sf)__A, __t,__s;
  __t  = __builtin_ia32_shufps(__a,__a,_MM_SHUFFLE (2,3,0,1));
  __s  = __builtin_ia32_addps(__a,__t);
  __t  = __builtin_ia32_movhlps(__s,__s);
  __s  = __builtin_ia32_addps(__s,__t);
  return __builtin_ia32_vec_ext_v4sf(__s,0);
#else // gcc
    __m128 __S;
    __S =  _mm_add_ps(__A,_mm_shuffle_ps(__A,__A,_MM_SHUFFLE (2,3,0,1)));
    return _mm_cvtss_f32(_mm_add_ps(__S,_mm_movehl_ps(__S,__S)));
#endif// gcc
  }
  /// horizontal maximum: max(A0,A1,A2,A3)
  /// extract the max of the four SPFP values from argument
  inline float _mm_getmax_ps(__m128 __A)
  {
#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
  __v4sf __a = (__v4sf)__A, __t,__s;
  __t  = __builtin_ia32_shufps(__a,__a,_MM_SHUFFLE (2,3,0,1));
  __s  = __builtin_ia32_maxps(__a,__t);
  __t  = __builtin_ia32_movhlps(__s,__s);
  __s  = __builtin_ia32_maxps(__s,__t);
  return __builtin_ia32_vec_ext_v4sf(__s,0);
#else // gcc
    __m128 __S;
    __S =  _mm_max_ps(__A,_mm_shuffle_ps(__A,__A,_MM_SHUFFLE (2,3,0,1)));
    return _mm_cvtss_f32(_mm_max_ps(__S,_mm_movehl_ps(__S,__S)));
#endif// gcc
  }

  /// horizontal minimum: min(A0,A1,A2,A3)
  /// extract the min of the four SPFP values from argument
  inline float _mm_getmin_ps(__m128 __A)
  {
#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
  __v4sf __a = (__v4sf)__A, __t,__s;
  __t  = __builtin_ia32_shufps(__a,__a,_MM_SHUFFLE (2,3,0,1));
  __s  = __builtin_ia32_minps(__a,__t);
  __t  = __builtin_ia32_movhlps(__s,__s);
  __s  = __builtin_ia32_minps(__s,__t);
  return __builtin_ia32_vec_ext_v4sf(__s,0);
#else // gcc
    __m128 __S;
    __S =  _mm_min_ps(__A,_mm_shuffle_ps(__A,__A,_MM_SHUFFLE (2,3,0,1)));
    return _mm_cvtss_f32(_mm_min_ps(__S,_mm_movehl_ps(__S,__S)));
#endif// gcc
  }

#ifdef __SSE2__
  /// horizontal minimum: min(A0,A1,A2,A3)
  /// extract the min of the four packed 32-bit integer values from argument
  inline int32 _mm_getmin_epi32(__m128i __A)
  {
    union WDutils__align16 {
      int   i[4];
      float x[4];
    } tmp;
    _mm_store_ps(tmp.x,(__m128)__A);
    int m=tmp.i[0];
    if(tmp.i[1]<m) m=tmp.i[1];
    if(tmp.i[2]<m) m=tmp.i[2];
    if(tmp.i[3]<m) m=tmp.i[3];
    return m;
  }
  /// horizontal sum: A0+A1+A2+A3
  /// extract the sum of the four packed 32-bit integer values from argument
  inline int32 _mm_getsum_epi32(__m128i __A)
  {
    __m128i __S;
    __S =  _mm_add_epi32(__A,_mm_shuffle_epi32(__A,_MM_SHUFFLE (2,3,0,1)));
    union { int32 i; float x; } tmp;
    tmp.x = _mm_cvtss_f32((__m128)_mm_add_epi32(__S,(__m128i)_mm_movehl_ps((__m128)__S,(__m128)__S)));
    return tmp.i;
  }
#endif// __SSE2__
#endif// __SSE__

  /// support for coding with SSE intrinsics, code using SSE intrinsics.
  namespace SSE {

    template<int N> struct __Aux;
    template<> struct __Aux<2> {
      template<typename __I> static __I top(__I i) { return (i+1) & ~1; }
      template<typename __I> static __I bot(__I i) { return i     & ~1; }
    };
    template<> struct __Aux<4> {
      template<typename __I> static __I top(__I i) { return (i+3) & ~3; }
      template<typename __I> static __I bot(__I i) { return i     & ~3; }
    };
    template<> struct __Aux<8> {
      template<typename __I> static __I top(__I i) { return (i+7) & ~7; }
      template<typename __I> static __I bot(__I i) { return i     & ~7; }
    };
    template<> struct __Aux<16> {
      template<typename __I> static __I top(__I i) { return (i+15)& ~15; }
      template<typename __I> static __I bot(__I i) { return i     & ~15; }
    };

    /// smallest multiple of N not less than i
    template<int N, typename __I> inline
    __I top(__I i) { return __Aux<N>::top(i); }
    /// largest multiple of N not greater than i
    template<int N, typename __I> inline
    __I bottom(__I i) { return __Aux<N>::bot(i); }

    ///
    /// simple array manipulations for unaligned arrays
    ///
    /// \note routines are up to 16/sizeof(T) times faster than simple code
    ///
    struct UnAligned {
      /// \name assign to each element of array
      //@{
      /// \code for(i=0; i!=n; ++i) f[i]=x; \endcode
      static void Ass(int*f, size_t n, int x);
      /// \code for(i=0; i!=n; ++i) f[i]=x; \endcode
      static void Ass(float*f, size_t n, float x);
      /// \code for(i=0; i!=n; ++i) f[i]=x; \endcode
      static void Ass(double*f, size_t n, double x);
      //
      /// \code for(size_t i=0; i!=n; ++i) a[i] = b[i]; \endcode
      template<typename T>
      static void Ass(T*a, size_t n, const T*b)
      { std::memcpy(a,b,n*sizeof(T)); }
      //
      /// \code for(i=0; i!=n; ++i) f[i]=-f[i]; \endcode
      static void Neg(int*f, size_t n);
      /// \code for(i=0; i!=n; ++i) f[i]=-f[i]; \endcode
      static void Neg(float*f, size_t n);
      /// \code for(i=0; i!=n; ++i) f[i]=-f[i]; \endcode
      static void Neg(double*f, size_t n);
      //
      /// \code for(i=0; i!=n; ++i) f[i]=0; \endcode
      template<typename T>
      static void Reset(T*f, size_t n)
      { Ass(f,n,T(0)); }
      //
      /// \code for(i=0; i!=n; ++i) f[i]+=x; \endcode
      static void Add(int*f, size_t n, int x);
      /// \code for(i=0; i!=n; ++i) f[i]+=x; \endcode
      static void Add(float*f, size_t n, float x);
      /// \code for(i=0; i!=n; ++i) f[i]+=x; \endcode
      static void Add(double*f, size_t n, double x);
      //
      /// \code for(i=0; i!=n; ++i) f[i]-=x; \endcode
      static void Sub(int*f, size_t n, int x);
      /// \code for(i=0; i!=n; ++i) f[i]-=x; \endcode
      static void Sub(float*f, size_t n, float x);
      /// \code for(i=0; i!=n; ++i) f[i]-=x; \endcode
      static void Sub(double*f, size_t n, double x);
      //
      /// \code for(i=0; i!=n; ++i) f[i]*=x; \endcode
      static void Mul(int*f, size_t n, int x);
      /// \code for(i=0; i!=n; ++i) f[i]*=x; \endcode
      static void Mul(float*f, size_t n, float x);
      /// \code for(i=0; i!=n; ++i) f[i]*=x; \endcode
      static void Mul(double*f, size_t n, double x);
      //
      /// \code for(i=0; i!=n; ++i) f[i]/=x; \endcode
      static void Div(float*f, size_t n, float x)
      { Mul(f,n,float(1.0/x)); }
      /// \code for(i=0; i!=n; ++i) f[i]/=x; \endcode
      static void Div(double*f, size_t n, double x)
      { Mul(f,n,1.0/x); }
      //
      /// \code for(i=0; i!=n; ++i) f[i]=x/f[i]; \endcode
      static void Inv(float*f, size_t n, float x);
      /// \code for(i=0; i!=n; ++i) f[i]=x/f[i]; \endcode
      static void Inv(double*f, size_t n, double x);
      //
      /// \code for(i=0; i!=n; ++i) f[i]=1/f[i]; \endcode
      template<typename T>
      static void Reciprocal(T*f, size_t n)
      { Inv(f,n,T(1)); }
      //
      /// \code for(i=0; i!=n; ++i) f[i]=std::sqrt(f[i]); \endcode
      static void Sqrt(float*f, size_t n);
      /// \code for(i=0; i!=n; ++i) f[i]=std::sqrt(f[i]); \endcode
      static void Sqrt(double*f, size_t n);
      //@}
      /// \name compute property of whole array
      //{@
      /// \code S=0; for(i=0; i!=n; ++i) S+=f[i]; return S; \endcode
      static int Sum(const int*f, size_t n);
      /// \code S=0; for(i=0; i!=n; ++i) S+=f[i]; return S; \endcode
      static float Sum(const float*f, size_t n);
      /// \code S=0; for(i=0; i!=n; ++i) S+=f[i]; return S; \endcode
      static double Sum(const double*f, size_t n);
      //
      /// \code S=0; for(i=0; i!=n; ++i) S+=f[i]*f[i]; return S; \endcode
      static int Norm(const int*f, size_t n);
      /// \code S=0; for(i=0; i!=n; ++i) S+=f[i]*f[i]; return S; \endcode
      static float Norm(const float*f, size_t n);
      /// \code S=0; for(i=0; i!=n; ++i) S+=f[i]*f[i]; return S; \endcode
      static double Norm(const double*f, size_t n);
      //@}
    };// class SSE::UnAligned

    ///
    /// simple array manipulations for aligned arrays
    ///
    /// \note all array arguments must be 16-byte aligned and array sizes
    ///       multiples of 16/sizeof(T) where T is the array type.
    /// \note routines are about 16/sizeof(T) times faster than simple code
    ///
    struct Aligned {
      /// \name assign to each element of array
      //@{
      /// \code for(i=0; i!=n; ++i) f[i]=x; \endcode
      static void Ass(int*f, size_t n, int x);
      /// \code for(i=0; i!=n; ++i) f[i]=x; \endcode
      static void Ass(float*f, size_t n, float x);
      /// \code for(i=0; i!=n; ++i) f[i]=x; \endcode
      static void Ass(double*f, size_t n, double x);
      //
      /// \code for(i=0; i!=n; ++i) a[i]=b[i]; \endcode
      /// \note this does not require 16-byte alignment
      template<typename T>
      static void Ass(T*a, size_t n, const T*b)
      { std::memcpy(a,b,n*sizeof(T)); }
      //
      /// \code for(i=0; i!=n; ++i) a[i]=w*b[i]; \endcode
      static void Ass(int*a, size_t n, int w, const int*b);
      /// \code for(i=0; i!=n; ++i) a[i]=w*b[i]; \endcode
      static void Ass(float*a, size_t n, float w, const float*b);
      /// \code for(i=0; i!=n; ++i) a[i]=w*b[i]; \endcode
      static void Ass(double*a, size_t n, double w, const double*b);
      //
      /// \code for(i=0; i!=n; ++i) f[i]=-f[i]; \endcode
      static void Neg(int*f, size_t n);
      /// \code for(i=0; i!=n; ++i) f[i]=-f[i]; \endcode
      static void Neg(float*f, size_t n);
      /// \code for(i=0; i!=n; ++i) f[i]=-f[i]; \endcode
      static void Neg(double*f, size_t n);
      //
      /// \code for(i=0; i!=n; ++i) f[i]=0; \endcode
      template<typename T>
      static void Reset(T*f, size_t n)
      { Ass(f,n,T(0)); }
      //
      /// \code for(i=0; i!=n; ++i) f[i]+=x; \endcode
      static void Add(int*f, size_t n, int x);
      /// \code for(i=0; i!=n; ++i) f[i]+=x; \endcode
      static void Add(float*f, size_t n, float x);
      /// \code for(i=0; i!=n; ++i) f[i]+=x; \endcode
      static void Add(double*f, size_t n, double x);
      //
      /// \code for(i=0; i!=n; ++i) a[i]+=b[i]; \endcode
      static void Add(int*a, size_t n, const int*b);
      /// \code for(i=0; i!=n; ++i) a[i]+=b[i]; \endcode
      static void Add(float*a, size_t n, const float*b);
      /// \code for(i=0; i!=n; ++i) a[i]+=b[i]; \endcode
      static void Add(double*a, size_t n, const double*b);
      //
      /// \code for(i=0; i!=n; ++i) a[i]+=w*b[i]; \endcode
      static void Add(int*a, size_t n, int w, const int*b);
      /// \code for(i=0; i!=n; ++i) a[i]+=w*b[i]; \endcode
      static void Add(float*a, size_t n, float w, const float*b);
      /// \code for(i=0; i!=n; ++i) a[i]+=w*b[i]; \endcode
      static void Add(double*a, size_t n, double w, const double*b);
      //
      /// \code for(i=0; i!=n; ++i) f[i]-=x; \endcode
      static void Sub(int*f, size_t n, int x);
      /// \code for(i=0; i!=n; ++i) f[i]-=x; \endcode
      static void Sub(float*f, size_t n, float x);
      /// \code for(i=0; i!=n; ++i) f[i]-=x; \endcode
      static void Sub(double*f, size_t n, double x);
      //
      /// \code for(i=0; i!=n; ++i) a[i]-=b[i]; \endcode
      static void Sub(int*a, size_t n, const int*b);
      /// \code for(i=0; i!=n; ++i) a[i]-=b[i]; \endcode
      static void Sub(float*a, size_t n, const float*b);
      /// \code for(i=0; i!=n; ++i) a[i]-=b[i]; \endcode
      static void Sub(double*a, size_t n, const double*b);
      //
      /// \code for(i=0; i!=n; ++i) a[i]-=w*b[i]; \endcode
      static void Sub(int*a, size_t n, int w, const int*b);
      /// \code for(i=0; i!=n; ++i) a[i]-=w*b[i]; \endcode
      static void Sub(float*a, size_t n, float w, const float*b);
      /// \code for(i=0; i!=n; ++i) a[i]-=w*b[i]; \endcode
      static void Sub(double*a, size_t n, double w, const double*b);
      //
      /// \code for(i=0; i!=n; ++i) f[i]*=x; \endcode
      static void Mul(int*f, size_t n, int x);
      /// \code for(i=0; i!=n; ++i) f[i]*=x; \endcode
      static void Mul(float*f, size_t n, float x);
      /// \code for(i=0; i!=n; ++i) f[i]*=x; \endcode
      static void Mul(double*f, size_t n, double x);
      //
      /// \code for(i=0; i!=n; ++i) a[i]*=b[i]; \endcode
      static void Mul(int*a, size_t n, const int*b);
      /// \code for(i=0; i!=n; ++i) a[i]*=b[i]; \endcode
      static void Mul(float*a, size_t n, const float*b);
      /// \code for(i=0; i!=n; ++i) a[i]*=b[i]; \endcode
      static void Mul(double*a, size_t n, const double*b);
      //
      /// \code for(i=0; i!=n; ++i) f[i]/=x; \endcode
      static void Div(float*f, size_t n, float x)
      { Mul(f,n,float(1.0/x)); }
      /// \code for(i=0; i!=n; ++i) f[i]/=x; \endcode
      static void Div(double*f, size_t n, double x)
      { Mul(f,n,1.0/x); }
      //
      /// \code for(i=0; i!=n; ++i) a[i]/=b[i]; \endcode
      static void Div(float*a, size_t n, const float*b);
      /// \code for(i=0; i!=n; ++i) a[i]/=b[i]; \endcode
      static void Div(double*a, size_t n, const double*b);
      //
      /// \code for(i=0; i!=n; ++i) f[i]=x/f[i]; \endcode
      static void Inv(float*f, size_t n, float x);
      /// \code for(i=0; i!=n; ++i) f[i]=x/f[i]; \endcode
      static void Inv(double*f, size_t n, double x);
      //
      /// \code for(i=0; i!=n; ++i) f[i]=x/b[i]; \endcode
      static void Inv(float*f, size_t n, float x, const float*b);
      /// \code for(i=0; i!=n; ++i) f[i]=x/b[i]; \endcode
      static void Inv(double*f, size_t n, double x, const double*b);
      //
      /// \code for(i=0; i!=n; ++i) f[i]=1/f[i]; \endcode
      template<typename T>
      static void Reciprocal(T*f, size_t n)
      { Inv(f,n,T(1)); }
      //
      /// \code for(i=0; i!=n; ++i) a[i]=1/b[i]; \endcode
      template<typename T>
      static void Reciprocal(T*a, size_t n, const T*b)
      { Inv(a,n,T(1),b); }
      //
      /// \code for(i=0; i!=n; ++i) f[i]=std::sqrt(f[i]); \endcode
      static void Sqrt(float*f, size_t n);
      /// \code for(i=0; i!=n; ++i) f[i]=std::sqrt(f[i]); \endcode
      static void Sqrt(double*f, size_t n);
      //
      /// \code for(i=0; i!=n; ++i) a[i]=std::sqrt(b[i]); \endcode
      static void Sqrt(float*a, size_t n, const float*b);
      /// \code for(i=0; i!=n; ++i) a[i]=std::sqrt(b[i]); \endcode
      static void Sqrt(double*a, size_t n, const double*b);
      //@}
      /// \name compute property of whole array
      //{@
      /// \code S=0; for(i=0; i!=n; ++i) S+=f[i]; return S; \endcode
      static int Sum(const int*f, size_t n);
      /// \code S=0; for(i=0; i!=n; ++i) S+=f[i]; return S; \endcode
      static float Sum(const float*f, size_t n);
      /// \code S=0; for(i=0; i!=n; ++i) S+=f[i]; return S; \endcode
      static double Sum(const double*f, size_t n);
      //
      /// \code S=0; for(i=0; i!=n; ++i) S+=f[i]*f[i]; return S; \endcode
      static int Norm(const int*f, size_t n);
      /// \code S=0; for(i=0; i!=n; ++i) S+=f[i]*f[i]; return S; \endcode
      static float Norm(const float*f, size_t n);
      /// \code S=0; for(i=0; i!=n; ++i) S+=f[i]*f[i]; return S; \endcode
      static double Norm(const double*f, size_t n);
      //
      /// \code S=0; for(i=0; i!=n; ++i) S+=a[i]*b[i]; return S; \endcode
      static int Dot(const int*a, size_t n, const int*b);
      /// \code S=0; for(i=0; i!=n; ++i) S+=a[i]*b[i]; return S; \endcode
      static float Dot(const float*a, size_t n, const float*b);
      /// \code S=0; for(i=0; i!=n; ++i) S+=a[i]*b[i]; return S; \endcode
      static double Dot(const double*a, size_t n, const double*b);
      //@}
    };// class SSE::Aligned

    /// contains some constants and code relevant for SSE coding
    /// instantinations for float and double
    template<typename __F> struct Traits;

    /// SSE::Traits<float>
    template<> struct Traits<float>
    {
      /// is SSE enabled for this type?
#ifdef __SSE__
      static const bool sse = true;
#else
      static const bool sse = false;
#endif
      /// alignment number: K floats align to 128 bytes
      static const int K=4;
      /// smallest multiple of K not less than n
      template<typename __I>
      static __I Top(__I n) { return top<K,__I>(n); }
      /// largest multiple of K not greater than n
      template<typename __I>
      static __I Bottom(__I n) { return bottom<K,__I>(n); }
      /// is an array of floats aligned to at least 4 bytes (=sizeof(float))?
      static bool is_minimum_aligned(float*f) {
	return (size_t(f) & 4)  == 0;
      }
      /// is an array aligned to 16 bytes
      static bool is_aligned(float*f) {
	return (size_t(f) & 16) == 0;
      }
    };

    /// SSE::Traits<int>
    template<> struct Traits<int>
    {
      /// is SSE enabled for this type?
#ifdef __SSE2__
      static const bool sse = true;
#else
      static const bool sse = false;
#endif
      /// alignment number: K floats align to 128 bytes
      static const int K=4;
      /// smallest multiple of K not less than n
      template<typename __I>
      static __I Top(__I n) { return top<K,__I>(n); }
      /// largest multiple of K not greater than n
      template<typename __I>
      static __I Bottom(__I n) { return bottom<K,__I>(n); }
      /// is an array of ints aligned to at least 4 bytes (=sizeof(int))?
      static bool is_minimum_aligned(int*f) {
	return (size_t(f) & 4)  == 0;
      }
      /// is an array aligned to 16 bytes
      static bool is_aligned(int*f) {
	return (size_t(f) & 16) == 0;
      }
    };

    /// SSE::Traits<double>
    template<> struct Traits<double>
    {
      /// is SSE enabled for this type?
#ifdef __SSE2__
      static const bool sse = true;
#else
      static const bool sse = false;
#endif
      /// alignment number: K doubles align to 128 bytes
      static const int K=2;
      /// smallest multiple of K not less than n
      template<typename __I>
      static __I Top(__I n) { return top<K,__I>(n); }
      /// largest multiple of K not greater than n
      template<typename __I>
      static __I Bottom(__I n) { return bottom<K,__I>(n); }
      /// is an array of floats aligned to at least 8 bytes (=sizeof(double))?
      static bool is_minimum_aligned(double*f) {
	return (size_t(f) & 8)  == 0;
      }
      /// is an array aligned to 16 bytes?
      static bool is_aligned(double*f) {
	return (size_t(f) & 16) == 0;
      }
    };
    /// smallest multiple of 16/sizeof(__F) not less than i
    template<typename __F, typename __I> inline
    __I Top(__I i) {
      return Traits<__F>::Top(i);
    }
    /// largest multiple of 16/sizeof(__F) not greater than i
    template<typename __F, typename __I> inline
    __I Bottom(__I i) {
      return Traits<__F>::Bottom(i);
    }
    ////////////////////////////////////////////////////////////////////////////
    /// An array of float or double supporting operations via SSE instructions
    /// \note not to be confused with WDutils::Array16 in memory.h
    template<typename _F>
    class Array16 {
      WDutilsStaticAssert( WDutilsSameType(_F,int) ||
			   WDutilsSameType(_F,float) ||
			   WDutilsSameType(_F,double) );
      /// copy ctor disabled to encourage references as return type etc.
      /// \note you may use  member \a assign() instead
      Array16(Array16 const&);
      /// \name data
      //@{
      const size_t _S; ///< # elements actually allocated
    protected:
      /// # elements requested
      /// \note under GCC we cannot use __N as this is #defined in
      ///       /usr/include/c++/4.3/x86_64-suse-linux/bits/c++config.h
      const size_t _N;
      _F    *const _A; ///< array with elements
      //@}
      /// resets elements with _N <= i < _S
      void reset_tail() const
      { for(size_t i=_N; i!=_S; ++i) _A[i] = _F(0); }
      /// check for size mismatch
      void check_size(Array16 const&B, const char*name) const WDutils_THROWING
      {
	if(B._N != _N)
	  WDutils_THROW("SSE::Array16<%s>::%s: size mismatch:%u vs %u\n",
			nameof(_F),name, unsigned(_N), unsigned(B._N));
      }
    public:
      /// default ctor
      Array16()
	: _S(0), _N(0), _A(0) {}
      /// ctor from given size
      explicit Array16(size_t n)
	: _S(Top<_F>(n)), _N(n), _A(WDutils_NEW16(_F,_S)) {}
      /// dtor
      ~Array16()
      { if(_A) WDutils_DEL16(_A); const_cast<_F*&>(_A)=0; }
      /// reset size
      /// \param[in] n new number of elements
      Array16&reset(size_t n);
      /// assign to another Array16
      /// \param[in] a another Array16 to become a copy of
      /// \note if sizes don't match, we reset ours first
      Array16&assign(Array16 const&a);
      /// number of elements
      size_t const&size() const { return _N; }
      /// const access
      _F const&operator[] (int i) const { return _A[i]; }
      /// non-const access
      _F&operator[] (int i) { return _A[i]; }
      /// direct access to data array
      _F*array() { return _A; }
      /// direct access to data array
      const _F*array() const { return _A; }
      /// \name unary operations
      //@{
      /// reset element-wise to 0
      /// \code for(int i=0; i!=size(); ++i) A[i] = 0; \endcode
      Array16&reset()
      { Aligned::Reset(_A,_S); return*this; }
      /// negate element-wise
      /// \code for(int i=0; i!=size(); ++i) A[i] = -A[i]; \endcode
      Array16&negate()
      { Aligned::Neg(_A,_S); return*this; }
      /// sum of all elements
      /// \code _F x(0); for(int i=0; i!=size(); ++i) x+=A[i]; return x;
      /// \endcode
      _F sum() const
      { reset_tail(); return Aligned::Sum(_A,_S); }
      /// sum of all elements squared
      /// \code _F x(0); for(int i=0; i!=size(); ++i) x+=A[i]*A[i]; return x;
      /// \endcode
      _F norm() const
      { reset_tail(); return Aligned::Norm(_A,_S); }
      //@}
      /// \name binary operations with scalar
      //@{
      /// assign element-wise to scalar
      /// \code for(int i=0; i!=size(); ++i) A[i] = x; \endcode
      Array16&operator=(_F x)
      { Aligned::Ass(_A,_S,x); return*this; }
      /// add a scalar to each element
      /// \code for(int i=0; i!=size(); ++i) A[i] += x; \endcode
      Array16&operator+=(_F x)
      { Aligned::Add(_A,_S,x); return*this; }
      /// subtract a scalar from each element
      /// \code for(int i=0; i!=size(); ++i) A[i] -= x; \endcode
      Array16&operator-=(_F x)
      { Aligned::Sub(_A,_S,x); return*this; }
      /// multiply each element by a scalar
      /// \code for(int i=0; i!=size(); ++i) A[i] *= x; \endcode
      Array16&operator*=(_F x)
      { Aligned::Mul(_A,_S,x); return*this; }
      /// divide each element by a scalar
      /// \code for(int i=0; i!=size(); ++i) A[i] *= x; \endcode
      Array16&operator/=(_F x) WDutils_THROWING
      { Aligned::Div(_A,_S,x); return*this; }
      //@}
      /// \name binary operations with another Array16
      //@{
      /// assign element-wise
      /// \code for(int i=0; i!=size(); ++i) A[i] = B[i]; \endcode
      /// \note it is an error if the number of elements do not match
      Array16&operator=(Array16 const&B) WDutils_THROWING;
      /// add element-wise
      /// \code for(int i=0; i!=size(); ++i) A[i] += B[i]; \endcode
      /// \note it is an error if the number of elements do not match
      Array16&operator+=(Array16 const&B) WDutils_THROWING
      {
	check_size(B,"operator+=(Array16&)");
	Aligned::Add(_A,_S,B._A);
	return*this;
      }
      /// subtract element-wise
      /// \code for(int i=0; i!=size(); ++i) A[i] += B[i]; \endcode
      /// \note it is an error if the number of elements do not match
      Array16&operator-=(Array16 const&B) WDutils_THROWING
      {
	check_size(B,"operator-=(Array16&)");
	Aligned::Sub(_A,_S,B._A);
	return*this;
      }
      /// dot product
      /// \code _F x(0); for(int i=0; i!=size(); ++i) x+=A[i]*B[i]; return x;
      /// \endcode
      /// \note it is an error if the number of elements do not match
      _F operator*(Array16 const&B) const WDutils_THROWING
      {
	check_size(B,"operator*(Array16&)");
	reset_tail();
	return Aligned::Dot(_A,_S,B._A);
      }
      //@}
      /// \name tertiary operations with scalar and another Array16
      //@{
      /// add weighted element-wise
      /// \code for(int i=0; i!=size(); ++i) A[i] += w*B[i]; \endcode
      Array16&addtimes(Array16 const&B, _F w) WDutils_THROWING
      {
	check_size(B,"addtimes(Array16&)");
	Aligned::Add(_A,_S,w,B._A);
	return*this;
      }
      /// subtract weighted element-wise
      /// \code for(int i=0; i!=size(); ++i) A[i] += w*B[i]; \endcode
      Array16&subtimes(Array16 const&B, _F w) WDutils_THROWING
      {
	check_size(B,"subtimes(Array16&)");
	Aligned::Sub(_A,_S,w,B._A);
	return*this;
      }
      //@}
    };
    /// filler object of size __S
    template<size_t __S> class Filler { char __F[__S]; };
    /// extend a given class to have size a multibple of 16-byte
    /// \note only default constructor possible for @a Base
    template<typename Base>
    class Extend16 : public Base
    {
      static const size_t Bytes = sizeof(Base);
      static const size_t Splus = Bytes & 15;
      static const size_t Added = Splus? 16-Splus : 0;
      Filler<Added> __Filler;
    };
    //
#ifdef __SSE__
    /// working horse actually implemented in sse.cc
    /// \note Not intended for consumption, use routine below instead.
    void __swap16(float*a, float*b, size_t n);
    /// swap a multiple of 16 bytes at 16-byte aligned memory.
    /// \param[in] a   16-byte aligned memory address
    /// \param[in] b   16-byte aligned memory address
    /// \param[in] n   multiple of 16: # bytes to swap at @a a and @a b
    /// \note For reasons of efficiency, the above conditions are not tested.
    ///       If either address is not 16-byte aligned, a run-time error
    ///       (segmentation fault) will result. If @a n is not a multiple of
    ///       16, the last @a n%16 bytes will not be swapped.
    inline void Swap16(void*a, void*b, size_t n)
    { __swap16(static_cast<float*>(a),static_cast<float*>(b),n>>2); }
    /// swaps 16-byte aligned objects with size a multiple of 16-bytes
    /// \param[in,out] x  pter to object, on return holds data of @a y
    /// \param[in,out] y  pter to object, on return holds data of @a x
    /// \note If the sizeof(Object) is not a multiple of 16, a compile-time
    ///       error occurs. However, if either @a x or @a y is not 16-byte
    ///       aligned a run-time error (segmentation fault) will result.
    template<typename Object>
    inline void SwapAligned(Object*a, Object*b)
    {
      WDutilsStaticAssert( sizeof(Object)%16 == 0 );
      __swap16(reinterpret_cast<float*>(a), reinterpret_cast<float*>(b),
	       sizeof(Object)>>2);
    }
#endif
  } // namespace SSE
} // namespace WDutils
//
#endif
