// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/sse.h
///
/// \brief  support for SSE coding and simple SSE supported code
///
/// \author Walter Dehnen
///
/// \date   2009-2012
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009-2012 Walter Dehnen
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

#if __cplusplus >= 201103L
# ifndef WDutils_included_type_traits
#  include <type_traits>
#  define WDutils_included_type_traits
# endif
#endif

#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
// with GCC use x86intrin.h
# ifndef _X86INTRIN_H_INCLUDED
extern "C" {
#  include <x86intrin.h>
}
# endif

#elif defined(__SSE__)
// use individual headers for intrinsics

# ifdef defined(__INTEL_COMPILER) && defined(_MM_MALLOC_H_INCLUDED)
# warning The intel compiler has seen GNU's _mm_malloc.h which declares _mm_malloc() and _mm_free() to have different linking than those declared in INTEL's xmmintrin.h header file, which we are going to include now. This may cause a compiler error, which can be prevented by ensuring that _mm_malloc.h is not explicitly included when using the intel compiler.
# endif

# ifndef WDutils_included_xmmintrin_h
extern "C" {
#    include <xmmintrin.h>
}
# endif // WDutils_included_xmmintrin_h

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
# endif   // __SSE2__
#endif    // __SSE__

#if __cplusplus < 201103L
# define noexcept
# define constexpr
#endif

//
#ifdef __SSE__
namespace WDutils {
#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
# define always_inline __attribute__((__always_inline__))
#else
# define always_inline
#endif
  ///
  /// \name SSE & AVX functions and operations
  //@{
  namespace meta {
    union float_and_int {
      float f; int32_t i; 
      float_and_int()
# if __cplusplus >= 201103L
      = default;
# else
      : {}
# endif
      float_and_int(int32_t k) : i(k) {}
    };
  }
# ifdef __SSE2__
#  define abs_mask_pi _mm_set1_epi32(0x7fffffff)
#  define sgn_mask_pi _mm_set1_epi32(0x80000000)
#  define abs_mask_ps _mm_castsi128_ps(abs_mask_pi)
#  define sgn_mask_ps _mm_castsi128_ps(sgn_mask_pi)
#  define abs_mask_pd _mm_castsi128_pd(_mm_set1_epi64x(0x7fffffffffffffffll))
#  define sgn_mask_pd _mm_castsi128_pd(_mm_set1_epi64x(0x8000000000000000ll))
# else
#  define abs_mask_ps _mm_set1_ps(WDutils::meta::float_and_int(0x7fffffff).f)
#  define sgn_mask_ps _mm_set1_ps(WDutils::meta::float_and_int(0x80000000).f)
# endif// __SSE2__
# ifdef __AVX__
#  define abs_mask_qi _mm256_set1_epi32(0x7fffffff)
#  define sgn_mask_qi _mm256_set1_epi32(0x80000000)
#  define abs_mask_qs _mm256_castsi256_ps(abs_mask_qi)
#  define sgn_mask_qs _mm256_castsi256_ps(sgn_mask_qi )
#  define abs_mask_qd						\
  _mm256_castsi256_pd(_mm256_set1_epi64x(0x7fffffffffffffffll))
#  define sgn_mask_qd						\
  _mm256_castsi256_pd(_mm256_set1_epi64x(0x8000000000000000ll))
# endif// __AVX__
#endif // __SSE__
  //
  ///
  /// generic support for coding with SSE/AVX intrinsics
  ///
  namespace SSE {
#if __cplusplus >= 201103L
    /// is @c VectorType a SSE or AVX floating-point vector type?
    template<typename VectorType>
    struct is_floating_vector_type
    {
      static const bool value = 0
# ifdef __SSE__
	|| std::is_same<VectorType,__m128>::value
# endif
# ifdef __SSE2__
	|| std::is_same<VectorType,__m128d>::value
# endif
# ifdef __AVX__
	|| std::is_same<VectorType,__m256>::value
	|| std::is_same<VectorType,__m256d>::value
# endif
	;
    };
    /// is @c VectorType a SSE or AVX integer vector type?
    template<typename VectorType>
    struct is_integer_vector_type
    {
      static const bool value = 0
# ifdef __SSE2__
	|| std::is_same<VectorType,__m128i>::value
# endif
# ifdef __AVX__
	|| std::is_same<VectorType,__m256i>::value
# endif
	;
    };
    /// is @c VectorType a SSE or AVX 32-bit integer or floating-point vector
    /// type?
    template<typename VectorType>
    struct is_supported_vector_type
    {
      static const bool value = 0
# ifdef __SSE__
	|| std::is_same<VectorType,__m128>::value
# endif
# ifdef __SSE2__
	|| std::is_same<VectorType,__m128d>::value
	|| std::is_same<VectorType,__m128i>::value
# endif
# ifdef __AVX__
	|| std::is_same<VectorType,__m256>::value
	|| std::is_same<VectorType,__m256d>::value
	|| std::is_same<VectorType,__m256i>::value
# endif
	;
    };
#endif// C++11
#ifdef __SSE__
    /// static info related to SSE/AVX types
    template<typename VectorType> struct static_type_info;
    /// __m128:  4 packed single-precision floating point numbers
    template<> struct static_type_info<__m128>
    {
      using vector_type = __m128;
      using single_type = float;
      static const int block_size = 4;
      static const int alignment = 16;
      static_type_info() = delete;
    };
# ifdef __SSE2__
    /// __m128d:  2 packed double-precision floating point numbers
    template<> struct static_type_info<__m128d>
    {
      using vector_type = __m128d;
      using single_type = double;
      static const int block_size = 2;
      static const int alignment = 16;
      static_type_info() = delete;
    };
    /// __m128i:  4 packed 32-bit signed integers
    /// \note gcc interprets __m128i differently as 2 int64_t. Our
    ///       interpretation appears the most useful for scientific computing.
    template<> struct static_type_info<__m128i>
    {
      using vector_type = __m128i;
      using single_type = int32_t;
      static const int block_size = 4;
      static const int alignment = 16;
      static_type_info() = delete;
    };
# endif
# ifdef __AVX__
    /// __m256:  8 packed single-precision floating point numbers
    template<> struct static_type_info<__m256>
    {
      using vector_type = __m256;
      using single_type = float;
      static const int block_size = 8;
      static const int alignment = 32;
      static_type_info() = delete;
    };
    /// __m256d:  4 packed double-precision floating point numbers
    template<> struct static_type_info<__m256d>
    {
      using vector_type = __m256d;
      using single_type = double;
      static const int block_size = 4;
      static const int alignment = 32;
      static_type_info() = delete;
    };
    /// __m256d:  8 packed 32-bit signed integers
    /// \note gcc interprets __m256i differently as 4 int64_t. Our
    ///       interpretation appears the most useful for scientific computing.
    template<> struct static_type_info<__m256i>
    {
      using vector_type = __m256i;
      using single_type = int32_t;
      static const int block_size = 8;
      static const int alignment = 32;
      static_type_info() = delete;
    };
# endif
    ///
    /// auxiliary code used to implement functionality in enclosing namespace
    ///
    namespace aux {
      ///
      /// auxiliary template, used to implement some SSE/AVX related functions
      /// which must be templates (because the intended vector type cannot be
      /// uniquely deduced from their arguments)
      ///
      template<typename VectorType> struct VecOps;
      // for __m128
      template<> struct VecOps<__m128> : public static_type_info<__m128>
      {
	static vector_type always_inline zero() noexcept
	{ return _mm_setzero_ps(); }
	static vector_type always_inline one() noexcept
	{ return _mm_set1_ps(1.f); }
	static vector_type always_inline set(const single_type x) noexcept
	{ return _mm_set1_ps(x); }
	static vector_type always_inline load(const single_type*p) noexcept
	{ return _mm_load_ps(p); }
	static vector_type always_inline loadu(const single_type*p) noexcept
	{ return _mm_loadu_ps(p); }
      };
# ifdef __SSE2__
      // for __m128d
      template<> struct VecOps<__m128d> : public static_type_info<__m128d>
      {
	static vector_type always_inline zero() noexcept
	{ return _mm_setzero_pd(); }
	static vector_type always_inline one() noexcept
	{ return _mm_set1_pd(1.0); }
	static vector_type always_inline set(const single_type x) noexcept
	{ return _mm_set1_pd(x); }
	static vector_type always_inline load(const single_type*p) noexcept
	{ return _mm_load_pd(p); }
	static vector_type always_inline loadu(const single_type*p) noexcept
	{ return _mm_loadu_pd(p); }
      };
      // for __m128i
      template<> struct VecOps<__m128i> : public static_type_info<__m128i>
      {
	static vector_type always_inline zero() noexcept
	{ return _mm_setzero_si128(); }
	static vector_type always_inline one() noexcept
	{ return _mm_set1_epi32(1); }
	static vector_type always_inline set(const single_type x) noexcept
	{ return _mm_set1_epi32(x); }
	static vector_type always_inline load(const single_type*p) noexcept
	{ return _mm_castps_si128(_mm_load_ps
				  (reinterpret_cast<const float*>(p))); }
	static vector_type always_inline loadu(const single_type*p) noexcept
	{ return _mm_castps_si128(_mm_loadu_ps
				  (reinterpret_cast<const float*>(p))); }
      };
# endif
# ifdef __AVX__
      // for __m256
      template<> struct VecOps<__m256> : public static_type_info<__m256>
      {
	static vector_type always_inline zero() noexcept
	{ return _mm256_setzero_ps(); }
	static vector_type always_inline one() noexcept
	{ return _mm256_set1_ps(1.f); }
	static vector_type always_inline set(const single_type x) noexcept
	{ return _mm256_set1_ps(x); }
	static vector_type always_inline load(const single_type*p) noexcept
	{ return _mm256_load_ps(p); }
	static vector_type always_inline loadu(const single_type*p) noexcept
	{ return _mm256_loadu_ps(p); }
      };
      // for __m256d
      template<> struct VecOps<__m256d> : public static_type_info<__m256d>
      {
	static vector_type always_inline zero() noexcept
	{ return _mm256_setzero_pd(); }
	static vector_type always_inline one() noexcept
	{ return _mm256_set1_pd(1.0); }
	static vector_type always_inline set(const single_type x) noexcept
	{ return _mm256_set1_pd(x); }
	static vector_type always_inline load(const single_type*p) noexcept
	{ return _mm256_load_pd(p); }
	static vector_type always_inline loadu(const single_type*p) noexcept
	{ return _mm256_loadu_pd(p); }
      };
      // for __m256i
      template<> struct VecOps<__m256i> : public static_type_info<__m256i>
      {
	static vector_type always_inline zero() noexcept
	{ return _mm256_setzero_si256(); }
	static vector_type always_inline one() noexcept
	{ return _mm256_set1_epi32(1); }
	static vector_type always_inline set(const single_type x) noexcept
	{ return _mm256_set1_epi32(x); }
	static vector_type always_inline load(const single_type*p) noexcept
	{ return _mm256_castps_si256(_mm256_load_ps
				     (reinterpret_cast<const float*>(p))); }
	static vector_type always_inline loadu(const single_type*p) noexcept
	{ return _mm256_castps_si256(_mm256_loadu_ps
				     (reinterpret_cast<const float*>(p))); }
     };
# endif
    } // namespace WDutils::SSE::aux
    ///
    /// \name SSE/AVX non-operator functions with unified name
    //@{

    ///
    /// addition
    ///
    inline __m128 always_inline add(__m128 x, __m128 y) noexcept
    { return _mm_add_ps(x,y); }
# ifdef __SSE2__
    inline __m128d always_inline add(__m128d x, __m128d y) noexcept
    { return _mm_add_pd(x,y); }
    inline __m128i always_inline add(__m128i x, __m128i y) noexcept
    { return _mm_add_epi32(x,y); }
# endif
# ifdef __AVX__
    inline __m256 always_inline add(__m256 x, __m256 y) noexcept
    { return _mm256_add_ps(x,y); }
    inline __m256d always_inline add(__m256d x, __m256d y) noexcept
    { return _mm256_add_pd(x,y); }
# endif
# ifdef __AVX2__
    inline __m256i always_inline add(__m256i x, __m256i y) noexcept
    { return _mm256_add_epi32(x,y); }
# endif

    ///
    /// subtraction
    ///
    inline __m128 always_inline sub(__m128 x, __m128 y) noexcept
    { return _mm_sub_ps(x,y); }
# ifdef __SSE2__
    inline __m128d always_inline sub(__m128d x, __m128d y) noexcept
    { return _mm_sub_pd(x,y); }
    inline __m128i always_inline sub(__m128i x, __m128i y) noexcept
    { return _mm_sub_epi32(x,y); }
# endif
# ifdef __AVX__
    inline __m256 always_inline sub(__m256 x, __m256 y) noexcept
    { return _mm256_sub_ps(x,y); }
    inline __m256d always_inline sub(__m256d x, __m256d y) noexcept
    { return _mm256_sub_pd(x,y); }
# endif
# ifdef __AVX2__
    inline __m256i always_inline sub(__m256i x, __m256i y) noexcept
    { return _mm256_sub_epi32(x,y); }
# endif

    ///
    /// negation
    ///
    inline __m128 always_inline neg(__m128 x) noexcept
    { return _mm_xor_ps(x,sgn_mask_ps); }
# ifdef __SSE2__
    inline __m128d always_inline neg(__m128d x) noexcept
    { return _mm_xor_pd(x,sgn_mask_pd); }
    inline __m128i always_inline neg(__m128i x) noexcept
    { return _mm_xor_si128(x,sgn_mask_pi); }
# endif
# ifdef __AVX__
    inline __m256 always_inline neg(__m256 x) noexcept
    { return _mm256_xor_ps(x,sgn_mask_qs); }
    inline __m256d always_inline neg(__m256d x) noexcept
    { return _mm256_xor_pd(x,sgn_mask_qd); }
    inline __m256i always_inline neg(__m256i x) noexcept
#  ifdef __AVX2__
    { return _mm256_xor_si256(x,sgn_mask_qi); }
#  else
    { return _mm256_castps_si256(_mm256_xor_ps(_mm256_castsi256_ps(x),
					       sgn_mask_qs)); }
#  endif
# endif

    ///
    /// multiplication
    ///
    inline __m128 always_inline mul(__m128 x, __m128 y) noexcept
    { return _mm_mul_ps(x,y); }
# ifdef __SSE2__
    inline __m128d always_inline mul(__m128d x, __m128d y) noexcept
    { return _mm_mul_pd(x,y); }
    inline __m128i always_inline mul(__m128i x, __m128i y) noexcept
    { return _mm_mul_epi32(x,y); }
# endif
# ifdef __AVX__
    inline __m256 always_inline mul(__m256 x, __m256 y) noexcept
    { return _mm256_mul_ps(x,y); }
    inline __m256d always_inline mul(__m256d x, __m256d y) noexcept
    { return _mm256_mul_pd(x,y); }
# endif
# ifdef __AVX2__
    inline __m256i always_inline mul(__m256i x, __m256i y) noexcept
    { return _mm256_mul_epi32(x,y); }
# endif

    ///
    /// division
    ///
    inline __m128 always_inline div(__m128 x, __m128 y) noexcept
    { return _mm_div_ps(x,y); }
# ifdef __SSE2__
    inline __m128d always_inline div(__m128d x, __m128d y) noexcept
    { return _mm_div_pd(x,y); }
# endif
# ifdef __AVX__
    inline __m256 always_inline div(__m256 x, __m256 y) noexcept
    { return _mm256_div_ps(x,y); }
    inline __m256d always_inline div(__m256d x, __m256d y) noexcept
    { return _mm256_div_pd(x,y); }
# endif

    ///
    /// reciprocal
    ///
    /// \note we are not using the SSE or AVX intrinsics for the reciprocal,
    ///       because its accuracy is not sufficient (see SSE documentation)
    ///
    inline __m128 always_inline rcp(__m128 x) noexcept
    { return _mm_div_ps(_mm_set1_ps(1.f),x); }
# ifdef __SSE2__
    inline __m128d always_inline rcp(__m128d x) noexcept
    { return _mm_div_pd(_mm_set1_pd(1.0),x); }
# endif
# ifdef __AVX__
    inline __m256 always_inline rcp(__m256 x) noexcept
    { return _mm256_div_ps(_mm256_set1_ps(1.f),x); }
    inline __m256d always_inline rcp(__m256d x) noexcept
    { return _mm256_div_pd(_mm256_set1_pd(1.0),x); }
# endif

    ///
    /// sqrt
    ///
    inline __m128 always_inline sqrt(__m128 x) noexcept
    { return _mm_sqrt_ps(x); }
# ifdef __SSE2__
    inline __m128d always_inline sqrt(__m128d x) noexcept
    { return _mm_sqrt_pd(x); }
# endif
# ifdef __AVX__
    inline __m256 always_inline sqrt(__m256 x) noexcept
    { return _mm256_sqrt_ps(x); }
    inline __m256d always_inline sqrt(__m256d x) noexcept
    { return _mm256_sqrt_pd(x); }
# endif

    ///
    /// maximum
    ///
    inline __m128 always_inline max(__m128 x, __m128 y) noexcept
    { return _mm_max_ps(x,y); }
# ifdef __SSE2__
    inline __m128d always_inline max(__m128d x, __m128d y) noexcept
    { return _mm_max_pd(x,y); }
# endif
# ifdef __SSE4_1__
    inline __m128i always_inline max(__m128i x, __m128i y) noexcept
    { return _mm_max_epi32(x,y); }
# endif
# ifdef __AVX__
    inline __m256 always_inline max(__m256 x, __m256 y) noexcept
    { return _mm256_max_ps(x,y); }
    inline __m256d always_inline max(__m256d x, __m256d y) noexcept
    { return _mm256_max_pd(x,y); }
# endif
# ifdef __AVX2__
    inline __m256i always_inline max(__m256i x, __m256i y) noexcept
    { return _mm256_max_epi32(x,y); }
# endif

    ///
    /// minimum
    ///
    inline __m128 always_inline min(__m128 x, __m128 y) noexcept
    { return _mm_min_ps(x,y); }
# ifdef __SSE2__
    inline __m128d always_inline min(__m128d x, __m128d y) noexcept
    { return _mm_min_pd(x,y); }
# endif
# ifdef __SSE4_1__
    inline __m128i always_inline min(__m128i x, __m128i y) noexcept
    { return _mm_min_epi32(x,y); }
# endif
# ifdef __AVX__
    inline __m256 always_inline min(__m256 x, __m256 y) noexcept
    { return _mm256_min_ps(x,y); }
    inline __m256d always_inline min(__m256d x, __m256d y) noexcept
    { return _mm256_min_pd(x,y); }
# endif
# ifdef __AVX2__
    inline __m256i always_inline min(__m256i x, __m256i y) noexcept
    { return _mm256_min_epi32(x,y); }
# endif

    ///
    /// bit-wise and
    ///
    inline __m128 always_inline bit_and(__m128 x, __m128 y) noexcept
    { return _mm_and_ps(x,y); }
# ifdef __SSE2__
    inline __m128d always_inline bit_and(__m128d x, __m128d y) noexcept
    { return _mm_and_pd(x,y); }
    inline __m128i always_inline bit_and(__m128i x, __m128i y) noexcept
    { return _mm_and_si128(x,y); }
# endif
# ifdef __AVX__
    inline __m256 always_inline bit_and(__m256 x, __m256 y) noexcept
    { return _mm256_and_ps(x,y); }
    inline __m256d always_inline bit_and(__m256d x, __m256d y) noexcept
    { return _mm256_and_pd(x,y); }
    inline __m256i always_inline bit_and(__m256i x, __m256i y) noexcept
#  ifdef __AVX2__
    { return _mm256_and_si256(x,y); }
#  else
    { return _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(x),
					       _mm256_castsi256_ps(y))); }
#  endif
# endif

    ///
    /// bit-wise or
    ///
    inline __m128 always_inline bit_or(__m128 x, __m128 y) noexcept
    { return _mm_or_ps(x,y); }
# ifdef __SSE2__
    inline __m128d always_inline bit_or(__m128d x, __m128d y) noexcept
    { return _mm_or_pd(x,y); }
    inline __m128i always_inline bit_or(__m128i x, __m128i y) noexcept
    { return _mm_or_si128(x,y); }
# endif
# ifdef __AVX__
    inline __m256 always_inline bit_or(__m256 x, __m256 y) noexcept
    { return _mm256_or_ps(x,y); }
    inline __m256d always_inline bit_or(__m256d x, __m256d y) noexcept
    { return _mm256_or_pd(x,y); }
    inline __m256i always_inline bit_or(__m256i x, __m256i y) noexcept
#  ifdef __AVX2__
    { return _mm256_or_si256(x,y); }
#  else
    { return _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(x),
					      _mm256_castsi256_ps(y))); }
#  endif
# endif

    ///
    /// bit-wise xor
    ///
    inline __m128 always_inline bit_xor(__m128 x, __m128 y) noexcept
    { return _mm_xor_ps(x,y); }
# ifdef __SSE2__
    inline __m128d always_inline bit_xor(__m128d x, __m128d y) noexcept
    { return _mm_xor_pd(x,y); }
    inline __m128i always_inline bit_xor(__m128i x, __m128i y) noexcept
    { return _mm_xor_si128(x,y); }
# endif
# ifdef __AVX__
    inline __m256 always_inline bit_xor(__m256 x, __m256 y) noexcept
    { return _mm256_xor_ps(x,y); }
    inline __m256d always_inline bit_xor(__m256d x, __m256d y) noexcept
    { return _mm256_xor_pd(x,y); }
    inline __m256i always_inline bit_xor(__m256i x, __m256i y) noexcept
#  ifdef __AVX2__
    { return _mm256_xor_si256(x,y); }
#  else
    { return _mm256_castps_si256(_mm256_xor_ps(_mm256_castsi256_ps(x),
					       _mm256_castsi256_ps(y))); }
#  endif
# endif

    ///
    /// abs(x)
    ///
    inline __m128 always_inline abs(__m128 x) noexcept
    { return _mm_and_ps(x,abs_mask_ps); }
# ifdef __SSE2__
    inline __m128d always_inline abs(__m128d x) noexcept
    { return _mm_and_pd(x,abs_mask_pd); }
    inline __m128i always_inline abs(__m128i x) noexcept
    { return _mm_and_si128(x,abs_mask_pi); }
# endif
# ifdef __AVX__
    inline __m256 always_inline abs(__m256 x) noexcept
    { return _mm256_and_ps(x,abs_mask_qs); }
    inline __m256d always_inline abs(__m256d x) noexcept
    { return _mm256_and_pd(x,abs_mask_qd); }
    inline __m256i always_inline abs(__m128i x) noexcept
#  ifdef __AVX2__
    { return _mm256_and_si256(x,abs_mask_qi); }
#  else
    { return _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(x),
					       abs_mask_qs)); }
#  endif
# endif

    ///
    /// -abs(x)
    ///
    inline __m128 always_inline negabs(__m128 x) noexcept
    { return _mm_or_ps(x,sgn_mask_ps); }
# ifdef __SSE2__
    inline __m128d always_inline negabs(__m128d x) noexcept
    { return _mm_or_pd(x,sgn_mask_pd); }
    inline __m128i always_inline negabs(__m128i x) noexcept
    { return _mm_or_si128(x,sgn_mask_pi); }
# endif
# ifdef __AVX__
    inline __m256 always_inline negabs(__m256 x) noexcept
    { return _mm256_or_ps(x,sgn_mask_qs); }
    inline __m256d always_inline negabs(__m256d x) noexcept
    { return _mm256_or_pd(x,sgn_mask_qd); }
    inline __m256i always_inline negabs(__m256i x) noexcept
#  ifdef __AVX2__
    { return _mm256_or_si256(x,sgn_mask_qi); }
#  else
    { return _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(x),
					      sgn_mask_qs)); }
#  endif
# endif

    ///
    /// abs(x-y)
    ///
    inline __m128 always_inline diff(__m128 x, __m128 y) noexcept
    { return _mm_and_ps(_mm_sub_ps(x,y),abs_mask_ps); }
# ifdef __SSE2__
    inline __m128d always_inline diff(__m128d x, __m128d y) noexcept
    { return _mm_and_pd(_mm_sub_pd(x,y),abs_mask_pd); }
    inline __m128i always_inline diff(__m128i x, __m128i y) noexcept
    { return abs(sub(x,y)); }
# endif
# ifdef __AVX__
    inline __m256 always_inline diff(__m256 x, __m256 y) noexcept
    { return _mm256_and_ps(_mm256_sub_ps(x,y),abs_mask_qs); }
    inline __m256d always_inline diff(__m256d x, __m256d y) noexcept
    { return _mm256_and_pd(_mm256_sub_pd(x,y),abs_mask_qd); }
    inline __m256i always_inline diff(__m256i x, __m256i y) noexcept
    { return abs(sub(x,y)); }
# endif

    ///
    /// just the sign bits
    ///
    inline __m128 always_inline signmask(__m128 x) noexcept
    { return _mm_and_ps(x,sgn_mask_ps); }
# ifdef __SSE2__
    inline __m128d always_inline signmask(__m128d x) noexcept
    { return _mm_and_pd(x,sgn_mask_pd); }
    inline __m128i always_inline signmask(__m128i x) noexcept
    { return _mm_and_si128(x,sgn_mask_pi); }
# endif
# ifdef __AVX__
    inline __m256 always_inline signmask(__m256 x) noexcept
    { return _mm256_and_ps(x,sgn_mask_qs); }
    inline __m256d always_inline signmask(__m256d x) noexcept
    { return _mm256_and_pd(x,sgn_mask_qd); }
    inline __m256i always_inline signmask(__m256i x) noexcept
#  ifdef __AVX2__
    { return _mm256_and_si256(x,sgn_mask_qi); }
#  else
    { return _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(x),
					       sgn_mask_qs)); }
#  endif
# endif
  
    ///
    /// sign(x)*y
    ///
    inline __m128 always_inline signmove(__m128 x, __m128 y) noexcept
    { return bit_or(signmask(x),abs(y)); }
# ifdef __SSE2__
    inline __m128d always_inline signmove(__m128d x, __m128d y) noexcept
    { return bit_or(signmask(x),abs(y)); }
    inline __m128i always_inline signmove(__m128i x, __m128i y) noexcept
    { return bit_or(signmask(x),abs(y)); }
# endif
# ifdef __AVX__
    inline __m256 always_inline signmove(__m256 x, __m256 y) noexcept
    { return bit_or(signmask(x),abs(y)); }
    inline __m256d always_inline signmove(__m256d x, __m256d y) noexcept
    { return bit_or(signmask(x),abs(y)); }
    inline __m256i always_inline signmove(__m256i x, __m256i y) noexcept
    { return bit_or(signmask(x),abs(y)); }
# endif

    ///
    /// equal
    ///
    inline __m128 always_inline cmpeq(__m128 x, __m128 y) noexcept
    { return _mm_cmpeq_ps(x,y); }
# ifdef __SSE2__
    inline __m128d always_inline cmpeq(__m128d x, __m128d y) noexcept
    { return _mm_cmpeq_pd(x,y); }
    inline __m128i always_inline cmpeq(__m128i x, __m128i y) noexcept
    { return _mm_cmpeq_epi32(x,y); }
# endif
# ifdef __AVX__
    inline __m256 always_inline cmpeq(__m256 x, __m256 y) noexcept
    { return _mm256_cmp_ps(x,y,_CMP_EQ_UQ); }
    inline __m256d always_inline cmpeq(__m256d x, __m256d y) noexcept
    { return _mm256_cmp_pd(x,y,_CMP_EQ_UQ); }
# endif
# ifdef __AVX2__
    inline __m256i always_inline cmpeq(__m256i x, __m256i y) noexcept
    { return _mm256_cmpeq_epi32(x,y); }
# endif

    ///
    /// not equal
    ///
    inline __m128 always_inline cmpneq(__m128 x, __m128 y) noexcept
    { return _mm_cmpneq_ps(x,y); }
# ifdef __SSE2__
    inline __m128d always_inline cmpneq(__m128d x, __m128d y) noexcept
    { return _mm_cmpneq_pd(x,y); }
# endif
# ifdef __AVX__
    inline __m256 always_inline cmpneq(__m256 x, __m256 y) noexcept
    { return _mm256_cmp_ps(x,y,_CMP_NEQ_UQ); }
    inline __m256d always_inline cmpneq(__m256d x, __m256d y) noexcept
    { return _mm256_cmp_pd(x,y,_CMP_NEQ_UQ); }
# endif

    ///
    /// less than
    ///
    inline __m128 always_inline cmplt(__m128 x, __m128 y) noexcept
    { return _mm_cmplt_ps(x,y); }
# ifdef __SSE2__
    inline __m128d always_inline cmplt(__m128d x, __m128d y) noexcept
    { return _mm_cmplt_pd(x,y); }
    inline __m128i always_inline cmplt(__m128i x, __m128i y) noexcept
    { return _mm_cmplt_epi32(x,y); }
# endif
# ifdef __AVX__
    inline __m256 always_inline cmplt(__m256 x, __m256 y) noexcept
    { return _mm256_cmp_ps(x,y,_CMP_LT_UQ); }
    inline __m256d always_inline cmplt(__m256d x, __m256d y) noexcept
    { return _mm256_cmp_pd(x,y,_CMP_LT_UQ); }
# endif
# ifdef __AVX2__
    inline __m256i always_inline cmpgt(__m256i x, __m256i y) noexcept
    { return _mm256_cmpgt_epi32(y,x); }
# endif

    ///
    /// less than or equal
    ///
    inline __m128 always_inline cmple(__m128 x, __m128 y) noexcept
    { return _mm_cmple_ps(x,y); }
# ifdef __SSE2__
    inline __m128d always_inline cmple(__m128d x, __m128d y) noexcept
    { return _mm_cmple_pd(x,y); }
# endif
# ifdef __AVX__
    inline __m256 always_inline cmple(__m256 x, __m256 y) noexcept
    { return _mm256_cmp_ps(x,y,_CMP_LE_UQ); }
    inline __m256d always_inline cmple(__m256d x, __m256d y) noexcept
    { return _mm256_cmp_pd(x,y,_CMP_LE_UQ); }
# endif

    ///
    /// greater than
    ///
    inline __m128 always_inline cmpgt(__m128 x, __m128 y) noexcept
    { return _mm_cmpgt_ps(x,y); }
# ifdef __SSE2__
    inline __m128d always_inline cmpgt(__m128d x, __m128d y) noexcept
    { return _mm_cmpgt_pd(x,y); }
    inline __m128i always_inline cmpgt(__m128i x, __m128i y) noexcept
    { return _mm_cmpgt_epi32(x,y); }
# endif
# ifdef __AVX__
    inline __m256 always_inline cmpgt(__m256 x, __m256 y) noexcept
    { return _mm256_cmp_ps(x,y,_CMP_GT_UQ); }
    inline __m256d always_inline cmpgt(__m256d x, __m256d y) noexcept
    { return _mm256_cmp_pd(x,y,_CMP_GT_UQ); }
# endif
# ifdef __AVX2__
    inline __m256i always_inline cmpgt(__m256i x, __m256i y) noexcept
    { return _mm256_cmpgt_epi32(x,y); }
# endif

    ///
    /// grater than or equal
    ///
    inline __m128 always_inline cmpge(__m128 x, __m128 y) noexcept
    { return _mm_cmpge_ps(x,y); }
# ifdef __SSE2__
    inline __m128d always_inline cmpge(__m128d x, __m128d y) noexcept
    { return _mm_cmpge_pd(x,y); }
# endif
# ifdef __AVX__
    inline __m256 always_inline cmpge(__m256 x, __m256 y) noexcept
    { return _mm256_cmp_ps(x,y,_CMP_GE_UQ); }
    inline __m256d always_inline cmpge(__m256d x, __m256d y) noexcept
    { return _mm256_cmp_pd(x,y,_CMP_GE_UQ); }
# endif

    ///
    /// integer whose bits are the sign bits of the operand
    ///
    inline int always_inline movemask(__m128 x) noexcept
    { return _mm_movemask_ps(x); }
# ifdef __SSE2__
    inline int always_inline movemask(__m128d x) noexcept
    { return _mm_movemask_pd(x); }
    inline int always_inline movemask(__m128i x) noexcept
    { return _mm_movemask_ps(_mm_castsi128_ps(x)); }
# endif
# ifdef __AVX__
    inline int always_inline movemask(__m256 x) noexcept
    { return _mm256_movemask_ps(x); }
    inline int always_inline movemask(__m256d x) noexcept
    { return _mm256_movemask_pd(x); }
    inline int always_inline movemask(__m256i x) noexcept
    { return _mm256_movemask_ps(_mm256_castsi256_ps(x)); }
# endif

    ///
    /// horizontal sum
    ///
    inline float always_inline sum(__m128 x) noexcept
    {
# if defined(__GNUC__) && !defined(__INTEL_COMPILER)
      __v4sf a,t,s;
      a = x;
      t = __builtin_ia32_shufps(a,a,_MM_SHUFFLE (2,3,0,1));
      s = __builtin_ia32_addps(a,t);
      t = __builtin_ia32_movhlps(s,s);
      s = __builtin_ia32_addps(s,t);
      return __builtin_ia32_vec_ext_v4sf(s,0);
# else
      __m128 s = _mm_add_ps(x,_mm_shuffle_ps(x,x,_MM_SHUFFLE (2,3,0,1)));
      return _mm_cvtss_f32(_mm_add_ps(s,_mm_movehl_ps(s,s)));
# endif
    }
# ifdef __SSE2__
    inline int always_inline sum(__m128i x) noexcept
    {
      __m128i s;
      s = _mm_add_epi32(x,_mm_shuffle_epi32(x,_MM_SHUFFLE (2,3,0,1)));
      WDutils::meta::float_and_int tmp;
//       tmp.x = _mm_cvtss_f32
// 	(reinterpret_cast<__m128>
// 	 (_mm_add_epi32(s,reinterpret_cast<__m128i>
// 			(_mm_movehl_ps(reinterpret_cast<__m128>(s),
// 				       reinterpret_cast<__m128>(s))))));
      tmp.f = _mm_cvtss_f32
	(_mm_castsi128_ps(_mm_add_epi32(s,_mm_castps_si128
					(_mm_movehl_ps(_mm_castsi128_ps(s),
						       _mm_castsi128_ps(s))))));
      return tmp.i;
    }
# endif

    ///
    /// horizontal sum in each element
    ///
    inline __m128 always_inline hadd(__m128 x) noexcept
    {
# if defined(__GNUC__) && !defined(__INTEL_COMPILER)
      __v4sf a,t,s;
      a = x;
      t = __builtin_ia32_shufps(a,a,_MM_SHUFFLE (2,3,0,1));
      s = __builtin_ia32_addps(a,t);
      t = __builtin_ia32_shufps(s,s,_MM_SHUFFLE (1,0,3,2));
      return __builtin_ia32_addps(s,t);
# else
      __m128 s = _mm_add_ps(x,_mm_shuffle_ps(x,x,_MM_SHUFFLE (2,3,0,1)));
      return _mm_add_ps(s,_mm_shuffle_ps(s,s,_MM_SHUFFLE (1,0,3,2)));
# endif
    }

    ///
    /// horizontal maximum
    ///
    inline float always_inline max(__m128 x) noexcept
    {
# if defined(__GNUC__) && !defined(__INTEL_COMPILER)
      __v4sf a,t,s;
      a = x;
      t = __builtin_ia32_shufps(a,a,_MM_SHUFFLE (2,3,0,1));
      s = __builtin_ia32_maxps(a,t);
      t = __builtin_ia32_movhlps(s,s);
      s = __builtin_ia32_maxps(s,t);
      return __builtin_ia32_vec_ext_v4sf(s,0);
# else
      __m128 s = _mm_max_ps(x,_mm_shuffle_ps(x,x,_MM_SHUFFLE (2,3,0,1)));
      return _mm_cvtss_f32(_mm_max_ps(s,_mm_movehl_ps(s,s)));
# endif
    }

    ///
    /// horizontal minimum
    ///
    inline float always_inline min(__m128 x) noexcept
    {
# if defined(__GNUC__) && !defined(__INTEL_COMPILER)
      __v4sf a,t,s;
      a = x;
      t = __builtin_ia32_shufps(a,a,_MM_SHUFFLE (2,3,0,1));
      s = __builtin_ia32_minps(a,t);
      t = __builtin_ia32_movhlps(s,s);
      s = __builtin_ia32_minps(s,t);
      return __builtin_ia32_vec_ext_v4sf(s,0);
# else
      __m128 s = _mm_min_ps(x,_mm_shuffle_ps(x,x,_MM_SHUFFLE (2,3,0,1)));
      return _mm_cvtss_f32(_mm_min_ps(s,_mm_movehl_ps(s,s)));
# endif
    }
# ifdef __SSE2__
    inline int always_inline min(__m128i x) noexcept
    // a poor man's implementation
    {
      union WDutils__align16 {
	int32_t i[4];
	float   x[4];
      } tmp;
      _mm_store_ps(tmp.x,_mm_castsi128_ps(x));
      int m=tmp.i[0];
      if(tmp.i[1]<m) m=tmp.i[1];
      if(tmp.i[2]<m) m=tmp.i[2];
      if(tmp.i[3]<m) m=tmp.i[3];
      return m;
    }
# endif
    ///
    /// return vector with elements equal to zero
    ///
    template<typename vec>
    inline vec always_inline zero() noexcept
    { return aux::VecOps<vec>::zero(); }

    ///
    /// return vector with all elements equal to one
    ///
    template<typename vec>
    inline vec always_inline one() noexcept
    { return aux::VecOps<vec>::one(); }

    ///
    /// return vector with all elements equal to same value
    ///
    template<typename vec>
    inline vec always_inline set(const typename aux::VecOps<vec>::single_type x)
      noexcept
    { return aux::VecOps<vec>::set(x); }

    ///
    /// return vector with elements equal to given values
    ///
    inline __m128 always_inline set(const float x, const float y,
				    const float z, const float w) noexcept
    { return _mm_set_ps(x,y,z,w); }
# ifdef __SSE2__
    inline __m128d always_inline set(const double x, const double y) noexcept
    { return _mm_set_pd(x,y); }
    inline __m128i always_inline set(const int x, const int y,
				     const int z, const int w) noexcept
    { return _mm_set_epi32(x,y,z,w); }
# endif
# ifdef __AVX__
    inline __m256 always_inline set(const float a, const float b,
				    const float d, const float d,
				    const float e, const float f,
				    const float g, const float h) noexcept
    { return _mm256_set_ps(a,b,c,d,e,f,g,h); }
    inline __m256d always_inline set(const double x, const double y,
				     const double z, const double w) noexcept
    { return _mm256_set_pd(x,y,z,w); }
    inline __m256i always_inline set(const int a, const int b,
				     const int d, const int d,
				     const int e, const int f,
				     const int g, const int h) noexcept
    { return _mm256_set_epi32(a,b,c,d,e,f,g,h); }
# endif

    ///
    /// return vector with elements equal to given values in reverse order
    ///
    inline __m128 always_inline setr(const float x, const float y,
				     const float z, const float w) noexcept
    { return _mm_setr_ps(x,y,z,w); }
# ifdef __SSE2__
    inline __m128d always_inline setr(const double x, const double y) noexcept
    { return _mm_setr_pd(x,y); }
    inline __m128i always_inline setr(const int x, const int y,
				      const int z, const int w) noexcept
    { return _mm_setr_epi32(x,y,z,w); }
# endif
# ifdef __AVX__
    inline __m256 always_inline setr(const float a, const float b,
				     const float d, const float d,
				     const float e, const float f,
				     const float g, const float h) noexcept
    { return _mm256_setr_ps(a,b,c,d,e,f,g,h); }
    inline __m256d always_inline setr(const double x, const double y,
				      const double z, const double w) noexcept
    { return _mm256_setr_pd(x,y,z,w); }
    inline __m256i always_inline setr(const int a, const int b,
				      const int d, const int d,
				      const int e, const int f,
				      const int g, const int h) noexcept
    { return _mm256_setr_epi32(a,b,c,d,e,f,g,h); }
# endif

    ///
    /// load vector from aligned memory
    ///
    template<typename vec>
    inline vec always_inline
    load(const typename aux::VecOps<vec>::single_type*p)
    { return aux::VecOps<vec>::load(p); }

    ///
    /// load vector from unaligned memory
    ///
    template<typename vec>
    inline vec always_inline
    loadu(const typename aux::VecOps<vec>::single_type*p)
    { return aux::VecOps<vec>::loadu(p); }

    ///
    /// \name store vector to aligned memory
    ///
    inline void always_inline store(float*p, __m128 x) noexcept
    { _mm_store_ps(p,x); }
# ifdef __SSE2__
    inline void always_inline store(double*p, __m128d x) noexcept
    { _mm_store_pd(p,x); }
    inline void always_inline store(int32_t*p, __m128i x) noexcept
    { _mm_store_ps(reinterpret_cast<float*>(p), _mm_castsi128_ps(x)); }
# endif
# ifdef __AVX__
    inline void always_inline store(float*p, __m256 x) noexcept
    { _mm256_store_ps(p,x); }
    inline void always_inline store(double*p, __m256d x) noexcept
    { _mm256_store_pd(p,x); }
    inline void always_inline store(int32_t*p, __m128i x) noexcept
    { _mm256_store_ps(reinterpret_cast<float*>(p), _mm256_castsi256_ps(x)); }
# endif

    ///
    /// \name store vector to unaligned memory
    ///
    inline void always_inline storeu(float*p, __m128 x) noexcept
    { _mm_storeu_ps(p,x); }
# ifdef __SSE2__
    inline void always_inline storeu(double*p, __m128d x) noexcept
    { _mm_storeu_pd(p,x); }
    inline void always_inline storeu(int*p, __m128i x) noexcept
    { _mm_storeu_ps(reinterpret_cast<float*>(p), _mm_castsi128_ps(x)); }
# endif
# ifdef __AVX__
    inline void always_inline storeu(float*p, __m256 x) noexcept
    { _mm256_storeu_ps(p,x); }
    inline void always_inline storeu(double*p, __m256d x) noexcept
    { _mm256_storeu_pd(p,x); }
    inline void always_inline storeu(int32_t*p, __m128i x) noexcept
    { _mm256_storeu_ps(reinterpret_cast<float*>(p), _mm256_castsi256_ps(x)); }
# endif
    ///
    /// conversions (not just casts) between vector types
    ///
# ifdef __SSE2__
    /// convert two pairs of double to four packed float
    inline __m128 always_inline convert(__m128d x, __m128d y) noexcept
    { return _mm_movelh_ps(_mm_cvtpd_ps(x),_mm_cvtpd_ps(y)); }
# endif
# ifdef __AVX__
    inline __m128 always_inline convert(__m256d x) noexcept
    { return _mm256_cvtpd_ps(x); }
    inline __m256 always_inline convert(__m256d x, __m256d x) noexcept
    { return _mm256_insertf128_ps(_mm256_castps128_ps256(_mm256_cvtpd_ps(x)),
				  _mm256_cvtpd_ps(y),1); }
    //@}
# endif
    //@}
#endif//__SSE__
  } // namespace SSE
  //
#ifdef __SSE__

# ifdef __INTEL_COMPILER
  ///
  /// \name SSE/AVX operators
  //@{

  ///
  /// x+=y
  ///
  inline __m128& always_inline operator+=(__m128&x, __m128 y) noexcept
  { return x = _mm_add_ps(x,y); }
#  ifdef __SSE2__
  inline __m128d& always_inline operator+=(__m128d&x, __m128d y) noexcept
  { return x = _mm_add_pd(x,y); }
  inline __m128i& always_inline operator+=(__m128i&x, __m128i y) noexcept
  { return x = _mm_add_epi32(x,y); }
#  endif
#  ifdef __AVX__
  inline __m256& always_inline operator+=(__m256&x, __m256 y) noexcept
  { return x = _mm256_add_ps(x,y); }
  inline __m256d& always_inline operator+=(__m256d&x, __m256d y) noexcept
  { return x = _mm256_add_pd(x,y); }
#  endif
# ifdef __AVX2__
  inline __m256i& always_inline operator+=(__m256i&x, __m256i y) noexcept
  { return x = _mm256_add_epi32(x,y); }
# endif

  ///
  /// x+y
  ///
  inline __m128 always_inline operator+(__m128 x, __m128 y) noexcept
  { return _mm_add_ps(x,y); }
#  ifdef __SSE2__
  inline __m128d always_inline operator+(__m128d x, __m128d y) noexcept
  { return _mm_add_pd(x,y); }
  inline __m128i always_inline operator+(__m128i x, __m128i y) noexcept
  { return _mm_add_epi32(x,y); }
#  endif
#  ifdef __AVX__
  inline __m256 always_inline operator+(__m256 x, __m256 y) noexcept
  { return _mm256_add_ps(x,y); }
  inline __m256d always_inline operator+(__m256d x, __m256d y) noexcept
  { return _mm256_add_pd(x,y); }
#  endif
# ifdef __AVX2__
  inline __m256i always_inline operator+(__m256i x, __m256i y) noexcept
  { return _mm256_add_epi32(x,y); }
# endif

  ///
  /// x-=y
  ///
  inline __m128& always_inline operator-=(__m128&x, __m128 y) noexcept
  { return x = _mm_sub_ps(x,y); }
#  ifdef __SSE2__
  inline __m128d& always_inline operator-=(__m128d&x, __m128d y) noexcept
  { return x = _mm_sub_pd(x,y); }
  inline __m128i& always_inline operator-=(__m128i&x, __m128i y) noexcept
  { return x = _mm_sub_epi32(x,y); }
#  endif
#  ifdef __AVX__
  inline __m256& always_inline operator-=(__m256&x, __m256 y) noexcept
  { return x = _mm256_sub_ps(x,y); }
  inline __m256d& always_inline operator-=(__m256d&x, __m256d y) noexcept
  { return x = _mm256_sub_pd(x,y); }
#  endif
#  ifdef __AVX2__
  inline __m256i& always_inline operator-(__m256i&x, __m256i y) noexcept
  { return x = _mm256_sub_epi32(x,y); }
#  endif

  ///
  /// x-y
  ///
  inline __m128 always_inline operator-(__m128 x, __m128 y) noexcept
  { return _mm_sub_ps(x,y); }
#  ifdef __SSE2__
  inline __m128d always_inline operator-(__m128d x, __m128d y) noexcept
  { return _mm_sub_pd(x,y); }
  inline __m128i always_inline operator-(__m128i x, __m128i y) noexcept
  { return _mm_sub_epi32(x,y); }
#  endif
#  ifdef __AVX__
  inline __m256 always_inline operator-(__m256 x, __m256 y) noexcept
  { return _mm256_sub_ps(x,y); }
  inline __m256d always_inline operator-(__m256d x, __m256d y) noexcept
  { return _mm256_sub_pd(x,y); }
#  endif
#  ifdef __AVX2__
  inline __m256i always_inline operator-(__m256i x, __m256i y) noexcept
  { return _mm256_sub_epi32(x,y); }
#  endif

  ///
  /// -x
  ///
  inline __m128 always_inline operator-(__m128 x) noexcept
  { return _mm_xor_ps(x,sgn_mask_ps); }
#  ifdef __SSE2__
  inline __m128d always_inline operator-(__m128d x) noexcept
  { return _mm_xor_pd(x,sgn_mask_pd); }
  inline __m128i always_inline operator-(__m128i x) noexcept
  { return _mm_xor_si128(x,sign_mask_pi); }
#  endif
#  ifdef __AVX__
  inline __m256 always_inline operator-(__m256 x) noexcept
  { return _mm256_xor_ps(x,sgn_mask_qs); }
  inline __m256d always_inline operator-(__m256d x) noexcept
  { return _mm256_xor_pd(x,sgn_mask_qd); }
  inline __m256i always_inline operator-(__m256i x) noexcept
#   ifdef __AVX2__
  { return _mm256_xor_si256(x,sgn_mask_qi); }
#   else
  { return _mm256_castps_si256(_mm256_xor_ps(_mm256_castsi256_ps(x),
					     sgn_mask_qs)); }
#   endif
#  endif

  ///
  /// x*=y
  ///
  inline __m128& always_inline operator*=(__m128&x, __m128 y) noexcept
  { return x = _mm_mul_ps(x,y); }
#  ifdef __SSE2__
  inline __m128d& always_inline operator*=(__m128d&x, __m128d y) noexcept
  { return x = _mm_mul_pd(x,y); }
  inline __m128i& always_inline operator*=(__m128i&x, __m128i y) noexcept
  { return x = _mm_mul_epi32(x,y); }
#  endif
#  ifdef __AVX__
  inline __m256& always_inline operator*=(__m256&x, __m256 y) noexcept
  { return x = _mm256_mul_ps(x,y); }
  inline __m256d& always_inline operator*=(__m256d&x, __m256d y) noexcept
  { return x = _mm256_mul_pd(x,y); }
#  endif
#  ifdef __AVX2__
  inline __m256i& always_inline operator*=(__m256i&x, __m256i y) noexcept
  { return x = _mm256_mul_epi32(x,y); }
#  endif

  ///
  /// x*y
  ///
  inline __m128 always_inline operator*(__m128 x, __m128 y) noexcept
  { return _mm_mul_ps(x,y); }
#  ifdef __SSE2__
  inline __m128d always_inline operator*(__m128d x, __m128d y) noexcept
  { return _mm_mul_pd(x,y); }
  inline __m128i always_inline operator*(__m128i x, __m128i y) noexcept
  { return _mm_mul_epi32(x,y); }
#  endif
#  ifdef __AVX__
  inline __m256 always_inline operator*(__m256 x, __m256 y) noexcept
  { return _mm256_mul_ps(x,y); }
  inline __m256d always_inline operator*(__m256d x, __m256d y) noexcept
  { return _mm256_mul_pd(x,y); }
#  endif
#  ifdef __AVX2__
  inline __m256i always_inline operator*(__m256i x, __m256i y) noexcept
  { return _mm256_mul_epi32(x,y); }
#  endif

  ///
  /// x/=y
  ///
  inline __m128& always_inline operator/=(__m128&x, __m128 y) noexcept
  { return x = _mm_div_ps(x,y); }
#  ifdef __SSE2__
  inline __m128d& always_inline operator/=(__m128d&x, __m128d y) noexcept
  { return x = _mm_div_pd(x,y); }
#  endif
#  ifdef __AVX__
  inline __m256& always_inline operator/=(__m256&x, __m256 y) noexcept
  { return x = _mm256_div_ps(x,y); }
  inline __m256d& always_inline operator/=(__m256d&x, __m256d y) noexcept
  { return x = _mm256_div_pd(x,y); }
#  endif

  ///
  /// x/y
  ///
  inline __m128 always_inline operator/(__m128 x, __m128 y) noexcept
  { return _mm_div_ps(x,y); }
#  ifdef __SSE2__
  inline __m128d always_inline operator/(__m128d x, __m128d y) noexcept
  { return _mm_div_pd(x,y); }
#  endif
#  ifdef __AVX__
  inline __m256 always_inline operator/(__m256 x, __m256 y) noexcept
  { return _mm256_div_ps(x,y); }
  inline __m256d always_inline operator/(__m256d x, __m256d y) noexcept
  { return _mm256_div_pd(x,y); }
#  endif

  ///
  /// x&=y
  ///
  inline __m128& always_inline operator&=(__m128&x, __m128 y) noexcept
  { return x = _mm_and_ps(x,y); }
#  ifdef __SSE2__
  inline __m128d& always_inline operator&=(__m128d&x, __m128d y) noexcept
  { return x = _mm_and_pd(x,y); }
  inline __m128i& always_inline operator&=(__m128i&x, __m128i y) noexcept
  { return x = _mm_and_epi32(x,y); }
#  endif
#  ifdef __AVX__
  inline __m256& always_inline operator&=(__m256&x, __m256 y) noexcept
  { return x = _mm256_and_ps(x,y); }
  inline __m256d& always_inline operator&=(__m256d&x, __m256d y) noexcept
  { return x = _mm256_and_pd(x,y); }
  inline __m256i&always_inline operator&=(__m256i&x, __m256d y) noexcept
#   ifdef __AVX2__
  { return x = _mm256_and_si256(x,y); }
#   else
  { return x = _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(x),
						 _mm256_castsi256_ps(y))); }
#   endif
#  endif

  ///
  /// x&y
  ///
  inline __m128 always_inline operator&(__m128 x, __m128 y) noexcept
  { return _mm_and_ps(x,y); }
#  ifdef __SSE2__
  inline __m128d always_inline operator&(__m128d x, __m128d y) noexcept
  { return _mm_and_pd(x,y); }
  inline __m128i always_inline operator&(__m128i x, __m128i y) noexcept
  { return _mm_and_epi32(x,y); }
#  endif
#  ifdef __AVX__
  inline __m256 always_inline operator&(__m256 x, __m256 y) noexcept
  { return _mm256_and_ps(x,y); }
  inline __m256d always_inline operator&(__m256d x, __m256d y) noexcept
  { return _mm256_and_pd(x,y); }
  inline __m256i always_inline operator&(__m256i x, __m256d y) noexcept
#   ifdef __AVX2__
  { return _mm256_and_si256(x,y); }
#   else
  { return _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(x),
					     _mm256_castsi256_ps(y))); }
#   endif
#  endif

  ///
  /// x|=y
  ///
  inline __m128& always_inline operator|=(__m128&x, __m128 y) noexcept
  { return x = _mm_or_ps(x,y); }
#  ifdef __SSE2__
  inline __m128d& always_inline operator|=(__m128d&x, __m128d y) noexcept
  { return x = _mm_or_pd(x,y); }
  inline __m128i& always_inline operator|=(__m128i&x, __m128i y) noexcept
  { return x = _mm_or_epi32(x,y); }
#  endif
#  ifdef __AVX__
  inline __m256& always_inline operator|=(__m256&x, __m256 y) noexcept
  { return x = _mm256_or_ps(x,y); }
  inline __m256d& always_inline operator|=(__m256d&x, __m256d y) noexcept
  { return x = _mm256_or_pd(x,y); }
  inline __m256i&always_inline operator|=(__m256i& x, __m256i y) noexcept
#   ifdef __AVX2__
  { return x = _mm256_or_si256(x,y); }
#   else
  { return x = _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(x),
						_mm256_castsi256_ps(y))); }
#   endif
#  endif

  ///
  /// x|y
  ///
  inline __m128 always_inline operator|(__m128 x, __m128 y) noexcept
  { return _mm_or_ps(x,y); }
#  ifdef __SSE2__
  inline __m128d always_inline operator|(__m128d x, __m128d y) noexcept
  { return _mm_or_pd(x,y); }
  inline __m128i always_inline operator|(__m128i x, __m128i y) noexcept
  { return _mm_or_epi32(x,y); }
#  endif
#  ifdef __AVX__
  inline __m256 always_inline operator|(__m256 x, __m256 y) noexcept
  { return _mm256_or_ps(x,y); }
  inline __m256d always_inline operator|(__m256d x, __m256d y) noexcept
  { return _mm256_or_pd(x,y); }
  inline __m256i always_inline operator|(__m256i x, __m256i y) noexcept
#   ifdef __AVX2__
  { return _mm256_or_si256(x,y); }
#   else
  { return _mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(x),
					    _mm256_castsi256_ps(y))); }
#   endif
#  endif

  ///
  /// x^=y
  ///
  inline __m128& always_inline operator^=(__m128&x, __m128 y) noexcept
  { return x = _mm_xor_ps(x,y); }
#  ifdef __SSE2__
  inline __m128d& always_inline operator^=(__m128d&x, __m128d y) noexcept
  { return x = _mm_xor_pd(x,y); }
  inline __m128i& always_inline operator^=(__m128i&x, __m128i y) noexcept
  { return x = _mm_xor_epi32(x,y); }
#  endif
#  ifdef __AVX__
  inline __m256& always_inline operator^=(__m256&x, __m256 y) noexcept
  { return x = _mm256_xor_ps(x,y); }
  inline __m256d& always_inline operator^=(__m256d&x, __m256d y) noexcept
  { return x = _mm256_xor_pd(x,y); }
  inline __m256i& always_inline operator^=(__m256i& x, __m256i y) noexcept
#   ifdef __AVX2__
  { return x = _mm256_xor_si256(x,y); }
#   else
  { return x = _mm256_castps_si256(_mm256_xor_ps(_mm256_castsi256_ps(x),
						 _mm256_castsi256_ps(y))); }
#   endif
#  endif

  ///
  /// x^y
  ///
  inline __m128 always_inline operator^(__m128 x, __m128 y) noexcept
  { return _mm_xor_ps(x,y); }
#  ifdef __SSE2__
  inline __m128d always_inline operator^(__m128d x, __m128d y) noexcept
  { return _mm_xor_pd(x,y); }
  inline __m128i always_inline operator^(__m128i x, __m128i y) noexcept
  { return _mm_xor_epi32(x,y); }
#  endif
#  ifdef __AVX__
  inline __m256 always_inline operator^(__m256 x, __m256 y) noexcept
  { return _mm256_xor_ps(x,y); }
  inline __m256d always_inline operator^(__m256d x, __m256d y) noexcept
  { return _mm256_xor_pd(x,y); }
  inline __m256i always_inline operator^(__m256i x, __m256i y) noexcept
#   ifdef __AVX2__
  { return _mm256_xor_si256(x,y); }
#   else
  { return _mm256_castps_si256(_mm256_xor_ps(_mm256_castsi256_ps(x),
					     _mm256_castsi256_ps(y))); }
#   endif
#  endif

  ///
  /// x < y
  ///
  inline __m128 always_inline operator<(__m128 x, __m128 y) noexcept
  { return _mm_cmplt_ps(x,y); }
#  ifdef __SSE2__
  inline __m128d always_inline operator<(__m128d x, __m128d y) noexcept
  { return _mm_cmplt_pd(x,y); }
  inline __m128i always_inline operator<(__m128i x, __m128i y) noexcept
  { return _mm_cmplt_epi32(x,y); }
#  endif
#  ifdef __AVX__
  inline __m256 always_inline operator<(__m256 x, __m256 y) noexcept
  { return _mm256_cmp_ps(x,y,_CMP_LT_UQ); }
  inline __m256d always_inline operator<(__m256d x, __m256d y) noexcept
  { return _mm256_cmp_pd(x,y,_CMP_LT_UQ); }
#  endif
#  ifdef __AVX2__
  inline __m256i always_inline operator<(__m256i x, __m256i y) noexcept
  { return _mm256_cmpgt_epi32(y,x); }
#  endif

  ///
  /// x <= y
  ///
  inline __m128 always_inline operator<=(__m128 x, __m128 y) noexcept
  { return _mm_cmple_ps(x,y); }
#  ifdef __SSE2__
  inline __m128d always_inline operator<=(__m128d x, __m128d y) noexcept
  { return _mm_cmple_pd(x,y); }
#  endif
#  ifdef __AVX__
  inline __m256 always_inline operator<=(__m256 x, __m256 y) noexcept
  { return _mm256_cmp_pd(x,y,_CMP_LE_UQ); }
  inline __m256d always_inline operator<=(__m256d x, __m256d y) noexcept
  { return _mm256_cmp_pd(x,y,_CMP_LE_UQ); }
#  endif

  ///
  /// x > y
  ///
  inline __m128 always_inline operator>(__m128 x, __m128 y) noexcept
  { return _mm_cmpgt_ps(x,y); }
#  ifdef __SSE2__
  inline __m128d always_inline operator>(__m128d x, __m128d y) noexcept
  { return _mm_cmpgt_pd(x,y); }
  inline __m128i always_inline operator>(__m128i x, __m128i y) noexcept
  { return _mm_cmpgt_epi32(x,y); }
#  endif
#  ifdef __AVX__
  inline __m256 always_inline operator>(__m256 x, __m256 y) noexcept
  { return _mm256_cmp_ps(x,y,_CMP_GT_UQ); }
  inline __m256d always_inline operator>(__m256d x, __m256d y) noexcept
  { return _mm256_cmp_pd(x,y,_CMP_GT_UQ); }
#  endif
#  ifdef __AVX2__
  inline __m256i always_inline operator>(__m256i x, __m256i y) noexcept
  { return _mm256_cmpgt_epi32(x,y); }
#  endif

  ///
  /// x >= y
  ///
  inline __m128 always_inline operator>=(__m128 x, __m128 y) noexcept
  { return _mm_cmpge_ps(x,y); }
#  ifdef __SSE2__
  inline __m128d always_inline operator>=(__m128d x, __m128d y) noexcept
  { return _mm_cmpge_pd(x,y); }
#  endif
#  ifdef __AVX__
  inline __m256 always_inline operator>=(__m256 x, __m256 y) noexcept
  { return _mm256_cmp_ps(x,y,_CMP_GE_UQ); }
  inline __m256d always_inline operator>=(__m256d x, __m256d y) noexcept
  { return _mm256_cmp_pd(x,y,_CMP_GE_UQ); }
#  endif
  
  //@}
# endif // __INTEL_COMPILER
# undef abs_mask_pi
# undef sgn_mask_pi
# undef abs_mask_ps
# undef sgn_mask_ps
# undef abs_mask_pd
# undef sgn_mask_pd
# undef abs_mask_qi
# undef sgn_mask_qi
# undef abs_mask_qs
# undef sgn_mask_qs
# undef abs_mask_qd
# undef sgn_mask_qd

  //
# ifdef __INTEL_COMPILER
  inline float xmm0(__m128 _A)
  {
    union { float f; int i; } tmp;
    tmp.i = _mm_extract_ps(_A,0);
    return tmp.f;
  }
  inline float xmm1(__m128 _A)
  {
    union { float f; int i; } tmp;
    tmp.i = _mm_extract_ps(_A,1);
    return tmp.f;
  }
  inline float xmm2(__m128 _A)
  {
    union { float f; int i; } tmp;
    tmp.i = _mm_extract_ps(_A,2);
    return tmp.f;
  }
  inline float xmm3(__m128 _A)
  {
    union { float f; int i; } tmp;
    tmp.i = _mm_extract_ps(_A,3);
    return tmp.f;
  }
#  ifdef __SSE2__
  inline double xmm0(__m128d _A)
  {
    union { double x; long long int i; } tmp;
    tmp.i = _mm_extract_ps(_A,0);
    return tmp.x;
  }
  inline double xmm1(__m128d _A)
  {
    union { double x; long long int i; } tmp;
    tmp.i = _mm_extract_ps(_A,1);
    return tmp.x;
  }
#  endif// __SSE2__
# elif defined(__GNUC__)  // __INTEL_COMPILER / __GNUC__
  inline float xmm0(__m128 _A)
  { return __builtin_ia32_vec_ext_v4sf(_A,0); }
  inline float xmm1(__m128 _A)
  { return __builtin_ia32_vec_ext_v4sf(_A,1); }
  inline float xmm2(__m128 _A)
  { return __builtin_ia32_vec_ext_v4sf(_A,2); }
  inline float xmm3(__m128 _A)
  { return __builtin_ia32_vec_ext_v4sf(_A,3); }
#  ifdef __SSE2__
  inline double xmm0(__m128d _A)
  { return __builtin_ia32_vec_ext_v2df(_A,0); }
  inline double xmm1(__m128d _A)
  { return __builtin_ia32_vec_ext_v2df(_A,1); }
#  endif// __SSE2__
# endif // __INTEL_COMPILER / __GNUC__
}
#endif // __SSE__
//
#undef noexcept
#undef constexpr
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
  namespace SSE {

    template<int N> struct _Aux;
    template<> struct _Aux<2> {
      template<typename _I> static _I top(_I i) { return (i+1) & ~1; }
      template<typename _I> static _I bot(_I i) { return i     & ~1; }
    };
    template<> struct _Aux<4> {
      template<typename _I> static _I top(_I i) { return (i+3) & ~3; }
      template<typename _I> static _I bot(_I i) { return i     & ~3; }
    };
    template<> struct _Aux<8> {
      template<typename _I> static _I top(_I i) { return (i+7) & ~7; }
      template<typename _I> static _I bot(_I i) { return i     & ~7; }
    };
    template<> struct _Aux<16> {
      template<typename _I> static _I top(_I i) { return (i+15)& ~15; }
      template<typename _I> static _I bot(_I i) { return i     & ~15; }
    };

    /// smallest multiple of N not less than i
    template<int N, typename _I> inline
    _I top(_I i) { return _Aux<N>::top(i); }
    /// largest multiple of N not greater than i
    template<int N, typename _I> inline
    _I bottom(_I i) { return _Aux<N>::bot(i); }

    ///
    /// simple array manipulations for unaligned arrays
    ///
    /// \note routines are up to 16/sizeof(T) times faster than simple code
    ///
    struct UnAligned {
#if __cplusplus >= 201103L
      UnAligned() = delete;
#endif
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
#if __cplusplus >= 201103L
      Aligned() = delete;
#endif
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
    template<typename _F> struct Traits;

    /// SSE::Traits<float>
    template<> struct Traits<float>
    {
#if __cplusplus >= 201103L
      Traits() = delete;
#endif
      /// is SSE enabled for this type?
#ifdef __SSE__
      static const bool sse = true;
#else
      static const bool sse = false;
#endif
      /// alignment number: K floats align to 128 bytes
      static const int K=4;
      /// smallest multiple of K not less than n
      template<typename _I>
      static _I Top(_I n) { return top<K,_I>(n); }
      /// largest multiple of K not greater than n
      template<typename _I>
      static _I Bottom(_I n) { return bottom<K,_I>(n); }
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
#if __cplusplus >= 201103L
      Traits() = delete;
#endif
      /// is SSE enabled for this type?
#ifdef __SSE2__
      static const bool sse = true;
#else
      static const bool sse = false;
#endif
      /// alignment number: K floats align to 128 bytes
      static const int K=4;
      /// smallest multiple of K not less than n
      template<typename _I>
      static _I Top(_I n) { return top<K,_I>(n); }
      /// largest multiple of K not greater than n
      template<typename _I>
      static _I Bottom(_I n) { return bottom<K,_I>(n); }
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
#if __cplusplus >= 201103L
      Traits() = delete;
#endif
      /// is SSE enabled for this type?
#ifdef __SSE2__
      static const bool sse = true;
#else
      static const bool sse = false;
#endif
      /// alignment number: K doubles align to 128 bytes
      static const int K=2;
      /// smallest multiple of K not less than n
      template<typename _I>
      static _I Top(_I n) { return top<K,_I>(n); }
      /// largest multiple of K not greater than n
      template<typename _I>
      static _I Bottom(_I n) { return bottom<K,_I>(n); }
      /// is an array of floats aligned to at least 8 bytes (=sizeof(double))?
      static bool is_minimum_aligned(double*f) {
	return (size_t(f) & 8)  == 0;
      }
      /// is an array aligned to 16 bytes?
      static bool is_aligned(double*f) {
	return (size_t(f) & 16) == 0;
      }
    };
    /// smallest multiple of 16/sizeof(_F) not less than i
    template<typename _F, typename _I> inline
    _I Top(_I i) {
      return Traits<_F>::Top(i);
    }
    /// largest multiple of 16/sizeof(_F) not greater than i
    template<typename _F, typename _I> inline
    _I Bottom(_I i) {
      return Traits<_F>::Bottom(i);
    }
    ////////////////////////////////////////////////////////////////////////////
    /// An array of float or double supporting operations via SSE instructions
    /// \note not to be confused with WDutils::Array16 in memory.h
    template<typename _F>
    class Array16 {
      WDutilsStaticAssert((is_same<_F,int>::value ||
			   is_same<_F,float>::value ||
			   is_same<_F,double>::value ));
      //  no copy ctor
      Array16(const Array16&) WDutilsCXX11Delete;
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
    /// filler object of size _S
    template<size_t _S> class Filler { char _F[_S]; };
    /// extend a given class to have size a multibple of 16-byte
    /// \note only default constructor possible for @a Base
    template<typename Base>
    class Extend16 : public Base
    {
      static const size_t Bytes = sizeof(Base);
      static const size_t Splus = Bytes & 15;
      static const size_t Added = Splus? 16-Splus : 0;
      Filler<Added> _Filler;
    };
    //
#ifdef __SSE__
    /// working horse actually implemented in sse.cc
    /// \note Not intended for consumption, use routine below instead.
    void _swap16(float*a, float*b, size_t n);
    /// swap a multiple of 16 bytes at 16-byte aligned memory.
    /// \param[in] a   16-byte aligned memory address
    /// \param[in] b   16-byte aligned memory address
    /// \param[in] n   multiple of 16: # bytes to swap at @a a and @a b
    /// \note For reasons of efficiency, the above conditions are not tested.
    ///       If either address is not 16-byte aligned, a run-time error
    ///       (segmentation fault) will result. If @a n is not a multiple of
    ///       16, the last @a n%16 bytes will not be swapped.
    inline void Swap16(void*a, void*b, size_t n)
    { _swap16(static_cast<float*>(a),static_cast<float*>(b),n>>2); }
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
      _swap16(reinterpret_cast<float*>(a), reinterpret_cast<float*>(b),
	      sizeof(Object)>>2);
    }
#endif
  } // namespace SSE
} // namespace WDutils
//
#endif
