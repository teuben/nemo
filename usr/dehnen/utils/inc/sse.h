// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/sse.h
///
/// \brief  support for SSE coding and simple SSE supported code
///
/// \author Walter Dehnen
///
/// \date   2009-2013
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009-2013 Walter Dehnen
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
// with GCC 4.4 use x86intrin.h
# if (__GNUC__ >= 4 && __GNUC_MINOR__ >= 4) || defined(__clang__)
extern "C" {
#  include <x86intrin.h>
}
# else
extern "C" {
#  ifdef __SSE__
#  include <xmmintrin.h>
#  endif
#  ifdef __SSE2__
#  include <emmintrin.h>
#  endif
#  if defined (__SSE4_2__) || defined (__SSE4_1__)
#  include <smmintrin.h>
#  endif
#  ifdef __AVX__
#  include <avxintrin.h>
#  endif
#  ifdef __AVX2__
#  include <avx2intrin.h>
#  endif
}
# endif
#elif defined(__SSE__)
// use individual headers for intrinsics

# ifdef defined(__INTEL_COMPILER) && defined(_MM_MALLOC_H_INCLUDED)
#  warning The intel compiler has seen _mm_malloc.h by GNU which declares _mm_malloc() and _mm_free() to have different linking than those declared in xmmintrin.h by INTEL, which we are going to include now. This may cause a compiler error, which can be prevented by ensuring that _mm_malloc.h is not explicitly included when using the intel compiler.
# endif

extern "C" {
# include <immintrin.h>
}
#endif    // __SSE__

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

#if __cplusplus < 201103L
# define noexcept
# define constexpr
# define static_assert(EXPR,MESG) WDutilsStaticAssert(EXPR)
#endif

#if   defined(__INTEL_COMPILER)
#  define alignas(K) __declspec(align(K))
#elif defined(__GNUC__)
#  define alignas(K) __attribute__ ((aligned(K)))
#elif __cplusplus < 201103L
#  error do not know how to enforce alignment with this compiler
#endif

//
namespace WDutils {

  /// is sizeof(T) a multiple of alignment or vice versa?
  template<typename T>
  constexpr bool aligns_at(size_t alignment) noexcept
  {
    return sizeof(T) >= alignment?
      (sizeof(T) % alignment) == 0 :
      (alignment % sizeof(T)) == 0 ;
  }

#ifdef __SSE__
# if defined(__clang__) || (defined(__GNUC__) && !defined(__INTEL_COMPILER))
#  define always_inline __attribute__((__always_inline__))
# else
#  define always_inline
# endif
  //
  namespace meta {
    union float_and_uint {
      float f; uint32_t i; 
      float_and_uint() WDutilsCXX11DefaultBody
      float_and_uint(uint32_t k) : i(k) {}
    };
    union double_and_uint {
      double d; uint64_t i; 
      double_and_uint() WDutilsCXX11DefaultBody
      double_and_uint(uint64_t k) : i(k) {}
    };
  }
  ///
  /// generic support for coding with SSE/AVX intrinsics
  ///
  namespace SSE {

    ///
    /// @a K packed floating-point number of type @c T
    ///
    /// The idea is to provide a unique and simple C++ interface to use the
    /// SSE and AVX instruction set with any compiler that supports the usual
    /// intrinsics (like _mm_mul_ps). The various member functions and
    /// operators are directly implemented in terms of these intrinsics, making
    /// for an efficient yet convenient SSE/AVX interface.
    ///

    template<int K, typename T> struct packed;

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdocumentation"
    ///
    /// @a packed_is_supported<a,b>::value=true if @a packed<a,b> is supported
    ///
#pragma clang diagnostic pop
    template<int _S, typename _T>
    struct packed_is_supported
    {
      static const bool value = 0
	|| (_S==4 && is_same<_T,float>::value)
#ifdef __SSE2__
	|| (_S==2 && is_same<_T,double>::value)
#endif
#ifdef __AVX__
	|| (_S==8 && is_same<_T,float>::value)
	|| (_S==4 && is_same<_T,double>::value)
#endif
	;
    };
    //
#if __cplusplus >= 201103L
    constexpr unsigned log2(unsigned n) noexcept
    {
      return n<2? 0 : 1+log2(n>>1);
    }
#endif

    ///
    /// base class for packed<K,T>
    ///
    template<int K> struct packed_base
    {
      static_assert(K>0 && 0==(K&(K-1)),"block_size must be power of 2 > 0");
      /// number of elements
      static const unsigned block_size = K;
      /// block_size = 1<<block_sft
#if __cplusplus >= 201103L && !defined(__INTEL_COMPILER)
      static const unsigned block_shft = log2(block_size);
#else
      static const unsigned block_shft = meta::Integer<block_size>::Log2;
#endif
      /// mask for obtaining sub-index within block
      static const unsigned block_trim = K-1;
      /// mask for obtaining aligned index
      static const unsigned block_mask =~block_trim;
      /// maximum signbits
      static const int max_signbits = (1<<block_size)-1;
      /// is given index aligned to block_index?
      constexpr static bool is_aligned(unsigned i) noexcept
      { return (i&block_trim)==0; }
      /// aligned index given an index
      constexpr static unsigned aligned_index(unsigned i) noexcept
      { return i & block_mask; }
      /// sub-index given an index
      constexpr static unsigned sub_index(unsigned i) noexcept
      { return i & block_trim; }
      /// block index given an index
      constexpr static unsigned block_index(unsigned i) noexcept
      { return i >> block_shft; }
      /// # blocks given # elements
      constexpr static unsigned num_blocks(unsigned n) noexcept
      { return (n+block_trim)>>block_shft; }
      /// # elements in full blocks, given # elements
      constexpr static unsigned blocked_num(unsigned n) noexcept
      { return (n+block_trim)&block_mask; }
    };
    //--------------------------------------------------------------------------
    ///
    /// 4 packed single-precision floating-point numbers
    ///
    //--------------------------------------------------------------------------
    template<> struct packed<4,float> : packed_base<4>
    {
      /// \name types, constants, and static methods
      //@
      using packed_base<4>::block_size;
      using packed_base<4>::block_shft;
      using packed_base<4>::block_trim;
      using packed_base<4>::block_mask;
      using packed_base<4>::max_signbits;
      using packed_base<4>::is_aligned;
      using packed_base<4>::aligned_index;
      using packed_base<4>::sub_index;
      using packed_base<4>::block_index;
      using packed_base<4>::num_blocks;
      using packed_base<4>::blocked_num;
      /// associated SSE/AVX vector type
      typedef __m128 data_type;
      /// associated element type
      typedef float element_type;
      /// equivalent array of elements
      typedef element_type element_block[block_size];
      /// required alignement (bytes)
      static const unsigned alignment = block_size*sizeof(element_type);
      /// aligned equivalent array of elements
      typedef element_block alignas(16) aligned_element_block;
      /// a packed with all elements equal to 0
      static packed always_inline zero() noexcept
      { return packed(_mm_setzero_ps()); }
      /// a packed with all elements equal to 1
      static packed always_inline one() noexcept
      { return packed(_mm_set1_ps(1.f)); }
      /// is given pointer appropriately aligned?
      static bool is_aligned(void*p) noexcept
      { return (size_t(p)&(alignment-1))==0; }
      /// offset (number of element_types) of pointer from alignment
      static size_t offset(element_type*p) noexcept
      { return (size_t(p)>>sizeof(element_type))&block_mask; }
      //@}

      /// \name construction and assignment
      //@{
      /// default ctor
      always_inline packed() WDutilsCXX11DefaultBody
#if __cplusplus >= 201103L
      /// copy ctor
      always_inline packed(packed const&) = default;
      /// copy operator
      packed& always_inline operator=(packed const&) = default;
#endif
      /// ctor from data_type
      explicit always_inline packed(data_type m) : _m(m) {}
      /// ctor from single integer value: set all element equal to single value
      explicit always_inline packed(int x) noexcept
      { _m = _mm_set1_ps(x); }
      /// ctor from single value: set all element equal to single value
      explicit always_inline packed(float x) noexcept
      { _m = _mm_set1_ps(x); }
      /// ctor from single value: set all element equal to single value
      explicit always_inline packed(double x) noexcept
      { _m = _mm_set1_ps(float(x)); }
      /// ctor from 4 values: set elements
      /// \note inverse order to _mm_set_ps(a,b,c,d)
      always_inline packed(element_type a, element_type b,
			   element_type c, element_type d) noexcept
      { _m = _mm_setr_ps(a,b,c,d); }
      /// set all elements equal to zero
      packed& always_inline set_zero() noexcept
      { _m = _mm_setzero_ps(); return*this; }
      /// set element equal to single value
      packed& always_inline set(element_type x) noexcept
      { _m = _mm_set1_ps(x); return*this; }
      /// set elements: [a,b,c,d]
      /// \note inverse order to _mm_set_ps(a,b,c,d)
      packed& always_inline set(element_type a, element_type b,
				element_type c, element_type d) noexcept
      { _m = _mm_setr_ps(a,b,c,d); return*this; }
      //@}

      /// \name data access and conversion
      //@{
      /// conversion to const vector type
      always_inline operator data_type const&() const noexcept
      { return _m; }
      /// direct data const access
      data_type const& always_inline data() const noexcept
      { return _m; }
      /// conversion to vector type
      always_inline operator data_type&() noexcept
      { return _m; }
      /// direct non-const data access
      data_type& always_inline data() noexcept
      { return _m; }
      /// constant element access, templated
      template<unsigned I>
      friend element_type at(packed) noexcept;
# ifdef __SSE2__
      /// upcast: convert lower two elements to  @c packed<2,double>
      friend packed<2,double> upcast_lo(packed) noexcept;
      /// upcast: convert upper two elements to  @c packed<2,double>
      friend packed<2,double> upcast_hi(packed) noexcept;
      /// downcast: convert 2 @c packed<2,double>  to  @c packed<4,float>
      friend packed downcast(packed<2,double>, packed<2,double>) noexcept;
# endif
# ifdef __AVX__
      /// downcast: convert  @c packed<4,double>  to  @c packed<4,float>
      friend packed downcast(packed<4,double>) noexcept;
      /// upcast:   convert  @c packed<4,float>   to  @c packed<4,double>  
      friend packed<4,double> upcast(packed) noexcept;
# endif
      //@}

      /// \name load and store
      //@{
      /// load from aligned memory location
      static packed always_inline load(const element_type*p) noexcept
      { return packed(_mm_load_ps(p)); }
      /// load from unaligned memory location
      static packed always_inline loadu(const element_type*p) noexcept
      { return packed(_mm_loadu_ps(p)); }
      /// load from aligned memory location, using template arg for alignment
      template<bool aligned> static typename enable_if< aligned, packed>::type
      always_inline load_t(const element_type*p) noexcept
      { return packed(_mm_load_ps(p)); }
      /// load from unaligned memory location, using template arg for alignment
      template<bool aligned> static typename enable_if<!aligned,packed>::type
      always_inline load_t(const element_type*p) noexcept
      { return packed(_mm_loadu_ps(p)); }
      /// store to aligned memory location
      void always_inline store(element_type*p) const noexcept
      { _mm_store_ps(p,_m); }
      /// store to unaligned memory location
      void always_inline storeu(element_type*p) const noexcept
      { _mm_storeu_ps(p,_m); }
      /// store to aligned memory location, using template arg for alignment
      template<bool aligned> typename enable_if< aligned>::type
      always_inline store_t(element_type*p) const noexcept
      { _mm_store_ps(p,_m); }
      /// store to unaligned memory location, using template arg for alignment
      template<bool aligned> typename enable_if<!aligned>::type
      always_inline store_t(element_type*p) const noexcept
      { _mm_storeu_ps(p,_m); }
      /// load aligned object with member data() returning const element_type*
      template<typename class_with_member_data>
      static packed always_inline pack(class_with_member_data const&a) noexcept
      { return load(a.data()); }
      /// load unaligned object with member data() returning const element_type*
      template<typename class_with_member_data>
      static packed always_inline packu(class_with_member_data const&a) noexcept
      { return loadu(a.data()); }
      /// load object with member data() returning const element_type*
      template<bool aligned, typename class_with_member_data>
      static packed always_inline pack_t(class_with_member_data const&a)
	noexcept
      { return load_t<aligned>(a.data()); }
      /// store to aligned object with member data() returning element_type*
      template<typename class_with_member_data>
      void always_inline unpack(class_with_member_data&a) const noexcept
      { store(a.data()); }
      /// store to aligned object with member data() returning element_type*
      template<typename class_with_member_data>
      void always_inline unpacku(class_with_member_data&a) const noexcept
      { storeu(a.data()); }
      /// store to aligned object with member data() returning element_type*
      template<bool aligned, typename class_with_member_data>
      void always_inline unpack_t(class_with_member_data&a) const noexcept
      { store_t<aligned>(a.data()); }
      /// store directly to aligned memory without polluting cashes
      void always_inline stream(element_type*p) const noexcept
      { _mm_stream_ps(p,_m); }
      //@}

      /// \name vectorised arithmetic operations
      //@{
      /// +=
      packed& always_inline operator+=(packed p) noexcept
      { _m = _mm_add_ps(p._m,_m); return*this; }
      /// +
      packed always_inline operator+ (packed p) const noexcept
      { return packed(_mm_add_ps(_m,p._m)); }
      /// -=
      packed& always_inline operator-=(packed p) noexcept
      { _m = _mm_sub_ps(_m,p._m); return*this; }
      /// -
      packed always_inline operator- (packed p) const noexcept
      { return packed(_mm_sub_ps(_m,p._m)); }
      /// unary -
      packed always_inline operator-() const noexcept
      { return packed(_mm_xor_ps(_m,sgn_mask())); }
      /// *=
      packed& always_inline operator*=(packed p) noexcept
      { _m = _mm_mul_ps(p._m,_m); return*this; }
      /// *
      packed always_inline operator* (packed p) const noexcept
      { return packed(_mm_mul_ps(_m,p._m)); }
      /// /=
      packed& always_inline operator/=(packed p) noexcept
      { _m = _mm_div_ps(_m,p._m); return*this; }
      /// /
      packed always_inline operator/ (packed p) const noexcept
      { return packed(_mm_div_ps(_m,p._m)); }
      /// sqrt
      friend packed always_inline sqrt(packed p) noexcept
      { return packed(_mm_sqrt_ps(p._m)); }
      /// square
      friend packed always_inline square(packed p) noexcept
      { return packed(_mm_mul_ps(p._m,p._m)); }
      /// reciprocal
      friend packed always_inline reciprocal(packed p) noexcept
      { return packed(_mm_div_ps(_mm_set1_ps(1.f),p._m)); }
      /// approximate reciprocal
      friend packed always_inline approximate_reciprocal(packed p) noexcept
      { return packed(_mm_rcp_ps(p._m)); }
      /// maximum of two packed
      friend packed always_inline max(packed a, packed b) noexcept
      { return packed(_mm_max_ps(a._m,b._m)); }
      /// minimum of two packed
      friend packed always_inline min(packed a, packed b) noexcept
      { return packed(_mm_min_ps(a._m,b._m)); }
      /// abs
      friend packed always_inline abs(packed p) noexcept
      { return packed(_mm_and_ps(p._m,abs_mask())); }
      /// -abs
      friend packed always_inline negabs(packed p) noexcept
      { return packed(_mm_or_ps(p._m,sgn_mask())); }
      /// abs(x-y)
      friend packed always_inline diff(packed a, packed b) noexcept
      { return packed(_mm_and_ps(_mm_sub_ps(a._m,b._m),abs_mask())); }
      /// x = abs(x-y)
      packed& always_inline make_diff(packed p) noexcept
      {
	_m = _mm_sub_ps(_m,p._m);
	_m = _mm_and_ps(_m,abs_mask());
	return*this;
      }
      /// sign(a)*b
      friend packed always_inline signmove(packed a, packed b) noexcept
      { return packed(_mm_or_ps(a.signmask(),_mm_and_ps(b._m,abs_mask()))); }
      //@}

      /// \name horizontal arithmetic operations
      //@{
      /// horizontal sum
      friend element_type always_inline sum(packed p) noexcept
      { 
	data_type s = _mm_add_ps(p._m,_mm_shuffle_ps(p._m,p._m,
						     _MM_SHUFFLE(2,3,0,1)));
	return _mm_cvtss_f32(_mm_add_ps(s,_mm_movehl_ps(s,s)));
      }
      /// sum of first two elements
      friend element_type always_inline sum2(packed p) noexcept
      {
	aligned_element_block q;
	_mm_store_ps(q,p._m);
	return q[0]+q[1];
      }
      /// sum of first three elements
      friend element_type always_inline sum3(packed p) noexcept
      {
	aligned_element_block q;
	_mm_store_ps(q,p._m);
	return q[0]+q[1]+q[2];
      }
      /// horizontal sum in each element
      friend packed always_inline hadd(packed p) noexcept
      {
	data_type s = _mm_add_ps(p._m,_mm_shuffle_ps(p._m,p._m,
						     _MM_SHUFFLE(2,3,0,1)));
	return packed(_mm_add_ps(s,_mm_shuffle_ps(s,s,_MM_SHUFFLE(1,0,3,2))));
      }
      /// horizontal maximum
      friend element_type always_inline max(packed p) noexcept
      {
	data_type s = _mm_max_ps(p._m,_mm_shuffle_ps(p._m,p._m,
						     _MM_SHUFFLE(2,3,0,1)));
	return _mm_cvtss_f32(_mm_max_ps(s,_mm_movehl_ps(s,s)));
      }
      /// horizontal minimum
      friend element_type always_inline min(packed p) noexcept
      {
	data_type s = _mm_min_ps(p._m,_mm_shuffle_ps(p._m,p._m,
						     _MM_SHUFFLE(2,3,0,1)));
	return _mm_cvtss_f32(_mm_min_ps(s,_mm_movehl_ps(s,s)));
      }
      //@}

      /// \name vectorised comparisons
      //@{
      /// <
      packed always_inline operator< (packed p) const noexcept
      { return packed(_mm_cmplt_ps(_m,p._m)); }
      /// <=
      packed always_inline operator<=(packed p) const noexcept
      { return packed(_mm_cmple_ps(_m,p._m)); }
      /// >
      packed always_inline operator> (packed p) const noexcept
      { return packed(_mm_cmpgt_ps(_m,p._m)); }
      /// >=
      packed always_inline operator>=(packed p) const noexcept
      { return packed(_mm_cmpge_ps(_m,p._m)); }
      /// ==
      packed always_inline operator==(packed p) const noexcept
      { return packed(_mm_cmpeq_ps(_m,p._m)); }
      /// !=
      packed always_inline operator!=(packed p) const noexcept
      { return packed(_mm_cmpneq_ps(_m,p._m)); }
      //@}

      /// \name vectorised logical operations
      /// &=
      packed& always_inline operator&=(packed p) noexcept
      { _m = _mm_and_ps(p._m,_m); return*this; }
      /// &
      packed always_inline operator& (packed p) const noexcept
      { return packed(_mm_and_ps(_m,p._m)); }
      /// |=
      packed& always_inline operator|=(packed p) noexcept
      { _m = _mm_or_ps(p._m,_m); return*this; }
      /// |
      packed always_inline operator| (packed p) const noexcept
      { return packed(_mm_or_ps(_m,p._m)); }
      /// ^=
      packed& always_inline operator^=(packed p) noexcept
      { _m = _mm_xor_ps(p._m,_m); return*this; }
      /// ^
      packed always_inline operator^ (packed p) const noexcept
      { return packed(_mm_xor_ps(_m,p._m)); }
      /// unary !
      packed always_inline operator! () const noexcept
      { return packed(_mm_xor_ps(_m,one_mask())); }
      /// and not: return this&!that
      packed always_inline andnot(packed p) const noexcept
      { return packed(_mm_andnot_ps(p._m,_m)); }
      //@}

      /// \name boolean or integer properties
      //@{
      /// return integer with bits equal to sign bits
      /// \note if @a p is the result of a boolean operation (comparison, or
      ///       bit-wise operations on comparison results), signbit(p) can be
      ///       used to in an if statement.
      friend always_inline int signbits(packed p) noexcept
      { return _mm_movemask_ps(p._m); }
      /// is any element negative
      friend always_inline bool has_negative(packed p) noexcept
      { return _mm_movemask_ps(p._m); }
      /// are all elements negative
      friend always_inline bool all_negative(packed p) noexcept
      { return _mm_movemask_ps(p._m) == max_signbits; }
      /// is any element non-zero
      friend always_inline bool has_non_zero(packed p) noexcept
      { return signbits(p!=packed::zero()); }
      /// are all elements non-zero
      friend always_inline bool all_non_zero(packed p) noexcept
      { return signbits(p!=packed::zero()) == max_signbits; }
      /// is any element zero
      friend always_inline bool has_zero(packed p) noexcept
      { return signbits(p==packed::zero()); }
      /// are all elements zero
      friend always_inline bool all_zero(packed p) noexcept
      { return signbits(p==packed::zero()) == max_signbits; }
      //@}

      /// \name miscellaneous
      //@{
      /// set all elements to Kth element of argument
      template<int K>
      friend packed single(packed p) noexcept;
      /// result = [b2,b3,a2,a3]
      friend packed always_inline movehl(packed a, packed b) noexcept
      { return packed(_mm_movehl_ps(a._m,b._m)); }
      /// result = []
      friend packed always_inline movelh(packed a, packed b) noexcept
      { return packed(_mm_movelh_ps(a._m,b._m)); }
      /// shuffle
      /// \note inverse order to _MM_SHUFFLE (here: first is first)
      template<int I0, int I1, int I2, int I3>
      friend packed shuffle(packed a, packed b) noexcept;
      /// blend two vectors depending on sign of third:  result = sign<0? x : y
      friend packed always_inline blend(packed sign, packed x, packed y)
	noexcept
      { 
# ifdef __SSE4_1__
	return packed(_mm_blendv_ps(y._m,x._m,sign._m));
# else
	// perhaps there is a better way using move and movemask?
        __m128 mask = _mm_cmplt_ps(sign._m,_mm_setzero_ps());
        return packed(_mm_or_ps(_mm_and_ps(mask,x._m),
				_mm_andnot_ps(mask,y._m)));
# endif // __SSE4_1__
      }
      /// combine two vectors depending on third:  result = mask? x : y
      /// \note All bits in mask[i] must be either 0 or 1.
      /// \note If we have SSE4, we can implement this using the blend
      ///       intrinsics, when only the highest (sign bit) in mask[i] is
      ///       required. However, this should not be relied upon. Use blend()
      ///       instead if you want to use that functionality explicitly.
      friend packed always_inline combine(packed mask, packed x, packed y)
	noexcept
      {
# ifdef __SSE4_1__
	return packed(_mm_blendv_ps(y._m,x._m,mask._m));
# else
	return packed(_mm_or_ps(_mm_and_ps(mask._m,x._m),
				_mm_andnot_ps(mask._m,y._m)));
# endif // __SSE4_1__
      }
      //@}
      static element_type always_inline nil_mask_elem() noexcept
      { return WDutils::meta::float_and_uint(0x0u).f; }
      static element_type always_inline all_mask_elem() noexcept
      { return WDutils::meta::float_and_uint(0xffffffffu).f; }
      static element_type always_inline sgn_mask_elem() noexcept
      { return WDutils::meta::float_and_uint(0x80000000u).f; }
      static element_type always_inline abs_mask_elem() noexcept
      { return WDutils::meta::float_and_uint(0x7fffffffu).f; }
    private:
      //
# ifdef __SSE2__
      static data_type always_inline nil_mask() noexcept
      { return _mm_castsi128_ps(_mm_set1_epi32(0x0)); }
      static data_type always_inline one_mask() noexcept
      { return _mm_castsi128_ps(_mm_set1_epi32(int(0xffffffff))); }
      static data_type always_inline sgn_mask() noexcept
      { return _mm_castsi128_ps(_mm_set1_epi32(int(0x80000000))); }
      static data_type always_inline abs_mask() noexcept
      { return _mm_castsi128_ps(_mm_set1_epi32(0x7fffffff)); }
# else
      static data_type always_inline nil_mask() noexcept
      { return _mm_set1_ps(WDutils::meta::float_and_uint(0x0u).f); }
      static data_type always_inline one_mask() noexcept
      { return _mm_set1_ps(WDutils::meta::float_and_uint(0xffffffffu).f); }
      static data_type always_inline sgn_mask() noexcept
      { return _mm_set1_ps(WDutils::meta::float_and_uint(0x80000000u).f); }
      static data_type always_inline abs_mask() noexcept
      { return _mm_set1_ps(WDutils::meta::float_and_uint(0x7fffffffu).f); }
# endif // __SSE2__
      /// just our sign bits
      data_type always_inline signmask() const noexcept
      { return _mm_and_ps(_m,sgn_mask()); }
      /// data
      data_type _m;
    };// SSE::packed<4,float>
    typedef packed<4,float> fvec4;
    //
    template<unsigned I> inline
    float always_inline at(fvec4 p) noexcept
    {
      static_assert(I < fvec4::block_size,"index out of range");
# if   defined(__clang__)
      return p._m[I];
# elif defined(__GNUC__) && !defined(__INTEL_COMPILER)
      return __builtin_ia32_vec_ext_v4sf(p._m,I);
# elif defined(__SSE4_1__)
      float tmp;
      _MM_EXTRACT_FLOAT(tmp,p.m,I);
      return tmp;
# else
      fvec4::aligned_element_block tmp;
      _mm_store_ps(tmp,p.m);
      return tmp[I];
# endif
    }
    //
    template<int K> inline
    fvec4 always_inline single(fvec4 p) noexcept
    {
      static_assert(K>=0 && K<fvec4::block_size,"K out of range");
      return fvec4(_mm_shuffle_ps(p._m,p._m,_MM_SHUFFLE(K,K,K,K)));
    }
    //
    template<int I0, int I1, int I2, int I3> inline
    fvec4 always_inline shuffle(fvec4 a, fvec4 b) noexcept
    {
      static_assert(I0>=0 && I0<fvec4::block_size &&
		    I1>=0 && I1<fvec4::block_size &&
		    I2>=0 && I2<fvec4::block_size &&
		    I3>=0 && I3<fvec4::block_size, "Is out of range");
      return fvec4(_mm_shuffle_ps(a._m,b._m,_MM_SHUFFLE(I3,I2,I1,I0)));
    }
  }
  //
  WDutils_TRAITS(SSE::fvec4,"fvec4");
  //
  namespace SSE {
# ifdef __AVX__
    //--------------------------------------------------------------------------
    ///
    /// 8 packed single-precision floating-point numbers
    ///
    //--------------------------------------------------------------------------
    template<> struct packed<8,float> : packed_base<8>
    {
      /// \name types, constants, and static methods
      //@
      using packed_base<8>::block_size;
      using packed_base<8>::block_shft;
      using packed_base<8>::block_trim;
      using packed_base<8>::block_mask;
      using packed_base<8>::max_signbits;
      using packed_base<8>::is_aligned;
      using packed_base<8>::aligned_index;
      using packed_base<8>::sub_index;
      using packed_base<8>::block_index;
      using packed_base<8>::num_blocks;
      using packed_base<8>::blocked_num;
      /// associated SSE/AVX vector type
      typedef __m256 data_type;
      /// associated element type
      typedef float element_type;
      /// equivalent array of elements
      typedef element_type element_block[block_size];
      /// required alignement (bytes)
      static const unsigned alignment = block_size*sizeof(element_type);
      /// aligned equivalent array of elements
      typedef element_block alignas(32) aligned_element_block;
      /// a packed with all elements equal to 0
      static packed always_inline zero() noexcept
      { return packed(_mm256_setzero_ps()); }
      /// a packed with all elements equal to 1
      static packed always_inline one() noexcept
      { return packed(_mm256_set1_ps(1.f)); }
      /// is given pointer appropriately aligned?
      static bool is_aligned(void*p) noexcept
      { return (size_t(p)&(alignment-1))==0; }
      /// offset (number of element_types) of pointer from alignment
      static size_t offset(element_type*p) noexcept
      { return (size_t(p)>>sizeof(element_type))&block_mask; }
      //@}

      /// \name construction and assignment
      //@{
      /// default ctor
      always_inline packed() WDutilsCXX11DefaultBody
#if __cplusplus >= 201103L
      /// copy ctor
      always_inline packed(packed const&) = default;
      /// copy operator
      packed& always_inline operator=(packed const&) = default;
#endif
      /// ctor from data_type
      explicit always_inline packed(data_type m) : _m(m) {}
      /// ctor from single integer value: set all element equal to single value
      explicit always_inline packed(int x) noexcept
      { _m = _mm256_set1_ps(x); }
      /// ctor from single value: set all element equal to single value
      explicit always_inline packed(float x) noexcept
      { _m = _mm256_set1_ps(x); }
      /// ctor from single value: set all element equal to single value
      explicit always_inline packed(double x) noexcept
      { _m = _mm256_set1_ps(float(x)); }
      /// ctor from 8 values: set elements
      /// \note inverse order to _mm256_set_ps(a,b,c,d,e,f,g,h)
      always_inline packed(element_type a, element_type b,
			   element_type c, element_type d,
			   element_type e, element_type f,
			   element_type g, element_type h) noexcept
      { _m = _mm256_setr_ps(a,b,c,d,e,f,g,h); }
      /// set all elements equal to zero
      packed& always_inline set_zero() noexcept
      { _m = _mm256_setzero_ps(); return*this; }
      /// set element equal to single value
      packed& always_inline set(float x) noexcept
      { _m = _mm256_set1_ps(x); return*this; }
      /// set element equal to single value
      packed& always_inline set(double x) noexcept
      { _m = _mm256_set1_ps(float(x)); return*this; }
      /// set elements
      /// \note inverse order to _mm256_set_ps(a,b,c,d,e,f,g,h)
      packed& always_inline set(element_type a, element_type b,
				element_type c, element_type d,
				element_type e, element_type f,
				element_type g, element_type h) noexcept
      { _m = _mm256_setr_ps(a,b,c,d,e,f,g,h); return*this; }
      /// ctor from fvec4
      always_inline packed(packed<4,float> a) noexcept
      { _m = _mm256_broadcast_ps(&(a.data())); }
//       /// ctor from 2 fvec4
//       always_inline packed(packed<4,float> a, packed<4,float> b) noexcept
      //@}

      /// \name data access and conversion
      //@{
      /// conversion to const vector type
      always_inline operator data_type const&() const noexcept
      { return _m; }
      /// direct data const access
      data_type const& always_inline data() const noexcept
      { return _m; }
      /// conversion to vector type
      always_inline operator data_type&() noexcept
      { return _m; }
      /// direct non-const data access
      data_type& always_inline data() noexcept
      { return _m; }
      /// obtain lower packed<4,float>
      friend packed<4,float> always_inline lower(packed p) noexcept
      { return packed<4,float>(_mm256_extractf128_ps(p._m,0)); }
      /// obtain upper packed<4,float>
      friend packed<4,float> always_inline upper(packed p) noexcept
      { return packed<4,float>(_mm256_extractf128_ps(p._m,1)); }
      /// obtain lower (I=0) or upper (I=1) packed<4,float>
      template<unsigned I>
      friend packed<4,float> extract(packed p) noexcept;
      /// constant element access, templated
      template<unsigned I>
      friend element_type at(packed p) noexcept;
      /// downcast: convert 2 @c packed<4,double> to  @c packed<8,float>
      friend packed downcast(packed<4,double>, packed<4,double>) noexcept;
      //@}

      /// \name load and store
      //@{
      /// load from aligned memory location
      static packed always_inline load(const element_type*p) noexcept
      { return packed(_mm256_load_ps(p)); }
      /// load from unaligned memory location
      static packed always_inline loadu(const element_type*p) noexcept
      { return packed(_mm256_loadu_ps(p)); }
      /// load from aligned memory location, using template arg for alignment
      template<bool aligned> static typename enable_if< aligned, packed>::type
      always_inline load_t(const element_type*p) noexcept
      { return packed(_mm256_load_ps(p)); }
      /// load from unaligned memory location, using template arg for alignment
      template<bool aligned> static typename enable_if<!aligned,packed>::type
      always_inline load_t(const element_type*p) noexcept
      { return packed(_mm256_loadu_ps(p)); }
      /// store to aligned memory location
      void always_inline store(element_type*p) const noexcept
      { _mm256_store_ps(p,_m); }
      /// store to unaligned memory location
      void always_inline storeu(element_type*p) const noexcept
      { _mm256_storeu_ps(p,_m); }
      /// store to aligned memory location, using template arg for alignment
      template<bool aligned> typename enable_if< aligned>::type
      always_inline store_t(element_type*p) const noexcept
      { _mm256_store_ps(p,_m); }
      /// store to unaligned memory location, using template arg for alignment
      template<bool aligned> typename enable_if<!aligned>::type
      always_inline store_t(element_type*p) const noexcept
      { _mm256_storeu_ps(p,_m); }
      /// load aligned object with member data() returning const element_type*
      template<typename class_with_member_data>
      static packed always_inline pack(class_with_member_data const&a) noexcept
      { return load(a.data()); }
      /// load unaligned object with member data() returning const element_type*
      template<typename class_with_member_data>
      static packed always_inline packu(class_with_member_data const&a) noexcept
      { return loadu(a.data()); }
      /// load object with member data() returning const element_type*
      template<bool aligned, typename class_with_member_data>
      static packed always_inline pack_t(class_with_member_data const&a)
	noexcept
      { return load_t<aligned>(a.data()); }
      /// store to aligned object with member data() returning element_type*
      template<typename class_with_member_data>
      void always_inline unpack(class_with_member_data&a) const noexcept
      { store(a.data()); }
      /// store to aligned object with member data() returning element_type*
      template<typename class_with_member_data>
      void always_inline unpacku(class_with_member_data&a) const noexcept
      { storeu(a.data()); }
      /// store to aligned object with member data() returning element_type*
      template<bool aligned, typename class_with_member_data>
      void always_inline unpack_t(class_with_member_data&a) const noexcept
      { store_t<aligned>(a.data()); }
      /// store directly to aligned memory without polluting cashes
      void always_inline stream(element_type*p) const noexcept
      { _mm256_stream_ps(p,_m); }
      //@}

      /// \name vectorised arithmetic operations
      //@{
      /// +=
      packed& always_inline operator+=(packed p) noexcept
      { _m = _mm256_add_ps(p._m,_m); return*this; }
      /// +
      packed always_inline operator+ (packed p) const noexcept
      { return packed(_mm256_add_ps(_m,p._m)); }
      /// -=
      packed& always_inline operator-=(packed p) noexcept
      { _m = _mm256_sub_ps(_m,p._m); return*this; }
      /// -
      packed always_inline operator- (packed p) const noexcept
      { return packed(_mm256_sub_ps(_m,p._m)); }
      /// unary -
      packed always_inline operator-() const noexcept
      { return packed(_mm256_xor_ps(_m,sgn_mask())); }
      /// *=
      packed& always_inline operator*=(packed p) noexcept
      { _m = _mm256_mul_ps(p._m,_m); return*this; }
      /// *
      packed always_inline operator* (packed p) const noexcept
      { return packed(_mm256_mul_ps(_m,p._m)); }
      /// /=
      packed& always_inline operator/=(packed p) noexcept
      { _m = _mm256_div_ps(_m,p._m); return*this; }
      /// /
      packed always_inline operator/ (packed p) const noexcept
      { return packed(_mm256_div_ps(_m,p._m)); }
      /// sqrt
      friend packed always_inline sqrt(packed p) noexcept
      { return packed(_mm256_sqrt_ps(p._m)); }
      /// square
      friend packed always_inline square(packed p) noexcept
      { return  packed(_mm256_mul_ps(p._m,p._m)); }
      /// reciprocal
      friend packed always_inline reciprocal(packed p) noexcept
      { return packed(_mm256_div_ps(_mm256_set1_ps(1.f),p._m)); }
      /// approximate reciprocal
      friend packed always_inline approximate_reciprocal(packed p) noexcept
      { return packed(_mm256_rcp_ps(p._m)); }
      /// maximum of two packed
      friend packed always_inline max(packed a, packed b) noexcept
      { return packed(_mm256_max_ps(a._m,b._m)); }
      /// minimum of two packed
      friend packed always_inline min(packed a, packed b) noexcept
      { return packed(_mm256_min_ps(a._m,b._m)); }
      /// abs
      friend packed always_inline abs(packed p) noexcept
      { return packed(_mm256_and_ps(p._m,abs_mask())); }
      /// -abs
      friend packed always_inline negabs(packed p) noexcept
      { return packed(_mm256_or_ps(p._m,sgn_mask())); }
      /// abs(x-y)
      friend packed always_inline diff(packed a, packed b) noexcept
      { return packed(_mm256_and_ps(_mm256_sub_ps(a._m,b._m),abs_mask())); }
      /// x = abs(x-y)
      packed& always_inline make_diff(packed p) noexcept
      {
	_m = _mm256_sub_ps(_m,p._m);
	_m = _mm256_and_ps(_m,abs_mask());
	return*this;
      }
      /// sign(a)*b
      friend packed always_inline signmove(packed a, packed b) noexcept
      { return packed(_mm256_or_ps(a.signmask(),
				   _mm256_and_ps(b._m,abs_mask()))); }
      //@}
      /// \name horizontal arithmetic operations
      //@{
      /// horizontal sum
      friend element_type always_inline sum(packed p) noexcept
      { return sum(lower(p)+upper(p)); }
#  if(0)
      /// horizontal sum in each element
      friend packed always_inline hadd(packed p) noexcept
#  endif
      /// horizontal maximum
      friend element_type always_inline max(packed p) noexcept
      { return max(max(lower(p),upper(p))); }
      /// horizontal minimum
      friend element_type always_inline min(packed p) noexcept
      { return min(min(lower(p),upper(p))); }
      //@}
      /// \name vectorised comparisons
      //@{
      /// <
      packed always_inline operator< (packed p) const noexcept
      { return packed(_mm256_cmp_ps(_m,p._m,_CMP_LT_OQ)); }
      /// <=
      packed always_inline operator<=(packed p) const noexcept
      { return packed(_mm256_cmp_ps(_m,p._m,_CMP_LE_OQ)); }
      /// >
      packed always_inline operator> (packed p) const noexcept
      { return packed(_mm256_cmp_ps(_m,p._m,_CMP_GT_OQ)); }
      /// >=
      packed always_inline operator>=(packed p) const noexcept
      { return packed(_mm256_cmp_ps(_m,p._m,_CMP_GE_OQ)); }
      /// ==
      packed always_inline operator==(packed p) const noexcept
      { return packed(_mm256_cmp_ps(_m,p._m,_CMP_EQ_UQ)); }
      /// !=
      packed always_inline operator!=(packed p) const noexcept
      { return packed(_mm256_cmp_ps(_m,p._m,_CMP_NEQ_UQ)); }
      //@}

      /// \name vectorised logical operations
      /// &=
      packed& always_inline operator&=(packed p) noexcept
      { _m = _mm256_and_ps(p._m,_m); return*this; }
      /// &
      packed always_inline operator& (packed p) const noexcept
      { return packed(_mm256_and_ps(_m,p._m)); }
      /// |=
      packed& always_inline operator|=(packed p) noexcept
      { _m = _mm256_or_ps(p._m,_m); return*this; }
      /// |
      packed always_inline operator| (packed p) const noexcept
      { return packed(_mm256_or_ps(_m,p._m)); }
      /// ^=
      packed& always_inline operator^=(packed p) noexcept
      { _m = _mm256_xor_ps(p._m,_m); return*this; }
      /// ^
      packed always_inline operator^ (packed p) const noexcept
      { return packed(_mm256_xor_ps(_m,p._m)); }
      /// unary !
      packed always_inline operator! () const noexcept
      { return packed(_mm256_xor_ps(_m,one_mask())); }
      /// and not: this=this&!that
      packed always_inline andnot(packed p) const noexcept
      { return packed(_mm256_andnot_ps(p._m,_m)); }
      //@}

      /// \name boolean or integer properties
      //@{
      /// return integer with bits equal to sign bits
      /// \note if @a p is the result of a boolean operation (comparison, or
      ///       bit-wise operations on comparison results), signbit(p) can be
      ///       used to in an if statement.
      friend always_inline int signbits(packed p) noexcept
      { return _mm256_movemask_ps(p._m); }
      /// is any element negative
      friend always_inline bool has_negative(packed p) noexcept
      { return _mm256_movemask_ps(p._m); }
      /// are all elements negative
      friend always_inline bool all_negative(packed p) noexcept
      { return _mm256_movemask_ps(p._m) == max_signbits; }
      /// is any element non-zero
      friend always_inline bool has_non_zero(packed p) noexcept
      { return signbits(p!=packed::zero()); }
      /// are all elements non-zero
      friend always_inline bool all_non_zero(packed p) noexcept
      { return signbits(p!=packed::zero()) == max_signbits; }
      /// is any element zero
      friend always_inline bool has_zero(packed p) noexcept
      { return signbits(p==packed::zero()); }
      /// are all elements zero
      friend always_inline bool all_zero(packed p) noexcept
      { return signbits(p==packed::zero()) == max_signbits; }
      //@}

      /// \name miscellaneous
      //@{
      /// blend two vectors depending on sign of third:  result = sign<0? x : y
      friend packed always_inline blend(packed sign, packed x, packed y)
	noexcept
      { return packed(_mm256_blendv_ps(y._m,x._m,sign._m)); }
      /// combine two vectors depending on third:  result = mask? x : y
      friend packed always_inline combine(packed mask, packed x, packed y)
	noexcept
      { return packed(_mm256_blendv_ps(y._m,x._m,mask._m)); }
      //@}
      static element_type always_inline nil_mask_elem() noexcept
      { return WDutils::meta::float_and_uint(0x0u).f; }
      static element_type always_inline all_mask_elem() noexcept
      { return WDutils::meta::float_and_uint(0xffffffffu).f; }
      static element_type always_inline sgn_mask_elem() noexcept
      { return WDutils::meta::float_and_uint(0x80000000u).f; }
      static element_type always_inline abs_mask_elem() noexcept
      { return WDutils::meta::float_and_uint(0x7fffffffu).f; }
    private:
      static data_type always_inline nil_mask() noexcept
      { return _mm256_castsi256_ps(_mm256_set1_epi32(0x0)); }
      static data_type always_inline one_mask() noexcept
      { return _mm256_castsi256_ps(_mm256_set1_epi32(int(0xffffffff))); }
      static data_type always_inline sgn_mask() noexcept
      { return _mm256_castsi256_ps(_mm256_set1_epi32(int(0x80000000))); }
      static data_type always_inline abs_mask() noexcept
      { return _mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff)); }
      /// just our sign bits
      data_type always_inline signmask() const noexcept
      { return _mm256_and_ps(_m,sgn_mask()); }
      /// data
      data_type _m;
    };// SSE::packed<8,float>
    //
    typedef packed<8,float> fvec8;
    //
    template<unsigned I> inline
    fvec4 always_inline extract(fvec8 p) noexcept
    {
      static_assert(I<2,"index out of range");
      return fvec4(_mm256_extractf128_ps(p._m,I));
    }
    //
    template<unsigned I> inline
    float always_inline at(fvec8 p) noexcept
    {
      static_assert(I<fvec8::block_size,"index out of range");
#  if   defined(__clang__)
      return p._m[I];
#  else
      return SSE::at<(I&3)>(extract<(I>>2)>(p));
#  endif
    }
  }
  //
  WDutils_TRAITS(SSE::fvec8,"fvec8");
  //
  namespace SSE {
# endif // __AVX__

# ifdef __SSE2__
    //--------------------------------------------------------------------------
    ///
    /// 2 packed double-precision floating-point numbers
    ///
    //--------------------------------------------------------------------------
    template<> struct packed<2,double> : packed_base<2>
    {
    public:
      /// \name types, constants, and static methods
      //@
      using packed_base<2>::block_size;
      using packed_base<2>::block_shft;
      using packed_base<2>::block_trim;
      using packed_base<2>::block_mask;
      using packed_base<2>::max_signbits;
      using packed_base<2>::is_aligned;
      using packed_base<2>::aligned_index;
      using packed_base<2>::sub_index;
      using packed_base<2>::block_index;
      using packed_base<2>::num_blocks;
      using packed_base<2>::blocked_num;
      /// associated element type
      typedef double element_type;
      /// associated SSE/AVX vector type
      typedef __m128d data_type;
      /// equivalent array of elements
      typedef element_type element_block[block_size];
      /// required alignement (bytes)
      static const unsigned alignment = block_size*sizeof(element_type);
      /// aligned equivalent array of elements
      typedef element_block alignas(16) aligned_element_block;
      /// a packed with all elements equal to 0
      static packed always_inline zero() noexcept
      { return packed(_mm_setzero_pd()); }
      /// a packed with all elements equal to 1
      static packed always_inline one() noexcept
      { return packed(_mm_set1_pd(1.0)); }
      /// is given pointer appropriately aligned?
      static bool is_aligned(void*p) noexcept
      { return (size_t(p)&(alignment-1))==0; }
      /// offset (number of element_types) of pointer from alignment
      static size_t offset(element_type*p) noexcept
      { return (size_t(p)>>sizeof(element_type))&block_mask; }
      //@}

      /// \name construction and assignment
      //@{
      /// default ctor
      always_inline packed() WDutilsCXX11DefaultBody
#if __cplusplus >= 201103L
      /// copy ctor
      always_inline packed(packed const&) = default;
      /// copy operator
      packed& always_inline operator=(packed const&) = default;
#endif
      /// ctor from data_type is private
      explicit always_inline packed(__m128d m) : _m(m) {}
      /// ctor from single integer value: set all element equal to single value
      explicit always_inline packed(int x) noexcept
      { _m = _mm_set1_pd(x); }
      /// ctor from single value: set all element equal to single value
      explicit always_inline packed(float x) noexcept
      { _m = _mm_set1_pd(double(x)); }
      /// ctor from single value: set all element equal to single value
      explicit always_inline packed(double x) noexcept
      { _m = _mm_set1_pd(x); }
      /// ctor from 2 values: set elements
      /// \note inverse order to _mm_set_pd(y,x)
      always_inline packed(element_type x, element_type y) noexcept
      { _m = _mm_set_pd(y,x); }
      /// set all elements equal to zero
      packed& always_inline set_zero() noexcept
      { _m = _mm_setzero_pd(); return*this; }
      /// set element equal to single value
      packed& always_inline set(element_type x) noexcept
      { _m = _mm_set1_pd(x); return*this; }
      /// set elements
      /// \note inverse order to _mm_set_pd(y,x)
      packed& always_inline set(element_type x, element_type y) noexcept
      { _m = _mm_set_pd(y,x); return*this; }
      //@}

      /// \name data access and conversion
      //@{
      /// conversion to const vector type
      always_inline operator data_type const&() const noexcept
      { return _m; }
      /// direct data const access
      data_type const& always_inline data() const noexcept
      { return _m; }
      /// conversion to vector type
      always_inline operator data_type&() noexcept
      { return _m; }
      /// direct non-const data access
      data_type& always_inline data() noexcept
      { return _m; }
      /// constant element access, templated
      template<unsigned I>
      friend element_type at(packed p) noexcept;
      /// downcast: convert two @c packed<2,double> to one @c packed<4,float>
      friend packed<4,float> always_inline downcast(packed a, packed b) noexcept
      { return packed<4,float>(_mm_movelh_ps(_mm_cvtpd_ps(a._m),
					     _mm_cvtpd_ps(b._m))); }
      /// upcast: convert lower two of packed<4,float> to packed<2,double>
      friend packed always_inline upcast_lo(packed<4,float> p) noexcept
      { return packed(_mm_cvtps_pd(p._m)); }
      /// upcast: convert upper two of packed<4,float> to packed<2,double>
      friend packed always_inline upcast_hi(packed<4,float> p) noexcept
      { return packed(_mm_cvtps_pd(_mm_movehl_ps(p._m,p._m))); }
      //@}

      /// \name load and store
      //@{
      /// load from aligned memory location
      static packed always_inline load(const element_type*p) noexcept
      { return packed(_mm_load_pd(p)); }
      /// load from unaligned memory location
      static packed always_inline loadu(const element_type*p) noexcept
      { return packed(_mm_loadu_pd(p)); }
      /// load from aligned memory location, using template arg for alignment
      template<bool aligned> static typename enable_if< aligned, packed>::type
      always_inline load_t(const element_type*p) noexcept
      { return packed(_mm_load_pd(p)); }
      /// load from unaligned memory location, using template arg for alignment
      template<bool aligned> static typename enable_if<!aligned,packed>::type
      always_inline load_t(const element_type*p) noexcept
      { return packed(_mm_loadu_pd(p)); }
      /// store to aligned memory location
      void always_inline store(element_type*p) const noexcept
      { _mm_store_pd(p,_m); }
      /// store to unaligned memory location
      void always_inline storeu(element_type*p) const noexcept
      { _mm_storeu_pd(p,_m); }
      /// store to aligned memory location, using template arg for alignment
      template<bool aligned> typename enable_if< aligned>::type
      always_inline store_t(element_type*p) const noexcept
      { _mm_store_pd(p,_m); }
      /// store to unaligned memory location, using template arg for alignment
      template<bool aligned> typename enable_if<!aligned>::type
      always_inline store_t(element_type*p) const noexcept
      { _mm_storeu_pd(p,_m); }
      /// load aligned object with member data() returning const element_type*
      template<typename class_with_member_data>
      static packed always_inline pack(class_with_member_data const&a) noexcept
      { return load(a.data()); }
      /// load unaligned object with member data() returning const element_type*
      template<typename class_with_member_data>
      static packed always_inline packu(class_with_member_data const&a) noexcept
      { return loadu(a.data()); }
      /// load object with member data() returning const element_type*
      template<bool aligned, typename class_with_member_data>
      static packed always_inline pack_t(class_with_member_data const&a)
	noexcept
      { return load_t<aligned>(a.data()); }
      /// store to aligned object with member data() returning element_type*
      template<typename class_with_member_data>
      void always_inline unpack(class_with_member_data&a) const noexcept
      { store(a.data()); }
      /// store to aligned object with member data() returning element_type*
      template<typename class_with_member_data>
      void always_inline unpacku(class_with_member_data&a) const noexcept
      { storeu(a.data()); }
      /// store to aligned object with member data() returning element_type*
      template<bool aligned, typename class_with_member_data>
      void always_inline unpack_t(class_with_member_data&a) const noexcept
      { store_t<aligned>(a.data()); }
      /// store directly to aligned memory without polluting cashes
      void always_inline stream(element_type*p) const noexcept
      { _mm_stream_pd(p,_m); }
      //@}

      /// \name vectorised arithmetic operations
      //@{
      /// +=
      packed& always_inline operator+=(packed p) noexcept
      { _m = _mm_add_pd(p._m,_m); return*this; }
      /// +
      packed always_inline operator+ (packed p) const noexcept
      { return packed(_mm_add_pd(_m,p._m)); }
      /// -=
      packed& always_inline operator-=(packed p) noexcept
      { _m = _mm_sub_pd(_m,p._m); return*this; }
      /// -
      packed always_inline operator- (packed p) const noexcept
      { return packed(_mm_sub_pd(_m,p._m)); }
      /// unary -
      packed always_inline operator-() const noexcept
      { return packed(_mm_xor_pd(_m,sgn_mask())); }
      /// *=
      packed& always_inline operator*=(packed p) noexcept
      { _m = _mm_mul_pd(p._m,_m); return*this; }
      /// *
      packed always_inline operator* (packed p) const noexcept
      { return packed(_mm_mul_pd(_m,p._m)); }
      /// /=
      packed& always_inline operator/=(packed p) noexcept
      { _m = _mm_div_pd(_m,p._m); return*this; }
      /// /
      packed always_inline operator/ (packed p) const noexcept
      { return packed(_mm_div_pd(_m,p._m)); }
      /// sqrt
      friend packed always_inline sqrt(packed p) noexcept
      { return packed(_mm_sqrt_pd(p._m)); }
      /// square
      friend packed always_inline square(packed p) noexcept
      { return  packed(_mm_mul_pd(p._m,p._m)); }
      /// reciprocal
      friend packed always_inline reciprocal(packed p) noexcept
      { return packed(_mm_div_pd(_mm_set1_pd(1.0),p._m)); }
      /// maximum of two packed
      friend packed always_inline max(packed a, packed b) noexcept
      { return packed(_mm_max_pd(a._m,b._m)); }
      /// minimum of two packed
      friend packed always_inline min(packed a, packed b) noexcept
      { return packed(_mm_min_pd(a._m,b._m)); }
      /// abs
      friend packed always_inline abs(packed p) noexcept
      { return packed(_mm_and_pd(p._m,abs_mask())); }
      /// -abs
      friend packed always_inline negabs(packed p) noexcept
      { return packed(_mm_or_pd(p._m,sgn_mask())); }
      /// abs(x-y)
      friend packed always_inline diff(packed a, packed b) noexcept
      { return packed(_mm_and_pd(_mm_sub_pd(a._m,b._m),abs_mask())); }
      /// x = abs(x-y)
      packed& always_inline make_diff(packed p) noexcept
      {
	_m = _mm_sub_pd(_m,p._m);
	_m = _mm_and_pd(_m,abs_mask());
	return*this;
      }
      /// sign(a)*b
      friend packed always_inline signmove(packed a, packed b) noexcept
      { return packed(_mm_or_pd(a.signmask(),_mm_and_pd(b._m,abs_mask()))); }
      //@}

      /// \name horizontal arithmetic operations
      //@{
      /// horizontal sum
      friend element_type always_inline sum(packed p) noexcept
      {
	aligned_element_block q;
	_mm_store_pd(q,p._m);
	return q[0]+q[1];
      }
      /// horizontal sum of first element
      friend element_type always_inline sum1(packed p) noexcept
      { return at<0>(p); }
#  if(0)
      /// horizontal sum in each element
      friend packed always_inline hadd(packed p) noexcept
#  endif
      /// horizontal maximum
      friend element_type always_inline max(packed p) noexcept
      {
	aligned_element_block q;
	_mm_store_pd(q,p._m);
	return q[0]>q[1]? q[0]:q[1];
      }
      /// horizontal minimum
      friend element_type always_inline min(packed p) noexcept
      {
	aligned_element_block q;
	_mm_store_pd(q,p._m);
	return q[0]<q[1]? q[0]:q[1];
      }
      //@}

      /// \name vectorised comparisons
      //@{
      /// <
      packed always_inline operator< (packed p) const noexcept
      { return packed(_mm_cmplt_pd(_m,p._m)); }
      /// <=
      packed always_inline operator<=(packed p) const noexcept
      { return packed(_mm_cmple_pd(_m,p._m)); }
      /// >
      packed always_inline operator> (packed p) const noexcept
      { return packed(_mm_cmpgt_pd(_m,p._m)); }
      /// >=
      packed always_inline operator>=(packed p) const noexcept
      { return packed(_mm_cmpge_pd(_m,p._m)); }
      /// ==
      packed always_inline operator==(packed p) const noexcept
      { return packed(_mm_cmpeq_pd(_m,p._m)); }
      /// !=
      packed always_inline operator!=(packed p) const noexcept
      { return packed(_mm_cmpneq_pd(_m,p._m)); }
      //@}

      /// \name vectorised logical operations
      /// &=
      packed&always_inline  operator&=(packed p) noexcept
      { _m = _mm_and_pd(p._m,_m); return*this; }
      /// &
      packed always_inline  operator& (packed p) const noexcept
      { return packed(_mm_and_pd(_m,p._m)); }
      /// |=
      packed&always_inline  operator|=(packed p) noexcept
      { _m = _mm_or_pd(p._m,_m); return*this; }
      /// |
      packed always_inline  operator| (packed p) const noexcept
      { return packed(_mm_or_pd(_m,p._m)); }
      /// ^=
      packed&always_inline  operator^=(packed p) noexcept
      { _m = _mm_xor_pd(p._m,_m); return*this; }
      /// ^
      packed always_inline  operator^ (packed p) const noexcept
      { return packed(_mm_xor_pd(_m,p._m)); }
      /// unary !
      packed always_inline operator! () const noexcept
      { return packed(_mm_xor_pd(_m,one_mask())); }
      /// and not: this=this&!that
      packed always_inline andnot(packed p) const noexcept
      { return packed(_mm_andnot_pd(p._m,_m)); }
      //@}

      /// \name boolean or integer properties
      //@{
      /// return integer with bits equal to sign bits
      /// \note if @a p is the result of a boolean operation (comparison, or
      ///       bit-wise operations on comparison results), signbit(p) can be
      ///       used to in an if statement.
      friend always_inline int signbits(packed p) noexcept
      { return _mm_movemask_pd(p._m); }
      /// is any element negative
      friend always_inline bool has_negative(packed p) noexcept
      { return _mm_movemask_pd(p._m); }
      /// are all elements negative
      friend always_inline bool all_negative(packed p) noexcept
      { return _mm_movemask_pd(p._m) == max_signbits; }
      /// is any element non-zero
      friend always_inline bool has_non_zero(packed p) noexcept
      { return signbits(p!=packed::zero()); }
      /// are all elements non-zero
      friend always_inline bool all_non_zero(packed p) noexcept
      { return signbits(p!=packed::zero()) == max_signbits; }
      /// is any element zero
      friend always_inline bool has_zero(packed p) noexcept
      { return signbits(p==packed::zero()); }
      /// are all elements zero
      friend always_inline bool all_zero(packed p) noexcept
      { return signbits(p==packed::zero()) == max_signbits; }
      //@}

      /// \name miscellaneous
      //@{
      /// set all elements to Kth element of argument
      template<int K>
      friend packed always_inline single(packed p) noexcept;
      /// blend two vectors depending on sign of third: result = sign<0? x : y
      friend packed always_inline blend(packed sign, packed x, packed y)
	noexcept
      { 
#  ifdef __SSE4_1__
	return packed(_mm_blendv_pd(y._m,x._m,sign._m));
#  else
	// perhapd there is a better way using move and movemask?
        __m128d mask = _mm_cmplt_pd(sign._m,_mm_setzero_pd());
        return packed(_mm_or_pd(_mm_and_pd(mask,x._m),
				_mm_andnot_pd(mask,y._m)));
#  endif
      }
      /// combine two vectors depending on third:  result = mask? x : y
      /// \note All bits in mask[i] must be either 0 or 1.
      /// \note If we have SSE4, we can implement this using the blend
      ///       intrinsics, when only the highest (sign bit) in mask[i] is
      ///       required. However, this should not be relied upon. Use blend()
      ///       instead if you want to use that functionality explicitly.
      friend packed always_inline combine(packed mask, packed x, packed y)
	noexcept
      {
#  ifdef __SSE4_1__
	return packed(_mm_blendv_pd(y._m,x._m,mask._m));
#  else
	return packed(_mm_or_pd(_mm_and_pd(mask._m,x._m),
				_mm_andnot_pd(mask._m,y._m)));
#  endif
      }
      //@}
      static element_type always_inline nil_mask_elem() noexcept
      { return WDutils::meta::double_and_uint(0x0lu).d; }
      static element_type always_inline all_mask_elem() noexcept
      { return WDutils::meta::double_and_uint(0xfffffffffffffffflu).d; }
      static element_type always_inline sgn_mask_elem() noexcept
      { return WDutils::meta::double_and_uint(0x8000000000000000lu).d; }
      static element_type always_inline abs_mask_elem() noexcept
      { return WDutils::meta::double_and_uint(0x7ffffffffffffffflu).d; }
    private:
      static __m128d always_inline nil_mask() noexcept
      { return _mm_castsi128_pd(_mm_set1_epi64x(0x0l)); }
      static __m128d always_inline one_mask() noexcept
      { return _mm_castsi128_pd(_mm_set1_epi64x(int64_t(0xffffffffffffffff))); }
      static __m128d always_inline sgn_mask() noexcept
      { return _mm_castsi128_pd(_mm_set1_epi64x(int64_t(0x8000000000000000))); }
      static __m128d always_inline abs_mask() noexcept
      { return _mm_castsi128_pd(_mm_set1_epi64x(0x7fffffffffffffff)); }
      /// just our sign bits
      data_type always_inline signmask() const noexcept
      { return _mm_and_pd(_m,sgn_mask()); }
      /// data
      data_type _m;
    };// SSE::packed<2,double>
    //
    typedef packed<2,double> dvec2;
    //
    template<unsigned I> inline
    double always_inline at(dvec2 p) noexcept
    {
      static_assert(I<dvec2::block_size,"index out of range");
# if   defined(__clang__)
      return p._m[I];
# elif defined(__GNUC__) && !defined(__INTEL_COMPILER)
      return __builtin_ia32_vec_ext_v2df(p._m,I);
# else
      dvec2::aligned_element_block q;
      _mm_store_pd(q,p._m);
      return q[I];
# endif
    }
    //
    template<int K> inline
    dvec2 always_inline single(dvec2 p) noexcept
    {
      static_assert(K>=0 && K<dvec2::block_size,"K out of range");
      return dvec2(_mm_shuffle_pd(p._m,p._m,_MM_SHUFFLE2(K,K)));
    }
  }
  //
  WDutils_TRAITS(SSE::dvec2,"dvec2");
  //
  namespace SSE {
# endif // __SSE2__

# ifdef __AVX__
    //--------------------------------------------------------------------------
    ///
    /// 4 packed double-precision floating-point numbers
    ///
    //--------------------------------------------------------------------------
    template<> struct packed<4,double> : packed_base<4>
    {
      /// \name types, constants, and static methods
      //@
      using packed_base<4>::block_size;
      using packed_base<4>::block_shft;
      using packed_base<4>::block_trim;
      using packed_base<4>::block_mask;
      using packed_base<4>::max_signbits;
      using packed_base<4>::is_aligned;
      using packed_base<4>::aligned_index;
      using packed_base<4>::sub_index;
      using packed_base<4>::block_index;
      using packed_base<4>::num_blocks;
      using packed_base<4>::blocked_num;
      /// associated element type
      typedef double element_type;
      /// associated SSE/AVX vector type
      typedef __m256d data_type;
      /// equivalent array of elements
      typedef element_type element_block[block_size];
      /// required alignement (bytes)
      static const unsigned alignment = block_size*sizeof(element_type);
      /// aligned equivalent array of elements
      typedef element_block alignas(32) aligned_element_block;
      /// a packed with all elements equal to 0
      static packed always_inline zero() noexcept
      { return packed(_mm256_setzero_pd()); }
      /// a packed with all elements equal to 1
      static packed always_inline one() noexcept
      { return packed(_mm256_set1_pd(1.0)); }
      /// is given pointer appropriately aligned?
      static bool is_aligned(void*p) noexcept
      { return (size_t(p)&(alignment-1))==0; }
      /// offset (number of element_types) of pointer from alignment
      static size_t offset(element_type*p) noexcept
      { return (size_t(p)>>sizeof(element_type))&block_mask; }
      //@}

      /// \name construction and assignment
      //@{
      /// default ctor
      always_inline packed() WDutilsCXX11DefaultBody
#if __cplusplus >= 201103L
      /// copy ctor
      always_inline packed(packed const&) = default;
      /// copy operator
      packed& always_inline operator=(packed const&) = default;
#endif
      /// ctor from data_type is private
      explicit always_inline packed(data_type m) : _m(m) {}
      /// ctor from single integer value: set all element equal to single value
      explicit always_inline packed(int x) noexcept
      { _m = _mm256_set1_pd(x); }
      /// ctor from single value: set all element equal to single value
      explicit always_inline packed(float x) noexcept
      { _m = _mm256_set1_pd(float(x)); }
      /// ctor from single value: set all element equal to single value
      explicit always_inline packed(double x) noexcept
      { _m = _mm256_set1_pd(x); }
      /// ctor from 4 values: set elements
      /// \note inverse ordre to _mm256_set_pd()
      always_inline packed(element_type a, element_type b,
			   element_type c, element_type d) noexcept
      { _m = _mm256_set_pd(d,c,b,a); }
      /// ctor from dvec2
      always_inline packed(packed<2,double> a) noexcept
      { _m = _mm256_broadcast_pd(&(a.data())); }
      /// set all elements equal to zero
      packed& always_inline set_zero() noexcept
      { _m = _mm256_setzero_pd(); return*this; }
      /// set element equal to single value
      packed& always_inline set(element_type x) noexcept
      { _m = _mm256_set1_pd(x); return*this; }
      /// set elements
      /// \note inverse ordre to _mm256_set_pd()
      packed& always_inline set(element_type a, element_type b,
				element_type c, element_type d) noexcept
      { _m = _mm256_set_pd(d,c,b,a); return*this; }
      //@}

      /// \name data access and conversion
      //@{
      /// conversion to const vector type
      always_inline operator data_type const&() const noexcept
      { return _m; }
      /// direct data const access
      data_type const& always_inline data() const noexcept
      { return _m; }
      /// conversion to vector type
      always_inline operator data_type&() noexcept
      { return _m; }
      /// direct non-const data access
      data_type& always_inline data() noexcept
      { return _m; }
      /// obtain lower (I=0) or upper (I=1) packed<2,double>
      template<unsigned I>
      friend packed<2,double> always_inline extract(packed p) noexcept;
      /// obtain lower packed<2,double>
      friend packed<2,double> always_inline lower(packed p) noexcept
      { return packed<2,double>(_mm256_extractf128_pd(p._m,0)); }
      /// obtain upper packed<2,double>
      friend packed<2,double> always_inline upper(packed p) noexcept
      { return packed<2,double>(_mm256_extractf128_pd(p._m,1)); }
      /// constant element access, templated
      template<unsigned I>
      friend element_type always_inline at(packed p) noexcept;
      /// downcast: convert 2 @c packed<4,double> to @c packed<8,float>
      friend packed<8,float> always_inline downcast(packed a, packed b) noexcept
      { return packed<8,float>(_mm256_insertf128_ps
			       (_mm256_castps128_ps256(_mm256_cvtpd_ps(a._m)),
				_mm256_cvtpd_ps(b._m),1)); }
      /// downcast: convert @c packed<4,double>  to  @c packed<4,float>
      friend packed<4,float> always_inline downcast(packed a) noexcept
      { return packed<4,float>(_mm256_cvtpd_ps(a._m)); }
      /// upcast:   convert  @c packed<4,float>  to  @c packed<4,double>  
      friend packed always_inline upcast(packed<4,float> a) noexcept
      { return packed(_mm256_cvtps_pd(a._m)); }
      //@}

      /// \name load and store
      //@{
      /// load from aligned memory location
      static packed always_inline load(const element_type*p) noexcept
      { return packed(_mm256_load_pd(p)); }
      /// load from unaligned memory location
      static packed always_inline loadu(const element_type*p) noexcept
      { return packed(_mm256_loadu_pd(p)); }
      /// load from aligned memory location, using template arg for alignment
      template<bool aligned> static typename enable_if< aligned, packed>::type
      always_inline load_t(const element_type*p) noexcept
      { return packed(_mm256_load_pd(p)); }
      /// load from unaligned memory location, using template arg for alignment
      template<bool aligned> static typename enable_if<!aligned,packed>::type
      always_inline load_t(const element_type*p) noexcept
      { return packed(_mm256_loadu_pd(p)); }
      /// store to aligned memory location
      void always_inline store(element_type*p) const noexcept
      { _mm256_store_pd(p,_m); }
      /// store to unaligned memory location
      void always_inline storeu(element_type*p) const noexcept
      { _mm256_storeu_pd(p,_m); }
      /// store to aligned memory location, using template arg for alignment
      template<bool aligned> typename enable_if< aligned>::type
      always_inline store_t(element_type*p) const noexcept
      { _mm256_store_pd(p,_m); }
      /// store to unaligned memory location, using template arg for alignment
      template<bool aligned> typename enable_if<!aligned>::type
      always_inline store_t(element_type*p) const noexcept
      { _mm256_storeu_pd(p,_m); }
      /// load aligned object with member data() returning const element_type*
      template<typename class_with_member_data>
      static packed always_inline pack(class_with_member_data const&a) noexcept
      { return load(a.data()); }
      /// load unaligned object with member data() returning const element_type*
      template<typename class_with_member_data>
      static packed always_inline packu(class_with_member_data const&a) noexcept
      { return loadu(a.data()); }
      /// load object with member data() returning const element_type*
      template<bool aligned, typename class_with_member_data>
      static packed always_inline pack_t(class_with_member_data const&a)
	noexcept
      { return load_t<aligned>(a.data()); }
      /// store to aligned object with member data() returning element_type*
      template<typename class_with_member_data>
      void always_inline unpack(class_with_member_data&a) const noexcept
      { store(a.data()); }
      /// store to aligned object with member data() returning element_type*
      template<typename class_with_member_data>
      void always_inline unpacku(class_with_member_data&a) const noexcept
      { storeu(a.data()); }
      /// store to aligned object with member data() returning element_type*
      template<bool aligned, typename class_with_member_data>
      void always_inline unpack_t(class_with_member_data&a) const noexcept
      { store_t<aligned>(a.data()); }
      /// store directly to aligned memory without polluting cashes
      void always_inline stream(element_type*p) const noexcept
      { _mm256_stream_pd(p,_m); }
      //@}

      /// \name vectorised arithmetic operations
      //@{
      /// +=
      packed& always_inline operator+=(packed p) noexcept
      { _m = _mm256_add_pd(p._m,_m); return*this; }
      /// +
      packed always_inline operator+ (packed p) const noexcept
      { return packed(_mm256_add_pd(_m,p._m)); }
      /// -=
      packed& always_inline operator-=(packed p) noexcept
      { _m = _mm256_sub_pd(_m,p._m); return*this; }
      /// -
      packed always_inline operator- (packed p) const noexcept
      { return packed(_mm256_sub_pd(_m,p._m)); }
      /// unary -
      packed always_inline operator-() const noexcept
      { return packed(_mm256_xor_pd(_m,sgn_mask())); }
      /// *=
      packed& always_inline operator*=(packed p) noexcept
      { _m = _mm256_mul_pd(p._m,_m); return*this; }
      /// *
      packed always_inline operator* (packed p) const noexcept
      { return packed(_mm256_mul_pd(_m,p._m)); }
      /// /=
      packed& always_inline operator/=(packed p) noexcept
      { _m = _mm256_div_pd(_m,p._m); return*this; }
      /// /
      packed always_inline operator/ (packed p) const noexcept
      { return packed(_mm256_div_pd(_m,p._m)); }
      /// sqrt
      friend packed always_inline sqrt(packed p) noexcept
      { return packed(_mm256_sqrt_pd(p._m)); }
      /// square
      friend packed always_inline square(packed p) noexcept
      { return  packed(_mm256_mul_pd(p._m,p._m)); }
      /// reciprocal
      friend packed always_inline reciprocal(packed p) noexcept
      { return packed(_mm256_div_pd(_mm256_set1_pd(1.0),p._m)); }
      /// maximum of two packed
      friend packed always_inline max(packed a, packed b) noexcept
      { return packed(_mm256_max_pd(a._m,b._m)); }
      /// minimum of two packed
      friend packed always_inline min(packed a, packed b) noexcept
      { return packed(_mm256_min_pd(a._m,b._m)); }
      /// abs
      friend packed always_inline abs(packed p) noexcept
      { return packed(_mm256_and_pd(p._m,abs_mask())); }
      /// -abs
      friend packed always_inline negabs(packed p) noexcept
      { return packed(_mm256_or_pd(p._m,sgn_mask())); }
      /// abs(x-y)
      friend packed always_inline diff(packed a, packed b) noexcept
      { return packed(_mm256_and_pd(_mm256_sub_pd(a._m,b._m),abs_mask())); }
      /// x = abs(x-y)
      packed& always_inline make_diff(packed p) noexcept
      {
	_m = _mm256_sub_pd(_m,p._m);
	_m = _mm256_and_pd(_m,abs_mask());
	return*this;
      }
      /// sign(a)*b
      friend packed always_inline signmove(packed a, packed b) noexcept
      { return packed(_mm256_or_pd(a.signmask(),
				   _mm256_and_pd(b._m,abs_mask()))); }
      //@}

      /// \name horizontal arithmetic operations
      //@{
      /// horizontal sum
      friend element_type always_inline sum(packed p) noexcept
      {
	aligned_element_block q;
	_mm256_store_pd(q,p._m);
	return q[0]+q[1]+q[2]+q[3];
      }
      /// sum of first three elements
      friend element_type always_inline sum3(packed p) noexcept
      {
	aligned_element_block q;
	_mm256_store_pd(q,p._m);
	return q[0]+q[1]+q[2];
      }
#if(0)
      /// horizontal sum in each element
      friend packed always_inline hadd(packed p) noexcept
      /// horizontal maximum
      friend element_type always_inline max(packed p) noexcept
      /// horizontal minimum
      friend element_type always_inline min(packed p) noexcept
#endif
      //@}

      /// \name vectorised comparisons
      //@{
      /// <
      packed always_inline operator< (packed p) const noexcept
      { return packed(_mm256_cmp_pd(_m,p._m,_CMP_LT_OQ)); }
      /// <=
      packed always_inline operator<=(packed p) const noexcept
      { return packed(_mm256_cmp_pd(_m,p._m,_CMP_LE_OQ)); }
      /// >
      packed always_inline operator> (packed p) const noexcept
      { return packed(_mm256_cmp_pd(_m,p._m,_CMP_GT_OQ)); }
      /// >=
      packed always_inline operator>=(packed p) const noexcept
      { return packed(_mm256_cmp_pd(_m,p._m,_CMP_GE_OQ)); }
      /// ==
      packed always_inline operator==(packed p) const noexcept
      { return packed(_mm256_cmp_pd(_m,p._m,_CMP_EQ_UQ)); }
      /// !=
      packed always_inline operator!=(packed p) const noexcept
      { return packed(_mm256_cmp_pd(_m,p._m,_CMP_NEQ_UQ)); }
      //@}

      /// \name vectorised logical operations
      /// &=
      packed& always_inline operator&=(packed p) noexcept
      { _m = _mm256_and_pd(p._m,_m); return*this; }
      /// &
      packed always_inline operator& (packed p) const noexcept
      { return packed(_mm256_and_pd(_m,p._m)); }
      /// |=
      packed& always_inline operator|=(packed p) noexcept
      { _m = _mm256_or_pd(p._m,_m); return*this; }
      /// |
      packed always_inline operator| (packed p) const noexcept
      { return packed(_mm256_or_pd(_m,p._m)); }
      /// ^=
      packed& always_inline operator^=(packed p) noexcept
      { _m = _mm256_xor_pd(p._m,_m); return*this; }
      /// ^
      packed always_inline operator^ (packed p) const noexcept
      { return packed(_mm256_xor_pd(_m,p._m)); }
      /// unary !
      packed always_inline operator! () const noexcept
      { return packed(_mm256_xor_pd(_m,one_mask())); }
      /// and not: this=this&!that
      packed always_inline andnot(packed p) const noexcept
      { return packed(_mm256_andnot_pd(p._m,_m)); }
      //@}

      /// \name boolean or integer properties
      //@{
      /// return integer with bits equal to sign bits
      /// \note if @a p is the result of a boolean operation (comparison, or
      ///       bit-wise operations on comparison results), signbit(p) can be
      ///       used to in an if statement.
      friend always_inline int signbits(packed p) noexcept
      { return _mm256_movemask_pd(p._m); }
      /// is any element negative
      friend always_inline bool has_negative(packed p) noexcept
      { return _mm256_movemask_pd(p._m); }
      /// are all elements negative
      friend always_inline bool all_negative(packed p) noexcept
      { return _mm256_movemask_pd(p._m) == max_signbits; }
      /// is any element non-zero
      friend always_inline bool has_non_zero(packed p) noexcept
      { return signbits(p!=packed::zero()); }
      /// are all elements non-zero
      friend always_inline bool all_non_zero(packed p) noexcept
      { return signbits(p!=packed::zero()) == max_signbits; }
      /// is any element zero
      friend always_inline bool has_zero(packed p) noexcept
      { return signbits(p==packed::zero()); }
      /// are all elements zero
      friend always_inline bool all_zero(packed p) noexcept
      { return signbits(p==packed::zero()) == max_signbits; }
      //@}

      /// \name miscellaneous
      //@{
      /// set all elements to Kth element of argument
      template<int K>
      friend packed always_inline single(packed p) noexcept;
      /// blend two vectors depending on sign of third:  result = sign<0? x : y
      friend packed blend(packed sign, packed x, packed y) noexcept
      { return packed(_mm256_blendv_pd(y._m,x._m,sign._m)); }
      /// combine two vectors depending on third:  result = mask? x : y
      friend packed combine(packed mask, packed x, packed y) noexcept
      { return packed(_mm256_blendv_pd(y._m,x._m,mask._m)); }
      //@}
      static element_type always_inline nil_mask_elem() noexcept
      { return WDutils::meta::double_and_uint(0x0lu).d; }
      static element_type always_inline all_mask_elem() noexcept
      { return WDutils::meta::double_and_uint(0xfffffffffffffffflu).d; }
      static element_type always_inline sgn_mask_elem() noexcept
      { return WDutils::meta::double_and_uint(0x8000000000000000lu).d; }
      static element_type always_inline abs_mask_elem() noexcept
      { return WDutils::meta::double_and_uint(0x7ffffffffffffffflu).d; }
    private:
      static data_type always_inline nil_mask() noexcept
      { return _mm256_castsi256_pd(_mm256_set1_epi64x(0x0l)); }
      static data_type always_inline one_mask() noexcept
      { return _mm256_castsi256_pd(_mm256_set1_epi64x(int64_t(0xffffffffffffffff))); }
      static data_type always_inline sgn_mask() noexcept
      { return _mm256_castsi256_pd(_mm256_set1_epi64x(int64_t(0x8000000000000000))); }
      static data_type always_inline abs_mask() noexcept
      { return _mm256_castsi256_pd(_mm256_set1_epi64x(0x7fffffffffffffff)); }
      /// just our sign bits
      data_type always_inline signmask() const noexcept
      { return _mm256_and_pd(_m,sgn_mask()); }
      /// data
      data_type _m;
    };// SSE::packed<4,double>
    //
    typedef packed<4,double> dvec4;
    //
    template<unsigned I> inline
    dvec2 always_inline extract(dvec4 p) noexcept
    {
      static_assert(I<2,"index out of range");
      return dvec2(_mm256_extractf128_pd(p._m,I));
    }
    //
    template<unsigned I> inline
    double always_inline at(dvec4 p) noexcept
    {
      static_assert(I<dvec4::block_size,"index out of range");
#  if   defined(__clang__)
      return p._m[I];
#  else
      return SSE::at<(I&1)>(extract<(I>>1)>(p));
#  endif
    }
    //
    template<int K> inline
    dvec4 always_inline single(dvec4 p) noexcept
    {
      static_assert(K>=0 && K<dvec4::block_size,"K out of range");
#  ifdef __AVX2__
#   warning more efficient implementation possible with AVX2
#  endif
      return dvec4(_mm256_permute_pd
		   (_mm256_permute2f128_pd(p._m,p._m,K&2?49:32),K&1?15:0));
    }
  }
  //
  WDutils_TRAITS(SSE::dvec4,"dvec4");
  //
  namespace SSE {
# endif // __AVX__
#endif  // __SSE__
    /// is_packed<X>::value is true if X is SSE::packed<K,T>
    template<typename vec>
    struct is_packed {
      static const bool value = false
#ifdef __SSE__
	|| is_same<fvec4,vec>::value
#endif
#ifdef __SSE2__
	|| is_same<dvec2,vec>::value
#endif
#ifdef __AVX__
	|| is_same<dvec4,vec>::value
	|| is_same<fvec8,vec>::value
#endif
	;
    };
    /// get_packed<X>::type is packed<K,X> with largest supported K
    template<typename T> struct get_packed;
#ifdef __SSE__
    template<> struct get_packed<float>
# ifdef __AVX__
    { typedef fvec8 type; };
# else
    { typedef fvec4 type; };
# endif
#endif
#ifdef __SSE2__
    template<> struct get_packed<double>
# ifdef __AVX__
    { typedef dvec4 type; };
# else
    { typedef dvec2 type; };
# endif
#endif
  } // namespace WDutils::SSE
  //
#ifdef __SSE__
# ifdef __INTEL_COMPILER
// #  ifndef __SSE4_1__
  inline float xmm0(__m128 _A)
  { float t; _MM_EXTRACT_FLOAT(t,_A,0); return t; }
  inline float xmm1(__m128 _A)
  { float t; _MM_EXTRACT_FLOAT(t,_A,1); return t; }
  inline float xmm2(__m128 _A)
  { float t; _MM_EXTRACT_FLOAT(t,_A,2); return t; }
  inline float xmm3(__m128 _A)
  { float t; _MM_EXTRACT_FLOAT(t,_A,3); return t; }
// #  endif// __SSE4_1__
#  ifdef __SSE2__
  inline double xmm0(__m128d _A)
  { double x; _mm_storel_pd(&x,_A); return x; }
  inline double xmm1(__m128d _A)
  { double x; _mm_storeh_pd(&x,_A); return x; }
#  endif// __SSE2__
# elif defined(__clang__)  // clang
  inline float xmm0(__m128 _A)
  { return _A[0]; }
  inline float xmm1(__m128 _A)
  { return _A[1]; }
  inline float xmm2(__m128 _A)
  { return _A[2]; }
  inline float xmm3(__m128 _A)
  { return _A[3]; }
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
# ifdef __AVX__
#  if defined(__clang__)  // clang
  inline float ymm0(__m256 _A)
  { return _A[0]; }
  inline float ymm1(__m256 _A)
  { return _A[1]; }
  inline float ymm2(__m256 _A)
  { return _A[2]; }
  inline float ymm3(__m256 _A)
  { return _A[3]; }
  inline float ymm4(__m256 _A)
  { return _A[4]; }
  inline float ymm5(__m256 _A)
  { return _A[5]; }
  inline float ymm6(__m256 _A)
  { return _A[6]; }
  inline float ymm7(__m256 _A)
  { return _A[7]; }
  //
  inline double ymm0(__m256d _A)
  { return _A[0]; }
  inline double ymm1(__m256d _A)
  { return _A[1]; }
  inline double ymm2(__m256d _A)
  { return _A[2]; }
  inline double ymm3(__m256d _A)
  { return _A[3]; }
#  else
  inline float ymm0(__m256 _A)
  { return xmm0(_mm256_extractf128_ps(_A,0)); }
  inline float ymm1(__m256 _A)
  { return xmm1(_mm256_extractf128_ps(_A,0)); }
  inline float ymm2(__m256 _A)
  { return xmm2(_mm256_extractf128_ps(_A,0)); }
  inline float ymm3(__m256 _A)
  { return xmm3(_mm256_extractf128_ps(_A,0)); }
  inline float ymm4(__m256 _A)
  { return xmm0(_mm256_extractf128_ps(_A,1)); }
  inline float ymm5(__m256 _A)
  { return xmm1(_mm256_extractf128_ps(_A,1)); }
  inline float ymm6(__m256 _A)
  { return xmm2(_mm256_extractf128_ps(_A,1)); }
  inline float ymm7(__m256 _A)
  { return xmm3(_mm256_extractf128_ps(_A,1)); }
  //
  inline double ymm0(__m256d _A)
  { return xmm0(_mm256_extractf128_pd(_A,0)); }
  inline double ymm1(__m256d _A)
  { return xmm1(_mm256_extractf128_pd(_A,0)); }
  inline double ymm2(__m256d _A)
  { return xmm0(_mm256_extractf128_pd(_A,1)); }
  inline double ymm3(__m256d _A)
  { return xmm1(_mm256_extractf128_pd(_A,1)); }
#  endif
# endif // __AVX__
#endif // __SSE__
//
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
      /// largest packed type
#if   defined(__AVX__)
      typedef packed<8,float> packed_type;
#elif defined(__SSE__)
      typedef packed<4,float> packed_type;
#endif
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
      /// largest packed type
#if   defined(__AVX__)
      typedef packed<4,double> packed_type;
#elif defined(__SSE__)
      typedef packed<2,double> packed_type;
#endif
    };
    ////////////////////////////////////////////////////////////////////////////
    namespace details {
      template<typename Functor>
      struct is_functor {
	static const bool value =
	  is_same<Functor,meta::assign    >::value ||
	  is_same<Functor,meta::add       >::value ||
	  is_same<Functor,meta::subtract  >::value ||
	  is_same<Functor,meta::multiply  >::value ||
	  is_same<Functor,meta::divide    >::value ||
	  is_same<Functor,meta::swap      >::value ||
	  is_same<Functor,meta::sqrt      >::value ||
	  is_same<Functor,meta::square    >::value ||
	  is_same<Functor,meta::reciprocal>::value ||
	  is_same<Functor,meta::maximum   >::value ||
	  is_same<Functor,meta::minimum   >::value;
      };
      // auxiliary for connect<>
      template<typename RealType>
      struct connector
      {
	// for(i=0; i!=N; ++i) Functor::operate(x[i], y[i]);
	template<typename Functor>
	static void connect(RealType*x, const RealType*y, unsigned n) noexcept
	{ for(unsigned i=0; i!=n; ++i) Functor::operate(x[i],y[i]); }
	// for(i=0; i!=N; ++i) Functor::operate(x[i], y[i]);
	static void swap(RealType*x, RealType*y, unsigned n) noexcept
	{ for(unsigned i=0; i!=n; ++i) meta::swap::operate(x[i],y[i]); }
	// for(i=0; i!=N; ++i) Functor::operate(x[i], y);
	template<typename Functor>
	static void foreach(RealType*x, const RealType y, unsigned n) noexcept
	{ for(unsigned i=0; i!=n; ++i) Functor::operate(x[i],y); }
      };
#ifdef __SSE__
      // auxiliary for struct connector<>
      template<typename real>
      struct connector_helper
      {
	//
	template<typename Functor>
	static always_inline void connect(real*x, real y) noexcept
	{ Functor::operate(*x,y); }
	//
	template<typename Functor, typename vec>
	static always_inline
	typename enable_if<!is_same<Functor,meta::assign>::value &&
	                   !is_same<Functor,meta::sqrt>::value &&
	                   !is_same<Functor,meta::square>::value &&
	                   !is_same<Functor,meta::reciprocal>::value &&
	                    is_packed<vec>::value>::type
	connect(real*x, vec const&vy) noexcept
	{
	  vec vx = vec::load(x);
	  Functor::operate(vx,vy);
	  vx.store(x);
	}
	//
	template<typename Functor, typename vec>
	static always_inline
	typename enable_if<is_same<Functor,meta::assign>::value>::type
	connect(real*x, vec const&vy) noexcept
	{ vy.store(x); }
	//
	template<typename Functor, typename vec>
	static always_inline
	typename enable_if<is_same<Functor,meta::sqrt>::value>::type
	connect(real*x, vec const&vy) noexcept
	{ sqrt(vy).store(x); }
	//
	template<typename Functor, typename vec>
	static always_inline
	typename enable_if<is_same<Functor,meta::square>::value>::type
	connect(real*x, vec const&vy) noexcept
	{ (vy*vy).store(x); }
	//
	template<typename Functor, typename vec>
	static always_inline
	typename enable_if<is_same<Functor,meta::reciprocal>::value>::type
	connect(real*x, vec const&vy) noexcept
	{
	  static vec one = vec::one();
	  (one/vy).store(x);
	}
	//
	template<unsigned block_size, bool y_aligned>
	static always_inline
	typename enable_if<block_size==1>::type
	swap(real*x, real*y) noexcept
	{ real tmp(*x); *x=*y; *y=tmp; }
	//
	template<unsigned block_size, bool y_aligned>
	static always_inline
	typename enable_if<block_size!=1>::type
	swap(real*x, real*y) noexcept
	{
	  typedef packed<block_size,real> vec;
	  vec vx = vec::load(x);
	  vec vy = vec::template load_t<y_aligned>(y);
	  vy.store(x);
	  vx.template store_t<y_aligned>(y);
	}
      };
#endif
#ifdef __SSE2__
      // templated array connection: SSE version for @c RealType = @c double
      template<>
      struct connector<double> : private connector_helper<double>
      {
	typedef connector_helper<double> helper;
	// for(i=0; i!=N; ++i) Functor::operate(x[i], y[i]);
	template<typename Functor>
	static void connect(double*x, const double*y, unsigned n) noexcept
	{
	  if(n && (size_t(x)&15)) {
	    helper::connect<Functor>(x,*y);
	    --n,++x,++y; 
	  }
# ifdef __AVX__
	  if(n>=2 && (size_t(x)&31)) {
	    helper::connect<Functor>(x,dvec2::loadu(y));
	    n-=2,x+=2,y+=2;
	  }
	  if(size_t(y)&31) {
	    for(; n>=4; n-=4,x+=4,y+=4)
	      helper::connect<Functor>(x,dvec4::loadu(y));
	    if(n>=2) {
	      helper::connect<Functor>(x,dvec2::loadu(y));
	      n-=2,x+=2,y+=2;
	    }
	  } else {
	    for(; n>=4; n-=4,x+=4,y+=4)
	      helper::connect<Functor>(x,dvec4::load(y));
	    if(n>=2) {
	      helper::connect<Functor>(x,dvec2::load(y));
	      n-=2,x+=2,y+=2;
	    }
	  }
# else
	  if(size_t(y)&15)
	    for(; n>=2; n-=2,x+=2,y+=2)
	      helper::connect<Functor>(x,dvec2::loadu(y));
	  else
	    for(; n>=2; n-=2,x+=2,y+=2)
	      helper::connect<Functor>(x,dvec2::load(y));
# endif
	  if(n)
	    helper::connect<Functor>(x,*y);
	}
	// for(i=0; i!=N; ++i) Functor::operate(x[i], y[i]);
	static void swap(double*x, double*y, unsigned n) noexcept
	{
	  if(n && (size_t(x)&15)) {
	    helper::swap<1,0>(x,y);
	    --n,++x,++y; 
	  }
# ifdef __AVX__
	  if(n>=2 && (size_t(x)&31)) {
	    helper::swap<2,0>(x,y);
	    n-=2,x+=2,y+=2;
	  }
	  if(size_t(y)&31) {
	    for(; n>=4; n-=4,x+=4,y+=4)
	      helper::swap<4,0>(x,y);
	    if(n>=2) {
	      helper::swap<2,0>(x,y);
	      n-=2,x+=2,y+=2;
	    }
	  } else {
	    for(; n>=4; n-=4,x+=4,y+=4)
	      helper::swap<4,1>(x,y);
	    if(n>=2) {
	      helper::swap<2,1>(x,y);
	      n-=2,x+=2,y+=2;
	    }
	  }
# else
	  if(size_t(y)&15)
	    for(; n>=2; n-=2,x+=2,y+=2)
	      helper::swap<2,0>(x,y);
	  else
	    for(; n>=2; n-=2,x+=2,y+=2)
	      helper::swap<2,1>(x,y);
# endif
	  if(n)
	    helper::swap<1,0>(x,y);
	}
	// for(i=0; i!=N; ++i) Functor::operate(x[i], y);
	template<typename Functor>
	static void foreach(double*x, const double y, unsigned n) noexcept
	{
	  if(n && (size_t(x)&15)) {
	    helper::connect<Functor>(x,y);
	    --n,++x; 
	  }
	  dvec2 y2(y);
# ifdef __AVX__
	  if(n>=2 && (size_t(x)&31)) {
	    helper::connect<Functor>(x,y2);
	    n-=2,x+=2;
	  }
	  dvec4 y4(y);
	  for(; n>=4; n-=4,x+=4)
	    helper::connect<Functor>(x,y4);
# endif
	  for(; n>=2; n-=2,x+=2)
	    helper::connect<Functor>(x,y2);
	  if(n)
	    helper::connect<Functor>(x,y);
	}
      };// struct SSE::details::connector<double>
#endif// __SSE2__
#ifdef __SSE__
      // templated array connection: SSE version for @c RealType = @c float
      template<>
      struct connector<float> : private connector_helper<float>
      {
	typedef connector_helper<float> helper;
	// for(i=0; i!=N; ++i) assignment_functor<float>(x[i], y[i]);
	template<typename Functor>
	static void connect(float*x, const float*y, unsigned n) noexcept
	{
	  for(; n && (size_t(x)&15); --n,++x,++y)
	    helper::connect<Functor>(x,*y);
# ifdef __AVX__
	  if(n>=4 && (size_t(x)&31)) {
	    helper::connect<Functor>(x,fvec4::loadu(y));
	    n-=4,x+=4,y+=4;
	  }
	  if(size_t(y)&31) {
	    for(; n>=8; n-=8,x+=8,y+=8)
	      helper::connect<Functor>(x,fvec8::loadu(y));
	    if(n>=4) {
	      helper::connect<Functor>(x,fvec4::loadu(y));
	      n-=4,x+=4,y+=4;
	    }
	  } else {
	    for(; n>=8; n-=8,x+=8,y+=8)
	      helper::connect<Functor>(x,fvec8::load(y));
	    if(n>=4) {
	      helper::connect<Functor>(x,fvec4::load(y));
	      n-=4,x+=4,y+=4;
	    }
	  }
# else
	  if(size_t(y)&15)
	    for(; n>=4; n-=4,x+=4,y+=4)
	      helper::connect<Functor>(x,fvec4::loadu(y));
	  else
	    for(; n>=4; n-=4,x+=4,y+=4)
	      helper::connect<Functor>(x,fvec4::load(y));
# endif
	  for(; n; --n,++x,++y)
	    helper::connect<Functor>(x,*y);
	}
	// for(i=0; i!=N; ++i) assignment_functor<float>(x[i], y[i]);
	static void swap(float*x, float*y, unsigned n) noexcept
	{
	  for(; n && (size_t(x)&15); --n,++x,++y)
	    helper::swap<1,0>(x,y);
# ifdef __AVX__
	  if(n>=4 && (size_t(x)&31)) {
	    helper::swap<4,0>(x,y);
	    n-=4,x+=4,y+=4;
	  }
	  if(size_t(y)&31) {
	    for(; n>=8; n-=8,x+=8,y+=8)
	      helper::swap<8,0>(x,y);
	    if(n>=4) {
	      helper::swap<4,0>(x,y);
	      n-=4,x+=4,y+=4;
	    }
	  } else {
	    for(; n>=8; n-=8,x+=8,y+=8)
	      helper::swap<8,1>(x,y);
	    if(n>=4) {
	      helper::swap<4,1>(x,y);
	      n-=4,x+=4,y+=4;
	    }
	  }
# else
	  if(size_t(y)&15)
	    for(; n>=4; n-=4,x+=4,y+=4)
	      helper::swap<4,0>(x,y);
	  else
	    for(; n>=4; n-=4,x+=4,y+=4)
	      helper::swap<4,1>(x,y);
# endif
	  for(; n; --n,++x,++y)
	    helper::swap<1,0>(x,y);
	}
	// for(i=0; i!=N; ++i) assignment_functor<float>(x[i], y);
	template<typename Functor>
	static void foreach(float*x, const float y, unsigned n) noexcept
	{
	  for(; n && (size_t(x)&15); --n,++x)
	    helper::connect<Functor>(x,y);
	  fvec4 y4(y);
# ifdef __AVX__
	  if(n>=4 && (size_t(x)&31)) {
	    helper::connect<Functor>(x,y4);
	    n-=4,x+=4;
	  }
	  fvec8 y8(y);
	  for(; n>=8; n-=8,x+=8)
	    helper::connect<Functor>(x,y8);
# endif
	  for(; n>=4; n-=4,x+=4)
	    helper::connect<Functor>(x,y4);
	  for(; n; --n,++x)
	    helper::connect<Functor>(x,y);
	}
      };// struct SSE::details::connector<float>
#endif// __SSE__
    } // namespace WDutils::SSE::details
    ///
    /// connect two arrays element wise:
    /// @code  for(i=0; i!=n; ++i) Functor::operate(x[i], y[i]); @endcode
    ///
    template<typename Functor, typename RealType>
    typename enable_if< details::is_functor<Functor>::value &&
		       !is_same<Functor,meta::swap >::value>::type
    connect(RealType*x, const RealType*y, unsigned n) noexcept
    { details::connector<RealType>::template connect<Functor>(x,y,n); }
    ///
    /// swap elements of two arrays:
    /// @code  for(i=0; i!=n; ++i) swap(x[i], y[i]); @endcode
    ///
    template<typename RealType>
    void swap(RealType*x, RealType*y, unsigned n) noexcept
    { details::connector<RealType>::swap(x,y,n); }
    //
    template<typename Functor, typename RealType>
    typename enable_if<is_same<Functor,meta::swap>::value>::type
    connect(RealType*x, RealType*y, unsigned n) noexcept
    { swap(x,y,n); }
    ///
    /// connect each array element with (the same) scalar:
    /// @code  for(i=0; i!=n; ++i) Functor::operate(x[i], y); @endcode
    ///
    template<typename Functor, typename RealType>
    typename enable_if< details::is_functor<Functor> ::value &&
                       !is_same<Functor,meta::swap>  ::value &&
                       !is_same<Functor,meta::divide>::value>::type
    foreach(RealType*x, const RealType y, unsigned n) noexcept
    { details::connector<RealType>::template foreach<Functor>(x,y,n); }
    //
    template<typename Functor, typename RealType>
    typename enable_if<is_same<Functor,meta::divide>::value>::type
    foreach(RealType*x, const RealType y, unsigned n) noexcept
    { foreach<meta::multiply>(x,RealType(1)/y,n); }
    //
    namespace details {
      // static templated array connection: non-SSE version
      template<typename RealType>
      struct static_connector
      {
	// for(i=0; i!=N; ++i) Functor::operate(x[i], y[i]);
	template<unsigned N, typename Functor>
	static typename enable_if<N>::type
	connect(RealType*x, const RealType*y) noexcept
	{
	  Functor::operate(*x,*y);
	  connect<N-1,Functor>(++x,++y);
	}
	template<unsigned N, typename Functor>
	static typename enable_if<N==0>::type
	always_inline connect(RealType*, const RealType*) noexcept {}
	// for(i=0; i!=N; ++i) swap(x[i], y[i]);
	template<unsigned N, typename Functor>
	static typename enable_if<N>::type
	swap(RealType*x, RealType*y) noexcept
	{
	  meta::swap::operate(*x,*y);
	  swap<N-1,Functor>(++x,++y);
	}
	template<unsigned N, typename Functor>
	static typename enable_if<N==0>::type
	always_inline swap(RealType*, RealType*) noexcept {}
	// for(i=0; i!=N; ++i) Functor::operate(x[i], y);
	template<unsigned N, typename Functor>
	static typename enable_if<N>::type
	foreach(RealType*x, const RealType y) noexcept
	{
	  Functor::operate(*x,y);
	  foreach<N-1,Functor>(++x,y);
	}
	template<unsigned N, typename Functor>
	static typename enable_if<N==0>::type
	always_inline foreach(RealType*, RealType) noexcept {}
      };
#ifdef __SSE2__
      // static templated array connection: SSE version for @c double
      template<>
      class static_connector<double> : connector_helper<double>
      {
	typedef connector_helper<double> helper;
# if __cplusplus >= 201103L
	//
	static constexpr unsigned half_block(unsigned block)
	{ return block>1? block/2 : 1; }
	//
	static constexpr unsigned block_size(unsigned N, unsigned max_block)
	{
	  return
	    N >= max_block? max_block : 
	    N >= half_block(max_block)? half_block(max_block) : 1;
	}
# else  //
	template<unsigned N, unsigned maxb>
	struct size_block {
	  static const unsigned half  = maxb>1? maxb/2 : 1;
	  static const unsigned block = N >= maxb? maxb : 
	                                N >= half? half : 1;
	};
#  define block_size(NUM,MAXBLOCK) size_block<NUM,MAXBLOCK>::block
# endif
	//
	template<typename Functor, unsigned N, unsigned max_block>
	static void always_inline connect_start(double*x, const double*y)
	  noexcept
	{
	  static const unsigned block = block_size(N,max_block);
	  static const size_t   align = sizeof(double)*block-1;
	  if(size_t(y)&align) connect_t<Functor,N,block,0>(x,y);
	  else                connect_t<Functor,N,block,1>(x,y);
	}
	//
	template<typename Functor, unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<N==0 && block==1>::type
	connect_t(double*, const double*) noexcept {}
	//
	template<typename Functor, unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<N==1 && block==1>::type
	connect_t(double*x, const double*y) noexcept
	{ helper::connect<Functor>(x,*y); }
	//
	template<typename Functor, unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<(N>1) && block==1>::type
	connect_t(double*x, const double*y) noexcept
	{
	  static_assert(N>1,"block size mismatch");
	  if(size_t(x)&15) {
	    helper::connect<Functor>(x,*y);
	    connect_start<Functor,N-1,2>(++x,++y);
	  } else
	    connect_start<Functor,N  ,2>(x,y);
	}
	//
	template<typename Functor, unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<block==2
# ifdef __AVX__
				  && (N<4)
# endif
	                         >::type
	connect_t(double*x, const double*y) noexcept
	{
	  static_assert(N>=block,"block size mismatch");
	  helper::connect<Functor>(x,dvec2::load_t<y_aligned>(y));
	  connect_t<Functor,N-2,block_size(N-2,2),y_aligned>(x+=2,y+=2);
	}
# ifdef __AVX__
	//
	template<typename Functor, unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<(N>=4) && block==2>::type
	connect_t(double*x, const double*y) noexcept
	{
	  static_assert(N>=block,"block size mismatch");
	  if(size_t(x)&31) {
	    helper::connect<Functor>(x,dvec2::load_t<y_aligned>(y));
	    connect_start<Functor,N-2,4>(x+=2,y+=2);
	  } else
	    connect_start<Functor,N  ,4>(x,y);
	}
	//
	template<typename Functor, unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<block==4>::type
	connect_t(double*x, const double*y) noexcept
	{
	  static_assert(N>=block,"block size mismatch");
	  helper::connect<Functor>(x,dvec4::load_t<y_aligned>(y));
	  connect_t<Functor,N-4,block_size(N-4,4),y_aligned>(x+=4,y+=4);
	}
# endif
      public:
	// for(i=0; i!=N; ++i) Functor::operator(x[i], y[i]);
	template<unsigned N, typename Functor>
	static void connect(double*x, const double*y) noexcept
	{ connect_t<Functor,N,1,0>(x,y); }

      private:
	//
	template<unsigned N, unsigned max_block>
	static void always_inline swap_start(double*x, double*y)
	  noexcept
	{
	  static const unsigned block = block_size(N,max_block);
	  static const size_t   align = sizeof(double)*block-1;
	  if(size_t(y)&align) swap_t<N,block,0>(x,y);
	  else                swap_t<N,block,1>(x,y);
	}
	//
	template<unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<N==0 && block==1>::type
	swap_t(double*, double*) noexcept {}
	//
	template<unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<N==1 && block==1>::type
	swap_t(double*x, double*y) noexcept
	{ helper::swap<1,0>(x,y); }
	//
	template<unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<(N>1) && block==1>::type
	swap_t(double*x, double*y) noexcept
	{
	  static_assert(N>1,"block size mismatch");
	  if(size_t(x)&15) {
	    helper::swap<1,0>(x,y);
	    swap_start<N-1,2>(++x,++y);
	  } else
	    swap_start<N  ,2>(x,y);
	}
	//
	template<unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<block==2
# ifdef __AVX__
				  && (N<4)
# endif
	                         >::type
	swap_t(double*x, double*y) noexcept
	{
	  static_assert(N>=block,"block size mismatch");
	  helper::swap<2,y_aligned>(x,y);
	  swap_t<N-2,block_size(N-2,2),y_aligned>(x+=2,y+=2);
	}
# ifdef __AVX__
	//
	template<unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<(N>=4) && block==2>::type
	swap_t(double*x, double*y) noexcept
	{
	  static_assert(N>=block,"block size mismatch");
	  if(size_t(x)&31) {
	    helper::swap<2,y_aligned>(x,y);
	    swap_start<N-2,4>(x+=2,y+=2);
	  } else
	    swap_start<N  ,4>(x,y);
	}
	//
	template<unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<block==4>::type
	swap_t(double*x, double*y) noexcept
	{
	  static_assert(N>=block,"block size mismatch");
	  helper::swap<4,y_aligned>(x,y);
	  swap_t<N-4,block_size(N-4,4),y_aligned>(x+=4,y+=4);
	}
# endif
      public:
	// for(i=0; i!=N; ++i) swap(x[i], y[i]);
	template<unsigned N>
	static void swap(double*x, double*y) noexcept
	{ swap_t<N,1,0>(x,y); }
	//
      private:
# ifdef __AVX__
#  define VECS_DECL dvec2, dvec4
# else
#  define VECS_DECL dvec2
# endif
	//
	template<typename Functor, unsigned N, unsigned block>
	static typename enable_if<N==0 && block==1>::type
	foreach_t(double*, double, VECS_DECL) noexcept {}
	//
	template<typename Functor, unsigned N, unsigned block>
	static typename enable_if<N==1 && block==1>::type
	foreach_t(double*x, const double y, VECS_DECL) noexcept
	{ helper::connect<Functor>(x,y); }
	//
#  undef  VECS_DECL
# ifdef __AVX__
#  define VECS_DECL const double y, const dvec2&y2, const dvec4&y4
#  define VECS_PASS y, y2, y4
# else
#  define VECS_DECL const double y, const dvec2&y2
#  define VECS_PASS y, y2
# endif
	//
	template<typename Functor, unsigned N, unsigned block>
	static typename enable_if<(N>1) && block==1>::type
	foreach_t(double*x, VECS_DECL) noexcept
	{
	  static_assert(N>1,"block size mismatch");
	  if(size_t(x)&15) {
	    helper::connect<Functor>(x,y);
	    foreach_t<Functor,N-1,block_size(N-1,2)>(++x, VECS_PASS);
	  } else
	    foreach_t<Functor,N,block_size(N,2)>(x, VECS_PASS);
	}
	//
	template<typename Functor, unsigned N, unsigned block>
	static typename enable_if<block==2
# ifdef __AVX__
				  && (N<4)
# endif
			       	   >::type
	foreach_t(double*x, VECS_DECL) noexcept
	{
	  static_assert(N>=block,"block size mismatch");
	  helper::connect<Functor>(x,y2);
	  foreach_t<Functor,N-2,block_size(N-2,2)>(x+=2, VECS_PASS);
	}
# ifdef __AVX__
	//
	template<typename Functor, unsigned N, unsigned block>
	static typename enable_if<(N>=4) && block==2>::type
	foreach_t(double*x, VECS_DECL) noexcept
	{
	  static_assert(N>=block,"block size mismatch");
	  if(size_t(x)&31) {
	    helper::connect<Functor>(x,y2);
	    foreach_t<Functor,N-2,block_size(N-2,4)>(x+=2, VECS_PASS);
	  } else
	    foreach_t<Functor,N,block_size(N,4)>(x, VECS_PASS);
	}
	//
	template<typename Functor, unsigned N, unsigned block>
	static typename enable_if<block==4>::type
	foreach_t(double*x, VECS_DECL) noexcept
	{
	  static_assert(N>=block,"block size mismatch");
	  helper::connect<Functor>(x,y4);
	  foreach_t<Functor,N-4,block_size(N-4,4)>(x+=4, VECS_PASS);
	}
# endif
# undef block_size
# undef VECS_DECL
# undef VECS_PASS
	//
	template<typename Functor, unsigned N>
	static typename enable_if<N==0>::type
	m_foreach(double*, double) noexcept {}
	//
	template<typename Functor, unsigned N>
	static typename enable_if<N==1>::type
	m_foreach(double*x, const double y) noexcept
	{ helper::connect<Functor>(x,y); }
# ifdef __AVX__
	//
	template<typename Functor, unsigned N>
	static typename enable_if<(N>=2 && N<4)>::type
	m_foreach(double*x, double y) noexcept
	{ foreach_t<Functor,N,1>(x,y,dvec2(y),dvec4()); }
	//
	template<typename Functor, unsigned N>
	static typename enable_if<(N>=4)>::type
	m_foreach(double*x, double y) noexcept
	{ foreach_t<Functor,N,1>(x,y,dvec2(y),dvec4(y)); }
#else   //
	template<typename Functor, unsigned N>
	static typename enable_if<(N>=2)>::type
	m_foreach(double*x, double y) noexcept
	{ foreach_t<Functor,N,1>(x,y,dvec2(y)); }
# endif
      public:
	// for(i=0; i!=N; ++i) Functor::operator(x[i], y);
	template<unsigned N, typename Functor>
	static void foreach(double*x, double y) noexcept
	{ m_foreach<Functor,N>(x,y); }
      };// struct SSE::details::static_connector<double>
#endif // __SSE2__
#ifdef __SSE__
      // static templated array connection: SSE version for @c float
      template<>
      struct static_connector<float> : private connector_helper<float>
      {
	typedef connector_helper<float> helper;
# if __cplusplus >= 201103L && !defined(__INTEL_COMPILER)
	//
	static constexpr unsigned smll_block(unsigned block)
	{ return block==8? 4 : 1; }
	//
	static constexpr unsigned block_size(unsigned N, unsigned max_block)
	{
	  return
	    N >= max_block? max_block  : 
	    N >= smll_block(max_block)? smll_block(max_block) : 1;
	}
# else  //
	template<unsigned N, unsigned maxb>
	struct size_block {
	  static const unsigned smll  = maxb==8? 4 : 1;
	  static const unsigned block = N >= maxb? maxb : 
	                                N >= smll? smll : 1;
	};
#  define block_size(NUM,MAXBLOCK) size_block<NUM,MAXBLOCK>::block
# endif
	//
	template<typename Functor, unsigned N, unsigned max_block>
	static void always_inline connect_start(float*x, const float*y) noexcept
	{
	  static const unsigned block = block_size(N,max_block);
	  static const size_t   align = sizeof(float)*block-1;
	  if(size_t(y)&align) connect_t<Functor,N,block,0>(x,y);
	  else                connect_t<Functor,N,block,1>(x,y);
	}
	//
	template<typename Functor, unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<N==0 && block==1>::type
	connect_t(float*, const float*) noexcept {}
	//
	template<typename Functor, unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<N==1 && block==1>::type
	connect_t(float*x, const float*y) noexcept
	{ helper::connect<Functor>(x,*y); }
	//
	template<typename Functor, unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<N==2 && block==1>::type
	connect_t(float*x, const float*y) noexcept
	{
	  helper::connect<Functor>(x++,*y++);
	  helper::connect<Functor>(x  ,*y  );
	}
	//
	template<typename Functor, unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<N==3 && block==1>::type
	connect_t(float*x, const float*y) noexcept
	{
	  helper::connect<Functor>(x++,*y++);
	  helper::connect<Functor>(x++,*y++);
	  helper::connect<Functor>(x  ,*y  );
	}
	//
	template<typename Functor, unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<(N>3) && block==1>::type
	connect_t(float*x, const float*y) noexcept
	{
	  static_assert(N>3,"block size mismatch");
	  switch((size_t(x)&15)) {
	  case 12:
	    helper::connect<Functor>(x++,*y++);
	    connect_start<Functor,N-1,4>(x,y);
	    break;
	  case 8:
	    helper::connect<Functor>(x++,*y++);
	    helper::connect<Functor>(x++,*y++);
	    connect_start<Functor,N-2,4>(x,y);
	    break;
	  case 4:
	    helper::connect<Functor>(x++,*y++);
	    helper::connect<Functor>(x++,*y++);
	    helper::connect<Functor>(x++,*y++);
	    connect_start<Functor,N-3,4>(x,y);
	    break;
	  default:
	    connect_start<Functor,N  ,4>(x,y);
	  }
	}
	//
	template<typename Functor, unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<block==4
# ifdef __AVX__
				  && (N<8)
# endif
			       	   >::type
	connect_t(float*x, const float*y) noexcept
	{
	  static_assert(N>=block,"block size mismatch");
	  helper::connect<Functor>(x,fvec4::load_t<y_aligned>(y));
	  connect_t<Functor,N-4,block_size(N-4,4),y_aligned>(x+=4,y+=4);
	}
# ifdef __AVX__
	//
	template<typename Functor, unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<(N>=8) && block==4>::type
	connect_t(float*x, const float*y) noexcept
	{
	  static_assert(N>=block,"block size mismatch");
	  if(size_t(x)&31) {
	    helper::connect<Functor>(x,fvec4::load_t<y_aligned>(y));
	    connect_start<Functor,N-4,8>(x+=4,y+=4);
	  } else
	    connect_start<Functor,N  ,8>(x,y);
	}
	//
	template<typename Functor, unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<block==8>::type
	connect_t(float*x, const float*y) noexcept
	{
	  static_assert(N>=block,"block size mismatch");
	  helper::connect<Functor>(x,fvec8::load_t<y_aligned>(y));
	  connect_t<Functor,N-8,block_size(N-8,8),y_aligned>(x+=8,y+=8);
	}
# endif
      public:
	/// unrolled for(i=0; i!=N; ++i) Functor::operator(x[i], y[i]);
	template<unsigned N, typename Functor>
	static void connect(float*x, const float*y) noexcept
	{ connect_t<Functor,N,1,0>(x,y); }
      private:
	//
	template<unsigned N, unsigned max_block>
	static void always_inline swap_start(float*x, float*y) noexcept
	{
	  static const unsigned block = block_size(N,max_block);
	  static const size_t   align = sizeof(float)*block-1;
	  if(size_t(y)&align) swap_t<N,block,0>(x,y);
	  else                swap_t<N,block,1>(x,y);
	}
	//
	template<unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<N==0 && block==1>::type
	swap_t(float*, float*) noexcept {}
	//
	template<unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<N==1 && block==1>::type
	swap_t(float*x, float*y) noexcept
	{ helper::swap<1,0>(x,y); }
	//
	template<unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<N==2 && block==1>::type
	swap_t(float*x, float*y) noexcept
	{
	  helper::swap<1,0>(x++,y++);
	  helper::swap<1,0>(x  ,y  );
	}
	//
	template<unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<N==3 && block==1>::type
	swap_t(float*x, float*y) noexcept
	{
	  helper::swap<1,0>(x++,y++);
	  helper::swap<1,0>(x++,y++);
	  helper::swap<1,0>(x  ,y  );
	}
	//
	template<unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<(N>3) && block==1>::type
	swap_t(float*x, float*y) noexcept
	{
	  static_assert(N>=4,"block size mismatch");
	  switch((size_t(x)&15)) {
	  case 12:
	    helper::swap<1,0>(x++,y++);
	    swap_start<N-1,4>(x,y);
	    break;
	  case 8:
	    helper::swap<1,0>(x++,y++);
	    helper::swap<1,0>(x++,y++);
	    swap_start<N-2,4>(x,y);
	    break;
	  case 4:
	    helper::swap<1,0>(x++,y++);
	    helper::swap<1,0>(x++,y++);
	    helper::swap<1,0>(x++,y++);
	    swap_start<N-3,4>(x,y);
	    break;
	  default:
	    swap_start<N  ,4>(x,y);
	  }
	}
	//
	template<unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<block==4
# ifdef __AVX__
				  && (N<8)
# endif
			       	   >::type
	swap_t(float*x, float*y) noexcept
	{
	  static_assert(N>=block,"block size mismatch");
	  helper::swap<4,y_aligned>(x,y);
	  swap_t<N-4,block_size(N-4,4),y_aligned>(x+=4,y+=4);
	}
# ifdef __AVX__
	//
	template<unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<(N>=8) && block==4>::type
	swap_t(float*x, float*y) noexcept
	{
	  static_assert(N>=block,"block size mismatch");
	  if(size_t(x)&31) {
	    helper::swap<4,y_aligned>(x,y);
	    swap_start<N-4,8>(x+=4,y+=4);
	  } else
	    swap_start<N  ,8>(x,y);
	}
	//
	template<unsigned N, unsigned block, bool y_aligned>
	static typename enable_if<block==8>::type
	swap_t(float*x, float*y) noexcept
	{
	  helper::swap<8,y_aligned>(x,y);
	  swap_t<N-8,block_size(N-8,8),y_aligned>(x+=8,y+=8);
	}
# endif
      public:
	/// unrolled for(i=0; i!=N; ++i) swap(x[i], y[i]);
	template<unsigned N>
	static void swap(float*x, float*y) noexcept
	{ swap_t<N,1,0>(x,y); }
	//
      private:
# ifdef __AVX__
#  define VECS_DECL fvec4, fvec8
# else
#  define VECS_DECL fvec4
# endif
	//
	template<typename Functor, unsigned N, unsigned block>
	static typename enable_if<N==0 && block==1>::type
	foreach_t(float*, float, VECS_DECL) noexcept {}
	//
	template<typename Functor, unsigned N, unsigned block>
	static typename enable_if<N==1 && block==1>::type
	foreach_t(float*x, const float y, VECS_DECL) noexcept
	{ helper::connect<Functor>(x,y); }
	//
	template<typename Functor, unsigned N, unsigned block>
	static typename enable_if<N==2 && block==1>::type
	foreach_t(float*x, const float y, VECS_DECL) noexcept
	{
	  helper::connect<Functor>(x++,y);
	  helper::connect<Functor>(x  ,y);
	}
	//
	template<typename Functor, unsigned N, unsigned block>
	static typename enable_if<N==3 && block==1>::type
	foreach_t(float*x, const float y, VECS_DECL) noexcept
	{
	  helper::connect<Functor>(x++,y);
	  helper::connect<Functor>(x++,y);
	  helper::connect<Functor>(x  ,y);
	}
	//
#  undef  VECS_DECL
# ifdef __AVX__
#  define VECS_DECL const float y, const fvec4&y4, const fvec8&y8
#  define VECS_PASS y, y4, y8
# else
#  define VECS_DECL const float y, const fvec4&y4
#  define VECS_PASS y, y4
# endif
	//
	template<typename Functor, unsigned N, unsigned block>
	static typename enable_if<(N>3) && block==1>::type
	foreach_t(float*x, VECS_DECL) noexcept
	{
	  static_assert(N>=4,"block size mismatch");
	  switch((size_t(x)&15)) {
	  case 12:
	    helper::connect<Functor>(x++,y);
	    foreach_t<Functor,N-1,block_size(N-1,4)>(x, VECS_PASS);
	    break;
	  case 8:
	    helper::connect<Functor>(x++,y);
	    helper::connect<Functor>(x++,y);
	    foreach_t<Functor,N-2,block_size(N-2,4)>(x, VECS_PASS);
	    break;
	  case 4:
	    helper::connect<Functor>(x++,y);
	    helper::connect<Functor>(x++,y);
	    helper::connect<Functor>(x++,y);
	    foreach_t<Functor,N-3,block_size(N-3,4)>(x, VECS_PASS);
	    break;
	  default:
	    foreach_t<Functor,N  ,block_size(N  ,4)>(x, VECS_PASS);
	  }
	}
	//
	template<typename Functor, unsigned N, unsigned block>
	static typename enable_if<block==4
# ifdef __AVX__
				  && (N<8)
# endif
			       	   >::type
	foreach_t(float*x, VECS_DECL) noexcept
	{
	  static_assert(N>=block,"block size mismatch");
	  helper::connect<Functor>(x,y4);
	  foreach_t<Functor,N-4,block_size(N-4,4)>(x+=4, VECS_PASS);
	}
# ifdef __AVX__
	//
	template<typename Functor, unsigned N, unsigned block>
	static typename enable_if<(N>=8) && block==4>::type
	foreach_t(float*x, VECS_DECL) noexcept
	{
	  static_assert(N>=block,"block size mismatch");
	  if(size_t(x)&31) {
	    helper::connect<Functor>(x,y4);
	    foreach_t<Functor,N-4,block_size(N-4,8)>(x+=4, VECS_PASS);
	  } else
	    foreach_t<Functor,N,  block_size(N,  8)>(x   , VECS_PASS);
	}
	//
	template<typename Functor, unsigned N, unsigned block>
	static typename enable_if<block==8>::type
	foreach_t(float*x, VECS_DECL) noexcept
	{
	  static_assert(N>=block,"block size mismatch");
	  helper::connect<Functor>(x,y8);
	  foreach_t<Functor,N-8,block_size(N-8,8)>(x+=8, VECS_PASS);
	}
# endif
# undef block_size
# undef VECS_DECL
# undef VECS_PASS
	//
	template<typename Functor, unsigned N>
	static typename enable_if<N==0>::type
	m_foreach(float*, float) noexcept {}
	//
	template<typename Functor, unsigned N>
	static typename enable_if<N==1>::type
	m_foreach(float*x, const float y) noexcept
	{ helper::connect<Functor>(x,y); }
	//
	template<typename Functor, unsigned N>
	static typename enable_if<N==2>::type
	m_foreach(float*x, const float y) noexcept
	{
	  helper::connect<Functor>(x++,y);
	  helper::connect<Functor>(x  ,y);
	}
	//
	template<typename Functor, unsigned N>
	static typename enable_if<N==3>::type
	m_foreach(float*x, const float y) noexcept
	{
	  helper::connect<Functor>(x++,y);
	  helper::connect<Functor>(x++,y);
	  helper::connect<Functor>(x  ,y);
	}
	//
# ifdef __AVX__
	template<unsigned N, typename Functor>
	static typename enable_if<(N>=4 && N<8)>::type
	m_foreach(float*x, float y) noexcept
	{ foreach_t<Functor,N,1>(x,y,fvec4(y),fvec8()); }
	//
	template<unsigned N, typename Functor>
	static typename enable_if<(N>=8)>::type
	m_foreach(float*x, float y) noexcept
	{ foreach_t<Functor,N,1>(x,y,fvec4(y),fvec8(y)); }
# else  //
	template<unsigned N, typename Functor>
	static typename enable_if<(N>=4)>::type
	m_foreach(float*x, float y) noexcept
	{ foreach_t<Functor,N,1>(x,y,fvec4(y)); }
# endif
      public:
	/// unrolled for(i=0; i!=N; ++i) Functor::operator(x[i], y);
	template<unsigned N, typename Functor>
	static void foreach(float*x, float y) noexcept
	{ m_foreach<N,Functor>(x,y); }
      };// struct SSE::details::static_connector<float>
#endif // __SSE__
    } // namespace WDutils::SSE::details
    ///
    /// connect two arrays element wise:
    /// @code  for(i=0; i!=n; ++i) Functor::operate(x[i], y[i]); @endcode
    ///
    template<typename Functor, unsigned N, typename RealType>
    typename enable_if< details::is_functor<Functor>::value &&
		       !is_same<Functor,meta::swap >::value>::type
    static_connect(RealType*x, const RealType*y) noexcept
    { details::static_connector<RealType>::template connect<N,Functor>(x,y); }
    /// copy array element wise
    template<unsigned N, typename RealType>
    void static_copy(RealType*x, const RealType*y) noexcept
    { static_connect<meta::assign,N>(x,y); }
    ///
    /// swap elements of two arrays:
    /// @code  for(i=0; i!=n; ++i) swap(x[i], y[i]); @endcode
    ///
    template<unsigned N, typename RealType>
    void static_swap(RealType*x, RealType*y) noexcept
    { details::static_connector<RealType>::template swap<N>(x,y); }
    //
    template<typename Functor, unsigned N, typename RealType>
    typename enable_if<is_same<Functor,meta::swap>::value>::type
    static_connect(RealType*x, RealType*y) noexcept
    { static_swap<N>(x,y); }
    ///
    /// connect each array element with (the same) scalar:
    /// @code  for(i=0; i!=n; ++i) Functor::operate(x[i], y); @endcode
    ///
    template<typename Functor, unsigned N, typename RealType>
    typename enable_if< details::is_functor<Functor >::value &&
                       !is_same<Functor,meta::swap  >::value &&
                       !is_same<Functor,meta::divide>::value>::type
    static_foreach(RealType*x, const RealType y) noexcept
    { details::static_connector<RealType>::template foreach<N,Functor>(x,y); }
    //
    template<typename Functor, unsigned N, typename RealType>
    typename enable_if<is_same<Functor,meta::divide>::value>::type
    static_foreach(RealType*x, const RealType y) noexcept
    { static_foreach<meta::multiply,N>(x,RealType(1)/y); }
#if(0)
    /// copy raw memory from srce to dest
    void*memory_copy(void*dest, const void*srce, size_t bytes)
    {
#ifdef __SSE__
      // if memory offset between dest & srce is not a multiple of 4, use memcpy
      int off_dest = size_t(dest)&31;
      int off_srce = size_t(srce)&31;
      if((off_dest-off_srce)&3)
#endif
	return std::memcpy(dest,srce,bytes);
#ifdef __SSE__
      // assert no undue overlap
# ifdef __AVX__
      WDutilsAssert(size_t(dest)+32    <= size_t(srce) &&
		    size_t(srce)+bytes <= size_t(dest));
# else
      WDutilsAssert(size_t(dest)+16    <= size_t(srce) &&
		    size_t(srce)+bytes <= size_t(dest));
# endif
      char      *p_dest = static_cast<char*> (dest);
      const char*p_srce = static_cast<const char*>(dest);
      // ensure alignment to  4 bytes
      for(; bytes>=4 && size_t(p_dest)&3; --bytes,++p_dest,++p_srce)
	*p_dest = *p_srce;
      // ensure alignment to 16 bytes
      for(; bytes>=16 && size_t(p_dest)&15; bytes-=4,p_dest+=4,p_srce+=4)
	reinterpret_cast<int32_t&>(*p_dest) =
	  reinterpret_cast<const int32_t&>(*p_srce);
# ifdef __AVX__
      // ensure alignment to 32 bytes
      if(bytes>=16 && size_t(p_dest)&31) {
	_mm_store_ps(reinterpret_cast<float*>(p_dest),
		     _mm_loadu_ps(reinterpret_cast<const float*>(p_srce)));
	bytes -=16;
	p_dest+=16;
	p_srce+=16;
      }
      // loop blocks of 32 bytes
      if(size_t(p_srce)&31) {
	for(; bytes>=32; bytes-=32,p_dest+=32,p_srce+=32)
	  _mm256_store_ps(reinterpret_cast<float*>(p_dest)
	  _mm256_loadu_ps(reinterpret_cast<const float*>(p_srce)));
	if(bytes>=16) {
	  _mm_store_ps(reinterpret_cast<float*>(p_dest),
	  _mm_loadu_ps(reinterpret_cast<const float*>(p_srce)));
	  bytes -=16;
	  p_dest+=16;
	  p_srce+=16;
	}
      } else {
	for(; bytes>=32; bytes-=32,p_dest+=32,p_srce+=32)
	  _mm256_store_ps(reinterpret_cast<float*>(p_dest)
	  _mm256_load_ps(reinterpret_cast<const float*>(p_srce)));
	if(bytes>=16) {
	  _mm_store_ps(reinterpret_cast<float*>(p_dest),
	  _mm_load_ps(reinterpret_cast<const float*>(p_srce)));
	  bytes -=16;
	  p_dest+=16;
	  p_srce+=16;
	}
      }
# else
      // loop blocks of 16 bytes
      if(size_t(p_srce)&15) {
	for(; bytes>=16; bytes-=16,p_dest+=16,p_srce+=16)
	  _mm_store_ps(reinterpret_cast<float*>(p_dest),
		       _mm_loadu_ps(reinterpret_cast<const float*>(p_srce)));
      } else {
	for(; bytes>=16; bytes-=16,p_dest+=16,p_srce+=16)
	  _mm_store_ps(reinterpret_cast<float*>(p_dest),
		       _mm_load_ps(reinterpret_cast<const float*>(p_srce)));
      }
# endif
      // do remaining blocks of 4 bytes
      for(; bytes>=4; bytes-=4,p_dest+=4,p_srce+=4)
	reinterpret_cast<int32_t&>(*p_dest) =
	  reinterpret_cast<const int32_t&>(*p_srce);
      // do remaining single bytes
      for(; bytes; --bytes,++p_dest,++p_srce)
	*p_dest = *p_srce;
#endif // __SSE__
      return dest;
    }
#endif
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
    /// extend Base by Add bytes
    template<typename Base, size_t Add>
    struct ExtendBy
    { class type : public Base { char _m_fill[Add]; }; };
    template<typename Base>
    struct ExtendBy<Base,0>
    { typedef Base type; };
    /// extension needed for alignment with K bytes
    template<typename Base, size_t K>
    class Extension
    {
      struct _tmp {
	WDutilsCXX11StaticAssert((K&(K-1))==0,"K not a power of two"); };
      static const size_t bytes = sizeof(Base);
      static const size_t splus = bytes & (K-1);
    public:
      static const size_t added = splus? K-splus : 0;
      typedef typename ExtendBy<Base,added>::type Extended;
      WDutilsStaticAssert((sizeof(Extended)&(K-1))==0);
    };
    /// extend a given class to have size a multibple of K bytes
    /// \note only default constructor possible for @a Base
#if __cplusplus >= 201103L
    template<typename Base, size_t K>
    using Extend = typename Extension<Base,K>::Extended;
    template<typename Base>
    using Extend16 = Extend<Base,16>;
#else
    template<typename Base, size_t K>
    struct Extend : public Extension<Base,K>::Extended {};
    template<typename Base>
    struct Extend16 : public Extension<Base,16>::Extended {};
#endif
    //
#ifdef __SSE__
    /// \note Not intended for consumption, use routine below instead.
    inline void _swap16(float*a, float*b, size_t n)
    {
      for(float*an=a+n; a<an; a+=4,b+=4) {
	__m128 t = _mm_load_ps(a);
	_mm_store_ps(a,_mm_load_ps(b));
	_mm_store_ps(b,t);
      }
    }
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
    /// \param[in,out] a  pter to object, on return holds data of @a y
    /// \param[in,out] b  pter to object, on return holds data of @a x
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
#undef noexcept
#undef constexpr
#undef static_assert
#undef always_inline
#undef alignas
//
#endif
