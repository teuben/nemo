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
extern "C" {
#  include <x86intrin.h>
}
#elif defined(__SSE__)
// use individual headers for intrinsics

# ifdef defined(__INTEL_COMPILER) && defined(_MM_MALLOC_H_INCLUDED)
#  warning The intel compiler has seen GNU's _mm_malloc.h which declares _mm_malloc() and _mm_free() to have different linking than those declared in INTEL's xmmintrin.h header file, which we are going to include now. This may cause a compiler error, which can be prevented by ensuring that _mm_malloc.h is not explicitly included when using the intel compiler.
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

#if(0)
  template<size_t A, typename T>
  struct alignment_traits
  {
    using element_type            = T;
    static const size_t alignment = A;
    static_assert((alignment&(alignment-1))==0,"alignment not a power of 2");
    static const bool   aligns       = aligns_at<element_type>(alignment);
    static const size_t element_size = sizeof(element_type);
    static const size_t block_size   =
      element_size>=alignment? 1:alignment/element_size;
    static const size_t block_trim   = block_size - 1;
    static const size_t block_mask   =~block_trim;
    static const bool   block_size_is_power_of_two = (block_size&block_trim)==0;
    constexpr static size_t num_blocks(size_t n) noexcept
    { return _num_blocks<block_size>(n); }
  private:
    template<size_t _block_size>
    constexpr static size_t _num_blocks(size_t n) noexcept
    { 
      static_assert(aligns,"type does not align");
      static_assert(block_size_is_power_of_two,"block size not a power of 2");
      return (n+block_trim)&block_mask;
    }
  };
#endif

#ifdef __SSE__
# if defined(__GNUC__) && !defined(__INTEL_COMPILER)
#  define always_inline __attribute__((__always_inline__))
# else
#  define always_inline
# endif
  //
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

    template<int K, typename T> class packed;

    ///
    /// packed_is_supported<a,b>::value is true if packed<a,b> is supported
    ///
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

    //--------------------------------------------------------------------------
    ///
    /// 4 packed single-precision floating-point numbers
    ///
    //--------------------------------------------------------------------------
    template<> struct packed<4,float>
    {
      /// \name types, constants, and static methods
      //@
      /// number of elements
      static const unsigned block_size = 4;
      /// associated SSE/AVX vector type
      typedef __m128 data_type;
      /// associated element type
      typedef float element_type;
      /// equivalent array of elements
      typedef element_type element_block[block_size];
      /// aligned equivalent array of elements
      typedef WDutils__align16 element_block aligned_element_block;
      /// block_size = 1<<block_sft
      static const unsigned block_shft = 2;
      /// mask for obtaining sub-index within block
      static const unsigned block_trim = 3;
      /// mask for obtaining aligned index
      static const unsigned block_mask =~block_trim;
      /// required alignement (bytes)
      static const unsigned alignment = block_size*sizeof(element_type);
      /// a packed with all elements equal to 0
      static packed always_inline zero() noexcept
      { return packed(_mm_setzero_ps()); }
      /// a packed with all elements equal to 1
      static packed always_inline one() noexcept
      { return packed(_mm_set1_ps(1.f)); }
      /// is given index aligned to block_index?
      constexpr static bool is_aligned(unsigned i) noexcept
      { return (i&block_trim)==0; }
      /// is given pointer appropriately aligned?
      constexpr static bool is_aligned(void*p) noexcept
      { return (size_t(p)&(alignment-1))==0; }
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
      //@}

      /// \name construction and assignment
      //@{
      /// default ctor
      always_inline packed() = default;
      /// move ctor
      always_inline packed(packed&&) = default;
      /// copy ctor
      always_inline packed(packed const&) = default;
#ifndef __INTEL_COMPILER
      /// move operator
      packed& always_inline operator=(packed&&) = default;
#endif
      /// copy operator
      packed& always_inline operator=(packed const&) = default;
      /// ctor from data_type
      explicit always_inline packed(data_type m) : _m(m) {}
      /// ctor from single value: set all element equal to single value
      explicit always_inline packed(element_type x) noexcept
      { _m = _mm_set1_ps(x); }
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
      friend element_type at(packed const&p) noexcept
      {
	static_assert(I<block_size,"index out of range");
# if   defined(__GNUC__) && !defined(__INTEL_COMPILER)
	return __builtin_ia32_vec_ext_v4sf(p._m,I);
# else
	union { float f; int32_t i; } tmp;
	tmp.i = _mm_extract_ps(p._m,I);
	return tmp.f;
# endif
      }
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
      /// load any aligned object that can be statically cast to const
      /// element_type*
      template<typename anything>
      static packed always_inline pack(anything const&a) noexcept
      { return load(static_cast<const element_type*>(a)); }
      /// load any unaligned object that can be statically cast to const
      /// element_type*
      template<typename anything>
      static packed always_inline packu(anything const&a) noexcept
      { return loadu(static_cast<const element_type*>(a)); }
      /// load from aligned memory location, using template arg for alignment
      template<bool aligned> static
      typename enable_if< aligned, packed>::type always_inline
      load_t(const element_type*p) noexcept
      { return packed(_mm_load_ps(p)); }
      /// load from unaligned memory location, using template arg for alignment
      template<bool aligned> static
      typename enable_if<!aligned,packed>::type always_inline
      load_t(const element_type*p) noexcept
      { return packed(_mm_loadu_ps(p)); }
      /// load any object that can be statically cast to const element_type*
      template<bool aligned, typename anything>
      static packed always_inline pack_t(anything const&a) noexcept
      { return load_t<aligned>(static_cast<const element_type*>(a)); }
      /// store to aligned memory location
      void always_inline store(element_type*p) const noexcept
      { _mm_store_ps(p,_m); }
      /// store to unaligned memory location
      void always_inline storeu(element_type*p) const noexcept
      { _mm_storeu_ps(p,_m); }
      /// store to any aligned object that can be statically cast to
      /// element_type*
      template<typename anything>
      void always_inline unpack(anything&a) const noexcept
      { store(static_cast<element_type*>(a)); }
      /// store to any unaligned object that can be statically cast to
      /// element_type*
      template<typename anything>
      void always_inline unpacku(anything&a) const noexcept
      { storeu(static_cast<element_type*>(a)); }
      /// store to aligned memory location, using template arg for alignment
      template<bool aligned>
      typename enable_if< aligned>::type always_inline
      store_t(element_type*p) const noexcept
      { _mm_store_ps(p,_m); }
      /// store to unaligned memory location, using template arg for alignment
      template<bool aligned>
      typename enable_if<!aligned>::type always_inline
      store_t(element_type*p) const noexcept
      { _mm_storeu_ps(p,_m); }
      /// store to any object that can be statically cast to element_type*
      template<bool aligned, typename anything>
      void always_inline unpack_t(anything&a) const noexcept
      { store_t<aligned>(static_cast<element_type*>(a)); }
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
	WDutils__align16 float q[4];
	_mm_store_ps(q,p._m);
	return q[0]+q[1];
      }
      /// sum of first three elements
      friend element_type always_inline sum3(packed p) noexcept
      {
	WDutils__align16 float q[4];
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

      /// \name miscellaneous
      //@{
      /// return integer with bits equal to sign bits
      friend always_inline int signbits(packed p) noexcept
      { return _mm_movemask_ps(p._m); }
      /// set all elements to Kth element of argument
      template<int K>
      friend packed always_inline single(packed p) noexcept
      {
	static_assert(K>=0 && K<block_size,"K out of range");
	return packed(_mm_shuffle_ps(p._m,p._m,_MM_SHUFFLE(K,K,K,K)));
      }
      /// result = [b2,b3,a2,a3]
      friend packed always_inline movehl(packed a, packed b) noexcept
      { return packed(_mm_movehl_ps(a._m,b._m)); }
      /// result = []
      friend packed always_inline movelh(packed a, packed b) noexcept
      { return packed(_mm_movelh_ps(a._m,b._m)); }
      /// shuffle
      /// \note inverse order to _MM_SHUFFLE (here: first is first)
      template<int I0, int I1, int I2, int I3>
      friend packed always_inline shuffle(packed a, packed b) noexcept
      {
	static_assert(I0>=0 && I0<block_size &&
		      I1>=0 && I1<block_size &&
		      I2>=0 && I2<block_size &&
		      I3>=0 && I3<block_size, "Is out of range");
	return packed(_mm_shuffle_ps(a._m,b._m,_MM_SHUFFLE(I3,I2,I1,I0)));
      }
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
      { return WDutils::meta::float_and_int(0x0).f; }
      static element_type always_inline all_mask_elem() noexcept
      { return WDutils::meta::float_and_int(0xffffffff).f; }
      static element_type always_inline sgn_mask_elem() noexcept
      { return WDutils::meta::float_and_int(0x80000000).f; }
      static element_type always_inline abs_mask_elem() noexcept
      { return WDutils::meta::float_and_int(0x7fffffff).f; }
    private:
      //
# ifdef __SSE2__
      static data_type always_inline nil_mask() noexcept
      { return _mm_castsi128_ps(_mm_set1_epi32(0x0)); }
      static data_type always_inline one_mask() noexcept
      { return _mm_castsi128_ps(_mm_set1_epi32(0xffffffff)); }
      static data_type always_inline sgn_mask() noexcept
      { return _mm_castsi128_ps(_mm_set1_epi32(0x80000000)); }
      static data_type always_inline abs_mask() noexcept
      { return _mm_castsi128_ps(_mm_set1_epi32(0x7fffffff)); }
# else
      static data_type always_inline nil_mask() noexcept
      { return _mm_set1_ps(WDutils::meta::float_and_int(0x0).f); }
      static data_type always_inline one_mask() noexcept
      { return _mm_set1_ps(WDutils::meta::float_and_int(0xffffffff).f); }
      static data_type always_inline sgn_mask() noexcept
      { return _mm_set1_ps(WDutils::meta::float_and_int(0x80000000).f); }
      static data_type always_inline abs_mask() noexcept
      { return _mm_set1_ps(WDutils::meta::float_and_int(0x7fffffff).f); }
# endif // __SSE2__
      /// just our sign bits
      data_type always_inline signmask() const noexcept
      { return _mm_and_ps(_m,sgn_mask()); }
      /// data
      data_type _m;
    };// SSE::packed<4,float>

    typedef packed<4,float> fvec4;

# ifdef __AVX__
    //--------------------------------------------------------------------------
    ///
    /// 8 packed single-precision floating-point numbers
    ///
    //--------------------------------------------------------------------------
    template<> struct packed<8,float>
    {
      /// \name types, constants, and static methods
      //@
      /// associated SSE/AVX vector type
      typedef __m256 data_type;
      /// number of elements
      static const unsigned block_size = 8;
      /// associated element type
      typedef float element_type;
      /// equivalent array of elements
      typedef element_type element_block[block_size];
      /// aligned equivalent array of elements
      typedef WDutils__align32 element_block aligned_element_block;
      /// block_size = 1<<block_sft
      static const unsigned block_shft = 3;
      /// mask for obtaining sub-index within block
      static const unsigned block_trim = 7;
      /// mask for obtaining aligned index
      static const unsigned block_mask =~block_trim;
      /// required alignement (bytes)
      static const unsigned alignment = block_size*sizeof(element_type);
      /// a packed with all elements equal to 0
      static packed always_inline zero() noexcept
      { return packed(_mm256_setzero_ps()); }
      /// a packed with all elements equal to 1
      static packed always_inline one() noexcept
      { return packed(_mm256_set1_ps(1.f)); }
      /// is given index aligned to block_index?
      constexpr static bool is_aligned(unsigned i) noexcept
      { return (i&block_trim)==0; }
      /// is given pointer appropriately aligned?
      constexpr static bool is_aligned(void*p) noexcept
      { return (size_t(p)&(alignment-1))==0; }
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
      //@}

      /// \name construction and assignment
      //@{
      /// default ctor
      always_inline packed() = default;
      /// move ctor
      always_inline packed(packed&&) = default;
      /// copy ctor
      always_inline packed(packed const&) = default;
#ifndef __INTEL_COMPILER
      /// move operator
      packed& always_inline operator=(packed&&) = default;
#endif
      /// copy operator
      packed& always_inline operator=(packed const&) = default;
      /// ctor from data_type
      explicit always_inline packed(data_type m) : _m(m) {}
      /// ctor from single value: set all element equal to single value
      explicit always_inline packed(element_type x) noexcept
      { _m = _mm256_set1_ps(x); }
      /// ctor from 4 values: set elements
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
      packed& always_inline set(element_type x) noexcept
      { _m = _mm256_set1_ps(x); return*this; }
      /// set elements
      /// \note inverse order to _mm256_set_ps(a,b,c,d,e,f,g,h)
      packed& always_inline set(element_type a, element_type b,
			   element_type c, element_type d,
			   element_type e, element_type f,
			   element_type g, element_type h) noexcept
      { _m = _mm256_setr_ps(a,b,c,d,e,f,g,h); return*this; }
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
      /// obtain packed<4,float> from elements [I,...,I+3]
      template<unsigned I>
      packed<4,float> always_inline extract() const noexcept
      {
	static_assert(I<block_size,"index out of range");
	return packed<4,float>(_mm256_extractf128_ps(_m,I));
      }
      /// constant element access, templated
      template<unsigned I>
      friend element_type at(packed const&p) noexcept
      {
	static_assert(I<block_size,"index out of range");
	return at<(I&3)>(extract<(I>>2)>());
      }
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
      template<bool aligned> static
      typename enable_if< aligned, packed>::type always_inline
      load_t(const element_type*p) noexcept
      { return packed(_mm256_load_ps(p)); }
      /// load any aligned object that can be statically cast to const
      /// element_type*
      template<typename anything>
      static packed always_inline pack(anything const&a) noexcept
      { return load(static_cast<const element_type*>(a)); }
      /// load any unaligned object that can be statically cast to const
      /// element_type*
      template<typename anything>
      static packed always_inline packu(anything const&a) noexcept
      { return loadu(static_cast<const element_type*>(a)); }
      /// load from unaligned memory location, using template arg for alignment
      template<bool aligned> static
      typename enable_if<!aligned,packed>::type always_inline
      load_t(const element_type*p) noexcept
      { return packed(_mm256_loadu_ps(p)); }
      /// load any object that can be statically cast to const element_type*
      template<bool aligned, typename anything>
      static packed always_inline pack_t(anything const&a) noexcept
      { return load_t<aligned>(static_cast<const element_type*>(a)); }
      /// store to aligned memory location
      void always_inline store(element_type*p) const noexcept
      { _mm256_store_ps(p,_m); }
      /// store to unaligned memory location
      void always_inline storeu(element_type*p) const noexcept
      { _mm256_storeu_ps(p,_m); }
      /// store to any aligned object that can be statically cast to
      /// element_type*
      template<typename anything>
      void always_inline unpack(anything&a) const noexcept
      { store(static_cast<element_type*>(a)); }
      /// store to any unaligned object that can be statically cast to
      /// element_type*
      template<typename anything>
      void always_inline unpacku(anything&a) const noexcept
      { storeu(static_cast<element_type*>(a)); }
      /// store to aligned memory location, using template arg for alignment
      template<bool aligned>
      typename enable_if< aligned>::type always_inline
      store_t(element_type*p) const noexcept
      { _mm256_store_ps(p,_m); }
      /// store to unaligned memory location, using template arg for alignment
      template<bool aligned>
      typename enable_if<!aligned>::type always_inline
      store_t(element_type*p) const noexcept
      { _mm256_storeu_ps(p,_m); }
      /// store to any object that can be statically cast to element_type*
      template<bool aligned, typename anything>
      void always_inline unpack_t(anything&a) const noexcept
      { store_t<aligned>(static_cast<element_type*>(a)); }
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
#  if(0)
      /// \name horizontal arithmetic operations
      //@{
      /// horizontal sum
      friend element_type always_inline sum(packed p) noexcept
      /// horizontal sum in each element
      friend packed always_inline hadd(packed p) noexcept
      /// horizontal maximum
      friend element_type always_inline max(packed p) noexcept
      /// horizontal minimum
      friend element_type always_inline min(packed p) noexcept
      //@}
#  endif
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

      /// \name miscellaneous
      //@{
      /// return integer with bits equal to sign bits
      friend always_inline int signbits(packed p) noexcept
      { return _mm256_movemask_ps(p._m); }
      /// blend two vectors depending on sign of third:  result = sign<0? x : y
      friend packed always_inline blend(packed sign, packed x, packed y)
	noexcept
      { return packed(_mm256_blendv_ps(y._m,x._m,sign._m)); }
      /// combine two vectors depending on third:  result = mask? x : y
      friend packed always_inline combine(packed mask, packed x, packed y)
	noexcept
      { return packed(_mm256_blendv_ps(y._m,x._m,mask._m)); }
      //@}
    private:
      static data_type always_inline nil_mask() noexcept
      { return _mm256_castsi256_ps(_mm256_set1_epi32(0x0)); }
      static data_type always_inline one_mask() noexcept
      { return _mm256_castsi256_ps(_mm256_set1_epi32(0xffffffff)); }
      static data_type always_inline sgn_mask() noexcept
      { return _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000)); }
      static data_type always_inline abs_mask() noexcept
      { return _mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff)); }
      /// just our sign bits
      data_type always_inline signmask() const noexcept
      { return _mm256_and_ps(_m,sgn_mask()); }
      /// data
      data_type _m;
    };// SSE::packed<8,float>
# endif // __AVX__

    typedef packed<8,float> fvec8;

# ifdef __SSE2__
    //--------------------------------------------------------------------------
    ///
    /// 2 packed double-precision floating-point numbers
    ///
    //--------------------------------------------------------------------------
    template<> struct packed<2,double>
    {
    public:
      /// \name types, constants, and static methods
      //@
      /// associated element type
      typedef double element_type;
      /// associated SSE/AVX vector type
      typedef __m128d data_type;
      /// number of element types hold
      static const unsigned block_size = 2;
      /// equivalent array of elements
      typedef element_type element_block[block_size];
      /// aligned equivalent array of elements
      typedef WDutils__align16 element_block aligned_element_block;
      /// block_size = 1<<block_sft
      static const unsigned block_shft = 1;
      /// mask for obtaining sub-index within block
      static const unsigned block_trim = 1;
      /// mask for obtaining aligned index
      static const unsigned block_mask =~block_trim;
      /// required alignement (bytes)
      static const unsigned alignment = block_size*sizeof(element_type);
      /// a packed with all elements equal to 0
      static packed always_inline zero() noexcept
      { return packed(_mm_setzero_pd()); }
      /// a packed with all elements equal to 1
      static packed always_inline one() noexcept
      { return packed(_mm_set1_pd(1.0)); }
      /// is given index aligned to block_index?
      constexpr static bool is_aligned(unsigned i) noexcept
      { return (i&block_trim)==0; }
      /// is given pointer appropriately aligned?
      constexpr static bool is_aligned(void*p) noexcept
      { return (size_t(p)&(alignment-1))==0; }
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
      //@}

      /// \name construction and assignment
      //@{
      /// default ctor
      always_inline packed() = default;
      /// move ctor
      always_inline packed(packed&&) = default;
      /// copy ctor
      always_inline packed(packed const&) = default;
#ifndef __INTEL_COMPILER
      /// move operator
      packed& always_inline operator=(packed&&) = default;
#endif
      /// copy operator
      packed& always_inline operator=(packed const&) = default;
      /// ctor from data_type is private
      explicit always_inline packed(__m128d m) : _m(m) {}
      /// ctor from single value: set all element equal to single value
      explicit always_inline packed(element_type x) noexcept
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
      friend element_type at(packed const&p) noexcept
      {
	static_assert(I<block_size,"index out of range");
# if   defined(__GNUC__) && !defined(__INTEL_COMPILER)
	return __builtin_ia32_vec_ext_v2df(p._m,I);
# else
	union { double d; int64_t i; } tmp;
	tmp.i = _mm_extract_pd(p._m,I);
	return tmp.d;
# endif
      }
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
      /// load any aligned object that can be statically cast to const
      /// element_type*
      template<typename anything>
      static packed always_inline pack(anything const&a) noexcept
      { return load(static_cast<const element_type*>(a)); }
      /// load any unaligned object that can be statically cast to const
      /// element_type*
      template<typename anything>
      static packed always_inline packu(anything const&a) noexcept
      { return loadu(static_cast<const element_type*>(a)); }
      /// load from aligned memory location, using template arg for alignment
      template<bool aligned> static
      typename enable_if< aligned, packed>::type always_inline
      load_t(const element_type*p) noexcept
      { return packed(_mm_load_pd(p)); }
      /// load from unaligned memory location, using template arg for alignment
      template<bool aligned> static
      typename enable_if<!aligned,packed>::type always_inline
      load_t(const element_type*p) noexcept
      { return packed(_mm_loadu_pd(p)); }
      /// load any object that can be statically cast to const element_type*
      template<bool aligned, typename anything>
      static packed always_inline pack_t(anything const&a) noexcept
      { return load_t<aligned>(static_cast<const element_type*>(a)); }
      /// store to aligned memory location
      void always_inline store(element_type*p) const noexcept
      { _mm_store_pd(p,_m); }
      /// store to unaligned memory location
      void always_inline storeu(element_type*p) const noexcept
      { _mm_storeu_pd(p,_m); }
      /// store to any aligned object that can be statically cast to
      /// element_type*
      template<typename anything>
      void always_inline unpack(anything&a) const noexcept
      { store(static_cast<element_type*>(a)); }
      /// store to any unaligned object that can be statically cast to
      /// element_type*
      template<typename anything>
      void always_inline unpacku(anything&a) const noexcept
      { storeu(static_cast<element_type*>(a)); }
      /// store to aligned memory location, using template arg for alignment
      template<bool aligned>
      typename enable_if< aligned>::type always_inline
      store_t(element_type*p) const noexcept
      { _mm_store_pd(p,_m); }
      /// store to unaligned memory location, using template arg for alignment
      template<bool aligned>
      typename enable_if<!aligned>::type always_inline
      store_t(element_type*p) const noexcept
      { _mm_storeu_pd(p,_m); }
      /// store to any object that can be statically cast to element_type*
      template<bool aligned, typename anything>
      void always_inline unpack_t(anything&a) const noexcept
      { store_t<aligned>(static_cast<element_type*>(a)); }
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
	WDutils__align16 double q[2];
	_mm_store_pd(q,p._m);
	return q[0]+q[1];
      }
      /// horizontal sum of first element
      friend element_type always_inline sum1(packed p) noexcept
      {
	WDutils__align16 double q[2];
	_mm_store_pd(q,p._m);
	return q[0];
      }
#  if(0)
      /// horizontal sum in each element
      friend packed always_inline hadd(packed p) noexcept
      /// horizontal maximum
      friend element_type always_inline max(packed p) noexcept
      /// horizontal minimum
      friend element_type always_inline min(packed p) noexcept
#  endif
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

      /// \name miscellaneous
      //@{
      /// return integer with bits equal to sign bits
      friend always_inline int signbits(packed p) noexcept
      { return _mm_movemask_pd(p._m); }
      /// set all elements to Kth element of argument
      template<int K>
      friend packed always_inline single(packed p) noexcept
      {
	static_assert(K>=0 && K<block_size,"K out of range");
	return packed(_mm_shuffle_pd(p._m,p._m,_MM_SHUFFLE2(K,K)));
      }
      /// blend two vectors depending on sign of third:  result = sign<0? x : y
      friend packed always_inline blend(packed sign, packed x, packed y)
	noexcept
      { 
#  ifdef __SSE4_1__
	return packed(_mm_blendv_pd(y._m,x._m,sign._m));
#  else
	// perhapd there is a better way using move and movemask?
        __m128 mask = _mm_cmplt_pd(sign._m,_mm_setzero_pd());
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
    private:
      static __m128d always_inline nil_mask() noexcept
      { return _mm_castsi128_pd(_mm_set1_epi64x(0x0)); }
      static __m128d always_inline one_mask() noexcept
      { return _mm_castsi128_pd(_mm_set1_epi64x(0xffffffffffffffff)); }
      static __m128d always_inline sgn_mask() noexcept
      { return _mm_castsi128_pd(_mm_set1_epi64x(0x8000000000000000)); }
      static __m128d always_inline abs_mask() noexcept
      { return _mm_castsi128_pd(_mm_set1_epi64x(0x7fffffffffffffff)); }
      /// just our sign bits
      data_type always_inline signmask() const noexcept
      { return _mm_and_pd(_m,sgn_mask()); }
      /// data
      data_type _m;
    };// SSE::packed<2,double>
# endif // __SSE2__

    typedef packed<2,double> dvec2;

# ifdef __AVX__
    //--------------------------------------------------------------------------
    ///
    /// 4 packed double-precision floating-point numbers
    ///
    //--------------------------------------------------------------------------
    template<> struct packed<4,double>
    {
      /// \name types, constants, and static methods
      //@
      /// associated element type
      typedef double element_type;
      /// associated SSE/AVX vector type
      typedef __m256d data_type;
      /// number of elements
      static const unsigned block_size = 4;
      /// equivalent array of elements
      typedef element_type element_block[block_size];
      /// aligned equivalent array of elements
      typedef WDutils__align32 element_block aligned_element_block;
      /// block_size = 1<<block_sft
      static const unsigned block_shft = 2;
      /// mask for obtaining sub-index within block
      static const unsigned block_trim = 3;
      /// mask for obtaining aligned index
      static const unsigned block_mask =~block_trim;
      /// required alignement (bytes)
      static const unsigned alignment = block_size*sizeof(element_type);
      /// a packed with all elements equal to 0
      static packed always_inline zero() noexcept
      { return packed(_mm256_setzero_pd()); }
      /// a packed with all elements equal to 1
      static packed always_inline one() noexcept
      { return packed(_mm256_set1_pd(1.0)); }
      /// is given index aligned to block_index?
      constexpr static bool is_aligned(unsigned i) noexcept
      { return (i&block_trim)==0; }
      /// is given pointer appropriately aligned?
      constexpr static bool is_aligned(void*p) noexcept
      { return (size_t(p)&(alignment-1))==0; }
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
      //@}

      /// \name construction and assignment
      //@{
      /// default ctor
      always_inline packed() = default;
      /// move ctor
      always_inline packed(packed&&) = default;
      /// copy ctor
      always_inline packed(packed const&) = default;
#ifndef __INTEL_COMPILER
      /// move operator
      packed& always_inline operator=(packed&&) = default;
#endif
      /// copy operator
      packed& always_inline operator=(packed const&) = default;
      /// ctor from data_type is private
      explicit always_inline packed(data_type m) : _m(m) {}
      /// ctor from single value: set all element equal to single value
      explicit always_inline packed(element_type x) noexcept
      { _m = _mm256_set1_pd(x); }
      /// ctor from 4 values: set elements
      /// \note inverse ordre to _mm256_set_pd()
      always_inline packed(element_type a, element_type b,
			   element_type c, element_type d) noexcept
      { _m = _mm256_set_pd(d,c,b,a); }
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
      /// obtain packed<2,double> from elements [I,I+1]
      template<unsigned I>
      packed<2,double> always_inline extract() const noexcept
      {
	static_assert(I<block_size,"index out of range");
	return packed<2,double>(_mm256_extractf128_pd(_m,I));
      }
      /// constant element access, templated
      template<unsigned I>
      friend element_type at(packed const&p) noexcept
      {
	static_assert(I<block_size,"index out of range");
	return at<(I&1)>(extract<(I>>1)>());
      }
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
      /// load any aligned object that can be statically cast to const
      /// element_type*
      template<typename anything>
      static packed always_inline pack(anything const&a) noexcept
      { return load(static_cast<const element_type*>(a)); }
      /// load any unaligned object that can be statically cast to const
      /// element_type*
      template<typename anything>
      static packed always_inline packu(anything const&a) noexcept
      { return loadu(static_cast<const element_type*>(a)); }
      /// load from aligned memory location, using template arg for alignment
      template<bool aligned> static
      typename enable_if< aligned, packed>::type always_inline
      load_t(const element_type*p) noexcept
      { return packed(_mm256_load_pd(p)); }
      /// load from unaligned memory location, using template arg for alignment
      template<bool aligned> static
      typename enable_if<!aligned,packed>::type always_inline
      load_t(const element_type*p) noexcept
      { return packed(_mm256_loadu_pd(p)); }
      /// load from any object that can be statically cast to const element*
      template<bool aligned, typename anything>
      static packed always_inline pack_t(anything const&a) noexcept
      { return load_t<aligned>(static_cast<const element_type*>(a)); }
      /// store to aligned memory location
      void always_inline store(element_type*p) const noexcept
      { _mm256_store_pd(p,_m); }
      /// store to unaligned memory location
      void always_inline storeu(element_type*p) const noexcept
      { _mm256_storeu_pd(p,_m); }
      /// store to any aligned object that can be statically cast to
      /// element_type*
      template<typename anything>
      void always_inline unpack(anything&a) const noexcept
      { store(static_cast<element_type*>(a)); }
      /// store to any unaligned object that can be statically cast to
      /// element_type*
      template<typename anything>
      void always_inline unpacku(anything&a) const noexcept
      { storeu(static_cast<element_type*>(a)); }
      /// store to aligned memory location, using template arg for alignment
      template<bool aligned>
      typename enable_if< aligned>::type always_inline
      store_t(element_type*p) const noexcept
      { _mm256_store_pd(p,_m); }
      /// store to unaligned memory location, using template arg for alignment
      template<bool aligned>
      typename enable_if<!aligned>::type always_inline
      store_t(element_type*p) const noexcept
      { _mm256_storeu_pd(p,_m); }
      /// store to any object that can be statically cast to element*
      template<bool aligned, typename anything>
      void always_inline unpack_t(anything&a) const noexcept
      { store_t<aligned>(static_cast<element_type*>(a)); }
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
	WDutils__align32 double q[4];
	_mm256_store_pd(q,p._m);
	return q[0]+q[1]+q[2]+q[3];
      }
      /// sum of first three elements
      friend element_type always_inline sum3(packed p) noexcept
      {
	WDutils__align32 double q[4];
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

      /// \name miscellaneous
      //@{
      /// return integer with bits equal to sign bits
      friend always_inline int signbits(packed p) noexcept
      { return _mm256_movemask_pd(p._m); }
      /// set all elements to Kth element of argument
      template<int K>
      friend packed always_inline single(packed p) noexcept
      {
	static_assert(K>=0 && K<block_size,"K out of range");
#ifdef __AVX2__
# warning more efficient implementation possible
#endif
	return packed(_mm256_permute_pd
		      (_mm256_permute2f128_pd(p._m,p._m,K&2? 49:32),K&1?15:0));
      }
      /// blend two vectors depending on sign of third:  result = sign<0? x : y
      friend packed blend(packed sign, packed x, packed y) noexcept
      { return packed(_mm256_blendv_pd(y._m,x._m,sign._m)); }
      /// combine two vectors depending on third:  result = mask? x : y
      friend packed combine(packed mask, packed x, packed y) noexcept
      { return packed(_mm256_blendv_pd(y._m,x._m,mask._m)); }
      //@}
    private:
      static data_type always_inline nil_mask() noexcept
      { return _mm256_castsi256_pd(_mm256_set1_epi64x(0x0)); }
      static data_type always_inline one_mask() noexcept
      { return _mm256_castsi256_pd(_mm256_set1_epi64x(0xffffffffffffffff)); }
      static data_type always_inline sgn_mask() noexcept
      { return _mm256_castsi256_pd(_mm256_set1_epi64x(0x8000000000000000)); }
      static data_type always_inline abs_mask() noexcept
      { return _mm256_castsi256_pd(_mm256_set1_epi64x(0x7fffffffffffffff)); }
      /// just our sign bits
      data_type always_inline signmask() const noexcept
      { return _mm256_and_pd(_m,sgn_mask()); }
      /// data
      data_type _m;
    };// SSE::packed<4,double>
    typedef packed<4,double> dvec4;

# endif // __AVX__
#endif  // __SSE__
  } // namespace WDutils::SSE
  //
#ifdef __SSE__
# ifdef __INTEL_COMPILER
  inline float xmm0(__m128 _A)
  { union { float f; int i; } t; t.i = _mm_extract_ps(_A,0); return t.f; }
  inline float xmm1(__m128 _A)
  { union { float f; int i; } t; t.i = _mm_extract_ps(_A,1); return t.f; }
  inline float xmm2(__m128 _A)
  { union { float f; int i; } t; t.i = _mm_extract_ps(_A,2); return t.f; }
  inline float xmm3(__m128 _A)
  { union { float f; int i; } t; t.i = _mm_extract_ps(_A,3); return t.f; }
#  ifdef __SSE2__
  inline double xmm0(__m128d _A)
  { double x; _mm_storel_pd(&x,_A); return x; }
  inline double xmm1(__m128d _A)
  { double x; _mm_storeh_pd(&x,_A); return x; }
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
# ifdef __AVX__
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
# endif // __AVX__
#endif // __SSE__
} // namespace WDutils
//
#undef noexcept
#undef constexpr
//

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
    /// filler object of size _S
    template<size_t _S> class Filler { char _F[_S]; };
    /// extend a given class to have size a multibple of K bytes
    /// \note only default constructor possible for @a Base
    template<typename Base
#if __cplusplus >= 201103L
	     , int K= 16
#endif
	     >
class
#if __cplusplus >= 201103L
    Extend
#else
    Extend16
#endif
     : public Base
    {
#if __cplusplus >= 201103L
      static_assert((K&(K-1))==0,"K not a power of two");
#else
      static const int K=16;
#endif
      static const size_t Bytes = sizeof(Base);
      static const size_t Splus = Bytes & (K-1);
      static const size_t Added = Splus? K-Splus : 0;
      Filler<Added> _Filler;
    };
#if __cplusplus >= 201103L
    template<typename Base>
    using Extend16 = Extend<Base,16>;
#endif
    //
#ifdef __SSE__
    /// working horse actually implemented in sse.cc
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
