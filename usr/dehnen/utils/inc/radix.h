// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/radix.h
///
/// \author Walter Dehnen
/// \date   2011
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but          
// WITHOUT ANY WARRANTY; without even the implied warranty of                   
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
// General Public License for more details.                                     
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 675
// Mass Ave, Cambridge, MA 02139, USA.
//
////////////////////////////////////////////////////////////////////////////////
#ifndef WDutils_included_radix_h
#define WDutils_included_radix_h

#ifndef WDutils_included_traits_h
#  include <traits.h>
#endif
#ifndef WDutils_included_memory_h
#  include <memory.h>
#endif
#ifndef WDutils_included_limits
#  define WDutils_included_limits
#  include <limits>
#endif
#ifdef WDutilsDevel
# ifndef WDutils_included_parallel_h
#  include <parallel.h>
# endif
#endif
////////////////////////////////////////////////////////////////////////////////
namespace WDutils {
  /// methods for stable radix sorting of integers and floating point numbers
  namespace Radix {
    template<int, typename traits> struct SortBits;
    /// stable radix sort of the lowest 8 bits
    template<typename traits> struct SortBits<8,traits>
    {
      typedef typename traits::input_type   input_type;
      typedef typename traits::integer_type integer_type;
      WDutilsStaticAssert(std::numeric_limits<integer_type>::is_integer);
      static void sort(unsigned N, input_type*X, input_type*Y)
      {
	uint32_t B0[256] = {0};
	static const integer_type m=0xff;
	for(uint32_t i=0; i!=N; ++i) {
	  register integer_type f=traits::integer(Y[i]=X[i]);
	  B0[f&m]++;
	}
	for(uint32_t i=0,s=0,t; i!=256; ++i) { t=B0[i]; B0[i]=s; s+=t; }
	for(uint32_t i=0; i!=N; ++i) X[B0[traits::integer(Y[i])&m]++]=Y[i];
      }
    };
    /// stable radix sort of the lowest 16 bits
    template<typename traits> struct SortBits<16,traits>
    {
      typedef typename traits::input_type   input_type;
      typedef typename traits::integer_type integer_type;
      WDutilsStaticAssert(std::numeric_limits<integer_type>::is_integer &&
			  sizeof(integer_type)>=2);
      static void sort(unsigned N, input_type*X, input_type*Y)
      {
	uint32_t B0[512]={0},*B1=B0+256;
	static const integer_type m=0xff;
	for(uint32_t i=0; i!=N; ++i) {
	  register integer_type f=traits::integer(X[i]);
	  B0[ f     &m] ++;
	  B1[(f>>=8)&m] ++;
	}
	for(uint32_t i=0,s=0,t; i!=256; ++i) { t=B0[i]; B0[i]=s; s+=t; }
	for(uint32_t i=0,s=0,t; i!=256; ++i) { t=B1[i]; B1[i]=s; s+=t; }
	for(uint32_t i=0; i!=N; ++i) Y[B0[traits::integer(X[i])   &m]++]=X[i];
	for(uint32_t i=0; i!=N; ++i) X[B1[traits::integer(Y[i])>>8&m]++]=Y[i];
      }
    };
    /// stable radix sort of the lowest 32 bits
    template<typename traits> struct SortBits<32,traits>
    {
      typedef typename traits::input_type   input_type;
      typedef typename traits::integer_type integer_type;
      WDutilsStaticAssert(std::numeric_limits<integer_type>::is_integer &&
			  sizeof(integer_type)>=4);
      static void sort(unsigned N, input_type*X, input_type*Y)
      {
	static const uint32_t  n0=0x800,n2=0x400;
	static const integer_type m0=0x7ff,m2=0x3ff;
	uint32_t B0[n0+n0+n2]={0},*B1=B0+n0,*B2=B1+n0;
	for(uint32_t i=0; i!=N; ++i) {
	  register integer_type f=traits::integer(Y[i]=X[i]);
	  B0[ f      &m0] ++;
	  B1[(f>>=11)&m0] ++;
	  B2[(f>>=11)&m2] ++;
	}
	for(uint32_t i=0,s=0,t; i!=n0; ++i) { t=B0[i]; B0[i]=s; s+=t; }
	for(uint32_t i=0,s=0,t; i!=n0; ++i) { t=B1[i]; B1[i]=s; s+=t; }
	for(uint32_t i=0,s=0,t; i!=n2; ++i) { t=B2[i]; B2[i]=s; s+=t; }
	for(uint32_t i=0; i!=N; ++i)
	  X[B0[traits::integer(Y[i])     &m0]++]=Y[i];
	for(uint32_t i=0; i!=N; ++i)
	  Y[B1[traits::integer(X[i])>>11 &m0]++]=X[i];
	for(uint32_t i=0; i!=N; ++i)
	  X[B2[traits::integer(Y[i])>>22 &m2]++]=Y[i];
      }
    };
    /// stable radix sort of the lowest 64 bits
    template<typename traits> struct SortBits<64,traits>
    {
      typedef typename traits::input_type   input_type;
      typedef typename traits::integer_type integer_type;
      WDutilsStaticAssert(std::numeric_limits<integer_type>::is_integer &&
			  sizeof(integer_type)>=8);
      static void sort(unsigned N, input_type*X, input_type*Y)
      {
	static const uint32_t  n0=0x800,n2=0x400,n=n0+n0+n2;
	static const integer_type m0=0x7ff,m2=0x3ff;
	uint32_t B0[n]={0},*B1=B0+n0,*B2=B1+n0;
	// sort the lower 32 bits
	for(uint32_t i=0; i!=N; ++i) {
	  register uint64_t f=traits::integer(X[i]);
	  B0[ f      &m0] ++;
	  B1[(f>>=11)&m0] ++;
	  B2[(f>>=11)&m2] ++;
	}
	for(uint32_t i=0,s=0,t; i!=n0; ++i) { t=B0[i]; B0[i]=s; s+=t; }
	for(uint32_t i=0,s=0,t; i!=n0; ++i) { t=B1[i]; B1[i]=s; s+=t; }
	for(uint32_t i=0,s=0,t; i!=n2; ++i) { t=B2[i]; B2[i]=s; s+=t; }
	for(uint32_t i=0; i!=N; ++i)
	  Y[B0[traits::integer(X[i])     &m0]++]=X[i];
	for(uint32_t i=0; i!=N; ++i)
	  X[B1[traits::integer(Y[i])>>11 &m0]++]=Y[i];
	for(uint32_t i=0; i!=N; ++i)
	  Y[B2[traits::integer(X[i])>>22 &m2]++]=X[i];
	// sort the upper 32 bits
	for(uint32_t i=0; i!=n; ++i) B0[i]=0;
	for(uint32_t i=0; i!=N; ++i) {
	  register uint64_t f=traits::integer(Y[i]);
	  B0[(f>>=32)&m0] ++;
	  B1[(f>>=11)&m0] ++;
	  B2[(f>>=11)&m2] ++;
	}
	for(uint32_t i=0,s=0,t; i!=n0; ++i) { t=B0[i]; B0[i]=s; s+=t; }
	for(uint32_t i=0,s=0,t; i!=n0; ++i) { t=B1[i]; B1[i]=s; s+=t; }
	for(uint32_t i=0,s=0,t; i!=n2; ++i) { t=B2[i]; B2[i]=s; s+=t; }
	for(uint32_t i=0; i!=N; ++i)
	  X[B0[traits::integer(Y[i])>>32 &m0]++]=Y[i];
	for(uint32_t i=0; i!=N; ++i)
	  Y[B1[traits::integer(X[i])>>43 &m0]++]=X[i];
	for(uint32_t i=0; i!=N; ++i)
	  X[B2[traits::integer(Y[i])>>54 &m2]++]=Y[i];
      }
    };

    /// traits for use in @c SortBits
    template<typename> struct SortTraits;

    /// for implementing @c SortTraits of integer types
    template<typename __I> struct SortTraitsInteger
    {
      WDutilsStaticAssert(std::numeric_limits<__I>::is_integer);
      typedef __I input_type;
      typedef __I integer_type;
      static integer_type integer(input_type x) { return x; }
    };

#define WDutilsRadixSortTraits(INT)	\
    template<> struct SortTraits<INT>	\
      : public SortTraitsInteger<INT> {}
    WDutilsRadixSortTraits(int8_t);
    WDutilsRadixSortTraits(int16_t);
    WDutilsRadixSortTraits(int32_t);
    WDutilsRadixSortTraits(int64_t);
    WDutilsRadixSortTraits(uint8_t);
    WDutilsRadixSortTraits(uint16_t);
    WDutilsRadixSortTraits(uint32_t);
    WDutilsRadixSortTraits(uint64_t);
#undef WDutilsRadixSortTraits
    /// SortTraits for floats
    /// \detail
    /// for implementing Sort of structs with a float member to be sorted on
    /// \note the order-preserving map from float to uint32_t is due to Michael
    ///       Herf (see http://www.stereopsis.com/radix.html)
    template<> struct SortTraits<float>
    {
      typedef float  input_type;
      typedef uint32_t integer_type;
      static integer_type integer(input_type const&x)
      { return forward(reinterpret_cast<integer_type const&>(x)); }
      static integer_type forward_map(integer_type f)
      { return -int32_t(f>>31) | 0x80000000; }
      static integer_type backward_map(integer_type f)
      { return ((f>>31) - 1) | 0x80000000; }
      static integer_type forward(integer_type f)
      { return f ^ forward_map(f); }
      static integer_type backward(integer_type f)
      { return f ^ backward_map(f); }
    };
    /// SortTraits for doubles
    /// \detail
    /// for implementing Sort of structs with a double member to be sorted on
    /// \note the order-preserving map from double to uint64_t is based on the
    ///       equivalent map between float and uint32_t due to Michael Herf see
    ///       http://www.stereopsis.com/radix.html
    template<> struct SortTraits<double>
    {
      typedef double input_type;
      typedef uint64_t integer_type;
      static integer_type integer(double const& x)
      { return forward(reinterpret_cast<integer_type const&>(x)); }
      static integer_type forward_map(integer_type f)
      { return -int64_t(f>>63) | 0x8000000000000000ll; }
      static integer_type backward_map(integer_type f)
      { return ((f>>63) - 1) | 0x8000000000000000ll; }
      static integer_type forward(integer_type f)
      { return f ^ forward_map(f); }
      static integer_type backward(integer_type f)
      { return f ^ backward_map(f); }
    };
    /// radix sort of any type with a SortTraits<>
    /// \param[in]      N  # elements
    /// \param[in,out]  X  array, sorted on return
    /// \param[in]      Y  auxiliary array of @a N elements
    /// \note radix sort provides stable sorting and costs O(N) time
    template<typename __T>
    inline void Sort(unsigned N, __T*X, __T*Y)
    { SortBits<8*sizeof(__T),SortTraits<__T> >::sort(N,X,Y); }
    /// radix sort of the lower K bits of any type with a SortTraits<>
    /// \param[in]      N  # elements
    /// \param[in,out]  X  array, sorted on return
    /// \param[in]      Y  auxiliary array of @a N elements
    /// \note radix sort provides stable sorting and costs O(N) time
    template<int K, typename __T>
    inline void SortLow(unsigned N, __T*X, __T*Y)
    { SortBits<K,SortTraits<__T> >::sort(N,X,Y); }
    /// radix sort of single-precision floating point numbers
    /// \param[in]      N  # elements
    /// \param[in,out]  X  array, sorted on return
    /// \param[in]      Y  auxiliary array of @a N @c floats
    /// \note radix sort provides stable sorting and costs O(N) time
    /// \note we map to 32-bit unsigned integers, sort them, and then map back.
    ///       the order-preserving map from float to uint32_t is due to Michael
    ///       Herf (see http://www.stereopsis.com/radix.html)
    void Sort(unsigned N, float*X, float*Y);
    /// radix sort of double-precision floating point numbers
    /// \param[in]      N  # elements
    /// \param[in,out]  X  array, sorted on return
    /// \param[in]      Y  auxiliary array of @a N @c doubles
    /// \note radix sort provides stable sorting and costs O(N) time
    /// \note we map to 64-bit unsigned integers, sort them, and then map back.
    ///       the order-preserving map from double to uint64_t is based on the
    ///       equivalent map between float and uint32_t due to Michael Herf see
    ///       http://www.stereopsis.com/radix.html
    void Sort(unsigned N, double*X, double*Y);
    /// radix sort
    /// \param[in]      N  # elements
    /// \param[in,out]  X  array, sorted on return
    /// \note radix sort provides stable sorting and costs O(N) time
    /// \note allocates and de-allocates an auxiliary array
    template<typename __T>
    inline void Sort(unsigned N, __T*X)
    {
      __T*Y=WDutils_NEW16(__T,N);
      Sort(N,X,Y);
      WDutils_DEL16(Y);
    }
#ifdef WDutilsDevel
    /// \name these parallel radix sorts scale badly
    //@{
    /// parallel radix sort of any type with a SortTraits<>
    /// \param[in]      N  # elements
    /// \param[in,out]  X  array, sorted on return
    /// \param[in]      Y  auxiliary array of @a N elements
    /// \param[in]   warn  warn about call from within parallel region?
    /// \note radix sort provides stable sorting and costs O(N) time
    template<typename __T>
    void PSort(unsigned N, __T*X, __T*Y, bool warn=1);
    /// parallel radix sort of any type with a SortTraits<>
    /// \param[in]      N  # elements
    /// \param[in,out]  X  array, sorted on return
    /// \param[in]   warn  warn about call from within parallel region?
    /// \note radix sort provides stable sorting and costs O(N) time
    template<typename __T>
    void PSort(unsigned N, __T*X, bool warn=1);
    /// parallel radix sort of the lower K bits of any type with a SortTraits<>
    /// \param[in]      N  # elements
    /// \param[in,out]  X  array, sorted on return
    /// \param[in]      Y  auxiliary array of @a N elements
    /// \param[in]   warn  warn about call from within parallel region?
    /// \note radix sort provides stable sorting and costs O(N) time
    template<int K, typename __T>
    void PSortLow(unsigned N, __T*X, __T*Y, bool warn=1);
    /// parallel radix sort of the lower K bits of any type with a SortTraits<>
    /// \param[in]      N  # elements
    /// \param[in,out]  X  array, sorted on return
    /// \param[in]   warn  warn about call from within parallel region?
    /// \note radix sort provides stable sorting and costs O(N) time
    template<int K, typename __T>
    void PSortLow(unsigned N, __T*X, bool warn=1);
    /// radix sort of single-precision floating point numbers in parallel
    //
    /// \warning the speed-up is usually not optimal, but still non-zero.
    //           the reason is presumably, that we have to do count the data
    //           3 times each in 3 separate loops, while in the serial version
    //           we can to this in just one loop.
    //
    /// \param[in]      N  # elements
    /// \param[in,out]  X  array, sorted on return
    /// \param[in]      Y  auxiliary array of @a N @c floats
    /// \param[in]      w  warn if called from within OMP parallel region?
    /// \note radix sort provides stable sorting and costs O(N) time
    /// \note we map to 32-bit unsigned integers, sort them, and then map back.
    ///       the order-preserving map from float to uint32_t is due to Michael
    ///       Herf (see http://www.stereopsis.com/radix.html)
    /// \note If called from within an OMP parallel region with >1 threads, we
    ///       assume that the data refer to a shared array and sort them. If you
    ///       instead want to sort private data, use Sort() instead.
    void PSort(unsigned N, float*X, float*Y, bool w=true);
    /// radix sort of single-precision floating point numbers in parallel
    //
    /// \warning the speed-up is usually not optimal, but still non-zero.
    //           the reason is presumably, that we have to do count the data
    //           3 times each in 3 separate loops, while in the serial version
    //           we can to this in just one loop.
    //
    /// \param[in]      N  # elements
    /// \param[in,out]  X  array, sorted on return
    /// \param[in]      w  warn if called from within OMP parallel region?
    /// \note radix sort provides stable sorting and costs O(N) time
    /// \note we map to 32-bit unsigned integers, sort them, and then map back.
    ///       the order-preserving map from float to uint32_t is due to Michael
    ///       Herf (see http://www.stereopsis.com/radix.html)
    /// \note allocates and de-allocates an auxiliary array
    /// \note If called from within an OMP parallel region with >1 threads, we
    ///       assume that the data refer to a shared array and sort them. If you
    ///       instead want to sort private data, use Sort() instead.
    void PSort(unsigned N, float*X, bool w=true);
    /// radix sort of double-precision floating point numbers in parallel
    //
    /// \warning this is not faster than serial sorting via Sort()
    //           the reason is presumably, that we have to do count the data
    //           6 times each in 6 separate loops, while in the serial version
    //           we can to this in just 2 loops.
    //
    /// \param[in]      N  # elements
    /// \param[in,out]  X  array, sorted on return
    /// \param[in]      Y  auxiliary array of @a N @c doubles
    /// \param[in]      w  warn if called from within OMP parallel region?
    /// \note radix sort provides stable sorting and costs O(N) time
    /// \note we map to 64-bit unsigned integers, sort them, and then map back.
    ///       the order-preserving map from double to uint64_t is based on the
    ///       equivalent map between float and uint32_t due to Michael Herf see
    ///       http://www.stereopsis.com/radix.html
    /// \note If called from within an OMP parallel region with >1 threads, we
    ///       assume that the data refer to a shared array and sort them. If you
    ///       instead want to sort private data, use Sort() instead.
    void PSort(unsigned N, double*X, double*Y, bool w=true);
    /// radix sort of double-precision floating point numbers in parallel
    //
    /// \warning this is not faster than serial sorting via Sort()
    //           the reason is presumably, that we have to do count the data
    //           6 times each in 6 separate loops, while in the serial version
    //           we can to this in just 2 loops.
    //
    /// \param[in]      N  # elements
    /// \param[in,out]  X  array, sorted on return
    /// \param[in]      w  warn if called from within OMP parallel region?
    /// \note radix sort provides stable sorting and costs O(N) time.
    /// \note we map to 64-bit unsigned integers, sort them, and then map back.
    ///       the order-preserving map from double to uint64_t is based on the
    ///       equivalent map between float and uint32_t due to Michael Herf see
    ///       http://www.stereopsis.com/radix.html
    /// \note allocates and de-allocates an auxiliary array
    /// \note If called from within an OMP parallel region with >1 threads, we
    ///       assume that the data refer to a shared array and sort them. If you
    ///       instead want to sort private data, use Sort() instead.
    void PSort(unsigned N, double*X, bool w=true);
    //@}
#endif
    /// radix sort of the first K bits
    /// \param[in]      N  # elements
    /// \param[in,out]  X  array, sorted on return
    /// \note radix sort provides stable sorting and costs O(N) time
    /// \note allocates and de-allocates an auxiliary array
    template<int K, typename __T>
    inline void SortLow(unsigned N, __T*X)
    {
      __T*Y=WDutils_NEW16(__T,N);
      SortLow<K>(N,X,Y);
      WDutils_DEL16(Y);
    }
  } // namespace WDutils::Radix
} // namespace WDutils {
#ifdef WDutilsDevel
# include <radixP.h>
#endif
//
#endif // WDutils_included_radix_h
