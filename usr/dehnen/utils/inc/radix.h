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

////////////////////////////////////////////////////////////////////////////////
namespace WDutils {
  template<int, typename traits> struct RadixSortBits;
  /// stable radix sort of the lower 8 of 8 or more bits
  template<typename traits> struct RadixSortBits<8,traits>
  {
    typedef typename traits::input_type   input_type;
    typedef typename traits::integer_type integer_type;
    WDutilsStaticAssert(std::numeric_limits<integer_type>::is_integer);
    static void sort(unsigned N, input_type*X, input_type*Y)
    {
      uint32 B0[256] = {0};
      static const integer_type m=0xff;
      for(uint32 i=0; i!=N; ++i) {
	register integer_type f=traits::integer(Y[i]=X[i]);
	B0[f&m]++;
      }
      for(uint32 i=0,s=0,t; i!=256; ++i) { t=B0[i]; B0[i]=s; s+=t; }
      for(uint32 i=0; i!=N; ++i) X[B0[traits::integer(Y[i])&m]++]=Y[i];
    }
  };
  /// stable radix sort of the lower 16 of 16 or more bits
  template<typename traits> struct RadixSortBits<16,traits>
  {
    typedef typename traits::input_type   input_type;
    typedef typename traits::integer_type integer_type;
    WDutilsStaticAssert(std::numeric_limits<integer_type>::is_integer &&
			sizeof(integer_type)>=2);
    static void sort(unsigned N, input_type*X, input_type*Y)
    {
      uint32 B0[512]={0},*B1=B0+256;
      static const integer_type m=0xff;
      for(uint32 i=0; i!=N; ++i) {
	register integer_type f=traits::integer(X[i]);
	B0[ f     &m] ++;
	B1[(f>>=8)&m] ++;
      }
      for(uint32 i=0,s=0,t; i!=256; ++i) { t=B0[i]; B0[i]=s; s+=t; }
      for(uint32 i=0,s=0,t; i!=256; ++i) { t=B1[i]; B1[i]=s; s+=t; }
      for(uint32 i=0; i!=N; ++i) Y[B0[traits::integer(X[i])   &m]++]=X[i];
      for(uint32 i=0; i!=N; ++i) X[B1[traits::integer(Y[i])>>8&m]++]=Y[i];
    }
  };
  /// stable radix sort of the lower 32 of 32 or more bits
  template<typename traits> struct RadixSortBits<32,traits>
  {
    typedef typename traits::input_type   input_type;
    typedef typename traits::integer_type integer_type;
    WDutilsStaticAssert(std::numeric_limits<integer_type>::is_integer &&
			sizeof(integer_type)>=4);
    static void sort(unsigned N, input_type*X, input_type*Y)
    {
      static const uint32  n0=0x800,n2=0x400;
      static const integer_type m0=0x7ff,m2=0x3ff;
      uint32 B0[n0+n0+n2]={0},*B1=B0+n0,*B2=B1+n0;
      for(uint32 i=0; i!=N; ++i) {
	register integer_type f=traits::integer(Y[i]=X[i]);
	B0[ f      &m0] ++;
	B1[(f>>=11)&m0] ++;
	B2[(f>>=11)&m2] ++;
      }
      for(uint32 i=0,s=0,t; i!=n0; ++i) { t=B0[i]; B0[i]=s; s+=t; }
      for(uint32 i=0,s=0,t; i!=n0; ++i) { t=B1[i]; B1[i]=s; s+=t; }
      for(uint32 i=0,s=0,t; i!=n2; ++i) { t=B2[i]; B2[i]=s; s+=t; }
      for(uint32 i=0; i!=N; ++i) X[B0[traits::integer(Y[i])     &m0]++]=Y[i];
      for(uint32 i=0; i!=N; ++i) Y[B1[traits::integer(X[i])>>11 &m0]++]=X[i];
      for(uint32 i=0; i!=N; ++i) X[B2[traits::integer(Y[i])>>22 &m2]++]=Y[i];
    }
  };
  /// stable radix sort of the lower 64 of 64 or more bits
  template<typename traits> struct RadixSortBits<64,traits>
  {
    typedef typename traits::input_type   input_type;
    typedef typename traits::integer_type integer_type;
    WDutilsStaticAssert(std::numeric_limits<integer_type>::is_integer &&
			sizeof(integer_type)>=8);
    static void sort(unsigned N, input_type*X, input_type*Y)
    {
      static const uint32  n0=0x800,n2=0x400,n=n0+n0+n2;
      static const integer_type m0=0x7ff,m2=0x3ff;
      uint32 B0[n]={0},*B1=B0+n0,*B2=B1+n0;
      // sort the lower 32 bits
      for(uint32 i=0; i!=N; ++i) {
	register uint64 f=traits::integer(X[i]);
	B0[ f      &m0] ++;
	B1[(f>>=11)&m0] ++;
	B2[(f>>=11)&m2] ++;
      }
      for(uint32 i=0,s=0,t; i!=n0; ++i) { t=B0[i]; B0[i]=s; s+=t; }
      for(uint32 i=0,s=0,t; i!=n0; ++i) { t=B1[i]; B1[i]=s; s+=t; }
      for(uint32 i=0,s=0,t; i!=n2; ++i) { t=B2[i]; B2[i]=s; s+=t; }
      for(uint32 i=0; i!=N; ++i) Y[B0[traits::integer(X[i])     &m0]++]=X[i];
      for(uint32 i=0; i!=N; ++i) X[B1[traits::integer(Y[i])>>11 &m0]++]=Y[i];
      for(uint32 i=0; i!=N; ++i) Y[B2[traits::integer(X[i])>>22 &m2]++]=X[i];
      // sort the upper 32 bits
      for(uint32 i=0; i!=n; ++i) B0[i]=0;
      for(uint32 i=0; i!=N; ++i) {
	register uint64 f=traits::integer(Y[i]);
	B0[(f>>=32)&m0] ++;
	B1[(f>>=11)&m0] ++;
	B2[(f>>=11)&m2] ++;
      }
      for(uint32 i=0,s=0,t; i!=n0; ++i) { t=B0[i]; B0[i]=s; s+=t; }
      for(uint32 i=0,s=0,t; i!=n0; ++i) { t=B1[i]; B1[i]=s; s+=t; }
      for(uint32 i=0,s=0,t; i!=n2; ++i) { t=B2[i]; B2[i]=s; s+=t; }
      for(uint32 i=0; i!=N; ++i) X[B0[traits::integer(Y[i])>>32 &m0]++]=Y[i];
      for(uint32 i=0; i!=N; ++i) Y[B1[traits::integer(X[i])>>43 &m0]++]=X[i];
      for(uint32 i=0; i!=N; ++i) X[B2[traits::integer(Y[i])>>54 &m2]++]=Y[i];
    }
  };

  /// traits for use in @c RadixSortBits
  template<typename> struct RadixSortTraits;

  /// for implementing @c RadixSortTraits of integer types
  template<typename __I> struct RadixSortTraitsInteger
  {
    WDutilsStaticAssert(std::numeric_limits<__I>::is_integer);
    typedef __I input_type;
    typedef __I integer_type;
    static integer_type integer(input_type x) { return x; }
  };

#define WDutilsRadixSortTraits(INT)		\
  template<> struct RadixSortTraits<INT>	\
    : public RadixSortTraitsInteger<INT> {}
  WDutilsRadixSortTraits(int8);
  WDutilsRadixSortTraits(int16);
  WDutilsRadixSortTraits(int32);
  WDutilsRadixSortTraits(int64);
  WDutilsRadixSortTraits(uint8);
  WDutilsRadixSortTraits(uint16);
  WDutilsRadixSortTraits(uint32);
  WDutilsRadixSortTraits(uint64);
#undef WDutilsRadixSortTraits
  /// RadixSortTraits for floats
  /// \note not used in RadixSort(float) below, but useful for implementing
  ///       RadixSort of structs with a float member to be sorted on
  /// \note the order-preserving map from float to uint32 is due to Michael
  ///       Herf (see http://www.stereopsis.com/radix.html)
  template<> struct RadixSortTraits<float>
  {
    typedef float  input_type;
    typedef uint32 integer_type;
    static uint32 forward(uint32 f)
    { return f ^ (-int32(f>>31) | 0x80000000); }
    static integer_type integer(input_type const&x)
    { return forward(reinterpret_cast<uint32 const&>(x)); }
  };
  /// RadixSortTraits for doubles
  /// \note not used in RadixSort(double) below, but useful for implementing
  ///       RadixSort of structs with a double member to be sorted on
  /// \note the order-preserving map from double to uint64 is based on the
  ///       equivalent map between float and uint32 due to Michael Herf see
  ///       http://www.stereopsis.com/radix.html
  template<> struct RadixSortTraits<double>
  {
    typedef double input_type;
    typedef uint64 integer_type;
    static uint64 forward(uint64 f)
    { return f ^ (-int64(f>>63) | 0x8000000000000000ll); }
    static uint64 integer(double const& x)
    { return forward(reinterpret_cast<uint64 const&>(x)); }
  };

  /// radix sort of any type with a RadixSortTraits<>
  /// \param[in]      N # elements
  /// \param[in,out]  X array, sorted on return
  /// \param[in]      Y auxiliary array of @a N elements
  /// \note radix sort provides stable sorting and costs O(N) time
  template<typename __T>
  inline void RadixSort(unsigned N, __T*X, __T*Y)
  { RadixSortBits<sizeof(__T),RadixSortTraits<__T> >::sort(N,X,Y); }

  /// radix sort of the lower K bits of any type with a RadixSortTraits<>
  /// \param[in]      N # elements
  /// \param[in,out]  X array, sorted on return
  /// \param[in]      Y auxiliary array of @a N elements
  /// \note radix sort provides stable sorting and costs O(N) time
  template<int K, typename __T>
  inline void RadixSortLow(unsigned N, __T*X, __T*Y)
  { RadixSortBits<K,RadixSortTraits<__T> >::sort(N,X,Y); }
  
  /// radix sort of single-precision floating point numbers
  /// \param[in]      N # elements
  /// \param[in,out]  X array, sorted on return
  /// \param[in]      Y auxiliary array of @a N @c floats
  /// \note radix sort provides stable sorting and costs O(N) time
  /// \note we map to 32-bit unsigned integers, sort them, and then map back
  /// \note the order-preserving map from float to uint32 is due to Michael
  ///       Herf (see http://www.stereopsis.com/radix.html)
  /// \note declared inline here rather than non-inline in some file radix.cc,
  ///       for otherwise the g++ compiler/linker produces code which cannot
  ///       resolve the related symbol (with dynamic linking).
  inline void RadixSort(unsigned N, float*X, float*Y)
  {
    static const uint32 n0=0x800,n2=0x400;
    static const uint32 m0=0x7ff,m2=0x3ff;
    uint32 B0[n0+n0+n2]={0},*B1=B0+n0,*B2=B1+n0;
    uint32*iX=reinterpret_cast<uint32*>(X);
    uint32*iY=reinterpret_cast<uint32*>(Y);
    for(uint32 i=0; i!=N; ++i) {
      register uint32 f=iX[i];
      iY[i]=f^=-int32(f>>31)|0x80000000;
      B0[f       &m0]++;
      B1[(f>>=11)&m0]++;
      B2[(f>>=11)&m2]++;
    }
    for(uint32 i=0,s=0,t; i!=n0; ++i) { t=B0[i]; B0[i]=s; s+=t; }
    for(uint32 i=0,s=0,t; i!=n0; ++i) { t=B1[i]; B1[i]=s; s+=t; }
    for(uint32 i=0,s=0,t; i!=n2; ++i) { t=B2[i]; B2[i]=s; s+=t; }
    for(uint32 i=0; i!=N; ++i) iX[B0[iY[i]     &m0]++] = iY[i];
    for(uint32 i=0; i!=N; ++i) iY[B1[iX[i]>>11 &m0]++] = iX[i];
    for(uint32 i=0; i!=N; ++i) {
      register uint32 f=iY[i];
      iX[B2[f>>22 &m2]++]=f^(((f>>31)-1)|0x80000000); 
    }
  }
  /// radix sort of double-precision floating point numbers
  /// \param[in]      N # elements
  /// \param[in,out]  X array, sorted on return
  /// \param[in]      Y auxiliary array of @a N @c doubles
  /// \note radix sort provides stable sorting and costs O(N) time
  /// \note we map to 64-bit unsigned integers, sort them, and then map back
  /// \note the order-preserving map from double to uint64 is based on the
  ///       equivalent map between float and uint32 due to Michael Herf see
  ///       http://www.stereopsis.com/radix.html
  /// \note declared inline here rather than non-inline in some file radix.cc,
  ///       for otherwise the g++ compiler/linker produces code which cannot
  ///       resolve the related symbol (with dynamic linking).
  inline void RadixSort(unsigned N, double*X, double*Y)
  {
    static const uint32 n0=0x800,n2=0x400, n=n0+n0+n2;
    static const uint64 m0=0x7ff,m2=0x3ff;
    uint32 B0[n]={0},*B1=B0+n0,*B2=B1+n0;
    uint64*iX=reinterpret_cast<uint64*>(X);
    uint64*iY=reinterpret_cast<uint64*>(Y);
    // sort the lower 32 bits
    for(uint32 i=0; i!=N; ++i) {
      register uint64 f=iX[i];
      iX[i]=f^=-int64(f>>63) | 0x8000000000000000ll;
      B0[f       &m0] ++;
      B1[(f>>=11)&m0] ++;
      B2[(f>>=11)&m2] ++;
    }
    for(uint32 i=0,s=0,t; i!=n0; ++i) { t=B0[i]; B0[i]=s; s+=t; }
    for(uint32 i=0,s=0,t; i!=n0; ++i) { t=B1[i]; B1[i]=s; s+=t; }
    for(uint32 i=0,s=0,t; i!=n2; ++i) { t=B2[i]; B2[i]=s; s+=t; }
    for(uint32 i=0; i!=N; ++i) iY[B0[iX[i]     &m0]++] = iX[i];
    for(uint32 i=0; i!=N; ++i) iX[B1[iY[i]>>11 &m0]++] = iY[i];
    for(uint32 i=0; i!=N; ++i) iY[B2[iX[i]>>22 &m2]++] = iX[i];
    // sort the upper 32 bits
    for(uint32 i=0; i!=n; ++i) B0[i]=0;
    for(uint32 i=0; i!=N; ++i) {
      register uint64 f=iY[i];
      B0[(f>>=32)&m0] ++;
      B1[(f>>=11)&m0] ++;
      B2[(f>>=11)&m2] ++;
    }
    for(uint32 i=0,s=0,t; i!=n0; ++i) { t=B0[i]; B0[i]=s; s+=t; }
    for(uint32 i=0,s=0,t; i!=n0; ++i) { t=B1[i]; B1[i]=s; s+=t; }
    for(uint32 i=0,s=0,t; i!=n2; ++i) { t=B2[i]; B2[i]=s; s+=t; }
    for(uint32 i=0; i!=N; ++i) iX[B0[iY[i]>>32 &m0]++] = iY[i];
    for(uint32 i=0; i!=N; ++i) iY[B1[iX[i]>>43 &m0]++] = iX[i];
    for(uint32 i=0; i!=N; ++i) {
      register uint64 f = iY[i];
      iX[B2[f    >>54 &m2]++] =  f ^ (((f>>63) - 1) | 0x8000000000000000ll);
    }
  }
  /// radix sort
  /// \param[in]      N # elements
  /// \param[in,out]  X array, sorted on return
  /// \note radix sort provides stable sorting and costs O(N) time
  /// \note allocates and de-allocates an auxiliary array
  template<typename __T>
  inline void RadixSort(unsigned N, __T*X)
  {
    __T*Y=WDutils_NEW16(__T,N);
    RadixSort(N,X,Y);
    WDutils_DEL16(Y);
  }
  /// radix sort of the first K bits
  /// \param[in]      N # elements
  /// \param[in,out]  X array, sorted on return
  /// \note radix sort provides stable sorting and costs O(N) time
  /// \note allocates and de-allocates an auxiliary array
  template<int K, typename __T>
  inline void RadixSortLow(unsigned N, __T*X)
  {
    __T*Y=WDutils_NEW16(__T,N);
    RadixSortLow<K>(N,X,Y);
    WDutils_DEL16(Y);
  }
} // namespace WDutils {
#endif // WDutils_included_radix_h
