// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/sse.h
///
/// \brief  support for SSE coding and simple SSE supported code
///
/// \author Walter Dehnen
///
/// \date   2009
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Walter Dehnen
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
#    include <xmmintrin.h>
#  endif

# ifdef __SSE2__
#  ifndef WDutils_included_emmintrin_h
#    define WDutils_included_emmintrin_h
#    include <emmintrin.h>
#  endif
# endif
#endif

#ifndef WDutils_included_exception_h
#  include <exception.h>
#endif

namespace WDutils {

  /// support for coding with SSE intrinsics, code using SSE intrinsics.
  namespace SSE {

    template<int N> struct __Aux {
      template<typename __I> static __I top(__I i) {
	__I m = i % N;
	return m? i+N-m : i;
      }
      template<typename __I> static __I bot(__I i) {
	return i - (i % N);
      }
    };
    template<> struct __Aux<2> {
      template<typename __I>
      static __I top(__I i) { return i&1? i+1:i; }
      template<typename __I>
      static __I bot(__I i) { return i&1? i-1:i; }
    };
    template<> struct __Aux<4> {
      template<typename __I>
      static __I top(__I i) {
	__I m = i & 3;
	return m? i+4-m : i;
      }
      template<typename __I> static __I bot(__I i) {
	return i & ~3;
      }
    };
    template<> struct __Aux<8> {
      template<typename __I>
      static __I top(__I i) {
	__I m = i & 7;
	return m? i+8-m : i;
      }
      template<typename __I> static __I bot(__I i) {
	return i & ~7;
      }
    };
    template<> struct __Aux<16> {
      template<typename __I>
      static __I top(__I i) {
	__I m = i & 15;
	return m? i+16-m : i;
      }
      template<typename __I> static __I bot(__I i) {
	return i & ~15;
      }
    };

    /// smallest multiple of N not less than i
    template<int N, typename __I> inline
    __I top(__I i) { return __Aux<N>::top(i); }
    /// largest multiple of N not greater than i
    template<int N, typename __I> inline
    __I bottom(__I i) { return __Aux<N>::bot(i); }
    /// \name simple manipulations applied to each element of an array
    //@{
    /// \code for(size_t i=0; i!=n; ++i) f[i] = x; \endcode
    /// \note up to 4 times faster than simple code
    void Assign(float*f, size_t n, float x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] = x; \endcode
    /// \note up to twice as fast than simple code
    void Assign(double*f, size_t n, double x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] = 0; \endcode
    /// \note up to 4 times faster than simple code
    inline void Reset(float*f, size_t n) { Assign(f,n,0.f); }
    /// \code for(size_t i=0; i!=n; ++i) f[i] = 0; \endcode
    /// \note up to twice as fast than simple code
    inline void Reset(double*f, size_t n) { Assign(f,n,0.0); }
    /// \code for(size_t i=0; i!=n; ++i) f[i] += x; \endcode
    /// \note up to 4 times faster than simple code
    void Add(float*f, size_t n, float x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] += x; \endcode
    /// \note up to twice as fast than simple code
    void Add(double*f, size_t n, double x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] -= x; \endcode
    /// \note up to 4 times faster than simple code
    void Subtract(float*f, size_t n, float x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] -= x; \endcode
    /// \note up to twice as fast than simple code
    void Subtract(double*f, size_t n, double x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] *= x; \endcode
    /// \note up to 4 times faster than simple code
    void Multiply(float*f, size_t n, float x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] *= x; \endcode
    /// \note up to twice as fast than simple code
    void Multiply(double*f, size_t n, double x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] /= x; \endcode
    /// \note up to 4 times faster than simple code
    inline void Divide(float*f, size_t n, float x) WDutils_THROWING { 
      if(x==0.f) WDutils_THROW("SSE::Divide() by 0\n");
      Multiply(f,n,1.f/x);
    }
    /// \code for(size_t i=0; i!=n; ++i) f[i] /= x; \endcode
    /// \note up to twice as fast than simple code
    inline void Divide(double*f, size_t n, double x) WDutils_THROWING { 
      if(x==0.0) WDutils_THROW("SSE::Divide() by 0\n");
      Multiply(f,n,1.0/x);
    }
    /// \code for(size_t i=0; i!=n; ++i) f[i] = x/f[i]; \endcode
    /// \note up to 4 times faster than simple code
    void Invert(float*f, size_t n, float x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] = x/f[i]; \endcode
    /// \note up to twice as fast than simple code
    void Invert(double*f, size_t n, double x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] = 1/f[i]; \endcode
    /// \note up to 4 times faster than simple code
    inline void Reciprocal(float*f, size_t n) { Invert(f,n,1.f); }
    /// \code for(size_t i=0; i!=n; ++i) f[i] = 1/f[i]; \endcode
    /// \note up to twice as fast than simple code
    inline void Reciprocal(double*f, size_t n) { Invert(f,n,1.0); }
    /// \code for(size_t i=0; i!=n; ++i) f[i] = std::sqrt(f[i]); \endcode
    /// \note up to 4 times faster than simple code
    void SquareRoot(float*f, size_t n);
    /// \code for(size_t i=0; i!=n; ++i) f[i] = std::sqrt(f[i]); \endcode
    /// \note up to twice as fast than simple code
    void SquareRoot(double*f, size_t n);
    //@}

    /// contains some constants and code relevant for SSE coding
    /// instantinations for float and double
    template<typename __F> struct Traits;

    /// SSE::Traits<float>
    template<> struct Traits<float>
    {
      /// is SSE enabled for this type?
#ifdef __SSE__
      const static bool sse = true;
#else
      const static bool sse = false;
#endif
      /// alignment number: K floats align to 128 bytes
      const static int K=4;
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

    /// SSE::Traits<double>
    template<> struct Traits<double>
    {
      /// is SSE enabled for this type?
#ifdef __SSE2__
      const static bool sse = true;
#else
      const static bool sse = false;
#endif
      /// alignment number: K doubles align to 128 bytes
      const static int K=2;
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
  } // namespace SSE  
} // namespace WDutils
//
#endif
