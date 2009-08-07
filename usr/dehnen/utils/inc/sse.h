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

#ifndef WDutils_included_meta_h
#  include <meta.h>
#endif
#ifndef WDutils_included_memory_h
#  include <memory.h>
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
    void Ass(float*f, size_t n, float x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] = x; \endcode
    /// \note up to twice as fast than simple code
    void Ass(double*f, size_t n, double x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] = x; \endcode
    /// \note up to 4 times faster than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 4
    void Ass16(float*f, size_t n, float x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] = x; \endcode
    /// \note up to twice as fast than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 2
    void Ass16(double*f, size_t n, double x);

    /// \code for(size_t i=0; i!=n; ++i) f[i] = -f[i]; \endcode
    /// \note up to 4 times faster than simple code
    void Neg(float*f, size_t n);
    /// \code for(size_t i=0; i!=n; ++i) f[i] = -f[i]; \endcode
    /// \note up to twice as fast than simple code
    void Neg(double*f, size_t n);
    /// \code for(size_t i=0; i!=n; ++i) f[i] = -f[i]; \endcode
    /// \note up to 4 times faster than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 4
    void Neg16(float*f, size_t n);
    /// \code for(size_t i=0; i!=n; ++i) f[i] = -f[i]; \endcode
    /// \note up to twice as fast than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 2
    void Neg16(double*f, size_t n);

    /// \code for(size_t i=0; i!=n; ++i) f[i] = 0; \endcode
    /// \note up to 4 times faster than simple code
    inline void Reset(float*f, size_t n) { Ass(f,n,0.f); }
    /// \code for(size_t i=0; i!=n; ++i) f[i] = 0; \endcode
    /// \note up to twice as fast than simple code
    inline void Reset(double*f, size_t n) { Ass(f,n,0.0); }
    /// \code for(size_t i=0; i!=n; ++i) f[i] = 0; \endcode
    /// \note up to 4 times faster than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 4
    inline void Reset16(float*f, size_t n) { Ass16(f,n,0.f); }
    /// \code for(size_t i=0; i!=n; ++i) f[i] = 0; \endcode
    /// \note up to twice as fast than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 2
    inline void Reset16(double*f, size_t n) { Ass16(f,n,0.0); }

    /// \code for(size_t i=0; i!=n; ++i) f[i] += x; \endcode
    /// \note up to 4 times faster than simple code
    void Add(float*f, size_t n, float x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] += x; \endcode
    /// \note up to twice as fast than simple code
    void Add(double*f, size_t n, double x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] += x; \endcode
    /// \note up to 4 times faster than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 4
    void Add16(float*f, size_t n, float x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] += x; \endcode
    /// \note up to twice as fast than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 2
    void Add16(double*f, size_t n, double x);
    /// \code for(size_t i=0; i!=n; ++i) a[i] += b[i]; \endcode
    /// \note up to 4 times faster than simple code
    /// \note assumes that a and b are 16byte aligned and n a multiple of 4
    void Add16(float*a, size_t n, const float*b);
    /// \code for(size_t i=0; i!=n; ++i) a[i] += b[i]; \endcode
    /// \note up to twice as fast than simple code
    /// \note assumes that a and b are 16byte aligned and n a multiple of 2
    void Add16(double*a, size_t n, const double*b);
    /// \code for(size_t i=0; i!=n; ++i) a[i] += w*b[i]; \endcode
    /// \note up to 4 times faster than simple code
    /// \note assumes that a and b are 16byte aligned and n a multiple of 4
    void Add16(float*a, size_t n, float w, const float*b);
    /// \code for(size_t i=0; i!=n; ++i) a[i] += w*b[i]; \endcode
    /// \note up to twice as fast than simple code
    /// \note assumes that a and b are 16byte aligned and n a multiple of 2
    void Add16(double*a, size_t n, double w, const double*b);

    /// \code for(size_t i=0; i!=n; ++i) f[i] -= x; \endcode
    /// \note up to 4 times faster than simple code
    void Sub(float*f, size_t n, float x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] -= x; \endcode
    /// \note up to twice as fast than simple code
    void Sub(double*f, size_t n, double x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] -= x; \endcode
    /// \note up to 4 times faster than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 4
    void Sub16(float*f, size_t n, float x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] -= x; \endcode
    /// \note up to twice as fast than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 2
    void Sub16(double*f, size_t n, double x);
    /// \code for(size_t i=0; i!=n; ++i) a[i] -= b[i]; \endcode
    /// \note up to 4 times faster than simple code
    /// \note assumes that a and b are 16byte aligned and n a multiple of 4
    void Sub16(float*a, size_t n, const float*b);
    /// \code for(size_t i=0; i!=n; ++i) a[i] -= b[i]; \endcode
    /// \note up to twice as fast than simple code
    /// \note assumes that a and b are 16byte aligned and n a multiple of 2
    void Sub16(double*a, size_t n, const double*b);
    /// \code for(size_t i=0; i!=n; ++i) a[i] -= w*b[i]; \endcode
    /// \note up to 4 times faster than simple code
    /// \note assumes that a and b are 16byte aligned and n a multiple of 4
    void Sub16(float*a, size_t n, float w, const float*b);
    /// \code for(size_t i=0; i!=n; ++i) a[i] -= w*b[i]; \endcode
    /// \note up to twice as fast than simple code
    /// \note assumes that a and b are 16byte aligned and n a multiple of 2
    void Sub16(double*a, size_t n, double w, const double*b);

    /// \code for(size_t i=0; i!=n; ++i) f[i] *= x; \endcode
    /// \note up to 4 times faster than simple code
    void Mul(float*f, size_t n, float x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] *= x; \endcode
    /// \note up to twice as fast than simple code
    void Mul(double*f, size_t n, double x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] *= x; \endcode
    /// \note up to 4 times faster than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 4
    void Mul16(float*f, size_t n, float x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] *= x; \endcode
    /// \note up to twice as fast than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 2
    void Mul16(double*f, size_t n, double x);
    /// \code for(size_t i=0; i!=n; ++i) a[i] *= b[i]; \endcode
    /// \note up to 4 times faster than simple code
    /// \note assumes that a and b are 16byte aligned and n a multiple of 4
    void Mul16(float*a, size_t n, const float*b);
    /// \code for(size_t i=0; i!=n; ++i) a[i] *= b[i]; \endcode
    /// \note up to twice as fast than simple code
    /// \note assumes that a and b are 16byte aligned and n a multiple of 2
    void Mul16(double*a, size_t n, const double*b);

    /// sum of all elements
    /// \code float S(0); for(size_t i=0; i!=n; ++i) S+=f[i]; return S;
    /// \endcode
    /// \note up to 4 times faster than simple code
    float Sum(const float*f, size_t n);
    /// sum of all elements
    /// \code double S(0); for(size_t i=0; i!=n; ++i) S+=f[i]; return S;
    /// \note up to twice as fast than simple code
    /// \endcode
    double Sum(const double*f, size_t n);
    /// sum of all elements
    /// \code float S(0); for(size_t i=0; i!=n; ++i) S+=f[i]; return S;
    /// \endcode
    /// \note up to 4 times faster than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 4
    float Sum16(const float*f, size_t n);
    /// sum of all elements
    /// \code double S(0); for(size_t i=0; i!=n; ++i) S+=f[i]; return S;
    /// \endcode
    /// \note up to twice as fast than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 2
    double Sum16(const double*f, size_t n);

    /// dot of all elements
    /// \code float S(0); for(size_t i=0; i!=n; ++i) S+=a[i]*b[i]; return S;
    /// \endcode
    /// \note up to 4 times faster than simple code
    /// \note assumes that a and b are 16byte aligned and n a multiple of 4
    float Dot16(const float*a, size_t n, const float*b);
    /// dot of all elements
    /// \code double S(0); for(size_t i=0; i!=n; ++i) S+=a[i]*b[i]; return S;
    /// \endcode
    /// \note up to twice as fast than simple code
    /// \note assumes that a and b are 16byte aligned and n a multiple of 2
    double Dot16(const double*a, size_t n, const double*b);

    /// \code for(size_t i=0; i!=n; ++i) f[i] /= x; \endcode
    /// \note up to 4 times faster than simple code
    inline void Div(float*f, size_t n, float x) WDutils_THROWING { 
      if(x==0.f) WDutils_THROW("SSE::Div() by 0\n");
      Mul(f,n,1.f/x);
    }
    /// \code for(size_t i=0; i!=n; ++i) f[i] /= x; \endcode
    /// \note up to twice as fast than simple code
    inline void Div(double*f, size_t n, double x) WDutils_THROWING { 
      if(x==0.0) WDutils_THROW("SSE::Div() by 0\n");
      Mul(f,n,1.0/x);
    }
    /// \code for(size_t i=0; i!=n; ++i) f[i] /= x; \endcode
    /// \note up to 4 times faster than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 4
    inline void Div16(float*f, size_t n, float x) WDutils_THROWING { 
      if(x==0.f) WDutils_THROW("SSE::Div16() by 0\n");
      Mul16(f,n,1.f/x);
    }
    /// \code for(size_t i=0; i!=n; ++i) f[i] /= x; \endcode
    /// \note up to twice as fast than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 2
    inline void Div16(double*f, size_t n, double x) WDutils_THROWING { 
      if(x==0.0) WDutils_THROW("SSE::Div16() by 0\n");
      Mul16(f,n,1.0/x);
    }
    /// \code for(size_t i=0; i!=n; ++i) a[i] /= b[i]; \endcode
    /// \note up to 4 times faster than simple code
    /// \note assumes that a and b are 16byte aligned and n a multiple of 4
    void Div16(float*a, size_t n, const float*b);
    /// \code for(size_t i=0; i!=n; ++i) a[i] /= b[i]; \endcode
    /// \note up to twice as fast than simple code
    /// \note assumes that a and b are 16byte aligned and n a multiple of 2
    void Div16(double*a, size_t n, const double*b);

    /// \code for(size_t i=0; i!=n; ++i) f[i] = x/f[i]; \endcode
    /// \note up to 4 times faster than simple code
    void Inv(float*f, size_t n, float x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] = x/f[i]; \endcode
    /// \note up to twice as fast than simple code
    void Inv(double*f, size_t n, double x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] = x/f[i]; \endcode
    /// \note up to 4 times faster than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 4
    void Inv16(float*f, size_t n, float x);
    /// \code for(size_t i=0; i!=n; ++i) f[i] = x/f[i]; \endcode
    /// \note up to twice as fast than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 2
    void Inv16(double*f, size_t n, double x);
    /// \code for(size_t i=0; i!=n; ++i) a[i] = x/b[i]; \endcode
    /// \note up to 4 times faster than simple code
    /// \note assumes that a and b are 16byte aligned and n a multiple of 4
    void Inv16(float*a, size_t n, float x, const float*b);
    /// \code for(size_t i=0; i!=n; ++i) a[i] = x/b[i]; \endcode
    /// \note up to twice as fast than simple code
    /// \note assumes that a and b are 16byte aligned and n a multiple of 2
    void Inv16(double*a, size_t n, double x, const double*b);

    /// \code for(size_t i=0; i!=n; ++i) f[i] = 1/f[i]; \endcode
    /// \note up to 4 times faster than simple code
    inline void Reciprocal(float*f, size_t n) { Inv(f,n,1.f); }
    /// \code for(size_t i=0; i!=n; ++i) f[i] = 1/f[i]; \endcode
    /// \note up to twice as fast than simple code
    inline void Reciprocal(double*f, size_t n) { Inv(f,n,1.0); }
    /// \code for(size_t i=0; i!=n; ++i) f[i] = 1/f[i]; \endcode
    /// \note up to 4 times faster than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 4
    inline void Reciprocal16(float*f, size_t n) { Inv16(f,n,1.f); }
    /// \code for(size_t i=0; i!=n; ++i) f[i] = 1/f[i]; \endcode
    /// \note up to twice as fast than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 2
    inline void Reciprocal16(double*f, size_t n) { Inv16(f,n,1.0); }
    /// \code for(size_t i=0; i!=n; ++i) a[i] = 1/b[i]; \endcode
    /// \note up to 4 times faster than simple code
    /// \note assumes that a and b are 16byte aligned and n a multiple of 4
    inline void Reciprocal16(float*a, size_t n, const float*b)
    { Inv16(a,n,1.f,b); }
    /// \code for(size_t i=0; i!=n; ++i) a[i] = 1/b[i]; \endcode
    /// \note up to twice as fast than simple code
    /// \note assumes that a and b are 16byte aligned and n a multiple of 2
    inline void Reciprocal16(double*a, size_t n, const double*b)
    { Inv16(a,n,1.0,b); }

    /// \code for(size_t i=0; i!=n; ++i) f[i] = std::sqrt(f[i]); \endcode
    /// \note up to 4 times faster than simple code
    void Sqrt(float*f, size_t n);
    /// \code for(size_t i=0; i!=n; ++i) f[i] = std::sqrt(f[i]); \endcode
    /// \note up to twice as fast than simple code
    void Sqrt(double*f, size_t n);
    /// \code for(size_t i=0; i!=n; ++i) f[i] = std::sqrt(f[i]); \endcode
    /// \note up to 4 times faster than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 4
    void Sqrt16(float*f, size_t n);
    /// \code for(size_t i=0; i!=n; ++i) f[i] = std::sqrt(f[i]); \endcode
    /// \note up to twice as fast than simple code
    /// \note assumes that f is 16byte aligned and n a multiple of 2
    void Sqrt16(double*f, size_t n);
    /// \code for(size_t i=0; i!=n; ++i) a[i] = std::sqrt(b[i]); \endcode
    /// \note up to 4 times faster than simple code
    /// \note assumes that a and b are 16byte aligned and n a multiple of 4
    void Sqrt16(float*a, size_t n, const float*b);
    /// \code for(size_t i=0; i!=n; ++i) a[i] = std::sqrt(b[i]); \endcode
    /// \note up to twice as fast than simple code
    /// \note assumes that a and b are 16byte aligned and n a multiple of 2
    void Sqrt16(double*a, size_t n, const double*b);
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

    ////////////////////////////////////////////////////////////////////////////
    /// An array of float or double supporting operations via SSE instructions
    template<typename _F>
    class Array16 {
      WDutilsStaticAssert( meta::TypeInfo<_F>::is_floating_point );
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
      { for(size_t i=_N; i!=_S; ++i) const_cast<_F&>(_A[i]) = _F(0); }
      /// check for size mismatch
      void check_size(Array16 const&B, const char*name) const WDutils_THROWING
      {
	if(B._N != _N)
	  WDutils_THROW("SSE::Array16<%s>::%s: size mismatch:%lu vs %lu\n",
			nameof(_F),name,_N,B._N);
      }
    public:
      /// default ctor
      Array16()
	: _S(0), _N(0), _A(0) {}
      /// ctor from given size
      explicit Array16(size_t n)
	: _S(Top<_F>(n)), _N(n), _A(new16<_F>(_S)) {}
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
      _F &operator[] (int i) { return _A[i]; }
      /// \name unary operations
      //@{
      /// reset element-wise to 0
      /// \code for(int i=0; i!=size(); ++i) A[i] = 0; \endcode
      Array16&reset()
      {
	Reset16(_A,_S);
	return*this;
      }
      /// negate element-wise
      /// \code for(int i=0; i!=size(); ++i) A[i] = -A[i]; \endcode
      Array16&negate()
      {
	Neg16(_A,_S);
	return*this;
      }
      //@}
      /// \name binary operations with scalar
      //@{
      /// assign element-wise to scalar
      /// \code for(int i=0; i!=size(); ++i) A[i] = x; \endcode
      Array16&operator=(_F x)
      {
	Ass16(_A,_S,x);
	return*this;
      }
      /// add a scalar to each element
      /// \code for(int i=0; i!=size(); ++i) A[i] += x; \endcode
      Array16&operator+=(_F x)
      {
	Add16(_A,_S,x);
	return*this;
      }
      /// subtract a scalar from each element
      /// \code for(int i=0; i!=size(); ++i) A[i] -= x; \endcode
      Array16&operator-=(_F x)
      {
	Sub16(_A,_S,x);
	return*this;
      }
      /// multiply each element by a scalar
      /// \code for(int i=0; i!=size(); ++i) A[i] *= x; \endcode
      Array16&operator*=(_F x)
      {
	Mul16(_A,_S,x);
	return*this;
      }
      /// divide each element by a scalar
      /// \code for(int i=0; i!=size(); ++i) A[i] *= x; \endcode
      Array16&operator/=(_F x) WDutils_THROWING
      {
	Div16(_A,_S,x);
	return*this;
      }
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
	Add16(_A,_S,B._A);
	return*this;
      }
      /// subtract element-wise
      /// \code for(int i=0; i!=size(); ++i) A[i] += B[i]; \endcode
      /// \note it is an error if the number of elements do not match
      Array16&operator-=(Array16 const&B) WDutils_THROWING
      {
	check_size(B,"operator-=(Array16&)");
	Sub16(_A,_S,B._A);
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
	return Dot16(_A,_S,B._A);
      }
      //@}
      /// \name tertiary operations with scalar and another Array16
      //@{
      /// add weighted element-wise
      /// \code for(int i=0; i!=size(); ++i) A[i] += w*B[i]; \endcode
      Array16&addtimes(Array16 const&B, _F w) WDutils_THROWING
      {
	check_size(B,"addtimes(Array16&)");
	Add16(_A,_S,w,B._A);
	return*this;
      }
      /// subtract weighted element-wise
      /// \code for(int i=0; i!=size(); ++i) A[i] += w*B[i]; \endcode
      Array16&subtimes(Array16 const&B, _F w) WDutils_THROWING
      {
	check_size(B,"subtimes(Array16&)");
	Sub16(_A,_S,w,B._A);
	return*this;
      }
      //@}
    };

  } // namespace SSE
} // namespace WDutils
//
#endif
