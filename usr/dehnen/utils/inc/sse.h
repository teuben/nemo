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
#ifndef WDutils_included_cstring
#  define WDutils_included_cstring
#  include <cstring>
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

    // auxiliary for routines below
    template<typename T> struct U16
    {
      WDutilsStaticAssert((// meta::TypeCompare<T,int>::identical    ||
			   meta::TypeCompare<T,float>::identical  ||
			   meta::TypeCompare<T,double>::identical ));
      static void Ass(T*, size_t, T);
      static void Neg(T*, size_t);
      static void Add(T*, size_t, T);
      static void Sub(T*, size_t, T);
      static void Mul(T*, size_t, T);
      static void Div(T*f, size_t n, T x) WDutils_THROWING
      { 
// 	WDutilsStaticAssert(meta::TypeInfo<T>::is_floating_point);
	if(x==T(0)) WDutils_THROW("SSE::Div() by 0\n");
	Mul(f,n,T(1)/x);
      }
      static void Inv(T*, size_t n, T x);
      static void Sqrt(T*f, size_t n);
      static T Sum(const T*, size_t);
    };
    // special case: T=int not all implemented (yet?)
    template<> struct U16<int>
    {
      typedef int T;
      static void Ass(T*, size_t, T);
//       static void Neg(T*, size_t);
      static void Add(T*, size_t, T);
//       static void Sub(T*, size_t, T);
//       static void Mul(T*, size_t, T);
//       static T Sum(const T*, size_t);
    };


    /// \name simple manipulations of each array element
    /// \note up to about 16/sizeof(T) times faster than simple code
    //@{

    /// assign same value to all elements
    /// \code for(size_t i=0; i!=n; ++i) f[i] = x; \endcode
    template<typename T>
    inline void Ass(int*f, size_t n, int x)
    { U16<T>::Ass(f,n,x); }

    /// assign element-wise to another array
    /// \code for(size_t i=0; i!=n; ++i) a[i] = b[i]; \endcode
    template<typename T>
    inline void Ass(T*a, size_t n, const T*b)
    { std::memcpy(a,b,n*sizeof(T)); }

    /// negate all elements
    /// \code for(size_t i=0; i!=n; ++i) f[i] = -f[i]; \endcode
    template<typename T>
    inline void Neg(T*f, size_t n)
    { U16<T>::Neg(f,n); }

    /// set eaqch element to zero
    /// \code for(size_t i=0; i!=n; ++i) f[i] = 0; \endcode
    template<typename T>
    inline void Reset(T*f, size_t n)
    { Ass(f,n,T(0)); }

    /// add same value to each elements
    /// \code for(size_t i=0; i!=n; ++i) f[i] += x; \endcode
    template<typename T>
    inline void Add(T*f, size_t n, T x)
    { U16<T>::Add(f,n,x); }

    /// subtract same value from each element
    /// \code for(size_t i=0; i!=n; ++i) f[i] -= x; \endcode
    template<typename T>
    inline void Sub(T*f, size_t n, T x)
    { U16<T>::Sub(f,n,x); }

    /// multiply each element with same value
    /// \code for(size_t i=0; i!=n; ++i) f[i] *= x; \endcode
    /// \note up to 4 times faster than simple code
    template<typename T>
    inline void Mul(T*f, size_t n, T x)
    { U16<T>::Mul(f,n,x); }

    /// sum of all elements
    /// \code float S(0); for(size_t i=0; i!=n; ++i) S+=f[i]; return S;
    /// \endcode
    template<typename T>
    inline T Sum(const T*f, size_t n)
    { return U16<T>::Sum(f,n); }

    /// divide each element by same value
    /// \code for(size_t i=0; i!=n; ++i) f[i] /= x; \endcode
    template<typename T>
    inline void Div(T*f, size_t n, T x)
    { U16<T>::Div(f,n,x); }

    /// replace each elements by its inverse times a constant
    /// \code for(size_t i=0; i!=n; ++i) f[i] = x/f[i]; \endcode
    template<typename T>
    inline void Inv(float*f, size_t n, float x)
    { U16<T>::Inv(f,n,x); }

    /// replace each element by its inverse
    /// \code for(size_t i=0; i!=n; ++i) f[i] = 1/f[i]; \endcode
    template<typename T>
    inline void Reciprocal(T*f, size_t n)
    { Inv(f,n,T(1)); }

    /// replace each element by its square root
    /// \code for(size_t i=0; i!=n; ++i) f[i] = std::sqrt(f[i]); \endcode
    template<typename T>
    inline void Sqrt(float*f, size_t n)
    { U16<T>::Sqrt(f,n); }
    //@}

    // auxiliary for struct Align below
    template<typename T> struct A16
    {
      WDutilsStaticAssert((// meta::TypeCompare<T,int>::identical    ||
			   meta::TypeCompare<T,float>::identical  ||
			   meta::TypeCompare<T,double>::identical));
      static void Ass(T*, size_t, T);
      static void Ass(T*, size_t, T, const T*);
      static void Neg(T*, size_t);
      static void Add(T*, size_t, T);
      static void Add(T*, size_t, const T*);
      static void Add(T*, size_t, T, const T*);
      static void Sub(T*, size_t, T);
      static void Sub(T*, size_t, const T*);
      static void Sub(T*, size_t, T, const T*);
      static void Mul(T*, size_t, T);
      static void Mul(T*, size_t, const T*);
      static void Div(T*f, size_t n, T x) WDutils_THROWING { 
// 	WDutilsStaticAssert(meta::TypeInfo<T>::is_floating_point);
	if(x==T(0)) WDutils_THROW("SSE::Aligned::Div() by 0\n");
	Mul(f,n,T(1)/x);
      }
      static void Div(T*, size_t, const T*);
      static void Inv(T*, size_t, T);
      static void Inv(T*, size_t, T, const T*);
      static void Sqrt(T*, size_t);
      static void Sqrt(T*, size_t, const T*);
      static T Sum(const T*, size_t);
      static T Dot(const T*, size_t, const T*);
    };
    // special case: T=int not all implemented (yet?)
    template<> struct A16<int>
    {
      typedef int T;
      static void Ass(T*, size_t, T);
//       static void Ass(T*, size_t, T, const T*);
//       static void Neg(T*, size_t);
      static void Add(T*, size_t, T);
      static void Add(T*, size_t, const T*);
//       static void Add(T*, size_t, T, const T*);
//       static void Sub(T*, size_t, T);
//       static void Sub(T*, size_t, const T*);
//       static void Sub(T*, size_t, T, const T*);
//       static void Mul(T*, size_t, T);
//       static void Mul(T*, size_t, const T*);
//       static void Div(T*, size_t, const T*);
//       static void Inv(T*, size_t, T);
//       static void Inv(T*, size_t, T, const T*);
//       static T Sum(const T*, size_t);
//       static T Dot(const T*, size_t, const T*);
    };

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
      /// assign each element to scalar
      /// \code for(size_t i=0; i!=n; ++i) f[i] = x; \endcode
      template<typename T>
      static void Ass(T*f, size_t n, T x)
      { A16<T>::Ass(f,n,x); }

      /// assign element-wise to another array
      /// \code for(size_t i=0; i!=n; ++i) a[i] = b[i]; \endcode
      /// \note this does not require 16-byte alignment
      template<typename T>
      static void Ass(T*a, size_t n, const T*b)
      { std::memcpy(a,b,n*sizeof(T)); }

      /// assign element-wise to weighted array
      /// \code for(size_t i=0; i!=n; ++i) a[i] = w*b[i]; \endcode
      template<typename T>
      static void Ass(T*a, size_t n, T w, const T*b)
      { A16<T>::Ass(a,n,w,b); }

      /// replace each element with its negative
      /// \code for(size_t i=0; i!=n; ++i) f[i] = -f[i]; \endcode
      template<typename T>
      static void Neg(T*f, size_t n)
      { A16<T>::Neg(f,n); }

      /// set each element to zero
      /// \code for(size_t i=0; i!=n; ++i) f[i] = 0; \endcode
      template<typename T>
      static void Reset(T*f, size_t n)
      { Aligned::Ass(f,n,T(0)); }

      /// add same scalar to each element
      /// \code for(size_t i=0; i!=n; ++i) f[i] += x; \endcode
      template<typename T>
      static void Add(T*f, size_t n, T x)
      { A16<T>::Add(f,n,x); }

      /// add another array element wise
      /// \code for(size_t i=0; i!=n; ++i) a[i] += b[i]; \endcode
      template<typename T>
      static void Add(T*a, size_t n, const T*b)
      { A16<T>::Add(a,n,b); }

      /// add weighted array element wise
      /// \code for(size_t i=0; i!=n; ++i) a[i] += w*b[i]; \endcode
      template<typename T>
      static void Add(T*a, size_t n, T w, const T*b)
      { A16<T>::Add(a,n,w,b); }

      /// subtract same scalar from each element
      /// \code for(size_t i=0; i!=n; ++i) f[i] -= x; \endcode
      template<typename T>
      static void Sub(T*f, size_t n, T x)
      { A16<T>::Sub(f,n,x); }

      /// subtract another array element wise
      /// \code for(size_t i=0; i!=n; ++i) a[i] -= b[i]; \endcode
      template<typename T>
      static void Sub(T*a, size_t n, const T*b)
      { A16<T>::Sub(a,n,b); }

      /// subtract weighted array element wise
      /// \code for(size_t i=0; i!=n; ++i) a[i] -= w*b[i]; \endcode
      template<typename T>
      static void Sub(T*a, size_t n, T w, const T*b)
      { A16<T>::Sub(a,n,w,b); }

      /// multiply each element by same scalar
      /// \code for(size_t i=0; i!=n; ++i) f[i] *= x; \endcode
      template<typename T>
      static void Mul(T*f, size_t n, T x)
      { A16<T>::Mul(f,n,x); }

      /// multiply element-wise with another array
      /// \code for(size_t i=0; i!=n; ++i) a[i] *= b[i]; \endcode
      template<typename T>
      static void Mul(T*a, size_t n, const T*b)
      { A16<T>::Mul(a,n,b); }

      /// divide each element by same scalar
      /// \code for(size_t i=0; i!=n; ++i) f[i] /= x; \endcode
      template<typename T>
      static void Div(T*f, size_t n, T x)
      { A16<T>::Div(f,n); }

      /// divide by elements of another array
      /// \code for(size_t i=0; i!=n; ++i) a[i] /= b[i]; \endcode
      template<typename T>
      static void Div(T*a, size_t n, const T*b)
      { A16<T>::Div(a,n,b); }

      /// replace each elements by its inverse times a constant
      /// \code for(size_t i=0; i!=n; ++i) f[i] = x/f[i]; \endcode
      template<typename T>
      static void Inv(T*f, size_t n, T x)
      { A16<T>::Inv(f,n,x); }

      /// set to constant divided by element of another array
      /// \code for(size_t i=0; i!=n; ++i) a[i] = x/b[i]; \endcode
      template<typename T>
      static void Inv(T*a, size_t n, T x, const T*b)
      { A16<T>::Inv(a,n,x,b); }

      /// replace each element by its inverse
      /// \code for(size_t i=0; i!=n; ++i) f[i] = 1/f[i]; \endcode
      template<typename T>
      static void Reciprocal(T*f, size_t n)
      { Inv(f,n,T(1)); }

      /// set each element to the inverse of another array
      /// \code for(size_t i=0; i!=n; ++i) a[i] = 1/b[i]; \endcode
      template<typename T>
      static void Reciprocal(T*a, size_t n, const T*b)
      { Inv(a,n,T(1),b); }

      /// replace each element by its square root
      /// \code for(size_t i=0; i!=n; ++i) f[i] = std::sqrt(f[i]); \endcode
      template<typename T>
      static void Sqrt(T*f, size_t n)
      { A16<T>::Sqrt(f,n); }

      /// set each element to the square root of another array
      /// \code for(size_t i=0; i!=n; ++i) a[i] = std::sqrt(b[i]); \endcode
      template<typename T>
      static void Sqrt(T*a, size_t n, const T*b)
      { A16<T>::Sqrt(a,n,b); }

      //@}
      /// \name compute property of whole array
      //{@

      /// sum of all elements
      /// \code float S(0); for(size_t i=0; i!=n; ++i) S+=f[i]; return S;
      /// \endcode
      template<typename T>
      static T Sum(const T*f, size_t n)
      { return A16<T>::Sum(f,n); }

      /// dot product between two arrays
      /// \code float S(0); for(size_t i=0; i!=n; ++i) S+=a[i]*b[i]; return S;
      /// \endcode
      template<typename T>
      static T Dot(const T*a, size_t n, const T*b)
      { return A16<T>::Dot(a,n,b); }

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

    /// SSE::Traits<int>
    template<> struct Traits<int>
    {
      /// is SSE enabled for this type?
#ifdef __SSE2__
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
      { for(size_t i=_N; i!=_S; ++i) _A[i] = _F(0); }
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
      _F&operator[] (int i) { return _A[i]; }
      /// \name unary operations
      //@{
      /// reset element-wise to 0
      /// \code for(int i=0; i!=size(); ++i) A[i] = 0; \endcode
      Array16&reset()
      {
	A16<_F>::Ass(_A,_S,_F(0));
	return*this;
      }
      /// negate element-wise
      /// \code for(int i=0; i!=size(); ++i) A[i] = -A[i]; \endcode
      Array16&negate()
      {
	A16<_F>::Neg(_A,_S);
	return*this;
      }
      /// sum of all elements
      /// \code _F x(0); for(int i=0; i!=size(); ++i) x+=A[i]; return x;
      /// \endcode
      _F sum() const
      {
	reset_tail();
	return A16<_F>::Sum(_A,_S);
      }
      //@}
      /// \name binary operations with scalar
      //@{
      /// assign element-wise to scalar
      /// \code for(int i=0; i!=size(); ++i) A[i] = x; \endcode
      Array16&operator=(_F x)
      {
	A16<_F>::Ass(_A,_S,x);
	return*this;
      }
      /// add a scalar to each element
      /// \code for(int i=0; i!=size(); ++i) A[i] += x; \endcode
      Array16&operator+=(_F x)
      {
	A16<_F>::Add(_A,_S,x);
	return*this;
      }
      /// subtract a scalar from each element
      /// \code for(int i=0; i!=size(); ++i) A[i] -= x; \endcode
      Array16&operator-=(_F x)
      {
	A16<_F>::Sub(_A,_S,x);
	return*this;
      }
      /// multiply each element by a scalar
      /// \code for(int i=0; i!=size(); ++i) A[i] *= x; \endcode
      Array16&operator*=(_F x)
      {
	A16<_F>::Mul(_A,_S,x);
	return*this;
      }
      /// divide each element by a scalar
      /// \code for(int i=0; i!=size(); ++i) A[i] *= x; \endcode
      Array16&operator/=(_F x) WDutils_THROWING
      {
	A16<_F>::Div(_A,_S,x);
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
	A16<_F>::Add(_A,_S,B._A);
	return*this;
      }
      /// subtract element-wise
      /// \code for(int i=0; i!=size(); ++i) A[i] += B[i]; \endcode
      /// \note it is an error if the number of elements do not match
      Array16&operator-=(Array16 const&B) WDutils_THROWING
      {
	check_size(B,"operator-=(Array16&)");
	A16<_F>::Sub(_A,_S,B._A);
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
	return A16<_F>::Dot(_A,_S,B._A);
      }
      //@}
      /// \name tertiary operations with scalar and another Array16
      //@{
      /// add weighted element-wise
      /// \code for(int i=0; i!=size(); ++i) A[i] += w*B[i]; \endcode
      Array16&addtimes(Array16 const&B, _F w) WDutils_THROWING
      {
	check_size(B,"addtimes(Array16&)");
	A16<_F>::Add(_A,_S,w,B._A);
	return*this;
      }
      /// subtract weighted element-wise
      /// \code for(int i=0; i!=size(); ++i) A[i] += w*B[i]; \endcode
      Array16&subtimes(Array16 const&B, _F w) WDutils_THROWING
      {
	check_size(B,"subtimes(Array16&)");
	A16<_F>::Sub(_A,_S,w,B._A);
	return*this;
      }
      //@}
    };

  } // namespace SSE
} // namespace WDutils
//
#endif
