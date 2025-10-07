// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file    utils/inc/tupel.h
///
/// \author  Walter Dehnen
///
/// \date    1996-2013
/// 
/// \brief   contains the definition of template class WDutils::tupel and
///	     all its members and friends
///
/// \version aug-2003: template metaprogramming to unroll loops automatically
/// \version nov-2003: non-standard #include files redundant; pseudo_tupel
/// \version nov-2004: added volume(), applied();
/// \version may-2005: removed tupel::add_times, ass_times, sub_times
/// \version jun-2005: removed pseudo_tupel, const_pseudo_tupel
/// \version jul-2005: added doxygen documentation, added minnorm()
/// \version mar-2006: removed reliance on friend namespace injection
/// \version jul-2006: added global function applied()
/// \version sep-2006: added member method reset()
/// \version mar-2007: made data protected
/// \version may-2008: output manipulator "print()" supported
/// \version oct-2008: added static assertion to N-specific constructors
/// \version sep-2009: formatted output preserves width and precision
/// \version nov-2010: namespace WDutils::Tuple, allows easy export
/// \version aug-2011: meta-templating improved
/// \version jun-2012: starting implementing some C++11 stuff
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 1996-2013  Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 675
// Mass Ave, Cambridge, MA 02139, USA.
// 
////////////////////////////////////////////////////////////////////////////////
#ifndef WDutils_included_tupel_h
#define WDutils_included_tupel_h

#ifndef WDutils_included_iostream
#  include <iostream>
#  define WDutils_included_iostream
#endif
#ifndef WDutils_included_cmath
#  include <cmath>
#  define WDutils_included_cmath
#endif
#ifndef WDutils_included_cstdio
#  include <cstdio>
#  define WDutils_included_cstdio
#endif
////////////////////////////////////////////////////////////////////////////////
#ifndef WDutils_included_exception_h
# if __cplusplus >= 201103L
#  define WDutilsCXX11StaticAssert(TEST,MSGS) static_assert(TEST,MSGS);
# else
namespace WDutils {
  template<bool> struct STATIC_ASSERTION_FAILURE;
  template<>     struct STATIC_ASSERTION_FAILURE<true> { enum { value = 1 }; };
#  define WDutilsStaticAssert(TEST)				\
  enum { _DUMMY = sizeof(WDutils::STATIC_ASSERTION_FAILURE<	\
    static_cast<bool>((TEST))>)					\
  }
}
# endif // __cplusplus
#endif  // WDutils_included_exception_h
////////////////////////////////////////////////////////////////////////////////
#ifndef WDutils_included_tupel_cc
#  include <tupel.cc>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace WDutils {

#if __cplusplus < 201103L && !defined(WDutils_included_meta_h)
  template<bool B, class T = void>
  struct enable_if {};
  template<class T>
  struct enable_if<true, T> { typedef T type; };

  template<class A, class B = void>
  struct is_same { static const bool value = false; };
  template<class A>
  struct is_same<A,A> { static const bool value = true; };

  template<class T> struct remove_const          { typedef T type; };
  template<class T> struct remove_const<const T> { typedef T type; };
 
  template<class T> struct remove_volatile             { typedef T type; };
  template<class T> struct remove_volatile<volatile T> { typedef T type; };

  template<class T>
  struct remove_cv {
    typedef typename 
    remove_volatile<typename remove_const<T>::type>::type type;
  };
  template<class T>
  struct is_floating_point {
    static const bool value = 
      is_same<float, typename remove_cv<T>::type>::value  ||
      is_same<double, typename remove_cv<T>::type>::value  ||
      is_same<long double, typename remove_cv<T>::type>::value;
  };
#endif

  /// namespace for class tupel<int,typename>
  /// \note To allow \code using namespace WDutils::Tuple; \endcode
  namespace Tuple {
    
    // /////////////////////////////////////////////////////////////////////////
    //
    // class WDutils::tupel
    //
    /// \brief
    /// Template: a tupel of N scalars of type X, held in an array X[N].
    ///
    /// Class tupel<N,X> in essence is just a wrapper around X[N]. It is
    /// designed to act as a what physicists call "vector" in N-dimensional
    /// Euclidian space, similar to std::valarray, but with the number of
    /// elements as template parameter. The member methods, such as assign,
    /// additions, etc, are coded using meta template programming in file
    /// tupel.cc, resulting in maximum efficency. For N=3, the operator ^
    /// between tupels refers to the usual vector cross product. In the falcON
    /// project, tupel<3,real> are used to represent position, velocity and
    /// all other vector quantities.
    ///
    // /////////////////////////////////////////////////////////////////////////
    template<int N, typename X> class tupel {
      // private/protected types and data
    private:
      typedef meta::taux<X,N-1> M;
    protected:
      X a[N];                  ///< data: an array of N elements of type X
    public:
      //------------------------------------------------------------------------
      /// \name static members and public types
      //@{
      typedef X        element_type;               ///< type of elements 
      static const int NDAT = N;                   ///< number of elements 
      static const int ORD  = 1;                   ///< tensor rank or order
      static       int size() { return N; }        ///< number of elements
      //@}
      //------------------------------------------------------------------------
      /// \name constructors
      //@{
      /// unitialized construction
      tupel() {}
      /// copy from another tupel
      tupel(tupel const&x)
      { M::v_as(a,x.a); }
      /// all elements equal a scalar
      explicit tupel(X x)
      { M::s_as(a,x); }
      /// set elements equal to elements of an array
      explicit tupel(const X* x)
      { M::v_as(a,x); }
      /// copy from a tupel with another element type
      template<typename S> tupel (tupel<N,S> const&x)
      { M::v_as(a,static_cast<const S*>(x)); }
      /// all elements equal a scalar of type different from element type
      template<typename S> explicit tupel(const S*x)
      { M::v_as(a,x); }
      /// from 2 elements (for N=2)
      tupel(X x0, X x1)
      {
	static_assert(N==2, "Static assertion failed");
	a[0]=x0;
	a[1]=x1;
      }
      /// from 3 elements (for N=3)
      tupel(X x0, X x1, X x2)
      {
	static_assert(N==3, "Static assertion failed");
	a[0]=x0;
	a[1]=x1;
	a[2]=x2;
      }
      /// from 4 elements (for N=4)
      tupel(X x0, X x1, X x2, X x3)
      {
	static_assert(N==4, "Static assertion failed");
	a[0]=x0;
	a[1]=x1;
	a[2]=x2;
	a[3]=x3;
      }
      /// from 5 elements (for N=5)
      tupel(X x0, X x1, X x2, X x3, X x4)
      {
	static_assert(N==5, "Static assertion failed");
	a[0]=x0;
	a[1]=x1;
	a[2]=x2;
	a[3]=x3;
	a[4]=x4;
      }
      /// from 6 elements (for N=6)
      tupel(X x0, X x1, X x2, X x3, X x4, X x5)
      { 
	static_assert(N==6, "Static assertion failed");
	a[0]=x0;
	a[1]=x1;
	a[2]=x2;
	a[3]=x3;
	a[4]=x4;
	a[5]=x5;
      }
      /// from 7 elements (for N=7)
      tupel(X x0, X x1, X x2, X x3, X x4, X x5, X x6)
      {
	static_assert(N==7, "Static assertion failed");
	a[0]=x0;
	a[1]=x1;
	a[2]=x2;
	a[3]=x3;
	a[4]=x4;
	a[5]=x5;
	a[6]=x6;
      }
      /// from 8 elements (for N=8)
      tupel(X x0, X x1, X x2, X x3, X x4, X x5, X x6, X x7)
      {
	static_assert(N==8, "Static assertion failed");
	a[0]=x0;
	a[1]=x1;
	a[2]=x2;
	a[3]=x3;
	a[4]=x4;
	a[5]=x5;
	a[6]=x6;
	a[7]=x7;
      }
      /// from 9 elements (for N=9)
      tupel(X x0, X x1, X x2, X x3, X x4, X x5, X x6, X x7, X x8) {
	static_assert(N==9, "Static assertion failed");
	a[0]=x0;
	a[1]=x1;
	a[2]=x2;
	a[3]=x3;
	a[4]=x4;
	a[5]=x5;
	a[6]=x6;
	a[7]=x7;
	a[8]=x8;
      }
      /// from 10 elements (for N=10)
      tupel(X x0, X x1, X x2, X x3, X x4, X x5, X x6, X x7, X x8, X x9) 
      {
	static_assert(N==10, "Static assertion failed");
	a[0]=x0;
	a[1]=x1;
	a[2]=x2;
	a[3]=x3;
	a[4]=x4;
	a[5]=x5;
	a[6]=x6;
	a[7]=x7;
	a[8]=x8;
	a[9]=x9;
      }
      //@}
      //------------------------------------------------------------------------
      /// \name element access
      //@{
      /// non-const element access via sub-script operator
      X      &operator[] (int i)       { return a[i]; }
      /// const element access via sub-script operator
      X const&operator[] (int i) const { return a[i]; }
#ifdef WD_TUPEL_FUNCOP
      X      &operator() (int i)       { return a[i]; }
      X const&operator() (int i) const { return a[i]; }
#endif
      //@}
      //------------------------------------------------------------------------
      /// \name unitary operators and methods
      //@{
      /// reset all elements to zero
      tupel&reset()
      {
	M::s_ze(a);
	return*this;
      }
      /// change sign for each element
      tupel&negate()
      {
	M::v_neg(a);
	return*this;
      }
      /// return negative of this
      tupel operator- () const
      {
	tupel y;
	M::v_nega(y.a,a);
	return y;
      }
      /// return norm := Sum element^2
      X norm() const
      { return M::v_norm(a); }
      /// return |x| := sqrt(norm())
      X abs() const 
      { return std::sqrt(norm()); }
      /// return minimum element
      X min() const
      { return M::v_min(a); }
      /// return maximum element
      X max() const
      { return M::v_max(a); }
      /// return maximum |element|
      X maxnorm() const
      { return M::v_amax(a); }
      /// return minimum |element|
      X minnorm() const 
      {  return M::v_amin(a); }
      /// return Product of elements
      X volume() const
      { return M::v_vol(a); }
      /// conversion to pointer to scalar
      operator       X*  ()       { return a; }
      /// conversion to const pointer to scalar
      operator const X*  () const { return a; }
      /// conversion to pointer to scalar
      X*data() { return a; }
      /// conversion to const pointer to scalar
      const X*data() const { return a; }
      //@}
      //------------------------------------------------------------------------
      /// \name binary operators with scalar
      //@{
      /// assign all elements to scalar
      template<typename S> tupel&operator= (S x)
      {
	M::s_as(a,x);
	return*this;
      }
      /// multiply all elements by scalar
      template<typename S> tupel&operator*= (S x)
      {
	M::s_ml(a,x);
	return*this;
      }
      /// divide all elements by scalar
      template<typename S> tupel&operator/= (S x)
      { return operator*=(S(1)/x); }
      /// return product with scalar
      template<typename S> tupel operator* (S x) const
      {
	tupel y;
	M::v_ast(y.a,a,x);
	return y;
      }
      /// return product with inverse of scalar
      template<typename S> tupel operator/ (S x) const
      { return operator*(S(1)/x); }
      /// are all elements equal to scalar?
      bool operator== (X x) const
      { return M::s_eq(a,x); }
      /// is any element un-equal to scalar?
      bool operator!= (X x) const
      { return M::s_neq(a,x); }
      //@}
      //------------------------------------------------------------------------
      /// \name binary operators with tupel<N,X>
      //@{
      /// set *this equal to x: set tupel::a[i] = x[i]
      tupel&operator= (tupel const&x)
      {
	M::v_as(a,x.a);
	return*this;
      }
      /// add x to *this: set tupel::a[i] += x[i]
      tupel&operator+= (tupel const&x)
      {
	M::v_ad(a,x.a);
	return*this;
      }
      /// subtract x from *this: set tupel::a[i] -= x[i]
      tupel&operator-= (tupel const&x)
      {
	M::v_su(a,x.a);
	return*this;
      }
      /// element-wise *= : set tupel::a[i] *= x[i]
      tupel&ass_mul(tupel const&x)
      {
	M::v_ml(a,x.a);
	return*this;
      }
      /// element-wise /= : set tupel::a[i] /= x[i]
      tupel&ass_div(tupel const&x) {
	M::v_dv(a,x.a);
	return*this;
      }
      /// is x the same as *this?
      bool operator==(tupel const&x) const
      { return M::v_eq(a,x.a); }
      /// is x not the same as *this?
      bool operator!=(tupel const&x) const
      { return M::v_neq(a,x.a); }
      /// sum: return *this + x
      tupel operator+ (tupel const&x) const
      {
	tupel y;
	M::v_sum(y.a,a,x.a);
	return y;
      }
      /// difference: return *this - x
      tupel operator- (tupel const&x) const
      {
	tupel y;
	M::v_dif(y.a,a,x.a);
	return y;
      }
      /// absolute difference: return |*this - x| element wise
      tupel abs_diff(tupel const&x) const
      {
	tupel y;
	M::v_adif(y.a,a,x.a);
	return y;
      }
      /// scalar (dot) product: return Sum tupel::a[i]*x[i]
      X operator* (tupel const&x) const
      { return M::v_dot(a,x.a); }
      /// difference squared: return (*this-x)^2 := Sum (tupel::a[i]-x[i])^2
      X dist_sq(tupel const&x) const
      { return M::v_diq(a,x.a); }
      /// distance: return |*this-x| := sqrt(Sum (tupel::a[i]-x[i])^2)
      X distance(tupel const&x) const
      { return std::sqrt(dist_sq (x)); }
      /// sum squared: return (this+x)^2 := Sum (tupel::a[i]+x[i])^2
      X sum_sq (tupel const&x) const
      { return M::v_suq(a,x.a); }
      /// update maximum element-wise: tupel::a[i] = max(tupel::a[i], x[i])
      tupel&up_max (tupel const&x)
      {
	M::v_uma(a,x.a);
	return*this;
      }
      /// update minimum element-wise: tupel::a[i] = min(tupel::a[i], x[i])
      tupel&up_min (tupel const&x)
      {
	M::v_umi(a,x.a);
	return*this;
      }
      /// update minimum and maximum element-wise:
      /// Min[i] = min(Min[i], tupel::a[i]); and
      /// Max[i] = max(Max[i], tupel::a[i]);
      void up_min_max(tupel &Min, tupel&Max) const
      { M::v_umia(Min.a,Max.a,a); }
      //@}
      //------------------------------------------------------------------------
      /// \name binary operators with tupel<N,S>
      //@{
      /// set *this equal to x: set tupel::a[i] = x[i]
      template<typename S> tupel&operator= (tupel<N,S> const&x)
      {
	M::v_as(a,static_cast<const S*>(x));
	return*this;
      }
      /// add x to *this: set tupel::a[i] += x[i]
      template<typename S> tupel&operator+= (tupel<N,S> const&x)
      {
	M::v_ad(a,static_cast<const S*>(x));
	return*this;
      }
      /// subtract x from *this: set tupel::a[i] -= x[i]
      template<typename S> tupel&operator-= (tupel<N,S> const&x)
      {
	M::v_su(a,static_cast<const S*>(x));
	return*this;
      }
      /// sum: return *this + x
      template<typename S> tupel operator+ (tupel<N,S> const&x) const
      {
	tupel y;
	M::v_sum(y.a,a,static_cast<const S*>(x));
	return y;
      }
      /// difference: return *this - x
      template<typename S> tupel operator- (tupel<N,S> const&x) const
      {
	tupel y;
	M::v_dif(y.a,a,static_cast<const S*>(x));
	return y;
      }
      /// absolute difference: return |*this - x| element wise
      template<typename S> tupel abs_diff(tupel<N,S> const&x) const
      {
	tupel y;
	M::v_adif(y.a,a,x.a);
	return y;
      }
      /// scalar (dot) product: return Sum tupel::a[i]*x[i]
      template<typename S> X operator* (tupel<N,S> const&x) const
      { return M::v_dot(a,static_cast<const S*>(x)); }
      /// difference squared: return (*this-x)^2 := Sum (tupel::a[i]-x[i])^2
      template<typename S> X dist_sq(tupel<N,S> const&x) const
      { return M::v_diq(a,static_cast<const S*>(x)); }
      /// distance: return |*this-x| := sqrt(Sum (tupel::a[i]-x[i])^2)
      template<typename S> X distance (tupel<N,S> const&x) const
      { return std::sqrt(dist_sq (x)); }
      /// sum squared: return (this+x)^2 := Sum (tupel::a[i]+x[i])^2
      template<typename S> X sum_sq (tupel<N,S> const&x) const
      { return M::v_suq(a,static_cast<const S*>(x)); }
      /// update minimum and maximum element-wise:
      /// \param[in,out] xmin  replaced by min(xmin[i], tupel::a[i]);
      /// \param[in,out] xmax  replaced by max(xmax[i], tupel::a[i]);
      template<typename S>
      void up_min_max(tupel<N,S>&xmin,  tupel<N,S>&xmax) const
      { M::v_umia(static_cast<S*>(xmin),static_cast<S*>(xmax),a); }
      //@}
      //------------------------------------------------------------------------
      /// \name miscellaneous
      //@{
      /// set elements equal to array elements
      tupel&copy(const X*x)
      {
	M::v_as(a,x);
	return*this;
      }
      /// set elements equal to array elements
      template<typename S> tupel&copy(const S*x)
      {
	M::v_as(a,x);
	return*this;
      }
      /// update maximum element-wise: tupel::a[i] = max(tupel::a[i], x_i + f)
      tupel&up_max(tupel const&x, X f)
      {
	M::v_umax(a,x.a,f);
	return*this;
      }
      /// update minimum element-wise: tupel::a[i] = min(tupel::a[i], x_i - f)
      tupel&up_min(tupel const&x, X f)
      {
	M::v_umix(a,x.a,f);
	return*this;
      }
      /// apply element-wise function call: tupel::a[i] = f(tupel::a[i])
      tupel&apply(X(*f)(X))
      {
	M::v_appl(a,a,f);
	return*this;
      }
      /// return element-wise function call: return f(tupel::a[i])
      tupel applied(X(*f)(X)) const
      {
	tupel y;
	M::v_appl(y.a,a,f);
	return y;
      }
      /// replace by element-wise function: tupel::a[i] = f(x[i])
      tupel&connect(tupel const&x, X(*f)(X))
      {
	M::v_appl(a,x.a,f);
	return*this;
      }
      /// replace by element-wise function: tupel::a[i] = f(x[i])
      tupel connected(tupel const&x, X(*f)(X))
      {
	tupel y(*this);
	M::v_appl(y.a,x.a,f);
	return y;
      }
      /// normalise: x[i] /= abs(x)
      tupel&normalise()
      {
	X n = norm();
	if(n) return operator /= (std::sqrt(n));
	else  return*this;
      }
      /// normalises: return x[i] / abs(x)
      tupel normalised() const
      {
	tupel y(*this);
	return y.normalize();
      }
      //
      tupel&normalize() { return normalise(); }
      tupel normalized() const { return normalised(); }
      /// is any element nan?
      bool isnan() const
      { return M::v_nan(a); }
      /// is any element inf?
      bool isinf() const
      { return M::v_inf(a); }
      //------------------------------------------------------------------------
    };
    // /////////////////////////////////////////////////////////////////////////
    /// \relates WDutils::tupel
    /// \name vector cross product in 2D and 3D
    //@{
    /// vector cross product for N=2: returns x[0]*y[1]- x[1]*y[0]
    template<typename X> inline
    X operator^ (tupel<2,X> const&x, tupel<2,X> const&y)
    { return x[0]*y[1] - x[1]*y[0]; }
    /// vector cross product for N=3
    template<typename X> inline
    tupel<3,X> operator^ (tupel<3,X> const&x, tupel<3,X> const&y)
    {
      return tupel<3,X>(x[1]*y[2] - x[2]*y[1],
			x[2]*y[0] - x[0]*y[2],
			x[0]*y[1] - x[1]*y[0]);
    }
    //@}
    // /////////////////////////////////////////////////////////////////////////
    /// \relates WDutils::Tuple::tupel
    /// \name formatted I/O
    //@{
    /// formatted output: space separated; keeps width and precision
    template<int N, typename X> inline
    std::ostream&operator<<(std::ostream&s, tupel<N,X> const&x)
    {
      meta::taux<X,N-1>::v_out(s,x);
      return s;
    }
    /// formatted input: read element-wise
    template<int N, typename X> inline
    std::istream&operator>>(std::istream&s, tupel<N,X> &x)
    {
      meta::taux<X,N-1>::v_in(s,x);
      return s;
    }
    //@}
    // /////////////////////////////////////////////////////////////////////////
    /// \relates WDutils::tupel
    /// \name functions taking tupel<> arguments
    //@{
    /// is y[i] == x for all i?
    template<int N, typename X> inline
    bool operator==(X x,tupel<N,X> const&y)
    { return y==x; }
    /// is y[i] != x for any i?
    template<int N, typename X> inline
    bool operator!=(X x,tupel<N,X> const&y)
    { return y!=x; }
    /// return maximum element of tupel
    template<int N, typename X> inline
    X max(tupel<N,X> const&x)
    { return x.max(); }
    /// return maximum |element| of tupel
    template<int N, typename X> inline
    X maxnorm(tupel<N,X> const&x)
    { return x.maxnorm(); }
    /// return minimum element of tupel
    template<int N, typename X> inline
    X min(tupel<N,X> const&x)
    { return x.min(); }
    /// return minimum |element| of tupel
    template<int N, typename X> inline
    X minnorm(tupel<N,X> const&x)
    { return x.minnorm(); }
    /// return norm of tupel: Sum x[i]^2
    template<int N, typename X> inline
    X norm(tupel<N,X> const&x)
    { return x.norm(); }
    /// return absolute value: sqrt(Sum x[i]^2)
    template<int N, typename X> inline
    X abs(tupel<N,X> const&x)
    { return x.abs(); }
    /// return Product of elements
    template<int N, typename X> inline
    X volume (tupel<N,X> const&x)
    { return x.volume(); }
    /// return absolute difference between two vectors
    template<int N, typename X, typename S> inline
    tupel<N,X> abs_diff(tupel<N,X> const&x, tupel<N,S> const&y)
    { return x.abs_diff(y); }
    /// return difference squared: (x-y)^2 := Sum (x[i]-y[i])^2
    template<int N, typename X, typename S> inline
    X dist_sq(tupel<N,X> const&x, tupel<N,S> const&y)
    { return x.dist_sq(y); }
    /// return distance: (x-y)^2 := Sum (x[i]-y[i])^2
    template<int N, typename X, typename S> inline
    X distance(tupel<N,X> const&x, tupel<N,S> const&y)
    { return x.distance(y); }
    /// return sum squared: (x-y)^2 := Sum (x[i]+y[i])^2
    template<int N, typename X, typename S> inline
    X sum_sq(tupel<N,X> const&x, tupel<N,S> const&y)
    { return x.sum_sq(y); }
    /// product with scalar: x[i] = y * v[i]
    template<int N, typename X, typename S> inline
    typename enable_if<is_floating_point<S>::value, tupel<N,X> >::type
    operator*(S y, tupel<N,X> const&v)
    { return v*y; }
    /// return element-wise function call: return f(tupel::a[i])
    template<int N, typename X> inline
    tupel<N,X> applied(tupel<N,X> const&x, X(*f)(X))
    { return x.applied(f); }
    /// replace by element-wise function: tupel::a[i] = f(x[i])
    template<int N, typename X> inline
    tupel<N,X> connected(tupel<N,X> const&x, tupel<N,X> const&y, X(*f)(X))
    { return x.connected(y,f); }
    /// is any element nan?
    template<int N, typename X> inline
    bool isnan(tupel<N,X> const&x)
    { return x.isnan(); }
    /// is any element inf?
    template<int N, typename X> inline
    bool isinf(tupel<N,X> const&x)
    { return x.isinf(); }
    /// normalise: x[i] /= abs(x)
    template<int N, typename X> inline
    tupel<N,X> &normalise(tupel<N,X> &x)
    { return x.normalise(); }
    /// normalises: x[i] = abs(x)
    template<int N, typename X> inline
    tupel<N,X> normalised(tupel<N,X> const&x)
    { return x.normalised(); }
    //
    template<int N, typename X> inline
    tupel<N,X> &normalize(tupel<N,X> &x)
    { return x.normalize(); }
    //
    template<int N, typename X> inline
    tupel<N,X> normalized(tupel<N,X> const&x)
    { return x.normalized(); }
    /// update maximum element-wise: x[i] = max(x[i], y[i])
    template<int N, typename X> inline
    void update_max(tupel<N,X>&x, tupel<N,X> const&y)
    { return x.up_max(y); }
    /// update minimum element-wise: x[i] = min(x[i], y[i])
    template<int N, typename X> inline
    void update_min(tupel<N,X>&x, tupel<N,X> const&y)
    { return x.up_min(y); }
    //@}
  } // namespace WDutils::Tuple
  using namespace Tuple;
  // ///////////////////////////////////////////////////////////////////////////
#ifdef WDutils_included_traits_h
  template<int N, typename T> struct traits< tupel<N,T> > {
    static const char *name () {
      static char _name[1024]={0};
      if(_name[0]==0)
	sprintf(_name,"tupel<%d,%s>",N,traits<T>::name());
      return _name;
    }
    static const char *names() {
      return name();
    }
    static const unsigned size = sizeof(tupel<N,T>);
  };
#endif
  // ///////////////////////////////////////////////////////////////////////////
#ifdef WDutils_included_io_h
  template<int N, typename X>
  inline smanip_fp_vec_width<X> print(tupel<N,X> const&x, int w, int p) {
    return smanip_fp_vec_width<X>(static_cast<const X*>(x),N,w,p);
  }
#endif
} // namespace WDutils {
////////////////////////////////////////////////////////////////////////////////
#ifndef WDutils_included_exception_h
# undef WDutilsCXX11Delete
# undef WDutilsCXX11Default
# undef WDutilsCXX11DefaultBody
#endif
////////////////////////////////////////////////////////////////////////////////
#endif// WDutils_included_tupel_h
