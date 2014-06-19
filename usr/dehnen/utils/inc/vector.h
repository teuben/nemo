// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file    utils/inc/vector.h
///
/// \author  Walter Dehnen
///
/// \date    2012-2013
/// 
/// \brief   contains the definition of class template WDutils::vector<int,T>
///
////////////////////////////////////////////////////////////////////////////////
///
/// \version jun-2012  implemented, based on tupel.h, fully C++11
/// \version jan-2013  support for vector wrappers and derived types
///
////////////////////////////////////////////////////////////////////////////////
#if __cplusplus < 201103L
#  error requiring C++11
#elif !defined(WDutils_included_vector_h)
#define WDutils_included_vector_h
////////////////////////////////////////////////////////////////////////////////
#ifndef WDutils_included_iostream
#  include <iostream>
#  define WDutils_included_iostream
#endif
#ifndef WDutils_included_cmath
#  include <cmath>
#  define WDutils_included_cmath
#endif
#ifndef WDutils_included_type_traits
#  include <type_traits>
#  define WDutils_included_type_traits
#endif
////////////////////////////////////////////////////////////////////////////////
namespace WDutils {
  /// a vector of N elements of type T, hold in a buffer
  template<int N, typename T> struct vector;
  /// a vector of N elements of type T, referred to by a pointer
  template<int N, typename T> struct vector_wrapper;
  /// a constant vector of N elements of type T, referred to by a pointer
  /// \note allowed_element_type<T>::value must be true
  template<int N, typename T> struct const_vector_wrapper;
  /// namespace for auxiliary functionality for vector<>
  namespace vector_details {
    ///
    /// auxiliary class template for loop unrolling
    ///
    /// \note First > Last is allowed, impying a reverse loop
    ///
    template<int First, int Last>
    struct unroll
    {
      static_assert(Last!=First,"logic error");
      typedef unroll<(First<Last? First+1:First-1), Last> next;
      static const int I=First;
      /// output
      template<typename X>
      static std::ostream&out(std::ostream&s, const X*x,
			      std::streamsize w, std::streamsize p)
      {
	s << x[I] << ' ';
	s.width(w);
	s.precision(p);
	return next::out(s,x,w,p);
      }
      /// const unary
      template<typename X, typename UnaryFunc>
      static void const_unary(const X*x, UnaryFunc f)
	noexcept(noexcept(f))
      { f(x[I]); next::const_unary(x,f); }
      /// unary
      template<typename X, typename UnaryFunc>
      static void unary(X*x, UnaryFunc f)
	noexcept(noexcept(f))
      { f(x[I]); next::unary(x,f); }
      /// collect unary
      template<typename X, typename Operator, typename UnaryFunc>
      static auto coll_unary(const X*x, Operator op, UnaryFunc f)
	noexcept(noexcept(f)) -> decltype(op(f(x[0]),f(x[1])))
      { return op(f(x[I]), next::coll_unary(x,op,f)); }
      /// const binary
      template<typename X, typename Y, typename BinaryFunc>
      static void const_binary(const X*x, const Y*y, BinaryFunc f)
	noexcept(noexcept(f))
      { f(x[I],y[I]); next::const_binary(x,y,f); }
      /// binary
      template<typename X, typename Y, typename BinaryFunc>
      static void binary(X*x, const Y*y, BinaryFunc f)
	noexcept(noexcept(f))
      { f(x[I],y[I]); next::binary(x,y,f); }
      /// collect binary
      template<typename X, typename Y, typename Operator, typename BinaryFunc>
      static auto coll_binary(const X*x, const Y*y, Operator op, BinaryFunc f)
	noexcept(noexcept(f) && noexcept(op)) -> decltype(f(x[I],y[I]))
      { return op(f(x[I],y[I]), next::coll_binary(x,y,op,f)); }
      /// tertiary
      template<typename X, typename Y, typename Z, typename TertiaryFunc>
      static void tertiary(X*x, const Y*y, const Z*z, TertiaryFunc f)
	noexcept(noexcept(f))
      { f(x[I],y[I],z[I]); next::tertiary(x,y,z,f); }
      /// constant tertiary
      template<typename X, typename Y, typename Z, typename TertiaryFunc>
      static void const_tertiary(const X*x, Y*y, Z*z, TertiaryFunc f)
	noexcept(noexcept(f))
      { f(x[I],y[I],z[I]); next::const_tertiary(x,y,z,f); }
      /// variadic set
      /// \note must not be used with reverse loop (First > Last)
      template<typename X, typename Y, typename... Args>
      static void set(X*x, Y a, Args... args) noexcept
      {
	static_assert(First<Last,"variadic set in reverse loop");
	static_assert(std::is_convertible<Y,X>::value,
		      "cannot convert argument to element type");
	x[I]=X(a); next::set(x,args...);
      }
    };
    // closing the loop: last element
    template<int I>
    struct unroll<I,I>
    {
      // output
      template<typename X>
      static std::ostream&out(std::ostream&s, const X*x,
			      std::streamsize, std::streamsize)
      { return s << x[I]; }
      // const unary
      template<typename X, typename UnaryFunc>
      static void const_unary(const X*x, UnaryFunc f)
	noexcept(noexcept(f))
      { f(x[I]); }
      // unary
      template<typename X, typename UnaryFunc>
      static void unary(X*x, UnaryFunc f)
	noexcept(noexcept(f))
      { f(x[I]); }
      // collect unary
      template<typename X, typename Operator, typename UnaryFunc>
      static auto coll_unary(const X*x, Operator, UnaryFunc f)
	noexcept(noexcept(f)) -> decltype(f(x[I]))
      { return f(x[I]); }
      // const binary
      template<typename X, typename Y, typename BinaryFunc>
      static void const_binary(const X*x, const Y*y, BinaryFunc f)
	noexcept(noexcept(f))
      { f(x[I],y[I]); }
      // binary
      template<typename X, typename Y, typename BinaryFunc>
      static void binary(X*x, const Y*y, BinaryFunc f)
	noexcept(noexcept(f))
      { f(x[I],y[I]); }
      // collect binary
      template<typename X, typename Y, typename Operator, typename BinaryFunc>
      static auto coll_binary(const X*x, const Y*y, Operator, BinaryFunc f)
	noexcept(noexcept(f)) -> decltype(f(x[I],y[I]))
      { return f(x[I],y[I]); }
      // tertiary
      template<typename X, typename Y, typename Z, typename TertiaryFunc>
      static void tertiary(X*x, const Y*y, const Z*z, TertiaryFunc f)
	noexcept(noexcept(f))
      { f(x[I],y[I],z[I]); }
      // constant tertiary
      template<typename X, typename Y, typename Z, typename TertiaryFunc>
      static void const_tertiary(const X*x, Y*y, Z*z, TertiaryFunc f)
	noexcept(noexcept(f))
      { f(x[I],y[I],z[I]); }
      // variadic set
      template<typename X, typename Y>
      static void set(X*x, Y a) noexcept
      { 
	static_assert(std::is_convertible<Y,X>::value,
		      "cannot convert first argument to element type");
	x[I] = X(a);
      }
    };
    //  square of argument
    template<typename X>
    inline auto square(X a) noexcept -> decltype(a*a) { return a*a; }
    //  forward decl
    template<int, typename, typename> struct non_const_vector_methods;
    //
    //  code idea pinched from http://stackoverflow.com/questions/6534041
    //
    namespace check
    {
      typedef char No[7];                        // nothing else has 7 bytes
      // check for existence of x*y
      template<typename left, typename right>
      No&operator*(const left&, const right&);
      template <typename left, typename right=left>
      struct mul
      {
	typedef decltype( (*static_cast<left *>(0)) *
			  (*static_cast<right*>(0)) ) _R;
	static const bool value = !std::is_same<_R,No>::value;
	typedef typename std::conditional<value,_R,void>::type type;
      };
      // check for existence of x*=y
      template<typename left, typename right>
      No&operator*=(left&, const right&);
      template <typename left, typename right=left>
      struct mul_ass
      {
	typedef decltype( (*static_cast<left *>(0)) *=
			  (*static_cast<right*>(0)) ) _R;
	static const bool value = !std::is_same<_R,No>::value;
	typedef typename std::conditional<value,_R,void>::type type;
      };
      // check for existence of x/y
      template<typename left, typename right>
      No&operator/(const left&, const right&);
      template <typename left, typename right=left>
      struct div
      {
	typedef decltype( (*static_cast<left *>(0)) /
			  (*static_cast<right*>(0)) ) _R;
	static const bool value = !std::is_same<_R,No>::value;
	typedef typename std::conditional<value,_R,void>::type type;
      };
      // check for existence of x/=y
      template<typename left, typename right>
      No&operator/=(left&, const right&);
      template <typename left, typename right=left>
      struct div_ass
      {
	typedef decltype( (*static_cast<left *>(0)) /=
			  (*static_cast<right*>(0)) ) _R;
	static const bool value = !std::is_same<_R,No>::value;
	typedef typename std::conditional<value,_R,void>::type type;
      };
      // check for existence of x+y
      template<typename left, typename right>
      No&operator+(const left&, const right&);
      template <typename left, typename right=left>
      struct add
      {
	typedef decltype( (*static_cast<left *>(0)) +
			  (*static_cast<right*>(0)) ) _R;
	static const bool value = !std::is_same<_R,No>::value;
	typedef typename std::conditional<value,_R,void>::type type;
      };
      // check for existence of x+=y
      template<typename left, typename right>
      No&operator+=(left&, const right&);
      template <typename left, typename right=left>
      struct add_ass
      {
	typedef decltype( (*static_cast<left *>(0)) +=
			  (*static_cast<right*>(0)) ) _R;
	static const bool value = !std::is_same<_R,No>::value;
	typedef typename std::conditional<value,_R,void>::type type;
      };
      // check for existence of x-y
      template<typename left, typename right>
      No&operator-(const left&, const right&);
      template <typename left, typename right=left>
      struct sub
      {
	typedef decltype( (*static_cast<left* >(0)) -
			  (*static_cast<right*>(0)) ) _R;
	static const bool value = !std::is_same<_R,No>::value;
	typedef typename std::conditional<value,_R,void>::type type;
      };
      // check for existence of x+=y
      template<typename left, typename right>
      No&operator-=(left&, const right&);
      template <typename left, typename right=left>
      struct sub_ass
      {
	typedef decltype( (*static_cast<left *>(0)) -=
			  (*static_cast<right*>(0)) ) _R;
	static const bool value = !std::is_same<_R,No>::value;
	typedef typename std::conditional<value,_R,void>::type type;
      };
      // check for existence of x<y
      template<typename left, typename right>
      No&operator<(const left&, const right&);
      template <typename left, typename right=left>
      struct lt
      {
	typedef decltype( (*static_cast<left* >(0)) <
			  (*static_cast<right*>(0)) ) _R;
	static const bool value = !std::is_same<_R,No>::value;
	typedef typename std::conditional<value,_R,void>::type type;
      };
      // check for existence of x>y
      template<typename left, typename right>
      No&operator<(const left&, const right&);
      template <typename left, typename right=left>
      struct gt
      {
	typedef decltype( (*static_cast<left* >(0)) >
			  (*static_cast<right*>(0)) ) _R;
	static const bool value = !std::is_same<_R,No>::value;
	typedef typename std::conditional<value,_R,void>::type type;
      };
      // check for existence of bool equal(x,y)
      template<typename left, typename right>
      No&equal(const left&, const right&);
      template <typename left, typename right=left>
      struct eq
      {
	typedef decltype( (*static_cast<left* >(0)) >
			  (*static_cast<right*>(0)) ) _R;
	static const bool value = !std::is_same<_R,No>::value;
	typedef typename std::conditional<value,_R,void>::type type;
      };
      // check for existence of -x
      template<typename operand>
      No&operator-(const operand&);
      template <typename operand>
      struct neg
      {
	typedef decltype( -(*static_cast<operand*>(0))) _R;
	static const bool value = !std::is_same<_R,No>::value;
	typedef typename std::conditional<value,_R,void>::type type;
      };
    }
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif
    ///
    /// specifies whether vector operations *,/,+,- with other_type and/or
    /// @c vector<other_type> are implemented
    ///
    /// \note we allow this if                                             \n
    ///       1) the respective operator is defined,                       \n
    ///       2) neither type is bool, and                                 \n
    ///       3) the non-element type is convertible to the element type.  \n
    ///       You may override this by providing a specialisation.
#ifdef __clang__
#  pragma clang diagnostic pop
#endif
    template<template<typename, typename> class bin_op,
	     typename element_type, typename other_type = element_type>
    struct allow_bin_op
    {
      static const bool value =
	! std::is_same<bool,other_type>::value              &&
	! std::is_same<bool,element_type>::value            &&
	( std::is_same<other_type,element_type>::value ||
	  std::is_convertible<other_type,element_type>::value ) &&
	bin_op<element_type, other_type>::value;
    };
    ///
    /// functionality for const vector operations (observers)
    ///
    /// \note there are also some stand-alone observer functions and operators
    template<int N, typename X, typename V>
    struct const_vector_base
    {
      static_assert(N>=0,"negative number of vector elements");
      static_assert(N!=0,"zero vector elements");
      /// total \# coefficients
      static const int num_elem = N;
      /// element type
      typedef X element_type;
      /// derived vector type
      typedef V vector_type;
      //
      template<typename other_type, typename other_vector_type>
      using other_const_vector =
	const_vector_base<num_elem,other_type,other_vector_type>;
      //
      typedef vector_details::unroll<0,N-1> do_unroll;
      /// const data access
      const element_type*data() const noexcept 
      { return static_cast<const vector_type*>(this)->buffer; }
      /// const element access: X[n]
      element_type const&operator[] (int n) const noexcept
      { return data()[n]; }
      /// const element access: X[n]
      element_type const&operator[] (unsigned n) const noexcept
      { return data()[n]; }
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wfloat-equal"
#endif
      /// equality
      /// \note you should not compare floating point numbers for equality
      bool operator==(element_type s) const noexcept
      {
	return do_unroll::
	  coll_unary(data(),
		     [] (bool x, bool y) noexcept -> bool { return x&&y; },
		     [s](X x)            noexcept -> bool { return x==s; } );
      }
      /// inequality
      /// \note you should not compare floating point numbers for inequality
      bool operator!=(element_type s) const noexcept
      { return ! operator==(s); }
      /// equality
      /// \note you should not compare floating point numbers for equality
      template<typename other_vector_type>
      bool operator==
      (other_const_vector<element_type,other_vector_type> const&other)
	const noexcept
      {
	return do_unroll::
	  coll_binary(data(),other.data(),
		      [](bool x, bool y) noexcept -> bool { return x&&y; },
		      [](X x, X y)       noexcept -> bool { return x==y; } );
      }
#ifdef __clang__
#  pragma clang diagnostic pop
#endif
      /// inequality
      /// \note you should not compare floating point numbers for inequality
      template<typename other_vector_type>
      bool operator!=
      (other_const_vector<element_type,other_vector_type> const&other)
	const noexcept
      { return ! operator==(other); }
      /// update minimum and maximum element-wise:
      ///  Min[i] = min(Min[i], (*this)[i]); and
      ///  Max[i] = max(Max[i], (*this)[i]);
      template<typename min_vector_type, typename max_vector_type>
      inline void up_min_max
      (non_const_vector_methods<N,X,min_vector_type>&Min,
       non_const_vector_methods<N,X,max_vector_type>&Max) const noexcept;
      /// write vector element-wise to std::ostream; keeps width and precision
      friend std::ostream&operator<<(std::ostream&s,
				     const_vector_base const&v)
      {
	auto w = s.width();
	auto p = s.precision();
	return do_unroll::out(s,v.data(),w,p);
      }
      /// return the range [I0..IN[ with constant access
      template<int I0, int IN>
      const_vector_wrapper<(IN-I0),X> range() const noexcept;
    };// struct vector_details::qconst_vector_base<>
    ///
    /// functionality for non-const vector operations (modifiers)
    ///
    template<int N, typename X, typename V>
    struct non_const_vector_methods
    {
      /// total \# coefficients
      static const int num_elem = N;
      /// element type
      typedef X element_type;
      /// derived vector type
      typedef V vector_type;
      //
      template<typename other_type, typename other_vector_type>
      using other_const_vector =
	const_vector_base<num_elem,other_type,other_vector_type>;
      //
      typedef vector_details::unroll<0,N-1> do_unroll;
      /// data access
      element_type*data() noexcept
      { return static_cast<vector_type*>(this)->buffer; }
      /// element access: X[n]
      element_type &operator[] (int n) noexcept
      { return data()[n]; }
      /// element access: X[n]
      element_type &operator[] (unsigned n) noexcept
      { return data()[n]; }
      /// copy from (possibly) larger vector
      template<int other_num_elem, typename other_type,
	       typename other_vector_type>
      typename std::enable_if<std::is_convertible<other_type,
						  element_type>::value,
			      vector_type&>::type
      copy(const_vector_base<other_num_elem,other_type,other_vector_type>
	   const&other) noexcept
      {
	static_assert(other_num_elem>=num_elem,"cannot copy < few elements");
	do_unroll::binary(data(), other.data(),
			  [] (X&x, other_type y) noexcept { x=X(y); } );
	return static_cast<vector_type&>(*this);
      }
      /// reset all elements to zero
      vector_type& reset() noexcept
      {
	do_unroll::unary(data(), [] (X&x) noexcept { x=X(0); } );
	return static_cast<vector_type&>(*this);
      }
      /// set all elements to the same scalar
      vector_type& set(element_type s) noexcept
      {
	do_unroll::unary(data(), [s] (X&x) noexcept { x=s; } );
	return static_cast<vector_type&>(*this);
      }
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif
      /// set elements to arguments
      /// \note will not compile if                                        \n
      ///       1) the number of arguments doesn't match N or              \n
      ///       2) the arguments aren't convertible to element_type.
#ifdef __clang__
#  pragma clang diagnostic pop
#endif
      template<typename other_type, typename... Args>
      vector_type&set(other_type a0, Args... args) noexcept
      { 
	static_assert(sizeof...(args)   < N, "too many arguments");
	static_assert(sizeof...(args)+2 > N, "too few arguments");
	do_unroll::set(data(),a0,args...);
	return static_cast<vector_type&>(*this);
      }
      /// change sign for each element
      vector_type&negate() noexcept
      { 
	static_assert(check::neg<X>::value,
		      "unary minus for elements not supported: "
		      "cannot implement vector::negate()");
	do_unroll::unary(data(), [] (X&x) noexcept { x=-x; });
	return static_cast<vector_type&>(*this);
      }
      /// multiply all elements by scalar
      vector_type& operator*=(element_type s) noexcept
      {
	do_unroll::unary(data(), [s] (X&x) noexcept { x*=s; } );
	return static_cast<vector_type&>(*this);
      }
      /// divide all elements by scalar
      vector_type&operator/= (element_type s) noexcept
      { return operator*=(X(1)/s); }
      /// add another vector: element-wise +=
      template<typename other_type, typename other_vector_type>
      vector_type&operator+=
      (other_const_vector<other_type,other_vector_type> const&other) noexcept
      {
	static_assert(allow_bin_op<check::add_ass,X,other_type>::value,
		      "add-and-assign of elements not supported: "
		      "cannot implement vector += vector");
	do_unroll::binary(data(), other.data(),
			  [] (X&x, other_type y) noexcept { x+=y; } );
	return static_cast<vector_type&>(*this);
      }
      /// subtract another vector: element-wise -=
      template<typename other_type, typename other_vector_type>
      vector_type&operator-=
      (other_const_vector<other_type,other_vector_type> const&other) noexcept
      {
	static_assert(allow_bin_op<check::sub_ass,X,other_type>::value,
		      "subtract-and-assign of elements not supported: "
		      "cannot implement vector -= vector");
	do_unroll::binary(data(), other.data(),
			  [] (X&x, other_type y) noexcept { x-=y; } );
	return static_cast<vector_type&>(*this);
      }
      /// element wise *=
      template<typename other_type, typename other_vector_type>
      vector_type&multiply_element_wise
      (other_const_vector<other_type,other_vector_type> const&other) noexcept
      {
	static_assert(allow_bin_op<check::mul_ass,X,other_type>::value,
		      "multiply-and-assign of elements not supported: "
		      "cannot implement vector::multiply_element_wise(vector)");
	do_unroll::binary(data(), other.data(),
			  [] (X&x, other_type y) noexcept { x*=y; } );
	return static_cast<vector_type&>(*this);
      }
      /// element wise /=
      template<typename other_type, typename other_vector_type>
      vector_type&divide_element_wise
      (other_const_vector<other_type,other_vector_type> const&other) noexcept
      {
	static_assert(allow_bin_op<check::mul_ass,X,other_type>::value,
		      "divide-and-assign of elements not supported: "
		      "cannot implement vector::divide_element_wise(vector)");
	do_unroll::binary(data(), other.data(),
			  [] (X&x, other_type y) noexcept { x/=y; } );
	return static_cast<vector_type&>(*this);
      }
      /// element-wise maximum: a[i] = max(a[i],other[i])
      template<typename other_vector_type>
      vector_type&up_max
      (other_const_vector<element_type,other_vector_type> const&other) noexcept
      {
	static_assert(allow_bin_op<check::gt,X>::value,
		      "'>' comparison of elements not supported: "
		      "cannot implement vector::up_max()");
	do_unroll::binary(data(), other.data(),
			  [] (X&x, X y) noexcept { if(y>x) x=y; } );
	return static_cast<vector_type&>(*this);
      }
      /// element-wise minimum: a[i] = min(a[i],other[i])
      template<typename other_vector_type>
      vector_type&up_min
      (other_const_vector<X,other_vector_type> const&other) noexcept
      {
	static_assert(allow_bin_op<check::lt,X>::value,
		      "'<' comparison of elements not supported: "
		      "cannot implement vector::up_min()");
	do_unroll::binary(data(), other.data(),
			  [] (X&x, X y) noexcept { if(y<x) x=y; } );
	return static_cast<vector_type&>(*this);
      }
      /// apply unary function element wise: vec[i] = f(vec[i])
      /// \note f(element_type) must be convertible to element_type
      template<typename UnaryFunction>
      vector_type&apply(UnaryFunction f) noexcept(noexcept(f))
      {
	static_assert(std::is_convertible<decltype(f(data()[0])),
		      element_type>::value,
		      "cannot apply function: "
		      "result not convertible to element type");
	do_unroll::unary(data(), [f] (X&x) noexcept(noexcept(f)) { x=f(x); } );
	return static_cast<vector_type&>(*this);
      }
      /// connect with another vector element wise: vec[i] = f(vec[i], other[i])
      /// \note f(element_type,other_type) must be convertible to element_type
      template<typename other_type, typename other_vector_type,
	       typename binary_function>
      vector_type&connect
      (other_const_vector<other_type,other_vector_type> const&other,
       binary_function f) noexcept(noexcept(f))
      {
	static_assert(std::is_convertible<decltype(f(data()[0],other[0])),
		      element_type>::value,
		      "cannot connect using binary function: "
		      "result not convertible to element type");
	do_unroll::binary(data(), other.data(), [f] (X&x, other_type y)
			  noexcept(noexcept(f))  { x=f(x,y); } );
	return static_cast<vector_type&>(*this);
      }
      /// normalise: vec[i] /= abs(vec)
      vector_type&normalise() noexcept
      {
	static_assert(allow_bin_op<check::mul,element_type>::value,
		      "multiplication of elements not suported: "
		      "cannot normalise");
	typedef typename check::mul<element_type>::type M;
	static_assert(allow_bin_op<check::add,M>::value,
		      "adding up of squared elements not supported: "
		      "cannot normalise");
	typedef decltype(M(0)*M(0)) T;
	auto n = do_unroll::
	  coll_unary(data(),
		     [] (M x, M y) noexcept -> T { return x+y; },
		     [] (X x)      noexcept -> M { return x*x; });
	if(n && n!=T(1)) return operator /= (std::sqrt(n));
	else return static_cast<vector_type&>(*this);
      }
      /// read vector element-wise from std::istream
      friend std::istream&operator>>(std::istream&s,
				     non_const_vector_methods&v) noexcept
      {
	do_unroll::unary(v.data(),[&s] (X&x) { s >> x; } );
	return s;
      }
      /// return the range [I0..IN[
      template<int I0, int IN>
      inline vector_wrapper<(IN-I0),X> range() noexcept;
    };// struct non_const_vector_methods
    ///
    /// functionality for const and non-const vector operations
    ///
    template<int N, typename X, typename V>
    struct vector_base :
      const_vector_base<N,X,V>, non_const_vector_methods<N,X,V>
    {
      typedef const_vector_base<N,X,V> cbase;
      typedef non_const_vector_methods<N,X,V> nbase;
      //
      typedef X element_type;
      typedef V vector_type;
      //
      using cbase::num_elem;
      using cbase::data;
      using cbase::operator[];
      using cbase::up_min_max;
      using cbase::range;
      //
      using nbase::data;
      using nbase::reset;
      using nbase::copy;
      using nbase::set;
      using nbase::negate;
      using nbase::multiply_element_wise;
      using nbase::divide_element_wise;
      using nbase::up_max;
      using nbase::up_min;
      using nbase::apply;
      using nbase::connect;
      using nbase::normalise;
      using nbase::operator[];
      using nbase::range;
      /// assignment from single value
      vector_type&operator=(element_type s) noexcept
      { return nbase::set(s); }
    };
  } // namespace WDutils::vector_details
  ///
  /// array of @a N elements of type @c X with full vector functionalities
  ///
  /// \note  @c allowed_element_type defaults to @c std::is_arithmetic.
  ///        @a allowed_element_type<X>::value must be true. 
  template<int N, typename X>
  struct vector :
    public vector_details::vector_base<N,X,vector<N,X> >
  {
    typedef X element_type;
    typedef vector_details::vector_base<N,X,vector<N,X> > vbase;
    //
    using vbase::num_elem;
    using vbase::data;
    using vbase::reset;
    using vbase::set;
    using vbase::negate;
    using vbase::multiply_element_wise;
    using vbase::divide_element_wise;
    using vbase::up_max;
    using vbase::up_min;
    using vbase::up_min_max;
    using vbase::apply;
    using vbase::connect;
    using vbase::normalise;
    using vbase::range;
    /// \name constructors and assignment operators
    //@{
    /// default constructor
    vector() = default;
    /// copy constructor
    vector(vector const&) = default;
    /// casting copy constructor
    template<typename other_type, typename other_vector_type,
	     class = typename
	     std::enable_if<std::is_convertible<other_type,
						element_type>::value
			    >::type >
    vector(vector_details::
	   const_vector_base<num_elem,other_type,other_vector_type> const&other)
    { vbase::copy(other); }
    /// assignment from another vector
    vector&operator=(vector const&) = default;
    /// casting assignment from vector of other type
    template<typename other_type, typename other_vector_type>
    typename std::enable_if<std::is_convertible<other_type,
						element_type>::value,
			    vector&>::type
    operator=(vector_details::
	      const_vector_base<num_elem,other_type,other_vector_type>
	      const&other)
    { return vbase::copy(other); }
    ///
    using vbase::operator=;
    /// construction from N values
    template<typename... Args>
    vector(element_type a0, Args... args)
    { vbase::set(a0, args...); }
    //@}
    /// the actual data
    //  NOTE  Xhere is no point in making buffer private, as access is granted
    //        anyway. Moreover it would invalidate the current template magic.
    element_type buffer[num_elem];
  private:
    using vbase::copy;
  };// class WDutils::vector<>
  ///
  /// interpreting a data buffer as vector with full functionalities
  ///
  template<int N, typename X>
  struct vector_wrapper :
    vector_details::vector_base<N,X,vector_wrapper<N,X> >
  {
    typedef X element_type;
    typedef vector_details::vector_base<N,X,vector_wrapper<N,X> > vbase;
    //
    using vbase::num_elem;
    using vbase::data;
    using vbase::reset;
    using vbase::set;
    using vbase::negate;
    using vbase::multiply_element_wise;
    using vbase::divide_element_wise;
    using vbase::up_max;
    using vbase::up_min;
    using vbase::up_min_max;
    using vbase::apply;
    using vbase::connect;
    using vbase::normalise;
    using vbase::range;
    /// \name constructors and assignment operators
    //@{
    /// no default constructor
    vector_wrapper() = delete;
    /// no copy constructor
    vector_wrapper(vector_wrapper const&) = delete;
    /// move constructor
    vector_wrapper(vector_wrapper &&) = default;
    /// constructor from pointer to data
    /// \note explicit not sensible
    vector_wrapper(element_type*b) noexcept
      : buffer(b) {}
    /// constructor from a another vector: just wrap other[0,N[
    template<typename other_vector_type>
    explicit vector_wrapper
    (vector_details::
     non_const_vector_methods<num_elem,element_type,other_vector_type>&other)
      noexcept : buffer(other.data()) {}
    //@}
    /// pointer to the actual data held somewhere else
    //  NOTE  There is no point in making buffer private, as access is granted
    //        anyway. Moreover it would invalidate the current template magic.
    element_type*const buffer;
  };// class WDutils::vector_wrapper<>
  ///
  /// interpreting a constant data buffer as vector with constant functionality
  ///
  template<int N, typename X>
  struct const_vector_wrapper :
    vector_details::const_vector_base<N,X,const_vector_wrapper<N,X>>
  {
    typedef
      vector_details::const_vector_base<N,X,const_vector_wrapper<N,X> > cbase;
    //
    typedef X element_type;
    //
    using cbase::num_elem;
    using cbase::data;
    using cbase::up_min_max;
    using cbase::range;
    /// \name constructors and assignment operators
    //@{
    /// no default constructor
    const_vector_wrapper() = delete;
    /// no copy constructor
    const_vector_wrapper(const_vector_wrapper const&) = delete;
    /// move constructor
    const_vector_wrapper(const_vector_wrapper &&) = default;
    /// constructor from pointer to data
    /// \note explicit not sensible
    const_vector_wrapper(const element_type*b) noexcept
      : buffer(b) {}
    /// constructor from a another vector: just wrap other[0,N[
    template<typename other_vector_type>
    explicit const_vector_wrapper
    (vector_details::const_vector_base<num_elem,element_type,other_vector_type>
     const &other) noexcept
      : buffer(other.data()) {}
    //@}
    /// pointer to the actual data held somewhere else
    //  NOTE  There is no point in making buffer private, as access is granted
    //        anyway. Moreover it would invalidate the current template magic.
    const element_type*const buffer;
  };// class WDutils::const_vector_wrapper<>
  //
  namespace vector_details {
    //
    // members of const_vector_base<>
    //
    // update minimum and maximum
    template<int N, typename X, typename V>
    template<typename min_vector_type, typename max_vector_type>
    inline void const_vector_base<N,X,V>::up_min_max
    (non_const_vector_methods<N,X,min_vector_type>&Min,
     non_const_vector_methods<N,X,max_vector_type>&Max) const noexcept
    {
      do_unroll::const_tertiary(data(), Min.data(), Max.data(),
				[] (X x, X&mi, X&ma) noexcept 
				{ if     (x<mi) mi=x;
				  else if(x>ma) ma=x; } );
    }
    // constant range
    template<int N, typename X, typename V>  template<int I0, int IN>
    inline const_vector_wrapper<(IN-I0),X>
    const_vector_base<N,X,V>::range() const noexcept
    {
      static_assert(0<=I0,"range includes negative indices");
      static_assert(I0<IN,"range's lowest > highest");
      static_assert(IN<=N,"range exceeds available");
      return { data()+I0 };
    }
    //
    // members of non_const_vector_methods<>
    //
    // non-constant range
    template<int N, typename X, typename V> template<int I0, int IN>
    inline vector_wrapper<(IN-I0),X>
    non_const_vector_methods<N,X,V>::range() noexcept
    { 
      static_assert(0<=I0,"range includes negative indices");
      static_assert(I0<IN,"range's lowest > highest");
      static_assert(IN<=N,"range exceeds available");
      return { data()+I0 };
    }
    //
    // non-member functions taking vector arguments and/or returning vector
    //
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif
    /// unary minus
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename V>
    inline vector<N,X> operator-(const_vector_base<N,X,V> const&v) noexcept
    {
      static_assert(check::neg<X>::value,
		    "unary minus for elements no supported: "
		    "cannot implement -vector");
      vector<N,X> ret;
      unroll<0,N-1>::binary(ret.data(), v.data(),
			    [] (X&r, X x) noexcept { r=-x; });
      return ret;
    }
    /// vector sum
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename Y, typename V, typename W>
    inline typename std::conditional<check::add<X,Y>::value,
				     vector<N, typename check::add<X,Y>::type>,
				     void>::type
    operator+(const_vector_base<N,X,V> const&v,
	      const_vector_base<N,Y,W> const&w) noexcept
    {
      static_assert(allow_bin_op<check::add,X,Y>::value,
		    "addition of elements not supported: "
		    "cannot implement vector+vector");
      typedef typename check::add<X,Y>::type R;
      vector<N,R> ret;
      unroll<0,N-1>::tertiary(ret.data(),v.data(),w.data(),
			      [] (R&r, X x, Y y) noexcept { r=x+y; } );
      return ret;
    }
    /// vector difference
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename Y, typename V, typename W>
    inline typename std::conditional<check::add<X,Y>::value,
				     vector<N, typename check::sub<X,Y>::type>,
				     void>::type
    operator-(const_vector_base<N,X,V> const&v,
	      const_vector_base<N,Y,W> const&w) noexcept
    {
      static_assert(allow_bin_op<check::sub,X,Y>::value,
		    "subtraction of elements not supported: "
		    "cannot implement vector-vector");
      typedef typename check::sub<X,Y>::type R;
      vector<N,R> ret;
      unroll<0,N-1>::tertiary(ret.data(),v.data(),w.data(),
			      [] (R&r, X x, Y y) noexcept { r=x-y; } );
      return ret;
    }
    /// vector times scalar
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename V, typename S>
    inline typename
    std::enable_if<allow_bin_op<check::mul,X,S>::value,
		   vector<N, typename check::mul<X,S>::type> >::type
    operator*(const_vector_base<N,X,V> const&v, S s) noexcept
    {
      typedef typename check::mul<X,S>::type R;
      vector<N,R> ret;
      unroll<0,N-1>::binary(ret.data(), v.data(),
			    [s](R&r, X x) noexcept { r=x*s; } );
      return ret;
    }
    /// scalar times vector
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename V, typename S>
    inline typename
    std::enable_if<allow_bin_op<check::mul,X,S>::value,
		   vector<N, typename check::mul<X,S>::type> >::type
    operator*(S s, const_vector_base<N,X,V> const&v) noexcept
    {
      typedef typename check::mul<S,X>::type R;
      vector<N,R> ret;
      unroll<0,N-1>::binary(ret.data(), v.data(),
			    [s](R&r, X x) noexcept { r=s*x; } );
      return ret;
    }
#ifdef WDutils_included_WDMath_h
    /// element-wise equality, save for floating point types
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename V, typename W>
    inline bool equal(const_vector_base<N,X,V> const&v,
		      const_vector_base<N,X,W> const&w) noexcept
    {
      static_assert(check::eq<X,X>::value,
		    "equal(element_type,element_type) not defined");
      static_assert(std::is_same<typename check::eq<X,X>::type, bool>::value,
		    "equal(element_type,element_type) not returning bool");
      return unroll<0,N-1>::
	coll_binary(v.data(), w.data(),
		    [] (bool x, bool y) noexcept->bool
		    { return x && y; },
		    [] (X    x, X    y) noexcept->bool
		    { return WDutils::equal(x,y); } );
    }
#endif
    /// vector dot product
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename Y, typename V, typename W>
    inline auto dot(const_vector_base<N,X,V> const&v,
		    const_vector_base<N,Y,W> const&w) noexcept
      -> decltype(v[0]*w[0]+v[1]*w[1])
    {
      static_assert(allow_bin_op<check::mul,X,Y>::value,
		    "element type multiplication not supported: "
		    "cannot implement vector dot product");
      typedef typename check::mul<X,Y>::type M;
      static_assert(allow_bin_op<check::add,M,M>::value,
		    "addition of product of elements not supported: "
		    "cannot implement vector dot product");
      typedef typename check::add<M,M>::type R;
      return unroll<0,N-1>::
	coll_binary(v.data(), w.data(),
		    [] (M x, M y) noexcept ->R { return x+y; },
		    [] (X x, Y y) noexcept ->M { return x*y; } );
    }
    /// vector dot product
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename Y, typename V, typename W>
    inline auto operator*(const_vector_base<N,X,V> const&v,
			  const_vector_base<N,Y,W> const&w) noexcept
      -> decltype(v[0]*w[0]+v[1]*w[1])
    { return dot(v,w); }
    /// vector division by scalar
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename V, typename S>
    inline auto operator/(const_vector_base<N,X,V> const&v, S s) noexcept
      -> decltype(v*(X(1)/s)) { return v * (X(1)/s); }
    /// vector cross product for N=2: returns x[0]*y[1] - x[1]*y[0]
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<typename X, typename Y, typename V, typename W>
    inline typename check::sub<typename check::mul<X,Y>::type>::type
    operator^(const_vector_base<2,X,V> const&x,
	      const_vector_base<2,Y,W> const&y) noexcept
    {
      static_assert(allow_bin_op<check::mul,X,Y>::value,
		    "multiplication of elements not supported: "
		    "cannot implement vector<2> ^ vector<2>");
      typedef typename check::mul<X,Y>::type M;
      static_assert(allow_bin_op<check::sub,M>::value,
		    "subtraction of element products not supported: "
		    "cannot implement vector<2> ^ vector<2>");
      return x[0]*y[1] - x[1]*y[0];
    }
    /// vector cross product for N=3
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<typename X, typename Y, typename V, typename W>
    inline vector<3, typename check::sub<typename check::mul<X,Y>::type>::type>
    operator^(const_vector_base<3,X,V> const&x,
	      const_vector_base<3,Y,W> const&y) noexcept
    {
      static_assert(allow_bin_op<check::mul,X,Y>::value,
		    "multiplication of elements not supported: "
		    "cannot implement vector<3> ^ vector<3>");
      typedef typename check::mul<X,Y>::type M;
      static_assert(allow_bin_op<check::sub,M>::value,
		    "subtraction of element products not supported: "
		    "cannot implement vector<3> ^ vector<3>");
      typedef typename check::sub<M>::type R;
      return vector<3,R>(x[1]*y[2] - x[2]*y[1],
			 x[2]*y[0] - x[0]*y[2],
			 x[0]*y[1] - x[1]*y[0]);
    }
    /// norm: sum of vector elements squared
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename V>
    inline typename check::add<typename check::mul<X>::type>::type
    norm(const_vector_base<N,X,V> const&v) noexcept
    {
      static_assert(allow_bin_op<check::mul,X>::value,
		    "element type multiplication not supported: "
		    "cannot implement vector norm");
      typedef typename check::mul<X>::type M;
      static_assert(allow_bin_op<check::add,M,M>::value,
		    "addition of product of elements not supported: "
		    "cannot implement vector norm");
      typedef typename check::add<M,M>::type R;
      return unroll<0,N-1>::
	coll_unary(v.data(),
		   [](M x, M y) noexcept -> R { return x+y; },
		   [](X x)      noexcept -> M { return x*x; });
    }
    /// absolute value: sqrt(norm(vector))
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename V>
    inline auto abs(const_vector_base<N,X,V> const&v)
      noexcept -> decltype(std::sqrt(norm(v)))
    {
      return std::sqrt(norm(v));
    }
    /// product of all elements
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename V>
    inline typename check::mul<X>::type
    volume(const_vector_base<N,X,V> const&v) noexcept
    {
      static_assert(allow_bin_op<check::mul,X>::value,
		    "multiplication of elements not supported: "
		    "cannot implement volume(vector)");
      typedef typename check::mul<X>::type M;
      static_assert(std::is_same<X,M>::value,
		    "element multiplication not type preserving: "
		    "cannot implement volume(vector)");
      return unroll<0,N-1>::
	coll_unary(v.data(),
		   [](X x, X y) noexcept ->X { return x*y; },
		   [](X x)      noexcept ->X { return x; });
    }
    /// minimum element of vector
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename V>
    inline X min(const_vector_base<N,X,V> const&v) noexcept
    {
      static_assert(allow_bin_op<check::lt,X>::value,
		    "'<' comparison of elements not supported: "
		    "cannot implement min(vector)");
      return unroll<0,N-1>::
	coll_unary(v.data(),
		   [](X x, X y) noexcept ->X { return x<y?x:y; },
		   [](X x)      noexcept ->X { return x; });
    }
    /// minimum |element| of vector
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename V>
    inline X min_abs(const_vector_base<N,X,V> const&v) noexcept
    {
      static_assert(std::is_same<X,decltype(std::abs(X(1)))>::value,
		    "std::abs(elements) not type preserving: "
		    "cannot implement min_abs(vector)");
      static_assert(allow_bin_op<check::lt,X>::value,
		    "'<' comparison of elements not supported: "
		    "cannot implement min_abs(vector)");
      return unroll<0,N-1>::
	coll_unary(v.data(),
		   [](X x, X y) noexcept ->X { return x<y?x:y; },
		   [](X x)      noexcept ->X { return std::abs(x); });
    }
    /// maximum element of vector
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename V>
    inline X max(const_vector_base<N,X,V> const&v) noexcept
    {
      static_assert(allow_bin_op<check::gt,X>::value,
		    "'>' comparison of elements not supported: "
		    "cannot implement max(vector)");
      return unroll<0,N-1>::
	coll_unary(v.data(),
		   [](X x, X y) noexcept -> X { return x>y?x:y; },
		   [](X x)      noexcept -> X { return x; });
    }
    /// maximum |element| of vector
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename V>
    inline X max_abs(const_vector_base<N,X,V> const&v) noexcept
    {
      static_assert(std::is_same<X,decltype(std::abs(X(1)))>::value,
		    "std::abs(elements) not type preserving: "
		    "cannot implement max_abs(vector)");
      static_assert(allow_bin_op<check::gt,X>::value,
		    "'>' comparison of elements not supported: "
		    "cannot implement max_abs(vector)");
      return unroll<0,N-1>::
	coll_unary(v.data(),
		   [](X x, X y) noexcept -> X { return x>y?x:y; },
		   [](X x) noexcept -> X { return std::abs(x); });
    }
    /// are all elements zero?
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename V>
    inline bool is_zero(const_vector_base<N,X,V> const&v) noexcept
    {
      return unroll<0,N-1>::
	coll_unary(v.data(),
		   [](bool x, bool y) noexcept -> bool { return x && y; },
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wfloat-equal"
#endif
		   [](X x) noexcept -> bool { return x == X(0); });
#ifdef __clang__
#  pragma clang diagnostic pop
#endif
    }
    /// is any element NaN?
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename V>
    bool isnan(const_vector_base<N,X,V> const&v) noexcept
    {
      return unroll<0,N-1>::
	coll_unary(v.data(),
		   [](bool x, bool y) noexcept -> bool { return x || y; },
		   [](X x) noexcept -> bool { return std::isnan(x); });
    }
    /// is any element inf?
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename V>
    bool isinf(const_vector_base<N,X,V> const&v) noexcept
    {
      return unroll<0,N-1>::
	coll_unary(v.data(),
		   [](bool x, bool y) noexcept -> bool { return x || y; },
		   [](X x) noexcept -> bool { return std::isinf(x); });
    }
    /// vector difference squared
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename Y, typename V, typename W>
    inline typename
    check::add<typename check::mul<typename check::sub<X,Y>::type>::type>::type
    dist_sq(const_vector_base<N,X,V> const&v,
	    const_vector_base<N,Y,W> const&w) noexcept
    {
      static_assert(allow_bin_op<check::sub,X,Y>::value,
		    "subtraction of elements not supported: "
		    "cannot implement dist_sq(vector,vector)");
      typedef typename check::sub<X,Y>::type S;
      static_assert(allow_bin_op<check::sub,X,Y>::value,
		    "multiplication of element differences not supported: "
		    "cannot implement dist_sq(vector,vector)");
      typedef typename check::mul<S,S>::type M;
      static_assert(allow_bin_op<check::sub,X,Y>::value,
		    "addition of squared element differences not supported: "
		    "cannot implement dist_sq(vector,vector)");
      typedef typename check::add<M,M>::type R;
      return unroll<0,N-1>::
	coll_binary(v.data(), w.data(),
		    [](M x, M y) noexcept ->R { return x+y; },
		    [](X x, Y y) noexcept ->M { return square(x-y); } );
    }
    /// vector distance = sqrt( norm(vector difference) )
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename Y, typename V, typename W>
    inline auto distance(const_vector_base<N,X,V> const&v,
			 const_vector_base<N,Y,W> const&w)
      noexcept -> decltype(std::sqrt(dist_sq(v,w)))
    { return std::sqrt(dist_sq(v,w)); }
    /// vector sum squared
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename Y, typename V, typename W>
    inline typename
    check::add<typename check::mul<typename check::add<X,Y>::type>::type>::type
    sum_sq(const_vector_base<N,X,V> const&v,
	   const_vector_base<N,Y,W> const&w) noexcept
    {
      static_assert(allow_bin_op<check::add,X,Y>::value,
		    "addition of elements not supported: "
		    "cannot implement sum_sq(vector,vector)");
      typedef typename check::sub<X,Y>::type S;
      static_assert(allow_bin_op<check::sub,X,Y>::value,
		    "multiplication of element differences not supported: "
		    "cannot implement sum_sq(vector,vector)");
      typedef typename check::mul<S,S>::type M;
      static_assert(allow_bin_op<check::sub,X,Y>::value,
		    "addition of squared element differences not supported: "
		    "cannot implement sum_sq(vector,vector)");
      typedef typename check::add<M,M>::type R;
      return unroll<0,N-1>::
	coll_binary(v.data(), w.data(),
		    [](M x, M y) noexcept ->R { return x+y; },
		    [](X x, Y y) noexcept ->M { return square(x-y); } );
    }
    /// apply unary function: R[i] = f(X[i]); f() doesn't return void
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename V, typename F>
    inline auto apply(const_vector_base<N,X,V> const&v, F f)
      noexcept(noexcept(f))
      -> typename std::enable_if<!std::is_void<decltype(f(v[0]))>::value,
				 vector<N,decltype(f(v[0]))> >::type
    {
      typedef decltype(f(v[0])) R;
      vector<N,R> ret;
      unroll<0,N-1>::binary(ret.data(), v.data(),
			    [f] (R&r, X x) noexcept { r=f(x); } );
      return ret;
    }
    /// apply binary function: R[i] = f(X[i],Y[i]); f() doesn't return void
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename V, typename Y, typename W,  typename F>
    inline auto apply(const_vector_base<N,X,V> const&v,
		      const_vector_base<N,Y,W> const&w, F f)
      noexcept(noexcept(f))
      -> typename std::enable_if<!std::is_void<decltype(f(v[0],w[0]))>::value,
				 vector<N,decltype(f(v[0],w[0]))> >::type
    {
      typedef decltype(f(v[0],w[0])) R;
      vector<N,R> ret;
      unroll<0,N-1>::
	tertiary(ret.data(), v.data(), w.data(),
		 [f] (R&r, X x, Y y) noexcept { r=f(x,y); } );
      return ret;
    }
    /// apply unary function: f(X[i]); f() returns void
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename V, typename F>
    inline auto apply(const_vector_base<N,X,V> const&v, F f)
      noexcept(noexcept(f))
      -> typename std::enable_if<std::is_void<decltype(f(v[0]))>::value>::type
    {
      unroll<0,N-1>::const_unary(v.data(),f);
    }
    /// apply binary function: f(X[i],Y[i]); f() returns void
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename V, typename Y, typename W,  typename F>
    inline auto apply(const_vector_base<N,X,V> const&v,
		      const_vector_base<N,Y,W> const&w, F f)
      noexcept(noexcept(f))
      -> typename std::enable_if<std::is_void<decltype(f(v[0],w[0]))>
				 ::value>::type
    {
      unroll<0,N-1>::const_binary(v.data(),w.data(),f);
    }
    /// vector absolute difference
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename V, typename Y, typename W>
    inline typename
    std::conditional<check::sub<X,Y>::value,
		     vector<N, typename check::sub<X,Y>::type>, void
		     >::type
    abs_diff(const_vector_base<N,X,V> const&v,
	     const_vector_base<N,Y,W> const&w)
    {
      static_assert(allow_bin_op<check::sub,X,Y>::value,
		    "subtraction of elements not supported: "
		    "cannot implement abs_diff(vector,vector)");
      typedef typename check::sub<X,Y>::type R;
      vector<N,R> ret;
      unroll<0,N-1>::
	tertiary(ret.data(), v.data(), w.data(),
		 [] (R&r, X x, Y y) noexcept { r=std::abs(x-y); } );
      return ret;
    }
    /// returns @a x/|x|, the unit vector in direction @a x
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int N, typename X, typename V>
    inline auto normalised(const_vector_base<N,X,V> const&v)
      noexcept -> decltype(v/std::sqrt(norm(v)))
    {
      auto n = norm(v);
      if(n>0) return v / (std::sqrt(n));
      else  return v;
    }
    /// range of elements
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int I0, int IN, int N, typename X, typename V>
    inline auto range(const_vector_base<N,X,V> const&vec)
      noexcept -> decltype(vec.template range<I0,IN>())
    { return vec.template range<I0,IN>(); }
    /// range of elements
    /// \relates WDutils::vector
    /// \relates WDutils::vector_wrapper
    /// \relates WDutils::const_vector_wrapper
    template<int I0, int IN, int N, typename X, typename V>
    inline auto range(non_const_vector_methods<N,X,V>&vec)
      noexcept -> decltype(vec.template range<I0,IN>())
    { return vec.template range<I0,IN>(); }
#ifdef __clang__
#  pragma clang diagnostic pop
#endif
} // namespace WDutils::vector_details
  //
  using vector_details::norm;
  using vector_details::abs;
  using vector_details::volume;
  using vector_details::dot;
  using vector_details::min;
  using vector_details::max;
  using vector_details::min_abs;
  using vector_details::max_abs;
  using vector_details::is_zero;
#ifdef WDutils_included_WDMath_h
  using vector_details::equal;
#endif
  using vector_details::isnan;
  using vector_details::isinf;
  using vector_details::dist_sq;
  using vector_details::distance;
  using vector_details::sum_sq;
  using vector_details::apply;
  using vector_details::abs_diff;
  using vector_details::normalised;
  using vector_details::operator+;
  using vector_details::operator-;
  using vector_details::operator*;
  using vector_details::operator/;
  using vector_details::operator^;
  using vector_details::range;
} // namespace WDutils
//
#endif
