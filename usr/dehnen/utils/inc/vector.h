// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file    utils/inc/vector.h
///
/// \author  Walter Dehnen
///
/// \date    2012
/// 
/// \brief   contains the definition of class template WDutils::vector<int,T>
///
////////////////////////////////////////////////////////////////////////////////
///
/// \version jun-2012  implemented, based on tupel.h, fully C++11
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
////////////////////////////////////////////////////////////////////////////////
namespace WDutils {
  /// namespace for auxiliary functionality for vector<>
  namespace Vector {
    ///
    /// auxiliary class template for loop unrolling in vector<>
    ///
    template<int N1, typename T, int I=0> struct unroll
    {
      static_assert(N1>=I,"logic error");
      typedef unroll<N1,T,I+1> next;
      // output
      static std::ostream&out(std::ostream&s, const T*x, unsigned w, unsigned p)
      {
	s << x[I] << ' ';
	s.width(w);
	s.precision(p);
	return next::out(s,x,w,p);
      }
      // unary
      template<typename UnaryFunc>
      static void unary(T*x, UnaryFunc f)
      { f(x[I]); next::unary(x,f); }
      // collect unary
      template<typename Operator, typename UnaryFunc>
      static auto coll_unary(const T*x, Operator op, UnaryFunc f)
	-> decltype(f(x[I]))
      { return op(f(x[I]), next::coll_unary(x,op,f)); }
      // binary
      template<typename Y, typename BinaryFunc>
      static void binary(T*x, const Y*y, BinaryFunc f)
      { f(x[I],y[I]); next::binary(x,y,f); }
      // collect binary
      template<typename Operator, typename BinaryFunc>
      static auto coll_binary(const T*x, const T*y, Operator op, BinaryFunc f)
	-> decltype(f(x[I],y[I]))
      { return op(f(x[I],y[I]), next::coll_binary(x,y,op,f)); }
      // tertiary
      template<typename TertiaryFunc>
      static void tertiary(T*x, const T*y, const T*z, TertiaryFunc f)
      { f(x[I],y[I],z[I]); next::tertiary(x,y,z,f); }
      // constant tertiary
      template<typename TertiaryFunc>
      static void const_tertiary(const T*x, T*y, T*z, TertiaryFunc f)
      { f(x[I],y[I],z[I]); next::const_tertiary(x,y,z,f); }
    };
    /// closing the loop: last element
    template<int I, typename T> struct unroll<I,T,I>
    {
      // output
      static std::ostream&out(std::ostream&s, const T*x, unsigned, unsigned)
      { return s << x[I]; }
      // unary
      template<typename UnaryFunc>
      static void unary(T*x, UnaryFunc f)
      { f(x[I]); }
      // collect unary
      template<typename Operator, typename UnaryFunc>
      static auto coll_unary(const T*x, Operator, UnaryFunc f)
	-> decltype(f(x[I]))
      { return f(x[I]); }
      // binary
      template<typename Y, typename BinaryFunc>
      static void binary(T*x, const Y*y, BinaryFunc f)
      { f(x[I],y[I]); }
      // collect binary
      template<typename Operator, typename BinaryFunc>
      static auto coll_binary(const T*x, const T*y, Operator, BinaryFunc f)
	-> decltype(f(x[I],y[I]))
      { return f(x[I],y[I]); }
      // tertiary
      template<typename TertiaryFunc>
      static void tertiary(T*x, const T*y, const T*z, TertiaryFunc f)
      { f(x[I],y[I],z[I]); }
      // constant tertiary
      template<typename TertiaryFunc>
      static void const_tertiary(const T*x, T*y, T*z, TertiaryFunc f)
      { f(x[I],y[I],z[I]); }
    };
    /// square of argument
    template<typename X>
    inline X square(X a) { return a*a; }
  } // namespace WDutils::Vector
  template<int, typename> class vector;
  ///
  /// array of @a N elements of type @c T
  ///
  /// \note  We minimise cross-type operations (i.e. between @c vector<3,float>
  ///        and @c vector<3,double>). If such functionality is required, one
  ///        should first raise the lower-precision operand (we do provide a
  ///        constructor and assignment from vector with another element type).
  template<int N, typename T> class vector
  {
  private:
    static_assert(N>0,"WDutils::vector<N>: N<=0");
    typedef Vector::unroll<N-1,T> unroll;
    /// data: an array of N elements of type T
    T a[N];
    //  friendship with vector of other type to allow assignment
    template <int,typename> friend class vector;
  public:
    /// type of elements 
    typedef T element_type;
    /// tensor rank or order
    static const int ORD  = 1;
    /// number of elements
    constexpr static int size()
    { return N; }
    /// \name construction and assignment
    //@{
    /// default ctor
    vector() = default;
    /// copy ctor
    vector(vector const&) = default;
    /// assignment from another vector
    vector&operator=(vector const&) = default;
    /// copy from vector of another element type
    template<typename U>
    vector(vector<N,U> const&other)
    { unroll::binary(a, other.a, [] (T&x, U y) { x=y; } ); }
    /// assignment from another vector of another element type
    template<typename U>
    vector&operator=(vector<N,U> const&other)
    {
      unroll::binary(a, other.a, [] (T&x, U y) { x=y; } );
      return*this;
    }
    /// ctor from a scalar: set all elements to scalar
    explicit vector(T s)
    { unroll::unary(a, [s] (T&x) { x=s; } ); }
    /// assignment to a scalar: set all elements to scalar
    vector&operator=(T s)
    {
      unroll::unary(a, [s] (T&x) { x=s; } );
      return*this;
    }
    /// ctor from 2 scalars for N=2
    vector(T a0, T a1)
    {
      static_assert(N==2,"vector<N>::vector(a0,a1): N!=2");
      a[0]=a0;
      a[1]=a1;
    }
    /// ctor from 3 scalars for N=3
    vector(T a0, T a1, T a2)
    {
      static_assert(N==3,"vector<N>::vector(a0,a1,a2): N!=3");
      a[0]=a0;
      a[1]=a1;
      a[2]=a2;
    }
    /// ctor from 4 scalars for N=4
    vector(T a0, T a1, T a2, T a3)
    {
      static_assert(N==4,"vector<N>::vector(a0,a1,a2,a3): N!=4");
      a[0]=a0;
      a[1]=a1;
      a[2]=a2;
      a[3]=a3;
    }
    /// ctor from 5 scalars for N=5
    vector(T a0, T a1, T a2, T a3, T a4)
    {
      static_assert(N==5,"vector<N>::vector(a0,a1,a2,a3,a4): N!=5");
      a[0]=a0;
      a[1]=a1;
      a[2]=a2;
      a[3]=a3;
      a[4]=a4;
    }
    /// ctor from 6 scalars for N=6
    vector(T a0, T a1, T a2, T a3, T a4, T a5)
    {
      static_assert(N==6,"vector<N>::vector(a0,a1,a2,a3,a4,a5): N!=6");
      a[0]=a0;
      a[1]=a1;
      a[2]=a2;
      a[3]=a3;
      a[4]=a4;
      a[5]=a5;
    }
    /// ctor from 7 scalars for N=7
    vector(T a0, T a1, T a2, T a3, T a4, T a5, T a6)
    {
      static_assert(N==7,"vector<N>::vector(a0,a1,a2,a3,a4,a5,a6): N!=7");
      a[0]=a0;
      a[1]=a1;
      a[2]=a2;
      a[3]=a3;
      a[4]=a4;
      a[5]=a5;
      a[6]=a6;
    }
    /// ctor from 8 scalars for N=8
    vector(T a0, T a1, T a2, T a3, T a4, T a5, T a6, T a7)
    {
      static_assert(N==8,"vector<N>::vector(a0,a1,a2,a3,a4,a5,a6,a7): N!=8");
      a[0]=a0;
      a[1]=a1;
      a[2]=a2;
      a[3]=a3;
      a[4]=a4;
      a[5]=a5;
      a[6]=a6;
      a[7]=a7;
    }
    /// ctor from 9 scalars for N=9
    vector(T a0, T a1, T a2, T a3, T a4, T a5, T a6, T a7, T a8)
    {
      static_assert(N==9,"vector<N>::vector(a0,a1,a2,a3,a4,a5,a6,a7,a8): N!=9");
      a[0]=a0;
      a[1]=a1;
      a[2]=a2;
      a[3]=a3;
      a[4]=a4;
      a[5]=a5;
      a[6]=a6;
      a[7]=a7;
      a[8]=a8;
    }
    /// ctor from 10 scalars for N=10
    vector(T a0, T a1, T a2, T a3, T a4, T a5, T a6, T a7, T a8, T a9)
    {
      static_assert(N==10,
		    "vector<N>::vector(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9): N!=10");
      a[0]=a0;
      a[1]=a1;
      a[2]=a2;
      a[3]=a3;
      a[4]=a4;
      a[5]=a5;
      a[6]=a6;
      a[7]=a7;
      a[8]=a8;
      a[9]=a9;
    }
    //@}
    /// \name element access
    //@{
    /// first element
    /// \note used as public const member in decltype()
    T const&first() const
    { return a[0]; }
    /// const element access
    T const&operator[] (int i) const
    { return a[i]; }
    /// non-const element access
    T&operator[] (int i)
    { return a[i]; }
    /// conversion to pointer to scalar
    operator T*() { return a; }
    /// conversion to const pointer to scalar
    operator const T*() const { return a; }
    /// conversion to pointer to scalar
    T*data() { return a; }
    /// conversion to const pointer to scalar
    const T*data() const { return a; }
    //@}
    /// \name unary operations
    //@{
    /// change sign for each element
    vector&negate()
    { 
      unroll::unary(a, [] (T&x) { x=-x; });
      return*this;
    }
    /// return negative of this
    vector operator-() const
    {
      vector v;
      unroll::binary(v.a, a, [] (T&x, T y) { x=-y; });
      return v;
    }
    /// norm: sum of elements squared
    T norm() const
    {
      return unroll::coll_unary(a,
				[] (T x, T y)->T { return x+y; },
				[] (T x)     ->T { return x*x; });
    }
    /// abs: sqrt(norm)
    /// \note returns appropriate type, e.g. @c double for @a abs(vect<N,int>)
    auto abs() const -> decltype(std::sqrt(this->first()))
    { return std::sqrt(norm()); }
    /// minimum element
    T min() const
    {
      return unroll::coll_unary(a,
				[] (T x, T y)->T { return x<y?x:y; },
				[] (T x)     ->T { return x; });
    }
    /// minimum |element|
    T minnorm() const
    {
      return unroll::coll_unary(a,
				[] (T x, T y)->T { return x<y?x:y; },
				[] (T x)     ->T { return std::abs(x); });
    }
    /// maximum element
    T max() const
    {
      return unroll::coll_unary(a,
				[] (T x, T y)->T { return x>y?x:y; },
				[] (T x)     ->T { return x; });
    }
    /// maximum |element|
    T maxnorm() const
    {
      return unroll::coll_unary(a,
				[] (T x, T y)->T { return x>y?x:y; },
				[] (T x)     ->T { return std::abs(x); });
    }
    /// product of elements
    T volume() const
    {
      return unroll::coll_unary(a,
				[] (T x, T y)->T { return x*y; },
				[] (T x)     ->T { return x; });
    }
    /// is any element NaN?
    bool isnan() const
    { 
      return
	unroll::coll_unary(a,
			   [] (bool x, bool y)->bool { return x || y; },
			   [] (T x)           ->bool { return std::isnan(x); });
    }
    /// is any element inf?
    bool isinf() const
    { 
      return
	unroll::coll_unary(a,
			   [] (bool x, bool y)->bool { return x || y; },
			   [] (T x)           ->bool { return std::isinf(x); });
    }
    //@}
    /// \name binary operations with scalar
    //@{
    /// multiply all elements by scalar
    vector& operator*=(T s)
    { unroll::unary(a, [s] (T&x) { x*=s; } ); return*this; }
    /// product with scalar
    vector operator*(T s) const
    {
      vector r;
      unroll::binary(r.a, a, [s] (T&z, T x) { z=x*s; } );
      return r;
    }
    /// divide all elements by scalar
    vector&operator/= (T s)
    { return operator*=(T(1)/s); }
    /// division by scalar
    vector operator/ (T s) const
    { return operator*(T(1)/s); }
    /// equality
    /// \note you should not compare floating point numbers for equality
    bool operator==(T s) const
    {
      return unroll::coll_unary(a,
				[]  (bool x, bool y)->bool { return x&&y; },
				[s] (T x) { return x==s; } );
    }
    /// inequality
    /// \note you should not compare floating point numbers for inequality
    bool operator!=(T s) const
    {
      return unroll::coll_unary(a,
				[]  (bool x, bool y)->bool { return x||y; },
				[s] (T x) { return x!=s; } );
    }
    //@}
    /// \name binary operations with other vector
    //@{
    /// equality
    /// \note you should not compare floating point numbers for equality
    bool operator==(vector const&other) const
    {
      return unroll::coll_binary(a,other.a,
				 [] (bool x, bool y)->bool { return x&&y; },
				 [] (T x, T y) { return x==y; } );
    }
    /// inequality
    /// \note you should not compare floating point numbers for inequality
    bool operator!=(vector const&other) const
    {
      return unroll::coll_binary(a,other.a,
				 [] (bool x, bool y)->bool { return x||y; },
				 [] (T x, T y) { return x!=y; } );
    }
    /// add another vector
    vector&operator+=(vector const&other)
    {
      unroll::binary(a, other.a, [] (T&x, T y) { x+=y; } );
      return*this;
    }
    /// subtract another vector
    vector&operator-=(vector const&other)
    {
      unroll::binary(a, other.a, [] (T&x, T y) { x-=y; } );
      return*this;
    }
    /// element wise *=
    vector&multiply_element_wise_with(vector const&other)
    {
      unroll::binary(a, other.a, [] (T&x, T y) { x*=y; } );
      return*this;
    }
    /// element wise /=
    vector&divide_element_wise_by(vector const&other)
    {
      unroll::binary(a, other.a, [] (T&x, T y) { x/=y; } );
      return*this;
    }
    /// vector sum
    vector operator+(vector const&other) const
    {
      vector r;
      unroll::tertiary(r.a,a,other.a, [] (T&s, T x, T y) { s = x+y; } );
      return r;
    }
    /// vector difference
    vector operator-(vector const&other) const
    {
      vector r;
      unroll::tertiary(r.a,a,other.a, [] (T&d, T x, T y) { d = x-y; } );
      return r;
    }
    /// absolute difference: return |*this - other| element wise
    vector abs_diff(vector const&other) const
    {
      vector r;
      unroll::tertiary(r.a,a,other.a, [] (T&d, T x, T y) { d=std::abs(x-y); } );
      return r;
    }
    /// vector difference squared
    T dist_sq(vector const&other) const
    {
      return
	unroll::coll_binary(a, other.a,
			    [] (T x, T y)->T { return x+y; },
			    [] (T x, T y)->T { return Vector::square(x-y); } );
    }
    /// vector sum squared
    T sum_sq(vector const&other) const
    {
      return
	unroll::coll_binary(a, other.a,
			    [] (T x, T y)->T { return x+y; },
			    [] (T x, T y)->T { return Vector::square(x+y); } );
    }
    /// vector dot product
    T operator*(vector const&other) const
    {
      return unroll::coll_binary(a, other.a,
				 [] (T x, T y)->T { return x+y; },
				 [] (T x, T y)->T { return x*y; } );
    }
    /// update maximum element-wise: a[i] = max(a[i],other[i])
    vector&up_max(vector const&other)
    {
      unroll::binary(a, other.a, [] (T&x, T y) { if(y>x) x=y; } );
      return*this;
    }
    /// update minimum element-wise: a[i] = max(a[i],other[i])
    vector&up_min(vector const&other)
    {
      unroll::binary(a, other.a, [] (T&x, T y) { if(y<x) x=y; } );
      return*this;
    }
    /// update minimum and maximum element-wise:
    ///  Min[i] = min(Min[i], (*this)[i]); and
    ///  Max[i] = max(Max[i], (*this)[i]);
    void up_min_max(vector&Min, vector&Max) const
    {
      unroll::const_tertiary(a,Min.a, Max.a,
			     [] (T x, T&mi, T&ma)
			     {
			       if     (x<mi) mi=x;
			       else if(x>ma) ma=x;
			     } );
    }
    //@}
    /// \name miscellaneous
    //@{
    /// apply unary function element wise: vec[i] = f(vec[i])
    template<typename UnaryFunction>
    vector&apply(UnaryFunction f)
    {
      unroll::unary(a, [f] (T&x) { x=f(x); } );
      return*this;
    }
    /// apply unary function element wise: ret[i] = f(vec[i])
    template<typename UnaryFunction>
    auto applied(UnaryFunction f) const
      -> vector<N,decltype(f(T{}))>
    {
      typedef decltype(f(T{})) Y;
      vector<N,Y> r;
      unroll::binary(r.a, a, [f] (Y&y, T x) { y=f(x); } );
      return r;
    }
    /// connect with another vector element wise: vec[i] = f(vec[i], other[i])
    template<typename BinaryFunction>
    vector&apply_binary(vector const&other, BinaryFunction f)
    {
      unroll::binary(a, other.a, [f] (T&x, T y) { x=f(x,y); } );
      return*this;
    }
    /// connect with another vector element wise: ret[i] = f(vec[i], other[i])
    template<typename BinaryFunction>
    auto applied_binary(vector const&other, BinaryFunction f) const
      -> vector<N,decltype(f(T{},T{}))>
    {
      typedef decltype(f(T{},T{})) Y;
      vector<N,Y> r;
      unroll::tertiary(r.a, a, other.a, [f] (Y&y, T x, T z) { y=f(x,z); } );
      return r;
    }
    /// normalize: vec[i] /= abs(vec)
    vector&normalize()
    {
      T n = norm();
      if(n) return operator /= (std::sqrt(n));
      else  return*this;
    }
    /// normalized: ret[i] = vec[i] / abs(vec)
    vector normalized() const
    {
      T n = norm();
      if(n) return this->operator/ (std::sqrt(n));
      else  return*this;
    }
    /// write vector element-wise to std::ostream; keeps width and precision
    friend std::ostream&operator<<(std::ostream&s, vector const&vec)
    {
      unsigned w = s.width();
      unsigned p = s.precision();
      return unroll::out(s,vec.a,w,p);
    }
    /// read vector element-wise from std::istream
    friend std::istream&operator>>(std::istream&s, vector&vec)
    {
      unroll::unary(vec.a,[&s] (T&x) { s >> x; } );
      return s;
    }
    //@}
  };// class WDutils::vector<N,T>
  /// \name properties of a vector
  //@{
  /// norm: sum of vector elements squared
  /// \relates WDutils::vector
  template<int N, typename T>
  inline T norm(vector<N,T> const&vec)
  { return vec.norm(); }
  /// absolute value: sqrt(norm(vector))
  /// \relates WDutils::vector
  template<int N, typename T>
  inline auto abs(vector<N,T> const&vec) -> decltype(vec.abs())
  { return vec.abs(); }
  /// product of elements
  /// \relates WDutils::vector
  template<int N, typename T>
  inline T volume(vector<N,T> const&vec)
  { return vec.volume(); }
  /// minimum element of vector
  /// \relates WDutils::vector
  template<int N, typename T>
  inline T min(vector<N,T> const&vec)
  { return vec.min(); }
  /// minimum |element| of vector
  /// \relates WDutils::vector
  template<int N, typename X> inline
  X minnorm(vector<N,X> const&x)
  { return x.minnorm(); }
  /// maximum element of vector
  /// \relates WDutils::vector
  template<int N, typename T>
  inline T max(vector<N,T> const&vec)
  { return vec.max(); }
  /// maximum |element| of vector
  /// \relates WDutils::vector
  template<int N, typename X> inline
  X maxnorm(vector<N,X> const&x)
  { return x.maxnorm(); }
  /// is any element NaN?
  /// \relates WDutils::vector
  template<int N, typename X> inline
  bool isnan(vector<N,X> const&x)
  { return x.isnan(); }
  /// is any element inf?
  /// \relates WDutils::vector
  template<int N, typename X> inline
  bool isinf(vector<N,X> const&x)
  { return x.isinf(); }
  /// returns @a x/|x|, the unit vector in direction @a x
  /// \relates WDutils::vector
  template<int N, typename X> inline
  vector<N,X> normalized(vector<N,X> const&x)
  { return x.normalized(); }
  //@}
  /// \name binary operations with vectors
  //@{
  /// scalar times vector
  /// \relates WDutils::vector
  template<int N, typename T>
  inline vector<N,T> operator*(T s, vector<N,T> const&vec)
  { return vec*s; }
  /// return absolute difference between two vectors
  /// \relates WDutils::vector
  template<int N, typename X>
  inline vector<N,X> abs_diff(vector<N,X> const&x, vector<N,X> const&y)
  { return x.abs_diff(y); }
  /// vector difference squared
  /// \relates WDutils::vector
  template<int N, typename T>
  inline T dist_sq(vector<N,T> const&x, vector<N,T> const&y)
  { return x.dist_sq(y); }
  /// vector sum squared
  /// \relates WDutils::vector
  template<int N, typename T>
  inline T sum_sq(vector<N,T> const&x, vector<N,T> const&y)
  { return x.sum_sq(y); }
  /// vector distance = sqrt( norm(vector difference) )
  /// \relates WDutils::vector
  template<int N, typename T>
  inline auto distance(vector<N,T> const&x, vector<N,T> const&y)
    -> decltype(std::sqrt(dist_sq(x,y)))
  { return std::sqrt(dist_sq(x,y)); }
  /// vector cross product for N=2: returns x[0]*y[1]- x[1]*y[0]
  /// \relates WDutils::vector
  template<typename T>
  inline T operator^ (vector<2,T> const&x, vector<2,T> const&y)
  { return x[0]*y[1] - x[1]*y[0]; }
  /// vector cross product for N=3
  /// \relates WDutils::vector
  template<typename T>
  inline vector<3,T> operator^ (vector<3,T> const&x, vector<3,T> const&y)
  {
    return vector<3,T>(x[1]*y[2] - x[2]*y[1],
		       x[2]*y[0] - x[0]*y[2],
		       x[0]*y[1] - x[1]*y[0]);
  }
  //@}
} // namespace WDutils
//
#endif // WDutils_included_vector_h
