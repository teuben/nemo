// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// tupel.h                                                                     |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 1996-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef included_tupel_h
#define included_tupel_h

#include <iostream>
#include <cmath>

namespace WD {

  template<int N, typename X>
  class tupel {

  public:
    // 0   static members and types
    typedef X element_type;
    static const int NDAT = N;
    static int size()
    {
      return N;
    }

  private:
    //    data: array of N elements
    //    essentially, tupel<N,X> is just a (very efficient) wrapper of this
    element_type a[N];

  public:
    // 1   constructors
    // 1.1 construction from single argument (scalar, array, or tupel)
    tupel()
    {
    }
    inline tupel(element_type const&);
    inline tupel(const element_type*);
    inline tupel(tupel const&);
#ifdef TUPEL_CROSS_TYPE
    template<typename scalar_type> tupel(tupel<N,scalar_type> const&);
#endif
    // 1.2 construction from list of elements (up to 10)
    inline tupel(element_type const&,
		 element_type const&);
    inline tupel(element_type const&,
		 element_type const&,
		 element_type const&);
    inline tupel(element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&);
    inline tupel(element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&);
    inline tupel(element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&);
    inline tupel(element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&);
    inline tupel(element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&);
    inline tupel(element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&);
    inline tupel(element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&,
		 element_type const&);
    // 2   element access
    inline element_type      &operator[]  (int i)
    {
      return a[i];
    }
    
    inline element_type const&operator[]  (int i) const
    {
      return a[i];
    }

    // 3   unitary operators
    // 3.1 negation and unatary minus
    inline tupel&negate();
    inline tupel operator-() const;
    // 3.2 minimum and maximum element, as well as maximum(abs(element))
    inline element_type min() const;
    inline element_type max() const;
    // 3.3 type conversion to pointer to element_type
    inline operator element_type* ()
    {
      return a;
    }

    inline operator const element_type* () const
    {
      return a;
    }

    // 4   binary operators with scalar
    // 4.1 with scalar of element_type
    inline tupel&operator=  (element_type const&);
    inline tupel&operator*= (element_type const&);
    inline tupel operator*  (element_type const&) const;
    inline tupel&operator/= (element_type const&);
    inline tupel operator/  (element_type const&) const;
    inline bool  operator== (element_type const&) const;
    inline bool  operator!= (element_type const&) const;
#ifdef TUPEL_CROSS_TYPE
    // 4.2 with scalar of type different from element_type
    template<typename scalar_type>
    inline tupel&operator= (scalar_type const&);
    template<typename scalar_type>
    inline tupel&operator*=(scalar_type const&);
    template<typename scalar_type>
    inline tupel operator* (scalar_type const&) const;
    template<typename scalar_type>
    inline tupel&operator/=(scalar_type const&);
    template<typename scalar_type>
    inline tupel operator/ (scalar_type const&) const;
#endif

    // 5   binary operators with tupel
    // 5.1 with tupel<N,element_type>
    inline tupel& operator=  (tupel const&);
    inline tupel& operator+= (tupel const&);
    inline tupel& operator-= (tupel const&);
    inline bool   operator== (tupel const&) const;
    inline bool   operator!= (tupel const&) const;
    inline tupel  operator+  (tupel const&) const;
    inline tupel  operator-  (tupel const&) const;
    //     add or subtract another tupel integer times
    template<int> inline tupel& add_int_times(tupel const&);
    template<int> inline tupel& sub_int_times(tupel const&);
    //     update the minimum/maximum in each element: a[i] = max(a[i],b[i])
    inline tupel& update_max (tupel const&);
    inline tupel& update_min (tupel const&);
#ifdef TUPEL_CROSS_TYPE
    // 5.2 with tupel<N,scalar_type>
    template<typename scalar_type>
    inline tupel& operator=  (tupel<N,scalar_type> const&);
    template<typename scalar_type>
    inline tupel& operator+= (tupel<N,scalar_type> const&);
    template<typename scalar_type>
    inline tupel& operator-= (tupel<N,scalar_type> const&);
    template<typename scalar_type>
    inline tupel  operator+  (tupel<N,scalar_type> const&) const;
    template<typename scalar_type>
    inline tupel  operator-  (tupel<N,scalar_type> const&) const;
    //     add or subtract another tupel integer times
    template<int, typename scalar_type>
    inline tupel& add_int_times(tupel<N,scalar_type> const&);
    template<int, typename scalar_type>
    inline tupel& sub_int_times(tupel<N,scalar_type> const&);
#endif

    // 6   miscellaneous
    // 6.1 copy from array (for security, operator= is not provided)
    inline tupel&copy(const element_type*);
#ifdef TUPEL_CROSS_TYPE
    template<typename scalar_type>
    inline tupel&copy(const scalar_type*);
#endif
    // 6.2 assign, add to, or subtract from: real times another tupel
    template<typename scalar_type>
    inline tupel&ass_times(tupel const&, scalar_type const&);
    template<typename scalar_type>
    inline tupel&add_times(tupel const&, scalar_type const&);
    template<typename scalar_type>
    inline tupel&sub_times(tupel const&, scalar_type const&);
    
  }; // class tupel<N,X>

  // 7  non-member methods that use tupel<>s
  // 7.1  formatted I/O
  template<int N, typename X> inline
  std::ostream&operator<< (std::ostream&s,tupel<N,X> const&);

  template<int N, typename X> inline
  std::istream&operator>> (std::istream&s,tupel<N,X>&);

  // 7.2 norm and absolute value: norm == sum_i x[i]*x[i]
  template<int N, typename X> inline
  X norm(tupel<N,X> const&);

  template<int N, typename X> inline
  X abs (tupel<N,X> const&);

  // 7.3 maximum(abs(element))
  template<int N, typename X> inline
  X maxnorm(tupel<N,X> const&);

  // 7.4 distance, squared distance, and squared sum 
  template<int N, typename X> inline
  X dist_sq (tupel<N,X> const&, tupel<N,X> const&);
  template<int N, typename X> inline
  X dist    (tupel<N,X> const&, tupel<N,X> const&);
  template<int N, typename X> inline
  X sum_sq  (tupel<N,X> const&, tupel<N,X> const&);
#ifdef TUPEL_CROSS_TYPE
  template<int N, typename X, typename Y> inline
  X dist_sq (tupel<N,X> const&, tupel<N,Y> const&);
  template<int N, typename X, typename Y> inline
  X dist    (tupel<N,X> const&, tupel<N,Y> const&);
  template<int N, typename X, typename Y> inline
  X sum_sq  (tupel<N,X> const&, tupel<N,Y> const&);
#endif

  // 7.5 vector dot product
  template<int N, typename X> inline
  X operator* (tupel<N,X> const&, tupel<N,X> const&);
#ifdef TUPEL_CROSS_TYPE
  template<int N, typename X, typename Y> inline
  X operator* (tupel<N,X> const&, tupel<N,Y> const&);
#endif

  // 7.6 vector cross product in 2D and 3D, coded as operator ^
  //     in 2D, we just return the scalar  a[0]*b[1]-a[1]*b[0],
  //     so that X^V gives Lz and (Lx,Ly,Lz) in 2D and 3D, respectively
  template<typename X> inline
  X operator^ (tupel<2,X> const&, tupel<2,X> const&);
  template<typename X> inline
  tupel<3,X> operator^ (tupel<3,X> const&, tupel<3,X> const&);
#ifdef TUPEL_CROSS_TYPE
  template<typename X, typename Y> inline
  X operator^ (tupel<2,X> const&, tupel<2,Y> const&);
  template<typename X, typename Y> inline
  tupel<3,X> operator^ (tupel<3,X> const&, tupel<3,Y> const&);
#endif

  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // What follows is the implementation of the declarations above.            //
  //                                                                          //
  // The implementation uses the template metaprogramming technique, which    //
  // makes the compiler to automatically unroll loops over the tupel index.   //
  // This make the implementation substantially faster than a plain naive     //
  // for loop over the index.                                                 //
  //////////////////////////////////////////////////////////////////////////////

  // in order to insulate the implementation details, we put them in a 
  // separate namespace
  namespace meta {

    // A.1  times<K>(x) returns K*x;
    //    optimized at compile time for K=-2,-1,0,1,2
    //    template struct times__<> to implement template method times<> below.
    template<int N, typename X>
    struct times__
    {
      static X act(X const&x)
      {
	return N * x;
      }
    };
    template<typename X>
    struct times__<2,X>
    {
      static X act(X const&x)
      {
	return x+x;
      }
    };
    template<typename X>
    struct times__<1,X>
    {
      static X act(X const&x)
      {
	return x;
      }
    };
    template<typename X>
    struct times__<0,X>
    {
      static X act(X const&x)
      {
	return X(0);
      }
    };
    template<typename X>
    struct times__<-1,X>
    {
      static X act(X const&x)
      {
	return -x;
      }
    };
    template<typename X>
    struct times__<-2,X>
    {
      static X act(X const&x)
      {
	return -x-x;
      }
    };

    //    implement times<>() using struct times__<>
    template<int N, typename X> inline
    X times(X const&x)
    { 
      return times__<N,X>::act(x);
    }

    // A.2  template method square(x) returns x*x (better than macro)
    template<typename X> inline
    X square(X const&x)
    {
      return x*x;
    }

    // A.3  abs<X>(x) to return |x|, optimized for float and double at compile
    //      template struct abs__ to implement abs<>()
    template<typename X> struct abs__
    {
      X act(X const&x)
      {
	return x < X(0) ? -x : x;
      }
    };
    template<> struct abs__<float>
    {
      float act(float const&x)
      {
	return fabs(x);
      }
    };
    template<> struct abs__<double>
    {
      double act(double const&x)
      {
	return fabs(x);
      }
    };
    
    //   implement abs<>(x) using struct abs__<>
    template<typename X> inline
    X abs(X const&x)
    {
      return abs__<X>::act(x);
    }
    

    // A.4  min and max from abritrary type allowing for operators < and >
    template<typename X> inline
    X min(X const&x,
	  X const&y)
    {
      return x<y? x : y;
    }

    template<typename X> inline
    X max(X const&x,
	  X const&y)
    {
      return x>y? x : y;
    }

    // A.5  class aux_tupel<N,X,I>
    //      its static methods will be used below to implement the member
    //      methods of class tupel<N,X>
    template<int N, typename X, int I=0> 
    struct aux_tupel
    {
      typedef X element_type;
      friend class tupel<N+1,element_type>;

      // typedef next to be the aux_tupel for the next element
      typedef aux_tupel<N,element_type,I+1> next;

      // s_as():  assign array elements to scalar
      template<typename scalar_type>
      static void s_as(element_type     *a,
		       const scalar_type&x)
      {
	a[I]=x;
	next::s_as(a,x);
      }

      // s_ml():  multiply array elements with scalar
      template<typename scalar_type>
      static void s_ml(element_type     *a,
		       const scalar_type&x)
      {
	a[I]*=x;
	next::s_ml(a,x);
      }

      // v_as():  assign array elements to array elements
      template<typename scalar_type>
      static void v_as(element_type     *a,
		       const scalar_type*b)
      {
	a[I]=b[I];
	next::v_as(a,b);
      }

      // v_ad():  add array elements to array elements
      template<typename scalar_type>
      static void v_ad(element_type     *a,
		       const scalar_type*b)
      {
	a[I]+=b[I];
	next::v_ad(a,b);
      }

      // v_su():  subtract array elements from array elements
      template<typename scalar_type>
      static void v_su(element_type     *a,
		       const scalar_type*b)
      {
	a[I]-=b[I];
	next::v_su(a,b);
      }

      // v_ast():  assign array elements to x times another array's elements
      template<typename scalar_type>
      static void v_ast(element_type      *a,
			const element_type*b,
			scalar_type  const&x)
      {
	a[I]=x*b[I];
	next::v_ast(a,b,x);
      }

      // v_adt():  add to array elements x times another array's elements
      template<typename scalar_type>
      static void v_adt(element_type      *a,
			const element_type*b,
		        scalar_type  const&x)
      {
	a[I]+=x*b[I];
	next::v_adt(a,b,x);
      }

      // v_adit():  add to array elements i times another array's elements
      template<int K>
      static void v_adit(element_type      *a,
			 const element_type*b)
      {
	a[I]+=times<K>(b[I]);
	next:: template v_adit<K>(a,b,x);
      }

      // v_sut():  subtract from array elements x times another array's elements
      template<typename scalar_type> 
      static void v_sut (element_type      *a,
			 const element_type*b,
			 scalar_type  const&x)
      {
	a[I]-=x*b[I];
	next::v_sut(a,b,x);
      }

      // v_suit(): subtract from array elements i times another array's elements
      template<int K> 
      static void v_suit (element_type      *a,
			  const element_type*b)
      {
	a[I]-=times<K>(b[I]);
	next:: template v_suit<K>(a,b,x);
      }

      // v_neg():  negate array elements
      static void v_neg (element_type*a)
      {
	a[I]=-a[I];
	next::v_neg(a);
      }
      
      // v_nega():  set array elements to minus another array's elements
      template<typename scalar_type>
      static void v_nega(element_type     *a,
			 const scalar_type*b)
      {
	a[I]=-b[I];
	next::v_nega(a,b);
      }

      // v_sum():  set array elements to sum of two other array's elements
      template<typename scalar_type> 
      static void v_sum(element_type      *a,
			const element_type*b,
			const scalar_type *c)
      {
	a[I]=b[I]+c[I];
	next::v_sum(a,b,c);
      }
      
      // v_diff():  set array elements to difference between
      //            two other array's elements
      template<typename scalar_type> 
      static void v_dif (element_type      *a,
			 const element_type*b,
			 const scalar_type *c)
      {
	a[I]=b[I]-c[I];
	next::v_dif(a,b,c);
      }

      // s_eq():  compare array elements with a scalar for equality
      static bool s_eq(const element_type*a,
		       const element_type&b)
      { 
	return a[I]==b && next::s_eq(a,b);
      }

      // s_neq():  compare array elements with a scalar for inequality
      static bool s_neq(const element_type*a,
			const element_type&b)
      { 
	return a[I]!=b || next::s_neq(a,b);
      }

      // v_eq():  compare two array's elements for equality
      static bool v_eq(const element_type*a,
		       const element_type*b)
      { 
	return a[I]==b[I] && next::v_eq(a,b);
      }

      // v_neq():  compare two array's elements for equality
      static bool v_neq(const element_type*a,
			const element_type*b)
      { 
	return a[I]!=b[I] || next::v_neq(a,b);
      }

      // v_min(): return mininum of elements of an array
      static element_type v_min(const element_type*a)
      {
	return min( a[I], next::v_min(a) );
      }

      // v_max(): return maxinum of elements of an array
      static element_type v_max(const element_type*a)
      {
	return max( a[I], next::v_max(a) );
      }

      // v_uma(): update maximum: a[I] = max(a[i],b[i])
      static void v_uma (element_type      *a,
			 const element_type*b)
      {
	if(b[I]>a[I]) a[I] = b[I];
	next::v_uma(a,b);
      }

      // v_umi(): update minimum: a[I] = min(a[i],b[i])
      static void v_umi (element_type      *a,
			 const element_type*b)
      {
	if(b[I]<a[I]) a[I] = b[I];
	next::v_umi(a,b);
      }

      // methods to be used by non-members of tupel<>
      // v_norm(): add squares of elements of array
      static element_type v_norm(const element_type*a)
      {
	return a[I]*a[I] + next::v_norm(a);
      }

      // v_amax(): return maxinum of absolute value of elements of an array
      static element_type v_amax(const element_type*a)
      {
	return max( abs(a[I]), next::v_amax(a) );
      }

      // v_dot(): add products of elements of two arrays
      template<typename scalar_type>
      static element_type v_dot(const element_type*a,
				const scalar_type *b)
      {
	return a[I]*b[I] + next::v_dot(a,b);
      }
      
      // v_diq(): add square of differences between elements of two arrays
      template<typename scalar_type>
      static element_type v_diq(const element_type*a,
				const scalar_type *b)
      {
	return square(a[I]-b[I]) + next::v_diq(a,b);
      }

      // v_suq(): add square of sum of elements of two arrays
      template<typename scalar_type>
      static element_type v_suq(const element_type*a,
				const scalar_type *b)
      {
	return square(a[I]+b[I]) + next::v_suq(a,b);
      }

      // v_out(): formatted output of elements of array
      static void v_out(std::ostream      &o,
			const element_type*a)
      {
	o<<a[I]<<' ';
	next::v_out(o,a);
      }

      // v_in(): ascii input of elements of array
      static void v_in(std::istream&i,
		       element_type*a)
      {
	i>>a[I];
	next::v_in(i,a);
      }

    }; // class aux_tupel<N,X,I>

    //    class aux_tupel<I,X,I>
    //    special case of I=N: last element
    template<int I, typename X> struct aux_tupel<I,X,I>
    {
      typedef X element_type;
      friend class tupel<I+1,element_type>;

      // s_as():  assign array elements to scalar
      template<typename scalar_type>
      static void s_as(element_type     *a,
		       const scalar_type&x)
      {
	a[I]=x;
      }

      // s_ml():  multiply array elements with scalar
      template<typename scalar_type>
      static void s_ml(element_type     *a,
		       const scalar_type&x)
      {
	a[I]*=x;
      }

      // v_as():  assign array elements to array elements
      template<typename scalar_type>
      static void v_as(element_type     *a,
		       const scalar_type*b)
      {
	a[I]=b[I];
      }

      // v_ad():  add array elements to array elements
      template<typename scalar_type>
      static void v_ad(element_type     *a,
		       const scalar_type*b)
      {
	a[I]+=b[I];
      }

      // v_su():  subtract array elements from array elements
      template<typename scalar_type>
      static void v_su(element_type     *a,
		       const scalar_type*b)
      {
	a[I]-=b[I];
      }

      // v_ast():  assign array elements to x times another array's elements
      template<typename scalar_type>
      static void v_ast(element_type      *a,
			const element_type*b,
			scalar_type  const&x)
      {
	a[I]=x*b[I];
      }

      // v_adt():  add to array elements x times another array's elements
      template<typename scalar_type>
      static void v_adt(element_type      *a,
			const element_type*b,
		        scalar_type  const&x)
      {
	a[I]+=x*b[I];
      }

      // v_adit():  add to array elements i times another array's elements
      template<int K>
      static void v_adit(element_type      *a,
			 const element_type*b)
      {
	a[I]+=times<K>(b[I]);
      }

      // v_sut():  subtract from array elements x times another array's elements
      template<typename scalar_type> 
      static void v_sut (element_type      *a,
			 const element_type*b,
			 scalar_type  const&x)
      {
	a[I]-=x*b[I];
      }

      // v_suit(): subtract from array elements i times another array's elements
      template<int K> 
      static void v_suit (element_type      *a,
			  const element_type*b)
      {
	a[I]-=times<K>(b[I]);
      }

      // v_neg():  negate array elements
      static void v_neg (element_type*a)
      {
	a[I]=-a[I];
      }
      
      // v_nega():  set array elements to minus another array's elements
      template<typename scalar_type>
      static void v_nega(element_type     *a,
			 const scalar_type*b)
      {
	a[I]=-b[I];
      }

      // v_sum():  set array elements to sum of two other array's elements
      template<typename scalar_type> 
      static void v_sum(element_type      *a,
			const element_type*b,
			const scalar_type *c)
      {
	a[I]=b[I]+c[I];
      }
      
      // v_diff():  set array elements to difference between
      //            two other array's elements
      template<typename scalar_type> 
      static void v_dif (element_type      *a,
			 const element_type*b,
			 const scalar_type *c)
      {
	a[I]=b[I]-c[I];
      }

      // s_eq():  compare array elements with a scalar for equality
      static bool s_eq(const element_type*a,
		       const element_type&b)
      { 
	return a[I]==b;
      }

      // s_neq():  compare array elements with a scalar for inequality
      static bool s_neq(const element_type*a,
			const element_type&b)
      { 
	return a[I]!=b;
      }

      // v_eq():  compare two array's elements for equality
      static bool v_eq(const element_type*a,
		       const element_type*b)
      { 
	return a[I]==b[I];
      }

      // v_neq():  compare two array's elements for equality
      static bool v_neq(const element_type*a,
			const element_type*b)
      { 
	return a[I]!=b[I];
      }

      // v_min(): return mininum of elements of an array
      static element_type v_min(const element_type*a)
      {
	return a[I];
      }

      // v_max(): return maxinum of elements of an array
      static element_type v_max(const element_type*a)
      {
	return a[I];
      }

      // v_uma(): update maximum: a[I] = max(a[i],b[i])
      static void v_uma (element_type      *a,
			 const element_type*b)
      {
	if(b[I]>a[I]) a[I] = b[I];
      }

      // v_umi(): update minimum: a[I] = min(a[i],b[i])
      static void v_umi (element_type      *a,
			 const element_type*b)
      {
	if(b[I]<a[I]) a[I] = b[I];
      }

      // methods to be used by non-members of tupel<>
      // v_amax(): return maxinum of absolute value of elements of an array
      static element_type v_amax(const element_type*a)
      {
	return abs(a[I]);
      }

      // v_norm(): add squares of elements of array
      static element_type v_norm(const element_type*a)
      {
	return a[I]*a[I];
      }

      // v_dot(): add products of elements of two arrays
      template<typename scalar_type>
      static element_type v_dot(const element_type*a,
				const scalar_type *b)
      {
	return a[I]*b[I];
      }
      
      // v_diq(): add square of differences between elements of two arrays
      template<typename scalar_type>
      static element_type v_diq(const element_type*a,
				const scalar_type *b)
      {
	return square(a[I]-b[I]);
      }


      // v_suq(): add square of sum of elements of two arrays
      template<typename scalar_type>
      static element_type v_suq(const element_type*a,
				const scalar_type *b)
      {
	return square(a[I]+b[I]);
      }

      // v_out(): formatted output of elements of array
      static void v_out(std::ostream      &o,
			const element_type*a)
      {
	o<<a[I];
      }

      // v_in(): ascii input of elements of array
      static void v_in(std::istream&i,
		       element_type*a)
      {
	i>>a[I];
      }

    }; // class aux_tupel<I,X,I>

  } // namespace meta

  // A.6  Now, finally, we can implement the members of tupel
  // A.6.1   constructors
  template<int N, typename X> inline
  tupel<N,X>::tupel(element_type const&x)
  {
    meta::aux_tupel<N-1,X>::s_as(a,x);
  }

  template<int N, typename X> inline
  tupel<N,X>::tupel(const element_type*x)
  {
    meta::aux_tupel<N-1,X>::v_as(a,x);
  }

  template<int N, typename X> inline
  tupel<N,X>::tupel(tupel const&x)
  {
    meta::aux_tupel<N-1,X>::v_as(a,x.a);
  }

#ifdef TUPEL_CROSS_TYPE
  template<int N, typename X> template<typename S>
  inline tupel<N,X>::tupel(tupel<N,S> const&x)
  {
    meta::aux_tupel<N-1,X>::v_as(a,static_cast<const S*>(x));
  }
#endif
  template<int N, typename X> inline
  tupel<N,X>::tupel(element_type const&x0,
		    element_type const&x1)
  {
    a[0]=x0;
    a[1]=x1;
  }
  template<int N, typename X> inline
  tupel<N,X>::tupel(element_type const&x0,
		    element_type const&x1,
		    element_type const&x2)
  {
    a[0]=x0;
    a[1]=x1;
    a[2]=x2;
  }
  template<int N, typename X> inline
  tupel<N,X>::tupel(element_type const&x0,
		    element_type const&x1,
		    element_type const&x2,
		    element_type const&x3)
  {
    a[0]=x0;
    a[1]=x1;
    a[2]=x2;
    a[3]=x3;
  }
  template<int N, typename X> inline
  tupel<N,X>::tupel(element_type const&x0,
		    element_type const&x1,
		    element_type const&x2,
		    element_type const&x3,
		    element_type const&x4)
  {
    a[0]=x0;
    a[1]=x1;
    a[2]=x2;
    a[3]=x3;
    a[4]=x4;
  }
  template<int N, typename X> inline
  tupel<N,X>::tupel(element_type const&x0,
		    element_type const&x1,
		    element_type const&x2,
		    element_type const&x3,
		    element_type const&x4,
		    element_type const&x5)
  { 
    a[0]=x0;
    a[1]=x1;
    a[2]=x2;
    a[3]=x3;
    a[4]=x4;
    a[5]=x5;
  }
  template<int N, typename X> inline
  tupel<N,X>::tupel(element_type const&x0,
		    element_type const&x1,
		    element_type const&x2,
		    element_type const&x3,
		    element_type const&x4,
		    element_type const&x5,
		    element_type const&x6)
  {
    a[0]=x0;
    a[1]=x1;
    a[2]=x2;
    a[3]=x3;
    a[4]=x4;
    a[5]=x5;
    a[6]=x6;
  }
  template<int N, typename X> inline
  tupel<N,X>::tupel(element_type const&x0,
		    element_type const&x1,
		    element_type const&x2,
		    element_type const&x3,
		    element_type const&x4,
		    element_type const&x5,
		    element_type const&x6,
		    element_type const&x7)
  {
    a[0]=x0;
    a[1]=x1;
    a[2]=x2;
    a[3]=x3;
    a[4]=x4;
    a[5]=x5;
    a[6]=x6;
    a[7]=x7;
  }
  template<int N, typename X> inline
  tupel<N,X>::tupel(element_type const&x0,
		    element_type const&x1,
		    element_type const&x2,
		    element_type const&x3,
		    element_type const&x4,
		    element_type const&x5,
		    element_type const&x6,
		    element_type const&x7,
		    element_type const&x8)
  {
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
  template<int N, typename X> inline
  tupel<N,X>::tupel(element_type const&x0,
		    element_type const&x1,
		    element_type const&x2,
		    element_type const&x3,
		    element_type const&x4,
		    element_type const&x5,
		    element_type const&x6,
		    element_type const&x7,
		    element_type const&x8,
		    element_type const&x9)
  {
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

  // A.6.2   unitary operators
  // A.6.2.1 negation and unatary minus
  template<int N, typename X> inline
  tupel<N,X>& tupel<N,X>::negate()
  {
    meta::aux_tupel<N-1,X>::v_neg(a);
    return*this;
  }

  template<int N, typename X> inline
  tupel<N,X> tupel<N,X>::operator-() const
  {
    tupel<N,X> x;
    meta::aux_tupel<N-1,X>::v_nega(x.a,a);
    return x;
  }
  // A.6.2.2 minimum and maximum element and max(abs(element))
  template<int N, typename X> inline
  X tupel<N,X>::min() const
  {
    return meta::aux_tupel<N-1,X>::v_min(a);
  }

  template<int N, typename X> inline
  X tupel<N,X>::max() const
  {
    return meta::aux_tupel<N-1,X>::v_max(a);
  }

  // A.6.3   binary operators with scalar
  // A.6.3.1 with scalar of element_type
  template<int N, typename X> inline
  tupel<N,X>&tupel<N,X>::operator=  (element_type const&x)
  {
    meta::aux_tupel<N-1,X>::s_as(a,x);
    return*this;
  }

  template<int N, typename X> inline
  tupel<N,X>&tupel<N,X>::operator*= (element_type const&x)
  {
    meta::aux_tupel<N-1,X>::s_ml(a,x);
    return*this;
  }

  template<int N, typename X> inline
  tupel<N,X> tupel<N,X>::operator*  (element_type const&x) const
  {
    tupel<N,X> y;
    meta::aux_tupel<N-1,X>::v_ast(x.a,a,x);
    return y;
  }

  template<int N, typename X> inline
  tupel<N,X>&tupel<N,X>::operator/= (element_type const&x)
  {
    return operator*= (X(1)/x);
  }

  template<int N, typename X> inline
  tupel<N,X> tupel<N,X>::operator/  (element_type const&x) const
  {
    return operator*(X(1)/x);
  }

  template<int N, typename X> inline
  bool tupel<N,X>::operator== (element_type const&x) const
  {
    return meta::aux_tupel<N-1,X>::s_eq(a,x);
  }

  template<int N, typename X> inline
  bool tupel<N,X>::operator!= (element_type const&x) const
  {
    return meta::aux_tupel<N-1,X>::s_neq(a,x);
  }

#ifdef TUPEL_CROSS_TYPE
  // A.6.3.2 with scalar of type different from element_type
  template<int N, typename X> template<typename scalar_type> inline
  tupel<N,X>&tupel<N,X>::operator=  (scalar_type const&x)
  {
    meta::aux_tupel<N-1,X>::s_as(a,x);
    return*this;
  }

  template<int N, typename X> template<typename scalar_type> inline
  tupel<N,X>&tupel<N,X>::operator*= (scalar_type const&x)
  {
    meta::aux_tupel<N-1,X>::s_ml(a,x);
    return*this;
  }

  template<int N, typename X> template<typename scalar_type> inline
  tupel<N,X> tupel<N,X>::operator*  (scalar_type const&x) const
  {
    tupel<N,X> y;
    meta::aux_tupel<N-1,X>::v_ast(x.a,a,x);
    return y;
  }

  template<int N, typename X> template<typename scalar_type> inline
  tupel<N,X>&tupel<N,X>::operator/= (scalar_type const&x)
  {
    return operator*=(scalar_type(1)/x);
  }

  template<int N, typename X> template<typename scalar_type> inline
  tupel<N,X> tupel<N,X>::operator/  (scalar_type const&x) const
  {
    return operator*(scalar_type(1)/x);
  }
#endif

  // A.6.4   binary operators with tupel
  // A.6.4.1 with tupel<N,element_type>
  template<int N, typename X> inline
  tupel<N,X>& tupel<N,X>::operator=  (tupel<N,X> const&v)
  {
    meta::aux_tupel<N-1,X>::v_as(a,v.a);
    return*this;
  }

  template<int N, typename X> inline
  tupel<N,X>& tupel<N,X>::operator+= (tupel<N,X> const&v)
  {
    meta::aux_tupel<N-1,X>::v_ad(a,v.a);
    return*this;
  }

  template<int N, typename X> inline
  tupel<N,X>& tupel<N,X>::operator-= (tupel<N,X> const&v)
  {
    meta::aux_tupel<N-1,X>::v_su(a,v.a);
    return*this;
  }

  template<int N, typename X> inline
  bool tupel<N,X>::operator== (tupel<N,X> const&v) const
  {
    return meta::aux_tupel<N-1,X>::v_eq(a,v.a);
  }

  template<int N, typename X> inline
  bool tupel<N,X>::operator!= (tupel<N,X> const&v) const
  {
    return meta::aux_tupel<N-1,X>::v_neq(a,v.a);
  }

  template<int N, typename X> inline
  tupel<N,X>  tupel<N,X>::operator+  (tupel<N,X> const&v) const
  {
    tupel<N,X> y; 
    meta::aux_tupel<N-1,X>::v_sum(y.a,a,v.a);
    return y;
  }

  template<int N, typename X> inline
  tupel<N,X>  tupel<N,X>::operator-  (tupel<N,X> const&v) const
  {
    tupel<N,X> y; 
    meta::aux_tupel<N-1,X>::v_dif(y.a,a,v.a);
    return y;
  }

  //     add or subtract another tupel integer times
  template<int N, typename X>  template<int K> inline
  tupel<N,X>& tupel<N,X>::add_int_times(tupel<N,X> const&v)
  {
    meta::aux_tupel<N-1,X>:: template v_adit<K>(a,v.a);
    return*this;
  }

  template<int N, typename X>  template<int K> inline
  tupel<N,X>& tupel<N,X>::sub_int_times(tupel<N,X> const&v)
  {
    meta::aux_tupel<N-1,X>:: template v_suit<K>(a,v.a);
    return*this;
  }

  //     update the minimum/maximum in each element: a[i] = max(a[i],b[i])
  template<int N, typename X> inline
  tupel<N,X>& tupel<N,X>::update_max (tupel<N,X> const&v)
  {
    meta::aux_tupel<N-1,X>::v_uma(a,v.a);
    return*this;
  }

  template<int N, typename X> inline
  tupel<N,X>& tupel<N,X>::update_min (tupel<N,X> const&v)
  {
    meta::aux_tupel<N-1,X>::v_umi(a,v.a);
    return*this;
  }
#ifdef TUPEL_CROSS_TYPE
  // A.6.4.2 with tupel<N,scalar_type>
  template<int N, typename X> template<typename scalar_type> inline
  tupel<N,X>& tupel<N,X>::operator=  (tupel<N,scalar_type> const&v)
  {
    meta::aux_tupel<N-1,X>::v_as(a,v.a);
    return*this;
  }

  template<int N, typename X> template<typename scalar_type> inline
  tupel<N,X>& tupel<N,X>::operator+= (tupel<N,scalar_type> const&v)
  {
    meta::aux_tupel<N-1,X>::v_ad(a,v.a);
    return*this;
  }

  template<int N, typename X> template<typename scalar_type> inline
  tupel<N,X>& tupel<N,X>::operator-= (tupel<N,scalar_type> const&v)
  {
    meta::aux_tupel<N-1,X>::v_su(a,v.a);
    return*this;
  }

  template<int N, typename X> template<typename scalar_type> inline
  tupel<N,X>  tupel<N,X>::operator+  (tupel<N,scalar_type> const&v) const
  {
    tupel<N,X> y; 
    meta::aux_tupel<N-1,X>::v_sum(y.a,a,v.a);
    return y;
  }

  template<int N, typename X> template<typename scalar_type> inline
  tupel<N,X>  tupel<N,X>::operator-  (tupel<N,scalar_type> const&v) const
  {
    tupel<N,X> y; 
    meta::aux_tupel<N-1,X>::v_dif(y.a,a,v.a);
    return y;
  }

  //     add or subtract another tupel integer times
  template<int N, typename X> template<int K, typename scalar_type> inline
  tupel<N,X>& tupel<N,X>::add_int_times(tupel<N,scalar_type> const&v)
  {
    meta::aux_tupel<N-1,X>:: template v_adit<K>(a,v.a);
    return*this;
  }

  template<int N, typename X> template<int K, typename scalar_type> inline
  tupel<N,X>& tupel<N,X>::sub_int_times(tupel<N,scalar_type> const&v)
  {
    meta::aux_tupel<N-1,X>:: template v_suit<K>(a,v.a);
    return*this;
  }
#endif // TUPEL_CROSS_TYPE

  // A.6.5   miscellaneous
  // A.6.5.1 copy from array (for security, operator= is not provided)
  template<int N, typename X> inline
  tupel<N,X>& tupel<N,X>::copy(const element_type*x)
  {
    meta::aux_tupel<N-1,X>::v_as(a,x);
    return*this;
  }
#ifdef TUPEL_CROSS_TYPE
  template<int N, typename X> template<typename scalar_type> inline
  tupel<N,X>& tupel<N,X>::copy(const scalar_type*x)
  {
    meta::aux_tupel<N-1,X>::v_as(a,x);
    return*this;
  }
#endif
  // A.6.5.2 assign, add to, or subtract from: real times another tupel
  template<int N, typename X> template<typename scalar_type> inline
  tupel<N,X>& tupel<N,X>::ass_times(tupel<N,X> const&v, scalar_type const&x)
  {
    meta::aux_tupel<N-1,X>::v_ast(a,v.a,x);
    return*this;
  }

  template<int N, typename X> template<typename scalar_type> inline
  tupel<N,X>& tupel<N,X>::add_times(tupel<N,X> const&v, scalar_type const&x)
  {
    meta::aux_tupel<N-1,X>::v_adt(a,v.a,x);
    return*this;
  }

  template<int N, typename X> template<typename scalar_type> inline
  tupel<N,X>& tupel<N,X>::sub_times(tupel<N,X> const&v, scalar_type const&x)
  {
    meta::aux_tupel<N-1,X>::v_sut(a,v.a,x);
    return*this;
  }

  // A.6.6   non-member methods for tupel<N,X>
  // A.6.6.1 formatted I/O
  template<int N, typename X> inline
  std::ostream&operator<< (std::ostream&s, tupel<N,X> const&v)
  {
    meta::aux_tupel<N-1,X>::v_out(s,static_cast<const X*>(v));
    return s;
  }

  template<int N, typename X> inline
  std::istream&operator>> (std::istream&s, tupel<N,X>&v)
  {
    meta::aux_tupel<N-1,X>::v_in(s,static_cast<X*>(v));
    return s;
  }

  // A.6.6.2 norm and absolute value: norm == sum_i x[i]*x[i]
  template<int N, typename X> inline
  X norm(tupel<N,X> const&v)
  {
    return meta::aux_tupel<N-1,X>::v_norm(static_cast<const X*>(v));
  }

  template<int N, typename X> inline
  X abs (tupel<N,X> const&v)
  {
    return sqrt(norm(v));
  }

  // A.6.6.3 max(abs(element))
  template<int N, typename X> inline
  X maxnorm(tupel<N,X> const&v)
  {
    return meta::aux_tupel<N-1,X>::v_amax(static_cast<const X*>(v));
  }

  // A.6.6.4 distance, squared distance, and squared sum 
  template<int N, typename X> inline
  X dist_sq (tupel<N,X> const&v, tupel<N,X> const&w)
  {
    return meta::aux_tupel<N-1,X>::v_diq(static_cast<const X*>(v),
					 static_cast<const X*>(w));
  }

  template<int N, typename X> inline
  X dist    (tupel<N,X> const&v, tupel<N,X> const&w)
  {
    return sqrt(dist_sq(v,w));
  }

  template<int N, typename X> inline
  X sum_sq  (tupel<N,X> const&v, tupel<N,X> const&w)
  {
    return meta::aux_tupel<N-1,X>::v_suq(static_cast<const X*>(v),
					 static_cast<const X*>(w));
  }
#ifdef TUPEL_CROSS_TYPE
  template<int N, typename X, typename Y> inline
  X dist_sq (tupel<N,X> const&v, tupel<N,Y> const&w)
  {
    return meta::aux_tupel<N-1,X>::v_diq(static_cast<const X*>(v),
					 static_cast<const Y*>(w));
  }

  template<int N, typename X, typename Y> inline
  X dist    (tupel<N,X> const&v, tupel<N,Y> const&w)
  {
    return sqrt(dist_sq(v,w));
  }

  template<int N, typename X, typename Y> inline
  X sum_sq  (tupel<N,X> const&v, tupel<N,Y> const&w)
  {
    return meta::aux_tupel<N-1,X>::v_suq(static_cast<const X*>(v),
					 static_cast<const Y*>(w));
  }
#endif
  // A.6.6.5 vector dot product
  template<int N, typename X> inline
  X operator* (tupel<N,X> const&v, tupel<N,X> const&w)
  {
    return meta::aux_tupel<N-1,X>::v_dot(static_cast<const X*>(v),
					 static_cast<const X*>(w));
  }
#ifdef TUPEL_CROSS_TYPE
  template<int N, typename X, typename Y> inline
  X operator* (tupel<N,X> const&v, tupel<N,Y> const&w)
  {
    return meta::aux_tupel<N-1,X>::v_dot(static_cast<const X*>(v),
					 static_cast<const Y*>(w));
  }
#endif

  // A.6.6.6 vector cross product in 2D and 3D, coded as operator ^
  template<typename X> inline
  X operator^ (tupel<2,X> const&v, tupel<2,X> const&w)
  {
    return v[0]*w[1] - v[1]*w[0];
  }
  template<typename X> inline
  tupel<3,X> operator^ (tupel<3,X> const&v, tupel<3,X> const&w)
  {
    return tupel<3,X>(v[1]*w[2] - v[2]*w[1],
		      v[2]*w[0] - v[0]*w[2],
		      v[0]*w[1] - v[1]*w[0]);
  }
#ifdef TUPEL_CROSS_TYPE
  template<typename X, typename Y> inline
  X operator^ (tupel<2,X> const&v, tupel<2,Y> const&w)
  {
    return v[0]*w[1] - v[1]*w[0];
  }
  template<typename X, typename Y> inline
  tupel<3,X> operator^ (tupel<3,X> const&v, tupel<3,Y> const&w)
  {
    return tupel<3,X>(v[1]*w[2] - v[2]*w[1],
		      v[2]*w[0] - v[0]*w[2],
		      v[0]*w[1] - v[1]*w[0]);
  }
#endif

} // namespace WD
#endif // included_tupel_h
