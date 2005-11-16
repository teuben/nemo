// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    inc/public/tupel.h                                                 
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    1996-2005                                                          
///                                                                             
/// \brief   contains the definition of template class falcON::tupel and        
///	     all its members.                                                   
///                                                                             
/// \version aug-2003: template metaprogramming to unroll loops automatically   
/// \version nov-2003: non-standard #include files redundant; pseudo_tupel      
/// \version nov-2004: added volume(), applied();                               
/// \version may-2005: removed tupel::add_times, ass_times, sub_times           
/// \version jun-2005: removed pseudo_tupel, const_pseudo_tupel                 
/// \version jul-2005: added doxygen documentation, added minnorm()             
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 1996-2005  Walter Dehnen                                       
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
#ifndef falcON_included_tupel_h
#define falcON_included_tupel_h

#ifndef falcON_included_iostream
#  include <iostream>
#  define falcON_included_iostream
#endif
#ifndef falcON_included_cmath
#  include <cmath>
#  define falcON_included_cmath
#endif
#ifndef falcON_included_tupel_cc
#  include <public/tupel.cc>
#endif

////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::tupel                                                        
  //                                                                            
  /// \brief                                                                    
  /// Template: a tupel of N scalars of type X, held in an array X[N].          
  /// Type tupel<3,real> is used in falcON to represent vectors (position,      
  /// velocity etc). All operations are inlined using metatemplate programming. 
  // ///////////////////////////////////////////////////////////////////////////
  template<int N, typename X> class tupel {
    // private types and data                                                   
  private:
    typedef meta::taux<X,N-1,0> M;
    X a[N];                          ///< data: an array of N elements of type X
  public:
    //--------------------------------------------------------------------------
    // static members and public types                                          
    typedef X        element_type;               ///< type of elements          
    typedef tupel    V;                          ///< type tupel                
    static const int NDAT = N;                   ///< number of elements        
    static const int ORD  = 1;                   ///< tensor rank or order      
    static       int size() { return N; }        ///< number of elements        
    //--------------------------------------------------------------------------
    // construction                                                             
    /// unitialized construction
    tupel() {}
    /// copy from another tupel
    tupel(tupel const&x) {
      M::v_as(a,x.a);
    }
    /// all elements equal a scalar
    explicit tupel(X const&x) {
      M::s_as(a,x);
    }
    /// set elements equal to elements of an array
    explicit tupel(const X* x) {
      M::v_as(a,x);
    }
    /// copy from a tupel with another element type
    template<typename S> tupel (tupel<N,S> const&x) {
      M::v_as(a,static_cast<const S*>(x));
    }
    /// all elements equal a scalar of type different from element type
    template<typename S> explicit tupel(const S*x) {
      M::v_as(a,x);
    }
    //--------------------------------------------------------------------------
    /// \name construction from list of elements (up to N=10)                   
    //{@
    /// from 2 elements (for N=2)
    tupel      (X const&x0, X const&x1) {
      a[0]=x0;
      a[1]=x1;
    }
    /// from 3 elements (for N=3)
    tupel      (X const&x0, X const&x1, X const&x2) {
      a[0]=x0;
      a[1]=x1;
      a[2]=x2;
    }
    /// from 4 elements (for N=4)
    tupel      (X const&x0, X const&x1, X const&x2, X const&x3) {
      a[0]=x0;
      a[1]=x1;
      a[2]=x2;
      a[3]=x3;
    }
    /// from 5 elements (for N=5)
    tupel      (X const&x0, X const&x1, X const&x2, X const&x3, X const&x4) {
      a[0]=x0;
      a[1]=x1;
      a[2]=x2;
      a[3]=x3;
      a[4]=x4;
    }
    /// from 6 elements (for N=6)
    tupel      (X const&x0, X const&x1, X const&x2, X const&x3, X const&x4,
		X const&x5) { 
      a[0]=x0;
      a[1]=x1;
      a[2]=x2;
      a[3]=x3;
      a[4]=x4;
      a[5]=x5;
    }
    /// from 7 elements (for N=7)
    tupel      (X const&x0, X const&x1, X const&x2, X const&x3, X const&x4,
		X const&x5, X const&x6) {
      a[0]=x0;
      a[1]=x1;
      a[2]=x2;
      a[3]=x3;
      a[4]=x4;
      a[5]=x5;
      a[6]=x6;
    }
    /// from 8 elements (for N=8)
    tupel      (X const&x0, X const&x1, X const&x2, X const&x3, X const&x4,
		X const&x5, X const&x6, X const&x7) {
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
    tupel      (X const&x0, X const&x1, X const&x2, X const&x3, X const&x4,
		X const&x5, X const&x6, X const&x7, X const&x8) {
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
    tupel      (X const&x0, X const&x1, X const&x2, X const&x3, X const&x4,
		X const&x5, X const&x6, X const&x7, X const&x8, X const&x9) {
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
    //--------------------------------------------------------------------------
    /// \name element access                                                    
    //@{
    /// non-const element access via sub-script operator
    X      &operator[] (int i)       { return a[i]; }
    /// const element access via sub-script operator
    X const&operator[] (int i) const { return a[i]; }
#ifdef falcON_TUPEL_FUNCOP
    X      &operator() (int i)       { return a[i]; }
    X const&operator() (int i) const { return a[i]; }
#endif
    //@}
    //--------------------------------------------------------------------------
    /// \name unitary operators and methods                                     
    //@{
    /// change sign for each element
    tupel&negate() {
      M::v_neg(a);
      return*this;
    }
    /// return negative of this
    tupel operator- () const {
      register tupel y;
      M::v_nega(y.a,a);
      return y;
    }
    /// return norm := Sum element^2
    X norm() const {
      return M::v_norm(a);
    }
    /// return |x| := sqrt(norm())
    X abs() const {
      return std::sqrt(norm());
    }
    /// return minimum element
    X min() const {
      return M::v_min(a);
    }
    /// return maximum element
    X max() const {
      return M::v_max(a);
    }
    /// return maximum |element|
    X maxnorm() const {
      return M::v_amax(a);
    }
    /// return minimum |element|
    X minnorm() const {
      return M::v_amin(a);
    }
    /// return Product of elements
    X volume() const {
      return M::v_vol(a);
    }
    /// conversion to pointer to scalar
    operator       X*  ()       { return a; }
    /// conversion to const pointer to scalar
    operator const X*  () const { return a; }
    //@}
    //--------------------------------------------------------------------------
    /// \name binary operators with scalar                                      
    //@{
    /// assign all elements to scalar
    template<typename S> tupel&operator= (S const&x) {
      M::s_as(a,x);
      return*this;
    }
    /// multiply all elements by scalar
    template<typename S> tupel&operator*= (S const&x) {
      M::s_ml(a,x);
      return*this;
    }
    /// divide all elements by scalar
    template<typename S> tupel&operator/= (S const&x) {
      return operator*=(S(1)/x);
    }
    /// return product with scalar
    template<typename S> tupel operator* (S const&x) const {
      register tupel y;
      M::v_ast(y.a,a,x);
      return y;
    }
    /// return product with inverse of scalar
    template<typename S> tupel operator/ (S const&x) const {
      return operator*(S(1)/x);
    }
    /// are all elements equal to scalar?
    bool operator== (X const&x) const {
      return M::s_eq(a,x);
    }
    /// is any element un-equal to scalar?
    bool operator!= (X const&x) const {
      return M::s_neq(a,x);
    }
    //@}
    //--------------------------------------------------------------------------
    /// \name binary operators with tupel<N,X>                                  
    //@{
    /// set *this equal to x: set tupel::a[i] = x[i]
    tupel&operator= (tupel const&x) {
      M::v_as(a,x.a);
      return*this;
    }
    /// add x to *this: set tupel::a[i] += x[i]
    tupel&operator+= (tupel const&x) {
      M::v_ad(a,x.a);
      return*this;
    }
    /// subtract x from *this: set tupel::a[i] -= x[i]
    tupel&operator-= (tupel const&x) {
      M::v_su(a,x.a);
      return*this;
    }
    /// element-wise *= : set tupel::a[i] *= x[i]
    tupel&ass_mul(tupel const&x) {
      M::v_ml(a,x.a);
      return*this;
    }
    /// element-wise /= : set tupel::a[i] /= x[i]
    tupel&ass_div(tupel const&x) {
      M::v_dv(a,x.a);
      return*this;
    }
    /// is x the same as *this?
    bool operator==(tupel const&x) const {
      return M::v_eq(a,x.a);
    }
    /// is x not the same as *this?
    bool operator!=(tupel const&x) const {
      return M::v_neq(a,x.a);
    }
    /// sum: return *this + x
    tupel operator+ (tupel const&x) const {
      register tupel y;
      M::v_sum(y.a,a,x.a);
      return y;
    }
    /// difference: return *this - x
    tupel operator- (tupel const&x) const {
      register tupel y;
      M::v_dif(y.a,a,x.a);
      return y;
    }
    /// scalar (dot) product: return Sum tupel::a[i]*x[i]
    X operator* (tupel const&x) const {
      return M::v_dot(a,x.a);
    }
    /// difference squared: return (*this-x)^2 := Sum (tupel::a[i]-x[i])^2
    X dist_sq(tupel const&x) const {
      return M::v_diq(a,x.a);
    }
    /// distance: return |*this-x| := sqrt(Sum (tupel::a[i]-x[i])^2)
    X dist(tupel const&x) const {
      return std::sqrt(dist_sq (x));
    }
    /// sum squared: return (this+x)^2 := Sum (tupel::a[i]+x[i])^2
    X sum_sq (tupel const&x) const {
      return M::v_suq(a,x.a);
    }
    /// update maximum element-wise: tupel::a[i] = max(tupel::a[i], x[i])
    tupel&up_max (tupel const&x) {
      M::v_uma(a,x.a);
      return*this;
    }
    /// update minimum element-wise: tupel::a[i] = min(tupel::a[i], x[i])
    tupel&up_min (tupel const&x) {
      M::v_umi(a,x.a);
      return*this;
    }
    /// update minimum and maximum element-wise:
    /// min[i] = min(min[i], tupel::a[i]); and
    /// max[i] = max(max[i], tupel::a[i]);
    void up_min_max(tupel &min, tupel&max) const {
      M::v_umia(min.a,max.a,a);
    }
    //@}
    //--------------------------------------------------------------------------
    /// \name binary operators with tupel<N,S>                                  
    //@{
    /// set *this equal to x: set tupel::a[i] = x[i]
    template<typename S> tupel&operator= (tupel<N,S> const&x) {
      M::v_as(a,static_cast<const S*>(x));
      return*this;
    }
    /// add x to *this: set tupel::a[i] += x[i]
    template<typename S> tupel&operator+= (tupel<N,S> const&x) {
      M::v_ad(a,static_cast<const S*>(x));
      return*this;
    }
    /// subtract x from *this: set tupel::a[i] -= x[i]
    template<typename S> tupel&operator-= (tupel<N,S> const&x) {
      M::v_su(a,static_cast<const S*>(x));
      return*this;
    }
    /// sum: return *this + x
    template<typename S> tupel operator+ (tupel<N,S> const&x) const {
      register tupel y;
      M::v_sum(y.a,a,static_cast<const S*>(x));
      return y;
    }
    /// difference: return *this - x
    template<typename S> tupel operator- (tupel<N,S> const&x) const {
      register tupel y;
      M::v_dif(y.a,a,static_cast<const S*>(x));
      return y;
    }
    /// scalar (dot) product: return Sum tupel::a[i]*x[i]
    template<typename S> X operator* (tupel<N,S> const&x) const {
      return M::v_dot(a,(const S*)x);
    }
    /// difference squared: return (*this-x)^2 := Sum (tupel::a[i]-x[i])^2
    template<typename S> X dist_sq(tupel<N,S> const&x) const {
      return M::v_diq(a,(const S*)x);
    }
    /// distance: return |*this-x| := sqrt(Sum (tupel::a[i]-x[i])^2)
    template<typename S> X dist (tupel<N,S> const&x) const {
      return std::sqrt(dist_sq (x));
    }
    /// sum squared: return (this+x)^2 := Sum (tupel::a[i]+x[i])^2
    template<typename S> X sum_sq (tupel<N,S> const&x) const {
      return M::v_suq(a,(const S*)x);
    }
    /// update minimum and maximum element-wise:
    /// min[i] = min(min[i], tupel::a[i]); and
    /// max[i] = max(max[i], tupel::a[i]);
    template<typename S> void up_min_max(tupel<N,S>&min, tupel<N,S>&max) const {
      M::v_umia((S*)min,(S*)max,a);
    }
    //@}
    //--------------------------------------------------------------------------
    /// \name miscellaneous                                                     
    //@{
    /// set elements equal to array elements
    tupel&copy(const X*x) {
      M::v_as(a,x);
      return*this;
    }
    /// set elements equal to array elements
    template<typename S> tupel&copy(const S*x) {
      M::v_as(a,x);
      return*this;
    }
    /// update maximum element-wise: tupel::a[i] = max(tupel::a[i], x_i + f)
    tupel&up_max(tupel const&x, X const&f) {
      M::v_umax(a,x.a,f);
      return*this;
    }
    /// update minimum element-wise: tupel::a[i] = min(tupel::a[i], x_i - f)
    tupel&up_min(tupel const&x, X const&f) {
      M::v_umix(a,x.a,f);
      return*this;
    }
    /// apply element-wise function call: tupel::a[i] = f(tupel::a[i])
    tupel&apply(X(*f)(X)) {
      M::v_appl(a,a,f);
      return*this;
    }
    /// return element-wise function call: return f(tupel::a[i])
    tupel applied(X(*f)(X)) const {
      register tupel y;
      M::v_appl(y.a,a,f);
      return y;
    }
    /// replace by element-wise function: tupel::a[i] = f(x[i])
    tupel&connect(tupel const&x, X(*f)(X)) {
      M::v_appl(a,x.a,f);
      return*this;
    }
    //@}
    //--------------------------------------------------------------------------
    /// \name formatted I/O                                                     
    //@{
    /// formatted output: space separated
    friend std::ostream&operator<<(std::ostream&s ,tupel const&x) {
      M::v_out(s,x.a);
      return s;
    }
    /// formatted input: read element-wise
    friend std::istream&operator>>(std::istream&s, tupel &x) {
      M::v_in (s,x.a);
      return s;
    }
    //@}
    //--------------------------------------------------------------------------
    /// \name functions taking tupel<> arguments                                
    //@{                                                                        
    /// is y[i] == x for all i?
    friend bool operator==(X const&x,tupel const&y) {
      return y==x;
    }
    /// is y[i] != x for any i?
    friend bool operator!=(X const&x,tupel const&y) {
      return y!=x;
    }
    /// return maximum element of tupel
    friend X max(tupel const&x) {
      return x.max();
    }
    /// return maximum |element| of tupel
    friend X maxnorm(tupel const&x) {
      return x.maxnorm();
    }
    /// return minimum element of tupel
    friend X min(tupel const&x) {
      return x.min();
    }
    /// return minimum |element| of tupel
    friend X minnorm(tupel const&x) {
      return x.minnorm();
    }
    /// return norm of tupel: Sum x[i]^2
    friend X norm(tupel const&x) {
      return x.norm();
    }
    /// return absolute value: sqrt(Sum x[i]^2)
    friend X abs(tupel const&x) {
      return x.abs();
    }
    /// return Product of elements
    friend X volume (tupel const&x) {
      return x.volume();
    }
    /// return difference squared: (x-y)^2 := Sum (x[i]-y[i])^2
    template<typename S> friend X dist_sq(tupel const&x, tupel<N,S> const&y) {
      return x.dist_sq(y);
    }
    /// return distance: (x-y)^2 := Sum (x[i]-y[i])^2
    template<typename S> friend X dist(tupel const&x, tupel<N,S> const&y) {
      return x.dist(y);
    }
    /// return sum squared: (x-y)^2 := Sum (x[i]+y[i])^2
    template<typename S> friend X sum_sq(tupel const&x, tupel<N,S> const&y) {
      return x.sum_sq(y);
    }
    /// product with scalar: x[i] = y * v[i]
    template<typename S> friend tupel operator*(S const&y, tupel const&v) {
      return v*y;
    }
    /// is any element nan?
    friend bool isnan(tupel const&x) {
      return meta::taux<X,N-1>::v_nan(x.a);
    }
    /// is any element inf?
    friend bool isinf(tupel const&x) {
      return meta::taux<X,N-1>::v_inf(x.a);
    }
    /// update maximum element-wise: x[i] = max(x[i], y[i])
    friend void update_max(tupel&x, tupel const&y) {
      return x.up_max(y);
    }
    /// update minimum element-wise: x[i] = min(x[i], y[i])
    friend void update_min(tupel&x, tupel const&y) {
      return x.up_min(y);
    }
    //@}
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  /// \relates falcON::tupel
  /// \name vector cross product in 2D and 3D
  //@{
  /// vector cross product for N=2: returns x[0]*y[1]- x[1]*y[0]
  template<typename X>
  inline X operator^ (tupel<2,X> const&x, tupel<2,X> const&y) {
    return x[0]*y[1] - x[1]*y[0];
  }
  /// vector cross product for N=3
  template<typename X>
  inline tupel<3,X> operator^ (tupel<3,X> const&x, tupel<3,X> const&y) {
    return tupel<3,X>(x[1]*y[2] - x[2]*y[1],
		      x[2]*y[0] - x[0]*y[2],
		      x[0]*y[1] - x[1]*y[0]);
  }
  //@}
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class traits< tupel<N,T> >                                               //
  //                                                                          //
  // ///////////////////////////////////////////////////////////////////////////
#ifdef falcON_included_traits_h
  template<int N, typename T> struct traits< tupel<N,T> > {
    static const char  *name() {
      return message("tupel<%d,%s>", N, traits<T>::name());
    }
    static const char  *names() {
      return message("tupel<%d,%s>s", N, traits<T>::name());
    }
    static const size_t size = sizeof(tupel<N,T>);
  };
#endif
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#endif// falcON_included_tupel_h
