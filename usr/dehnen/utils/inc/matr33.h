// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/matr33.h
///
/// \author Walter Dehnen
///                                                                             
/// \date   2010-2012
///
////////////////////////////////////////////////////////////////////////////////
///
/// \version June-2012 WD  using vector_class to allow vector.h and tupel.h
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2010-2012 Walter Dehnen
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
#ifndef WDutils_included_matr33_h
#define WDutils_included_matr33_h

#ifndef WDutils_included_exception_h
#  include <exception.h>
#endif
#ifndef WDutils_included_inline_h
#  include <inline.h>
#endif

#if __cpluspluc < 201103L
#  define noexcept
#endif
namespace WDutils {

  template<typename X> class SymmMatrix33;

  /// a general 3x3 matrix
  /// \note I just started this, will be adding functionality as needed.
  template<typename X>
  class Matrix33 {
    static const int N=9;
    union {
      X M[3][3];
      X A[N];
    };
    friend class SymmMatrix33<X>;
    Matrix33(X a0, X a1, X a2, X a3, X a4, X a5, X a6, X a7, X a8) noexcept
    { A[0]=a0; A[1]=a1; A[2]=a2;
      A[3]=a3; A[4]=a4; A[5]=a5;
      A[6]=a6; A[7]=a7; A[8]=a8; }
  public:
    typedef SymmMatrix33<X> SymmMatrix;  ///< associated symmetric matrix type
    /// default ctor
    Matrix33() noexcept {}
    /// copy ctor
    Matrix33(Matrix33 const&m)  noexcept
    { for(int i=0; i!=N; ++i) A[i] =m.A[i]; }
    /// ctor from symmetric matrix
    inline Matrix33(SymmMatrix const&m) noexcept;
    /// ctor from scalar: set all elements to scalar
    explicit Matrix33(X x) noexcept
    { for(int i=0; i!=N; ++i) A[i] =x; }
    /// \name operations with scalar
    //@{
    /// product with scalar
    Matrix33&operator*=(X x) noexcept
    { for(int i=0; i!=N; ++i) A[i]*=x; return*this; }
    /// divide by scalar
    Matrix33&operator/=(X x) noexcept
    { return operator*=(1/x); }
    //}@
    /// \name inter matrix operations
    //@{
    /// assign to matrix
    Matrix33&operator=(Matrix33 const&m) noexcept
    { for(int i=0; i!=N; ++i) A[i] =m.A[i]; return*this; }
    /// assign to symmetric matrix
    inline Matrix33&operator =(SymmMatrix const&m) noexcept;
    /// add matrix
    Matrix33&operator+=(Matrix33 const&m) noexcept
    { for(int i=0; i!=N; ++i) A[i]+=m.A[i]; return*this; }
    /// add symmetric matrix
    inline Matrix33&operator+=(SymmMatrix const&m) noexcept;
    /// sum of two matrices
    Matrix33 operator+(Matrix33 const&m) const noexcept
    { Matrix33 r; for(int i=0; i!=N; ++i) r.A[i]=A[i]+m.A[i]; return r; }
    /// sum of matrix and symmetric matrix
    inline Matrix33 operator+ (SymmMatrix const&m) const;
    /// subtract matrix
    Matrix33&operator-=(Matrix33 const&m) noexcept
    { for(int i=0; i!=N; ++i) A[i]-=m.A[i]; return*this; }
    /// subtract symmetric matrix
    inline Matrix33&operator-=(SymmMatrix const&m) noexcept;
    /// difference between two matrices
    Matrix33 operator-(Matrix33 const&m) const noexcept
    { Matrix33 r; for(int i=0; i!=N; ++i) r.A[i]=A[i]-m.A[i]; return r; }
    /// difference between matrix and symmetric matrix
    inline Matrix33 operator- (SymmMatrix const&m) const noexcept;
    /// square of matrix
    Matrix33 square() const noexcept
    {
      Matrix33 m;
      X t01  = M[0][1]*M[1][0];
      X t02  = M[0][2]*M[2][0];
      X t12  = M[1][2]*M[2][1];
      m.A[0] = M[0][0]*M[0][0]+      t01      +      t02;
      m.A[1] = M[0][0]*M[0][1]+M[0][1]*M[1][1]+M[0][2]*M[2][1];
      m.A[1] = M[0][0]*M[0][2]+M[0][1]*M[1][2]+M[0][2]*M[2][2];
      m.A[1] = M[1][0]*M[0][0]+M[1][1]*M[1][0]+M[1][2]*M[2][0];
      m.A[1] =       t01      +M[1][1]*M[1][1]+      t12;
      m.A[1] = M[1][0]*M[0][2]+M[1][1]*M[1][2]+M[1][2]*M[2][2];
      m.A[1] = M[2][0]*M[0][0]+M[2][1]*M[1][0]+M[2][2]*M[2][0];
      m.A[1] = M[2][0]*M[0][1]+M[2][1]*M[1][1]+M[2][2]*M[2][1];
      m.A[1] =       t02      +      t12      +M[2][2]*M[2][2];
      return m;
    }
    /// product of two matrices
    Matrix33 operator*(Matrix33 const&m) const noexcept
    {
      return Matrix33(M[0][0]*m[0][0]+M[0][1]*m[1][0]+M[0][2]*m[2][0],
		      M[0][0]*m[0][1]+M[0][1]*m[1][1]+M[0][2]*m[2][1],
		      M[0][0]*m[0][2]+M[0][1]*m[1][2]+M[0][2]*m[2][2],
		      M[1][0]*m[0][0]+M[1][1]*m[1][0]+M[1][2]*m[2][0],
		      M[1][0]*m[0][1]+M[1][1]*m[1][1]+M[1][2]*m[2][1],
		      M[1][0]*m[0][2]+M[1][1]*m[1][2]+M[1][2]*m[2][2],
		      M[2][0]*m[0][0]+M[2][1]*m[1][0]+M[2][2]*m[2][0],
		      M[2][0]*m[0][1]+M[2][1]*m[1][1]+M[2][2]*m[2][1],
		      M[2][0]*m[0][2]+M[2][1]*m[1][2]+M[2][2]*m[2][2]);
    }
    /// product of matrix with symmetric matrix
    inline Matrix33 operator*(SymmMatrix const&m) const noexcept;
    /// symmetric part of matrix
    inline SymmMatrix symm_part() const noexcept;
    /// inverse of matrix
    Matrix33 inverse() const noexcept
    {
      Matrix33 m;
      m[0][0] = M[1][1]*M[2][2]-M[1][2]*M[2][1];
      m[0][1] = M[0][2]*M[2][1]-M[0][1]*M[2][2];
      m[0][2] = M[0][1]*M[1][2]-M[0][2]*M[1][1];
      m[1][0] = M[1][2]*M[2][0]-M[1][0]*M[2][2];
      m[1][1] = M[0][0]*M[2][2]-M[0][2]*M[2][0];
      m[1][2] = M[0][2]*M[1][0]-M[0][0]*M[1][2];
      m[2][0] = M[1][0]*M[2][1]-M[1][1]*M[2][0];
      m[2][1] = M[0][1]*M[2][0]-M[0][0]*M[2][1];
      m[2][2] = M[0][0]*M[1][1]-M[0][1]*M[1][0];
      return m/= M[0][0]*m[0][0]-M[1][0]*m[0][1]+M[2][0]*m[0][2];
    }
    /// remove trace-part of matrix
    Matrix33&traceless() noexcept
    {
      X tmp = X(0.333333333333333333) * trace();
      A[0] -= tmp;
      A[4] -= tmp;
      A[8] -= tmp;
      return*this;
    }
    //@}
    /// \name predicates
    //@{
    /// determinant of matrix
    X det() const noexcept
    { 
      return 
	M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1]) -
	M[1][0]*(M[0][1]*M[2][2]-M[2][1]*M[0][2]) +
	M[2][0]*(M[0][1]*M[1][2]-M[1][1]*M[0][2]);
    }
    /// trace = Mkk
    /// \note norm = 1/3 trace^2 + shear^2 + vorticity^2
    X trace() const noexcept
    { return M[0][0]+M[1][1]+M[2][2]; }
    /// trace of M^2
    X trace_of_square() const noexcept
    {
      return   M[0][0]*M[0][0] + M[1][1]*M[1][1] + M[2][2]*M[2][2]+
	twice (M[0][1]*M[1][0] + M[0][2]*M[2][0] + M[1][2]*M[2][1]);
    }
    /// trace of A*B
    /// \note trace(A*B) = trace(B*A)
    X trace_of_prod(Matrix33 const&m) const noexcept
    {
      return
	M[0][0]*m[0][0]+M[0][1]*m[1][0]+M[0][2]*m[2][0] +
	M[1][0]*m[0][1]+M[1][1]*m[1][1]+M[1][2]*m[2][1] +
	M[2][0]*m[0][2]+M[2][1]*m[1][2]+M[2][2]*m[2][2] ;
    }
    /// trace of A*B
    X trace_of_prod(SymmMatrix const&m) const noexcept;
    /// norm = trace(M*M^t) = Mij*Mij
    /// \note norm = 1/3 trace^2 + shear^2 + vorticity^2
    X norm() const noexcept
    { X n(0); for(int i=0; i!=N; ++i) n+=A[i]*A[i]; return n; }
    /// shear: sqrt(norm(traceless symmetric part of matrix))
    /// \note norm = 1/3 trace^2 + shear^2 + vorticity^2
    X shear() const noexcept
    {
      X tmp = X(0.3333333333333333333) * trace();
      return std::sqrt(WDutils::square(A[0]-tmp) +
		       WDutils::square(A[4]-tmp) +
		       WDutils::square(A[8]-tmp) +
		       X(0.5)*(WDutils::square(A[1]+A[3])+
			       WDutils::square(A[2]+A[6])+
			       WDutils::square(A[5]+A[7])));
    }
    /// vorticity-squared = norm(anti-symmetric part of matrix)
    X vorticity_sq() const noexcept
    {
      return X(0.5)*(WDutils::square(A[1]-A[3])+
		     WDutils::square(A[2]-A[6])+
		     WDutils::square(A[5]-A[7]));
    }
    /// vorticity: sqrt(norm(anti-symmetric part of matrix))
    /// \note norm = 1/3 trace^2 + shear^2 + vorticity^2
    X vorticity() const noexcept
    { return std::sqrt(vorticity_sq()); }
    //@}
    /// \name data access
    //@{
    /// const element access
    const X* operator[](int i) const noexcept
    { return M[i]; }
    /// non-const element access
    X* operator[](int i) noexcept
    { return M[i]; }
    /// const access to internal representation
    X const&operator()(int i) const noexcept
    { return A[i]; }
    /// non-const access to internal representation
    X&operator()(int i) noexcept
    { return A[i]; }
    //@}
    /// \name operations with vector
    //@{
    /// set to dyadic vector product
    template<template<int,typename> class vector_class>
    Matrix33&dyadic(vector_class<3,X> const&x,
		    vector_class<3,X> const&y) noexcept
    {
      for(int i=0; i!=3; ++i)
	for(int j=0; j!=3; ++j)
	  M[i][j]=x[i]*y[j];
      return*this;
    }
    /// set to dyadic vector product
    template<template<int,typename> class vector_class>
    Matrix33&dyadic(vector_class<3,X> const&x) noexcept
    {
      A[0]     =x[0]*x[0];
      A[1]=A[3]=x[0]*x[1];
      A[2]=A[6]=x[0]*x[2];
      A[4]     =x[1]*x[1];
      A[5]=A[7]=x[1]*x[2];
      A[8]     =x[2]*x[2];
      return*this;
    }
    /// add dyadic vector product
    template<template<int,typename> class vector_class>
    Matrix33&add_dyadic(vector_class<3,X> const&x,
			vector_class<3,X> const&y) noexcept
    {
      for(int i=0; i!=3; ++i)
	for(int j=0; j!=3; ++j)
	  M[i][j]+=x[i]*y[j];
      return*this;
    }
    /// add dyadic vector product times a scalar
    template<template<int,typename> class vector_class>
    Matrix33&add_dyadic(X fac, vector_class<3,X> const&x, 
			       vector_class<3,X> const&y) noexcept
    {
      for(int i=0; i!=3; ++i) {
	X tmp = fac*x[i];
	for(int j=0; j!=3; ++j)
	  M[i][j]+=tmp*y[j];
      }
      return*this;
    }
    /// product with vector from right
    template<template<int,typename> class vector_class>
    vector_class<3,X> operator*(vector_class<3,X> const&x) const noexcept
    {
      return vector_class<3,X>(M[0][0]*x[0] + M[0][1]*x[1] + M[0][2]*x[2],
			       M[1][0]*x[0] + M[1][1]*x[1] + M[1][2]*x[2],
			       M[2][0]*x[0] + M[2][1]*x[1] + M[2][2]*x[2]);
    }
    //@}
  };// class Matrix33

  //
  /// a symmetric 3x3 matrix
  /// \note I just started this, will be adding functionality as needed.
  template<typename X>
  class SymmMatrix33 {
    static const int N = 6;
    X A[N];
    friend class Matrix33<X>;
    SymmMatrix33(X a0, X a1, X a2, X a3, X a4, X a5)
    { A[0]=a0; A[1]=a1; A[2]=a2; A[3]=a3; A[4]=a4; A[5]=a5; }
  public:
    typedef Matrix33<X> Matrix;     ///< associated general matrix type
    /// default ctor
    SymmMatrix33() noexcept {}
    /// copy ctor
    SymmMatrix33(SymmMatrix33 const&M)  noexcept
    { for(int i=0; i!=N; ++i) A[i]=M.A[i]; }
    /// ctor from scalar: set all elements to scalar
    explicit SymmMatrix33(X x) noexcept
    { for(int i=0; i!=N; ++i) A[i] =x; }
    /// \name operations with scalar
    //@{
    /// assign all elements to scalar
    SymmMatrix33&operator=(X x) noexcept
    { for(int i=0; i!=N; ++i) A[i]=x; return*this; }
    /// multiply with scalar
    SymmMatrix33&operator*=(X x) noexcept
    { for(int i=0; i!=N; ++i) A[i]*=x; return*this; }
    /// divide by scalar
    SymmMatrix33&operator/=(X x) noexcept
    { return operator*=(1/x); }
    //}@
    /// \name inter matrix operations
    //@{
    /// assign to symmetric matrix
    SymmMatrix33&operator=(SymmMatrix33 const&m) noexcept
    { for(int i=0; i!=N; ++i) A[i]=m.A[i]; return*this; }
    /// add symmetric matrix
    SymmMatrix33&operator+=(SymmMatrix33 const&m) noexcept
    { for(int i=0; i!=N; ++i) A[i]+=m.A[i]; return*this; }
    /// sum of two symmetric matrices
    SymmMatrix33 operator+(SymmMatrix33 const&m) const noexcept
    { SymmMatrix33 r; for(int i=0; i!=N; ++i) r[i]=A[i]+m.A[i]; return r; }
    /// subtract symmetric matrix
    SymmMatrix33&operator-=(SymmMatrix33 const&m) noexcept
    { for(int i=0; i!=N; ++i) A[i]-=m.A[i]; return*this; }
    /// difference between two symmetric matrices
    SymmMatrix33 operator-(SymmMatrix33 const&m) const noexcept
    { SymmMatrix33 r; for(int i=0; i!=N; ++i) r[i]=A[i]-m.A[i]; return r; }
    /// product with symmetric matrix
    /// \note the produce of two symmetric matrices is in general not symmetric
    Matrix operator*(SymmMatrix33 const&m) const noexcept
    {
      register X t11=A[1]*m.A[1], t22=A[2]*m.A[2], t44=A[3]*m.A[3];
      return Matrix(A[0]*m.A[0] + t11         + t22,
		    A[0]*m.A[1] + A[1]*m.A[3] + A[2]*m.A[4],
		    A[0]*m.A[2] + A[1]*m.A[4] + A[2]*m.A[5],
		    A[1]*m.A[0] + A[3]*m.A[1] + A[4]*m.A[2],
		    t11         + A[3]*m.A[3] + t44,
		    A[1]*m.A[2] + A[3]*m.A[4] + A[4]*m.A[5],
		    A[2]*m.A[0] + A[4]*m.A[1] + A[5]*m.A[2],
		    A[2]*m.A[1] + A[4]*m.A[3] + A[5]*m.A[4],
		    t22         + t44         + A[5]*m.A[5]);
    }
    /// product with general matrix
    Matrix operator*(Matrix const&m) const noexcept
    {
      return Matrix(A[0]*m[0][0] + A[1]*m[1][0] + A[2]*m[2][0],
		    A[0]*m[0][1] + A[1]*m[1][1] + A[2]*m[2][1],
		    A[0]*m[0][2] + A[1]*m[1][2] + A[2]*m[2][2],
		    A[1]*m[0][0] + A[3]*m[1][0] + A[4]*m[2][0],
		    A[1]*m[0][1] + A[3]*m[1][1] + A[4]*m[2][1],
		    A[1]*m[0][2] + A[3]*m[1][2] + A[4]*m[2][2],
		    A[2]*m[0][0] + A[4]*m[1][0] + A[5]*m[2][0],
		    A[2]*m[0][1] + A[4]*m[1][1] + A[5]*m[2][1],
		    A[2]*m[0][2] + A[4]*m[1][2] + A[5]*m[2][2]);
    }
    /// unit matrix
    SymmMatrix33&unity() noexcept
    { A[0]=A[3]=A[5]=X(1); A[1]=A[2]=A[4]=X(0); return*this; }
    /// remove trace-part of matrix
    SymmMatrix33&traceless() noexcept
    {
      X tmp = X(0.333333333333333333) * trace();
      A[0] -= tmp;
      A[3] -= tmp;
      A[5] -= tmp;
      return*this;
    }
    /// invert symmetrix matrix
    /// \param[out]  M inverse of @c *this
    /// \note the inverse of a symmetric matrix is itself symmetric
    void invert(SymmMatrix33&M) const noexcept
    {
      M(0) = A[3]*A[5] - A[4]*A[4];
      M(1) = A[2]*A[4] - A[1]*A[5];
      M(2) = A[1]*A[4] - A[2]*A[3];
      M(3) = A[0]*A[5] - A[2]*A[2];
      M(4) = A[1]*A[2] - A[0]*A[4];
      M(5) = A[0]*A[3] - A[1]*A[1];
      M   /= A[0]*M(0) - A[1]*M(1) + A[2]*M(2);
    }
    /// inverse of symmetrix matrix
    /// \note the inverse of a symmetric matrix is itself symmetric
    SymmMatrix33 inverse() const noexcept
    {
      SymmMatrix33 M;
      invert(M);
      return M;
    }
    //@}
    /// \name predicates
    //@{
    /// determinant of matrix
    X det() const noexcept
    { 
      return
	A[0]*(A[3]*A[5]-A[4]*A[4])
	+A[1]*(twice(A[2]*A[4])-A[1]*A[5])
	-A[3]*A[2]*A[2];
    }
    /// trace = Mkk
    /// \note norm = 1/3 trace^2 + shear^2
    X trace() const noexcept
    { return A[0]+A[3]+A[5]; }
    /// trace of M^2 (same as norm for symmetric matrix)
    X trace_of_square() const noexcept
    { return norm(); }
    /// trace of A*B
    X trace_of_prod(Matrix const&m) const noexcept
    {
      return
	A[0]* m[0][0]          +A[3]* m[1][1]          +A[5]* m[2][2]         +
	A[1]*(m[1][0]+m[0][1]) +A[2]*(m[2][0]+m[0][2]) +A[4]*(m[2][1]+m[1][2]);
    }
    /// trace of A*B
    X trace_of_prod(SymmMatrix33 const&m) const noexcept
    {
      return  A[0]*m.A[0] + A[3]*m.A[3] + A[5]*m.A[5] +
	twice(A[1]*m.A[1] + A[2]*m.A[2] + A[4]*m.A[4]);
    }
    /// norm = trace(M*M^t) = Mij*Mij
    /// \note norm = 1/3 trace^2 + shear^2
    X norm() const noexcept
    { 
      return  A[0]*A[0] +A[3]*A[3] +A[5]*A[5] +
	twice(A[1]*A[1] +A[2]*A[2] +A[4]*A[4]);
    }
    /// shear: sqrt(norm(traceless part of symmetric matrix))
    /// \note norm = 1/3 trace^2 + shear^2
    X shear() const noexcept
    {
      X tmp = X(0.3333333333333333333) * trace();
      return std::sqrt(square(A[0]-tmp) +
		       square(A[3]-tmp) +
		       square(A[5]-tmp) +
		       twice(A[1]*A[1]+A[2]*A[2]+A[4]*A[4]));
    }
    //@}
    /// \name data access
    //@{
    /// const element access
    X const&operator()(int i) const noexcept
    { return A[i]; }
    /// non-const element access
    X &operator()(int i) noexcept
    { return A[i]; }
    //@}
    /// \name operations with vector
    //@{
    /// set to dyadic vector self-product
    template<template<int,typename> class vector_class>
    SymmMatrix33&dyadic(vector_class<3,X> const&x) noexcept
    {
      A[0]=x[0]*x[0];
      A[1]=x[0]*x[1];
      A[2]=x[0]*x[2];
      A[3]=x[1]*x[1];
      A[4]=x[1]*x[2];
      A[5]=x[2]*x[2];
      return*this;
    }
    /// set to symmetrised dyadic vector product
    template<template<int,typename> class vector_class>
    SymmMatrix33&dyadic_symm(vector_class<3,X> const&x,
			     vector_class<3,X> const&y) noexcept
    {
      A[0]=x[0]*y[0];
      A[1]=X(0.5)*(x[0]*y[1]+x[1]*y[0]);
      A[2]=X(0.5)*(x[0]*y[2]+x[2]*y[0]);
      A[3]=x[1]*y[1];
      A[4]=X(0.5)*(x[1]*y[2]+x[2]*y[1]);
      A[5]=x[2]*y[2];
      return*this;
    }
    /// add dyadic vector product
    template<template<int,typename> class vector_class>
    SymmMatrix33&add_dyadic(vector_class<3,X> const&x) noexcept
    {
      A[0]+=x[0]*x[0];
      A[1]+=x[0]*x[1];
      A[2]+=x[0]*x[2];
      A[3]+=x[1]*x[1];
      A[4]+=x[1]*x[2];
      A[5]+=x[2]*x[2];
      return*this;
    }
    /// add symmetrised dyadic vector product
    template<template<int,typename> class vector_class>
    SymmMatrix33&add_dyadic_symm(vector_class<3,X> const&x,
				 vector_class<3,X> const&y) noexcept
    {
      A[0]+=x[0]*y[0];
      A[1]+=X(0.5)*(x[0]*y[1]+x[1]*y[0]);
      A[2]+=X(0.5)*(x[0]*y[2]+x[2]*y[0]);
      A[3]+=x[1]*y[1];
      A[4]+=X(0.5)*(x[1]*y[2]+x[2]*y[1]);
      A[5]+=x[2]*y[2];
      return*this;
    }
    /// add dyadic vector product times scalar
    template<template<int,typename> class vector_class>
    SymmMatrix33&add_dyadic(X fac, vector_class<3,X> const&x) noexcept
    {
      X tmp=fac*x[0];
      A[0]+=tmp*x[0];
      A[1]+=tmp*x[1];
      A[2]+=tmp*x[2];
      tmp  =fac*x[1];
      A[3]+=tmp*x[1];
      A[4]+=tmp*x[2];
      A[5]+=fac*x[2]*x[2];
      return*this;
    }
    /// add symmetrised dyadic vector product times scalar
    template<template<int,typename> class vector_class>
    SymmMatrix33&add_dyadic_symm(X fac, vector_class<3,X> const&x,
				        vector_class<3,X> const&y) noexcept
    {
      A[0]+=fac*x[0]*y[0];
      A[1]+=fac*X(0.5)*(x[0]*y[1]+x[1]*y[0]);
      A[2]+=fac*X(0.5)*(x[0]*y[2]+x[2]*y[0]);
      A[3]+=fac*x[1]*y[1];
      A[4]+=fac*X(0.5)*(x[1]*y[2]+x[2]*y[1]);
      A[5]+=fac*x[2]*y[2];
      return*this;
    }
    /// compute to moment of intertia tensor
    template<template<int,typename> class vector_class>
    SymmMatrix33&moment_of_intertia(X mass, vector_class<3,X> const&x)
      noexcept 
    {
      X xx = x[0]*x[0],
	yy = x[1]*x[1],
	zz = x[2]*x[2];
      A[0] = mass * (yy+zz);
      A[1] =-mass * x[0]*x[1];
      A[2] =-mass * x[0]*x[2];
      A[3] = mass * (xx+zz);
      A[4] =-mass * x[1]*x[2];
      A[5] = mass * (xx+yy);
      return*this;
    }
    /// add moment of intertia tensor
    template<template<int,typename> class vector_class>
    SymmMatrix33&add_moment_of_intertia(X mass, vector_class<3,X> const&x)
      noexcept 
    {
      X xx = x[0]*x[0],
	yy = x[1]*x[1],
	zz = x[2]*x[2];
      A[0] += mass * (yy+zz);
      A[1] -= mass * x[0]*x[1];
      A[2] -= mass * x[0]*x[2];
      A[3] += mass * (xx+zz);
      A[4] -= mass * x[1]*x[2];
      A[5] += mass * (xx+yy);
      return*this;
    }
    /// subtract moment of intertia tensor
    template<template<int,typename> class vector_class>
    SymmMatrix33&sub_moment_of_intertia(X mass, vector_class<3,X> const&x)
      noexcept 
    {
      X xx = x[0]*x[0],
	yy = x[1]*x[1],
	zz = x[2]*x[2];
      A[0] -= mass * (yy+zz);
      A[1] += mass * x[0]*x[1];
      A[2] += mass * x[0]*x[2];
      A[3] -= mass * (xx+zz);
      A[4] += mass * x[1]*x[2];
      A[5] -= mass * (xx+yy);
      return*this;
    }
    /// product with vector
    template<template<int,typename> class vector_class>
    vector_class<3,X> operator*(vector_class<3,X> const&x) const noexcept
    {
      return vector_class<3,X>(A[0]*x[0] + A[1]*x[1] + A[2]*x[2],
			       A[1]*x[0] + A[3]*x[1] + A[4]*x[2],
			       A[2]*x[0] + A[4]*x[1] + A[5]*x[2]);
    }
    //@}

  };// class SymmMatrix33
  //
  template<typename X>
  inline Matrix33<X>::Matrix33(SymmMatrix const&m) noexcept
  {
    A[0] =m.A[0]; A[1] =m.A[1]; A[2] =m.A[2];
    A[3] =m.A[1]; A[4] =m.A[3]; A[5] =m.A[4];
    A[6] =m.A[2]; A[7] =m.A[4]; A[8] =m.A[5];
  }
  //
  template<typename X>
  inline Matrix33<X>&Matrix33<X>::operator =(SymmMatrix const&m) noexcept
  {
    A[0] =m.A[0]; A[1] =m.A[1]; A[2] =m.A[2];
    A[3] =m.A[1]; A[4] =m.A[3]; A[5] =m.A[4];
    A[6] =m.A[2]; A[7] =m.A[4]; A[8] =m.A[5];
    return*this;
  }
  //
  template<typename X>
  inline Matrix33<X>&Matrix33<X>::operator+=(SymmMatrix const&m) noexcept
  {
    A[0]+=m.A[0]; A[1]+=m.A[1]; A[2]+=m.A[2];
    A[3]+=m.A[1]; A[4]+=m.A[3]; A[5]+=m.A[4];
    A[6]+=m.A[2]; A[7]+=m.A[4]; A[8]+=m.A[5];
    return*this;
  }
  //
  template<typename X>
  inline Matrix33<X>&Matrix33<X>::operator-=(SymmMatrix const&m) noexcept
  {
    A[0]-=m.A[0]; A[1]-=m.A[1]; A[2]-=m.A[2];
    A[3]-=m.A[1]; A[4]-=m.A[3]; A[5]-=m.A[4];
    A[6]-=m.A[2]; A[7]-=m.A[4]; A[8]-=m.A[5];
    return*this;
  }
  //
  template<typename X>
  inline Matrix33<X> Matrix33<X>::operator+ (SymmMatrix const&m) const noexcept
  {
    Matrix33<X> r;
    r[0]=A[0]+m.A[0]; r[1]=A[1]+m.A[1]; r[2]=A[2]+m.A[2];
    r[3]=A[3]+m.A[1]; r[4]=A[4]+m.A[3]; r[5]=A[5]+m.A[4];
    r[6]=A[6]+m.A[2]; r[7]=A[7]+m.A[4]; r[8]=A[8]+m.A[5];
    return r;
  }
  //
  template<typename X>
  inline Matrix33<X> Matrix33<X>::operator- (SymmMatrix const&m) const noexcept
  {
    Matrix33<X> r;
    r[0]=A[0]-m.A[0]; r[1]=A[1]-m.A[1]; r[2]=A[2]-m.A[2];
    r[3]=A[3]-m.A[1]; r[4]=A[4]-m.A[3]; r[5]=A[5]-m.A[4];
    r[6]=A[6]-m.A[2]; r[7]=A[7]-m.A[4]; r[8]=A[8]-m.A[5];
    return r;
  }
  //
  template<typename X>
  inline SymmMatrix33<X> Matrix33<X>::symm_part() const noexcept
  {
    X half(0.5);
    return SymmMatrix(A[0], half*(A[1]+A[3]), half*(A[2]+A[6]),
		      A[3], half*(A[5]+A[7]),
		      A[8]);
  }
  //
  template<typename X>
  inline Matrix33<X> Matrix33<X>::operator*(SymmMatrix const&m) const noexcept
  {
    return Matrix33(M[0][0]*m.A[0]+M[0][1]*m.A[1]+M[0][2]*m.A[2],
		    M[0][0]*m.A[1]+M[0][1]*m.A[3]+M[0][2]*m.A[4],
		    M[0][0]*m.A[2]+M[0][1]*m.A[4]+M[0][2]*m.A[5],
		    M[1][0]*m.A[0]+M[1][1]*m.A[1]+M[1][2]*m.A[2],
		    M[1][0]*m.A[1]+M[1][1]*m.A[3]+M[1][2]*m.A[4],
		    M[1][0]*m.A[2]+M[1][1]*m.A[4]+M[1][2]*m.A[5],
		    M[2][0]*m.A[0]+M[2][1]*m.A[1]+M[2][2]*m.A[2],
		    M[2][0]*m.A[1]+M[2][1]*m.A[3]+M[2][2]*m.A[4],
		    M[2][0]*m.A[2]+M[2][1]*m.A[4]+M[2][2]*m.A[5]);
  }
  //
  template<typename X>
  inline X Matrix33<X>::trace_of_prod(SymmMatrix const&m) const noexcept
  {
    return M[0][0]*m.A[0] + M[1][1]*m.A[3] + M[2][2]*m.A[5] +
      (M[0][1]+M[1][0])*m.A[1] +
      (M[0][2]+M[2][0])*m.A[2] +
      (M[1][2]+M[2][1])*m.A[4];
  }
  /// \name functionality of Matrix33 and SymmMatrix33
  /// \related Matrix33
  /// \related SymmMatrix33
  //@{
  /// product of vector times matrix
  template<typename X, template<int, typename> class vector_class>
  inline
  vector_class<3,X> operator*(vector_class<3,X> const&x, Matrix33<X> const&M)
    noexcept
  {
    return vector_class<3,X>(M[0][0]*x[0] + M[1][0]*x[1] + M[2][0]*x[2],
			     M[0][1]*x[0] + M[1][1]*x[1] + M[2][1]*x[2],
			     M[0][2]*x[0] + M[1][2]*x[1] + M[2][2]*x[2]);
  }
  /// trace of matrix
  template<typename X> inline X trace(Matrix33<X> const&M) noexcept
  { return M.trace(); }
  /// trace of symmetric matrix
  template<typename X> inline X trace(SymmMatrix33<X> const&M) noexcept
  { return M.trace(); }
  /// trace of square of matrix
  template<typename X> inline X trace_of_square(Matrix33<X> const&M) noexcept
  { return M.trace_of_square(); }
  /// trace of square of symmetric matrix (same as norm)
  template<typename X> inline X trace_of_square(SymmMatrix33<X> const&M)
    noexcept
  { return M.norm(); }
  /// trace of matrix product
  template<typename X> inline X trace_of_prod(Matrix33<X> const&A,
					      Matrix33<X> const&B) noexcept
  { return A.trace_of_prod(B); }
  /// trace of matrix product
  template<typename X> inline X trace_of_prod(Matrix33<X> const&A,
					      SymmMatrix33<X> const&B) noexcept
  { return A.trace_of_prod(B); }
  /// trace of matrix product
  template<typename X> inline X trace_of_prod(SymmMatrix33<X> const&A,
					      SymmMatrix33<X> const&B) noexcept
  { return A.trace_of_prod(B); }
  /// trace of matrix product
  template<typename X> inline X trace_of_prod(SymmMatrix33<X> const&A,
					      Matrix33<X> const&B) noexcept
  { return A.trace_of_prod(B); }
  /// inverse of matrix
  template<typename X> inline Matrix33<X> inverse(Matrix33<X> const&M) noexcept
  { return M.inverse(); }
  /// inverse of symmetric matrix
  template<typename X> inline SymmMatrix33<X> inverse(SymmMatrix33<X> const&M)
    noexcept
  { return M.inverse(); }
  /// determinant of matrix
  template<typename X> inline X det(Matrix33<X> const&M) noexcept
  { return M.det(); }
  /// determinant of symmetric matrix
  template<typename X> inline X det(SymmMatrix33<X> const&M) noexcept
  { return M.det(); }
  /// norm of matrix
  template<typename X> inline X norm(Matrix33<X> const&M) noexcept
  { return M.norm(); }
  /// norm of symmetric matrix
  template<typename X> inline X norm(SymmMatrix33<X> const&M) noexcept
  { return M.norm(); }
  /// shear of matrix
  template<typename X> inline X shear(Matrix33<X> const&M) noexcept
  { return M.shear(); }
  /// shear of symmetric matrix
  template<typename X> inline X shear(SymmMatrix33<X> const&M) noexcept
  { return M.shear(); }
  /// vorticity of matrix
  template<typename X> inline X vorticity(Matrix33<X> const&M) noexcept
  { return M.vorticity(); }
  /// vorticity-squared of matrix
  template<typename X> inline X vorticity_sq(Matrix33<X> const&M) noexcept
  { return M.vorticity_sq(); }
  /// unit matrix as symmetric matrix
  template<typename X> inline SymmMatrix33<X> unity() noexcept
  { SymmMatrix33<X> M; return M.unity(); }
  /// product of vector times symmetric matrix
  template<typename X, template<int,typename> class vector_class>
  inline vector_class<3,X>
  operator*(vector_class<3,X> const&x, SymmMatrix33<X> const&M) noexcept
  { return M*x; }
  /// symmetric part of 3x3 matrix
  template<typename X> inline SymmMatrix33<X>
  symm_part(Matrix33<X> const&M) noexcept
  { return M.symm_part(); }
  //@}

  /// construct rotation matrix
  /// \param[in] u   unit vector: rotation axis
  /// \param[in] a   rotation angle [radian] about @a axis
  /// \return matrix for rotating 3D vectors by angle @a around axis @a u
  /// \note We assume but don't check that @a u is a unit vector
  template<typename X, template<int,typename> class vector_class>
  inline Matrix33<X> rotation_matrix_na(vector_class<3,X> const&u, X a)
    noexcept
  {
    Matrix33<X> R; R.dyadic(u);
    X c   = std::cos(a);
    X s   = std::sin(a);
    R    *= 1-c;
    R(0) += c;
    R(1) -= u[2]*s;
    R(2) += u[1]*s;
    R(3) += u[2]*s;
    R(4) += c;
    R(5) -= u[0]*s;
    R(6) -= u[1]*s;
    R(7) += u[0]*s;
    R(8) += c;
    return R;
  }
  /// construct rotation matrix
  /// \param[in] axis   rotation axis, not necessarily normalised
  /// \param[in] angle  rotation angle [radian] about @a axis
  /// \return matrix for rotating 3D vectors by @a angle around @a axis
  template<typename X, template<int,typename> class vector_class>
  inline Matrix33<X> rotation_matrix(vector_class<3,X> const&axis, X angle)
    noexcept
  { return rotation_matrix_na(axis.normalized(),angle); }
} // namespace WDutils
//
#undef noexcept 
//
#endif // WDutils_included_symm33_h
