// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/symm33.h
///
/// \author Walter Dehnen
///                                                                             
/// \date   2010
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2010 Walter Dehnen
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
#ifndef WDutils_included_symm33_h
#define WDutils_included_symm33_h

#ifndef WDutils_included_exception_h
#  include <exception.h>
#endif
#ifndef WDutils_included_tupel_h
#  include <tupel.h>
#endif
#ifndef WDutils_included_sse_h
#  include <sse.h>
#endif

namespace WDutils {

  /// functionality for symmetric 3x3 matrices
  ///
  /// \note I just started this, will be adding functionality as needed.
  template<typename __T>
  struct Symm33 {
    typedef __T           real;          ///<  scalar type
    typedef tupel<3,real> vector;        ///<  associated vector type
    typedef real          matrix[6];     ///<  symmetric 3x3 matrix
    /// assign matrix
    static void assign(const matrix&M, matrix&R)
    { R[0]=M[0]; R[1]=M[1]; R[2]=M[2]; R[3]=M[3]; R[4]=M[4]; R[5]=M[5]; }
    /// add matrix
    static void add(const matrix&M, matrix&R)
    { R[0]+=M[0]; R[1]+=M[1]; R[2]+=M[2]; R[3]+=M[3]; R[4]+=M[4]; R[5]+=M[5]; }
    /// subtract matrix
    static void sub(const matrix&M, matrix&R)
    { R[0]-=M[0]; R[1]-=M[1]; R[2]-=M[2]; R[3]-=M[3]; R[4]-=M[4]; R[5]-=M[5]; }
    /// multiply matrix with scalar
    static void prod(matrix&M, real fac)
    { M[0]*=fac; M[1]*=fac; M[2]*=fac; M[3]*=fac; M[4]*=fac; M[5]*=fac; }
    /// product of matrix with scalar
    static void prod(const matrix&M, real f, matrix&R)
    { R[0]=f*M[0]; R[1]=f*M[1]; R[2]=f*M[2];
      R[3]=f*M[3]; R[4]=f*M[4]; R[5]=f*M[5]; }
    /// matrix inverse
    //  20*,8+,1/
    static void invert(const matrix&M, matrix&Mi)
    {
      real tmp;
      Mi[0] = M[3]*M[5]-M[4]*M[4];
      tmp   = M[2]*M[4];
      Mi[1] = tmp      -M[1]*M[5];
      Mi[2] = M[1]*M[4]-M[2]*M[3];
      Mi[3] = M[0]*M[5]-tmp;
      tmp   = M[3]*Mi[3];
      Mi[4] = M[1]*M[2]-M[0]*M[4];
      tmp  += M[4]*Mi[4];
      Mi[5] = M[0]*M[3]-M[1]*M[1];
      tmp  += M[5]*Mi[5];
      tmp   = real(1)/tmp;
      Mi[0]*= tmp;
      Mi[1]*= tmp;
      Mi[2]*= tmp;
      Mi[3]*= tmp;
      Mi[4]*= tmp;
      Mi[5]*= tmp;
    }
    /// symmetrised matrix product
    /// \note In general, the product of two symmetric matrices is not
    ///       symmetric. Here, we return the symmetric part of that product.
    //  27*,18+
    static void prod(const matrix&A, const matrix&B, matrix&AB)
    {
      real t11=A[1]*B[1], t22=A[2]*B[2], t44=A[4]*B[4];
      AB[0] = A[0]*B[0]+t11      +t22;
      AB[3] = t11      +A[3]*B[3]+t44;
      AB[5] = t22      +t44      +A[5]*B[5];
      AB[1] = real(0.5)*(A[0]*B[1]+A[1]*B[3]+A[2]*B[4] +
			 A[1]*B[0]+A[3]*B[1]+A[4]*B[2] );
      AB[2] = real(0.5)*(A[0]*B[2]+A[1]*B[4]+A[2]*B[5] +
			 A[2]*B[0]+A[4]*B[1]+A[5]*B[2] );
      AB[4] = real(0.5)*(A[1]*B[2]+A[3]*B[4]+A[4]*B[5] +
			 A[2]*B[1]+A[4]*B[3]+A[5]*B[4] );
    }
    /// add symmetrised matrix product
    /// \note In general, the product of two symmetric matrices is not
    ///       symmetric. Here, we add the symmetric part of that product.
    //  27*,24+
    static void add_prod(const matrix&A, const matrix&B, matrix&R)
    {
      real t11=A[1]*B[1], t22=A[2]*B[2], t44=A[4]*B[4];
      R[0] += A[0]*B[0]+t11      +t22;
      R[3] += t11      +A[3]*B[3]+t44;
      R[5] += t22      +t44      +A[5]*B[5];
      R[1] += real(0.5)*(A[0]*B[1]+A[1]*B[3]+A[2]*B[4] +
			 A[1]*B[0]+A[3]*B[1]+A[4]*B[2] );
      R[2] += real(0.5)*(A[0]*B[2]+A[1]*B[4]+A[2]*B[5] +
			 A[2]*B[0]+A[4]*B[1]+A[5]*B[2] );
      R[4] += real(0.5)*(A[1]*B[2]+A[3]*B[4]+A[4]*B[5] +
			 A[2]*B[1]+A[4]*B[3]+A[5]*B[4] );
    }
    /// square of symmetric matrix
    //  15*,12+
    static void square(const matrix&A, matrix&AA)
    {
      real t11=A[1]*A[1], t22=A[2]*A[2], t44=A[4]*A[4];
      AA[0] = A[0]*A[0]+t11      +t22;
      AA[3] = t11      +A[3]*A[3]+t44;
      AA[5] = t22      +t44      +A[5]*A[5];
      AA[1] = A[0]*A[1]+A[1]*A[3]+A[2]*A[4];
      AA[2] = A[0]*A[2]+A[1]*A[4]+A[2]*A[5];
      AA[4] = A[1]*A[2]+A[3]*A[4]+A[4]*A[5];
    }
    /// matrix times vector
    static void prod(const matrix&A, vector const&x, vector&xA)
    {
      xA[0] = A[0]*x[0]+A[1]*x[1]+A[2]*x[2];
      xA[1] = A[1]*x[0]+A[3]*x[1]+A[4]*x[2];
      xA[2] = A[2]*x[0]+A[4]*x[1]+A[5]*x[2];
    }
    /// trace of matrix
    static real trace(const matrix&A)
    { return A[0]+A[3]+A[5]; }
    /// remove trace carrying part from matrix
    static void traceless(matrix&A)
    {
      real tmp = real(0.33333333333333333) * trace(A);
      A[0] -= tmp;
      A[3] -= tmp;
      A[5] -= tmp;
    }
    /// matric norm
    static real norm(const matrix&A)
    {
      return  A[0]*A[0]+A[3]*A[3]+A[5]*A[5]+
	twice(A[1]*A[1]+A[2]*A[2]+A[4]*A[4]);
    }
    /// return element (i,j)
    /// \note this is only provided for completeness
    static real const&element(const matrix&M, int i, int j)
    { return M[K(i,j)]; }
    /// return element (i,j)
    /// \note this is only provided for completeness
    static real &element(matrix&M, int i, int j)
    { return M[K(i,j)]; }
  private:
    static int K(int i, int j)
    {
      const static int Ktab[3][3] ={{0,1,2},{1,3,4},{2,4,5}};
      return Ktab[i][j];
    }
  };// struct Symm33
#if(0)  // still under construction

#ifdef __SSE__
  /// specialisation for single precision in the case of SSE support
  ///
  /// \note I just started this, will be adding functionality as needed.
  template<>
  struct Symm33<float> {
    /// type used for symmetrix 3x3 matrix.
    /// \note the last 2 number are not used.
    typedef WDutils__align16 float matrix[8];



    /// return element (i,j)
    /// \note this is only provided for completeness
    static float const&element(const matrix&m, int i, int j)
    { return M[Symm33<int>::K(i,j)]; }
    /// return element (i,j)
    /// \note this is only provided for completeness
    static float&element(matrix&m, int i, int j)
    { return M[Symm33<int>::K(i,j)]; }
  };// struct Symm33<float>
#endif  // __SSE__

#ifdef __SSE2__
  /// specialisation for double precision in the case of SSE support
  ///
  /// \note I just started this, will be adding functionality as needed.
  template<>
  struct Symm33<double> {
    /// type used for symmetrix 3x3 matrix.
    typedef WDutils__align16 double matrix[6];



    /// return element (i,j)
    /// \note this is only provided for completeness
    static double const&element(const matrix&m, int i, int j)
    { return M[Symm33<int>::K(i,j)]; }
    /// return element (i,j)
    /// \note this is only provided for completeness
    static double&element(matrix&m, int i, int j)
    { return M[Symm33<int>::K(i,j)]; }
  };// struct Symm33<float>
#endif  // __SSE__
#endif


} // namespace WDutils
//
#endif // WDutils_included_symm33_h
