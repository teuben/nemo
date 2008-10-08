// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    utils/inc/numerics.h                                               
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    1994-2008                                                          
///                                                                             
/// \todo    finish doxygen documentation                                       
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 1994-2008  Walter Dehnen                                       
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
//                                                                              
// CONTENTS                                                                     
//                                                                              
// sorting                                              v0.2 & v0.3             
// find position in ordered table                       v0.0                    
// polynomial interpolation on grids                    v0.0                    
// root finding                                         v0.0                    
// bracket a minimum                                    v0.0                    
// function minimization                                v0.0                    
// Burlisch-Stoer integration of 1D real integrals      v0.0                    
// Runge-Kutte 4th order integrator                     v0.0                    
// Legendre polynomials and their derivatives           v0.1                    
// cubic spline (#included from walter/spln.h)          v0.1                    
// Gauss-Legendre integration: points & weights         v0.1                    
// eigenvalues and vectors of symmetric matrices        v0.4                    
//                                                                              
////////////////////////////////////////////////////////////////////////////////
#ifndef WDutils_included_numerics_h
#define WDutils_included_numerics_h

#ifndef WDutils_included_iostream
#  include <iostream>
#  define WDutils_included_iostream
#endif
#ifndef WDutils_included_cstdlib
#  include <cstdlib>
#  define WDutils_included_cstdlib
#endif
#ifndef WDutils_included_inline_h
#  include <inline.h>
#endif
#ifndef WDutils_included_memory_h
#  include <memory.h>
#endif
#ifndef WDutils_included_tupel_h
#  include <tupel.h>
#endif

////////////////////////////////////////////////////////////////////////////////
namespace WDutils {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// \name finding position in ordered table                                   
  //@{                                                                          
  //----------------------------------------------------------------------------
  /// find position in ordered table (see NR)                                   
  ///                                                                           
  /// \param T (template parameter) type with < and > operators, usually a
  ///          scalar
  /// \param[in] xarr  array of T, must be ordered (ascending or descending)
  /// \param[in] n     size of array xarr
  /// \param[in] x     value to find position fo
  /// \param[in] j     initial guess for position
  /// \return position jlo such that xarr[jlo] <= x < xarr[jlo+1]
  ///
  /// The ordered table xarr is hunted for jlo such that
  ///    xarr[jlo] <= x < xarr[jlo+1].
  /// For an ascendingly ordered array, we return -1 if x < xarr[0], n-1 if
  /// x == xarr[n-1], and n if x > xarr[n--1].
  template<typename T>
  int hunt(const T*xarr, int n, T x, int j) {
    int  jlo=j,jhi,l=n-1;
    bool ascnd=(xarr[l]>xarr[0]);
    if(!ascnd && xarr[l]==xarr[0] ) return -1;	    // x_0 = x_l                
    if( ascnd && x<xarr[0] || !ascnd && x>xarr[0] ) return -1;
    if( ascnd && x>xarr[l] || !ascnd && x<xarr[l] ) return  n;

    if(jlo<0 || jlo>l) {                            // input guess not useful,  
      jlo = -1;                                     //   go to bisection below  
      jhi = n;
    } else {
      int inc = 1;
      if(x>=xarr[jlo] == ascnd) {                   // hunt upward              
	if(jlo == l) return (x==xarr[l])? l : n;
	jhi = jlo+1;
	while(x>=xarr[jhi] == ascnd) {              // not done hunting         
	  jlo =jhi;
	  inc+=inc;                                 // so double the increment  
	  jhi =jlo+inc;
	  if(jhi>l) {                               // off end of table         
	    jhi=n;
	    break;
	  }
	}
      } else {                                      // hunt downward            
	if(jlo == 0) return ascnd? -1 : 0;
	jhi = jlo;
	jlo-= 1;
	while(x<xarr[jlo] == ascnd) {               // not done hunting         
	  jhi = jlo;
	  inc+= inc;                                // so double the increment  
	  jlo = jhi-inc;
	  if(jlo < 0) {                             // off end of table         
	    jlo = 0;
	    break;
	  }
	}
      }
    }
    while (jhi-jlo != 1) {                          // bisection phase          
      int jm=(jhi+jlo) >> 1;
      if(x>=xarr[jm] == ascnd) jlo=jm;
      else jhi=jm;
    }
    return jlo;
  }
  //----------------------------------------------------------------------------
  /// find position in ordered table (see NR)
  ///
  /// \param T (template param) type with < and > operators, usually a scalar
  /// \param[in,out] k  jlo such that xarr[jlo] <= x < xarr[jlo+1]
  /// \param[in]     x  array of T, must be ordered (ascending or descending) 
  /// \param[in]     n  size of array xarr
  /// \param[in]     xi value to find position for
  /// \return        position jlo such that xarr[jlo] <= x < xarr[jlo+1]
  ///
  /// If the original value for k already gives the position, we return.
  /// Otherwise, we guess k from linear interpolation and then invoke hunt().
  /// If x is not in range, we throw an error.
  template<typename T>
  inline void find(int&k, int n, const T*x, T xi) WDutils_THROWING
  {
    if(k<0 || k>=n-1 || x[k]>xi || x[k+1]<xi) {
      k = int( (xi-x[0]) / (x[n-1]-x[0]) * (n-1) );
      k = hunt(x,n,xi,k);
      if(k<0 || k>=n) 
	WDutils_THROW("find(): x=%f out of range [%f,%f]\n", xi,x[0],x[n-1]);
    }
  }
  //@}
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// \name sorting and related                                                 
  //{@                                                                          
  //----------------------------------------------------------------------------
  /// the numbers 0 to n-1 are ordered in ascending order of A[i]
  /// using the heap-sort algorithm; based on a routine given in NR.
  /// \note a routine HeapSort() can be found in file heap.h.
  /// \param sortable type of values to be sorted
  /// \param sortit class with operator[](int) returning sortable
  /// \param A (input) array or values to be sorted
  /// \param n (input) number of elements
  /// \param indx (output) index table so that A[indx[i]] is sorted
  template<typename sortable, class sortit>
  void HeapIndex(const sortit&A, int n, int*indx)
  {
    if(n<=0) return;
    if(n==1) { indx[0]=0; return; }
    for(int j=0; j!=n; ++j) indx[j] = j;
    int l = n>>1;
    int ir= n-1;
    sortable q;
    for(;;) {
      int indxt;
      if(l>0)
	q = A[indxt=indx[--l]];
      else {
	q = A[indxt=indx[ir]];
	indx[ir] = indx[0];
	if(--ir == 0) {
	  indx[0] = indxt;
	  return;
	}
      }
      int i = l;
      int j = (l<<1) + 1;
      while(j<=ir) {
	if(j < ir && A[indx[j]] < A[indx[j+1]] ) j++;
	if(q < A[indx[j]] ) {
	  indx[i] = indx[j];
	  j+= 1+(i=j);
	} else
	  j = ir+1;
      }
      indx[i] = indxt;
    }
  }
  //----------------------------------------------------------------------------
  /// the numbers 0 to n-1 are ordered in ascending order of A[i]
  /// using the heap-sort algorithm; based on a routine given in NR.
  /// \note a routine HeapSort() can be found in file heap.h.
  /// \param sortable type of values to be sorted
  /// \param A (input) array or values to be sorted
  /// \param n (input) number of elements
  /// \param indx (output) index table so that A[indx[i]] is sorted
  template<typename sortable>
  void HeapIndex(const sortable*A, int n, int*indx)
    // based on a routine given in NR                                           
    // the numbers 0 to n-1 are ordered in ascending order of A[i]              
  {
    if(n<=0) return;
    if(n==1) { indx[0]=0; return; }
    for(int j=0; j!=n; ++j) indx[j] = j;
    int l = n>>1;
    int ir= n-1;
    for(;;) {
      const sortable*q;
      int indxt;
      if(l>0)
	q = A+(indxt=indx[--l]);
      else {
	q = A+(indxt=indx[ir]);
	indx[ir] = indx[0];
	if(--ir == 0) {
	  indx[0] = indxt;
	  return;
	}
      }
      int i = l;
      int j = (l<<1) + 1;
      while(j<=ir) {
	if(j  < ir && A[indx[j]] < A[indx[j+1]] ) j++;
	if(*q < A[indx[j]] ) {
	  indx[i] = indx[j];
	  j+= 1+(i=j);
	} else
	  j = ir+1;
      }
      indx[i] = indxt;
    }
  }
  //----------------------------------------------------------------------------
  /// the numbers 0 to n-1 are ordered in ascending order of A[i]
  /// using the heap-sort algorithm; based on a routine given in NR
  /// \param sortable type of values to be sorted
  /// \param A (input) Array or values to be sorted
  /// \param I (output) index table so that A[I[i]] is sorted
  template<typename sortable>
  void HeapIndex(Array<sortable,1>const&A, Array<int,1>&I) WDutils_THROWING {
    if(A.size() != I.size()) WDutils_THROW("size mismatch in HeapIndex()\n");
    HeapIndex(A.array(),A.size(),I.array());
  }
  //----------------------------------------------------------------------------
  /// given an array of values, produce a table of their ranks
  template<typename sortable>
  void HeapRank(const sortable*A, int n, int*rank)
  {
    if(n<=0) return;
    int*indx = WDutils_NEW(int,n);
    HeapIndex(A,n,indx);
    for(int r=0; r!=n; ++r) rank[indx[r]] = r;
    WDutils_DEL_A(indx);
  }
  //----------------------------------------------------------------------------
  /// given an array of values, produce a table of their ranks
  template<typename sortable>
  void HeapRank(Array<sortable,1>const&A, Array<int,1>&I) WDutils_THROWING {
    if(A.size() != I.size()) WDutils_THROW("size mismatch in HeapRank()\n");
    HeapRank(A.array(),A.size(),I.array());
  }
  //----------------------------------------------------------------------------
  /// find percentiles of a 1D distribution of values                           
  /// instantinations for float and double                                      
  template<typename scalar> class FindPercentile {
    void *DATA;
    /// copy constructor disabled; use references instead of copies
    FindPercentile(const FindPercentile&);
    void setup(const scalar*, int, const scalar*);
  public:
    /// ctor: setup sort tree
    /// \param[in] F  array with elements
    /// \param[in] N  number of elements
    /// \param[in] W  (optional) array of weights
    FindPercentile(const scalar*F, int N, const scalar*W=0) : DATA(0) {
      setup(F,N,W);
    }
    /// ctor: setup sort tree
    /// \param[in] F  Array with elements
    FindPercentile(Array<scalar,1>const&F) : DATA(0) {
      setup(F.array(),F.size(),0);
    }
    /// ctor: setup sort tree
    /// \param[in] F  Array<> with elements
    /// \param[in] W  Array<> of weights
    FindPercentile(Array<scalar,1>const&F,
		   Array<scalar,1>const&W) WDutils_THROWING : DATA(0) {
      if(F.size() != W.size())
	WDutils_THROW("size mismatch in FindPercentile\n");
      setup(F.array(),F.size(),W.array());
    }
    /// dtor: de-allocate sorttree
    ~FindPercentile();
    /// \param f percentile note that P in [0,1]
    /// \return value corresponding to cumulative weight f*W_total
    scalar Percentile(scalar f) const;
    /// \param r rank in [0,N]
    /// \return value corresponding to rank r
    scalar Percentile(int r) const;
  };
  //----------------------------------------------------------------------------
  /// just find a single percentile                                             
  /// \return percentile value: r values are less than this
  /// \param r rank
  /// \param F array with values
  /// \param N size of array
  template<typename scalar>
  scalar percentile(int r, const scalar*F, int N) {
    FindPercentile<scalar> FP(F,N);
    return FP.Percentile(r);
  }
  //----------------------------------------------------------------------------
  /// just find a single percentile                                             
  /// \return percentile value: a fraction f of the total weight is at y<F
  /// \param f fraction
  /// \param F array with values
  /// \param N size of array
  /// \param W array with weights
  template<typename scalar>
  scalar percentile(scalar f, const scalar*F, int N, const scalar*W)
  {
    FindPercentile<scalar> FP(F,N,W);
    return FP.Percentile(f);
  }
  //----------------------------------------------------------------------------
  /// find index given ranks
  /// instantinations for float and double
  template<typename scalar> class FindRank {
    void *DATA;
    /// copy constructor disabled; use references instead of copies
    FindRank(const FindRank&);
    void setup(const scalar*, unsigned);
  public:
    /// ctor: setup sort tree
    /// \param[in] F  array with elements
    /// \param[in] N  number of elements
    FindRank(const scalar*F, unsigned N) : DATA(0) {
      setup(F,N);
    }
    /// ctor: setup sort tree
    /// \param[in] F  Array with elements
    FindRank(Array<scalar,1>const&F) : DATA(0) {
      setup(F.array(),F.size());
    }
    /// dtor: de-allocate sorttree
    ~FindRank();
    /// find index given rank
    /// \param[in] rank  rank in [0,N-1] for which index is required
    /// \return    index of element with that rank.
    unsigned Index(unsigned rank) const;
    /// find value given rank
    /// \param[in] rank  rank in [0,N-1] for which value is required
    /// \return    value of element with that rank
    scalar Value(unsigned rank) const;
  };
  //@}
  //----------------------------------------------------------------------------
  /// class to support reporting file and line number on error
  class FileLineFind {
  protected:
    const char* file;
    const int   line;
    /// find the index of the first point to be used in interpolation
    /// \param[in] n size of table == total number of points
    /// \param[in] m number of points used for interpolation
    /// \param[in] x table of points
    /// \param[in] xi point to be interpolated at
    /// \param[in,out] j index of first point to be used in interpolation
    /// \return number of points required in interpolation (1 or min{m,n})
    template<typename X>
    static int find(int&j, int n, int m, const X*x, X xi)
    {
      int M=m<n? m:n;
      j = int( (xi-x[0]) / (x[n-1]-x[0]) * (n-1) );
      j = hunt(x,n,xi,j) - (M+1)/2 + 1;
      if(j>=0 && j<n && x[j]==xi) M = 1;
      else if(j<0)		  j = 0;
      else if(j>n-M)		  j = n-M;
      return M;
    }
    FileLineFind() : file(0), line(0) {}
    /// constructor: take file and line
    FileLineFind(const char*f, int l) : file(f), line(l) {}
  };
  //----------------------------------------------------------------------------
  /// \name polynomial interpolation in 1D
  //@{
  /// supporting macro Polev.
  class PolynomialEvaluation : private FileLineFind {
    /// polynomial interpolation using n values; adapted from NR
    /// \param[in] x  array of points
    /// \param[in] y  array of values
    /// \param[in] n  number of points
    /// \param[in] P  auxialiary array of size n
    /// \param[in] xi position to be interpolated at
    /// \return y(xi) as interpolated
    template<typename X, typename Y>
    Y polint(int n, const X*x, const Y*y, Y*P, X xi) WDutils_THROWING
    {
      for(int i=0;i!=n;++i)
	P[i]=y[i];
      for(int m=1;m!=n;++m)
	for(int i=0;i<n-m;++i) {
	  if(x[i]==x[i+m]) {
	    if(file)
	      WDutils_THROWN("[%s:%d]: x's not distinct in Polev()",file,line);
	    else
	      WDutils_THROW ("x's not distinct in polev()");
	  }
	  P[i]= ( (xi-x[i+m])*P[i] + (x[i]-xi)*P[i+1] ) / (x[i] - x[i+m]);
	}
      return P[0];
    }
    //..........................................................................
  public:
    /// default constructor
    PolynomialEvaluation() : FileLineFind() {}
    /// constructor: take file and line
    PolynomialEvaluation(const char*f, int l) : FileLineFind(f,l) {}
    /// polynomial interpolation using m of n values
    /// \param n  total size of arrays
    /// \param m  number of points used in interpolation
    /// \param xi position to find function value at
    /// \param x  array of points
    /// \param y  array of values
    /// \note Together with the macro polev this implements the function
    /// polev() for given m.
    template<typename X, typename Y>
    Y operator()(X xi, const X*x, const Y*y, int n, int m) WDutils_THROWING
    {
      int j;
      Y*P = WDutils_NEW(Y,n);
      Y yi= find(j,n,m,x,xi)==1? y[j] : polint(m,x+j,y+j,P,xi);
      WDutils_DEL_A(P);
      return yi;
    }
    /// polynomial interpolation using 4 of n values
    /// \param n  total size of arrays
    /// \param xi position to find function value at
    /// \param x  array of points
    /// \param y  array of values
    /// \note Together with the macro polev this implements the function
    /// polev() for m=4
    template<typename X, typename Y> inline
    Y operator()(X xi, const X*x, const Y*y, int n) WDutils_THROWING
    {
      int j;
      Y P[4];
      return find(j,n,4,x,xi)==1? y[j] : polint(4,x+j,y+j,P,xi);
    }
    /// polynomial interpolation using 4 of n values, taking Array<T> arguments
    /// \param xi position to find function value at
    /// \param x  array of points
    /// \param y  array of values
    /// \note Together with the macro polev this implements the function
    /// polev() for m=4 and Array<> arguments
    template<typename X, typename Y>
    Y operator()(X xi, const Array<X,1>&x, const Array<Y,1>&y) WDutils_THROWING
    {
      if(x.size() != y.size())
	if(file)
	  WDutils_THROWN("[%s:%d]: Array size mismatch in Polev()",file,line);
	else
	  WDutils_THROW ("Array size mismatch in polev()");
      return operator()(xi,x.array(),y.array(),x.size());
    }
  };
  /// macro Polev: implements functions Polev() like polev().
  /// The idea is to implement the "functions" Polev() via a macro such that
  /// on error the file and line of the call to Polev() can be reported. \n
  /// The trick is simple: the macro expands code like
  /// \code Polev(xi,x,y,n); \endcode into
  /// \code PolynomialEvaluation(__FILE__,__LINE__)(xi,x,y,n); \endcode
  /// The first argument list invokes the constructor and the second the
  /// operator() members of class PolynomialEvaluation.\n
  /// \note We provide ordinary functions polev() within the namespace WDutils
  ///       with the same semantics as Polev(). The only difference being (i)
  ///       the error reporting (with file and line number in Polev()) and (ii)
  ///       the fact that Polev as a macro is in the global namespace
#define Polev WDutils::PolynomialEvaluation(__FILE__,__LINE__)
  /// polynomial interpolation using m of n values
  /// \param[in] xi position to find function value at
  /// \param[in] x  array of points
  /// \param[in] y  array of values
  /// \param[in] n  total size of arrays
  /// \param[in] m  number of points used in interpolation
  template<typename X, typename Y>
  inline Y polev(X xi, const X*x, const Y*y, int n, int m) WDutils_THROWING {
    return PolynomialEvaluation(0,0)(xi,x,y,n,m);
  }
  /// polynomial interpolation using 4 of n values
  /// \param[in] xi position to find function value at
  /// \param[in] x  array of points
  /// \param[in] y  array of values
  /// \param[in] n  total size of arrays
  template<typename X, typename Y>
  inline Y polev(X xi, const X*x, const Y*y, int n) WDutils_THROWING {
    return PolynomialEvaluation(0,0)(xi,x,y,n);
  }
  /// polynomial interpolation using 4 of n values, taking Array<T> arguments
  /// \param[in] xi position to find function value at
  /// \param[in] x  array of points
  /// \param[in] y  array of values
  template<typename X, typename Y>
  inline Y polev(X xi, const Array<X,1>&x, const Array<Y,1>&y) WDutils_THROWING{
    return PolynomialEvaluation()(xi,x,y);
  }
  //----------------------------------------------------------------------------
  /// like polev(), but no extrapolation; gives boundary values instead
  template<typename X, typename Y>
  inline Y ipolev(X xi, const X*x, const Y*y, int n, int m)
  {
    if(xi < x[0]   && x[0]   < x[n-1]) return y[0];
    if(xi > x[0]   && x[0]   > x[n-1]) return y[0];
    if(xi > x[n-1] && x[n-1] > x[0]  ) return y[n-1];
    if(xi < x[n-1] && x[n-1] < x[0]  ) return y[n-1];
    return polev(xi,x,y,n,m);
  }
  //----------------------------------------------------------------------------
  /// like polev(), but no extrapolation; gives boundary values instead
  template<typename X, typename Y>
  inline Y ipolev(X xi, const X*x, const Y*y, int n)
  {
    if(xi < x[0]   && x[0]   < x[n-1]) return y[0];
    if(xi > x[0]   && x[0]   > x[n-1]) return y[0];
    if(xi > x[n-1] && x[n-1] > x[0]  ) return y[n-1];
    if(xi < x[n-1] && x[n-1] < x[0]  ) return y[n-1];
    return polev(xi,x,y,n);
  }
  //@}
  //----------------------------------------------------------------------------
  /// \name root finding
  //@{
  /// encodes NR rtsafe
  class RootSafe: private FileLineFind {
  public:
    RootSafe() : FileLineFind() {}
    RootSafe(const char*f, int l) : FileLineFind(f,l) {}
    template <typename X>
    X operator() (void(*func)(X,X&,X&), X x1, X x2, X xacc)
      throw(WDutils::exception)
    {
      const int maxit=100;
      X xh,xl,dx,dxo,f,df,fh,fl,rts,temp;
      func(x1,fl,df);
      func(x2,fh,df);
      if(fl*fh >= 0.) {
	if(file)
	  throw exception("[%s.%d]: root must be bracketed in Rtsafe()",
			  file,line);
	else
	  throw exception("root must be bracketed rtsafe()");
      }
      if(fl<0.) {
	xl  = x1;
	xh  = x2;
      } else {
	xh  = x1;
	xl  = x2;
	temp= fl;
	fl  = fh;
	fh  = temp;
      }
      rts = 0.5*(x1+x2);
      dxo = abs(x2-x1);
      dx  = dxo;
      func(rts,f,df);
      for(int j=0; j!=maxit; ++j) {
	if((((rts-xh)*df-f)*((rts-xl)*df-f)>= 0.) || (abs(2.*f)>abs(dxo*df))) {
	  dxo = dx;
	  dx  = 0.5*(xh-xl);
	  rts = xl+dx;
	  if(xl==rts) return rts;
	} else {
	  dxo = dx;
	  dx  = f/df;
	  temp=rts;
	  rts-= dx;
	  if(temp==rts) return rts;
	}
	if(abs(dx)<xacc) return rts;
	func(rts,f,df);
	if(f<0.) {
	  xl  = rts;
	  fl  = f;
	} else {
	  xh  = rts;
	  fh  = f;
	}
      }
      if(file)
	throw exception("[%s.%d]: "
			"maximum number of %d iterations exceeded in Rtsafe()",
			file,line,maxit);
      else
	throw exception("maximum number of %d iterations exceeded in rtsafe()",
			maxit);
      return rts;
    }
  };
  /// like NR rtsafe
  /// \return root: x value at which f(x)=0
  /// \param[in] func function returning f & df/dx (args 2&3) given x (1st arg)
  /// \param[in] x1   left boundary of interval
  /// \param[in] x2   right boundary of interval
  /// \param[in] xacc desired accuracy for root
  /// \note f(x1)*f(x2) \b must not be positive on input
  template <typename X>
  X rtsafe(void(*func)(X,X&,X&), X x1, X x2, X xacc) throw(WDutils::exception)
  {
    return RootSafe()(func,x1,x2,xacc);
  }
  /// macro with same functionality as function rtsafe() above, except that on
  /// error file and line number of the call are reported.
#define Rtsafe WDutils::RootSafe(__FILE__,__LINE__)
  //@}
  //----------------------------------------------------------------------------
  // bracketing a minimum                                                       
  //----------------------------------------------------------------------------
  template<typename scalar_type>
  void minbrak(scalar_type& ax,                    // I/O ax                    
	       scalar_type& bx,                    // I/O bx                    
	       scalar_type& cx,                    // O:  cx                    
	       scalar_type(*f)(const scalar_type)) // I: f(x)                   
  {
    const scalar_type GLIMIT=100.0, GOLD=1.618034, TINY=1.e-20;
    scalar_type ulim,u,r,q,fu,fa=f(ax),fb=f(bx),fc;
    if(fb<fa) { swap(ax,bx); swap(fa,fb); }
    cx = bx+GOLD*(bx-ax);
    fc = f(cx);
    while(fb > fc) {
      r    = (bx-ax)*(fb-fc);
      q    = (bx-cx)*(fb-fa);
      u    = bx-((bx-cx)*q-(bx-ax)*r)/(2*sign(max(abs(q-r),TINY),q-r));
      ulim = bx+GLIMIT*(cx-bx);
      if((bx-u)*(u-cx) > 0.0) {                    // try parabolic u in [b,c]  
	fu = f(u);
	if(fu < fc) {
	  ax = bx;
	  bx = u;
	  return;
	} else if(fu>fb) {
	  cx = u;
	  return;
	}
	u = cx + GOLD * (cx-bx);
	fu= f(u);
      } else if((cx-u)*(u-ulim) > 0.0) {
	fu = f(u);
	if(fu < fc) {
#define SHFT(a,b,c,d) (a)=(b); (b)=(c); (c)=(d);
	  SHFT(bx,cx,u,cx+GOLD*(cx-bx))
	  SHFT(fb,fc,fu,f(u))
	}
      } else if((u-ulim)*(ulim-cx) >= 0.0) {
	u  = ulim;
	fu = f(u);
      } else {
	u  = cx + GOLD * (cx-bx);
	fu = f(u);
      }
      SHFT(ax,bx,cx,u)
      SHFT(fa,fb,fc,fu)
#undef  SHFT
    }
  }
  //----------------------------------------------------------------------------
  // function minimization                                                      
  //----------------------------------------------------------------------------
  template<typename scalar_type> scalar_type
  brent(                                           // R: f_min = f(x_min)       
	scalar_type ax,                            // I: ax   these bracket the 
	scalar_type bx,                            // I: bx   minimum, e.g. out-
	scalar_type cx,                            // I: cx   put of minbrak()  
	scalar_type(*f)(scalar_type),              // I: function f(x)          
	scalar_type tol,                           // I: accuracy wanted        
	scalar_type& xmin)                         // O: x_min                  
  {
    const int   itmax=100;
    const scalar_type cgold=0.381966, zeps=1.e-10;
    int         iter;
    scalar_type a,b,d=0.,e=0.,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
    a = (ax<cx) ? ax:cx;
    b = (ax>cx) ? ax:cx;
    x = w = v = bx;
    fw= fv= fx= f(x);
    for(iter=0; iter!=itmax; ++iter) {
      xm   = 0.5*(a+b);
      tol1 = tol*abs(x)+zeps;
      tol2 = 2.0* tol1;
      if(abs(x-xm) <= (tol2-0.5*(b-a))) {
	xmin = x;
	return fx;
      }
      if(abs(e) > tol1) {
	r = (x-w) * (fx-fv);
	q = (x-v) * (fx-fw);
	p = (x-v) * q - (x-w) * r;
	q = 2 * (q-r);
	if(q>0.) p =-p;
	q     = abs(q);
	etemp = e;
	e     = d;
	if(abs(p)>=abs(0.5*q*etemp) || p<=q*(a-x) || p>=q*(b-x))
	  d = cgold * (e = (x>=xm) ? a-x : b-x);
	else {
	  d = p / q;
	  u = x + d;
	  if ((u-a)<tol2 || (b-u) < tol2) d=sign(tol1,xm-x);
	}
      } else
	d = cgold * (e= (x>=xm) ? a-x : b-x);
      u  = (abs(d)>=tol1) ? x+d : x+sign(tol1,d);
      fu = f(u);
      if(fu<=fx) {
	if(u>=x) a=x;
	else     b=x;
#define SHFT(a,b,c,d) (a)=(b); (b)=(c); (c)=(d);
	SHFT(v,w,x,u)
	SHFT(fv,fw,fx,fu)
#undef SHFT
      } else {
	if(u<x) a=u;
	else    b=u;
	if(fu<=fw || w==x) {
	  v = w;
	  w = u;
	  fv= fw;
	  fw= fu;
	} else if (fu<=fv || v==x || v==w) {
	  v = u;
	  fv= fu;
	}
      }
    }
    WDutils_Error("in brent(): exceeding iterations");
    xmin   =x;
    return fx;
  }
  // ///////////////////////////////////////////////////////////////////////////
  //
  /// Burlisch-Stoer integration of 1D real integrals
  //
  /// \return approximated value for integral
  /// \param[in]  func  pointer to function to be integrated
  /// \param[in]  a     lower boundary of integration interval
  /// \param[in]  b     upper boundary of integration interval
  /// \param[in]  eps   desired relative accuracy
  /// \param[out] err   (optional) actual relative error of the return value
  /// \param[in]  abort (optional) abort if exceeding maximum of iterations?
  /// \param[in]  miter (optional) maximum number of iterations
  ///                                                                           
  /// Quadrature program using the Bulirsch sequence and rational extrapolation.
  /// The algorithm is puplished in Bulirsch & Stoer, Num. Math. 9, 271-278
  /// (1967), where a routine in ALGOL is given. This is a straightforward
  /// translation into C++.
  ///
  /// \warning Do not use this routine for integrating low order polynomials (up
  /// to fourth order) or periodic functions with period equal to the interval
  /// of integration or linear combinations of both.
  //
  // ///////////////////////////////////////////////////////////////////////////
  double qbulir(double(*func)(double),
		double  a,
		double  b,
		double  eps,
		double* err  =0,
                bool    abort=true,
		int     miter=25);
  //----------------------------------------------------------------------------
  // Runge-Kutta 4th order integrator for ODEs                                  
  //----------------------------------------------------------------------------
  template<typename scalar_type, typename vector_type>
  inline
  vector_type rk4(vector_type const&y,
		  vector_type const&dy0,
		  scalar_type x,
		  scalar_type h,
		  vector_type(*derivs)(scalar_type, vector_type const&))
  {
    const scalar_type hh=0.5*h, h6=h/6., xh=x+hh;
    vector_type yt,dym,dyt;
    dyt = derivs(xh,y+hh*dy0);
    dym = derivs(xh,y+hh*dyt);
    yt  = y+h*dym;
    dym+= dyt;
    dyt = derivs(x+h,yt);
    return y+h6*(dy0+dyt+dym+dym);
  }
  //----------------------------------------------------------------------------
  template<typename scalar_type, typename vector_type>
  inline
  vector_type rk4(vector_type const&y,
		  scalar_type x,
		  scalar_type h,
		  vector_type(*derivs)(scalar_type, vector_type const&))
  {
    return rk4(y,derivs(x,y),x,h,derivs);
  }
  //----------------------------------------------------------------------------
  // Legendre polynomials and their derivatives                                 
  //----------------------------------------------------------------------------
  template<typename S, int N>
  void LegendrePeven(S*p, double x) 
  {
    // based on a routine from J.J. Binney                                      
    // evaluates even Legendre Polys up to l=2*(N-1) at x                       
    double x2=x*x;
    p[0] = 1.;
    p[1] = 1.5*x2-0.5;
    for(int n=2; n<N; n++) {
      int l  = 2*(n-1);
      int l2 = l+l;
      p[n] = - p[n-2] * l*(l-1)/double((l2+1)*(l2-1))
	     + p[n-1] * (x2 - (l2*l+l2-1)/double((l2-1)*(l2+3)) );
      p[n]*= (l2+1)*(l2+3) / double((l+1)*(l+2));
    }
  }
  //----------------------------------------------------------------------------
  template<typename S, int N>
  void dLegendrePeven(S*p, S*d, double x)
    // based on a routine from J.J. Binney                                      
    // evaluates even Legendre Polys and its derivs up to l=2*(N-1) at x        
  {
    double x2=x*x;
    p[0] = 1.;
    d[0] = 0.;
    p[1] = 1.5*x2-0.5;
    d[1] = 1.5;
    for(int n=2; n<N; ++n) {
      int l  = 2*(n-1);
      int l2 = l+l;
      p[n] = - p[n-2] * l*(l-1)/double((l2+1)*(l2-1))
	     + p[n-1] * (x2 - (l2*l+l2-1)/double((l2-1)*(l2+3)) );
      p[n]*= (l2+1)*(l2+3) / double((l+1)*(l+2));
      d[n] = - d[n-2] * l*(l-1)/double((l2+1)*(l2-1))
	     + d[n-1] * (x2 - (l2*l+l2-1)/double((l2-1)*(l2+3)) )
	     + p[n-1];
      d[n]*= (l2+1)*(l2+3) / double((l+1)*(l+2));
    }
    x2 = x+x;
    for(int n=0; n<N; n++)
      d[n] *= x2;
  }
  //----------------------------------------------------------------------------
  template<typename S, int N> inline
  void LegendrePeven(tupel<N,S>& p, double x) {
    return LegendrePeven<S,N>(p,x);
  }
  //----------------------------------------------------------------------------
  template<typename S, int N> inline
  void dLegendrePeven(tupel<N,S>& p, tupel<N,S>& d, double x) {
    return dLegendrePeven<S,N>(p,d,x);
  }
  //----------------------------------------------------------------------------
  // Gauss-Legendre integration: points & weights                               
  //----------------------------------------------------------------------------
  void GaussLegendre(double*, double*, unsigned);
  // ///////////////////////////////////////////////////////////////////////////
  //
  /// \name eigensystem of symmetric matrix using Jacobi transformation
  /// \note not fully tested
  //@{
  // ---------------------------------------------------------------------------
  /// Eigen values and vectors for symmetric matrix using Jacobi transformation
  /// \param      N (template parameter) size of matrix
  /// \param      X (template parameter) scalar type (either float or double)
  /// \param[in]  M  symmetric matrix
  /// \param[out] V  matrix with Eigenvectors of M
  /// \param[out] D  vector with Eigenvalues of M
  /// \param[out] R  number of rotations required
  ///
  /// The eigen values and vectors of a symmetric matrix are computed using
  /// Jacobi transformations (see NR section 11.1). This is not efficient for
  /// efficient for large N, hence we have coded N as template parameter.
  /// --------------------------------------------------------------------------
  template<int N, typename X>
  void EigenSymJacobi(const X M[N][N], X V[N][N], X D[N], int&R)
  {
    const X   zero    = X(0);
    const X   half    = X(0.5);
    const X   one     = X(1);
    const int MaxIter = sizeof(X)*13;
    // copy M to A, set V to unity, copy diagonal of A to B and D
    X A[N][N], B[N], Z[N];
    for(int ip=0; ip!=N; ++ip) {
      for(int iq=0; iq!=N; ++iq) {
	A[ip][iq] = M[ip][iq];
	V[ip][iq] = zero;
      }
      V[ip][ip] = one;
      B[ip]     = A[ip][ip];
      D[ip]     = A[ip][ip];
      Z[ip]     = zero;
    }
    // perform iteration
    R = 0;
    for(int iter=0; iter!=MaxIter; ++iter) {
      X sm(zero);
      for(int ip=0; ip!=N-1; ++ip)
	for(int iq=ip+1; iq!=N; ++iq)
	  sm += abs(A[ip][iq]);
      if(sm == zero) return;
      X tresh = iter<3? sm/X(5*N*N) : zero;
      for(int ip=0; ip!=N-1; ++ip) {
	for(int iq=ip+1; iq!=N; ++iq) {
	  X g = 100 * abs(A[ip][iq]);
	  if(iter > 3                    &&
	     abs(D[ip])+g == abs(D[ip])  &&
	     abs(D[iq])+g == abs(D[iq])     )
	    A[ip][iq] = zero;
	  else if(abs(A[ip][iq]) > tresh) {
	    X h = D[iq]-D[ip], t;
	    if(abs(h)+g == abs(h))
	      t = A[ip][iq]/h;
	    else {
	      X theta = half*h/A[ip][iq];
	      t = one/(abs(theta)+sqrt(one+theta*theta));
	      if(theta < zero) t = -t;
	    }
	    X c   = one/sqrt(one+t*t);
	    X s   = t*c;
	    X tau = s/(one+c);
	    h     = t*A[ip][iq];
	    Z[ip] -= h;
	    Z[iq] += h;
	    D[ip] -= h;
	    D[iq] += h;
	    A[ip][iq] = zero;
#define Rotate(M,i,j,k,l)			\
  {						\
    g = M[i][j];				\
    h = M[k][l];				\
    M[i][j] = g-s*(h+g*tau);			\
    M[k][l] = h+s*(g-h*tau);			\
  }
	    for(int j=0;    j!=ip; ++j) Rotate(A,j,ip,j,iq);
	    for(int j=ip+1; j!=iq; ++j) Rotate(A,ip,j,j,iq);
	    for(int j=iq+1; j!=N;  ++j) Rotate(A,ip,j,iq,j);
	    for(int j=0;    j!=N;  ++j) Rotate(V,j,ip,j,iq);
#undef Rotate
	    ++R;
	  }
	}
      }
      for(int ip=0; ip!=N; ++ip) {
	B[ip] += Z[ip];
	D[ip]  = B[ip];
	Z[ip]  = zero;
      }
    }
    WDutils_THROW("EigenSymJacobi(): number iteration exceeds %d\n",MaxIter);
  }
  // ---------------------------------------------------------------------------
  /// function template sorting the eigenvalues & vectors by straight insertion
  ///
  /// \param         N (template parameter) size of matrix
  /// \param         X (template parameter) scalar type (float or double)
  /// \param[in,out] V  matrix with Eigenvectors
  /// \param[in,out] D  vector with Eigenvalues
  // ---------------------------------------------------------------------------
  template<int N, typename X>
  void EigenSort(X V[N][N], X D[N])
  {
    for(int k,i=0; i!=N-1; ++i) {
      X p = D[k=i];
      for(int j=i+1; j!=N; ++j)
	if(D[j] >= p) p=D[k=j];
      if(k!=i) {
	D[k] = D[i];
	D[i] = p;
	for(int j=0; j!=N; ++j) {
	  p       = V[j][i];
	  V[j][i] = V[j][k];
	  V[j][k] = p;
	}
      }
    }
  }
  // ---------------------------------------------------------------------------
  /// Sorted eigen values and vectors of symmetric matrix with Jacobi transform
  /// \param      N  (template parameter) size of matrix
  /// \param      X  (template parameter) scalar type (float or double)
  /// \param[in]  M  symmetric matrix
  /// \param[out] V  matrix with Eigenvectors, sorted
  /// \param[out] D  vector with Eigenvalues = columns, sorted
  /// \param[out] R  number of rotations required
  /// This routine simply combines EigenSymJacobi() and EigenSort().
  template<int N, typename X>
  void EigenSymJacobiSorted(const X M[N][N], X V[N][N], X D[N], int&R)
  {
    EigenSymJacobi(M,V,D,R);
    EigenSort(V,D);
  }
  //@}                                                                          
  // ---------------------------------------------------------------------------
  /// replaces matrix M with its transposed
  template<int N, typename X>
  void Transpose(X M[N][N])
  {
    for(int i=1; i!=N; ++i)
      for(int j=0; j!=i; ++j)
	swap(M[i][j], M[j][i]);
  }
  // ---------------------------------------------------------------------------
  /// \name eigensystem of symmetric matrix using Householder transformation
  /// \note not fully tested
  //@{
  // ---------------------------------------------------------------------------
  /// reduce real symmetric matrix to tridiagonal form
  /// \param T (template parameter) prepare for eigenvector extraction?
  /// \param X (template parameter) scalar type (float or double)
  /// \param[in]     N size of matrix
  /// \param[in,out] A on input: real symmetric matrix,\n
  ///                  on putput: input required by \a EigenSystemTridiagonal()
  /// \param[out]    D diagonal elements of tridiagonal form
  /// \param[out]    E off-diagonal elements of tridiagonal form
  ///
  /// Householder reduction of a real symmetric matrix
  /// \a A[0..\a N -1][0..\a N -1]. On output, \a A is replaced by an orthogonal
  /// matrix effecting the transformation. \a D[0..\a N -1] returns the diagonal
  /// elements of the tridiagonal matrix, and \a E[0..\a N -1] the off-diagonal
  /// elements, with \a E[0]=0.
  /// If \a T is set to false, \a A has no sensible meaning on output. See NR
  /// section 11.2 for details.
  /// \warning not tested, may be buggy
  template<bool T, typename X> void HouseholderReduction(int N, X**A, X*D, X*E);
  // ---------------------------------------------------------------------------
  /// compute eigensystem of tridiagonal symmetric matrix
  ///
  /// \param X (template parameter) only X=float and X=double are implemented
  /// \param[in]     N  size of matrix
  /// \param[in,out] D  on input: diagonal elements; on output: eigenvalues
  /// \param[in,out] E  on input: off-diagonal elements; on output: destroyed
  /// \param[in,out] Z  on input: see below; on output: EVs corresponding to D
  ///
  /// QL algorithm with implicit shifts to determine the eigenvalues and eigen-
  /// vectors of a real symmetric tridiagonal matrix, or of a real symmetric
  /// matrix previously reduced by \a HouseholderReduction(). In the first case,
  /// \a Z on input must be the unit matrix. In the second case, \a Z on input
  /// must be the matrix returned by \a HouseholderReduction(). For details, see
  /// NR section 11.3.
  /// \warning not tested, may be buggy
  template<typename X> void EigenSystemTridiagonal(int N, X*D, X*E, X**Z);
  // ---------------------------------------------------------------------------
  /// compute eigenvalues of tridiagonal symmetric matrix
  ///
  /// \param X (template parameter) only X=float and X=double are implemented
  /// \param[in]     N   size of matrix
  /// \param[in,out] D  on input: diagonal elements; on output: eigenvalues
  /// \param[in,out] E  on input: off-diagonal elements; on output: destroyed
  ///
  /// QL algorithm with implicit shifts to determine the eigenvalues of a real
  /// symmetric tridiagonal matrix, or of a real symmetric matrix previously
  /// reduced by \a HouseholderReduction(). For details, see NR section 11.3.
  /// \warning not tested, may be buggy
  template<typename X> void EigenValuesTridiagonal(int N, X*D, X*E);
  // ---------------------------------------------------------------------------
  /// eigensystem of symmetric matrix, replaces original matrix
  ///
  /// \param EIGENVECTORS (template parameter) want eigenvectors?
  /// \param X (template parameter) scalar type (either float or double)
  /// \param[in]     N  size of matrix
  /// \param[in,out] A  on input: symmetric matrix; on output: eigenvectors
  /// \param[out]    D  vector with eigenvalues 
  /// \param[out]    V  (optional) vector with eigenvectors
  ///
  /// This combines HouseholderReduction() and EigenSystemTridiagonal() or
  /// EigenValuesTridiagonal() for \a N known at run time. If \a N is known at
  /// compile time, use \a EigenSymmetricFixed() below.
  /// \warning not tested, may be buggy
  template<bool EIGENVECTORS, typename X>
  void EigenSymmetricReplace(int N, X**A, X*D)
  {
    X*E = WDutils_NEW(X,N);
    HouseholderReduction<EIGENVECTORS>(N,A,D,E);
    if(EIGENVECTORS) EigenSystemTridiagonal(N,D,E,A);
    else             EigenValuesTridiagonal(N,D,E);
    falcON_DEL_A(E);
  }
  // ---------------------------------------------------------------------------
  /// eigensystem of symmetric matrix, keeps original matrix
  ///
  /// \param X (template parameter) scalar type (either float or double)
  /// \param[in]  N  size of matrix
  /// \param[in]  M  symmetric matrix
  /// \param[out] D  vector with eigenvalues
  /// \param[out] V  (optional) vector with eigenvectors
  ///
  /// This combines HouseholderReduction() and EigenSystemTridiagonal() or
  /// EigenValuesTridiagonal() for \a N known at run time. If \a N is known at
  /// compile time, use \a EigenSymmetricFixed() below.
  /// \warning not tested, may be buggy
  template<typename X>
  void EigenSymmetricKeep(int N, const X**M, X*D, X**V=0)
  {
    X*E = WDutils_NEW(X,N);
    if(V) {
      for(int i=0; i!=N; ++i)
	for(int j=0; j!=N; ++j)
	  V[i][j] = M[i][j];
      HouseholderReduction<1>(N,V,D,E);
      EigenSystemTridiagonal(N,D,E,V);
    } else {
      V = WDutils_NEW(X*,N);
      X* m = WDutils_NEW(X,N*N);
      for(int i=0; i!=N; ++i, m+=N) {
	V[i] = m;
	for(int j=0; j!=N; ++j)
	  V[i][j] = M[i][j];
      }
      HouseholderReduction<0>(N,V,D,E);
      EigenValuesTridiagonal(N,D,E);
      falcON_DEL_A(V[0]);
      falcON_DEL_A(V);
    }
    falcON_DEL_A(E);
  }
  // ---------------------------------------------------------------------------
  /// eigensystem of symmetric matrix, keeps original matrix, N template param
  ///
  /// \param      N  (template parameter) size of matrix
  /// \param      X  (template parameter) scalar type (float or double)
  /// \param[in]  M  matrix
  /// \param[out] D  vector with Eigenvalues
  /// \param[out] V  (optional) matrix with Eigenvectors
  ///
  /// This combines HouseholderReduction() and EigenSystemTridiagonal() for
  /// fixed \a N (known as template parameter at compile time).
  /// \warning not tested, may be buggy
  template<int N, typename X>
  void EigenSymmetricFixed(const X M[N][N], X D[N], X**V=0)
  {
    X E[N];
    if(V) {
      for(int i=0; i!=N; ++i)
	for(int j=0; j!=N; ++j)
	  V[i][j] = M[i][j];
      HouseholderReduction<1>(N,V,D,E);
      EigenSystemTridiagonal (N,D,E,V);

    } else {
      X Z[N][N];
      for(int i=0; i!=N; ++i)
	for(int j=0; j!=N; ++j)
	  Z[i][j] = M[i][j];
      HouseholderReduction<0>(N,Z,D,E);
      EigenValuesTridiagonal (N,D,E);
    }
  }
  //@}
  //////////////////////////////////////////////////////////////////////////////
} // namespace WDutils {
////////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_numerics_h
