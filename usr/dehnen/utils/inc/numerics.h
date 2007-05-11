// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    utils/inc/numerics.h                                               
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    1994-2007                                                          
///                                                                             
/// \todo    finish doxygen documentation                                       
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 1994-2007  Walter Dehnen                                       
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
  /// \param xarr (input) array of T, must be ordered (ascending or descending) 
  /// \param n    (input) size of array xarr                                    
  /// \param x    (input) value to find position for                            
  /// \param j    (input) initial guess for position                            
  /// \return     position jlo such that xarr[jlo] <= x < xarr[jlo+1]           
  ///                                                                           
  /// The ordered table xarr is hunted for jlo such that                        
  ///    xarr[jlo] <= x < xarr[jlo+1].                                          
  /// For an ascendingly ordered array, we return -1 if x < xarr[0], n-1 if     
  /// x == xarr[n-1], and n if x > xarr[n--1].                                  
  //                                                                            
  //----------------------------------------------------------------------------
  template<typename T>
  int hunt(const T*xarr, const int n, const T x, const int j) {
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
  /// \param T (template parameter) type with < and > operators, usually a      
  ///          scalar                                                           
  /// \param k    (in/out) jlo such that xarr[jlo] <= x < xarr[jlo+1]           
  /// \param x    (input) array of T, must be ordered (ascending or descending) 
  /// \param n    (input) size of array xarr                                    
  /// \param xi   (input) value to find position for                            
  /// \param j    (input) initial guess for position                            
  /// \return     position jlo such that xarr[jlo] <= x < xarr[jlo+1]           
  ///                                                                           
  /// If the original value for k already gives the position, we return.        
  /// Otherwise, we guess k from linear interpolation and then invoke hunt().   
  /// If x is not in range, we throw an error.                                  
  //                                                                            
  template<typename T>
  inline void find(int& k, const int n, T const *x, const T xi)
  {
    if(k<0 || k>=n-1 || x[k]>xi || x[k+1]<xi) {
      k = int( (xi-x[0]) / (x[n-1]-x[0]) * (n-1) );
      k = hunt(x,n,xi,k);
      if(k<0 || k>=n) 
	error("[%s.%d]: in %s: x=%f out of range [%f,%f]\n",
	      __FILE__,__LINE__,"find()",xi,x[0],x[n-1]);
    }
  }
  //@}
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// \name sorting and related                                                 
  //{@                                                                          
  //----------------------------------------------------------------------------
  /// the numbers 0 to n-1 are ordered in ascending order of A[i]
  /// using the heap-sort algorithm; based on a routine given in NR
  /// \param sortable type of values to be sorted
  /// \param sortit class with operator[](int) returning sortable
  /// \param A (input) array or values to be sorted
  /// \param n (input) number of elements
  /// \param indx (output) index table so that A[indx[i]] is sorted
  template<typename sortable, class sortit>
  void HeapIndex(const sortit&A, const int n, int *indx)
  {
    if(n<=0) return;
    if(n==1) { indx[0]=0; return; }
    register int      l,j,ir,indxt,i;
    register sortable q;
    for(j=0; j!=n; ++j) indx[j] = j;
    l = n>>1;
    ir= n-1;
    for(;;) {
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
      i = l;
      j = (l<<1) + 1;
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
  /// using the heap-sort algorithm; based on a routine given in NR
  /// \param sortable type of values to be sorted
  /// \param A (input) array or values to be sorted
  /// \param n (input) number of elements
  /// \param indx (output) index table so that A[indx[i]] is sorted
  template<typename sortable>
  void HeapIndex(const sortable*A, const int n, int *indx)
    // based on a routine given in NR                                           
    // the numbers 0 to n-1 are ordered in ascending order of A[i]              
  {
    if(n<=0) return;
    if(n==1) { indx[0]=0; return; }
    register int      l,j,ir,indxt,i;
    register sortable q;
    for(j=0; j!=n; ++j) indx[j] = j;
    l = n>>1;
    ir= n-1;
    for(;;) {
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
      i = l;
      j = (l<<1) + 1;
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
  /// using the heap-sort algorithm; based on a routine given in NR
  /// \param sortable type of values to be sorted
  /// \param A (input) Array or values to be sorted
  /// \param indx (output) index table so that A[indx[i]] is sorted
  template<typename sortable>
  void HeapIndex(Array<sortable,1>const&A, Array<int,1>const&I) {
    if(A.size() != I.size()) error("size mismatch in HeapIndex()\n");
    HeapIndex(A.array(),A.size(),I.array());
  }
  //----------------------------------------------------------------------------
  /// obtain approximate function values for given ranks.
  ///
  /// Given a table of approximately evenly distributed ranks R[i], make an 
  /// approximate table of values Y[i] so that rank(Y[i]) ~= R[i].
  /// For size of table of ranks m << n, this is considerably faster than  
  /// using a full sort.
  ///
  /// \param sortable type of values to be sorted
  /// \param sortit class with operator[](int) returning sortable
  /// \param A (input) objects to be sorted
  /// \param n (input) number of objects
  /// \param R (input) table of ranks
  /// \param Y (output) table of approximate values at given ranks
  /// \param m (input) size of tables R and Y, usually m << n
  /// \param indx (output) index table so that A[indx[i]] is sorted
  /// \param Amin (optional input) minimum value
  /// \param Amax (optional input) maximum value
  template<typename sortable, typename sortit> inline
  void set_value_table(sortit         const&A,
		       int            const&n,
		       const sortable*const&R,
		       sortable      *const&Y,
		       int            const&m,
		       const sortable* Amin= 0,
		       const sortable* Amax= 0);
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
    /// \param F (input) array with elements
    /// \param N (input) number of elements
    /// \param W (optional input) array of weights
    FindPercentile(const scalar*F, int N, const scalar*W=0) : DATA(0) {
      setup(F,N,W);
    }
    /// ctor: setup sort tree
    /// \param F (input) Array with elements
    FindPercentile(Array<scalar,1>const&F) : DATA(0) {
      setup(F.array(),F.size(),0);
    }
    /// ctor: setup sort tree
    /// \param F (input) Array<> with elements
    /// \param W (input) Array<> of weights
    FindPercentile(Array<scalar,1>const&F,
		   Array<scalar,1>const&W) : DATA(0) {
      if(F.size() != W.size()) error("size mismatch in FindPercentile\n");
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
  //@}
  //----------------------------------------------------------------------------
  // polynomial interpolation on grids                                          
  //----------------------------------------------------------------------------
  template<typename scalar_type>
  inline int find_for_polev(int& j, int n, int m,
			    const scalar_type *x, scalar_type xi)
  {
    int M=m;
    j = int( (xi-x[0]) / (x[n-1]-x[0]) * (n-1) );
    j = hunt(x,n,xi,j) - (m+1)/2 + 1;
    if(j>=0 && j<n && x[j]==xi)	 M = 1; 	    // no interpolation required
    else if(j<0)		 j = 0;
    else if(j>n-M)		 j = n-M;
    return M;
  }
  //----------------------------------------------------------------------------
  template<int m, typename scalar_type>
  inline int find_for_polev_T(int& j, int n,
			      const scalar_type *x, scalar_type xi)
  {
    int M=m;
    j = int( (xi-x[0]) / (x[n-1]-x[0]) * (n-1) );
    j = hunt(x,n,xi,j) - (m+1)/2 + 1;
    if(j>=0 && j<n && x[j]==xi)	 M = 1; 	    // no interpolation required
    else if(j<0)		 j = 0;
    else if(j>n-M)		 j = n-M;
    return M;
  }
  //----------------------------------------------------------------------------
  template<typename scalar_type, typename num_type>
  inline num_type polint(const scalar_type *xa,
			 const num_type *ya,
			 int n,
			 scalar_type x)
  {
    // polynom interpolation using n values
    num_type y,*P=WDutils_NEW(num_type,n);
    for(int i=0;i!=n;++i) P[i]=ya[i];
    for(int m=1;m!=n;++m)
      for(int i=0;i<n-m;++i) {
	if(xa[i]==xa[i+m]) WDutils_ErrorF("x's not distinct","polev()");
	P[i]= ( (x-xa[i+m])*P[i] + (xa[i]-x)*P[i+1] ) / (xa[i] - xa[i+m]);
      }
    y = P[0];
    WDutils_DEL_A(P);
    return y;
  }
  //----------------------------------------------------------------------------
  template<int n, typename scalar_type, typename num_type>
  inline num_type polint_T(const scalar_type*xa,
			   const num_type   *ya,
			   scalar_type x) 
  {
    // polynom interpolation using n values
    num_type y, P[n];
    for(int i=0;i!=n;++i) P[i]=ya[i];
    for(int m=1;m!=n;++m)
      for(int i=0;i<n-m;++i) {
	if(xa[i]==xa[i+m]) WDutils_ErrorF("x's not distinct","polev()");
	P[i]= ( (x-xa[i+m])*P[i] + (xa[i]-x)*P[i+1] ) / (xa[i] - xa[i+m]);
      }
    y = P[0];    
    return y;
  }
  //----------------------------------------------------------------------------
  // polynomial interpolation using m of n values                               
  template<typename scalar_type, typename num_type>
  inline num_type polev(scalar_type x,
			const scalar_type*xarr,
			const num_type*yarr,
			int n, int m)
  {
    // given the arrays xarr, yarr, polev returns y(x) using m of n values.
    int j, M;
    switch(m) {
    case 2 : 
      M = find_for_polev_T<2>(j,n,xarr,x);
      return polint_T<2>(xarr+j, yarr+j, x);
    case 3 : 
      M = find_for_polev_T<3>(j,n,xarr,x);
      return polint_T<3>(xarr+j, yarr+j, x);
    case 4 : 
      M = find_for_polev_T<4>(j,n,xarr,x);
      return polint_T<4>(xarr+j, yarr+j, x);
    case 5 : 
      M = find_for_polev_T<5>(j,n,xarr,x);
      return polint_T<5>(xarr+j, yarr+j, x);
    case 6 : 
      M = find_for_polev_T<6>(j,n,xarr,x);
      return polint_T<6>(xarr+j, yarr+j, x);
    default: 
      M = find_for_polev(j,n,m,xarr,x);
      return polint(xarr+j, yarr+j, M, x);
    }
  }
  //----------------------------------------------------------------------------
  // polynomial interpolation using 4 of n values                               
  template<typename scalar_type, typename num_type>
  inline num_type polev(scalar_type x,
			const scalar_type*xarr,
			const num_type*yarr,
			int n)
  {
    // given the arrays xarr, yarr, polev returns y(x) using 4 of n values.
    int j, M=find_for_polev_T<4>(j,n,xarr,x);
    return polint_T<4>(xarr+j, yarr+j, x);
  }
  //----------------------------------------------------------------------------
  // like polev, but taking Array<T>
  template<typename scalar_type, typename num_type>
  inline num_type polev(scalar_type x,
			const Array<scalar_type,1>&xarr,
			const Array<num_type,1>&yarr)
  {
    if(xarr.size() != yarr.size()) error("size mismatch in polev()\n");
    return polev(x,xarr.array(),yarr.array(),xarr.size());
  }
  //----------------------------------------------------------------------------
  template<typename scalar_type, typename num_type>
  inline num_type ipolev(scalar_type x,
			 const scalar_type*xarr,
			 const num_type   *yarr,
			 int n, int m)
  {
    // like polev, but no extrapolation. gives boundary values instead          
    if(x < xarr[0]   && xarr[0]   < xarr[n-1]) return yarr[0];
    if(x > xarr[0]   && xarr[0]   > xarr[n-1]) return yarr[0];
    if(x > xarr[n-1] && xarr[n-1] > xarr[0]  ) return yarr[n-1];
    if(x < xarr[n-1] && xarr[n-1] < xarr[0]  ) return yarr[n-1];
    return polev<scalar_type,num_type>(x,xarr,yarr,n,m);
  }
  //----------------------------------------------------------------------------
  template<typename scalar_type, typename num_type>
  inline num_type ipolev(scalar_type x,
			 const scalar_type*xarr,
			 const num_type*yarr,
			 int n)
  {
    // like polev, but no extrapolation. gives boundary values instead          
    if(x < xarr[0]   && xarr[0]   < xarr[n-1]) return yarr[0];
    if(x > xarr[0]   && xarr[0]   > xarr[n-1]) return yarr[0];
    if(x > xarr[n-1] && xarr[n-1] > xarr[0]  ) return yarr[n-1];
    if(x < xarr[n-1] && xarr[n-1] < xarr[0]  ) return yarr[n-1];
    return polev<scalar_type,num_type>(x,xarr,yarr,n);
  }
  //----------------------------------------------------------------------------
  // root finding                                                               
  //----------------------------------------------------------------------------
  template <typename scalar_type>
  scalar_type rtsafe(void(*func)(scalar_type,scalar_type&,scalar_type&),
		     scalar_type x1,
		     scalar_type x2,
		     scalar_type xacc)
  {
    const int            maxit=100;
    register scalar_type xh,xl,dx,dxo,f,df,fh,fl,rts,swap,temp;
    func(x1,fl,df);
    func(x2,fh,df);
    if(fl*fh >= 0.) WDutils_ErrorF("root must be bracketed","rtsafe()");
    if(fl<0.) {
      xl  = x1;
      xh  = x2;
    } else {
      xh  = x1;
      xl  = x2;
      swap= fl;
      fl  = fh;
      fh  = swap;
    }
    rts = 0.5*(x1+x2);
    dxo = abs(x2-x1);
    dx  = dxo;
    func(rts,f,df);
    for (register int j=0; j!=maxit; ++j) {
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
    WDutils_WarningF("max number of iterations exceeded","rtsafe()");
    return rts;
  }
  //----------------------------------------------------------------------------
  // bracketing a minimum                                                       
  //----------------------------------------------------------------------------
  template<typename scalar_type>
  void minbrak(scalar_type& ax,                    // I/O ax                    
	       scalar_type& bx,                    // I/O bx                    
	       scalar_type& cx,                    // O:  cx                    
	       scalar_type(*f)(const scalar_type)) // I: f(x)                   
  {
    const    scalar_type GLIMIT=100.0, GOLD=1.618034, TINY=1.e-20;
    register scalar_type ulim,u,r,q,fu,fa=f(ax),fb=f(bx),fc;
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
	const scalar_type ax,                      // I: ax   these bracket the 
	const scalar_type bx,                      // I: bx   minimum, e.g. out-
	const scalar_type cx,                      // I: cx   put of minbrak()  
	scalar_type(*f)(const scalar_type),        // I: function f(x)          
	const scalar_type tol,                     // I: accuracy wanted        
	scalar_type& xmin)                         // O: x_min                  
  {
    const    int         itmax=100;
    const    scalar_type cgold=0.381966, zeps=1.e-10;
    register int         iter;
    register scalar_type a,b,d=0.,e=0.,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,
                         u,v,w,x,xm;
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
    WDutils_ErrorF("exceeding iterations","brent()");
    xmin   =x;
    return fx;
  }
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// Burlisch-Stoer integration of 1D real integrals                           
  //                                                                            
  /// \return approximated value for integral                                   
  /// \param  func  (input) pointer to function to be integrated                
  /// \param  a     (input) lower boundary of integration interval              
  /// \param  b     (input) upper boundary of integration interval              
  /// \param  eps   (input) desired relative accuracy                           
  /// \param  err   (output, optional) actual relative error of the return value
  /// \param  abort (input, optional) abort if exceeding maximum of iterations? 
  /// \param  miter (input, optional) maximum number of iterations              
  ///                                                                           
  /// Quadrature program using the Bulirsch sequence and rational extrapolation.
  /// The algorithm is puplished in Bulirsch & Stoer, Num. Math. 9, 271-278     
  /// (1967), where a routine in ALGOL is given. This is a straightforward      
  /// translation into C++.                                                     
  ///                                                                           
  /// \note CAUTION:                                                            
  /// Do not use this routine for integrating low order polynomials (up to      
  /// fourth order) or periodic functions with period equal to the interval of  
  /// integration or linear combinations of both.                               
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
  vector_type rk4(const vector_type y,
		  const vector_type dy0,
		  const scalar_type x,
		  const scalar_type h,
		  vector_type(*derivs)(scalar_type const&, vector_type const&))
  {
    const    scalar_type hh=0.5*h, h6=h/6., xh=x+hh;
    register vector_type yt,dym,dyt;
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
  vector_type rk4(const vector_type y,
		  const scalar_type x,
		  const scalar_type h,
		  vector_type(*derivs)(scalar_type const&, vector_type const&))
  {
    return rk4(y,derivs(x,y),x,h,derivs);
  }
  //----------------------------------------------------------------------------
  // Legendre polynomials and their derivatives                                 
  //----------------------------------------------------------------------------
  template<typename S, int N>
  void LegendrePeven(S*p, const double x) 
  {
    // based on a routine from J.J. Binney                                      
    // evaluates even Legendre Polys up to l=2*(N-1) at x                       
    register int    n,l,l2;
    register double x2=x*x;
    p[0] = 1.;
    p[1] = 1.5*x2-0.5;
    for(n=2; n<N; n++) {
      l    = 2*(n-1);
      l2   = l+l;
      p[n] = - p[n-2] * l*(l-1)/double((l2+1)*(l2-1))
	     + p[n-1] * (x2 - (l2*l+l2-1)/double((l2-1)*(l2+3)) );
      p[n]*= (l2+1)*(l2+3) / double((l+1)*(l+2));
    }
  }
  //----------------------------------------------------------------------------
  template<typename S, int N>
  void dLegendrePeven(S*p, S*d, const double x)
    // based on a routine from J.J. Binney                                      
    // evaluates even Legendre Polys and its derivs up to l=2*(N-1) at x        
  {
    register int    n,l,l2;
    register double x2=x*x;
    p[0] = 1.;
    d[0] = 0.;
    p[1] = 1.5*x2-0.5;
    d[1] = 1.5;
    for(n=2; n<N; n++) {
      l    = 2*(n-1);
      l2   = l+l;
      p[n] = - p[n-2] * l*(l-1)/double((l2+1)*(l2-1))
	     + p[n-1] * (x2 - (l2*l+l2-1)/double((l2-1)*(l2+3)) );
      p[n]*= (l2+1)*(l2+3) / double((l+1)*(l+2));
      d[n] = - d[n-2] * l*(l-1)/double((l2+1)*(l2-1))
	     + d[n-1] * (x2 - (l2*l+l2-1)/double((l2-1)*(l2+3)) )
	     + p[n-1];
      d[n]*= (l2+1)*(l2+3) / double((l+1)*(l+2));
    }
    x2 = x+x;
    for(n=0; n<N; n++)
      d[n] *= x2;
  }
  //----------------------------------------------------------------------------
  template<typename S, int N> inline
  void LegendrePeven(tupel<N,S>& p, const double x) {
    return LegendrePeven<S,N>(p,x);
  }
  //----------------------------------------------------------------------------
  template<typename S, int N> inline
  void dLegendrePeven(tupel<N,S>& p, tupel<N,S>& d, const double x) {
    return dLegendrePeven<S,N>(p,d,x);
  }
  //----------------------------------------------------------------------------
  // Gauss-Legendre integration: points & weights                               
  //----------------------------------------------------------------------------
  void GaussLegendre(double*, double*, const unsigned);
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// \name eigensystem of symmetric matrix                                     
  //@{                                                                          
  //                                                                            
  /// \name eigensystem of symmetric matrix using Jacobi transformation         
  //@{                                                                          
  //                                                                            
  // ---------------------------------------------------------------------------
  /// Eigen values and vector for symmetric matrix using Jacobi transformation  
  ///                                                                           
  /// \param N (template parameter) size of matrix                              
  /// \param X (template parameter) scalar type (either float or double)        
  /// \param M (input)  matrix                                                  
  /// \param V (output) matrix with Eigenvectors                                
  /// \param D (output) vector with Eigenvalues                                 
  /// \param R (output) number of rotations required                            
  ///                                                                           
  /// The eigen values and vectors of a symmetric matrix are computed using     
  /// Jacobi transformations (see NR section 11.1). This is not efficient for   
  /// efficient for large N, hence we have coded N as template parameter.       
  // ---------------------------------------------------------------------------
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
    error("EigenSymJacobi(): number iteration exceeds %d\n",MaxIter);
  }
  // ---------------------------------------------------------------------------
  /// function template sorting the eigenvalues & vectors                       
  ///                                                                           
  /// \param N (template parameter) size of matrix                              
  /// \param X (template parameter) scalar type (either float or double)        
  /// \param V (input/output)  matrix with Eigenvectors                         
  /// \param D (input/output)  vector with Eigenvalues                          
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
  /// Sorted eigen values and vector of symmetric matrix with Jacobi transform  
  ///                                                                           
  /// \param N (template parameter) size of matrix                              
  /// \param X (template parameter) scalar type (either float or double)        
  /// \param M (input)  matrix                                                  
  /// \param V (output) matrix with Eigenvectors                                
  /// \param D (output) vector with Eigenvalues                                 
  /// \param R (output) number of rotations required                            
  ///                                                                           
  /// This combines EigenSymJacobi() and EigenSort() via                        
  /// \code                                                                     
  ///   EigenSymJacobi(M,V,D,R);                                                
  ///   EigenSort(V,D);                                                         
  /// \endcode                                                                  
  // ---------------------------------------------------------------------------
  template<int N, typename X>
  void EigenSymJacobiSorted(const X M[N][N], X V[N][N], X D[N], int&R)
  {
    EigenSymJacobi(M,V,D,R);
    EigenSort(V,D);
  }
  //@}                                                                          
  //                                                                            
  /// \name eigensystem of symmetric matrix using Householder transformation    
  //@{                                                                          
  //                                                                            
  // ---------------------------------------------------------------------------
  /// reduce real symmetric matrix to tridiagonal form                          
  ///                                                                           
  /// \param T (template parameter) prepare for eigenvector extraction?         
  /// \param X (template parameter) only \a X=float and \a X=double are         
  ///          implemented                                                      
  /// \param N (input)  size of matrix                                          
  /// \param A on input: real symmetric matrix,                                 
  ///          on putput: input required by \a EigenSystemTridiagonal()         
  /// \param D (output) diagonal elements of tridiagonal form                   
  /// \param E (output) off-diagonal elements of tridiagonal form               
  ///                                                                           
  /// Householder reduction of a real symmetric matrix                          
  /// \a A[0..\a N -1][0..\a N -1]. On output, \a A is replaced by an orthogonal
  /// matrix effecting the transformation. \a D[0..\a N -1] returns the diagonal
  /// elements of the tridiagonal matrix, and \a E[0..\a N -1] the off-diagonal 
  /// elements, with \a E[0]=0.                                                 
  /// If \a T is set to false, \a A has no sensible meaning on output. See NR   
  /// section 11.2 for details.                                                 
  // ---------------------------------------------------------------------------
  template<bool T, typename X> void HouseholderReduction(int N, X**A, X*D, X*E);
  // ---------------------------------------------------------------------------
  /// compute eigensystem of tridiagonal symmetric matrix                       
  ///                                                                           
  /// \param X (template parameter) only X=float and X=double are implemented   
  /// \param N (input)  size of matrix                                          
  /// \param D on input: diagonal elements; on output: eigenvalues              
  /// \param E on input: off-diagonal elements; on output: destroyed            
  /// \param Z on input: see below; on output: eigenvectors corresponding to D  
  ///                                                                           
  /// QL algorithm with implicit shifts to determine the eigenvalues and eigen- 
  /// vectors of a real symmetric tridiagonal matrix, or of a real symmetric    
  /// matrix previously reduced by \a HouseholderReduction(). In the first case,
  /// \a Z on input must be the unit matrix. In the second case, \a Z on input  
  /// must be the matrix returned by \a HouseholderReduction(). For details,    
  /// see NR section 11.3.                                                      
  // ---------------------------------------------------------------------------
  template<typename X> void EigenSystemTridiagonal(int N, X*D, X*E, X**Z);
  // ---------------------------------------------------------------------------
  /// compute eigenvalues of tridiagonal symmetric matrix                       
  ///                                                                           
  /// \param X (template parameter) only X=float and X=double are implemented   
  /// \param N (input) size of matrix                                           
  /// \param D on input: diagonal elements; on output: eigenvalues              
  /// \param E on input: off-diagonal elements; on output: destroyed            
  ///                                                                           
  /// QL algorithm with implicit shifts to determine the eigenvalues of a real  
  /// symmetric tridiagonal matrix, or of a real symmetric matrix previously    
  /// reduced by \a HouseholderReduction(). For details, see NR section 11.3.   
  // ---------------------------------------------------------------------------
  template<typename X> void EigenValuesTridiagonal(int N, X*D, X*E);
  // ---------------------------------------------------------------------------
  /// eigensystem of symmetric matrix, replaces original matrix                 
  ///                                                                           
  /// \param EIGENVECTORS (template parameter) want eigenvectors?               
  /// \param X (template parameter) scalar type (either float or double)        
  /// \param N (input)  size of matrix                                          
  /// \param A on input: symmetric matrix; on output: eigenvectors              
  /// \param D (output) vector with eigenvalues                                 
  /// \param V (output, optional) vector with eigenvectors                      
  ///                                                                           
  /// This combines HouseholderReduction() and EigenSystemTridiagonal() or      
  /// EigenValuesTridiagonal() for \a N known at run time. If \a N is known at  
  /// compile time, use \a EigenSymmetricFixed() below.                         
  // ---------------------------------------------------------------------------
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
  /// \param N (input)  size of matrix                                          
  /// \param M (input)  symmetric matrix                                        
  /// \param D (output) vector with eigenvalues                                 
  /// \param V (output, optional) vector with eigenvectors                      
  ///                                                                           
  /// This combines HouseholderReduction() and EigenSystemTridiagonal() or      
  /// EigenValuesTridiagonal() for \a N known at run time. If \a N is known at  
  /// compile time, use \a EigenSymmetricFixed() below.                         
  // ---------------------------------------------------------------------------
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
  /// \param N (template parameter) size of matrix                              
  /// \param X (template parameter) scalar type (either float or double)        
  /// \param M (input)  matrix                                                  
  /// \param D (output) vector with Eigenvalues                                 
  /// \param V (output, optional) matrix with Eigenvectors                      
  ///                                                                           
  /// This combines HouseholderReduction() and EigenSystemTridiagonal() for     
  /// fixed \a N (known as template parameter at compile time).                 
  // ---------------------------------------------------------------------------
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
//                                                                            //
// inline functions and implementation details                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace {
  using namespace WDutils;
#define PLEAF static_cast<leaf*>
#define PRNGE static_cast<range*>
  template<typename scalar> class sorttree {
  public:
    //--------------------------------------------------------------------------
    // a node just has a pointer to its lower sub-node                          
    //--------------------------------------------------------------------------
    struct node { node* LO;
      node() : LO(0) {} };
    //--------------------------------------------------------------------------
    // a leaf is a node holding a single scalar X                               
    //--------------------------------------------------------------------------
    struct leaf  : public node { scalar  X; };
    //--------------------------------------------------------------------------
    // a range is a node with an additional pointer to its upper sub-node and   
    // a center C and size S. it contains leafs with X in [C-S, C+S[.           
    //--------------------------------------------------------------------------
    // a twig range has UP=0 and holds in LO a linked list of its leaf.         
    //    NOTE that the leafs may be deleted using delete_leafs(). After call   
    //         to this method, the pointer LO becomes meaningless.              
    // a branch range has UP and LO pointing to its sub-ranges                  
    //--------------------------------------------------------------------------
    struct range : public node { scalar C,S; int N; range *UP;
      scalar lo () const { return C-S; }
      scalar up () const { return C+S; }
      scalar dia() const { return S+S; }
    };
  private:
    typedef block_alloc<range> ranger;             // allocator for new ranges  
    //--------------------------------------------------------------------------
    const int Ncrit;                               // split ranges with N>Ncrit 
    ranger    RG;                                  // allocator for new ranges  
    range    *ROOT;                                // the root range            
    leaf     *L0;                                  // the arrays of leafs       
    //--------------------------------------------------------------------------
    int           A,I,NI;                          // for making value table    
    const scalar* RANK;                            // table: ranks              
    scalar      * VALUE;                           // table: values (to be made)
    //--------------------------------------------------------------------------
    range* split(range* const&RNGE) {
      // split the twig range RNGE into two sub-ranges                          
      register scalar H=scalar(0.5)*RNGE->S;
      register range* LOWR = RG.new_element();
      LOWR->C = RNGE->C - H;
      LOWR->S = H;
      LOWR->N = 0;
      LOWR->UP= 0;
      register range* UPPR = RG.new_element();
      UPPR->C = RNGE->C + H;
      UPPR->S = H;
      UPPR->N = 0;
      UPPR->UP= 0;
      register leaf *LEAF,*NEXT;
      for(LEAF=PLEAF(RNGE->LO); LEAF; LEAF=NEXT) {
	NEXT = PLEAF(LEAF->LO);
	if(LEAF->X < RNGE->C) {
	  LEAF->LO = LOWR->LO;
	  LOWR->LO= LEAF;
	  LOWR->N++;
	} else {
	  LEAF->LO = UPPR->LO;
	  UPPR->LO= LEAF;
	  UPPR->N++;
	}
      }
      RNGE->LO = LOWR;
      RNGE->UP = UPPR;
      if(LOWR->N == 0) return UPPR;
      if(UPPR->N == 0) return LOWR;
      return 0;
    }
    //--------------------------------------------------------------------------
    void addleaf(register range* R, leaf* const&L) {
      // add single leaf L to range R (usually the root range).                 
      // if the twig range into which L falls overspills, it is splitted        
      for(;;) {
	if(R->UP) {
	  R->N++;
	  R = (L->X < R->C) ? PRNGE(R->LO) : R->UP ;
	} else {
	  L->LO = R->LO;
	  R->LO = L;
	  R->N++;
	  if(R->N > Ncrit) while( R = split(R) );
	  return;
	}
      }
    }
    //--------------------------------------------------------------------------
    void treewalk(range* const&R) {
      // create a value table VALUE, given the rank table RANK using linear     
      // interpolation of ranks between the values of the ranges borders        
      if(R->UP) {
	treewalk(PRNGE(R->LO));
	treewalk(R->UP);
      } else if(R->N) {
	register int nA = A+R->N;
	if(RANK[I] >= A && RANK[I] <= nA && I < NI) {
	  VALUE[I]  = R->lo() + (RANK[I]-A)*R->dia()/scalar(R->N);
	  I++;
	}
	A = nA;
      }
    }
    //--------------------------------------------------------------------------
    void merge(range* our, const range* their) {
      // merge two ranges which must have identical C,S                         
      // NOTE see comments with merge() below                                   
      if(their->UP)                                // 1   theirs is branch      
	if(our->UP) {                              // 1.1 ours  is branch       
	  merge(our->UP,their->UP);                //    merge upper sub-ranges 
	  merge(PRNGE(our->LO),PRNGE(their->LO));  //    merge lower sub-ranges 
	} else {                                   // 1.2   ours  is twig       
	  register leaf *LEAF=PLEAF(our->LO), *NEXT;
	  our->LO = their->LO;                     //    make our range         
	  our->UP = their->UP;                     //    a copy of their range  
	  our->N  = their->N;                      //    also copy number       
	  for(; LEAF; LEAF=NEXT) {                 //    loop our leafs         
	    NEXT = PLEAF(LEAF->LO);                //      get next in list     
	    addleaf(our,LEAF);                     //      add them to range    
	  }
	}
      else {                                       // 2   theirs  is twig       
	register leaf *LEAF=PLEAF(their->LO), *NEXT;
	for(; LEAF; LEAF=NEXT) {                   //    loop their leafs       
	  NEXT = PLEAF(LEAF->LO);                  //      get next in list     
	  addleaf(our,LEAF);                       //      add them to ours     
	}
      }
    }
    //--------------------------------------------------------------------------
    // construction: build a tree with twigs of Ncrit leafs                     
    //--------------------------------------------------------------------------
  public:
  template<typename sortit>
    sorttree(sortit const&A,
	     int    const&n,
	     int    const&nc,
	     scalar const&xmin,
	     scalar const&xmax) :
      Ncrit ( max(1,nc) ),
      RG    ( max(10,n/Ncrit) ),
      ROOT  ( RG.new_element() ),
      L0    ( WDutils_NEW(leaf,n) )
    {
      ROOT->C = scalar(0.5)*( ceil(xmax) + floor(xmin) );
      ROOT->S = scalar(0.5)*( ceil(xmax) - floor(xmin) );
      ROOT->N = 0;
      ROOT->UP= 0;
      for(register int i=0; i!=n; ++i) {
	L0[i].X = A[i];
	addleaf(ROOT,L0+i);
      }
    }
    //--------------------------------------------------------------------------
    // delete leafs                                                             
    //--------------------------------------------------------------------------
    void delete_leafs()
    {
      if(L0) { WDutils_DEL_A(L0); L0=0; }
    }
    //--------------------------------------------------------------------------
    // destruction: delete leafs and implicitly call ~ranger                    
    //--------------------------------------------------------------------------
    ~sorttree()
    {
      if(L0) WDutils_DEL_A(L0);
    }
    //--------------------------------------------------------------------------
    // give const access to our root                                            
    //--------------------------------------------------------------------------
    const range* const&root() { return ROOT; }
    //--------------------------------------------------------------------------
    // make a (approximate) value table given a rank table                      
    //--------------------------------------------------------------------------
    void make_value_table(const scalar *const&rank,
			  scalar *      const&value,
			  int           const&n) {
      A    = 0;
      I    = 0;
      NI   = n;
      RANK = rank;
      VALUE= value;
      treewalk(ROOT);
    }
    //--------------------------------------------------------------------------
    // merge with another sorttree                                              
    // NOTE The merged tree will mix leafs and ranges allocated in this and     
    //      that. Thus, any operation on the merged tree will cause an run-time 
    //      error (Segmentation fault), if that sorttree has been deleted after 
    //      the merger.                                                         
    // NOTE While it is, of course, possible to design code without this feature
    //      it would be less efficient.                                         
    //--------------------------------------------------------------------------
    void merge(const sorttree* const that) WDutils_THROWING {
      if(that->Ncrit != this->Ncrit)
	WDutils_THROW("sorttree::merge(): Ncrit differs: cannot merge");
      if(that->L0==0 || this->L0==0)
	WDutils_THROW("sorttree::merge(): leafs deleted: cannot merge");
      if(that->ROOT->C != this->ROOT->C ||
	 that->ROOT->S != this->ROOT->S)
	WDutils_THROW("sorttree::merge(): ranges differ: cannot merge");
      Ncrit;
      merge(ROOT,that->ROOT);
    }
  };
#undef PLEAF
#undef PRNGE
} // namespace {
////////////////////////////////////////////////////////////////////////////////
namespace WDutils {
  template<typename sortable, typename sortit> inline
  void set_value_table(sortit         const&A,      // I: objects to be sorted  
		       int            const&nA,     // I: # objects             
		       const sortable*const&rank,   // I: table of ranks        
		       sortable      *const&value,  // O: table of values       
		       int            const&nt,     // I: size of these tables  
		       const sortable*xmin,         //[I: minimum value]        
		       const sortable*xmax)         //[I: maximum value]        
  {
    if(xmin && xmax) {
      sorttree<sortable> ST(A,nA,nA/(3*nt),*xmin,*xmax);
      ST.make_value_table(rank,value,nt);
      if(0    == rank[0]   ) value[0]    = *xmin;
      if(nA-1 == rank[nt-1]) value[nt-1] = *xmax;
    } else {
      register sortable Xmin=A[0], Xmax=Xmin;
      for(register int i=1; i!=nA; ++i) {
	if(A[i] < Xmin) Xmin = A[i];
	if(A[i] > Xmax) Xmax = A[i];
      }
      sorttree<sortable> ST(A,nA,nA/(3*nt),Xmin,Xmax);
      ST.make_value_table(rank,value,nt);
      if(0    == rank[0]   ) value[0]    = Xmin;
      if(nA-1 == rank[nt-1]) value[nt-1] = Xmax;
    }
  }
} // namespace WDutils {
////////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_numerics_h
