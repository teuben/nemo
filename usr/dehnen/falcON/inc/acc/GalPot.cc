// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// GalPot.cc                                                                   |
//                                                                             |
// Copyright (C) 1996-2007 Walter Dehnen                                       |
//                                                                             |
// This program is free software; you can redistribute it and/or modify        |
// it under the terms of the GNU General Public License as published by        |
// the Free Software Foundation; either version 2 of the License, or (at       |
// your option) any later version.                                             |
//                                                                             |
// This program is distributed in the hope that it will be useful, but         |
// WITHOUT ANY WARRANTY; without even the implied warranty of                  |
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU           |
// General Public License for more details.                                    |
//                                                                             |
// You should have received a copy of the GNU General Public License           |
// along with this program; if not, write to the Free Software                 |
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                   |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// source code for class GalaxyPotential and dependencies                      |
//                                                                             |
// TO BE COMPILED AND LINKED TO THE ENDUSER                                    |
//                                                                             |
// Version 0.0    15. July      1997                                           |
// Version 0.1    24. March     1998                                           |
// Version 0.2    22. September 1998                                           |
// Version 0.3    07. June      2001                                           |
// Version 0.4    22. April     2002                                           |
// Version 0.5    05. December  2002                                           |
// Version 0.6    05. February  2003                                           |
// Version 0.7    23. September 2004  fixed "find(): x out of range" error     |
// Version 0.8    24. June      2005  explicit construction of tupel           |
// Version 0.9    06. November  2007  consistent with GalPot package           |
//                                                                             |
//-----------------------------------------------------------------------------+
#define  GalPot_cc
#ifndef GalPot_h
#  include "GalPot.h"	                          // make sure GalPot.h is known
#endif
namespace GalPot {                                  // v0.4                     
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //  Part 1: auxiliary functions and classes that are not already defined    //
  //          in GalPot.h & GalPot_pre.h                                      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  const double
    Pi  = 3.14159265358979323846264338328,
    TPi = 2*Pi,
    FPi = 4*Pi;

  inline int    sign  (double x) { return x<0.? -1 : x>0? 1:0; }
  inline double square(double x) { return x*x; }
  inline double cube  (double x) { return x*x*x; }

  inline void SwallowRestofLine(std::istream& from)
  {
    char c;
    do from.get(c); while( from.good() && c!='\n');
  }

  template <class ALLOCTYPE>
  inline int Alloc1D(ALLOCTYPE* &A, const int N)
  {
    A = new ALLOCTYPE[N]; 
    if(!A) return 1;
    return 0;
  }

  template <class ALLOCTYPE>
  inline int Alloc2D(ALLOCTYPE** &A, const int N[2])
  {
    register int i, iN1;
    A    = new ALLOCTYPE*[N[0]];	 if(!A) return 1;
    A[0] = new ALLOCTYPE [N[0]*N[1]];  if(!A[0]) return 1;
    for(i=1, iN1=N[1]; i<N[0]; i++,iN1+=N[1])
      A[i] = A[0] + iN1;
    return 0;
  }

  template <class ALLOCTYPE>
  inline void Free1D(ALLOCTYPE* A)
  { 
    delete[] A; 
  }

  template <class ALLOCTYPE>
  inline void Free2D(ALLOCTYPE** A)
  {
    delete[] A[0];
    delete[] A;
  }

  inline void error(const char* msgs)
  {
    std::cerr<< "GalPot ERROR: " << msgs << std::endl;
    std::exit(1);
  }

  inline void warning(const char* msgs)
  {
    std::cerr<< "GalPot WARNING: " << msgs << std::endl;
  }

  template<class S>
  int hunt(const S*xarr, const int n, const S x, const int j)
    // hunts the ordered table xarr for jlo such that xarr[jlo]<=x<xarr[jlo+1]
    // on input j provides a guess for the final value of jlo.
    // for an ascendingly ordered array, we return
    //  -1 for         x < x[0]
    //   i for x[i] <= x < x[i+1]  if  0<=i<n
    // n-1 for         x == x[n-1]
    // n   for         x >  x[n-1]
  {
    int jm,jlo=j,jhi,l=n-1;
    int ascnd=(xarr[l]>xarr[0]);

    if(!ascnd && xarr[l]==xarr[0] ) return -1;  // x_0 = x_l
    if((ascnd && x<xarr[0]) || (!ascnd && x>xarr[0]) ) return -1;
    if((ascnd && x>xarr[l]) || (!ascnd && x<xarr[l]) ) return  n;

    if(jlo<0 || jlo>l) {                   // input guess not useful,
      jlo = -1;                            //    go to bisection below
      jhi = n;
    } else {
      int inc = 1;
      if((x>=xarr[jlo]) == ascnd) {          // hunt upward
	if(jlo == l) return (x==xarr[l])? l : n;
	jhi = jlo+1;
	while((x>=xarr[jhi]) == ascnd) {     // not done hunting
	  jlo =jhi;
	  inc+=inc;                        // so double the increment
	  jhi =jlo+inc;
	  if(jhi>l) {                      // off end of table
	    jhi=n;
	    break;
	  }
	}
      } else {                             // hunt downward
	if(jlo == 0) return -1;
	jhi = jlo;
	jlo-= 1;
	while((x<xarr[jlo]) == ascnd) {      // not done hunting
	  jhi = jlo;
	  inc+= inc;                       // so double the increment
	  jlo = jhi-inc;
	  if(jlo < 0) {                    // off end of table
	    jlo = 0;
	    break;
	  }
	}
      }
    }
    while (jhi-jlo != 1) {                 // bisection phase
      jm=(jhi+jlo) >> 1;
      if((x>=xarr[jm]) == ascnd) jlo=jm;
      else jhi=jm;
    }
    return jlo;
  }

  template<class S>
  inline void find(int& klo, const int n, S *x, const S xi)
    // for an ascendingly ordered array, we return
    //  -1 for         x < x[0]
    //   i for x[i] <= x < x[i+1]  if  0<=i<n
    // n-2 for         x == x[n-1]
    // n   for         x >  x[n-1]
  {
    if(klo<0 || klo>=n-1 || x[klo]>xi || x[klo+1]<xi) {
      klo = int( (xi-x[0]) / (x[n-1]-x[0]) * (n-1) );
      klo = hunt(x,n,xi,klo);
      if(klo<0 || klo>=n) error("find(): x out of range");
      if(klo == n-1) klo = n-2;                         // v.07
    }
  }

  template<class S, class T>
  void spline(S*  const x,	// input:   table of points
	      T*  const y,	// input:   table of function values
	      int const n,	// input:   size of above tables
	      T   const yp1,	// input:   y'(x[0])
	      T   const ypn,	// input:   y'(x[n-1])
	      T*        y2,	// output:  table of y''
	      int const nat1=0,	// input:   natural spline at x[0] ?
	      int const natn=0)	// input:   natural spline at x[n-1] ?
  { 
    const S  zero=0., half=0.5, one=1., two=2., three=3., six=6.;
    register int i;
    register S   qn,p,sig,dx,dx1,dx2;
    T   un,dy,dy1;
    T *u = new T[n-1];
    S *v = new S[n-1];
    dx   = x[1] - x[0];
    dy   = y[1] - y[0];
    if(nat1)
      u[0] = v[0] = zero;
    else {
      v[0] =-half;
      u[0] = three/dx * (dy/dx-yp1);
    }
    for(i=1; i<n-1; i++) {
      dx1  = x[i+1]-x[i];
      dx2  = x[i+1]-x[i-1];
      dy1  = y[i+1]-y[i];
      sig  = dx/dx2;
      p    = sig*v[i-1]+two;
      v[i] = (sig-one)/p;
      u[i] = ( six*(dy1/dx1-dy/dx)/dx2 - sig*u[i-1] ) / p;
      dx   = dx1;
      dy   = dy1;
    }
    if(natn)
      un = qn = zero;
    else {
      qn = half;
      un = three/dx * (ypn - dy/dx);
    }
    y2[n-1]= (un-qn*u[n-2]) / (qn*v[n-2]+one);
    for(i=n-2; i>=0; i--)
      y2[i] = v[i]*y2[i+1] + u[i];
    delete[] u;
    delete[] v;
  }

  template<class S, class T>
  void splinTarr(
		 S   const xl,		// input:   xl < x
		 S   const xh,		// input:   x < xh
		 S   const xi,		// input:   xlo <= xi < xhi
		 T*  const yl,	   	// input:   y_k(xl)
		 T*  const yh,		// input:   y_k(xh)
		 T*  const y2l,		// input:   y2_k(xl)
		 T*  const y2h,		// input:   y2_k(xh)
		 int const K,		// input:   k=0,...,K-1
		 T*  y,			// output:  y_k(xi)
		 T*  dy=0,		// output:  dy_k/dx(xi)     if dy  != 0
		 T*  d2y=0)		// output:  d^2y_k/d^2x(xi) if d2y != 0
  {
    // computes K splines simultaneously
    const    S   zero=0., one=1., three=3., six=6.;
    register T   *Y=y, *Yl=yl, *Yh=yh, *Y2l=y2l, *Y2h=y2h, *YK=y+K;
    register S   h,h6,hh,A,B,Aq,Bq,Ap,Bp;
    if((h=xh-xl)==zero) error("splinTarr bad X input");
    h6=h/six;
    hh=h*h6;
    A =(xh-xi)/h;  Aq=A*A;  Ap=(Aq-one)*A*hh;
    B =one-A;      Bq=B*B;  Bp=(Bq-one)*B*hh;
    if(d2y) {
      register T *dY=dy, *d2Y=d2y;
      register S Au=h6*(three*Aq-one), Bu=h6*(three*Bq-one);
      for(; Y<YK; Y++,Yl++,Yh++,Y2l++,Y2h++,dY++,d2Y++) {
	*Y   = A**Yl+B**Yh+Ap**Y2l+Bp**Y2h;
	*dY  = (*Yh-*Yl)/h + Bu**Y2h-Au**Y2l;
	*d2Y = A**Y2l + B**Y2h;
      }
    } else if(dy) {
      register T *dY=dy;
      register S Au=h6*(three*Aq-one), Bu=h6*(three*Bq-one);
      for(; Y<YK; Y++,Yl++,Yh++,Y2l++,Y2h++,dY++) {
	*Y   = A**Yl+B**Yh+Ap**Y2l+Bp**Y2h;
	*dY  = (*Yh-*Yl)/h + Bu**Y2h-Au**Y2l;
      }
    } else {
      for(; Y<YK; Y++,Yl++,Yh++,Y2l++,Y2h++)
	*Y   = A**Yl+B**Yh+Ap**Y2l+Bp**Y2h;
    }
  }

  template<class S, class T>
  void Pspline(
	       S*  const x,		// input:   table of points x
	       T*  const Y,		// input:   table of y(x)
	       T*  const Y1,		// input:   table of dy/dx
	       int const n,		// input:   size of these tables
	       T *Y3)			// output:  table of d^3y/dx^3
  // Given dy/dx on a grid, d^3y/dx^3 is computed such that the (unique) 
  // polynomials of 5th order between two adjacent grid points that give y,dy/dx
  // and d^3y/dx^3 on the grid are continuous in d^2y/dx^2, i.e. give the same
  // value at the grid points. At the grid boundaries  d^3y/dx^3=0  is adopted.
  {
    const S      zero=0.,one=1.,three=3.,seven=7.,ten=10.,twelve=12.; 
    register int i;
    register S   p,sig,dx,dx1,dx2;
    register T   dy=Y[1]-Y[0], dy1=dy;
    S *v = new S[n-1];
    dx   = x[1]-x[0];
    Y3[0]= v[0] = zero;
    for(i=1; i<n-1; i++) {
      dx1  = x[i+1]-x[i];
      dx2  = x[i+1]-x[i-1];
      dy1  = Y[i+1]-Y[i];
      sig  = dx/dx2;
      p    = sig*v[i-1]-three;
      v[i] = (sig-one)/p;
      Y3[i]= twelve*(   seven*Y1[i]*dx2/(dx*dx1) 
			+ three*(Y1[i-1]/dx+Y1[i+1]/dx1)
		        - ten*(dy/(dx*dx) + dy1/(dx1*dx1))  ) / dx2;
      Y3[i]= (Y3[i] - sig*Y3[i-1] ) / p;
      dx   = dx1;
      dy   = dy1;
    }
    Y3[n-1] = zero;
    for(i=n-2; i>=0; i--)
      Y3[i] += v[i]*Y3[i+1];
    delete[] v;
  }

  template<class S, class T>
  T PsplinT(				// return:  y(xi)
	    const S& x0,		// input:   x0 <= x
	    const S& x1,		// input:   x  <= x1 
	    const T& Y0,		// input:   Y(x0)
	    const T& Y1,		// input:   Y(x1)
	    const T& Y10,		// input:   dY(x0)
	    const T& Y11,		// input:   dY(x1)
	    const T& Y30,		// input:   d3Y(x0)
	    const T& Y31,		// input:   d3Y(x0)
	    const S& xi,		// input:   x-value where y is wanted
	    T* dYi=0,			// output:  dy/dx(xi)     if dy  != 0
	    T* d2Yi=0)			// output:  d^2y/d^2x(xi) if d2y != 0
  {
    const    S zero=0.,one=1.,two=2.,five=5.,six=6.,nine=9.,fe=48.;
    register S h,hi,hf, A,B,C,D,Aq,Bq;
    if((h=x1-x0)==zero) error("PsplinT bad X input");
    hi = one/h;
    hf = h*h;
    A  = hi*(x1-xi); Aq = A*A;
    B  = one-A;      Bq = B*B;
    C  = h*Aq*B;
    D  =-h*Bq*A;
    register T 	t1 = hi*(Y1-Y0),
      C2 = Y10-t1,
      C3 = Y11-t1,
      t2 = six*(Y10+Y11-t1-t1)/hf,
      C4 = Y30-t2,
      C5 = Y31-t2;
    hf/= fe;
    register T 	Yi = A*Y0+ B*Y1+ C*C2+ D*C3+ 
      hf*(C*(Aq+Aq-A-one)*C4+ D*(Bq+Bq-B-one)*C5);
    if(dYi) {
      register S BAA=B-A-A, ABB=A-B-B;
      hf  += hf;
      *dYi = t1 + (A*ABB)*C2 + (B*BAA)*C3
	+ hf*A*B*((one+A-five*Aq)*C4+ (one+B-five*Bq)*C5);
      if(d2Yi) {
	*d2Yi = BAA*C2 - ABB*C3;
	*d2Yi+= *d2Yi  + hf * ( (two*Aq*(nine*B-A)-one) * C4
				+(two*Bq*(B-nine*A)+one) * C5 );
	*d2Yi*= hi;
      }
    }
    return Yi;
  }

  template<class S, class T>
  void PsplinTarr(
		  const S& xl,		// input:   xl <= x
		  const S& xh,		// input:   x  <= xh 
		  const S& xi,		// input:   x-value where y is wanted
		  T* const yl,		// input:   Y_k(xl)
		  T* const yh,		// input:   Y_k(xh)
		  T* const y1l,		// input:   dY_k(xl)
		  T* const y1h,		// input:   dY_k(xh)
		  T* const y3l,		// input:   d3Y_k(xl)
		  T* const y3h,		// input:   d3Y_k(xh)
		  const int K,		// input:   k=0,...,K-1
		  T* y,			// output:  y_k(xi)
		  T* dy=0,		// output:  dy_k/dx(xi)     if dy  != 0
		  T* d2y=0)		// output:  d^2y_k/d^2x(xi) if d2y != 0
  {
    const    S zero=0.,one=1.,two=2.,five=5.,six=6.,nine=9.,fe=48.;
    register S h,hi,hq,hf, A,B,C,D,E,F,Aq,Bq;
    register T C2=*yl,C3=C2,C4=C2,C5=C2, t1=C2,t2=C2;
    register T *Y=y,*Yl=yl,*Yh=yh,*Y1l=y1l,*Y1h=y1h,*Y3l=y3l,*Y3h=y3h,*YK=y+K;
    if((h=xh-xl)==zero) error("PsplinTarr(): bad X input");
    hi = one/h;
    hq = h*h;
    hf = hq/fe;
    A  = hi*(xh-xi); Aq = A*A;
    B  = one-A;      Bq = B*B;
    C  = h*Aq*B;
    D  =-h*Bq*A;
    E  = hf*C*(Aq+Aq-A-one);
    F  = hf*D*(Bq+Bq-B-one);
    if(d2y) {
      register S hf2= hf+hf, BAA=B-A-A, ABB=A-B-B,
	AB=A*B, ABh=hf2*AB, Cp=Aq-AB-AB, Dp=Bq-AB-AB,
	Ep=ABh*(one+A-five*Aq), Fp=ABh*(one+B-five*Bq),
	Epp=hf2*(two*Aq*(nine*B-A)-one), Fpp=hf2*(two*Bq*(B-nine*A)+one);
      register T *dY=dy, *d2Y=d2y;
      for(; Y<YK; Y++,Yl++,Yh++,Y1l++,Y1h++,Y3l++,Y3h++,dY++,d2Y++) {
	t1   = hi*(*Yh-*Yl);
	C2   = *Y1l-t1;
	C3   = *Y1h-t1;
	t2   = six*(*Y1l+*Y1h-t1-t1)/hq;
	C4   = *Y3l-t2;
	C5   = *Y3h-t2;
	*Y   = A**Yl + B**Yh + C*C2 + D*C3+ E*C4 + F*C5;
	*dY  = t1+ Cp*C2 + Dp*C3 + Ep*C4 + Fp*C5;
	*d2Y = BAA*C2 - ABB*C3;
	*d2Y+= *d2Y + Epp*C4 + Fpp*C5;
	*d2Y*= hi;
      }
    } else if(dy) {
      register S AB=A*B, ABh=(hf+hf)*AB, Cp=Aq-AB-AB, Dp=Bq-AB-AB,
	Ep=ABh*(one+A-five*Aq), Fp=ABh*(one+B-five*Bq);
      register T *dY=dy;
      for(; Y<YK; Y++,Yl++,Yh++,Y1l++,Y1h++,Y3l++,Y3h++,dY++) {
	t1  = hi*(*Yh-*Yl);
	C2  = *Y1l-t1;
	C3  = *Y1h-t1;
	t2  = six*(*Y1l+*Y1h-t1-t1)/hq;
	C4  = *Y3l-t2;
	C5  = *Y3h-t2;
	*Y  = A**Yl + B**Yh + C*C2 + D*C3+ E*C4 + F*C5;
	*dY = t1+ Cp*C2 + Dp*C3 + Ep*C4 + Fp*C5;
      }
    } else {
      for(; Y<YK; Y++,Yl++,Yh++,Y1l++,Y1h++,Y3l++,Y3h++) {
	t1  = hi*(*Yh-*Yl);
	C2  = *Y1l-t1;
	C3  = *Y1h-t1;
	t2  = six*(*Y1l+*Y1h-t1-t1)/hq;
	C4  = *Y3l-t2;
	C5  = *Y3h-t2;
	*Y  = A**Yl+ B**Yh+ C*C2+ D*C3+ E*C4+ F*C5;
      }
    }
  }

  template<class S, class T>
  void Pspline2D(
		 S*  const x[2],  // input:  tables of points x0, x1
		 T** const y[3],  // input:  tables of y, dy/dx0, dy/dx1
		 int const n[2],  // input:  sizes of above tables: n[0],n[1]
		 T** a[4])        // output: tables: coeffs a[0],a[1],a[2],a[3]
  {
    // 2D Pspline with natural boundary conditions
    register   T   z=y[0][0][0];
    z = 0.;
    register int i,j;
    T *t = new T[n[0]];
    T *t1= new T[n[0]];
    T *t3= new T[n[0]];
    // 1. for each x1 do 1D Pspline for y in x0
    for(j=0; j<n[1]; j++) {
      for(i=0; i<n[0]; i++) {
	t[i]  = y[0][i][j];		// y
	t1[i] = y[1][i][j];		// dy/dx0
      }
      Pspline(x[0],t,t1,n[0],t3);
      for(i=0; i<n[0]; i++)
	a[0][i][j] = t3[i];		// d^3y/dx0^3
    }
    // 2. for each x0 do 
    // 1D Pspline for y and splines for dy/dx0, d^3y/dx0^3 in x1
    for(i=0; i<n[0]; i++) {
      Pspline(x[1],y[0][i],y[2][i],n[1],a[1][i]);
      spline (x[1],y[1][i],n[1],z,z,a[2][i],0,1);
      spline (x[1],a[0][i],n[1],z,z,a[3][i],0,1);
    }
    delete[] t;
    delete[] t1;
    delete[] t3;
  }

  template<class S, class T>
  T Psplev2D(		      // return:  y(x0i,x1i)
	     S*  const x[2],  // input:   tables of points x0, x1
	     T** const y[3],  // input:   tables of y, dy/dx0, dy/dx1
	     T** const a[4],  // input:   tables of coeffs a[0],a[1],a[2],a[3]
	     int const n[2],  // input:   sizes of above tables: n[0],n[1]
	     S   const xi[2], // input:   (x0,x1)-value where y is wanted
	     T*  d1=0,	      // output:  gradient of y   if d1 != 0
	     T** d2=0)	      // output:  d^2y/dxi/dxj    if d2 != 0
  {
    static int l0=0, l1=0;
    find(l0,n[0],x[0],xi[0]);
    find(l1,n[1],x[1],xi[1]);
    register int k0=l0+1, k1=l1+1;

    T fl[2] ={y[0][l0][l1], y[0][k0][l1]}, fh[2] ={y[0][l0][k1], y[0][k0][k1]},
      f1l[2]={y[2][l0][l1], y[2][k0][l1]}, f1h[2]={y[2][l0][k1], y[2][k0][k1]},
      f3l[2]={a[1][l0][l1], a[1][k0][l1]}, f3h[2]={a[1][l0][k1], a[1][k0][k1]},
      flo[4]={y[1][l0][l1], y[1][k0][l1], a[0][l0][l1], a[0][k0][l1]},
      fhi[4]={y[1][l0][k1], y[1][k0][k1], a[0][l0][k1], a[0][k0][k1]},
      f2l[4]={a[2][l0][l1], a[2][k0][l1], a[3][l0][l1], a[3][k0][l1]},
      f2h[4]={a[2][l0][k1], a[2][k0][k1], a[3][l0][k1], a[3][k0][k1]};
    if(d2) {
      T F[2], G[4], dF[2], dG[4], d2F[2], d2G[4];
      PsplinTarr (x[1][l1],x[1][k1],xi[1],fl,fh,f1l,f1h,f3l,f3h,2,F,dF,d2F);
      splinTarr  (x[1][l1],x[1][k1],xi[1],flo,fhi,f2l,f2h,4,G,dG,d2G);
      d2[1][1]=PsplinT(x[0][l0],x[0][k0],d2F[0],d2F[1],
		       d2G[0],d2G[1],d2G[2],d2G[3],xi[0]);
      d1[1]   =PsplinT(x[0][l0],x[0][k0],dF[0],dF[1],
		       dG[0],dG[1],dG[2],dG[3],xi[0],d2[1]);
      d2[0][1]=d2[1][0];
      return   PsplinT(x[0][l0],x[0][k0],F[0],F[1],
		       G[0],G[1],G[2],G[3],xi[0],d1,d2[0]);
    } else if(d1) {
      T F[2], G[4], dF[2], dG[4];
      PsplinTarr(x[1][l1],x[1][k1],xi[1],fl,fh,f1l,f1h,f3l,f3h,2,F,dF);
      splinTarr (x[1][l1],x[1][k1],xi[1],flo,fhi,f2l,f2h,4,G,dG);
      d1[1] =PsplinT(x[0][l0],x[0][k0],dF[0],dF[1],
		     dG[0],dG[1],dG[2],dG[3],xi[0]);
      return PsplinT(x[0][l0],x[0][k0],F[0],F[1],
		     G[0],G[1],G[2],G[3],xi[0],d1);
    }
    T F[2], G[4];
    PsplinTarr(x[1][l1],x[1][k1],xi[1],fl,fh,f1l,f1h,f3l,f3h,2,F);
    splinTarr (x[1][l1],x[1][k1],xi[1],flo,fhi,f2l,f2h,4,G);
    return PsplinT(x[0][l0],x[0][k0],F[0],F[1],G[0],G[1],G[2],G[3],xi[0]);
  }

  template<class S, int N>
  void LegendrePeven(tupel<N,S>& p, const double x)
    // based on a routine from J.J. Binney
    // evaluates even Legendre Polys up to l=2*(N-1) at x
  {
    register int    n,l,l2;
    register double x2=x*x;
    p[0] = 1.;
    p[1] = 1.5*x2-0.5;
    for(n=2; n<N; n++) {
      l    = 2*(n-1);
      l2   = 2*l;
      p[n] = - p[n-2] * l*(l-1)/double((l2+1)*(l2-1))
	+ p[n-1] * (x2 - (l2*l+l2-1)/double((l2-1)*(l2+3)) );
      p[n]*= (l2+1)*(l2+3) / double((l+1)*(l+2));
    }
  }

  template<class S, int N>
  void dLegendrePeven(tupel<N,S>& p, tupel<N,S>& d, const double x)
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
      l2   = 2*l;
      p[n] = - p[n-2] * l*(l-1)/double((l2+1)*(l2-1))
	+ p[n-1] * (x2 - (l2*l+l2-1)/double((l2-1)*(l2+3)) );
      p[n]*= (l2+1)*(l2+3) / double((l+1)*(l+2));
      d[n] = - d[n-2] * l*(l-1)/double((l2+1)*(l2-1))
	+ d[n-1] * (x2 - (l2*l+l2-1)/double((l2-1)*(l2+3)) )
	+ p[n-1];
      d[n]*= (l2+1)*(l2+3) / double((l+1)*(l+2));
    }
    x2 = 2*x;
    for(n=0; n<N; n++)
      d[n] *= x2;
  }

  template<class FUNC>
  typename FUNC::ValType qbulir(FUNC const            &func,
				typename FUNC::ValType a,
				typename FUNC::ValType b,
				typename FUNC::ValType eps_)
  {
    typename FUNC::ValType ba=b-a;
    if(ba==0.) return 0.;

    int                    i,n=2,nn=3,mx=25,m,mr, bo,bu=0,odd=1;
    typename FUNC::ValType d[7],dt[7];
    typename FUNC::ValType c,d1,ddt,den,e,eps,err,eta=1.e-7,
                           gr,hm,nt,sm,t,t1,t2,t2a,ta,tab=0.,tb,v=0.,w;

    while(eta+1. != 1.) eta *= 0.5;
    eta  *= 2.;                    // eta = actual computing accuracy

    eps   = max(eps_,eta);
    sm    = 0.;
    gr    = 0.;
    t1    = 0.;
    t2    = 0.5*(func(a)+func(b));
    t2a   = t2;
    tb    = abs(t2a);
    c     = t2*ba;
    dt[0] = c;

    for(m=1;m<=mx;m++) {            // iterate over the refinements
      bo = (m>=7);
      hm = ba/n;
      if(odd) {
	for(i=1;i<=n;i+=2) {
	  w  = func(a+i*hm);
	  t2+= w;
	  tb+= abs(w);
	}
	nt   = t2;
	tab  = tb * abs(hm);
	d[1] = 16./9.;
	d[3] = 64./9.;
	d[5] = 256./9.;
      } else {
	for(i=1;i<=n;i+=6) {
	  w  = i*hm;
	  t1+= func(a+w) + func(b-w);
	}
	nt   = t1+t2a;
	t2a  = t2;
	d[1] = 2.25;
	d[3] = 9.;
	d[5] = 36.;
      }
      ddt   = dt[0];
      t     = nt*hm;
      dt[0] = t;
      nt    = dt[0];
      if(bo) {
	mr   = 6;
	d[6] = 64.;
	w    = 144.;
      } else {
	mr   = m;
	d[m] = n*n;
	w    = d[m];
      }
      for(i=1;i<=mr;i++) {
	d1  = d[i]*ddt;
	den = d1-nt;
	e   = nt-ddt;
	if(den != 0) {
	  e /= den;
	  v  = nt*e;
	  nt = d1*e;
	  t += v;
	} else {
	  nt = 0.;
	  v  = 0.;
	}
	ddt   = dt[i];
	dt[i] = v;
      }
      ta = c;
      c  = t;
      if(!bo) t -= v;
      v   = t-ta;
      t  += v;
      err = abs(v);
      if(ta<t) {
	d1 = ta;
	ta = t;
	t  = d1;
      }
      bo = bo || (ta<gr && t>sm);
      if(bu && bo && err < eps*tab*w) break;
      gr   = ta;
      sm   = t;
      odd  = !odd;
      i    = n;
      n    = nn;
      nn   = i+i;
      bu   = bo;
      d[2] = 4.;
      d[4] = 16.;
    }
    v = tab*eta;
    if(m==mx) std::cerr << " qbulir exceeding maximum of iterations\n";
    return c;
  }

  template<class C, class S=double> struct Adaptor {
    typedef S ValType;
    typedef const C* c_pter_C;
    typedef S (C::*f_pter_C)(S) const;
    const c_pter_C o;
    const f_pter_C f;
    Adaptor(c_pter_C __o, f_pter_C __f) : o(__o), f(__f) {}
    S operator() (S x) const { return (o->*f)(x); }
  };

  void GaussLegendre(double *x, double *w, int n)
  {
    register double eps=1.e-10;
    for(register double ep1=1.0+eps; 1.!=ep1; eps*=0.5, ep1=1.0+eps) ;
    eps  *=2.;                       // eps = actual computing accuracy
    register int j,i,m=(n+1)/2;
    register double z1,z,pp,p3,p2,p1;
    for (i=0;i<m;i++) {
      z=cos(Pi*(i+0.75)/(n+0.5));
      do {
	p1 = 1.0;
	p2 = 0.0;
	for(j=0;j<n;j++) {
	  p3 = p2;
	  p2 = p1;
	  p1 = ( (2*j+1)*z*p2 - j*p3 ) / double(j+1);
	}
	pp = n * (z*p1-p2) / (z*z-1.0);
	z1 = z;
	z  = z1 - p1 / pp;
      } while (abs(z-z1)>eps);
      x[i]     =-z;
      x[n-1-i] = z;
      w[i]     = 2. / ((1.0-z*z)*pp*pp);
      w[n-1-i] = w[i];
    }
  }

  double I0(const double x)
  {
    register double ax=abs(x),y;
    if(ax < 3.75) {
      y = x/3.75;
      y*= y;
      return 1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
		+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
    }
    y = 3.75/ax;
    return (exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
	   +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
	   +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
	   +y*0.392377e-2))))))));
  }

  double I1(const double x)
  {
    register double ans,ax=abs(x),y;
    if(ax < 3.75) {
      y = x/3.75;
      y*= y;
      ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
		 +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
    } else {
      y=3.75/ax;
      ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
		      -y*0.420059e-2));
      ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
		    +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
      ans *= (exp(ax)/sqrt(ax));
    }
    return x < 0.0 ? -ans : ans;
  }

  double K0(const double x)
  {
    if(x<0.) error("negative argument in K0(x)");
    register double y;
    if(x <= 2.) {
      y = x*x/4.;
      return (-log(x/2.0)*I0(x))+(-0.57721566+y*(0.42278420
	    +y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
	    +y*(0.10750e-3+y*0.74e-5))))));
    }
    y = 2./x;
    return (exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
	  +y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
	  +y*(-0.251540e-2+y*0.53208e-3))))));
  }

  double K1(const double x)
  {
    if(x<0.) error("negative argument in K1(x)");
    register double y;
    if(x <= 2.) {
      y=x*x/4.0;
      return (log(x/2.0)*I1(x))+(1.0/x)*(1.0+y*(0.15443144
	    +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
	    +y*(-0.110404e-2+y*(-0.4686e-4)))))));
    }
    y = 2./x;
    return (exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
	  +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
	  +y*(0.325614e-2+y*(-0.68245e-3)))))));
  }

  double Kn(const int n, const double x)
  {
    if(n<0)  error("negative n in Kn(x)");
    if(x<0.) error("negative argument in Kn(x)");
    if(n==0) return K0(x);
    if(n==1) return K1(x);
    register int j;
    register double bk,bkm,bkp,tox;

    tox = 2./x;
    bkm = K0(x);
    bk  = K1(x);
    for(j=1; j<n; j++) {
      bkp = bkm+j*tox*bk;
      bkm = bk;
      bk  = bkp;
    }
    return bk;
  }
}                                               // namespace GalPot; v0.4       
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//         Part 2: member functions of classes defined in GalPot.h            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
using namespace GalPot;                         // v0.4                         
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class DiskAnsatz                                                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//   rho(R,z) = Sig(R) * h(z)                                                 //
//                                                                            //
// with                                                                       //
//                  R0   R            R                                       //
//   Sig(R) = exp(- -- - -- - eps cos(--))                                    //
//                  R    Rd           Rd                                      //
//                                                                            //
// and                                                                        //
//                                                                            //
//   h(z)   = delta(z)                   for  d=0                             //
//            (1/2 d)  * exp(-|z/d|)     for  d>0                             //
//            (1/4|d|) * sech^2(|z/2d|)  for  d<0                             //
//                                                                            //
// The potential part returned by operator() amounts to                       //
//                                                                            //
//                                                                            //
//   Phi(r,z) = Sig(r) * H(z)                                                 //
//                                                                            //
// where r = sqrt(R^2+z^2) is the spherical polar radius and H(z) is the      //
// second integral of h(z).                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

void DiskAnsatz::setup(const DiskPar& p)
{
  S0 = abs(p(0));
  Rd = abs(p(1));
  zd = abs(p(2));
  R0 = abs(p(3));
  eps= p(4);
  if(zd==0.)  thin       = 1;
  else    thin       = 0;
  if(R0==0.)  hollow     = 0;
  else    hollow     = 1;
  if(p(2)<0.) isothermal = 1;
  else    isothermal = 0;
  zdoRd = zd/Rd;
  R0oRd = R0/Rd;
  Rd2   = Rd*Rd;
  fac   = TPi*Grav*S0;
}

double DiskAnsatz::SurfaceDensity(double R) const
{
  if(hollow && R==0.) return 0.;
  const double y=R/Rd;
  return eps? S0*exp(-R0/R-y+eps*cos(y)) : S0*exp(-R0/R-y);
}

inline double DiskAnsatz::mass_integrand(double y) const
{
  if(y<=0. || y>=1.) return 0.;
  const double y1=1.-y, x=y/y1;
  return eps? exp(-R0oRd/x-x+eps*cos(x))*x/(y1*y1) : exp(-R0oRd/x-x)*x/(y1*y1);
}

double DiskAnsatz::mass(double R) const
{
  const double F=TPi*S0*Rd2;
  if(R<=0.) {         // give total mass
    if(eps)
      return F*qbulir(Adaptor<DiskAnsatz>(this,&DiskAnsatz::mass_integrand),
		      0.,1.,1.e-6);
    if(hollow) return FPi*S0*R0*Rd*Kn(2,2.*sqrt(R0/Rd));
    return F;
  }                   // give mass inside R (but integrate z from -oo to oo)
  return F*qbulir(Adaptor<DiskAnsatz>(this,&DiskAnsatz::mass_integrand),
		  0.,R/(Rd+R),1.e-6);
}

double DiskAnsatz::Density(double R, double z) const
{
  double rhR;
  if(eps)         rhR= (hollow && R==0.)? 0. : exp(-R0/R-R/Rd+eps*cos(R/Rd));
  else if(hollow) rhR= (R==0.)? 0. : exp(-R0/R-R/Rd);
  else            rhR= exp(-R/Rd);
  if(thin)                            // vertically thin disk: return Sigma
    return S0*rhR;
  else if(isothermal) {               // vertically isothermal disk
    double x=abs(z/zd), ex=exp(-x), ex1=1.+ex;
    return S0*rhR * ex/(ex1*ex1*zd);
  }                                   // vertically exponential disk
  double x=abs(z/zd), ex=exp(-x);
  return S0*rhR * 0.5*ex/zd;
}

double DiskAnsatz::Residual(double r, double st, double ct) const
// gives aimed Laplace(Phi_multipole)
{
  if(ct==0. || S0==0.) return 0;
  register double R=r*st, z=r*ct, g,gp,gpp, F,f,fp,fpp;
  // deal with the vertical part
  if(thin) {
    g   = abs(z);
    gp  = sign(z);
    gpp = 0.;
  } else if(isothermal) {
    register double x,sh1;
    x   = abs(z/zd);
    gpp = exp(-x);
    sh1 = 1.+gpp;
    g   = 2*zd*(0.5*x+log(0.5*sh1));
    gp  = sign(z)*(1.-gpp)/sh1;
    gpp/= 0.5*sh1*sh1*zd;
  } else {
    register double x;
    x   = abs(z/zd);
    gpp = exp(-x);
    g   = zd*(gpp-1+x);
    gp  = sign(z)*(1.-gpp);
    gpp/= zd;
  }
  // deal with the radial part
  if(hollow && r==0.) F=f=fp=fpp=0.;
  else if(eps) {
    if(hollow) {
      register double rq=r*r,cr=cos(r/Rd),sr=sin(r/Rd);
      F   = (R==0)? 0. : exp(-R0/R-R/Rd+eps*cos(R/Rd));
      f   = exp(-R0/r-r/Rd+eps*cr);
      fp  = R0/rq-(1.+eps*sr)/Rd;
      fpp = (fp*fp-2.*R0/(rq*r)-eps*cr/Rd2)*f;
      fp *= f;
    } else {
      register double cr=cos(r/Rd),sr=sin(r/Rd);
      F   = exp(-R/Rd+eps*cos(R/Rd));
      f   = exp(-r/Rd+eps*cr);
      fp  = -(1.+eps*sr)/Rd;
      fpp = (fp*fp-eps*cr/Rd2)*f;
      fp *= f;
    }
  } else {
    if(hollow) {
      register double rq=r*r;
      F   = (R==0)? 0. : exp(-R0/R-R/Rd);
      f   = exp(-R0/r-r/Rd);
      fp  = R0/rq-1./Rd;
      fpp = (fp*fp-2.*R0/(rq*r))*f;
      fp *= f;
    } else {
      F   = exp(-R/Rd);
      f   = exp(-r/Rd);
      fp  =-f/Rd;
      fpp = f/Rd2;
    }
  }

  return fac * ((F-f)*gpp - 2*fp*(g+z*gp)/r - fpp*g);
}

double DiskAnsatz::operator() (double R, double z, double r, double* dP) const
{
  if(S0==0.) {
    if(dP) dP[0]=dP[1]=0.;
    return 0.;
  }
  register double g,f;
  if(dP) {
    register double gp,fp;
    if(thin) {
      g  = abs(z);
      gp = sign(z);
    } else if(isothermal) {
      register double x,ex,sh1;
      x  = abs(z/zd);
      ex = exp(-x);
      sh1= 1.+ex;
      g  = 2*zd*(0.5*x+log(0.5*sh1));
      gp = sign(z)*(1.-ex)/sh1;
    } else {
      register double x,ex;
      x  = abs(z/zd);
      ex = exp(-x);
      g  = zd*(ex-1+x);
      gp = sign(z)*(1.-ex);
    }

    if(hollow && r==0.) f=fp=0.;
    else if(eps) {
      if(hollow) {
	register double rq=r*r,cr=cos(r/Rd),sr=sin(r/Rd);
	f   = exp(-R0/r-r/Rd+eps*cr);
	fp  = (R0/rq-(1.+eps*sr)/Rd)*f;
      } else {
	register double cr=cos(r/Rd),sr=sin(r/Rd);
	f   = exp(-r/Rd+eps*cr);
	fp  = -(1.+eps*sr)*f/Rd;
      }
    } else {
      if(hollow) {
	register double rq=r*r;
	f   = exp(-R0/r-r/Rd);
	fp  = (R0/rq-1./Rd)*f;
      } else {
	f  = exp(-r/Rd);
	fp =-f/Rd;
      }
    }
    dP[0] = fac *  R/r * fp * g;
    dP[1] = fac * (z/r * fp * g + f * gp);
  } else {
    if(thin) {
      g  =abs(z);
    } else if(isothermal) {
      register double x,ex,sh1;
      x  = abs(z/zd);
      ex = exp(-x);
      sh1= 1.+ex;
      g  = 2*zd*(0.5*x+log(0.5*sh1));
    } else {
      register double x,ex;
      x  = abs(z/zd);
      ex = exp(-x);
      g  = zd*(ex-1+x);
    }

    if(hollow && r==0.) f=0.;
    else if(eps) {
      if(hollow) f = exp(-R0/r-r/Rd+eps*cos(r/Rd));
      else       f = exp(-r/Rd+eps*cos(r/Rd));
    } else {
      if(hollow) f = exp(-R0/r-r/Rd);
      else       f = exp(-r/Rd);
    }
  }

  return fac*f*g;
}

double DiskAnsatz::Laplace(double R, double z) const
{
  if(S0==0.) return 0;
  register double r=hypot(R,z), g,gp,gpp, F,f,fp,fpp;
  // deal with the vertical part
  if(thin) {
    g  =abs(z);
    gp =sign(z);
    gpp=0.;
  } else if(isothermal) {
    register double x,sh1;
    x   = abs(z/zd);
    gpp = exp(-x);
    sh1 = 1.+gpp;
    g   = 2*zd*(0.5*x+log(0.5*sh1));
    gp  = sign(z)*(1.-gpp)/sh1;
    gpp/= 0.5*sh1*sh1*zd;
  } else {
    register double x;
    x   = abs(z/zd);
    gpp = exp(-x);
    g   = zd*(gpp-1+x);
    gp  = sign(z)*(1.-gpp);
    gpp/= zd;
  }
  // deal with the radial part
  if(hollow && r==0.) F=f=fp=fpp=0.;
  else if(eps) {
    if(hollow) {
      register double rq=r*r,cr=cos(r/Rd),sr=sin(r/Rd);
      F   = (R==0)? 0. : exp(-R0/R-R/Rd+eps*cos(R/Rd));
      f   = exp(-R0/r-r/Rd+eps*cr);
      fp  = R0/rq-(1.+eps*sr)/Rd;
      fpp = (fp*fp-2.*R0/(rq*r)-eps*cr/Rd2)*f;
      fp *= f;
    } else {
      register double cr=cos(r/Rd),sr=sin(r/Rd);
      F   = exp(-R/Rd+eps*cos(R/Rd));
      f   = exp(-r/Rd+eps*cr);
      fp  = -(1.+eps*sr)/Rd;
      fpp = (fp*fp-eps*cr/Rd2)*f;
      fp *= f;
    }
  } else {
    if(hollow) {
      register double rq=r*r;
      F   = (R==0)? 0. : exp(-R0/R-R/Rd);
      f   = exp(-R0/r-r/Rd);
      fp  = R0/rq-1./Rd;
      fpp = (fp*fp-2.*R0/(rq*r))*f;
      fp *= f;
    } else {
      F   = exp(-R/Rd);
      f   = exp(-r/Rd);
      fp  =-f/Rd;
      fpp = f/Rd2;
    }
  }
  return fac * (f*gpp + 2*fp*(g+z*gp)/r + fpp*g);
}

Frequs DiskAnsatz::kapnuom(double R) const
// returns dPhi/dR, d^2Phi/dR^2, d^2Phi/dz^2 at z=0
{
  if(S0==0.) return Frequs(0.);
  register double er, gpp;
  if(hollow) er = (R==0)? 0. : exp(-R0/R-R/Rd+eps*cos(R/Rd));
  else   er = exp(-R/Rd+eps*cos(R/Rd));
  if(thin) {                          // vertically thin disk
    warning("KapNuOm(Phi) involves delta(z) at z=0");
    gpp=0.;
  } else if(isothermal)               // vertically isothermal disk
    gpp=0.5/zd;
  else                                // vertically exponential disk
    gpp=1./zd;
  Frequs om;
  om[0] = 0.;
  om[1] = fac*er*gpp;
  om[2] = 0.;
  return om;
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class Disks                                                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

void Disks::reset(int N, const DiskPar* p)
{
  if(D) delete[] D;
  nd  = N;
  D   = new DiskAnsatz[nd];
  Dup = D+nd;
  for(int i=0; i!=nd; ++i) D[i].setup(p[i]);
}

Disks::Disks(std::istream& from)
{
  if(!from) error("Trying to construct Disks from a closed std::istream");
  DiskPar P;
  from >> nd;
  SwallowRestofLine(from);
  D   = new DiskAnsatz[nd];
  Dup = D+nd;
  for(DiskAnsatz *p=D; p!=Dup; ++p) {
    from >> P;
    SwallowRestofLine(from);
    p->setup(P);
  }
  nemo_dprintf(4,"Disks: read %d parameters\n",nd);
}

Disks::Disks(const Disks& DS) : nd(DS.nd)
{
  D   = new DiskAnsatz[nd];
  Dup = D+nd;
  for(int i=0; i!=nd; ++i)
    D[i].setup(DS.Parameter(i));
}

Disks::Disks(int N, const DiskPar* p) : nd(N)
{
  D   = new DiskAnsatz[nd];
  Dup = D+nd;
  for(int i=0; i!=nd; ++i) D[i].setup(p[i]);
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class SpheroidDensity                                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

void SpheroidDensity::setup(const SphrPar& d)
{
  rh0 = d(0);
  q   = d(1);
  gam = d(2);
  bet = d(3);
  r0  = d(4);
  rcut= d(5);

  if(rcut<=0.) rci = 0.;
  else     rci = 1./rcut;
  beg = bet-gam;
  qi  = 1./q;
  r0i = 1./r0;
}

double SpheroidDensity::Density(double R, double z) const
{
  register double m = hypot(R,z*qi), m0=m*r0i, rho=rh0;
  if(gam==0.5)   rho /= sqrt(m0);
  else if(gam==1.)    rho /= m0;
  else if(gam==2.)    rho /= m0*m0;
  else if(gam!=0.)    rho /= pow(m0,gam);
  m0 += 1;
  if(beg==1.)    rho /= m0;
  else if(beg==2.)    rho /= m0*m0;
  else if(beg==3.)    rho /= m0*m0*m0;
  else                rho /= pow(m0,beg);
  if(rci)             rho *= exp(-square(m*rci));
  return rho;
}

double SpheroidDensity::mass_integrand(double y) const
{
  if(rci) {
    double y1=1.-y, m=r0*y/y1;
    return pow(y,2.-gam) * pow(y1,bet-4.) * exp(-square(m*rci));
  }
  return pow(y,2-gam) * pow(1.-y,bet-4.);
}

double SpheroidDensity::mass(double m) const
{
  return FPi*q*rh0*cube(r0)*
    qbulir(Adaptor<SpheroidDensity>(this,&SpheroidDensity::mass_integrand),
	   0.,m/(m+r0),1.e-6);
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class Spheroids                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

void Spheroids::reset(int N, const SphrPar* p)
{
  if(S) delete[] S;
  ns  = N;
  S   = new SpheroidDensity[ns];
  Sup = S+ns;
  for(int i=0; i<ns; i++) S[i].setup(p[i]);
}

Spheroids::Spheroids(std::istream& from)
{
  if(!from) error("Trying to construct Spheroids from a closed istream");
  SphrPar P;
  from >> ns;
  SwallowRestofLine(from);
  S   = new SpheroidDensity[ns];
  Sup = S+ns;
  for(SpheroidDensity *p=S; p!=Sup; ++p) {
    from >> P;
    SwallowRestofLine(from);
    p->setup(P);
  }
  nemo_dprintf(4,"Spheroids: read %d parameters\n",ns);
}

Spheroids::Spheroids(const Spheroids& SP) : ns(SP.ns)
{
  S   = new SpheroidDensity[ns];
  Sup = S+ns;
  for(int i=0; i!=ns; ++i) 
    S[i].setup(SP.Parameter(i));
}

Spheroids::Spheroids(int N, const SphrPar* p) : ns(N)
{
  S   = new SpheroidDensity[ns];
  Sup = S+ns;
  for(int i=0; i!=ns; ++i) S[i].setup(p[i]);
}

double Spheroids::beta() const
{
  double b=1.e3;
  for(SpheroidDensity *p=S; p!=Sup; ++p)
    b = min(b, p->outer_power());
  return (b==1.e3)? -1 : b;
}

double Spheroids::gamma() const
{
  double g=0.;
  for(SpheroidDensity *p=S; p!=Sup; ++p)
    g = max(g, p->inner_power());
  return g;
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class Multipole                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace {
  const int N =LMAX/2+1;   // number of multipoles
  const int N2=3*N/2;	   // number of grid point for cos[theta] in [0,1]
  const int N4=5*N/2; 	   // number of points used to integrate over cos[theta]

  typedef tupel<N,double> DBN;
}

void Multipole::AllocArrays()
{
  if(LR) {
    lLc = new double[K[0]];
    d2R = new double[K[0]];
    d2L = new double[K[0]];
  }
  logr = new double[K[0]];
  X[0] = logr;
  X[1] = new double[K[1]];
  Alloc2D(Y[0],K); Alloc2D(Y[1],K); Alloc2D(Y[2],K);
  Alloc2D(Z[0],K); Alloc2D(Z[1],K); Alloc2D(Z[2],K); Alloc2D(Z[3],K);
}

Multipole::Multipole(int         Kk,
                     double      ri,
                     double      ra,
                     double      g,
                     double      b,
                     PotResidual*PR,
                     int         lr)
{
  nemo_dprintf(4,"Multipole::Multipole() ... \n");
  LR   = lr;
  K[0] = Kk;
  K[1] = N2;
  AllocArrays();
  setup(ri,ra,g,b,PR);
  nemo_dprintf(4," done Multipole::Multipole()\n");
}

void Multipole::reset(double      ri,
		      double      ra,
		      double      g,
		      double      b,
		      PotResidual*PR,
		      int         lr)
{
  LR = lr;
  setup(ri,ra,g,b,PR);
}

void Multipole::setup(double      ri,
		      double      ra,
                      double      g,
		      double      b,
                      PotResidual*PR)
{
  Rmin = ri;
  Rmax = ra;
  gamma= g; 
  beta = b;
  lRmin= log(Rmin);
  lRmax= log(Rmax); 
  g2   = 2.-gamma;

  const    DBN    Zero=DBN(0.);
  const    double half=0.5, three=3., sixth=1./6.,
    dlr =(lRmax-lRmin)/double(K[0]-1);
  register int    i,l,k,ll,lli1;
  register double dx,dx2,xl_ll,xh_ll,risq,ril2,dP;
  DBN    A[4],P2l,dP2l,EX;
  //
  // 0  check for inconsistencies in input
  //
  nemo_dprintf(5,"Multipole::setup(): 0\n");
  if(beta>0. && beta<3.) {
    warning("beta in ]0,3[ unsuitable for Multipole expansion; "
	    "we will take beta=3.2\n");
    beta=3.2;
  }
  //
  // 1  compute expansion of the density
  //
  nemo_dprintf(5,"Multipole::setup(): 1\n");
  double
    *ct   = new double[N4],
    *st   = new double[N4],
    *wi   = new double[N4],
    *r    = new double[K[0]];
  DBN
    *W    = new DBN   [N4],
    *rhol = new DBN   [K[0]],
    *rhl2 = new DBN   [K[0]];
  //
  // 1.1 set points and weights for integration over cos(theta)
  //
  nemo_dprintf(6,"Multipole::setup(): 1.1\n");
  GaussLegendre(ct,wi,N4);
  nemo_dprintf(8,"Multipole::setup(): N4=%d\n",N4);
  for(i=0; i<N4; i++) {
    ct[i] = 0.5 * (ct[i]+1.);
    st[i] = sqrt(1.-ct[i]*ct[i]);
    wi[i] = 0.5 * wi[i];
    LegendrePeven(W[i],ct[i]);
    W[i] *= wi[i];
  }
  //
  // 1.2 integrate over cos(theta)
  //
  nemo_dprintf(6,"Multipole::setup(): 1.2\n");
  for(k=0; k<K[0]; k++) {
    logr[k] = k<K[0]-1? lRmin+dlr*k : lRmax;                 // v0.7
    r[k]    = exp(logr[k]);
    rhol[k] = 0.;
    for(i=0; i<N4; i++)
      rhol[k] += W[i] * PR->Residual(r[k],st[i],ct[i]);
  }
  delete[] ct;
  delete[] st;
  delete[] wi;
  delete[] W;
  //
  // 1.3 establish spline in r needed for integration
  //
  nemo_dprintf(6,"Multipole::setup(): 1.3\n");
  spline(r,rhol,K[0],(-gamma/r[0])*rhol[0],Zero,rhl2,0,1);
  //
  // 2. compute potential's expansion
  //
  nemo_dprintf(5,"Multipole::setup(): 2\n");
  DBN
    *P1   = new DBN[K[0]],
    *P2   = new DBN[K[0]],
    *Phil = new DBN[K[0]],
    *dPhl = new DBN[K[0]];
  //
  // 2.1 set P1[k][l] r[k]^(-1-2l) = Int[rho_2l(x,l) x^(2l+2), {x,0,r[k]}]
  //
  //     for r < Rmin we take  rho_2l proportional r^-gamma
  //
  nemo_dprintf(6,"Multipole::setup(): 2.1\n");
  risq  = Rmin*Rmin;
  for(l=0; l<N; l++) {
    P1[0][l] = rhol[0][l] * risq / double(2*l+3-gamma);
    EX[l]    = exp(-(1+2*l)*dlr);
  }
  for(k=0; k<K[0]-1; k++) {
    dx   = r[k+1]-r[k];
    dx2  = dx*dx;
    A[0] = r[k+1]*rhol[k] - r[k]*rhol[k+1] + sixth*r[k]*r[k+1] *
      ( (r[k+1]+dx)*rhl2[k] - (r[k]-dx)*rhl2[k+1] );
    A[1] = rhol[k+1]-rhol[k]
      + sixth * ( (dx2-three*r[k+1]*r[k+1]) * rhl2[k]
		  -(dx2-three*r[k]*r[k])     * rhl2[k+1] );
    A[2] = half  * (r[k+1]*rhl2[k] - r[k]*rhl2[k+1]);
    A[3] = sixth * (rhl2[k+1]-rhl2[k]);
    for(l=0,ll=2; l<N; l++,ll+=2) {
      xl_ll = r[k]*EX(l);
      xh_ll = r[k+1];
      for(i=0,lli1=ll+1,dP=0.; i<4; i++,lli1++) {
	xl_ll*= r[k];
	xh_ll*= r[k+1];
	dP   += A[i](l) * (xh_ll - xl_ll) / lli1;
      }
      P1[k+1][l] = EX(l) * P1[k](l) + dP / dx;
    }
  }
  //
  // 2.2 set P2[k][l] = r[k]^(2l) Int[rho_2l(x,l) x^(1-2l), {x,r[k],Infinity}]
  //
  //     for r > Rmax we take  rho_2l proportional r^-beta if beta>0
  //                                  = 0                  if beta<=0
  //
  nemo_dprintf(6,"Multipole::setup(): 2.2\n");
  if(beta>0.) {
    risq  = Rmax*Rmax;
    for(l=0; l<N; l++) {
      P2[K[0]-1][l] = rhol[K[0]-1][l] * risq / double(beta+2*l-2);
      EX[l] = exp(-2*l*dlr);
    }
  } else {
    P2[K[0]-1] = 0.;
    for(l=0; l<N; l++)
      EX[l] = exp(-2*l*dlr);
  }
  for(k=K[0]-2; k>=0; k--) {
    risq = r[k]*r[k];
    dx   = r[k+1]-r[k];
    dx2  = dx*dx;
    A[0] = r[k+1]*rhol[k] - r[k]*rhol[k+1] + sixth*r[k]*r[k+1] *
      ( (r[k+1]+dx)*rhl2[k] - (r[k]-dx)*rhl2[k+1] );
    A[1] = rhol[k+1]-rhol[k]
      + sixth * ( (dx2-three*r[k+1]*r[k+1]) * rhl2[k]
		  -(dx2-three*r[k]*r[k])     * rhl2[k+1] );
    A[2] = half  * (r[k+1]*rhl2[k] - r[k]*rhl2[k+1]);
    A[3] = sixth * (rhl2[k+1]-rhl2[k]);
    for(l=0,ll=1,ril2=1.; l<N; l++,ll-=2,ril2*=risq) {
      xl_ll = r[k];
      xh_ll = r[k+1]*EX(l);
      for(i=0,lli1=ll+1,dP=0.; i<4; i++,lli1++) {
	xl_ll *= r[k];
	xh_ll *= r[k+1];
	if(lli1) dP += A[i](l) * (xh_ll - xl_ll) / lli1;
	else     dP += A[i](l) * ril2 * dlr;
      }
      P2[k][l] = EX(l) * P2[k+1](l) + dP / dx;
    }
  }
  //
  // 2.3 put together the Phi_2l(r) and dPhi_2l(r)/dlog[r]
  //
  nemo_dprintf(6,"Multipole::setup(): 2.3\n");
  for(k=0; k<K[0]; k++)
    for(l=ll=0; l<N; l++,ll+=2) {
      Phil[k][l] =-P1[k](l) - P2[k](l);                   // Phi_2l
      dPhl[k][l] = (ll+1)*P1[k](l) - ll*P2[k](l);         // dPhi_2l/dlogr
    }
  if(gamma<2)
    Phi0 = Phil[0](0) - dPhl[0](0) / g2;
  delete[] r;
  delete[] rhol;
  delete[] rhl2;
  delete[] P1;
  delete[] P2;
  //
  // 3. establish L_circ(R) on the logarithmic grid
  //
  nemo_dprintf(5,"Multipole::setup(): 3\n");
  if(LR) {
    tg3 = 2./(3.-gamma);
    g3h = 0.5*(3.-gamma);
    lLc = new double[K[0]];
    d2R = new double[K[0]];
    d2L = new double[K[0]];
    LegendrePeven(P2l,0.);
    for(k=0; k<K[0]; k++)
      lLc[k] = 0.5 * ( 2*logr[k] + log(dPhl[k]*P2l) );
    spline(lLc,logr,K[0],tg3,2.,d2R,0,0);
    spline(logr,lLc,K[0],g3h,.5,d2L,0,0);
    lzmin = lLc[0];
    lzmax = lLc[K[0]-1];
  }
  //
  // 4.  Put potential and its derivatives on a 2D grid in log[r] & cos[theta]
  //
  nemo_dprintf(5,"Multipole::setup(): 4\n");
  //
  // 4.1 set linear grid in theta
  //
  nemo_dprintf(6,"Multipole::setup(): 4.1\n");
  for(i=0; i<N2; i++) 
    X[1][i] = double(i) / double(N2-1);
  //
  // 4.2 set dPhi/dlogr & dPhi/dcos[theta] 
  //
  nemo_dprintf(6,"Multipole::setup(): 4.2\n");
  for(i=0; i<N2; i++) {
    dLegendrePeven(P2l,dP2l,X[1][i]);
    for(k=0; k<K[0]; k++) {
      Y[0][k][i] = Phil[k] * P2l;			// Phi
      Y[1][k][i] = dPhl[k] * P2l;			// d Phi / d logR
      Y[2][k][i] = Phil[k] * dP2l;		// d Phi / d cos(theta)
    }
  }
  delete[] Phil;
  delete[] dPhl;
  //
  // 4.3 establish 2D Pspline of Phi in log[r] & cos[theta]
  //
  nemo_dprintf(6,"Multipole::setup(): 4.3\n");
  Pspline2D(X,Y,K,Z);
  nemo_dprintf(5,"Multipole::setup(): done\n");
}

Multipole::~Multipole()
{
  if(LR) {
    delete[] lLc;
    delete[] d2R;
    delete[] d2L;
  }
  delete[] X[0]; delete[] X[1];
  Free2D(Y[0]); Free2D(Y[1]); Free2D(Y[2]);
  Free2D(Z[0]); Free2D(Z[1]); Free2D(Z[2]); Free2D(Z[3]);
}

double Multipole::operator() (double r, double ct, double st, double* dP) const
{
  double Xi[2];
  double lr=log(r), Phi;
  Xi[0] = min(lRmax,max(lRmin,lr));
  Xi[1] = abs(ct);
  Phi   = Psplev2D(X,Y,Z,K,Xi,dP);
  if(dP) dP[1]*= sign(ct);
  if(lr < lRmin) {
    if(g2>0.) {
      Phi = (Phi-Phi0)*exp(g2*(lr-Xi[0]));
      if(dP) dP[0] = g2*Phi;
      Phi+= Phi0;
    } else if(g2==0.) {
      if(dP) dP[0] = Phi/lRmin;
      Phi*= lr/lRmin;
    } else {
      Phi*= exp(g2*(lr-Xi[0]));
      if(dP) dP[0] = g2*Phi;
    }
  } else if(lr > lRmax) {
    Phi *= Rmax/r;
    if(dP) dP[0] =-Phi;
  }
  if(dP) {
    register double temp;
    dP[0]/= r;
    dP[1]*=-st/r;
    temp  = ct*dP[0] - st*dP[1];
    dP[0] = st*dP[0] + ct*dP[1];
    dP[1] = temp;
  }
  return Phi;
}

double Multipole::vcsquare(double R, double &dvcqdR) const
{
  const int n2[2]={2,2};
  double Xi[2], dP[2], **d2P;
  double lr=log(R), Phi;
  Alloc2D(d2P,n2);
  Xi[0] = min(lRmax,max(lRmin,lr));
  Xi[1] = 0.;
  Phi   = Psplev2D(X,Y,Z,K,Xi,dP,d2P);
  if(lr < lRmin) {
    if(g2>0.) {
      dP[0]     = g2*(Phi-Phi0)*exp(g2*(lr-Xi[0]));
      d2P[0][0] = g2*dP[0];
    } else if(g2==0.) {
      dP[0]     = Phi/lRmin;
      d2P[0][0] = 0.;
    } else {
      dP[0]     = g2*Phi*exp(g2*(lr-Xi[0]));
      d2P[0][0] = g2*dP[0];
    }
  } else if(lr > lRmax) {
    dP[0]     =-Phi*Rmax/R;
    d2P[0][0] =-dP[0];
  }
  dvcqdR = d2P[0][0] / R;     // dvc^2 / dR = (1/R) d^2 Phi/ d(lnR)^2
  Free2D(d2P);
  return dP[0];               // vc^2       = dPhi/dlnR
}

double Multipole::Laplace(double r, double ct) const
{
  const int m[2]={2,2};
  double lr=log(r), Phi, Lap;
  double Xi[2], *dP, **d2P;
  Alloc1D(dP,2);
  Alloc2D(d2P,m);
  Xi[0] = min(lRmax,max(lRmin,lr));
  Xi[1] = abs(ct);
  Phi   = Psplev2D(X,Y,Z,K,Xi,dP,d2P);
  dP[1]*=sign(ct);
  if(lr < lRmin) {
    if(g2>0.) {
      dP[0]     = g2*(Phi-Phi0)*exp(g2*(lr-Xi[0]));
      d2P[0][0] = g2*dP[0];
    } else if(g2==0.) {
      dP[0]     = Phi/lRmin;
      d2P[0][0] = 0.;
    } else {
      dP[0]     = g2*Phi*exp(g2*(lr-Xi[0]));
      d2P[0][0] = g2*dP[0];
    }
  } else if(lr > lRmax) {
    dP[0]     =-Phi*Rmax/r;
    d2P[0][0] =-dP[0];
  }
  Lap = ( dP[0]+d2P[0][0] + (1.-ct*ct)*d2P[1][1]-2.*ct*dP[1] ) / (r*r);
  Free1D(dP);
  Free2D(d2P);
  return Lap;
}

Frequs Multipole::kapnuom(double R) const
// returns dPhi/dR, d^2Phi/dR^2/ d^2Phi/dz^2
{
  const int m[2]={2,2};
  double lr=log(R), Rq=R*R, Phi;
  double Xi[2], *dP, **d2P;
  Alloc1D(dP,2);
  Alloc2D(d2P,m);
  Xi[0] = min(lRmax,max(lRmin,lr));
  Xi[1] = 0.;
  Phi   = Psplev2D(X,Y,Z,K,Xi,dP,d2P);
  if(lr < lRmin) {
    if(g2>0.) {
      dP[0]     = g2*(Phi-Phi0)*exp(g2*(lr-Xi[0]));
      d2P[0][0] = g2*dP[0];
    } else if(g2==0.) {
      dP[0]     = Phi/lRmin;
      d2P[0][0] = 0.;
    } else {
      dP[0]     = g2*Phi*exp(g2*(lr-Xi[0]));
      d2P[0][0] = g2*dP[0];
    }
  } else if(lr > lRmax) {
    dP[0]     =-Phi*Rmax/R;
    d2P[0][0] =-dP[0];
  }
  Frequs om;
  om[0] = (d2P[0][0]-dP[0]) / Rq;
  om[1] = (dP[0] + d2P[1][1]) / Rq;
  om[2] = dP[0] / R;
  Free1D(dP);
  Free2D(d2P);
  return om; 
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class GalaxyPotential                                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

double GalaxyPotential::operator() (double R, double z) const
{
  register double r  =hypot(R,z), pot=M(r,z/r,R/r);
  for(register DiskAnsatz *p=D; p<Dup; p++) pot+= (*p)(R,z,r);
  return pot;
}

double GalaxyPotential::operator() (double R, double z,
                                    double&dR, double&dz) const
{
  double d[2];
  register double r  =hypot(R,z), pot=M(r,z/r,R/r,d);
  dR = d[0]; dz = d[1];
  for(register DiskAnsatz *p=D; p<Dup; p++) {
    pot += (*p)(R,z,r,d);
    dR  += d[0];
    dz  += d[1];
  }
  return pot;
}

void GalaxyPotential::OortConstants(double R, double &A, double &B) const
{
  double vc, dvc;
  vc  = sqrt(M.vcsquare(R,dvc));
  dvc/= 2.*vc;
  A   = 0.5 * (vc/R - dvc);
  B   =-0.5 * (vc/R + dvc);
}

double GalaxyPotential::Laplace(double R, double z) const
{
  register double r=hypot(R,z), L=M.Laplace(r,z/r);
  register DiskAnsatz *p=D;
  for(; p<Dup; p++) L += p->Laplace(R,z);
  return L;
}

Frequs GalaxyPotential::KapNuOm  (double R) const
{
  Frequs Om = M.kapnuom(R);
  for(register DiskAnsatz *p=D; p<Dup; p++) Om += p->kapnuom(R);
  Om[2]/= R;
  Om[0]+= 3*Om(2);
  Om.apply(&sqrt);
  return Om;                  // omega, kappa, nu
}

// end of file GalPot.cc ///////////////////////////////////////////////////////
