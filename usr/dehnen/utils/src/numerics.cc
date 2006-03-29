// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    src/exception.cc                                                   
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    1994-2005                                                          
///                                                                             
/// \todo    add doxygen documentation                                          
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 1994-2005  Walter Dehnen                                       
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
#include <numerics.h>
#include <Pi.h>
#include <cmath>

//------------------------------------------------------------------------------
// Burlisch-Stoer integration of 1D real integrals                              
//------------------------------------------------------------------------------
// Quadrature program using the Bulirsch sequence and rational extrapolation.   
// The algorithm is puplished in Bulirsch & Stoer, Num. Math. 9, 271-278 (1967),
// where a routine in ALGOL is given. This routine is a straightforward         
// translation into C++.                                                        
// CAUTION:                                                                     
// Do not use this routine for integrating low order polynomials (up to fourth  
// order) or periodic functions with period equal to the interval of integration
// or linear combinations of both.                                              
// INPUT:  func   pointer to function to be integrated.                         
//         a,b    lower and upper boundaries of the integration interval;       
//         eps    desired relativ accuracy;                                     
// OUTPUT: return approximated value for the integral;                          
//         err    actual relative error of the return value.                    
// -----------------------------------------------------------------------------
double WDutils::qbulir(double(*func)(double),
		       double  a,
		       double  b, 
		       double  eps_,
		       double *erro,
		       bool    abort,
		       int     mx)
{
  register double ba=b-a;
  if(ba==0.) return 0.;
  register int    i,n=2,nn=3,m,mr, bo,bu=0,odd=1;
  register double c,d1,ddt,den,e,eps,eta=1.e-7,gr,hm,nt,err,
                  sm,t,t1,t2,t2a,ta,tab=0.,tb,v=0.,w;
  double          d[7],dt[7];

  while(eta+1. != 1.) eta *=0.5;
  eta  *=2.;                       // eta = actual computing accuracy
  eps   = max(eps_,eta);
  sm    = 0.;
  gr    = 0.;
  t1    = 0.;
  t2    = 0.5*((*func)(a)+(*func)(b));
  t2a   = t2;
  tb    = abs(t2a);
  c     = t2*ba;
  dt[0] = c;
  for(m=0; m!=mx; ++m) {           // iterate over the refinements
    bo = (m>=7);
    hm = ba/n;
    if(odd) {
      for(i=1;i<=n;i+=2) {
	w  = (*func)(a+i*hm);
	t2+= w;
	tb+= abs(w);
      }
      nt  = t2;
      tab = tb * abs(hm);
      d[1]=16./9.;
      d[3]=64./9.;
      d[5]=256./9.;
    } else {
      for(i=1;i<=n;i+=6) {
	w  = i*hm;
	t1+= (*func)(a+w) + (*func)(b-w);
      }
      nt  = t1+t2a;
      t2a =t2;
      d[1]=9./4.;
      d[3]=9.;
      d[5]=36.;
    }
    ddt  =dt[0];
    t    =nt*hm;
    dt[0]=t;
    nt   =dt[0];
    if(bo) {
      mr  =6;
      d[6]=64.;
      w   =144.;
    } else {
      mr  =m;
      d[m]=n*n;
      w   =d[m];
    }
    for(i=1;i<=mr;i++) {
      d1 =d[i]*ddt;
      den=d1-nt;
      e  =nt-ddt;
      if(den != 0.) {
	e /= den;
	v  = nt*e;
	nt = d1*e;
	t += v;
      } else {
	nt = 0.;
	v  = 0.;
      }
      ddt  = dt[i];
      dt[i]= v;
    }
    ta = c;
    c  = t;
    if(!bo) t -= v;
    v  = t-ta;
    t += v;
    err= abs(v);
    if(ta<t) {
      d1 = ta;
      ta = t;
      t  = d1;
    }
    bo = bo || (ta<gr && t>sm);
    if(bu && bo && err < eps*tab*w) break;
    gr = ta;
    sm = t;
    odd= !odd;
    i  = n;
    n  = nn;
    nn = i+i;
    bu = bo;
    d[2]=4.;
    d[4]=16.;
  }
  v = tab*eta;
  if(err<v) err = v;
  if(erro) *erro = err/(tab*w);
  if(m==mx) {
    if(abort) WDutils_ErrorF("max number of iterations exceeded","qbulir()");
    else      WDutils_WarningF("max number of iterations exceeded","qbulir()");
  }
  return c;
}
//------------------------------------------------------------------------------
// Gauss-Legendre integration: points & weights                                 
//------------------------------------------------------------------------------
void WDutils::GaussLegendre(double *x, double *w, const unsigned n)
{
  register double eps=1.e-10;
  for(register double ep1=1.0+eps; 1.!=ep1; eps*=0.5, ep1=1.0+eps);
//   register double eps;
//   for(eps=1.e-10; (eps+1.)!=1.; eps*=0.5);
  eps *=2.;
  register int    i,m=(n+1)/2;
  register double z1,z,pp,p3,p2,p1;
  for (i=0;i<m;i++) {
    z=std::cos(Pi*(i+0.75)/(n+0.5));
    do {
      p1 = 1.0;
      p2 = 0.0;
      for(register unsigned j=0; j!=n; ++j) {
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
//------------------------------------------------------------------------------
