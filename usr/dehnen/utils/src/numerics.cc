// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file    src/numerics.cc
///
/// \author  Walter Dehnen
///
/// \date    1994-2008,2010
///
/// \todo    add doxygen documentation
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 1994-2008  Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 675
// Mass Ave, Cambridge, MA 02139, USA.
//
////////////////////////////////////////////////////////////////////////////////
#include <numerics.h>
#include <WDMath.h>
#ifndef WDutils_included_cstring
#  include <cstring>
#  define WDutils_included_cstring
#endif
#ifndef WDutils_included_limits
#  include <limits>
#  define WDutils_included_limits
#endif

#ifdef __INTEL_COMPILER
#pragma warning (disable:383) /* reference to temporary used */
#pragma warning (disable:424) /* extra ";" ignored */
#pragma warning (disable:981) /* operands are evaluated in unspecified order */
#endif

////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////
double WDutils::qbulir(double(*func)(double),
		       double  a,
		       double  b, 
		       double  eps,
		       double *erro,
		       bool    abort,
		       int     mx)
{
  // get actual computing precision and adjust accuracy
  const double eta = std::numeric_limits<double>::epsilon();

  if(eps<eta) eps=eta;
  // get integration interval and return zero if it vanishes
  double ba=b-a;
  if(ba<eta) return 0.;
  // initialise some variables
  bool   bo,bu=0,odd=1;
  int    i,m,n=2,nn=3,mr;
  double d[7],dt[7];
  double ddt,hm,nt,err(0),t,ta,tab=(0),v=(0),w(0),sm(0),gr(0),t1(0),
    t2 =0.5*(func(a)+func(b)),
    t2a=t2,
    tb =abs(t2a),
    c  =t2*ba;
  dt[0]=c;
  // iterate over the refinements
  for(m=1; m<=mx; ++m) {
    bo = m>=7;
    hm = ba/n;
    if(odd) {
      for(i=1; i<=n; i+=2) {
	w = func(a+i*hm);
	t2 += w;
	tb += abs(w);
      }
      nt   = t2;
      tab  = tb * abs(hm);
      d[1] = 1.77777777777777777778; // 16/9
      d[3] = 7.11111111111111111111; // 64/9
      d[5] = 28.4444444444444444444; // 256.9
    } else {
      for(i=1; i<=n; i+=6) {
	w = i*hm;
	t1+= func(a+w) + func(b-w);
      }
      nt   = t1+t2a;
      t2a  = t2;
      d[1] = 2.25; // 9/4
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
      double d1  = d[i]*ddt;
      double den = d1-nt;
//    double e   = nt-ddt;
      if(abs(den)>eta) {
	double e = (nt-ddt)/den;
//	e /= den;
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
    v  = t-ta;
    t += v;
    err= abs(v);
    if(ta<t) {
      double d1 = ta;
      ta = t;
      t  = d1;
    }
    bo  = bo || (ta<gr && t>sm);
    if(bu && bo && err < eps*tab*w) break;
    gr  = ta;
    sm  = t;
    odd = !odd;
    i   = n;
    n   = nn;
    nn  = i+i;
    bu  = bo;
    d[2]= 4.;
    d[4]= 16.;
  }
  v = tab*eta;
  if(err<v) err = v;
  if(erro) *erro = err/(tab*w);
  if(m==mx) {
    if(abort) WDutils_Error("in qbulir(): max number of iterations exceeded");
    else    WDutils_Warning("in qbulir(): max number of iterations exceeded");
  }
  return c;
}
////////////////////////////////////////////////////////////////////////////////
// Gauss-Legendre integration: points & weights                                 
//------------------------------------------------------------------------------
void WDutils::GaussLegendre(double *x, double *w, const unsigned n)
{
  const double eps=std::numeric_limits<double>::epsilon();
  const int m=(n+1)/2;
  for(int i=0; i!=m; ++i) {
    double z1,pp,z=std::cos(Pi*(i+0.75)/(n+0.5));
    do {
      double p1 = 1.0;
      double p2 = 0.0;
      for(unsigned j=0; j!=n; ++j) {
	double p3 = p2;
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
////////////////////////////////////////////////////////////////////////////////
// Householder reduction of real symmetric matrix                               
//------------------------------------------------------------------------------
namespace WDutils {
  template<bool EIGENVECTORS, typename X>
  void HouseholderReduction(int n, X**a, X*d, X*e)
  { 
    const X zero(0), one(1);
    for(int i=n-1; i; --i) {
      X h = zero;
      if(i>1) {
	X sc = zero;
	for(int k=0; k!=i; ++k)
	  sc += abs(a[i][k]);
	if(iszero(sc))
// 	if(sc == zero)
	  e[i] = a[i][i-1];
	else {
	  X in = one/sc;
	  for(int k=0; k!=i; ++k) {
	    a[i][k] *= in;
	    h += square(a[i][k]);
	  }
	  X f = a[i][i-1];
	  X g = f>=zero? -sqrt(h) : sqrt(h);
	  e[i] = sc*g;
	  h -= f*g;
	  a[i][i-1] = f-g;
	  f = zero;
	  in = one/h;
	  for(int j=0; j!=i; ++j) {
	    if(EIGENVECTORS)
	      a[j][i] = a[i][j] * in;
	    g = zero;
	    for(int k=0; k<=j; ++k)
	      g += a[j][k]*a[i][k];
	    for(int k=j+1; k!=i; ++k)
	      f += a[k][j]*a[i][k];
	    e[j] = g*in;
	    f += e[j]*a[i][j];
	  }
	  X hh = f/(h+h);
	  for(int j=0; j!=i; ++j) {
	    f = a[i][j];
	    e[j] = g = e[j]-hh*f;
	    for(int k=0; k<=j; ++k)
	      a[j][k] -= f*e[k]+g*a[i][k];
	  }
	}
      } else
	e[i] = a[i][i];
      d[i] = h;
    }
    if(EIGENVECTORS)
      d[0] = zero;
    e[0] = zero;
    if(EIGENVECTORS)
      for(int i=0; i!=n; ++i) {
	if(!iszero(d[i])) {
// 	if(d[i]) {
	  for(int j=0; j!=i; ++j) {
	    X g = zero;
	    for(int k=0; k!=i; ++k)
	      g += a[i][k]*a[k][j];
	    for(int k=0; k!=i; ++k)
	      a[k][j] -= g*a[k][i];
	  }
	}
	d[i] = a[i][i];
	a[i][i] = one;
	for(int j=0; j!=i; ++j)
	  a[j][i] = a[i][j] = zero;
      }
    else
      for(int i=0; i!=n; ++i)
	d[i] = a[i][i];
  }
  //----------------------------------------------------------------------------
  template void HouseholderReduction<1,float >(int, float **, float *, float *);
  template void HouseholderReduction<1,double>(int, double**, double*, double*);
  template void HouseholderReduction<0,float >(int, float **, float *, float *);
  template void HouseholderReduction<0,double>(int, double**, double*, double*);
  //----------------------------------------------------------------------------
  template<typename X> void EigenSystemTridiagonal(int n, X*d, X*e, X**z)
  {
    const X zero(0), one(1);
    for(int i=1; i!=n; ++i)
      e[i-1] = e[i];
    e[n-1] = zero;
    for(int l=0; l!=n; ++l) {
      int iter=0, m;
      do {
	for(m=l; m!=n-1; ++m) {
	  X dd = abs(d[m])+abs(d[m+1]);
	  if(insignificant(e[m],dd)) break;
//	  if(abs(e[m])+dd == dd) break;
	}
	if(m != l) {
	  if(iter++ == 30) WDutils_Error("in EigenSystemTridiagonal(): "
					 "max number of iterations exceeded");
	  X g = (d[l+1]-d[l])/twice(e[l]);
	  X r = hypot(g,one);
	  g   = d[m]-d[l]+e[l]/(g+sign(r,g));
	  X s = one;
	  X c = one;
	  X p = zero;
	  int i=m-2;
	  for(; i>=0; --i) {
	    X f = s*e[i];
	    X b = c*e[i];
	    r = hypot(f,g);
	    e[i+1] = r;
	    if(iszero(r)) {
// 	    if(r==zero) {
	      d[i+1] -= p;
	      e[m] = zero;
	      break;
	    }
	    s = f/r;
	    c = g/r;
	    g = d[i+1]-p;
	    r = (d[i]-g)*s+twice(c*b);
	    p = s*r;
	    d[i+1] = g+p;
	    g = c*r-b;
	    for(int k=0; k!=n; ++k) {
	      f = z[k][i+1];
	      z[k][i+1] = s*z[k][i]+c*f;
	      z[k][i]   = c*z[k][i]-s*f;
	    }
	  }
	  if(iszero(r) && i>=0) continue;
// 	  if(r==zero && i>=0) continue;
	}
      } while(m!=l);
    }
  }
  //----------------------------------------------------------------------------
  template void EigenSystemTridiagonal<float >(int, float *, float *, float **);
  template void EigenSystemTridiagonal<double>(int, double*, double*, double**);
  //----------------------------------------------------------------------------
  template<typename X> void EigenValuesTridiagonal(int n, X*d, X*e)
  {
    const X zero(0), one(1);
    for(int i=1; i!=n; ++i)
      e[i-1] = e[i];
    e[n-1] = zero;
    for(int l=0; l!=n; ++l) {
      int iter=0, m;
      do {
	for(m=l; m!=n-1; ++m) {
	  X dd = abs(d[m])+abs(d[m+1]);
	  if(insignificant(e[m],dd)) break;
// 	  if(abs(e[m])+dd == dd) break;
	}
	if(m != l) {
	  if(iter++ == 30) WDutils_Error("in EigenValuesTridiagonal(): "
					 "max number of iterations exceeded");
	  X g = (d[l+1]-d[l])/twice(e[l]);
	  X r = hypot(g,one);
	  g   = d[m]-d[l]+e[l]/(g+sign(r,g));
	  X s = one;
	  X c = one;
	  X p = zero;
	  int i=m-2;
	  for(; i>=0; --i) {
	    X f = s*e[i];
	    X b = c*e[i];
	    r = hypot(f,g);
	    e[i+1] = r;
	    if(iszero(r)) {
// 	    if(r==zero) {
	      d[i+1] -= p;
	      e[m] = zero;
	      break;
	    }
	    s = f/r;
	    c = g/r;
	    g = d[i+1]-p;
	    r = (d[i]-g)*s+twice(c*b);
	    p = s*r;
	    d[i+1] = g+p;
	    g = c*r-b;
	  }
	  if(iszero(r) && i>=0) continue;
// 	  if(r==zero && i>=0) continue;
	}
      } while(m!=l);
    }
  }
  //----------------------------------------------------------------------------
  template void EigenValuesTridiagonal<float >(int, float *, float *);
  template void EigenValuesTridiagonal<double>(int, double*, double*);
} // namespace WDutils
////////////////////////////////////////////////////////////////////////////////
//
// finding percentiles in rank or cumulative weight
//
////////////////////////////////////////////////////////////////////////////////
namespace {
  using namespace WDutils;
  //////////////////////////////////////////////////////////////////////////////
  /// auxiliary for FindPercentile
  template<typename scalar>
  class Ranker {
    friend class FindPercentile<scalar>;
    /// represents a single value
    struct point {
      scalar   X;       ///< value
      scalar   W;       ///< weight
      unsigned I;       ///< index
    };
    /// represents a range of points
    struct range {
      unsigned N;       ///< number of points in range
      unsigned R;       ///< rank of first point
      scalar   W;       ///< cumulative weight at first point
      range   *S;       ///< pter to left sub-range
      range() {}
      range(int n) : N(n), R(0), W(0), S(0) {}
    };
    scalar             SumW;        ///< total weight
    point             *P;           ///< table of points
    range              Root;        ///< root range
    block_alloc<range> RangeAlloc;  ///< allocator for ranges
    /// swap two points
    static void swap(point*a, point*b)
    {
      point p;
      memcpy(&p, a, sizeof(point));
      memcpy( a, b, sizeof(point));
      memcpy( b,&p, sizeof(point));
    }
    /// split a list of points
    /// \param[in]  p list of points to split
    /// \param[in]  n size of list
    /// \param[in]  s split point
    /// \param[out] w total weight of left part
    /// \return     number of points in left part
    /// \note If s is in the range of positions, the split will always generate
    ///       non-empty sub-ranges.
    static unsigned splitlist(point*p, unsigned np, scalar s, scalar&w)
    {
      point*l,*u,*n=p+np;
      w = 0;
      // find l first element which is not smaller than S
      for(l=p;   l!=n && l->X<s; ++l)            w += l->W;
      if(l==n) return np;
      // find u first element beyond l which is not greater than S
      for(u=l+1; u!=n && u->X>s; ++u) {}
      while(u!=n) {
	swap(l,u);                               w += l->W;
	for(++l;            l!=n && l->X<s; ++l) w += l->W;
	for(u=max(u+1,l+1); u!=n && u->X>s; ++u) {}
      }
      return l - p;
    }
    /// split a range.
    /// we split at the element in mid-index.
    void split(range*A) WDutils_THROWING;
    /// ctor: create points and set root
    Ranker(const scalar*, unsigned, const scalar*, unsigned) WDutils_THROWING;
    /// ctor: create points and set root
    Ranker(unsigned, scalar(*)(unsigned), unsigned) WDutils_THROWING;
    /// ctor: create points and set root
    Ranker(unsigned, void(*)(unsigned,scalar&,scalar&),unsigned)
    WDutils_THROWING;
  public:
    /// dtor: delete points
    ~Ranker() WDutils_THROWING { WDutils_DEL_A(P); }
    /// is a range pointer valid?
    bool is_valid(const range*r)
    {
      return RangeAlloc.is_element(r);
    }
    /// rank points and ranges at given rank
    /// \param[in] r rank
    /// \return range of size 1 with position of given rank
    /// The sort tree is refined such that at rank R the range has size 1 and
    /// the point or rank R has internal index R.
    const range*RankR(unsigned r) WDutils_THROWING
    {
      range*A=&Root;
      if(r>=A->N) WDutils_THROW("FindPercentile<%s>::FindRank: "
				"r=%d >= N=%d\n",nameof(scalar),r,A->N);
      while(A->N>1) {
	if(A->S==0) split(A);
	A = r<A->S[1].R? A->S : A->S+1;
      }
      return A;
    }
    /// rank points and ranges at given cumulative weight
    /// \param[in] w cumulative weight
    /// \return range of size 1 with cumulative weight just below W
    /// The sort tree is refined such that at cumulative weight W the range has
    /// size 1 and the point rank and index match.
    const range*RankW(scalar w) WDutils_THROWING
    {
      if(w>SumW) WDutils_THROW("FindPercentile<%s>::FindCumulativeWeight: "
			       "w=%f >= Wtot=%f\n",nameof(scalar),w,SumW);
      range*A=&Root;
      while(A->N>1) {
	if(A->S==0) split(A);
	A = w<A->S[1].W? A->S : A->S+1;
      }
      return A;
    }
  };// class FindPercentile<>
} // namespace {
namespace WDutils {
  WDutils_TRAITS(Ranker<float>::point,"Ranker<float>::point");
  WDutils_TRAITS(Ranker<float>::range,"Ranker<float>::range");
  WDutils_TRAITS(Ranker<double>::point,"Ranker<double>::point");
  WDutils_TRAITS(Ranker<double>::range,"Ranker<double>::range");
}
//
namespace {
  // Ranker::Ranker
  template<typename T>
  Ranker<T>::Ranker(const T*a, unsigned n, const T*w, unsigned k)
    WDutils_THROWING :
    SumW       ( 0 ),
    P          ( WDutils_NEW(point,n) ),
    Root       ( n ),
    RangeAlloc ( k? 4*k*int(1+log(double(n))) : 10*int(1+log(double(n))) )
  {
    for(unsigned i=0; i!=n; ++i) {
      P[i].X = a[i];
      P[i].I = i;
      P[i].W = w? w[i] : T(1);
      if(P[i].W<=0)
	WDutils_THROW("FindPercentile: weight #%d = %f <= 0\n",i,P[i].W);
      SumW += P[i].W;
    }
  }
  //
  template<typename T>
  Ranker<T>::Ranker(unsigned n, T(*F)(unsigned), unsigned k)
    WDutils_THROWING :
    SumW       ( 0 ),
    P          ( WDutils_NEW(point,n) ),
    Root       ( n ),
    RangeAlloc ( k? 4*k*int(1+log(double(n))) : 10*int(1+log(double(n))) )
  {
    for(unsigned i=0; i!=n; ++i) {
      P[i].X = F(i);
      P[i].I = i;
      P[i].W = 1;
      SumW += P[i].W;
    }
  }
  //
  template<typename T>
  Ranker<T>::Ranker(unsigned n, void(*F)(unsigned, T&, T&), unsigned k)
    WDutils_THROWING :
    SumW       ( 0 ),
    P          ( WDutils_NEW(point,n) ),
    Root       ( n ),
    RangeAlloc ( k? 4*k*int(1+log(double(n))) : 10*int(1+log(double(n))) )
  {
    for(unsigned i=0; i!=n; ++i) {
      P[i].I = i;
      F(i,P[i].X,P[i].W);
      if(P[i].W<=0)
	WDutils_THROW("FindPercentile: weight #%d = %f <= 0\n",i,P[i].W);
      SumW += P[i].W;
    }
  }
  // Ranker::split
  template<typename T>
  void Ranker<T>::split(range*A) WDutils_THROWING
  {
    if(A->N >= 2) {
      unsigned L;
      T        W;
      if(A->N > 2) {
	L = splitlist(P+A->R, A->N, P[A->R+A->N/2].X, W);
      } else {
	if(P[A->R].X > P[A->R+1].X) swap(P+A->R, P+A->R+1);
	W    = P[A->R].W;
	L    = 1;
      }
      A->S = RangeAlloc.new_elements(2);
      A->S[0].R   = A->R;    A->S[1].R = A->R+L;
      A->S[0].N   = L;       A->S[1].N = A->N-L;
      A->S[0].W   = A->W;    A->S[1].W = A->W+W;
      A->S[0].S   = 0;       A->S[1].S = 0;
    } else 
      WDutils_THROW("FindPercentile: cannot split range with N=%d<2\n",A->N);
  }
} // namespace {
////////////////////////////////////////////////////////////////////////////////
namespace WDutils {
#define RANKER         static_cast<Ranker<T>*>(DATA)
#define Range(HANDLE)  static_cast<const typename Ranker<T>::range*>(HANDLE)
  //
  template<typename T>
  void FindPercentile<T>::setup(const T*X, unsigned N,
				const T*W, unsigned K)
    WDutils_THROWING
  {
    if(DATA) WDutils_THROW("FindPercentile<%s>::setup(): DATA=%p != 0\n",
			   nameof(T),DATA);
    DATA = new Ranker<T>(X,N,W,K);
  }
  //
  template<typename T>
  void FindPercentile<T>::setup(unsigned N, T(*X)(unsigned),
				unsigned K)
    WDutils_THROWING
  {
    if(DATA) WDutils_THROW("FindPercentile<%s>::setup(): DATA=%p != 0\n",
			   nameof(T),DATA);
    DATA = new Ranker<T>(N,X,K);
  }
  //
  template<typename T>
  void FindPercentile<T>::setup(unsigned N,
				void(*F)(unsigned, T&, T&),
				unsigned K)
    WDutils_THROWING
  {
    if(DATA) WDutils_THROW("FindPercentile<%s>::setup(): DATA=%p != 0\n",
			   nameof(T),DATA);
    DATA = new Ranker<T>(N,F,K);
  }
  //
  template<typename T>
  T FindPercentile<T>::TotalWeight() const
  {
    return RANKER->SumW;
  }
  //
  template<typename T>
  unsigned FindPercentile<T>::TotalNumber() const
  {
    return RANKER->Root.N;
  }
  //
  template<typename T> typename FindPercentile<T>::
  handle FindPercentile<T>::FindRank(unsigned r) const WDutils_THROWING
  {
    return RANKER->RankR(r);
  }
  //
  template<typename T> typename FindPercentile<T>::
  handle FindPercentile<T>::FindCumulativeWeight(T w) const WDutils_THROWING
  {
    return RANKER->RankW(w);
  }
  //
  template<typename T>
  unsigned FindPercentile<T>::Index(handle h, bool c) const WDutils_THROWING
  {
    if(c && (!RANKER->is_valid(Range(h)) || Range(h)->N!=1) )
      WDutils_THROW("FindPercentile<%s>::Index(): invalid handle\n",nameof(T));
    return RANKER->P[Range(h)->R].I;
  }
  //
  template<typename T>
  T FindPercentile<T>::Weight(handle h, bool c) const WDutils_THROWING
  {
    if(c && (!RANKER->is_valid(Range(h)) || Range(h)->N!=1) )
      WDutils_THROW("FindPercentile<%s>::Weight(): invalid handle\n",nameof(T));
    return RANKER->P[Range(h)->R].W;
  }
  //
  template<typename T>
  T FindPercentile<T>::Position(handle h, bool c) const WDutils_THROWING
  {
    if(c && (!RANKER->is_valid(Range(h)) || Range(h)->N!=1) )
      WDutils_THROW("FindPercentile<%s>::Position(): invalid handle\n",
		    nameof(T));
    return RANKER->P[Range(h)->R].X;
  }
  //
  template<typename T>
  unsigned FindPercentile<T>::Rank(handle h, bool c) const WDutils_THROWING
  {
    if(c && (!RANKER->is_valid(Range(h)) || Range(h)->N!=1) )
      WDutils_THROW("FindPercentile<%s>::Rank(): invalid handle\n",nameof(T));
    return Range(h)->R;
  }
  //
  template<typename T>
  T FindPercentile<T>::CumulativeWeight(handle h, bool c) const WDutils_THROWING
  {
    if(c && (!RANKER->is_valid(Range(h)) || Range(h)->N!=1) )
      WDutils_THROW("FindPercentile<%s>::Weight(): invalid handle\n",nameof(T));
    return Range(h)->W;
  }
  //
  template<typename T>
  FindPercentile<T>::~FindPercentile()
  {
    WDutils_DEL_O(static_cast<Ranker<T>*>(DATA));
    DATA = 0;
  }
#undef Range
#undef RANKER
  //
  template class FindPercentile<float>;
  template class FindPercentile<double>;
} // namespace WDutils
////////////////////////////////////////////////////////////////////////////////
