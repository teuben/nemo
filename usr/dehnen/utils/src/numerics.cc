// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    src/numerics.cc                                                    
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    1994-2007                                                          
///                                                                             
/// \todo    add doxygen documentation                                          
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
#include <numerics.h>
#include <WDMath.h>
#ifndef WDutils_included_cstring
#  include <cstring>
#  define WDutils_included_cstring
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
  static double eta=0.;
  if(eta == 0.) {
    eta = 1.e-7;
    while(eta+1. != 1.) eta *=0.5;
    eta  *=2.;                       // eta = actual computing accuracy
  }
  if(eps<eta) eps=eta;
  // get integration interval and return zero if it vanishes
  double ba=b-a;
  if(ba==0.) return 0.;
  // initialise some variables
  bool   bo,bu=0,odd=1;
  int    i,m,n=2,nn=3,mr;
  double d[7],dt[7];
  double ddt,hm,nt,err,t,ta,tab=(0),v=(0),w,sm(0),gr(0),t1(0),
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
      if(den) {
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
	if(sc == zero)
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
	if(d[i]) {
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
	  if(abs(e[m])+dd == dd) break;
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
	    if(r==zero) {
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
	  if(r==zero && i>=0) continue;
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
	  if(abs(e[m])+dd == dd) break;
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
	    if(r==zero) {
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
	  if(r==zero && i>=0) continue;
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

  template<typename scalar>
  class PercentileFinder {

  public:

    struct point {
      scalar q,w;
      point() {}
      point(point const&P) : q(P.q), w(P.w) {}
      point& operator=(point const&P) { q=P.q; w=P.w; return*this; }
    }; // struct point

    struct range {
      range   *r[2]; ///< lower and upper sub-range (if range is split)
      scalar   q[2]; ///< lower and upper value of q
      scalar   w[2]; ///< lower and upper cumulative weight
      int      n[2]; ///< lower and upper cumulative rank-1
      point   *P;    ///< table of points to be used

      /// default ctor: unsplit range
      range() { r[0]=0; r[1]=0; }

      point*begin() const { return P+n[0]; }
      point*end() const { return P+n[1]; }
      bool is_final() const { return r[0]==0; }
      bool contains(scalar W) const { return W>w[0] && W<=w[1]; }
      bool contains(int    R) const { return R>n[0] && R<=n[1]; }
      void dump() const {
	std::cerr
	  << "### WDutils debug info: "
	  << "PercentileFinder::range has:"
	  << " q="<<q[0]<<','<<q[1]
	  << " w="<<w[0]<<','<<w[1]
	  << " n="<<n[0]<<','<<n[1]<<'\n';
      }
    }; // struct range

    friend class WDutils::traits<point>;
    friend class WDutils::traits<range>;

  private:
    unsigned           Ntot;
    point             *PointsA, *PointsB;
    range             *Root;
    block_alloc<range> RangeAlloc;
    scalar             Wtot;

    void split(range*);

    template<typename T> range* findrange(T X) {
      register range* R = Root;
      while(R->q[1] > R->q[0]) {
	if(R->is_final()) split(R);
	R = R->r[0]->contains(X) ? R->r[0] : R->r[1];
      }
      return R;
    }

    template<typename T> double findQ(T const&X) {
      return findrange(X)->begin()->q;
    }

  public:
    /// ctor
    /// \param F array with values to be sorted
    /// \param N number of values
    /// \param W array with weights
    PercentileFinder(const scalar*F, int N, const scalar*W=0) WDutils_THROWING;

    /// dtor: de-allocate points and destruct RangeAlloc
    ~PercentileFinder()
    {
      WDutils_DEL_A(PointsA);
      WDutils_DEL_A(PointsB);
    }

    /// find value containing a given fraction of the total weight
    /// \param f fraction of the total weight
    /// \return value so that the cumulative weight equals f*Wtot
    scalar FindValue(scalar f) WDutils_THROWING {
      if(f<0. || f>1.) WDutils_THROW("error in percentile finding\n");
      return
	f==0 ? Root->q[0] :
	f==1 ? Root->q[1] : findQ(f*Wtot);
    }

    /// value corresponding to given rank
    /// \param r rank
    /// \return value at rank r
    scalar FindValue(int r) WDutils_THROWING {
      if(r<0. || r>Ntot) WDutils_THROW("error in percentile finding\n");
      return
	r==0    ? Root->q[0] :
	r==Ntot ? Root->q[1] : findQ(r);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  // does not compile with gcc 4.3.1, which seems a compiler bug
#if defined(__GNUC__) && ( __GNUC__ < 4 || __GNUC_MINOR__ < 3)
}
namespace WDutils {
  template<typename scalar>
  struct traits<typename PercentileFinder<scalar>::point> {
    static const char* name() { return "PercentileFinder::point"; }
  };
  template<typename scalar>
  struct traits<typename PercentileFinder<scalar>::range> {
    static const char* name() { return "PercentileFinder::range"; }
  };
}
namespace {
#endif
  //////////////////////////////////////////////////////////////////////////////

  template<typename scalar>
  PercentileFinder<scalar>::PercentileFinder(const scalar*F, 
					     int          N,
					     const scalar*W) WDutils_THROWING :
    Ntot       ( N ),
    PointsA    ( WDutils_NEW(point,Ntot) ),
    PointsB    ( WDutils_NEW(point,Ntot) ),
    RangeAlloc ( int(10+log(double(Ntot))) ),
    Wtot       ( 0. )
  {
    Root = RangeAlloc.new_element();
    scalar t, q[2] = { F[0], F[0] };
    for(int n=0; n!=N; ++n) {
      scalar f = F[n];
      update_min(q[0],f);
      update_max(q[1],f);
      PointsA[n].q = f;
      scalar w = W? W[n] : 1;
      if(w <= 0) WDutils_THROW("PercentileFinder: weight=%f <= 0\n",w);
      PointsA[n].w = w;
      Wtot += w;
    }
    Root->q[0] = q[0];
    Root->q[1] = q[1];
    Root->w[0] = 0.;
    Root->w[1] = Wtot;
    Root->n[0] = 0;
    Root->n[1] = Ntot;
    Root->P    = PointsA;
    if(debug(8)) Root->dump();
  }

  template<typename scalar>
  void PercentileFinder<scalar>::split(range*R) {
    point *P = R->P == PointsA ? PointsB : PointsA;
    scalar s = 0.5 * (R->q[0] + R->q[1]);
    scalar w = 0, q[2];
    int    l = R->n[0], h = R->n[1];
    q[0] = R->q[0];
    q[1] = R->q[1];
    for(point*p=R->begin(); p!=R->end(); ++p)
      if(p->q < s) {
	w     += p->w;
	P[l++] = *p;
	update_max(q[0],p->q);
      } else {
	P[--h] = *p;
	update_min(q[1],p->q);
      }
    R->r[0]       = RangeAlloc.new_element();
    R->r[0]->q[0] = R->q[0];
    R->r[0]->q[1] = q[0];
    R->r[0]->w[0] = R->w[0];
    R->r[0]->w[1] = R->w[0] + w;
    R->r[0]->n[0] = R->n[0];
    R->r[0]->n[1] = l;
    R->r[0]->P    = P;
    if(debug(8)) R->r[0]->dump();
      
    R->r[1]       = RangeAlloc.new_element();
    R->r[1]->q[0] = q[1];
    R->r[1]->q[1] = R->q[1];
    R->r[1]->w[0] = R->w[0] + w;
    R->r[1]->w[1] = R->w[1];
    R->r[1]->n[0] = h;
    R->r[1]->n[1] = R->n[1];
    R->r[1]->P    = P;
    if(debug(8)) R->r[1]->dump();
  }
} // namespace {
///////////////////////////////////////////////////////////////////////////////
namespace WDutils {
  template<typename scalar> void
  FindPercentile<scalar>::setup(const scalar*F, int N, const scalar*W)
  {
    if(DATA) WDutils_DEL_O(static_cast<PercentileFinder<scalar>*>(DATA));
    DATA = new PercentileFinder<scalar>(F,N,W);
  }
  template<typename scalar>
  FindPercentile<scalar>::~FindPercentile() {
    WDutils_DEL_O(static_cast<PercentileFinder<scalar>*>(DATA));
  }
  template<typename scalar>
  scalar FindPercentile<scalar>::Percentile(scalar f) const {
    return static_cast<PercentileFinder<scalar>*>(DATA)->FindValue(f);
  }
  template<typename scalar>
  scalar FindPercentile<scalar>::Percentile(int r) const {
    return static_cast<PercentileFinder<scalar>*>(DATA)->FindValue(r);
  }
  template class FindPercentile<double>;
  template class FindPercentile<float>;
} // namespace WDutils
////////////////////////////////////////////////////////////////////////////////
//
// finding index given rank
//
////////////////////////////////////////////////////////////////////////////////
namespace {

  /// auxiliary for FindIndex
  template<typename scalar>
  class Ranker {
    /// represents a single value
    struct point {
      scalar   X;       ///< value
      unsigned I;       ///< index
    };
    /// represents a range of points
    struct range {
      unsigned R;       ///< rank of first point
      unsigned N;       ///< number of points in range
      range   *S;       ///< pter to left sub-range
      range() : S(0) {}
    };
    unsigned           N;           ///< # points
    point             *P;           ///< table of points
    block_alloc<range> RangeAlloc;  ///< allocator for ranges
    range              Root;        ///< root range
    /// swap two points
    static void swap(point*a, point*b)
    {
      point p;
      memcpy(&p, a, sizeof(point));
      memcpy( a, b, sizeof(point));
      memcpy( b,&p, sizeof(point));
    }
    /// split a list of points
    /// \param[in] P list of points to split
    /// \param[in] N size of list
    /// \param[in] S split point
    /// \return      number of points with value < S (can be null)
    static unsigned splitlist(point*p, unsigned np, scalar S)
    {
      point*l,*u,*n=p+np;
      // find first element not smaller than S
      for(l=p;   l!=n && l->X<S; ++l);
      if(l==n) return np;
      // find first element not greater than S
      for(u=l+1; u!=n && u->X>S; ++u);
      while(u!=n) {
	swap(l,u);
	for(++l; l!=n && l->X<S; ++l);
	for(++u; u!=n && u->X>S; ++u);
      }
      return l - p;
    }
    /// split a range.
    /// we split at the element in mid-index.
    void split(range*A);
  public:
    /// ctor: create points and set root
    Ranker(const scalar*A, unsigned n, unsigned k=0) :
      N          ( n ),
      P          ( WDutils_NEW(point,N) ),
      RangeAlloc ( k? 4*k*int(1+log(double(N))) : 10*int(1+log(double(N))) )
    {
      Root.R   = 0;
      Root.N   = N;
      for(unsigned i=0; i!=N; ++i) {
	P[i].X = A[i];
	P[i].I = i;
      }
    }
    /// dtor: delete points
    ~Ranker() { WDutils_DEL_A(P); }
    /// find index for given rank
    unsigned Index(unsigned rank)
    {
      if(rank>=N) WDutils_THROW("FindRank::Index(): rank=%d >= N=$%d\n",rank,N);
      range*A=&Root;
      unsigned original_range_count = RangeAlloc.N_used();
      while(A->N>1) {
	if(A->S==0) split(A);
	A = rank<A->S[1].R? A->S : A->S+1;
      }
      DebugInfo(6,"Ranker::Index(): required %d new ranges\n",
		RangeAlloc.N_used() - original_range_count);
      return P[A->R].I;
    }
    /// find value for given rank
    scalar Value(unsigned rank)
    {
      if(rank>=N) WDutils_THROW("FindRank::Value(): rank=%d >= N=$%d\n",rank,N);
      range*A=&Root;
      unsigned original_range_count = RangeAlloc.N_used();
      while(A->N>1) {
	if(A->S==0) split(A);
	A = rank<A->S[1].R? A->S : A->S+1;
      }
      DebugInfo(6,"Ranker::Value(%u) = %f: required %d new ranges\n",
		rank, P[A->R].X, RangeAlloc.N_used() - original_range_count);
      return P[A->R].X;
    }
  };
}
namespace WDutils {
WDutils_TRAITS(::Ranker<float>::range,"<anonymous>::Ranker<float>::range");
WDutils_TRAITS(::Ranker<float>::point,"<anonymous>::Ranker<float>::point");
WDutils_TRAITS(::Ranker<double>::range,"<anonymous>::Ranker<double>::range");
WDutils_TRAITS(::Ranker<double>::point,"<anonymous>::Ranker<double>::point");
}
namespace {
  template<typename T>
  void Ranker<T>::split(range*A)
  {
    if(A->N >= 2) {
      unsigned L=1;
      if(A->N > 2)
	L = splitlist(P+A->R, A->N, P[A->R+A->N/2].X);
      else if(P[A->R].X > P[A->R+1].X)
	swap(P+A->R, P+A->R+1);
      A->S = RangeAlloc.new_elements(2);
      A->S[0].R   = A->R;    A->S[1].R   = A->R+L;
      A->S[0].N   = L;       A->S[1].N   = A->N-L;
    } else 
      WDutils_Error("FindIndex: trying to split range with %d points",A->N);
  }
} // namespace {
////////////////////////////////////////////////////////////////////////////////
namespace WDutils {
  template<typename scalar>
  void FindRank<scalar>::setup(const scalar*F, unsigned N,unsigned K) 
  {
    if(DATA) WDutils_DEL_O(static_cast<Ranker<scalar>*>(DATA));
    DATA = new Ranker<scalar>(F,N,K);
  }
  template<typename scalar>
  unsigned FindRank<scalar>::Index(unsigned rank) const
  {
    return static_cast<Ranker<scalar>*>(DATA)->Index(rank);
  }
  template<typename scalar>
  scalar FindRank<scalar>::Value(unsigned rank) const
  {
    return static_cast<Ranker<scalar>*>(DATA)->Value(rank);
  }
  template<typename scalar>
  FindRank<scalar>::~FindRank()
  {
    WDutils_DEL_O(static_cast<Ranker<scalar>*>(DATA));
    DATA = 0;
  }
  template class FindRank<double>;
  template class FindRank<float>;
} // namespace WDutils
////////////////////////////////////////////////////////////////////////////////
