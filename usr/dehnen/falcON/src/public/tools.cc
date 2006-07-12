// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/public/tools.cc                                                 
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2002-2006                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2002-2006 Walter Dehnen                                        
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
#include <public/tools.h>

using namespace falcON;
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// implementing falcON::centre_of_mass()                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
vect falcON::centre_of_mass(const bodies*B) {
  register vect_d X(0.);
  register double W(0.);
  LoopSubsetBodies(B,b) {
    W += mass(b);
    X += mass(b) * pos(b);
  }
  return W == zero? vect(zero) : vect(X/W);
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// implementing falcON::find_centre()                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace { using namespace falcON;
  inline double weight(real const&m, real const&p, unsigned const&alpha)
  {
    switch(alpha) {
    case 0u: return m;
    case 1u: return m * abs(p);
    case 2u: return m * power<2>(p);
    case 3u: return m * power<3>(abs(p));
    default: return m * pow(abs(p),alpha);
    }
  }
  inline double weight(bodies::iterator const&Bi, unsigned const&alpha)
  {
    return alpha? weight(mass(Bi), pot(Bi), alpha) : mass(Bi);
  }
}
////////////////////////////////////////////////////////////////////////////////
void falcON::find_centre_alpha(const bodies*B,
			       real         f,
			       unsigned     alpha,
			       unsigned     Nmin,
			       vect        &xc,
			       real        &rc,
			       bool         ini,
			       vect        *vc,
			       real        *rhc) falcON_THROWING
{
  unsigned nb = B->N_subset();
  if(nb < Nmin)
    falcON_THROW("find_centre(): # bodies=%d < Nmin=%d",nb,Nmin);
  const real faq = square(f<one? f : 1./f);        // reduction factor squared  
  vect     Xc;                                     // initial/current position  
  real     Rq;                                     // initial/current radius^2  
  double   W;                                      // Sum w_i                   
  vect_d   X;                                      // Sum w_i x_i               
  unsigned N;                                      // counter: # bodies         
  // 1. find starting position & radius for iteration                           
  if(ini) {                                        // IF use input as initial   
    Xc = xc;                                       //   set Xc  to given value  
    Rq = square(rc);                               //   set R^2 to given value  
    N  = 0u;                                       //   reset initial N         
  } else {                                         // ELSE (compute initial)    
    W = 0.;                                        //   reset Sum w_i           
    X = 0.;                                        //   reset Sum w_i x_i       
    N = 0u;                                        //   reset # bodies          
    register double Q(0.);                         //   reset Sum w_i x^2_i     
    LoopSubsetBodies(B,b) {                        //   LOOP bodies             
      register double w = weight(b,alpha);         //     w_i = m_i*p_i^alpha   
      W+= w;                                       //     Sum w_i               
      X+= w * pos(b);                              //     Sum w_i x_i           
      Q+= w * norm(pos(b));                        //     Sum w_i x_i^2         
    }                                              //   END LOOP                
    Xc = X/W;                                      //   Xc  = <x>               
    Rq = Q/W - norm(X);                            //   R^2 = <(x-<x>)^2>       
//     // TEST
//     std::cerr<<" 0:"<<Xc<<"  "<<sqrt(Rq)<<'\n';
//     // tensor_set
  }                                                // ENDIF                     
  // 2. iterate: exclude bodies at |X-X0| >= R, reducing R by fac each iteration
  int  I(0);                                       // counter: iterations       
  bool SHRINK(1);                                  // adjust OR shrink radius   
  for(; N != Nmin; ++I) {                          // DO until BREAK            
    if(I) {
      Xc  = X/W;                                   //   set new center position 
      Rq *= SHRINK? faq : pow(Nmin/double(N),0.6666666666666667);
    }
    W = 0.;                                        //   reset Sum w_i           
    X = 0.;                                        //   reset Sum w_i x_i       
    N = 0u;                                        //   reset # bodies          
    LoopSubsetBodies(B,b)                          //   LOOP bodies             
      if(dist_sq(pos(b),Xc) < Rq) {                //     IF |X-Xc| < R         
	register double w = weight(b,alpha);       //       w_i = m_i*p_i^alpha 
	W+= w;                                     //       Sum w_i             
	X+= w * pos(b);                            //       Sum w_i x_i         
	N++;                                       //       count               
      }                                            //   END LOOP                
    if(I && SHRINK && N < Nmin) SHRINK = 0;        //   IF (N < Nmin): adjust   
//     // TEST
//     std::cerr<<std::setw(2)<<I<<": X="<<Xc<<" R="<<sqrt(Rq)<<" N="<<N<<'\n';
//     // tensor_set
  }                                                // END DO                    
  xc = Xc;                                         // set output value: position
  rc = sqrt(Rq);                                   // set output value: radius  
  // 4. estimate velocity and density of center                                 
  if(vc || rhc) {                                  // IF V OR rho wanted        
    register double M(0.), W(0.);                  //   Sum m_i, Sum w_i        
    register vect_d V(0.);                         //   Sum w_i v_i             
    LoopSubsetBodies(B,b)                          //   LOOP bodies             
      if(dist_sq(pos(b),Xc) < Rq) {                //     IF |X-Xc| < R         
	M+= mass(b);                               //       Sum m_i             
	register double w = weight(b,alpha);       //       w_i = m_i*p_i^alpha 
	W+= w;                                     //       Sum w_i             
	V+= w * vel(b);                            //       Sum w_i v_i         
      }                                            //   END LOOP                
    if(vc)  *vc  = V/W;                            //   set velocity            
    if(rhc) *rhc = 3.*M/(FPi*Rq*rc);               //   set density estimate    
  }
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// implementing falcON::estimate_density_peak                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace { using namespace falcON;
  //----------------------------------------------------------------------------
  class DensLeaf : public OctTree::Leaf {
    real      &wght()       { return SCAL; }
    real const&wght() const { return SCAL; }
  public:
    friend real const&wght  (const DensLeaf*const&L) { return L->wght(); }
    template<typename bodies_type>
    void set_weight(const bodies_type*const&B, unsigned const&A) {
      wght() = A ?
	weight(B->mass(mybody()), B->pot(mybody()), A) :
	B->mass(mybody());
    }
  };
  //----------------------------------------------------------------------------
  class DensCell : public OctTree::Cell {
    real const&wght() const { return AUX1.SCAL; }
  public:
    typedef DensLeaf leaf_type;
    real&wght() { return AUX1.SCAL; }
    friend real const&wght  (const DensCell*const&C) {return C->wght(); }
  };
  //----------------------------------------------------------------------------
  typedef OctTree::CellIter<DensCell> DensCell_iter;
} // namespace {
////////////////////////////////////////////////////////////////////////////////
falcON_TRAITS(DensLeaf,"DensLeaf");
falcON_TRAITS(DensCell,"DensCell");
////////////////////////////////////////////////////////////////////////////////
void falcON::estimate_density_peak(OctTree *TREE,
				   unsigned alpha,
				   unsigned Nmin,
				   vect    &X0,
				   real    &H)
{
  // 1. set weights in leafs
  LoopLeafs(DensLeaf,TREE,Li)
    Li->set_weight(TREE->my_bodies(),alpha);
  // 2. pass weights up the tree & find cell with maximum weight density
  DensCell_iter max_cell;
  real          rho_max = zero;
  LoopCellsUp(DensCell_iter,TREE,Ci) {
    register real w=0.;
    LoopLeafKids(DensCell_iter,Ci,l) w += wght(l);
    LoopCellKids(DensCell_iter,Ci,c) w += wght(c);
    Ci->wght() = w;
    w /= power<3>(times<2>(radius(Ci)));
    if(w > rho_max) {
      rho_max  = w;
      max_cell = Ci;
    }
  }
  // 3. loop leafs in cell with maximum weight density to find center
  vect xw(zero);
  LoopAllLeafs(DensCell_iter,max_cell,l) xw += wght(l) * pos(l);
  X0 = xw / wght(max_cell);
  // 4. finally use the number density to estimate H
  H  = times<2>(radius(max_cell)) *
       pow(double(Nmin)/double(number(max_cell))/FPit, 0.3333333333333333);
}
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// implementing falcON::find_lagrange_rad()                                     
//                                                                              
////////////////////////////////////////////////////////////////////////////////
namespace {

  class RadiusFinder {

#if defined(__GNUC__) && (__GNUC__ < 3 || __GNUC__ == 3 && __GNUC_MINOR__ < 4)
    // patch for fix a bug with gcc version < 3.4
  public:
#endif

    struct point { 
      double q, m;           // q=r^2, mass

      point() {}
      point          (point const&P) : q(P.q), m(P.m) {}
      point&operator=(point const&P) { q=P.q; m=P.m; return *this; }

    };

    struct range {
      range   *r[2];         // lower & upper sub-range (if range is split)
      double   q[2];         // lower & upper q=r^2
      double   m[2];         // lower & upper cumulative mass
      int      n[2];         //(lower & upper cumulative rank) - 1
      point   *P;            // use this table of points.

      range() {              // default constructor: unsplit range
	r[0]=0;
	r[1]=0;
      }

      point*begin_points()               const { return P + n[0]; }
      point*end_points  ()               const { return P + n[1]; }
      bool  is_final    ()               const { return r[0] == 0; }
      bool  contains    (double const&M) const { return M>m[0] && M<=m[1]; }
      bool  contains    (int    const&R) const { return R>n[1] && R<=n[1]; }

    };

    friend class falcON::traits<point>;
    friend class falcON::traits<range>;

  private:
    unsigned           Ntot;
    point             *PointsA, *PointsB;
    range             *Root;
    block_alloc<range> RangeAlloc;
    double             Mtot;

    void split(range*);

    template<typename T> range* findrange(T const&X) {
      register range* R = Root;
      while(R->q[1] > R->q[0]) {
	if(R->is_final()) split(R);
	R = R->r[0]->contains(X) ? R->r[0] : R->r[1];
      }
      return R;
    }

    template<typename T> double findQ(T const&X) {
      return findrange(X)->begin_points()->q;
    }

  public:

    // constructor: allocate PointsA,B; inititialize PointsA & Root
    RadiusFinder(const bodies *    ,        // set of bodies
		 int            = 1,        // # M[r] anticipated
		 const vect   * = 0);       //[centre offset]

    double FindLagrangeRadius(double const&M) {
      // find radius containing the fraction M of the total mass
      if(M > 1.) warning("M/Mtot > 1 -> Lagrange radius = oo");
      return
	M<=0. ? 0. :
	M>=1. ? sqrt(Root->q[1]) : sqrt(findQ(M*Mtot));
    }
#if(0)
    double FindRankRadius(int const&R) {
      // find radius of the particle with radial rank R
      if(R >=Ntot) warning("rank >= N -> rank radius = oo");
      return
	R<0       ? 0. :
	R==0      ? sqrt(Root->q[0]) :
	R+1>=Ntot ? sqrt(Root->q[1]) : sqrt(findQ(R-1));
    }

    double FindRankRadiusAndMass(int const&Rank,
				 double   &Mcum) {
      // find radius of the particle X with radial rank R
      // also return mass within that radius, counting particle X half
      if(Rank > Ntot) warning("rank > N -> rank radius = oo");
      range *R = findrange(Rank);
      Mcum = 0.5*(R->m[0]+R->m[1]);
      return sqrt(R->begin_points()->q);
    }
#endif

    ~RadiusFinder() {
      falcON_DEL_A(PointsA);
      falcON_DEL_A(PointsB);
    }

  };
} // namespace {
////////////////////////////////////////////////////////////////////////////////
falcON_TRAITS(RadiusFinder::point,"RadiusFinder::point");
falcON_TRAITS(RadiusFinder::range,"RadiusFinder::range");
////////////////////////////////////////////////////////////////////////////////
namespace {
  void RadiusFinder::split(range*R) {    // split range
    point *P = R->P == PointsA ? PointsB : PointsA;
    const    double s = R->q[0] == 0.? 0.1*R->q[1] : sqrt(R->q[0] * R->q[1]);
    register double m = 0., q[2];
    register int    l = R->n[0], h = R->n[1];
    q[0] = R->q[0];
    q[1] = R->q[1];
    for(register point* p=R->begin_points(); p!=R->end_points(); ++p)
      if(p->q < s) {
	m     += p->m;
	P[l++] = *p;
	update_max(q[0],p->q);
      } else {
	P[--h] = *p;
	update_min(q[1],p->q);
      }
    
    R->r[0]       = RangeAlloc.new_element();
    R->r[0]->q[0] = R->q[0];
    R->r[0]->q[1] = q[0];
    R->r[0]->m[0] = R->m[0];
    R->r[0]->m[1] = R->m[0] + m;
    R->r[0]->n[0] = R->n[0];
    R->r[0]->n[1] = l;
    R->r[0]->P    = P;
      
    R->r[1]       = RangeAlloc.new_element();
    R->r[1]->q[0] = q[1];
    R->r[1]->q[1] = R->q[1];
    R->r[1]->m[0] = R->m[0] + m;
    R->r[1]->m[1] = R->m[1];
    R->r[1]->n[0] = h;
    R->r[1]->n[1] = R->n[1];
    R->r[1]->P    = P;

  }

  RadiusFinder::RadiusFinder(const bodies*B,
			     int          N,
			     const vect  *c) :
    Ntot       ( B->N_subset() ),
    PointsA    ( falcON_NEW(point,Ntot) ),
    PointsB    ( falcON_NEW(point,Ntot) ),
    RangeAlloc ( int(10+N*log(double(Ntot))) ),
    Mtot       ( 0. )
  {
    Root = RangeAlloc.new_element();
    register double t, q[2];
    q[1] = q[0] = c?
      dist_sq(pos(B->begin_all_bodies()),*c) :
      norm   (pos(B->begin_all_bodies()));
    LoopSubsetBodies(B,b) {
      t = c? dist_sq(pos(b),*c) : norm(pos(b));
#ifdef DEBUG
      if(std::isnan(t)) error("body position contains nan\n");
#endif
      update_min(q[0],t);
      update_max(q[1],t);
      PointsA[bodyindex(b)].q = t;
      t    = mass(b);
#ifdef DEBUG
      if(std::isnan(t)) error("body mass is nan\n");
#endif
      PointsA[bodyindex(b)].m = t;
      Mtot+= t;
    }
    Root->q[0] = q[0];
    Root->q[1] = q[1];
    Root->m[0] = 0.;
    Root->m[1] = Mtot;
    Root->n[0] = 0;
    Root->n[1] = Ntot;
    Root->P    = PointsA;
  }
}
////////////////////////////////////////////////////////////////////////////////
void falcON::find_lagrange_rad(const bodies*B,
			       unsigned     K,
			       const double*M,
			       double      *R,
			       const vect  *C)
{
  RadiusFinder RF(B,K,C);
  for(int k=0; k!=K; ++k)
    R[k] = RF.FindLagrangeRadius(M[k]);
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// implementing falcON::find_density_centre()                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace {
  using namespace falcON;
  //////////////////////////////////////////////////////////////////////////////
  template<int N=2> struct ferrers;
  template<> struct ferrers<3> {
    static double norm() { return FPi / 19.6875; }
    static void   diff(double const&m, double const&xq, double D[3]) {
      register double t = 1.-xq, d=m*t, d2=d*t, d3=d2*t;
      D[2] = 24*d;
      D[1] =-6*d2;
      D[0] = d3;
    }
    static void   diff1(double const&m, double const&xq, double D[2]) {
      register double t = 1.-xq, d2=m*t*t, d3=d2*t;
      D[1] =-6*d2;
      D[0] = d3;
    }
#if(0)
    static double d1(double const&m, double const&xq) {
      return -6*m*square(1.-xq);
    }
#endif
  };
}
////////////////////////////////////////////////////////////////////////////////
namespace {
  using namespace falcON;

  void gr(const bodies*B,                          // I: set of bodies          
	  vect_d const&x,                          // I: trial position         
	  double const&r,                          // I: trial radius           
	  unsigned    &n,                          // O: N(|r-x|<h)             
	  double      &rho,                        // O: rho_h(r)               
	  vect_d      &g,                          // O: - drho/dr              
	  vect_d      &v)                          // O: velocity (if v known)  
  {
    const double rq = r*r, irq=1./rq;
    n   = 0;
    rho = 0.;
    v   = 0.;
    g   = 0.;
    if(B->have_vel()) {
      LoopSubsetBodies(B,b) {
	register vect_d R(x); R -= pos(b);
	register double Rq = norm(R);
	if(Rq < rq) {
	  register double D[2];
	  ferrers<3>::diff1(mass(b),Rq*irq,D);
	  rho+= D[0];
	  v  += D[0] * vel(b);
	  g  += D[1] * R;
	  ++n;
	}
      }
      v /= rho;
    } else {
      LoopSubsetBodies(B,b) {
	register vect_d R(x); R -= pos(b);
	register double Rq = norm(R);
	if(Rq < rq) {
	  register double D[2];
	  ferrers<3>::diff1(mass(b),Rq*irq,D);
	  rho+= D[0];
	  g  += D[1] * R;
	  ++n;
	}
      }
    }
    register double tmp = 1. / ( ferrers<3>::norm()*power<3>(r) );
    rho *= tmp;
    tmp *= irq;
    g   *= tmp;
  }

  void di(const bodies*B,                          // I: set of bodies          
	  vect_d const&x,                          // I: trial position         
	  double const&r,                          // I: trial radius           
	  vect_d const&h,                          // I: trial direction        
	  double      &d1,                         // O: directional deriv ->h  
	  double      &d2)                         // O: 2nd ---                
  {
    const double rq = r*r, irq=1./rq;
    const double hqirq = norm(h) * irq;
    register double _d1 = 0.;
    register double _d2 = 0.;
    LoopSubsetBodies(B,b) {
      register vect_d R(x); R -= pos(b);
      register double tmp = norm(R);
      if(tmp < rq) {
	register double D[3];
	ferrers<3>::diff(mass(b),tmp*irq,D);
	tmp  = h * R * irq;
	_d1 += tmp * D[1];
	_d2 += hqirq  * D[1] +  tmp * tmp * D[2];
      }
    }
    register double tmp = 1. / ( ferrers<3>::norm()*power<3>(r) );
    d1 = _d1 * tmp;
    d2 = _d2 * tmp;
  }

} // namespace {

bool falcON::find_density_centre(const bodies*B,
				 unsigned     N,
				 vect        &xc,
				 real        &hc,
				 vect        *vc,
				 real        *rhc)
{
  // we use a conjugate gradient method, whereby approximating the line
  // maximisation by the 1st and 2nd directional derivatives.
  // Near the maximum, the 2nd derivative should be negative.
  // If it isn't, we cannot use it for line maximisation, but simply go
  // in the direction of h at most the size of the current radius far.
  // 
  // It seems that the algorithm converges, even if initially set off by more
  // than the initial radius (which already was 3 times the expectation).

//   // TEST
//   {
//     std::cerr<<" initial xc="<<xc<<" hc="<<hc<<'\n'
// 	     <<" do you want to change the initial xc,hc (1/0)? ";
//     bool want;
//     std::cin >> want;
//     if(want) {
//       std::cerr<<" give xc"; std::cin>>xc;
//       std::cerr<<" give hc"; std::cin>>hc;
//     }
//   }
//   // TSET

  const unsigned nb = B->N_subset();
  if(nb < N)
    falcON_THROW("find_density_centre(): N=% < Ncen = %\n",nb,N);

  const int max_i = 100;
  unsigned          n,no;
  double            rh,r(hc),ro,dr,d1,d2;
  vect_d            x(xc), g, go, h, v;
  if(r <= 0.) r=0.1;
  // initialize
  gr(B,x,r,n,rh,g,v);
  while(n==0) {
    r += r;
    gr(B,x,r,n,rh,g,v);
  }
  h = g;
//   // TEST
//   std::cerr<<" i: x="<<x<<" r="<<r<<" rh="<<rh<<" g="<<g<<" n="<<n<<'\n';
//   // TSET
  // iterate using cg method
  int i=0;
  for(; i < max_i; ++i) {
    di(B,x,r,h,d1,d2);
    dr = ro-r;
    ro = r;
    if(i && (no-N)*(n-N)<=0 && dr!=0)
      r += dr * double(N-n)/double(no-n);
    else if(n!=N)
      r *= 0.7+0.3*cbrt(double(N)/double(n));
    register vect_d dx(h);
    if(d2*r >=-abs(d1))    // 2nd derivate not really negative -> cannot use it
      dx /= d1;
    else                   // 2nd derivative < 0 as it should be near maximum  
      dx *= -d1/d2;
    dx*= r / sqrt(norm(dx)+r*r);   // restrict total size of step to r        
    x += dx;
    no = n;
    go = g;
    gr(B,x,r,n,rh,g,v);
    while(n==0) {
      r += r;
      gr(B,x,r,n,rh,g,v);
    }
//     // TEST
//     std::cerr<<std::setw(2)<<i
// 	     <<": x="<<x<<" r="<<r<<" rh="<<rh<<" g="<<g<<" h="<<h
// 	     <<" d1="<<d1
// 	     <<" d2="<<d2
// 	     <<" dx="<<dx
// 	     <<" n="<<n<<'\n';
//     // TSET
    if(abs(g)*r<1.e-8*rh && n==N) break;
    h  = g + h * (((g-go)*g)/(go*go));
  }
  xc = x;
  hc = r;
  if(rhc) *rhc = rh;
  if(vc)  *vc  = v;
  return i < max_i;
}
////////////////////////////////////////////////////////////////////////////////
#ifdef falcON_PROPER
#  include <proper/tools.cc>
#endif



