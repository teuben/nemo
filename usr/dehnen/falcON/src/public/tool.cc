// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// tool.cc                                                                     |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <public/tool.h>
#include <public/Pi.h>

using namespace nbdy;
////////////////////////////////////////////////////////////////////////////////
vect nbdy::centre_of_mass(const bodies*const&B) {
  register vect_d X=0.;
  register double W=0.;
  LoopBodies(bodies,B,Bi) {
    W += mass(Bi);
    X += mass(Bi) * pos(Bi);
  }
  return vect(X/W);
}
////////////////////////////////////////////////////////////////////////////////
namespace { using namespace nbdy; using nbdy::uint;
  inline double weight(real const&m, real const&p, uint const&alpha)
  {
    switch(alpha) {
    case 0u: return m;
    case 1u: return m * abs(p);
    case 2u: return m * power<2>(p);
    case 3u: return m * power<3>(abs(p));
    default: return m * pow(abs(p),alpha);
    }
  }
  inline double weight(bodies::iterator const&Bi, uint const&alpha)
  {
    return alpha? weight(mass(Bi), pot(Bi), alpha) : mass(Bi);
  }
}
////////////////////////////////////////////////////////////////////////////////
void nbdy::find_centre(const bodies*const&B,       // I  : bodies               
		       real         const&f,       // I  : reduction factor     
		       uint         const&alpha,   // I  : weigth=m*|Phi|^alpha 
		       uint         const&Nmin,    // I  : min #bodies in center
		       vect              &xc,      // I/O: center position      
		       real              &rc,      // I/O: centre radius        
		       bool         const&ini,     //[I  : use input as initial 
		       vect              *vc,      //[O  : center velocity]     
		       real              *rhc,     //[O  : estimate: rho(xc)]   
		       uint         const&b0,      //[I  : begin of bodies]     
		       uint         const&bn)      //[I  : end   of bodies]     
{
  const    real faq = square(f<one? f : 1./f);     // reduction factor squared  
  vect Xc;                                         // initial/current position  
  real Rq;                                         // initial/current radius^2  
  double W;                                        // Sum w_i                   
  vect_d X;                                        // Sum w_i x_i               
  uint   N;                                        // counter: # bodies         
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
    LoopBodiesRange(bodies,B,Bi,b0,bn) {           //   LOOP range of bodies    
      register double w = weight(Bi,alpha);        //     w_i = m_i*p_i^alpha   
      W+= w;                                       //     Sum w_i               
      X+= w * pos(Bi);                             //     Sum w_i x_i           
      Q+= w * norm(pos(Bi));                       //     Sum w_i x_i^2         
    }                                              //   END LOOP                
    Xc = X/W;                                      //   Xc  = <x>               
    Rq = Q/W - norm(X);                            //   R^2 = <(x-<x>)^2>       
//     // TEST
//     std::cerr<<" 0:"<<Xc<<"  "<<sqrt(Rq)<<'\n';
//     // TSET
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
    LoopBodiesRange(bodies,B,Bi,b0,bn)             //   LOOP range of bodies    
      if(dist_sq(pos(Bi),Xc) < Rq) {               //     IF |X-Xc| < R         
	register double w = weight(Bi,alpha);      //       w_i = m_i*p_i^alpha 
	W+= w;                                     //       Sum w_i             
	X+= w * pos(Bi);                           //       Sum w_i x_i         
	N++;                                       //       count               
      }                                            //   END LOOP                
    if(I && SHRINK && N < Nmin) SHRINK = 0;        //   IF (N < Nmin): adjust   
//     // TEST
//     std::cerr<<std::setw(2)<<I<<": X="<<Xc<<" R="<<sqrt(Rq)<<" N="<<N<<'\n';
//     // TSET
  }                                                // END DO                    
  xc = Xc;                                         // set output value: position
  rc = sqrt(Rq);                                   // set output value: radius  
  // 4. estimate velocity and density of center                                 
  if(vc || rhc) {                                  // IF V OR rho wanted        
    register double M(0.), W(0.);                  //   Sum m_i, Sum w_i        
    register vect_d V(0.);                         //   Sum w_i v_i             
    LoopBodiesRange(bodies,B,Bi,b0,bn)             //   LOOP range of bodies    
      if(dist_sq(pos(Bi),Xc) < Rq) {               //     IF |X-Xc| < R         
	M+= mass(Bi);                              //       Sum m_i             
	register double w = weight(Bi,alpha);      //       w_i = m_i*p_i^alpha 
	W+= w;                                     //       Sum w_i             
	V+= w * vel(Bi);                           //       Sum w_i v_i         
      }                                            //   END LOOP                
    if(vc)  *vc  = V/W;                            //   set velocity            
    if(rhc) *rhc = 3.*M/(FPi*Rq*rc);               //   set density estimate    
  }
}
////////////////////////////////////////////////////////////////////////////////
namespace { using namespace nbdy;
  template<int N=2> struct ferrers;
  template<> struct ferrers<2> {
    static double norm() { return FPi / 13.125; }
    static void   diff(double const&m, double const&xq, double D[3]) {
      D[2] = 8. * m;
      register double tmp = 1 - xq;
      D[1] =-times<4>(tmp) * m;
      D[0] = power<2>(tmp) * m;
    }
  };
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
    static double d1(double const&m, double const&xq) {
      return -6*m*square(1.-xq);
    }
  };
}
////////////////////////////////////////////////////////////////////////////////
#if (0)
void nbdy::find_centre(const bodies*const&B,       // I  : bodies               
		       uint         const&Nmin,    // I  : min #bodies in center
		       vect              &xc,      // I/O: center position      
		       real              &hc,      // I/O: centre radius        
		       vect              *vc,      //[O  : center velocity]     
		       real              *rhc,     //[O  : estimate: rho(xc)]   
		       uint         const&b0,      //[I  : begin of bodies]     
		       uint         const&bn)      //[I  : end   of bodies]     
{
  const int    miter  = 100;
  const double FACMAX = 1.;
  int    N;
  double RHO = 0., RH  = 0., d2RH, FAC=0.5;
  vect_d dRH = 1., Xc=xc, Vc, dXO, dX;
  for(int iter=0; iter<miter; ++iter) {
    N   = 0;
    RH  = 0.;
    dRH = 0.;
    d2RH= 0.;
    Vc  = 0.;
    const double hq=hc*hc, ihq=1./hq;
    LoopBodiesRange(bodies,B,Bi,b0,bn) {
      register vect_d R(Xc); R -= pos(Bi);
      register double Rq = norm(R);
      if(Rq < hq) {
	register double D[3];
	ferrers<2>::diff(mass(Bi),Rq*ihq,D);
	RH   += D[0];
	dRH  += R *= D[1];
	d2RH += times<3>(hq*D[1]) + Rq*D[2];
	Vc.add_times(vel(Bi),D[0]);
	++N;
      }
    }
    Vc  /= RH;
    RH  /= ferrers<2>::norm() * power<3>(hc);
    dRH /= ferrers<2>::norm() * power<5>(hc);
    d2RH/= ferrers<2>::norm() * power<7>(hc);
    dX   =-dRH / d2RH;
    // TEST
    std::cerr<<std::setw(2)<<iter<<": X="<<xc<<" H="<<hc
	     <<" N="<<N<<" drho="<<dRH<<" d2rho="<<d2RH
	     <<" rho="<<RH<<" dX="<<dX<<'\n';
    // TSET
    if(N == 0 || RH < RHO) {
      // new position wasn't better
      Xc  = xc;                         // goto old position
      RH  = RHO;                        // and old density 
      dX  = dXO;                        // and old dX
      FAC*= 0.5;                        // shrink factor
      if(N==0) hc += hc;                // IF N=0 increase hc
    } else if(RH > RHO) { 
      // found a position with higher density
      xc  = Xc;                         // remember it     
      dXO = dX;                         // remember dX
      RHO = RH;                         // remember density 
      if(rhc) *rhc = RH;                // update output
      if(vc)  *vc  = Vc;                // update output
      if(abs(dRH)*hc < 1.e-3*RH  && N == Nmin) return;
                                        // converged? then return
      hc *= 0.7+0.3*pow(double(Nmin)/double(N),0.33333333333333333);
                                        // adapt radius
      FAC = min(FACMAX, 1.1*FAC);       // adjust factor
    }
    Xc += FAC * dX * hc/sqrt(hc*hc+norm(dX)); // new guess for X_c
  }
  return;
}
#endif
////////////////////////////////////////////////////////////////////////////////
namespace {
  using namespace nbdy;
  void gr(const bodies*const&B,                    // I: bodies                 
	  uint         const&b0,                   // I: begin of bodies        
	  uint         const&bn,                   // I: end   of bodies        
	  vect_d       const&x,                    // I: trial position         
	  double       const&r,                    // I: trial radius           
	  uint              &n,                    // O: N(|r-x|<h)             
	  double            &rho,                  // O: rho_h(r)               
	  vect_d            &g,                    // O: - drho/dr              
	  vect_d            &v)                    // O: velocity               
  {
    const double rq = r*r, irq=1./rq;
    n   = 0;
    rho = 0.;
    v   = 0.;
    g   = 0.;
    LoopBodiesRange(bodies,B,Bi,b0,bn) {
      register vect_d R(x); R -= pos(Bi);
      register double Rq = norm(R);
      if(Rq < rq) {
	register double D[2];
	ferrers<3>::diff1(mass(Bi),Rq*irq,D);
	rho  += D[0];
	v.add_times(vel(Bi),D[0]);
	g.add_times(R,D[1]);
	++n;
      }
    }
    register double tmp = 1. / ( ferrers<3>::norm()*power<3>(r) );
    v   /= rho;
    rho *= tmp;
    tmp *= irq;
    g   *= tmp;
  }

  void di(const bodies*const&B,                    // I: bodies                 
	  uint         const&b0,                   // I: begin of bodies        
	  uint         const&bn,                   // I: end   of bodies        
	  vect_d       const&x,                    // I: trial position         
	  double       const&r,                    // I: trial radius           
	  vect_d       const&h,                    // I: trial direction        
	  double            &d1,                   // O: directional deriv ->h  
	  double            &d2)                   // O: 2nd ---                
  {
    const double rq = r*r, irq=1./rq;
    const double hqirq = norm(h) * irq;
    register double _d1 = 0.;
    register double _d2 = 0.;
    LoopBodiesRange(bodies,B,Bi,b0,bn) {
      register vect_d R(x); R -= pos(Bi);
      register double tmp = norm(R);
      if(tmp < rq) {
	register double D[3];
	ferrers<3>::diff(mass(Bi),tmp*irq,D);
	tmp  = h * R * irq;
	_d1 += tmp * D[1];
	_d2 += hqirq  * D[1] +  tmp * tmp * D[2];
      }
    }
    register double tmp = 1. / ( ferrers<3>::norm()*power<3>(r) );
    d1 = _d1 * tmp;
    d2 = _d2 * tmp;
  }
}
void nbdy::find_centre(const bodies*const&B,       // I  : bodies               
		       uint         const&N,       // I  : #bodies in center    
		       vect              &xc,      // I/O: center position      
		       real              &hc,      // I/O: centre radius        
		       vect              *vc,      //[O  : center velocity]     
		       real              *rhc,     //[O  : estimate: rho(xc)]   
		       uint         const&b0,      //[I  : begin of bodies]     
		       uint         const&bn)      //[I  : end   of bodies]     
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
//     vect_d dX;
//     std::cerr<<" give offset "<<std::endl;
//     std::cin >> dX;
//     xc += dX;
//   }
//   // TSET
  const int max_i = 100;
  uint              n;
  double            rh,r(hc),d1,d2;
  vect_d            x(xc), g, go, h, v;
  // initialize
  gr(B,b0,bn,x,r,n,rh,g,v);
  while(n==0) {
    r += r;
    gr(B,b0,bn,x,r,n,rh,g,v);
  }
  h = g;
//   // TEST
//   std::cerr<<" i: x="<<x<<" r="<<r<<" rh="<<rh<<" g="<<g<<" n="<<n<<'\n';
//   // TSET
  // iterate using cg method
  for(int i=0; i < max_i; ++i) {
    go = g;
    di(B,b0,bn,x,r,h,d1,d2);
    register vect_d dx(h);
    if(d2*r >=-abs(d1))    // 2nd derivate not really negative -> cannot use it
      dx *= r*sign(d1) / sqrt(norm(h)+square(r*d1));
    else                   // 2nd derivative < 0 as it should be near maximum  
      dx *= -d1/d2;
    x += dx;
    r *= 0.7+0.3*cbrt(double(N)/double(n));
    gr(B,b0,bn,x,r,n,rh,g,v);
    while(n==0) {
      r += r;
      gr(B,b0,bn,x,r,n,rh,g,v);
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
}
////////////////////////////////////////////////////////////////////////////////
namespace { using namespace nbdy;
  //----------------------------------------------------------------------------
  class dens_leaf : public basic_leaf {
    real      &wght()       { return SCAL; }
    real const&wght() const { return SCAL; }
  public:
    friend real const&wght  (const dens_leaf*const&L) { return L->wght(); }
    template<typename bodies_type>
    void set_weight(const bodies_type*const&B, uint const&A) {
      wght() = A ?
	weight(B->mass(mybody()), B->pot(mybody()), A) :
	B->mass(mybody());
    }
  };
  //----------------------------------------------------------------------------
  struct UpdateWeights {
    const oct_tree*tree;
    const uint     alfa;
    UpdateWeights(const oct_tree*const&t,
		  uint           const&a) : tree(t), alfa(a) {}
    template<typename bodies_type>
    uint operator() (const bodies_type*const&B) const {
      LoopLeafs(dens_leaf,tree,Li)
	Li->set_weight(B,alfa);
      return 0u;
    }
  };    
  //----------------------------------------------------------------------------
  class dens_cell : public basic_cell {
    real const&wght() const { return AUX1.SCAL; }
  public:
    typedef dens_leaf leaf_type;
    real&wght() { return AUX1.SCAL; }
    friend real const&wght  (const dens_cell*const&C) {return C->wght(); }
  };
  //----------------------------------------------------------------------------
  typedef oct_tree::CellIter<dens_cell> dens_cell_iter;
}
////////////////////////////////////////////////////////////////////////////////
void nbdy::estimate_density_peak(oct_tree*const&TREE,
				 uint     const&alpha,
				 uint     const&Nmin,
				 vect          &X0,
				 real          &H)
{
  // 1. set weights in leafs
  TREE->UseBodies(UpdateWeights(TREE,alpha));
  // 2. pass weights up the tree & find cell with maximum weight density
  dens_cell_iter max_cell;
  real           rho_max = zero;
  LoopCellsUp(dens_cell_iter,TREE,Ci) {
    register real w=0.;
    LoopLeafKids(dens_cell_iter,Ci,l) w += wght(l);
    LoopCellKids(dens_cell_iter,Ci,c) w += wght(c);
    Ci->wght() = w;
    w /= power<3>(times<2>(radius(Ci)));
    if(w > rho_max) {
      rho_max  = w;
      max_cell = Ci;
    }
  }
  // 3. loop leafs in cell with maximum weight density to find center
  vect xw(zero);
  LoopAllLeafs(dens_cell_iter,max_cell,l) xw += wght(l) * pos(l);
  X0 = xw / wght(max_cell);
  // 4. finally use the number density to estimate H
  H  = times<2>(radius(max_cell)) *
       pow(double(Nmin)/double(number(max_cell))/FPit, 0.3333333333333333);
}
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// implementing nbdy::find_lagrange_rad()                                       
//                                                                              
////////////////////////////////////////////////////////////////////////////////
#include <public/memo.h>

namespace {

  class RadiusFinder {

    struct point { 
      double q, m;           // q=r^2, mass

      point() {}
      point          (point const&P) : q(P.q), m(P.m) {}
      point&operator=(point const&P) { q=P.q; m=P.m; return *this; }

    } *PointsA, *PointsB;

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

    } *Root;

    block_alloc<range> RangeAlloc;
    double             Mtot;
    int                Ntot;

    void split(range*);

    template<typename T> double findQ(T const&X) {
      register range* R = Root;
      while(R->q[1] > R->q[0]) {
	if(R->is_final()) split(R);
	R = R->r[0]->contains(X) ? R->r[0] : R->r[1];
      }
      return R->begin_points()->q;
    }

  public:

    RadiusFinder(const bodies*const&, int const& = 1);
    // constructor: allocate PointsA,B; inititialize PointsA & Root

    double FindLagrangeRadius(double const&M) {
      if(M > 1.) warning("M/Mtot > 1 -> Lagrange radius = oo");
      return
	M<=0. ? 0. :
	M>=1. ? sqrt(Root->q[1]) : sqrt(findQ(M*Mtot));
    }

    double FindRankRadius(int const&R) {
      if(R > Ntot) warning("rank > N -> rank radius = oo");
      return
	R<1     ? 0. :
	R==1    ? sqrt(Root->q[0]) :
	R>=Ntot ? sqrt(Root->q[1]) : sqrt(findQ(R-1));
    }

    ~RadiusFinder() {
      delete[] PointsA;
      delete[] PointsB;
    }

  };

  void RadiusFinder::split(range*R) {    // split range
    point *P = R->P == PointsA ? PointsB : PointsA;
    const    double s = R->q[0] == 0.? 0.1*R->q[1] : sqrt(R->q[0] * R->q[1]);
    register double m = 0., q[2];
    register int    l = R->n[0], h=R->n[1];
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

  RadiusFinder::RadiusFinder(const bodies*const&B,
			     int          const&N) :
    PointsA    ( new point[B->N_bodies()] ),
    PointsB    ( new point[B->N_bodies()] ),
    RangeAlloc ( int(10+N*log(double(B->N_bodies()))) ),
    Mtot       ( 0. ),
    Ntot       ( B->N_bodies() ) {
    Root = RangeAlloc.new_element();
    register double t, q[2];
    q[1] = q[0] = norm(B->pos(0));
    for(int b=0; b!=B->N_bodies(); ++b) {
      t = norm(B->pos(b));
      if(std::isnan(t)) error("body position contains nan\n");
      update_min(q[0],t);
      update_max(q[1],t);
      PointsA[b].q = t;
      t            = B->mass(b);
      if(std::isnan(t)) error("body mass is nan\n");
      PointsA[b].m = t;
      Mtot        += t;
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
void nbdy::find_lagrange_rad(const bodies*const&B,
		                   int    const&N,
			     const double*const&M,
		                   double*const&R)
{
  RadiusFinder RF(B,N);
  for(int n=0; n!=N; ++n)
    R[n] = RF.FindLagrangeRadius(M[n]);
}




