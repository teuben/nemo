// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// tools.cc                                                                    |
//                                                                             |
// Copyright (C) 2002-2005  Walter Dehnen                                      |
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
#include <public/tools.h>

using namespace falcON;
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// implementing falcON::centre_of_mass()                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
vect falcON::centre_of_mass(const bodies*const&B) {
  register vect_d X(0.);
  register double W(0.);
  LoopAllBodies(B,Bi) {
    W += mass(Bi);
    X += mass(Bi) * pos(Bi);
  }
  return vect(X/W);
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
void falcON::find_centre(const bodies*const&B,     // I  : bodies               
			 real         const&f,     // I  : reduction factor     
			 unsigned     const&alpha, // I  : weigth=m*|Phi|^alpha 
			 unsigned     const&Nmin,  // I  : min #bodies in center
			 vect              &xc,    // I/O: center position      
			 real              &rc,    // I/O: centre radius        
			 bool               ini,   //[I  : use input as initial 
			 vect              *vc,    //[O  : center velocity]     
			 real              *rhc,   //[O  : estimate: rho(xc)]   
			 unsigned           b0,    //[I  : first body]          
			 unsigned           nb)    //[I  : # bodies]            
  falcON_THROWING
{
  if(b0 >= B->N_bodies())
    falcON_THROW("find_centre(): first body >= N");
  else if(0 == nb)
    nb = B->N_bodies() - b0;
  else if(b0 + nb > B->N_bodies()) {
    warning("find_centre(): cannot use %d, but only %d bodies",
	    nb, B->N_bodies()-b0);
    nb = B->N_bodies() - b0;
  }
  if(nb < Nmin)
    falcON_THROW("find_centre(): # bodies=%d < Nmin=%d",nb,Nmin);
  const real faq = square(f<one? f : 1./f);        // reduction factor squared  
  const body first(B->bodyNo(b0));                 // first body considered     
  const body end(B->bodyNo(b0+nb));                // end of bodies             
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
    for(body Bi(first); Bi!=end; ++Bi) {           //   LOOP bodies             
      register double w = weight(Bi,alpha);        //     w_i = m_i*p_i^alpha   
      W+= w;                                       //     Sum w_i               
      X+= w * pos(Bi);                             //     Sum w_i x_i           
      Q+= w * norm(pos(Bi));                       //     Sum w_i x_i^2         
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
    for(body Bi(first); Bi!=end; ++Bi)             //   LOOP bodies             
      if(dist_sq(pos(Bi),Xc) < Rq) {               //     IF |X-Xc| < R         
	register double w = weight(Bi,alpha);      //       w_i = m_i*p_i^alpha 
	W+= w;                                     //       Sum w_i             
	X+= w * pos(Bi);                           //       Sum w_i x_i         
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
    for(body Bi(first); Bi!=end; ++Bi)             //   LOOP bodies             
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
falcON_TRAITS(DensLeaf,"DensLeaf","DensLeafs");
falcON_TRAITS(DensCell,"DensCell","DensCells");
////////////////////////////////////////////////////////////////////////////////
void falcON::estimate_density_peak(OctTree *const&TREE,
				   unsigned const&alpha,
				   unsigned const&Nmin,
				   vect         &X0,
				   real         &H)
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
    RadiusFinder(body  const &,
		 unsigned     ,
		 int           = 1,        // # M[r] anticipated
		 const vect  * = 0);       //[centre offset]

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
      delete[] PointsA;
      delete[] PointsB;
    }

  };
} // namespace {
////////////////////////////////////////////////////////////////////////////////
falcON_TRAITS(RadiusFinder::point,
	      "RadiusFinder::point","RadiusFinder::points");
falcON_TRAITS(RadiusFinder::range,
	      "RadiusFinder::range","RadiusFinder::range");
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

  RadiusFinder::RadiusFinder(body  const &b0,
			     unsigned     nb,
			     int          Nrad,
			     const vect  *c) :
    Ntot       ( nb ),
    PointsA    ( falcON_NEW(point,nb) ),
    PointsB    ( falcON_NEW(point,nb) ),
    RangeAlloc ( int(10+Nrad*log(double(nb))) ),
    Mtot       ( 0. )
  {
    Root = RangeAlloc.new_element();
    register double t, q[2];
    q[1] = q[0] = c? dist_sq(pos(b0),*c) : norm(pos(b0));
    const body bn(b0,nb);
    for(body b(b0); b!=bn; ++b) {
      t = c? dist_sq(pos(b),*c) : norm(pos(b));
      if(std::isnan(t)) error("body position contains nan\n");
      update_min(q[0],t);
      update_max(q[1],t);
      PointsA[bodyindex(b)].q = t;
      t    = mass(b);
      if(std::isnan(t)) error("body mass is nan\n");
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
void falcON::find_lagrange_rad(body  const &b,
			       int          K,
			       const double*M,
			       double      *R,
			       const vect  *C,
			       unsigned     N)
{
  const bodies*B = b.my_bodies();
  if(N==0 || bodyindex(b) + N > B->N_bodies())
    N = B->N_bodies() - bodyindex(b);
  RadiusFinder RF(b,N,K,C);
  for(int k=0; k!=K; ++k)
    R[k] = RF.FindLagrangeRadius(M[k]);
}
////////////////////////////////////////////////////////////////////////////////
#ifdef falcON_PROPER
#  include <proper/tools.cc>
#endif



