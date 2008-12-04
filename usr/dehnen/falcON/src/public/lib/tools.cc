// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/public/tools.cc                                                 
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2002-2008                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2002-2008 Walter Dehnen                                        
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
  double   W(0.);                                  // Sum w_i                   
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
  unsigned          n,no(0);
  double            rh,r(hc),ro(0),dr,d1,d2;
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



