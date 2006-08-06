// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// center.h                                                                    |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2004                                               |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// defines class centering that may be used to implement centered external     |
// acceleration fields.                                                        |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_center_h
#define falcON_included_center_h

#ifndef falcON_included_algorithm
#  include <algorithm>
#  define falcON_included_algorithm
#endif
#ifndef falcON_included_tupel_h
#  include <utils/tupel.h>
#endif
#ifndef falcON_included_inline_h
#  include <utils/inline.h>
#endif
#ifndef falcON_included_iomanip
#  include <iomanip>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace {
  using WDutils::tupel;
  using WDutils::abs;
  using WDutils::norm;
  using std::abs;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class centering                                                          //
  //                                                                          //
  // finds the global density maximum of a given body distribution as the     //
  // position where                                                           //
  //                                                                          //
  //               1                                                          //
  //   rho(xc) := ---  Sum           m_i W(|x_i-x_c|)                         //
  //              h^3  |x_i-xc| < h                                           //
  //                                                                          //
  // becomes maximal with h set such that the number of bodies contributing   //
  // equals Nmin. The kernel W(r) is taken to be the Ferrers n=3 sphere.      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class centering {
  private:
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
    mutable double H;                              // size of center            
    mutable double C[3];                           // position of center        
    int            Nmin;                           // # bodies in center        
//     // TEST
//     mutable bool   TEST;
//     // TSET
    //--------------------------------------------------------------------------
    // kernel                                                                   
    //--------------------------------------------------------------------------
    struct kernel {
      template<int NDIM> static double norm() { 
	const double Pi = 3.14159265358979323846264338328;
	return NDIM == 2? 0.25*Pi : 4*Pi / 19.6875;
      }
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
	return -6*m*WDutils::square(1.-xq);
      }
    };
    //--------------------------------------------------------------------------
    // compute gradient of density w.r.t. center position                       
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar> static
    void gr(const scalar            *M,            // I: masses                 
	    const tupel<NDIM,scalar>*X,            // I: positions              
	    int                const&N,            // I: size of arrays         
	    tupel<NDIM,double> const&x,            // I: trial position         
	    double             const&r,            // I: trial radius           
	    int                     &n,            // O: N(|r-x|<h)             
	    double                  &rho,          // O: rho_h(r)               
	    tupel<NDIM,double>      &g)            // O: - drho/dr              
    {
      const double rq = r*r, irq=1./rq;
      n   = 0;
      rho = 0.;
      g   = 0.;
      for(int b=0; b!=N; ++b) {
	register tupel<NDIM,double> R(x); R -= X[b];
	register double Rq = norm(R);
	if(Rq < rq) {
	  register double D[2];
	  kernel::diff1(M[b],Rq*irq,D);
	  rho  += D[0];
	  g    += R * D[1];
	  ++n;
	}
      }
      register double tmp = 1. / ( kernel::norm<NDIM>() *r*r*r );
      rho *= tmp;
      tmp *= irq;
      g   *= tmp;
    }
    //--------------------------------------------------------------------------
    // compute 1st and 2nd derivative of rho in given direction                 
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar> static
    void di(const scalar            *M,            // I: masses                 
	    const tupel<NDIM,scalar>*X,            // I: positions              
	    int                const&N,            // I: size of arrays         
	    tupel<NDIM,double> const&x,            // I: trial position         
	    double             const&r,            // I: trial radius           
	    tupel<NDIM,double> const&h,            // I: trial direction        
	    double                  &d1,           // O: directional deriv ->h  
	    double                  &d2)           // O: 2nd ---                
    {
      const double rq = r*r, irq=1./rq, hqirq = norm(h) * irq;
      d1 = 0.;
      d2 = 0.;
      for(int b=0; b!=N; ++b) {
	register tupel<NDIM,double> R(x); R -= X[b];
	register double tmp = norm(R);
	if(tmp < rq) {
	  register double D[3];
	  kernel::diff(M[b],tmp*irq,D);
	  tmp = (h*R) * irq;
	  d1 += tmp * D[1];
	  d2 += hqirq  * D[1] +  tmp * tmp * D[2];
	}
      }
      register double tmp = 1. / ( kernel::norm<NDIM>() *r*r*r );
      d1 *= tmp;
      d2 *= tmp;
    }
    //--------------------------------------------------------------------------
    // compute the center of mass and the rms radius of all bodies              
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar> static
    void centre_of_mass(const scalar            *M,// I: masses                 
			const tupel<NDIM,scalar>*X,// I: positions              
			int                const&N,// I: size of arrays         
			tupel<NDIM,double>      &c,// O: centre of mass         
			double                  &r)// O: rms radius             
    {
      r        = 0.;
      c        = 0.;
      double m = 0.;
      for(int b=0; b!=N; ++b) {
	c += M[b] * X[b];
	r += M[b] * norm(X[b]);
	m += M[b];
      }
//       // TEST
//       std::cerr<<"\n centre_of_mass(): M="<<m<<'\n';
//       // TSET
      m  = 1./m;
      c *= m;
      r *= m;
      r -= norm(c);
      r  = std::sqrt(r);
    }
    //--------------------------------------------------------------------------
    // compute the center of mass of those bodies within r of position x        
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar> static
    int  centre_of_mass(                           // R: number satisfying crit 
			const scalar            *M,// I: masses                 
			const tupel<NDIM,scalar>*X,// I: positions              
			int                const&N,// I: size of arrays         
			tupel<NDIM,double> const&x,// I: consider only |X-x|    
			double             const&r,// I:                < r     
			tupel<NDIM,double>      &c)// O: centre of mass         
    {
      const double rq=r*r;
      c        = 0.;
      int    n = 0;
      double m = 0.;
      for(int b=0; b!=N; ++b)
	if(rq > x.dist_sq(X[b])) {
	  c += M[b] * X[b];
	  m += M[b];
	  ++n;
	}
      c /= m;
      return n;
    }
    //--------------------------------------------------------------------------
    // find the maximimum of rho(xc) given a trial position, using the gc method
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar>
    bool iterate(                                  // R: was okay?              
		 const scalar            *M,       // I: masses                 
		 const tupel<NDIM,scalar>*X,       // I: positions              
		 int                const&N) const // I: size of arrays         
    {
      // We use a conjugate gradient method, whereby approximating the line     
      // maximisation by the 1st and 2nd directional derivatives. Near the      
      // maximum, the 2nd derivative should be negative. If it isn't, we cannot 
      // use it for line maximisation, but simply go in the direction of h at   
      // most the size of the current radius far.                               
      //                                                                        
      // It seems that the algorithm converges, even if initially set off by    
      // more than the initial radius.                                          
      
      const int max_i = 200;
      int       n,nr,no;
      double    rh,r(H),rr,ro,dr,d1,d2;
      tupel<NDIM,double> x(C),g,go,h;
      // initialize
      gr<NDIM,scalar>(M,X,N,x,r,n,rh,g);
      while(n==0) {
	r += r;
	gr<NDIM,scalar>(M,X,N,x,r,n,rh,g);
      }
      h = g;
//       // TEST
//       std::ofstream TEST("GrowRidigBar.test");
//       TEST<<"\n centering::iterate():\n"
// 	  <<" i: x="<<x<<" r="<<r<<" rh="<<rh<<" g="<<g<<" n="<<n<<'\n';
//       // TSET
      // iterate using cg method
      int i=0, ir, io;
      for(; i < max_i; ++i) {
	di<NDIM,scalar>(M,X,N,x,r,h,d1,d2);
	dr = ro-r;
	rr = r;
// 	// TEST
// 	bool inter=i>1 && i-io<4 && (no-Nmin)*(n-Nmin)<0 && dr!=0.;
// 	TEST<<" ro="<<ro<<" no="<<no
// 	    <<" r="<<r<<" n="<<n<<" inter="<<inter<<'\n';
// 	// TSET
	if(i>1 && i-io<4 && (no-Nmin)*(n-Nmin)<0 && dr!=0.)
	  r += dr * double(Nmin-n)/double(no-n);
	else if(n!=Nmin)
	  r *= 0.6+0.4*cbrt(double(Nmin)/double(n));
	register tupel<NDIM,double> dx(h);
	if(d2*r >=-abs(d1)) dx /= d1;
	else                dx *=-d1/d2;
	dx*= r/sqrt(norm(dx)+r*r);
	x += dx;
	nr = n;
	ir = i;
	go = g;
	gr<NDIM,scalar>(M,X,N,x,r,n,rh,g);
	while(n==0) {
	  r += r;
	  gr<NDIM,scalar>(M,X,N,x,r,n,rh,g);
	}
	if(i==0 || (Nmin-n)*(Nmin-nr) < 0 ) {
	  io = ir;
	  ro = rr;
	  no = nr;
	}
// 	// TEST
// 	TEST<<std::setw(2)<<i
// 	    <<": x="<<x<<" r="<<r<<" rh="<<rh<<" g="<<g<<" h="<<h
// 	    <<" d1="<<d1
// 	    <<" d2="<<d2
// 	    <<" dx="<<dx
// 	    <<" n="<<n<<' '<<inter<<'\n';
// 	// TSET
	if(abs(dx) < 1.e-6*r && abs(Nmin-n)<2) break;
// 	if(abs(g)*r<1.e-8*rh && abs(Nmin-n)<2) break;
	if(i==0 || i % 50)
	  h = g + h * (((g-go)*g)/(go*go));
	else
	  h = g;
      }
      if(i < max_i) {
	H    = r;
	C[0] = x[0];
	C[1] = x[1];
	C[2] = NDIM > 2? x[2] : 0.;
	return true;
      }
//       // TEST
//       TEST.flush();
//       error(" centering did not converge");
//       // TSET
      return false;
    }
    //--------------------------------------------------------------------------
    // update: (re-) compute the center position                                
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar>
    bool _update(                                  // R: was okay?              
		 const scalar            *M,       // I: masses                 
		 const tupel<NDIM,scalar>*X,       // I: positions              
		 int                const&N) const // I: size of arrays         
    {
      if(H == 0.) {                                // IF first time centering   
	tupel<NDIM,double> x;                      //   we need a trial position
	centre_of_mass<NDIM,scalar>(M,X,N,x,H);    //   center of mass          
	tupel<NDIM,double> c=x;
	double h = H;
	while(1) {
	  int n = centre_of_mass<NDIM,scalar>(M,X,N,c,h,x);
	  if(4*Nmin > n) break;
	  H    = h;
	  C[0] = c[0];
	  C[1] = c[1];
	  C[2] = NDIM > 2? c[2] : 0.;
	  c    = x;
	  h   *= 0.95;
	}
      }                                            // ENDIF                     
      return iterate<NDIM,scalar>(M,X,N);          // iterate from trial        
    }
  protected:
    //--------------------------------------------------------------------------
    centering() : H(0.) 
//       // TEST
// 		  , TEST(0)
//       // TSET
    {}
    //--------------------------------------------------------------------------
//     // TEST
//     bool& test() const { return TEST; }
//     // TSET
    //--------------------------------------------------------------------------
    // instead of construction                                                  
    //--------------------------------------------------------------------------
    void init(int n) {
      Nmin = n;
    }
    //--------------------------------------------------------------------------
    // update: (re-) compute the center position                                
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar>
    bool update(                                   // converged?                
		const scalar*const&M,              // I: masses                 
		const scalar*const&X,              // I: positions              
		int          const&N) const        // I: size of arrays         
    {
      return _update(M,static_cast<const tupel<NDIM,scalar>*>
		     (static_cast<const void*>(X)), N);
    }
    //--------------------------------------------------------------------------
    // return number of bodies in center                                        
    //--------------------------------------------------------------------------
    int const&N_min () const { return Nmin; }
    //--------------------------------------------------------------------------
    // return center position found with last update.                           
    //--------------------------------------------------------------------------
    const double* centre() const { return C; }
  };
} // namespace {
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_center_h
