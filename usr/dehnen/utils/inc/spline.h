// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    inc/spline.h                                                       
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    1994-2005                                                          
///                                                                             
/// \todo    finish doxygen documentation                                       
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
#ifndef WDutils_included_spline_h
#define WDutils_included_spline_h

#ifndef WDutils_included_memory_h
#  include <memory.h>
#endif
#ifndef WDutils_included_numerics_h
#  include <numerics.h>
#endif

////////////////////////////////////////////////////////////////////////////////
namespace WDutils {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class WDutils::cubic_splines                                               
  //                                                                            
  /// this is a mere collection of (static) methods                             
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  template<class scalar_type, class table_type = scalar_type>
  class cubic_splines {
  public:
    //--------------------------------------------------------------------------
    /// construct cubic spline
    /// \param n  (input) size of tables
    /// \param x  (input) table of points
    /// \param y  (input) table of function values
    /// \param y2 (output) table of y''
    /// \param yp1 (optional input) pter to y' at x[0], use y''=0 if not given
    /// \param ypn (optional input) pter to y' at x[n-1], use y''=0 if not given
    static void construct(int n,
			  const scalar_type*x,
			  const table_type *y,
			  table_type *const y2,
			  const table_type *yp1=0,
			  const table_type *ypn=0)
    { 
      const scalar_type
	half =scalar_type(0.5),
	one  =scalar_type(1),
	two  =scalar_type(2),
	three=scalar_type(3),
	six  =scalar_type(6);
      register scalar_type  qn,p,sig,dx,dx1,dx2;
      register table_type   un,dy,dy1;
      table_type  *u = WDutils_NEW(table_type,n-1);
      scalar_type *v = WDutils_NEW(scalar_type,n-1);
      dx = x[1] - x[0];
      dy = y[1] - y[0];
      if(yp1) {
	v[0] =-half;
	u[0] = three/dx * (dy/dx- *yp1);
      } else
	u[0] = v[0] = table_type(0);
      for(register int i=1; i!=n-1; ++i) {
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
      if(ypn) {
	qn = half;
	un = three/dx * (*ypn - dy/dx);
      } else
	un = qn = table_type(0);
      y2[n-1]= (un-qn*u[n-2]) / (qn*v[n-2]+one);
      for(register int i=n-2; i>=0; --i)
        y2[i] = v[i]*y2[i+1] + u[i];
      WDutils_DEL_A(u);
      WDutils_DEL_A(v);
    }
    //--------------------------------------------------------------------------
    /// construct cubic spline                                                  
    /// \param n  (input) size of tables
    /// \param x  (input) table of points
    /// \param y  (input) table of function values
    /// \param y2 (output) table of y''
    /// \param yp1 (input) y' at x[0]
    /// \param ypn (input) y' at x[n-1]
    /// \param y2 (output) table of y''(x)
    /// \param n1 (optional input) use natural boundary at x[0]?
    /// \param nn (optional input) use natural boundary at x[n-1]?
    static void construct(int n,
			  const scalar_type*x,
			  const table_type *y,
			  table_type const&yp1,
			  table_type const&ypn,
			  table_type*const y2,
			  bool n1=false,
			  bool nn=false)
    {
      construct(n,x,y,y2, n1? &yp1 : 0, nn? &ypn : 0);
    }
    //--------------------------------------------------------------------------
    /// evaluate spline at x=xi with xl <= xi <= xh and yl=y(xl) etc...
    /// \param xi (input) x-value where spline is desired
    /// \param xl (input) x-point left to xi: xl <= xi
    /// \param xh (input) x-point right to xi: xi <= xh
    /// \param yl (input) y(xl)
    /// \param yh (input) y(xh)
    /// \param y2l (input) y''(xl)
    /// \param y2h (input) y''(xh)
    /// \param yi  (output) y(xi)
    /// \param dyi (optional output) y'(xi)
    /// \param d2y1 (optional output) y''(xi)
    static void evaluate(scalar_type const&xi,
			 scalar_type const&xl,
			 scalar_type const&xh,
			 table_type  const&yl,
			 table_type  const&yh,
			 table_type  const&y2l,
			 table_type  const&y2h,
			 table_type       *yi,
			 table_type       *dyi = 0,
			 table_type       *d2yi= 0)
    {
      const scalar_type
	zero =scalar_type(0),
	one  =scalar_type(1),
	three=scalar_type(3),
	six  =scalar_type(6);
      register scalar_type h=xh-xl,h6,A,B;
      if(h==zero) WDutils_Error("cubic spline: bad x input");
      h6 = h / six;
      A  = (xh-xi) / h;
      B  = one - A;
      if(dyi) {
        register scalar_type Aq=A*A, Bq=B*B;
        *dyi = (yh-yl)/h + h6*((three*Bq-one)*y2h-(three*Aq-one)*y2l);
        if(d2yi) *d2yi = A*y2l + B*y2h;
        *yi = A*yl +B*yh +((Aq-one)*A*y2l+(Bq-one)*B*y2h)*(h*h6);
      }
      *yi = A*yl +B*yh +(((A*A-one)*A)*y2l+((B*B-one)*B)*y2h)*(h*h6);
    }
    //--------------------------------------------------------------------------
    /// evaluate spline at x=xi                                                 
    /// \param xi (input) x-value where spline is desired
    /// \param x  (input) array with x[0] <= xi <= x[1]
    /// \param y  (input) array with y(x)
    /// \param y2 (input) array with y''(x)
    /// \param yi  (output) y(xi)
    /// \param dyi (optional output) y'(xi)
    /// \param d2y1 (optional output) y''(xi)
    static table_type evaluate(scalar_type const &xi,
			       const scalar_type *x,
			       const table_type  *y,
			       const table_type  *y2,
			       table_type  *dyi = 0,
			       table_type  *d2yi= 0)
    {
      register table_type yi;
      evaluate(xi, x[0],x[1], y[0],y[1], y2[0],y2[1], &yi,dyi,d2yi);
      return yi;
    }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class WDutils::spline                                                      
  //                                                                            
  /// a cubic spline for arbitrary x- and y-type                                
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  template<class scalar_type, class table_type = scalar_type>
  class spline : private cubic_splines<scalar_type,table_type>
  {
  private:
    //--------------------------------------------------------------------------
    /// \name data (private)
    //@{
    const int          n;                         ///< size of tables
    const scalar_type *x;                         ///< ordered table of points
    const table_type  *y;                         ///< table of function values
    table_type  *const y2;                        ///< table of y''
    mutable int        lo;                        ///< x[lo] <= xi <= x[lo+1]
    //@}
  public:
    //--------------------------------------------------------------------------
    /// \name construction: establish a cubic spline
    //{@
    /// constructor from arrays
    /// \param _n size of tables
    /// \param _x table of points
    /// \param _y table of function values
    /// \param yp1 (optional) first derivative at _x[0]
    /// \param ypn (optional) first derivative at _x[_n-1]
    /// If null pointers are given for the first derivative of y at _x[0] and
    /// _x[_n-1] (which is the default), natural boundary conditions are used.
    spline(int _n,
	   const scalar_type*_x,
	   const table_type *_y,
	   const table_type *yp1=0,
	   const table_type *ypn=0)
      : n(_n), x(_x), y(_y), y2(WDutils_NEW(table_type,n)), lo(0)
    {
      construct(n,x,y,y2,yp1,ypn);
    }
    //--------------------------------------------------------------------------
    /// constructor from Array<>s
    /// \param _x Array<> with points
    /// \param _y Array<> with function values
    /// \param yp1 (optional) first derivative at _x[0]
    /// \param ypn (optional) first derivative at _x[_n-1]
    /// If null pointers are given for the first derivative of y at _x[0] and
    /// _x[_n-1] (which is the default), natural boundary conditions are used.
    spline(Array<scalar_type,1> const&_x,
	   Array<table_type ,1> const&_y,
	   const table_type *yp1=0,
	   const table_type *ypn=0)
      : n(_x.size()), x(_x.array()), y(_y.array()),
	y2(WDutils_NEW(table_type,n)), lo(0)
    {
      if(_x.size() != _y.size())
	error("size mismatch in spline constructions\n");
      construct(n,x,y,y2,yp1,ypn);
    }
    //@}
    //--------------------------------------------------------------------------
    /// destructor
    ~spline() { WDutils_DEL_A(y2); }
    //--------------------------------------------------------------------------
    /// spline evaluation
    /// \return y(xi)
    /// \param xi point at which to compute spline
    /// \param dy (optional output) y'(xi)
    /// \param d2y (optional output) y''(xi)
    table_type operator() (scalar_type xi,
			   table_type *dy = 0,
			   table_type *d2y= 0) const
    {
      find(lo,n,x,xi);
      if(lo==n-1) --lo;
      return evaluate(xi,x+lo,y+lo,y2+lo,dy,d2y);
    }
    //--------------------------------------------------------------------------
    /// \name direct constant access to data
    //{@
    /// size of tables
    int         const& N()            const { return n; }
    /// points: x_i
    scalar_type const& X(const int i) const { return x[i]; } 
    /// values: y_i
    table_type  const& Y(const int i) const { return y[i]; } 
    /// first point
    scalar_type const& first_X()      const { return x[0]; } 
    /// first value
    table_type  const& first_Y()      const { return y[0]; } 
    /// last point
    scalar_type const& last_X ()      const { return x[n-1]; } 
    /// last value
    table_type  const& last_Y ()      const { return y[n-1]; } 
    //@}
  };
  //////////////////////////////////////////////////////////////////////////////
  template<class scalar_type, class table_type>
  struct traits<spline<scalar_type,table_type> > {
    static const char  *name () {
      return message("spline<%s,%s>",nameof(scalar_type),nameof(table_type));
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class WDutils::penta_splines                                               
  //                                                                            
  // this is a mere collection of (static) methods                              
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class penta_splines {
  public:
    //--------------------------------------------------------------------------
    // construct 5th order spline                                               
    //--------------------------------------------------------------------------
    template<class scalar_type, class table_type>
    static void construct(int         const n,	  // I:  size of tables         
			  const scalar_type*x,	  // I:  table: points          
			  const table_type *y,    // I:  table: f(x)            
			  const table_type *y1,   // I:  table: df/dx           
			  table_type *const y3)   // O:  table: d^3f/dx^3       
    {
      scalar_type p,sig,dx,dx1,dx2;
      table_type  dy=y[1]-y[0], dy1=dy;
      scalar_type *v = WDutils_NEW(scalar_type,n-1);
      dx   = x[1]-x[0];
      y3[0]= v[0] = 0;
      for(int i=1; i<n-1; i++) {
	dx1  = x[i+1]-x[i];
	dx2  = x[i+1]-x[i-1];
	dy1  = y[i+1]-y[i];
	sig  = dx/dx2;
	p    = sig*v[i-1]-3;
	v[i] = (sig-1)/p;
	y3[i]= times<12>( times< 7>( y1[i]*dx2/(dx*dx1) ) +
			  times< 3>( y1[i-1]/dx+y1[i+1]/dx1 ) -
			  times<10>( dy/(dx*dx) + dy1/(dx1*dx1) ) ) / dx2;
	y3[i]= (y3[i] - sig*y3[i-1] ) / p;
	dx   = dx1;
	dy   = dy1;
      }
      y3[n-1] = table_type(0);
      for(int i=n-2; i>=0; i--)
	y3[i] += v[i]*y3[i+1];
      WDutils_DEL_A(v);
    }
    //--------------------------------------------------------------------------
    // evaluate spline at x=xi with xl <= xi <= xh and yl=y(xl) etc...          
    //--------------------------------------------------------------------------
    template<class scalar_type, class table_type>
    static void evaluate(scalar_type const &xi,
			 scalar_type const &x0,
			 scalar_type const &x1,
			 table_type  const &y0,
			 table_type  const &y1,
			 table_type  const &y10,
			 table_type  const &y11,
			 table_type  const &y30,
			 table_type  const &y31,
			 table_type  *yi,
			 table_type  *dyi = 0,
			 table_type  *d2yi= 0)
    {
      const scalar_type ife=1./48.;
      scalar_type h =x1-x0;
      if(h==0) error("penta_splines::evaluate(): bad x input");
      scalar_type
	hi = scalar_type(1)/h,
	hf = h*h,
	A  = hi*(x1-xi), Aq = A*A,
	B  = 1-A,        Bq = B*B,
	C  = h*Aq*B,
	D  =-h*Bq*A;
      table_type
 	t1 = hi*(y1-y0),
	C2 = y10-t1,
	C3 = y11-t1,
	t2 = times<6>(y10+y11-t1-t1)/hf,
	C4 = y30-t2,
	C5 = y31-t2;
      hf *= ife;
      *yi= A*y0+ B*y1+ C*C2+ D*C3+ hf*(C*(Aq+Aq-A-1)*C4+ D*(Bq+Bq-B-1)*C5);
      if(dyi) {
	scalar_type BAA=B-A-A, ABB=A-B-B;
	hf  += hf;
	*dyi = t1 + (A*ABB)*C2 + (B*BAA)*C3
	     + hf*A*B*((1+A-times<5>(Aq))*C4+ (1+B-times<5>(Bq))*C5);
	if(d2yi) {
	  *d2yi = BAA*C2 - ABB*C3;
	  *d2yi+= *d2yi  + hf * ( (twice(Aq)*(times<9>(B)-A)-1) * C4 +
				  (twice(Bq)*(B-times<9>(A))+1) * C5 );
	  *d2yi*= hi;
	}
      }
    }
    //--------------------------------------------------------------------------
    // evaluate penta spline at x=xi                                            
    //--------------------------------------------------------------------------
    template<class scalar_type, class table_type>
    static table_type evaluate(scalar_type const &xi,
			       const scalar_type *x,
			       const table_type  *y,
			       const table_type  *y1,
			       const table_type  *y3,
			       table_type  *dyi = 0,
			       table_type  *d2yi= 0)
    {
      table_type yi;
      evaluate(xi,x[0],x[1],y[0],y[1],y1[0],y1[1],y3[0],y3[1],&yi,dyi,d2yi);
      return yi;
    }
    //--------------------------------------------------------------------------
    // evaluate penta spline at x=xi                                            
    //--------------------------------------------------------------------------
    template<class scalar_type, class table_type>
    static table_type eval(scalar_type const &xi,
			   const scalar_type *x,
			   const table_type  *y,
			   const table_type  *y1,
			   const table_type  *y3,
			   table_type  *dyi = 0,
			   table_type  *d2yi= 0)
    {
      table_type yi;
      evaluate(xi,x[0],x[1],y[0],y[1],y1[0],y1[1],y3[0],y3[1],&yi,dyi,d2yi);
      return yi;
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class WDutils::Pspline                                                   //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<class scalar_type, class table_type = scalar_type>
  class Pspline : private penta_splines
  {
  private:
    //--------------------------------------------------------------------------
    // data of class WDutils::spline                                            
    //--------------------------------------------------------------------------
    const int          n;                         // size of tables             
    const scalar_type *x;                         // ordered table of points    
    const table_type  *y;                         // table of function values   
    const table_type  *y1;                        // table of 1st derivatives   
    table_type  *const y3;                        // table of 3rd derivatives   
    mutable int        lo;                        // x[lo] <= xi <= x[lo+1]     
  public:
    //--------------------------------------------------------------------------
    // constructor: establish a cubic spline.                                   
    //              if zero pointers are given for the first derivative of y    
    //              at x[0] and x[n-1], natural boundary conditions are used    
    //--------------------------------------------------------------------------
    Pspline(int         const _n,                 // I: size of tables          
	    const scalar_type*_x,                 // I: ordered table of points 
	    const table_type *_y,                 // I: table of y              
	    const table_type *_y1)                // I: table of y'             
      : n(_n), x(_x), y(_y), y1(_y1), y3(WDutils_NEW(table_type,n)), lo(0)
    {
      construct(n,x,y,y1,y3);
    }
    //--------------------------------------------------------------------------
    // destructor                                                               
    //--------------------------------------------------------------------------
    ~Pspline() { WDutils_DEL_A(y3); }
    //--------------------------------------------------------------------------
    // spline evaluation                                                        
    //--------------------------------------------------------------------------
    table_type operator() (                       // R: y(xi)                   
			   const scalar_type xi,  // I: xi                      
			   table_type *dy = 0,    //[O]: dy/dx     at x=xi      
			   table_type *d2y= 0)    //[O]: d^2y/dx^2 at x=xi      
      const
    {
      find(lo,n,x,xi);
      if(lo==n-1) --lo;
      return evaluate(xi,x+lo,y+lo,y1+lo,y3+lo,dy,d2y);
    }
    //--------------------------------------------------------------------------
    // access to tables                                                         
    //--------------------------------------------------------------------------
    int         const& N       ()      const { return n; }
    scalar_type const& X       (int i) const { return x [i]; } 
    table_type  const& Y       (int i) const { return y [i]; } 
    table_type  const& Y1      (int i) const { return y1[i]; } 
    scalar_type const& first_X ()      const { return x [0]; } 
    table_type  const& first_Y ()      const { return y [0]; } 
    table_type  const& first_Y1()      const { return y1[0]; } 
    scalar_type const& last_X  ()      const { return x [n-1]; } 
    table_type  const& last_Y  ()      const { return y [n-1]; } 
    table_type  const& last_Y1 ()      const { return y1[n-1]; } 
  };
}
////////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_spln_h
