// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// vect.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 1999-2003                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// last change 28/07/03: happy gcc 3.3                                         |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_vect_h
#define falcON_included_vect_h

#ifndef falcON_included_iostream
#  include <iostream>
#  define falcON_included_iostream
#endif

#ifndef falcON_included_inln_h
#  include <public/inln.h>
#endif

#ifndef falcON_NDIM
#  error falcON_NDIM not defined in vect.h
#endif
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
#if falcON_NDIM == 2

#define CI(op,b,ox) A[0] op b[0] ox  A[1] op b[1]
#define CX(op,b,ox) A[0] op b    ox  A[1] op b
#define FI(op,b)    A[0] op b[0];    A[1] op b[1]
#define FX(op,b)    A[0] op b;       A[1] op b
#define DO(macro)   macro(0); macro(1)

#else

#define CI(op,b,ox) A[0] op b[0] ox  A[1] op b[1] ox A[2] op b[2]
#define CX(op,b,ox) A[0] op b    ox  A[1] op b    ox A[2] op b
#define FI(op,b)    A[0] op b[0];    A[1] op b[1];   A[2] op b[2]
#define FX(op,b)    A[0] op b;       A[1] op b;      A[2] op b
#define DO(macro)   macro(0); macro(1); macro(2)

#endif

#define OP      operator
#define TS      template<typename SCAL>
#define   VECT        vector<SCAL>
#define c_SCAL  const SCAL&
#define c_VECT  const vector<SCAL>&
#define c_EVEC  const sym1  <SCAL>&
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // some forward declarations                                                //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  TS class vector;
  TS class sym1;
  TS class sym2;
  TS class sym3;
  TS class sym4;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::vector                                                       //
  //                                                                          //
  //--------------------------------------------------------------------------//
  //                                                                          //
  // a tupel with nbdy::Ndim elements                                         //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename SCALAR>                    // scalar type                   
  class vector {
  public:
    static const int NDAT = Ndim;
    //--------------------------------------------------------------------------
    // types                                                                    
    //--------------------------------------------------------------------------
  public:
    typedef SCALAR               element_type; // type of elements              
  private:
    typedef element_type         scal;         // type for const returns        
    typedef vector               vect;         // type for const returns        
    typedef int          const&c_indx;         // type for const args           
    typedef scal         const&c_scal;         // type for const args           
    typedef vect         const&c_vect;         // type for const args           
    typedef sym1<scal>   const&c_sym1;         // type for const args           
    //--------------------------------------------------------------------------
    // friendships                                                              
    //--------------------------------------------------------------------------
    friend class sym1<element_type>;
    friend class sym2<element_type>;
    friend class sym3<element_type>;
    friend class sym4<element_type>;
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
  private:
    element_type A[Ndim];
    //--------------------------------------------------------------------------
    // construction                                                             
    //--------------------------------------------------------------------------
  public:
       vector ()         {}
       vector (c_scal x) { FX(=,x); }
       vector (c_vect x) { FI(=,x.A); }
       vector (c_sym1 x) { FI(=,x.A); }
    TS vector (c_VECT x) { FI(=,x); }
    TS vector (c_EVEC x) { FI(=,x); }
    //--------------------------------------------------------------------------
    // non-const methods                                                        
    //--------------------------------------------------------------------------
       scal* pointer()                   { return A; }
       scal &OP[]   (c_indx i)           { return A[i]; }
       vect &OP=    (c_scal x)           { FX(=,x); return *this; }
       vect &OP=    (const scal* x)      { FI(=,x); return *this; }
       vect &OP=    (c_vect x)           { FI(=,x.A); return *this; }
       vect &OP=    (c_sym1 x)           { FI(=,x.A); return *this; }
    TS vect &OP=    (c_VECT x)           { FI(=,x); return *this; }
    TS vect &OP=    (c_EVEC x)           { FI(=,x); return *this; }
       vect &OP+=   (c_vect x)           { FI(+=,x.A); return *this; }
       vect &OP+=   (c_sym1 x)           { FI(+=,x.A); return *this; }
    TS vect &OP+=   (c_VECT x)           { FI(+=,x); return *this; }
    TS vect &OP+=   (c_EVEC x)           { FI(+=,x); return *this; }
       vect &OP-=   (c_vect x)           { FI(-=,x.A); return *this; }
       vect &OP-=   (c_sym1 x)           { FI(-=,x.A); return *this; }
    TS vect &OP-=   (c_VECT x)           { FI(-=,x); return *this; }
    TS vect &OP-=   (c_EVEC x)           { FI(-=,x); return *this; }
    TS vect &OP*=   (c_SCAL x)           { FX(*=,x); return *this; }
    TS vect &OP/=   (c_SCAL x)           { return OP*=(scal(1)/x); }
       vect &negate ()                   { FI(=-,A); return *this; }
    TS vect &add_mul(c_VECT y, c_scal x) { FI(+=,x*y.A); return *this; }
    TS vect &sub_mul(c_VECT y, c_scal x) { FI(-=,x*y.A); return *this; }
    TS vect &add_mul(c_EVEC y, c_scal x) { FI(+=,x*y.A); return *this; }
    TS vect &sub_mul(c_EVEC y, c_scal x) { FI(-=,x*y.A); return *this; }
    //--------------------------------------------------------------------------
#define JOB(i) A[i] += y.A[i]+y.A[i]
    TS vect &add_twice(c_VECT y)         { DO(JOB); return *this; }
    TS vect &add_twice(c_EVEC y)         { DO(JOB); return *this; }
#undef  JOB
#define JOB(i) A[i] -= y.A[i]+y.A[i]
    TS vect &sub_twice(c_VECT y)         { DO(JOB); return *this; }
    TS vect &sub_twice(c_EVEC y)         { DO(JOB); return *this; }
#undef  JOB
    //--------------------------------------------------------------------------
#define JOB(i) A[i]=f(A[i])
    vect &apply  (scal(*f)(c_scal))   { DO(JOB); return *this; }
#undef  JOB
    //--------------------------------------------------------------------------
#define JOB(i) A[i]=f(A[i],x.A[i])
    TS vect &connect(c_VECT x, scal(*f)(c_scal, c_scal))
    { DO(JOB); return *this; }
#undef  JOB
    //--------------------------------------------------------------------------
    // const methods                                                            
    //--------------------------------------------------------------------------
    const scal* const_pointer() const { return A; }
    scal const&OP[](c_indx i) const { return A[i]; }
       vect OP*    (c_scal x) const { register vect y(*this); y*=x; return y; }
    friend
       vect OP*    (c_scal x,
		    c_vect y)       { return y*x; }
       vect OP/    (c_scal x) const { register vect y(*this); y*=1/x;return y; }
    TS vect OP+    (c_VECT x) const { register vect y(*this); y+=x; return y; }
    TS vect OP+    (c_EVEC x) const { register vect y(*this); y+=x; return y; }
    TS vect OP-    (c_VECT x) const { register vect y(*this); y-=x; return y; }
    TS vect OP-    (c_EVEC x) const { register vect y(*this); y-=x; return y; }
       vect OP-    ()         const { register vect y(0); y-= *this; return y; }
       scal norm   ()         const { return CI(*,A,+); }
    friend
       scal norm   (c_vect x)       { return x.norm(); }
    TS scal OP*    (c_VECT x) const { return CI(*,x.A,+); }
    TS scal OP*    (c_EVEC x) const { return CI(*,x.A,+); }
    //--------------------------------------------------------------------------
#if falcON_NDIM == 2
    scal maxnorm   ()         const { return max(abs(A[0]), abs(A[1])); }
    scal min       ()         const { return min(A[0],A[1]); }
    scal max       ()         const { return max(A[0],A[1]); }
    TS scal dist_sq(c_VECT x) const { return square(A[0]-x.A[0])
				           + square(A[1]-x.A[1]); }
    TS scal dist_sq(c_EVEC x) const { return square(A[0]-x.A[0])
				           + square(A[1]-x.A[1]); }
    TS scal sum_sq (c_VECT x) const { return square(A[0]+x.A[0])
				           + square(A[1]+x.A[1]); }
    TS scal sum_sq (c_EVEC x) const { return square(A[0]+x.A[0])
				           + square(A[1]+x.A[1]); }
    TS scal OP^    (c_VECT x) const { return A[0]*x[1] - A[1]*x[0]; }
    TS scal OP^    (c_EVEC x) const { return A[0]*x[1] - A[1]*x[0]; }
#else
    //==========================================================================
    scal maxnorm   ()         const { return max(abs(A[0]), abs(A[1]),
						 abs(A[2])); }
    scal min       ()         const { return min(A[0],A[1],A[2]); }
    scal max       ()         const { return max(A[0],A[1],A[2]); }
    TS scal dist_sq(c_VECT x) const { return square(A[0]-x.A[0])
				           + square(A[1]-x.A[1])
				           + square(A[2]-x.A[2]); }
    TS scal dist_sq(c_EVEC x) const { return square(A[0]-x.A[0])
				           + square(A[1]-x.A[1])
				           + square(A[2]-x.A[2]); }
    TS scal sum_sq (c_VECT x) const { return square(A[0]+x.A[0])
				           + square(A[1]+x.A[1])
				           + square(A[2]+x.A[2]); }
    TS scal sum_sq (c_EVEC x) const { return square(A[0]+x.A[0])
				           + square(A[1]+x.A[1])
				           + square(A[2]+x.A[2]); }
    TS vect OP^    (c_VECT x) const {
      register vect z;
      z.A[0] = A[1]*x[2] - A[2]*x[1];
      z.A[1] = A[2]*x[0] - A[0]*x[2];
      z.A[2] = A[0]*x[1] - A[1]*x[0];
      return z;
    }
    TS vect OP^    (c_EVEC x) const {
      register vect z;
      z.A[0] = A[1]*x[2] - A[2]*x[1];
      z.A[1] = A[2]*x[0] - A[0]*x[2];
      z.A[2] = A[0]*x[1] - A[1]*x[0];
      return z;
    }
#endif
    TS scal dist(c_VECT x) const { return sqrt(dist_sq(x)); }
    TS scal dist(c_EVEC x) const { return sqrt(dist_sq(x)); }
    //--------------------------------------------------------------------------
    TS friend scal dist_sq(c_vect x, c_VECT y) { return x.dist_sq(y); }
    TS friend scal dist_sq(c_vect x, c_EVEC y) { return x.dist_sq(y); }
    TS friend scal sum_sq (c_vect x, c_VECT y) { return x.sum_sq(y); }
    TS friend scal sum_sq (c_vect x, c_EVEC y) { return x.sum_sq(y); }
    TS friend scal dist   (c_vect x, c_VECT y) { return x.dist(y); }
    TS friend scal dist   (c_vect x, c_EVEC y) { return x.dist(y); }
    //--------------------------------------------------------------------------
       bool OP==(c_scal x) const { return CX(==,x,&&); }
    TS bool OP==(c_VECT x) const { return CI(==,x.A,&&); }
    TS bool OP==(c_EVEC x) const { return CI(==,x.A,&&); }
       bool OP!=(c_scal x) const { return CX(!=,x,||); }
    TS bool OP!=(c_VECT x) const { return CI(!=,x.A,||); }
    TS bool OP!=(c_EVEC x) const { return CI(!=,x.A,||); }
    //--------------------------------------------------------------------------
    // ascii I/O                                                                
    //--------------------------------------------------------------------------
    friend std::ostream& OP<< (std::ostream&s, c_vect x)
    {
      return write_array(s,x.A,NDAT);
    }
    //--------------------------------------------------------------------------
    friend std::istream& OP>> (std::istream&s, vect &x)
    {
      return read_array(s,x.A,NDAT);
    }
  };
}
#undef DO
#undef CI
#undef CX
#undef FI
#undef FX
#undef OP
#undef TS
#undef VECT
#undef c_SCAL
#undef c_EVEC
#undef c_VECT
////////////////////////////////////////////////////////////////////////////////
#endif	// falcON_included_vect_h
