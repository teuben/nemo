// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// vect.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 1999-2002                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef included_vect_h
#define included_vect_h

#ifndef included_iostream
#  include <iostream>
#  define included_iostream
#endif

#ifndef included_inln_h
#  include <public/inln.h>
#endif

#ifndef NDIM
#  error NDIM not defined in __FILE__
#endif
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
#if NDIM == 2

#define CI(op,b,ox) A[0] op b[0] ox  A[1] op b[1]
#define CX(op,b,ox) A[0] op b    ox  A[1] op b
#define FI(op,b)    A[0] op b[0];    A[1] op b[1]
#define FX(op,b)    A[0] op b;       A[1] op b

#else

#define CI(op,b,ox) A[0] op b[0] ox  A[1] op b[1] ox A[2] op b[2]
#define CX(op,b,ox) A[0] op b    ox  A[1] op b    ox A[2] op b
#define FI(op,b)    A[0] op b[0];    A[1] op b[1];   A[2] op b[2]
#define FX(op,b)    A[0] op b;       A[1] op b;      A[2] op b

#endif

#define OP      operator
#define TS      template<typename SCALAR>
#define   VECT        vector<SCALAR>
#define c_SCAL  const SCALAR&
#define c_VECT  const vector<SCALAR>&
#define c_EVEC  const sym1  <SCALAR>&
  //////////////////////////////////////////////////////////////////////////////
  // some forward declarations                                                  
  //////////////////////////////////////////////////////////////////////////////
  template<typename REAL> class vector;
  template<typename REAL> class sym1;
  template<typename REAL> class sym2;
  template<typename REAL> class sym3;
  template<typename REAL> class sym4;
  template<typename REAL> std::ostream& OP<< (std::ostream&, const vector<REAL>&);
  template<typename REAL> std::istream& OP>> (std::istream&,       vector<REAL>&);
  //////////////////////////////////////////////////////////////////////////////
  // class nbdy::vector                                                         
  //////////////////////////////////////////////////////////////////////////////
  template<typename REAL>                      // scalar type                   
  class vector {                               // still over templated          
  public:
    static const int NDAT = NDIM;
    //--------------------------------------------------------------------------
    // friendships                                                              
    //--------------------------------------------------------------------------
    friend class sym1<REAL>;
    friend class sym2<REAL>;
    friend class sym3<REAL>;
    friend class sym4<REAL>;
    friend std::ostream& OP<< <>(std::ostream&, const vector<REAL>&);
    friend std::istream& OP>> <>(std::istream&,       vector<REAL>&);
    //--------------------------------------------------------------------------
    // types                                                                    
    //--------------------------------------------------------------------------
  public:
    typedef REAL                 element_type; // type of elements              
  private:
    typedef REAL                 real;         // type for const returns        
    typedef vector               vect;         // type for const returns        
    typedef const int          c_indx;         // type for const args           
    typedef const REAL        &c_real;         // type for const args           
    typedef const vect        &c_vect;         // type for const args           
    typedef const sym1<REAL>  &c_sym1;         // type for const args           
    //--------------------------------------------------------------------------
    // data                                                                     
    //--------------------------------------------------------------------------
  private:
    REAL A[NDIM];
    //--------------------------------------------------------------------------
    // construction                                                             
    //--------------------------------------------------------------------------
  public:
       vector ()         {}
       vector (c_real x) { FX(=,x); }
       vector (c_vect x) { FI(=,x.A); }
       vector (c_sym1 x) { FI(=,x.A); }
    TS vector (c_VECT x) { FI(=,x); }
    TS vector (c_EVEC x) { FI(=,x); }
    //--------------------------------------------------------------------------
    // non-const methods                                                        
    //--------------------------------------------------------------------------
       real* pointer()                   { return A; }
       real &OP[]   (c_indx i)           { return A[i]; }
       vect &OP=    (c_real x)           { FX(=,x); return *this; }
       vect &OP=    (const real* x)      { FI(=,x); return *this; }
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
    TS vect &OP/=   (c_SCAL x)           { return OP*=(real(1)/x); }
       vect &negate ()                   { FI(=-,A); return *this; }
    TS vect &add_mul(c_VECT y, c_real x) { FI(+=,x*y.A); return *this; }
    TS vect &sub_mul(c_VECT y, c_real x) { FI(-=,x*y.A); return *this; }
    TS vect &add_mul(c_EVEC y, c_real x) { FI(+=,x*y.A); return *this; }
    TS vect &sub_mul(c_EVEC y, c_real x) { FI(-=,x*y.A); return *this; }
    //--------------------------------------------------------------------------
    TS vect &add_twice(c_VECT y)         { A[0] += y.A[0]+y.A[0];
                                           A[1] += y.A[1]+y.A[1];
#if NDIM==3
                                           A[2] += y.A[2]+y.A[2];
#endif
					   return *this; }
    TS vect &add_twice(c_EVEC y)         { A[0] += y.A[0]+y.A[0];
                                           A[1] += y.A[1]+y.A[1];
#if NDIM==3
                                           A[2] += y.A[2]+y.A[2];
#endif
					   return *this; }
    TS vect &sub_twice(c_VECT y)         { A[0] -= y.A[0]+y.A[0];
                                           A[1] -= y.A[1]+y.A[1];
#if NDIM==3
                                           A[2] -= y.A[2]+y.A[2];
#endif
					   return *this; }
    TS vect &sub_twice(c_EVEC y)         { A[0] -= y.A[0]+y.A[0];
                                           A[1] -= y.A[1]+y.A[1];
#if NDIM==3
                                           A[2] -= y.A[2]+y.A[2];
#endif
					   return *this; }
    //--------------------------------------------------------------------------
    vect &apply  (real(*f)(c_real))   { 
      A[0]=f(A[0]); A[1]=f(A[1]); 
#if NDIM==3
      A[2]=f(A[2]);
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    TS vect &connect(c_VECT x, real(*f)(c_real, c_real)) {
      A[0]=f(A[0],x.A[0]); A[1]=f(A[1],x.A[1]);
#if NDIM==3
      A[2]=f(A[2],x.A[2]);
#endif
      return *this;
    }
    //--------------------------------------------------------------------------
    // const methods                                                            
    //--------------------------------------------------------------------------
    const real* const_pointer() const { return A; }
    real const&OP[](c_indx i) const { return A[i]; }
       vect OP*    (c_real x) const { register vect y(*this); y*=x; return y; }
       vect OP/    (c_real x) const { register vect y(*this); y*=1/x;return y; }
    TS vect OP+    (c_VECT x) const { register vect y(*this); y+=x; return y; }
    TS vect OP+    (c_EVEC x) const { register vect y(*this); y+=x; return y; }
    TS vect OP-    (c_VECT x) const { register vect y(*this); y-=x; return y; }
    TS vect OP-    (c_EVEC x) const { register vect y(*this); y-=x; return y; }
       vect OP-    ()         const { register vect y(0); y-= *this; return y; }
       real norm   ()         const { return CI(*,A,+); }
    TS real OP*    (c_VECT x) const { return CI(*,x.A,+); }
    TS real OP*    (c_EVEC x) const { return CI(*,x.A,+); }
    //--------------------------------------------------------------------------
#if NDIM == 2
    real maxnorm   ()         const { return max(abs(A[0]), abs(A[1])); }
    real min       ()         const { return min(A[0],A[1]); }
    real max       ()         const { return max(A[0],A[1]); }
    TS real dist_sq(c_VECT x) const { return square(A[0]-x.A[0])
				           + square(A[1]-x.A[1]); }
    TS real dist_sq(c_EVEC x) const { return square(A[0]-x.A[0])
				           + square(A[1]-x.A[1]); }
    TS real sum_sq (c_VECT x) const { return square(A[0]+x.A[0])
				           + square(A[1]+x.A[1]); }
    TS real sum_sq (c_EVEC x) const { return square(A[0]+x.A[0])
				           + square(A[1]+x.A[1]); }
    TS real OP^    (c_VECT x) const { return A[0]*x.A[1] - A[1]*x.A[0]; }
    TS real OP^    (c_EVEC x) const { return A[0]*x.A[1] - A[1]*x.A[0]; }
#else
    //==========================================================================
    real maxnorm   ()         const { return max(abs(A[0]), abs(A[1]),
						 abs(A[2])); }
    real min       ()         const { return min(A[0],A[1],A[2]); }
    real max       ()         const { return max(A[0],A[1],A[2]); }
    TS real dist_sq(c_VECT x) const { return square(A[0]-x.A[0])
				           + square(A[1]-x.A[1])
				           + square(A[2]-x.A[2]); }
    TS real dist_sq(c_EVEC x) const { return square(A[0]-x.A[0])
				           + square(A[1]-x.A[1])
				           + square(A[2]-x.A[2]); }
    TS real sum_sq (c_VECT x) const { return square(A[0]+x.A[0])
				           + square(A[1]+x.A[1])
				           + square(A[2]+x.A[2]); }
    TS real sum_sq (c_EVEC x) const { return square(A[0]+x.A[0])
				           + square(A[1]+x.A[1])
				           + square(A[2]+x.A[2]); }
    TS vect OP^    (c_VECT x) const {
      register vect z;
      z.A[0] = A[1]*x.A[2] - A[2]*x.A[1];
      z.A[1] = A[2]*x.A[0] - A[0]*x.A[2];
      z.A[2] = A[0]*x.A[1] - A[1]*x.A[0];
      return z;
    }
    TS vect OP^    (c_EVEC x) const {
      register vect z;
      z.A[0] = A[1]*x.A[2] - A[2]*x.A[1];
      z.A[1] = A[2]*x.A[0] - A[0]*x.A[2];
      z.A[2] = A[0]*x.A[1] - A[1]*x.A[0];
      return z;
    }
#endif
    //--------------------------------------------------------------------------
       bool OP==(c_real x) const { return CX(==,x,&&); }
    TS bool OP==(c_VECT x) const { return CI(==,x.A,&&); }
    TS bool OP==(c_EVEC x) const { return CI(==,x.A,&&); }
       bool OP!=(c_real x) const { return CX(!=,x,||); }
    TS bool OP!=(c_VECT x) const { return CI(!=,x.A,||); }
    TS bool OP!=(c_EVEC x) const { return CI(!=,x.A,||); }
  };
  //////////////////////////////////////////////////////////////////////////////
  // related functions                                                          
  //////////////////////////////////////////////////////////////////////////////
  template<typename REAL> inline
  vector<REAL> OP*(const REAL&x, const vector<REAL>& y) { return y*x; }
  //----------------------------------------------------------------------------
  template<typename REAL> inline
  REAL norm(const vector<REAL>& x) { return x.norm(); }
  //----------------------------------------------------------------------------
  template<typename REAL> inline
  REAL dist_sq(const vector<REAL>&x,const vector<REAL>&y){ return x.dist_sq(y); }
  //----------------------------------------------------------------------------
  template<typename REAL> inline
  REAL sum_sq(const vector<REAL>&x,const vector<REAL>&y){ return x.sum_sq(y); }
  //----------------------------------------------------------------------------
  template<typename REAL> inline
  std::ostream& OP<<(std::ostream&s, const vector<REAL>& x) {
    return write_array(s,x.A,NDIM);
  }
  //----------------------------------------------------------------------------
  template<typename REAL> inline
  std::istream& OP>>(std::istream&s, vector<REAL>& x) {
    return read_array(s,x.A,NDIM);
  }
  //----------------------------------------------------------------------------
}
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
#endif	// included_vect_h
