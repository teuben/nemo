// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// tupl.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 1996-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// defines:                                                                    |
//                                                                             |
// class tupel<N,scalar>                                                       |
// class pseudo_tupel<N,scalar>                                                |
// class const_pseudo_tupel<N,scalar>                                          |
//                                                                             |
// new design aug-2003: template metaprogramming to unroll loops automatically |
//            nov-2003: non-standard #include files redundant; pseudo_tupel    |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_tupl_h
#define falcON_included_tupl_h

#ifndef falcON_included_iostream
#  include <iostream>
#  define falcON_included_iostream
#endif
#ifndef falcON_included_frst_h
#  include <public/frst.h>
#endif
#ifndef falcON_included_meta_h
#  include <public/meta.h>
#endif

////////////////////////////////////////////////////////////////////////////////
#define i_      inline
#define c_      const
#define op      operator
#define tS      template<typename S>
#define tP      template<typename P>
#define rt      return*this
#define cS      const S
#define tX      template<typename X>
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  template<int, typename> class tupel;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::pseudo_tupel                                                 //
  //                                                                          //
  // encodes the numbers A[d][I] with d=0...N-1 like a tupel                  //
  // used to access non-const data in nbdy::abodies                           //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
#define W       tupel<N,S>
#define cW      const tupel<N,S>
  template<int N, typename X> class pseudo_tupel {
    typedef X      *pX;
    typedef const X cX;
    //--------------------------------------------------------------------------
    pX const *a;
    const int I;
  public:
    //--------------------------------------------------------------------------
    pseudo_tupel(pX const*const&A, int const&i) : a(A), I(i) {}
    const X&op[](int const&d) const { return a[d][I]; }
          X&op[](int const&d)       { return a[d][I]; }
    //--------------------------------------------------------------------------
    // operations with scalar                                                   
    //--------------------------------------------------------------------------
    tS pseudo_tupel& op = (cS&x) { meta::taux<S,N-1,0>::s_app(*this,x); rt; }
    tS pseudo_tupel& op+= (cS&x) { meta::taux<S,N-1,0>::s_apa(*this,x); rt; }
    tS pseudo_tupel& op-= (cS&x) { meta::taux<S,N-1,0>::s_aps(*this,x); rt; }
    //--------------------------------------------------------------------------
    // operations with tupel<N,S>                                               
    //--------------------------------------------------------------------------
    tS pseudo_tupel& op = (cW&x) {
      meta::taux<S,N-1,0>::v_app(*this,static_cast<cS*>(x)); rt; }
    tS pseudo_tupel& op+= (cW&x) {
      meta::taux<S,N-1,0>::v_apa(*this,static_cast<cS*>(x)); rt; }
    tS pseudo_tupel& op-= (cW&x) {
      meta::taux<S,N-1,0>::v_aps(*this,static_cast<cS*>(x)); rt; }
    tS pseudo_tupel&ass_times(cW&x, cS&f) {
      meta::taux<S,N-1,0>::v_appt(*this,static_cast<cS*>(x),f); rt; }
    tS pseudo_tupel&add_times(cW&x, cS&f) {
      meta::taux<S,N-1,0>::v_apat(*this,static_cast<cS*>(x),f); rt; }
    tS pseudo_tupel&sub_times(cW&x, cS&f) {
      meta::taux<S,N-1,0>::v_apst(*this,static_cast<cS*>(x),f); rt; }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::const_pseudo_tupel                                           //
  //                                                                          //
  // encodes the numbers A[d][I] with d=0...N-1 like a tupel                  //
  // used to access const data in nbdy::abodies                               //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<int N, typename X> class const_pseudo_tupel {
    typedef const X *pX;
    pX const *a;
    const int I;
  public:
    const_pseudo_tupel(pX const*const&A, int const&i) : a(A), I(i) {}
    const X&op[](int const&d) const { return a[d][I]; }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::tupel                                                        //
  //                                                                          //
  // a tupel of N scalars of type X, held in an array X[N]                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<int N, typename X> class tupel {
    // private types and data                                                   
  private:
    typedef const int               cI;
    typedef const X                 cX;
    typedef const X                *pX;
    typedef const tupel             cV;
    typedef pseudo_tupel<N,X>        P;
    typedef const_pseudo_tupel<N,X> cP;
    typedef meta::taux<X,N-1,0>      M;
    X a[N];
  public:
    // static members and public types                                          
    typedef X        element_type;
    typedef tupel    V;
    static const int NDAT = N;
    static const int ORD  = 1;
    static           int size() { return N; }
    // construction                                                             
       tupel      ()        {}
       tupel      (cX&x)    { M::s_as(a,x); }
       tupel      (pX x)    { M::v_as(a,x); }
       tupel      (cV&x)    { M::v_as(a,x.a); }
       tupel      (cP&x)    { M::v_ass(a,x); }
    tS tupel      (cW&x)    { M::v_as(a,(c_ S*)x); }
    tS tupel      (cS*x)    { M::v_as(a,x); }
    // construction from list of elements (up to 10)                            
       tupel      (cX&x0,cX&x1) {
	 a[0]=x0;a[1]=x1; }
       tupel      (cX&x0,cX&x1,cX&x2) {
	 a[0]=x0;a[1]=x1;a[2]=x2; }
       tupel      (cX&x0,cX&x1,cX&x2,cX&x3) {
	 a[0]=x0;a[1]=x1;a[2]=x2;a[3]=x3; }
       tupel      (cX&x0,cX&x1,cX&x2,cX&x3,cX&x4) {
	 a[0]=x0;a[1]=x1;a[2]=x2;a[3]=x3;a[4]=x4; }
       tupel      (cX&x0,cX&x1,cX&x2,cX&x3,cX&x4,cX&x5) { 
	 a[0]=x0;a[1]=x1;a[2]=x2;a[3]=x3;a[4]=x4;a[5]=x5; }
       tupel      (cX&x0,cX&x1,cX&x2,cX&x3,cX&x4,cX&x5,cX&x6) {
	 a[0]=x0;a[1]=x1;a[2]=x2;a[3]=x3;a[4]=x4;a[5]=x5;a[6]=x6; }
       tupel      (cX&x0,cX&x1,cX&x2,cX&x3,cX&x4,cX&x5,cX&x6,cX&x7) {
	 a[0]=x0;a[1]=x1;a[2]=x2;a[3]=x3;a[4]=x4;a[5]=x5;a[6]=x6;a[7]=x7; }
       tupel      (cX&x0,cX&x1,cX&x2,cX&x3,cX&x4,cX&x5,cX&x6,cX&x7,cX&x8) {
	 a[0]=x0;a[1]=x1;a[2]=x2;a[3]=x3;a[4]=x4;a[5]=x5;a[6]=x6;a[7]=x7;
	 a[8]=x8; }
       tupel      (cX&x0,cX&x1,cX&x2,cX&x3,cX&x4,cX&x5,cX&x6,cX&x7,cX&x8,cX&x9){
	 a[0]=x0;a[1]=x1;a[2]=x2;a[3]=x3;a[4]=x4;a[5]=x5;a[6]=x6;a[7]=x7;
	 a[8]=x8;a[9]=x9; } 
    // element access                                                           
       X   &op[]  (cI&i)    { return a[i]; }
       X c_&op[]  (cI&i) c_ { return a[i]; }
    // unitary operators                                                        
       V&negate   ()        { M::v_neg(a); rt; }
       V op-	  ()     c_ { register V y; M::v_nega(y.a,a); return y; }
       X norm     ()     c_ { return M::v_norm(a); }
       X abs      ()     c_ { return sqrt(norm()); }
       X min      ()     c_ { return M::v_min(a); }
       X max      ()     c_ { return M::v_max(a); }
       X maxnorm  ()     c_ { return M::v_amax(a); }
       op    X*   ()        { return a; }
       op c_ X*   ()     c_ { return a; }
    // binary operators with pseudo_tupel<N,X> & const_pseudo_tupel<N,X>        
       V&op=      (cP&x)    { M::v_ass(a,x); rt; }
       V&op=      (cP c_&x) { M::v_ass(a,x); rt; }
       V&op+=     (cP&x)    { M::v_asa(a,x); rt; }
       V&op+=     (cP c_&x) { M::v_asa(a,x); rt; }
       V&op-=     (cP&x)    { M::v_asu(a,x); rt; }
       V&op-=     (cP c_&x) { M::v_asu(a,x); rt; }
       V&op=      ( P&x)    { M::v_ass(a,x); rt; }
       V&op=      ( P c_&x) { M::v_ass(a,x); rt; }
       V&op+=     ( P&x)    { M::v_asa(a,x); rt; }
       V&op+=     ( P c_&x) { M::v_asa(a,x); rt; }
       V&op-=     ( P&x)    { M::v_asu(a,x); rt; }
       V&op-=     ( P c_&x) { M::v_asu(a,x); rt; }
    // binary operators with scalar S                                           
    tS V&op=      (cS&x)    { M::s_as(a,x); rt; }
    tS V&op*=     (cS&x)    { M::s_ml(a,x); rt; }
    tS V&op/=     (cS&x)    { return op*=(S(1)/x); }
    tS V op*      (cS&x) c_ { register V y; M::v_ast(y.a,a,x); return y; }
    tS V op/      (cS&x) c_ { return op*(S(1)/x); }
       bool  op== (cX&x) c_ { return M::s_eq(a,x); }
       bool  op!= (cX&x) c_ { return M::s_neq(a,x); }
    // binary operators with tupel<N,X>                                         
       V&op=      (cV&x)    { M::v_as(a,x.a); rt; }
       V&op+=     (cV&x)    { M::v_ad(a,x.a); rt; }
       V&add_twice(cV&x)    { meta::tauxZ<2,X,N-1,0>::v_ad(a,x.a); rt; }
       V&op-=     (cV&x)    { M::v_su(a,x.a); rt; }
       V&sub_twice(cV&x)    { meta::tauxZ<2,X,N-1,0>::v_su(a,x.a); rt; }
       bool  op== (cV&x) c_ { return M::v_eq(a,x.a); }
       bool  op!= (cV&x) c_ { return M::v_neq(a,x.a); }
       V op+      (cV&x) c_ { register V y; M::v_sum(y.a,a,x.a); return y; }
       V op-      (cV&x) c_ { register V y; M::v_dif(y.a,a,x.a); return y; }
       X op*      (cV&x) c_ { return M::v_dot(a,x.a); }
       X dist_sq  (cV&x) c_ { return M::v_diq(a,x.a); }
       X dist     (cV&x) c_ { return sqrt(dist_sq (x)); }
       X sum_sq   (cV&x) c_ { return M::v_suq(a,x.a);}
       V&up_max   (cV&x)    { M::v_uma(a,x.a); rt; }
       V&up_min   (cV&x)    { M::v_umi(a,x.a); rt; }
    // binary operators with tupel<N,S>                                         
    tS V&op=      (cW&x)    { M::v_as(a,(cS*)x); rt; }
    tS V&op+=     (cW&x)    { M::v_ad(a,(cS*)x); rt; }
    tS V&add_twice(cW&x)    { meta::tauxZ<2,X,N-1,0>::v_ad(a,(cS*)x); rt; }
    tS V&op-=     (cW&x)    { M::v_su(a,(cS*)x); rt; }
    tS V&sub_twice(cW&x)    { meta::tauxZ<2,X,N-1,0>::v_su(a,(cS*)x); rt; }
    tS V op+      (cW&x) c_ { register V y; M::v_sum(y.a,a,(cS*)x); return y; }
    tS V op-      (cW&x) c_ { register V y; M::v_dif(y.a,a,(cS*)x); return y; }
    tS X op*      (cW&x) c_ { return M::v_dot(a,(cS*)x); }
    tS X dist_sq  (cW&x) c_ { return M::v_diq(a,(cS*)x); }
    tS X dist     (cW&x) c_ { return sqrt(dist_sq (x)); }
    tS X sum_sq   (cW&x) c_ { return M::v_suq(a,(cS*)x);}
    // miscellaneous                                                            
    tS V&ass_times(cV&x,cS&f) { M::v_ast(a,x.a,f); rt; }
    tS V&add_times(cV&x,cS&f) { M::v_adt(a,x.a,f); rt; }
    tS V&sub_times(cV&x,cS&f) { M::v_sut(a,x.a,f); rt; }
    V&up_max      (cV&x,cX&f) { M::v_umax(a,x.a,f); rt; }
    V&up_min      (cV&x,cX&f) { M::v_umix(a,x.a,f); rt; }
    // formatted I/O                                                            
    friend std::ostream&op<<(std::ostream&s,cV&x) { M::v_out(s,x.a); return s; }
    friend std::istream&op>>(std::istream&s,V &x) { M::v_in (s,x.a); return s; }
  };
#undef  W
#undef  cW
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // non-member functions, taking tupel<> arguments                           //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
#define tNX    template<int N, typename X>
#define tNXY   template<int N, typename X, typename Y>
  // norm, dist_sq, dist, and sum_sq of tupel(s)                                
  tNX  i_ X max    (tupel<N,X> c_&x)                  { return x.max(); }
  tNX  i_ X maxnorm(tupel<N,X> c_&x)                  { return x.maxnorm(); }
  tNX  i_ X min    (tupel<N,X> c_&x)                  { return x.min(); }
  tNX  i_ X norm   (tupel<N,X> c_&x)                  { return x.norm(); }
  tNX  i_ X abs    (tupel<N,X> c_&x)                  { return x.abs(); }
  tNXY i_ X dist_sq(tupel<N,X> c_&x, tupel<N,Y> c_&y) { return x.dist_sq(y); }
  tNXY i_ X dist   (tupel<N,X> c_&x, tupel<N,Y> c_&y) { return x.dist(y); }
  tNXY i_ X sum_sq (tupel<N,X> c_&x, tupel<N,Y> c_&y) { return x.sum_sq(y); }
  // product scalar * tupel                                                     
  tNXY i_ tupel<N,X> op*(Y c_&y, tupel<N,X> c_&v) { return v*y; }
  // cross product for 2D and 3D, coded as operator^                            
  tX   i_ X          op^(tupel<2,X> c_&x, tupel<2,X> c_&y) {
    return meta::taux<X,0>::v_cr2(static_cast<c_ X*>(x),
				  static_cast<c_ X*>(y));
  }
  tX   i_ tupel<3,X> op^(tupel<3,X> c_&x, tupel<3,X> c_&y) {
    register tupel<3,X> z;
    meta::taux<X,0>::v_cr3(static_cast<   X*>(z),
			   static_cast<c_ X*>(x), 
			   static_cast<c_ X*>(y));
    return z;
  }
#undef  tX
#undef  tNX
#undef  tNXY
  //////////////////////////////////////////////////////////////////////////////
#undef  cS
#undef  tS
#undef  op
#undef  c_
#undef  i_
#undef  rt
////////////////////////////////////////////////////////////////////////////////
}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
#endif	                                           // falcON_included_tupl_h    
