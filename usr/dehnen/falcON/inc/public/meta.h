// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// meta.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2003-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_meta_h
#define falcON_included_meta_h

#ifndef falcON_included_inln_h
#  include <public/inln.h>
#endif
////////////////////////////////////////////////////////////////////////////////
#define c_      const
#define sv      static void
#define sb      static bool
#define sX      static X
#define cS      const S
#define tX      template<typename X>
#define tS      template<typename S>
#define tP      template<typename P>
////////////////////////////////////////////////////////////////////////////////
namespace meta {
  //////////////////////////////////////////////////////////////////////////////
  using nbdy::square;
  using nbdy::min;
  using nbdy::max;
  using nbdy::update_min;
  using nbdy::update_max;
  using nbdy::times;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class meta::taux<scalar,N,I>                                             //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename X, int N, int I=0> class taux {
    typedef const X       cX;
    typedef taux<X,N,I+1> M;
  public:
    tS sv s_as  ( X*a,cS&x)      { a[I] =x; M::s_as  (a,x); }
       sv s_ze  ( X*a)           { a[I] =X(0); M::s_ze(a); }
    tS sv s_ml  ( X*a,cS&x)      { a[I]*=x; M::s_ml(a,x); }
    tS sv v_as  ( X*a,cS*b)      { a[I] =b[I]; M::v_as(a,b); }
    tS sv v_ad  ( X*a,cS*b)      { a[I]+=b[I]; M::v_ad(a,b); }
    tS sv v_su  ( X*a,cS*b)      { a[I]-=b[I]; M::v_su(a,b); }
    tS sv v_ast ( X*a,cX*b,cS&x) { a[I] =x*b[I]; M::v_ast(a,b,x); }
    tS sv v_adt ( X*a,cX*b,cS&x) { a[I]+=x*b[I]; M::v_adt(a,b,x); }
    tS sv v_sut ( X*a,cX*b,cS&x) { a[I]-=x*b[I]; M::v_sut(a,b,x); }
       sv v_neg ( X*a)           { a[I]=-a[I]; M::v_neg(a); }
    tS sv v_nega( X*a,cS*b)      { a[I]=-b[I]; M::v_nega(a,b); }
    tS sv v_sum ( X*a,cX*b,cS*c) { a[I] =b[I]+c[I]; M::v_sum(a,b,c); }
    tS sv v_dif ( X*a,cX*b,cS*c) { a[I] =b[I]-c[I]; M::v_dif(a,b,c); }
       sb s_eq  (cX*a,cX&b)      { return a[I]==b && M::s_eq(a,b); }
       sb s_neq (cX*a,cX&b)      { return a[I]!=b || M::s_neq(a,b); }
       sb v_eq  (cX*a,cX*b)      { return a[I]==b[I] && M::v_eq(a,b); }
       sb v_neq (cX*a,cX*b)      { return a[I]!=b[I] || M::v_neq(a,b); }
       sX v_norm(cX*a)           { return a[I] *a[I] + M::v_norm(a); }
    tS sX v_dot (cX*a,cS*b)      { return a[I] *b[I] + M::v_dot(a,b); }
    tS sX v_diq (cX*a,cS*b)      { return square(a[I]-b[I])+M::v_diq(a,b); }
    tS sX v_suq (cX*a,cS*b)      { return square(a[I]+b[I])+M::v_suq(a,b); }
       sX v_min (cX*a)           { return min(a[I], M::v_min(a)); }
       sX v_max (cX*a)           { return max(a[I], M::v_max(a)); }
       sX v_amax(cX*a)           { return max(abs(a[I]), M::v_amax(a)); }
       sv v_uma ( X*a,cX*b)      { update_max(a[I],b[I]); M::v_uma(a,b); } 
       sv v_umi ( X*a,cX*b)      { update_min(a[I],b[I]); M::v_umi(a,b); } 
       sv v_umax( X*a,cX*b,cX&x) { update_max(a[I],b[I],x); M::v_umax(a,b,x);}
       sv v_umix( X*a,cX*b,cX&x) { update_min(a[I],b[I],x); M::v_umix(a,b,x);}
       sv v_out (std::ostream&o, cX*a) { o<<a[I]<<' '; M::v_out(o,a); }
       sv v_in  (std::istream&i,  X*a) { i>>a[I]; M::v_in(i,a); }
    tP sv s_app(P &p, cX&x)      { p[I] = x; M::s_app(p,x); }
    tP sv s_apa(P &p, cX&x)      { p[I]+= x; M::s_apa(p,x); }
    tP sv s_aps(P &p, cX&x)      { p[I]-= x; M::s_aps(p,x); }
    tP sv v_ass(X*a,P c_&p)      { a[I] =p[I]; M::v_ass(a,p); }
    tP sv v_asa(X*a,P c_&p)      { a[I]+=p[I]; M::v_asa(a,p); }
    tP sv v_asu(X*a,P c_&p)      { a[I]-=p[I]; M::v_asu(a,p); }
    tP sv v_app(P &p, cX*a)      { p[I] =a[I]; M::v_app(p,a); }
    tP sv v_apa(P &p, cX*a)      { p[I]+=a[I]; M::v_apa(p,a); }
    tP sv v_aps(P &p, cX*a)      { p[I]-=a[I]; M::v_aps(p,a); }
    tP sv v_asst(X*a,P c_&p,cX&x){ a[I] = x*p[I]; M::v_asst(a,p,x); }
    tP sv v_asat(X*a,P c_&p,cX&x){ a[I]+= x*p[I]; M::v_asat(a,p,x); }
    tP sv v_asut(X*a,P c_&p,cX&x){ a[I]-= x*p[I]; M::v_asut(a,p,x); }
    tP sv v_appt(P &p, cX*a,cX&x){ p[I] = x*a[I]; M::v_appt(p,a,x); }
    tP sv v_apat(P &p, cX*a,cX&x){ p[I]+= x*a[I]; M::v_apat(p,a,x); }
    tP sv v_apst(P &p, cX*a,cX&x){ p[I]-= x*a[I]; M::v_apst(p,a,x); }
  };
  //----------------------------------------------------------------------------
  template<typename X, int I> class taux<X,I,I> {
    typedef const X      cX;
  public:
    tS sv s_as  ( X*a,cS&x)      { a[I] =x; }
       sv s_ze  ( X*a)           { a[I] =X(0); }
    tS sv s_ml  ( X*a,cS&x)      { a[I]*=x; }
    tS sv v_as  ( X*a,cS*b)      { a[I] =b[I]; }
    tS sv v_ad  ( X*a,cS*b)      { a[I]+=b[I]; }
    tS sv v_su  ( X*a,cS*b)      { a[I]-=b[I]; }
    tS sv v_ast ( X*a,cX*b,cS&x) { a[I] =x*b[I]; }
    tS sv v_adt ( X*a,cX*b,cS&x) { a[I]+=x*b[I]; }
    tS sv v_sut ( X*a,cX*b,cS&x) { a[I]-=x*b[I]; }
       sv v_neg ( X*a)           { a[I]=-a[I]; }
    tS sv v_nega( X*a,cS*b)      { a[I]=-b[I]; }
    tS sv v_sum ( X*a,cX*b,cS*c) { a[I] =b[I]+c[I]; }
    tS sv v_dif ( X*a,cX*b,cS*c) { a[I] =b[I]-c[I]; }
       sb s_eq  (cX*a,cX&b)      { return a[I]==b; }
       sb s_neq (cX*a,cX&b)      { return a[I]!=b; }
       sb v_eq  (cX*a,cX*b)      { return a[I]==b[I]; }
       sb v_neq (cX*a,cX*b)      { return a[I]!=b[I]; }
       sX v_norm(cX*a)           { return a[I] *a[I]; }
    tS sX v_dot (cX*a,cS*b)      { return a[I] *b[I]; }
    tS sX v_diq (cX*a,cS*b)      { return square(a[I]-b[I]); }
    tS sX v_suq (cX*a,cS*b)      { return square(a[I]+b[I]); }
       sX v_min (cX*a)           { return a[I]; }
       sX v_max (cX*a)           { return a[I]; }
       sX v_amax(cX*a)           { return abs(a[I]); }
       sv v_uma ( X*a,cX*b)      { update_max(a[I],b[I]); } 
       sv v_umi ( X*a,cX*b)      { update_min(a[I],b[I]); } 
       sv v_umax( X*a,cX*b,cX&x) { update_max(a[I],b[I],x); }
       sv v_umix( X*a,cX*b,cX&x) { update_min(a[I],b[I],x); }
       sv v_out (std::ostream&o, cX*a) { o<<a[I]; }
       sv v_in  (std::istream&i,  X*a) { i>>a[I]; }
    tP sv s_app (P &p, cX&x)     { p[I] = x; }
    tP sv s_apa (P &p, cX&x)     { p[I]+= x; }
    tP sv s_aps (P &p, cX&x)     { p[I]-= x; }
    tP sv v_ass (X*a,P c_&p)     { a[I] = p[I]; }
    tP sv v_asa (X*a,P c_&p)     { a[I]+= p[I]; }
    tP sv v_asu (X*a,P c_&p)     { a[I]-= p[I]; }
    tP sv v_app (P &p, cX*a)     { p[I] = a[I]; }
    tP sv v_apa (P &p, cX*a)     { p[I]+= a[I]; }
    tP sv v_aps (P &p, cX*a)     { p[I]-= a[I]; }
    tP sv v_asst(X*a,P c_&p,cX&x){ a[I] = x*p[I]; }
    tP sv v_asat(X*a,P c_&p,cX&x){ a[I]+= x*p[I]; }
    tP sv v_asut(X*a,P c_&p,cX&x){ a[I]-= x*p[I]; }
    tP sv v_appt(P &p, cX*a,cX&x){ p[I] = x*a[I]; }
    tP sv v_apat(P &p, cX*a,cX&x){ p[I]+= x*a[I]; }
    tP sv v_apst(P &p, cX*a,cX&x){ p[I]-= x*a[I]; }
  };
  //----------------------------------------------------------------------------
  template<typename X> class taux<X,0,0> {
    typedef const X      cX;
  public:
    sX v_cr2  (cX*a, cX*b)       { return a[0]*b[1] - a[1]*b[0]; }
    sv v_cr3  ( X*z, cX*a, cX*b) {
      z[0] = a[1]*b[2] - a[2]*b[1];
      z[1] = a[2]*b[0] - a[0]*b[2];
      z[2] = a[0]*b[1] - a[1]*b[0];
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class meta::tauxZ<Z,scalar,N,I>                                          //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<int Z, typename X, int N, int I=0> class tauxZ {
    typedef const X         cX;
    typedef tauxZ<Z,X,N,I+1> M;
  public:
    tS sv v_as  ( X*a,cS*b)      { a[I] =times<Z>(b[I]); M::v_as(a,b); }
    tS sv v_ad  ( X*a,cS*b)      { a[I]+=times<Z>(b[I]); M::v_ad(a,b); }
    tS sv v_su  ( X*a,cS*b)      { a[I]-=times<Z>(b[I]); M::v_su(a,b); }
    tS sv v_ast ( X*a,cX*b,cS&x) { a[I] =times<Z>(x*b[I]); M::v_ast(a,b,x);}
    tS sv v_adt ( X*a,cX*b,cS&x) { a[I]+=times<Z>(x*b[I]); M::v_adt(a,b,x);}
    tS sv v_sut ( X*a,cX*b,cS&x) { a[I]-=times<Z>(x*b[I]); M::v_sut(a,b,x);}
  };
  //----------------------------------------------------------------------------
  template<int Z, typename X, int I> class tauxZ<Z,X,I,I> {
    typedef const X      cX;
  public:
    tS sv v_as  ( X*a,cS*b)      { a[I] =times<Z>(b[I]); }
    tS sv v_ad  ( X*a,cS*b)      { a[I]+=times<Z>(b[I]); }
    tS sv v_su  ( X*a,cS*b)      { a[I]-=times<Z>(b[I]); }
    tS sv v_ast ( X*a,cX*b,cS&x) { a[I] =times<Z>(x*b[I]); }
    tS sv v_adt ( X*a,cX*b,cS&x) { a[I]+=times<Z>(x*b[I]); }
    tS sv v_sut ( X*a,cX*b,cS&x) { a[I]-=times<Z>(x*b[I]); }
  };
  //////////////////////////////////////////////////////////////////////////////
  const double iN[16] = { 1.,                      // 1/N                       
			  1.,
			  0.5,
			  0.3333333333333333333333333,
			  0.25,
			  0.2,
			  0.1666666666666666666666667,
			  0.1428571428571428571428571,
			  0.125,
			  0.1111111111111111111111111,
			  0.1,
			  0.09090909090909090909090909,
			  0.08333333333333333333333333,
			  0.076923076923076923076923076,
			  0.071428571428571428571428571,
			  0.066666666666666666666666667 };
  //////////////////////////////////////////////////////////////////////////////
  const double iF[16] = { 1.,                      // 1/N!                      
			  1.,
			  0.5,
			  0.166666666666666666666667,
			  0.416666666666666666666667e-1,
			  0.833333333333333333333333e-2,
			  0.138888888888888888888889e-2,
			  0.198412698412698412698413e-3,
			  0.248015873015873015873016e-4,
			  0.275573192239858906525573e-5,
			  0.275573192239858906525573e-6,
			  0.250521083854417187750521e-7,
			  0.208767569878680989792101e-8,
			  0.160590438368216145993924e-9,
			  0.114707455977297247138517e-10,
			  0.764716373181981647590113e-12 };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // struct ONE<N>                                                            //
  //                                                                          //
  // F: factorial N!!                                                         //
  // G: N!!                                                                   //
  // H: (2*N-1)!!                                                             //
  // P: 2^N                                                                   //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<int N> struct ONE {
    enum {
      F  = N*ONE<N-1>::F,
      G  = N*ONE<N-2>::G,
      H  = (N+N-1)*ONE<N-1>::H,
      P  = 4*ONE<N-2>::P
    }; };
  template<> struct ONE<1> { enum { F=1, G=1, H=1, P=2 }; };
  template<> struct ONE<0> { enum { F=1, G=1, H=1, P=1 }; };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // struct TWO<L,M>                                                          //
  //                                                                          //
  // B: binomial (L,M)                                                        //
  // I: index of (L,M)                                                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<int L, int M> struct TWO {
    enum { B = TWO<L-1,M-1>::B + TWO<L-1,M>::B,
	   I = (L*(L+1))/2+M
    }; };
  template<int M> struct TWO<0,M> { enum { B=1, I=0 }; };
  template<int L> struct TWO<L,0> { enum { B=1, I=(L*(L+1))/2 }; };
  template<int L> struct TWO<L,L> { enum { B=1, I=(L*(L+3))/2 }; };
  template<>      struct TWO<0,0> { enum { B=1, I=0 }; };
  //////////////////////////////////////////////////////////////////////////////
}                                                  // END: namespace meta       
////////////////////////////////////////////////////////////////////////////////
#undef c_
#undef sv
#undef sb
#undef sX
#undef cS
#undef tX
#undef tS
#undef tP
////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_included_meta_h    
