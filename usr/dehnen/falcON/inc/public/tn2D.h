// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// tn2D.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 1999-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// symmetric tensors of order K in 2 dimensions                                |
//                                                                             |
// A tensor of order K in D dimensions has K indices which each can have the D |
// values 0,...,D-1. Thus in 2D, these are 0 or 2. A tensor is symmetric       |
// if it is unchanged by any permutation of its indices. Thus, for a symmetric |
// tensor the ordering of the indices does not matter, only the number of      |
// indices with given value. This allows us to greatly reduce the indexing by  |
// introducing a set of "super indices":                                       |
//                                                                             |
//     k_i = number of indices i_0,...,i_K that are equal to i.                |
//                                                                             |
// Since                                                                       |
//                                                                             |
//     Sum_{i=0}^{D-1} k_i = K,                                                |
//                                                                             |
// there is some redundancy, if we are going to keep the order K as a useful   |
// number.                                                                     |
// In 2D, kx,ky may be replaced by an alternative set of "super indices":      |
//                                                                             |
//     K  = kx + ky  :  number of indices;                                     |
//     L  =      ky  :  number of indices with value >= 1;                     |
//                                                                             |
// Conversely                                                                  |
//                                                                             |
//     kx = K - L,                                                             |
//     ky = L.                                                                 |
//                                                                             |
// At given K, we have                                                         |
//                                                                             |
//     0 <= L <= K.                                                            |
//                                                                             |
// This implies that at given order K, there are                               |
//                                                                             |
//     K+1                                                                     |
//                                                                             |
// independent tensor elements (for a general asymmetric tensor it's 2^K).     |
//                                                                             |
// The number of possible permutations of the set of indices (i1,...,iK) is    |
// given by                                                                    |
//                                                                             |
//     P(K,L) = binomial(K,L)                                                  |
//            = K! / (kx! * ky!)                                               |
//                                                                             |
// There are numerous possible operations that can be done with symmetric      |
// tensors. Here, we support only those of immediate interest for falcON.      |
//                                                                             |
// 1. elementary operations and connections with scalar                        |
// 2. the K-fold outer product of a vector with itself                         |
// 3. the symmetrized outer product between tensors: T_K = T_K1 o T_K2, K=K1+K2|
// 4. the symmetrized outer product between a tensor and the (symmetrized)     |
//    N-fold outer product of delta_ij                                         |
// 5. the complete inner product: T_K = T_K1 * T_K2, K=K1-K2, K1 > K2          |
//    which sums over all indices of T_K2.                                     |
//                                                                             |
// In particular, we do NOT support intermediate products, such as for example |
// the ordinary matrix product, which may be considered an inner product over  |
// one index and an outer product over the second.                             |
// Note, however, that such an intermediate type of product may be equivalent  |
// to first a symmetrized outer product with the delta tensor of some order    |
// followed by an complete inner product.                                      |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_tn2D_h
#define falcON_included_tn2D_h

#ifndef falcON_included_tupl_h
#  include <public/tupl.h>
#endif

namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::symt2D<K,X>                                                  //
  //                                                                          //
  // a symmetric 2D tensor of order K with elements of type X                 //
  //                                                                          //
  // NOTES:                                                                   //
  //        for K=0, use a scalar instead                                     //
  //        for K=1, use a tupel<2,X> instead                                 //
  //        for K=2, we define a specialisation below                         //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<int K, typename X=real> class symt2D {
    symt2D(symt2D const&);                         // no copy constructor       
    // public static members and type                                           
  public:
    typedef X element_type;
    static const int NDAT = K+1;                   // # elements                
    static const int ORD  = K;                     // order of tensor           
    // private types                                                            
  private:
    typedef symt2D           T;                    // tensor                    
    typedef const T         cT;                    // const tensor              
    typedef meta::taux<X,K>  M;                    // meta looping l=0...K      
    // static methods (type conversions to X* and const X*)                     
    template<class A>
    static X*      pX(A      &x) { return static_cast<      X*>(x); }
    template<class A>
    static const X*pX(A const&x) { return static_cast<const X*>(x); }
    // datum                                                                    
    X a[NDAT];
  public:
    // construction                                                             
    symt2D()          {}
    symt2D(X const&x) { M::s_as(a,x); }
    // element access                                                           
    T      &operator[] (int i)       { return a[i]; }
    T const&operator[] (int i) const { return a[i]; }
//     T      &operator() (int L)       { return a[L]; }
//     T const&operator() (int L) const { return a[L]; }
    // unitary operators                                                        
    T&negate         ()               { M::v_neg(a); return*this; }
    T&set_zero       ()               { M::s_ze (a); return*this; }
    operator       X*()               { return a; }
    operator const X*() const         { return a; }
    // binary operators with scalar                                             
    T&operator=      (X const&s)      { M::s_as(a,s); return*this; }
    T&operator*=     (X const&s)      { M::s_ml(a,s); return*this; }
    T&operator/=     (X const&s)      { return operator*=(X(1)/s); }
    bool operator==  (X const&s) const{ return M::s_eq(a,s); }
    bool operator!=  (X const&s) const{ return M::s_neq(a,s); }
    // binary operators with symt2D<K,X>                                        
    T&operator=      (cT&t)           { M::v_as (a,t.a);   return*this; }
    T&operator+=     (cT&t)           { M::v_ad (a,t.a);   return*this; }
    T&operator-=     (cT&t)           { M::v_su (a,t.a);   return*this; }
    T&ass_times      (cT&t,X const&s) { M::v_ast(a,t.a,s); return*this; }
    T&add_times      (cT&t,X const&s) { M::v_adt(a,t.a,s); return*this; }
    T&sub_times      (cT&t,X const&s) { M::v_sut(a,t.a,s); return*this; }
    // unit matrix                                                              
    T&add_unity      ()               { a[0]+=X(1); a[K]+=X(1); return*this;}
    T&sub_unity      ()               { a[0]-=X(1); a[K]-=X(1); return*this;}
    T&ass_unity      ()               { set_zero(); return add_unity(X(1)); }
    T&add_unity      (X const&s)      { a[0]+=s; a[K]+=s; return*this; }
    T&sub_unity      (X const&s)      { a[0]-=s; a[K]-=s; return*this; }
    T&ass_unity      (X const&s)      { set_zero(); return add_unity(s); }
    // outer product of vector: T_{i}   := V_i1 * ... * V_ik                    
    // we need to know          T_{i-1}                                         
    T&ass_out_prd   (symt2D<K-1,X> const&, tupel<2,X> const&);
    T&add_out_prd   (symt2D<K-1,X> const&, tupel<2,X> const&);
    T&sub_out_prd   (symt2D<K-1,X> const&, tupel<2,X> const&);
    T&ass_out_prd   (symt2D<K-1,X> const&, tupel<2,X> const&, X const&);
    T&add_out_prd   (symt2D<K-1,X> const&, tupel<2,X> const&, X const&);
    T&sub_out_prd   (symt2D<K-1,X> const&, tupel<2,X> const&, X const&);
    // symmetric outer product                                                  
    //                                                                          
    // consider the outer product A between two symmetric tensors B and C of    
    // orders K_B and K_C, respectively. Clearly, the straightforward product   
    // will not be symmetric, but we can symmetrize it by summing over all      
    // Binom(K_A,K_B) ways to distribute the indices between B and C (note that 
    // K_A = K_B+K_C).                                                          
    // we =/+=/-= symmetric outer product of two lower-order tensors            
    template<int Q> T&ass_sym_prd(const symt2D<K-Q,X>&, const symt2D<Q,X>&);
    template<int Q> T&add_sym_prd(const symt2D<K-Q,X>&, const symt2D<Q,X>&);
    template<int Q> T&sub_sym_prd(const symt2D<K-Q,X>&, const symt2D<Q,X>&);
    // we =/+=/-= symmetric outer product of vector and tensor<K-1>             
    T&ass_sym_prd(symt2D<K-1,X> const&, tupel<2,X> const&);
    T&add_sym_prd(symt2D<K-1,X> const&, tupel<2,X> const&);
    T&sub_sym_prd(symt2D<K-1,X> const&, tupel<2,X> const&);
    T&ass_sym_prd(tupel<2,X> const&, symt2D<K-1,X> const&);
    T&add_sym_prd(tupel<2,X> const&, symt2D<K-1,X> const&);
    T&sub_sym_prd(tupel<2,X> const&, symt2D<K-1,X> const&);
    // we =/+=/-= s times the symmetric outer product of two lower-order tensors
    template<int Q> T&ass_sym_prd(const symt2D<K-Q,X>&,
				  const symt2D<Q  ,X>&, X const&);
    template<int Q> T&add_sym_prd(const symt2D<K-Q,X>&,
				  const symt2D<Q  ,X>&, X const&);
    template<int Q> T&sub_sym_prd(const symt2D<K-Q,X>&,
				  const symt2D<Q  ,X>&, X const&);
    // we =/+=/-= s times the symmetric outer product of vector and tensor<K-1> 
    T&ass_sym_prd(symt2D<K-1,X> const&, tupel<2,X> const&, X const&);
    T&add_sym_prd(symt2D<K-1,X> const&, tupel<2,X> const&, X const&);
    T&sub_sym_prd(symt2D<K-1,X> const&, tupel<2,X> const&, X const&);
    T&ass_sym_prd(tupel<2,X> const&, symt2D<K-1,X> const&, X const&);
    T&add_sym_prd(tupel<2,X> const&, symt2D<K-1,X> const&, X const&);
    T&sub_sym_prd(tupel<2,X> const&, symt2D<K-1,X> const&, X const&);
    // symmetric outer product of tensor with delta^(n)                         
    template<int Q> T&ass_del_prd(symt2D<Q,X> const&);
    template<int Q> T&ass_del_prd(symt2D<Q,X> const&, X const&);
    template<int Q> T&add_del_prd(symt2D<Q,X> const&);
    template<int Q> T&add_del_prd(symt2D<Q,X> const&, X const&);
    template<int Q> T&sub_del_prd(symt2D<Q,X> const&);
    template<int Q> T&sub_del_prd(symt2D<Q,X> const&, X const&);
    // symmetric outer product of vector with delta^(n)                         
    T&ass_del_prd(tupel<2,X> const&);
    T&ass_del_prd(tupel<2,X> const&, X const&);
    T&add_del_prd(tupel<2,X> const&);
    T&add_del_prd(tupel<2,X> const&, X const&);
    T&sub_del_prd(tupel<2,X> const&);
    T&sub_del_prd(tupel<2,X> const&, X const&);
    // symmetric outer product of delta^(n)                                     
    T&ass_del_prd();
    T&ass_del_prd(X const&);
    T&add_del_prd();
    T&add_del_prd(X const&);
    T&sub_del_prd();
    T&sub_del_prd(X const&);
    // complete inner product:                                                  
    // A{K-K2}[l-l2,m-m2] = A{K}[l,m] * C{K2}[l2,m2] * P{K2}[l2,m2]             
    void ass_inn_prd(tupel<2,X>   const&, symt2D<K-1,X>&) const;
    template<int Q>
    void ass_inn_prd(symt2D<Q,X>  const&, symt2D<K-Q,X>&) const;
    void ass_inn_prd(symt2D<K-1,X>const&, tupel<2,X>   &) const;
    void ass_inn_prd(symt2D<K,X>  const&, X            &) const;

    void add_inn_prd(tupel<2,X>   const&, symt2D<K-1,X>&) const;
    template<int Q>
    void add_inn_prd(symt2D<Q,X>  const&, symt2D<K-Q,X>&) const;
    void add_inn_prd(symt2D<K-1,X>const&, tupel<2,X>   &) const;
    void add_inn_prd(symt2D<K,X>  const&, X            &) const;

    void sub_inn_prd(tupel<2,X>   const&, symt2D<K-1,X>&) const;
    template<int Q>
    void sub_inn_prd(symt2D<Q,X>  const&, symt2D<K-Q,X>&) const;
    void sub_inn_prd(symt2D<K-1,X>const&, tupel<2,X>   &) const;
    void sub_inn_prd(symt2D<K,X>  const&, X            &) const;
    X    ass_inn_prd(cT&p) const;
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::symt2D<2,X>                                                  //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename X> class symt2D<2,X> {
    symt2D(symt2D const&);                         // no copy constructor       
    // public static members and type                                           
  public:
    typedef X element_type;
    static const int NDAT = 3;
    static const int K    = 2;
    // private types                                                            
  private:
    typedef symt3D           T;                    // tensor                    
    typedef const T         cT;                    // const tensor              
    typedef meta::taux<X,2>  M;                    // meta looping (l,m) to K   
    // static methods (type conversions to X* and const X*)                     
    template<class A>
    static X*      pX(A      &x) { return static_cast<      X*>(x); }
    template<class A>
    static const X*pX(A const&x) { return static_cast<const X*>(x); }
    // datum                                                                    
    X a[NDAT];
  public:
    // construction                                                             
    symt2D()          {}
    symt2D(X const&x) { M::s_as(a,x); }
    // element access                                                           
    T      &operator[] (int i)       { return a[i]; }
    T const&operator[] (int i) const { return a[i]; }
//     T      &operator() (int L)       { return a[L]; }
//     T const&operator() (int L) const { return a[L]; }
    // unitary operators                                                        
    T  &negate   ()                     { M::v_neg(a); return*this; }
    T  &set_zero ()                     { M::s_ze (a); return*this; }
    operator    X*     ()               { return a; }
    operator const X*  () const         { return a; }
    // binary operators with scalar                                             
    T  &operator=      (X const&s)      { M::s_as(a,s); return*this; }
    T  &operator*=     (X const&s)      { M::s_ml(a,s); return*this; }
    T  &operator/=     (X const&s)      { return operator*=(X(1)/s); }
    bool operator==    (X const&s) const{ return M::s_eq(a,s); }
    bool operator!=    (X const&s) const{ return M::s_neq(a,s); }
    // binary operators with symt2D<2,X>                                        
    T  &operator=      (cT&t)           { M::v_as (a,t.a);   return*this; }
    T  &operator+=     (cT&t)           { M::v_ad (a,t.a);   return*this; }
    T  &operator-=     (cT&t)           { M::v_su (a,t.a);   return*this; }
    T  &ass_times      (cT&t,X const&s) { M::v_ast(a,t.a,s); return*this; }
    T  &add_times      (cT&t,X const&s) { M::v_adt(a,t.a,s); return*this; }
    T  &sub_times      (cT&t,X const&s) { M::v_sut(a,t.a,s); return*this; }
    // unit matrix                                                              
    T  &add_unity      ()               { a[0]+=X(1); a[K]+=X(1); return*this; }
    T  &sub_unity      ()               { a[0]-=X(1); a[K]-=X(1); return*this; }
    T  &ass_unity      ()               { set_zero(); return add_unity(); }
    T  &add_unity      (X const&s)      { a[0]+=s; a[K]+=s; return*this; }
    T  &sub_unity      (X const&s)      { a[0]-=s; a[K]-=s; return*this; }
    T  &ass_unity      (X const&s)      { set_zero(); return add_unity(s); }
    // outer product of vector: T_{ij} := V_i * V_j                             
    T&ass_out_prd(tupel<2,X> const&);
    T&add_out_prd(tupel<2,X> const&);
    T&sub_out_prd(tupel<2,X> const&);
    T&ass_out_prd(tupel<2,X> const&, X const&);
    T&add_out_prd(tupel<2,X> const&, X const&);
    T&sub_out_prd(tupel<2,X> const&, X const&);
    // symmetric outer product                                                  
    // we =/+=/-= symmetric outer product of two vectors                        
    T&ass_sym_prd(tupel<2,X> const&, tupel<2,X> const&);
    T&add_sym_prd(tupel<2,X> const&, tupel<2,X> const&);
    T&sub_sym_prd(tupel<2,X> const&, tupel<2,X> const&);
    // we =/+=/-= s times the symmetric outer product of two vectors            
    T&ass_sym_prd(tupel<2,X> const&, tupel<2,X> const&, X const&);
    T&add_sym_prd(tupel<2,X> const&, tupel<2,X> const&, X const&);
    T&sub_sym_prd(tupel<2,X> const&, tupel<2,X> const&, X const&);
    // symmetric outer product of delta^(n)                                     
    T&ass_del_prd();
    T&ass_del_prd(X const&);
    T&add_del_prd();
    T&add_del_prd(X const&);
    T&sub_del_prd();
    T&sub_del_prd(X const&);
    // complete inner product                                                   
    void ass_inn_prd(tupel<2,X>  const&, tupel<2,X>&) const;
    void ass_inn_prd(symt2D<2,X> const&, X         &) const;
    void add_inn_prd(tupel<2,X>  const&, tupel<2,X>&) const;
    void add_inn_prd(symt2D<2,X> const&, X         &) const;
    void sub_inn_prd(tupel<2,X>  const&, tupel<2,X>&) const;
    void sub_inn_prd(symt2D<2,X> const&, X         &) const;
    X    ass_inn_prd(symt2D<2,X> const&) const;
    // trace                                                                    
    X trace() const { return a[0] + a[2]; }
    friend X trace (cT&t) { return t.trace(); }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::symt2D<1,X>                                                  //
  //                                                                          //
  // NOT to be used, use tupel<2,X> instead.                                  //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename X> class symt2D<1,X> {};
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // formatted I/O                                                            //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<int N, typename X> inline
  std::ostream&operator<< (std::ostream&s, symt2D<N,X> const&x) {
    meta::taux<X,symt2D<N,X>::KK>::v_out(s,static_cast<const X*>(x));
    return s;
  }
  template<int N, typename X> inline
  std::istream&operator>> (std::istream&s, symt2D<N,X>&x) {
    meta::taux<X,symt2D<N,X>::KK>::v_in(s,static_cast<X*>(x));
    return s;
  }
  //////////////////////////////////////////////////////////////////////////////
}                                                  // END: namespace nbdy       
#include <public/tn2D.cc>
////////////////////////////////////////////////////////////////////////////////
#endif	                                           // falcON_included_tn2D_h    
