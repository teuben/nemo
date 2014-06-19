// -*- C++ -*-                                                                 
////////////////////////////////////////////////////////////////////////////////
//
// tensor.h
//
// Copyright (C) 1999-2006,2009,2012 Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 675
// Mass Ave, Cambridge, MA 02139, USA.
//
////////////////////////////////////////////////////////////////////////////////
//
// symmetric tensors of order K in 3 dimensions
//
// A tensor of order K in D dimensions has K indices which each can have the D
// values 0,...,D-1. Thus in 3D, these are 0,1, or 2. A tensor is symmetric if
// it is unchanged by any permutation of its indices. Thus, for a symmetric
// tensor the ordering of the indices does not matter, only the number of
// indices with given value. This allows us to greatly reduce the indexing by
// introducing a set of "super indices":
//
//     k_j := number of indices i_0,...,i_K that are equal to j.
//
// Since
//
//     Sum_{j=0}^{D-1} k_j = K,
//
// there is some redundancy, if we are going to keep the order K as a useful
// number.
// In 3D, kx,ky,kz may be replaced by an alternative set of super indices:
//
//     K  = kx + ky + kz  :  number of indices;
//     L  =      ky + kz  :  number of indices with value >= 1;
//     M  =           kz  :  number of indices with value >= 2.
//
// Conversely
//
//     kx = K - L,
//     ky = L - M,
//     kz = M.
//
// At given K, we have
//
//     0 <= L <= K;
//     0 <= M <= L.
//
// This implies that at given order K, there are
//
//     (K+1)*(K+2)/2
//
// independent tensor elements (for a general asymmetric tensor it's 3^K).
//
// The number of possible permutations of the set of indices (i1,...,iK) is
// given by
//
//     P(K,L,M) = binomial(K,L) * binomial(L,M)
//              = K! / (kx! * ky! * kz!)
//
// There are numerous possible operations that can be done with symmetric
// tensors. Here, we support only those of immediate interest for falcON.
//
// 1. elementary operations and connections with scalar
// 2. the K-fold outer product of a vector with itself
// 3. the symmetrized outer product between tensors: T_K = T_K1 o T_K2,
//    K=K1+K2
// 4. the symmetrized outer product between a tensor and the (symmetrized)
//    N-fold outer product of delta_ij
// 5. the complete inner product: T_K = T_K1 * T_K2, K=K1-K2, K1 > K2 which
//    sums over all indices of T_K2.
//
// In particular, we do NOT support intermediate products, such as for example
// the ordinary matrix product, which may be considered an inner product over
// one index and an outer product over the second.  Note, however, that such
// an intermediate type of product may be equivalent to first a symmetrized
// outer product with the delta tensor of some order followed by an complete
// inner product.
//
////////////////////////////////////////////////////////////////////////////////
//
// RULES
//
// 0. Notations
//                                                       (n)
// 0.1 n-fold outer product of vector with itself:      X
//
// 0.2 symmetric outer product between tensors:         X o Y
//
//
// 1. Outer Product of a Sum of Tensors
//
// The outer self-product of a sum of tensors is equal to
//
//         (n)      n   (k)    (n-k)
//   (A + B)   = Sum   A    o B
//                 k=0
//
// The usual binomial coefficient is absorbed into the symmetric product.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef falcON_included_tensor_h
#define falcON_included_tensor_h

#ifndef falcON_included_utils_h
#  include <public/utils.h>
#endif
#ifndef falcON_included_tupel_h
#  include <utils/tupel.h>                       // for WDutils::meta::taux<>
#endif
// #ifndef falcON_included_tupel_cc
// #  include <utils/tupel.cc>                      // for WDutils::meta::taux<>
// #endif


namespace falcON {
  //
  // class falcON::symt3D<K,X>
  //
  // a symmetric 3D tensor of order K with elements of type X
  //
  // NOTES:
  //        for K=0, use a scalar instead
  //        for K=1, use a falcONVec<3,X> instead
  //        for K=2, we define a specialisation below
  //
  template<int K, typename X=real> class symt3D {
    symt3D(symt3D const&);                         // no copy constructor
    // public static members and type
  public:
    typedef X element_type;
    static const int NDAT = ((K+1)*(K+2))/2;       // # elements
    static const int K0   = (K*(K+1))/2;           // index of element [1,1]
    static const int KK   = NDAT-1;                // index of element [2,2]
    static const int ORD  = K;                     // order of tensor
    // private types
  private:
    typedef symt3D           T;                    // tensor
    typedef const T         cT;                    // const tensor
    typedef meta::taux<X,KK> M;                    // meta looping (l,m) to K
    // static methods (type conversions to X* and const X*)
    template<class A>
    static X*      pX(A      &x) { return static_cast<      X*>(x); }
    template<class A>
    static const X*pX(A const&x) { return static_cast<const X*>(x); }
    // datum
    X a[NDAT];
  public:
    // construction
    symt3D()                   {}
    explicit symt3D(X const&x) { M::s_as(a,x); }
    // element access
    X      &operator[] (int i)              { return a[i]; }
    X const&operator[] (int i)        const { return a[i]; }
    X      &operator() (int l, int m)       { return a[(l*(l+1))/2+m]; }
    X const&operator() (int l, int m) const { return a[(l*(l+1))/2+m]; }
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
    // binary operators with symt3D<K,X>
    T&operator=      (cT&t)           { M::v_as (a,t.a);   return*this; }
    T&operator+=     (cT&t)           { M::v_ad (a,t.a);   return*this; }
    T&operator-=     (cT&t)           { M::v_su (a,t.a);   return*this; }
    T&ass_times      (cT&t,X const&s) { M::v_ast(a,t.a,s); return*this; }
    T&add_times      (cT&t,X const&s) { M::v_adt(a,t.a,s); return*this; }
    T&sub_times      (cT&t,X const&s) { M::v_sut(a,t.a,s); return*this; }
    // unit matrix
    T&add_unity      () { a[0]+=X(1);a[K0]+=X(1);a[KK]+=X(1);return*this;}
    T&sub_unity      () { a[0]-=X(1);a[K0]-=X(1);a[KK]-=X(1);return*this;}
    T&ass_unity      () { set_zero(); return add_unity(X(1)); }
    T&add_unity      (X const&s) { a[0]+=s; a[K0]+=s; a[KK]+=s; return*this; }
    T&sub_unity      (X const&s) { a[0]-=s; a[K0]-=s; a[KK]-=s; return*this; }
    T&ass_unity      (X const&s) { set_zero(); return add_unity(s); }
    // outer product of vector: T_{i}   := V_i1 * ... * V_ik
    // we need to know          T_{i-1}
    T&ass_out_prd   (symt3D<K-1,X> const&, falcONVec<3,X> const&);
    T&add_out_prd   (symt3D<K-1,X> const&, falcONVec<3,X> const&);
    T&sub_out_prd   (symt3D<K-1,X> const&, falcONVec<3,X> const&);
    T&ass_out_prd   (symt3D<K-1,X> const&, falcONVec<3,X> const&, X const&);
    T&add_out_prd   (symt3D<K-1,X> const&, falcONVec<3,X> const&, X const&);
    T&sub_out_prd   (symt3D<K-1,X> const&, falcONVec<3,X> const&, X const&);
    // symmetric outer product
    //
    // consider the outer product A between two symmetric tensors B and C of
    // orders K_B and K_C, respectively. Clearly, the straightforward product
    // will not be symmetric, but we can symmetrize it by summing over all
    // Binom(K_A,K_B) ways to distribute the indices between B and C (note that
    // K_A = K_B+K_C).
    //
    // we =/+=/-= symmetric outer product of two lower-order tensors
    template<int Q> T&ass_sym_prd(const symt3D<K-Q,X>&, const symt3D<Q,X>&);
    template<int Q> T&add_sym_prd(const symt3D<K-Q,X>&, const symt3D<Q,X>&);
    template<int Q> T&sub_sym_prd(const symt3D<K-Q,X>&, const symt3D<Q,X>&);
    // we =/+=/-= symmetric outer product of vector and tensor<K-1>
    T&ass_sym_prd(symt3D<K-1,X> const&, falcONVec<3,X> const&);
    T&add_sym_prd(symt3D<K-1,X> const&, falcONVec<3,X> const&);
    T&sub_sym_prd(symt3D<K-1,X> const&, falcONVec<3,X> const&);
    T&ass_sym_prd(falcONVec<3,X> const&, symt3D<K-1,X> const&);
    T&add_sym_prd(falcONVec<3,X> const&, symt3D<K-1,X> const&);
    T&sub_sym_prd(falcONVec<3,X> const&, symt3D<K-1,X> const&);
    // we =/+=/-= s times the symmetric outer product of two lower-order tensors
    template<int Q> T&ass_sym_prd(const symt3D<K-Q,X>&,
				  const symt3D<Q  ,X>&, X const&);
    template<int Q> T&add_sym_prd(const symt3D<K-Q,X>&,
				  const symt3D<Q  ,X>&, X const&);
    template<int Q> T&sub_sym_prd(const symt3D<K-Q,X>&,
				  const symt3D<Q  ,X>&, X const&);
    // we =/+=/-= s times the symmetric outer product of vector and tensor<K-1> 
    T&ass_sym_prd(symt3D<K-1,X> const&, falcONVec<3,X> const&, X const&);
    T&add_sym_prd(symt3D<K-1,X> const&, falcONVec<3,X> const&, X const&);
    T&sub_sym_prd(symt3D<K-1,X> const&, falcONVec<3,X> const&, X const&);
    T&ass_sym_prd(falcONVec<3,X> const&, symt3D<K-1,X> const&, X const&);
    T&add_sym_prd(falcONVec<3,X> const&, symt3D<K-1,X> const&, X const&);
    T&sub_sym_prd(falcONVec<3,X> const&, symt3D<K-1,X> const&, X const&);
    // symmetric outer product of tensor with delta^(n)
    template<int Q> T&ass_del_prd(symt3D<Q,X> const&);
    template<int Q> T&ass_del_prd(symt3D<Q,X> const&, X const&);
    template<int Q> T&add_del_prd(symt3D<Q,X> const&);
    template<int Q> T&add_del_prd(symt3D<Q,X> const&, X const&);
    template<int Q> T&sub_del_prd(symt3D<Q,X> const&);
    template<int Q> T&sub_del_prd(symt3D<Q,X> const&, X const&);
    // symmetric outer product of vector with delta^(n)
    T&ass_del_prd(falcONVec<3,X> const&);
    T&ass_del_prd(falcONVec<3,X> const&, X const&);
    T&add_del_prd(falcONVec<3,X> const&);
    T&add_del_prd(falcONVec<3,X> const&, X const&);
    T&sub_del_prd(falcONVec<3,X> const&);
    T&sub_del_prd(falcONVec<3,X> const&, X const&);
    // symmetric outer product of delta^(n)
    T&ass_del_prd();
    T&ass_del_prd(X const&);
    T&add_del_prd();
    T&add_del_prd(X const&);
    T&sub_del_prd();
    T&sub_del_prd(X const&);
    // complete inner product:
    // A{K-K2}[l-l2,m-m2] = A{K}[l,m] * C{K2}[l2,m2] * P{K2}[l2,m2]
    void ass_inn_prd(falcONVec<3,X> const&, symt3D<K-1,X>&) const;
    template<int Q>
    void ass_inn_prd(symt3D<Q,X> const&, symt3D<K-Q,X>&) const;
    void ass_inn_prd(symt3D<K-1,X> const&, falcONVec<3,X>&) const;
    void ass_inn_prd(symt3D<K,X> const&, X&) const;

    void add_inn_prd(falcONVec<3,X> const&, symt3D<K-1,X>&) const;
    template<int Q>
    void add_inn_prd(symt3D<Q,X> const&, symt3D<K-Q,X>&) const;
    void add_inn_prd(symt3D<K-1,X> const&, falcONVec<3,X>&) const;
    void add_inn_prd(symt3D<K,X> const&, X&) const;

    void sub_inn_prd(falcONVec<3,X> const&, symt3D<K-1,X>&) const;
    template<int Q>
    void sub_inn_prd(symt3D<Q,X> const&, symt3D<K-Q,X>&) const;
    void sub_inn_prd(symt3D<K-1,X> const&, falcONVec<3,X>&) const;
    void sub_inn_prd(symt3D<K,X> const&, X&) const;
    X    ass_inn_prd(cT&p) const;
  };
  //
  // class falcON::symt3D<2,X>
  //
  template<typename X> class symt3D<2,X> {
#if __cplusplus >= 201103L
    symt3D(symt3D const&) = delete;                // no copy constructor
#else
    symt3D(symt3D const&);                         // no copy constructor
#endif
    // public static members and type
  public:
    typedef X element_type;
    static const int NDAT = 6;
    static const int K0   = 3;
    static const int KK   = 5;
    // private types
  private:
    typedef symt3D           T;                    // tensor
    typedef const T         cT;                    // const tensor
    typedef meta::taux<X,5>  M;                    // meta looping (l,m) to K
    // static methods (type conversions to X* and const X*)
    template<class A>
    static X*      pX(A      &x) { return static_cast<      X*>(x); }
    template<class A>
    static const X*pX(A const&x) { return static_cast<const X*>(x); }
    // datum
    X a[NDAT];
  public:
    // construction
    symt3D() {}
    explicit symt3D(X const&x) { M::s_as(a,x); }
    // element access
    X      &operator[] (int i)              { return a[i]; }
    X const&operator[] (int i)        const { return a[i]; }
    X      &operator() (int l, int m)       { return a[(l*(l+1))/2+m]; }
    X const&operator() (int l, int m) const { return a[(l*(l+1))/2+m]; }
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
    // binary operators with symt3D<2,X>
    T  &operator=      (cT&t)           { M::v_as (a,t.a);   return*this; }
    T  &operator+=     (cT&t)           { M::v_ad (a,t.a);   return*this; }
    T  &operator-=     (cT&t)           { M::v_su (a,t.a);   return*this; }
    T  &ass_times      (cT&t,X const&s) { M::v_ast(a,t.a,s); return*this; }
    T  &add_times      (cT&t,X const&s) { M::v_adt(a,t.a,s); return*this; }
    T  &sub_times      (cT&t,X const&s) { M::v_sut(a,t.a,s); return*this; }
    // unit matrix
    T  &add_unity      () {
      a[0]+=X(1); a[K0]+=X(1); a[KK]+=X(1); return*this; }
    T  &sub_unity      () {
      a[0]-=X(1); a[K0]-=X(1); a[KK]-=X(1); return*this; }
    T  &ass_unity      () {
      set_zero(); return add_unity(); }
    T  &add_unity      (X const&s) {
      a[0]+=s; a[K0]+=s; a[KK]+=s; return*this; }
    T  &sub_unity      (X const&s) {
      a[0]-=s; a[K0]-=s; a[KK]-=s; return*this; }
    T  &ass_unity      (X const&s) {
      set_zero(); return add_unity(s); }
    // outer product of vector: T_{ij} := V_i * V_j
    T&ass_out_prd(falcONVec<3,X> const&);
    T&add_out_prd(falcONVec<3,X> const&);
    T&sub_out_prd(falcONVec<3,X> const&);
    T&ass_out_prd(falcONVec<3,X> const&, X const&);
    T&add_out_prd(falcONVec<3,X> const&, X const&);
    T&sub_out_prd(falcONVec<3,X> const&, X const&);
    // symmetric outer product
    // we =/+=/-= symmetric outer product of two vectors
    T&ass_sym_prd(falcONVec<3,X> const&, falcONVec<3,X> const&);
    T&add_sym_prd(falcONVec<3,X> const&, falcONVec<3,X> const&);
    T&sub_sym_prd(falcONVec<3,X> const&, falcONVec<3,X> const&);
    // we =/+=/-= s times the symmetric outer product of two vectors
    T&ass_sym_prd(falcONVec<3,X> const&, falcONVec<3,X> const&, X const&);
    T&add_sym_prd(falcONVec<3,X> const&, falcONVec<3,X> const&, X const&);
    T&sub_sym_prd(falcONVec<3,X> const&, falcONVec<3,X> const&, X const&);
    // symmetric outer product of delta^(n)
    T&ass_del_prd();
    T&ass_del_prd(X const&);
    T&add_del_prd();
    T&add_del_prd(X const&);
    T&sub_del_prd();
    T&sub_del_prd(X const&);
    // complete inner product
    void ass_inn_prd(falcONVec<3,X>  const&, falcONVec<3,X>&) const;
    void ass_inn_prd(symt3D<2,X> const&, X         &) const;
    void add_inn_prd(falcONVec<3,X>  const&, falcONVec<3,X>&) const;
    void add_inn_prd(symt3D<2,X> const&, X         &) const;
    void sub_inn_prd(falcONVec<3,X>  const&, falcONVec<3,X>&) const;
    void sub_inn_prd(symt3D<2,X> const&, X         &) const;
    X    ass_inn_prd(symt3D<2,X> const&) const;
    // trace
    X trace() const { return a[0] + a[3] + a[5]; }
    friend X trace (cT&t) { return t.trace(); }
  };
  //
  // class falcON::symt3D<1,X>
  //
  // NOT to be used, use falcONVec<3,X> instead.
  //
  template<typename X> class symt3D<1,X> {};
  //
  // formatted I/O
  //
  template<int N, typename X> inline
  std::ostream&operator<< (std::ostream&s, symt3D<N,X> const&x) {
    meta::taux<X,symt3D<N,X>::KK>::v_out(s,static_cast<const X*>(x));
    return s;
  }
  template<int N, typename X> inline
  std::istream&operator>> (std::istream&s, symt3D<N,X>&x) {
    meta::taux<X,symt3D<N,X>::KK>::v_in(s,static_cast<X*>(x));
    return s;
  }
} // namespace falcON
namespace WDutils {
  //
  template<int N, typename T> struct traits< falcON::symt3D<N,T> > {
    static const char  *name() {
      return message("symt3D<%d,%s>", N, traits<T>::name());
    }
  };
} // namespace WDutils
#include <public/tensor.cc>
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_tensor_h
