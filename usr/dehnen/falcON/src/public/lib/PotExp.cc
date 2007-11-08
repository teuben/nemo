// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
/// \file src/public/lib/PotExp.cc                                             |
//                                                                             |
// Copyright (C) 1994-1996, 2004-2007 Walter Dehnen                            |
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
#include <public/PotExp.h>
#include <public/simd.h>
#include <cstring>
#include <iomanip>
#include <rpc/rpc.h>
#ifdef   falcON_NEMO
#  ifdef   MAX               /* not only Peter Teuben, but also somebody */
#    undef MAX               /* else had the splendid idea to define a   */
#  endif                     /* C-macro "MAX".                           */
#  ifdef   MIN               /* not only Peter Teuben, but also somebody */
#    undef MIN               /* else had the splendid idea to define a   */
#  endif                     /* C-macro "MIN".                           */
#  include <stdinc.h>
#  include <cstring>
#endif

namespace falcON {
  bool is_appended(const char*name, char c, char*copy);
}

using namespace falcON;
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// auxiliary methods for implementing class falcON::PotExp                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // types                                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  using falcON::tupel;
  using falcON::fvec4;
  typedef PotExp::scalar   scalar;
  typedef PotExp::symmetry symmetry;
  typedef PotExp::Anlm     Anlm;
  //----------------------------------------------------------------------------
  //                                                                            
  // class CasRec                                                               
  //                                                                            
  // holds the numbers a_m with m=-L...L                                        
  //                                                                            
  // this type is actually not really needed; we keep it for reference.         
  //                                                                            
  // NOTE                                                                       
  //    Depending on the adopted symmetry, we ASSUME certain elements to be     
  //    zero or trivial copies, but do NOT ENFORCE this. In particular:         
  //    1. spherical or cylindrical symmetry: only the m=0 element is used.     
  //    2. triaxial symmetry: only elements with even m>0 are used and it is    
  //       implicitly assumed that A[-m]=A[m], i.e. only  cos(m phi)  terms.    
  //    3. reflexion symmetry: only elements with even m are used.              
  //                                                                            
  class CasRec {
    friend class PotExpAccess;
    int     L;                                     // m in [-L,L]               
    scalar *D, *A;                                 // pointer to memory         
  public:
    // construction & destruction                                               
    CasRec(int l) : L(l), D( falcON_NEW(scalar,l+l+1) ), A(l+D) {}
    ~CasRec() { falcON_DEL_A(D); }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // element access                                                           
    //    NOTE: that elements ASSUMED to be zero may have any value!            
    scalar      &operator() (int m)       { return A[m]; }
    scalar const&operator() (int m) const { return A[m]; }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    int const& mmax() const { return L; }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // set to cas(m*phi) := cos(m*phi) + sin(m*phi)                             
    template<symmetry>
    CasRec&set(scalar const&, scalar const&);      // I: cos(phi), sin(phi)     
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    template<symmetry> CasRec&reset   ();              // A[m] = 0              
    template<symmetry> CasRec&assign  (scalar const&); // A[m] = x              
    template<symmetry> CasRec&multiply(scalar const&); // A[m]*= x              
    template<symmetry> CasRec&divide  (scalar const&); // A[m]/= x              
    template<symmetry> CasRec&copy    (CasRec const&); // A[m] = B[m]           
    template<symmetry> CasRec&add     (CasRec const&); // A[m]+= B[m]           
    template<symmetry> CasRec&sub     (CasRec const&); // A[m]-= B[m]           
    template<symmetry> scalar dot     (CasRec const&) const; // dot product     
    template<symmetry> CasRec&addtimes(CasRec const&,
				       scalar const&); // A[m]+= x*B[m]         
    template<symmetry> CasRec&subtimes(CasRec const&,
				       scalar const&); // A[m]-= x*B[m]         
  }; // class CasRec
  //----------------------------------------------------------------------------
  //                                                                            
  // class YlmRec                                                               
  //                                                                            
  // holds the numbers a_lm with l=0...L, m=-l...l                              
  //                                                                            
  // NOTE                                                                       
  //    Depending on the adopted symmetry, we ASSUME certain elements to be     
  //    zero or trivial copies, but do NOT ENFORCE this. In particular:         
  //    1. spherical symmetry: only the l=0,m=0 element is used.                
  //    2. cylindrical symmetry: only the even l, m=0 elements are used.        
  //    3. triaxial symmetry: only elements with (l,m) even are considered      
  //       and it is implicitly assumed that A[-m]=A[m]; only m>0 are used.     
  //    4. reflexion symmetry: all elements with even (l,m) are used.           
  class YlmRec {
    friend class PotExpAccess;
    int     L,L1,L1Q;                              // L, L+1, (L+1)^2           
    scalar *A;
    scalar      &e (int l, int m)       { return A[l*(l+1)+m]; }
    scalar const&e (int l, int m) const { return A[l*(l+1)+m]; }
  public:
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // construction & destruction                                               
    YlmRec(int l)
      : L(l),   L1(L+1), L1Q(L1*L1), A(falcON_NEW(scalar,L1Q)) {}
    YlmRec(YlmRec const&a)
      : L(a.L), L1(L+1), L1Q(L1*L1), A(falcON_NEW(scalar,L1Q)) {
      memcpy(A,a.A,L1Q*sizeof(scalar));
    }
    ~YlmRec() { falcON_DEL_A(A); }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // element access                                                           
    //    NOTE: that elements ASSUMED to be zero may have any value!            
    scalar      &operator() (int l, int m)       { return e(l,m); }
    scalar const&operator() (int l, int m) const { return e(l,m); }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    int const& lmax() const { return L; }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void table_print(symmetry, std::ostream&, int=6) const;
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  private: // only accessible to PotExp in file pexp.cc                         
    // set Ylm to non-normalized scalar-valued spherical harmonics, i.e.        
    //    Y(l,m;theta,phi) = P(l,|m|; cos[theta]) * cas(m*phi)                  
    template<symmetry>
    YlmRec&set(scalar const&ct, scalar const&st,   // I: cos(the),sin(the)      
	       scalar const&cp, scalar const&sp);  // I: cos(phi),sin(phi)      
    template<symmetry>
    YlmRec&set(YlmRec&, YlmRec&,                   // O: dY/dthe, dY/dphi       
	       scalar const&ct, scalar const&st,   // I: cos(the),sin(the)      
	       scalar const&cp, scalar const&sp);  // I: cos(phi),sin(phi)      
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // set Ylm to N_lm multiplied by 4Pi := (2l+1) (l-|m|)! / (l+|m|)!          
    YlmRec& Nlm();
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    template<symmetry> YlmRec&reset   ();              // A[l,m] = 0            
    template<symmetry> YlmRec&assign  (scalar const&); // A[l,m] = x            
    template<symmetry> YlmRec&multiply(scalar const&); // A[l,m]*= x            
    template<symmetry> YlmRec&divide  (scalar const&); // A[l,m]/= x            
    template<symmetry> YlmRec&copy    (YlmRec const&); // A[l,m] = B[l,m]       
    template<symmetry> YlmRec&add     (YlmRec const&); // A[l,m]+= B[l,m]       
    template<symmetry> YlmRec&sub     (YlmRec const&); // A[l,m]-= B[l,m]       
    template<symmetry> scalar dot     (YlmRec const&) const; // dot product     
    template<symmetry> YlmRec&addtimes(YlmRec const&,
				       scalar const&); // A[l,m]+= x*B[l,m]     
    template<symmetry> YlmRec&subtimes(YlmRec const&,
				       scalar const&); // A[l,m]-= x*B[l,m]     
  }; // class YlmRec
  //----------------------------------------------------------------------------
  //                                                                            
  // sub-type AnlRec                                                            
  //                                                                            
  // holds the numbers a_nl with n=0...N, l=0...L                               
  //                                                                            
  // NOTE                                                                       
  //    Depending on the adopted symmetry, we ASSUME certain elements to be     
  //    zero or trivial copies, but do NOT ENFORCE this. In particular:         
  //    1. spherical symmetry: only the l=0 elements are used.                  
  //    2. lower symmetry: only elements with even l are used.                  
  //    3. no symmetry: all elements are used                                   
  class AnlRec {
    friend class PotExpAccess;
    int     N1,L1;
    scalar *A;
  public:
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // construction & destruction                                               
    AnlRec(int n, int l)
      : N1(n+1), L1(l+1), A(falcON_NEW(scalar,N1*L1)) {}
    AnlRec(AnlRec const&a)
      : N1(a.N1), L1(a.L1), A(falcON_NEW(scalar,N1*L1)) {
      memcpy(A,a.A,N1*L1*sizeof(scalar));
    }
    ~AnlRec() { falcON_DEL_A(A); }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // element access                                                           
    scalar      &operator() (int n, int l)       { return A[n*L1+l]; }
    scalar const&operator() (int n, int l) const { return A[n*L1+l]; }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void table_print(symmetry, std::ostream&, int=6) const;
  private: // only accessible to PotExp in file pexp.cc                         
    template<symmetry> AnlRec&reset   ();              // A[n,l] = 0            
    template<symmetry> AnlRec&assign  (scalar const&); // A[n,l] = x            
    template<symmetry> AnlRec&multiply(scalar const&); // A[n,l]*= x            
    template<symmetry> AnlRec&divide  (scalar const&); // A[n,l]/= x            
    template<symmetry> AnlRec&copy    (AnlRec const&); // A[n,l] = B[n,l]       
    template<symmetry> AnlRec&add     (AnlRec const&); // A[n,l]+= B[n,l]       
    template<symmetry> AnlRec&sub     (AnlRec const&); // A[n,l]-= B[n,l]       
    template<symmetry> scalar dot     (AnlRec const&) const; // dot product     
    template<symmetry> AnlRec&addtimes(AnlRec const&,
				       scalar const&); // A[n,l]+= x*B[n,l]     
    template<symmetry> AnlRec&subtimes(AnlRec const&,
				       scalar const&); // A[n,l]-= x*B[n,l]     
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // static variables to hold info about alpha used by static methods           
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  scalar AL=1, iAL=1, AL1=2;
  inline void setAL(scalar const&al) {
    AL = al;
    iAL= 1/al;
    AL1= 1+al;
  }
} // namespace {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // variables to hold info about scale radius used by methods in this file and 
  // in file PotExpSSE.cc                                                       
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
namespace falcON { namespace P {
  scalar R0, IR0;
#ifdef falcON_SSE
  void setIR04(scalar);
#endif
  inline void setR0(scalar r0) {
    R0  = r0;
    IR0 = 1/r0;
#ifdef falcON_SSE
    setIR04(IR0);
#endif
  }
} } // namespace falcON::P
namespace {
  using namespace falcON::P;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // auxiliary inline functions                                                 
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  inline symmetry operator & (symmetry a, symmetry b) {
    return PotExp::symmetry(int(a) & int(b));
  }
  //----------------------------------------------------------------------------
  inline symmetry operator | (symmetry a, symmetry b) {
    return PotExp::symmetry(int(a) | int(b));
  }
  //----------------------------------------------------------------------------
  inline scalar lambda(int l) {
    return AL*(l+l+1)+0.5;
  }
  //----------------------------------------------------------------------------
  inline void SetXiFi(scalar&xi, scalar&fi, scalar const&r)
  {
    // sets                                       
    //          s-1         2                     
    //   xi := ----- = 1 - ----                   
    //          s+1        s+1                    
    // and                                        
    //            1                               
    //   fi := -------                            
    //         (s+1)^a                            
    // with                                       
    //         (1/a)                              
    //   s := r                                   
    //                                            
    if       (AL == 0.5) {
      fi = 1/(r*r+1);
      xi = 1-fi-fi;
      fi = sqrt(fi);
    } else if(AL == 1.0) {
      fi = 1/(r+1);
      xi = 1-fi-fi;
    } else if(AL == 2.0) {
      fi = 1/(sqrt(r)+1);
      xi = 1-fi-fi;
      fi*= fi;
    } else {
      fi = 1/(pow(r,iAL)+1);
      xi = 1-fi-fi;
      fi = pow(fi,AL);
    }
  }
  //----------------------------------------------------------------------------
  inline void SetXiFi(scalar&xi, scalar&dxi,
		      scalar&fi, scalar&dfi,
		      scalar const&r)
  {
    // sets                                       
    //          s-1         2                     
    //   xi := ----- = 1 - ----,                  
    //          s+1        s+1                    
    // and                                        
    //            1                               
    //   fi := -------                            
    //         (s+1)^a                            
    //                                            
    // as well as                                 
    //                                            
    //          d xi     2 s/r                    
    //   dxi := ---- = ---------,                 
    //          d r    a (s+1)^2                  
    // and                                        
    //          d fi       s/r                    
    //   dfi := ---- = - -----------              
    //          d r      (s+1)^(a+1)              
    // with                                       
    //         (1/a)                              
    //   s := r                                   
    //                                            
    if       (AL == 0.5) {
      dfi = 1/(r*r+1);                             //  1/(s+1)                  
      fi  = sqrt(dfi);                             //  1/(s+1)^(1/2)            
      xi  = 1-dfi-dfi;                             //  1-2/(s+1) = (s-1)/(s+1)  
      dxi = 4*dfi;                                 //  4/(s+1)                  
      dfi*= r;                                     //  r/(s+1)                  
      dxi*= dfi;                                   //  4*r/(s+1)^2              
      dfi*=-fi;                                    // -r/(s+1)^(3/2)            
    } else if(AL == 1.0) {
      fi  = 1/(r+1);                               //  1/(r+1)                  
      xi  = 1-fi-fi;                               //  1-2/(r+1) = (r-1)/(r+1)  
      dxi = fi*fi;                                 //  1/(r+1)^2                
      dfi =-dxi;                                   // -1/(r+1)^2                
      dxi+= dxi;                                   //  2/(r+1)^2                
    } else if(AL == 2.0) {
      const scalar s=sqrt(r);                      //  s := r^(1/a) = sqrt(r)   
      fi  = 1/(s+1);                               //  1/(s+1)                  
      xi  = 1-fi-fi;                               //  1-2/(s+1) = (s-1)/(s+1)  
      dfi = fi;                                    //  1/(s+1)                  
      fi *= fi;                                    //  1/(s+1)^2                
      dxi = fi/s;                                  //  1/s/[s+1]^2=s/r/[s+1]^2  
      dfi*=-dxi;                                   // -s/r/[s+1]^3              
    } else {
      const scalar s=pow(r,iAL);                   //  s = r^(1/a)              
      fi  = 1/(s+1);                               //  1/(s+1)                  
      xi  = 1-fi-fi;                               //  1-2/(r+1) = (r-1)/(r+1)  
      dfi = s*fi/r;                                //  s/r/(s+1)                
      dxi = iAL*(dfi+dfi)*fi;                      //  2*s/r/(s+1)^2/a          
      fi  = pow(fi,AL);                            //  1/(s+1)^a                
      dfi*=-fi;                                    // -s/r/(s+1)(1+a)           
    }
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // PotExpAccess                                                             //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class PotExpAccess {
  protected:
    static int const &L  (CasRec const&C) { return C.L; }
    static int const &L  (YlmRec const&Y) { return Y.L; }
    static int const &L1 (YlmRec const&Y) { return Y.L1; }
    static int const &L1Q(YlmRec const&Y) { return Y.L1Q; }
    static int const &N1 (AnlRec const&A) { return A.N1; }
    static int const &L1 (AnlRec const&A) { return A.L1; }
    static int const &N  (Anlm   const&A) { return A.N; }
    static int const &N1 (Anlm   const&A) { return A.N1; }
    static int const &L  (Anlm   const&A) { return A.L; }
    static int const &L1 (Anlm   const&A) { return A.L1; }
    static int const &L1Q(Anlm   const&A) { return A.L1Q; }

    static scalar      &e(CasRec      &C, int i) { return C.A[i]; }
    static scalar const&e(CasRec const&C, int i) { return C.A[i]; }
    static scalar      &e(YlmRec      &C, int i) { return C.A[i]; }
    static scalar const&e(YlmRec const&C, int i) { return C.A[i]; }
    static scalar      &e(AnlRec      &C, int i) { return C.A[i]; }
    static scalar const&e(AnlRec const&C, int i) { return C.A[i]; }
    static scalar      &e(Anlm        &C, int i) { return C.A[i]; }
    static scalar const&e(Anlm   const&C, int i) { return C.A[i]; }

    static scalar      *p(CasRec      &C) { return C.A; }
    static scalar const*p(CasRec const&C) { return C.A; }
    static scalar      *p(YlmRec      &C) { return C.A; }
    static scalar const*p(YlmRec const&C) { return C.A; }
    static scalar      *p(AnlRec      &C) { return C.A; }
    static scalar const*p(AnlRec const&C) { return C.A; }
    static scalar      *p(Anlm        &C) { return C.A; }
    static scalar const*p(Anlm   const&C) { return C.A; }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // symmetry dependent operations                                            //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<symmetry> struct AUX;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // 1. the general case: no symmetry                                         //
  //                                                                          //
  // n=0...L; l=0...L; m=-l...l                                               //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<> struct AUX<PotExp::none> : private PotExpAccess {
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(CasRec&A, CasRec const&B, scalar const&x) {
      for(int m=-L(A); m<=L(A); ++m)
	Connector::op(e(A,m),e(B,m),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(YlmRec&A, YlmRec const&B, scalar const&x) {
      for(int i=0; i!=L1Q(A); ++i)
	Connector::op(e(A,i),e(B,i),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(AnlRec&A, AnlRec const&B, scalar const&x) {
      for(int i=0; i!=N1(A)*L1(A); ++i)
	Connector::op(e(A,i),e(B,i),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(Anlm&A, Anlm const&B, scalar const&x) {
      for(int i=0; i!=N1(A)*L1Q(A); ++i)
	Connector::op(e(A,i),e(B,i),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(Anlm&A, AnlRec const&P, YlmRec const&Y,
			scalar const&x) {
      register       scalar*Ai = p(A);
      register const scalar*Pi = p(P);
      for(int n=0; n!=N1(A); ++n) {
	register const scalar*Yi = p(Y);
	for(int l=0; l<L1(A); ++l,++Pi) {
	  const scalar Pnl = *Pi;
	  for(int m=-l; m<=l; ++m,++Ai,++Yi)
	    Connector::op(*Ai, Pnl* *Yi, x);
	}
      }
    }
    //--------------------------------------------------------------------------
    static scalar Dot(CasRec const&A, CasRec const&B) {
      register scalar x=0.;
      for(int m=-L(A); m<=L(A); ++m)
	x+=e(A,m)*e(B,m);
      return x;
    }
    //--------------------------------------------------------------------------
    static scalar Dot(YlmRec const&A, YlmRec const&B) {
      register scalar x=0.;
      for(int i=0; i!=L1Q(A); ++i)
	x += e(A,i)*e(B,i);
      return x;
    }
    //--------------------------------------------------------------------------
    static scalar Dot(AnlRec const&A, AnlRec const&B) {
      register scalar x=0.;
      for(int i=0; i!=N1(A)*L1(A); ++i)
	x += e(A,i)*e(B,i);
      return x;
    }
    //--------------------------------------------------------------------------
    static scalar Dot(Anlm const&A, Anlm const&B) {
      register scalar x=0.;
      for(int i=0; i!=N1(A)*L1Q(A); ++i)
	x += e(A,i)*e(B,i);
      return x;
    }
    //--------------------------------------------------------------------------
    static scalar Dot(Anlm const&A, AnlRec const&P, YlmRec const&Y) {
      // P = Sum_nlm A_nlm * P_nl * Y_lm                                        
      scalar       x=0;
      const scalar*Ai = p(A);
      const scalar*Pi = p(P);
      for(int n=0; n!=N1(A); ++n) {
	const scalar*Yi = p(Y);
	for(int l=0; l<L1(A); ++l,++Pi) {
	  register scalar y=0;
	  for(int m=-l; m<=l; ++m,++Ai,++Yi)
	    y += *Ai * *Yi;
	  x += *Pi * y;
	}
      }
      return x;
    }
    //--------------------------------------------------------------------------
    template<typename T>
    static scalar Dot(tupel<3,T>&dx,
		      Anlm   const&A, AnlRec const&P, AnlRec const&R, 
		      YlmRec const&Y, YlmRec const&T, YlmRec const&Q) {
      // P = Sum_nlm A_nlm * P_nl * Y_lm                                        
      // dP/d(r,the,phi)                                                        
      scalar x=0,xr=0,xt=0,xp=0;
      const scalar*Ai = p(A);
      const scalar*Pi = p(P);
      const scalar*Ri = p(R);
      for(int n=0; n!=N1(A); ++n) {
	const scalar*Yi = p(Y);
	const scalar*Ti = p(T);
	const scalar*Qi = p(Q);
	for(int l=0; l<L1(A); ++l,++Pi,++Ri) {
	  register scalar y=0,yt=0,yp=0;
	  for(int m=-l; m<=l; ++m,++Ai,++Yi,++Ti,++Qi) {
	    y  += *Ai * *Yi;
	    yt += *Ai * *Ti;
	    yp += *Ai * *Qi;
	  }
	  x  += *Pi * y;
	  xr += *Ri * y;
	  xt += *Pi * yt;
	  xp += *Pi * yp;
	}
      }
      dx[0] = xr;
      dx[1] = xt;
      dx[2] = xp;
      return x;
    }
    //--------------------------------------------------------------------------
    static void SetCas(CasRec&A, scalar const&c, scalar const&s) {
      e(A,0) = 1;
      for(int m=0,m1=1; m<L(A); ++m,++m1) {
	e(A, m1) = c * e(A, m) + s * e(A,-m);
	e(A,-m1) = c * e(A,-m) - s * e(A, m);
      }
    }
    //--------------------------------------------------------------------------
    static void SetPlm(YlmRec&Y, scalar const&ct, scalar const&st) {
      e(Y,0) = 1;                                  // P(0,0) = 1                
      // 1. compute P(l,l) using the GR 8.731.4                                 
      for(int l = 0,                               // LOOP l=0...L-1            
	    tlp = 1,                               //   2*l+1                   
	    i   = 0,                               //   i(l  ,l  ) = l*(l+2)    
	    ip  = 3;                               //   i(l+1,l+1) = (l+1)*(l+3)
	  l < L(Y);                                //   l+1 to run to L         
	  ++l,                                     //   l          <- l+1       
	    tlp+= 2,                               //   2*l+1      <- 2*l+1 + 2 
	    i   = ip,                              //   i(l  ,l  ) <- i(l+1,l+1)
	    ip += tlp+2)                           //   i(l+1,l+1) <- i(l+2,l+2)
	e(Y,ip) = -tlp * st * e(Y,i);              //   GR 8.731.4              
      // 2. compute P(l,m) using GR 8.731.4(1)                                  
      for(int m = 0; m<L(Y); ++m)                  // LOOP m=0...L,1            
	for(int l = m,                             //   LOOP l=m...L-1,1        
	      tlp = m+m+1,                         //     2*l+1                 
	      lpm = m+m,                           //     l+m                   
	      lmp = 1,                             //     l-m+1                 
	      im  = m*m,                           //     i(l-1,m)=l*(l-1)+m    
	      i   = m*(m+2),                       //     i(l  ,m)=l*(l+1)+m    
	      ip  = m+(m+1)*(m+2);                 //     i(l+1,m)=(l+1)*(l+2)+m
	    l < L(Y);                              //     l+1 to run to L       
	    ++l,                                   //     l        <- l+1       
	      tlp+= 2,                             //     2*l+1    <- 2*l+1 + 2 
	      ++lpm,                               //     l+m      <- l+1+m     
	      ++lmp,                               //     l-m+1    <- l-m-1+1   
	      im  = i,                             //     i(l-1,m) <- i(l  ,m)  
	      i   = ip,                            //     i(l  ,m) <- i(l+1,m)  
	      ip += tlp+1)                         //     i(l+1,m) <- i(l+2,m)  
	  if(l==m)                                 //     IF l==m:              
	    e(Y,ip) = tlp*ct*e(Y,i);               //       set P(m+1,m)        
	  else                                     //     ELSE (l>m):           
	    e(Y,ip) =(tlp*ct*e(Y,i) - lpm*e(Y,im)) / scalar(lmp);
    }
    //--------------------------------------------------------------------------
    static void SetPlm(YlmRec&Y, YlmRec&T, scalar const&ct, scalar const&st) {
      e(Y,0) = 1;                                  //  P(0,0)      = 1          
      e(T,0) = 0;                                  // dP(0,0)/dthe = 0          
      for(int l = 0,                               // LOOP l=0...L-1            
	    tlp = 1,                               //   2*l+1                   
	    i   = 0,                               //   i(l  ,l  ) = l*(l+2)    
	    ip  = 3;                               //   i(l+1,l+1) = (l+1)*(l+3)
	  l < L(Y);                                //   l+1 to run to L         
	  ++l,                                     //   l          <- l+1       
	    tlp+= 2,                               //   2*l+1      <- 2*l+1 + 2 
	    i   = ip,                              //   i(l  ,l  ) <- i(l+1,l+1)
	    ip += tlp+2) {                         //   i(l+1,l+1) <- i(l+2,l+2)
	e(Y,ip) = -tlp* st*e(Y,i);                 //   GR 8.731.4              
	e(T,ip) = -tlp*(st*e(T,i)+ct*e(Y,i));      //   dP(l+1,l+1)/dthe        
      }                                            // END LOOP                  
      for(int m = 0; m<L(Y); ++m)                  // LOOP m=0...L,1            
	for(int l = m,                             //   LOOP l=m...L-1          
	      tlp = m+m+1,                         //     2*l+1                 
	      lpm = m+m,                           //     l+m                   
	      lmp = 1,                             //     l-m+1                 
	      im  = m*m,                           //     i(l-1,m)=l*(l-1)+m    
	      i   = m*(m+2),                       //     i(l  ,m)=l*(l+1)+m    
	      ip  = m+(m+1)*(m+2);                 //     i(l+1,m)=(l+1)*(l+2)+m
	    l < L(Y);                              //     l+1 to run to L       
	    ++l,                                   //     l        <- l+1       
	      tlp+= 2,                             //     2*l+1    <- 2*l+1 + 2 
	      ++lpm,                               //     l+m      <- l+1+m     
	      ++lmp,                               //     l-m+1    <- l-m-1+1   
	      im  = i,                             //     i(l-1,m) <- i(l  ,m)  
	      i   = ip,                            //     i(l  ,m) <- i(l+1,m)  
	      ip += tlp+1) {                       //     i(l+1,m) <- i(l+2,m)  
	  if(l==m) {                               //     IF l==m:              
	    e(Y,ip) = tlp* ct*e(Y,i);              //        P(m+1,m)           
	    e(T,ip) = tlp*(ct*e(T,i)-st*e(Y,i));   //       dP(m+1,m)/dthe      
	  } else {                                 //     ELSE (l>m):           
	    const scalar ilmp = 1/scalar(lmp);     //       1/(l-m+1)           
	    e(Y,ip)  = ( tlp* ct*e(Y,i)          -lpm*e(Y,im) ) * ilmp;
	    e(T,ip)  = ( tlp*(ct*e(T,i)-st*e(Y,i)) -lpm*e(T,im) ) * ilmp;
	  }
	}                                          // END LOOPS(m,l)            
    }
    //--------------------------------------------------------------------------
    static void SetYlm(YlmRec&Y,
		       scalar const&ct, scalar const&st,
		       scalar const&cp, scalar const&sp) {
      SetPlm(Y,ct,st);                             // set Plm(theta)            
      register scalar _Cp=1.,Cp;                   // cas( m*phi)               
      register scalar _Cm=1.,Cm;                   // cas(-m*phi)               
      for(int m=1; m<L1(Y); ++m) {                 // LOOP m=1...L              
	Cp = cp * _Cp + sp * _Cm;                  //   C( m) := cas( m*phi)    
	Cm = cp * _Cm - sp * _Cp;                  //   C(-m) := cas(-m*phi)    
	for(int l = m,                             //   LOOP l=m...L            
	      im  = m*m,                           //     i(l,-m) = l*(l+1)-m   
	      ip  = im+m+m;                        //     i(l, m) = l*(l+1)+m   
	    l < L1(Y);                             //     l to run to L         
	    ++l,                                   //     l       <- l+1        
	      im += l+l,                           //     i(l,-m) <- i(l+1,-m)  
	      ip += l+l) {                         //     i(l, m) <- i(l+1, m)  
	  e(Y,im) = Cm * e(Y,ip);                  //     Y(l,-m)=P(l,m)*C(-m)  
	  e(Y,ip)*= Cp;                            //     Y(l, m)=P(l,m)*C( m)  
	}                                          //   END LOOP(l)             
	_Cm = Cm;                                  //   cas[-(m-2)) <- C(-m)    
	_Cp = Cp;                                  //   cas[  m-2 ) <- C( m)    
      }                                            // END LOOP(m)               
    }
    //--------------------------------------------------------------------------
    static void SetYlm(YlmRec&Y, YlmRec&T, YlmRec&P,
		       scalar const&ct, scalar const&st,
		       scalar const&cp, scalar const&sp) {
      SetPlm(Y,T,ct,st);                           // set Plm(the), dPlm/dthe   
      register scalar _Cp=1.,Cp;                   // cas( m*phi)               
      register scalar _Cm=1.,Cm;                   // cas(-m*phi)               
      for(int l=0; l<L1(Y); ++l) e(P,l*(l+1))=0;   // dY(l,m=0)/dp=0            
      for(int m=1; m<L1(Y); ++m) {                 // LOOP m=1...L              
	Cp = cp * _Cp + sp * _Cm;                  //   C( m) := cas( m*phi)    
	Cm = cp * _Cm - sp * _Cp;                  //   C(-m) := cas(-m*phi)    
	for(int l = m,                             //   LOOP l=m...L            
	      im  = m*m,                           //     i(l,-m) = l*(l+1)-m   
	      ip  = im+m+m;                        //     i(l, m) = l*(l+1)+m   
	    l < L1(Y);                             //     l to run to L         
	    ++l,                                   //     l       <- l+1        
	      im += l+l,                           //     i(l,-m) <- i(l+1,-m)  
	      ip += l+l) {                         //     i(l, m) <- i(l+1, m)  
	  e(P,im) = -m * Cp * e(Y,ip);             //     dY(l,-m)/dp=-m*P*C(+m)
	  e(P,ip) =  m * Cm * e(Y,ip);             //     dY(l,+m)/dp= m*P*C(-m)
	  e(Y,im) = Cm * e(Y,ip);                  //      Y(l,-m)   =   P*C(-m)
	  e(Y,ip)*= Cp;                            //      Y(l,+m)   =   P*C(+m)
	  e(T,im) = Cm * e(T,ip);                  //     dY(l,-m)/dt= P_t*C(-m)
	  e(T,ip)*= Cp;                            //     dY(l,+m)/dt= P_t*C(+m)
	}                                          //   END LOOP(l)             
	_Cm = Cm;                                  //   C(-(m-2)) <- C(-m)      
	_Cp = Cp;                                  //   C(  m-2 ) <- C( m)      
      }                                            // END LOOP(m)               
    }
    //--------------------------------------------------------------------------
    static void SetPsi(AnlRec&P, scalar const&r, scalar const&GM) {
      // routine checked (08/07/04) against MAPLE                               
      //                                                                        
      //                            l                          (1/a)            
      //                     G*M * r           (2*l+1)*a+1/2  r      - 1        
      // sets P(n,l) = ---------------------  G              (----------)       
      //                 (1/a)     (2*l+1)*a   n               (1/a)            
      //               (r       + 1)                          r      + 1        
      //                                                                        
      // where G denote the Gegenbauer polynomials.                             
      // see also Zhao (1996, MNRAS 278, 488)                                   
      register scalar xi;                          // (r^(1/a)-1)/(r^(1/a)+1)   
      register scalar fi;                          // 1/(1+r^(1/a))^a           
      SetXiFi(xi,fi,r);
      e(P,0) = GM*fi;                              // P(0,0) = G*M/(1+r^(1/a))^a
      fi  *= r*fi;                                 // fi = r / (1+r^(1/a))^(2a) 
      for(int l=0,lp=1; lp<L1(P); ++l,++lp)        // LOOP lp=dl...L,1          
	e(P,lp) = fi * e(P,l);                     //   P(0,l+1) = fi*P(0,l)    
      if(N1(P)==1) return;                         // no terms n>0 --> DONE     
      register scalar tlm = twice(lambda(0));      // 2*lambda(l=0)             
      const    scalar Dtl = 4*AL;                  // 2*(lambda(l+1)-lambda(l)) 
      const    scalar xi2 = xi+xi;                 // 2*xi                      
      for(int l=0,ln=L1(P);                        // LOOP l=0...L,1            
	  l<L1(P);                                 //   l+1 runs to L           
	  ++l,++ln,tlm+=Dtl) {                     //   increment stuff         
	e(P,ln) = tlm * xi * e(P,l);               //   P(1,l) = 2*lam*xi*P(0,l)
	register scalar tlmn2xi = (tlm+2)*xi;      //   2*(lam+n)*xi            
	register scalar tlm1n   = tlm;             //   2*lam+n-1               
	for(int n1= 2,                             //   LOOP n+1=2...N          
	      im  = l,                             //     i(n-1,l)              
	      i   = im+L1(P),                      //     i(n  ,l)              
	      ip  = i +L1(P);                      //     i(n+1,l)              
	    n1 < N1(P);                            //     n+1 runs til N        
	    ++n1,                                  //     n        <- n+1       
	      im  = i,                             //     i(n-1,l) <- i(n  ,l)  
	      i   = ip,                            //     i(n  ,l) <- i(n+1,l)  
	      ip += L1(P),                         //     i(n+1,l) += L+1       
	      tlmn2xi += xi2,                      //     2*(lam+n)*xi          
	      ++tlm1n)                             //     2*lam+n-1             
	  e(P,ip) =(tlmn2xi*e(P,i) - tlm1n*e(P,im))/scalar(n1); // GR 8.933.1   
      }                                            // END LOOPS (l,n)           
    }
    //--------------------------------------------------------------------------
    static void SetPsi(AnlRec&P, AnlRec&D, scalar const&r) {
      // sets Psi_nl(r) and d Psi_nl(r) / dr                                    
      register scalar xi,dxi;                      // (r^(1/a)-1)/(r^(1/a)+1)   
      register scalar fi,dfi;                      // 1/(1+r^(1/a))^a           
      SetXiFi(xi,dxi,fi,dfi,r);
      e(P,0) = fi;                                 // P(0,0) = 1/(1+r^(1/a))^a  
      e(D,0) =dfi;                                 // D(0,0) = dP(0,0)/dr       
      dfi *= twice(r*fi); fi  *= fi;               // fi  -> r*fi^2             
      dfi += fi;          fi  *= r;                // dfi -> fi^2 + 2*r*fi*dfi  
      for(int l=0,lp=1; lp<L1(P); ++l,++lp) {      // LOOP lp=dl...L,1          
	e(P,lp) = fi*e(P,l);                       //   P(0,l+1) = fi*P(0,l)    
	e(D,lp) = fi*e(D,l) + dfi*e(P,l);          //   D(0,l+1) = dP(0,l+1)/dr 
      }                                            // END LOOP(l)               
      if(N1(P)==1) return;                         // no terms n>0 --> DONE     
      register scalar tlm  = twice(lambda(0));     // 2*lambda(l=0)             
      const    scalar Dtl  = 4*AL;                 // 2*(lambda(l+1)-lambda(l)) 
      const    scalar xi2  = xi+xi;                // 2*xi                      
      const    scalar dxi2 = dxi+dxi;              // d(2*xi)/dr                
      for(int l=0,ln=L1(P);                        // LOOP l=0...L,1            
	  l<L1(P);                                 //   l+1 runs to L           
	  ++l,++ln,tlm+=Dtl) {                     //   increment stuff         
	e(P,ln) = tlm * xi*e(P,l);                 //   P(1,l) = 2*lam*xi*P(0,l)
	e(D,ln) = tlm *(xi*e(D,l)+dxi*e(P,l));     //   D(1,l) = dP(1,l)/dr     
	register scalar tlmn2xi  = (tlm+2)*xi;     //   2*(lam+n)*xi            
	register scalar dtlmn2xi = (tlm+2)*dxi;    //   d(2*(lam+n)*xi)/dr      
	register scalar tlm1n    = tlm;            //   2*lam+n-1               
	for(int n1= 2,                             //   LOOP n+1=2...N          
	      im  = l,                             //     i(n-1,l)              
	      i   = im+L1(P),                      //     i(n  ,l)              
	      ip  = i +L1(P);                      //     i(n+1,l)              
	    n1 < N1(P);                            //     n+1 runs til N        
	    ++n1,                                  //     n        <- n+1       
	      im  = i,                             //     i(n-1,l) <- i(n  ,l)  
	      i   = ip,                            //     i(n  ,l) <- i(n+1,l)  
	      ip += L1(P),                         //     i(n+1,l) += L+1       
	      tlmn2xi += xi2,                      //     2*(lam+n)*xi          
	      dtlmn2xi+= dxi2,                     //     2*(lam+n)*dxi         
	      ++tlm1n) {                           //     2*lam+n-1             
	  const scalar in1 = 1/scalar(n1);         //     see GR 8.933.1        
	  e(P,ip) =(tlmn2xi*e(P,i)                 - tlm1n*e(P,im)) * in1;
	  e(D,ip) =(tlmn2xi*e(D,i) + dtlmn2xi*e(P,i) - tlm1n*e(D,im)) * in1;
	}                                          //   END LOOP(n)             
      }                                            // END LOOP(l)               
    }
    //--------------------------------------------------------------------------
    static void SetKnl (AnlRec&K) {
      // computes the normalisation constants Knl
      const scalar iA1 = iAL / AL1;
      for(int l=0; l!=L1(K); ++l) {
	const scalar lm = lambda(l);
	for(int n=0; n!=N1(K); ++n)
	  K(n,l) = iA1 * (square(n+lm) - 0.25);
      }
    }
    //--------------------------------------------------------------------------
    static void SetAnl (AnlRec&A, AnlRec const&K) {
      // computes the A_nl = 1/I_nl/r0^2 given K_nl
      const scalar ln2 = log(scalar(2));
      for(int l=0; l!=L1(A); ++l) {
	const scalar lm = lambda(l);
	A(0,l) = lm * exp(4*lm*ln2+2*lgamma(lm)-lgamma(lm+lm));
	for(int n=0,n1=1; n1<N1(A); ++n,++n1)
	  A(n1,l) = n1*(n1+lm)/((n+lm)*(n+lm+lm)) * A(n,l);
      }
      const scalar t1 = square(IR0) / (TPi * AL1);
      for(int n=0; n!=N1(A); ++n)
	for(int l=0; l!=L1(A); ++l)
	  A(n,l) *= t1 / K(n,l);
    }
    //--------------------------------------------------------------------------
    static void SetNlm (YlmRec&Y) {
      // computing the N_lm multiplied by 4Pi
      // i.e. on return Y(l,m) holds (2l+1) (l-|m|)! / (l+|m|)!
      Y(0,0) = 1.;
      for(int l=1; l<L1(Y); ++l) {
	Y(l,0) = Y(l-1,0) + 2.;
	for(int m=0,m1=1; m<l; ++m,++m1)
	  Y(l,-m1) = Y(l,m1) = Y(l,m) / scalar((l+m1)*(l-m));
      }
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // 2. reflexion symmetry                                                    //
  //                                                                          //
  // n=0...L; l=0...L,2; m=-l...l,2                                           //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<> struct AUX<PotExp::reflexion> : private PotExpAccess {
    //--------------------------------------------------------------------------
    template<class Connector> static
    void Connect(CasRec&A, CasRec const&B, scalar const&x) {
      for(int m=-L(A); m<=L(A); m+=2)
	Connector::op(e(A,m),e(B,m),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector> static
    void Connect(YlmRec&A, YlmRec const&B, scalar const&x) {
      for(int l=0,i=0; l<L1(A); l+=2,i+=4*l-2)
	for(int m=-l; m<=l; m+=2)
	  Connector::op(e(A,i+m),e(B,i+m),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(AnlRec&A, AnlRec const&B, scalar const&x)
    {
      for(int n=0,j=0; n!=N1(A); ++n,j+=L1(A))
	for(int l=0,i=j; l<L1(A); l+=2,i+=2)
	  Connector::op(e(A,i),e(B,i),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(Anlm&A, Anlm const&B, scalar const&x)
    {
      for(int n=0,j=0; n!=N1(A); ++n,j+=L1Q(A))
	for(int l=0; l<L1(A); l+=2)
	  for(int m=-l,i=j+l*l; m<=l; m+=2,i+=2)
	    Connector::op(e(A,i),e(B,i),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(Anlm&A, AnlRec const&P, YlmRec const&Y,
			scalar const&x) {
      scalar      *An = p(A);
      const scalar*Pn = p(P);
      for(int n=0; n!=N1(A); ++n, An+=L1Q(A), Pn+=L1(A)) {
	scalar      *Al = An;
	const scalar*Pl = Pn;
	const scalar*Yl = p(Y);
	for(int l=0; l<L1(A); l+=2,Al+=4*l-2,Pl+=2,Yl+=4*l-2) {
	  const scalar Pnl = *Pl;
	  for(int m=-l; m<=l; m+=2)
	    Connector::op(Al[m], Pnl*Yl[m], x);
	}
      }
    }
    //--------------------------------------------------------------------------
    static scalar Dot(CasRec const&A, CasRec const&B) {
      register scalar x=0.;
      for(int m=-L(A); m<=L(A); m+=2)
	x += e(A,m)*e(B,m);
      return x;
    }
    //--------------------------------------------------------------------------
    static scalar Dot(YlmRec const&A, YlmRec const&B) {
      register scalar x=0.;
      for(int l=0,i=0; l<L1(A); l+=2,i+=4*l-2)
	for(int m=-l; m<=l; m+=2)
	  x += e(A,i+m) * e(B,i+m);
      return x;
    }
    //--------------------------------------------------------------------------
    static scalar Dot(AnlRec const&A, AnlRec const&B) {
      register scalar x=0.;
      for(int n=0,j=0; n!=N1(A); ++n,j+=L1(A))
	for(int l=0,i=j; l<L1(A); l+=2,i+=2)
	  x += e(A,i) * e(B,i);
      return x;
    }
    //--------------------------------------------------------------------------
    static scalar Dot(Anlm const&A, Anlm const&B) {
      register scalar x=0.;
      for(int n=0,j=0; n!=N1(A); ++n,j+=L1Q(A))
	for(int l=0; l<L1(A); l+=2)
	  for(int m=-l,i=j+l*l; m<=l; m+=2,i+=2)
	    x += e(A,i) * e(B,i);
      return x;
    }
    //--------------------------------------------------------------------------
    static scalar Dot(Anlm const&A, AnlRec const&P, YlmRec const&Y) {
      // P = Sum_nlm A_nlm * P_nl * Y_lm                                        
      scalar       x=0;
      const scalar*An = p(A);
      const scalar*Pn = p(P);
      for(int n=0; n!=N1(A); ++n, An+=L1Q(A), Pn+=L1(A)) {
	const scalar*Al = An;
	const scalar*Pi = Pn;
	const scalar*Yl = p(Y);
	for(int l=0; l<L1(A); l+=2,Al+=4*l-2,Pi+=2,Yl+=4*l-2) {
	  register scalar y=0;
	  for(int m=-l; m<=l; m+=2)
	    y += Al[m] * Yl[m];
	  x += *Pi * y;
	}
      }
      return x;
    }
    //--------------------------------------------------------------------------
    template<typename T>
    static scalar Dot(tupel<3,T>&dx,
		      Anlm   const&A, AnlRec const&P, AnlRec const&R, 
		      YlmRec const&Y, YlmRec const&T, YlmRec const&Q) {
      // P = Sum_nlm A_nlm * P_nl * Y_lm                                        
      // dP/d(r,the,phi)                                                        
      scalar x=0,xr=0,xt=0,xp=0;
      const scalar*An = p(A);
      const scalar*Pn = p(P);
      const scalar*Rn = p(R);
      for(int n=0; n!=N1(A); ++n,An+=L1Q(A),Pn+=L1(A),Rn+=L1(A)) {
	const scalar*Al = An;
	const scalar*Pi = Pn;
	const scalar*Ri = Rn;
	const scalar*Yl = p(Y);
	const scalar*Tl = p(T);
	const scalar*Ql = p(Q);
	for(int l=0; l<L1(A);
	    l+=2,Al+=4*l-2,Pi+=2,Ri+=2,Yl+=4*l-2,Tl+=4*l-2,Ql+=4*l-2) {
	  register scalar y=0,yt=0,yp=0;
	  for(int m=-l; m<=l; m+=2) {
	    y  += Al[m] * Yl[m];
	    yt += Al[m] * Tl[m];
	    yp += Al[m] * Ql[m];
	  }
	  x  += *Pi * y;
	  xr += *Ri * y;
	  xt += *Pi * yt;
	  xp += *Pi * yp;
	}
      }
      dx[0] = xr;
      dx[1] = xt;
      dx[2] = xp;
      return x;
    }
    //--------------------------------------------------------------------------
    static void SetCas(CasRec&A, scalar const&c, scalar const&s) {
      e(A,0) = 1;
      const scalar cc=c*c-s*s, ss=2*c*s;
      for(int m=0,m2=2; m2<=L(A); m+=2,m2+=2) {
	e(A, m2) = cc * e(A, m) + ss * e(A,-m);
	e(A,-m2) = cc * e(A,-m) - ss * e(A, m);
      }
    }
    //--------------------------------------------------------------------------
    static void SetPlm(YlmRec&Y, scalar const&ct, scalar const&st) {
      e(Y,0) = 1;                                  // P(0,0) = 1                
      // 1. compute P(l,l) using the GR 8.731.4                                 
      for(int l = 0,                               // LOOP l=0...L-1            
	    tlp = 1,                               //   2*l+1                   
	    i   = 0,                               //   i(l  ,l  ) = l*(l+2)    
	    ip  = 3;                               //   i(l+1,l+1) = (l+1)*(l+3)
	  l < L(Y);                                //   l+1 to run to L         
	  ++l,                                     //   l          <- l+1       
	    tlp+= 2,                               //   2*l+1      <- 2*l+1 + 2 
	    i   = ip,                              //   i(l  ,l  ) <- i(l+1,l+1)
	    ip += tlp+2)                           //   i(l+1,l+1) <- i(l+2,l+2)
	e(Y,ip) = -tlp*st*e(Y,i);                  //   GR 8.731.4              
      // 2. compute P(l,m) using GR 8.731.4(1)                                  
      for(int m = 0; m<L(Y); m+=2)                 // LOOP m=0...L,2            
	for(int l = m,                             //   LOOP l=m...L-1,1        
	      tlp = m+m+1,                         //     2*l+1                 
	      lpm = m+m,                           //     l+m                   
	      lmp = 1,                             //     l-m+1                 
	      im  = m*m,                           //     i(l-1,m)=l*(l-1)+m    
	      i   = m*(m+2),                       //     i(l  ,m)=l*(l+1)+m    
	      ip  = m+(m+1)*(m+2);                 //     i(l+1,m)=(l+1)*(l+2)+m
	    l < L(Y);                              //     l+1 to run to L       
	    ++l,                                   //     l        <- l+1       
	      tlp+= 2,                             //     2*l+1    <- 2*l+1 + 2 
	      ++lpm,                               //     l+m      <- l+1+m     
	      ++lmp,                               //     l-m+1    <- l-m-1+1   
	      im  = i,                             //     i(l-1,m) <- i(l  ,m)  
	      i   = ip,                            //     i(l  ,m) <- i(l+1,m)  
	      ip += tlp+1)                         //     i(l+1,m) <- i(l+2,m)  
	  if(l==m)                                 //     IF l==m:              
	    e(Y,ip) = tlp*ct*e(Y,i);               //       set P(m+1,m)        
	  else                                     //     ELSE (l>m):           
	    e(Y,ip) =(tlp*ct*e(Y,i)-lpm*e(Y,im)) / scalar(lmp);
    }
    //--------------------------------------------------------------------------
    static void SetPlm(YlmRec&Y, YlmRec&T, scalar const&ct, scalar const&st) {
      e(Y,0) = 1;                                  //  P(0,0)      = 1          
      e(T,0) = 0;                                  // dP(0,0)/dthe = 0          
      for(int l = 0,                               // LOOP l=0...L-1            
	    tlp = 1,                               //   2*l+1                   
	    i   = 0,                               //   i(l  ,l  ) = l*(l+2)    
	    ip  = 3;                               //   i(l+1,l+1) = (l+1)*(l+3)
	  l < L(Y);                                //   l+1 to run to L         
	  ++l,                                     //   l          <- l+1       
	    tlp+= 2,                               //   2*l+1      <- 2*l+1 + 2 
	    i   = ip,                              //   i(l  ,l  ) <- i(l+1,l+1)
	    ip += tlp+2) {                         //   i(l+1,l+1) <- i(l+2,l+2)
	e(Y,ip) = -tlp* st*e(Y,i);                 //    P(l+1,l+1)             
	e(T,ip) = -tlp*(st*e(T,i)+ct*e(Y,i));      //   dP(l+1,l+1)/dthe        
      }                                            // END LOOP                  
      for(int m = 0; m<L(Y); m+=2)                 // LOOP m=0...L,2            
	for(int l = m,                             //   LOOP l=m...L-1          
	      tlp = m+m+1,                         //     2*l+1                 
	      lpm = m+m,                           //     l+m                   
	      lmp = 1,                             //     l-m+1                 
	      im  = m*m,                           //     i(l-1,m)=l*(l-1)+m    
	      i   = m*(m+2),                       //     i(l  ,m)=l*(l+1)+m    
	      ip  = m+(m+1)*(m+2);                 //     i(l+1,m)=(l+1)*(l+2)+m
	    l < L(Y);                              //     l+1 to run to L       
	    ++l,                                   //     l        <- l+1       
	      tlp+= 2,                             //     2*l+1    <- 2*l+1 + 2 
	      ++lpm,                               //     l+m      <- l+1+m     
	      ++lmp,                               //     l-m+1    <- l-m-1+1   
	      im  = i,                             //     i(l-1,m) <- i(l  ,m)  
	      i   = ip,                            //     i(l  ,m) <- i(l+1,m)  
	      ip += tlp+1) {                       //     i(l+1,m) <- i(l+2,m)  
	  if(l==m) {                               //     IF l==m:              
	    e(Y,ip)  = tlp* ct*e(Y,i);             //       P(m+1,m)            
	    e(T,ip)  = tlp*(ct*e(T,i)-st*e(Y,i));  //       dP(m+1,m)/dthe      
	  } else {                                 //     ELSE (l>m):           
	    const scalar ilmp = 1/scalar(lmp);     //       1/(l-m+1)           
	    e(Y,ip) = (tlp* ct*e(Y,i)          -lpm*e(Y,im)) * ilmp;
	    e(T,ip) = (tlp*(ct*e(T,i)-st*e(Y,i)) -lpm*e(T,im)) * ilmp;
	  }
	}                                          // END LOOPS(m,l)            
    }
    //--------------------------------------------------------------------------
    static void SetYlm(YlmRec&Y,
		       scalar const&ct, scalar const&st,
		       scalar const&cp, scalar const&sp) {
      SetPlm(Y,ct,st);                             // set Plm(theta)            
      register scalar _Cp=1.,Cp;                   // cas( m*phi)               
      register scalar _Cm=1.,Cm;                   // cas(-m*phi)               
      const scalar cc=cp*cp-sp*sp, ss=2*cp*sp;
      for(int m=2; m<L1(Y); m+=2) {                // LOOP m=2...L,2            
	Cp = cc * _Cp + ss * _Cm;                  //   C( m) := cas( m*phi)    
	Cm = cc * _Cm - ss * _Cp;                  //   C(-m) := cas(-m*phi)    
	for(int l = m,                             //   LOOP l=m...L,2          
	      di  = 2*(m+m-1),                     //     4*l-2 =i(l,m)-i(l-2,m)
	      im  = m*m,                           //     i(l,-m) = l*(l+1)-m   
	      ip  = im+m+m;                        //     i(l, m) = l*(l+1)+m   
	    l < L1(Y);                             //     l to run to L         
	    l  += 2,                               //     l       <- l+2        
	      di += 8,                             //     4*l-2   <- 4*l-2 + 8  
	      im += di,                            //     i(l,-m) <- i(l,-m) +di
	      ip += di) {                          //     i(l, m) <- i(l, m) +di
	  e(Y,im)  = Cm * e(Y,ip);                 //     Y(l,-m)=P(l,m)*cas[-m)
	  e(Y,ip) *= Cp;                           //     Y(l, m)=P(l,m)*cas[ m)
	}                                          //   END LOOP(l)             
	_Cm = Cm;                                  //   cas[-(m-2)) <- cas[-m)  
	_Cp = Cp;                                  //   cas[  m-2 ) <- cas[ m)  
      }                                            // END LOOP(m)               
    }
    //--------------------------------------------------------------------------
    static void SetYlm(YlmRec&Y, YlmRec&T, YlmRec&P,
		       scalar const&ct, scalar const&st,
		       scalar const&cp, scalar const&sp) {
      SetPlm(Y,T,ct,st);                           // set Plm(the), dPlm/dthe   
      register scalar _Cp=1.,Cp;                   // cas( m*phi)               
      register scalar _Cm=1.,Cm;                   // cas(-m*phi)               
      const scalar cc=cp*cp-sp*sp, ss=2*cp*sp;
      for(int l=0; l<L1(Y); l+=2) e(P,l*(l+1))=0;  // dY(l,m=0)/dp=0            
      for(int m=2; m<L1(Y); m+=2) {                // LOOP m=2...L,2            
	Cp = cc * _Cp + ss * _Cm;                  //   C( m) := cas( m*phi)    
	Cm = cc * _Cm - ss * _Cp;                  //   C(-m) := cas(-m*phi)    
	for(int l = m,                             //   LOOP l=m...L,2          
	      di  = 2*(m+m-1),                     //     4*l-2 =i(l,m)-i(l-2,m)
	      im  = m*m,                           //     i(l,-m) = l*(l+1)-m   
	      ip  = im+m+m;                        //     i(l, m) = l*(l+1)+m   
	    l < L1(Y);                             //     l to run to L         
	    l  += 2,                               //     l       <- l+2        
	      di += 8,                             //     4*l-2   <- 4*l-2 + 8  
	      im += di,                            //     i(l,-m) <- i(l,-m) +di
	      ip += di) {                          //     i(l, m) <- i(l, m) +di
	  e(P,im)  = -m * Cp * e(Y,ip);            //     dY(l,-m)/dp=-m*P*C(+m)
	  e(P,ip)  =  m * Cm * e(Y,ip);            //     dY(l,+m)/dp=+m*P*C(-m)
	  e(Y,im)  = Cm * e(Y,ip);                 //      Y(l,-m)   =   P*C(-m)
	  e(Y,ip) *= Cp;                           //      Y(l,+m)   =   P*C(+m)
	  e(T,im)  = Cm * e(T,ip);                 //     dY(l,-m)/dt= P_t*C(-m)
	  e(T,ip) *= Cp;                           //     dY(l,+m)/dt= P_t*C(+m)
	}                                          //   END LOOP(l)             
	_Cm = Cm;                                  //   C(-(m-2)) <- C(-m)      
	_Cp = Cp;                                  //   C(  m-2 ) <- C( m)      
      }                                            // END LOOP(m)               
    }
    //--------------------------------------------------------------------------
    static void SetPsi (AnlRec&P, scalar const&r, scalar const&GM) {
      // routine checked (08/07/04) against MAPLE                               
      //                                                                        
      //                            l                          (1/a)            
      //                     G*M * r           (2*l+1)*a+1/2  r      - 1        
      // sets P(n,l) = ---------------------  G              (----------)       
      //                 (1/a)     (2*l+1)*a   n               (1/a)            
      //               (r       + 1)                          r      + 1        
      //                                                                        
      // where G denote the Gegenbauer polynomials.                             
      // see also Zhao (1996, MNRAS 278, 488)                                   
      register scalar xi;                          // (r^(1/a)-1)/(r^(1/a)+1)   
      register scalar fi;                          // 1/(1+r^(1/a))^a           
      SetXiFi(xi,fi,r);
      e(P,0) = GM*fi;                              // P(0,0) = G*M/(1+r^(1/a))^a
      fi    *= r*fi;                               // fi = r / (1+r^(1/a))^(2a) 
      fi    *= fi;                                 // fi^2                      
      for(int l=0,lp=2; lp<L1(P); l+=2,lp+=2)      // LOOP lp=dl...L,2          
	e(P,lp) = fi * e(P,l);                     //   P(0,l+2) = fi^2*P(0,l)  
      if(N1(P)==1) return;                         // no terms n>0 --> DONE     
      register scalar tlm = twice(lambda(0));      // 2*lambda(l=0)             
      const    scalar Dtl = 8*AL;                  // 2*(lambda(l+dl)-lambda(l))
      const    scalar xi2 = xi+xi;                 // 2*xi                      
      for(int l=0,ln=L1(P);                        // LOOP l=0...L,1            
	  l<L1(P);                                 //   l+1 runs to L           
	  l+=2,ln+=2,tlm+=Dtl) {                   //   increment stuff         
	e(P,ln) = tlm * xi * e(P,l);               //   P(1,l) = 2*lam*xi*P(0,l)
	register scalar tlmn2xi = (tlm+2)*xi;      //   2*(lam+n)*xi            
	register scalar tlm1n   = tlm;             //   2*lam+n-1               
	for(int n1= 2,                             //   LOOP n+1=2...N          
	      im  = l,                             //     i(n-1,l)              
	      i   = im+L1(P),                      //     i(n  ,l)              
	      ip  = i +L1(P);                      //     i(n+1,l)              
	    n1 < N1(P);                            //     n+1 runs til N        
	    ++n1,                                  //     n        <- n+1       
	      im  = i,                             //     i(n-1,l) <- i(n  ,l)  
	      i   = ip,                            //     i(n  ,l) <- i(n+1,l)  
	      ip += L1(P),                         //     i(n+1,l) += L+1       
	      tlmn2xi += xi2,                      //     2*(lam+n)*xi          
	      ++tlm1n)                             //     2*lam+n-1             
	  e(P,ip) = ( tlmn2xi    * e(P,i) - tlm1n * e(P,im) ) / scalar(n1);
	//P(n+1,0)= ( 2*(lam+n)*xi * P(n,l)-(2*lam+n-1) * P(n-1,l) ) / (n+1)    
	// see GR 8.933.1                                                       
      }                                            // END LOOPS (l,n)           
    }
    //--------------------------------------------------------------------------
    static void SetPsi(AnlRec&P, AnlRec&D, scalar const&r) {
      // sets Psi_nl(r) and d Psi_nl(r) / dr                                    
      register scalar xi,dxi;                      // (r^[1/a]-1)/(r^[1/a]+1)   
      register scalar fi,dfi;                      // 1/(1+r^(1/a))^a           
      SetXiFi(xi,dxi,fi,dfi,r);
      e(P,0) = fi;                                 // P(0,0) = 1/(1+r^(1/a))^a  
      e(D,0) =dfi;                                 // D(0,0) = dP(0,0)/dr       
      dfi *= twice(r*fi); fi  *= fi;               // fi  -> r*fi^2             
      dfi += fi;          fi  *= r;                // dfi -> fi^2 + 2*r*fi*dfi  
      dfi*= fi+fi;                                 //   fi ->fi^2               
      fi *= fi;                                    //   dfi->2*fi*dfi           
      for(int l=0,lp=2; lp<L1(P); l+=2,lp+=2) {    // LOOP lp=dl...L,2          
	e(P,lp) = fi*e(P,l);                       //   P(0,l+2) = fi^2*P(0,l)  
	e(D,lp) = fi*e(D,l) + dfi*e(P,l);          //   D(0,l+2) = dP(0,l+2)/dr 
      }                                            // END LOOP(l)               
      if(N1(P)==1) return;                         // no terms n>0 --> DONE     
      register scalar tlm  = twice(lambda(0));     // 2*lambda(l=0)             
      const    scalar Dtl  = 8*AL;                 // 2*(lambda(l+dl)-lambda(l))
      const    scalar xi2  = xi+xi;                // 2*xi                      
      const    scalar dxi2 = dxi+dxi;              // d(2*xi)/dr                
      for(int l=0,ln=L1(P);                        // LOOP l=0...L,1            
	  l<L1(P);                                 //   l+1 runs to L           
	  l+=2,ln+=2,tlm+=Dtl) {                   //   increment stuff         
	e(P,ln) = tlm* xi*e(P,l);                  //   P(1,l) = 2*lam*xi*P(0,l)
	e(D,ln) = tlm*(xi*e(D,l)+dxi*e(P,l));      //   D(1,l) = dP(1,l)/dr     
	register scalar tlmn2xi  = (tlm+2)*xi;     //   2*(lam+n)*xi            
	register scalar dtlmn2xi = (tlm+2)*dxi;    //   d(2*(lam+n)*xi)/dr      
	register scalar tlm1n    = tlm;            //   2*lam+n-1               
	for(int n1= 2,                             //   LOOP n+1=2...N          
	      im  = l,                             //     i(n-1,l)              
	      i   = im+L1(P),                      //     i(n  ,l)              
	      ip  = i +L1(P);                      //     i(n+1,l)              
	    n1 < N1(P);                            //     n+1 runs til N        
	    ++n1,                                  //     n        <- n+1       
	      im  = i,                             //     i(n-1,l) <- i(n  ,l)  
	      i   = ip,                            //     i(n  ,l) <- i(n+1,l)  
	      ip += L1(P),                         //     i(n+1,l) += L+1       
	      tlmn2xi += xi2,                      //     2*(lam+n)*xi          
	      dtlmn2xi+= dxi2,                     //     2*(lam+n)*dxi         
	      ++tlm1n) {                           //     2*lam+n-1             
	  const scalar in1 = 1/scalar(n1);         //     see GR 8.933.1        
	  e(P,ip) =(tlmn2xi*e(P,i)                 - tlm1n*e(P,im)) * in1;
	  e(D,ip) =(tlmn2xi*e(D,i) + dtlmn2xi*e(P,i) - tlm1n*e(D,im)) * in1;
	}                                          //   END LOOP(n)             
      }                                            // END LOOP(l)               
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // 3. triaxial symmetry                                                     //
  //                                                                          //
  // n=0...L; l=0...L,2; m=0...l,2; assuming A[n,l,-m) = A[n,l,m)             //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<> struct AUX<PotExp::triaxial> : private PotExpAccess {
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(CasRec&A, CasRec const&B, scalar const&x) {
      for(int m=0; m<=L(A); m+=2)
	Connector::op(e(A,m),e(B,m),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(YlmRec&A, YlmRec const&B, scalar const&x) {
      for(int l=0,i=0; l<L1(A); l+=2,i+=4*l-2)
	for(int m=0; m<=l; m+=2)
	  Connector::op(e(A,i+m),e(B,i+m),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(AnlRec&A, AnlRec const&B, scalar const&x) {
      AUX<PotExp::reflexion>:: template Connect<Connector>(A,B,x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(Anlm&A, Anlm const&B, scalar const&x) {
      for(int n=0,j=0; n!=N1(A); ++n,j+=L1Q(A))
	for(int l=0; l<L1(A); l+=2)
	  for(int m=0,i=j+l*l+l; m<=l; m+=2,i+=2)
	    Connector::op(e(A,i),e(B,i),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(Anlm&A, AnlRec const&P, YlmRec const&Y,
			scalar const&x) {
      scalar      *An = p(A);
      const scalar*Pn = p(P);
      for(int n=0; n!=N1(A); ++n, An+=L1Q(A), Pn+=L1(A)) {
	scalar      *Al = An;
	const scalar*Pi = Pn;
	const scalar*Yl = p(Y);
	for(int l=0; l<L1(A); l+=2,Al+=4*l-2,Pi+=2,Yl+=4*l-2) {
	  const scalar Pnl = *Pi;
	  for(int m=0; m<=l; m+=2)
	    Connector::op(Al[m], Pnl*Yl[m], x);
	}
      }
    }
    //--------------------------------------------------------------------------
    static scalar Dot(CasRec const&A, CasRec const&B) {
      register scalar x=e(A,0)*e(B,0);
      for(int m=2; m<=L(A); m+=2)
	x+= twice(e(A,m)*e(B,m));
      return x;
    }
    //--------------------------------------------------------------------------
    static scalar Dot(YlmRec const&A, YlmRec const&B) {
      register scalar x=0.;
      for(int l=0,i=0; l<L1(A); l+=2,i+=4*l-2) {
	x += e(A,i) * e(B,i);
	for(int m=2; m<=l; m+=2)
	  x += twice(e(B,i+m) * e(B,i+m));
      }
      return x;
    }
    //--------------------------------------------------------------------------
    static scalar Dot(AnlRec const&A, AnlRec const&B) {
      return AUX<PotExp::reflexion>::Dot(A,B);
    }
    //--------------------------------------------------------------------------
    static scalar Dot(Anlm const&A, Anlm const&B) {
      register scalar x=0.;
      for(int n=0,j=0; n!=N1(A); ++n,j+=L1Q(A))
	for(int l=0; l<L1(A); l+=2)
	  for(int m=0,i=j+l*l+l; m<=l; m+=2,i+=2)
	    x += m==0? e(A,i)*e(B,i) : twice(e(A,i)*e(B,i));
      return x;
    }
    //--------------------------------------------------------------------------
    static scalar Dot(Anlm const&A, AnlRec const&P, YlmRec const&Y) {
      // P = Sum_nlm A_nlm * P_nl * Y_lm                                        
      scalar       x=0;
      const scalar*An = p(A);
      const scalar*Pn = p(P);
      for(int n=0; n!=N1(A); ++n, An+=L1Q(A), Pn+=L1(A)) {
	const scalar*Al = An;
	const scalar*Pi = Pn;
	const scalar*Yl = p(Y);
	for(int l=0; l<L1(A); l+=2,Al+=4*l-2,Pi+=2,Yl+=4*l-2) {
	  register scalar y=0;
	  for(int m=0; m<=l; m+=2)
	    y += (m==0)? Al[m]*Yl[m] : twice(Al[m]*Yl[m]);
	  x += *Pi * y;
	}
      }
      return x;
    }
    //--------------------------------------------------------------------------
    template<typename T>
    static scalar Dot(tupel<3,T>&dx,
		      Anlm   const&A, AnlRec const&P, AnlRec const&R, 
		      YlmRec const&Y, YlmRec const&T, YlmRec const&Q) {
      // P = Sum_nlm A_nlm * P_nl * Y_lm                                        
      // dP/d(r,the,phi)                                                        
      scalar x=0,xr=0,xt=0,xp=0;
      const scalar*An = p(A);
      const scalar*Pn = p(P);
      const scalar*Rn = p(R);
      for(int n=0; n!=N1(A); ++n,An+=L1Q(A),Pn+=L1(A),Rn+=L1(A)) {
	const scalar*Al = An;
	const scalar*Pi = Pn;
	const scalar*Ri = Rn;
	const scalar*Yl = p(Y);
	const scalar*Tl = p(T);
	const scalar*Ql = p(Q);
	for(int l=0; l<L1(A);
	    l+=2,Al+=4*l-2,Pi+=2,Ri+=2,Yl+=4*l-2,Tl+=4*l-2,Ql+=4*l-2) {
	  register scalar y=0,yt=0,yp=0;
	  for(int m=0; m<=l; m+=2)
	    if(m==0) {
	      y  += Al[m] * Yl[m];
	      yt += Al[m] * Tl[m];
	    } else {
	      y  += twice(Al[m] * Yl[m]);
	      yt += twice(Al[m] * Tl[m]);
	      yp += twice(Al[m] * Ql[m]);
	    }
	  x  += *Pi * y;
	  xr += *Ri * y;
	  xt += *Pi * yt;
	  xp += *Pi * yp;
	}
      }
      dx[0] = xr;
      dx[1] = xt;
      dx[2] = xp;
      return x;
    }
    //--------------------------------------------------------------------------
    static void SetCas(CasRec&A, scalar const&c, scalar const&s) {
      e(A,0) = 1;
      const scalar cc=c*c-s*s, ss=2*c*s;
      for(int m=0,m2=2; m2<=L(A); m+=2,m2+=2) {
	e(A, m2) = cc * e(A, m) + ss * e(A,-m);
	e(A,-m2) = cc * e(A,-m) - ss * e(A, m);
      }
      for(int m=2; m<=L(A); m+=2)
	e(A,m) = scalar(0.5) * (e(A,m) + e(A,-m));
    }
    //--------------------------------------------------------------------------
    static void SetYlm(YlmRec&Y,
		       scalar const&ct, scalar const&st,
		       scalar const&cp, scalar const&sp) {
      AUX<PotExp::reflexion>::SetPlm(Y,ct,st);     // set Plm(theta)            
      register scalar _Cp=1.,Cp;                   // cas( m*phi)               
      register scalar _Cm=1.,Cm;                   // cas(-m*phi)               
      const scalar cc=cp*cp-sp*sp, ss=2*cp*sp;
      for(int m=2; m<L1(Y); m+=2) {                // LOOP m=2...L,2            
	Cp = cc * _Cp + ss * _Cm;                  //   cas[ m) := cas( m*phi)  
	Cm = cc * _Cm - ss * _Cp;                  //   cas[-m) := cas(-m*phi)  
	for(int l = m,                             //   LOOP l=m...L,2          
	      di  = 2*(m+m-1),                     //     4*l-2 =i(l,m)-i(l-2,m)
	      ip  = m*(m+2);                       //     i(l,m)= l*(l+1)+m     
	    l < L1(Y);                             //     l to run to L         
	    l  += 2,                               //     l       <- l+2        
	      di += 8,                             //     4*l-2   <- 4*l-2 + 8  
	      ip += di) {                          //     i(l, m) <- i(l, m)+di 
	  e(Y,ip) *= 0.5 * (Cp+Cm);                //     Y(l, m) = P(l,m)*cos  
	}                                          //   END LOOP(l)             
	_Cm = Cm;                                  //   cas[-(m-2)) <- cas[-m)  
	_Cp = Cp;                                  //   cas[  m-2 ) <- cas[ m)  
      }                                            // END LOOP(m)               
    }
    //--------------------------------------------------------------------------
    static void SetYlm(YlmRec&Y, YlmRec&T, YlmRec&P,
		       scalar const&ct, scalar const&st,
		       scalar const&cp, scalar const&sp) {
      AUX<PotExp::reflexion>::SetPlm(Y,T,ct,st);   // set Plm(the), dPlm/dthe   
      register scalar _Cp=1.,Cp;                   // cas( m*phi)               
      register scalar _Cm=1.,Cm;                   // cas(-m*phi)               
      const scalar cc=cp*cp-sp*sp, ss=2*cp*sp;
      for(int l=0; l<L1(Y); l+=2) e(P,l*(l+1))=0;  // dY(l,m=0)/dp=0            
      for(int m=2; m<L1(Y); m+=2) {                // LOOP m=2...L,2            
	Cp = cc * _Cp + ss * _Cm;                  //   cas[ m) := cas( m*phi)  
	Cm = cc * _Cm - ss * _Cp;                  //   cas[-m) := cas(-m*phi)  
	for(int l = m,                             //   LOOP l=m...L,2          
	      di  = 2*(m+m-1),                     //     4*l-2 =i(l,m)-i(l-2,m)
	      ip  = m*(m+2);                       //     i(l,m)= l*(l+1)+m     
	    l < L1(Y);                             //     l to run to L         
	    l  += 2,                               //     l       <- l+2        
	      di += 8,                             //     4*l-2   <- 4*l-2 + 8  
	      ip += di) {                          //     i(l, m) <- i(l, m)+di 
	  const scalar Cos=0.5*(Cp+Cm);            //     cos[m)                
	  e(P,ip)  = e(Y,ip)*0.5*m*(Cm-Cp);        //     dY/dphi =-P  *m*sin[m)
	  e(Y,ip) *= Cos;                          //      Y      = P    *cos[m)
	  e(T,ip) *= Cos;                          //     dY/dth  =dP/dth*cos[m)
	}                                          //   END LOOP(l)             
	_Cm = Cm;                                  //   cas[-(m-2)) <- cas[-m)  
	_Cp = Cp;                                  //   cas[  m-2 ) <- cas[ m)  
      }                                            // END LOOP(m)               
    }
    //--------------------------------------------------------------------------
    static void SetPsi (AnlRec&P, scalar const&r, scalar const&GM) {
      AUX<PotExp::reflexion>::SetPsi(P,r,GM);
    }
    //--------------------------------------------------------------------------
    static void SetPsi (AnlRec&P, AnlRec&D, scalar const&r) {
      AUX<PotExp::reflexion>::SetPsi(P,D,r);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // 4. cylindrical symmetry                                                  //
  //                                                                          //
  // n=0...L; l=0...L,2; m=0                                                  //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<> struct AUX<PotExp::cylindrical> : private PotExpAccess {
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(CasRec&A, CasRec const&B, scalar const&x) {
      Connector::op(e(A,0),e(B,0),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(YlmRec&A, YlmRec const&B, scalar const&x) {
      for(int l=0,i=0; l<L1(A); l+=2,i+=4*l-2)
	Connector::op(e(A,i),e(B,i),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(AnlRec&A, AnlRec const&B, scalar const&x) {
      AUX<PotExp::reflexion>::template Connect<Connector>(A,B,x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(Anlm&A, Anlm const&B, scalar const&x) {
      for(int n=0,j=0; n!=N1(A); ++n,j+=L1Q(A))
	for(int l=0,i=j; l<L1(A); l+=2, i+=4*l-2)
	  Connector::op(e(A,i),e(B,i),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(Anlm&A, AnlRec const&P, YlmRec const&Y,
			scalar const&x) {
      scalar      *An = p(A);
      const scalar*Pn = p(P);
      for(int n=0; n!=N1(A); ++n, An+=L1Q(A), Pn+=L1(A)) {
	scalar      *Al = An;
	const scalar*Pi = Pn;
	const scalar*Yl = p(Y);
	for(int l=0; l<L1(A); l+=2,Al+=4*l-2,Pi+=2,Yl+=4*l-2)
	  Connector::op(*Al, *Pi * *Yl, x);
      }
    }
    //--------------------------------------------------------------------------
    static scalar Dot(CasRec const&A, CasRec const&B) {
      return e(A,0) * e(B,0);
    }
    //--------------------------------------------------------------------------
    static scalar Dot(YlmRec const&A, YlmRec const&B) {
      register scalar x=0.;
      for(int l=0,i=0; l<L1(A); l+=2,i+=4*l-2)
	x += e(A,i) * e(B,i);
      return x;
    }
    //--------------------------------------------------------------------------
    static scalar Dot(AnlRec const&A, AnlRec const&B) {
      return AUX<PotExp::reflexion>::Dot(A,B);
    }
    //--------------------------------------------------------------------------
    static scalar Dot(Anlm const&A, Anlm const&B) {
      register scalar x=0.;
      for(int n=0,j=0; n!=N1(A); ++n,j+=L1Q(A))
	for(int l=0,i=j; l<L1(A); l+=2, i+=4*l-2)
	  x += e(A,i)*e(B,i);
      return x;
    }
    //--------------------------------------------------------------------------
    static scalar Dot(Anlm const&A, AnlRec const&P, YlmRec const&Y) {
      // P = Sum_nlm A_nlm * P_nl * Y_lm                                        
      scalar       x  = 0;
      const scalar*An = p(A);
      const scalar*Pn = p(P);
      for(int n=0; n!=N1(A); ++n, An+=L1Q(A), Pn+=L1(A)) {
	const scalar*Al = An;
	const scalar*Pi = Pn;
	const scalar*Yl = p(Y);
	for(int l=0; l<L1(A); l+=2,Al+=4*l-2,Pi+=2,Yl+=4*l-2)
	  x += *Al * *Pi * *Yl;
      }
      return x;
    }
    //--------------------------------------------------------------------------
    template<typename T>
    static scalar Dot(tupel<3,T>&dx,
		      Anlm   const&A, AnlRec const&P, AnlRec const&R, 
		      YlmRec const&Y, YlmRec const&T, YlmRec const&) {
      // P = Sum_nlm A_nlm * P_nl * Y_lm                                        
      // dP/d(r,the,phi)                                                        
      scalar x=0,xr=0,xt=0;
      const scalar*An = p(A);
      const scalar*Pn = p(P);
      const scalar*Rn = p(R);
      for(int n=0; n!=N1(A); ++n,An+=L1Q(A),Pn+=L1(A),Rn+=L1(A)) {
	const scalar*Al = An;
	const scalar*Pi = Pn;
	const scalar*Ri = Rn;
	const scalar*Yl = p(Y);
	const scalar*Tl = p(T);
	for(int l=0; l<L1(A); l+=2,Al+=4*l-2,Pi+=2,Ri+=2,Yl+=4*l-2,Tl+=4*l-2) {
	  x  += *Al * *Pi * *Yl;
	  xr += *Al * *Ri * *Yl;
	  xt += *Al * *Pi * *Tl;
	}
      }
      dx[0] = xr;
      dx[1] = xt;
      dx[2] = 0;
      return x;
    }
    //--------------------------------------------------------------------------
    static void SetCas(CasRec&A, scalar const&, scalar const&) {
      e(A,0) = 1;
    }
    //--------------------------------------------------------------------------
    static void SetYlm(YlmRec&Y,
		       scalar const&ct, scalar const&st,
		       scalar const&  , scalar const&) {
      e(Y,0) = 1;                                  // Y(0,0) = 1                
      e(Y,2) = ct;                                 // Y(1,0) = cos(theta)       
      for(int l = 1,                               // LOOP l+1=1...L            
	    lp  = 2,                               //   (l+1)                   
	    tlp = 3,                               //   2*l+1                   
	    i   = 2,                               //   i(l  ,0)  = l*(l+1)     
	    ip  = 6,                               //   i(l+1,0)  = (l+1)*(l+2) 
	    im  = 0;                               //   i(l-1,0)  = (l-1)*l     
	  l < L(Y);                                //   l+1 to run to L         
	  ++l,                                     //   l        <- l+1         
	    ++lp,                                  //   l+1      <- l+2         
	    tlp+= 2,                               //   2*l+1    <- 2*l+1 + 2   
	    im  = i,                               //   i(l-1,0) <- i(l  ,0)    
	    i   = ip,                              //   i(l  ,0) <- i(l+1,0)    
	    ip += lp+lp)                           //   i(l+1,0) <- i(l+2,0)    
	e(Y,ip) = (  tlp * ct * e(Y,i) - l * e(Y,im) ) / scalar(lp);
    }
    //--------------------------------------------------------------------------
    static void SetYlm(YlmRec&Y, YlmRec&T, YlmRec&P,
		       scalar const&ct, scalar const&st,
		       scalar const&  , scalar const&) {
      e(Y,0) = 1;                                  // Y(0,0) = 1                
      e(T,0) = 0;                                  // dY(0,0)/dthe = 0          
      e(P,0) = 0;                                  // dY(0,0)/dphi = 0          
      e(Y,2) = ct;                                 // Y(1,0)       = cos(the)   
      e(T,2) =-st;                                 // dY(1,0)/dthe = sin(the)   
      e(P,2) =  0;                                 // dY(1,0)/dphi = 0          
      for(int l = 1,                               // LOOP l+1=1...L            
	    lp  = 2,                               //   (l+1)                   
	    tlp = 3,                               //   2*l+1                   
	    i   = 2,                               //   i(l  ,0)  = l*(l+1)     
	    ip  = 6,                               //   i(l+1,0)  = (l+1)*(l+2) 
	    im  = 0;                               //   i(l-1,0)  = (l-1)*l     
	  l < L(Y);                                //   l+1 to run to L         
	  ++l,                                     //   l        <- l+1         
	    ++lp,                                  //   l+1      <- l+2         
	    tlp+= 2,                               //   2*l+1    <- 2*l+1 + 2   
	    im  = i,                               //   i(l-1,0) <- i(l  ,0)    
	    i   = ip,                              //   i(l  ,0) <- i(l+1,0)    
	    ip += lp+lp) {                         //   i(l+1,0) <- i(l+2,0)    
	const scalar ilp = 1/scalar(lp);
	e(Y,ip) = (tlp* ct*e(Y,i)          - l*e(Y,im)) * ilp;
	e(T,ip) = (tlp*(ct*e(T,i)-st*e(Y,i)) - l*e(T,im)) * ilp;
	e(P,ip) = 0.;
      }                                            // END LOOP(l)               
    }
    //--------------------------------------------------------------------------
    static void SetPsi (AnlRec&P, scalar const&r, scalar const&GM) {
      AUX<PotExp::reflexion>::SetPsi(P,r,GM);
    }
    //--------------------------------------------------------------------------
    static void SetPsi (AnlRec&P, AnlRec&D, scalar const&r) {
      AUX<PotExp::reflexion>::SetPsi(P,D,r);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // 5. spherical symmetry                                                    //
  //                                                                          //
  // n=0...L; l=0; m=0                                                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<> struct AUX<PotExp::spherical> : private PotExpAccess {
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(CasRec&A, CasRec const&B, scalar const&x) {
      Connector::op(e(A,0),e(B,0),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(YlmRec&A, YlmRec const&B, scalar const&x) {
      Connector::op(e(A,0),e(B,0),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(AnlRec&A, AnlRec const&B, scalar const&x) {
      for(int n=0,i=0; n!=N1(A); ++n,i+=L1(A))
	Connector::op(e(A,i),e(B,i),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(Anlm&A, Anlm const&B, scalar const&x) {
      for(int n=0,i=0; n!=N1(A); ++n,i+=L1Q(A))
	Connector::op(e(A,i),e(B,i),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(Anlm&A, AnlRec const&P, YlmRec const&Y,
			scalar const&x) {
      scalar      *An = p(A);
      const scalar*Pn = p(P);
      if(e(Y,0) == scalar(1))
	for(int n=0; n!=N1(A); ++n, An+=L1Q(A), Pn+=L1(A))
	  Connector::op(*An, *Pn, x);
      else
	for(int n=0; n!=N1(A); ++n, An+=L1Q(A), Pn+=L1(A))
	  Connector::op(*An, *Pn * e(Y,0), x);
    }
    //--------------------------------------------------------------------------
    static scalar Dot(CasRec const&A, CasRec const&B) {
      return e(A,0) * e(B,0);
    }
    //--------------------------------------------------------------------------
    static scalar Dot(YlmRec const&A, YlmRec const&B) {
      return e(A,0) * e(B,0);
    }
    //--------------------------------------------------------------------------
    static scalar Dot(AnlRec const&A, AnlRec const&B) {
      register scalar x=0;
      for(int n=0,i=0; n!=N1(A); ++n,i+=L1(A))
	x += e(A,i)*e(B,i);
      return x;
    }
    //--------------------------------------------------------------------------
    static scalar Dot(Anlm const&A, Anlm const&B) {
      register scalar x=0;
      for(int n=0,i=0; n!=N1(A); ++n,i+=L1Q(A))
	x += e(A,i)*e(B,i);
      return x;
    }
    //--------------------------------------------------------------------------
    static scalar Dot(Anlm const&A, AnlRec const&P, YlmRec const&Y) {
      // P = Sum_nlm A_nlm * P_nl * Y_lm                                        
      scalar       x  = 0;
      const scalar*An = p(A);
      const scalar*Pn = p(P);
      if(e(Y,0) == scalar(1))
	for(int n=0; n!=N1(A); ++n, An+=L1Q(A), Pn+=L1(A))
	  x += *An * *Pn;
      else
	for(int n=0; n!=N1(A); ++n, An+=L1Q(A), Pn+=L1(A))
	  x += *An * *Pn * e(Y,0);
      return x;
    }
    //--------------------------------------------------------------------------
    template<typename T>
    static scalar Dot(tupel<3,T>&dx,
		      Anlm   const&A, AnlRec const&P, AnlRec const&R, 
		      YlmRec const&Y, YlmRec const&, YlmRec const&) {
      // P = Sum_nlm A_nlm * P_nl * Y_lm                                        
      // dP/d(r,the,phi)                                                        
      scalar x=0,xr=0;
      const scalar*An = p(A);
      const scalar*Pn = p(P);
      const scalar*Rn = p(R);
      if(e(Y,0) == scalar(1))
	for(int n=0; n!=N1(A); ++n,An+=L1Q(A),Pn+=L1(A),Rn+=L1(A)) {
	  x  += *An * *Pn;
	  xr += *An * *Rn;
	}
      else
	for(int n=0; n!=N1(A); ++n,An+=L1Q(A),Pn+=L1(A),Rn+=L1(A)) {
	  x  += *An * *Pn * e(Y,0);
	  xr += *An * *Rn * e(Y,0);
	}
      dx[0] = xr;
      dx[1] = 0;
      dx[2] = 0;
      return x;
    }
    //--------------------------------------------------------------------------
    static void SetCas(CasRec&A, scalar const&, scalar const&) {
      e(A,0) = 1;
    }
    //--------------------------------------------------------------------------
    static void SetYlm(YlmRec&Y,
		       scalar const&, scalar const&,
		       scalar const&, scalar const&) {
      e(Y,0) = 1;                                  // Y(0,0) = 1                
    }
    //--------------------------------------------------------------------------
    static void SetYlm(YlmRec&Y, YlmRec&T, YlmRec&P,
		       scalar const&, scalar const&,
		       scalar const&, scalar const&) {
      e(Y,0) = 1;                                  // Y(0,0) = 1                
      e(T,0) = 0;                                  // dY(0,0)/dthe = 0          
      e(P,0) = 0;                                  // dY(0,0)/dphi = 0          
    }
    //--------------------------------------------------------------------------
    static void SetPsi (AnlRec&P, scalar const&r, scalar const&GM) {
      // routine checked (08/07/04) against MAPLE                               
      //                                                                        
      //                            l                          (1/a)            
      //                     G*M * r           (2*l+1)*a+1/2  r      - 1        
      // sets P(n,l) = ---------------------  G              (----------)       
      //                 (1/a)     (2*l+1)*a   n               (1/a)            
      //               (r       + 1)                          r      + 1        
      //                                                                        
      // where G denote the Gegenbauer polynomials.                             
      // see also Zhao (1996, MNRAS 278, 488)                                   
      register scalar xi;                          // (r^(1/a)-1)/(r^(1/a)+1)   
      register scalar fi;                          // 1/(1+r^(1/a))^a           
      SetXiFi(xi,fi,r);
      e(P,0) = GM*fi;                              // P(0,0) = G*M/(1+r^(1/a))^a
      if(N1(P)==1) return;                         // no terms n>0 --> DONE     
      const scalar tlm = twice(lambda(0));         // 2*lambda(l=0)             
      const scalar xi2 = xi+xi;                    // 2*xi                      
      e(P,L1(P)) = tlm * xi * e(P,0);              // P(1,0) = 2*lam*xi*P(0,0)  
      register scalar tlmn2xi = (tlm+2)*xi;        // 2*(lam+n)*xi              
      register scalar tlm1n   = tlm;               // 2*lam+n-1                 
      for(int n1= 2,                               // LOOP n+1=2...N            
	    im  = 0,                               //   i(n-1,l)                
	    i   = im+L1(P),                        //   i(n  ,l)                
	    ip  = i +L1(P);                        //   i(n+1,l)                
	  n1 < N1(P);                              //   n+1 runs til N          
	  ++n1,                                    //   n        <- n+1         
	    im  = i,                               //   i(n-1,l) <- i(n  ,l)    
	    i   = ip,                              //   i(n  ,l) <- i(n+1,l)    
	    ip += L1(P),                           //   i(n+1,l) += L+1         
	    tlmn2xi += xi2,                        //   2*(lam+n)*xi            
	    ++tlm1n)                               //   2*lam+n-1               
	e(P,ip) = (tlmn2xi*e(P,i)-tlm1n*e(P,im))/scalar(n1); // see GR 8.933.1  
    }
    //--------------------------------------------------------------------------
    static void SetPsi (AnlRec&P, AnlRec&D, scalar const&r) {
      // sets Psi_nl(r) and dPsi_nl/dr                                          
      register scalar xi,dxi;                      // (r^(1/a)-1)/(r^(1/a)+1)   
      register scalar fi,dfi;                      // 1/(1+r^(1/a))^a           
      SetXiFi(xi,dxi,fi,dfi,r);
      e(P,0) =  fi;                                // P(0,0) = 1/(1+r^(1/a))^a  
      e(D,0) = dfi;                                // D(0,0) = dP(0,0)/dr       
      if(N1(P)==1) return;                         // no terms n>0 --> DONE     
      const scalar tlm  = twice(lambda(0));        // 2*lambda(l=0)             
      const scalar xi2  = xi+xi;                   // 2*xi                      
      const scalar dxi2 = dxi+dxi;                 // 2*dxi/dr                  
      e(P,L1(P)) = tlm* xi*e(P,0);                 // P(1,0) = 2*lam*xi*P(0,0)  
      e(D,L1(P)) = tlm*(xi*e(D,0)+dxi*e(P,0));     // D(1,0) = dP(1,0)/dr       
      register scalar tlmn2xi  = (tlm+2)*xi;       // 2*(lam+n)*xi              
      register scalar dtlmn2xi = (tlm+2)*dxi;      // 2*(lam+n)*dxi             
      register scalar tlm1n    = tlm;              // 2*lam+n-1                 
      for(int n1= 2,                               // LOOP n+1=2...N            
	    im  = 0,                               //   i(n-1,l)                
	    i   = im+L1(P),                        //   i(n  ,l)                
	    ip  = i +L1(P);                        //   i(n+1,l)                
	  n1 < N1(P);                              //   n+1 runs til N          
	  ++n1,                                    //   n        <- n+1         
	    im  = i,                               //   i(n-1,l) <- i(n  ,l)    
	    i   = ip,                              //   i(n  ,l) <- i(n+1,l)    
	    ip += L1(P),                           //   i(n+1,l) += L+1         
	    tlmn2xi += xi2,                        //   2*(lam+n)*xi            
	    dtlmn2xi+= dxi2,                       //   2*(lam+n)*dxi           
	    ++tlm1n) {                             //   2*lam+n-1               
	const scalar in1 = 1/scalar(n1);           //   see GR 8.933.1          
	e(P,ip) = (tlmn2xi*e(P,i)              -tlm1n*e(P,im)) * in1;
	e(D,ip) = (tlmn2xi*e(D,i)+dtlmn2xi*e(P,i)-tlm1n*e(D,im)) * in1;
      }
    }
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // function templates                                                       //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<symmetry S> inline
  void SetCas (CasRec&Y, scalar const&c, scalar const&s) {
    AUX<S>::SetCas(Y,c,s);
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline
  void SetYlm (YlmRec&Y,
	       scalar const&ct, scalar const&st,
	       scalar const&cp, scalar const&sp) {
    AUX<S>::SetYlm(Y,ct,st,cp,sp);
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline
  void SetYlm (YlmRec&Y, YlmRec&Yt, YlmRec&Yp,
	       scalar const&ct, scalar const&st,
	       scalar const&cp, scalar const&sp) {
    AUX<S>::SetYlm(Y,Yt,Yp,ct,st,cp,sp);
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline
  void SetPsi (AnlRec&P, scalar const&r, scalar const&GM) {
    AUX<S>::SetPsi(P,r,GM);
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline
  void SetPsi (AnlRec&P, AnlRec&Pr, scalar const&r) {
    AUX<S>::SetPsi(P,Pr,r);
  }
  //----------------------------------------------------------------------------
  template<symmetry S, typename T> inline
  scalar EvalG(tupel<3,T>&d,
	       Anlm   const&C, AnlRec const&P,  AnlRec const&Pr,
	       YlmRec const&Y, YlmRec const&Yt, YlmRec const&Yp) {
    return R0 * AUX<S>::Dot(d,C,P,Pr,Y,Yt,Yp);
  }
  //----------------------------------------------------------------------------
  template<symmetry S, typename T> inline
  scalar EvalP(Anlm const&C, AnlRec const&P, YlmRec const&Y) {
    return R0 * AUX<S>::Dot(C,P,Y);
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // symmetry as run-time argument                                            //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  void SetCas(symmetry sym, CasRec&C, scalar const&c, scalar const&s) {
    if(sym & PotExp::zrot)
      return SetCas<PotExp::cylindrical>(C,c,s);
    if(sym & PotExp::axes)
      return SetCas<PotExp::triaxial>   (C,c,s);
    if(sym & PotExp::pint)
      return SetCas<PotExp::reflexion>  (C,c,s);
    return   SetCas<PotExp::none>       (C,c,s);
  }
  //----------------------------------------------------------------------------
  void SetYlm(symmetry sym, YlmRec&Y, 
	      scalar const&ct, scalar const&st,
	      scalar const&cp, scalar const&sp) {
    if(sym & PotExp::arot)
      return SetYlm<PotExp::spherical>  (Y,ct,st,cp,sp);
    if(sym & PotExp::zrot)
      return SetYlm<PotExp::cylindrical>(Y,ct,st,cp,sp);
    if(sym & PotExp::axes)
      return SetYlm<PotExp::triaxial>   (Y,ct,st,cp,sp);
    if(sym & PotExp::pint)
      return SetYlm<PotExp::reflexion>  (Y,ct,st,cp,sp);
    return   SetYlm<PotExp::none>       (Y,ct,st,cp,sp);
  }
  //----------------------------------------------------------------------------
  void SetYlm(symmetry sym, YlmRec&Y, YlmRec&T, YlmRec&P,
	      scalar const&ct, scalar const&st,
	      scalar const&cp, scalar const&sp) {
    if(sym & PotExp::arot)
      return SetYlm<PotExp::spherical>  (Y,T,P,ct,st,cp,sp);
    if(sym & PotExp::zrot)
      return SetYlm<PotExp::cylindrical>(Y,T,P,ct,st,cp,sp);
    if(sym & PotExp::axes)
      return SetYlm<PotExp::triaxial>   (Y,T,P,ct,st,cp,sp);
    if(sym & PotExp::pint)
      return SetYlm<PotExp::reflexion>  (Y,T,P,ct,st,cp,sp);
    return   SetYlm<PotExp::none>       (Y,T,P,ct,st,cp,sp);
  }
  //----------------------------------------------------------------------------
  void SetPsi(symmetry sym, AnlRec&P, scalar const&r, scalar const&m) {
    if(sym & PotExp::arot)
      return SetPsi<PotExp::spherical>  (P,r,m);
    if(sym & PotExp::zrot)
      return SetPsi<PotExp::cylindrical>(P,r,m);
    if(sym & PotExp::axes)
      return SetPsi<PotExp::triaxial>   (P,r,m);
    if(sym & PotExp::pint)
      return SetPsi<PotExp::reflexion>  (P,r,m);
    return   SetPsi<PotExp::none>       (P,r,m);
  }
  //----------------------------------------------------------------------------
  void SetPsi(symmetry sym, AnlRec&P, AnlRec&Pr, scalar const&r) {
    if(sym & PotExp::arot)
      return SetPsi<PotExp::spherical>  (P,Pr,r);
    if(sym & PotExp::zrot)
      return SetPsi<PotExp::cylindrical>(P,Pr,r);
    if(sym & PotExp::axes)
      return SetPsi<PotExp::triaxial>   (P,Pr,r);
    if(sym & PotExp::pint)
      return SetPsi<PotExp::reflexion>  (P,Pr,r);
    return   SetPsi<PotExp::none>       (P,Pr,r);
  }
  //////////////////////////////////////////////////////////////////////////////
  inline void SetKnl (AnlRec&K) {
    AUX<PotExp::none>::SetKnl(K);
  }
  //----------------------------------------------------------------------------
  inline void SetAnl (AnlRec&A, AnlRec const&K) {
    AUX<PotExp::none>::SetAnl(A,K);
  }
  //----------------------------------------------------------------------------
  inline void SetNlm (YlmRec&Y) {
    AUX<PotExp::none>::SetNlm(Y);
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // structs used as template parameter for AUX<>::Connect<>                  //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  struct __setX { template<typename X> static void op (X&a,X const& ,X const&x)
    { a =x; } };
  struct __mulX { template<typename X> static void op (X&a,X const& ,X const&x)
    { a*=x; } };
  struct __setB { template<typename X> static void op (X&a,X const&b,X const&)
    { a =b; } };
  struct __mulB { template<typename X> static void op (X&a,X const&b,X const&)
    { a*=b; } };
  struct __addB { template<typename X> static void op (X&a,X const&b,X const&)
    { a+=b; } };
  struct __subB { template<typename X> static void op (X&a,X const&b,X const&)
    { a-=b; } };
  struct __setT { template<typename X> static void op (X&a,X const&b,X const&x)
    { a =x*b; } };
  struct __addT { template<typename X> static void op (X&a,X const&b,X const&x)
    { a+=x*b; } };
  struct __subT { template<typename X> static void op (X&a,X const&b,X const&x)
    { a-=x*b; } };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // Spherical(): computing spherical coordinates from Cartesian              //
  // Cartesian(): computing cartesian acceleration from spherical coordinates //
  //              and (a_r, a_theta, a_phi)                                   //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename X, typename Y>
  inline void Spherical(X&rd, X&ct, X&st, X&cp, X&sp,
			tupel<3,Y> const& x) {
    register X
      A  = x[0]*x[0]+x[1]*x[1],                    // R^2= x^2 + y^2
      B  = sqrt(A);                                // R  = sqrt(x^2 + y^2)
    A   += x[2]*x[2];                              // r^2= x^2 + y^2 + z^2
    A    = sqrt(A);                                // r  = sqrt(x^2 + y^2 + z^2)
    rd   = A*IR0;                                  // x  = r/r0
    if(B) {                                        // IF R!=0
      A  = X(1)/A;                                 //   1/r
      ct = x[2]*A;                                 //   cos(the) = z/r
      st = B*A;                                    //   sin(the) = R/r
      B  = X(1)/B;                                 //   1/R
      cp = x[0]*B;                                 //   cos(phi) = x/R
      sp = x[1]*B;                                 //   sin(phi) = y/R
    } else {                                       // ELSE R==0 
      ct = x[2]>=X(0)? X(1) : X(-1);
      st = X(0);
      cp = X(1);
      sp = X(0);
    }
  }
  //----------------------------------------------------------------------------
  // given (ar,ath,aph) and (cos[the],sin[the],cos[phi],sin[phi]), we compute   
  // in place (ax,ay,az), erasing (ar,ath,aph)                                  
  template<typename X, typename Y>
  inline void Cartesian(tupel<3,Y> &a,
			X const&ct, X const&st, X const&cp, X const&sp) {
    register Y
      am = a[0] * st + a[1] * ct,                  // am = ar*st + at*ct        
      az = a[0] * ct - a[1] * st;                  // az = ar*ct - at*st        
    a[0] = am   * cp - a[2] * sp;                  // ax = am*cp - ap*sp        
    a[1] = am   * sp + a[2] * cp;                  // ay = am*sp + ap*cp        
    a[2] = az;
  }
} // namespace {
////////////////////////////////////////////////////////////////////////////////
namespace falcON { namespace P {
#ifdef falcON_SSE
  void Spherical4(fvec4&, fvec4&, fvec4&, fvec4&, fvec4&,
		  const tupel<3,double>*);
  void Spherical4(fvec4&, fvec4&, fvec4&, fvec4&, fvec4&,
		  const tupel<3,float>*);
  void Cartesian4(tupel<3,double>*,
		  fvec4 const&, fvec4 const&, fvec4 const&, fvec4 const&);
  void Cartesian4(tupel<3,float>*,
		  fvec4 const&, fvec4 const&, fvec4 const&, fvec4 const&);
  //----------------------------------------------------------------------------
#else
  template<typename T>
  void Spherical4(fvec4&rd,
		  fvec4&ct, fvec4&st,
		  fvec4&cp, fvec4&sp, const tupel<3,T> *X) {
    Spherical(rd[0],ct[0],st[0],cp[0],sp[0], X[0]);
    Spherical(rd[1],ct[1],st[1],cp[1],sp[1], X[1]);
    Spherical(rd[2],ct[2],st[2],cp[2],sp[2], X[2]);
    Spherical(rd[3],ct[3],st[3],cp[3],sp[3], X[3]);
  }
  //----------------------------------------------------------------------------
  template<typename X>
  void Cartesian4(tupel<3,X>*a,
		  fvec4 const&ct, fvec4 const&st,
		  fvec4 const&cp, fvec4 const&sp) {
    Cartesian(a[0],ct[0],st[0],cp[0],sp[0]);
    Cartesian(a[1],ct[1],st[1],cp[1],sp[1]);
    Cartesian(a[2],ct[2],st[2],cp[2],sp[2]);
    Cartesian(a[3],ct[3],st[3],cp[3],sp[3]);
  }
#endif
} }
using namespace falcON::P;
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class CasRec                                                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
template<symmetry S> inline CasRec &CasRec::set(scalar const&cp,
						scalar const&sp) {
  AUX<S>::SetCas(*this,cp,sp);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline CasRec &CasRec::assign(scalar const&x) {
  AUX<S>::template Connect<__setX>(*this,*this,x);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline CasRec &CasRec::reset() {
  return assign<S>(scalar(0));
}
//------------------------------------------------------------------------------
template<symmetry S> inline CasRec &CasRec::multiply(scalar const&x) {
  AUX<S>::template Connect<__mulX>(*this,*this,x);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline CasRec &CasRec::divide(scalar const&x) {
  return multiply<S>(scalar(1)/x);
}
//------------------------------------------------------------------------------
template<symmetry S> inline CasRec &CasRec::copy(CasRec const&B) {
  AUX<S>::template Connect<__setB>(*this,B,scalar(0));
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline CasRec &CasRec::add(CasRec const&B) {
  AUX<S>::template Connect<__addB>(*this,B,scalar(0));
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline CasRec &CasRec::addtimes(CasRec const&B,
						     scalar const&x) {
  AUX<S>::template Connect<__addT>(*this,B,x);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline CasRec &CasRec::sub(CasRec const&B) {
  AUX<S>::template Connect<__subB>(*this,B,scalar(0));
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline CasRec &CasRec::subtimes(CasRec const&B,
						     scalar const&x) {
  AUX<S>::template Connect<__subT>(*this,B,x);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline scalar CasRec::dot(CasRec const&B) const {
  return AUX<S>::Dot(*this,B);
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class YlmRec                                                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
template<symmetry S> inline
YlmRec &YlmRec::set(scalar const&ct, scalar const&st,
		    scalar const&cp, scalar const&sp) {
  AUX<S>::SetYlm(*this,ct,st,cp,sp);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline
YlmRec &YlmRec::set(YlmRec&T, YlmRec&P,
		    scalar const&ct, scalar const&st,
		    scalar const&cp, scalar const&sp) {
  AUX<S>::SetYlm(*this,T,P,ct,st,cp,sp);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline YlmRec &YlmRec::assign(scalar const&x) {
  AUX<S>::template Connect<__setX>(*this,*this,x);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline YlmRec &YlmRec::reset() {
  return assign<S>(scalar(0));
}
//------------------------------------------------------------------------------
template<symmetry S> inline YlmRec &YlmRec::multiply(scalar const&x) {
  AUX<S>::template Connect<__mulX>(*this,*this,x);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline YlmRec &YlmRec::divide(scalar const&x) {
  return multiply<S>(scalar(1)/x);
}
//------------------------------------------------------------------------------
template<symmetry S> inline YlmRec &YlmRec::copy(YlmRec const&B) {
  AUX<S>::template Connect<__setB>(*this,B,scalar(0));
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline YlmRec &YlmRec::add(YlmRec const&B) {
  AUX<S>::template Connect<__addB>(*this,B,scalar(0));
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline YlmRec &YlmRec::addtimes(YlmRec const&B,
						     scalar const&x) {
  AUX<S>::template Connect<__addT>(*this,B,x);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline YlmRec &YlmRec::sub(YlmRec const&B) {
  AUX<S>::template Connect<__subB>(*this,B,scalar(0));
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline YlmRec &YlmRec::subtimes(YlmRec const&B,
						     scalar const&x) {
  AUX<S>::template Connect<__subT>(*this,B,x);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline scalar YlmRec::dot(YlmRec const&B) const {
  return AUX<S>::Dot(*this,B);
}
//------------------------------------------------------------------------------
inline YlmRec& YlmRec::Nlm() {
  SetNlm(*this);
  return *this;
}
//------------------------------------------------------------------------------
void YlmRec::table_print(symmetry     s,
			 std::ostream&o,
			 int          p) const
{
  int w = p + 6;
  o << "# l   m   C\n"
    << "# -----------------\n";
  int dl = (s & PotExp::pint)? 2 : 1;
  int lu = (s & PotExp::arot)? 0 : L;
  for(int l=0; l<=lu; l+=dl) {
    if(l) o << "#\n";
    int ml = (s & PotExp::axes)? 0 : -l;
    int mu = (s & PotExp::zrot)? 0 :  l;
    int dm = (s & PotExp::pint)? 2 :  1;
    for(int m=ml; m<=mu; m+=dm) {
      o << ' ' << std::setw(2) << l
	<< ' ' << std::setw(3) << m << "  "<< A[l*(l+1)+m] <<'\n';
    }
  }
  o.flush();
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class AnlRec                                                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
template<symmetry S> inline AnlRec &AnlRec::assign(scalar const&x) {
  AUX<S>::template Connect<__setX>(*this,*this,x);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline AnlRec &AnlRec::reset() {
  return assign<S>(scalar(0));
}
//------------------------------------------------------------------------------
template<symmetry S> inline AnlRec &AnlRec::multiply(scalar const&x) {
  AUX<S>::template Connect<__mulX>(*this,*this,x);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline AnlRec &AnlRec::divide(scalar const&x) {
  return multiply<S>(scalar(1)/x);
}
//------------------------------------------------------------------------------
template<symmetry S> inline AnlRec &AnlRec::copy(AnlRec const&B) {
  AUX<S>::template Connect<__setB>(*this,B,scalar(0));
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline AnlRec &AnlRec::add(AnlRec const&B) {
  AUX<S>::template Connect<__addB>(*this,B,scalar(0));
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline AnlRec &AnlRec::addtimes(AnlRec const&B,
						     scalar const&x) {
  AUX<S>::template Connect<__addT>(*this,B,x);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline AnlRec &AnlRec::sub(AnlRec const&B) {
  AUX<S>::template Connect<__subB>(*this,B,scalar(0));
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline AnlRec &AnlRec::subtimes(AnlRec const&B,
						     scalar const&x) {
  AUX<S>::template Connect<__subT>(*this,B,x);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline scalar AnlRec::dot(AnlRec const&B) const {
  return AUX<S>::Dot(*this,B);
}
//------------------------------------------------------------------------------
void AnlRec::table_print(symmetry     s,
			 std::ostream&o,
			 int          p) const
{
  int w = p + 6;
  o << "# l";
  for(int n=0; n!=N1; ++n) {
    for(int i=0; i!=p; ++i) o << ' ';
    o << "C(n="<<std::setw(2)<<n<<')';
  }
  o << '\n';
  o << "# ------";
  for(int n=0; n!=N1; ++n)
    for(int i=0; i!=w+1; ++i) o << '-';
  o << "-\n";
  int dl = (s & PotExp::pint)? 2 : 1;
  int lu = (s & PotExp::arot)? 0 : L1-1;
  for(int l=0; l<=lu; l+=dl) {
    if(l) o << "#\n";
    o << ' ' << std::setw(2) << l << "  ";
    for(int n=0; n!=N1; ++n)
      o << ' ' << std::setprecision(p) << std::setw(w) << A[n*L1+l];
    o <<'\n';
  }
  o.flush();
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::PotExp::Anlm  and  related                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace {
  template<symmetry S> inline Anlm &Anlm_assign(Anlm        &A,
					      scalar const&x) {
    AUX<S>::template Connect<__setX>(A,A,x);
    return A;
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline Anlm &Anlm_reset(Anlm&A) {
    return Anlm_assign<S>(A,scalar(0));
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline Anlm &Anlm_multiply(Anlm        &A,
						  scalar const&x) {
    AUX<S>::template Connect<__mulX>(A,A,x);
    return A;
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline Anlm &Anlm_divide(Anlm        &A,
						scalar const&x) {
    return Anlm_multiply<S>(A,scalar(1)/x);
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline Anlm &Anlm_copy(Anlm      &A,
					      Anlm const&B) {
    AUX<S>::template Connect<__setB>(A,B,scalar(0));
    return A;
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline Anlm &Anlm_add(Anlm      &A,
					     Anlm const&B) {
    AUX<S>::template Connect<__addB>(A,B,scalar(0));
    return A;
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline Anlm &Anlm_subtract(Anlm      &A,
						  Anlm const&B) {
    AUX<S>::template Connect<__subB>(A,B,scalar(0));
    return A;
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline Anlm &Anlm_multiply(Anlm      &A,
						  Anlm const&B) {
    AUX<S>::template Connect<__mulB>(A,B,scalar(0));
    return A;
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline Anlm &Anlm_addtimes(Anlm        &A,
						  Anlm   const&B,
						  scalar const&x) {
    AUX<S>::template Connect<__addT>(A,B,x);
    return A;
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline Anlm &Anlm_subtimes(Anlm        &A,
						  Anlm   const&B,
						  scalar const&x) {
    AUX<S>::template Connect<__subT>(A,B,x);
    return A;
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline scalar Anlm_dot(Anlm const&A,
					      Anlm const&B) {
    return AUX<S>::Dot(A,B);
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline Anlm &Anlm_assign(Anlm        &A,
						AnlRec const&P,
						YlmRec const&Y) {
    AUX<S>::template Connect<__setB>(A,P,Y,scalar(0));
    return A;
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline Anlm &Anlm_add(Anlm        &A,
					     AnlRec const&P,
					     YlmRec const&Y) {
    AUX<S>::template Connect<__addB>(A,P,Y,scalar(0));
    return A;
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline Anlm &Anlm_subtract(Anlm        &A,
						  AnlRec const&P,
						  YlmRec const&Y) {
    AUX<S>::template Connect<__subB>(A,P,Y,scalar(0));
    return A;
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline Anlm &Anlm_addtimes(Anlm        &A,
						  AnlRec const&P,
						  YlmRec const&Y,
						  scalar const&x) {
    AUX<S>::template Connect<__addT>(A,P,Y,x);
    return A;
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline Anlm &Anlm_subtimes(Anlm        &A,
						  AnlRec const&P,
						  YlmRec const&Y,
						  scalar const&x) {
    AUX<S>::template Connect<__subT>(A,P,Y,x);
    return A;
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline Anlm &Anlm_multiply(Anlm        &A,
						  AnlRec const&P,
						  YlmRec const&Y) {
    AUX<S>::template Connect<__mulB>(A,P,Y,scalar(0));
    return A;
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline scalar Anlm_dot(Anlm   const&A,
					      AnlRec const&P,
					      YlmRec const&Y) {
    return AUX<S>::Dot(A,P,Y);
  }
  //----------------------------------------------------------------------------
  template<symmetry S, typename T>
  inline scalar Anlm_dot(Anlm const  &A,
			 tupel<3,T>  &d,
			 AnlRec const&P,
			 AnlRec const&Pr, 
			 YlmRec const&Y,
			 YlmRec const&Yt,
			 YlmRec const&Yp) {
    return AUX<S>::Dot(d,A,P,Pr,Y,Yt,Yp);
  }
} // namespace {
////////////////////////////////////////////////////////////////////////////////
Anlm &Anlm::reset(symmetry S) {
  switch(S) {
  case spherical:   return Anlm_reset<spherical>  (*this);
  case cylindrical: return Anlm_reset<cylindrical>(*this);
  case triaxial:    return Anlm_reset<triaxial>   (*this);
  case reflexion:   return Anlm_reset<reflexion>  (*this);
  default:          return Anlm_reset<none>       (*this);
  }
}
//------------------------------------------------------------------------------
Anlm&Anlm::assign  (scalar const&x, symmetry S) {
  switch(S) {
  case spherical:   return Anlm_assign<spherical>  (*this,x);
  case cylindrical: return Anlm_assign<cylindrical>(*this,x);
  case triaxial:    return Anlm_assign<triaxial>   (*this,x);
  case reflexion:   return Anlm_assign<reflexion>  (*this,x);
  default:          return Anlm_assign<none>       (*this,x);
  }
}
//------------------------------------------------------------------------------
Anlm&Anlm::multiply  (scalar const&x, symmetry S) {
  switch(S) {
  case spherical:   return Anlm_multiply<spherical>  (*this,x);
  case cylindrical: return Anlm_multiply<cylindrical>(*this,x);
  case triaxial:    return Anlm_multiply<triaxial>   (*this,x);
  case reflexion:   return Anlm_multiply<reflexion>  (*this,x);
  default:          return Anlm_multiply<none>       (*this,x);
  }
}
//------------------------------------------------------------------------------
Anlm&Anlm::divide  (scalar const&x, symmetry S) {
  switch(S) {
  case spherical:   return Anlm_divide<spherical>  (*this,x);
  case cylindrical: return Anlm_divide<cylindrical>(*this,x);
  case triaxial:    return Anlm_divide<triaxial>   (*this,x);
  case reflexion:   return Anlm_divide<reflexion>  (*this,x);
  default:          return Anlm_divide<none>       (*this,x);
  }
}
//------------------------------------------------------------------------------
Anlm&Anlm::copy    (Anlm   const&A, symmetry S) {
  switch(S) {
  case spherical:   return Anlm_copy<spherical>  (*this,A);
  case cylindrical: return Anlm_copy<cylindrical>(*this,A);
  case triaxial:    return Anlm_copy<triaxial>   (*this,A);
  case reflexion:   return Anlm_copy<reflexion>  (*this,A);
  default:          return Anlm_copy<none>       (*this,A);
  }
}
//------------------------------------------------------------------------------
Anlm&Anlm::add    (Anlm   const&A, symmetry S) {
  switch(S) {
  case spherical:   return Anlm_add<spherical>  (*this,A);
  case cylindrical: return Anlm_add<cylindrical>(*this,A);
  case triaxial:    return Anlm_add<triaxial>   (*this,A);
  case reflexion:   return Anlm_add<reflexion>  (*this,A);
  default:          return Anlm_add<none>       (*this,A);
  }
}
//------------------------------------------------------------------------------
Anlm&Anlm::subtract(Anlm   const&A, symmetry S) {
  switch(S) {
  case spherical:   return Anlm_subtract<spherical>  (*this,A);
  case cylindrical: return Anlm_subtract<cylindrical>(*this,A);
  case triaxial:    return Anlm_subtract<triaxial>   (*this,A);
  case reflexion:   return Anlm_subtract<reflexion>  (*this,A);
  default:          return Anlm_subtract<none>       (*this,A);
  }
}
//------------------------------------------------------------------------------
Anlm&Anlm::multiply(Anlm   const&A, symmetry S) {
  switch(S) {
  case spherical:   return Anlm_multiply<spherical>  (*this,A);
  case cylindrical: return Anlm_multiply<cylindrical>(*this,A);
  case triaxial:    return Anlm_multiply<triaxial>   (*this,A);
  case reflexion:   return Anlm_multiply<reflexion>  (*this,A);
  default:          return Anlm_multiply<none>       (*this,A);
  }
}
//------------------------------------------------------------------------------
scalar Anlm::dot    (Anlm   const&A, symmetry S) const {
  switch(S) {
  case spherical:   return Anlm_dot<spherical>  (*this,A);
  case cylindrical: return Anlm_dot<cylindrical>(*this,A);
  case triaxial:    return Anlm_dot<triaxial>   (*this,A);
  case reflexion:   return Anlm_dot<reflexion>  (*this,A);
  default:          return Anlm_dot<none>       (*this,A);
  }
}
//------------------------------------------------------------------------------
Anlm&Anlm::addtimes    (Anlm   const&A, scalar const&x, symmetry S) {
  switch(S) {
  case spherical:   return Anlm_addtimes<spherical>  (*this,A,x);
  case cylindrical: return Anlm_addtimes<cylindrical>(*this,A,x);
  case triaxial:    return Anlm_addtimes<triaxial>   (*this,A,x);
  case reflexion:   return Anlm_addtimes<reflexion>  (*this,A,x);
  default:          return Anlm_addtimes<none>       (*this,A,x);
  }
}
//------------------------------------------------------------------------------
Anlm&Anlm::subtimes    (Anlm   const&A, scalar const&x, symmetry S) {
  switch(S) {
  case spherical:   return Anlm_subtimes<spherical>  (*this,A,x);
  case cylindrical: return Anlm_subtimes<cylindrical>(*this,A,x);
  case triaxial:    return Anlm_subtimes<triaxial>   (*this,A,x);
  case reflexion:   return Anlm_subtimes<reflexion>  (*this,A,x);
  default:          return Anlm_subtimes<none>       (*this,A,x);
  }
}
//------------------------------------------------------------------------------
void Anlm::print(symmetry s, std::ostream&o, int p) const {
  int w = p + 6;
  for(int n=0; n!=N1; ++n) {
    o << '\n';
    int dl = (s & PotExp::pint)? 2 : 1;
    for(int l=0; l!=L1; l+=dl) {
      int ml = (s & PotExp::axes)? 0 : -l;
      int mu = (s & PotExp::zrot)? 0 :  l;
      int dm = (s & PotExp::pint)? 2 :  1;
      o << '\n';
      for(int m=ml; m<=mu; m+=dm)
	o << " C(" << std::setw(1) << n
	  << ','   << std::setw(1) << l
	  << ','   << std::setw(2) << m
	  << ") =" << std::setprecision(p) << std::setw(w) << A[n*L1Q+l*(l+1)+m]
	  << '\n';
    }
  }
  o.flush();
}
//------------------------------------------------------------------------------
void Anlm::table_print(symmetry     s,
		       std::ostream&o,
		       int          p) const
{
  int w = p + 6;
  o << "# l   m";
  for(int n=0; n!=N1; ++n) {
    for(int i=0; i!=p; ++i) o << ' ';
    o << "C(n="<<std::setw(2)<<n<<')';
  }
  o << '\n';
  o << "# ------";
  for(int n=0; n!=N1; ++n)
    for(int i=0; i!=w+1; ++i) o << '-';
  o << "-\n";
  int dl = (s & PotExp::pint)? 2 : 1;
  int lu = (s & PotExp::arot)? 0 : L;
  for(int l=0; l<=lu; l+=dl) {
    if(l) o << "#\n";
    int ml = (s & PotExp::axes)? 0 : -l;
    int mu = (s & PotExp::zrot)? 0 :  l;
    int dm = (s & PotExp::pint)? 2 :  1;
    for(int m=ml; m<=mu; m+=dm) {
      o << ' ' << std::setw(2) << l
	<< ' ' << std::setw(3) << m << "  ";
      for(int n=0; n!=N1; ++n)
	o << ' ' << std::setprecision(p) << std::setw(w) << A[n*L1Q+l*(l+1)+m];
      o <<'\n';
    }
  }
  o.flush();
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::PotExp                                                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // CBlock used for computing the coefficients                               //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename T>
  class CBlock : private PotExpAccess {
    typedef tupel<3,T> V;
    //--------------------------------------------------------------------------
    T            M[4];                             // mass                      
    V            X[4];                             // (x,y,z)                   
    fvec4        rd,ct,st,cp,sp;                   // (r/r0,cos/sin(the/phi))   
    int          K;                                // counter                   
    Anlm        &C;                                // coeffs to add to          
    AnlRec       Psi;                              // for Psi_nl(r)             
    YlmRec       Ylm;                              // for Y_lm(the,phi)         
    //--------------------------------------------------------------------------
    void load(T const&m, V const&x) {              // load a body into buffer   
      M[K] = m;                                    //   remember its mass and   
      X[K] = x;                                    //   its position            
      ++K;                                         //   increment counter       
    }
    //--------------------------------------------------------------------------
    template<symmetry SYM>                         // symmetry at compile time  
    void flush() {                                 // flush the buffer          
      Spherical4(rd,ct,st,cp,sp,X);                //   compute spherical coords
      for(int k=0; k!=K; ++k) {                    //   LOOP buffer             
	SetPsi<SYM>(Psi,rd[k],M[k]);               //     set Psi_nl(r_i)       
	SetYlm<SYM>(Ylm,ct[k],st[k],cp[k],sp[k]);  //     set Y_lm(the_i,phi_i) 
	Anlm_add<SYM>(C,Psi,Ylm);                  //     add to C_nlm          
      }                                            //   END LOOP                
      K = 0;                                       //   reset the counter       
    }
    //--------------------------------------------------------------------------
    template<symmetry SYM>                         // symmetry at compile time  
    void AddCoeffs(                                // add to C due to bodies    
		   int       n,                    // I: # bodies               
		   const   T*m,                    // I: body masses            
		   const   V*x,                    // I: body positions         
		   const int*f,                    // I: body flags             
		   int       mark) {               // I: source mark            
      K = 0;                                       //   reset buffer            
      if(mark && f) {                              //   IF only marked bodies   
	for(int i=0; i!=n; ++i)                    //     LOOP bodies           
	  if(f[i]&mark && m[i] != 0) {             //       IF marked & massive 
	    load(m[i],x[i]);                       //         load into buffer  
	    if(K==4) flush<SYM>();                 //         flush full buffer 
	  }                                        //     END LOOP              
      } else {                                     //   ELSE (all bodies)       
	for(int i=0; i!=n; ++i)                    //     LOOP bodies           
	  if(m[i] != 0) {                          //       IF massive          
	    load(m[i],x[i]);                       //         load into buffer  
	    if(K==4) flush<SYM>();                 //         flush full buffer 
	  }                                        //     END LOOP              
      }                                            //   ENDIF                   
      if(K) flush<SYM>();                          //   flush non-empty buffer  
    }
  public:
    //--------------------------------------------------------------------------
    CBlock(Anlm&c) :                               // constructor               
      K   (0),                                     //   reset counter           
      C   (c),                                     //   get coefficients        
      Psi (C.nmax(),C.lmax()),                     //   init AnlRec for Psi_nl  
      Ylm (C.lmax())                               //   init YlmRec for Y_lm    
    {}
    //--------------------------------------------------------------------------
    void AddCoeffs(                                // add to C due to bodies    
		   symmetry,                       // I: symmetry at run-time   
		   int,                            // I: # bodies               
		   const   T*,                     // I: body masses            
		   const   V*,                     // I: body positions         
		   const int*,                     // I: body flags             
		   int);                           // I: source indicator       
  }; // class CBlock<T>                                                         
  //----------------------------------------------------------------------------
  template<typename T>
  void CBlock<T>::AddCoeffs(symmetry  s,
			    int       n,
			    const   T*m,
			    const   V*x,
			    const int*f,
			    int       mark) {
    if     (s & PotExp::arot) AddCoeffs<PotExp::spherical  >(n,m,x,f,mark);
    else if(s & PotExp::zrot) AddCoeffs<PotExp::cylindrical>(n,m,x,f,mark);
    else if(s & PotExp::axes) AddCoeffs<PotExp::triaxial   >(n,m,x,f,mark);
    else if(s & PotExp::pint) AddCoeffs<PotExp::reflexion  >(n,m,x,f,mark);
    else                      AddCoeffs<PotExp::none       >(n,m,x,f,mark);
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // GBlock used for computing gravity from the coefficients                  //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename T>
  class GBlock {
    typedef tupel<3,T> V;
    //--------------------------------------------------------------------------
    int          I[4];                             // body index                
    T            P[4];                             // potential                 
    V            X[4];                             // (x,y,z)                   
    V            A[4];                             // (ar,ath,aph) / (ax,ay,az) 
    fvec4        rd,ct,st,cp,sp;                   // (r/r0,cos/sin(the/phi))   
    int          K;                                // counter                   
    Anlm   const&C;                                // coeffs to use             
    AnlRec       Psi,dPr;                          // for Psi_nl(r)             
    YlmRec       Ylm,dYt,dYp;                      // for Y_lm(the,phi)         
    //--------------------------------------------------------------------------
    void load(int i, V const&x) {
      I[K] = i;
      X[K] = x;
      ++K;
    }
    //--------------------------------------------------------------------------
    template<symmetry SYM>                         // symmetry at compile time  
    void flush(                                    // flush buffer              
	       T  *p,                              // O: potentials             
	       V  *a,                              // O: accelerations          
	       int add) {                          // I: add or assign?         
      Spherical4(rd,ct,st,cp,sp,X);                //   set spherical coords    
      for(int k=0; k!=K; ++k) {                    //   LOOP buffer             
	SetPsi<SYM>(Psi,dPr,rd[k]);                //     set Psi, dPsi/dr      
	SetYlm<SYM>(Ylm,dYt,dYp,
		    ct[k],st[k],cp[k],sp[k]);      //     set Ylm, dYlm/d(th,ph)
	P[k]=EvalG<SYM>(A[k],C,Psi,dPr,Ylm,dYt,dYp); //   set Pot, (ar,ath,aph) 
      }                                            //   END LOOP                
      Cartesian4(A,ct,st,cp,sp);                   //   set (ax,ay,az)          
      if(add&1) for(int k=0; k!=K; ++k) p[I[k]]-= P[k]; // add potential        
      else      for(int k=0; k!=K; ++k) p[I[k]] =-P[k]; // OR assign it         
      if(add&2) for(int k=0; k!=K; ++k) a[I[k]]+= A[k]; // add acceleration     
      else      for(int k=0; k!=K; ++k) a[I[k]] = A[k]; // OR assign it         
      K = 0;                                       // reset counter             
    }
    //--------------------------------------------------------------------------
    template<symmetry SYM>                         // symmetry at compile time  
    void AddGravity(                               // add gravity to single body
		    V const&x,                     // I: positions              
		    T      &p,                     // O: potential              
		    V      &a,                     // O: acceleration           
		    int     add) {                 // I: add or assign?         
      T Rd,Ct,St,Cp,Sp;                            //   spherical coords        
      Spherical(Rd,Ct,St,Cp,Sp,x);                 //   set spherical coords    
      SetPsi<SYM>(Psi,dPr,Rd);                     //   set Psi, dPsi/dr        
      SetYlm<SYM>(Ylm,dYt,dYp,Ct,St,Cp,Sp);        //   set Ylm, dYlm/d(th,ph)  
      V Ac;
      if(add&1)
	p-= EvalG<SYM>(Ac,C,Psi,dPr,Ylm,dYt,dYp);  //   set P, dP/d(r,the,phi)  
      else
	p =-EvalG<SYM>(Ac,C,Psi,dPr,Ylm,dYt,dYp);
      Cartesian(Ac,Ct,St,Cp,Sp);                   //   set (ax,ay,az)          
      if(add&2)
	a+= Ac;
      else
	a = Ac;
    }
    //--------------------------------------------------------------------------
    template<symmetry SYM>                         // symmetry at compile time  
    void AddGravity(                               // add gravity to bodies     
		    int       n,                   // I: # bodies               
		    const   V*x,                   // I: positions              
		    T        *p,                   // O: potentials             
		    V        *a,                   // O: accelerations          
		    const int*f,                   // I: body flags             
		    int       add) {               // I: add or assign?         
      K = 0;                                       //   reset buffer            
      if(f) {                                      //   IF only flagged bodies  
	for(int i=0; i!=n; ++i) if(f[i] & 1) {     //     LOOP flagged bodies   
	  load(i,x[i]);                            //       load into buffer    
	  if(K==4) flush<SYM>(p,a,add);            //       flush full buffer   
	}                                          //     END LOOP              
      } else {                                     //   ELSE: all bodies        
	for(int i=0; i!=n; ++i) {                  //     LOOP all bodies       
	  load(i,x[i]);                            //       load into buffer    
	  if(K==4) flush<SYM>(p,a,add);            //       flush full buffer   
	}                                          //     END LOOP              
      }                                            //   ENDIF                   
      if(K) flush<SYM>(p,a,add);                   //   flush non-empty buffer  
    }
  public:
    //--------------------------------------------------------------------------
    GBlock(Anlm const&c) :                         // constructor               
      K   (0),                                     //   reset counter           
      C   (c),                                     //   get coefficients        
      Psi (C.nmax(),C.lmax()),                     //   init AnlRec for Psi_nl  
      dPr (C.nmax(),C.lmax()),                     //   init AnlRec for dPsi/dr 
      Ylm (C.lmax()),                              //   init YlmRec for Y_lm    
      dYt (C.lmax()),                              //   init YlmRec for dY/dthe 
      dYp (C.lmax())                               //   init YlmRec for dY/dphi 
    {}
    //--------------------------------------------------------------------------
    void AddGravity(                               // add gravity set of bodies 
		    symmetry,                      // I: symmetry at run time   
		    int,                           // I: # bodies               
		    const   V*,                    // I: positions              
		    T        *,                    // O: potentials             
		    V        *,                    // O: accelerations          
		    const int*,                    // I: body flags             
		    int);                          // I: add or assign?         
    //--------------------------------------------------------------------------
    void AddGravity(                               // add gravity to single body
		    symmetry,                      // I: symmetry at run time   
		    const   V&,                    // I: position               
		    T        &,                    // O: potential              
		    V        &,                    // O: acceleration           
		    int);                          // I: add or assign?         
  }; // class GBlock                                                            
  //----------------------------------------------------------------------------
  template<typename T>
  void GBlock<T>::AddGravity(
			     symmetry  s,          // I: symmetry at run time   
			     int       n,          // I: # bodies               
			     const   V*x,          // I: positions              
			     T        *p,          // O: potentials             
			     V        *a,          // O: accelerations          
			     const int*f,          // I: body flags             
			     int       add) {      // I: add or assign?         
    if     (s & PotExp::arot) AddGravity<PotExp::spherical  >(n,x,p,a,f,add);
    else if(s & PotExp::zrot) AddGravity<PotExp::cylindrical>(n,x,p,a,f,add);
    else if(s & PotExp::axes) AddGravity<PotExp::triaxial   >(n,x,p,a,f,add);
    else if(s & PotExp::pint) AddGravity<PotExp::reflexion  >(n,x,p,a,f,add);
    else                      AddGravity<PotExp::none       >(n,x,p,a,f,add);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  void GBlock<T>::AddGravity(
			     symmetry  s,          // I: symmetry at run time   
			     const   V&x,          // I: positions              
			     T        &p,          // O: potentials             
			     V        &a,          // O: accelerations          
			     int       add) {      // I: add or assign?         
    if     (s & PotExp::arot) AddGravity<PotExp::spherical  >(x,p,a,add);
    else if(s & PotExp::zrot) AddGravity<PotExp::cylindrical>(x,p,a,add);
    else if(s & PotExp::axes) AddGravity<PotExp::triaxial   >(x,p,a,add);
    else if(s & PotExp::pint) AddGravity<PotExp::reflexion  >(x,p,a,add);
    else                      AddGravity<PotExp::none       >(x,p,a,add);
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // PBlock used for computing potentials from the coefficients               //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename T>
  class PBlock {
    typedef tupel<3,T> V;
    //--------------------------------------------------------------------------
    int          I[4];                             // body index                
    T            P[4];                             // potential                 
    V            X[4];                             // (x,y,z)                   
    fvec4        rd,ct,st,cp,sp;                   // (r/r0,cos/sin(the/phi))   
    int          K;                                // counter                   
    Anlm   const&C;                                // coeffs to use             
    AnlRec       Psi;                              // for Psi_nl(r)             
    YlmRec       Ylm;                              // for Y_lm(the,phi)         
    //--------------------------------------------------------------------------
    void load(int i, V const&x) {
      I[K] = i;
      X[K] = x;
      ++K;
    }
    //--------------------------------------------------------------------------
    template<symmetry SYM>                         // symmetry at compile time  
    void flush(                                    // flush buffer              
	       T  *p,                              // O: potentials             
	       int add) {                          // I: add or assign?         
      Spherical4(rd,ct,st,cp,sp,X);                //   set spherical coords    
      for(int k=0; k!=K; ++k) {                    //   LOOP buffer             
	SetPsi<SYM>(Psi,rd[k],1);                  //     set Psi               
	SetYlm<SYM>(Ylm,ct[k],st[k],cp[k],sp[k]);  //     set Ylm               
	P[k]=EvalP<SYM,T>(C,Psi,Ylm);              //     set Pot               
      }                                            //   END LOOP                
      if(add&1) for(int k=0; k!=K; ++k) p[I[k]]-= P[k]; // add potential        
      else      for(int k=0; k!=K; ++k) p[I[k]] =-P[k]; // OR assign it         
      K = 0;                                       // reset counter             
    }
    //--------------------------------------------------------------------------
    template<symmetry SYM>                         // symmetry at compile time  
    T Potential(                                   // r: Potential of body      
		V const&x) {                       // I: positions              
      T Rd,Ct,St,Cp,Sp;                            //   spherical coords        
      Spherical(Rd,Ct,St,Cp,Sp,x);                 //   set spherical coords    
      SetPsi<SYM>(Psi,Rd,1);                       //   set Psi                 
      SetYlm<SYM>(Ylm,Ct,St,Cp,Sp);                //   set Ylm                 
      return EvalP<SYM,T>(C,Psi,Ylm);              //   compute P               
    }
    //--------------------------------------------------------------------------
    template<symmetry SYM>                         // symmetry at compile time  
    void AddPotential(                             // add gravity to bodies     
		      int       n,                 // I: # bodies               
		      const   V*x,                 // I: positions              
		      T        *p,                 // O: potentials             
		      const int*f,                 // I: body flags             
		      int       add) {             // I: add or assign?         
      K = 0;                                       //   reset buffer            
      if(f) {                                      //   IF only flagged bodies  
	for(int i=0; i!=n; ++i) if(f[i] & 1) {     //     LOOP flagged bodies   
	  load(i,x[i]);                            //       load into buffer    
	  if(K==4) flush<SYM>(p,add);              //       flush full buffer   
	}                                          //     END LOOP              
      } else {                                     //   ELSE: all bodies        
	for(int i=0; i!=n; ++i) {                  //     LOOP all bodies       
	  load(i,x[i]);                            //       load into buffer    
	  if(K==4) flush<SYM>(p,add);              //       flush full buffer   
	}                                          //     END LOOP              
      }                                            //   ENDIF                   
      if(K) flush<SYM>(p,add);                     //   flush non-empty buffer  
    }
  public:
    //--------------------------------------------------------------------------
    PBlock(Anlm const&c) :                         // constructor               
      K   (0),                                     //   reset counter           
      C   (c),                                     //   get coefficients        
      Psi (C.nmax(),C.lmax()),                     //   init AnlRec for Psi_nl  
      Ylm (C.lmax())                               //   init YlmRec for Y_lm    
    {}
    //--------------------------------------------------------------------------
    void AddPotential(                             // add gravity set of bodies 
		      symmetry,                    // I: symmetry at run time   
		      int,                         // I: # bodies               
		      const   V*,                  // I: positions              
		      T        *,                  // O: potentials             
		      const int*,                  // I: body flags             
		      int);                        // I: add or assign?         
    //--------------------------------------------------------------------------
    T Potential(                                   // R: potential: single body 
		symmetry,                          // I: symmetry at run time   
		const   V&);                       // I: position               
  }; // class PBlock                                                            
  //----------------------------------------------------------------------------
  template<typename T>
  void PBlock<T>::AddPotential(
			       symmetry  s,        // I: symmetry at run time   
			       int       n,        // I: # bodies               
			       const   V*x,        // I: positions              
			       T        *p,        // O: potentials             
			       const int*f,        // I: body flags             
			       int       add) {    // I: add or assign?         
    if     (s & PotExp::arot) AddPotential<PotExp::spherical  >(n,x,p,f,add);
    else if(s & PotExp::zrot) AddPotential<PotExp::cylindrical>(n,x,p,f,add);
    else if(s & PotExp::axes) AddPotential<PotExp::triaxial   >(n,x,p,f,add);
    else if(s & PotExp::pint) AddPotential<PotExp::reflexion  >(n,x,p,f,add);
    else                      AddPotential<PotExp::none       >(n,x,p,f,add);
  }
  //----------------------------------------------------------------------------
  template<typename T>
  T PBlock<T>::Potential(                          // R: potential: single body 
			 symmetry  s,              // I: symmetry at run time   
			 const   V&x) {            // I: positions              
    if     (s & PotExp::arot) Potential<PotExp::spherical  >(x);
    else if(s & PotExp::zrot) Potential<PotExp::cylindrical>(x);
    else if(s & PotExp::axes) Potential<PotExp::triaxial   >(x);
    else if(s & PotExp::pint) Potential<PotExp::reflexion  >(x);
    else                      Potential<PotExp::none       >(x);
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //  SelfGrav used for computing Self-Gravity                                //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename T>
  class SelfGrav {
    typedef tupel<3,T> V;
    //--------------------------------------------------------------------------
    int          N,N4;                             // # bodies, size of arrays  
    fvec4       *RD,*CT,*ST,*CP,*SP;               // arrays for spherical coord
    Anlm        &C;                                // coeffs to use             
    AnlRec       Psi,dPr;                          // for Psi_nl(r)             
    YlmRec       Ylm,dYt,dYp;                      // for Y_lm(the,phi)         
    //--------------------------------------------------------------------------
  public:
    SelfGrav(Anlm&c, int n, const V*x) :
      N   (n),
      N4  ((n+3)/4),
      RD  (new fvec4[N4]), // MUST NOT USE falcON_NEW            
      CT  (new fvec4[N4]), // because this is fvec4::new[], which
      ST  (new fvec4[N4]), // calls falcON::malloc16, which      
      CP  (new fvec4[N4]), // calls falcON_DEL() anyway!         
      SP  (new fvec4[N4]), //                                    
      C   (c),                                     //   get coefficients        
      Psi (C.nmax(),C.lmax()),                     //   init AnlRec for Psi_nl  
      dPr (C.nmax(),C.lmax()),                     //   init AnlRec for dPsi/dr 
      Ylm (C.lmax()),                              //   init YlmRec for Y_lm    
      dYt (C.lmax()),                              //   init YlmRec for dY/dthe 
      dYp (C.lmax())                               //   init YlmRec for dY/dphi 
    {
      // compute spherical coordinates for ALL bodies                           
      V X[4];
      for(int i=0,i4=0; i4!=N4; ++i4) {
	for(int k=0; i!=N && k!=4; ++i,++k)
	  X[k] = x[i];
	Spherical4(RD[i4],CT[i4],ST[i4],CP[i4],SP[i4],X);
      }
    }
    //--------------------------------------------------------------------------
    ~SelfGrav() {
      delete[] RD;  // MUST NOT USE falcON_DEL_A        
      delete[] CT;  // because this is fvec4::delete[], 
      delete[] ST;  // which calls falcON::free16,      
      delete[] CP;  // which calls falcON_DEL_A anyway! 
      delete[] SP;  //                                  
    }
    //--------------------------------------------------------------------------
    // add up C_nlm from all bodies or all bodies marked                        
    template<symmetry SYM>
    void AddCoef(const T  *m,
		 const int*f,
		 int       k) {
      if(f && k) {
	for(int i=0,i4=0; i4!=N4; ++i4)
	  for(int k=0; i!=N && k!=4; ++i,++k) if(f[i] & k && m[i]) {
	    SetPsi<SYM>(Psi,RD[i4][k],m[i]);
	    SetYlm<SYM>(Ylm,CT[i4][k],ST[i4][k],CP[i4][k],SP[i4][k]);
	    Anlm_add<SYM>(C,Psi,Ylm);
	  }
      } else {
	for(int i=0,i4=0; i4!=N4; ++i4)
	  for(int k=0; i!=N && k!=4; ++i,++k) if(m[i]) {
	    SetPsi<SYM>(Psi,RD[i4][k],m[i]);
	    SetYlm<SYM>(Ylm,CT[i4][k],ST[i4][k],CP[i4][k],SP[i4][k]);
	    Anlm_add<SYM>(C,Psi,Ylm);
	  }
      }
    }
    //--------------------------------------------------------------------------
    // set gravity due to coefficients for all or all active bodies             
    template<symmetry SYM>
    void SetGrav(T        *p,
		 V        *a,
		 const int*f,
		 int       add) {
      int I[4];
      T   P[4];
      V   A[4];
      if(f) {
	for(int i=0,i4=0,K; i4!=N4; ++i4) {
	  for(int k=0; i!=N && k!=4; ++i,++k) if(f[i] & 1) {
	    I[K] = i;
	    SetPsi<SYM>(Psi,dPr,RD[i4][K]);
	    SetYlm<SYM>(Ylm,dYt,dYp,CT[i4][K],ST[i4][K],CP[i4][K],SP[i4][K]);
	    P[K] = EvalG<SYM>(A[K],C,Psi,dPr,Ylm,dYt,dYp);
	  }
	  Cartesian4(A,CT[i4],ST[i4],CP[i4],SP[i4]);
	  if(add&1) for(int k=0; k!=K; ++k) p[I[k]]-= P[k];
	  else      for(int k=0; k!=K; ++k) p[I[k]] =-P[k];
	  if(add&2) for(int k=0; k!=K; ++k) a[I[k]]+= A[k];
	  else      for(int k=0; k!=K; ++k) a[I[k]] = A[k];
	}
      } else {
	for(int i=0,i4=0,K; i4!=N4; ++i4) {
	  for(K=0; i!=N && K!=4; ++i,++K) {
	    I[K] = i;
	    SetPsi<SYM>(Psi,dPr,RD[i4][K]);
	    SetYlm<SYM>(Ylm,dYt,dYp,CT[i4][K],ST[i4][K],CP[i4][K],SP[i4][K]);
	    P[K] = EvalG<SYM>(A[K],C,Psi,dPr,Ylm,dYt,dYp);
	  }
	  Cartesian4(A,CT[i4],ST[i4],CP[i4],SP[i4]);
	  if(add&1) for(int k=0; k!=K; ++k) p[I[k]]-= P[k];
	  else      for(int k=0; k!=K; ++k) p[I[k]] =-P[k];
	  if(add&2) for(int k=0; k!=K; ++k) a[I[k]]+= A[k];
	  else      for(int k=0; k!=K; ++k) a[I[k]] = A[k];
	}
      }
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  template<symmetry SYM>
  inline void normalize(Anlm&C, Anlm const&K, scalar const&G) {
    Anlm_multiply<SYM>(C,K);
    if(G != 1)
      Anlm_multiply<SYM>(C,G);
  }
} // namespace {
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class PotExp                                                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
PotExp::PotExp(scalar   a,                         // parameter alpha           
	       scalar   r,                         // scale radius              
	       int      n,                         // maximum n in expansion    
	       int      l,                         // maximum l in expansion    
	       symmetry s) :                       // symmetry used             
  SYM  (l==0? spherical : s),
  N    (n),
  L    (l),
  AL   (a),
  R0   (r),
  Knlm (n,l),
  STATE(0)
{
  if(AL<0.5) {
    snprintf(ERR,256,"PotExp: alpha<0.5\n");
    STATE |= 2;
    return;
  }
  if(N<0) {
    snprintf(ERR,256,"PotExp: nmax must be >= 0\n");
    STATE |= 2;
    return;
  }
  if(L&1) {
    snprintf(ERR,256,"PotExp: lmax must be even\n");
    STATE |= 2;
    return;
  }
  if(L==0 && s != spherical) {
    snprintf(WARN,256,"PotExp: assuming spherical symmetry since lmax=0\n");
    STATE |= 1;
  }
  setAL (AL);
  setR0 (R0);
  YlmRec Nlm(l);
  AnlRec Knl(n,l), Anl(n,l);
  SetNlm(Nlm);
  SetKnl(Knl);
  SetAnl(Anl,Knl);
  Anlm_assign<none>(Knlm,Anl,Nlm);
}
//------------------------------------------------------------------------------
#define CHECKMISMATCH(FUNC)						\
if(N != C.nmax() || L != C.lmax()) {					\
  if(N  !=C.nmax())							\
    if(L!=C.lmax())							\
      snprintf(ERR,256,"PotExp::%s(): max n&l mismatch\n",FUNC);	\
    else								\
      snprintf(ERR,256,"PotExp::%s(): max n mismatch\n",FUNC);		\
  else if(L!=C.lmax())							\
    snprintf(ERR,256,"PotExp::%s(): max l mismatch\n",FUNC);		\
  STATE |= 2;								\
  return;								\
}
//------------------------------------------------------------------------------
template<typename T>                               // T: double or float        
void PotExp::AddCoeffs  (Anlm            &C,       // O: C_nlm coefficients     
			 int              n,       // I: # bodies               
			 const T         *m,       // I: masses                 
			 const tupel<3,T>*x,       // I: positions              
			 const int       *f,       // I: flags        see Note  
			 int              k) const //[I: k]           see Note  
  //                                                                            
  // computes                                                                   
  //                                                                            
  //          1     N                                                           
  // C_nlm = ----  Sum m_i * Psi_nl(r_i) * Y_lm(theta_i,phi_i)                  
  //         r0^2  i=1                                                          
  //                                                                            
  // where the sum includes either all (mass-carrying) bodies if k==0, or       
  // all bodies whose flag contains (at least one of) the bits in mark.         
  //                                                                            
{
  CHECKMISMATCH("AddCoeffs");
  setAL(AL);
  setR0(R0);
  CBlock<T> B4(C);
  B4.AddCoeffs(SYM,n,m,x,f,k);
}
//------------------------------------------------------------------------------
template void PotExp::
AddCoeffs<float>(Anlm&, int, const float*,
		 const tupel<3,float>*, const int*, int) const;
template void PotExp::
AddCoeffs<double>(Anlm&, int, const double*,
		  const tupel<3,double>*, const int*, int) const;
//------------------------------------------------------------------------------
void PotExp::Normalize(Anlm&C, scalar const&G) const {
  //                                                                            
  // sets  C_nlm *= G * K_nlm                                                   
  //                                                                            
  CHECKMISMATCH("Normalize");
  if     (SYM & arot) ::normalize<spherical  >(C,Knlm,G);
  else if(SYM & zrot) ::normalize<cylindrical>(C,Knlm,G);
  else if(SYM & axes) ::normalize<triaxial   >(C,Knlm,G);
  else if(SYM & pint) ::normalize<reflexion  >(C,Knlm,G);
  else                ::normalize<none       >(C,Knlm,G);
}
//------------------------------------------------------------------------------
template<typename T>                               // T: double or float        
void PotExp::SetGravity (Anlm const      &C,       // I: C_nlm coefficients     
			 int              n,       // I: # bodies               
			 const tupel<3,T>*x,       // I: positions              
			 T               *p,       // O: potentials             
			 tupel<3,T>      *a,       // O: accelrations           
			 const int       *f,       // I: body flags   see Note 1
			 int              d) const // I: add?         see Note 2
  //                                                                            
  // computes for all bodies or for all active bodies                           
  //                                                                            
  // P(x_i) = - Sum   C_nlm * Psi_nl(r_i) * Y_lm(theta_i,phi_i)                 
  //           n,l,m                                                            
  //                                                                            
  // and -dP/dx(x_i) in cartesian coordinates                                   
  //                                                                            
  // potential    is added if add&1                                             
  // acceleration is added if add&2                                             
  //                                                                            
{
  CHECKMISMATCH("SetGravity");
  setAL(AL);
  setR0(R0);
  GBlock<T> B4(C);
  B4.AddGravity(SYM,n,x,p,a,f,d);
}
//------------------------------------------------------------------------------
template void PotExp::
SetGravity<float>(Anlm const&, int, const tupel<3,float>*,
		  float*, tupel<3,float>*, const int*, int) const;
template void PotExp::
SetGravity<double>(Anlm const&, int, const tupel<3,double>*,
		   double*, tupel<3,double>*, const int*, int) const;
//------------------------------------------------------------------------------
template<typename T>                               // T: double or float        
void PotExp::SetGravity (Anlm       const&C,       // I: C_nlm coefficients     
			 tupel<3,T> const&x,       // I: position               
			 T               &p,       // O: potential              
			 tupel<3,T>      &a,       // O: acceleration           
			 int              d) const // I: add?         see Note 2
  //                                                                            
  //                                                                            
  // computes for a single body:                                                
  //                                                                            
  // P(x_i) = - Sum   C_nlm * Psi_nl(r_i) * Y_lm(theta_i,phi_i)                 
  //           n,l,m                                                            
  //                                                                            
  // and -dP/dx(x_i) in cartesian coordinates                                   
  //                                                                            
  // potential    is added if add&1                                             
  // acceleration is added if add&2                                             
  //                                                                            
{
  CHECKMISMATCH("SetGravity");
  setAL(AL);
  setR0(R0);
  GBlock<T> B4(C);
  B4.AddGravity(SYM,x,p,a,d);
}
//------------------------------------------------------------------------------
template void PotExp::SetGravity<float>(Anlm const&, tupel<3,float> const&,
					float&, tupel<3,float>&, int) const;
template void PotExp::SetGravity<double>(Anlm const&, tupel<3,double> const&,
					 double&, tupel<3,double>&, int) const;
//------------------------------------------------------------------------------
template<typename T>                               // T: double or float        
void PotExp::SetPotential(Anlm const      &C,      // I: C_nlm coefficients     
			  int              n,      // I: # bodies               
			  const tupel<3,T>*x,      // I: positions              
			  T               *p,      // O: potentials             
			  const int       *f,      // I: body flags   see Note 1
			  int            add) const// I: add?         see Note 2
  //                                                                            
  // computes for all bodies or for all active bodies                           
  //                                                                            
  // P(x_i) = - Sum   C_nlm * Psi_nl(r_i) * Y_lm(theta_i,phi_i)                 
  //           n,l,m                                                            
  //                                                                            
  // potential    is added if add&1                                             
  //                                                                            
{
  CHECKMISMATCH("SetPotential");
  setAL(AL);
  setR0(R0);
  PBlock<T> B4(C);
  B4.AddPotential(SYM,n,x,p,f,add);
}
//------------------------------------------------------------------------------
template void PotExp::
SetPotential<float>(Anlm const&, int, const tupel<3,float>*,
		    float*, const int*, int) const;
template void PotExp::
SetPotential<double>(Anlm const&, int, const tupel<3,double>*,
		     double*, const int*, int) const;
//------------------------------------------------------------------------------
template<typename T>                               // T: double or float        
void PotExp::SelfGravity(Anlm            &C,       // O: normalized C_nlm       
			 int              n,       // I: # bodies               
			 const T         *m,       // I: masses                 
			 const tupel<3,T>*x,       // I: positions              
			 T               *p,       // O: potentials             
			 tupel<3,T>      *a,       // O: accelrations           
			 const int       *f,       // I: flags        Notes 1&2 
			 int              k,       // I: mark         see Note 1
			 bool             l,       // I: all          see Note 2
			 int              d,       // I: add?         see Note 3
			 scalar           G) const //[I: const of Gravity]      
{
  CHECKMISMATCH("AddCoeffs");
  setAL(AL);
  setR0(R0);
  C.reset();
  SelfGrav<T> SG(C,n,x);
  if       (SYM & arot) {
    SG.template AddCoef<spherical>(m,f,k);
    ::normalize<spherical>(C,Knlm,G);
    SG.template SetGrav<spherical>(p, a, l?0:f, d);
  } else if(SYM & zrot) {
    SG.template AddCoef<cylindrical>(m,f,k);
    ::normalize<cylindrical>(C,Knlm,G);
    SG.template SetGrav<cylindrical>(p, a, l?0:f, d);
  } else if(SYM & axes) {
    SG.template AddCoef<triaxial>(m,f,k);
    ::normalize<triaxial>(C,Knlm,G);
    SG.template SetGrav<triaxial>(p, a, l?0:f, d);
  } else if(SYM & pint) {
    SG.template AddCoef<reflexion>(m,f,k);
    ::normalize<reflexion>(C,Knlm,G);
    SG.template SetGrav<reflexion>(p, a, l?0:f, d); 
  } else {
    SG.template AddCoef<none>(m,f,k);
    ::normalize<none>(C,Knlm,G);
    SG.template SetGrav<none>(p, a, l?0:f, d);
  }
}
//------------------------------------------------------------------------------
template void PotExp::
SelfGravity<float>(Anlm&, int, const float*, const tupel<3,float>*,
		   float*, tupel<3,float>*, const int*,
		   int, bool, int, scalar) const;
template void PotExp::
SelfGravity<double>(Anlm&, int, const double*, const tupel<3,double>*,
		    double*, tupel<3,double>*, const int*,
		    int, bool, int, scalar) const;
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::AnlmIO                                                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#define TRY_XDR(x,func)							\
  if(!(x)) falcON_THROW("AnlmIO::%s(): XDR operation failed",func);
//------------------------------------------------------------------------------
void AnlmIO::open_for_write(const char*file_name) falcON_THROWING
{
  // open file and connect xdr stream to it                                     
  if(open != closed)
    falcON_THROW("AnlmIO::open_for_write(): already open");
#ifdef falcON_NEMO
  char fname[1024];
  if(is_appended(file_name,'!',fname))
    file = stropen(fname,"w!");
  else
    file = stropen(file_name,"w");
#else
  file = fopen(file_name,"w");
#endif
  if(!file)
    falcON_THROW("cannot open file \"%s\" for writing",file_name);
  if(xdrs == 0) xdrs = new XDR;
  xdrstdio_create( xdrs, file, XDR_ENCODE );
  // create header                                                              
  char header[10] = "anlmfile", *p=header;
  // write identifier
  TRY_XDR(xdr_string(xdrs, &p, 10),"open_for_write");
  open = writing;
}
//------------------------------------------------------------------------------
void AnlmIO::open_for_read(const char*file_name) falcON_THROWING
{
  // open file and connect xdr stream to it                                     
  if(open != closed)
    falcON_THROW("AnlmIO::open_for_read(): already open");
#ifdef falcON_NEMO
  file = stropen(file_name,"r");
#else
  file = fopen(file_name,"r");
#endif
  if(!file)
    falcON_THROW("cannot open file \"%s\" for reading",file_name);
  if(xdrs == 0) xdrs = new XDR;
  xdrstdio_create(xdrs, file, XDR_DECODE);
  // read header                                                                
  char header[10], shead[10] = "anlmfile", *p=header;
  // read identifier
  TRY_XDR(xdr_string(xdrs, &p, 10),"open_for_read");
  if( strcmp(header,shead) )
    falcON_THROWING("file \"%s\" is not a XDR *anlmfile*, "
		    "cannot open for reading",file_name);
  open = reading;
}
//------------------------------------------------------------------------------
void AnlmIO::read(PotExp::symmetry&sym, double&a, double&r,
		  PotExp::Anlm&A, double &t)
  falcON_THROWING
{
  if(open != reading)
    falcON_THROW("AnlmIO::read(): stream not opened for reading");
  if(feof(file)) falcON_THROW("AnlmIO::read(): seen end of file\n");
  if(ferror(file)) falcON_THROW("AnlmIO::read(): I/O error\n");
  int n,l,s;
  TRY_XDR(xdr_double(xdrs,&t),"read");                 // read time
  TRY_XDR(xdr_double(xdrs,&a),"read");                 // read alpha
  TRY_XDR(xdr_double(xdrs,&r),"read");                 // read scale
  TRY_XDR(xdr_int   (xdrs,&s),"read");                 // read symmetry
  TRY_XDR(xdr_int   (xdrs,&n),"read");                 // read max n
  TRY_XDR(xdr_int   (xdrs,&l),"read");                 // read max l
  sym = PotExp::symmetry(s);
  A.reset(n,l);
  Anlm::scalar* aN=A.A+(n+1)*square(l+1);
  for(Anlm::scalar*a=A.A; a!=aN; ++a)
    TRY_XDR(xdr_double(xdrs,a),"read");
}
//------------------------------------------------------------------------------
void AnlmIO::write(PotExp::symmetry sym, double a, double r,
		   PotExp::Anlm const& A, double t)
  falcON_THROWING
{
  if(open != writing)
    falcON_THROW("AnlmIO::write(): stream not opened for writing");
  int s = sym;
  TRY_XDR(xdr_double(xdrs,&t),"write");                           // write time
  TRY_XDR(xdr_double(xdrs,&a),"write");                           // write alpha
  TRY_XDR(xdr_double(xdrs,&r),"write");                           // write scale
  TRY_XDR(xdr_int   (xdrs,&s),"write");                           // write symm
  TRY_XDR(xdr_int   (xdrs,const_cast<int*>(&(A.nmax()))),"write");// write max n
  TRY_XDR(xdr_int   (xdrs,const_cast<int*>(&(A.lmax()))),"write");// write max l
  Anlm::scalar* aN=A.A+(A.nmax()+1)*square(A.lmax()+1);
  for(const Anlm::scalar*a=A.A; a!=aN; ++a)
    TRY_XDR(xdr_double(xdrs,const_cast<double*>(a)),"write");
}
//------------------------------------------------------------------------------
bool AnlmIO::is_good() const {
  return !feof(file) && !ferror(file);
}
//------------------------------------------------------------------------------
void AnlmIO::close()
{
  if(open) {
    xdr_destroy (xdrs);
    falcON_DEL_O(xdrs);
#ifdef falcON_NEMO
    strclose( file );
#else
    fclose  ( file );
#endif
  }
  open = closed;
  xdrs = 0;
  file = 0;
}
#undef TRY_XDR
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#if (0)
////////////////////////////////////////////////////////////////////////////////
#include <iomanip>
#include <iostream>
#include <ctime>
using namespace std;
int main()
{
  double al;
  int _sym,nmax,lmax;
  cout<<" symmetry: 0 -> none\n"
      <<"           1 -> reflexion\n"
      <<"           2 -> triaxial\n"
      <<"           3 -> cylindrical\n"
      <<"           4 -> spherical\n"
      <<" symmetry    = ";
  cin >> _sym;
  cout<<" alpha       = ";
  cin >> al;
  cout<<" lmax (even) = ";
  cin >> lmax;
  cout<<" nmax        = ";
  cin >> nmax;
  symmetry sym(_sym==4? PotExp::spherical   :
	       _sym==3? PotExp::cylindrical :
	       _sym==2? PotExp::triaxial    :
	       _sym==1? PotExp::reflexion   : PotExp::none);

  PotExp PE(al,1.,nmax,lmax,sym );
  cout<<" PotExp::symmetry: \""<<PE.symmetry_name()<<"\""<<endl;

#if(0) // test SetYlm
  cout<<"\n testing SetYlm()\n\n";
  YlmRec Y00(lmax),Yl0(lmax),Yh0(lmax),Y0l(lmax),Y0h(lmax);
  YlmRec Y(lmax),Yt(lmax),Yp(lmax);
  for(;;) {
    const double d=1.e-4, td=d+d;
    double the,phi;
    cout<<" cos(theta) in [-1,1] = "; cin>>the; the = acos(the);
    cout<<" phi                  = "; cin>>phi;
    SetYlm(sym,Y00,cos(the  ),sin(the  ),cos(phi  ),sin(phi  ));
    SetYlm(sym,Y0l,cos(the  ),sin(the  ),cos(phi-d),sin(phi-d));
    SetYlm(sym,Y0h,cos(the  ),sin(the  ),cos(phi+d),sin(phi+d));
    SetYlm(sym,Yl0,cos(the-d),sin(the-d),cos(phi  ),sin(phi  ));
    SetYlm(sym,Yh0,cos(the+d),sin(the+d),cos(phi  ),sin(phi  ));
    SetYlm(sym,Y,Yt,Yp,cos(the  ),sin(the  ),cos(phi  ),sin(phi  ));

#define PRINTY(__L,__M)							\
    {									\
      scalar								\
        y  =Y00(__L,__M),						\
	yt =(Yh0(__L,__M)-Yl0(__L,__M))/td,				\
	yp =(Y0h(__L,__M)-Y0l(__L,__M))/td,				\
	yy =Y(__L,__M),							\
	yyt=Yt(__L,__M),						\
	yyp=Yp(__L,__M);						\
      cout<< " (" << setw( 1) <<  (__L) << ',' << setw( 2) << (__M)	\
	  << "): Y="  << setw(13) << y					\
	  << " dY/dt="<< setw(13) << yt					\
	  << " dY/dp="<< setw(13) << yp					\
	  << "\n        "						\
	  << " Y="    << setw(13) << yy					\
	  << " dY/dt="<< setw(13) << yyt				\
	  << " dY/dp="<< setw(13) << yyp				\
	  << (abs(y -yy )>d*abs(y) ? " PROBLEMS: Y  " : " ")		\
	  << (abs(yt-yyt)>d*abs(yt)? " PROBLEMS: Yt " : " ")		\
	  << (abs(yp-yyp)>d*abs(yp)? " PROBLEMS: Yp " : " ") << endl;	\
    }

    cout<<'\n';
    PRINTY(0,0)
    if(sym & PotExp::arot) continue;
    int dl = (sym & PotExp::pint)? 2 : 1;
    for(int l=dl; l<=lmax; l+=dl) {
      cout<<'\n';
      if     (sym & PotExp::zrot)
	PRINTY(l,0)
      else if(sym & PotExp::axes)
	for(int m=0; m<=l; m+=2)
	  PRINTY(l,m)
      else if(sym & PotExp::pint)
	for(int m=-l; m<=l; m+=2)
	  PRINTY(l,m)
      else
	for(int m=-l; m<=l; ++m)
	  PRINTY(l,m)
    }
  }
#endif

#if(0) // test PotExp::SetPsi
  cout<<"\n testing SetPsi()\n\n";
  AnlRec P0(nmax,lmax), Pl(nmax,lmax), Ph(nmax,lmax);
  AnlRec PP(nmax,lmax), Pr(nmax,lmax);
  for(;;) {
    const double d=1.e-4, td=d+d;
    double r;
    cout<<" r/r0 = "; cin>>r;
    SetPsi(sym,P0,r,1);
    SetPsi(sym,Pl,r-d,1);
    SetPsi(sym,Ph,r+d,1);
    SetPsi(sym,PP,Pr,r);
#define PRINTP(__N,__L)						\
    cout<< " (" << (__N) << ',' << (__L)			\
	<< "): P="  << setw(13) <<  P0(__N,__L)			\
	<< " dP/dr="<< setw(13) << (Ph(__N,__L)-Pl(__N,__L))/td	\
	<< "\n       "						\
	<< " P="    << setw(13) <<  PP(__N,__L)			\
	<< " dP/dr="<< setw(13) <<  Pr(__N,__L)	<<endl;

    int dl = (sym & PotExp::pint)? 2 : 1;
    for(int n=0; n<=nmax; ++n) {
      cout<<'\n';
      PRINTP(n,0);
      if(! sym & PotExp::arot)
	for(int l=dl; l<=lmax; l+=dl)
	  PRINTP(n,l)
    }
  }
#endif

#if(0) // test Add(Anlm)

  cout<<"\n testing Ass(Anlm)\n\n";
  AnlRec P(nmax,lmax);
  YlmRec Y(lmax);
  Anlm   A(PE);

  for(;;) {
    double r,ct,st,phi,cp,sp;
    cout<<" r/r0 = "; cin>>r;
    cout<<" cos(theta) in [-1,1] = "; cin>>ct; st=sqrt(1.-ct*ct);
    cout<<" phi                  = "; cin>>phi;
    A = PotExp::scalar(0);
    PE.SetPsi(P,r,1);
    PE.SetYlm(Y,ct,st,cos(phi),sin(phi));
    PE.Add(A,P,Y);
    
#define PRINTA(__N,__L,__M)						\
    cout<< " A("  << setw( 1) << (__N)					\
        << ','    << setw( 1) << (__L)					\
	<< ','    << setw( 2) << (__M)					\
	<< ") ="  << setw(12) <<  A(__N,__L,__M)			\
	<< " P("  << setw( 1) << (__N)					\
	<< ','    << setw( 1) << (__L)					\
	<< ")*Y(" << setw( 1) << (__L)					\
	<< ','    << setw( 2) << (__M)					\
	<< ") ="  << setw(12) << P(__N,__L) * Y(__L,__M)		\
        <<((abs(A(__N,__L,__M) - P(__N,__L) * Y(__L,__M))>1.e-8)?	\
           "        PROBLEM\n" : "\n");

    cout<<'\n';
    for(int n=0; n<=nmax; ++n) {
      cout<<'\n';
      PRINTA(n,0,0)
      if(sym & PotExp::arot) continue;
      int dl = (sym & PotExp::pint)? 2 : 1;
      for(int l=dl; l<=lmax; l+=dl) {
	cout<<'\n';
	if     (sym & PotExp::zrot)
	  PRINTA(n,l,0)
	else if(sym & PotExp::axes)
	  for(int m=0; m<=l; m+=2)
	    PRINTA(n,l,m)
	else if(sym & PotExp::pint)
	  for(int m=-l; m<=l; m+=2)
	    PRINTA(n,l,m)
        else
	  for(int m=-l; m<=l; ++m)
	    PRINTA(n,l,m)
      }
    }
  }
#endif
  
#if(0) // test PotExp::Spherical4()
  int N;
  cout<<" testing PotExp::Spherical4()\n\n"
      <<" N = "; cin>>N;

  const int N4 = 4*N;
  vect  Y[4];
  vect *X = new vect[N4];
  sphc *S = new PotExp::sphc[N4], S1(0.), S2(0.);S[0]
  for(int i=0; i!=4; ++i) {
    cout << " X["<<i<<"] = ";
    cin  >> Y[i];
  }
  for(int n=0; n!=N; ++n)
    for(int i=0; i!=4; ++i)
      X[4*n+i] = Y[i];

  register clock_t cpu0;

  cpu0 = clock();
  for(int i=0; i!=N4; ++i)
    Spherical(S[i],X[i]);
  register real t1 = (clock() - cpu0)/real(CLOCKS_PER_SEC);
  for(int i=0; i!=N4; ++i) S1 += S[i];
  cout<<" PotExp::Spherical():\n"
      <<"   time      = "<<t1<<" sec\n"
      <<"   Sum_i S_i = "<<S1<<endl;

  cpu0 = clock();
  for(int i=0; i!=N4; i+=4)
    Spherical4(S+i,X+i);
  register real t2 = (clock() - cpu0)/real(CLOCKS_PER_SEC);
  for(int i=0; i!=N4; ++i) S2 += S[i];
  cout<<" PotExp::Spherical4():\n"
      <<"   time      = "<<t2<<" sec\n"
      <<"   Sum_i S_i = "<<S2<<endl;

//   PotExp::Spherical4(S,X);
//   for(int i=0; i!=4; ++i)
//     cout<<" S["<<i<<"] = "<<Q[i]<<"\n      = "<<S[i]<<endl;
#endif

#if(0) // test PotExp::Cartesian4()
  cout<<" testing PotExp::Cartesian4()\n\n";
  vect X[4], A[4], B[4];
  float Rd,Ct,St,Cp,Sp;
  fvec4 rd,ct,st,cp,sp;
  for(;;) {
    cout << " (x ,y ,z ) = "; cin >> X[0];
    cout << " (ar,at,ap) = "; cin >> A[0];
    Spherical (Rd,Ct,St,Cp,Sp,X[0]);
    Spherical4(rd,ct,st,cp,sp,X);
    cout<<" S = "<<Rd   <<' '<<Ct   <<' '<<St   <<' '<<Cp   <<' '<<Sp   <<'\n'
	<<"   = "<<rd[0]<<' '<<ct[0]<<' '<<st[0]<<' '<<cp[0]<<' '<<sp[0]<<'\n';
    B[0] = A[0];
    Cartesian (B[0],Ct,St,Cp,Sp);
    Cartesian4(A,ct,st,cp,sp);
    cout<<" A = "<<B[0]<<'\n'
	<<"   = "<<A[0]<<'\n';
  }
#endif
}
#endif // TESTING
