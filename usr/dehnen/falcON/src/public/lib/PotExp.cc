// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
/// \file src/public/lib/PotExp.cc                                             |
//                                                                             |
// Copyright (C) 1994-1996, 2004-2009 Walter Dehnen                            |
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
#undef MAX
#undef MIN
#ifdef   falcON_NEMO
#  include <public/nemo++.h>
#  include <cstring>
#endif

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
  using   falcON::tupel;
  using   falcON::fvec4;
  typedef PotExp::scalar   scalar;
  typedef PotExp::symmetry symmetry;
  typedef PotExp::Anlm     Anlm;
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
#if 0 // not used
    scalar const&e (int l, int m) const { return A[l*(l+1)+m]; }
#endif
  public:
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // construction & destruction                                               
    YlmRec(int l)
      : L(l),   L1(L+1), L1Q(L1*L1), A(falcON_NEW(scalar,L1Q)) {}
#if 0 // not used
    YlmRec(YlmRec const&a)
      : L(a.L), L1(L+1), L1Q(L1*L1), A(falcON_NEW(scalar,L1Q)) {
      memcpy(A,a.A,L1Q*sizeof(scalar));
    }
#endif
    ~YlmRec() { falcON_DEL_A(A); }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // element access                                                           
    //    NOTE: that elements ASSUMED to be zero may have any value!            
    scalar      &operator() (int l, int m)       { return e(l,m); }
#if 0 // not used
    scalar const&operator() (int l, int m) const { return e(l,m); }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    int const& lmax() const { return L; }
#endif
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void table_print(symmetry, std::ostream&, int=6) const;
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  private: // only accessible to PotExp in file pexp.cc                         
    // set Ylm to non-normalized scalar-valued spherical harmonics, i.e.        
    //    Y(l,m;theta,phi) = P(l,|m|; cos[theta]) * cas(m*phi)                  
    template<symmetry>
    YlmRec&set(scalar ct, scalar st,               // I: cos(the),sin(the)      
	       scalar cp, scalar sp);              // I: cos(phi),sin(phi)      
    template<symmetry>
    YlmRec&set(YlmRec&, YlmRec&,                   // O: dY/dthe, dY/dphi       
	       scalar ct, scalar st,               // I: cos(the),sin(the)      
	       scalar cp, scalar sp);              // I: cos(phi),sin(phi)      
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#if 0 // not used
    // set Ylm to N_lm multiplied by 4Pi := (2l+1) (l-|m|)! / (l+|m|)!          
    YlmRec& Nlm();
#endif
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    template<symmetry> YlmRec&reset   ();              // A[l,m] = 0            
    template<symmetry> YlmRec&assign  (scalar);        // A[l,m] = x            
    template<symmetry> YlmRec&multiply(scalar);        // A[l,m]*= x            
    template<symmetry> YlmRec&divide  (scalar);        // A[l,m]/= x            
    template<symmetry> YlmRec&copy    (YlmRec const&); // A[l,m] = B[l,m]       
    template<symmetry> YlmRec&add     (YlmRec const&); // A[l,m]+= B[l,m]       
    template<symmetry> YlmRec&sub     (YlmRec const&); // A[l,m]-= B[l,m]       
    template<symmetry> scalar dot     (YlmRec const&) const; // dot product     
    template<symmetry> YlmRec&addtimes(YlmRec const&,
				       scalar);        // A[l,m]+= x*B[l,m]     
    template<symmetry> YlmRec&subtimes(YlmRec const&,
				       scalar);        // A[l,m]-= x*B[l,m]     
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
#if 0 // not used
    AnlRec(AnlRec const&a)
      : N1(a.N1), L1(a.L1), A(falcON_NEW(scalar,N1*L1)) {
      memcpy(A,a.A,N1*L1*sizeof(scalar));
    }
#endif
    ~AnlRec() { falcON_DEL_A(A); }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // element access                                                           
    scalar      &operator() (int n, int l)       { return A[n*L1+l]; }
    scalar const&operator() (int n, int l) const { return A[n*L1+l]; }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void table_print(symmetry, std::ostream&, int=6) const;
  private: // only accessible to PotExp in file pexp.cc                         
    template<symmetry> AnlRec&reset   ();              // A[n,l] = 0            
    template<symmetry> AnlRec&assign  (scalar);        // A[n,l] = x            
    template<symmetry> AnlRec&multiply(scalar);        // A[n,l]*= x            
    template<symmetry> AnlRec&divide  (scalar);        // A[n,l]/= x            
    template<symmetry> AnlRec&copy    (AnlRec const&); // A[n,l] = B[n,l]       
    template<symmetry> AnlRec&add     (AnlRec const&); // A[n,l]+= B[n,l]       
    template<symmetry> AnlRec&sub     (AnlRec const&); // A[n,l]-= B[n,l]       
    template<symmetry> scalar dot     (AnlRec const&) const; // dot product     
    template<symmetry> AnlRec&addtimes(AnlRec const&,
				       scalar);        // A[n,l]+= x*B[n,l]     
    template<symmetry> AnlRec&subtimes(AnlRec const&,
				       scalar);        // A[n,l]-= x*B[n,l]     
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // static variables to hold info about alpha used by static methods           
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  scalar AL=1, iAL=1, AL1=2;
  inline void setAL(scalar al) {
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
  using falcON::P::R0;
  using falcON::P::IR0;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // auxiliary inline functions                                                 
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  inline symmetry operator & (symmetry a, symmetry b) {
    return PotExp::symmetry(int(a) & int(b));
  }
  //----------------------------------------------------------------------------
#if 0 // not used
  inline symmetry operator | (symmetry a, symmetry b) {
    return PotExp::symmetry(int(a) | int(b));
  }
#endif
  //----------------------------------------------------------------------------
  inline scalar lambda(int l) {
    return AL*(l+l+1)+0.5;
  }
  //----------------------------------------------------------------------------
  inline void SetXiFi(scalar&xi, scalar&fi, scalar r)
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
		      scalar r)
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
    static int const &L  (YlmRec const&Y) { return Y.L; }
    static int const &L1 (YlmRec const&Y) { return Y.L1; }
    static int const &L1Q(YlmRec const&Y) { return Y.L1Q; }
    static int const &N1 (AnlRec const&A) { return A.N1; }
    static int const &L1 (AnlRec const&A) { return A.L1; }
#if 0 // not used
    static int const &N  (Anlm   const&A) { return A.N; }
    static int const &L  (Anlm   const&A) { return A.L; }
#endif
    static int const &N1 (Anlm   const&A) { return A.N1; }
    static int const &L1 (Anlm   const&A) { return A.L1; }
    static int const &L1Q(Anlm   const&A) { return A.L1Q; }

    static scalar      &e(YlmRec      &C, int i) { return C.A[i]; }
    static scalar const&e(YlmRec const&C, int i) { return C.A[i]; }
    static scalar      &e(AnlRec      &C, int i) { return C.A[i]; }
    static scalar const&e(AnlRec const&C, int i) { return C.A[i]; }
    static scalar      &e(Anlm        &C, int i) { return C.A[i]; }
    static scalar const&e(Anlm   const&C, int i) { return C.A[i]; }

#if 0 // not used
    static scalar      *p(YlmRec      &C) { return C.A; }
    static scalar      *p(AnlRec      &C) { return C.A; }
#endif
    static scalar const*p(YlmRec const&C) { return C.A; }
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
    static void Connect(YlmRec&A, YlmRec const&B, scalar x) {
      for(int i=0; i!=L1Q(A); ++i)
	Connector::op(e(A,i),e(B,i),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(AnlRec&A, AnlRec const&B, scalar x) {
      for(int i=0; i!=N1(A)*L1(A); ++i)
	Connector::op(e(A,i),e(B,i),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(Anlm&A, Anlm const&B, scalar x) {
      for(int i=0; i!=N1(A)*L1Q(A); ++i)
	Connector::op(e(A,i),e(B,i),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(Anlm&A, AnlRec const&P, YlmRec const&Y, scalar x) {
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
#if 0 // not used
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
#endif
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
    template<typename _T>
    static scalar Dot(tupel<3,_T>&dx,
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
    static void SetPlm(YlmRec&Y, scalar ct, scalar st) {
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
    static void SetPlm(YlmRec&Y, YlmRec&T, scalar ct, scalar st) {
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
		       scalar ct, scalar st, scalar cp, scalar sp) {
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
		       scalar ct, scalar st, scalar cp, scalar sp) {
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
    static void SetPsi(AnlRec&P, scalar r, scalar GM) {
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
    static void SetPsi(AnlRec&P, AnlRec&D, scalar r) {
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
    void Connect(YlmRec&A, YlmRec const&B, scalar x) {
      for(int l=0,i=0; l<L1(A); l+=2,i+=4*l-2)
	for(int m=-l; m<=l; m+=2)
	  Connector::op(e(A,i+m),e(B,i+m),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(AnlRec&A, AnlRec const&B, scalar x)
    {
      for(int n=0,j=0; n!=N1(A); ++n,j+=L1(A))
	for(int l=0,i=j; l<L1(A); l+=2,i+=2)
	  Connector::op(e(A,i),e(B,i),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(Anlm&A, Anlm const&B, scalar x)
    {
      for(int n=0,j=0; n!=N1(A); ++n,j+=L1Q(A))
	for(int l=0; l<L1(A); l+=2)
	  for(int m=-l,i=j+l*l; m<=l; m+=2,i+=2)
	    Connector::op(e(A,i),e(B,i),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(Anlm&A, AnlRec const&P, YlmRec const&Y, scalar x) {
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
#if 0 // not used
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
#endif
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
    template<typename _T>
    static scalar Dot(tupel<3,_T>&dx,
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
    static void SetPlm(YlmRec&Y, scalar ct, scalar st) {
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
    static void SetPlm(YlmRec&Y, YlmRec&T, scalar ct, scalar st) {
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
		       scalar ct, scalar st, scalar cp, scalar sp) {
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
		       scalar ct, scalar st, scalar cp, scalar sp) {
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
    static void SetPsi (AnlRec&P, scalar r, scalar GM) {
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
    static void SetPsi(AnlRec&P, AnlRec&D, scalar r) {
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
    static void Connect(YlmRec&A, YlmRec const&B, scalar x) {
      for(int l=0,i=0; l<L1(A); l+=2,i+=4*l-2)
	for(int m=0; m<=l; m+=2)
	  Connector::op(e(A,i+m),e(B,i+m),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(AnlRec&A, AnlRec const&B, scalar x) {
      AUX<PotExp::reflexion>:: template Connect<Connector>(A,B,x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(Anlm&A, Anlm const&B, scalar x) {
      for(int n=0,j=0; n!=N1(A); ++n,j+=L1Q(A))
	for(int l=0; l<L1(A); l+=2)
	  for(int m=0,i=j+l*l+l; m<=l; m+=2,i+=2)
	    Connector::op(e(A,i),e(B,i),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(Anlm&A, AnlRec const&P, YlmRec const&Y, scalar x) {
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
#if 0 // not used
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
#endif
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
    template<typename _T>
    static scalar Dot(tupel<3,_T>&dx,
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
    static void SetYlm(YlmRec&Y,
		       scalar ct, scalar st, scalar cp, scalar sp) {
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
		       scalar ct, scalar st, scalar cp, scalar sp) {
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
    static void SetPsi (AnlRec&P, scalar r, scalar GM) {
      AUX<PotExp::reflexion>::SetPsi(P,r,GM);
    }
    //--------------------------------------------------------------------------
    static void SetPsi (AnlRec&P, AnlRec&D, scalar r) {
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
    static void Connect(YlmRec&A, YlmRec const&B, scalar x) {
      for(int l=0,i=0; l<L1(A); l+=2,i+=4*l-2)
	Connector::op(e(A,i),e(B,i),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(AnlRec&A, AnlRec const&B, scalar x) {
      AUX<PotExp::reflexion>::template Connect<Connector>(A,B,x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(Anlm&A, Anlm const&B, scalar x) {
      for(int n=0,j=0; n!=N1(A); ++n,j+=L1Q(A))
	for(int l=0,i=j; l<L1(A); l+=2, i+=4*l-2)
	  Connector::op(e(A,i),e(B,i),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(Anlm&A, AnlRec const&P, YlmRec const&Y, scalar x) {
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
#if 0 // not used
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
#endif
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
    template<typename _T>
    static scalar Dot(tupel<3,_T>&dx,
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
    static void SetYlm(YlmRec&Y,
		       scalar ct, scalar st, scalar, scalar) {
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
		       scalar ct, scalar st, scalar, scalar) {
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
    static void SetPsi (AnlRec&P, scalar r, scalar GM) {
      AUX<PotExp::reflexion>::SetPsi(P,r,GM);
    }
    //--------------------------------------------------------------------------
    static void SetPsi (AnlRec&P, AnlRec&D, scalar r) {
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
    static void Connect(YlmRec&A, YlmRec const&B, scalar x) {
      Connector::op(e(A,0),e(B,0),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(AnlRec&A, AnlRec const&B, scalar x) {
      for(int n=0,i=0; n!=N1(A); ++n,i+=L1(A))
	Connector::op(e(A,i),e(B,i),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(Anlm&A, Anlm const&B, scalar x) {
      for(int n=0,i=0; n!=N1(A); ++n,i+=L1Q(A))
	Connector::op(e(A,i),e(B,i),x);
    }
    //--------------------------------------------------------------------------
    template<class Connector>
    static void Connect(Anlm&A, AnlRec const&P, YlmRec const&Y, scalar x) {
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
#if 0 // not used
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
#endif
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
    template<typename _T>
    static scalar Dot(tupel<3,_T>&dx,
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
    static void SetYlm(YlmRec&Y, scalar, scalar, scalar, scalar) {
      e(Y,0) = 1;                                  // Y(0,0) = 1                
    }
    //--------------------------------------------------------------------------
    static void SetYlm(YlmRec&Y, YlmRec&T, YlmRec&P,
		       scalar, scalar, scalar, scalar) {
      e(Y,0) = 1;                                  // Y(0,0) = 1                
      e(T,0) = 0;                                  // dY(0,0)/dthe = 0          
      e(P,0) = 0;                                  // dY(0,0)/dphi = 0          
    }
    //--------------------------------------------------------------------------
    static void SetPsi (AnlRec&P, scalar r, scalar GM) {
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
    static void SetPsi (AnlRec&P, AnlRec&D, scalar r) {
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
  void SetYlm (YlmRec&Y, scalar ct, scalar st, scalar cp, scalar sp) {
    AUX<S>::SetYlm(Y,ct,st,cp,sp);
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline
  void SetYlm (YlmRec&Y, YlmRec&Yt, YlmRec&Yp,
	       scalar ct, scalar st, scalar cp, scalar sp) {
    AUX<S>::SetYlm(Y,Yt,Yp,ct,st,cp,sp);
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline
  void SetPsi (AnlRec&P, scalar r, scalar GM) {
    AUX<S>::SetPsi(P,r,GM);
  }
  //----------------------------------------------------------------------------
  template<symmetry S> inline
  void SetPsi (AnlRec&P, AnlRec&Pr, scalar r) {
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
  //----------------------------------------------------------------------------
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
  scalar(*fu)(scalar);
  scalar(*fb)(scalar,scalar);
  scalar(*ft)(scalar,scalar,scalar);
  struct __neg  { static void op(scalar&a, scalar  , scalar )  {a = -a; } };
  struct __setX { static void op(scalar&a, scalar  , scalar x) {a =  x; } };
  struct __mulX { static void op(scalar&a, scalar  , scalar x) {a*=  x; } };
  struct __setB { static void op(scalar&a, scalar b, scalar )  {a =  b; } };
  struct __mulB { static void op(scalar&a, scalar b, scalar )  {a*=  b; } };
  struct __addB { static void op(scalar&a, scalar b, scalar )  {a+=  b; } };
  struct __subB { static void op(scalar&a, scalar b, scalar )  {a-=  b; } };
#if 0 // not used
  struct __setT { static void op(scalar&a, scalar b, scalar x) {a =x*b; } };
#endif
  struct __addT { static void op(scalar&a, scalar b, scalar x) {a+=x*b; } };
  struct __subT { static void op(scalar&a, scalar b, scalar x) {a-=x*b; } };
  struct __una  { static void op(scalar&a, scalar  , scalar  ) {a =fu(a);  } };
  struct __binX { static void op(scalar&a, scalar  , scalar x) {a =fb(a,x);} };
  struct __binB { static void op(scalar&a, scalar b, scalar  ) {a =fb(a,b);} };
  struct __tert { static void op(scalar&a, scalar b, scalar x) {a =ft(a,b,x);}};
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
  inline void Cartesian(tupel<3,Y> &a, X rd, X ct, X st, X cp, X sp) {
    if(rd) {
      register X ir = IR0/rd;
      a[1] *= ir;
      if(st) a[2] *= ir/st;
      else   a[2]  = Y(0);
    } else {
      a[1] = Y(0);
      a[2] = Y(0);
    }
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
  void Cartesian4(tupel<3,double>*, fvec4 const&,
		  fvec4 const&, fvec4 const&, fvec4 const&, fvec4 const&);
  void Cartesian4(tupel<3,float>*, fvec4 const&,
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
  void Cartesian4(tupel<3,X>*a, fvec4 const&rd,
		  fvec4 const&ct, fvec4 const&st,
		  fvec4 const&cp, fvec4 const&sp) {
    Cartesian(a[0],rd[0],ct[0],st[0],cp[0],sp[0]);
    Cartesian(a[1],rd[1],ct[1],st[1],cp[1],sp[1]);
    Cartesian(a[2],rd[2],ct[2],st[2],cp[2],sp[2]);
    Cartesian(a[3],rd[3],ct[3],st[3],cp[3],sp[3]);
  }
#endif
} }
using namespace falcON::P;
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class YlmRec                                                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
template<symmetry S> inline
YlmRec &YlmRec::set(scalar ct, scalar st, scalar cp, scalar sp) {
  AUX<S>::SetYlm(*this,ct,st,cp,sp);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline
YlmRec &YlmRec::set(YlmRec&T, YlmRec&P,
		    scalar ct, scalar st, scalar cp, scalar sp) {
  AUX<S>::SetYlm(*this,T,P,ct,st,cp,sp);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline YlmRec &YlmRec::assign(scalar x) {
  AUX<S>::template Connect<__setX>(*this,*this,x);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline YlmRec &YlmRec::reset() {
  return assign<S>(scalar(0));
}
//------------------------------------------------------------------------------
template<symmetry S> inline YlmRec &YlmRec::multiply(scalar x) {
  AUX<S>::template Connect<__mulX>(*this,*this,x);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline YlmRec &YlmRec::divide(scalar x) {
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
template<symmetry S> inline YlmRec &YlmRec::addtimes(YlmRec const&B, scalar x) {
  AUX<S>::template Connect<__addT>(*this,B,x);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline YlmRec &YlmRec::sub(YlmRec const&B) {
  AUX<S>::template Connect<__subB>(*this,B,scalar(0));
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline YlmRec &YlmRec::subtimes(YlmRec const&B, scalar x) {
  AUX<S>::template Connect<__subT>(*this,B,x);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline scalar YlmRec::dot(YlmRec const&B) const {
  return AUX<S>::Dot(*this,B);
}
//------------------------------------------------------------------------------
#if 0 // not used
inline YlmRec& YlmRec::Nlm() {
  SetNlm(*this);
  return *this;
}
#endif
//------------------------------------------------------------------------------
void YlmRec::table_print(symmetry     s,
			 std::ostream&o,
			 int          p) const
{
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
template<symmetry S> inline AnlRec &AnlRec::assign(scalar x) {
  AUX<S>::template Connect<__setX>(*this,*this,x);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline AnlRec &AnlRec::reset() {
  return assign<S>(scalar(0));
}
//------------------------------------------------------------------------------
template<symmetry S> inline AnlRec &AnlRec::multiply(scalar x) {
  AUX<S>::template Connect<__mulX>(*this,*this,x);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline AnlRec &AnlRec::divide(scalar x) {
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
template<symmetry S> inline AnlRec &AnlRec::addtimes(AnlRec const&B, scalar x) {
  AUX<S>::template Connect<__addT>(*this,B,x);
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline AnlRec &AnlRec::sub(AnlRec const&B) {
  AUX<S>::template Connect<__subB>(*this,B,scalar(0));
  return *this;
}
//------------------------------------------------------------------------------
template<symmetry S> inline AnlRec &AnlRec::subtimes(AnlRec const&B, scalar x) {
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
#define CONNECT(C,B,X)							   \
  switch(S) {								   \
  case spherical:   AUX<spherical>  :: Connect<C>(*this,B,X); return*this; \
  case cylindrical: AUX<cylindrical>:: Connect<C>(*this,B,X); return*this; \
  case triaxial:    AUX<triaxial>   :: Connect<C>(*this,B,X); return*this; \
  case reflexion:   AUX<reflexion>  :: Connect<C>(*this,B,X); return*this; \
  default:	    AUX<none>       :: Connect<C>(*this,B,X); return*this; \
  }
//----------------------------------------------------------------------------
Anlm&Anlm::assign(scalar x, symmetry S) {
  CONNECT(__setX,*this,x);
}
Anlm&Anlm::negate(symmetry S) {
  CONNECT(__neg ,*this,scalar(0));
}
Anlm&Anlm::multiply(scalar x, symmetry S) {
  CONNECT(__mulX,*this,x);
}
Anlm&Anlm::copy(Anlm const&B, symmetry S) {
  CONNECT(__setB,B,scalar(0));
}
Anlm&Anlm::add(Anlm const&B, symmetry S) {
  CONNECT(__addB,B,scalar(0));
}
Anlm&Anlm::subtract(Anlm const&B, symmetry S) {
  CONNECT(__subB,B,scalar(0));
}
Anlm&Anlm::multiply(Anlm const&B, symmetry S) {
  CONNECT(__mulB,B,scalar(0));
}
Anlm&Anlm::addtimes(Anlm const&B, scalar x, symmetry S) {
  CONNECT(__addT,B,x);
}
Anlm&Anlm::subtimes(Anlm const&B, scalar x, symmetry S) {
  CONNECT(__subT,B,x);
}
Anlm&Anlm::unary(scalar(*f)(scalar), symmetry S) {
  ::fu = f;
  CONNECT(__una,*this,scalar(0));
}
Anlm&Anlm::binary(scalar(*f)(scalar,scalar), scalar x, symmetry S) {
  ::fb = f;
  CONNECT(__binX,*this,x);
}
Anlm&Anlm::binary(scalar(*f)(scalar,scalar), Anlm const&B, symmetry S) {
  ::fb = f;
  CONNECT(__binB,B,scalar(0));
}
Anlm&Anlm::tertiary(scalar(*f)(scalar,scalar,scalar), Anlm const&B, scalar x,
		    symmetry S) {
  ::ft = f;
  CONNECT(__tert,B,x);
}
#undef CONNECT
//------------------------------------------------------------------------------
scalar Anlm::dot(Anlm const&B, symmetry S) const {
  switch(S) {
  case spherical:   return AUX<spherical>  ::Dot(*this,B);
  case cylindrical: return AUX<cylindrical>::Dot(*this,B);
  case triaxial:    return AUX<triaxial>   ::Dot(*this,B);
  case reflexion:   return AUX<reflexion>  ::Dot(*this,B);
  default:          return AUX<none>       ::Dot(*this,B);
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
//------------------------------------------------------------------------------
#ifdef falcON_MPI
Anlm& Anlm::global_sum(Anlm const&C, const MPI::Communicator*Comm)
{
  reset(C.nmax(),C.lmax());
  COMMUN(Comm)->AllReduce<MPI::Sum>(C.A,A,N1*L1Q);
  return*this;
}
Anlm& Anlm::global_sum(const MPI::Communicator*Comm)
{
  COMMUN(Comm)->AllReduceInPlace<MPI::Sum>(A,N1*L1Q);
  return*this;
}
#endif
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
    void load(T m, V const&x) {                    // load a body into buffer   
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
	AUX<SYM>::template Connect<__addB>(C,Psi,Ylm,scalar(0));
	                                           //     add to C_nlm          
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
    //--------------------------------------------------------------------------
    template<symmetry SYM>                         // symmetry at compile time  
    static void AddCoeffs(int n, int m,
			  const V *x,
			  const T**y,
			  Anlm    *C)
    {
      PotExp::scalar rd,ct,st,cp,sp;
      AnlRec         Psi(C->nmax(),C->lmax());
      YlmRec         Ylm(C->lmax());
      for(int i=0; i!=n; ++i) {
	Spherical(rd,ct,st,cp,sp,x[i]);
	SetPsi<SYM>(Psi,rd,T(1));
	SetYlm<SYM>(Ylm,ct,st,cp,sp);
	for(int j=0; j!=m; ++m)
	  AUX<SYM>::template Connect<__addT>(C[j],Psi,Ylm,y[i][j]);
      }
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
    //--------------------------------------------------------------------------
    static void AddCoeffs(symmetry,
			  int, int,
			  const V*,
			  const T**,
			  Anlm*);
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
  //----------------------------------------------------------------------------
  template<typename T>
  void CBlock<T>::AddCoeffs(symmetry s, int n, int m,
			    const V*x, const T**y, Anlm*C) {
    if     (s & PotExp::arot) AddCoeffs<PotExp::spherical  >(n,m,x,y,C);
    else if(s & PotExp::zrot) AddCoeffs<PotExp::cylindrical>(n,m,x,y,C);
    else if(s & PotExp::axes) AddCoeffs<PotExp::triaxial   >(n,m,x,y,C);
    else if(s & PotExp::pint) AddCoeffs<PotExp::reflexion  >(n,m,x,y,C);
    else                      AddCoeffs<PotExp::none       >(n,m,x,y,C);
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
      Cartesian4(A,rd,ct,st,cp,sp);                //   set (ax,ay,az)          
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
      Cartesian(Ac,Rd,Ct,St,Cp,Sp);                //   set (ax,ay,az)          
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
#if 0 // not used
    T Potential(                                   // R: potential: single body 
		symmetry,                          // I: symmetry at run time   
		const   V&);                       // I: position               
#endif
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
#if 0 // not used
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
#endif
  //////////////////////////////////////////////////////////////////////////////
  template<symmetry SYM>
  inline void normalize(Anlm&C, Anlm const&K, scalar G) {
    AUX<SYM>::template Connect<__mulB>(C,K,scalar(0));
    if(G != 1)
      AUX<SYM>::template Connect<__mulX>(C,C,G);
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
    SNprintf(ERR,256,"PotExp: alpha<0.5\n");
    STATE |= 2;
    return;
  }
  if(N<0) {
    SNprintf(ERR,256,"PotExp: nmax must be >= 0\n");
    STATE |= 2;
    return;
  }
  if(L&1) {
    SNprintf(ERR,256,"PotExp: lmax must be even\n");
    STATE |= 2;
    return;
  }
  if(L==0 && s != spherical) {
    SNprintf(WARN,256,"PotExp: assuming spherical symmetry since lmax=0\n");
    STATE |= 1;
  }
  setAL (AL);
  setR0 (R0);
  YlmRec Nlm(l);
  AnlRec Knl(n,l), Anl(n,l);
  SetNlm(Nlm);
  SetKnl(Knl);
  SetAnl(Anl,Knl);
  AUX<none>::Connect<__setB>(Knlm,Anl,Nlm,scalar(0));
}
//------------------------------------------------------------------------------
#define CHECKMISMATCH(FUNC,COEFF)					\
  if(N != (COEFF).nmax()) {						\
    if(L != (COEFF).lmax())						\
      SNprintf(ERR,256,"PotExp::%s(): Anlm have (n,l)_max=(%d,%d), "	\
	       "expected (%d,%d)\n",					\
	       FUNC,(COEFF).nmax(),(COEFF).lmax(),N,L);			\
    else								\
      SNprintf(ERR,256,"PotExp::%s(): Anlm have n_max=%d, "		\
	       "expected %d\n",	 FUNC,(COEFF).nmax(),N);		\
    STATE |= 2;								\
    return;								\
  } else if(L != (COEFF).lmax()) {					\
    SNprintf(ERR,256,"PotExp::%s(): Anlm have l_max=%d, "		\
	     "expected %d\n", FUNC,(COEFF).lmax(),L);			\
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
  CHECKMISMATCH("AddCoeffs",C);
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
template<typename T>                               // T: double or float        
void PotExp::AddCoeffs  (int              m,       // I: # quantities           
			 Anlm            *C,       // O: C_nlm coefficients     
			 int              n,       // I: # bodies               
			 const tupel<3,T>*x,       // I: positions              
			 const T        **y) const // I: quantities per body    
  //                                                                            
  // computes                                                                   
  //                                                                            
  //             1     N                                                        
  // C[k]_nlm = ----  Sum y_ik * Psi_nl(r_i) * Y_lm(theta_i,phi_i)              
  //            r0^2  i=1                                                       
  //                                                                            
  // where the sum includes either all (mass-carrying) bodies if k==0, or       
  // all bodies whose flag contains (at least one of) the bits in mark.         
  //                                                                            
{
  for(int j=0; j!=m; ++j) { CHECKMISMATCH("AddCoeffs",C[j]); }
  setAL(AL);
  setR0(R0);
  CBlock<T>::AddCoeffs(SYM,n,m,x,y,C);
}
//------------------------------------------------------------------------------
template void PotExp::
AddCoeffs<float>(int, Anlm*, int, const tupel<3,float>*, const float**) const;
template void PotExp::
AddCoeffs<double>(int, Anlm*, int,const tupel<3,double>*,const double**) const;
//------------------------------------------------------------------------------
void PotExp::Normalize(Anlm&C, scalar G) const {
  //                                                                            
  // sets  C_nlm *= G * K_nlm                                                   
  //                                                                            
  CHECKMISMATCH("Normalize",C);
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
  CHECKMISMATCH("SetGravity",C);
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
  CHECKMISMATCH("SetGravity",C);
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
  CHECKMISMATCH("SetPotential",C);
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
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class falcON::AnlmIO                                                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#define XDRS static_cast<XDR*>(xdrs)
#define TRY_XDR(x,func,name)						\
  if(!(x)) falcON_THROW("AnlmIO::%s(): XDR operation \"%s\" failed",func,name);
//------------------------------------------------------------------------------
void AnlmIO::open_for_write(const char*file_name) falcON_THROWING
{
  DebugInfo(6,"AnlmIO::open_for_write(\"%s\")\n",file_name);
  // open file and connect xdr stream to it                                     
  if(open != closed)
    falcON_THROW("AnlmIO::open_for_write(): already open");
  file = fopen(file_name,"w");
  if(!file)
    falcON_THROW("cannot open file \"%s\" for writing",file_name);
  if(xdrs == 0) xdrs = new XDR;
  xdrstdio_create( XDRS, file, XDR_ENCODE );
  // create header                                                              
  char header[10] = "anlmfile", *p=header;
  // write identifier
  TRY_XDR(xdr_string(XDRS, &p, 10),"open_for_write","writing header");
  open = writing;
}
//------------------------------------------------------------------------------
void AnlmIO::open_for_read(const char*file_name) falcON_THROWING
{
  DebugInfo(6,"AnlmIO::open_for_read(\"%s\")\n",file_name);
  // open file and connect xdr stream to it                                     
  if(open != closed)
    falcON_THROW("AnlmIO::open_for_read(): already open");
  file = fopen(file_name,"r");
  if(!file)
    falcON_THROW("cannot open file \"%s\" for reading",file_name);
  if(xdrs == 0) xdrs = new XDR;
  xdrstdio_create(XDRS, file, XDR_DECODE);
  // read header                                                                
  char header[10], shead[10] = "anlmfile", *p=header;
  // read identifier
  TRY_XDR(xdr_string(XDRS, &p, 10),"open_for_read","reading header");
  if( strcmp(header,shead) )
    falcON_THROW("file \"%s\" is not a XDR *anlmfile*, "
		 "cannot open for reading",file_name);
  open = reading;
}
//------------------------------------------------------------------------------
void AnlmIO::write(PotExp::symmetry sym, double a, double r,
		   PotExp::Anlm const& A, double t)
  falcON_THROWING
{
  if(open != writing)
    falcON_THROW("AnlmIO::write(): stream not opened for writing");
  int s = sym;
  TRY_XDR(xdr_double(XDRS,&t),"write","writing time");
  TRY_XDR(xdr_double(XDRS,&a),"write","writing alpha");
  TRY_XDR(xdr_double(XDRS,&r),"write","writing scale");
  TRY_XDR(xdr_int(XDRS,&s),"write","writing symmetry");
  TRY_XDR(xdr_int(XDRS,const_cast<int*>(&(A.nmax()))),"write","writing nmax");
  TRY_XDR(xdr_int(XDRS,const_cast<int*>(&(A.lmax()))),"write","writing lmax");
  PotExp::scalar* aN=A.A+(A.nmax()+1)*square(A.lmax()+1);
  for(const PotExp::scalar*a=A.A; a!=aN; ++a)
    TRY_XDR(xdr_double(XDRS,const_cast<double*>(a)),"write","writing Anlm");
}
#undef TRY_XDR
//------------------------------------------------------------------------------
#define TRY_XDR(x,func,name) if(!(x)) return false;
bool AnlmIO::read(PotExp::symmetry&sym, double&a, double&r,
		  PotExp::Anlm&A, double &t)
  falcON_THROWING
{
  if(open != reading)
    falcON_THROW("AnlmIO::read(): stream not opened for reading");
  if(feof(file)) falcON_THROW("AnlmIO::read(): seen end of file\n");
  if(ferror(file)) falcON_THROW("AnlmIO::read(): I/O error\n");
  int n,l,s;
  TRY_XDR(xdr_double(XDRS,&t),"read","reading time");
  TRY_XDR(xdr_double(XDRS,&a),"read","reading alpha");
  TRY_XDR(xdr_double(XDRS,&r),"read","reading scale");
  TRY_XDR(xdr_int(XDRS,&s),"read","reading symmetry");
  TRY_XDR(xdr_int(XDRS,&n),"read","reading nmax");
  TRY_XDR(xdr_int(XDRS,&l),"read","reading lmax");
  sym = PotExp::symmetry(s);
  A.reset(n,l);
  PotExp::scalar* aN=A.A+(n+1)*square(l+1);
  for(PotExp::scalar*a=A.A; a!=aN; ++a)
    TRY_XDR(xdr_double(XDRS,a),"read","reading Anlm");
  return true;
}
#undef TRY_XDR
//------------------------------------------------------------------------------
bool AnlmIO::is_good() const {
  return !feof(file) && !ferror(file);
}
//------------------------------------------------------------------------------
void AnlmIO::close()
{
  if(open) {
    xdr_destroy (XDRS);
    falcON_DEL_O(XDRS);
    fclose(file);
  }
  open = closed;
  xdrs = 0;
  file = 0;
}
#undef XDRS
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
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  void SetYlm(double*a, int l, tupel<3,double> const&x)
  {
    YlmRec Y(l);
    double rd,ct,st,cp,sp;
    Spherical(rd,ct,st,cp,sp,x);
    ::SetYlm<PotExp::none>(Y,ct,st,cp,sp);
    double*end=a+square(l+1);
    for(double*b=&(Y(0,0)); a!=end; ++a,++b) *a = *b;
  }
}
