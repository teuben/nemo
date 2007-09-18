// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   inc/public/PotExp.h                                                 
///                                                                             
/// \brief  potential expansion in to basis functions of Zhao (1996) type.      
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   1994-1996, 2004, 2007                                               
///                                                                             
/// \note                                                                       
/// The code and its implementation in PotExp.cc are supposed to be included in 
/// an implementation of an external acceleration field. For this reason, usage 
/// of embedded falcON features, such as bodies, exceptions, error() has been   
/// avoided.\n                                                                  
/// \note                                                                       
/// For a more customized version, which directly deals with bodies, see file   
/// inc/proper/pexp.h.                                                          
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 1994-2007  Walter Dehnen, Paul McMillan                        
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
#ifndef falcON_included_PotExp_h
#define falcON_included_PotExp_h

#ifndef falcON_included_iostream
#  include <iostream>
#  define falcON_included_iostream
#endif
#ifndef falcON_included_basic_h
#  include <public/basic.h>
#endif

////////////////////////////////////////////////////////////////////////////////
namespace { class PotExpAccess; }                  // forward declaration       
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  class PotExp;                                    // forward declaration       
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::PotExp                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class PotExp {
    friend class PotExpAccess;
  public:
    //--------------------------------------------------------------------------
    // type used internally for real-valued numbers                             
    typedef double scalar;
    //--------------------------------------------------------------------------
    // public enum symmetry                                                     
    //                                                                          
    // The potential expansion code may enforce some type of symmetry.          
    // There are four degrees of symmetry ranging from reflexion to spherical   
    // symmetry. Each degree of symmetry implies all lower degrees.             
    enum symmetry {
      // basic symmetrizing properties for the A(l,m)                           
      none       = 0,                      // all terms used                    
      pint       = 1,                      // A(l,m)=0  IF l or m are odd       
      axes       = 2,                      // A(l,m)=A(l,-m)=[A(l,m)+A_(l,-m)]/2
      zrot       = 4,                      // A(l,m)=0  IF m!=0                 
      arot       = 8,                      // A(l,m)=0  IF l!=0                 
      // degrees of symmetry:                                                   
      reflexion   = pint,                  // reflexion symmetry wrt. origin    
      triaxial    = reflexion | axes,      // reflexion symmetry wrt. x,y,z axes
      cylindrical = triaxial  | zrot,      // axi symmetry                      
      spherical   = cylindrical | arot     // spherical symmetry                
    };
    //--------------------------------------------------------------------------
    static const char* name_of_sym(symmetry s) {
      switch(s) {
      case spherical:   return "spherical";
      case cylindrical: return "cylindrical";
      case triaxial:    return "triaxial";
      case reflexion:   return "reflexion";
      case none:        return "none";
      default:          return "ERRORNEOUS";
      }
    }
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // nested class Anlm                                                      //
    //                                                                        //
    // holds the numbers a_nlm with n=0...N, l=0...L, m=-l...l                //
    //                                                                        //
    // NOTE                                                                   //
    //    Depending on the adopted symmetry, we ASSUME certain elements to be //
    //    zero or trivial copies, but do NOT ENFORCE this. In particular:     //
    //    1. spherical symmetry: only the l=0,m=0 elements are used.          //
    //    2. cylindrical symmetry: only the even l, m=0 elements are used.    //
    //    3. triaxial symmetry: only elements with (l,m) even are considered  //
    //       and it is implicitly assumed that A[-m]=A[m]; only m>0 are used. //
    //    4. reflexion symmetry: all elements with even (l,m) are used.       //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    class Anlm {
      friend class PotExp;
      friend class ::PotExpAccess;
      //------------------------------------------------------------------------
      // type used internally for real-valued numbers                           
      typedef double scalar;
      //------------------------------------------------------------------------
      // data                                                                   
      const int N,L, N1,L1,L1Q;
      scalar   *A;
      //------------------------------------------------------------------------
      // construction & destruction                                             
      Anlm (int n, int l) :
	N(n), L(l),
	N1(N+1), L1(L+1), L1Q(L1*L1), A(falcON_NEW(scalar,N1*L1Q)) {}
    public:
      template<typename POTEXP>
      Anlm (POTEXP const&P) :
	N(P.Nmax()), L(P.Lmax()),
	N1(N+1), L1(L+1), L1Q(L1*L1), A(falcON_NEW(scalar,N1*L1Q)) {}
      Anlm (Anlm const&a) :
	N(a.N), L(a.L),
	N1(N+1), L1(L+1), L1Q(L1*L1), A(falcON_NEW(scalar,N1*L1Q)) {
	memcpy(A,a.A,N1*L1Q*sizeof(scalar));
      }
      ~Anlm() { delete[] A; }
      //------------------------------------------------------------------------
      int const&nmax() const { return N; }
      int const&lmax() const { return L; }
      //------------------------------------------------------------------------
      // element access                                                         
      //    NOTE: that elements ASSUMED to be zero may have any value!          
      scalar      &operator() (int n, int l, int m) {
	return A[n*L1Q+l*(l+1)+m];
      }
      scalar const&operator() (int n, int l, int m) const {
	return A[n*L1Q+l*(l+1)+m];
      }
      //------------------------------------------------------------------------
      Anlm&reset   (symmetry=none);                // A[n,l,m] = 0              
      Anlm&assign  (scalar const&, symmetry=none); // A[n,l,m] = x              
      Anlm&multiply(scalar const&, symmetry=none); // A[n,l,m]*= x              
      Anlm&divide  (scalar const&, symmetry=none); // A[n,l,m]/= x              
      Anlm&copy    (Anlm   const&, symmetry=none); // A[n,l,m] = B[n,l,m]       
      Anlm&add     (Anlm   const&, symmetry=none); // A[n,l,m]+= B[n,l,m]       
      Anlm&subtract(Anlm   const&, symmetry=none); // A[n,l,m]-= B[n,l,m]       
      Anlm&multiply(Anlm   const&, symmetry=none); // A[n,l,m]*= B[n,l,m]       
      scalar dot   (Anlm   const&, symmetry=none)  // dot product               
	const;
      Anlm&addtimes(Anlm   const&,                 // A[n,l,m]+= x*B[n,l,m]     
		    scalar const&, symmetry=none);
      Anlm&subtimes(Anlm   const&,                 // A[n,l,m]-= x*B[n,l,m]     
		    scalar const&, symmetry=none);
      Anlm&apply   (scalar(*f)(scalar)) {          // A[n,l,m] = f(A[n,l,m])    
	scalar* const AN = A+(N1*L1Q);
	for(scalar*a=A; a!=AN; ++a)
	  *a = f(*a);
	return *this;
      }
      //------------------------------------------------------------------------
      // make nice print outs                                                   
      void print(                                  // nice formated print       
		 symmetry,                         // I: symmetry to assume     
		 std::ostream&,                    // I: ostream to write to    
		 int = 6) const;                   //[I: precision]             
      void table_print(                            // print table               
		       symmetry,                   // I: symmetry to assume     
		       std::ostream&,              // I: ostream to write to    
		       int = 6) const;             //[I: precision]             
    }; // class Anlm
  private:
    //--------------------------------------------------------------------------
    // data members                                                             
    const symmetry  SYM;                   // symmetry used by this PotExp      
    const int       N;                     // max N in expansion                
    const int       L;                     // max L in expansion, MUST be even  
    const scalar    AL, R0;                // alpha, scale radius               
    Anlm            Knlm;                  // normalization constants           
    mutable char    STATE;                 // 0: okay, 1: warning 2: error      
    mutable char    ERR[64],WARN[64];      // error & warning message           
    //--------------------------------------------------------------------------
    // const data access                                                        
  protected:
    Anlm        const& K           () const { return Knlm; }
  public:
    symmetry    const& Symmetry    () const { return SYM; }
    scalar      const& alpha       () const { return AL; }
    scalar      const& scale       () const { return R0; }
    int         const& Nmax        () const { return N; }
    int         const& Lmax        () const { return L; }
    bool               is_okay     () const { return STATE == 0; }
    bool               has_warning () const { return STATE & 1; }
    bool               has_error   () const { return STATE & 2; }
    const char*        error_msg   () const { return ERR; }
    const char*        warning_msg () const { STATE &= ~1; return WARN; }
    //--------------------------------------------------------------------------
    PotExp(scalar,                         // parameter alpha                   
	   scalar,                         // scale radius                      
	   int,                            // maximum n in expansion            
	   int,                            // maximum l in expansion            
	   symmetry = reflexion);          // symmetry used                     
    //--------------------------------------------------------------------------
    const char* symmetry_name() const {
      return name_of_sym(SYM);
    }
    //--------------------------------------------------------------------------
    // add C_nlm due to a set of bodies                                         
    // Note     If the pointer to body flags is non-zero AND the last argument  
    //          'k' is non-zero too, only the contributions to C_nlm of bodies  
    //          for which (flag & k) is true are added.                         
    template<typename T>                           // T: double or float        
    void AddCoeffs  (Anlm            &,            // O: C_nlm coefficients     
		     int              ,            // I: # bodies               
		     const T         *,            // I: masses                 
		     const tupel<3,T>*,            // I: positions              
		     const int       *,            // I: flags        see Note  
		     int             = 0) const;   //[I: k]           see Note  
    //--------------------------------------------------------------------------
    // multiply C_nlm with G*A_nl*N_lm                                          
    void Normalize  (Anlm          &,              // I/O: C_nlm coefficients   
		     scalar const  & = 1) const;   //[I: constant of Gravity]   
    //--------------------------------------------------------------------------
    // compute Pot & Acc due to a set of coefficients for a set of bodies       
    // Notes                                                                    
    //       1. If the pointer to body flags is 0, gravity is computed for all  
    //          bodies, otherwise only for those for which (flag & 1) is true.  
    //       2. If the 1st bit of the last argument is set, potentials are      
    //          added, otherwise assigned.                                      
    //          Likewise, if the second bit of that arguement  is set,          
    //          accelerations are added, otherwise assigned. Thus, a value      
    //          of 0 means both potentials and acclerations get assigned        
    template<typename T>                           // T: double or float        
    void SetGravity (Anlm const      &,            // I: C_nlm coefficients     
		     int              ,            // I: # bodies               
		     const tupel<3,T>*,            // I: positions              
		     T               *,            // O: potentials             
		     tupel<3,T>      *,            // O: accelerations          
		     const int       *,            // I: body flags   see Note 1
		     int              ) const;     // I: add?         see Note 2
    //--------------------------------------------------------------------------
    // compute Pot & Acc due to a set of coefficients for a single body         
    template<typename T>                           // T: double or float        
    void SetGravity (Anlm       const&,            // I: C_nlm coefficients     
		     tupel<3,T> const&,            // I: position               
		     T               &,            // O: potential              
		     tupel<3,T>      &,            // O: acceleration           
		     int              ) const;     // I: add?         see Note 2
    //--------------------------------------------------------------------------
    // compute Pot due to a set of coefficients for a set of bodies             
    // Notes                                                                    
    //       1. If the pointer to body flags is 0, gravity is computed for all  
    //          bodies, otherwise only for those for which (flag & 1) is true.  
    //       2. If the 1st bit of the last argument is set, potentials are      
    //          added, otherwise assigned.                                      
    template<typename T>                           // T: double or float        
    void SetPotential(Anlm const      &,           // I: C_nlm coefficients     
		      int              ,           // I: # bodies               
		      const tupel<3,T>*,           // I: positions              
		      T               *,           // O: potentials             
		      const int       *,           // I: body flags   see Note 1
		      int              ) const;    // I: add?         see Note 2
    //--------------------------------------------------------------------------
    // self-gravity:                                                            
    // 1. compute C_nlm from a set of bodies                                    
    // 2. compute gravity from the C_nlm for the same bodies                    
    // 3. return the normalized coefficients                                    
    //                                                                          
    // WARNING NOTE                                                             
    // I wrote this code in the presumption that is will be faster than the     
    // sequence of AddCoeffs(), Normalize(), and SetGravity(), because the      
    // spherical coordinates are computed once only (when the coefficients are  
    // computed) and remembered for when gravity is computed from the coeffs.   
    // However, it turns out that in fact this code is slower! Presumably, this 
    // is because the computation of the spherical coordinates is in fact faster
    // than loading them from memory.                                           
    // Moreover, this routines is not suitable for a body data layout in blocks 
    // (though it could be amended to fit this).                                
    //                                                                          
    // Notes                                                                    
    //       1. If the pointer to body flags is non-zero AND the argument 'mark'
    //          is non-zero too, only the contributions to C_nlm of bodies      
    //          for which (flag & k) is true are added.                         
    //       2. If the pointer to body flags is 0 or the argument 'all' is true,
    //          gravity is computed for all bodies, otherwise only for those    
    //          for which (flag & 1) is true.                                   
    //       3. If the 1st bit of the argument 'add' is set, potentials are     
    //          added, otherwise assigned.                                      
    //          Likewise, if the second bit of that arguement  is set,          
    //          accelerations are added, otherwise assigned. Thus, a value      
    //          of 0 means both potentials and acclerations get assigned        
    template<typename T>                           // T: double or float        
    void SelfGravity(Anlm            &,            // O: normalized C_nlm       
		     int              ,            // I: # bodies               
		     const T         *,            // I: masses                 
		     const tupel<3,T>*,            // I: positions              
		     T               *,            // O: potentials             
		     tupel<3,T>      *,            // O: accelerations          
		     const int       *,            // I: flags        Notes 1&2 
		     int              ,            // I: mark         see Note 1
		     bool             ,            // I: all          see Note 2
		     int              ,            // I: add?         see Note 3
		     scalar          =1) const;    //[I: const of Gravity]      
    //--------------------------------------------------------------------------
  };// class PotExp
} // namespace falcON
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_PotExp_h
