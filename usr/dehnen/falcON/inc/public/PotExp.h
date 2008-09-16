// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   inc/public/PotExp.h                                                 
///                                                                             
/// \brief  potential expansion in to basis functions of Zhao (1996) type.      
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   1996, 2004-2008                                                     
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
// Copyright (C) 1994-2008  Walter Dehnen, Paul McMillan                        
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
#ifndef falcON_included_cstdio
# include <cstdio>
# define falcON_included_cstdio
#endif
#ifndef falcON_included_cstring
# include <cstring>
# define falcON_included_cstring
#endif
#ifdef falcON_MPI
# ifndef falcON_included_mpi_falcON_h
#  include <parallel/mpi_falcON.h>
# endif
#endif

// /////////////////////////////////////////////////////////////////////////////
namespace { class PotExpAccess; }                  // forward declaration       
// /////////////////////////////////////////////////////////////////////////////
namespace falcON {
  class AnlmIO;                                    // forward declaration       
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::PotExp                                                       
  //                                                                            
  /*!
    \brief Expansion of gravitational potential in basis functions, based on
    spherical harmonics (Zhao 1996).                                          

    The potential expansion is based on potential-density pairs
    \f[
    4\pi\,\rho_{nlm}(\mbox{\boldmath$x$}) =
    - \Delta\Psi_{nlm}(\mbox{\boldmath$x$}),
    \f]
    which satisfy the bi-orthogonality relation
    \f[
    \int\mathrm{d}^3\!\mbox{\boldmath$x$}\, \rho_{nlm}(\mbox{\boldmath$x$})\,
    \Psi_{ijk}(\mbox{\boldmath$x$})=
    \frac{\delta_{ni}\,\delta_{lj}\,\delta_{mk}}{N_{nlm}}
    \f]
    with some normalisation constants \f$N_{nlm}\f$.
    Thus, given \f$N\f$ mass points with masses \f$m_i\f$ and positions
    \f$\mathbf{x}_i\f$, we can obtain the coefficients \f$A_{nlm}\f$ for
    the potential expansion as
    \f[
    A_{nlm} = G\,N_{nlm} \sum_{i=1}^N m_i \Psi_{nlm}(\mbox{\boldmath$x$}_i)
    \f]
    with \f$G\f$ Newton's constant of gravity.
    Conversely, the gravitational potential and force are then given as
    \f{eqnarray}
    \Phi(\mbox{\boldmath$x$}) &=& - \sum_{nlm} A_{nlm}\,
    \Psi_{nlm}(\mbox{\boldmath$x$}) \\
    \mbox{\boldmath$F$}(\mbox{\boldmath$x$}) &=& \sum_{nlm} A_{nlm}\,
    \mbox{\boldmath$\nabla$}\Psi_{nlm}(\mbox{\boldmath$x$}).
    \f}
    The lowest order basis function has functional form (with
    \f$ r=|\mbox{\boldmath$x$}|\f$)
    \f[
    \Psi_{000} = \frac{1}{(r^{1/\alpha} + r_s^{1/\alpha})^{\alpha}}.
    \f]
    Thus, for \f$\alpha=1/2\f$ we get the Clutton-Brock (1972) expansion based
    on the Plummer sphere as lowest order, while for \f$\alpha=1\f$ we get
    the Hernquist & Ostriker (1992) expansion based on a Hernquist sphere as
    lowest order.
    \note
    The code here is not explicitly using bodies. For a frontend that
    does directly take bodies arguments, see classes PotentialExpansion
    and NPotExpansion (both only with the proprietary version).
  */
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class PotExp {
    friend class PotExpAccess;
  public:
    //--------------------------------------------------------------------------
    /// type used internally for real-valued numbers                            
    typedef double scalar;
    //--------------------------------------------------------------------------
    /// type to encode symmetry properties                                      
    //                                                                          
    /// The potential expansion code may enforce some type of symmetry.         
    /// There are four degrees of symmetry ranging from reflexion to spherical  
    /// symmetry. Each degree of symmetry implies all lower degrees.            
    enum symmetry {
      // basic symmetrizing properties for the A(l,m)                           
      none       = 0,                    ///< all terms used                    
      pint       = 1,                    ///< A(l,m)=0  IF l or m are odd       
      axes       = 2,                    ///< A(l,m)=A(l,-m)=[A(l,m)+A_(l,-m)]/2
      zrot       = 4,                    ///< A(l,m)=0  IF m!=0                 
      arot       = 8,                    ///< A(l,m)=0  IF l!=0                 
      // degrees of symmetry:                                                   
      reflexion   = pint,                ///< reflexion symmetry wrt. origin    
      triaxial    = reflexion | axes,    ///< reflexion symmetry wrt. x,y,z axes
      cylindrical = triaxial  | zrot,    ///< axial symmetry                    
      spherical   = cylindrical | arot   ///< spherical symmetry                
    };
    //--------------------------------------------------------------------------
    /// given a symmetry, return a name
    static const char* name_of_sym(symmetry s) {
      switch(s) {
      case spherical:   return "spherical";
      case cylindrical: return "cylindrical";
      case triaxial:    return "triaxial";
      case reflexion:   return "reflexion";
      case none:        return "none";
      default: {
	bool empty = true;
	std::string name("Errorneous: ");
	if(s & pint) { name += empty? "pint":"|pint"; empty=false; }
	if(s & axes) { name += empty? "axes":"|axes"; empty=false; }
	if(s & zrot) { name += empty? "zrot":"|zrot"; empty=false; }
	if(s & arot) { name += empty? "arot":"|arot"; empty=false; }
	return name.c_str();
      }
      }
    }
    // /////////////////////////////////////////////////////////////////////////
    //                                                                          
    // nested class Anlm                                                        
    //                                                                          
    /// holds the coefficients a_nlm with n in [0,N], l in [0,L], m in [-l,l].  
    /// Depending on the adopted symmetry, we ASSUME certain elements to be zero
    /// or trivial copies, but do NOT ENFORCE this. In particular:            \n
    ///   1. spherical symmetry: only the l=0,m=0 elements are used.          \n
    ///   2. cylindrical symmetry: only the even l, m=0 elements are used.    \n
    ///   3. triaxial symmetry: only elements with (l,m) even are considered    
    ///      and it is implicitly assumed that A[-m]=A[m]; only m>0 are used. \n
    ///   4. reflexion symmetry: all elements with even (l,m) are used.       \n
    //                                                                          
    // /////////////////////////////////////////////////////////////////////////
    class Anlm {
      friend class PotExp;
      friend class ::PotExpAccess;
      friend class AnlmIO;
      //------------------------------------------------------------------------
      // data                                                                   
      const int N;            ///< maximum value for n
      const int L;            ///< maximum value for l
      const int N1,L1,L1Q;    ///< nmax+1, lmax+1, (lmax+1)^2
      scalar   *A;            ///< pointer to coefficient data
      //------------------------------------------------------------------------
    public:
      /// like destruction followed by construction
      /// \param nmax maximum for n
      /// \param lmax maximum for l
      void reset(int nmax, int lmax) {
	if(nmax!=N || lmax!=L) {
	  falcON_DEL_A(A);
	  *(const_cast<int*>(&N)) = nmax;    // dirty trick to change const data
	  *(const_cast<int*>(&L)) = lmax;
	  *(const_cast<int*>(&N1)) = N+1;
	  *(const_cast<int*>(&L1)) = L+1;
	  *(const_cast<int*>(&L1Q)) = L1*L1;
	  A = falcON_NEW(scalar,N1*L1Q);
	}
      }
      //------------------------------------------------------------------------
      /// construction from given n_max & l_max, includes default constructor
      Anlm (int nmax=0, int lmax=0) :
	N(nmax), L(lmax),
	N1(N+1), L1(L+1), L1Q(L1*L1), A(falcON_NEW(scalar,N1*L1Q)) {}
      /// construction form a potential expansion
      template<typename POTEXP>
      explicit Anlm (POTEXP const&P) :
	N(P.Nmax()), L(P.Lmax()),
	N1(N+1), L1(L+1), L1Q(L1*L1), A(falcON_NEW(scalar,N1*L1Q)) {}
      /// copy constructor
      explicit Anlm (Anlm const&a) :
	N(a.N), L(a.L),
	N1(N+1), L1(L+1), L1Q(L1*L1), A(falcON_NEW(scalar,N1*L1Q)) {
	memcpy(A,a.A,N1*L1Q*sizeof(scalar));
      }
      /// destructor: delete data array
      ~Anlm() { falcON_DEL_A(A); }
      //------------------------------------------------------------------------
      /// return maximum possible n
      int const&nmax() const { return N; }
      /// return maximum possible l
      int const&lmax() const { return L; }
      //------------------------------------------------------------------------
      /// non-const element access, given n,l,m
      /// \note: that elements ASSUMED to be zero may have any value!           
      scalar      &operator() (int n, int l, int m) {
	return A[n*L1Q+l*(l+1)+m];
      }
      /// const element access, given n,l,m
      /// \note: that elements ASSUMED to be zero may have any value!           
      scalar const&operator() (int n, int l, int m) const {
	return A[n*L1Q+l*(l+1)+m];
      }
      //------------------------------------------------------------------------
      /// assign to scalar: A[n,l,m] = x
      /// \param x scalar to assign to
      /// \param s if not none: only relevant coefficients are dealt with
      Anlm&assign(scalar x, symmetry s=none);
      /// reset: A[n,l,m] = 0
      /// \param s if not none: only relevant coefficients are dealt with
      Anlm&reset(symmetry s=none) { return assign(scalar(0),s); }
      /// negate: A[n,l,m] =-A[n,l,m]
      /// \param s if not none: only relevant coefficients are dealt with
      Anlm&negate(symmetry s=none);
      /// multiply by scalar: A[n,l,m]*= x
      /// \param x scalar to multiply with
      /// \param s if not none: only relevant coefficients are dealt with
      Anlm&multiply(scalar x, symmetry s=none);
      /// divide by scalar: A[n,l,m]/= x
      /// \param x scalar to divide by
      /// \param s if not none: only relevant coefficients are dealt with
      Anlm&divide(scalar x, symmetry s=none) { return multiply(1/x,s); }
      /// assign element wise: A[n,l,m] = B[n,l,m]
      /// \param B Anlm to copy
      /// \param s if not none: only relevant coefficients are dealt with
      Anlm&copy(Anlm const&B, symmetry s=none);
      /// add element wise: A[n,l,m] += B[n,l,m]
      /// \param B Anlm to add
      /// \param s if not none: only relevant coefficients are dealt with
      Anlm&add(Anlm const&B, symmetry s=none);
      /// subtract element wise: A[n,l,m] -= B[n,l,m]
      /// \param B Anlm to subtract
      /// \param s if not none: only relevant coefficients are dealt with
      Anlm&subtract(Anlm const&B, symmetry s=none);
      /// multiply element wise: A[n,l,m] *= B[n,l,m]
      /// \param B Anlm to multiply element-wise with
      /// \param s if not none: only relevant coefficients are dealt with
      Anlm&multiply(Anlm const&B, symmetry s=none);
      /// dot product: return Sum_nlm A[n,l,m]*B[n,l,m]
      /// \param B Anlm to perform dot product with
      /// \param s if not none: only relevant coefficients are dealt with
      scalar dot(Anlm const&B, symmetry s=none) const;
      /// add element wise times scalar: A[n,l,m] += x*B[n,l,m]
      /// \param B Anlm operand
      /// \param x scalar operand
      /// \param s if not none: only relevant coefficients are dealt with
      Anlm&addtimes(Anlm const&B, scalar x, symmetry s=none);
      /// subtract element wise times scalar: A[n,l,m] -= x*B[n,l,m]
      /// \param B Anlm operand
      /// \param x scalar operand
      /// \param s if not none: only relevant coefficients are dealt with
      Anlm&subtimes(Anlm const&B, scalar x, symmetry s=none);
      /// general unary operation: A[n,l,m] = f(A[n,l,m])
      /// \param f function defining operation
      /// \param s if not none: only relevant coefficients are dealt with
      Anlm&unary(scalar(*f)(scalar), symmetry s=none);
      /// general binary operation with Anlm: A[n,l,m] = f(A[n,l,m],B[n,l,m])
      /// \param f function defining operation
      /// \param B Anlm operand
      /// \param s if not none: only relevant coefficients are dealt with
      Anlm&binary(scalar(*f)(scalar,scalar), Anlm const&B, symmetry s=none);
      /// general binary operation with scalar: A[n,l,m] = f(A[n,l,m],x)
      /// \param f function defining operation
      /// \param x scalar operand
      /// \param s if not none: only relevant coefficients are dealt with
      Anlm&binary(scalar(*f)(scalar,scalar), scalar x, symmetry s=none);
      /// general tertiary operation: A[n,l,m] = f(A[n,l,m],B[n,l,m],x)
      /// \param f function defining operation
      /// \param B Anlm operand
      /// \param x scalar operand
      /// \param s if not none: only relevant coefficients are dealt with
      Anlm&tertiary(scalar(*f)(scalar,scalar,scalar), Anlm const&B, scalar x,
		    symmetry s=none);
      /// apply element wise operation: A[n,l,m] = f(A[n,l,m])
      /// \param f function to apply to each element
      /// \note this is identical to unary(f,none)
      Anlm&apply(scalar(*f)(scalar)) {
	scalar* const AN = A+(N1*L1Q);
	for(scalar*a=A; a!=AN; ++a) *a = f(*a);
	return *this;
      }
      //------------------------------------------------------------------------
#ifdef falcON_MPI
      /// set to sum of A_nlm over all MPI processes
      Anlm&global_sum(Anlm const&A, const MPI::Communicator*C);
#endif
      //------------------------------------------------------------------------
      /// make formated nice print
      /// \param sym symmetry to assume (suppress printing of ignored terms)
      /// \param out where to print to
      /// \param prec precision for printing scalars
      void print(symmetry sym, std::ostream&out, int prec = 6) const;
      /// print coefficients as table
      /// \param sym symmetry to assume (suppress printing of ignored terms)
      /// \param out where to print to
      /// \param prec precision for printing scalars
      void table_print(symmetry sym, std::ostream&out, int prec = 6) const;
    };// sub-class Anlm
  private:
    //--------------------------------------------------------------------------
    // data members                                                             
    const symmetry  SYM;                   ///< symmetry used by this PotExp    
    const int       N;                     ///< max n in expansion              
    const int       L;                     ///< max l in expansion              
    const scalar    AL, R0;                ///< alpha, scale radius             
    Anlm            Knlm;                  ///< normalisation constants         
    mutable char    STATE;                 ///< 0: okay, 1: warning 2: error    
    mutable char    ERR[256],WARN[256];    ///< error & warning message         
    //--------------------------------------------------------------------------
    // const data access                                                        
  protected:
    Anlm        const& K           () const { return Knlm; }
  public:
    /// return assumed symmetry
    symmetry    const& Symmetry    () const { return SYM; }
    /// return parameter alpha
    scalar      const& alpha       () const { return AL; }
    /// return scale radius
    scalar      const& scale       () const { return R0; }
    /// return maximum n
    int         const& Nmax        () const { return N; }
    /// return maximum n
    int         const& Lmax        () const { return L; }
    /// has an error or warning occured?
    bool               is_okay     () const { return STATE == 0; }
    /// has a warning occured?
    bool               has_warning () const { return STATE & 1; }
    /// has an error occured?
    bool               has_error   () const { return STATE & 2; }
    /// error message, if any
    const char*        error_msg   () const { return ERR; }
    /// warning message, if any
    const char*        warning_msg () const { STATE &= ~1; return WARN; }
    //--------------------------------------------------------------------------
    /// constructor
    /// \param alpha parameter \f$\alpha\f$
    /// \param scale scale radius \f$r_s\f$
    /// \param nmax maximum n in expansion
    /// \param lmax maximum l in expansion
    /// \param sym symmetry to be assumed, default: reflexion symmetry
    PotExp(scalar alpha, scalar scale,
	   int nmax, int lmax, symmetry sym=reflexion);
    //--------------------------------------------------------------------------
    /// return name of symmetry used
    const char* symmetry_name() const {
      return name_of_sym(SYM);
    }
    //--------------------------------------------------------------------------
    /*!
    \brief add coefficients due to a set of bodies.

    This amounts to the operation
    \f[
      A_{nlm} += \sum_i m_i \Psi_{nlm}(\mbox{\boldmath$x$}_i).
    \f]
    \note The \f$A_{nlm}\f$ are not yet normalised, this must be done
          via PotExp::Normalize().
    \param T type of scalar used for masses & positions (float or double)
    \param A coefficients \f$A_{nlm}\f$ to add to
    \param n number of bodies
    \param m array with body masses
    \param x array with body positions
    \param f array with body flags (can be null)
    \param k flag (see note below)
    \note If flags are provided and the last argument is non-zero, only
          the contributions of bodies with (f[i] & k) is true are added.
    */
    template<typename T>
    void AddCoeffs(Anlm&A, int n, const T*m, const tupel<3,T>*x,
	           const int*f, int k=0) const;
    //--------------------------------------------------------------------------
    /*!
    \brief normalise the coefficients so that gravity can be computed.

    This amounts to \f$ A_{nlm} \to G\,A_{nlm} N_{nlm}\f$.
     \param A coefficients \f$A_{nlm}\f$ to be normalised
     \param G Newton's gravitation constant
    */
    void Normalize(Anlm&A, scalar G=1) const;
    //--------------------------------------------------------------------------
    /*!
    \brief compute potential and accelerations due to a set of coefficients 
           for a set of bodies.
	       
    This amounts to the operations
    \f{eqnarray*}
    \Phi(\mbox{\boldmath$x$}) &=& - \sum_{nlm} A_{nlm}\,
    \Psi_{nlm}(\mbox{\boldmath$x$}) \\
    \mbox{\boldmath$F$}(\mbox{\boldmath$x$}) &=& \sum_{nlm} A_{nlm}\,
    \mbox{\boldmath$\nabla$}\Psi_{nlm}(\mbox{\boldmath$x$}).
    \f}
    \note The coefficients must have been normalised via PotExp::Normalize().
    \param T type of scalar used for positions, potential, accelerations (float or double)
    \param A coefficients \f$A_{nlm}\f$ to be used
    \param n number of bodies
    \param x array with body positions
    \param P array with body potentials to be assigned/added to
    \param F array with body accelerations to be assigned/added to
    \param f array with body flags (can be null)
    \param add flag (see note below)
    \note If flags are provided gravity is computed only for bodies for 
          which f[i]&1 is true.
    \note If the 1st bit of the last argument is set, potentials are      
          added, otherwise assigned.\n
	  Likewise, if the second bit of that arguement is set,
	  accelerations are added, otherwise assigned. Thus, a value
	  of 0 means both potentials and acclerations get assigned.
    */
    template<typename T>
    void SetGravity (Anlm const&A, int n, const tupel<3,T>*x, T*P, tupel<3,T>*F,
		     const int*f, int add) const;
    //--------------------------------------------------------------------------
    /// compute potential and accelerations due to a set of coefficients 
    /// for a single body.
    template<typename T>
    void SetGravity (Anlm const&A, tupel<3,T> const&x, T&P, tupel<3,T>&F,
		     int add) const;
    //--------------------------------------------------------------------------
    /*!
    \brief compute potential due to a set of coefficients for a set of bodies.
	       
    This amounts to the operations
    \f[
    \Phi(\mbox{\boldmath$x$}) = - \sum_{nlm} A_{nlm}\,
    \Psi_{nlm}(\mbox{\boldmath$x$}) \\
    \f]
    \note The coefficients must have been normalised via PotExp::Normalize().
    \param T type of scalar used for positions, potential (float or double)
    \param A coefficients \f$A_{nlm}\f$ to be used
    \param n number of bodies
    \param x array with body positions
    \param P array with body potentials to be assigned/added to
    \param f array with body flags (can be null)
    \param add flag (see note below)
    \note If flags are provided gravity is computed only for bodies for 
          which f[i]&1 is true.
    \note If the 1st bit of the last argument is set, potentials are      
          added, otherwise assigned.
    */
    template<typename T>
    void SetPotential(Anlm const&A, int n, const tupel<3,T>*x, T*P,
		      const int*f, int add) const;
    //--------------------------------------------------------------------------
    /// self-gravity:\n
    /// 1 compute \f$A_{nlm}\f$ from a set of bodies\n
    /// 2 compute gravity from the \f$A_{nlm}\f$ for the same bodies\n
    /// 3 return the normalized coefficients\n
    ///                                                                         
    /// \warning
    /// I wrote this code in the presumption that it will be faster than the
    /// sequence of AddCoeffs(), Normalize(), and SetGravity(), because the
    /// spherical coordinates are computed once only (when the coefficients are
    /// computed) and remembered for when gravity is computed from the coeffs.
    /// However, it turns out that in fact this code is slower! Presumably, this
    /// is because the computation of the spherical coordinates is in fact
    /// faster than loading them from memory.\n
    /// Moreover, this routine is not suitable for a body data layout in blocks
    /// (though it could be amended to fit this).
    ///
    /// \param T type of scalar used for positions, potential (float or double)
    /// \param A coefficients \f$A_{nlm}\f$; on output: normalized coeffs
    /// \param n number of bodies
    /// \param m array with body masses
    /// \param x array with body positions
    /// \param P array with body potentials to be assigned/added to
    /// \param F array with body accelerations to be assigned/added to
    /// \param f array with body flags (can be null)
    /// \param k flag (see note below)
    /// \param all compute gravity for all (or only those with f[i]&1 == true)?
    /// \param add flag (see note below)
    /// \param G Newton's gravitation constant
    /// \note If flags are provided and 'k' is non-zero, only the contributions 
    ///       of bodies with (f[i] & k) is true are added to the \f$A_{nlm}\f$.
    /// \note If the 1st bit of the last argument is set, potentials are      
    ///       added, otherwise assigned.\n
    ///       Likewise, if the second bit of that arguement  is set,
    ///       accelerations are added, otherwise assigned. Thus, a value
    ///       of 0 means both potentials and acclerations get assigned.
    template<typename T>
    void SelfGravity(Anlm&A, int n, const T*m, const tupel<3,T>*x,
		     T*P, tupel<3,T>*F, const int*f,
		     int k, bool all, int add, scalar G=1) const;
    //--------------------------------------------------------------------------
  protected:
    template<typename T>
    void AddCoeffs(int m, Anlm*A, int n, const tupel<3,T>*x, const T**y) const;
  };// class PotExp
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class AnlmIO                                                               
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class AnlmIO {
    AnlmIO (AnlmIO const&);
    AnlmIO&operator=(AnlmIO const&);
    //--------------------------------------------------------------------------
    enum { closed=0, writing=1, reading=2 };
    int  open;  ///< state
    void*xdrs;  ///< xdr stream
    FILE*file;  ///< C file stream
    //--------------------------------------------------------------------------
  protected:
    AnlmIO() : open(0), xdrs(0), file(0) {}
   ~AnlmIO() { close(); }
    //--------------------------------------------------------------------------
    void open_for_read (const char*) falcON_THROWING;
    void open_for_write(const char*) falcON_THROWING;
    void open_for_read (std::string const&f) falcON_THROWING
    { open_for_read(f.c_str()); }
    void open_for_write(std::string const&f) falcON_THROWING
    { open_for_write(f.c_str()); }
    //--------------------------------------------------------------------------
    /// \param sym symmetry
    /// \param al  alpha
    /// \param rs  scale
    /// \param A   coeffs
    /// \param t   time
    void write(PotExp::symmetry sym, double al, double rs,
	       PotExp::Anlm const& A, double t) falcON_THROWING;
    /// \param sym symmetry
    /// \param al  alpha
    /// \param rs  scale
    /// \param A   coeffs
    /// \param t   time
    /// \return was read operation successful?
    bool read (PotExp::symmetry&sym, double&al, double&rs,
	       PotExp::Anlm&A, double&t) falcON_THROWING;
  public:
    bool is_open() const { return open!=closed; }
    bool is_good() const;
    void close();
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class AnlmI                                                                
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class AnlmI : public AnlmIO {
    AnlmI (AnlmI const&);
    AnlmI&operator=(AnlmI const&);
    //--------------------------------------------------------------------------
  public:
    AnlmI(const char*f) falcON_THROWING { open_for_read(f); }
    AnlmI(std::string const&f) falcON_THROWING { open_for_read(f.c_str()); }
    void open(const char*f) {
      if(is_open()) close();
      open_for_read(f);
    }
    void open(std::string const&f) {
      if(is_open()) close();
      open_for_read(f);
    }
    AnlmIO::read;
    operator bool() const { return is_open(); }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class AnlmO                                                                
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class AnlmO : public AnlmIO {
    AnlmO (AnlmO const&);
    AnlmO&operator=(AnlmO const&);
    //--------------------------------------------------------------------------
  public:
    AnlmO(const char*f) falcON_THROWING { open_for_write(f); }
    AnlmO(std::string const&f) falcON_THROWING { open_for_write(f.c_str()); }
    void open(const char*f) {
      if(is_open()) close();
      open_for_write(f);
    }
    void open(std::string const&f) {
      if(is_open()) close();
      open_for_write(f);
    }
    AnlmIO::write;
    operator bool() const { return is_open(); }
  };
} // namespace falcON
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_PotExp_h
