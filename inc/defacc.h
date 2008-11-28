/*
 *******************************************************************************
 *                                                                              
 * defacc.h                                                                     
 *                                                                              
 * Copyright Walter Dehnen, 2003-2008                                           
 * e-mail:   walter.dehnen@astro.le.ac.uk                                       
 * address:  Department of Physics and Astronomy, University of Leicester       
 *           University Road, Leicester LE1 7RH, United Kingdom                 
 *                                                                              
 *******************************************************************************
 *                                                                              
 *  declaration of functions that any implementation of an external,            
 *  dynamically loadable acceleration must define. See also acceleration.h      
 *                                                                              
 *******************************************************************************
 *                                                                              
 *  to ensure that your implementation of external acceleration is consistent,  
 *  include this header file into your C or C++ source code.                    
 *                                                                              
 *******************************************************************************
 * version 0.0  17/06/2004  WD                                                  
 * version 1.0  22/06/2004  WD  support for acc & pot via C++ & macros          
 * version 1.1  25/06/2004  WD  somewhat improved documentation                 
 * version 2.0  18/08/2004  WD  seperate accelertion and potential              
 *                              no array allocation in case of adding           
 *                              StaticSphericalModel, _NO_AUX_DEFACC            
 * version 3.0  24/08/2004  WD  no re-initialisation but independent acc fields;
 *                              allows superpositions using GrowCombined        
 * version 3.1  17/09/2004  WD  somewhat improved documentation                 
 * version 3.2  12/11/2004  WD  catching std::bad_alloc and report error        
 * version 3.3  04/05/2007  WD  added global scope resolution operators         
 * version 3.4  12/06/2008  WD  #including stdinc.h                            
 * version 3.5  11/09/2008  WD  no inclusion of NEMO header files               
 *                                                                              
 *******************************************************************************
 *                                                                              
 * Contents:                                                                    
 *                                                                              
 * 1 Declaration of C-linkable routines required by implementations of external 
 *   acceleration fields and, optionally, external potential.                   
 *                                                                              
 * C++ only:                                                                    
 *                                                                              
 * 2 Definition of helper templates for manipulating arrays or scalar like      
 *   NDIM dimensional vectores.                                                 
 *                                                                              
 * 3 Definitions of templates and macros that ease the implementation of the    
 *   C-linkable routines in C++.                                                
 *                                                                              
 *******************************************************************************
 */
#ifndef _defacc_h
#define _defacc_h

#include <acceleration.h>                   /* for definition of acc_pter     */

#ifdef __cplusplus
# include <cstring>
  extern "C" {
#else
# include <string.h>
#endif

#if(0)
#include <stdinc.h>                         /* for error, warning, dprintf    */
#endif

/*
 *******************************************************************************
 *                                                                              
 * 1 Declaration of C-linkable routines required by implementations of external 
 *   acceleration fields and, optionally, external potential.                   
 *                                                                              
 *******************************************************************************
 * 
 * functions to be defined in #including file.                                  
 *                                                                              
 *  Notes:                                                                      
 *  1 the 4th argument to iniacceleration() is boolean and returns whether      
 *    masses are required as input for acceleration().                          
 *  2 the 5th argument to iniacceleration() is boolean and returns whether      
 *    velocities are required as input for acceleration().                      
 *    Velocities may be used to compute friction forces, such as the drag       
 *    a gaseous disk is generating on stars crossing it.                        
 *  3 arrays are passed to acceleration() as pointer to void. They must be      
 *    either all of type float or all of type double as indicated by the last   
 *    argument being 'f' or 'd', respectively.                                  
 *  4 arrays of vector quantities are in the order x0,y0,z0, x1,y1,z1, ...      
 *  5 if the pointer to flags is NULL, all bodies are supposed to be active,    
 *    otherwise only those for which (f[i] & 1) is true.                        
 *  6 the argument "indicator" of acceleration() indicates whether the          
 *    accelerations and potential shall be assigned or added.                   
 *    If bit 0 is set, the potential    is added, otherwise assigned,           
 *    If bit 1 is set, the acceleration is added, otherwise assigned.           
 *    So, 0 means both are assigned.                                            
 *                                                                              
 */

void iniacceleration(                 /* return: void                         */
		     const double*,   /* input:  array with parameters        */
		     int,             /* input:  number of parameters         */
		     const char*,     /* input:  data file name (may be 0     */
		     acc_pter*,       /* output: pter to acceleration()       */
		     bool*,           /* output: acceleration() needs masses? */
		     bool*);          /* output: acceleration() needs vel's?  */

/*
 * NOTE 1
 * the routine pointed to by the pointer returned from iniacceleration() must   
 * also be defined in the same file that defines iniacceleration() so that it   
 * is loaded at the same time.
 */

/*
 * Declaration of NEMO functions that may be used in the implementation of the  
 * above. They are resolved in the NEMO library (this way, we avoid #including  
 * nemo.h) 
 *                                                                              
 * NOTE (added with version 3.5)                                                
 * Nemo defines many functions to take pointer arguments when a const pointer   
 * should have been used, for instance void error(char*,...).                   
 * With g++ 4.3.1 this causes code like                                         
 *     error("this is an error message");                                       
 * to generate the warning message                                              
 *     warning: deprecated conversion from string constant to ‘char*’          
 * indicating that this will become an error in future versions. This implies   
 * that we cannot simply include nemo header files into C++ code as easily as   
 * we used to.                                                                 
 * Instead, here we use a trick: declaring these routines taking const pointers.
 * This works, as both have the same C-linkage. However, the compiler will      
 * complain if it sees this and the nemo definitions simultaneously.            
 *******************************************************************************
 */

#ifndef _stdinc_h
  extern void warning(const char *, ...);   /* differs from stdinc.h */
  extern void error  (const char *, ...);   /* differs from stdinc.h */
  extern bool nemo_debug(int);
  typedef int (*dprintf_pter)(int, const char*, ...);
  extern dprintf_pter get_dprintf(const char*, int);
# define nemo_dprintf  get_dprintf(__FILE__,__LINE__)
#endif

#ifdef POT_DEF
/*
 * functions that may be defined in #including file.                            
 *
 * NOTE 2
 * here we use the old NEMO style of only one routine per potential type, not   
 * allowing for superpositions of, say, several Miyamoto-Nagai disks.           
 */
void inipotential    (                /* return: void                         */
		      const int*,     /* input:  number of parameters         */
		      const double*,  /* input:  array with parameters        */
		      const char*);   /* input:  data file name (may be NULL) */
void potential_double(const int*,     /* input:  number of dimensions         */
		      const double*,  /* input:  position                     */
		      double*,        /* output: acceleration                 */
		      double*,        /* output: potential                    */
		      const double*); /* input:  time                         */
void potential_float (const int*,     /* input:  number of dimensions         */
		      const float*,   /* input:  position                     */
		      float*,         /* output: acceleration                 */
		      float*,         /* output: potential                    */
		      const float*);  /* input:  time                         */
#endif

#ifdef __cplusplus
}
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// 2 Definition of helper templates for manipulating arrays or scalar like    //
//   NDIM dimensional vectores.                                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace {
  //
  // auxiliary template classes used by the helper functions below
  //
  template<int N, int I=0> struct Dims {
    template<typename X>
    static X square(X x) { return x*x; }
    typedef Dims<N,I+1> D;
    template<typename X>
    static X norm(const X*a) { return a[I]*a[I] + D::norm(a); }
    template<typename X, typename Y>
    static void set (X*a, Y x) { a[I] =x; D::set(a,x); }
    template<typename X, typename Y>
    static void ass (X*a, const Y*b) { a[I]=b[I]; D::ass(a,b); }
    template<typename X, typename Y, typename Z>
    static void asssum (X*a, const Y*b, const Z*c) { a[I] =b[I]+c[I]; 
						     D::asssum(a,b,c); }
    template<typename X, typename Y, typename Z>
    static void assdif (X*a, const Y*b, const Z*c) { a[I] =b[I]-c[I];
						     D::assdif(a,b,c); }
    template<typename X, typename Y>
    static void add (X*a, const Y*b) { a[I]+=b[I]; D::add(a,b); }
    template<typename X, typename Y>
    static void sub (X*a, const Y*b) { a[I]-=b[I]; D::sub(a,b); }
    template<typename X, typename Y>
    static void asstimes(X*a,const Y*b,Y x) { a[I] =x*b[I]; D::asstimes(a,b,x);}
    template<typename X, typename Y>
    static void addtimes(X*a,const Y*b,Y x) { a[I]+=x*b[I]; D::addtimes(a,b,x);}
    template<typename X, typename Y>
    static void mul (X*a, Y x) { a[I]*=x; D::mul(a,x); }
    template<typename X, typename Y>
    static X    dot (const X*a, const Y*b) { return a[I]*b[I] + D::dot(a,b); }
    template<typename X, typename Y>
    static X    distsq(const X*a, const Y*b) { return square(a[I]-b[I])
						     + D::distsq(a,b); }
  };
  template<int I> struct Dims<I,I> {
    template<typename X>
    static X square(X x) { return x*x; }
    template<typename X>
    static X norm(const X*a) { return a[I]*a[I]; }
    template<typename X, typename Y>
    static void set (X*a, Y x) { a[I] =x; }
    template<typename X, typename Y>
    static void ass (X*a, const Y*b) { a[I] =b[I]; }
    template<typename X, typename Y, typename Z>
    static void asssum (X*a, const Y*b, const Z*c) { a[I] =b[I]+c[I]; }
    template<typename X, typename Y, typename Z>
    static void assdif (X*a, const Y*b, const Z*c) { a[I] =b[I]-c[I]; }
    template<typename X, typename Y>
    static void add (X*a, const Y*b) { a[I]+=b[I]; }
    template<typename X, typename Y>
    static void sub (X*a, const Y*b) { a[I]-=b[I]; }
    template<typename X, typename Y>
    static void asstimes(X*a,const Y*b,Y x) { a[I] =x*b[I]; }
    template<typename X, typename Y>
    static void addtimes(X*a,const Y*b,Y x) { a[I]+=x*b[I]; }
    template<typename X, typename Y>
    static void mul (X*a, Y x) { a[I]*=x; }
    template<typename X, typename Y>
    static X    dot (const X*a, const Y*b) { return a[I]*b[I]; }
    template<typename X, typename Y>
    static X    distsq(const X*a, const Y*b) { return square(a[I]-b[I]); }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // helper template functions useful for manipulating pointers to scalar     //
  // as N dimensional vectors.                                                //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  // set vector elements equal to scalar:  Forall_i: a_i = x                    
  template<int N, typename X, typename Y>
  inline void v_set(X*a, Y x) { Dims<N-1>::set(a,x); }
  //----------------------------------------------------------------------------
  // set vector equal to vector: Forall_i: a_i = b_i                            
  template<int N, typename X, typename Y>
  inline void v_ass(X*a, const Y*b) { Dims<N-1>::ass(a,b); }
  //----------------------------------------------------------------------------
  // set vector equal to sum of two vectors: Forall_i a_i = b_i + c_i           
  template<int N, typename X, typename Y, typename Z>
  inline void v_asssum(X*a, const Y*b, const Z*c) { Dims<N-1>::asssum(a,b,c); }
  //----------------------------------------------------------------------------
  // set vector equal to difference between two vectors: Forall_i a_i = b_i-c_i 
  template<int N, typename X, typename Y, typename Z>
  inline void v_assdif(X*a, const Y*b, const Z*c) { Dims<N-1>::assdif(a,b,c); }
  //----------------------------------------------------------------------------
  // add vector to vector: Forall_i a_i = a_i + b_i                             
  template<int N, typename X, typename Y>
  inline void v_add(X*a, const Y*b) { Dims<N-1>::add(a,b); }
  //----------------------------------------------------------------------------
  // subtract vector from vector: Forall_i a_i = a_i - b_i                      
  template<int N, typename X, typename Y>
  inline void v_sub(X*a, const Y*b) { Dims<N-1>::sub(a,b); }
  //----------------------------------------------------------------------------
  // set vector equal to vector times scalar: Forall_i: a_i = b_i * x           
  template<int N, typename X, typename Y>
  inline void v_asstimes(X*a, const Y*b,Y x) { Dims<N-1>::asstimes(a,b,x); }
  //----------------------------------------------------------------------------
  // add vector times scalar to vector: Forall_i: a_i = a_i + b_i * x           
  template<int N, typename X, typename Y>
  inline void v_addtimes(X*a, const Y*b,Y x) { Dims<N-1>::addtimes(a,b,x); }
  //----------------------------------------------------------------------------
  // multiply vector by scalar: Forall_i a_i = a_i * x                          
  template<int N, typename X, typename Y>
  inline void v_mul(X*a, Y x) { Dims<N-1>::mul(a,x); }
  //----------------------------------------------------------------------------
  // return norm(vector) = a^2 = Sum_i a_i^2                                    
  template<int N, typename X>
  inline X v_norm(const X*a) { return Dims<N-1>::norm(a); }
  //----------------------------------------------------------------------------
  // return vector dot product = a.b = Sum_i a_i * b_i                          
  template<int N, typename X, typename Y>
  inline X v_dot(const X*a, const Y*b) { return Dims<N-1>::dot(a,b); }
  //----------------------------------------------------------------------------
  // return distance squared between two vectors = |a-b|^2 = Sum_i (a_i-b_i)^2  
  template<int N, typename X, typename Y>
  inline X v_distsq(const X*a, const Y*b) { return Dims<N-1>::distsq(a,b); }
  //----------------------------------------------------------------------------
} // namespace {
#ifndef __NO_AUX_DEFACC
# include <cstring>
# include <new>
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// 3 Definitions of templates and macros that ease the implementation of the  //
//   C-linkable routines in C++.                                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// To use these macros, define a class with the following minimum properties. //
//                                                                            //
// class NAME {                                                               //
// public:                                                                    //
//   static const char* name();                                               //
//   NAME(const double*pars,                                                  //
//        int          npar,                                                  //
//        const char  *file);                                                 //
//   bool NeedMass() const;                                                   //
//   bool NeedVels() const;                                                   //
//   template<int NDIM, typename scalar>                                      //
//   inline void set_time(double       time,                                  //
//                        int          nbod,                                  //
//                        const scalar*mas, // pter may be 0 if NeedMass()==0 //
// 	                  const scalar*pos,                                   //
// 	                  const scalar*vel) // pter may be 0 if NeedVels()==0 //
//                        const;                                              //
//   template<int NDIM, typename scalar>                                      //
//   inline void acc(const scalar*mas,      // pter may be 0 if NeedMass()==0 //
// 	             const scalar*pos,                                        //
// 	             const scalar*vel,      // pter may be 0 if NeedVels()==0 //
// 	             scalar      &pot,                                        //
// 	             scalar      *acc) const;                                 //
// };                                                                         //
//                                                                            //
//                                                                            //
// NOTES  1. In implementing this class, in particular the template over the  //
//           number of dimensions, you may use some helper templates above.   //
//        2. For a static spherical potential, we provide below an implemen-  //
//           tation employing, as template argument, a simpler class for you  //
//           to define. It's recommended to use for these type of potentials. //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
namespace {

  // class AccInstall
  template<class Acceleration>
  class AccInstall {
    //--------------------------------------------------------------------------
    // data                                                                     
    //                                                                          
    // NOTE  we remember the parameters and data file only for the sake of      
    //       NEMO external potentials (not accelerations) so that we can warn   
    //       about re-initialisation with different parameters                  
    double            *Pars;    // current parameter
    int                Npar;    // number of current parameter
    char              *File;    // current data file
    const Acceleration Acc;     // our acceleration
    double             Time;    // last/actual simulation time
    bool               First;   // true only after initialization
    //--------------------------------------------------------------------------
    // private methods                                                          
    void remember(const double*pars, int npar, const char*file) {
      // remember parameters and data file
      if(Pars) delete[] Pars;
      if(pars && npar>0) {
	try {
	  Pars = new double[npar];
	} catch(std::bad_alloc E) {
	  error("[%s:%d]: caught std::bad_alloc\n",__FILE__,__LINE__);
	}
	for(int n=0; n!=npar; ++n) Pars[n] = pars[n];
      } else
	Pars = 0;
      Npar = npar;
      if(File) delete[] File;
      if(file) {
	size_t n = strlen(file)+1;
	try {
	  File = new char[n];
	} catch(std::bad_alloc E) {
	  error("[%s:%d]: caught std::bad_alloc\n",__FILE__,__LINE__);
	}
	strncpy(File,file,n);
      } else
	File = 0;
    }
    //--------------------------------------------------------------------------
    template<typename scalar>
    void set_time(int          ndim,
		  double       time,
		  int          nbod,
		  const scalar*mas,
		  const scalar*pos,
		  const scalar*vel)
      // if time has changed or this is the first call, Acc::set_time()
    {
      if(First || Time != time) {
	Time  = time;
	First = false;
	switch(ndim) {
	case 2: 
	  Acc.template set_time<2>(Time,nbod,mas,pos,vel);
	  break;
	case 3: 
	  Acc.template set_time<3>(Time,nbod,mas,pos,vel);
	  break;
	default:
	  ::error("acceleration: ndim=%d not supported", ndim);
	}
      }
    }
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar>
    void acc_T(double       t,
	       int          nb,
	       const scalar*m,
	       const scalar*x,
	       const scalar*v,
	       const int   *f,
	       scalar      *pot,
	       scalar      *acc,
	       int          add)
      // compute gravity and add/assign pot/acc
    {
      set_time<scalar>(NDIM,t,nb,m,x,v);
      if(add & 1) 
	if(add & 2) {
	  // add both potential and acceleration
	  for(int n=0,nn=0; n!=nb; ++n,nn+=NDIM)
	    if(f==0 || f[n] & 1) {
	      register scalar P,A[NDIM];
	      Acc.template acc<NDIM>(m+n, x+nn, v+nn, P, A);
	      pot[n] += P;
	      v_add<NDIM>(acc+nn,A);
	    }
	} else {
	  // add potential, assign acceleration
	  for(int n=0,nn=0; n!=nb; ++n,nn+=NDIM)
	    if(f==0 || f[n] & 1) {
	      register scalar P;
	      Acc.template acc<NDIM>(m+n, x+nn, v+nn, P, acc+nn);
	      pot[n] += P;
	    }
      } else {
        if(add & 2) {
	  // assign potential, add acceleration
	  for(int n=0,nn=0; n!=nb; ++n,nn+=NDIM)
	    if(f==0 || f[n] & 1) {
	      register scalar A[NDIM];
	      Acc.template acc<NDIM>(m+n, x+nn, v+nn, pot[n], A);
	      v_add<NDIM>(acc+nn,A);
	    }
        } else {
	  // assign both potential and acceleration
	  for(int n=0,nn=0; n!=nb; ++n,nn+=NDIM)
	    if(f==0 || f[n] & 1) {
	      Acc.template acc<NDIM>(m+n, x+nn, v+nn, pot[n], acc+nn);
	    }
	}
      }
    }
    //--------------------------------------------------------------------------
    // public methods                                                           
  public:
    static const char* name()
    {
      return Acceleration::name();
    }
    //--------------------------------------------------------------------------
    AccInstall(const double*pars,
	       int          npar,
	       const char  *file,
	       bool        *need_mass,
	       bool        *need_vels) :
      Pars  ( 0 ),
      Npar  ( 0 ),
      File  ( 0 ),
      Acc   ( pars, npar, file ), 
      First ( true)
    {
      nemo_dprintf(4,"AccInstall() npar=%d, file=%s\n",npar,file);
      remember(pars, npar, file);
      if(need_mass)
	*need_mass = Acc.NeedMass();
      else if(Acc.NeedMass())
	::error("inipotential: cannot use \"%s\", since masses are needed",
		Acceleration::name());
      if(need_vels)
	*need_vels = Acc.NeedVels();
      else if(Acc.NeedVels())
	::error("inipotential: cannot use \"%s\", since velocities are needed",
		Acceleration::name());
    }
    //--------------------------------------------------------------------------
    bool differ(const double*pars, int npar, const char*file)
    {
      if(File==0 && file!=0) return true;
      if(File!=0 && file==0) return true;
      if(File && file && strcmp(File,file)) return true;
      if(npar != Npar) return true;
      for(int n=0; n!=npar; ++n)
	if(pars[n] != Pars[n]) return true;
      return false;
    }
    //--------------------------------------------------------------------------
    void acc(int        nd,
	     double     t,
	     int        nb,
	     const void*m,
	     const void*x,
	     const void*v,
	     const int *f,
	     void      *p,
	     void      *a,
	     int        add,
	     char       type)
    {
      switch(nd) {
      case 2:
	switch(type) {
	case 'f': return acc_T<2>(t,nb,
				  static_cast<const float*>(m),
				  static_cast<const float*>(x),
				  static_cast<const float*>(v),
				  f,
				  static_cast<float*>(p),
				  static_cast<float*>(a),
				  add);
	case 'd': return acc_T<2>(t,nb,
				  static_cast<const double*>(m),
				  static_cast<const double*>(x),
				  static_cast<const double*>(v),
				  f,
				  static_cast<double*>(p),
				  static_cast<double*>(a),
				  add);
	default: ::error("acceleration \"%s\": unknown type ('%s')",
			 Acceleration::name(),type);
	} break;
      case 3:
	switch(type) {
	case 'f': return acc_T<3>(t,nb,
				  static_cast<const float*>(m),
				  static_cast<const float*>(x),
				  static_cast<const float*>(v),
				  f,
				  static_cast<float*>(p),
				  static_cast<float*>(a),
				  add);
	case 'd': return acc_T<3>(t,nb,
				  static_cast<const double*>(m),
				  static_cast<const double*>(x),
				  static_cast<const double*>(v),
				  f,
				  static_cast<double*>(p),
				  static_cast<double*>(a),
				  add);
	default: ::error("acceleration \"%s\": unknown type ('%s')",
			 Acceleration::name(),type);
	} break;
      default: ::error("acceleration \"%s\": ndim=%d unsupported",
		       Acceleration::name(),nd);
      }
    }
    //--------------------------------------------------------------------------
#ifdef POT_DEF
    template<typename scalar>
    void pot_T(const int   *ndim,
	       const scalar*pos,
	       scalar      *acc,
	       scalar      *pot,
	       const scalar*time)
    {
      set_time(*ndim,*time,0,
	       static_cast<const scalar*>(0),
	       static_cast<const scalar*>(0),
	       static_cast<const scalar*>(0));
      switch(*ndim) {
      case 2: Acc.template acc<2>(static_cast<const scalar*>(0), pos, 
				  static_cast<const scalar*>(0), *pot, acc);
	break;
      case 3: Acc.template acc<3>(static_cast<const scalar*>(0), pos, 
				  static_cast<const scalar*>(0), *pot, acc);
	break;
      default: ::error("potential \"%s\": ndim=%d not supported",
		       Acceleration::name(),ndim);
      }
    }
#endif
  }; // class AccInstall<>
  //////////////////////////////////////////////////////////////////////////////
} // namespace {

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// macro __DEF_ACC(NAME)                                                      //
//                                                                            //
// defines the externally linkable routine iniacceleration() from a class     //
// of name NAME that satisfies the above scheme.                              //
//                                                                            //
// Every new call to iniacceleration() will issue a NEW instantination of     //
// the acceleration field (with different parameters hopefully!). We allow    //
// up to 10 instantinations.                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#define __DEF__ACC__NO(NUM)						\
void acceleration##NUM(int        d,					\
		       double     t,					\
		       int        n,					\
		       const void*m,					\
		       const void*x,					\
		       const void*v,					\
		       const int *f,					\
		       void      *p,					\
		       void      *a,					\
		       int        i,					\
		       char       y)					\
{ (MyAcc[NUM])->acc(d,t,n,m,x,v,f,p,a,i,y); }

#define __DEF__ACC(NAME)						\
namespace {								\
  const int                AccMax= 10;					\
  typedef AccInstall<NAME> MyAccInstall;				\
  MyAccInstall            *MyAcc[AccMax] = {0};				\
  int                      AccN  = 0;					\
__DEF__ACC__NO(0)							\
__DEF__ACC__NO(1)							\
__DEF__ACC__NO(2)							\
__DEF__ACC__NO(3)							\
__DEF__ACC__NO(4)							\
__DEF__ACC__NO(5)							\
__DEF__ACC__NO(6)							\
__DEF__ACC__NO(7)							\
__DEF__ACC__NO(8)							\
__DEF__ACC__NO(9)							\
  acc_pter Accs[AccMax] = {&acceleration0,				\
			   &acceleration1,				\
			   &acceleration2,				\
			   &acceleration3,				\
			   &acceleration4,				\
			   &acceleration5,				\
			   &acceleration6,				\
			   &acceleration7,				\
			   &acceleration8,				\
			   &acceleration9};				\
}									\
void iniacceleration(const double*pars,					\
		     int          npar,					\
		     const char  *file,					\
                     acc_pter    *accel,				\
		     bool        *need_m,				\
		     bool        *need_v)				\
{									\
  nemo_dprintf(4,"iniacceleration() called\n");				\
  if(AccN == AccMax) {							\
    ::warning("iniacceleration(): request to initialize "		\
	      "more than %d accelerations of type \"%s\"",		\
	      AccMax, MyAccInstall::name());				\
    *accel = 0;								\
    return;								\
  }									\
  try {									\
    MyAcc[AccN] = new MyAccInstall(pars,npar,file,need_m,need_v);	\
  } catch(std::bad_alloc E) {						\
    ::error("[%s:%d]: caught std::bad_alloc\n",__FILE__,__LINE__);	\
  }									\
  *accel = Accs[AccN++];						\
}

#ifdef POT_DEF
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// macro __DEF_POT(NAME)                                                      //
//                                                                            //
// defines the externally linkable routines inipotential(), potential_float() //
// and potential_double() from a class of name NAME that satisfies the above  //
// scheme.                                                                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

# define __DEF__POT(NAME)						\
namespace {								\
  typedef AccInstall<NAME> MyPotInstall;				\
  MyPotInstall            *MyPot;					\
}									\
void inipotential(const int   *npar,					\
		  const double*pars,					\
		  const char  *file)					\
{									\
  if(MyPot) {								\
    if(MyPot->differ(pars,*npar,file))					\
      ::warning("inipotential(): re-initializing \"%s\" "		\
		"with different parameters or data file",		\
		MyPotInstall::name());					\
    else								\
      ::warning("inipotential(): re-initializing \"%s\" "		\
		"with identical parameters and data file",		\
		MyPotInstall::name());					\
    delete MyPot;							\
  }									\
  try {									\
    MyPot = new MyPotInstall(pars,*npar,file,0,0);			\
  } catch(std::bad_alloc E) {						\
    ::error("[%s:%d]: caught std::bad_alloc\n",__FILE__,__LINE__);	\
  }									\
}									\
void potential_double(const int   *ndim,				\
		      const double*pos,					\
		      double      *acc,					\
		      double      *pot,					\
		      const double*time)				\
{									\
  MyPot->pot_T(ndim,pos,acc,pot,time);					\
}									\
void potential_float(const int  *ndim,					\
		     const float*pos,					\
		     float      *acc,					\
		     float      *pot,					\
		     const float*time)					\
{									\
  MyPot->pot_T(ndim,pos,acc,pot,time);					\
}

#else  // POT_DEF
# define __DEF__POT
#endif // POT_DEF
////////////////////////////////////////////////////////////////////////////////
namespace {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class implementing the template argument described above                 //
  // for a static spherical potential                                         //
  //                                                                          //
  // template argument class StaticSphericalModel must satisfy the following: //
  // class NAME {                                                             //
  // public:                                                                  //
  //   static const char name();                                              //
  //   NAME(const double*pars,                                                //
  //        int          npar,                                                //
  //        const char  *file);                                               //
  //   template<typename scalar>                                              //
  //   inline void potacc(scalar const&radius_squared,                        //
  //                      scalar      &potential,                             //
  //                      scalar      &minus_dpotdr_over_r) const;            //
  // };                                                                       //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<class StaticSphericalModel>
  class SphericalPot : private StaticSphericalModel {
  public:
    static const char* name() { return StaticSphericalModel::name(); }
    //--------------------------------------------------------------------------
    SphericalPot(const double*pars,
		 int          npar,
		 const char  *file) : StaticSphericalModel(pars,npar,file) {}
    //--------------------------------------------------------------------------
    bool NeedMass() const { return false; }
    bool NeedVels() const { return false; }
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar>
    void set_time(double       time,
		  int          nbod,
		  const scalar*mas,
		  const scalar*pos,
		  const scalar*vel) const {}
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar>
    inline void acc(const scalar*,
		    const scalar*pos,
		    const scalar*,
		    scalar      &pot,
		    scalar      *acc) const
    {
      register scalar dpdr_over_r;
      StaticSphericalModel::potacc(v_norm<NDIM>(pos), pot, dpdr_over_r);
      v_asstimes<NDIM>(acc, pos, dpdr_over_r);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
} // namespace {
////////////////////////////////////////////////////////////////////////////////
#endif // __NO_AUX_DEFACC
#endif // __cplusplus
#endif // _defacc_h
