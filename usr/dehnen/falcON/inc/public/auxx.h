// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// auxx.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2002                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// defines                                                                     |
//                                                                             |
// some global macros                                                          |
// some global constants                                                       |
// some global typedefs                                                        |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef included_auxx_h
#define included_auxx_h

#ifndef included_iostream
#  include <iostream>
#  define included_iostream
#endif
#ifndef included_cmath
#  include <cmath>
#  define included_cmath
#endif
#ifndef included_cstdlib
#  include <cstdlib>
#  define included_cstdlib
#endif

#if (defined(ALLOW_MPI) && defined(PARALLEL) )     // have MPI and want it?     
#  ifndef USE_MPI
#    define USE_MPI                                //   then use it             
#  endif
#else                                              // else                      
#  undef  USE_MPI                                  //   don't use it            
#endif

#ifndef included_exit_h
#  include <public/exit.h>
#endif
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// 1. INTERNAL PRECISION                                                        
//                                                                              
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {

#if defined(DOUBLE_SINGLE) || defined(DOUBLE_DOUBLE)

# undef REAL_IS_FLOAT
  typedef double real;
# ifdef linux
  using ::cbrt;
# endif

#else  // defined(DOUBLE_SINGLE) || defined(DOUBLE_DOUBLE)

# define REAL_IS_FLOAT
  typedef float real;
# ifdef linux
  inline real cbrt(real x) { return ::cbrtf(x); }
# endif

#endif  // defined(DOUBLE_SINGLE) || defined(DOUBLE_DOUBLE)

  using std::sqrt;
  using std::exp;
  using std::log;
  using std::pow;
#ifndef linux
  inline real cbrt(real x) { return std::pow(x,0.333333333333333333333); }
#endif
  inline real sqrt(unsigned x) { return sqrt(real(x)); }
  inline real pow (float  x, unsigned i) { return std::pow(x, int(i)); }
  inline real pow (double x, unsigned i) { return std::pow(x, int(i)); }

}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// 2. DATA TRANSFER VIA BODIES OR ARRAYS & PRECISION USED IN ARRAY ARGS         
//                                                                              
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
#if defined(SINGLE_DOUBLE) || defined(DOUBLE_DOUBLE)
  typedef double areal;
#else
  typedef float areal;
#endif
}
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// 3. DIMENSIONALITY                                                            
//                                                                              
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// 3.1 Two dimensions                                                           
//                                                                              
#ifdef TWODIMENSIONAL

#  define NDIM 2                         // # dimensions                        
#  define NSUB 4                         // 2^NDIM                              

#else
//                                                                              
// 3.2 Three dimensions                                                         
//                                                                              

#  define NDIM 3                         // # dimensions                        
#  define NSUB 8                         // 2^NDIM                              

#endif

#ifndef included_vect_h
#  include <public/vect.h>
#endif

#ifndef included_tens_h
#  include <public/tens.h>
#endif

namespace nbdy {
  typedef class vector<real>  vect;      // plain NDIM dimensional vector       
  typedef class sym1  <real>  ten1;      // symmetric tensor of order 1 = vector
  typedef class sym2  <real>  ten2;      // symmetric tensor of order 2 = matrix
  typedef class sym3  <real>  ten3;      // symmetric tensor of order 3         
  typedef class sym4  <real>  ten4;      // symmetric tensor of order 4         
}

#ifdef TWODIMENSIONAL
namespace nbdy {
  typedef real                   amom;
}
#else
namespace nbdy {
  typedef vect                   amom;
}
#endif

////////////////////////////////////////////////////////////////////////////////
//                                                                              
// 4. USEFUL CONSTANTS AND TYPEDEFS                                             
//                                                                              
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  typedef unsigned short   indx;              // use only in non-register memory
  typedef unsigned int     uint;
  //----------------------------------------------------------------------------
  const   real zero    = 0.,
               sixth   = 0.166666666666666666666667,
               fifth   = 0.2,
               quarter = 0.25,
               third   = 0.333333333333333333333333,
               half    = 0.5,
               one     = 1.,
               two     = 2.,
               three   = 3.,
               four    = 4.,
               six     = 6.,
               ten     = 10.,
               if2     = 0.5,                           // 1 / 2!              
               if3     = 0.166666666666666666666667,    // 1 / 3!              
               if4     = 0.0416666666666666666666667,   // 1 / 4!              
               if5     = 0.00833333333333333333333333,  // 1 / 5!              
               if6     = 0.00138888888888888888888889;  // 1 / 6!              

}
////////////////////////////////////////////////////////////////////////////////
#endif                                              // included_auxx_h          
