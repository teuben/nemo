// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// auxx.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2003                                          |
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
#ifndef falcON_included_auxx_h
#define falcON_included_auxx_h

#ifndef falcON_included_iostream
#  include <iostream>
#  define falcON_included_iostream
#endif
#ifndef falcON_included_cmath
#  include <cmath>
#  define falcON_included_cmath
#endif
#ifndef falcON_included_cstdlib
#  include <cstdlib>
#  define falcON_included_cstdlib
#endif

#if defined(falcON_MPI) && defined(falcON_PARALLEL)
                                                   // have MPI and want it?     
#  ifndef falcON_USE_MPI
#    define falcON_USE_MPI                         //   then use it             
#  endif
#else                                              // else                      
#  undef  falcON_USE_MPI                           //   don't use it            
#endif

#ifndef falcON_included_exit_h
#  include <public/exit.h>
#endif
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// 1. INTERNAL PRECISION                                                        
//                                                                              
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {

#if defined(falcON_DOUBLE_SINGLE) || defined(falcON_DOUBLE_DOUBLE)

# undef falcON_REAL_IS_FLOAT
  typedef double real;
# ifdef linux
  using ::cbrt;
# endif

#else

# define falcON_REAL_IS_FLOAT
  typedef float real;
# ifdef linux
  inline real cbrt(real x) { return ::cbrtf(x); }
# endif

#endif

  using std::sqrt;
  using std::exp;
  using std::log;
  using std::pow;

#ifndef linux
  inline real cbrt(real x)
  { 
    return real( std::pow( double(x), 0.333333333333333333333 ) );
  }
#endif

#if defined(falcON_REAL_IS_FLOAT) && defined (__GNUC__) && (__GNUC__ < 3)
#  ifdef linux
  inline real sqrt(real x) { return sqrtf(x); }
#  else
  inline real sqrt(real x) { return sqrt(double(x)); }
#  endif
#endif
  inline real sqrt(unsigned x)           { return sqrt(real(x)); }
  inline real pow (float  x, unsigned i) { return std::pow(x, int(i)); }
  inline real pow (double x, unsigned i) { return std::pow(x, int(i)); }

}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// 2. DATA TRANSFER VIA BODIES OR ARRAYS & PRECISION USED IN ARRAY ARGS         
//                                                                              
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
#if defined(falcON_SINGLE_DOUBLE) || defined(falcON_DOUBLE_DOUBLE)
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
#ifndef falcON_NDIM
#  define falcON_NDIM 3
#endif
//                                                                              
// 3.1 # dimensions                                                             
//                                                                              
#if falcON_NDIM == 2
#  define falcON_NSUB 4                  // 2^NDIM                              
#elif falcON_NDIM == 3
#  define falcON_NSUB 8                  // 2^NDIM                              
#else
#  error falcON_NDIM neither 2 nor 3
#endif

namespace nbdy {
  const int Ndim = falcON_NDIM;
  const int Nsub = falcON_NSUB;
}

//                                                                              
// 3.2 vector and tensors                                                       
//                                                                              

#ifndef falcON_included_vect_h
#  include <public/vect.h>
#endif

#ifndef falcON_included_tens_h
#  include <public/tens.h>
#endif

namespace nbdy {
  typedef class vector<real>   vect;     // plain Ndim dimensional vector       
  typedef class vector<double> vect_d;   // plain Ndim dimensional vector       
  typedef class sym1  <real>   ten1;     // symmetric tensor of order 1 = vector
  typedef class sym2  <real>   ten2;     // symmetric tensor of order 2 = matrix
  typedef class sym3  <real>   ten3;     // symmetric tensor of order 3         
  typedef class sym4  <real>   ten4;     // symmetric tensor of order 4         

#if falcON_NDIM == 2
  typedef real                   amom;
  typedef double                 amom_d;
#else
  typedef vect                   amom;
  typedef vector<double>         amom_d;
#endif

}

////////////////////////////////////////////////////////////////////////////////
//                                                                              
// 4. USEFUL CONSTANTS AND TYPEDEFS                                             
//                                                                              
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  typedef unsigned short   indx;              // use only in non-register memory
  typedef unsigned int     uint;
  //----------------------------------------------------------------------------
  const   real zero       = 0.,
               sixth      = 0.166666666666666666666667,
               fifth      = 0.2,
               quarter    = 0.25,
               third      = 0.333333333333333333333333,
               half       = 0.5,
               one        = 1.,
               threehalfs = 1.5,
               two        = 2.,
               three      = 3.,
               four       = 4.,
               six        = 6.,
               eight      = 8.,
               ten        = 10.,
               twelve     = 12.,
               sixten     = 16.,
               if2        = 0.5,                         // 1 / 2!              
               if3        = 0.166666666666666666666667,  // 1 / 3!              
               if4        = 0.0416666666666666666666667, // 1 / 4!              
               if5        = 0.00833333333333333333333333,// 1 / 5!              
               if6        = 0.00138888888888888888888889;// 1 / 6!              

}
////////////////////////////////////////////////////////////////////////////////
#endif                                              // falcON_included_auxx_h   
