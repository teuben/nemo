// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// auxx.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// defines                                                                     |
//                                                                             |
// some global macros                                                          |
// some global constants                                                       |
// some external constants                                                     |
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
#ifndef falcON_included_frst_h
#  include <public/frst.h>
#endif

////////////////////////////////////////////////////////////////////////////////
//                                                                              
// general falcON macros                                                        
//                                                                              
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
// MPI parallelism                                                              
//------------------------------------------------------------------------------
#if defined(falcON_MPI) && defined(falcON_PARALLEL)
                                                   // have MPI and want it?     
#  ifndef falcON_USE_MPI
#    define falcON_USE_MPI                         //   then use it             
#  endif
#else                                              // else                      
#  undef  falcON_USE_MPI                           //   don't use it            
#endif

//------------------------------------------------------------------------------
// SPH code?                                                                    
//------------------------------------------------------------------------------
#if defined(falcON_SPH) && !defined(falcON_INDI)
#  warning 'falcON_SPH' #defined but not 'falcON_INDI'
#  warning will #define 'falcON_INDI' now
#  define falcON_INDI
#endif

//------------------------------------------------------------------------------
// adaptive individual softening?                                               
//------------------------------------------------------------------------------
#if defined(falcON_ADAP) && !defined(falcON_INDI)
#  warning
#  warning 'falcON_ADAP' #defined but not 'falcON_INDI' 
#  warning we will #undef 'falcON_ADAP'
#  warning
#endif
//------------------------------------------------------------------------------
// expansion order used in gravity                                              
//------------------------------------------------------------------------------
#if defined(Walter) && !defined(falcON_included_ordr_h)
#  include <walter/ordr.h>
#endif

#ifndef falcON_ORDER
#  define falcON_ORDER 3
#endif
//------------------------------------------------------------------------------
// falcON_REAL_IS_FLOAT  and type nbdy::real                                    
//------------------------------------------------------------------------------
namespace nbdy {
#if defined(falcON_DOUBLE_SINGLE) || defined(falcON_DOUBLE_DOUBLE)
#  undef falcON_REAL_IS_FLOAT
  typedef double real;
#else
#  define falcON_REAL_IS_FLOAT
  typedef float real;
#endif
}
//------------------------------------------------------------------------------
// data type used for arrays (only important for FORTRAN & C export)            
//------------------------------------------------------------------------------
namespace nbdy {
#if defined(falcON_SINGLE_DOUBLE) || defined(falcON_DOUBLE_DOUBLE)
  typedef double areal;
#else
  typedef float areal;
#endif
}
//------------------------------------------------------------------------------
// SSE code?                                                                    
//------------------------------------------------------------------------------
#if defined(falcON_REAL_IS_FLOAT) && defined(falcON_SSE)
#  define  falcON_SSE_CODE
#else
#  undef   falcON_SSE_CODE
#endif
//------------------------------------------------------------------------------
// Dimensionality                                                               
//------------------------------------------------------------------------------
#ifndef falcON_NDIM
#  define falcON_NDIM 3
#endif
namespace nbdy {
  const int Ndim = falcON_NDIM;                    // # dimensions              
  const int Nsub = 1<<Ndim;                        // 2^(# dimensions)          
}
//------------------------------------------------------------------------------
// macros used to test timings                                                  
//------------------------------------------------------------------------------
#ifdef TEST_TIMING
#  include <ctime>
#  define SET_I        std::clock_t __C0_TIMING = std::clock();
#  define SET_T(TEXT)  std::cerr<<TEXT					\
                                <<(std::clock()-__C0_TIMING)/		\
                                   double(CLOCKS_PER_SEC)<<'\n';	\
                       __C0_TIMING = std::clock();
#else
#  define SET_I        {}
#  define SET_T(TEXT)  {}
#endif

////////////////////////////////////////////////////////////////////////////////
//                                                                              
// elementary falcON functions and types                                        
//                                                                              
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
// falcON error & warnings                                                      
//------------------------------------------------------------------------------
#ifndef falcON_included_exit_h
#  include <public/exit.h>
#endif

//------------------------------------------------------------------------------
//  vectors and tensors                                                         
//------------------------------------------------------------------------------
#ifndef falcON_included_inln_h
#  include <public/inln.h>
#endif
#ifndef falcON_included_tupl_h
#  include <public/tupl.h>
#endif
#ifndef falcON_included_tn3D_h
#  include <public/tn3D.h>
#endif

namespace nbdy {
  typedef tupel<Ndim,real  > vect;                 // plain Ndim dim vector     
  typedef tupel<Ndim,double> vect_d;               // plain Ndim dim vector     
  typedef pseudo_tupel      <Ndim,areal>   ps_vect;
  typedef const_pseudo_tupel<Ndim,areal> c_ps_vect;

#if falcON_NDIM == 2
  typedef real               amom;
  typedef double             amom_d;
#else
  typedef vect               amom;
  typedef vect_d             amom_d;
#endif
}
//------------------------------------------------------------------------------
// useful constants and typedefs                                                
//------------------------------------------------------------------------------
namespace nbdy {
  typedef unsigned short   indx;                   // use only in non-register  
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
               twelve     = 12.;
}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_included_auxx_h    
