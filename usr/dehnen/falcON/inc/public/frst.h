// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// frst.h                                                                      |
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
// elementary macros not specific for N-body code.                             |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_frst_h
#define falcON_included_frst_h

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

////////////////////////////////////////////////////////////////////////////////
//                                                                              
// check for known compiler                                                     
//                                                                              
////////////////////////////////////////////////////////////////////////////////

#if !defined(__GNUC__) && !defined (__PGCC__) && !defined (__INTEL_COMPILER)
#  warning " "
#  warning " falcON: you are using an unknown compiler"
#  warning " compilation may fail or produce buggy code"
#  warning " "
#endif

#if defined(__PGCC__) || (defined(__GNUC__) && __GNUC__ < 3)
#  define falcON_non_standard_math
// #  warning " you are using a known non-standard C++ compiler; compilation may fail or produce buggy code"
#endif

////////////////////////////////////////////////////////////////////////////////
//                                                                              
// elementary math functions for real                                           
//                                                                              
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {

  using std::sqrt;
  using std::exp;
  using std::log;
  using std::pow;

#if defined(falcON_non_standard_math)
#  ifdef linux
  inline float sqrt(float x)          { return ::sqrtf(x); }
  inline float exp (float x)          { return ::expf (x); }
  inline float log (float x)          { return ::logf (x); }
  inline float pow (float x, float y) { return ::powf (x,y); }
#  else
  inline float sqrt(float x)          { return sqrt(double(x)); }
  inline float exp (float x)          { return exp (double(x)); }
  inline float log (float x)          { return log (double(x)); }
  inline float pow (float x, float y) { return pow(double(x),double(y)); }
#  endif
#endif

#ifdef linux
  inline float cbrt(float x)          { return ::cbrtf(x); }
  using ::cbrt;
#else
  inline float cbrt(float x)          { 
    return float( std::pow( double(x), 0.333333333333333333333 ) );
  }
  inline double cbrt(double x)        { 
    return std::pow( x, 0.333333333333333333333 );
  }
#endif

  inline double sqrt(unsigned x) { return sqrt(double(x)); }
}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// memory alignment, in particular to 16 bytes                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
// macro falcON__align16                                                        
//     - forces the corresponding variable/type to be 16-byte aligned           
//     - works with icc (icpc) and gcc (g++) [versions > 3]                     
//     - use it like "struct falcON__align16 name { ... };                      
//------------------------------------------------------------------------------
#if defined (__INTEL_COMPILER)
#  define falcON__align16 __declspec(align(16)) 
#elif defined (__GNUC__) && __GNUC__ > 2
#  define falcON__align16 __attribute__ ((aligned(16)))
#else
#  define falcON__align16
#endif
//------------------------------------------------------------------------------
namespace nbdy {
  //---------------------------------------------------------------------------+
  //                                                                           |
  // methods to check memory alignment                                         |
  //                                                                           |
  //---------------------------------------------------------------------------+
  inline bool is_aligned(const void*p, int al)
  {
    return size_t(p) % al == 0;
  }
  //----------------------------------------------------------------------------
  inline bool is_aligned16(const void*p)
  {
    return size_t(p) % 16 == 0;
  }
  //---------------------------------------------------------------------------+
  //                                                                           |
  // methods malloc16() and free16()                                           |
  //                                                                           |
  //     - allocate and free memory that is 16-byte aligned                    |
  //     - free16() must be used to free memory allocated by malloc16()        |
  //     - free16() must not be used to free memory allocated otherwise.       |
  //                                                                           |
  //---------------------------------------------------------------------------+
  inline void* malloc16(                           // R: memory address         
                        int const&n)               // I: # bytes                
  {
    char *p = new char[n+16+sizeof(void*)];        // alloc: (n+16)b + pter     
    char *q = p + 16;                              // go 16b up                 
    size_t off = size_t(q) % 16;                   // offset from 16b alignment 
    if(off) q += 16-off;                           // IF offset, shift          
    *((void**)(q-sizeof(void*))) = p;              // remember allocation point 
    return (void*)q;                               // return aligned address    
  }
  //----------------------------------------------------------------------------
  inline void  free16  (void*const&q)              // I: memory address         
  {
    delete[] (char*)( *( (void**) ( ( (char*)q )-sizeof(void*) ) ) );
  }
  //----------------------------------------------------------------------------
  template<typename T> inline T* new16(int const&n)
  {
    return static_cast<T*>(malloc16(n * sizeof(T)));
  }
  //----------------------------------------------------------------------------
  template<typename T> inline void delete16(T*const&q)
  {
    free16(static_cast<void*>(q));
  }
}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_included_frst_h    
