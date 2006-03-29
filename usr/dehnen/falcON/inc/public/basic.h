// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    inc/public/basic.h                                                 
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    2002-2006                                                          
///                                                                             
/// \brief   Declarations of some basic functions needed in most source code    
///                                                                             
///          Contents:                                                          
///          \li falcON::error, falcON::warning, falcON::message, and           
///              falcON::debug_info                                             
///          \li exception handling (falcON::exception)                         
///          \li falcON::compile_info and falcON::run_info                      
///          \li memory allocation and de-allocation support                    
///          \li 16byte memory alignment support                                
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2002-2006  Walter Dehnen                                       
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
#ifndef falcON_included_basic_h
#define falcON_included_basic_h 1

#ifndef falcON_included_iostream
#  include <iostream>
#  define falcON_included_iostream
#endif
#ifndef falcON_included_string
#  include <string>
#  define falcON_included_string
#endif
#ifndef falcON_included_cstdlib
#  include <cstdlib>
#  define falcON_included_cstdlib
#endif
#ifndef falcON_included_cmath
#  include <cmath>
#  define falcON_included_cmath
#endif
#ifndef falcON_included_utils_h
#  include <public/utils.h>
#endif

namespace falcON {
  //----------------------------------------------------------------------------
  /// Returns the falcON base directory (reads the enviroment variable FALCON)
  inline const char* directory() { return getenv("FALCON"); }
  //----------------------------------------------------------------------------
  /// Returns the falcON library directory
  const char* libdir();
  //----------------------------------------------------------------------------
  /// Print an error message to stderr and abort.
  /// Uses a printf() style format string as first argument, further arguments
  /// must match format, exactly as in printf, which will be called.
  /// \param fmt gives the format in C printf() style
  void error(const char*fmt, ... );
  //----------------------------------------------------------------------------
  /// Print a warning message to stderr.
  /// Uses a printf() style format string as first argument, further arguments
  /// must match format, exactly as in printf, which will be called.
  /// \param fmt gives the format in C printf() style
  void warning(const char*fmt, ...);
  //----------------------------------------------------------------------------
  /// Print debugging information to stderr.                                    
  /// Debugging information is printed if the debugging level exceeds the       
  /// debugging depth (1st argument). The debugging level is found from         
  /// falcON::run_info::debug_level() and initialized at start-up of main().    
  /// Uses a printf() style format string as first argument, further arguments  
  /// must match format, exactly as in printf, which will be called.            
  /// \param deb gives the debugging depth: info is printed only if it is       
  ///            less or equal to debugging level                               
  /// \param fmt gives the format in C printf() style                           
  void debug_info(int deb, const char*fmt, ... );
  //----------------------------------------------------------------------------
  /// \name macros controling the usage of throw exception vs error             
  //@{                                                                          
#if 1
  /// use instead of <tt> throw(falcON::exception) </tt> after function
  /// declaration
#  define falcON_THROWING 
  /// use "falcON_THROW(fmt, data)" instead of "error(fmt, data)" or "throw
  /// falcON::exception(fmt, data)"
#  define falcON_THROW         falcON::error
  /// use "falcON_RETHROW(E)" to re-throw a caught exception "E"
#  define falcON_RETHROW(E)    falcON::error(text(E))
#else
#  define falcON_THROWING      throw(falcON::exception)
#  define falcON_THROW         throw falcON::exception
#  define falcON_RETHROW(E)    throw E
#endif
  //@}
  //============================================================================
  //
  // mechanism to avoid status mismatch
  //
  //============================================================================
  enum Status {
    public_version  = 0,
    proper_version  = 1,
    debug_version   = 2,
    nemo_version    = 4,
    sph_version     = 8,
    mpi_version     = 16,
    real_is_double  = 32 };
  //----------------------------------------------------------------------------
  inline Status CurrentStatus() {
    int status = public_version;
#ifdef falcON_DEBUG
    status |= debugg_version;
#endif
#ifdef falcON_PROPER
    status |= proper_version;
#endif
#ifdef falcON_NEMO
    status |= nemo_version;
#endif
#ifdef DEBUG
    status |= debug_version;
#endif
#ifdef falcON_SPH
    status |= sph_version;
#endif
#ifdef falcON_MPI
    status |= mpi_version;
#endif
#if defined(falcON_DOUBLE)
    status |= real_is_double;
#endif
    return Status(status);
  }
  Status LibraryStatus();
  void CheckAgainstLibrary(falcON::Status) falcON_THROWING;
  //============================================================================
  ///
  /// support for code information; initialized in ::main() at start-up
  ///
  namespace compile_info {
    bool const&is_set();              ///< is info provided?
    const char*version();             ///< version number of main
    const char*origin();              ///< version origin (author)
    const char*time();                ///< time of compilation
    const char*compiler();            ///< compiler used
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // useful macros for simple error messages                                    
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
#define falcON_Error(MSGS)					\
  error("[%s.%d]: %s",__FILE__,__LINE__,MSGS)
  //----------------------------------------------------------------------------
#define falcON_Except(MSGS)					\
  falcON_THROW("[%s.%d]: %s",__FILE__,__LINE__,MSGS)
  //----------------------------------------------------------------------------
#define falcON_Warning(MSGS)					\
  warning("[%s.%d]: %s",__FILE__,__LINE__,MSGS)
  //----------------------------------------------------------------------------
#define falcON_ErrorF(MSGS,FUNC)				\
  error("[%s.%d]: in %s: %s",__FILE__,__LINE__,FUNC,MSGS)
  //----------------------------------------------------------------------------
#define falcON_ExceptF(MSGS,FUNC)					\
  falcON_THROW("[%s.%d]: in %s: %s",__FILE__,__LINE__,FUNC,MSGS)
  //----------------------------------------------------------------------------
#define falcON_WarningF(MSGS,FUNC)				\
  warning("[%s.%d]: in %s: %s",__FILE__,__LINE__,FUNC,MSGS)
  //////////////////////////////////////////////////////////////////////////////
  ///                                                                           
  /// C MACRO to be used for array allocation                                   
  ///                                                                           
  /// Calling WDutils::NewArray<TYPE>(), which in case of an error generates an 
  /// error message detailing the source file and line of the call. In case     
  /// the debugging level exceeds 10, we always print debugging information     
  /// about memory allocation.                                                  
  ///                                                                           
  /// \param  TYPE name of the element type                                     
  /// \param  SIZE number of elements                                           
#define falcON_NEW(TYPE,SIZE)					\
  WDutils::NewArray<TYPE>(SIZE,__FILE__,__LINE__,"falcON")
  //////////////////////////////////////////////////////////////////////////////
  ///                                                                           
  /// C MACRO to be used for array de-allocation                                
  ///                                                                           
  /// Calling WDutils::DelArray<TYPE>(), which in case of an error generates an 
  /// error message detailing the source file and line of the call. In case     
  /// the debugging level exceeds 10, we always print debugging information     
  /// about memory de-allocation.                                               
  ///                                                                           
  /// \param P  pointer to be de-allocated                                      
#define falcON_DEL_A(P) WDutils::DelArray(P,__FILE__,__LINE__,"falcON")
  //////////////////////////////////////////////////////////////////////////////
  ///                                                                           
  /// C MACRO to be used for object de-allocation                               
  ///                                                                           
  /// Calling falcON::DelObject<TYPE>(), which in case of an error generates    
  /// an error message detailing the source file and line of the call.          
  ///                                                                           
  /// \param P  pointer to object to be de-allocated                            
#define falcON_DEL_O(P) DelObject(P,__FILE__,__LINE__,"falcON")
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // reporting actions to narrow code position of mysterious aborts           //
  //                                                                          //
  // instead of functions begin_report() and end_report(), we use an object   //
  // struct report, so that the implicit call to the destructor can be used   //
  // to replace a call to end_report().                                       //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
#ifdef falcON_PROPER
#  ifndef falcON_RepAction
#    define falcON_RepAction 0                     // used to be 1, new Sep 2004
#  endif
  struct report {                                  // for reporting actions     
    static void open_file(const char*,             // open file for report      
			  const char* =0);         // I: comment/history        
    static const char*file_name();                 // return file's name        
    static void close_file();                      // close file for report     
    static void info(const char*, ...);            // to be called anywhere     
    report(const char*, ...);                      // to be called explicitly   
   ~report();                                      // to be called implicitly   
  };
  //////////////////////////////////////////////////////////////////////////////
#else
  // trivial implementation to avoid macros in application files                
  struct report { report(const char*, ...) {}
                  static void info(const char*, ...) {} };
#endif
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // check for known compiler                                                 //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////

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

  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // elementary math functions for float                                      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////

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
  using ::cbrt;
  inline float cbrt(float x)          { return ::cbrtf(x); }
#else
  inline float cbrt(float x)          { 
    return float( std::pow( double(x), 0.333333333333333333333 ) );
  }
  inline double cbrt(double x)        { 
    return std::pow( x, 0.333333333333333333333 );
  }
#endif
  //============================================================================
  //
  /// \defgroup  Mem16  memory alignment to 16 bytes
  /// \name memory alignment to 16 bytes
  //@{

  //----------------------------------------------------------------------------
  /// Macro enforcing memory alignment to 16 bytes
  /// \ingroup Mem16
  ///
  /// Forces the corresponding variable/type to be 16-byte aligned; Works with
  /// icc (icpc) and gcc (g++) [versions > 3]; Use it like \code 
  ///    struct falcON__align16 name { ... };              \endcode
#if defined (__INTEL_COMPILER)
#  define falcON__align16 __declspec(align(16)) 
#elif defined (__GNUC__) && __GNUC__ > 2
#  define falcON__align16 __attribute__ ((aligned(16)))
#else
#  define falcON__align16
#endif
  //----------------------------------------------------------------------------
  //
  // methods to check memory alignment
  //
  //----------------------------------------------------------------------------
  /// is a given memory address aligned?
  /// \ingroup Mem16
  /// \param p  memory address to be tested
  /// \param al alignemt to a bytes will be tested
  inline bool is_aligned(const void*p, int al)
  {
    return size_t(p) % al == 0;
  }
  //----------------------------------------------------------------------------
  /// is a given memory address aligned to a 16 bytes memory location?
  /// \param p  memory address to be tested
  inline bool is_aligned16(const void*p)
  {
    return size_t(p) % 16 == 0;
  }
  //----------------------------------------------------------------------------
  ///
  /// Allocate memory at a address aligned to a 16 byte memory location
  /// \ingroup Mem16
  ///
  /// Will allocate slightly more memory than required to ensure we can find
  /// the required amount at a 16b memory location. To de-allocate, you \b must
  /// use falcON::free16(), otherwise an error will occur!
  ///
  /// \return   a newly allocated memory address at a 16 byte memory location
  /// \param n  number of bytes to allocate
  /// \version  debugged 02-09-2004 WD
  inline void* malloc16(size_t n) falcON_THROWING
  {
    // linear memory model:                                                 
    // ^    = 16byte alignment points                                       
    // S    = sizeof(void*) (assumed 4 in this sketch)                      
    // PPPP = memory where p is stored (needed in deletion)                 
    // def0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef01
    //    ^               ^               ^               ^               ^ 
    //        |-S-|       |  ->  at least n bytes                           
    //        |   |       |                                                 
    //        p   q   ->  q                                                 
    //    |--off--|-16-off|                                                 
    //                |-S-|                                                 
    //                PPPP|                                                 
    //                                                                      
    // the original allocation gave p, we return the shifted q and remember 
    // the original allocation address at PPPP.                             
    char *p = falcON_NEW(char,n+16+sizeof(void*)); // alloc: (n+16)b + pter     
    char *q = p + sizeof(void*);                   // go sizeof pointer up      
    size_t off = size_t(q) % 16;                   // offset from 16b alignment 
    if(off) q += 16-off;                           // IF offset, shift          
    *((void**)(q-sizeof(void*))) = p;              // remember allocation point 
    return static_cast<void*>(q);                  // return aligned address    
  }
  //----------------------------------------------------------------------------
  ///
  /// De-allocate memory previously allocated with falcON::malloc16()
  /// \ingroup Mem16
  ///
  /// This routine \b must be used to properly de-allocate memory that has been
  /// previously allocated by falcON::malloc16(); other de-allocation will
  /// inevitably result in a run-time \b error!
  ///
  /// \param q  pointer previously allocated by falcON::malloc16()
  inline void  free16  (void*q) falcON_THROWING
  {
    falcON_DEL_A( (char*)( *( (void**) ( ( (char*)q )-sizeof(void*) ) ) ) );
  }
  //----------------------------------------------------------------------------
  ///
  /// Allocate memory at a address alignged to at a 16 byte memory location
  /// \ingroup Mem16
  ///
  /// Will allocate slightly more memory than required to ensure we can find
  /// the required amount at a 16b memory location. To de-allocate, you \b must
  /// use falcON::delete16<>(), otherwise an error will occur!
  ///
  /// \return   a newly allocated memory address at a 16 byte memory location
  /// \param T  (template parameter) type of objects to allocate
  /// \param n  number of objects to allocate
  template<typename T> inline T* new16(size_t n) falcON_THROWING
  {
    return static_cast<T*>(malloc16(n * sizeof(T)));
  }
  //----------------------------------------------------------------------------
  ///
  /// De-allocate memory previously allocated with falcON::new16().
  /// \ingroup Mem16
  ///
  /// This routine \b must be used to properly de-allocate memory that has been
  /// previously allocated by falcON::new16(); other de-allocation will
  /// inevitably result in an error.
  ///
  /// \param T  (template parameter) type of objects q points to.
  /// \param q  pointer previously allocated by falcON::new16()
  template<typename T> inline void delete16(T* q)falcON_THROWING 
  {
    free16(static_cast<void*>(q));
  }
  //@}
  //////////////////////////////////////////////////////////////////////////////
} // namespace falcON
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_basic_h
