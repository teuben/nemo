// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    inc/public/basic.h                                                 
///                                                                             
/// \brief   Declarations of some basic functions needed in most source code    
///                                                                             
///          Contents:                                                          
///          \li exception handling (falcON::exception)                         
///          \li falcON::compile_info                                           
///          \li memory allocation and de-allocation support                    
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    2002-2008                                                          
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2002-2008  Walter Dehnen                                       
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
#ifndef falcON_included_types_h
#  include <public/types.h>
#endif

namespace falcON {
  //----------------------------------------------------------------------------
  /// Returns the falcON base directory (reads the enviroment variable FALCON)
  inline const char* directory() { return getenv("FALCON"); }
  //----------------------------------------------------------------------------
  /// Returns the falcON library directory
  const char* libdir();
  //----------------------------------------------------------------------------
  /// to handle error output
  struct Error : public WDutils::Error
  {
    /// default constructor
    Error()
      : WDutils::Error("falcON") {}
    /// constructor: get file name & line number
    Error(const char*f, int l)
      : WDutils::Error(f,l,"falcON") {}
  };
  //----------------------------------------------------------------------------
  /// to handle warning output
  struct Warning : public WDutils::Warning
  {
    /// default constructor
    Warning()
      : WDutils::Warning("falcON") {}
    /// constructor: get file name & line number
    Warning(const char*f, int l)
      : WDutils::Warning(f,l,"falcON") {}
  };
  //----------------------------------------------------------------------------
  /// to handle debug output
  struct DebugInformation : public WDutils::DebugInformation
  {
    /// default constructor
    DebugInformation()
      : WDutils::DebugInformation("falcON ") {}
    /// constructor: get file name & line number
    DebugInformation(const char*f, int l)
      : WDutils::DebugInformation(f,l,"falcON ") {}
  };
  //----------------------------------------------------------------------------
  /// \name macros controling the usage of throw exception vs error             
  //@{                                                                          
  /// error message (and exit), reporting [file:line]
  /// use "falcON_Error(fmt, data)" instead of "error(fmt, data)"
#define falcON_Error           falcON::Error(__FILE__,__LINE__)
  /// error message (and exit), NOT reporting [file:line]
  /// use "falcON_Error(fmt, data)" instead of "error(fmt, data)"
#define falcON_ErrorN          falcON::Error()
  /// print warning message, reporting [file:line]
  /// use "falcON_Warning(fmt, data)" instead of "warning(fmt, data)"
#define falcON_Warning	       falcON::Warning(__FILE__,__LINE__)
  /// print warning message, NOT reporting [file:line]
  /// use "falcON_Warning(fmt, data)" instead of "warning(fmt, data)"
#define falcON_WarningN	       falcON::Warning()
  //----------------------------------------------------------------------------
#ifdef WDutils_EXCEPTIONS
#  define falcON_EXCEPTIONS
#  define falcON_THROWING      throw(falcON::exception)
#  define falcON_THROWER       throw WDutils::Thrower
#  define falcON_THROWN        throw falcON::exception
#  define falcON_RETHROW(E)    WDutils_RETHROW(E)
#else
#  undef falcON_EXCEPTIONS
#  define falcON_THROWER       falcON::Error
  /// use instead of <tt> throw(falcON::exception) </tt> after function
  /// declaration
#  define falcON_THROWING 
  /// instead of throwing an exception: error 
  /// use "WDutils_THROWN(fmt, data)" instead of "error(fmt, data)" or "throw
  /// WDutils::exception(fmt, data)"
#  define falcON_THROWN	       falcON_ErrorN
  /// use "falcON_RETHROW(E)" to re-throw a caught exception "E"
#  define falcON_RETHROW(E)    falcON_Error  (text(E))
#endif
  /// instead of throwing an exception: error with [file:line]
  /// use "falcON_THROW(fmt, data)" instead of "error(fmt, data)" or "throw
  /// falcON::exception(fmt, data)"
#  define falcON_THROW         falcON_THROWER(__FILE__,__LINE__)
  //@}
#undef  DebugInfo
#define DebugInfo    falcON::DebugInformation(__FILE__,__LINE__)
#undef  DebugInfoN
#define DebugInfoN   falcON::DebugInformation()
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
    real_is_double  = 16 };
  //----------------------------------------------------------------------------
  inline Status CurrentStatus() {
    int status = public_version;
#ifdef falcON_PROPER
    status |= proper_version;
#endif
#ifdef DEBUG
    status |= debug_version;
#endif
#ifdef falcON_NEMO
    status |= nemo_version;
#endif
#ifdef falcON_SPH
    status |= sph_version;
#endif
#ifdef falcON_DOUBLE
    status |= real_is_double;
#endif
    return Status(status);
  }
  Status LibraryStatus();
  void CheckAgainstLibrary(falcON::Status,const char*) falcON_THROWING;
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
  WDutils::NewArray<TYPE>(SIZE,__FILE__,__LINE__,"falcON ")
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
#define falcON_DEL_A(P) WDutils::DelArray(P,__FILE__,__LINE__,"falcON ")
  //////////////////////////////////////////////////////////////////////////////
  ///                                                                           
  /// C MACRO to be used for object de-allocation                               
  ///                                                                           
  /// Calling falcON::DelObject<TYPE>(), which in case of an error generates    
  /// an error message detailing the source file and line of the call.          
  ///                                                                           
  /// \param P  pointer to object to be de-allocated                            
#define falcON_DEL_O(P) DelObject(P,__FILE__,__LINE__,"falcON ")
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

#if !defined(__GNUC__) && !defined (__PGCC__) && \
    !defined (__INTEL_COMPILER)	&& !defined (__PATHCC__)
#  warning " "
#  warning " falcON: you are using an unknown compiler"
#  warning " compilation may fail or produce buggy code"
#  warning " "
#endif

#if defined(__PGCC__) || (defined(__GNUC__) && __GNUC__ < 3)
#  define falcON_non_standard_math
// #  warning " you are using a known non-standard C++ compiler; compilation may fail or produce buggy code"
#endif
} // namespace falcON
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_basic_h
