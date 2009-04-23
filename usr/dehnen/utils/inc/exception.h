// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/exception.h
///
/// \author Walter Dehnen
///                                                                             
/// \date   2000-2008
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2008 Walter Dehnen
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
#ifndef WDutils_included_exception_h
#define WDutils_included_exception_h

#ifndef WDutils_included_string
# define WDutils_included_string
# include <string>
#endif
// /////////////////////////////////////////////////////////////////////////////
//                                                                              
//  WDutils                                                                     
//                                                                              
/// generally useful code of Walter Dehnen, used in project falcON,             
/// public under the GNU public licence                                         
///                                                                             
// /////////////////////////////////////////////////////////////////////////////
namespace WDutils {
  // ///////////////////////////////////////////////////////////////////////////
  //
  //  macro for compile-time assertion, stolen from the boost library
  //
  //----------------------------------------------------------------------------
  template<bool> struct STATIC_ASSERTION_FAILURE;
  template<>     struct STATIC_ASSERTION_FAILURE<true> { enum { value = 1 }; };
  /// \brief macro for compile-time assertion
  ///
  /// \code
  ///   WDutilsStaticAssert(constant expression);
  /// \endcode
  /// will cause a compiler error if the expression evaluates to false. This
  /// relies on sizeof() an incomplete type causing an error, though
  /// "STATIC_ASSERTION_FAILURE" along with the line of the actual
  /// instantination causing the error will also appear in the
  /// compiler-generated error message.
#define WDutilsStaticAssert(TEST)				\
  enum { __DUMMY = sizeof(WDutils::STATIC_ASSERTION_FAILURE<	\
    static_cast<bool>(TEST)>)					\
  }
  // ///////////////////////////////////////////////////////////////////////////
  //
  //  WDutils::RunInfo
  //
  /// provides information about the running process
  ///
  /// only one object exists, the static RunInfo::Info
  ///
  // ///////////////////////////////////////////////////////////////////////////
  class RunInfo {
  private:
    bool __host_known;
    bool __user_known;
    bool __pid_known;
    bool __name_known;
    bool __cmd_known;
    bool __is_mpi_proc;
    char __time   [100];
    char __host   [100];
    char __user   [100];
    char __pid     [20];
    char __name   [100];
    char __cmd   [1024];
    int  __debug;
    int  __mpi_proc;
    int  __mpi_size;
    RunInfo();
    static RunInfo Info;
  public:
    /// reset the debugging level
    static void set_debug_level(int d) { Info.__debug = d; }
    /// provide info about MPI
    static void set_mpi_proc(int p, int s) {
      Info.__is_mpi_proc = 1;
      Info.__mpi_proc=p;
      Info.__mpi_size=s;
    }
    /// is host name known?
    static bool const&host_known() { return Info.__host_known; }
    /// is user name known?
    static bool const&user_known() { return Info.__user_known; }
    /// is user pid known?
    static bool const&pid_known() { return Info.__pid_known; }
    /// is name of the running program known?
    static bool const&name_known() { return Info.__name_known; }
    /// is command line is known?
    static bool const&cmd_known() { return Info.__cmd_known; }
    /// is this process part of an MPI run?
    static bool const&is_mpi_proc() { return Info.__is_mpi_proc; }
    /// string with full time of run
    static const char*time() { return Info.__time; }
    /// string with host nam
    static const char*host() { return Info.__host; }
    /// string with user name
    static const char*user() { return Info.__user; }
    /// string with user pid
    static const char*pid() { return Info.__pid; }
    /// string with name of the running program
    static const char*name() { return Info.__name; }
    /// string with command line
    static const char*cmd() { return Info.__cmd; }
    /// return debugging level
    static const int &debug_level() { return Info.__debug; }
    /// return our rank within MPI, if we are part of a MPI run
    static const int &mpi_proc() { return Info.__mpi_proc; }
    /// return our size of MPI::World, if we are part of a MPI run
    static const int &mpi_size() { return Info.__mpi_size; }
    /// return true if debug level >= given debug depth
    static bool debug(int depth) { return Info.__debug >= depth; }
    /// print a log-file header
    static void header(std::ostream&out);
  };
  /// is debugging level exceeded by debugging depth (argument)?                
  /// \param d debugging depth
  inline bool debug(int d) {
    return RunInfo::debug(d);
  }
  // ///////////////////////////////////////////////////////////////////////////
  /// \name print debugging information to stderr, reporting [file:line]        
  //@{                                                                          
  /// to be used for reporting debug info
  struct DebugInformation {
    const char*file, *lib;    ///< file and library name
    const int  line;          ///< line number
    /// constructor: get library name
    DebugInformation(const char*__lib = "WDutils")
      : file(0), lib(__lib), line(0) {}
    /// constructor: get file name, line number, and library name
    DebugInformation(const char*__file, int __line, const char*__lib= "WDutils")
      : file(__file), lib(__lib), line(__line) {}
    /// print info message to stderr, report [file:line] if known.
    /// \param[in] fmt debug info message (C-type format string)
    /// \param[in] ... data to be formated
    void operator() (const char*fmt, ...) const
#ifdef __GNUC__
      __attribute__ ((format (printf, 2, 3)))
#endif
      ;
    /// print info message to stderr, report [file:line] if known.
    /// \param[in] lev level: only report if less than debug_level()
    /// \param[in] fmt debug info message (C-type format string)
    /// \param[in] ... data to be formated
    void operator() (int lev, const char*fmt, ...) const
#ifdef __GNUC__
      __attribute__ ((format (printf, 3, 4)))
#endif
      ;
  };
  /// print debug info to stderr and report [file:line]).
  /// use like NEMO's debug_info(), i.e. with EXACTLY the same syntax:
  /// \code
  /// void DebugInfo(int debug_level, const char*format, ...);
  /// void DebugInfo(const char*format, ...); 
  /// \endcode
#define DebugInfo  WDutils::DebugInformation(__FILE__,__LINE__)
  /// print debug info to stderr (without reporting [file:line]).
  /// use like NEMO's debug_info(), i.e. with EXACTLY the same syntax.
  /// \code
  /// void DebugInfoN(int debug_level, const char*format, ...);
  /// void DebugInfoN(const char*format, ...); 
  /// \endcode
#define DebugInfoN WDutils::DebugInformation()
  //@}
  // ///////////////////////////////////////////////////////////////////////////
  /// \name exception treatment                                                 
  //@{                                                                          
  /// simple exception with error message
  struct exception : protected std::string {
    /// copy constructor
    exception(exception const&e) : std::string(e) {}
    /// construction from C-style format string + data.
    /// Uses a printf() style format string as first argument, further arguments
    /// must match format, exactly as in printf, which will be called.
    /// \param[in] fmt gives the format in C printf() style
    explicit exception(const char*fmt, ...);
    /// return error message 
    const char*text() const { return c_str(); }
  };
  /// return error message given an exception
  inline const char*text(exception const&e) { return e.text(); }
  /// for generating exceptions
  struct Thrower {
    const char*file;          ///< file name
    const int  line;          ///< line number
    /// default constructor: set data to NULL
    Thrower() : file(0), line(0) {}
    /// constructor: get file name, and line number
    Thrower(const char*__file, int __line) : file(__file), line(__line) {}
    /// generate an exception
    /// \param[in] fmt  gives the format in C printf() style
    exception operator()(const char*fmt, ...) const
#ifdef __GNUC__
      __attribute__ ((format (printf, 2, 3)))
#endif
      ;
  };
  //@}
  // ///////////////////////////////////////////////////////////////////////////
  /// C++ wrapper around a C string.                                            
  /// Construction from C-type format string + data;                            
  /// Type conversion to const char*                                            
  /// Useful for generating a C-style string containing formatted data.         
  class message {
    message(message const&);                       // no copy constructor       
    static const size_t size = 1024;
    char __text[size];
  public:
    /// Generate a string from format + data.
    /// Uses a printf() style format string as first argument, further arguments
    /// must match format, exactly as in printf, which will be called.
    /// \param fmt gives the format in C printf() style
    explicit message(const char* fmt, ...) throw(exception);
    /// conversion to C-style string
    operator const char*() const { return __text; }
    /// return C-style string
    const char* text() const { return __text; }
  };
  // ///////////////////////////////////////////////////////////////////////////
  /// \name error treatment (alternative to throwing an exception)              
  //@{                                                                          
  /// to be used for error reporting
  struct Error {
    const char*file, *lib;    ///< file and library name
    const int  line;          ///< line number
    /// constructor: get library name
    Error(const char*__lib = "WDutils")
      : file(0), lib(__lib), line(0) {}
    /// constructor: get file name, line number, and library name
    Error(const char*__file, int __line, const char*__lib = "WDutils")
      : file(__file), lib(__lib), line(__line) {}
    /// print error message to stderr, report [file:line] if known.
    void operator() (const char*fmt, ...) const
#ifdef __GNUC__
      __attribute__ ((noreturn, format (printf, 2, 3)))
#endif
      ;
  };
  /// print error message to stderr, reporting [file:line], and exit.
  /// use like NEMO's error(), i.e. with the same syntax:
  /// \code
  /// void WDutils_Error(const char*format, ...);
  /// \endcode
#define WDutils_Error      WDutils::Error(__FILE__,__LINE__)
  /// print error message to stderr and exit.
  /// use like NEMO's error(), i.e. with the same syntax:
  /// \code
  /// void WDutils_ErrorN(const char*format, ...);
  /// \endcode
#define WDutils_ErrorN     WDutils::Error()
  //@}
  // ///////////////////////////////////////////////////////////////////////////
  /// \name warning treatment                                                   
  //@{                                                                          
  /// to be used for warning reporting
  struct Warning {
    const char*file, *lib;    ///< file and library name
    const int  line;          ///< line number
    /// constructor: get library name
    Warning(const char*__lib = "WDutils")
      : file(0), lib(__lib), line(0) {}
    /// constructor: get file name, line number, and library name
    Warning(const char*__file, int __line, const char*__lib = "WDutils")
      : file(__file), lib(__lib), line(__line) {}
    /// print error message to stderr, report [file:line] if known.
    void operator() (const char*fmt, ...) const
#ifdef __GNUC__
      __attribute__ ((format (printf, 2, 3)))
#endif
      ;
  };
  /// print warning message to stderr, reporting [file:line].
  /// use like NEMO's warning(), i.e. with the same syntax:
  /// \code
  /// void WDutils_Warning(const char*format, ...);
  /// \endcode
#define WDutils_Warning	     WDutils::Warning(__FILE__,__LINE__)
  /// print warning message to stderr (without reporting [file:line]).
  /// use like NEMO's warning(), i.e. with the same syntax:
  /// \code
  /// void WDutils_WarningN(const char*format, ...);
  /// \endcode
#define WDutils_WarningN     WDutils::Warning()
  //@}
  // ///////////////////////////////////////////////////////////////////////////
  /// \name macros and code controling the usage of throw exception vs error    
  //@{                                                                          
#if 0
#  undef  WDutils_EXCEPTIONS
  /// use instead of <tt> throw(WDutils::exception) </tt> after function
  /// declaration
#  define WDutils_THROWING
  /// instead of throwing an exception: error with [file:line]
  /// use "WDutils_THROW(fmt, data)" instead of "error(fmt, data)" or "throw
  /// WDutils::exception(fmt, data)"
#  define WDutils_THROW	        WDutils_Error
  /// instead of throwing an exception: error 
  /// use "WDutils_THROW(fmt, data)" instead of "error(fmt, data)" or "throw
  /// WDutils::exception(fmt, data)"
#  define WDutils_THROWN	WDutils_ErrorN
#  define WDutils_THROWER       WDutils::Error
  //----------------------------------------------------------------------------
  /// use "WDutils_RETHROW(E)" to re-throw a caught exception "E"
#  define WDutils_RETHROW(E)    WDutils_Error  (text(E))
#else
#  define WDutils_EXCEPTIONS
#  define WDutils_THROWING      throw(WDutils::exception)
#  define WDutils_THROWER       throw WDutils::Thrower
#  define WDutils_THROWN        throw WDutils::exception
#  define WDutils_RETHROW(E)    throw E
#endif
#define WDutils_THROW  WDutils_THROWER(__FILE__,__LINE__)
  //@}
  // ///////////////////////////////////////////////////////////////////////////
  /// a safer snprintf.                                                         
  /// If the string size is too small or if a formatting error occurs,          
  /// an exception is thrown, which also gives source file name and line number.
  /// Otherwise the behaviour is identical to ANSI C99 snprintf().              
  /// \return bytes written                                                     
  /// \param str string to write into                                           
  /// \param size maximum number of bytes to write, including trailing 0.      
  /// \param fmt format string                                                  
  int snprintf(char*str, size_t size, const char* fmt, ... ) WDutils_THROWING;
  // ///////////////////////////////////////////////////////////////////////////
  struct snprintf__ {
    const char* file;
    const int   line;
    snprintf__(const char*f, int l) : file(f), line(l) {}
    int operator() (char*str, size_t size, const char* fmt, ... )
      WDutils_THROWING
#ifdef __GNUC__
      __attribute__ ((format (printf, 4, 5)))
#endif
      ;
  };
  /// macro WDutilsSNprintf to be used like snprintf.
  /// reports source file and line on error.
#define SNprintf WDutils::snprintf__(__FILE__,__LINE__)

  // ///////////////////////////////////////////////////////////////////////////
} // namespace WDutils
// /////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_exception_h
