// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/exception.h
///
/// \author Walter Dehnen
///                                                                             
/// \date   2000-2012
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2012 Walter Dehnen
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

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || (__cplusplus >= 201103L)
#  define WDutilsCXX11
#else
#  undef  WDutilsCXX11
#endif

#if __cplusplus >= 201103L
#  define WDutilsCXX11Delete      = delete
#  define WDutilsCXX11Default     = default
#  define WDutilsCXX11DefaultBody = default;
#else
#  define WDutilsCXX11Delete
#  define WDutilsCXX11Default
#  define WDutilsCXX11DefaultBody {}
#endif

#ifndef WDutils_included_string
#  define WDutils_included_string
#  include <string>
#endif
#ifndef WDutils_included_limits
#  define WDutils_included_limits
#  include <limits>
#endif
#ifndef WDutils_included_cstdlib
#  define WDutils_included_cstdlib
#  include <cstdlib>
#endif
#ifdef WDutilsCXX11
# ifndef WDutils_included_type_traits
#  include <type_traits>
# endif
#endif

//                                                                              
//  WDutils                                                                     
//                                                                              
/// generally useful code of Walter Dehnen, used in project falcON,             
/// public under the GNU public licence                                         
///                                                                             
namespace WDutils {
#ifdef WDutilsCXX11
  using std::is_same;
#else
  template<typename __T1, typename __T2> struct is_same
  { static const bool value = false; };
  template<typename __T> struct is_same<__T,__T>
  { static const bool value = true; };
#endif
  /// static type information: an extension of std::numeric_limits<>
  //
  namespace meta {
    /// \note The only additional member indicates a floating point type; @c
    ///       std::numeric_limits<>::is_integer cannot be used instead as it's
    ///       @c false for all non-fundamental types.
    template<typename __T> struct TypeInfo :
      public std::numeric_limits<__T> {
      static const bool is_floating_point = false;
    };
    template<> struct TypeInfo<float> :
      public std::numeric_limits<float> {
      static const bool is_floating_point = true;
    };
    template<> struct TypeInfo<double> :
      public std::numeric_limits<double> {
      static const bool is_floating_point = true;
    };
    template<> struct TypeInfo<long double> :
      public std::numeric_limits<long double> {
      static const bool is_floating_point = true;
    };
  }
  using meta::TypeInfo;
  /// provides information about the running process
  /// \note only one object exists, the static RunInfo::Info
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
    int  __pid_num;
    int  __debug;
    int  __mpi_proc;
    int  __mpi_size;
    int  __omp_size;
    int  __omp_proc;
#if defined(__unix) || defined(__DARWIN_UNIX03)
    long long __sec, __usec;
#elif defined(WIN32)
    long long __timecount;
    double __timetick;
#endif
    /// default ctor
    RunInfo();
    /// static Info
    static RunInfo Info;
    //  no copy ctor and no operator=
    RunInfo           (const RunInfo&) WDutilsCXX11Delete ;
    RunInfo& operator=(const RunInfo&) WDutilsCXX11Delete ;
  public:
    /// reset the debugging level
    static void set_debug_level(int d)
    { Info.__debug = d; }
    /// provide info about MPI
    static void set_mpi_proc(int p, int s)
    {
      Info.__is_mpi_proc = 1;
      Info.__mpi_proc=p;
      Info.__mpi_size=s;
    }
    /// maximum # processors available for openMP
    static int max_omp_proc()
    { return Info.__omp_proc; }
    /// # openMP threads to be used, may exceed @a max_omp_proc()
    /// \note defaults to max_omp_proc, implying openMP is used if available
    static int omp_threads()
    { return Info.__omp_size; }
    /// set # openMP threads
    /// \note If @a arg[0] == 't', we set # threads to # processors.\n
    ///       If @a arg[0] == 'f', we set # threads to 1 (no openMP).\n
    ///       Otherwise, we try to convert @a arg to an integer number and
    ///       take that. This may exceed the # processors.
    static void set_omp(const char*arg);
    /// set # openMP threads
    static void set_omp(int n);
    /// shall openMP parallelism be used?
    static bool use_omp()
    { return Info.__omp_size > 1; }
    /// is host name known?
    static bool host_known()
    { return Info.__host_known; }
    /// is user name known?
    static bool user_known()
    { return Info.__user_known; }
    /// is user pid known?
    static bool pid_known()
    { return Info.__pid_known; }
    /// is name of the running program known?
    static bool name_known()
    { return Info.__name_known; }
    /// is command line is known?
    static bool cmd_known()
    { return Info.__cmd_known; }
    /// string with full time of run
    static const char*time()
    { return Info.__time; }
    /// string with host nam
    static const char*host()
    { return Info.__host; }
    /// string with user name
    static const char*user()
    { return Info.__user; }
    /// string with process id
    static const char*pid()
    { return Info.__pid; }
    /// numerical valud of user process id
    static int const&pid_num()
    { return Info.__pid_num; }
    /// string with name of the running program
    static const char*name()
    { return Info.__name; }
    /// string with command line
    static const char*cmd()
    { return Info.__cmd; }
    /// return debugging level
    static int debug_level()
    { return Info.__debug; }
    /// is this process part of an MPI run?
    static bool is_mpi_proc()
    { return Info.__is_mpi_proc; }
    /// return our rank within MPI, if we are part of a MPI run
    static int mpi_proc()
    { return Info.__mpi_proc; }
    /// return our size of MPI::World, if we are part of a MPI run
    static int mpi_size()
    { return Info.__mpi_size; }
    /// return true if debug level >= given debug depth
    static bool debug(int depth)
    { return Info.__debug >= depth; }
    /// print a log-file header
    static void header(std::ostream&out);
#if defined(__unix) || defined(__DARWIN_UNIX03) || defined(WIN32)
    /// time in seconds since start of the program (or more accurately since
    /// construction of the RunInfo object)
    /// \return time in seconds
    static double WallClock();
#endif
#if defined(__unix) || defined(__DARWIN_UNIX03)
    /// \param[out] sec   full seconds since start of the program
    /// \param[out] usec  micro seconds since start of the program
    static void WallClock(unsigned&sec, unsigned&usec);
#endif
  };
  /// is debugging level exceeded by debugging depth (argument)?
  /// \relates DebugInformation
  /// \param d debugging depth
  inline bool debug(int d)
  { return RunInfo::debug(d); }
  /* 
     Version 2.4 and later of GCC defines a magical variable
     `__PRETTY_FUNCTION__' which contains the name of the function currently
     being defined.  This is broken in G++ before version 2.6.  C9x has a
     similar variable called __func__, but the GCC macro is preferrable since
     it demangles C++ function names.
  */
#ifdef __GNUC__
#  if (__GNUC__ > 3) || ((__GNUC__ == 2) && (__GNUC_MINOR__ >= 6))
#    define WDutilsThisFunction	__PRETTY_FUNCTION__
#  elif defined __STDC_VERSION__ && __STDC_VERSION__ >= 199901L
#    define WDutilsThisFunction	__func__
#  endif
#else
#  define WDutilsThisFunction	0
#endif
  //
  /// \name print debugging information to stderr, reporting [file:line]        
  //@{
  /// to implement DebugInformation, Error, and Warning
  template<typename ReportTraits>
  struct Reporting {
    const char*library;
    const char*file,*func;      ///< names: file, function
    const int  line;            ///< line number
    /// constructor: get library and function name
    explicit Reporting(const char*__lib)
      : library(__lib), file(0), func(0), line(0) {}
    /// constructor: get library and function name
    explicit Reporting(const char*__func, const char*__lib)
      : library(__lib), file(0), func(__func), line(0) {}
    /// constructor: get file name, line number, and library name
    Reporting(const char*__file, int __line, const char*__lib)
      : library(__lib), file(__file), func(0), line(__line) {}
    /// constructor: get file name, func name, line number, and library name
    Reporting(const char*__func, const char*__file, int __line,
	      const char*__lib)
      : library(__lib), file(__file), func(__func), line(__line) {}
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
  private:
    //  no copy ctor and no operator=
    Reporting           (const Reporting&) WDutilsCXX11Delete;
    Reporting& operator=(const Reporting&) WDutilsCXX11Delete;
  };
  /// traits for DebugInformation
  struct DebugInfoTraits {
    static bool condition(int lev) { return RunInfo::debug(lev); }
    static const char*issue() { return "Debug Info"; }
    static void after() {}
  };
  typedef Reporting<DebugInfoTraits> DebugInformation;
  /// print debug info to stderr and report [file:line]func:
  /// use like NEMO's debug_info(), i.e. with EXACTLY the same syntax:
  /// \code
  ///   void DebugInfo(int debug_level, const char*format, ...);
  ///   void DebugInfo(const char*format, ...); 
  /// \endcode
#define DebugInfo \
  WDutils::DebugInformation(WDutilsThisFunction,__FILE__,__LINE__,"WDutils")
  /// print debug info to stderr and report func:
  /// use like NEMO's debug_info(), i.e. with EXACTLY the same syntax:
  /// \code
  ///   void DebugInfoF(int debug_level, const char*format, ...);
  ///   void DebugInfoF(const char*format, ...); 
  /// \endcode
#define DebugInfoF WDutils::DebugInformation(WDutilsThisFunction,"WDutils")
  /// print debug info to stderr (without reporting [file:line]).
  /// use like NEMO's debug_info(), i.e. with EXACTLY the same syntax.
  /// \code
  ///   void DebugInfoN(int debug_level, const char*format, ...);
  ///   void DebugInfoN(const char*format, ...); 
  /// \endcode
#define DebugInfoN WDutils::DebugInformation("WDutils")
  /// traits for Error
  struct ErrorTraits {
    static bool condition(int) { return true; }
    static const char*issue() { return "Error"; }
    static void after() { std::exit(1); }
  };
  typedef Reporting<ErrorTraits> Error;
  /// print error message to stderr, reporting [file:line]func, and exit.
  /// use like NEMO's error(), i.e. with the same syntax:
  /// \code
  ///   void WDutils_Error(const char*format, ...);
  /// \endcode
#define WDutils_Error \
  WDutils::Error(WDutilsThisFunction,__FILE__,__LINE__,"WDutils")
  /// print error message to stderr, reporting [file:line]func, and exit.
  /// use like NEMO's error(), i.e. with the same syntax:
  /// \code
  ///   void WDutils_ErrorF(const char*format, ...);
  /// \endcode
#define WDutils_ErrorF WDutils::Error(WDutilsThisFunction,"WDutils")
  /// print error message to stderr and exit.
  /// use like NEMO's error(), i.e. with the same syntax:
  /// \code
  ///   void WDutils_ErrorN(const char*format, ...);
  /// \endcode
#define WDutils_ErrorN WDutils::Error("WDutils")
  /// traits for Warning
  struct WarningTraits {
    static bool condition(int) { return true; }
    static const char*issue() { return "Warning"; }
    static void after() {}
  };
  typedef Reporting<WarningTraits> Warning;
  /// print warning message to stderr, reporting [file:line]func
  /// use like NEMO's warning(), i.e. with the same syntax:
  /// \code
  ///   void WDutils_Warning(const char*format, ...);
  /// \endcode
#define WDutils_Warning \
  WDutils::Warning(WDutilsThisFunction,__FILE__,__LINE__,"WDutils")
  /// print warning message to stderr, reporting [file:line]func
  /// use like NEMO's warning(), i.e. with the same syntax:
  /// \code
  ///   void WDutils_WarningF(const char*format, ...);
  /// \endcode
#define WDutils_WarningF WDutils::Warning(WDutilsThisFunction,"WDutils")
  /// print warning message to stderr
  /// use like NEMO's warning(), i.e. with the same syntax:
  /// \code
  ///   void WDutils_WarningN(const char*format, ...);
  /// \endcode
#define WDutils_WarningN WDutils::Warning("WDutils")
  /// \name exception treatment                                                 
  //@{                                                                          
  /// simple exception with error message
  struct exception : protected std::string {
    /// copy constructor
#if __cplusplus >= 201103L
    exception(exception const&) = default;
#else
    exception(exception const&e)
      : std::string(e) {}
#endif
    /// construction from C-style format string + data.
    /// Uses a printf() style format string as first argument, further arguments
    /// must match format, exactly as in printf, which will be called.
    /// \param[in] fmt gives the format in C printf() style
    explicit exception(const char*fmt, ...);
    /// return error message 
    const char*text() const
    { return c_str(); }
  };
  /// return error message given an exception
  inline const char*text(exception const&e)
  { return e.text(); }
  //
  class ThrowGuard;
  /// for generating exceptions
  class Thrower {
    friend class ThrowGuard;
    typedef void(*handler)(const char*,const char*, int, const char*);
    static handler  InsteadOfThrow;  ///< make error if OMP::IsParallel()
    const  char    *file,*func;      ///< file name, function name
    const  int      line;            ///< line number
    //  no copy ctor and no operator=
    Thrower           (const Thrower&) WDutilsCXX11Delete;
    Thrower& operator=(const Thrower&) WDutilsCXX11Delete;
  public:
    /// default constructor: set data to NULL
    Thrower()
      : file(0), func(0), line(0) {}
    /// constructor: get function name
    explicit Thrower(const char*__func)
      : file(0), func(__func), line(0) {}
    /// constructor: get file name, and line number
    Thrower(const char*__file, int __line)
      : file(__file), func(0), line(__line) {}
    /// constructor: get file & function name, and line number
    Thrower(const char*__func, const char*__file, int __line)
      : file(__file), func(__func), line(__line) {}
    /// generate an exception; for usage in WDutils_THROW
    /// \param[in] fmt  gives the format in C printf() style
    exception operator()(const char*fmt, ...) const
#ifdef __GNUC__
      __attribute__ ((format (printf, 2, 3)))
#endif
      ;
    /// generate an exception; for usage in WDutilsAssert
    /// \param[in] expr  boolean expression: throw exception if false
    exception operator()(bool expr) const;
  };
  /// method invoking an error, suitable as @a Thrower::handler
  inline void MakeError(const char*func, const char*file, int line,
			const char*mess)
  { Error(func,file,line,"WDutils")(mess); }
  /// guard against throwing an exception inside an openMP parallel region
  struct ThrowGuard
  {
    /// ctor: guard against throwing of @c WDutils::exception via @c
    /// WDutils::Thrower inside omp parallel regions
    ThrowGuard() : OldHandler(Thrower::InsteadOfThrow)
    { Thrower::InsteadOfThrow = WDutils::MakeError; }
    /// dtor: replace handler with original
    ~ThrowGuard()
    { Thrower::InsteadOfThrow = OldHandler; }
  protected:
    explicit ThrowGuard(Thrower::handler H)
      : OldHandler(Thrower::InsteadOfThrow)
    { Thrower::InsteadOfThrow = H; }
  private:
    const Thrower::handler OldHandler;
    //  no copy ctor and no operator=
    ThrowGuard           (const ThrowGuard&) WDutilsCXX11Delete;
    ThrowGuard& operator=(const ThrowGuard&) WDutilsCXX11Delete;
  };
  //@}
  //
  /// C++ wrapper around a C string.                                            
  /// Construction from C-type format string + data;                            
  /// Type conversion to const char*                                            
  /// Useful for generating a C-style string containing formatted data.         
  class message {
    static const size_t size = 1024;
    char __text[size];
    //  no copy ctor and no operator=
    message           (const message&) WDutilsCXX11Delete;
    message& operator=(const message&) WDutilsCXX11Delete;
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
  /// \name macros and code controling the usage of throw exception vs error    
  //@{                                                                          
#if 0
#  undef  WDutils_EXCEPTIONS
  /// use instead of <tt> throw(WDutils::exception) </tt> after function
  /// declaration
#  define WDutils_THROWING
  /// instead of throwing an exception: error 
  /// use "WDutils_THROW(fmt, data)" instead of "error(fmt, data)" or "throw
  /// WDutils::exception(fmt, data)"
#  define WDutils_THROWN	WDutils_ErrorN
  /// instead of throwing an exception: error with [file:line]
  /// use "WDutils_THROW(fmt, data)" instead of "error(fmt, data)" or "throw
  /// WDutils::exception(fmt, data)"
#  define WDutils_THROWER       WDutils::Error
  //----------------------------------------------------------------------------
  /// use "WDutils_RETHROW(E)" to re-throw a caught exception "E"
#  define WDutils_RETHROW(E)    WDutils_Error (text(E))
#else // 0/1
#  define WDutils_EXCEPTIONS
  /// use instead of <tt> throw(WDutils::exception) </tt> after function
  /// declaration
#  define WDutils_THROWING      throw(WDutils::exception)
#  define WDutils_THROWER       throw WDutils::Thrower
#  define WDutils_THROWN        throw WDutils::exception
#  define WDutils_RETHROW(E)    throw E
#endif
  /// use to report an error like <tt> WDutils_THROW("x=%f<0",x); </tt>
#define WDutils_THROW  \
  WDutils_THROWER(WDutilsThisFunction,__FILE__,__LINE__)
  /// use to report an error like <tt> WDutils_THROW("x=%f<0",x); </tt>
#define WDutils_THROWF WDutils_THROWER(WDutilsThisFunction)
  //@}
  //
  //  macro for compile-time assertion, stolen from the boost library
  //
#ifdef WDutilsCXX11
  /// \brief macro for compile-time assertion
#  define WDutilsStaticAssert(TEST) static_assert(TEST,#TEST);
#else
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
    static_cast<bool>((TEST))>)					\
  }
#endif
  /// \name assertion which throws an excpetion rather than abort
  //@{
  // is NDEBUG is defined, do nothing
#ifdef  NDEBUG
# define WDutilsAssert(expr)   (static_cast<void>(0))
# define WDutilsAssertE(expr)  (static_cast<void>(0))
#else
  /// throws exception with "assertion failed" message
# ifdef __GNUC__
  inline
  void AssertFail(const char*, const char*, const char*, unsigned)
    WDutils_THROWING __attribute__ ((__noreturn__));
# endif
  inline
  void AssertFail(const char*assertion, const char*func, 
		  const char*file, unsigned line) WDutils_THROWING
  { WDutils_THROWER(func,file,line)("Assertion \"%s\" failed",assertion); }
  //
# ifdef __GNUC__
  inline
  void AssertFailE(const char*, const char*, const char*, unsigned)
    __attribute__ ((__noreturn__));
# endif
  inline
  void AssertFailE(const char*assertion, const char*func,
		   const char*file, unsigned line)
  { 
    WDutils::Error(func,file,line,"WDutils")
      ("Assertion \"%s\" failed",assertion);
    std::exit(1);
  }
  /// use instead of assert(): throws an exception
# define WDutilsAssert(expr)						\
  ((expr)								\
  ? static_cast<void>(0)						\
  : WDutils::AssertFail(__STRING(expr),WDutilsThisFunction,__FILE__,__LINE__))
  /// almost identical to assert()
# define WDutilsAssertE(expr)						\
  ((expr)								\
  ? static_cast<void>(0)						\
  : WDutils::AssertFailE(__STRING(expr),WDutilsThisFunction,__FILE__,__LINE__))
#endif // NDEBUG
  //@}
  //
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
  /// macro to be used like snprintf.
  /// reports source file and line on error.
#define SNprintf WDutils::snprintf__(__FILE__,__LINE__)

  // ///////////////////////////////////////////////////////////////////////////
} // namespace WDutils
// /////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_exception_h
