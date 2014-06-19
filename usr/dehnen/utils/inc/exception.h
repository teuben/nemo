// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/exception.h
///
/// \author Walter Dehnen
///                                                                             
/// \date   2000-14
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-14 Walter Dehnen
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
#if __cplusplus >= 201103L
#  ifndef WDutils_included_iostream
#    define WDutils_included_iostream
#    include <iostream>
#  endif
#  ifndef WDutils_included_sstream
#    define WDutils_included_sstream
#    include <sstream>
#  endif
#endif
#ifndef WDutils_included_limits
#  define WDutils_included_limits
#  include <limits>
#endif
#ifndef WDutils_included_cstdlib
#  define WDutils_included_cstdlib
#  include <cstdlib>
#endif
#ifndef WDutils_included_stdexcept
#  define WDutils_included_stdexcept
#  include <stdexcept>
#endif
#if __cplusplus >= 201103L
#  ifndef WDutils_included_type_traits
#    include <type_traits>
#  endif
#endif
#if __cplusplus >= 201103L && defined(_OPENMP) && \
  !defined(WDutils_included_omp_h)
#  include <omp.h>
#  define WDutils_included_omp_h
#endif

//                                                                              
//  WDutils                                                                     
//                                                                              
/// generally useful code of Walter Dehnen, used in project falcON,             
/// public under the GNU public licence                                         
///                                                                             
namespace WDutils {

  /// provides information about the running process
  /// \note only one object exists, the static RunInfo::Info
  class RunInfo 
  {
  private:
    bool _m_host_known;
    bool _m_user_known;
    bool _m_pid_known;
    bool _m_name_known;
    bool _m_cmd_known;
    bool _m_is_mpi_proc;
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wunused-private-field"
#endif
    bool _m_dummy_bool[2]; // padding to 8 bytes
#ifdef __clang__
#  pragma clang diagnostic pop
#endif
    char _m_time   [104];
    char _m_host   [104];
    char _m_user   [104];
    char _m_pid     [24];
    char _m_name   [104];
    char _m_cmd   [1024];
    int       _m_pid_num,_m_debug;
    int       _m_mpi_proc,_m_mpi_size;
    int       _m_omp_proc,_m_omp_size;
    unsigned  _m_tbb_proc,_m_tbb_size;
    void     *_m_tbb_init;
#if defined(__unix) || defined(__DARWIN_UNIX03)
    long long _m_sec, _m_usec;
#elif defined(WIN32)
    long long _m_timecount;
    double    _m_timetick;
#endif
    /// default constructor
    RunInfo();
    /// destructor
    ~RunInfo();
    /// static Info
    static RunInfo Info;
    //  no copy ctor and no operator=
    RunInfo           (const RunInfo&) WDutilsCXX11Delete ;
    RunInfo& operator=(const RunInfo&) WDutilsCXX11Delete ;
  public:
    /// reset the debugging level
    static void set_debug_level(int d)
    { Info._m_debug = d; }
    /// provide info about MPI
    static void set_mpi_proc(int p, int s)
    {
      Info._m_is_mpi_proc = 1;
      Info._m_mpi_proc=p;
      Info._m_mpi_size=s;
    }

    /// \name openMP stuff
    //@{
    /// set \# openMP threads
    /// \note If @a arg[0] == 't', we set \# threads to \# processors.
    ///       If @a arg[0] == 'f', we set \# threads to 1 (no openMP).
    ///       Otherwise, we try to convert @a arg to an integer number and
    ///       take that. This may exceed the \# processors.
    static void set_omp(const char*arg);
    /// set number of openMP threads
    static void set_omp(int n);
    /// maximum \# processors available for openMP
    static int max_omp_proc()
    { return Info._m_omp_proc; }
    /// number of openMP threads to be used, may exceed @a max_omp_proc()
    /// \note defaults to max_omp_proc, implying openMP is used if available
    static int omp_threads()
    { return Info._m_omp_size; }
    /// shall openMP parallelism be used?
    static bool use_omp()
    { return Info._m_omp_size > 1; }
    //@}

    /// \name TBB stuff
    //@{
#ifdef WDutilsTBB
#  define _END_FUNCTION_
#else
#  define _END_FUNCTION_ {}
#endif
    /// set \# TBB threads
    /// \note If @a arg[0] == 't', we set \# threads to automatic
    ///       If @a arg[0] == 'f', we set \# threads to 1.
    ///       Otherwise, we try to convert @a arg to an integer number and
    ///       take that. This may exceed the \# processors.
    static void set_tbb(const char*) _END_FUNCTION_;
    /// set \# TBB threads, if arg=0, uses tbb::task_scheduler_init::automatic
    /// \note use n=1 to de-activate tbb parallelism
    static void set_tbb(unsigned =0) _END_FUNCTION_;
#undef _END_FUNCTION_
    /// use tbb
    static bool use_tbb()
    { return Info._m_tbb_size > 1; }
    /// maximum \# processors available for TBB
    static unsigned max_tbb_proc()
    { return Info._m_tbb_proc; }
    /// \# TBB threads to be used, may exceed @a max_tbb_proc()
    static unsigned tbb_threads()
    { return Info._m_tbb_size; }
    /// is TBB parallelism activated through set_tbb()?
    /// \note you can use TBB parallelism even if tbb_is_active() returns false
    static bool tbb_is_active()
    { return tbb_threads() > 0; }
    //@}
#if __cplusplus >= 201103L && defined(WDutilsDevel)
    /// return a unique short id for the current std::thread
    static unsigned thread_id();
#endif
    /// is host name known?
    static bool host_known()
    { return Info._m_host_known; }
    /// is user name known?
    static bool user_known()
    { return Info._m_user_known; }
    /// is user pid known?
    static bool pid_known()
    { return Info._m_pid_known; }
    /// is name of the running program known?
    static bool name_known()
    { return Info._m_name_known; }
    /// is command line is known?
    static bool cmd_known()
    { return Info._m_cmd_known; }
    /// string with full time of run
    static const char*time()
    { return Info._m_time; }
    /// string with host nam
    static const char*host()
    { return Info._m_host; }
    /// string with user name
    static const char*user()
    { return Info._m_user; }
    /// string with process id
    static const char*pid()
    { return Info._m_pid; }
    /// numerical valud of user process id
    static int const&pid_num()
    { return Info._m_pid_num; }
    /// string with name of the running program
    static const char*name()
    { return Info._m_name; }
    /// string with command line
    static const char*cmd()
    { return Info._m_cmd; }
    /// return debugging level
    static int debug_level()
    { return Info._m_debug; }
    /// is this process part of an MPI run?
    static bool is_mpi_proc()
    { return Info._m_is_mpi_proc; }
    /// return our rank within MPI, if we are part of a MPI run
    static int mpi_proc()
    { return Info._m_mpi_proc; }
    /// return our size of MPI::World, if we are part of a MPI run
    static int mpi_size()
    { return Info._m_mpi_size; }
    /// return true if debug level >= given debug depth
    static bool debug(int depth)
    { return Info._m_debug >= depth; }
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
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif
  /// is debugging level exceeded by debugging depth (argument)?
  /// \relates DebugInformation
  /// \param d debugging depth
#ifdef __clang__
#  pragma clang diagnostic pop
#endif
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
#if __cplusplus >= 201103L
  inline std::ostringstream&make_ostr(std::ostringstream&ostr) noexcept
  {
    return ostr;
  }
  template<class T, class... R>
  inline std::ostringstream&make_ostr(std::ostringstream&ostr,
				      T const&x, R&&... r)
  {
    ostr<<x;
    return make_ostr(ostr, std::forward<R>(r)...);
  }
  /// make a C++ string from any number of arguments using C++ formatted output
  template<class... Args>
  inline std::string make_string(Args&&... args)
  {
    std::ostringstream ostr;
    make_ostr(ostr,std::forward<Args>(args)...);
    return std::move(ostr.str());
  }
  inline std::string make_string(const char*text)
  {
    return {text};
  }
  inline std::string make_string(std::string const&text)
  {
    return {text};
  }
#endif // C++11
  //
  /// \name print debugging information to stderr, reporting [file:line]        
  //@{
  /// to implement DebugInformation, Error, and Warning
  template<typename ReportTraits>
  struct Reporting {
    const char    *library;
    const char    *file,*func;      ///< names: file, function
    const unsigned line;            ///< line number
    const unsigned flag;            ///< currently only: write thread id?
    /// constructor: get library and function name
    explicit Reporting(const char*_m_lib, unsigned _m_flag=1)
      : library(_m_lib), file(0), func(0), line(0), flag(_m_flag) {}
    /// constructor: get library and function name
    explicit Reporting(const char*_m_func, const char*_m_lib,
		       unsigned _m_flag=1)
      : library(_m_lib), file(0), func(_m_func), line(0), flag(_m_flag) {}
    /// constructor: get file name, line number, and library name
    Reporting(const char*_m_file, unsigned _m_line, const char*_m_lib,
	      unsigned _m_flag=1)
      : library(_m_lib), file(_m_file), func(0), line(_m_line), flag(_m_flag) {}
    /// constructor: get file name, func name, line number, and library name
    Reporting(const char*_m_func, const char*_m_file, unsigned _m_line,
	      const char*_m_lib, unsigned _m_flag=1)
      : library(_m_lib), file(_m_file), func(_m_func), line(_m_line)
      , flag(_m_flag) {}
    /// print info message to stderr, report [file:line] if known.
    /// \param[in] fmt debug info message (C-type format string)
    /// \note  data to be formated must obey the usual C-style formatting rules
    void operator() (const char*fmt, ...) const
#ifdef __GNUC__
      __attribute__ ((format (printf, 2, 3)))
#endif
      ;
    /// print info message to stderr, report [file:line] if known.
    /// \param[in] lev level: only report if less than debug_level()
    /// \param[in] fmt debug info message (C-type format string)
    /// \note  data to be formated must obey the usual C-style formatting rules
    void operator() (int lev, const char*fmt, ...) const
#ifdef __GNUC__
      __attribute__ ((format (printf, 3, 4)))
#endif
      ;
    /// report message constructed from any number of arguments
#if __cplusplus >= 201103L
    template<class... Args>
    void report(Args&&... args) const
    {
      std::ostringstream ostr;
      print_header(ostr);
      make_ostr(ostr,std::forward<Args>(args)...);
      std::clog << ostr.str() << std::endl;      
    }
    template<class... Args>
    void debug(int level, Args&&... args) const
    {
      if(ReportTraits::condition(level)) {
	std::ostringstream ostr;
	print_header(ostr,level);
	make_ostr(ostr,std::forward<Args>(args)...);
	std::clog << ostr.str() << std::endl;      
      }
    }
#endif
  private:
#if __cplusplus >= 201103L
    /// print header of report() output
    void print_header(std::ostringstream&, int=0) const;
#endif
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
  WDutils::DebugInformation(__FILE__,__LINE__,"WDutils")
  /// print debug info to stderr and report func:
  /// use like NEMO's debug_info(), i.e. with EXACTLY the same syntax:
  /// \code
  ///   void DebugInfoF(int debug_level, const char*format, ...);
  ///   void DebugInfoF(const char*format, ...); 
  /// \endcode
#define DebugInfoF \
  WDutils::DebugInformation(WDutilsThisFunction,__FILE__,__LINE__,"WDutils")
  /// print debug info to stderr (without reporting [file:line]).
  /// use like NEMO's debug_info(), i.e. with EXACTLY the same syntax.
  /// \code
  ///   void DebugInfoN(int debug_level, const char*format, ...);
  ///   void DebugInfoN(const char*format, ...); 
  /// \endcode
#define DebugInfoN  WDutils::DebugInformation("WDutils")
#define DebugInfoN0 WDutils::DebugInformation("WDutils",0u)
#if __cplusplus >= 201103L
#  define DebugInfo11  \
  WDutils::DebugInformation(__FILE__,__LINE__,"WDutils").debug
#  define DebugInfo11F \
  WDutils::DebugInformation(WDutilsThisFunction,__FILE__,__LINE__,\
			    "WDutils").debug
#  define DebugInfo11N \
  WDutils::DebugInformation("WDutils").debug
#endif
  /// traits for Error
  struct ErrorTraits {
    static bool condition(int) { return true; }
    static const char*issue() { return "Error"; }
    static __attribute__((noreturn)) void after() { std::terminate(); }
  };
  typedef Reporting<ErrorTraits> Error;
  /// print error message to stderr, reporting [file:line]func, and exit.
  /// use like NEMO's error(), i.e. with the same syntax:
  /// \code
  ///   void WDutils_Error(const char*format, ...);
  /// \endcode
#define WDutils_Error \
  WDutils::Error(__FILE__,__LINE__,"WDutils")
  /// print error message to stderr, reporting [file:line]func, and exit.
  /// use like NEMO's error(), i.e. with the same syntax:
  /// \code
  ///   void WDutils_ErrorF(const char*format, ...);
  /// \endcode
#define WDutils_ErrorF \
  WDutils::Error(WDutilsThisFunction,__FILE__,__LINE__,"WDutils")
  /// print error message to stderr and exit.
  /// use like NEMO's error(), i.e. with the same syntax:
  /// \code
  ///   void WDutils_ErrorN(const char*format, ...);
  /// \endcode
#define WDutils_ErrorN WDutils::Error("WDutils")
#if __cplusplus >= 201103L
#  define WDutilsError11  \
  WDutils::Error(__FILE__,__LINE__,"WDutils").report
#  define WDutilsError11F \
  WDutils::Error(WDutilsThisFunction,__FILE__,__LINE__,"WDutils").report
#  define WDutilsError11N \
  WDutils::Error("WDutils").report
#endif
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
#define WDutils_Warning WDutils::Warning(__FILE__,__LINE__,"WDutils")
  /// print warning message to stderr, reporting [file:line]func
  /// use like NEMO's warning(), i.e. with the same syntax:
  /// \code
  ///   void WDutils_WarningF(const char*format, ...);
  /// \endcode
#define WDutils_WarningF \
  WDutils::Warning(WDutilsThisFunction,__FILE__,__LINE__,"WDutils")
  /// print warning message to stderr
  /// use like NEMO's warning(), i.e. with the same syntax:
  /// \code
  ///   void WDutils_WarningN(const char*format, ...);
  /// \endcode
#define WDutils_WarningN WDutils::Warning("WDutils")
#if __cplusplus >= 201103L
#  define WDutilsWarning11  \
  WDutils::Warning(__FILE__,__LINE__,"WDutils").report
#  define WDutilsWarning11F \
  WDutils::Warning(WDutilsThisFunction,__FILE__,__LINE__,"WDutils").report
#  define WDutilsWarning11N \
  WDutils::Warning("WDutils").report
#endif
  /// \name exception treatment                                                 
  //@{                                                                          
  /// simple exception with error message
  struct exception : std::runtime_error
  {
    /// construction from std::runtime_error
    explicit exception(std::runtime_error const&e)
      : std::runtime_error(e) {}
    /// construction from string
    explicit exception(std::string const&s)
      : std::runtime_error(s) {}
    /// construction from C-style format string + data.
    /// Uses a printf() style format string as first argument, further arguments
    /// must match format, exactly as in printf, which will be called.
    /// \param[in] fmt gives the format in C printf() style
    explicit exception(const char*fmt, ...);
    ///
    using std::runtime_error::what;
  };
  /// make an exception from any number of arguments
#if __cplusplus >= 201103L
  template<class... Args>
  exception make_exception(Args&&... args)
  {
    return exception(make_string(std::forward<Args>(args)...));
  }
#endif
  /// return error message given an exception
  inline const char*text(exception const&e)
  { return e.what(); }
  //
  struct ThrowGuard;
  /// for generating exceptions
  class Thrower {
    friend struct ThrowGuard;
    typedef void(*handler)(const char*, unsigned, const char*);
    static handler  InsteadOfThrow;  ///< make error if OMP::IsParallel()
    const  char    *file,*func;      ///< file name, function name
    const  unsigned line;            ///< line number
    //  no copy ctor and no operator=
    Thrower           (const Thrower&) WDutilsCXX11Delete;
    Thrower& operator=(const Thrower&) WDutilsCXX11Delete;
  public:
    static handler instead_of_throw() { return InsteadOfThrow; }
    /// default constructor: set data to NULL
    Thrower()
      : file(0), func(0), line(0) {}
    /// constructor: get function name
    explicit Thrower(const char*_m_func)
      : file(0), func(_m_func), line(0) {}
    /// constructor: get file name, and line number
    Thrower(const char*_m_file, unsigned _m_line)
      : file(_m_file), func(0), line(_m_line) {}
    /// constructor: get file & function name, and line number
    Thrower(const char*_m_func, const char*_m_file, unsigned _m_line)
      : file(_m_file), func(_m_func), line(_m_line) {}
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
#if __cplusplus >= 201103L
    template<class... Args>
    exception throw_it(Args&&... args) const
    {
      const bool error = false
#  ifdef _OPENMP
	|| (omp_get_level() && InsteadOfThrow)
#  endif
	;
      std::ostringstream ostr;
      if(!error && file)
	ostr << '[' << file << ':' << line << "] ";
      if(func)
	ostr << "in " << func;
      make_ostr(ostr,std::forward<Args>(args)...);
#  ifdef _OPENMP
      if(error)
	InsteadOfThrow(file,line,ostr.str().c_str());
#  endif
      return exception(ostr.str());
    }
#endif
  };
  /// method invoking an error, suitable as @a Thrower::handler
  inline void MakeError(const char*file, unsigned line, const char*mess)
  {
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wformat-security"
#endif
    Error(file,line,"WDutils")(mess);
#ifdef __clang__
#  pragma clang diagnostic pop
#endif
  }
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
    char _m_text[size];
    //  no copy ctor and no operator=
    message           (const message&) WDutilsCXX11Delete;
    message& operator=(const message&) WDutilsCXX11Delete;
  public:
    /// Generate a string from format + data.
    /// Uses a printf() style format string as first argument, further arguments
    /// must match format, exactly as in printf, which will be called.
    /// \param fmt gives the format in C printf() style
    explicit message(const char* fmt, ...);
    /// conversion to C-style string
    operator const char*() const { return _m_text; }
    /// return C-style string
    const char* text() const { return _m_text; }
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
#  define WDutils_RETHROW(E)    WDutils_Error (E.what())
#else // 0/1
#  define WDutils_EXCEPTIONS
  /// use instead of <tt> throw(WDutils::exception) </tt> after function
  /// declaration
#  define WDutils_THROWING      throw(WDutils::exception)
#  define WDutils_THROWER       throw WDutils::Thrower
#  define WDutils_THROWN        throw WDutils::exception
#  define WDutils_RETHROW(E)    throw E
#  if __cplusplus >= 201103L
#    define WDutilsThrow11N     throw make_exception
#  endif
#endif
  /// use to report an error like <tt> WDutils_THROW("x=%f<0",x); </tt>
#define WDutils_THROW  \
  WDutils_THROWER(WDutilsThisFunction,__FILE__,__LINE__)
  /// use to report an error like <tt> WDutils_THROW("x=%f<0",x); </tt>
#define WDutils_THROWF WDutils_THROWER(WDutilsThisFunction)
#if __cplusplus >= 201103L
#  define WDutilsThrow11  \
  WDutils_THROWER(WDutilsThisFunction,__FILE__,__LINE__).throw_it
#  define WDutilsThrow11F \
  WDutils_THROWER(WDutilsThisFunction).throw_it
#endif
  //@}
  //
  //  macro for compile-time assertion, stolen from the boost library
  //
#if __cplusplus >= 201103L
  /// \brief macro for compile-time assertion
#  define WDutilsCXX11StaticAssert(TEST,MSGS) static_assert(TEST,MSGS);
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
#  define WDutilsStaticAssert(TEST)				\
  enum { __DUMMY = sizeof(WDutils::STATIC_ASSERTION_FAILURE<	\
    static_cast<bool>((TEST))>)					\
  }
#  define WDutilsCXX11StaticAssert(TEST,MSGS)	\
  WDutilsStaticAssert(TEST)
#endif
  /// \name assertion which throws an excpetion rather than abort
  //@{
  // is NDEBUG is defined, do nothing
#ifdef  NDEBUG
#  define WDutilsAssert(expr)          (static_cast<void>(0))
#  define WDutilsAssertE(expr)         (static_cast<void>(0))
#  define WDutilsAssertIf(cond,expr)   (static_cast<void>(0))
#  define WDutilsAssertEIf(cond,expr)  (static_cast<void>(0))
#else
  /// throws exception with "assertion failed" message
#  ifdef __GNUC__
  inline
  void AssertFail(const char*, const char*, unsigned)
    __attribute__ ((__noreturn__));
#  endif
  inline
  void AssertFail(const char*assertion, const char*file, unsigned line)
  { WDutils_THROWER(file,line)("Assertion \"%s\" failed",assertion); }
  //
#  ifdef __GNUC__
  inline
  void AssertFailE(const char*, const char*, unsigned)
    __attribute__ ((__noreturn__));
#  endif
  inline
  void AssertFailE(const char*assertion, const char*file, unsigned line)
  { 
    WDutils::Error(file,line,"WDutils")
      ("Assertion \"%s\" failed",assertion);
    std::exit(1);
  }
  /// use instead of assert(): throws an exception
#  define WDutilsAssert(expr)						\
  ((expr)								\
  ? static_cast<void>(0)						\
  : WDutils::AssertFail(__STRING(expr),__FILE__,__LINE__))
  /// almost identical to assert()
#  define WDutilsAssertE(expr)						\
  ((expr)								\
  ? static_cast<void>(0)						\
  : WDutils::AssertFailE(__STRING(expr),__FILE__,__LINE__))
  //
  template<bool Condition> struct AssertIf;
  template<> struct AssertIf<true>
  { static bool test(bool expr) { return expr; } };
  template<> struct AssertIf<false>
  { static bool test(bool) { return true; } };
  /// use instead of @a if(cond) WDutilsAssert(expr) with cond a constexpr
#  define WDutilsAssertIf(cond,expr)					\
    ((AssertIf<cond>(expr))						\
  ? static_cast<void>(0)						\
  : WDutils::AssertFail(__STRING(expr),__FILE__,__LINE__))
  /// use instead of @a if(cond) WDutilsAssertE(expr) with cond a constexpr
#  define WDutilsAssertEIf(cond,expr)					\
    ((AssertIf<cond>(expr))						\
  ? static_cast<void>(0)						\
  : WDutils::AssertFailE(__STRING(expr),__FILE__,__LINE__))
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
  int snprintf(char*str, size_t size, const char* fmt, ... );
  // ///////////////////////////////////////////////////////////////////////////
  struct snprintf__ {
    const char*    file;
    const unsigned line;
    snprintf__(const char*f, unsigned l) : file(f), line(l) {}
    int operator() (char*str, size_t size, const char* fmt, ... )
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
