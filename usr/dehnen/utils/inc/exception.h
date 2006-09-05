// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    utils/inc/exception.h                                              
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    2000-2006                                                          
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2000-2006  Walter Dehnen                                       
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
    char __time    [20];
    char __host   [100];
    char __user   [100];
    char __pid     [10];
    char __name   [100];
    char __cmd   [1024];
    int  __debug;
    RunInfo();
    static RunInfo Info;
  public:
    /// reset the debugging level
    static void set_debug_level(int d) { Info.__debug = d; }
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
    /// return true if debug level >= given debug depth
    static bool debug(int depth) { return Info.__debug >= depth; }
  };
  /// is debugging level exceeded by debugging depth (argument)?                
  /// \param d debugging depth
  inline bool debug(int d) {
    return RunInfo::debug(d);
  }
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// abort with error signal (1st argument)
  void exit(int = 1);
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// Print an header followed by a message to stderr.                          
  //                                                                            
  /// We use a template for va_list, as we don't want to include the C header   
  /// stdarg.h. The only implementation is for VA_LIST = va_list.               
  /// This routine is used by error(), warning(), and debug_info().             
  /// \note The macros va_start() and va_end MUST be called outside of printerr.
  template<typename VA_LIST>
  void printerr(const char*header,
		const char*fmt,
		VA_LIST   &list,
		bool       name = true);

  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// Print an error message to stderr and abort.                               
  //                                                                            
  /// Uses a printf() style format string as first argument, further arguments  
  /// must match format, exactly as in printf, which will be called.            
  /// \param fmt gives the format in C printf() style                           
  void error(const char*fmt, ... );

  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// Print a warning message to stderr.                                        
  //                                                                            
  /// Uses a printf() style format string as first argument, further arguments  
  /// must match format, exactly as in printf, which will be called.            
  /// \param fmt gives the format in C printf() style                           
  void warning(const char*fmt, ... );

  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// Print debugging information to stderr.                                    
  //                                                                            
  /// Debugging information is printed if the debugging level exceeds the       
  /// debugging depth (1st argument). The debugging level is found from         
  /// falcON::run_info::debug_level() and initialized at start-up of main().    
  /// Uses a printf() style format string as first argument, further arguments  
  /// must match format, exactly as in printf, which will be called.            
  /// \param deb gives the debugging depth: info is printed only if it is       
  ///            less or equal to debugging level                               
  /// \param fmt gives the format in C printf() style                           
  void debug_info(int deb, const char*fmt, ... );

  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  //  WDutils::message                                                          
  //                                                                            
  /// C++ wrapper around a C string.                                            
  ///                                                                           
  /// Construction from C-type format string + data;                            
  /// Type conversion to const char*                                            
  /// Useful for generating a C-style string containing formatted data.         
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class message {
    message(message const&);                       // no copy constructor       
    static const int size = 1024;
    char __text[size];
  public:
    /// Generate a string from format + data.
    /// Uses a printf() style format string as first argument, further arguments
    /// must match format, exactly as in printf, which will be called.
    /// \param fmt gives the format in C printf() style
    explicit message(const char* fmt, ...);
    /// conversion to C-style string
    operator const char*() const { return __text; }
    /// return C-style string
    const char* text() const { return __text; }
  };

  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  //  WDutils::exception                                                        
  //                                                                            
  /// base class for exceptions in falcON; derived from std::string             
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  struct exception : protected std::string {
    /// copy constructor
    exception(exception const&e) : std::string(e) {}
    /// construction from C-style format string + data.
    /// Uses a printf() style format string as first argument, further arguments
    /// must match format, exactly as in printf, which will be called.
    /// \param fmt gives the format in C printf() style
    explicit exception(const char*fmt, ...);
    /// return C-style string
    const char*text() const { return c_str(); }
    /// return C-style string
    friend const char*text(exception const&e) { return e.text(); }
  };

  // ///////////////////////////////////////////////////////////////////////////
  ///                                                                           
  /// \name macros controling the usage of throw exception vs error             
  //@{                                                                          
  // ///////////////////////////////////////////////////////////////////////////
#if 1
  //----------------------------------------------------------------------------
  /// use instead of <tt> throw(WDutils::exception) </tt> after function
  /// declaration
#  define WDutils_THROWING 
  //----------------------------------------------------------------------------
  /// use "WDutils_THROW(fmt, data)" instead of "error(fmt, data)" or "throw
  /// WDutils::exception(fmt, data)"
#  define WDutils_THROW         WDutils::error
  //----------------------------------------------------------------------------
  /// use "WDutils_RETHROW(E)" to re-throw a caught exception "E"
#  define WDutils_RETHROW(E)    WDutils::error(text(E))
#else
#  define WDutils_THROWING      throw(WDutils::exception)
#  define WDutils_THROW         throw WDutils::exception
#  define WDutils_RETHROW(E)    throw E
#endif
  //@}
  // ///////////////////////////////////////////////////////////////////////////
  ///                                                                           
  /// \name useful macros for simple error messages                             
  //@{                                                                          
  // ///////////////////////////////////////////////////////////////////////////
  /// error with information about file and line
#define WDutils_Error(MSGS)					\
  WDutils::error("[%s.%d]: %s",__FILE__,__LINE__,MSGS)
  //----------------------------------------------------------------------------
  /// throw exception with information about file and line
#define WDutils_Except(MSGS)					\
  WDutils_THROW("[%s.%d]: %s",__FILE__,__LINE__,MSGS)
  //----------------------------------------------------------------------------
  /// warning with information about file and line
#define WDutils_Warning(MSGS)					\
  WDutils::warning("[%s.%d]: %s",__FILE__,__LINE__,MSGS)
  //----------------------------------------------------------------------------
  /// error with information about file, line, and function (to be given)
#define WDutils_ErrorF(MSGS,FUNC)				\
  WDutils::error("[%s.%d]: in %s: %s",__FILE__,__LINE__,FUNC,MSGS)
  //----------------------------------------------------------------------------
  /// throw exception with information about file, line, and func (to be given)
#define WDutils_ExceptF(MSGS,FUNC)					\
  WDutils_THROW("[%s.%d]: in %s: %s",__FILE__,__LINE__,FUNC,MSGS)
  //----------------------------------------------------------------------------
  /// warning with information about file, line, and function (to be given)
#define WDutils_WarningF(MSGS,FUNC)				\
  WDutils::warning("[%s.%d]: in %s: %s",__FILE__,__LINE__,FUNC,MSGS)
  //@}
  // ///////////////////////////////////////////////////////////////////////////
} // namespace WDutils
// /////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_exception_h
