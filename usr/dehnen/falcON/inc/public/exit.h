// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// exit.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_exit_h
#define falcON_included_exit_h 1

#ifndef falcON_included_cstdlib
#  include <cstdlib>
#  define falcON_included_cstdlib
#endif

namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  void set_name(const char*);
  //----------------------------------------------------------------------------
  void exit    (const int);                        // I: error signal           
  //----------------------------------------------------------------------------
  void error   (const char*,                       // I: error message format   
                ...         );                     //[I: parameters]            
  //----------------------------------------------------------------------------
  void warning (const char*,                       // I: warning message format 
                ...         );                     //[I: parameters]            
  //----------------------------------------------------------------------------
  template<typename TYPE>
  TYPE* PointerCheck(TYPE*pter, const char*file, int const&line)
  {
    if (pter == 0 ) error("[%s.%d]: cannot allocate memory",file,line);
    return pter;
  }
  //////////////////////////////////////////////////////////////////////////////
  // useful macros for simple error messages                                    
  //////////////////////////////////////////////////////////////////////////////
#define falcON_New(TYPE,SIZE)					\
  PointerCheck(new TYPE[SIZE],__FILE__,__LINE__)
  //----------------------------------------------------------------------------
#define falcON_New16(TYPE,SIZE)					\
  PointerCheck(new16<TYPE>(SIZE),__FILE__,__LINE__)
  //----------------------------------------------------------------------------
#define falcON_Memory(POINTER)					\
  PointerCheck(POINTER,__FILE__,__LINE__)
  //----------------------------------------------------------------------------
#define falcON_Error(MSGS)					\
  error("[%s.%d]: %s",__FILE__,__LINE__,MSGS)
  //----------------------------------------------------------------------------
#define falcON_Warning(MSGS)					\
  warning("[%s.%d]: %s",__FILE__,__LINE__,MSGS)
  //----------------------------------------------------------------------------
#define falcON_ErrorF(MSGS,FUNC)				\
  error("[%s.%d]: in %s: %s",__FILE__,__LINE__,FUNC,MSGS)
  //----------------------------------------------------------------------------
#define falcON_WarningF(MSGS,FUNC)				\
  warning("[%s.%d]: in %s: %s",__FILE__,__LINE__,FUNC,MSGS)
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
#    define falcON_RepAction 1
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
  // support for code information                                             //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  namespace compile_info {                         // info for general usage    
    bool const&is_set();                           //   is info provided?       
    const char*version();                          //   version number          
    const char*origin();                           //   version origin          
    const char*time();                             //   time of compilation     
    const char*compiler();                         //   compiler used           
  }                                                // implemented in main.h     
  namespace run_info {                             // info for general usage    
    bool const&host_known();
    bool const&user_known();
    bool const&pid_known();
    const char*time();                             //   time of run             
    const char*host();                             //   host of run             
    const char*user();                             //   user of run             
    const char*pid();                              //   pid of run              
    void       init();                             //   set info up             
  }                                                // implemented in exit.cc    
}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
#endif                                             // falcON_included_exit_h    
