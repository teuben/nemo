// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// exit.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002                                               |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef included_exit_h
#define included_exit_h 1

namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  void exit   (const int);                         // I: error signal           
  //----------------------------------------------------------------------------
  void error  (const char*,                        // I: error message format   
               ...         );                      //[I: parameters]            
  //----------------------------------------------------------------------------
  void warning(const char*,                        // I: warning message format 
               ...         );                      //[I: parameters]            
  //////////////////////////////////////////////////////////////////////////////
}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
// useful macros for simple error messages                                      
////////////////////////////////////////////////////////////////////////////////
#define MemoryCheck(POINTER)						       \
if((POINTER) == 0)							       \
  error("[%s.%d]: cannot allocate memory", __FILE__, __LINE__)
//------------------------------------------------------------------------------
#define NbdyError(MSGS)                                                 \
  error("[%s.%d]: %s",__FILE__,__LINE__,MSGS);
//------------------------------------------------------------------------------
#define NbdyWarning(MSGS)                                               \
  warning("[%s.%d]: %s",__FILE__,__LINE__,MSGS);
//------------------------------------------------------------------------------
#define NbdyErrorF(MSGS,FUNC)					\
  {								\
    if(FUNC)							\
      error("[%s.%d]: in %s: %s",__FILE__,__LINE__,FUNC,MSGS);	\
    else							\
      error("[%s.%d]: %s",__FILE__,__LINE__,MSGS);		\
  }
//------------------------------------------------------------------------------
#define NbdyWarningF(MSGS,FUNC)						\
  {									\
    if(FUNC)								\
      warning("[%s.%d]: in %s: %s",__FILE__,__LINE__,FUNC,MSGS);	\
    else								\
      warning("[%s.%d]: %s",__FILE__,__LINE__,MSGS);			\
  }
////////////////////////////////////////////////////////////////////////////////
#endif                                             // included_exit_h           
