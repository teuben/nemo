// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// ionl.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 1994-2001                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef included_ioln_h
#define included_ioln_h

#ifndef included_exit_h
#  include <public/exit.h>
#endif
#ifndef included_iostream
#  include <iostream>
#  define included_iostream
#endif
#ifndef included_fstream
#  include <fstream>
#  define included_fstream
#endif
#ifndef included_cstdio
#  include <cstdio>
#  define included_cstdio
#endif
//------------------------------------------------------------------------------
namespace nbdy {
  //----------------------------------------------------------------------------
  // inline functions for opening file I/O and, if not successfull either       
  // - issue a warning and return:               open()                         
  // - issue an error:                           open_error()                   
  // both exist for opening ofstreams, ofstreams for appending, and ifstreams,  
  // as well as general fstreams.                                               
  //----------------------------------------------------------------------------
  inline bool open      (                           // R:   success?            
			 std::ofstream& S,          // I/O: stream to open      
			 const char* file,          // I:   file name           
			 std::ios::openmode mode = std::ios::out |
			                           std::ios::trunc)
  {
    S.open(file,mode);
    if(! S.is_open() ) {
      warning("cannot open file \"%s\" for output",file);
      return false;
    }
    return true;
  }
  //----------------------------------------------------------------------------
  inline void open_error(                           // no success -> abort      
			 std::ofstream& S,          // I/O: stream to open      
			 const char* file,          // I: file name             
			 std::ios::openmode mode = std::ios::out |
			                           std::ios::trunc)
  {
    S.open(file,mode);
    if(! S.is_open() ) error("cannot open file \"%s\" for output",file);
  }
  //////////////////////////////////////////////////////////////////////////////
  // open_to_append returns an integer with the following meaning:              
  // 0  couldn't open a file                                                    
  // 1  file did already exist and has been opened for appending                
  // 2  file could not be opened for appending, but opened as normal            
  inline int open_to_append(                        // 0/1/                     
			    std::ofstream& S,       // I/O: stream to open      
			    const char* file)       // I: file name             
  {
    S.open(file,std::ios::out | std::ios::app);
    if(S.is_open() ) return 1;
    S.open(file,std::ios::out);
    if(S.is_open() ) return 2;
    warning("cannot open file \"%s\" for appending",file);
    return 0;
  }
  //----------------------------------------------------------------------------
  inline int open_to_append_error(                  // no success -> abort      
				   std::ofstream& S,// I/O: stream to open      
				   const char* file)// I: file name             
  {
    S.open(file,std::ios::out | std::ios::app);
    if(S.is_open() ) return 1;
    S.open(file,std::ios::out);
    if(S.is_open() ) return 2;
    error("cannot open file \"%s\" for appending",file);
  }
  //////////////////////////////////////////////////////////////////////////////
  inline bool open(std::ifstream& S, const char* file,
		   std::ios::openmode mode=std::ios::in)
  {
    S.open(file,mode);
    if(! S.is_open() ) {
      warning("cannot open file \"%s\" for input",file);
      return 0;
    }
    return 1;
  }
  //----------------------------------------------------------------------------
  inline void open_error(std::ifstream& S, const char* file,
			 std::ios::openmode mode=std::ios::in)
  {
    S.open(file,mode);
    if(! S.is_open() ) error("cannot open file \"%s\" for input",file);
  }
  //////////////////////////////////////////////////////////////////////////////
  inline bool open(std::fstream& S, const char* file,
		   std::ios::openmode mode)
  {
    S.open(file,mode);
    if(! S.is_open() ) {
      warning("cannot open file \"%s\"",file);
      return 0;
    }
    return 1;
  }
  //----------------------------------------------------------------------------
  inline void open_error(std::fstream& S, const char* file,
			 std::ios::openmode mode)
  {
    S.open(file,mode);
    if(! S.is_open() ) error("cannot open file \"%s\"",file);
  }
  //----------------------------------------------------------------------------
  // read all characters until '\n' (inclusive)                                 
  //----------------------------------------------------------------------------
  inline void SwallowRestofLine(std::istream& from)
  {
    char c;
    do from.get(c); while( from.good() && c !='\n');
  }
  //----------------------------------------------------------------------------
  // return "st", "nd", "rd", "th" given i to make a ordered number             
  //----------------------------------------------------------------------------
  inline const char* stndrdth(const int i)
  {
    register int ia= (i<0)? -i:i;
    switch( ia % 100 ) {
    case 11: 
    case 12: 
    case 13: return "th";
    default:
      switch( ia % 10 ) {
      case 1:  return "st";
      case 2:  return "nd"; 
      case 3:  return "rd";
      default: return "th";
      }   
    }
  }
  //----------------------------------------------------------------------------
  // if the next char in istream equals a specific one, eat it                  
  //----------------------------------------------------------------------------
  inline std::istream& skip_char(std::istream& i, const char c) {
    char x;
    i>>x;
    if(x==c) return i;
    return i.putback(x);
  }
  //----------------------------------------------------------------------------
  // if the next char in istream equals a specific one, eat the rest of the line
  //----------------------------------------------------------------------------
  inline std::istream& skip_line(std::istream& i, const char c) {
    char x;
    i>>x;
    if(x==c) { SwallowRestofLine(i); return i; }
    else       return i.putback(x);
  }
  //----------------------------------------------------------------------------
  // like skip_line, but returns a boolean                                      
  //----------------------------------------------------------------------------
  inline bool eat_line(std::istream& i, const char c) {
    char x;
    i>>x;
    if(x==c) { SwallowRestofLine(i); return true; }
    else     { i.putback(x);         return false; }
  }
}
//-----------------------------------------------------------------------------+
// here are some alternatives for std::setw, std::setprecision, which on some  |
// occasions do not compile properly with gcc 2.95.3                           |
//---------------------------------------------------------------------------- +
#if defined (__GNUC__) && (__GNUC__ < 3)
// actually, with gcc version 3.2, these crash (segmentation fault on setw)     
namespace nbdy {
  template<typename TP> class Smanip {
    TP (std::ios::*_f)(TP);
    TP _a;
  public:
    inline Smanip(TP (std::ios::*f)(TP), TP a) : _f(f), _a(a) {}
    inline std::istream& operator()(std::istream&i) const
    { (i.*_f)(_a); return i; }
    inline std::ostream& operator()(std::ostream&o) const
    { (o.*_f)(_a); return o; }
  };
  template<typename TP> inline 
  std::istream& operator>>(std::istream&i, Smanip<TP> const &m)
  { return m(i); }
  template<typename TP> inline
  std::ostream& operator<<(std::ostream&o, Smanip<TP> const &m)
  { return m(o); }

  inline Smanip<int> setw(int a)
         { return Smanip<int>( &std::ios::width ,a); }
  inline Smanip<int> setprecision(int a)
         { return Smanip<int>( &std::ios::precision ,a); }
  inline Smanip<std::ios::fmtflags> setiosflags(std::ios::fmtflags a)
         { return Smanip<std::ios::fmtflags>( &std::ios::setf ,a); }
}
#else // use std manipulators
#ifndef included_iomanip
#  include <iomanip>
#  define included_iomanip
#endif
namespace nbdy {
  using std::setw;
  using std::setprecision;
  using std::setiosflags;
}
#endif
//------------------------------------------------------------------------------
#endif // included_ioln_h
