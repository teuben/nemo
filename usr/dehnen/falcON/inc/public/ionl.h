// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// ionl.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 1994-2003                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_ioln_h
#define falcON_included_ioln_h

#ifndef falcON_included_exit_h
#  include <public/exit.h>
#endif
#ifndef falcON_included_iostream
#  include <iostream>
#  define falcON_included_iostream
#endif
#ifndef falcON_included_fstream
#  include <fstream>
#  define falcON_included_fstream
#endif
#ifndef falcON_included_cstdio
#  include <cstdio>
#  define falcON_included_cstdio
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
    return 0;
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
    register int ia = (i<0)? -i:i;
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
  //----------------------------------------------------------------------------
  // I/O of arrays whose size is known at compile time                          
  //----------------------------------------------------------------------------
  template<int N, int I=0> struct meta_io {
  template<typename X> static void write(std::ostream&s, const X*x)
    { s<<x[I]<<' '; meta_io<N,I+1>::write(s,x); }
  template<typename X> static void read (std::istream&s,       X*x)
    { s>>x[I];      meta_io<N,I+1>::read (s,x); }
  };
  //----------------------------------------------------------------------------
  template<int N> struct meta_io<N,N> {
  template<typename X> static void write(std::ostream&s, const X*x) { s<<x[N]; }
  template<typename X> static void read (std::istream&s,       X*x) { s>>x[N]; }
  };
  //----------------------------------------------------------------------------
  template<int N, typename X> inline
  std::ostream& write_arr(std::ostream&s, const X* x)
  { 
    meta_io<N-1,0>::write(s,x);
    return s;
  }
  //----------------------------------------------------------------------------
  template<int N, typename X> inline
  std::istream& read_arr(std::istream&s, X* x)
  { 
    char c=0;
    s >> c;
    if(c == '(') {
      meta_io<N-1,0>::read(s,x);
      s >> c;
      if(c != ')') s.clear(std::ios::badbit);
    } else {
      s.putback(c);
      meta_io<N-1,0>::read(s,x);
    }
    return s;
  }
  //----------------------------------------------------------------------------
  // I/O of arrays whose size is known at run time                              
  //----------------------------------------------------------------------------
  template<typename S> inline
  std::ostream& write_array(std::ostream&s, const S* x, unsigned N)
  {
    s << x[0];
    for(register unsigned i=1; i!=N; ++i) s<<" "<<x[i];
    return s;
  }
  //----------------------------------------------------------------------------
  template<typename S> inline
  std::istream& read_array(std::istream&s, S* x, unsigned N) {
    register S y[N];
    char c=0;
    s >> c;
    if(c == '(') {
      for(register unsigned i=0; i!=N; ++i) s >> y[i];
      s >> c;
      if(c != ')') s.clear(std::ios::badbit);
    } else {
      s.putback(c);
      for(register unsigned i=0; i!=N; ++i) s >> y[i];
    }
    for(register unsigned i=0; i!=N; ++i) x[i] = y[i];
    return s;
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
#ifndef falcON_included_iomanip
#  include <iomanip>
#  define falcON_included_iomanip
#endif
namespace nbdy {
  using std::setw;
  using std::setprecision;
  using std::setiosflags;
}
#endif
//------------------------------------------------------------------------------
#endif // falcON_included_ioln_h
