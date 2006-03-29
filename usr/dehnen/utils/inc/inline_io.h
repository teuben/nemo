// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    inc/inline_io.h                                                    
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    2000-2005                                                          
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2000-2005  Walter Dehnen                                       
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
#ifndef WDutils_included_inline_io_h
#define WDutils_included_inline_io_h

#ifndef WDutils_included_basic_h
#  include <exception.h>
#endif
#ifndef WDutils_included_iostream
#  include <iostream>
#  define WDutils_included_iostream
#endif
#ifndef WDutils_included_fstream
#  include <fstream>
#  define WDutils_included_fstream
#endif
#ifndef WDutils_included_cstdio
#  include <cstdio>
#  define WDutils_included_cstdio
#endif
//------------------------------------------------------------------------------
namespace WDutils {
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
  inline std::istream& eatl(std::istream& from)
  {
    SwallowRestofLine(from);
    return from;
  }
  //----------------------------------------------------------------------------
  // return "st", "nd", "rd", "th" given i to make a ordered number             
  //----------------------------------------------------------------------------
  inline const char* stndrdth(int i)
  {
    if(i<0) i=-i;
    switch( i % 100 ) {
    case 11: 
    case 12: 
    case 13: return "th";
    default:
      switch( i % 10 ) {
      case 1:  return "st";
      case 2:  return "nd"; 
      case 3:  return "rd";
      default: return "th";
      }   
    }
  }
  //----------------------------------------------------------------------------
  // return " " or nothing (zero pointer) for a being positive or negative      
  //----------------------------------------------------------------------------
  template<typename scalar>
  inline const char* neg_space(scalar const&a)
  {
    return a < scalar(0)? "" : " ";
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
  //----------------------------------------------------------------------------
} // namespace WDutils

#ifndef WDutils_included_iomanip
# define WDutils_included_iomanip
# include <iomanip>
#endif
namespace WDutils {
  using std::setw;
  using std::setprecision;
  using std::setiosflags;
}
//------------------------------------------------------------------------------
#endif // WDutils_included_ioln_h
