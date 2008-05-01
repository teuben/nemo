// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    utils/inc/inline_io.h                                              
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    2000-2008                                                          
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2000-2008  Walter Dehnen                                       
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
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// \name inline functions for opening file I/O                               
  //@{                                                                          
  // ///////////////////////////////////////////////////////////////////////////
  /// try to open an output file.
  /// if not successful, we issue a warning and return false
  /// \return successfully opened?
  /// \param S ofstream to open
  /// \param file name of file to open
  /// \param mode (optional) open mode, default: output & truncation
  inline bool open(std::ofstream& S, const char* file,
		   std::ios::openmode mode = std::ios::out | std::ios::trunc)
  {
    S.open(file,mode);
    if(! S.is_open() ) {
      warning("cannot open file \"%s\" for output",file);
      return false;
    }
    return true;
  }
  /// try to open an output file.
  /// if not successful, we issue a fatal error
  /// \param S ofstream to open
  /// \param file name of file to open
  /// \param mode (optional) open mode, default: output & truncation
  inline
  void open_error(std::ofstream& S, const char* file,
		  std::ios::openmode mode = std::ios::out | std::ios::trunc)
  {
    S.open(file,mode);
    if(! S.is_open() ) error("cannot open file \"%s\" for output",file);
  }
  /// try to open an output file for appending.
  /// if not successful, we issue a warning and return 0\n
  /// if file did already exist and was opened for appending, return 1\n
  /// if did not exist and but a new one has been opened, return 2
  /// \return see above
  /// \param S ofstream to open
  /// \param file name of file to open
  inline int open_to_append(std::ofstream& S, const char* file)
  {
    S.open(file,std::ios::out | std::ios::app);
    if(S.is_open() ) return 1;
    S.open(file,std::ios::out);
    if(S.is_open() ) return 2;
    warning("cannot open file \"%s\" for appending",file);
    return 0;
  }
  /// try to open an output file for appending.
  /// if not successful, we issue a fatal error
  /// if file did already exist and was opened for appending, return 1\n
  /// if did not exist and but a new one has been opened, return 2
  /// \return see above
  /// \param S ofstream to open
  /// \param file name of file to open
  inline int open_to_append_error(std::ofstream& S, const char* file)
  {
    S.open(file,std::ios::out | std::ios::app);
    if(S.is_open() ) return 1;
    S.open(file,std::ios::out);
    if(S.is_open() ) return 2;
    error("cannot open file \"%s\" for appending",file);
    return 0;
  }
  /// try to open an input file.
  /// if not successful, we issue a warning and return false
  /// \return successfully opened?
  /// \param S ifstream to open
  /// \param file name of file to open
  /// \param mode (optoinal) opening mode, default: input
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
  /// try to open an input file.
  /// if not successful, we issue a fatal error
  /// \param S ifstream to open
  /// \param file name of file to open
  /// \param mode (optoinal) opening mode, default: input
  inline void open_error(std::ifstream& S, const char* file,
			 std::ios::openmode mode=std::ios::in)
  {
    S.open(file,mode);
    if(! S.is_open() ) error("cannot open file \"%s\" for input",file);
  }
  /// try to open a file for in-, output, or both.
  /// if not successful, we issue a warning an return false
  /// \return successfully opened?
  /// \param S fstream to open
  /// \param file name of file to open
  /// \param mode opening mode
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
  /// try to open a file for in-, output, or both.
  /// if not successful, we issue a fatal warning
  /// \param S fstream to open
  /// \param file name of file to open
  /// \param mode opening mode
  inline void open_error(std::fstream& S, const char* file,
			 std::ios::openmode mode)
  {
    S.open(file,mode);
    if(! S.is_open() ) error("cannot open file \"%s\"",file);
  }
  //@}
  //----------------------------------------------------------------------------
  /// swallow the rest of the current line.
  /// reads all character up to and including the next '\n'.
  /// \return istream read (this allows to put this routine in a >> >> sequence)
  /// \param from istream to read from
  inline std::istream&SwallowRestofLine(std::istream& from)
  {
    char c;
    do from.get(c); while( from.good() && c !='\n');
    return from;
  }
  //----------------------------------------------------------------------------
  /// return shorthand for making ordered number.
  /// given an integer, return "st", "nd", "rd" or "th" such that appended to
  /// the integer it makes correct ordered number, e.g. "3rd", "8th", "101st".
  /// \return "st", "nd", "rd" or "th"
  /// \param i integer
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
  /// for a>=0 return " ", and for a<0 nothing (zero pointer), such that output
  /// like cout << neg_space(a) << a looks aligned whatever sign a has.
  template<typename scalar>
  inline const char* neg_space(scalar const&a)
  {
    return a < scalar(0)? 0 : " ";
  }
  //----------------------------------------------------------------------------
  /// eat the next char in an istream, if it equals a specific one
  /// \return istream used
  /// \param i istream to read from
  /// \param c character to eat
  inline std::istream& skip_char(std::istream& i, const char c) {
    char x;
    i>>x;
    if(x==c) return i;
    return i.putback(x);
  }
  //----------------------------------------------------------------------------
  /// eat the next entire line, if next character equals a specific one
  /// \return istream used
  /// \param i istream to read from
  /// \param c character to match if the entire line is to be skipped
  inline std::istream& skip_line(std::istream& i, const char c) {
    char x;
    i>>x;
    if(x==c) { SwallowRestofLine(i); return i; }
    else       return i.putback(x);
  }
  //----------------------------------------------------------------------------
  /// eat the next entire line, if next character equals a specific one
  /// \return has line been skipped?
  /// \param i istream to read from
  /// \param c character to match if the entire line is to be skipped
  inline bool eat_line(std::istream& i, const char c) {
    char x;
    i>>x;
    if(x==c) { SwallowRestofLine(i); return true; }
    else     { i.putback(x);         return false; }
  }
  // ///////////////////////////////////////////////////////////////////////////
  /// \name I/O of arrays whose size is known at compile time
  //@{
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
  /// write and array in the format "a1 a2 a3 ... aN".
  /// \param N template parameter: size of array
  /// \param X template parameter: type of array elements
  /// \return ostream used
  /// \param s ostream to write to
  /// \param x pointer to first element
  template<int N, typename X> inline
  std::ostream& write_arr(std::ostream&s, const X* x)
  { 
    meta_io<N-1,0>::write(s,x);
    return s;
  }
  //----------------------------------------------------------------------------
  /// read an array from a space-separated format.
  /// \param N template parameter: size of array
  /// \param X template parameter: type of array elements
  /// \return istream used
  /// \param s istream to read from
  /// \param x pointer to first element
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
  //@}
  // ///////////////////////////////////////////////////////////////////////////
  /// \name I/O of arrays whose size is known at run time
  //@{
  /// write and array in the format "a1 a2 a3 ... aN".
  /// \param X template parameter: type of array elements
  /// \return ostream used
  /// \param s ostream to write to
  /// \param x pointer to first element
  /// \param N size of array
  template<typename X> inline
  std::ostream& write_array(std::ostream&s, const X* x, unsigned N)
  {
    s << x[0];
    for(register unsigned i=1; i!=N; ++i) s<<" "<<x[i];
    return s;
  }
  //----------------------------------------------------------------------------
  /// read an array from a space-separated format.
  /// \param X template parameter: type of array elements
  /// \return istream used
  /// \param s istream to read from
  /// \param x pointer to first element
  /// \param N size of array
  template<typename X> inline
  std::istream& read_array(std::istream&s, X* x, unsigned N) {
    register X y[N];
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
  //}@
  //----------------------------------------------------------------------------
  template<typename X>
  struct smanip_fp_width {
    X   x;
    int p,w;
    int width(double l) { // given precision, what is minimum width
      int il = 1+int(l);
      int fw = l<0? 3+p-il : il>=p? il : p+1;
      int ew = p+5;
      il = fw<ew? fw:ew;
      return x<0? il+1 : il;
    }
    smanip_fp_width(X __x, int __w, int __p) : x(__x), p(__p), w(__w)
    {
      if(x == 0) return;
      double l=std::log10(std::abs(x));
      w =std::max(w,width(l));       // minimum width to achieve
      for(++p; width(l)<=w; ++p);    // try for more precision
      --p;
    }
  };
  template<typename X>
  inline std::ostream& operator<<(std::ostream&o, smanip_fp_width<X> const&m) {
    int ow = o.width(m.w);
    int op = o.precision(m.p);
    o << m.x;
    o.width(ow);
    o.precision(op);
    return o;
  }
  /// manipulator: write a floating point number with minimum width but maximum
  /// precision.
  /// The floating point number \a x is written out with the maximum precision
  /// possible in the \a w characters wide field. However, we will at least 
  /// write it with precision \a p, even if this means overrunning the width.
  template<typename X>
  inline smanip_fp_width<X> print(X x, int w, int p) {
    return smanip_fp_width<X>(x,w,p);
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
#endif // WDutils_included_inline_io_h
