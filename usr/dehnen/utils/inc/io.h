// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/io.h
///
/// \brief  contains declarations of classes WDutils::input, WDutils::output,
///         as well as WDutils::FortranIRec and WDutils::FortranORec
///
/// \author Walter Dehnen
/// \date   2000-2010
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2010 Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but          
// WITHOUT ANY WARRANTY; without even the implied warranty of                   
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
// General Public License for more details.                                     
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 675
// Mass Ave, Cambridge, MA 02139, USA.
//
////////////////////////////////////////////////////////////////////////////////
#ifndef WDutils_included_io_h
#define WDutils_included_io_h

#ifndef WDutils_included_iostream
#  include <iostream>
#  define Wdutils_included_iostream
#endif
#ifndef WDutils_included_fstream
#  include <fstream>
#  define Wdutils_included_fstream
#endif
#ifndef WDutils_included_cstdio
#  include <cstdio>
#  define WDutils_included_cstdio
#endif
#ifndef WDutils_included_cstring
#  include <cstring>
#  define WDutils_included_cstring
#endif
#ifndef WDutils_included_cmath
#  include <cmath>
#  define Wdutils_included_cmath
#endif
#ifndef WDutils_included_string
#  include <string>
#  define WDutils_included_string
#endif
#ifndef WDutils_included_exception_h
#  include <exception.h>
#endif
#ifndef WDutils_included_traits_h
#  include <traits.h>
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
  /// \param[in,out] S ofstream to open
  /// \param[in]     file name of file to open
  /// \param[in]     mode (optional) open mode, default: output & truncation
  inline bool open(std::ofstream& S, const char* file,
		   std::ios::openmode mode = std::ios::out | std::ios::trunc)
  {
    S.open(file,mode);
    if(! S.is_open() ) {
      WDutils_Warning("cannot open file \"%s\" for output",file);
      return false;
    }
    return true;
  }
  /// try to open an output file.
  /// if not successful, we issue a fatal error
  /// \param[in,out] S ofstream to open
  /// \param[in]     file name of file to open
  /// \param[in]     mode (optional) open mode, default: output & truncation
  inline
  void open_error(std::ofstream& S, const char* file,
		  std::ios::openmode mode = std::ios::out | std::ios::trunc)
  {
    S.open(file,mode);
    if(! S.is_open() ) WDutils_Error("cannot open file \"%s\" for output",file);
  }
  /// try to open an output file for appending.
  /// if not successful, we issue a warning and return 0\n
  /// if file did already exist and was opened for appending, return 1\n
  /// if did not exist and but a new one has been opened, return 2
  /// \return see above
  /// \param[in,out] S ofstream to open
  /// \param[in]     file name of file to open
  inline int open_to_append(std::ofstream& S, const char* file)
  {
    S.open(file,std::ios::out | std::ios::app);
    if(S.is_open() ) return 1;
    S.open(file,std::ios::out);
    if(S.is_open() ) return 2;
    WDutils_Warning("cannot open file \"%s\" for appending",file);
    return 0;
  }
  /// try to open an output file for appending.
  /// if not successful, we issue a fatal error
  /// if file did already exist and was opened for appending, return 1\n
  /// if did not exist and but a new one has been opened, return 2
  /// \return see above
  /// \param[in,out] S ofstream to open
  /// \param[in]     file name of file to open
  inline int open_to_append_error(std::ofstream& S, const char* file)
  {
    S.open(file,std::ios::out | std::ios::app);
    if(S.is_open() ) return 1;
    S.open(file,std::ios::out);
    if(S.is_open() ) return 2;
    WDutils_Error("cannot open file \"%s\" for appending",file);
    return 0;
  }
  /// try to open an input file.
  /// if not successful, we issue a warning and return false
  /// \return successfully opened?
  /// \param[in,out] S ifstream to open
  /// \param[in]     file name of file to open
  /// \param[in]     mode (optoinal) opening mode, default: input
  inline bool open(std::ifstream& S, const char* file,
		   std::ios::openmode mode=std::ios::in)
  {
    S.open(file,mode);
    if(! S.is_open() ) {
      WDutils_Warning("cannot open file \"%s\" for input",file);
      return 0;
    }
    return 1;
  }
  /// try to open an input file.
  /// if not successful, we issue a fatal error
  /// \param[in,out] S ifstream to open
  /// \param[in]     file name of file to open
  /// \param[in]     mode (optoinal) opening mode, default: input
  inline void open_error(std::ifstream& S, const char* file,
			 std::ios::openmode mode=std::ios::in)
  {
    S.open(file,mode);
    if(! S.is_open() ) WDutils_Error("cannot open file \"%s\" for input",file);
  }
  /// try to open a file for in-, output, or both.
  /// if not successful, we issue a warning an return false
  /// \return successfully opened?
  /// \param[in,out] S fstream to open
  /// \param[in]     file name of file to open
  /// \param[in]     mode opening mode
  inline bool open(std::fstream& S, const char* file,
		   std::ios::openmode mode)
  {
    S.open(file,mode);
    if(! S.is_open() ) {
      WDutils_Warning("cannot open file \"%s\"",file);
      return 0;
    }
    return 1;
  }
  /// try to open a file for in-, output, or both.
  /// if not successful, we issue a fatal warning
  /// \param[in,out] S fstream to open
  /// \param[in]     file name of file to open
  /// \param[in]     mode opening mode
  inline void open_error(std::fstream& S, const char* file,
			 std::ios::openmode mode)
  {
    S.open(file,mode);
    if(! S.is_open() ) WDutils_Error("cannot open file \"%s\"",file);
  }
  //@}
  //----------------------------------------------------------------------------
  /// swallow the rest of the current line.
  /// reads all characters up to and including the next '\n'.
  /// \return istream read (this allows to put this routine in a >> >> sequence)
  /// \param[in,out] from istream to read from
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
  /// \param[in] i integer
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
  /// \param[in,out] i istream to read from
  /// \param[in]     c character to eat
  inline std::istream& skip_char(std::istream& i, const char c) {
    char x;
    i>>x;
    if(x==c) return i;
    return i.putback(x);
  }
  //----------------------------------------------------------------------------
  /// eat the next entire line, if next character equals a specific one
  /// \return istream used
  /// \param[in,out] i istream to read from
  /// \param[in]     c character to match if the entire line is to be skipped
  inline std::istream& skip_line(std::istream& i, const char c) {
    char x;
    i>>x;
    if(x==c) { SwallowRestofLine(i); return i; }
    else       return i.putback(x);
  }
  //----------------------------------------------------------------------------
  /// eat the next entire line, if next character equals a specific one
  /// \return has line been skipped?
  /// \param[in,out] i istream to read from
  /// \param[in]     c character to match if the entire line is to be skipped
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
  //@}
  // ///////////////////////////////////////////////////////////////////////////
  /// \name I/O of arrays whose size is known at run time
  //@{
  /// write and array in the format "a1 a2 a3 ... aN".
  /// \return ostream used
  /// \param[in,out] s ostream to write to
  /// \param[in]     x pointer to first element
  /// \param[in]     N size of array
  template<typename X> inline
  std::ostream& write_array(std::ostream&s, const X* x, unsigned N)
  {
    s << x[0];
    for(register unsigned i=1; i!=N; ++i) s<<" "<<x[i];
    return s;
  }
  //----------------------------------------------------------------------------
  /// read an array from a space-separated format.
  /// \return istream used
  /// \param[in,out] s istream to read from
  /// \param[in]     x pointer to first element
  /// \param[in]     N size of array
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
    int p,w,s;
    int width(double l) { // given precision, what is minimum width
      int il = 1+int(l);
      int fw = l<0? 3+p-il : il>=p? il : p+1;
      int ew = p+5;
      il = fw<ew? fw:ew;
      return x<0? il+1 : il;
    }
    smanip_fp_width(X __x, int __w, int __p, int __s=0)
      : x(__x), p(__p), w(__w), s(__s)
    {
      if(s) {
	if(x>0 && w>1) --w;
	else           s=0;
      }
      if(x == 0) return;
      double l=std::log10(std::abs(x));
      w =std::max(w,width(l));       // minimum width to achieve
      for(++p; width(l)<=w; ++p) ;   // try for more precision
      --p;
    }
  };
  template<typename X>
  inline std::ostream& operator<<(std::ostream&o, smanip_fp_width<X> const&m) {
    if(m.s)
      o << ' ';
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
  inline smanip_fp_width<X> print(X x, int w, int p, int s=0) {
    return smanip_fp_width<X>(x,w,p,s);
  }
  //----------------------------------------------------------------------------
  template<typename X>
  struct smanip_fp_vec_width {
    const X*x;
    int     n,p,w;
    smanip_fp_vec_width(const X*__x, int __n, int __w, int __p)
      : x(__x), n(__n), p(__p), w(__w) {}
  };
  template<typename X>
  inline std::ostream& operator<<(std::ostream&o,
				  smanip_fp_vec_width<X> const&m) {
    if(m.n) {
      o << smanip_fp_width<X>(m.x[0],m.w,m.p);
      for(int i=1; i<m.n; ++i)
	o << ' ' << smanip_fp_width<X>(m.x[i],m.w,m.p);
    }
    return o;
  }
  /// manipulator: write an array of floating point numbers, each with minimum
  /// width but maximum precision.
  template<typename X>
  inline smanip_fp_vec_width<X> print(const X*x, int n, int w, int p) {
    return smanip_fp_vec_width<X>(x,n,w,p);
  }
#ifdef WDutils_included_tupel_h
  template<int N, typename X>
  inline smanip_fp_vec_width<X> print(tupel<N,X> const&x, int w, int p) {
    return smanip_fp_vec_width<X>(static_cast<const X*>(x),N,w,p);
  }
#endif
  // ///////////////////////////////////////////////////////////////////////////
  //
  // WDutils::FileSize()
  //
  /// computes the file size by opening the file for reading and seeking.
  /// taken from www.codeproject.com/file/filesize.asp
  /// \note the reported file size may be too small (according to above source) 
  /// \param  sFileName name of file 
  /// \return size of file sFileName in bytes
  size_t FileSize(const char*sFileName);
  // ///////////////////////////////////////////////////////////////////////////
  class FortranIRec;
  class FortranORec;
  // ///////////////////////////////////////////////////////////////////////////
  //
  // class WDutils::iofile
  //
  // ///////////////////////////////////////////////////////////////////////////
  class iofile {
  protected:
    static const int FNAME_MAX_SIZE = 256;
    char             FNAME[FNAME_MAX_SIZE];
    const char      *FILE;
    iofile() : FILE(0) {}
    //
    // some compilers (icpc 11.1 for instance) complain about no virtual
    // destructor. However, when providing one, we get the most bizarre of
    // run-time errors with gcc: a segmentation fault in using output with
    // std::cout. (essentially, output::OUT refers to a different std::cout
    // object than the executable for some reason currently beyond my
    // comprehension). It affects, amongst others things gyrfalcON.
    //
    //     virtual~iofile() {}
    //
    void setfile(const char*fname) {
      DebugInfo(12,"iofile::setfile(%s): FILE=%p\n",fname,(void*)(FILE));
      if(fname && fname[0]) {
	strncpy(FNAME,fname,FNAME_MAX_SIZE);
	FILE = FNAME;
      } else
	FILE = 0;
      DebugInfo(12,"iofile::setfile(%s): FILE=%p = %s\n",fname,
		(void*)(FILE),FILE);
    }
  public:
    /// give file name, if any
    const char*const&file() const { return FILE; }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //
  // class WDutils::output
  //
  /// wrapper around std::ostream with some additional features
  ///
  /// On opening a file, a filename "-" is interpreted as \c stdout, but will
  /// be opened only if no other output or nemo_out is writing to \c stdout.
  /// Similarly, a filename "." is interpreted as sink, i.e. nothing will ever
  /// be written. The operator << is defined (as template) to map \c
  /// std::ostream::operator<<(), for other operations on \c std::ostream, use
  /// the member method stream().
  ///
  // ///////////////////////////////////////////////////////////////////////////
  class output : public iofile {
    std::ostream *OUT;
    bool          APPENDING;
    FortranORec  *FREC;
    friend class FortranORec;
    //--------------------------------------------------------------------------
    output           (output const&); // not implemented
    output& operator=(output const&); // not implemented
    //--------------------------------------------------------------------------
    void __open (bool);
  public:
    /// \name const boolean information
    //@{
    /// ready for output?
    bool   is_open     () const { return OUT != 0; }
    /// appending output?
    bool   is_appending() const { return APPENDING; }
    /// writes to file?
    bool   is_file     () const { return OUT && OUT != &std::cout; }
    /// writes to stdout?
    bool   is_stdout   () const { return OUT == &std::cout; }
    /// ready for output?
    operator bool      () const { return OUT != 0; }
    //@}
    //--------------------------------------------------------------------------
    /// \name type conversion to ostream
    //@{
    /// conversion to std::ostream&
    operator std::ostream      & ()       { return *OUT; }
    /// conversion to std::ostream const&
    operator std::ostream const& () const { return *OUT; }
    /// conversion to std::ostream&
    std::ostream      & stream   ()       { return *OUT; }
    /// conversion to std::ostream const&
    std::ostream const& stream   () const { return *OUT; }
    //@}
    //--------------------------------------------------------------------------
    /// \name destruction and closing
    //@{
    /// close any open files; \c stdout is freed for other output
    void close();
    /// like close()
    ~output() { close(); }
    //@}
    //--------------------------------------------------------------------------
    /// \name construction and opening
    //@{
    /// construction from nothing: nothing is opened
    output() : OUT(0) , APPENDING(false), FREC(0) {}
    /// construction from file name and potential option for appending.
    /// If \a fname equals 0 or ".", nothing is opened.  If \a fname equals "-",
    /// and no other output or nemo_out is opened to \c stdout, we map to \c
    /// stdout.  Otherwise, a file of name \a fname is created for output. An
    /// existing file of the same name is deleted unless \a append is true, in
    /// which case, we append to that existing file.
    explicit
    output(const char*fname, bool append=0) : APPENDING(false), FREC(0) {
      DebugInfo(4,"output(%s,%d)\n",fname,append);
      setfile(fname);
      __open(append);
    }
    /// construction from file name and potential option for appending.
    /// same as above, except for type of first argument
    explicit
    output(std::string const&fname, bool append=0) : APPENDING(false), FREC(0) {
      setfile(fname.c_str());
      __open(append);
    }
    /// close possible old stream, then proceed as in construction
    void open(const char*fname, bool append = 0) {
      close();
      setfile(fname);
      __open(append);
    }
    /// close possible old stream, then proceed as in construction
    void open(std::string const&fname, bool append = 0) {
      close();
      setfile(fname.c_str());
      __open(append);
    }
    /// open file with name made from \a format string and \a tag.
    ///
    /// A new file name is created from the C-style \a format string and the
    /// data \a tag provided via code equivalent to
    /// \code
    /// sprintf(filename, format, tag);
    /// \endcode
    /// If this file name differs from the current, the old file is closed and
    /// the new one opened. The idea is to provide the possibility of numbered
    /// output files as in the following code snippet
    /// \code
    /// output out;
    /// for(int i=0; i!=20; ++i) {
    ///   out.reopen("file%02d.dat",i);
    ///   out << i << std::endl;
    /// }
    /// \endcode
    /// creating the files \c file00.dat, \c file01.dat, ... \c file19.dat
    /// \return     has a new file been opened (and the old closed)?
    /// \param[in]  format C-style format string for file to open
    /// \param[in]  tag    datum needed in generating file to open
    /// \param[in]  append (optional) append (or overwrite) existing file?
    template<typename T> 
    bool reopen(const char*format, T const&tag, bool append=0)
    {
      char FNEW[FNAME_MAX_SIZE];
      SNprintf(FNEW,FNAME_MAX_SIZE,format,tag);
      if(OUT==0 || strcmp(FNEW,FNAME)) {
	open(FNEW,append);
	return true;
      } else
	return false;
    }
    //@}
    /// does file with name generated from \a format plus data \a tag exist?
    /// \param[in]  format C-style format string for file name
    /// \param[in]  tag    datum needed in generating file to open
    /// \return     true if file with name exists
    /// \note For non-unix systems, we try to construct an std::ifstream, which
    ///       essentially tests for existence and read permission.
    template<typename T>
    static bool file_exists(const char*format, T const&tag)
    {
      char FNEW[FNAME_MAX_SIZE];
      SNprintf(FNEW,FNAME_MAX_SIZE,format,tag);
#ifdef unix
      // simply use the access system call. F_OK asks for existence only.
      return 0 == access(FNEW, F_OK);
#else
      // in fact, this test is different, as existence and read permission are
      // required to open an ifstream.
      std::ifstream IN(FNEW);
      return IN.is_open();
#endif
    }
    //--------------------------------------------------------------------------
    /// flush() output
    void flush() { if(OUT) OUT->flush(); }
    //--------------------------------------------------------------------------
    /// \name formated output using operator <<
    //@{
    /// output of a single datum \e x
    template<typename X>
    output& operator<< (X const&x) {
      if(OUT) (*OUT) << x;
      return*this;
    }
    /// output of a character string
    output& operator<< (const char*str) {
      if(OUT) (*OUT) << str;
      return*this;
    }
    /// output of a manipulator \e m
    output& operator<< (std::ostream& (*p)(std::ostream&)) {
      if(OUT) (*OUT) << p;
      return*this;
    }
    /// output of a manipulator \e m
    output& operator<< (std::ios_base& (*p)(std::ios_base&)) {
      if(OUT) (*OUT) << p;
      return*this;
    }
    //@}
    //--------------------------------------------------------------------------
    /// unformatted output
    void write(const char*a, size_t n) {
      if(OUT) OUT->write(a,n);
    }
    //--------------------------------------------------------------------------
    /// call if opening any output to stdout from NEMO main
    /// will not allow more than one open output to stdout.
    static void open_std () WDutils_THROWING;
    /// call if closing any output to stdout from NEMO main
    static void close_std();
  };
  // ///////////////////////////////////////////////////////////////////////////
  //
  // class WDutils::input
  //
  /// wrapper around std::istream with some additional features
  ///
  /// On opening a file, a filename "-" is interpreted as \c stdin, but will
  /// be opened only if no other input or nemo_in is reading from \c stdin.
  /// The operator >> is defined (as template) to map
  /// std::istream::operator>>(), for other operations on std::istream, use
  /// the member method stream().
  ///
  // ///////////////////////////////////////////////////////////////////////////
  class input : public iofile {
    std::istream *IN;
    FortranIRec  *FREC;
    friend class FortranIRec;
    //--------------------------------------------------------------------------
    input           (input const&); // not implemented
    input& operator=(input const&); // not implemented
    //--------------------------------------------------------------------------
    void __open();
  public:
    /// \name const boolean information
    //@{
    /// ready for input?
    bool   is_open () const { return IN != 0; }
    /// reads from file?
    bool   is_file () const { return IN && IN != &std::cin; }
    /// reads from stdin?
    bool   is_stdin() const { return IN == &std::cin; }
    /// ready for input?
    operator bool  () const { return IN != 0; }
    //@}
    //--------------------------------------------------------------------------
    /// \name type conversion to istream
    //@{
    /// conversion to std::istream&
    operator std::istream      & ()       { return *IN; }
    /// conversion to std::istream const&
    operator std::istream const& () const { return *IN; }
    /// conversion to std::istream&
    std::istream      & stream   ()       { return *IN; }
    /// conversion to std::istream const&
    std::istream const& stream   () const { return *IN; }
    //@}
    //--------------------------------------------------------------------------
    /// \name destruction and closing
    //@{
    /// if \c stdin: clear \c stdin, otherwise delete ifstream
    void close();
    /// like close()
    ~input() { close(); }
    //@}
    //--------------------------------------------------------------------------
    /// \name construction and opening
    //@{
    /// construction from nothing: nothing is opened
    input() : IN(0), FREC(0) {}
    /// construction from file name.
    /// If \e fname equals "-", and no other input or nemo_in is opened to \c
    /// stdin, we map to \c stdin.  Otherwise, a file of name \e fname is
    /// opened for input.
    explicit
    input(const char*fname) : FREC(0) {
      setfile(fname);
      __open();
    }
    /// construction from file name.
    /// If \e fname equals "-", and no other input or nemo_in is opened to \c
    /// stdin, we map to \c stdin.  Otherwise, a file of name \e fname is
    /// opened for input.
    explicit
    input(std::string const&fname) : FREC(0) {
      setfile(fname.c_str());
      __open();
    }
    /// close possible old stream, then proceed as in construction
    void open(const char*fname) {
      close();
      setfile(fname);
      __open();
    }
    //@}
    /// close possible old stream, then proceed as in construction
    void open(std::string const&fname) {
      close();
      setfile(fname.c_str());
      __open();
    }
    //--------------------------------------------------------------------------
    /// \name formated input using operator >>
    //@{
    /// input of single datum \e x
    template<typename X>
    input& operator>> (X&x) {
      if(IN) (*IN) >> x;
      return*this;
    }
    /// input of manipulator \e p
    input& operator>> (std::istream& (*p)(std::istream&)) {
      if(IN) (*IN) >> p;
      return*this;
    }
    /// input of manipulator \e p
    input& operator>> (std::ios_base& (*p)(std::ios_base&)) {
      if(IN) (*IN) >> p;
      return*this;
    }
    //@}
    //--------------------------------------------------------------------------
    /// unformatted input
    void read(char*a, size_t n) {
      if(IN) IN->read(a,n);
    }
    //--------------------------------------------------------------------------
    /// call if opening any input to stdin from NEMO main
    /// will not allow more than one open input from stdin.
    static void open_std () WDutils_THROWING;
    /// call if closing any output to stdout from NEMO main
    static void close_std();
  };
  // ///////////////////////////////////////////////////////////////////////////
  //
  // class WDutils::FortranIRec
  //
  /// represents a record for unformatted FORTRAN style input;
  /// used for reading gadget data files
  ///
  // ///////////////////////////////////////////////////////////////////////////
  class FortranIRec {
  private:
    input           &IN;                  // related input stream
    const unsigned   HSZE;                // size of header: 4 or 8
    const bool       SWAP;                // swap bytes for headers
    unsigned         SIZE;                // size (bytes) of record
    mutable unsigned READ;                // number of bytes already read
    //--------------------------------------------------------------------------
    FortranIRec           (FortranIRec const&); // not implemented
    FortranIRec& operator=(FortranIRec const&); // not implemented
    unsigned read_size() throw(WDutils::exception);
  public:
    //--------------------------------------------------------------------------
    /// constructor: read buffer with size information
    /// \param[in,out] in  WDutils::input to read from
    /// \param[in]     rec (optional) size of Fortran record header: 4 or 8
    /// \param[in]     bswap (optional) swap bytes for size information?
    FortranIRec(input&in, unsigned rec=4, bool bswap=0)
      throw(WDutils::exception);
    //--------------------------------------------------------------------------
    /// close: same as destruction
    void close() throw(WDutils::exception);
    //--------------------------------------------------------------------------
    /// destructor: read to end of record, read end buffer
    ~FortranIRec() throw(WDutils::exception) { close(); }
    //--------------------------------------------------------------------------
    /// read some bytes
    ///
    /// If more bytes are wanted than left in the record, only those left
    /// in the record will be read and a warning be issued.
    ///
    /// \return     number of bytes actually read
    /// \param[out] buf buffer to read into
    /// \param[in]  n   number of bytes to read
    unsigned read_bytes(char*buf, unsigned n) throw(WDutils::exception);
    //--------------------------------------------------------------------------
    /// read some data of any type
    ///
    /// If more data are wanted than left in the record, only those left
    /// in the record will be read and a warning be issued.
    ///
    /// \return     number of data actually read
    /// \param[out] buf buffer to read into
    /// \param[in]  n   number of data to read
    template<typename T>
    unsigned read(T*buf, unsigned n) throw(WDutils::exception) {
      if(READ+n*sizeof(T) > SIZE) {
	WDutils_Warning("FortranIRec::read(): cannot read %d, but only %d %s\n",
			n, (SIZE-READ)/sizeof(T), nameof(T));
	n = (SIZE-READ)/sizeof(T);
      }
      if(n) read_bytes(static_cast<char*>
		      (static_cast<void*>(buf)), sizeof(T)*n);
      return n;
    }
    //--------------------------------------------------------------------------
    /// skip some bytes
    ///
    /// \param n   number of bytes to skip
    void skip_bytes(unsigned n);
    //--------------------------------------------------------------------------
    /// read a single FORTRAN record in one go (you have to know its size!)
    ///
    /// \param[in]  in  input stream to read from
    /// \param[out] buf data buffer to read into
    /// \param[in]  n   number of data of type T to read
    /// \param[in]  rec size of FORTRAN record header; must be 4 or 8
    template<typename T>
    static void Read(input &in, T*buf, unsigned n, unsigned rec=4)
      throw(WDutils::exception)
    {
      FortranIRec FIR(in,rec);
      if( sizeof(T) * n > FIR.size() )
	throw exception("ReadFortranRecord(): cannot read %d %s: "
			"only %d bytes in record (required are %d)\n",
			n,nameof(T),FIR.size(),sizeof(T)*n);
      if( sizeof(T) * n < FIR.size() )
	WDutils_Warning("ReadFortranRecord(): "
			"reading %d %s: only %d of %d in record\n",
		       n,nameof(T),sizeof(T)*n,FIR.size());
      FIR.read(buf,n);
    }
    //--------------------------------------------------------------------------
    /// information on number of bytes already read
    unsigned const&bytes_read() const { return READ; }
    //--------------------------------------------------------------------------
    /// information on number of bytes yet to be read
    unsigned bytes_unread() const { return SIZE-READ; }
    //--------------------------------------------------------------------------
    /// information on total size of record
    unsigned const&size() const { return SIZE; }
  };// class FortranIRec
  // ///////////////////////////////////////////////////////////////////////////
  //
  // class WDutils::FortranORec
  //
  /// represents a record for unformatted FORTRAN style output;
  /// used for writing gadget data files
  ///
  // ///////////////////////////////////////////////////////////////////////////
  class FortranORec {
  private:
    output          &OUT;                 // related output stream
    const unsigned   HSZE;                // size of header: 4 or 8
    unsigned         SIZE;                // size (bytes) of record
    mutable unsigned WRITTEN;             // number of bytes already written
    //--------------------------------------------------------------------------
    FortranORec           (FortranORec const&); // not implemented
    FortranORec& operator=(FortranORec const&); // not implemented
    void write_size() throw(WDutils::exception);
  public:
    //--------------------------------------------------------------------------
    /// constructor: write buffer with size information
    /// \param[in,out] out  output stream to write to
    /// \param[in]     size size (in bytes) of record
    /// \param[in]     rec  (optional) size of header must be 4 or 8
    FortranORec(output&out, unsigned size, unsigned rec=4)
      throw(WDutils::exception);
    //--------------------------------------------------------------------------
    /// destructor: write to end of record, write end buffer
    ~FortranORec() throw(WDutils::exception) { close(); }
    //--------------------------------------------------------------------------
    /// close: same as destruction
    void close() throw(WDutils::exception);
    //--------------------------------------------------------------------------
    /// write some bytes
    ///
    /// If more bytes are to be written than left for the record, we only
    /// write as many as the record allows but issue a warning.
    ///
    /// \return number of bytes actually written
    /// \param buf buffer to write
    /// \param n   number of bytes to write
    unsigned write_bytes(const char*buf, unsigned n) throw(WDutils::exception);
    //--------------------------------------------------------------------------
    /// fill some bytes with a given value
    ///
    /// \param[in]   n   number of bytes to fill
    /// \param[in]   val value to fill them with
    void fill_bytes(unsigned n, char val=0);
    //--------------------------------------------------------------------------
    /// write some data of any type
    ///
    /// If more data are to be written than left for the record, we only
    /// write as many as the record allows but issue a warning.
    ///
    /// \return    number of data actually written
    /// \param[in] buf buffer to write
    /// \param[in] n   number of data to write
    template<typename T>
    unsigned write(const T*buf, unsigned n) throw(WDutils::exception) {
      if(WRITTEN + n*sizeof(T) > SIZE) {
	WDutils_Warning("FortranORec::write(): "
			"cannot write %u, but only %lu  %s\n",
			n, (unsigned long)((SIZE-WRITTEN)/sizeof(T)),
			nameof(T));
	n = (SIZE-WRITTEN)/sizeof(T);
      }
      if(n) write_bytes(static_cast<const char*>
		       (static_cast<const void*>(buf)), sizeof(T)*n);
      return n;
    }
    //--------------------------------------------------------------------------
    /// write a single FORTRAN record in one go
    ///
    /// \param[in,out] out output stream to write to
    /// \param[in]     buf data buffer to write from
    /// \param[in]     n   number of data of type T to write
    /// \param[in]     rec size of FORTRAN record header; must be 4 or 8
    template<typename T>
    static void Write(output&out, const T*buf, unsigned n, unsigned rec=4)
      throw(WDutils::exception)
    {
      FortranORec FOR(out, sizeof(T)*n, rec);
      FOR.write(buf,n);
    }
    //--------------------------------------------------------------------------
    /// information on number of bytes already written
    unsigned const&bytes_written() const { return WRITTEN; }
    //--------------------------------------------------------------------------
    /// information on number of bytes yet to be written
    unsigned bytes_free() const { return SIZE-WRITTEN; }
    //--------------------------------------------------------------------------
    /// information on total size of record
    unsigned const&size() const { return SIZE; }
  };// class FortranORec
  // ///////////////////////////////////////////////////////////////////////////
  //
  // support for byte-swapping
  //
  // ///////////////////////////////////////////////////////////////////////////
  namespace {
    inline void swap_char(char&A, char&B) { char T(A); A=B; B=T; }
    template<int B> struct __bswap;
    template<> struct __bswap<1> {
      static void swap(void*, unsigned) {}
    };
    template<> struct __bswap<2> {
      static void swap(void*vdat, unsigned cnt) {
	char*dat = static_cast<char*>(vdat);
	for(; cnt; --cnt,dat+=2) {
	  swap_char(dat[0],dat[1]);
	}
      }
    };
    template<> struct __bswap<4> {
      static void swap(void*vdat, unsigned cnt) {
	char*dat = static_cast<char*>(vdat);
	for(; cnt; --cnt, dat+=4) {
	  swap_char(dat[0],dat[3]);
	  swap_char(dat[1],dat[2]);
	}
      }
    };
    template<> struct __bswap<8> {
      static void swap(void*vdat, unsigned cnt) {
	char*dat = static_cast<char*>(vdat);
	for(; cnt; --cnt,dat+=8) {
	  swap_char(dat[0],dat[7]);
	  swap_char(dat[1],dat[6]);
	  swap_char(dat[2],dat[5]);
	  swap_char(dat[3],dat[4]);
	}
      }
    };
    template<> struct __bswap<16> {
      static void swap(void*vdat, unsigned cnt) {
	char*dat = static_cast<char*>(vdat);
	for(; cnt; --cnt,dat+=16) {
	  swap_char(dat[0],dat[15]);
	  swap_char(dat[1],dat[14]);
	  swap_char(dat[2],dat[13]);
	  swap_char(dat[3],dat[12]);
	  swap_char(dat[4],dat[11]);
	  swap_char(dat[5],dat[10]);
	  swap_char(dat[6],dat[ 9]);
	  swap_char(dat[7],dat[ 8]);
	}
      }
    };
  }
  // ///////////////////////////////////////////////////////////////////////////
  /// swap the bytes for one object of any type
  ///
  /// in order for this to work, the sizeof() the type must be 1,2,4,8, or 16
  ///
  /// \param[in,out] bdat element to swap bytes for
  template<typename B> inline
  void swap_bytes(B&bdat) {
    __bswap<sizeof(B)>::swap(static_cast<void*>(&bdat), 1);
  }
  /// swap the bytes of elements of any type
  ///
  /// in order for this to work, the sizeof() the type must be 1,2,4,8, or 16
  ///
  /// \param[in,out] bdat first element to swap bytes for
  /// \param[in]     cnt  number of elments to swap bytes for
  template<typename B> inline
  void swap_bytes(B*bdat, unsigned cnt) {
    __bswap<sizeof(B)>::swap(static_cast<void*>(bdat), cnt);
  }
  /// swap the bytes of elements of unknown type but known size
  ///
  /// \param[in,out] vdat pointer to first element
  /// \param[in]     len  size of the elements, must be 1,2,4,8, or 16
  /// \param[in]     cnt  number of elments to swap bytes for
  inline
  void swap_bytes(void*vdat, unsigned len, unsigned cnt) WDutils_THROWING {
    switch(len) {
    case  1: return __bswap< 1>::swap(vdat,cnt);  // does nothing
    case  2: return __bswap< 2>::swap(vdat,cnt);
    case  4: return __bswap< 4>::swap(vdat,cnt);
    case  8: return __bswap< 8>::swap(vdat,cnt);
    case 16: return __bswap<16>::swap(vdat,cnt);
    default: WDutils_THROW("swap_bytes(): sizeof(type)=%d: not supported\n",
			   len);
    }
  }
  // ///////////////////////////////////////////////////////////////////////////
  WDutils_TRAITS(WDutils::output,"output");
  WDutils_TRAITS(WDutils::input,"input");
  WDutils_TRAITS(WDutils::FortranIRec,"FortranIRec");
  WDutils_TRAITS(WDutils::FortranORec,"FortranORec");
} // namespace WDutils {
// /////////////////////////////////////////////////////////////////////////////
#ifndef WDutils_included_iomanip
# define WDutils_included_iomanip
# include <iomanip>
#endif
namespace WDutils {
  using std::setw;
  using std::setprecision;
  using std::setiosflags;
}
////////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_io_h
