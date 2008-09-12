// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   inc/public/io.h                                                     
///                                                                             
/// \brief  contains declarations of classes falcON::input and falcON::output,  
///         as well as NEMO I/O support                                         
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2000-2008                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2000-2008 Walter Dehnen                                        
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
#ifndef falcON_included_io_h
#define falcON_included_io_h

#ifndef falcON_included_iostream
#  include <iostream>
#  define falcON_included_iostream
#endif
#ifndef falcON_included_cstdio
#  include <cstdio>
#  define falcON_included_cstdio
#endif
#ifndef falcON_included_cstring
#  include <cstring>
#  define falcON_included_cstring
#endif
#ifndef falcON_included_string
#  include <string>
#  define falcON_included_string
#endif
#ifndef falcON_included_basic_h
#  include <public/basic.h>
#endif

//------------------------------------------------------------------------------
namespace falcON {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // falcON::FileSize()                                                         
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
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
  // class falcON::iofile                                                       
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class iofile {
  protected:
    static const int FNAME_MAX_SIZE = 256;
    char             FNAME[FNAME_MAX_SIZE];
    const char      *FILE;
    iofile() : FILE(0) {}
    void setfile(const char*file) {
      if(file && file[0]) {
	strncpy(FNAME,file,FNAME_MAX_SIZE);
	FILE = FNAME;
      } else
	FILE = 0;
    }
  public:
    /// give file name, if any
    const char*const&file() const { return FILE; }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::output                                                       
  //                                                                            
  /// wrapper around std::ostream with some additional features                 
  ///                                                                           
  /// On opening a file, a filename "-" is interpreted as \c stdout, but will be
  /// opened only if no other output or nemo_out is writing to \c stdout.       
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
    bool   is_file     () const {  return OUT && OUT != &std::cout; }
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
    /// If \e file equals 0 or ".", nothing is opened.  If \e file equals "-",
    /// and no other output or nemo_out is opened to \c stdout, we map to
    /// \c stdout.  Otherwise, a file of name \e file is created for output. An
    /// existing file of the same name is deleted unless \e append is true, in
    /// which case, we append to that existing file.
    explicit
    output(const char*file, bool append=0) : FREC(0), APPENDING(false) {
      setfile(file);
      __open(append);
    }
    /// construction from file and potential option for appending.
    /// If \e file equals "-", and no other output or nemo_out is opened to 
    /// \c stdout, we map to \c stdout.  Otherwise, a file of name \e file is
    /// created for output. An existing file of the same name is deleted unless
    /// \e append is true, in which case, we append to that existing file.
    explicit
    output(std::string const&file, bool append=0) : FREC(0), APPENDING(false) {
      setfile(file.c_str());
      __open(append);
    }
    /// close possible old stream, then proceed as in construction         
    void open(const char*file, bool append = 0) {
      close();
      setfile(file);
      __open(append);
    }
    /// close possible old stream, then proceed as in construction         
    void open(std::string const&file, bool append = 0) {
      close();
      setfile(file.c_str());
      __open(append);
    }
    /// open file with name made from \e format string and \e tag.
    /// A new file name is created from the C-style \e format string and the
    /// data \e tag provided via \code sprintf(filename, format, tag) \endcode
    /// If this file name differs from the current, the old file is closed and
    /// the new one opened. The idea is to provide the possibility of numbered
    /// output files as in the following code \code
    /// output out;
    /// for(int i=0; i!=20; ++i) {
    ///   out.re_open("file%02d.dat",i);
    ///   out << i << std::endl;
    /// } \endcode creating the files \c file00.dat, \c file01.dat, ...
    /// \c file19.dat
    /// \return whether a new file has been opened (and the old closed)
    /// \param  format (input) C-style format string for file to open
    /// \param  tag    (input) datum needed in generating file to open
    /// \param  append (input, optional) append (or overwrite) existing file?
    template<typename T> 
    bool reopen(const char*format, T const&tag, bool append=0) {
      char FNEW[FNAME_MAX_SIZE];
      SNprintf(FNEW,FNAME_MAX_SIZE,format,tag);
      if(OUT==0 || strcmp(FNEW,FNAME)) {
	open(FNEW,append);
	return true;
      } else
	return false;
    }
    //@}
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
    output& operator<< (std::ostream::__ios_type &(*p)
                              (std::ostream::__ios_type &)) {
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
    static void open_std () falcON_THROWING;
    /// call if closing any output to stdout from NEMO main
    static void close_std();
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::input                                                        
  //                                                                            
  /// wrapper around std::istream with some additional features                 
  ///                                                                           
  /// On opening a file, a filename "-" is interpreted as \c stdin, but will be 
  /// opened only if no other input or nemo_in is reading from \c stdin.  The   
  /// operator >> is defined (as template) to map std::istream::operator>>(),   
  /// for other operations on std::istream, use the member method stream().     
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
    /// If \e file equals "-", and no other input or nemo_in is opened to \c
    /// stdin, we map to \c stdin.  Otherwise, a file of name \e file is opened
    /// for input.
    explicit
    input(const char*file) : FREC(0) {
      setfile(file);
      __open();
    }
    /// construction from file name.
    /// If \e file equals "-", and no other input or nemo_in is opened to \c
    /// stdin, we map to \c stdin.  Otherwise, a file of name \e file is opened
    /// for input.
    explicit
    input(std::string const&file) : FREC(0) {
      setfile(file.c_str());
      __open();
    }
    /// close possible old stream, then proceed as in construction
    void open(const char*file) {
      close();
      setfile(file);
      __open();
    }
    //@}
    /// close possible old stream, then proceed as in construction
    void open(std::string const&file) {
      close();
      setfile(file.c_str());
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
    /// input of manipulator \e m
    input& operator>> (std::istream& (*p)(std::istream&)) {
      if(IN) (*IN) >> p;
      return*this;
    }
    /// input of manipulator \e m
    input& operator>> (std::istream::__ios_type &(*p)
			      (std::istream::__ios_type &)) {
      if(IN) (*IN) >> p;
      return*this;
    }
    /// input of manipulator \e m
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
    static void open_std () falcON_THROWING;
    /// call if closing any output to stdout from NEMO main
    static void close_std();
  };
  // ///////////////////////////////////////////////////////////////////////////
  /// truncated string by one trailing character
  /// if the last character in name equals c then\n
  /// - we copy name into copy, except for character c\n
  /// - we return true\n
  /// otherwise\n
  /// - we return false
  bool is_appended(const char*name, char c, char*copy);
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::FortranIRec                                                  
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
    unsigned read_size() throw(falcON::exception);
  public:
    //--------------------------------------------------------------------------
    /// constructor: read buffer with size information
    /// \param in  falcON::input to read from
    /// \param rec (optional) size of Fortran record header: 4 or 8
    /// \param bswap (optional) swap bytes for size information?
    FortranIRec(input&in, unsigned rec=4, bool bswap=0)
      throw(falcON::exception);
    //--------------------------------------------------------------------------
    /// close: same as destruction
    void close() throw(falcON::exception);
    //--------------------------------------------------------------------------
    /// destructor: read to end of record, read end buffer
    ~FortranIRec() throw(falcON::exception) { close(); }
    //--------------------------------------------------------------------------
    /// read some bytes
    ///
    /// If more bytes are wanted than left in the record, only those left
    /// in the record will be read and a warning be issued.
    ///
    /// \return    number of bytes actually read
    /// \param buf buffer to read into
    /// \param n   number of bytes to read
    unsigned read_bytes(char*buf, unsigned n) throw(falcON::exception);
    //--------------------------------------------------------------------------
    /// read some data of any type
    ///
    /// If more data are wanted than left in the record, only those left
    /// in the record will be read and a warning be issued.
    ///
    /// \param T   data type
    /// \return    number of data actually read
    /// \param buf buffer to read into
    /// \param n   number of data to read
    template<typename T>
    unsigned read(T*buf, unsigned n) throw(falcON::exception) {
      if(READ+n*sizeof(T) > SIZE) {
	falcON_Warning("FortranIRec::read(): cannot read %d, but only %d %s\n",
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
    /// \param T   type of data to read
    /// \param in  input stream to read from
    /// \param buf data buffer to read into
    /// \param n   number of data of type T to read
    /// \param rec size of FORTRAN record header; must be 4 or 8
    template<typename T>
    static void Read(input &in, T*buf, unsigned n, unsigned rec=4)
      throw(falcON::exception)
    {
      FortranIRec FIR(in,rec);
      if( sizeof(T) * n > FIR.size() )
	throw exception("ReadFortranRecord(): cannot read %d %s: "
			"only %d bytes in record (required are %d)\n",
			n,nameof(T),FIR.size(),sizeof(T)*n);
      if( sizeof(T) * n < FIR.size() )
	falcON_Warning("ReadFortranRecord(): "
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
  // class falcON::FortranORec                                                  
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
    void write_size() throw(falcON::exception);
  public:
    //--------------------------------------------------------------------------
    /// constructor: write buffer with size information
    /// \param out output stream to write to
    /// \param size size (in bytes) of record
    /// \param rec (optional) size of Fortran record header must be 4 or 8
    FortranORec(output&out, unsigned size, unsigned rec=4)
      throw(falcON::exception);
    //--------------------------------------------------------------------------
    /// destructor: write to end of record, write end buffer
    ~FortranORec() throw(falcON::exception) { close(); }
    //--------------------------------------------------------------------------
    /// close: same as destruction
    void close() throw(falcON::exception);
    //--------------------------------------------------------------------------
    /// write some bytes
    ///
    /// If more bytes are to be written than left for the record, we only
    /// write as many as the record allows but issue a warning.
    ///
    /// \return number of bytes actually written
    /// \param buf buffer to write
    /// \param n   number of bytes to write
    unsigned write_bytes(const char*buf, unsigned n) throw(falcON::exception);
    //--------------------------------------------------------------------------
    /// fill some bytes with a given value
    ///
    /// \param n   number of bytes to fill
    /// \param val value to fill them with
    void fill_bytes(unsigned n, char val=0);
    //--------------------------------------------------------------------------
    /// write some data of any type
    ///
    /// If more data are to be written than left for the record, we only
    /// write as many as the record allows but issue a warning.
    ///
    /// \param T   data type
    /// \return    number of data actually written
    /// \param buf buffer to write
    /// \param n   number of data to write
    template<typename T>
    unsigned write(const T*buf, unsigned n) throw(falcON::exception) {
      if(WRITTEN + n*sizeof(T) > SIZE) {
	falcON_Warning("FortranORec::write(): "
		       "cannot write %d, but only %d %s\n",
		       n, (SIZE-WRITTEN)/sizeof(T), nameof(T));
	n = (SIZE-WRITTEN)/sizeof(T);
      }
      if(n) write_bytes(static_cast<const char*>
		       (static_cast<const void*>(buf)), sizeof(T)*n);
      return n;
    }
    //--------------------------------------------------------------------------
    /// write a single FORTRAN record in one go
    ///
    /// \param T   type of data to write
    /// \param out output stream to write to
    /// \param buf data buffer to write from
    /// \param n   number of data of type T to write
    /// \param rec size of FORTRAN record header; must be 4 or 8
    template<typename T>
    static void Write(output&out, const T*buf, unsigned n, unsigned rec=4)
      throw(falcON::exception)
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
    template<int B> struct __bswap {};
    template<> struct __bswap<1> {
      static void swap(void*vdat, unsigned cnt) { return; }
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
  /// \param bdat element to swap bytes for
  template<typename B> inline
  void swap_bytes(B&bdat) {
    __bswap<sizeof(B)>::swap(static_cast<void*>(&bdat), 1);
  }
  /// swap the bytes of elements of any type
  ///
  /// in order for this to work, the sizeof() the type must be 1,2,4,8, or 16
  ///
  /// \param bdat first element to swap bytes for
  /// \param cnt  number of elments to swap bytes for
  template<typename B> inline
  void swap_bytes(B*bdat, unsigned cnt) {
    __bswap<sizeof(B)>::swap(static_cast<void*>(bdat), cnt);
  }
  /// swap the bytes of elements of unknown type but known size
  ///
  /// \param vdat pointer to first element
  /// \param len  size of the elements, must be 1,2,4,8, or 16
  /// \param cnt  number of elments to swap bytes for
  inline
  void swap_bytes(void*vdat, unsigned len, unsigned cnt) falcON_THROWING {
    switch(len) {
    case  1: return;
    case  2: return __bswap< 2>::swap(vdat,cnt);
    case  4: return __bswap< 4>::swap(vdat,cnt);
    case  8: return __bswap< 8>::swap(vdat,cnt);
    case 16: return __bswap<16>::swap(vdat,cnt);
    default: falcON_THROW("swap_bytes(): sizeof(type)=%d: not supported\n",len);
    }
  }
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// structure modelled after gadget/allvars.h                                 
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  struct GadgetHeader {
    int          npart[6];          ///< # particles per type in this file
    double       masstab[6];        ///< if non-zero: mass of particle of type
    double       time;              ///< simulation time of snapshot
    double       redshift;          ///< redshift of snapshot
    int          flag_sfr;
    int          flag_feedback;
    unsigned int npartTotal[6];     ///< # particles per type in whole snapshot
    int          flag_cooling;
    int          num_files;         ///< # file for this snapshot
    double       BoxSize;
    double       Omega0;
    double       OmegaLambda;
    double       HubbleParam;
    int          flag_stellarage;
    int          flag_metals;
    unsigned int npartTotalHighWord[6];
    int          flag_entropy_instead_u;
    char         fill[60];          ///< to get sizeof(GadgetHeader)=256
    //--------------------------------------------------------------------------
    /// default constructor: set all data to 0
    GadgetHeader();
    //--------------------------------------------------------------------------
    /// try to read a GadgetHeader from an input file
    ///
    /// If the size of the Fortran record == 256 == sizeof(GadgetHeader), we
    /// read the header and return true.\n
    /// If the size of the Fortran record == byte_swapped(256), then we assume
    /// the file is of different endianess. We read the header, byte-swap it
    /// and return true.\n
    /// Otherwise, the data are not consistent with a GadgetHeader, so we return
    /// false.
    /// \return have read successfully
    /// \param  in   input stream to read from
    /// \param  rec  size of Fortran record header (must be 4 or 8)
    /// \param  swap (output) need byte-swap?
    bool Read(input& in, unsigned rec, bool& swap)
      throw(falcON::exception);
    //--------------------------------------------------------------------------
    /// check whether two GadgetHeaders could possibly come from different data
    /// files for the same snapshot
    bool mismatch(GadgetHeader const&H) const;
    //--------------------------------------------------------------------------
    /// on some ICs, npartTotal[] = 0. Here we remedy for this error
    void check_simple_npart_error();
    //--------------------------------------------------------------------------
    /// dump all the header data
    void dump(std::ostream&out) const;
  };
  // ///////////////////////////////////////////////////////////////////////////
} // namespace falcON {
// /////////////////////////////////////////////////////////////////////////////
falcON_TRAITS(falcON::output,"output");
falcON_TRAITS(falcON::input,"input");
falcON_TRAITS(falcON::FortranIRec,"FortranIRec");
falcON_TRAITS(falcON::FortranORec,"FortranORec");
falcON_TRAITS(falcON::GadgetHeader,"GadgetHeader");
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_io_h
