// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   inc/public/io.h                                                     
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2000-2006                                                           
///                                                                             
/// \brief  contains declarations of classes falcON::input and falcON::output,  
///         as well as NEMO I/O support                                         
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2000-2005 Walter Dehnen                                        
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

#ifndef falcON_included_cstdio
#  include <cstdio>
#  define falcON_included_cstdio
#endif
#ifndef falcON_included_cstring
#  include <cstring>
#  define falcON_included_cstring
#endif
#ifndef falcON_included_types_h
#  include <public/fields.h>
#endif

#ifdef  falcON_NEMO
//------------------------------------------------------------------------------
namespace falcON {
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
  class output {
    static const int FNAME_MAX_SIZE = 200;
    std::ostream *OUT;
    char          FNAME[FNAME_MAX_SIZE];
    char         *FILE;
    bool          APPENDING;
    //--------------------------------------------------------------------------
    void __open (bool);
    void __close();
  public:
    /// \name const boolean information                                         
    //@{
    bool   is_open     () const { return OUT != 0; }   ///< ready for output?
    bool   is_appending() const { return APPENDING; }  ///< appending output?
    operator bool      () const { return OUT != 0; }   ///< ready for output?
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
    /// \name construction and opening                                          
    //@{
    /// construction from nothing: nothing is opened
    output() : OUT(0) , FILE(0), APPENDING(false) {}
    /// construction from file name and potential option for appending.
    /// If \e file equals 0 or ".", nothing is opened.  If \e file equals "-",
    /// and no other output or nemo_out is opened to \c stdout, we map to
    /// \c stdout.  Otherwise, a file of name \e file is created for output. An
    /// existing file of the same name is deleted unless \e append is true, in
    /// which case, we append to that existing file.
    explicit
    output(const char*file, bool append = 0) : APPENDING(false) {
      strncpy(FNAME,file,FNAME_MAX_SIZE);
      FILE = file? FNAME : 0;
      __open(append);
    }
    /// close possible old stream, then proceed as in construction         
    void open(const char*file, bool append = 0) {
      __close();
      strncpy(FNAME,file,FNAME_MAX_SIZE);
      FILE = file? FNAME : 0;
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
      snprintf(FNEW,FNAME_MAX_SIZE,format,tag);
      if(OUT==0 || strcmp(FNEW,FNAME)) {
	open(FNEW,append);
	return true;
      } else
	return false;
    }
    //@}
    //--------------------------------------------------------------------------
    /// \name destruction and closing                                           
    //@{
    /// close any open files; \c stdout is freed for other output
    void close() {
      __close();
      OUT = 0;
    }
    /// like close()
    ~output() { __close(); }
    //@}
    //--------------------------------------------------------------------------
    /// return name of output file if any
    const char* file() const { return FILE; }
    /// flush() output
    void flush() { if(OUT) OUT->flush(); }
    //--------------------------------------------------------------------------
    /// \name formated output using operator <<
    //@{
    /// output of a single datum \e x
    template<typename X>
    std::ostream& operator<< (X const&x) {
      return (*OUT) << x;
    }
    /// output of a manipulator \e m
    std::ostream& operator<< (std::ostream& (*p)(std::ostream&)) {
      return (*OUT) << p;
    }
    /// output of a manipulator \e m
    std::ostream& operator<< (std::ostream::__ios_type &(*p)
                              (std::ostream::__ios_type &)) {
      return (*OUT) << p;
    }
    /// output of a manipulator \e m
    std::ostream& operator<< (std::ios_base& (*p)(std::ios_base&)) {
      return (*OUT) << p;
    }
    //@}
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
  class input {
    std::istream  *IN;
    void __open(const char*);
    void __close();
  public:
    /// \name const boolean information                                         
    //@{
    bool   is_open() const { return IN != 0; }         ///< ready for input?
    operator bool () const { return IN != 0; }         ///< ready for input?
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
    /// \name construction and opening                                          
    //@{
    /// construction from nothing: nothing is opened
    input() : IN(0) {}
    /// construction from file name.
    /// If \e file equals "-", and no other input or nemo_in is opened to \c
    /// stdin, we map to \c stdin.  Otherwise, a file of name \e file is opened
    /// for input.
    explicit
    input(const char*file) { __open(file); }
    /// close possible old stream, then proceed as in construction
    void open(const char*file) {
      __close();
      __open(file);
    }
    //@}
    //--------------------------------------------------------------------------
    /// \name destruction and closing
    //@{
    /// if \c stdin: clear \c stdin, otherwise delete ifstream
    void close() {
      __close();
      IN = 0;
    }
    /// like close()
    ~input() { __close(); }
    //@}
    //--------------------------------------------------------------------------
    /// \name formated input using operator >>
    //@{
    /// input of single datum \e x
    template<typename X>
    std::istream& operator>> (X&x) {
      return (*IN) >> x;
    }
    /// input of manipulator \e m
    std::istream& operator>> (std::istream& (*p)(std::istream&)) {
      return (*IN) >> p;
    }
    /// input of manipulator \e m
    std::istream& operator>> (std::istream::__ios_type &(*p)
                              (std::istream::__ios_type &)) {
      return (*IN) >> p;
    }
    /// input of manipulator \e m
    std::istream& operator>> (std::ios_base& (*p)(std::ios_base&)) {
      return (*IN) >> p;
    }
    //@}
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::nemo_io                                                      
  //                                                                            
  /// base class for \a nemo_in and \a nemo_out                                 
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class nemo_io {
    //--------------------------------------------------------------------------
    //                                                                          
    // enums to identify the                                                    
    // - type of basic data                                                     
    // - type of data structure                                                 
    // - meaning of data                                                        
    //                                                                          
    //--------------------------------------------------------------------------
  public:
    //--------------------------------------------------------------------------
    //                                                                          
    // enum DataType                                                            
    //                                                                          
    /// used to identify different nemo I/O able data types                     
    //                                                                          
    //--------------------------------------------------------------------------
    enum DataType {
      Null    = 0,
      Short   = 1,
      Integer = 2,
      Long    = 3,
      Single  = 4,
      Double  = 5,
#ifdef falcON_REAL_IS_FLOAT
      Real    = Single,
      NotReal = Double
#else
      Real    = Double,
      NotReal = Single
#endif
    };
    //--------------------------------------------------------------------------
    /// returns the name of a Datatype as string
    static const char* type_name(DataType t) {
      switch(t) {
      case Null: return "Null";
      case Short: return "Short";
      case Integer: return "Integer";
      case Single: return "Single";
      case Double: return "Double";
      default: return "Unknown";
      }
    }
    //--------------------------------------------------------------------------
    /// two DataTypes are coercing if they are of different floating point type
    static bool coercing(DataType t1, DataType t2) {
      return 
	(t1 == Double && t2 == Single) || (t2 == Double && t1 == Single);
    }
    //--------------------------------------------------------------------------
    /// different Datatypes are mismatched, unless they are coercing
    static bool mismatch(DataType t1, DataType t2) {
      return t1!=t2 && !coercing(t1,t2);
    }
    //--------------------------------------------------------------------------
    //                                                                          
    // enum Field                                                               
    //                                                                          
    /// identifiers for the various data fields supported by nemo_io            
    //                                                                          
    //--------------------------------------------------------------------------
    enum Field {
      null    = 0,
      mass    = 1 <<  0,
      pos     = 1 <<  1,
      vel     = 1 <<  2,
      eps     = 1 <<  3,
      key     = 1 <<  4,
      tau     = 1 <<  5,
      pot     = 1 <<  6,
      acc     = 1 <<  7,
      jerk    = 1 <<  8,
      dens    = 1 <<  9,
      aux     = 1 << 10,
      lev     = 1 << 11,
      num     = 1 << 12,
      posvel  = 1 << 13,
      SPHh    = 1 << 14,
      SPHnum  = 1 << 15,
      SPHu    = 1 << 16,
      SPHudot = 1 << 17,
      SPHurad = 1 << 18,
      SPHentr = 1 << 19,
      SPHdens = 1 << 20,
      SPHhdot = 1 << 21,
      SPHfact = 1 << 22,
      SPHcs   = 1 << 23
    };
    //--------------------------------------------------------------------------
    /// get the nemo_io::DataType used for a given Field                        
    static DataType type(Field f) {
      switch(f) {
      case mass   : return Real;
      case pos    : return Real;
      case vel    : return Real;
      case eps    : return Real;
      case key    : return Integer;
      case tau    : return Real;
      case pot    : return Real;
      case acc    : return Real;
      case jerk   : return Real;
      case dens   : return Real;
      case aux    : return Real;
      case lev    : return Short;
      case num    : return Integer;
      case posvel : return Real;
      case SPHh   : return Real;
      case SPHnum : return Integer;
      case SPHu   : return Real;
      case SPHudot: return Real;
      case SPHurad: return Real;
      case SPHentr: return Real;
      case SPHdens: return Real;
      case SPHhdot: return Real;
      case SPHfact: return Real;
      case SPHcs  : return Real;
      default     : return Null;
      }
    }
    //--------------------------------------------------------------------------
    /// is a Field SPH?                                                         
    static bool is_sph(Field f) {
      switch(f) {
      case SPHh:
      case SPHnum:
      case SPHu:
      case SPHudot:
      case SPHurad:
      case SPHentr:
      case SPHdens:
      case SPHhdot:
      case SPHfact:
      case SPHcs:   return true;
      default:      return false;
      }
    }
    //--------------------------------------------------------------------------
    /// what \a Field corresponds to a given \a fieldbit                        
    static Field field(fieldbit f)
    {
      switch(value(f)) {
      case fieldbit::m: return mass;
      case fieldbit::x: return pos;
      case fieldbit::v: return vel;
      case fieldbit::e: return eps;
      case fieldbit::k: return key;
      case fieldbit::t: return tau;
      case fieldbit::p: return pot;
      case fieldbit::a: return acc;
      case fieldbit::j: return jerk;
      case fieldbit::r: return dens;
      case fieldbit::y: return aux;
      case fieldbit::l: return lev;
      case fieldbit::n: return num;
      case fieldbit::H: return SPHh;
      case fieldbit::N: return SPHnum;
      case fieldbit::U: return SPHu;
      case fieldbit::I: return SPHudot;
      case fieldbit::E: return SPHurad;
      case fieldbit::S: return SPHentr;
      case fieldbit::R: return SPHdens;
      case fieldbit::J: return SPHhdot;
      case fieldbit::F: return SPHfact;
      case fieldbit::C: return SPHcs;
      default         : return null;
      }
    }
    //--------------------------------------------------------------------------
    // what \a fieldbit corresponds to a given \a Field                         
    static fieldbit bit(Field F) falcON_THROWING
    {
      switch(F) {
      case mass:    return fieldbit::m;
      case pos:     return fieldbit::x;
      case vel:     return fieldbit::v;
      case eps:     return fieldbit::e;
      case key:     return fieldbit::k;
      case tau:     return fieldbit::t;
      case pot:     return fieldbit::p;
      case acc:     return fieldbit::a;
      case jerk:    return fieldbit::j;
      case dens:    return fieldbit::r;
      case aux:     return fieldbit::y;
      case lev:     return fieldbit::l;
      case num:     return fieldbit::n;
      case SPHh:    return fieldbit::H;
      case SPHnum:  return fieldbit::N;
      case SPHu:    return fieldbit::U;
      case SPHudot: return fieldbit::I;
      case SPHurad: return fieldbit::E;
      case SPHentr: return fieldbit::S;
      case SPHdens: return fieldbit::R;
      case SPHhdot: return fieldbit::J;
      case SPHfact: return fieldbit::F;
      case SPHcs:   return fieldbit::C;
      default: falcON_THROW("unaccountable nemo_io::Field\n");
      }
    }
    //--------------------------------------------------------------------------
    /// given a collection of \a Field s, return the corresponding fieldset     
    static fieldset fields(nemo_io::Field F)
    {
      fieldset data;
      for(fieldbit f; f; ++f)
	if(field(f) & F) data |= fieldset(f);
      return data;
    }
    //--------------------------------------------------------------------------
    //                                                                          
    // private data of class nemo_io                                            
    //                                                                          
    //--------------------------------------------------------------------------
  protected:
    void *STREAM;
    //--------------------------------------------------------------------------
    //                                                                          
    // protected member methods                                                 
    //                                                                          
    //--------------------------------------------------------------------------
    nemo_io&open    (const char*, const char*);    ///< open new file, close old
    void    close   ();                            ///< close open file         
    nemo_io() : STREAM (0) {}
    ~nemo_io        () { close(); }                ///< close file & free memory
  public:
    bool    is_open () const { return STREAM!=0; } ///< ready for output ?      
    operator bool   () const { return is_open(); } ///< ready for output ?      
  };// class nemo_io
  class snap_in;
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::nemo_in                                                      
  //                                                                            
  /// represents a NEMO input stream, derived from \a nemo_io                   
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class nemo_in : public nemo_io {
    friend class snap_in;
    //--------------------------------------------------------------------------
    // data are all private                                                     
    //--------------------------------------------------------------------------
  private:
    mutable snap_in *SNAP_IN;                      // if non-zero: open snapshot
    bool             IS_PIPE;                      // are we a pipe?            
    nemo_in(nemo_in const&);                       // not implemented           
    void*const& stream() const { return STREAM; }  // our nemo I/O stream       
    //--------------------------------------------------------------------------
    // public member methods                                                    
    //--------------------------------------------------------------------------
  public:
    /// construction from nothing: do nothing yet
    nemo_in() {}
    //--------------------------------------------------------------------------
    /// construction from file name: open file for NEMO input
    ///
    /// If \a file equals "-", set pipe, ie. read from stdin instead from a file
    ///
    /// \param file (input) name of file to be opened
    /// \param mode (input, optional) open mode, see nemo documentation
    explicit nemo_in(const char*file,
		     const char*mode="r");
    //--------------------------------------------------------------------------
    /// destruction: close stream; no \a snap_in must be opened
    ~nemo_in() falcON_THROWING;
    //--------------------------------------------------------------------------
    /// close any old file and open new file
    ///
    /// \param file (input) name of file to be opened
    nemo_in&open(const char*file) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// close open input stream
    void close() falcON_THROWING;
    //--------------------------------------------------------------------------
    /// can a snapshot be opened?
    ///
    /// \return true if a \a snap_in can be constructed from *this
    bool has_snapshot() const;
    //--------------------------------------------------------------------------
    /// are we a pipe?
    bool const&is_pipe() const { return IS_PIPE; }
  };
  class data_in;
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::snap_in                                                      
  //                                                                            
  /// represents a snapshot on a NEMO input stream, \a nemo_in                  
  ///                                                                           
  /// At any time, only one \a snap_in can exist for any \a nemo_in.            
  /// Construction of a second \a snap_in from the same \a nemo_in will cause   
  /// a fatal error. A \a \snap_in provides information on which data are       
  /// contained in the snapshot and allows to construct \a data_in for          
  /// reading them.                                                             
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class snap_in {
    friend class data_in;
    //--------------------------------------------------------------------------
    // data are all private                                                     
    //--------------------------------------------------------------------------
  private:
    mutable data_in *DATA_IN;                      // if non-zero: open snapshot
    nemo_in const   &INPUT;                        // our input stream          
    mutable int      FIELDS_READ;                  // fields read already       
    bool             HAS_NSPH, HAS_TIME;           // do we have these data?    
    unsigned         NBOD, NSPH;                   // # bodies, # SPH bodies    
    double           TIME;                         // simulations time          
    //--------------------------------------------------------------------------
    void*const& stream() const {                   // our nemo stream           
      return INPUT.stream();
    }
    //--------------------------------------------------------------------------
    // public member methods                                                    
    //--------------------------------------------------------------------------
  public:
    /// only constructor: open a snapshot set from a NEMO input stream
    ///
    /// \param in (input) NEMO input stream
    explicit snap_in(nemo_in const&in) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// destruction: close snapshot set in NEMO input stream
    ~snap_in() falcON_THROWING;
    //--------------------------------------------------------------------------
    /// do we have # SPH bodies?
    bool has_Nsph() const { return HAS_NSPH; }
    /// do we have simul time?  
    bool has_time() const { return HAS_TIME; }
    /// return # bodies         
    unsigned const&Nbod() const { return NBOD; }
    /// return # SPH bodies     
    unsigned const&Nsph() const { return NSPH; }
    /// return the number of data expected for a given field
    unsigned       N(nemo_io::Field F) const {
      return nemo_io::is_sph(F)? NSPH : NBOD;
    }
    /// return the simulation time of the snapshot
    double   const&time() const { return TIME; }
    //--------------------------------------------------------------------------
    /// do we have data for the field \e f ?
    bool has(nemo_io::Field f) const;
    /// have data for the field \e f already been read ?
    bool has_been_read(nemo_io::Field f) const {
      return FIELDS_READ & f;
    }
    /// return: which fields have already been read from this snapshot
    int  const&fields_read() const { 
      return FIELDS_READ;
    }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::data_in                                                      
  //                                                                            
  /// represents the data of one field from one snapshot of a NEMO input stream 
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class data_in {
    //--------------------------------------------------------------------------
    // data are all private                                                     
    //--------------------------------------------------------------------------
  private:
    snap_in const          &INPUT;                 // our snapshot              
    const nemo_io::Field    FIELD;                 // sort of data we are for   
    unsigned                NREAD;                 // how many have been read?  
    unsigned                NTOT;                  // how many are there?       
    nemo_io::DataType       TYPE;                  // which data type?          
    unsigned                SUBN;                  // how many items per datum  
    static const int        NDIM = Ndim;           // # spatial dimensions      
    //--------------------------------------------------------------------------
    // public member methods                                                    
    //--------------------------------------------------------------------------
  public:
    //--------------------------------------------------------------------------
    /// construction: create a valid data input with the expected number of     
    ///               bodies and no type mismatch; otherwise throw an exception 
    ///
    /// \param s (input) input snapshot to be read from
    /// \param f (input) data field to be read
    data_in(snap_in const &s,
	    nemo_io::Field f) falcON_THROWING;
    /// destruction: close data set
    ~data_in();
    /// which data field are we reading?
    nemo_io::Field     const&field()      const { return FIELD; }
    /// which data type are we reading?
    nemo_io::DataType  const&type()       const { return TYPE; }
    /// to we need to coerce between float and double?
    bool must_coerce() const { 
      return nemo_io::coercing(TYPE, nemo_io::type(FIELD));
    }
    /// # items per datum (i.e. 1 for scalars, 3 for vectors)
    unsigned           const&sub_N()      const { return SUBN; }
    /// total # data present in input snapshot
    unsigned           const&N()          const { return NTOT; }
    /// # data already read from input
    unsigned           const&N_read()     const { return NREAD; }
    /// # data yet to be read from input
    unsigned                 N_unread()   const { 
      return NTOT>NREAD? NTOT-NREAD : 0u;
    }
    //--------------------------------------------------------------------------
    /// read all data into array given
    /// \param data (input) pointer to data array, must have sufficient memory
    void read(void*data);
    /// read only \e n data into array given
    /// \param data (input) pointer to data array, must have sufficient memory
    /// \param n (input) number of bodies to read data for
    void read(void*data, unsigned n);
  };
  class snap_out;
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::nemo_out                                                     
  //                                                                            
  /// represents a NEMO output stream, derived from \a nemo_io                  
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class nemo_out : public nemo_io {
    friend class snap_out;
    //--------------------------------------------------------------------------
    // data are all private                                                     
    //--------------------------------------------------------------------------
  private:
    mutable snap_out *SNAP_OUT;                    // if non-zero open snapshot 
    bool              IS_PIPE, IS_SINK;            // are we a pipe or sink?    
    nemo_out(nemo_out const&);                     // not implemented           
    void*const& stream() const { return STREAM; }  // our nemo I/O stream       
    //--------------------------------------------------------------------------
    // public member methods                                                    
    //--------------------------------------------------------------------------
  public:
    /// construction from nothing: do nothing yet
    nemo_out() {}
    //--------------------------------------------------------------------------
    /// open a nemo_out, special options by filename                            
    ///                                                                         
    /// If \e file equals "-", and no output or other nemo_out writes to
    /// \c stdout, then output will be made to \c stdout.
    ///
    /// If \e file equals ".", no output will be made at all (sink).
    ///
    /// Otherwise, if no file \e file exists already, a new file of that name
    /// is created for output. If such a file exists already, then usual NEMO
    /// behaviour is to abort with an error message. Here, we extend that in
    /// the following way: If the last character of \e file is '!', an
    /// existing file of the same name (except for the trailing '!') will be
    /// overwritten. Similarly, if the last character of \e file is '@', an
    /// existing file of the same name (except for the trailing '!') will be
    /// appended to. If no file of that name exists, we proceed as usual and
    /// open a file of name \e file (omitting any trailing '!' or '@').
    ///
    /// \return *this
    /// \param file (input) name of file to be opened, see detailled comments
    /// \param append (input) append existing file anyway?
    nemo_out&open(const char* file,
		  bool append = false) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// close open output stream
    void close() falcON_THROWING;
    //--------------------------------------------------------------------------
    /// open a nemo_out, special options by filename, see nemo_out::open()      
    ///                                                                         
    /// \param file (input) name of file to be opened, see detailled comments
    /// \param append (input) append existing file anyway?
    explicit nemo_out(const char*file,
		      bool       append = false) falcON_THROWING
    : IS_PIPE(0), IS_SINK(0), SNAP_OUT(0) { open(file,append); }
    //--------------------------------------------------------------------------
    /// close open output stream
    ~nemo_out() falcON_THROWING { close(); }
    //--------------------------------------------------------------------------
    /// are we a pipe?
    bool const&is_pipe() const { return IS_PIPE; }
    /// are we a sink?
    bool const&is_sink() const { return IS_SINK; }
    /// are we open and not a sink?
    bool is_writing   () const { return is_open() && !is_sink(); }
    /// convertion to bool: true if is_writing()
    operator bool     () const { return is_writing(); }
  };
  class data_out;
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::snap_out                                                     
  //                                                                            
  /// represents a snapshot on a NEMO output stream, \a nemo_out                
  ///                                                                           
  /// At any time, only one \a snap_out can exist for any \a nemo_out.          
  /// Construction of a second \a snap_out from the sane \a nemo_out will cause 
  /// a fatal error. A \a \snap_out allows to construct \a data_out for         
  /// writing data.                                                             
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class snap_out {
    friend class data_out;
    //--------------------------------------------------------------------------
    // data are all private                                                     
    //--------------------------------------------------------------------------
  private:
    nemo_out const   &OUTPUT;                      // our output stream         
    mutable data_out *DATA_OUT;                    // if non-zero: open data_out
    mutable int       FIELDS_WRITTEN;              // data already written out  
    unsigned          NBOD, NSPH;                  // # bodies, # SPH bodies    
    //--------------------------------------------------------------------------
    void*const& stream() const {                   // our nemo stream           
      return OUTPUT.stream();
    }
    //--------------------------------------------------------------------------
    // public member methods                                                    
    //--------------------------------------------------------------------------
  public:
    /// construction: open NEMO snapshot set
    ///
    /// \param out  (input) NEMO output stream
    /// \param Nbod (input) total # bodies
    /// \param Nsph (input) # SPH bodies (must come first)
    /// \param time (input) simulation time
    snap_out(nemo_out const&out,
	     unsigned       Nbod,
	     unsigned       Nsph,
	     double         time) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// destruction: close snapshot set in NEMO output stream
    ~snap_out() falcON_THROWING;
    //--------------------------------------------------------------------------
    /// return total # bodies
    unsigned const&Nbod() const { return NBOD; }
    /// return # SPH bodies
    unsigned const&Nsph() const { return NSPH; }
    /// have data for given field been written already?
    bool has_been_written(nemo_io::Field f) const {
      return FIELDS_WRITTEN & f;
    }
    /// for which fields have data been written yet?
    int const&fields_written() const {
      return FIELDS_WRITTEN;
    }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::data_out                                                     
  //                                                                            
  /// represents the data of one field from one snapshot of a NEMO output stream
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class data_out {
    //--------------------------------------------------------------------------
    // data are all private                                                     
    //--------------------------------------------------------------------------
  private:
    snap_out const        &OUTPUT;                 // our snapshot              
    const nemo_io::Field   FIELD;                  // sort of data we are for   
    unsigned               NWRITTEN;               // how many have been written
    const unsigned         NTOT;                   // how many in total         
    nemo_io::DataType      TYPE;                   // which data type?          
    unsigned               SUBN;                   // how many items per datum  
    static const int       NDIM = Ndim;            // # spatial dimensions      
    //--------------------------------------------------------------------------
  public:
    /// construction: open a NEMO data set
    ///
    /// \param s (input) our \a snap_out snapshot set to write to
    /// \param f (input) the data field to write out with this
    data_out(snap_out const&s,
	     nemo_io::Field f) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// destruction: close data set
    ~data_out();
    //--------------------------------------------------------------------------
    /// return out data field
    nemo_io::Field     const&field()      const { return FIELD; }
    /// return total # bodies to write (=Nbod or =Nsph)
    unsigned           const&N()          const { return NTOT; }
    /// return data type to be written out
    nemo_io::DataType  const&data_type()  const { return TYPE; }
    /// return # bodies for which data have been written out already
    unsigned           const&N_written()  const { return NWRITTEN; }
    /// return # bodies for which data have yet to be written.
    unsigned                 N_free()     const { 
      return NTOT>NWRITTEN? NTOT-NWRITTEN : 0u;
    }
    /// write all data from array given (which must hold enough data)
    ///
    /// \param data (input) pointer to data to be written
    void write(const void*data);
    /// write data for \e n bodies from array given (must hold enough data)
    ///
    /// \param data (input) pointer to data to be written
    /// \param n    (input) # bodies for which to write data
    void write(const void*data, unsigned n);
  };
  //////////////////////////////////////////////////////////////////////////////
  /// is a given time in a given time range?
  ///
  /// \param t (input) simulation time
  /// \param r (input) time range as string, using NEMO notation
  bool time_in_range(double t, const char*r);
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_NEMO
#endif // falcON_included_io_h
