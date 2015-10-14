// -*- C++ -*-
// /////////////////////////////////////////////////////////////////////////////
//
/// \file    inc/public/nemo++.h
//
/// \brief   provides an interfact to any NEMO library functionality
/// \author  Walter Dehnen
/// \date    2008-2009
//
// /////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008-2009  Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 675
// Mass Ave, Cambridge, MA 02139, USA.
//
// /////////////////////////////////////////////////////////////////////////////
#ifndef falcON_included_nemopp_h
#define falcON_included_nemopp_h 1

#ifdef falcON_NEMO

#ifndef __cplusplus
#  error "nemo++.h must only be used in C++"
#endif
#ifndef falcON_included_fields_h
#  include <public/fields.h>
#endif
#ifndef falcON_included_cstdio
#  include <cstdio>
#  define falcON_included_cstdio
#endif

namespace falcON {

  // ///////////////////////////////////////////////////////////////////////////
  /// debugging level
  /// \return depth <= debugging level
  int nemo_debug_level();

  // ///////////////////////////////////////////////////////////////////////////
  /// \name expression parsing
  //@{
  /// parse an expression into an array of values, see "man 3 herinp"
  /// \param T (template parameter) type of value
  /// \return     number of array elements set, < 0 if parse error
  /// \param[in]  e expression
  /// \param[out] x array with values
  /// \param[in]  n size of array
  /// \note implementations exist for int,unsigned,long,bool,float,double
  /// \note implemented in nemo++.cc using NEMO
  template<typename T>
  int nemoinp(const char*e, T*x, int n);

  /// parse an expression of the format d:m:s.s[,d:m:s.s[, ... ]]
  /// \return     number of array elements set, < 0 if parse error
  /// \param[in]  e expression of format d:m:s.s[,d:m:s.s[, ... ]]
  /// \param[out] x array with values
  /// \param[in]  n size of array
  /// \note implementations exist for int,double,float,long,bool
  /// \note implemented in nemo++.cc using NEMO's nemoinpx
  int nemoinpx(const char*e, double*x, int n);
  //@}

  // ///////////////////////////////////////////////////////////////////////////
  /// \name parameter input (you don't need to include NEMO's getparam.h)
  //@{
  
  /// initialise nemo run-time parameter
  /// \param[in] argv argument list given to executable
  /// \param[in] defv list specifying parameters
  /// \note implemented in nemo++.cc via NEMO's initparam()
  void initparam(const char**argv, const char**defv);

  /// finalise parameter
  /// \note implemented in nemo++.cc via NEMO's finiparam()
  void finiparam();

  /// does a parameter of given name a value?
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \return true if p has a value associated with it
  /// \note implemented in nemo++.cc via NEMO's hasvalue()
  bool hasvalue(const char*p);

  /// get character string
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \return value associated with p
  /// \note implemented in nemo++.cc via NEMO's getparam()
  const char*getparam(const char*p);

  /// get character string, default to empty string
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \return value associated with p if any, otherwise 0
  inline const char* getparam_z(const char*p) {
    return hasvalue(p)? getparam(p) : 0;
  }

  /// get bool
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \return value associated with p
  /// \note implemented in nemo++.cc via NEMO's getbparam()
  bool getbparam(const char*p);

  /// get bool, default to false
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \return value associated with p if any, otherwise false
  inline bool getbparam_z(const char*p) {
    return hasvalue(p)? getbparam(p) : false;
  }

  /// get int
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \return value associated with p
  /// \note implemented in nemo++.cc via NEMO's getiparam()
  int getiparam(const char*p);

  /// get int, default to 0
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \return value associated with p if any, otherwise 0
  inline int getiparam_z(const char*p) {
    return hasvalue(p)? getiparam(p) : 0;
  }

  /// get unsigned
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \return value associated with p
  inline unsigned getuparam(const char* p) {
    int i = getiparam(p);
    if(i<0) falcON_THROW("parameter %s requires positive integer, got %d\n",
			 p,i);
    return unsigned(i);
  }

  /// get unsigned, default to 0u
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \return value associated with p if any, otherwise 0u
  inline unsigned getuparam_z(const char*p) {
    return hasvalue(p)? getuparam(p) : 0u;
  }

  /// get long
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \return value associated with p
  /// \note implemented in nemo++.cc via NEMO's getlparam()
  long getlparam(const char*p);

  /// get long, default to 0
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \return value associated with p if any, otherwise 0
  /// \note implemented in nemo++.cc via NEMO's getlparam()
  inline long getlparam_z(const char*p) {
    return hasvalue(p)? getlparam(p) : 0;
  }

  /// get double
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \return value associated with p
  /// \note implemented in nemo++.cc via NEMO's getdparam()
  double getdparam(const char*p);

  /// get double, default to 0.0
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \return value associated with p if any, otherwise 0
  /// \note implemented in nemo++.cc via NEMO's getdparam()
  inline double getdparam_z(const char*p) {
    return hasvalue(p)? getdparam(p) : 0.;
  }

  /// get float
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \return value associated with p
  inline float getfparam(const char*p) {
    return static_cast<float>(getdparam(p));
  }

  /// get float, default to 0.0f
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \return value associated with p if any, otherwise 0
  /// \note implemented in nemo++.cc via NEMO's getdparam()
  inline float getfparam_z(const char*p) {
    return hasvalue(p)? getfparam(p) : 0.f;
  }

  /// get real
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \return value associated with p
  inline real getrparam(const char*p)  {
    return static_cast<real>(getdparam(p));
  }

  /// get real, default to 0
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \return value associated with p if any, otherwise 0
  inline real getrparam_z(const char*p) {
    return hasvalue(p)? getrparam(p) : zero;
  }

  /// get a fieldset
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \return value associated with p
  inline fieldset getioparam(const char*p) {
    return fieldset(getparam(p));
  }

  /// get a fieldset, default to fieldset::empty
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \return value associated with p if any, fieldset::empty
  inline fieldset getioparam_z(const char*p) {
      return hasvalue(p)? getioparam(p) : fieldset::empty; 
  }

  /// get a fieldset, default to fieldset::all
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \return value associated with p if any, fieldset::all
  inline fieldset getioparam_a(const char*p) {
      return hasvalue(p)? getioparam(p) : fieldset::all; 
  }

  /// get vect
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \note implemented in nemo++.cc using NEMO
  vect getvparam(const char*p);

  /// get vect, but allow single value, i.e. (0,0,0) if value of p is "0".
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \note implemented in nemo++.cc using NEMO
  vect getvrparam(const char*p);

  /// read vect from nemo parameter into 2nd argument, return pointer
  /// \return on success: pointer to vect, otherwise: NULL
  /// \param[in]  p parameter name as specified in 2nd argument to initparam()
  /// \param[out] x if returning not NULL, this vector is assigned to
  /// \note implemented in nemo++.cc using NEMO
  vect* getvparam_z(const char* p, vect&x);

  /// read vect from nemo param into 2nd arg, allow single value, return pter
  /// \return on success: pointer to vect, otherwise: NULL
  /// \param[in]  p parameter name as specified in 2nd argument to initparam()
  /// \param[out] x if returning not NULL, this vector is assigned to
  /// \note implemented in nemo++.cc using NEMO
  /// \note implemented in nemo++.cc using NEMO
  vect* getvrparam_z(const char* p, vect&x);

  /// read array of type T
  /// \param[in]  p parameter name as specified in 2nd argument to initparam()
  /// \param[out] a  array
  /// \param[in]  m  physical size of array
  /// \return     number of elements actually read.
  /// \note implemented in nemo++.cc using NEMO
  /// \note implementations for int,unsigned,bool,double,float
  template<typename T>
  inline int getaparam(const char*p, T*a, int m) {
    return nemoinp(getparam(p),a,m);
  }

  /// read array of type T, but allow for non-existence
  /// \param[in]  p parameter name as specified in 2nd argument to initparam()
  /// \param[out] a  array
  /// \param[in]  m  physical size of array
  /// \return     number of elements actually read.
  /// \note implemented in nemo++.cc using NEMO
  /// \note implementations for int,unsigned,bool,double,float
  template<typename T>
  int getaparam_z(const char*p, T*a, int m) {
    if(!hasvalue(p)) {
      for(int i=0; i!=m; ++i) a[i] = T(0);
      return 0;
    } else 
      return nemoinp(getparam(p),a,m);
  }

#ifdef falcON_included_PotExp_h
  /// get PotExp::symmetry
  /// \param[in] p parameter name as specified in 2nd argument to initparam()
  /// \return value associated with p
  inline PotExp::symmetry getsymparam(const char*p) {
    int _sym (getiparam(p));
    return
      _sym==4? PotExp::spherical   :
      _sym==3? PotExp::cylindrical :
      _sym==2? PotExp::triaxial    :
      _sym==1? PotExp::reflexion   : PotExp::none;
  }
#endif

  //@}
  // ///////////////////////////////////////////////////////////////////////////
  /// \name support for dynamically loading code
  //@{

  /// get full name of a file
  /// \param[in] file file name without path
  /// \return full file name, including full path
  /// \note implemented in nemo++.cc using NEMO's fullname()
  const char*fullname(const char*file);

  /// search for file in path
  /// \param[in] path path to search for (colon separated list)
  /// \param[in] file name of file to search for
  /// \return full file name of first match found in path, if any
  /// \note implemented in nemo++.cc using NEMO's pathfind()
  const char* pathfind(const char*path, const char*file);

  /// load all symbols from program
  /// \param[in] p name of program, usually this one
  /// \note implemented in nemo++.cc using NEMO's mysymbols()
  void mysymbols(const char*p);

  /// load all local symbols
  inline void localsymbols() {
    mysymbols(getparam("argv0"));
  }

  /// load all symbols from an .so file
  /// \param[in] f name of .so file
  /// \note implemented in nemo++.cc using NEMO's loadobj()
  void loadobj(const char*f);

  /// find function with C-linkage
  /// \param[in] f name of function with C-linkage to find
  /// \return pointer to function
  /// \note implemented in nemo++.cc using NEMO's findfn()
  void(*findfn(const char*f))();

  /// translates a generic symbol name.
  /// this routine translates a generic symbol name into the one which can be
  /// found in the symbol table of the host loader
  /// \param[in,out] f name of symbol with C-linkage
  /// \note implemented in nemo++.cc using NEMO's mapsys()
  void mapsys(char*f);

  //@}
} // namespace falcON

namespace falcON {
  // ///////////////////////////////////////////////////////////////////////////
  //
  // support for nemo file I/O
  //
  // ///////////////////////////////////////////////////////////////////////////
  /// base class for \a nemo_in and \a nemo_out                                 
  class nemo_io {
  public:
    /// used to identify different nemo I/O able data types                     
    enum DataType {
      Null    = 0,
      Byte    = 1,
      Short   = 2,
      Integer = 3,
      Long    = 4,
      Single  = 5,
      Double  = 6,
#ifdef falcON_REAL_IS_FLOAT
      Real    = Single,
      NotReal = Double
#else
      Real    = Double,
      NotReal = Single
#endif
    };
    /// returns the name of a Datatype as string recognised by NEMO methods
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
    /// two DataTypes are coercing if they are of different floating point type
    static bool coercing(DataType t1, DataType t2) {
      return 
	(t1 == Double && t2 == Single) || (t2 == Double && t1 == Single);
    }
    /// different Datatypes are mismatched, unless they are coercing
    static bool mismatch(DataType t1, DataType t2) {
      return t1!=t2 && !coercing(t1,t2);
    }
    /// identifiers for the various data fields supported by nemo_io
    enum Field {
      null    = 0,
      mass    = 1 <<  0,
      pos     = 1 <<  1,
      vel     = 1 <<  2,
      eps     = 1 <<  3,
      key     = 1 <<  4,
      step    = 1 <<  5,
      pot     = 1 <<  6,
      acc     = 1 <<  7,
      jerk    = 1 <<  8,
      dens    = 1 <<  9,
      aux     = 1 << 10,
      zet     = 1 << 11,
      lev     = 1 << 12,
      num     = 1 << 13,
      posvel  = 1 << 14,
      phden   = 1 << 15,
      torb    = 1 << 16,
      Size    = 1 << 17,
      Gasnum  = 1 << 18,
      Uin     = 1 << 19,
      Uindot  = 1 << 20,
      Uinrad  = 1 << 21,
      Entr    = 1 << 22,
      Gasdens = 1 << 23,
      Sizedot = 1 << 24,
      Sphfact = 1 << 25,
      Csound  = 1 << 26,
      AlphaAV = 1 << 27,
      DivV    = 1 << 28,
      MolWght = 1 << 29,
      Spin    = 1 << 30
    };
    /// get the nemo_io::DataType used for a given Field
    static DataType type(Field f) {
      switch(f) {
      case mass   : return Real;
      case pos    : return Real;
      case vel    : return Real;
      case eps    : return Real;
      case key    : return Integer;
      case step   : return Real;
      case pot    : return Real;
      case acc    : return Real;
      case jerk   : return Real;
      case dens   : return Real;
      case aux    : return Real;
      case zet    : return Real;
      case lev    : return Short;
      case num    : return Integer;
      case posvel : return Real;
      case phden  : return Real;
      case torb   : return Real;
      case Size   : return Real;
      case Gasnum : return Integer;
      case Uin    : return Real;
      case Uindot : return Real;
      case Uinrad : return Real;
      case Entr   : return Real;
      case Gasdens: return Real;
      case Sizedot: return Real;
      case Sphfact: return Real;
      case Csound : return Real;
      case AlphaAV: return Real;
      case DivV   : return Real;
      case MolWght: return Real;
      case Spin   : return Real;
      default     : return Null;
      }
    }
    /// what \a Field corresponds to a given \a fieldbit
    static Field field(fieldbit f) {
      switch(value(f)) {
      case fieldbit::m: return mass;
      case fieldbit::x: return pos;
      case fieldbit::v: return vel;
      case fieldbit::e: return eps;
      case fieldbit::k: return key;
      case fieldbit::s: return step;
      case fieldbit::p: return pot;
      case fieldbit::q: return pot;
      case fieldbit::a: return acc;
      case fieldbit::j: return jerk;
      case fieldbit::r: return dens;
      case fieldbit::y: return aux;
      case fieldbit::z: return zet;
      case fieldbit::l: return lev;
      case fieldbit::n: return num;
      case fieldbit::d: return phden;
      case fieldbit::t: return torb;
      case fieldbit::H: return Size;
      case fieldbit::N: return Gasnum;
      case fieldbit::U: return Uin;
      case fieldbit::I: return Uindot;
      case fieldbit::E: return Uinrad;
      case fieldbit::K: return Entr;
      case fieldbit::R: return Gasdens;
      case fieldbit::A: return AlphaAV;
      case fieldbit::D: return DivV;
      case fieldbit::J: return Sizedot;
      case fieldbit::F: return Sphfact;
      case fieldbit::C: return Csound;
      case fieldbit::M: return MolWght;
      case fieldbit::S: return Spin;
      default         : return null;
      }
    }
    /// what \a fieldbit corresponds to a given \a Field
    static fieldbit bit(Field F) {
      switch(F) {
      case mass:    return fieldbit::m;
      case pos:     return fieldbit::x;
      case vel:     return fieldbit::v;
      case eps:     return fieldbit::e;
      case key:     return fieldbit::k;
      case step:    return fieldbit::s;
      case pot:     return fieldbit::p;
      case acc:     return fieldbit::a;
      case jerk:    return fieldbit::j;
      case dens:    return fieldbit::r;
      case aux:     return fieldbit::y;
      case zet:     return fieldbit::z;
      case lev:     return fieldbit::l;
      case num:     return fieldbit::n;
      case phden:   return fieldbit::d;
      case torb:    return fieldbit::t;
      case Size:    return fieldbit::H;
      case Gasnum:  return fieldbit::N;
      case Uin:     return fieldbit::U;
      case Uindot:  return fieldbit::I;
      case Uinrad:  return fieldbit::E;
      case Entr:    return fieldbit::K;
      case Gasdens: return fieldbit::R;
      case Sizedot: return fieldbit::J;
      case Sphfact: return fieldbit::F;
      case Csound:  return fieldbit::C;
      case AlphaAV: return fieldbit::A;
      case DivV:    return fieldbit::D;
      case MolWght: return fieldbit::M;
      case Spin:    return fieldbit::S;
      default: falcON_Warning("unaccountable nemo_io::Field\n");
	return fieldbit::invalid;
      }
    }
    /// given a collection of \a Field s, return the corresponding fieldset     
    static fieldset fields(nemo_io::Field F)
    {
      fieldset data;
      for(fieldbit f; f; ++f)
	if(field(f) & F) data |= fieldset(f);
      return data;
    }
    //--------------------------------------------------------------------------
  protected:
    std::FILE *STREAM;              ///< file handle
    bool IN, OUT, PIPE, SINK;       ///< input?, output?, pipe?, sink?
    /// close stream if open, then open new file
    /// \param[in] file  name of file to open, "-": stdin/stdout, "." sink
    /// \param[in] mode  r,w,w!,s for read, write, overwrite, scratch
    void open(const char*file, const char*mode);
    /// close stream if open
    void close();
    /// default ctor: set stream to (nil)
    nemo_io() : STREAM(0), IN(0), OUT(0), PIPE(0), SINK(0) {}
    /// dtor: close stream
    ~nemo_io() { close(); }
  public:
    /// are we a pipe?
    bool const&is_pipe() const { return PIPE; }
    /// are we a sink?
    bool const&is_sink() const { return SINK; }
    /// are we ready for I/O?
    bool is_open() const { return STREAM!=0; }
    /// are we ready for input?
    bool is_reading() const { return IN; }
    /// are we ready for ouput?
    bool is_writing() const { return OUT && !SINK; }
  };// class nemo_io
  class snap_in;
  // ///////////////////////////////////////////////////////////////////////////
  /// represents a NEMO input stream, derived from \a nemo_io
  class nemo_in : public nemo_io {
    friend class snap_in;
    //--------------------------------------------------------------------------
  private:
    mutable snap_in *SNAP;                  ///< open snapshot, if any
    nemo_in(nemo_in const&);                // not implemented
    //--------------------------------------------------------------------------
  public:
    /// close open input stream (and close an open snap_in, if any)
    void close();
    /// open new file (and close any old input stream)
    /// \param[in] file name of file to be opened, may be (nil)
    nemo_in&open(const char*file=0) {
      close();
      nemo_io::open(file,"r");
      return *this;
    }
    /// default ctor and ctor from file name: open file for NEMO input.
    /// If \a file equals "-", set pipe, ie. read from stdin instead from a file
    /// \param[in] file name of file to be opened, may be (nil)
    explicit nemo_in(const char*file=0) : SNAP(0) {
      nemo_io::open(file,"r");
    }
    /// dtor: close()
    ~nemo_in() { close(); }
    /// can a snapshot be opened?
    /// \return true if a \a snap_in can be constructed from *this
    bool has_snapshot() const;
    /// conversion to bool: are we ready for input?
    operator bool() const { return is_reading(); }
  };// class nemo_in
  class data_in;
  // ///////////////////////////////////////////////////////////////////////////
  /// represents a snapshot on a NEMO input stream, \a nemo_in
  ///
  /// At any time, only one \a snap_in can exist for any \a nemo_in.
  /// Construction of a second \a snap_in from the same \a nemo_in will cause
  /// a fatal error. A \a snap_in provides information on which data are
  /// contained in the snapshot and allows to construct \a data_in for reading
  /// them.
  class snap_in {
    friend class data_in;
    //--------------------------------------------------------------------------
  private:
    nemo_in const   &INPUT;                        // our input          
    mutable data_in *DATA;                         // if non-zero: open snapshot
    mutable int      FIELDS_READ;                  // fields read already       
    bool             HAS_TIME;                     // have simulation time?     
    unsigned         NTOT, NBOD[bodytype::NUM];    // # bodies, # bodies / type 
    double           TIME;                         // simulations time          
    //--------------------------------------------------------------------------
    std::FILE* stream() const { return INPUT.STREAM; }
    //--------------------------------------------------------------------------
  public:
    /// ctor: open a snapshot set from a NEMO input stream
    /// \param[in] in NEMO input stream
    explicit snap_in(nemo_in const&in);
    //--------------------------------------------------------------------------
    /// dtor: close snapshot set in NEMO input stream
    ~snap_in();
    //--------------------------------------------------------------------------
    /// do we have simul time?  
    bool has_time() const { return HAS_TIME; }
    /// return total # bodies
    unsigned const&Ntot() const { return NTOT; }
    /// return array with # bodies per type
    const unsigned*Nbod() const { return NBOD; }
    /// return # bodies per given type
    unsigned const&Nbod(bodytype t) const { return NBOD[t]; }
    /// return the number of data expected for a given field
    unsigned N(nemo_io::Field F) const {
      unsigned n(0u);
      fieldbit b= F==nemo_io::posvel? fieldbit(fieldbit::x) : nemo_io::bit(F);
      for(bodytype t; t; ++t)
	if(t.allows(b)) n += NBOD[t];
      return n;
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
  };// class snap_in
  // ///////////////////////////////////////////////////////////////////////////
  /// represents the data of one field from one snapshot of a NEMO input stream
  class data_in {
    //--------------------------------------------------------------------------
  private:
    snap_in const          &INPUT;                 // our snapshot              
    const nemo_io::Field    FIELD;                 // sort of data we are for   
    unsigned                NREAD;                 // how many have been read?  
    unsigned                NTOT;                  // # data actually on file   
    nemo_io::DataType       TYPE;                  // which data type?          
    unsigned                SUBN;                  // how many items per datum  
    static const int        NDIM = Ndim;           // # spatial dimensions      
  public:
    /// ctor: create a valid data input with the expected number of bodies and
    ///       no type mismatch; otherwise throw an exception
    /// \param[in] s  input snapshot to be read from
    /// \param[in] f data field to be read
    data_in(snap_in const &s, nemo_io::Field f);
    /// dtor: close data set
    ~data_in();
    /// which data field are we reading?
    nemo_io::Field     const&field()      const { return FIELD; }
    /// which data type are we reading?
    nemo_io::DataType  const&type()       const { return TYPE; }
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
    /// read data into array given
    /// \param[in] data  pointer to data array, must have sufficient memory
    /// \param[in] n     number of bodies to read data for; default: all left
    /// \note We coerce data on the fly (convert floating point to real)
    void read(void*data, unsigned n=0);
    /// read position and velocity from phases
    /// \note must only be called when field() == nemo_io::posvel
    /// \param[in] pos   array for position data, only read if non-null
    /// \param[in] vel   array for velocity data, only read if non-null
    /// \param[in] n     number of bodies to read data for; default: all left
    /// \note We coerce data on the fly (convert floating point to real)
    void read_phases(void*pos, void*vel, unsigned n=0);
  };// class data_in
  class snap_out;
  // ///////////////////////////////////////////////////////////////////////////
  /// represents a NEMO output stream, derived from \a nemo_io
  class nemo_out : public nemo_io {
    friend class snap_out;
    //--------------------------------------------------------------------------
  private:
    mutable snap_out *SNAP;                    // if non-zero open snapshot 
    bool              IS_PIPE, IS_SINK;            // are we a pipe or sink?    
    nemo_out(nemo_out const&);                     // not implemented           
  public:
    /// close open output stream
    void close();
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
    /// existing file of the same name (except for the trailing '@') will be
    /// appended to. If no file of that name exists, we proceed as usual and
    /// open a file of name \e file (omitting any trailing '!' or '@').
    ///
    /// \return *this
    /// \param[in] file name of file to be opened, see detailled comments
    /// \param[in] append (optional) append existing file anyway?
    nemo_out&open(const char*file=0, bool append=false);
    /// default ctor and ctor from file name: open file for NEMO output.
    /// If \a file equals "-", set pipe, ie. write to stdout instead to a file
    /// \param[in] file   name of file to be opened, may be (nil)
    /// \note             see detailled comments for open()
    /// \param[in] append (optional) append existing file anyway?
    explicit nemo_out(const char*file=0, bool append=false)
      : SNAP(0) { open(file,append); }
    /// close open output stream
    ~nemo_out() { close(); }
    /// convertion to bool: true if is_writing()
    operator bool() const { return is_writing(); }
  };// class nemo_out
  class data_out;
  // ///////////////////////////////////////////////////////////////////////////
  /// represents a snapshot on a NEMO output stream, \a nemo_out
  ///
  /// At any time, only one \a snap_out can exist for any \a nemo_out.
  /// Construction of a second \a snap_out from the sane \a nemo_out will
  /// cause a fatal error. A \a snap_out allows to construct \a data_out for
  /// writing data.
  class snap_out {
    friend class data_out;
  private:
    nemo_out const   &OUTPUT;                      // our output stream         
    mutable data_out *DATA;                        // if non-zero: open data_out
    mutable int       FIELDS_WRITTEN;              // data already written out  
    unsigned          NTOT, NBOD[bodytype::NUM];   // # bodies, # bodies / type 
    //
    std::FILE* stream() const { return OUTPUT.STREAM; }
  public:
    /// ctor: open NEMO snapshot set
    /// \param[in] out   NEMO output stream
    /// \param[in] Nbod  total # bodies per bodytype
    /// \param[in] time  simulation time
    snap_out(nemo_out const&out, const unsigned Nbod[bodytype::NUM],
	     double time);
    /// dtor: close snapshot set in NEMO output stream
    ~snap_out();
    /// return total # bodies
    unsigned const&Ntot() const { return NTOT; }
    /// return # bodies per type
    unsigned const&Nbod(bodytype t) const { return NBOD[t]; }
    /// return the number of data expected for a given field
    unsigned N(nemo_io::Field F) const {
      unsigned n(0u);
      fieldbit b= F==nemo_io::posvel? fieldbit(fieldbit::x) : nemo_io::bit(F);
      for(bodytype t; t; ++t)
	if(t.allows(b)) n += NBOD[t];
      return n;
    }
    /// have data for given field been written already?
    bool has_been_written(nemo_io::Field f) const {
      return FIELDS_WRITTEN & f;
    }
    /// for which fields have data been written yet?
    int const&fields_written() const {
      return FIELDS_WRITTEN;
    }
  };// class snap_out
  // ///////////////////////////////////////////////////////////////////////////
  /// represents the data of one field from one snapshot of a NEMO output stream
  class data_out {
  private:
    snap_out const        &OUTPUT;                 // our snapshot              
    const nemo_io::Field   FIELD;                  // sort of data we are for   
    unsigned               NWRITTEN;               // how many have been written
    const unsigned         NTOT;                   // how many in total         
    nemo_io::DataType      TYPE;                   // which data type?          
    unsigned               SUBN;                   // how many items per datum  
    static const int       NDIM = Ndim;            // # spatial dimensions      
  public:
    /// ctor: open a NEMO data set
    /// \param[in] s  our \a snap_out snapshot set to write to
    /// \param[in] f  the data field to write out with this
    data_out(snap_out const&s, nemo_io::Field f);
    /// dtor: close data set
    ~data_out();
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
    /// \param[in] data  pointer to data to be written
    void write(const void*data);
    /// write data for \e n bodies from array given (must hold enough data)
    ///
    /// \param[in] data  pointer to data to be written
    /// \param[in] n     # bodies for which to write data
    void write(const void*data, unsigned n);
  };
  //////////////////////////////////////////////////////////////////////////////
  /// is a given time in a given time range?
  /// \param[in] t (input) simulation time
  /// \param[in] r (input) time range as string, using NEMO notation
  bool time_in_range(double t, const char*r);
} // namespace falcON
#ifdef falcON_NEMO
falcON_TRAITS(falcON::nemo_io,"nemo_io");
falcON_TRAITS(falcON::nemo_in,"nemo_in");
falcON_TRAITS(falcON::snap_in,"snap_in");
falcON_TRAITS(falcON::data_in,"data_in");
falcON_TRAITS(falcON::nemo_out,"nemo_out");
falcON_TRAITS(falcON::snap_out,"snap_out");
falcON_TRAITS(falcON::data_out,"data_out");
#endif
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_NEMO
#endif // falcON_included_nemopp_h
