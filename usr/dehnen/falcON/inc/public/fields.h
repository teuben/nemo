// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file    inc/public/fields.h
///
/// \brief   body flags and support for body data management
///
///          Contents:
///          \li falcON::flags: holding boolean information about each body
///          \li falcON::fieldbit: representing a single body data field
///          \li falcON::fieldset: representing a set of body data fields
///          \li template falcON::field_traits 
///          \li further constants and macros
///
/// \author  Walter Dehnen
///#
/// \date    2002-2008
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2002-2008  Walter Dehnen
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
////////////////////////////////////////////////////////////////////////////////
#ifndef falcON_included_fields_h
#define falcON_included_fields_h 1

#ifndef falcON_included_cstring
#  include <cstring>
#  define falcON_included_cstring
#endif
#ifndef falcON_included_basic_h
#  include <public/basic.h>
#endif

////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  // ///////////////////////////////////////////////////////////////////////////
  //
  // class falcON::flags
  //
  /// An integer, the bits of which are used as flags for body properties;
  /// each body, as well as tree cells and leafs have a flag.
  ///
  // ///////////////////////////////////////////////////////////////////////////
  class flags {
  public:
    /// enumeration holding the individual flags
    enum single {
      empty         = 0,        ///< no flag set
      active        = 1 <<  0,  ///< a body/tree node is active
      remove        = 1 <<  1,  ///< body to be removed, see bodies::remove()
      sph           = 1 <<  2,  ///< body is SPH particle
      sink          = 1 <<  3,  ///< body is a sink particle
      newbody       = 1 <<  4,  ///< body is new
      adjust_step   = 1 <<  5,  ///< help with time integration
      not_longer    = 1 <<  6,  ///< help with time integration
      not_shorter   = 1 <<  7,  ///< help with time integration
      sph_special   = 1 <<  8,  ///< help with time integration for SPH
      interacting   = 1 <<  9,  ///< help with SPH estimation
      subtree       = 1 << 10,  ///< help with sub-tree building
      marked        = 1 << 11,  ///< auxiliary for several purposes
      sticky        = 1 << 12,  ///< body is a sticky particle
      ignore        = 1 << 15,  ///< flag used to define subsets of bodies
      // flags for cells only
      all_active    = 1 << 20,  ///< all leafs in cell are active
      all_sph       = 1 << 21,  ///< all leafs in cell are SPH particles
      all_sink      = 1 << 21,  ///< all leafs in cell are SPH particles
      all_sticky    = 1 << 23,  ///< all leafs in cell are sticky particles
      subtree_cell  = 1 << 24
    };
    /// enumeration holding some combinations of flags
    enum combined {
      leaf_flags    = active|sph|sink|sticky|sph_special|interacting,
      body_flags    = active|sph|sink|sticky|sph_special|remove,
      subtree_flags = subtree|subtree_cell,
      activ_flags   = active|all_active
    };
    //--------------------------------------------------------------------------
    /// data: 32 bits -> up to 32 flags can be supported
    typedef int value_type;
  private:
    value_type val;
    //--------------------------------------------------------------------------
    /// \name constructors and assignment
    //@{
    /// construction from integer: set equal
    explicit flags(int i) : val(i) {}
  public:
    /// unitialized construction: all flags are unset
    flags() : val(empty) {}
    /// copy constructor
    flags(flags const&f) : val(f.val) {}
    /// construction from single flag
    flags(single f) : val(f) {}
    /// construction from combined flag
    flags(combined f) : val(f) {}
    /// assign from flags
    flags& operator=(flags const&f) {
      val = f.val;
      return *this;
    }
    /// assign from single flag
    flags& operator=(single f) {
      val = f;
      return *this;
    }
    /// assign from combined flag
    flags& operator=(combined f) {
      val = f;
      return *this;
    }
    //@}
    //--------------------------------------------------------------------------
    /// \name operations
    //@{
    /// add another set of flags
    flags& operator|= (flags const&f) {
      val |= f.val;
      return *this;
    }
    /// add a specific flag
    flags& operator|= (single b) {
      val |= b;
      return *this;
    }
    /// add a combined flag
    flags& operator|= (combined b) {
      val |= b;
      return *this;
    }
    /// combine sets of flags
    flags operator| (flags const&f) const {
      return flags( val | f.val );
    }
    /// combine set of flags with a single flag
    flags operator| (single b) const {
      return flags( val | b );
    }
    /// combine set of flags with a combined flag
    flags operator| (combined b) const {
      return flags( val | b );
    }
    /// combine single flag with set of flags
    friend flags operator| (single b, flags const&f);
    /// combine single flag with set of flags
    friend flags operator| (combined b, flags const&f);
    /// combine single flags
    friend flags operator| (single b, single c);
    /// combine combined with single flags
    friend flags operator| (combined b, single c);
    /// combine single with combined flags
    friend flags operator| (single b, combined c);
    /// combine combined flags
    friend flags operator| (combined b, combined c);
    //@}
    //--------------------------------------------------------------------------
    /// \name general boolean methods
    //@{
    bool operator== (flags const&f) const { return val == f.val; }
    bool operator!= (flags const&f) const { return val != f.val; }
    bool are_set    (flags const&f) const { return val & f.val; }
    bool is_set     (single        b) const { return val & b; }
    //@}
    //--------------------------------------------------------------------------
    /// \name specific boolean methods taking reference to flags
    //@{
    friend bool is_active  (flags const&f);
    friend bool to_remove  (flags const&f);
    friend bool is_sph     (flags const&f);
    friend bool is_sink    (flags const&f);
    friend bool is_sticky  (flags const&f);
    friend bool is_new     (flags const&f);
    friend bool in_subtree (flags const&f);
    friend bool is_marked  (flags const&f);
    friend bool al_active  (flags const&f);
    friend bool al_sph     (flags const&f);
    friend bool al_sink    (flags const&f);
    friend bool al_sticky  (flags const&f);
    friend bool in_subset  (flags const&f);
    //@}
    //--------------------------------------------------------------------------
    /// \name specific boolean methods taking pointer to flags
    //@{
    friend bool is_active  (const flags*f);
    friend bool to_remove  (const flags*f);
    friend bool is_sph     (const flags*f);
    friend bool is_sink    (const flags*f);
    friend bool is_sticky  (const flags*f);
    friend bool is_new     (const flags*f);
    friend bool in_subtree (const flags*f);
    friend bool is_marked  (const flags*f);
    friend bool al_active  (const flags*f);
    friend bool al_sph     (const flags*f);
    friend bool al_sink    (const flags*f);
    friend bool al_sticky  (const flags*f);
    friend bool in_subset  (const flags*f);
    //@}
    //--------------------------------------------------------------------------
    /// \name non-const methods                                                 
    //@{
    /// reset all bits to zero
    void reset() { val = 0; }
    /// set equal to part of another set of flags
    void set_to_part(flags const&f, int mask) { val = f.val & mask;}
    /// add a single flag
    void add(single f) { (*this) |= f; }
    /// add another set of flags
    void add(flags const&f) { (*this) |= f; }
    /// add part of another set of flags
    void add_part(flags const&f, int mask) { val |= f.val & mask; }
    /// add part of another set of flags
    void add_part(flags const&f, single mask) { val |= f.val & mask; }
    /// set part equal to corresponding part of another set of flags
    void set_part(flags const&f, int mask) { 
      val &= ~mask;
      val |= f.val & mask;
    }
    /// reset a single flag
    void un_set(single f) { val &= ~ int(f); }
    /// reset all flags in a given set of flags
    void un_set(flags const&f) { val &= ~(f.val); }
    //@}
    //--------------------------------------------------------------------------
    /// type conversion to bool
    operator bool () const { return val != 0; }
    //--------------------------------------------------------------------------
    /// \name formatted I/O                                                     
    //@{
    /// formatted output: simply write integer value
    friend std::ostream& operator<< (std::ostream&, flags const&);
    /// formatted input: read as integer
    friend std::istream& operator>> (std::istream&, flags&);
    //@}
  };// class flags
  // ///////////////////////////////////////////////////////////////////////////
  //
  // inline definitions of friends of class flags
  // also serve to inject these functions into namespace falcON
  //
  // ///////////////////////////////////////////////////////////////////////////
  inline flags operator| (flags::single b, flags const&f) {
    return flags( f.val | int(b) );
  }
  inline flags operator| (flags::combined b, flags const&f) {
    return flags( f.val | int(b) );
  }
  inline flags operator| (flags::single b, flags::single c) {
    return flags( int(b) | int(c) );
  }
  inline flags operator| (flags::combined b, flags::single c) {
    return flags( int(b) | int(c) );
  }
  inline flags operator| (flags::single b, flags::combined c) {
    return flags( int(b) | int(c) );
  }
  inline flags operator| (flags::combined b, flags::combined c) {
    return flags( int(b) | int(c) );
  }
  inline bool is_active (flags const&f) { return f.is_set(flags::active); }
  inline bool to_remove (flags const&f) { return f.is_set(flags::remove); }
  inline bool is_sph    (flags const&f) { return f.is_set(flags::sph); }
  inline bool is_sink   (flags const&f) { return f.is_set(flags::sink); }
  inline bool is_sticky (flags const&f) { return f.is_set(flags::sticky); }
  inline bool is_new    (flags const&f) { return f.is_set(flags::newbody); }
  inline bool in_subtree(flags const&f) { return f.is_set(flags::subtree); }
  inline bool is_marked (flags const&f) { return f.is_set(flags::marked); }
  inline bool al_active (flags const&f) { return f.is_set(flags::all_active); }
  inline bool al_sph    (flags const&f) { return f.is_set(flags::all_sph); }
  inline bool al_sink   (flags const&f) { return f.is_set(flags::all_sink); }
  inline bool al_sticky (flags const&f) { return f.is_set(flags::all_sticky); }
  inline bool ignore    (flags const&f) { return f.is_set(flags::ignore); }
  inline bool in_subset (flags const&f) { return !f.is_set(flags::ignore); }
  inline bool is_active (const flags*f) { return f->is_set(flags::active); }
  inline bool to_remove (const flags*f) { return f->is_set(flags::remove); }
  inline bool is_sph    (const flags*f) { return f->is_set(flags::sph); }
  inline bool is_sink   (const flags*f) { return f->is_set(flags::sink); }
  inline bool is_sticky (const flags*f) { return f->is_set(flags::sticky); }
  inline bool is_new    (const flags*f) { return f->is_set(flags::newbody); }
  inline bool in_subtree(const flags*f) { return f->is_set(flags::subtree); }
  inline bool is_marked (const flags*f) { return f->is_set(flags::marked); }
  inline bool al_active (const flags*f) { return f->is_set(flags::all_active); }
  inline bool al_sph    (const flags*f) { return f->is_set(flags::all_sph); }
  inline bool al_sink   (const flags*f) { return f->is_set(flags::all_sink); }
  inline bool al_sticky (const flags*f) { return f->is_set(flags::all_sticky); }
  inline bool ignore    (const flags*f) { return f->is_set(flags::ignore); }
  inline bool in_subset (const flags*f) { return !f->is_set(flags::ignore); }
  inline std::ostream& operator<< (std::ostream& o, flags const&f) {
    return o << f.val;
  }
  inline std::istream& operator>> (std::istream& i, flags&f) {
    return i >> f.val;
  }
  //////////////////////////////////////////////////////////////////////////////
  //
  // N-body data type management
  //
  // There are many different N-body data (mass, position, etc), each with its
  // unique data type and I/O requirements etc. Below are several types and
  // constants designed for managing this.
  //
  // class fieldbit
  // represents individual body data types as integer, i.e. 0 corresponds to
  // mass, 1 to positions, etc. fieldbit can be used as a java-style iterator
  // over allowed data types.
  //
  // class fieldset
  // represents a collection of body data types as integer each bit
  // corresponding to the data whose fieldbit matches its position. fieldset
  // are used to indicate which body data should be, e.g., allocated, read,
  // written, or copied, etc.
  //
  //////////////////////////////////////////////////////////////////////////////
  //
  // some constants related to N-body data
  //
  //////////////////////////////////////////////////////////////////////////////
  namespace BodyData {
    const int         KSTD    = 0;                 // std data start here       
    const int         KSPH    = 22;                // sph data start here       
    const int         KSINK   = 36;                // sink data start here      
    const int         NQUANT  = 37;                // total # data              
    //--------------------------------------------------------------------------
    /// array with the one-char data tags used as enum names
    /// \note not (to be) referred to outside of this file
    const char* const SQUANT  ="mxvuefkspqajriyzlnchdtHNUYIEKRADJFCMS";
    // letters taken are marked here
    //             abcdefghijklmnopqrstuvwxyz
    // lower case  x xxxx xxxxxxxxxxxxxxx xxx
    // upper case  x xxxx xxxx xx   xx x   x 
    //--------------------------------------------------------------------------
    /// array with full-length human readable names for the N-body data
    /// \note not (to be) referred to outside of this file
    const char* const QFULLNAME[NQUANT] = 
      { "mass", "position", "velocity", "predicted velocity",
	"softening length", "body flags", "key", "time step", "potential",
	"external potential", "acceleration", "jerk", "density",
	"auxiliary integer", "auxiliary scalar", "auxiliary vector",
	"time-step level", "number of partners", "node number",
	"Peano-Hilbert key", "phase-space density", "orbital period",
	"SPH smoothing length", "number of SPH partners", "U_internal",
	"U_predicted", "(dU/dt)_total", "(dU/dt)_external", "entropy function",
	"gas density", "alpha_visc", "div(v)", "dlog(h)/dt", "factor",
	"sound speed", "molecular weight", "spin"
      };
    /// array with function names for the N-body data
    /// \note not (to be) referred to outside of this file
    const char* const QFUNCNAME[NQUANT] = 
      { "mass", "pos", "vel", "vprd", "eps", "flag", "key", "step", "pot",
	"pex", "acc", "jerk", "rho", "integ", "aux", "zet", "level", "num",
	"node", "peano", "phden", "torb",
	"size", "snum", "uin", "uprd", "udot", "udex", "entr", "srho", "alfa",
	"divv", "hdot", "fact", "csnd", "molw", "spin"
      };
    /// array with five-chacter names for the N-body data
    /// \note not (to be) referred to outside of this file
    const char* const QFIVENAME[NQUANT] = 
      { "mass ", "pos  ", "vel  ", "vprd ", "eps  ", "flag ", "key  ", "step ",
	"pot  ", "pex  ", "acc  ", "jerk ", "rho  ", "integ", "aux  ", "zet  ",
	"level", "num  ", "node ", "peano", "phden", "torb ",
	"size ", "snum ", "uin  ", "uprd ", "udot ", "udex ", "entr ", "srho ",
	"alfa ", "divv ", "hdot ", "fact ", "csnd ", "molw ", "spin "
      };
    /// array with the sizeof() the N-body data
    const size_t ZQUANT[NQUANT] = {
      //            source properties: 8
      sizeof(real),           ///< mass
      sizeof(vect),           ///< position
      sizeof(vect),           ///< velocity
      sizeof(vect),           ///< velocity_predicted
      sizeof(real),           ///< softening length
      sizeof(flags),          ///< flags
      sizeof(unsigned),       ///< key
      sizeof(double),         ///< time step
      //            non-source properties: 12
      sizeof(real),           ///< potential
      sizeof(real),           ///< external potential
      sizeof(vect),           ///< acceleration
      sizeof(vect),           ///< jerk
      sizeof(real),           ///< density (mass)
      sizeof(int),            ///< auxiliary integer
      sizeof(real),           ///< auxiliary scalar
      sizeof(vect),           ///< auxiliary vector
      sizeof(indx),           ///< level
      sizeof(unsigned),       ///< number of partners
      sizeof(indx),           ///< node number
      sizeof(peanokey),       ///< Peano-Hilbert key
      sizeof(real),           ///< phase-space density
      sizeof(real),           ///< orbital period
      //            SPH properties: 13
      sizeof(real),           ///< size
      sizeof(unsigned),       ///< number of SPH partners
      sizeof(real),           ///< U_internal
      sizeof(real),           ///< U_predicted
      sizeof(real),           ///< (dU/dt)_total
      sizeof(real),           ///< (dU/dt)_external
      sizeof(real),           ///< entropy function
      sizeof(real),           ///< gas-density
      sizeof(real),           ///< alpha_visc
      sizeof(real),           ///< div(v)
      sizeof(real),           ///< dlogh/dt
      sizeof(real),           ///< f_i
      sizeof(real),           ///< sound speed
      sizeof(real),           ///< sound speed
      //            SINK properties: 1
      sizeof(vect)            ///< spin
    };
  } // namespace BodyData
  // ///////////////////////////////////////////////////////////////////////////
  //
  // class falcON::fieldbit
  //
  /// Represents a single individual body data field
  ///
  /// Each body data field corresponds to another integer: 0 means mass, 1
  /// positions, etc. In any application, \b never use these integer values,
  /// but \b always the enumeration value in enum fieldbit::bits, because
  /// these are guaranteed to be preserved during code changes. fieldbit can
  /// be used as a java-style iterator over data fields.
  ///
  // ///////////////////////////////////////////////////////////////////////////
  class fieldbit {
  public:
    //--------------------------------------------------------------------------
    /// type of internal representation
    typedef int value_type;
    /// \name static members
    //@{
    /// total # data
    static const int   NQUANT= BodyData::NQUANT;
    /// letter for body datum
    static const char& SNAME(int i) { return BodyData::SQUANT[i]; }
    /// 5-character function name for body dataum
    static const char* FIVENAME(int i) { return BodyData::QFIVENAME[i]; }
    /// function name for body dataum
    static const char* FUNCNAME(int i) { return BodyData::QFUNCNAME[i]; }
    /// full name for body dataum
    static const char* FULLNAME(int i) { return BodyData::QFULLNAME[i]; }
    /// size of type for body datum
    static size_t const& SIZE(int i) { return BodyData::ZQUANT[i]; }
    //@}
    //--------------------------------------------------------------------------
  private:
    value_type val;
  public:
    //--------------------------------------------------------------------------
    /// specifies the field - value relation
    ///
    /// value names are single character, lower case corresponds to general
    /// body data, upper case to SPH body data.
    enum bits { 
      //            source properties: 8
      m       = 0,            ///< mass
      x       = 1,            ///< position
      v       = 2,            ///< velocity
      u       = 3,            ///< predicted velocity
      V       = u,            ///< predicted velocity (again)
      e       = 4,            ///< individual gravitational softening length
      f       = 5,            ///< falcON::flags
      k       = 6,            ///< integer key or identifier
      s       = 7,            ///< time step (not used in falcON)
      //            non-source properties: 10
      p       = 8,            ///< N-body gravitational potential
      q       = 9,            ///< external gravitational potential
      a       = 10,           ///< accelaration (gravity and SPH)
      j       = 11,           ///< jerk (not used in falcON)
      r       = 12,           ///< mass-density
      i       = 13,           ///< auxiliary integer datum
      y       = 14,           ///< auxiliary scalar datum
      z       = 15,           ///< auxiliary vector datum
      l       = 16,           ///< time-step level (short)
      n       = 17,           ///< number of neighbours
      c       = 18,           ///< node number (short unsigned)
      h       = 19,           ///< Peano-Hilbert key
      d       = 20,           ///< phase-space density
      t       = 21,           ///< orbital period
      //            SPH properties: 12
      H       = 22,           ///< SPH smoothing length
      N       = 23,           ///< number of SPH interaction partners
      U       = 24,           ///< internal energy U
      Y       = 25,           ///< predicted U
      I       = 26,           ///< (dU/dt)_total
      E       = 27,           ///< (dU/dt)_external
      K       = 28,           ///< entropy function
      R       = 29,           ///< gas-density rho
      A       = 30,           ///< alpha_visc
      D       = 31,           ///< div(v)
      J       = 32,           ///< dlnh/dt
      F       = 33,           ///< factor f_i
      C       = 34,           ///< sound speed
      M       = 35,           ///< molecular weight
      S       = 36,           ///< spin vector
      invalid = NQUANT        ///< not corresponding to any field
    };
    //--------------------------------------------------------------------------
    /// \name construction and assignment
    //@{
    /// unitialized: mass
    fieldbit()                    : val(0) {}
    /// copy constructor
    fieldbit(fieldbit const   &b) : val(b.val) {}
    /// from fieldbit::bits
    fieldbit(bits              b) : val(b) {}
    /// from integer: use value
    explicit fieldbit(int      i) : val(i) {}
    /// from unsigned: use value
    explicit fieldbit(unsigned i) : val(i) {}
    /// from char: match against \b name of enum values in fieldbit::bits
    explicit fieldbit(char     l) : val(invalid) {
      for(int v=0; v!=invalid; ++v)
	if(l == BodyData::SQUANT[v] ) {
	  val = v;
	  return;
	}
    }
    /// copy operator
    fieldbit &operator= (fieldbit const&b) { val=b.val; return *this; }
    //@}
    //--------------------------------------------------------------------------
    /// \name support for looping over data fields
    //@{
    /// prefix increment: get next field
    fieldbit &operator++()                 { ++val; return *this; }
    /// conversion to bool: are we still a valid field
    operator bool       ()           const { return val<invalid; }
    /// conversion to value_type
    operator value_type ()           const { return val; }
    /// equality
    bool operator ==    (fieldbit f) const { return val == f.val; }
    /// equality
    bool operator ==    (bits     b) const { return val == b; }
    /// unequality
    bool operator !=    (fieldbit f) const { return val != f.val; }
    /// unequality
    bool operator !=    (bits     b) const { return val != b; }
    //@}
    //--------------------------------------------------------------------------
    /// the integer value
    /// \note taking const ref instead of plain value, because otherwise gcc
    /// version 4.1.2 issues a nonsense warning (obviously a compiler BUG)
    friend value_type const&value(fieldbit const&f);
    /// the sizeof() the datum corresponding to the data field
    friend size_t const&size(fieldbit f);
    /// the single character abbreviating the data field, same as
    /// the name of the corresponding enum value in fieldbit::bits
    friend char const&letter(fieldbit f);
    /// the full name of the data field
    friend const char*fullname(fieldbit f);
    /// the 5-character function name of the data field
    friend const char*fivename(fieldbit f);
    /// the function name of the data field
    friend const char*funcname(fieldbit f);
    /// is the data field a SPH quantity?
    friend bool is_sph(fieldbit f);
    /// is the data field supported by NEMO I/O?
    friend bool is_nemo(fieldbit);
    /// is the data field an integer?
    friend bool is_integer(fieldbit);
    /// is the data field a floating point scalar?
    friend bool is_scalar(fieldbit);
    /// is the data field a vector?
    friend bool is_vector(fieldbit);
    /// is the data field a vector?
    friend bool is_vector(bits);
    /// print datum of type known only at run time
    friend std::ostream& print_field(std::ostream&, const void*, fieldbit);
  };// class fieldbit
  // ///////////////////////////////////////////////////////////////////////////
  //
  // inline definitions of some friends of class fieldbit also serve to inject
  // these functions into namespace falcON
  //
  // ///////////////////////////////////////////////////////////////////////////
  inline fieldbit::value_type const&value (fieldbit const&f) { return f.val; }
  inline size_t const&size  (fieldbit f) { return BodyData::ZQUANT[f.val]; }
  inline char const&letter  (fieldbit f) { return BodyData::SQUANT[f.val]; }
  inline const char*fullname(fieldbit f) { return BodyData::QFULLNAME[f.val]; }
  inline const char*funcname(fieldbit f) { return BodyData::QFUNCNAME[f.val]; }
  inline const char*fivename(fieldbit f) { return BodyData::QFIVENAME[f.val]; }
  inline bool       is_sph  (fieldbit f) {
    return f.val >= BodyData::KSPH && f.val < BodyData::KSINK;
  }
  // ///////////////////////////////////////////////////////////////////////////
  //
  // class falcON::fieldset
  //
  /// Represents a set of body data fields
  ///                                                                           
  /// falcON::fieldset represents a set of body data fields. Each field
  /// corresponds to a unique bit in a 64bit integer (i.e. up to 64 different
  /// fields are possible with this implementation). The position of the bit
  /// is identical to the value of the corresponding fieldbit::bits. Thus, a
  /// set containing only, say, the position has the bit fieldbit::x set and
  /// all others zero; this corresponds to the integer value 1\<\<fieldbit::x.
  /// A combination of several fields is obtained by bit-wise or. Other
  /// bit-wise operations allow the manipulation of individual fields.
  ///
  // ///////////////////////////////////////////////////////////////////////////
  class fieldset {
    //--------------------------------------------------------------------------
  public:
    /// type used in internal representation
    typedef uint64 value_type;
  private:
    value_type val;
  public:
    const static value_type nil=0, one=1;
    /// specifies single-field sets and some combined sets
    ///
    /// Single-field sets have the same name as the corresponding fieldbit::bits
    /// value. They are useful as building blocks for combined sets.
    enum bits {
      o       = nil,                 ///< an empty set
      // STD & SPH:
      m       = one << fieldbit::m,  ///< just masses
      x       = one << fieldbit::x,  ///< just positions
      v       = one << fieldbit::v,  ///< just velocities
      u       = one << fieldbit::u,  ///< just predicted velocities
      e       = one << fieldbit::e,  ///< just individual softening lengths
      f       = one << fieldbit::f,  ///< just falcON::flags
      k       = one << fieldbit::k,  ///< just integer keys
      s       = one << fieldbit::s,  ///< just time steps
      p       = one << fieldbit::p,  ///< just N-body gravitational potentials
      q       = one << fieldbit::q,  ///< just external gravitational potentials
      a       = one << fieldbit::a,  ///< just accelarations
      j       = one << fieldbit::j,  ///< just jerks
      r       = one << fieldbit::r,  ///< just mass-densities
      i       = one << fieldbit::i,  ///< just auxiliary integer data
      y       = one << fieldbit::y,  ///< just auxiliary scalar data
      z       = one << fieldbit::z,  ///< just auxiliary vector data
      l       = one << fieldbit::l,  ///< just time-step levels
      n       = one << fieldbit::n,  ///< just numbers of neighbours
      c       = one << fieldbit::c,  ///< just node numbers
      h       = one << fieldbit::h,  ///< just Peano-Hilbert keys
      d       = one << fieldbit::d,  ///< just phase-space densities
      t       = one << fieldbit::t,  ///< just orbital periods
      // SPH only:
      H       = one << fieldbit::H,  ///< just SPH smoothing lengths
      N       = one << fieldbit::N,  ///< just number of SPH iaction partners
      U       = one << fieldbit::U,  ///< just internal energies U
      Y       = one << fieldbit::Y,  ///< just predicted U
      I       = one << fieldbit::I,  ///< just (dU/dt)_total
      E       = one << fieldbit::E,  ///< just (dU/dt)_external
      K       = one << fieldbit::K,  ///< just entropy functions
      R       = one << fieldbit::R,  ///< just gas-densities rho
      A       = one << fieldbit::A,  ///< just alpha_visc
      D       = one << fieldbit::D,  ///< just div(v)
      V       = one << fieldbit::V,  ///< just predicted velocities
      J       = one << fieldbit::J,  ///< just dlnh/dt
      F       = one << fieldbit::F,  ///< just factors f_i
      C       = one << fieldbit::C,  ///< just sound speeds
      M       = one << fieldbit::M,  ///< just molecular weights
      // SINK only:
      S       = one << fieldbit::S,  ///< just spin vectors
      /// default SPH quantities
      sphdef  = H|R|V|J|F|C,
      /// all SPH quantities
      sphmax  = sphdef|N|U|Y|I|E|K|A|D|M,
      /// SPH quantities supported by NEMO I/O
      sphnemo = H|N|U|I|E|K|R|A|D|J|F|C|M,
      /// all floating point scalar SPH quantities
      sphscal = H|U|Y|I|E|K|R|A|D|J|F|C|M,
      /// N-body and external gravitational potentials
      potent  = p|q,
      /// phases = {x,v}
      phases  = x|v,
      /// masses & phases
      basic   = phases|m,
      /// masses, phases, potentials, forces, flags
      gravity = basic|p|a|f,
      /// non-SPH source properties
      source  = m|x|v|u|e|f|k|c|s,
      /// non-SPH non-source properties
      nonsource = p|q|a|j|r|i|y|z|l|n|d|h|t,
      /// all integer-type quantities
      integers= f|k|c|l|i|n|h|N,
      /// all floating point scalar quantities
      scalars = m|e|s|p|q|r|y|d|t|sphscal,
      /// all vector quantities
      vectors = x|v|a|j|u|z|S,
      /// all quantities supported by NEMO Input
      nemoin  = m|x|v|e|k|s|p|a|r|y|z|l|n|d|t|sphnemo|S,
      /// all quantities supported by NEMO Output
      nemo    = nemoin | q,
      /// all quantities at all
      all     = source|nonsource|sphmax,
      /// all SPH quantities (again)
      SPH     = sphmax,
      /// all standard quantities
      STD     = source|nonsource,
      /// all sink quantities
      SINK    = S,
      /// all non-SPH quantities
      nonSPH  = all & ~SPH,
      /// an empty set (again)
      empty   = o
    };
    //--------------------------------------------------------------------------
    /// \name constructors and assignment                                       
    //@{
    /// default constructor: empty set
    fieldset() : val(0) {}
    /// copy constructor
    fieldset(fieldset const&b) : val(b.val) {}
    /// from fieldset::bits
    fieldset(bits b) : val(b) {}
    /// from fieldbit::bits
    fieldset(fieldbit::bits b) : val(one << b) {}
    /// from value_type: take value
    explicit fieldset(value_type i) : val(i) {}
    /// from falcON::fieldbit: single-field set
    explicit fieldset(fieldbit f) : val(one << falcON::value(f)) {}
    /// from string: match each character against a single-field set's name, 
    /// corresponding to the name of the fieldset::bits value; for example,
    /// fieldset("mxv") yields a fieldset with masses, positions and velocities
    explicit fieldset(const char*c) : val(0) {
      for(fieldbit f; f; ++f)
	if(std::strchr(c,letter(f))) val |= one << falcON::value(f);
    }
    /// copy assignment operator
    fieldset&operator= (fieldset b) { val  = b.val; return *this; }
    //@}
    //--------------------------------------------------------------------------
    /// \name operations                                                        
    //@{
    /// combine with another set
    fieldset&operator|= (fieldset b) { val |= b.val; return *this; }
    /// add: combine with another set
    fieldset&operator+= (fieldset b) { val |= b.val; return *this; }
    /// subtract: delete fields from other set
    fieldset&operator-= (fieldset b) { val &= ~(b.val); return *this; }
    /// subtract: delete fields from other set
    fieldset&operator-= (bits b) { val &= ~b; return *this; }
    /// difference: fields in this set and not in other
    fieldset operator- (fieldset b) const { return fieldset(val & ~(b.val)); }
    /// difference: fields in this set and not in other
    fieldset operator- (bits b) const { return fieldset(val & ~b); }
    /// overlap: fields present in both set
    fieldset&operator&= (fieldset b) { val &= b.val; return *this; }
    /// overlap: fields present in both set
    fieldset&operator&= (bits b) { val &= b; return *this; }
    /// combination of two sets
    fieldset operator| (fieldset b) const { return fieldset(val | b.val); }
    /// combination of two sets
    fieldset operator| (bits b) const { return fieldset(val | b); }
    /// combination of two sets
    friend fieldset operator| (bits, bits);
    /// combination of two sets
    friend fieldset operator| (bits, fieldset);
    /// cross section of two sets
    fieldset operator& (fieldset b) const { return fieldset(val & b.val); }
    /// cross section of two sets
    fieldset operator& (bits b) const { return fieldset(val & b); }
    /// cross section of two sets
    friend fieldset operator& (bits, fieldset);
    /// is set empty?
    bool is_empty() const { return val == nil; }
    /// is set empty?
    operator bool() const { return val != nil; }
    /// are two sets identical?
    bool operator== (fieldset b) const { return val == b.val; }
    /// are two sets identical?
    bool operator== (bits b) const { return val == b; }
    /// are two sets identical?
    friend bool     operator== (bits b, fieldset d) { return d == b; }
    /// are two sets not identical?
    bool operator!= (fieldset b) const { return val != b.val; }
    /// are two sets not identical?
    bool operator!= (bits b) const { return val != b; }
    /// are two sets not identical?
    friend bool operator!= (bits, fieldset);
    /// complementary set
    fieldset operator~ () const { return fieldset(~val); }
    //@}
    //--------------------------------------------------------------------------
    /// \name non-operator member methods and friends                           
    //@{
    /// do two sets intersect?
    bool intersect(fieldset b) const { return val & b.val; }
    /// do two sets intersect?
    bool intersect(bits b) const { return val & b; }
    /// does this set contain \b all of set \c b ?
    bool contain(fieldset b) const { return (val & b.val)==b.val; }
    /// does this set contain \b all of set \c b ?
    bool contain(bits b) const { return (val & b)==b; }
    /// does this set contain field \c f ?
    bool contain(fieldbit f) const { return val & (one << falcON::value(f)); }
    /// does this set contain field \c f ?
    bool contain(fieldbit::bits f) const { return contain(fieldbit(f)); }
    /// fields in \c b but not in this set
    fieldset missing(fieldset b) const { return fieldset(b.val^val) & b; }
    /// fields in \c b but not in this set
    fieldset missing(bits b) const { return fieldset(b^val) & b; }
    /// return the integer value (combination of all bits)
    value_type &value() { return val; }
    /// return the integer value (combination of all bits)
    friend value_type const&value(fieldset const&);
    //@}
    //--------------------------------------------------------------------------
    // description                                                              
    char* make_word(char* w) const {
      if(val) {
	char *l = w;
	for(fieldbit f; f; ++f)
	  if(contain(f)) *(l++) = letter(f);
	*l = 0;
      } else {
	w[0] = 'o';
	w[1] = 0;
      }
      return w;
    }
    //
    class wlist {
      char W[BodyData::NQUANT+1];
    public:
      explicit wlist(const fieldset*d) { d->make_word(W); }
      operator const char* () const { return W; }
      friend const char*word(wlist const&);
    };
    const char* word() const { return wlist(this); }
    /// return a string of the letters corresponding to the fields in the set;
    /// useful for ascii output
    friend const char* word(fieldset);
    /// return a string of the letters corresponding to the fields in the set;
    /// useful for ascii output
    friend const char* word(const fieldset*);
    //--------------------------------------------------------------------------
    /// the number of bytes required to hold all fields for a single body
    size_t bytes() const {
      size_t n = 0;
      for(fieldbit f; f; ++f)
	if(contain(f)) n += size(f);
      return n;
    }
    /// the smallest # bytes to hold any one of the fields for a single body
    size_t maxsize() const {
      size_t n = 0;
      for(fieldbit f; f; ++f)
	if(contain(f)) update_max(n,size(f));
      return n;
    }
    /// formatted output: just give string of the letters corresponding to the
    /// fields in the set
    friend std::ostream& operator<< (std::ostream&, const fieldset&);
    /// formatted input: expect same format as used in output
    friend std::istream& operator>> (std::istream&, fieldset&);
  };
  // ///////////////////////////////////////////////////////////////////////////
  //
  // inline definitions of friends of class fieldset and fieldset::wlist also
  // serve to inject these functions into namespace falcON
  //
  // ///////////////////////////////////////////////////////////////////////////
  inline fieldset operator| (fieldset::bits a, fieldset::bits b) {
    return fieldset(fieldset::value_type(a) | fieldset::value_type(b));
  }
  inline fieldset operator| (fieldset::bits b, fieldset d) { return d | b; }
  inline fieldset operator& (fieldset::bits b, fieldset d) { return d & b; }
  inline bool operator!= (fieldset::bits b, fieldset d) { return d != b; }
  inline fieldset::value_type const&value(fieldset const&b) { return b.val; }
  inline const char*word(fieldset::wlist const&w) { return w.W; }
  inline const char*word(fieldset d) { return d.word(); }
  inline const char*word(const fieldset*d) { return d->word(); }
  inline std::ostream& operator<< (std::ostream&s, const fieldset&b) {
    if(value(b)) {
      for(fieldbit f; f; ++f)
	if(b.contain(f)) s << letter(f);
    } else
      s << 'o';
    return s;
  }
  inline std::istream& operator>> (std::istream&s, fieldset&b) {
    char c[32];
    s >> c;
    b = fieldset(c);
    return s;
  }
  inline bool is_integer(fieldbit f) { return fieldset::integers & 1<<f.val; }
  inline bool is_scalar (fieldbit f) { return fieldset::scalars & 1<<f.val; }
  inline bool is_vector (fieldbit f) { return fieldset::vectors & 1<<f.val; }
  inline bool is_vector (fieldbit::bits b) { return fieldset::vectors & 1<<b; }
  // ///////////////////////////////////////////////////////////////////////////
  //
  // inline definitions of a friend of class fieldbit; also serves to inject
  // this function into namespace falcON
  //
  // ///////////////////////////////////////////////////////////////////////////
  inline bool is_nemo(fieldbit f) { return fieldset::nemo & (1 << value(f)); }
  // ///////////////////////////////////////////////////////////////////////////
  //
  // data related to bodytype management
  //
  // ///////////////////////////////////////////////////////////////////////////
  const unsigned BT_NUM = 3;
  const fieldset BT_DATA[BT_NUM] = {
    fieldset::STD | fieldset::SPH | fieldset::SINK,
    fieldset::STD | fieldset::SPH,
    fieldset::STD
  };
  // ///////////////////////////////////////////////////////////////////////////
  //
  //  class falcON::bodytype
  //
  /// Represents the type of bodies
  ///
  /// Currently only two types are supported: SPH (gas particles) and STD
  /// (non-SPH bodies). bodytype can be used like a jave-style iterator to
  /// loop over all bodytypes.
  ///
  // ///////////////////////////////////////////////////////////////////////////
  class bodytype {
  public:
    /// types of bodies
    enum bits {
      sink = 0, ///< sink: sink particles
      gas  = 1, ///< gas: SPH particles
      std  = 2  ///< standard: non-SPH bodies
    };
    //--------------------------------------------------------------------------
  private:
    char val;
    //--------------------------------------------------------------------------
  public:
    /// unitialized: gas
    bodytype()                 : val(0) {}
    /// from bodytype::bits
    bodytype(bits b)           : val(b) {}
    /// copy contructor
    bodytype(bodytype const&b) : val(b.val) {}
    /// copy assignment operator
    bodytype& operator= (bodytype const&b) { val=b.val; return *this; }
    /// assignment operator from bodytype::bits
    bodytype& operator= (bits b) { val=b; return *this; }
    /// \name support for looping over bodytypes
    //@{
    /// prefix++: get next type
    bodytype& operator++()       { ++val; return *this; }
    /// is this a valid bodytype?
    operator bool       () const { return val < BT_NUM; }
    /// equality
    bool operator ==    (bodytype t) const { return val == t.val; }
    /// equality
    bool operator ==    (bits     b) const { return val == b; }
    /// inequality
    bool operator !=    (bodytype t) const { return val != t.val; }
    /// inequality
    bool operator !=    (bits     b) const { return val != b; }
    /// less than
    bool operator <     (bodytype t) const { return val <  t.val; }
    /// less than
    bool operator <     (bits     b) const { return val <  b; }
    /// less than or equal
    bool operator <=    (bodytype t) const { return val <= t.val; }
    /// less than or equal
    bool operator <=    (bits     b) const { return val <= b; }
    /// greater than
    bool operator >     (bodytype t) const { return val >  t.val; }
    /// greater than
    bool operator >     (bits     b) const { return val >  b; }
    /// greater than or equal
    bool operator >=    (bodytype t) const { return val >= t.val; }
    /// greater than or equal
    bool operator >=    (bits     b) const { return val >= b; }
    //@}
    /// conversion to int: return No of type
    operator int        () const { return val; }
    /// conversion to long: return No of type
    operator long       () const { return val; }
    /// is bodytype STD
    bool is_std() const { return val == std; }
    /// is bodytype SPH
    bool is_sph() const { return val == gas; }
    /// is bodytype SINK?
    bool is_sink() const { return val == sink; }
    /// return the full set of fields bodies of this type may hold
    fieldset allows() const { return BT_DATA[val]; }
    /// can bodies of this bodytype hold field \c f?
    bool allows(fieldbit f) const {
      return BT_DATA[val].contain(f);
    }
    /// return a name like "sink", "gas", or "std"
    const char* name() const {
      switch(val) {
      case sink: return "sink";
      case gas:  return "gas";
      default:   return "std";
      }
    }
    /// return the integer value
    char const&value() const { return val; }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //
  //  class falcON::bodytypes
  //
  /// Represents a set of different bodytypes
  ///
  // ///////////////////////////////////////////////////////////////////////////
  class bodytypes {
    /// type used in internal representation
    typedef uint8 value_type;
  private:
    value_type val;
    explicit bodytypes(int i) : val(i) {}
  public:
    const static value_type nil=0, one=1;
    /// types of bodytypes
    enum bits {
      none    = nil,                     ///< no bodies at all
      sink    = one << bodytype::sink,   ///< only sink bodies
      gas     = one << bodytype::gas,    ///< only gas bodies
      std     = one << bodytype::std,    ///< only std bodies
      nonsink = gas|std,               ///< all but sink bodies
      all     = sink|gas|std           ///< all types of bodies
    };
    //--------------------------------------------------------------------------
    /// \name constructors and assignment                                       
    //@{
    /// default ctor: no types
    bodytypes() : val(0) {}
    /// from bodytypes::bits
    bodytypes(bits b) : val(b) {}
    /// from bodytype
    explicit bodytypes(bodytype t) : val(one<<int(t)) {}
    /// from bodytype::bits
    explicit bodytypes(bodytype::bits t) : val(one<<t) {}
    /// copy constructor
    bodytypes(bodytypes const&b) : val(b.val) {}
    /// copy assignment
    bodytypes& operator=(bodytypes const&b) { val = b.val; return*this; }
    //@}
    //--------------------------------------------------------------------------
    /// \name operations                                                        
    //@{
    /// combine with other bodytypes
    bodytypes&operator|= (bodytypes b) { val |= b.val; return *this; }
    /// combination of bodytypes
    bodytypes operator| (bodytypes b) const { return bodytypes(val | b.val); }
    /// combination of bodytypes
    bodytypes operator| (bits b) const { return bodytypes(val | b); }
    /// cross section of bodytypes
    bodytypes operator& (bodytypes b) const { return bodytypes(val & b.val); }
    /// cross section of bodytypes
    bodytypes operator& (bits b) const { return bodytypes(val & b); }
    /// no bodytypes at all?
    operator bool() const { return val != 0; }
    /// are bodytypes identical?
    bool operator== (bodytypes b) const { return val == b.val; }
    /// are bodytypes identical?
    bool operator== (bits b) const { return val == b; }
    /// are bodytypes not identical?
    bool operator!= (bodytypes b) const { return val != b.val; }
    /// are bodytypes not identical?
    bool operator!= (bits b) const { return val != b; }
    //@}
    //--------------------------------------------------------------------------
    /// \name non-operator member methods                                       
    //@{
    /// is a certain type within our range?
    bool contain(bodytype t) const { return val & one<<int(t); }
    /// are all types b within our range?
    bool contain(bodytypes b) const { return (val & b.val) == b.val; }
    /// return the integer value (combination of all bits)
    int value() const { return val; }
    //@}
  };
  // ///////////////////////////////////////////////////////////////////////////
  //
  // struct field_zero<type>
  //
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T> struct field_zero {
    static void set_zero(T&x) { x = T(0); }
  };
  template<> struct field_zero<flags> {
    static void set_zero(flags&x) { x = flags::empty; }
  };
  template<typename T, int N> struct field_zero< tupel<N,T> > {
    static void set_zero(tupel<N,T>&x) { x = T(0); }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //
  // struct field_type<type>
  //
  // ///////////////////////////////////////////////////////////////////////////
  template<typename TYPE> struct field_type
  : public field_zero<TYPE> {
    //
    field_zero<TYPE>::set_zero;
    //
    static TYPE*array(void*D) {
      return static_cast<TYPE*>(D);}
    //
    static const TYPE*c_array(const void*D) {
      return static_cast<const TYPE*>(D);}
  };
  //////////////////////////////////////////////////////////////////////////////
  //
  // struct field_traits<int>  int = fieldbit
  //
  // used in templated members of class bodies below.
  // define all related types; aid nemo I/O
  //
  //////////////////////////////////////////////////////////////////////////////
  /// \brief Template over body data field
  /// provides information about data type and name etc for data field.
  template<int I> struct field_traits;
#  define DefFieldTraits(BIT,TYPE)					\
  template<> struct field_traits< BIT > :				\
  public field_type<TYPE> {						\
    typedef TYPE type;							\
    static const bool   is_sph = BIT >= BodyData::KSPH			\
                              && BIT < BodyData::KSINK;			\
    static const size_t size   = sizeof(type);                          \
    static const char&word     () { return BodyData::SQUANT[BIT]; }	\
    static const char*fullname () { return BodyData::QFULLNAME[BIT]; }	\
    static const char*funcname () { return BodyData::QFUNCNAME[BIT]; }	\
    static const char*fivename () { return BodyData::QFIVENAME[BIT]; }	\
  };
  //----------------------------------------------------------------------------
  DefFieldTraits(fieldbit::m, real);               // mass                      
  DefFieldTraits(fieldbit::x, vect);               // position                  
  DefFieldTraits(fieldbit::v, vect);               // velocity                  
  DefFieldTraits(fieldbit::u, vect);               // predicted velocity        
  DefFieldTraits(fieldbit::e, real);               // softening length          
  DefFieldTraits(fieldbit::f, flags)               // body flags                
  DefFieldTraits(fieldbit::k, unsigned);           // body key                  
  DefFieldTraits(fieldbit::s, double);             // time step                 
  DefFieldTraits(fieldbit::p, real);               // internal potential        
  DefFieldTraits(fieldbit::q, real);               // external potential        
  DefFieldTraits(fieldbit::a, vect);               // acceleration              
  DefFieldTraits(fieldbit::j, vect);               // jerk                      
  DefFieldTraits(fieldbit::r, real);               // mass density              
  DefFieldTraits(fieldbit::i, int);                // auxiliary integer         
  DefFieldTraits(fieldbit::y, real);               // auxiliary scalar          
  DefFieldTraits(fieldbit::z, vect);               // auxiliary vector          
  DefFieldTraits(fieldbit::l, indx);               // time-step level           
  DefFieldTraits(fieldbit::n, unsigned);           // # neighbours              
  DefFieldTraits(fieldbit::c, indx);               // node number               
  DefFieldTraits(fieldbit::h, peanokey);           // Peano-Hilbert key         
  DefFieldTraits(fieldbit::d, real);               // phase-space density       
  DefFieldTraits(fieldbit::t, real);               // orbital periods           
  DefFieldTraits(fieldbit::H, real);               // SPH: smoothing length h   
  DefFieldTraits(fieldbit::N, unsigned);           // SPH: # neighbours         
  DefFieldTraits(fieldbit::U, real);               // SPH: internal energy U    
  DefFieldTraits(fieldbit::Y, real);               // SPH: predicted U_in       
  DefFieldTraits(fieldbit::I, real);               // SPH: (dU/dt)_total        
  DefFieldTraits(fieldbit::E, real);               // SPH: (dU/dt)_external     
  DefFieldTraits(fieldbit::K, real);               // SPH: entropy function     
  DefFieldTraits(fieldbit::R, real);               // SPH: gas density          
  DefFieldTraits(fieldbit::A, real);               // SPH: alpha_visc           
  DefFieldTraits(fieldbit::D, real);               // SPH: div(v)               
  DefFieldTraits(fieldbit::J, real);               // SPH: dh/dt                
  DefFieldTraits(fieldbit::F, real);               // SPH: f_i                  
  DefFieldTraits(fieldbit::C, real);               // SPH: sound speed          
  DefFieldTraits(fieldbit::M, real);               // SPH: molecular weights    
  DefFieldTraits(fieldbit::S, vect);               // SINK: spins               
#undef DefFieldTraits
  //////////////////////////////////////////////////////////////////////////////
  //
  // macro for defining named functions
  //
  // named functions MUST not be defined in any other way to ensure that any
  // changes in the data type layout etc. are automatically reflected.
  //
  //////////////////////////////////////////////////////////////////////////////
#define DEF_NAMED_STD(MACRO)			\
  MACRO(fieldbit::m,mass);			\
  MACRO(fieldbit::x,pos);			\
  MACRO(fieldbit::v,vel);			\
  MACRO(fieldbit::u,vprd);			\
  MACRO(fieldbit::e,eps);			\
  MACRO(fieldbit::f,flag);			\
  MACRO(fieldbit::k,key);			\
  MACRO(fieldbit::s,step);			\
  MACRO(fieldbit::p,pot);			\
  MACRO(fieldbit::q,pex);			\
  MACRO(fieldbit::a,acc);			\
  MACRO(fieldbit::j,jerk);			\
  MACRO(fieldbit::r,rho);			\
  MACRO(fieldbit::i,integ);			\
  MACRO(fieldbit::y,aux);			\
  MACRO(fieldbit::z,zet);			\
  MACRO(fieldbit::l,level);			\
  MACRO(fieldbit::n,num);			\
  MACRO(fieldbit::c,node);			\
  MACRO(fieldbit::h,peano);			\
  MACRO(fieldbit::d,phden);			\
  MACRO(fieldbit::t,torb);
#define DEF_NAMED_SPH(MACRO)			\
  MACRO(fieldbit::H,size);			\
  MACRO(fieldbit::N,snum);			\
  MACRO(fieldbit::U,uin);			\
  MACRO(fieldbit::Y,uprd);			\
  MACRO(fieldbit::I,udot);			\
  MACRO(fieldbit::E,udex);			\
  MACRO(fieldbit::K,entr);			\
  MACRO(fieldbit::R,srho);			\
  MACRO(fieldbit::A,alfa);			\
  MACRO(fieldbit::D,divv);			\
  MACRO(fieldbit::J,hdot);			\
  MACRO(fieldbit::F,fact);			\
  MACRO(fieldbit::C,csnd);			\
  MACRO(fieldbit::M,molw);
#define DEF_NAMED_SINK(MACRO)			\
  MACRO(fieldbit::S,spin);
#define DEF_NAMED(MACRO)			\
  DEF_NAMED_STD(MACRO)				\
  DEF_NAMED_SPH(MACRO)				\
  DEF_NAMED_SINK(MACRO)
#define DEF_VECTORS(MACRO)			\
  MACRO(fieldbit::x,pos);			\
  MACRO(fieldbit::v,vel);			\
  MACRO(fieldbit::u,vprd);			\
  MACRO(fieldbit::a,acc);			\
  MACRO(fieldbit::j,jerk);			\
  MACRO(fieldbit::S,spin);
  // ///////////////////////////////////////////////////////////////////////////
  inline std::ostream& print_field(std::ostream&s, const void*x, fieldbit f)
  {
    switch(value(f)) {
#define PRINT(BIT,NAME)							\
  case BIT: return s <<							\
      *(static_cast<const field_traits<BIT>::type *>(x))
    DEF_NAMED(PRINT)
#undef PRINT
    default:
      falcON_THROWN("falcON::print(fieldbit %c): unknown field\n",
		    letter(f));
      return s;
    }
  }
  //////////////////////////////////////////////////////////////////////////////
  //
  // templates for looping fields
  //
  //////////////////////////////////////////////////////////////////////////////
  template< template<int> class LOOP, int BIT, int BITEND >
  struct LoopFields {
    template< typename T > static void loop(T&x) { 
      LOOP<BIT>::act(x);
      LoopFields<LOOP, BIT+1, BITEND>::loop(x);
    }
    template< typename T > static void const_loop(T const&x) { 
      LOOP<BIT>::act(x);
      LoopFields<LOOP, BIT+1, BITEND>::const_loop(x);
    }
  };
  template< template<int> class LOOP, int BIT >
  struct LoopFields<LOOP, BIT, BIT>  {
    template< typename T > static void loop(T&) {}
    template< typename T > static void const_loop(T const&) {}
  };
  //----------------------------------------------------------------------------
  template< template<int> class LOOP>
  struct LoopAllFields {
    template< typename T > static void loop(T&x) {
      LoopFields<LOOP, 0, BodyData::NQUANT>::loop(x);
    }
    template< typename T > static void const_loop(T const&x) {
      LoopFields<LOOP, 0, BodyData::NQUANT>::const_loop(x);
    }
  };
  //----------------------------------------------------------------------------
  template< template<int> class LOOP>
  struct LoopSINKFields {
    template< typename T > static void loop(T&x) {
      LoopFields<LOOP, BodyData::KSINK, BodyData::NQUANT>::loop(x);
    }
    template< typename T > static void const_loop(T const&x) {
      LoopFields<LOOP, BodyData::KSINK, BodyData::NQUANT>::const_loop(x);
    }
  };
  //----------------------------------------------------------------------------
  template< template<int> class LOOP>
  struct LoopSPHFields {
    template< typename T > static void loop(T&x) {
      LoopFields<LOOP, BodyData::KSPH, BodyData::KSINK>::loop(x);
    }
    template< typename T > static void const_loop(T const&x) {
      LoopFields<LOOP, BodyData::KSPH, BodyData::KSINK>::const_loop(x);
    }
  };
  //----------------------------------------------------------------------------
  template< template<int> class LOOP>
  struct LoopSTDFields {
    template< typename T > static void loop(T&x) {
      LoopFields<LOOP, 0, BodyData::KSPH>::loop(x);
    }
    template< typename T > static void const_loop(T const&x) {
      LoopFields<LOOP, 0, BodyData::KSPH>::const_loop(x);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
} // namespace falcON {
falcON_TRAITS(falcON::flags,"flags");
falcON_TRAITS(falcON::fieldbit,"fieldbit");
falcON_TRAITS(falcON::fieldset,"fieldset");
falcON_TRAITS(falcON::bodytype,"bodytype");
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_fields_h
