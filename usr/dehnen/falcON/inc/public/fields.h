// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    inc/public/fields.h                                                
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    2002-2005                                                          
///                                                                             
/// \brief   body flags and support for body data management                    
///                                                                             
///          Contents:                                                          
///          \li falcON::flag: holding boolean information about each body      
///          \li falcON::fieldbit: representing a single body data field        
///          \li falcON::fieldset: representing a set of body data fields       
///          \li template falcON::field_traits                                  
///          \li further contants and macros                                    
///                                                                             
////////////////////////////////////////////////////////////////////////////////
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
#ifndef falcON_included_fields_h
#define falcON_included_fields_h 1

#ifndef falcON_included_types_h
#  include <public/types.h>
#endif

////////////////////////////////////////////////////////////////////////////////
namespace falcON {

  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::flag                                                         
  //                                                                            
  /// An integer, the bits of which are used as flags for body properties;      
  /// each body, as well as tree cells and leafs have a flag.                   
  // ///////////////////////////////////////////////////////////////////////////
  class flag {
    //--------------------------------------------------------------------------
    // type of class flag                                                       
    //--------------------------------------------------------------------------
  public:
    /// enumeration holding some often used flag values.                        
    enum {
      ACTIVE        =  1<<0,  ///< a body/tree node is active                   
      REMOVE        =  1<<1,  ///< body to be removed; see bodies::remove()     
      SPH           =  1<<2,  ///< a body is SPH, tree node contains SPH bodies 
      STICKY        =  1<<3,  ///< a body is sticky particle, tree node contains
      AL_ACTIVE     =  1<<4,  ///< all children of a cell are active            
      AL_SPH        =  1<<5,  ///< all children of a cell are SPH               
      AL_STICKY     =  1<<6,  ///< all children of a cell are sticky            
      SPH_SPECIAL   =  1<<7,  ///< help with time integration for SPH           
      ADJUST_STEP   =  1<<8,  ///< help with time integration                   
      // compounds of flags                                                     
      LEAF_FLAGS    =  ACTIVE|SPH|STICKY,     ///< sum of all leaf flags        
      BODY_FLAGS    =  LEAF_FLAGS|REMOVE,     ///< sum of all body flags        
      ACTIVE_FLAGS  =  ACTIVE|AL_ACTIVE,      ///< sum of all active flags      
      SPH_FLAGS     =  SPH   |AL_SPH,         ///< sum of all sph  flags        
      STICKY_FLAGS  =  STICKY|AL_STICKY       ///< sum of all sticky flags      
    };
    //--------------------------------------------------------------------------
    // data of class flag                                                       
    //--------------------------------------------------------------------------
  private:
    /// data: 32 bits -> up to 32 flags can be supported
    int FLAG;
    //--------------------------------------------------------------------------
    /// \name constructors                                                      
    //@{
  public:
    /// without argument: all bits are unset
    flag           ()             : FLAG(0) {}
    /// set equal to integer
    flag           (int  F)       : FLAG(F) {}
    /// copy constructor
    flag           (flag const&F) : FLAG(F.FLAG) {}
    /// copy operator
    flag& operator=(flag const&F) { FLAG = F.FLAG; return *this; }
    //@}
    //--------------------------------------------------------------------------
    /// \name non-const methods                                                 
    //@{
    /// reset all bits to zero
    void reset       ()                         { FLAG = 0; }
    /// set equal to another flag
    void set_to      (const flag*F)             { FLAG = F->FLAG; }
    /// set equal to another flag
    void set_to      (const flag&F)             { FLAG = F.FLAG; }
    /// set equal to part of another flag
    void set_to_part (const flag*F,const int&P) { FLAG = F->FLAG&P;}
    /// set equal to part of another flag
    void set_to_part (const flag&F,const int&P) { FLAG = F.FLAG&P;}
    /// add all bits of another flag
    void add         (const flag*F)             { FLAG|= F->FLAG; }
    /// add all bits of another flag
    void add         (const flag&F)             { FLAG|= F.FLAG; }
    /// add part of bits of another flag
    void add_part    (const flag*F,const int&P) { FLAG|= F->FLAG&P;}
    /// add part of bits of another flag
    void add_part    (const flag&F,const int&P) { FLAG|= F.FLAG&P;}
    /// set part of us equal to corresponding part of another flag
    void set_part    (const flag*F,const int&P) { FLAG = FLAG&~P | F->FLAG&P; }
    /// set part of us equal to corresponding part of another flag
    void set_part    (const flag&F,const int&P) { FLAG = FLAG&~P | F.FLAG&P; }
    /// reset all bits which are set in another flag
    void un_set      (const flag*F)             { FLAG&= ~F->FLAG; }
    /// reset all bits which are set in another flag
    void un_set      (const flag&F)             { FLAG&= ~F.FLAG; }
    //@}
    //--------------------------------------------------------------------------
    // boolean information via members                                          
    //--------------------------------------------------------------------------
    /// is any bit set in the argument set here?
    bool is_set(int T) const { return FLAG & T; }
    //--------------------------------------------------------------------------
    /// \name boolean information via friends                                   
    //@{
    /// is bit indicating activity set?
    friend bool is_active (const flag*F) { return F->is_set(ACTIVE); }
    /// is bit indicating activity set?
    friend bool is_active (flag const&F) { return F.is_set(ACTIVE); }
    /// is bit indicating activity of all children set?
    friend bool al_active (const flag*F) { return F->is_set(AL_ACTIVE); }
    /// is bit indicating activity of all children set?
    friend bool al_active (flag const&F) { return F.is_set(AL_ACTIVE); }
    /// is bit indicating SPH set?
    friend bool is_sph    (const flag*F) { return F->is_set(SPH); }
    /// is bit indicating SPH set?
    friend bool is_sph    (flag const&F) { return F.is_set(SPH); }
    /// is bit indicating SPH of all children set?
    friend bool al_sph    (const flag*F) { return F->is_set(AL_SPH); }
    /// is bit indicating SPH of all children set?
    friend bool al_sph    (flag const&F) { return F.is_set(AL_SPH); }
    /// is bit indicating sticky set?
    friend bool is_sticky (const flag*F) { return F->is_set(STICKY); }
    /// is bit indicating sticky set?
    friend bool is_sticky (flag const&F) { return F.is_set(STICKY); }
    /// is bit indicating sticky of all children set?
    friend bool al_sticky (const flag*F) { return F->is_set(AL_STICKY); }
    /// is bit indicating sticky of all children set?
    friend bool al_sticky (flag const&F) { return F.is_set(AL_STICKY); }
    /// is any bit set in argument set here?
    friend bool is_set    (const flag*F, int const&T) { return F->is_set(T); }
    /// is any bit set in argument set here?
    friend bool is_set    (flag const&F, int const&T) { return F.is_set(T); }
    //--------------------------------------------------------------------------
    //@}
    //--------------------------------------------------------------------------
    // conversion to int                                                        
    //--------------------------------------------------------------------------
    /// conversion to const int
    operator const int& () const { return FLAG; }
    /// conversion to int
    operator       int& ()       { return FLAG; }
    //--------------------------------------------------------------------------
    /// \name formatted I/O                                                     
    //@{
    /// formatted output: simply write integer value
    friend std::ostream& operator<< (std::ostream& o, const flag& F) {
      return o<<F.FLAG;
    }
    /// formatted input: read as integer
    friend std::istream& operator>> (std::istream& i, flag& F) {
      return i>>F.FLAG;
    }
    //@}
    //--------------------------------------------------------------------------
  };
  falcON_TRAITS(flag,"flag");

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // N-body data type management                                              //
  //                                                                          //
  // There are many different N-body data (mass, position, etc), each with    //
  // its unique data type and I/O requirements etc. Below are several types   //
  // and constant designed for managing this.                                 //
  //                                                                          //
  // class fieldbit                                                           //
  // represents individual body data types as integer, i.e. 0 corresponds to  //
  // mass, 1 to positions, etc. fieldbit can be used as a java-style iterator //
  // over allowes data types.                                                 //
  //                                                                          //
  // class fieldset                                                           //
  // represents a collection of body data types as integer each bit           //
  // corresponding to the data whose fieldbit matches its position. fieldset  //
  // are used to indicate which body data should be, e.g., allocated, read,   //
  // written, or copies, etc.                                                 //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // some global constants related to N-body data                             //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////

  const int         BD_NOSPH   = 18;               // # non-sph data            
  const int         BD_NQUANT  = 30;               // total # data              
  const char* const BD_SQUANT  ="mxvwefktpqajryzlnhHNUYIESRDJFC";
                                                   // one-char data tags        
  //----------------------------------------------------------------------------
  // array with human readable names for the N-body data                        
  const char* const BD_QNAME[BD_NQUANT] = 
    { "mass", "position", "velocity", "predicted velocity", "softening length",
      "body flag", "key", "time step", "potential", "external potential",
      "acceleration", "jerk", "density",  "auxiliary scalar", 
      "auxiliary vector", "time-step level", "number of partners",
      "Peano-Hilbert key", "SPH smoothing length",
      "number of SPH partners", "U_internal", "U_predicted", "(dU/dt)_total",
      "(dU/dt)_external", "entropy", "gas density", "d(gas density)/dt",
      "dlog(h)/dt", "factor", "sound speed"
    };
  //----------------------------------------------------------------------------
  // array with the sizeof() the N-body data                                    
  const size_t BD_ZQUANT[BD_NQUANT] = {
    //            source properties: 8
    sizeof(real),      // mass
    sizeof(vect),      // position
    sizeof(vect),      // velocity
    sizeof(vect),      // velocity_predicted
    sizeof(real),      // softening length
    sizeof(flag),      // flag
    sizeof(int),       // key
    sizeof(real),      // time step
    //            sink properties: 10
    sizeof(real),      // potential
    sizeof(real),      // external potential
    sizeof(vect),      // acceleration
    sizeof(vect),      // jerk
    sizeof(real),      // density (mass)
    sizeof(real),      // auxiliary scalar
    sizeof(vect),      // auxiliary scalar
    sizeof(indx),      // level
    sizeof(unsigned),  // number of partners
    sizeof(peanokey),  // Peano-Hilbert key
    //            SPH properties: 12
    sizeof(real),      // size
    sizeof(unsigned),  // number of SPH partners
    sizeof(real),      // U_internal
    sizeof(real),      // U_predicted
    sizeof(real),      // (dU/dt)_total
    sizeof(real),      // (dU/dt)_external
    sizeof(real),      // entropy
    sizeof(real),      // gas-density
    sizeof(real),      // dgas-density/dt
    sizeof(real),      // dlogh/dt
    sizeof(real),      // f_i
    sizeof(real)       // sound speed
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::fieldbit                                                     
  //                                                                            
  /// Represents a single individual body data field                            
  ///                                                                           
  /// Each body data field corresponds to another integer: 0 means mass, 1      
  /// positions, etc. In any application, \b never use these integer values, but
  /// \b always the enumeration value in enum fieldbit::bits, because these are 
  /// guaranteed to be preserved during code changes. fieldbit can be used as a 
  /// java-style iterator over data fields.                                     
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class fieldbit {
    int val;
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
      w       = 3,            ///< predicted velocity
      e       = 4,            ///< individual gravitational softening length
      f       = 5,            ///< falcON::flag
      k       = 6,            ///< integer key
      t       = 7,            ///< time step (not used in falcON)
      //            sink properties: 10
      p       = 8,            ///< N-body gravitational potential
      q       = 9,            ///< external gravitational potential
      a       = 10,           ///< accelaration (gravity and SPH)
      j       = 11,           ///< jerk (not used in falcON)
      r       = 12,           ///< mass-density
      y       = 13,           ///< auxiliary scalar datum
      z       = 14,           ///< auxiliary vector datum
      l       = 15,           ///< time-step level (short)
      n       = 16,           ///< number of neighbours
      h       = 17,           ///< Peano-Hilbert key
      //            SPH properties: 12
      H       = 18,           ///< SPH smoothing length
      N       = 19,           ///< number of SPH interaction partners
      U       = 20,           ///< internal energy U
      Y       = 21,           ///< predicted U
      I       = 22,           ///< (dU/dt)_total
      E       = 23,           ///< (dU/dt)_external
      S       = 24,           ///< entropy/entropy function (currently not used)
      R       = 25,           ///< gas-density rho
      D       = 26,           ///< drho/dt
      V       = w,            ///< predicted velocity (again)
      J       = 27,           ///< dh/dt or dlnh/dt
      F       = 28,           ///< factor f_i
      C       = 29,           ///< sound speed
      invalid = BD_NQUANT     ///< not corresponding to any field
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
	if(l == BD_SQUANT[v] ) {
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
    friend int         const&value   (fieldbit f) { return f.val; }
    /// the sizeof() the datum corresponding to the data field
    friend size_t      const&size    (fieldbit f) { return BD_ZQUANT[f.val]; }
    /// the single character abbreviating the data field, same as
    /// the name of the corresponding enum value in fieldbit::bits
    friend char        const&letter  (fieldbit f) { return BD_SQUANT[f.val]; }
    /// the name of the data field
    friend const char* const&name    (fieldbit f) { return BD_QNAME [f.val]; }
    /// is the data field a SPH quantity?
    friend bool              is_sph  (fieldbit f) { return f.val >= BD_NOSPH; }
    /// is the data field supported by NEMO I/O?
    friend bool              is_nemo (fieldbit);
  };
  falcON_TRAITS(fieldbit,"fieldbit");

  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::fieldset                                                     
  //                                                                            
  /// Represents a set of body data fields                                      
  ///                                                                           
  /// falcON::fieldset represents a set of body data fields. Each field         
  /// corresponds to a unique bit in a 32bit integer (i.e. up to 32 different   
  /// fields are possible with this implementation). The position of the bit    
  /// is identical to the value of the corresponding fieldbit::bits. Thus, a    
  /// set containing only, say, the position has the bit fieldbit::x set and    
  /// all others zero; this corresponds to the integer value 1\<\<fieldbit::x.  
  /// A combination of several fields is obtained by bit-wise or. Other bit-wise
  /// operations allow the manipulation of individual fields.                   
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class fieldset {
    //--------------------------------------------------------------------------
    unsigned val;
  public:
    //--------------------------------------------------------------------------
    /// specifies single-field sets and some combined sets
    ///
    /// Single-field sets have the same name as the corresponding fieldbit::bits
    /// value. They are useful as building blocks for combined sets.
    enum bits {
      o       = 0,                  ///< an empty set
      m       = 1 << fieldbit::m,   ///< just masses
      x       = 1 << fieldbit::x,   ///< just positions
      v       = 1 << fieldbit::v,   ///< just velocities
      w       = 1 << fieldbit::w,   ///< just predicted velocities
      e       = 1 << fieldbit::e,   ///< just individual softening lengths
      f       = 1 << fieldbit::f,   ///< just falcON::flags
      k       = 1 << fieldbit::k,   ///< just integer keys
      t       = 1 << fieldbit::t,   ///< just time steps
      p       = 1 << fieldbit::p,   ///< just N-body gravitational potentials
      q       = 1 << fieldbit::q,   ///< just external gravitational potentials
      a       = 1 << fieldbit::a,   ///< just accelarations
      j       = 1 << fieldbit::j,   ///< just jerks
      r       = 1 << fieldbit::r,   ///< just mass-densities
      y       = 1 << fieldbit::y,   ///< just auxiliary scalar data
      z       = 1 << fieldbit::z,   ///< just auxiliary vector data
      l       = 1 << fieldbit::l,   ///< just time-step levels
      n       = 1 << fieldbit::n,   ///< just numbers of neighbours
      h       = 1 << fieldbit::h,   ///< just Peano-Hilbert keys
      H       = 1 << fieldbit::H,   ///< just SPH smoothing lengths
      N       = 1 << fieldbit::N,   ///< just number of SPH interaction partners
      U       = 1 << fieldbit::U,   ///< just internal energies U
      Y       = 1 << fieldbit::Y,   ///< just predicted U
      I       = 1 << fieldbit::I,   ///< just (dU/dt)_total
      E       = 1 << fieldbit::E,   ///< just (dU/dt)_external
      S       = 1 << fieldbit::S,   ///< just entropies/entropy functions
      R       = 1 << fieldbit::R,   ///< just gas-densities rho
      D       = 1 << fieldbit::D,   ///< just drho/dt
      V       = 1 << fieldbit::V,   ///< just predicted velocities
      J       = 1 << fieldbit::J,   ///< just dh/dt or dlnh/dt
      F       = 1 << fieldbit::F,   ///< just factors f_i
      C       = 1 << fieldbit::C,   ///< just sound speeds
      /// default SPH quantities
      sphdef  = H|R|V|J|F|C,
      /// all SPH quantities
      sphmax  = sphdef|N|U|Y|I|E|S|D,
      /// SPH quantities supported by NEMO I/O
      sphnemo = H|N|U|I|E|S|R|J|F|C,
      /// all floating point scalar SPH quantities
      sphscal = H|U|Y|I|E|S|R|D|J|F|C,
      /// N-body and external gravitational potentials
      potent  = p|q,
      /// phases = {x,v}
      phases  = x|v,
      /// masses & phases
      basic   = phases|m,
      /// masses, phases, potentials, forces, flags
      gravity = basic|p|a|f,
      /// non-SPH source properties
      source  = m|x|v|w|e|f|k|t,
      /// non-SPH sink properties
      sink    = p|q|a|j|r|y|z|l|n|h,
      /// all integer-type quantities
      integers= f|k|l|n|h|N,
      /// all floating point scalar quantities
      scalars = m|e|t|p|q|r|y|sphscal,
      /// all vector quantities
      vectors = x|v|a|j|w|z,
      /// all quantities supported by NEMO Input
      nemoin  = m|x|v|e|f|k|t|p|a|r|y|l|n|sphnemo,
      /// all quantities supported by NEMO Output
      nemo    = nemoin | q,
      /// all quantities at all
      all     = source|sink|sphmax,
      /// all SPH quantities (again)
      SPH     = sphmax,
      /// all standard quantities
      STD     = source|sink,
      /// all non-SPH quantities
      nonSPH  = all & ~SPH,
      /// an empty set (again)
      empty   = o
    };
    //==========================================================================
    /// \name constructors and assignment                                       
    //@{
    /// unitialized: empty set
    fieldset()                    : val(0) {}
    /// copy constructor
    fieldset(fieldset const   &b) : val(b.val) {}
    /// from fieldset::bits
    fieldset(bits              b) : val(b) {}
    /// from integer: take value
    explicit fieldset(int      i) : val(i) {}
    /// from unsigned: take value
    explicit fieldset(unsigned i) : val(i) {}
    /// from falcON::fieldbit: single-field set
    explicit fieldset(fieldbit f) : val(1 << falcON::value(f)) {}
    /// from char: single-field set, match against \b name of enum values in
    /// fieldset::bits; for example, fieldset('m') gives the single-field set
    /// with just masses
    explicit fieldset(char     c) : val(0) {
      for(fieldbit f; f; ++f)
	if(c == letter(f) ) {
	  val |= 1 << falcON::value(f);
	  return;
	}
    }
    /// from string: match each character against a single-field set's name, 
    /// corresponding to the name of the fieldset::bits value; for example,
    /// fieldset("mxv") yields a fieldset with masses, positions and velocities
    explicit fieldset(const char    *c) : val(0) {
      for(fieldbit f; f; ++f)
	if(std::strchr(c,letter(f))) val |= 1 << falcON::value(f);
    }
    /// copy assignment operator
    fieldset&operator = (fieldset b) { val  = b.val; return *this; }
    //@}
    //==========================================================================
    /// \name operations                                                        
    //@{
    /// combine with another set
    fieldset&operator|= (fieldset b) { val |= b.val; return *this; }
    /// add: combine with another set
    fieldset&operator+= (fieldset b) { val |= b.val; return *this; }
    /// subtract: delete fields from other set
    fieldset&operator-= (fieldset b) { val &= ~(b.val); return *this; }
    /// subtract: delete fields from other set
    fieldset&operator-= (bits     b) { val &= ~b; return *this; }
    /// difference: fields in this set and not in other
    fieldset operator-  (fieldset b) { return fieldset(val & ~(b.val)); }
    /// overlap: fields present in both set
    fieldset&operator&= (fieldset b) { val &= b.val; return *this; }
    /// overlap: fields present in both set
    fieldset&operator&= (bits     b) { val &= b; return *this; }
    /// combination of two sets
    fieldset operator|  (fieldset b) const { return fieldset(val | b.val); }
    /// combination of two sets
    fieldset operator|  (bits     b) const { return fieldset(val | b); }
    /// combination of two sets
    friend fieldset operator|  (bits b, fieldset d) { return d | b; }
    /// cross section of two sets
    fieldset operator&  (fieldset b) const { return fieldset(val & b.val); }
    /// cross section of two sets
    fieldset operator&  (bits     b) const { return fieldset(val & b); }
    /// cross section of two sets
    friend fieldset operator&  (bits b, fieldset d) { return d & b; }
    /// is set empty?
    bool     is_empty   ()           const { return val != 0; }
    /// is set empty?
    operator bool       ()           const { return val != 0; }
    /// are two sets identical?
    bool     operator== (fieldset b) const { return val == b.val; }
    /// are two sets identical?
    bool     operator== (bits     b) const { return val == b; }
    /// are two sets identical?
    friend bool     operator== (bits b, fieldset d) { return d == b; }
    /// are two sets not identical?
    bool     operator!= (fieldset b) const { return val != b.val; }
    /// are two sets not identical?
    bool     operator!= (bits     b) const { return val != b; }
    /// are two sets not identical?
    friend bool     operator!= (bits b, fieldset d) { return d != b; }
    /// complementary set
    fieldset operator~  ()           const { return fieldset(~val); }
    //@}
    //==========================================================================
    // non-operator member methods and friends                                  
    //==========================================================================
    /// do two sets intersect?
    bool      intersect  (fieldset b) const { return val & b.val; }
    /// do two sets intersect?
    bool      intersect  (bits     b) const { return val & b; }
    /// does this set contain \b all of set \c b ?
    bool      contain    (fieldset b) const { return (val & b.val)==b.val; }
    /// does this set contain \b all of set \c b ?
    bool      contain    (bits     b) const { return (val & b)==b; }
    /// does this set contain field \c f ?
    bool      contain    (fieldbit f) const { 
      return val & (1 << falcON::value(f)); }
    /// does this set contain field \c f ?
    bool      contain    (fieldbit::bits f) const {
      return contain(fieldbit(f)); }
    /// fields in \c b but not in this set
    fieldset  missing    (fieldset b) const { return fieldset(b.val^val) & b; }
    /// fields in \c b but not in this set
    fieldset  missing    (bits     b) const { return fieldset(b^val) & b; }
    /// return the integer value (combination of all bits)
    unsigned const&value      ()           const { return val; }
    /// return the integer value (combination of all bits)
    friend 
    unsigned const&value      (fieldset const&b) { return b.val; }
    //==========================================================================
    // description                                                              
    //==========================================================================
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
    //--------------------------------------------------------------------------
    class wlist {
      char W[BD_NQUANT+1];
    public:
      explicit wlist(const fieldset*d) { d->make_word(W); }
      operator const char* () const { return W; }
    };
    const char* word() const { return wlist(this); }
    /// return a string of the letters corresponding to the fields in the set;
    /// useful for ascii output
    friend const char* word(      fieldset d) { return d. word(); }
    /// return a string of the letters corresponding to the fields in the set;
    /// useful for ascii output
    friend const char* word(const fieldset*d) { return d->word(); }
    //==========================================================================
    /// the number of bytes required to hold all fields for a single body
    inline size_t  bytes       ()            const {
      size_t n = 0;
      for(fieldbit f; f; ++f)
	if(contain(f)) n += size(f);
      return n;
    }
    //==========================================================================
    /// formatted output: just give string of the letters corresponding to the
    /// fields in the set
    friend std::ostream& operator<< (std::ostream&s, const fieldset&b) {
      if(b.value()) {
	for(fieldbit f; f; ++f)
	  if(b.contain(f)) s << letter(f);
      } else
	s << 'o';
      return s;
    }
    //==========================================================================
    /// formatted input: expect same format as used in output
    friend std::istream& operator>> (std::istream&s, fieldset&b) {
      char c[32];
      s >> c;
      b = fieldset(c);
      return s;
    }
  };
  falcON_TRAITS(fieldset,"fieldset");
  //////////////////////////////////////////////////////////////////////////////
  inline bool is_nemo(fieldbit f) {
    return fieldset::nemo & 1 << value(f);
  }
  // ///////////////////////////////////////////////////////////////////////////
  const unsigned BT_NUM = 2;
  const unsigned BT_DATA[BT_NUM] = { fieldset::STD | fieldset::SPH,
				     fieldset::STD };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  //  class falcON::bodytype                                                    
  //                                                                            
  /// Represents the type of bodies                                             
  ///                                                                           
  /// Currently only two types are supported: SPH (gas particles) and STD       
  /// (non-SPH bodies). bodytype can be used like a jave-style iterator to loop 
  /// over all bodytypes.                                                       
  // ///////////////////////////////////////////////////////////////////////////
  class bodytype {
  public:
    /// types of bodies
    enum bits {
      gas = 0, ///< gas: SPH particles
      std = 1  ///< standard: non-SPH bodies
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
    //--------------------------------------------------------------------------
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
    /// is bodytype supporting SPH?
    bool is_sph() const { return val == gas; }
    /// return the full set of fields bodies of this type may hold
    fieldset allows() const { return fieldset(BT_DATA[val]); }
    /// can bodies of this bodytype hold field \c f?
    bool allows(fieldbit f) const {
      return fieldset(BT_DATA[val]).contain(f);
    }
    /// return a name like "std" or "gas"
    const char* name() const {
      switch(val) {
      case gas: return "gas";
      default:  return "std";
      }
    }
  };
  falcON_TRAITS(bodytype,"bodytype");

  // ///////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // struct field_zero<type>                                                  //
  //                                                                          //
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T> struct field_zero {
    static void set_zero(T&x) { x = T(0); }
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
    //--------------------------------------------------------------------------
    field_zero<TYPE>::set_zero;
    //--------------------------------------------------------------------------
    static TYPE*array(void*D) {
      return static_cast<TYPE*>(D);}
    //--------------------------------------------------------------------------
    static const TYPE*c_array(const void*D) {
      return static_cast<const TYPE*>(D);}
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // struct field_traits<int>  int = fieldbit                                 //
  //                                                                          //
  // used in templated members of class bodies below.                         //
  // define all related types; aid nemo I/O                                   //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  /// \brief Template over body data field
  /// provides information about data type, and name etc or data field.
  template<int I> struct field_traits;
#  define DefFieldTraits(BIT,TYPE)				\
  template<> struct field_traits< BIT > :			\
  public field_type<TYPE> {					\
    typedef TYPE type;						\
    static const bool   is_sph = BIT >= BD_NOSPH;		\
    static const char   word() { return BD_SQUANT[BIT]; }	\
    static const char  *name() { return BD_QNAME[BIT]; }	\
  };
  //----------------------------------------------------------------------------
  DefFieldTraits( 0, real);                        // mass                      
  DefFieldTraits( 1, vect);                        // position                  
  DefFieldTraits( 2, vect);                        // velocity                  
  DefFieldTraits( 3, vect);                        // predicted velocity        
  DefFieldTraits( 4, real);                        // softening length          
  DefFieldTraits( 5, flag);                        // body flag                 
  DefFieldTraits( 6, int );                        // body key                  
  DefFieldTraits( 7, real);                        // time step                 
  DefFieldTraits( 8, real);                        // internal potential        
  DefFieldTraits( 9, real);                        // external potential        
  DefFieldTraits(10, vect);                        // acceleration              
  DefFieldTraits(11, vect);                        // jerk                      
  DefFieldTraits(12, real);                        // mass density              
  DefFieldTraits(13, real);                        // auxiliary scalar          
  DefFieldTraits(14, vect);                        // auxiliary vector          
  DefFieldTraits(15, indx);                        // time-step level           
  DefFieldTraits(16, unsigned);                    // # neighbours              
  DefFieldTraits(17, peanokey);                    // Peano-Hilbert key         
  DefFieldTraits(18, real);                        // SPH: smoothing length h   
  DefFieldTraits(19, unsigned);                    // SPH: # neighbours         
  DefFieldTraits(20, real);                        // SPH: internal energy U    
  DefFieldTraits(21, real);                        // SPH: predicted U_in       
  DefFieldTraits(22, real);                        // SPH: (dU/dt)_total        
  DefFieldTraits(23, real);                        // SPH: (dU/dt)_external     
  DefFieldTraits(24, real);                        // SPH: entropy              
  DefFieldTraits(25, real);                        // SPH: gas density          
  DefFieldTraits(26, real);                        // SPH: d(gas density)/dt    
  DefFieldTraits(27, real);                        // SPH: dh/dt                
  DefFieldTraits(28, real);                        // SPH: f_i                  
  DefFieldTraits(29, real);                        // SPH: sound speed          
#undef DefFieldTraits
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // macro for defining named functions                                       //
  //                                                                          //
  // named functions MUST not be defined in any other way to ensure that any  //
  // changes in the data type layout etc. are automatically reflected.        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
#define DEF_NAMED_NONSPH(MACRO)			\
  MACRO(fieldbit::m,mass);			\
  MACRO(fieldbit::x,pos);			\
  MACRO(fieldbit::v,vel);			\
  MACRO(fieldbit::w,vprd);			\
  MACRO(fieldbit::e,eps);			\
  MACRO(fieldbit::f,flg);			\
  MACRO(fieldbit::k,key);			\
  MACRO(fieldbit::t,tau);			\
  MACRO(fieldbit::p,pot);			\
  MACRO(fieldbit::q,pex);			\
  MACRO(fieldbit::a,acc);			\
  MACRO(fieldbit::j,jerk);			\
  MACRO(fieldbit::r,rho);			\
  MACRO(fieldbit::y,aux);			\
  MACRO(fieldbit::z,zet);			\
  MACRO(fieldbit::l,level);			\
  MACRO(fieldbit::n,num);			\
  MACRO(fieldbit::h,peano);
#define DEF_NAMED_SPH(MACRO)			\
  MACRO(fieldbit::H,size);			\
  MACRO(fieldbit::N,snum);			\
  MACRO(fieldbit::U,uin);			\
  MACRO(fieldbit::Y,uprd);			\
  MACRO(fieldbit::I,udot);			\
  MACRO(fieldbit::E,udex);			\
  MACRO(fieldbit::S,entr);			\
  MACRO(fieldbit::R,srho);			\
  MACRO(fieldbit::D,drho);			\
  MACRO(fieldbit::J,hdot);			\
  MACRO(fieldbit::F,fact);			\
  MACRO(fieldbit::C,csnd);
#define DEF_NAMED(MACRO)			\
  DEF_NAMED_NONSPH(MACRO)			\
  DEF_NAMED_SPH(MACRO)
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // templates for looping fields                                             //
  //                                                                          //
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
      LoopFields<LOOP, 0, BD_NQUANT>::loop(x);
    }
    template< typename T > static void const_loop(T const&x) {
      LoopFields<LOOP, 0, BD_NQUANT>::const_loop(x);
    }
  };
  //----------------------------------------------------------------------------
  template< template<int> class LOOP>
  struct LoopSPHFields {
    template< typename T > static void loop(T&x) {
      LoopFields<LOOP, BD_NOSPH, BD_NQUANT>::loop(x);
    }
    template< typename T > static void const_loop(T const&x) {
      LoopFields<LOOP, BD_NOSPH, BD_NQUANT>::const_loop(x);
    }
  };
  //----------------------------------------------------------------------------
  template< template<int> class LOOP>
  struct LoopSTDFields {
    template< typename T > static void loop(T&x) {
      LoopFields<LOOP, 0, BD_NOSPH>::loop(x);
    }
    template< typename T > static void const_loop(T const&x) {
      LoopFields<LOOP, 0, BD_NOSPH>::const_loop(x);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_fields_h
