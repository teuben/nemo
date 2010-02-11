// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   inc/body.h
///
/// \brief  contains declarations of class falcON::bodies,
///	    class falcON::snapshot and some useful macros;
///         currently problems with doxygen documentation.
///                                                                             
/// \author Walter Dehnen
/// \date   2000-2009
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2009 Walter Dehnen
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
#ifndef falcON_included_body_h
#define falcON_included_body_h 1

#ifndef falcON_included_fields_h
#  include <public/fields.h>
#endif
////////////////////////////////////////////////////////////////////////////////
/// All public code for this project is in namespace falcON
namespace falcON {
#ifdef falcON_NEMO
  class nemo_in;
  class snap_in;
  class data_in;
  class nemo_out;
  class snap_out;
  class data_out;
#endif
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // meta programmed sum over array                                             
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  template<int I, int N1> struct __Sum {
    template<typename T> static T sum(const T*X) {
      return X[I]+__Sum<I+1,N1>::sum(X); }
  };
  template<int I> struct __Sum<I,I> {
    template<typename T> static T sum(const T*X) { 
      return X[I]; }
  };
  template<int N, typename T> T sum(const T X[N]) {
    return N? __Sum<0,N-1>::sum(X) : T(0);
  }

  class ebodies;                                   // declared in forcesC.cc    
  class BodyFileter;                               // declared in bodyfunc.h    
  class ParallelSnapshot;                          // parallel/snapshot.h       
  class forces;                                    // forces.h                  
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  //  class falcON::bodies                                                      
  //                                                                            
  /// Holds and manages the body data.                                          
  ///                                                                           
  /// Body data are hold in a linked list of bodies::block s, which store the   
  /// data in a 'structure of array' (SoA) form. Sub-type bodies::iterator      
  /// provides access to body data (via its members and friends) mimicking an   
  /// 'array of structure' (AoS) layout. Sub-type bodies::index provides an     
  /// alternative data access via members of bodies. class falcON::bodies also  
  /// supports data I/O and adding/removing of data fields and bodies.          
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////

  class bodies {
    friend class ParallelSnapshot;
    //==========================================================================
    //                                                                          
    // static data                                                              
    //                                                                          
    //==========================================================================
  public:
    enum {
      DefSTD      = fieldset::gravity,
      DefSPH      = fieldset::sphdef,
      DefaultBits = DefSTD | DefSPH
    };
    class iterator;
    //==========================================================================
    //                                                                          
    // class bodies::block                                                      
    //                                                                          
    /// Holds a block of body data.                                             
    ///                                                                         
    /// This class is not for direct public usage. A block holds the data of    
    /// up to bodies::index::max_bodies bodies of the same bodytype.            
    ///                                                                         
    //==========================================================================
    class block {
      friend class bodies;
      friend class bodies::iterator;
      friend class ParallelSnapshot;
      //------------------------------------------------------------------------
      const bodytype     TYPE;                     // type of bodies hold       
      unsigned           NALL;                     // # data                    
      unsigned           NBOD;                     // # bodies hold             
      unsigned           NO;                       // this==bodies::BLOCK[NO]   
      unsigned           FIRST;                    // total index of first body
      unsigned           LOCALFIRST;               // local index of first body
      void              *DATA[BodyData::NQUANT];   // pointers to body data     
      block             *NEXT;                     // blocks: in linked list    
      const bodies      *BODS;                     // pointer back to bodies    
      //------------------------------------------------------------------------
      // data access                                                            
      //------------------------------------------------------------------------
    public:
      void      *data_void(fieldbit f)       { return DATA[value(f)]; }
      const void*data_void(fieldbit f) const { return DATA[value(f)]; }
      //------------------------------------------------------------------------
      template<int BIT>
      typename field_traits<BIT>::type* data() const {
	return field_traits<BIT>::array(DATA[BIT]);
      }
      //------------------------------------------------------------------------
      template<int BIT>
      const typename field_traits<BIT>::type* const_data() const {
	return field_traits<BIT>::c_array(DATA[BIT]);
      }
      //------------------------------------------------------------------------
      template<int BIT>
      typename field_traits<BIT>::type&datum(int i) const {
#if defined(DEBUG) || defined(EBUG)
	if(0 == DATA[BIT])
	  falcON_THROW("trying to access non-allocated \"%s\" data in bodies\n",
		       field_traits<BIT>::fullname());
	if(i < 0 || i >= NBOD)
	  falcON_THROW("index out of range in \"%s\" data access of bodies\n",
		       field_traits<BIT>::fullname());
#endif
	return data<BIT>()[i];
      }
      //------------------------------------------------------------------------
      template<int BIT>
      typename field_traits<BIT>::type const&const_datum(int i) const {
#if defined(DEBUG) || defined(EBUG)
	if(0 == DATA[BIT])
	  falcON_THROW("trying to access non-allocated \"%s\" data in bodies\n",
		       field_traits<BIT>::fullname());
	if(i < 0 || i >= NBOD)
	  falcON_THROW("index out of range in \"%s\" data access of bodies\n",
		       field_traits<BIT>::fullname());
#endif
	return const_data<BIT>()[i];
      }
      //------------------------------------------------------------------------
      flags      &flag(int i)       { return datum<fieldbit::f>(i); }
      flags const&flag(int i) const { return const_datum<fieldbit::f>(i); }
      //------------------------------------------------------------------------
      void set_first(unsigned f) { FIRST=f; LOCALFIRST=f; }
      void set_first(unsigned f, unsigned l) { FIRST=f; LOCALFIRST=l; }
      //------------------------------------------------------------------------
      bool has_field(fieldbit f) const { return DATA[value(f)] != 0; }
      //------------------------------------------------------------------------
      bool               is_sph     () const { return TYPE.is_sph(); }
      bool               is_empty   () const { return NBOD == 0; }
      unsigned     const&N_alloc    () const { return NALL; }
      unsigned     const&N_bodies   () const { return NBOD; }
      unsigned           N_free     () const { return NALL - NBOD; }
      unsigned     const&my_No      () const { return NO; }
      const bodies*const&my_bodies  () const { return BODS; }
      bodytype     const&type       () const { return TYPE; }
      unsigned     const&first      () const { return FIRST; }
      unsigned     const&localfirst () const { return LOCALFIRST; }
      unsigned           localend   () const { return LOCALFIRST+NBOD; }
      //------------------------------------------------------------------------
      block*const&next             () const { return NEXT; }
      block*      next_of_same_type() const {
	return (NEXT && NEXT->TYPE == TYPE)? NEXT : 0;
      }
      //------------------------------------------------------------------------
      static const block*next(const block*&B) {
	if(B) B=B->NEXT;
	return B;
      }
      static const block*next_of_same_type(const block*&B) {
	if(B && B->NEXT && B->NEXT->TYPE == B->TYPE) B = B->NEXT;
	return B;
      }
      //  31/07/08 WD: allowance for empty blocks
      /// replace \a B with the first non-empty block from \a B, if any.
      static const block*non_empty(const block*&B) {
	while(B && B->NBOD==0) B=B->NEXT;
	return B;
      }
      /// replace \a B with the next non-empty block from \a B, if any.
      static const block*non_empty_of_same_type(const block*&B) {
	while(B && B->NBOD==0) B=B->next_of_same_type();
	return B;
      }
      /// return the first non-empty block from this, if any.
      const block*first_non_empty() const {
	const block*B=this;
	return non_empty(B);
      }
      /// return the next non-empty block from this, if any.
      const block*next_non_empty() const {
	const block*B=NEXT;
	return non_empty(B);
      }
      /// return the first non-empty block from this of the same type, if any.
      const block*first_non_empty_of_same_type() const {
	const block*B=this;
	return non_empty_of_same_type(B);
      }
      /// return the next non-empty block of the same type, if any.
      const block*next_non_empty_of_same_type() const {
	const block*B=next_of_same_type();
	return non_empty_of_same_type(B);
      }
      //------------------------------------------------------------------------
      // flag manipulations                                                     
      //------------------------------------------------------------------------
      void reset_flags() const;
      void flag_all_as_active() const falcON_THROWING;
      //------------------------------------------------------------------------
      // construction & data allocation                                         
      //------------------------------------------------------------------------
    private:
      void set_data_void(fieldbit f, void*D) {
	if(0 != D  &&  0 != DATA[value(f)])
	  falcON_Warning("over writing pointer to allocated memory");
	DATA[value(f)] = D;
      }
      //------------------------------------------------------------------------
      // copy a single body within the block                                    
      fieldset copy_body  (                           // R: data copied         
			   unsigned,                  // I: from this index     
			   unsigned,                  // I: to   this index     
			   fieldset = fieldset::all); //[I: copy these data]    
      //------------------------------------------------------------------------
      // copy a number of bodies across blocks                                  
      fieldset copy_bodies(                           // R: data copied         
			   const block*,              // I: from this block     
			   unsigned,                  // I:  and this index     
			   unsigned,                  // I: to   our index here 
			   unsigned,                  // I: these many bodies   
			   fieldset = fieldset::all)  //[I: copy these data]    
	falcON_THROWING;
      //------------------------------------------------------------------------
      block();                                     // not implemented           
      block(block const&);                         // not implemented           
      //------------------------------------------------------------------------
      void clone(block*);
      //------------------------------------------------------------------------
    public:
      block(unsigned,                              // I: our No                 
	    unsigned,                              // I: # data to allocate     
	    unsigned,                              // I: # bodies <= ----       
	    unsigned,                              // I: first body's index     
	    bodytype,                              // I: type of bodies to hold 
	    fieldset,                              // I: data to allocate       
	    bodies *) falcON_THROWING;             // I: ^ back to my bodies    
      //------------------------------------------------------------------------
      void link      (block*n) {
	NEXT = n;
      }
      //------------------------------------------------------------------------
      void reset_data(fieldset ) const falcON_THROWING;
      void add_field (fieldbit ) falcON_THROWING;
      void del_field (fieldbit ) falcON_THROWING;
      void add_fields(fieldset ) falcON_THROWING;
      void del_fields(fieldset ) falcON_THROWING;
      void set_fields(fieldset ) falcON_THROWING;
      void remove    (unsigned&) falcON_THROWING;
      void swap_bytes(fieldbit ) falcON_THROWING;
      ~block() falcON_THROWING;
      //------------------------------------------------------------------------
      void skip(unsigned&, flags) const falcON_THROWING;
      //------------------------------------------------------------------------
      // copy up to NALL bodies                                                 
      // - we copy only bodies of the same type as hold here                    
      // - we only copy bodies whose flag matches last argument                 
      // - if the block copied is finished, we take its NEXT, starting at i=0   
      // - on return the block pointer is either NULL or together with i they   
      //   give the first body which was not copied because NALL was exceeded   
      // - NBOD is set to the bodies copied                                     
      fieldset copy(                               // R: data copied            
		    const block*&,                 // I: block[s] to copy from  
		    unsigned    &,                 // I: 1st body to copy       
		    fieldset = fieldset::all,      //[I: which data to copy]    
		    flags    = flags::empty)       //[I: flags for copying]     
	falcON_THROWING;
#ifdef falcON_NEMO
      //------------------------------------------------------------------------
      // NEMO I/O support                                                       
    protected:
      void read_data   (data_in &,unsigned,unsigned)          falcON_THROWING;
      void read_posvel (data_in &,unsigned,unsigned,fieldset) falcON_THROWING;
      void write_data  (data_out&,unsigned,unsigned) const    falcON_THROWING;
      void write_potpex(data_out&,unsigned,unsigned) const    falcON_THROWING;
#endif
      //------------------------------------------------------------------------
      // Gadget I/O support                                                     
#ifdef falcON_REAL_IS_FLOAT
      void read_Fortran (FortranIRec&, fieldbit, unsigned, unsigned, bool)
	falcON_THROWING;
      void write_Fortran(FortranORec&, fieldbit, unsigned, unsigned) const
	falcON_THROWING;
#endif
    };
    //==========================================================================
    //                                                                          
    // sub-type bodies::index                                                   
    //                                                                          
    /// An integer-like identifier for access to body data.                     
    /// class bodies and bodies::block are friends and can access the members   
    /// no() and in() to find the body data                                     
    ///                                                                         
    /// \note                                                                   
    /// An index is not preserved by methods that change the number of bodies   
    ///                                                                         
    //==========================================================================
  public:
    class index {
      friend class bodies;
      friend class bodies::block;
      friend class ParallelSnapshot;
      //------------------------------------------------------------------------
      /// useful constants, determine interpretation of index::I
      enum
      {
	block_bits = 8,                 ///< number of bits used for blocks No. 
	index_bits = 24,                ///< number of bits used for sub-index. 
	max_blocks = 1 << block_bits,   ///< maximum number of blocks.          
	max_bodies = 1 << index_bits,   ///< maximum number of bodies per block.
	index_mask = max_bodies - 1     ///< mask for extracting sub-index.     
      };
      //------------------------------------------------------------------------
      /// data are stored in 32bits.
      unsigned I;
    public:
      //------------------------------------------------------------------------
      /// return No of block, encoded in highest index::block_bits bits.
      unsigned no() const { return I >> index_bits; }
      //------------------------------------------------------------------------
      /// return sub-index, encoded in lowest index::index_bits bits.
      unsigned in() const { return I &  index_mask; }
      //------------------------------------------------------------------------
      /// return numerical value
      unsigned const&value() const { return I; }
      //------------------------------------------------------------------------
      /// default constructor: initialized
      index() {}
      //------------------------------------------------------------------------
      /// copy constructor.
      index(index const&i) : I(i.I) {}
      //------------------------------------------------------------------------
      /// constructor from block No (b) and index (n) within block.
      index(unsigned b, unsigned n) : I (n | b<<index_bits) {}
      //------------------------------------------------------------------------
      /// make a copy.
      index& operator= (index const&i) { I=i.I; return *this; }
      //------------------------------------------------------------------------
      /// output in format "no:index".
      friend std::ostream& operator<<(std::ostream&, const index&);
      //------------------------------------------------------------------------
      /// equality
      bool operator== (index const&i) const { return I == i.I; }
      //------------------------------------------------------------------------
      /// inequality
      bool operator!= (index const&i) const { return I != i.I; }
    };
    //==========================================================================
    //                                                                          
    // sub-type bodies::TimeSteps                                               
    //                                                                          
    //==========================================================================
    class TimeSteps {
      int      KMAX;                     ///< tau_max = 0.5^kmax
      unsigned NSTEPS, HIGHEST;          ///< # levels, highest level
      double  *TAU,*TAUQ,*TAUH;          ///< individual tau, tau^2, tau/2
    public:
      //------------------------------------------------------------------------
      TimeSteps(int km, unsigned ns)
	: KMAX   ( km ),
	  NSTEPS ( ns ),
	  HIGHEST( ns? ns-1 : 0 ),
	  TAU    ( NSTEPS? falcON_NEW(double, NSTEPS) : 0 ),
	  TAUQ   ( NSTEPS? falcON_NEW(double, NSTEPS) : 0 ),
	  TAUH   ( NSTEPS? falcON_NEW(double, NSTEPS) : 0 )
      {
	if(NSTEPS>0) {
	  TAU [0] = pow(0.5,KMAX);
	  TAUH[0] = 0.5*TAU[0];
	  TAUQ[0] = square(TAU[0]);
	  for(unsigned n=1; n!=NSTEPS; ++n) {
	    TAU [n] = TAUH[n-1];
	    TAUH[n] = 0.5*TAU[n];
	    TAUQ[n] = square(TAU[n]);
	  }
	} else 
	  falcON_Error("bodies::TimeSteps: ns=%d < 1\n",NSTEPS);
      }
      //------------------------------------------------------------------------
      ~TimeSteps() {
	if(TAU)  { falcON_DEL_A(TAU ); TAU =0; }
	if(TAUQ) { falcON_DEL_A(TAUQ); TAUQ=0; }
	if(TAUH) { falcON_DEL_A(TAUH); TAUH=0; }
      }
      //------------------------------------------------------------------------
      bool is_leap_frog() const { return NSTEPS == 1; }
      int highest_level() const { return HIGHEST; }
      const int&kmax() const { return KMAX; }
      unsigned const&Nsteps() const { return NSTEPS; }
      /// shortest time step
      const double&tau_min() const { return TAU[HIGHEST]; }
      /// longest time step
      const double&tau_max() const { return TAU[0]; }
      /// time step of level l
      double const&tau(int l) const { return TAU[l]; }
      /// half time step of level l
      double const&tauh(int l) const { return TAUH[l]; }
      /// time step squared, given of level l )
      double const&tauq(int l) const { return TAUQ[l]; }
      /// return time step of body 
      inline double const&tau(iterator const&) const;
      /// return time step squared of body 
      inline double const&tauq(iterator const&) const;
      /// return half time step of body 
      inline double const&tauh(iterator const&) const;
      /// return table of time steps
      const double*tau () const { return TAU; }
      /// return table of half time steps
      const double*tauq() const { return TAUQ; }
      /// return table of time steps squared
      const double*tauh() const { return TAUH; }
    };
    //==========================================================================
    //                                                                          
    // data of class bodies                                                     
    //                                                                          
    //==========================================================================
  private:
    unsigned         NALL[bodytype::NUM];          // # bodies allocated        
    unsigned         NBOD[bodytype::NUM];          // # bodies in use           
    unsigned         NNEW[bodytype::NUM];          // # bodies created new      
    unsigned         NDEL[bodytype::NUM];          // # bodies removed          
    unsigned         NTOT;                         // total # bodies in use     
    fieldset         BITS;                         // body data allocated       
    unsigned         NBLK;                         // # blocks in use           
    block           *BLOCK[index::max_blocks];     // table: blocks             
    // NOTE BLOCK[] may be filled randomly
    block           *TYPES[bodytype::NUM];         // table: bodies per bodytype
    block           *FIRST;                        // first block of bodies     
    const TimeSteps *TSTEPS;                       // time steps                
    mutable bool     SRCC;                         // source data changed?      
    mutable bool     SPHC;                         // SPH data changed?         
    const bool       C_FORTRAN;                    // we are used for C/FORTRAN 
    mutable
    const forces    *FORCES;                       // forces, if any            
    //--------------------------------------------------------------------------
  public:
    void set_forces(const forces*f) const { FORCES=f; }
    const forces*const&Forces() const { return FORCES; }
    block *const&first_block() const { return FIRST; }
    //==========================================================================
    // \name Functions providing information about the number of bodies
    /// # bodies allocated for a given bodytype.
    unsigned const&N_alloc (bodytype t) const
    { return NALL[int(t)]; }
    /// total # bodies allocated.
    unsigned       N_alloc ()           const
    { return sum<bodytype::NUM>(NALL); }
    /// array with # bodies in use per bodytype.
    const unsigned*N_bodies_per_type()  const
    { return NBOD; }
    /// # bodies in use for a given bodytype.
    unsigned const&N_bodies(bodytype t) const
    { return NBOD[int(t)]; }
    /// # bodies in use for given set of bodytypes.
    unsigned N_bodies(bodytypes T) const {
      unsigned N(0u);
      for(bodytype t; t; ++t)
	if(T.contain(t))
	  N += NBOD[int(t)];
      return N;
    }
    /// total # bodies in use.
    unsigned const&N_bodies()           const { return NTOT; }
    /// # bodies not flagged to be ignored
    unsigned       N_subset()           const;
    /// # SPH bodies in use.
    unsigned const&N_sph   ()           const { return NBOD[bodytype::gas]; }
    /// # standard bodies in use.
    unsigned const&N_std   ()           const { return NBOD[bodytype::std]; }
    /// # standard bodies in use.
    unsigned const&N_sink  ()           const { return NBOD[bodytype::sink]; }
    /// # bodies allocated but not used for a given bodytype.
    unsigned       N_free  (bodytype t) const { return N_alloc(t)-N_bodies(t); }
    /// total # bodies allocated but not used.
    unsigned       N_free  ()           const { return N_alloc() -N_bodies(); }
    //--------------------------------------------------------------------------
    bool const&srce_data_changed() const { return SRCC; }
    bool const&sph_data_changed() const { return SPHC; }
    //==========================================================================
    // \name Functions providing information about body data hold
    /// sum of all body data hold.
    fieldset const&all_data  ()           const { return  BITS; }
    /// do we have none of these data?
    bool           have_not  (fieldset b) const { return !BITS.intersect(b); }
    /// do we have all of these data?
    bool           have_all  (fieldset b) const { return  BITS.contain(b); }
    /// do we have some of these data?
    bool           have_some (fieldset b) const { return  BITS.intersect(b); }
    /// do we have this datum?
    bool           have      (fieldbit f) const { return  BITS.contain(f); }
    /// do we hold the body datum indicated?
#define HaveField(BIT,NAME)						\
    bool have_##NAME() const { return BITS.contain(BIT); }
    DEF_NAMED(HaveField)
#undef HaveField
    //==========================================================================
    // \name non-const data access using bodies::index
    /// return lvalue (can be changed)
    template<int BIT>
    typename field_traits<BIT>::type&datum(index i) const {
      return BLOCK[i.no()]->datum<BIT>(i.in());
    }
#define NonConstAccess(BIT,NAME)				\
    field_traits<BIT>::type      &    NAME    (index i) const {	\
      return datum<BIT>(i);					\
    }
    DEF_NAMED(NonConstAccess)
#undef NonConstAccess
    //==========================================================================
    // \name const data access using bodies::index
    /// return const datum
    template<int BIT>
    typename field_traits<BIT>::type const&const_datum(index i) const {
      return BLOCK[i.no()]->const_datum<BIT>(i.in());
    }
#define ConstAccess(BIT,NAME)					\
    field_traits<BIT>::type const&c_##NAME    (index i) const {	\
      return const_datum<BIT>(i);				\
    }
    DEF_NAMED(ConstAccess)
#undef ConstAccess
    //==========================================================================
    // \name further methods using bodies::index
    /// return true if the named or indicated datum is hold
    bool has_field(index i, fieldbit f) const {
      return BLOCK[i.no()] && BLOCK[i.no()]->has_field(f);
    }
#define HasDatum(Bit,Name)					\
    bool has_##Name(index i) const {				\
      return BLOCK[i.no()] && BLOCK[i.no()]->has_field(Bit);	\
    }
    DEF_NAMED(HasDatum)
#undef HasDatum
    /// is the index a valid one?
    bool is_valid(index i) const {
      return 0 != BLOCK[i.no()]  &&  i.in() < BLOCK[i.no()]->N_bodies();
    }
    /// index to the first of all bodies
    index first() const { return index(FIRST->my_No(),0); }
    /// running index of body (between 0 and N-1).
    /// \note in MPI parallel code this is a \b global index
    unsigned bodyindex(index i) const {
      return BLOCK[i.no()]->first() + i.in();
    }
    /// running index of body (between 0 and N-1).
    /// \note in MPI parallel code this is a \b local index, starting at 0
    unsigned localindex(index i) const {
      return BLOCK[i.no()]->localfirst() + i.in();
    }
    /// comparison between indices, also works in parallel across nodes
    bool is_less(index a, index b) const {
      return  (a.no() == b.no() &&  a.in() < b.in())
	||    BLOCK[a.no()]->FIRST < BLOCK[b.no()]->FIRST;
    }
    /// bodytype of body.
    bodytype const&type(index i) const {
      return BLOCK[i.no()]->type();
    }
    //==========================================================================
    // \name support for time steps (including access via bodies::index)
    /// do we have time step data?
    bool have_steps() const { return TSTEPS!=0; }
    /// let the user provide time step data
    void set_steps(const TimeSteps*T) { TSTEPS = T; }
    /// give TimeSteps 
    const TimeSteps*steps() const { return TSTEPS; }
    /// time step given a time step level
    double const&tau (indx l) const { return TSTEPS->tau(l); }
    /// time step squared given a time step level
    double const&tauq(indx l) const { return TSTEPS->tauq(l); }
    /// half time step given a time step level
    double const&tauh(indx l) const { return TSTEPS->tauh(l); }
    /// time step of body given its index
    /// \note if leapfrog: return global time step
    double const&tau(index i) const {
      return have_level()? tau (level(i)) : tau(0);
    }
    /// time step squared of body given its index
    /// \note if leapfrog: return global time step squared
    double const&tauq(index i) const {
      return have_level()? tauq(level(i)) : tauq(0);
    }
    /// half time step of body given its index
    /// \note if leapfrog: return half global time step
    double const&tauh(index i) const {
      return have_level()? tauh(level(i)) : tauh(0);
    }
    //==========================================================================
    //                                                                          
    // sub-type bodies::iterator                                                
    //                                                                          
    /// C++ style iterator over bodies, provides access to body data.           
    ///                                                                         
    /// Use iterator for easiest access to data of bodies; used a lot in falcON;
    /// Obtain object of type iterator from member of bodies;                   
    /// non-const data access via members, const data access via friends.       
    ///                                                                         
    /// The following two examples assign the position of the body referred to  
    /// iterator \c i1 to the body referred to by iterator \c i2                
    /// \code                                                                   
    ///    i2.datum<fieldbit::x>() = const_dat<fieldbit::x>(i1);                
    ///    i2.pos() = pos(i1);                                                  
    /// \endcode                                                                
    /// The first method allows templated manipulations, while the second is    
    /// more human readable.                                                    
    ///                                                                         
    //==========================================================================
    class iterator {
      friend class bodies;
      const block* B;                              ///< our block
      unsigned     K;                              ///< our index within block
      //------------------------------------------------------------------------
      /// construction from pointer to block and index within block
      explicit 
      iterator(const block*b,
	       unsigned    k=0) : B(b), K(k) {}
      //------------------------------------------------------------------------
      //  31/07/08 WD: allow for empty blocks
      void next_block() {
	block::non_empty(B=B->next());             //   get next non-empty block
	K = 0;                                     //   set index to 0          
      }
      //------------------------------------------------------------------------
    public:
#if defined(DEBUG) || defined(EBUG)
# define CheckInvalid(NAME)					\
	if(!is_valid())						\
	  falcON_THROW("body::%s() called on invalid body",NAME);
#else
# define CheckInvalid(NAME)
#endif

      /// return pointer to my bodies::block
      const block *const&my_block () const {
	return B;
      }
      /// return index within block
      unsigned const&my_subindex() const {
	return K;
      }
      /// return const pointer to my falcON::bodies
      const bodies*const&my_bodies() const {
	CheckInvalid("my_bodies");
	return B->my_bodies();
      }
      /// return number of my bodies::block
      unsigned const&my_block_No () const {
	CheckInvalid("my_block");
	return B->my_No();
      }
      /// return running body index (between 0, N-1)
      /// \note  in parallel code, the index is global and hence unique
      unsigned my_index() const {
	CheckInvalid("my_index");
	return B->first() + K;
      }
      /// return running body index (between 0, N-1)
      /// \note in parallel code, the index is local only
      unsigned my_localindex() const {
	CheckInvalid("my_index");
	return B->localfirst() + K;
      }
      /// friend return const pointer to falcON::bodies of iterator
      friend const bodies*const&bodies_of(iterator const&);
      /// friend return pointer to bodies::block of iterator
      friend const block *const&block_of(iterator const&);
      /// friend return number of bodies::block of iterator
      friend unsigned const&block_No(iterator const&);
      /// friend returning sub-index within block
      friend unsigned const&subindex(iterator const&);
      /// friend returning global running body index (between 0, N_global-1)
      friend unsigned bodyindex(iterator const&);
      /// friend returning local running body index (between 0, N_local-1)
      friend unsigned localindex(iterator const&);
      /// conversion to bodies::index
      operator index() const {
	CheckInvalid("operator index");
	return index(B->my_No(), K);
      }
      /// are we refering to a valid body?
      bool is_valid() const {
	return B != 0;
      }
      /// conversion to bool: are we referring to a valid body?
      operator bool() const {
	return is_valid();
      }
      /// are we identical to another iterator?
      bool operator== (iterator const&i) const {
	return K==i.K && B==i.B;
      }
      /// are we different from another iterator?
      bool operator!= (iterator const&i) const {
	return K!=i.K || B!=i.B;
      }
      /// are we before another iterator? (assuming the same bodies)
      bool operator<  (iterator const&i) const {
	CheckInvalid("operator<");
	return B==i.B? K< i.K : B->first()<i.B->first();
      }
      /// are we before or equal to another iterator? (assuming the same bodies)
      bool operator<= (iterator const&i) const {
	CheckInvalid("operator<=");
	return B==i.B? K<=i.K : B->first()<i.B->first();
      }
      /// are we after another iterator? (assuming the same bodies)
      bool operator>  (iterator const&i) const {
	CheckInvalid("operator>");
	return B==i.B? K> i.K : B->first()>i.B->first();
      }
      /// are we after or equal to another iterator? (assuming the same bodies)
      bool operator>= (iterator const&i) const {
	CheckInvalid("operator>=");
	return B==i.B? K>=i.K : B->first()>i.B->first();
      }
      //------------------------------------------------------------------------
#define HasDatum(Bit,Name) friend bool has_##Name(iterator const&);
      /// has our body the datum indicated?
      friend bool has_field(iterator const&, fieldbit);
      DEF_NAMED(HasDatum)
#undef HasDatum
      //------------------------------------------------------------------------
      /// pre-fix: get next body (note: there is no post-fix operator++)
      iterator& operator++() {
	CheckInvalid("operator++");
	if(++K == B->N_bodies()) next_block();
	return *this;
      }
      /// jump k bodies forward (or to last valid body)
      iterator& operator+=(unsigned k) {
	while(is_valid() && k) {
	  unsigned s = min(B->N_bodies()-K, k);
	  K += s;
	  k -= s;
	  if(K >= B->N_bodies()) next_block();
	}
	return *this;
      }
      /// return *this += k
      iterator operator+ (unsigned k) const {
	return iterator(*this) += k;
      }
      //------------------------------------------------------------------------
      /// default constructor makes an invalid body
      iterator() : B(0), K(0) {}
      /// copy constructor
      iterator(iterator const&i) : B(i.B), K(i.K) {}
      /// iterator offset by \e offset from \e i
      iterator(iterator const&i,
	       unsigned       offset) : B(i.B), K(i.K) {
	*this += offset;
      }
      /// assignment: make exact copy
      iterator& operator=(iterator const&i) {
	B = i.B;
	K = i.K;
	return *this;
      }
      //------------------------------------------------------------------------
      /// return lvalue (can be changed)
      template<int BIT>
      typename field_traits<BIT>::type &datum() {
#if defined(DEBUG) || defined(EBUG)
	if(!is_valid())
	  falcON_THROW("body::datum<%c>() called on invalid body",
		       field_traits<BIT>::word());
#endif
	return B-> template datum<BIT>(K);
      }
      //------------------------------------------------------------------------
#define NonConstAccess(Bit,Name)			\
      field_traits<Bit>::type&Name() {			\
	CheckInvalid(field_traits<Bit>::funcname());	\
        return B-> datum<Bit>(K);			\
      }
      DEF_NAMED(NonConstAccess)
#undef NonConstAccess
      //------------------------------------------------------------------------
      /// return const datum
      template<int BIT>
      typename field_traits<BIT>::type const&const_dat() const {
#if defined(DEBUG) || defined(EBUG)
	if(! is_valid() )
	  falcON_THROW("body::const_dat<%c>() called on invalid body",
		       field_traits<BIT>::word());
#endif
	return B-> template const_datum<BIT>(K);
      }
#define ConstAccess(Bit,Name)						\
      friend field_traits<Bit>::type const&Name(iterator const&);
      DEF_NAMED(ConstAccess)
#undef ConstAccess
      /// return angular momentum
      friend vect angmom(iterator const&);
      /// return type of body
      friend bodytype const&type(iterator const&);
      /// return time step of body 
      friend double const&tau(iterator const&);
      /// return time step squared of body 
      friend double const&tauq(iterator const&);
      /// return half time step of body 
      friend double const&tauh(iterator const&);
      //------------------------------------------------------------------------
      /// conversion to flags
      operator const flags&() const {
	return const_dat<fieldbit::f>();
      }
      /// is the single flag f set?
      bool flag_is_set(flags::single f) const {
	return const_dat<fieldbit::f>().is_set(f);
      }
      /// are all flags in f set?
      bool flags_are_set(flags const&f) const {
	return const_dat<fieldbit::f>().are_set(f);
      }
      /// flag this body as active
      void flag_as_active    () { flag().add(flags::active); }
      /// flag this body as inactive
      void unflag_active     () { flag().un_set(flags::active); }
      /// flag this body as not new
      void unflag_new        () { flag().un_set(flags::newbody); }
      void flag_as_sticky    () { flag().add(flags::sticky); }
      void unflag_sticky     () { flag().un_set(flags::sticky); }
      /// flag this body for removal by bodies::remove()
      void flag_for_removal  () { flag().add(flags::remove); }
      /// flag this body for usage with longer time step
      void allow_longer      () { flag().un_set(flags::not_longer); }
      /// flag this body not to use longer time step
      void forbid_longer     () { flag().add(flags::not_longer); }
      /// flag this body for usage with shorter time step
      void allow_shorter     () { flag().un_set(flags::not_shorter); }
      /// flag this body not to use shorter time step
      void forbid_shorter    () { flag().add(flags::not_shorter); }
      /// flag this body as being marked
      void mark              () { flag().add(flags::marked); }
      /// flag this body as not being marked
      void unmark            () { flag().un_set(flags::marked); }
      /// flag this body as being in subset
      void into_subset       () { flag().un_set(flags::ignore); }
      /// flag this body as not being in subset
      void outof_subset      () { flag().add(flags::ignore); }
      /// flag this body as being specially SPH marked
      void mark_SPH_special  () { flag().add(flags::sph_special); }
      /// flag this body as not being specially SPH marked
      void unmark_SPH_special() { flag().un_set(flags::sph_special); }
      //------------------------------------------------------------------------
      /// is body active?
      bool is_active     () const {
	return falcON::is_active(const_dat<fieldbit::f>()); }
      /// is body to be removed?
      bool to_remove     () const {
	return falcON::to_remove(const_dat<fieldbit::f>()); }
      /// is body SPH particle?
      bool is_sph        () const {
// 	return falcON::is_sph(const_dat<fieldbit::f>()); }
	return B->type().is_sph(); }
      /// is body STD particle?
      bool is_std        () const {
	return B->type().is_std(); }
      /// is body SINK particle?
      bool is_sink       () const {
	return B->type().is_sink(); }
      /// is body sticky particle?
      bool is_sticky     () const {
	return falcON::is_sticky(const_dat<fieldbit::f>()); }
      /// is body new?
      bool is_new        () const {
	return falcON::is_new(const_dat<fieldbit::f>()); }
      /// is this body flagged as being marked?
      bool is_marked     () const {
	return falcON::is_marked(const_dat<fieldbit::f>()); }
      /// is this body not flagged as being ignored?
      bool in_subset     () const {
	return falcON::in_subset(const_dat<fieldbit::f>()); }
      /// has this body non-zero mass?
      bool is_source     () const;
      /// is this body allowed to use a longer time step?
      bool may_go_longer () const {
	return !flag_is_set(flags::not_longer);}
      /// is this body allowed to use a shorter time step?
      bool may_go_shorter() const {
	return !flag_is_set(flags::not_shorter); }
      //------------------------------------------------------------------------
      /// increment until next body in subset
      iterator& next_in_subset();
      /// increment until next active body
      iterator& next_active();
      //------------------------------------------------------------------------
      /// friend: is body active?
      friend bool is_active(const iterator&);
      /// is body to be removed?
      friend bool to_remove(const iterator&);
      /// is body SPH particle?
      friend bool is_sph(const iterator&);
      /// is body STD particle?
      friend bool is_std(const iterator&);
      /// is body SINK particle?
      friend bool is_sink(const iterator&);
      /// is body sticky particle?
      friend bool is_sticky(const iterator&);
      /// is body new?
      friend bool is_new(const iterator&);
      /// friend: is this body flagged as being marked?
      friend bool is_marked(const iterator&);
      /// is this body not flagged as being ignored?
      friend bool in_subset(const iterator&);
      /// friend: has this body mass?
      friend bool is_source(const iterator&);
      /// friend: is this body allowed to use a longer time step?
      friend bool may_go_longer(const iterator&);
      /// friend: is this body allowed to use a shorter time step?
      friend bool may_go_shorter(const iterator&);
      //------------------------------------------------------------------------
      // I/O                                                                    
      /// formatted output: write bodyindex
      friend std::ostream& operator<<(std::ostream&, const iterator&);
      //------------------------------------------------------------------------
#ifdef falcON_REAL_IS_FLOAT
      iterator& read_Fortran (FortranIRec&, fieldbit, unsigned, bool=0)
	falcON_THROWING;
      iterator& write_Fortran(FortranORec&, fieldbit, unsigned) falcON_THROWING;
#endif
#ifdef falcON_NEMO
      iterator& read_data (data_in &, unsigned =0) falcON_THROWING;
      iterator& write_data(data_out&, unsigned =0) falcON_THROWING;
    private:
      iterator& write_potpex(data_out&, unsigned =0) falcON_THROWING;
      iterator& read_posvel (data_in &, fieldset, unsigned =0) falcON_THROWING;
#endif
      //------------------------------------------------------------------------
    };// class iterator
#undef CheckInvalid
    //==========================================================================
    // \name generate iterators                                                
    // \detail use these methods to obtain begin and end iterators when
    // looping over all or a subset of all bodies
    /// begin of all bodies
    //  31/07/08 WD: allow for empty FIRST block
    iterator begin_all_bodies() const {
      return iterator(FIRST->first_non_empty());
    }
    /// end of all bodies == invalid iterator
    iterator end_all_bodies  () const { return iterator(0); }
    /// begin of bodies of given bodytype, if any.
    /// if no bodies of given type exist, invalid iterator is returned
    //  31/07/08 WD: allow for empty TYPES[t] block
    iterator begin_typed_bodies(bodytype t) const {
      return iterator(TYPES[t]->first_non_empty());
    }
    /// begin of active bodies
    //  31/07/08 WD: allow for empty FIRST block
    iterator begin_active_bodies() const {
      iterator B=begin_all_bodies();
      while(B && !is_active(B)) ++B;
      return B;
    }
    /// begin of bodies in subset
    //  04/11/08 WD: generated; required for macro LoopSubsetBodies
    iterator begin_subset_bodies() const {
      iterator B=begin_all_bodies();
      if(have_flag()) while(B && !in_subset(B)) ++B;
      return B;
    }
    /// end of bodies of given bodytype, if any.
    /// if no bodies of given type exist, invalid iterator is returned
    /// otherwise the next other type body, or an invalid iterator
    //  31/07/08 WD: allow for empty TYPES[t] block
    iterator end_typed_bodies(bodytype t) const {
      if(TYPES[t]==0) return iterator(0); 
      for(++t; t; ++t)
	if(TYPES[t]) return iterator(TYPES[t]->first_non_empty());
      return iterator(0);
    }
    /// body of a given bodies::index
    iterator bodyIn(index i) const {
      const block*p = BLOCK[i.no()];
      if(p==0) return iterator(0);
      const unsigned in = i.in();
      return in < p->N_bodies()? iterator(p,in) : iterator(0);
    }
    /// body of a given local running index in [0, N_bodies()[
    /// \note this makes sense only locally
    iterator bodyNo(unsigned i) const {
      const block*p = FIRST;
      while(p && i >= p->localend()) p = p->next();
      return p? iterator(p, i - p->localfirst()) : iterator(0);
    }
    /// an invalid body (useful as default argument)
    static iterator bodyNil() { return iterator(0); }
    //==========================================================================
    // \name construction, destruction, and management of bodies and body data
    ///
    /// default constructor
    explicit 
    bodies(fieldset Bd=fieldset(DefaultBits)) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// constructor
    ///
    /// \param[in] N  array with number of bodies per bodytype
    /// \param[in] Bd body data fields to allocate (default: mxvapfHRVJFC)
    explicit 
    bodies(const unsigned N[bodytype::NUM], fieldset Bd=fieldset(DefaultBits))
      falcON_THROWING;
    //--------------------------------------------------------------------------
    /// copy constructor
    ///
    /// \param[in] B  bodies to copy
    /// \param[in] Bd body data fields to copy
    /// \param[in] C  if non-zero: copy only bodies which have this flag set
    /// \param[in] T  only copy bodies with type contained in T
    bodies(bodies const&B, fieldset Bd=fieldset::all, flags C=flags::empty,
	   bodytypes T=bodytypes::all) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// destructor: delete all data
    ~bodies() falcON_THROWING;
    //--------------------------------------------------------------------------
    /// Add an unsupported body data field
    ///
    /// If the body data field \a f is not yet supported, then we allocate body
    /// data indicated by \a f for all bodies which can have it (i.e. a SPH data
    /// field is only allocated for SPH bodies).
    void add_field(fieldbit f) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// Add unsupported body data fields
    ///
    /// Those body data fields in \a b which are not yet supported are allocated
    /// for all bodies which can have it (i.e. SPH data fields are only
    /// allocated for SPH bodies)
    void add_fields(fieldset b) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// Delete a supported body data field
    ///
    /// If the body data field \a f is presently allocated, it is deleted
    void del_field(fieldbit f) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// Delete supported body data fields
    ///
    /// Those body data fields in \a b which are presently allocated will be
    /// deleted. Fields not allocated are not affected
    void del_fields(fieldset b) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// Resets N, data: equivalent to destructor followed by constructor
    ///
    /// \param[in] N  array with number of bodies per bodytype
    /// \param[in] Bd body data fields to allocate (default: mxvapfHRVJFC)
    void reset(const unsigned N[bodytype::NUM],
	       fieldset Bd= fieldset(DefaultBits)) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// Resets N, keeps data the same (N[] = bodies per bodytype)
    ///
    /// \param[in] N number of bodies per bodytype
    void resetN(const unsigned N[bodytype::NUM]) falcON_THROWING
    { reset(N,BITS); }
    //--------------------------------------------------------------------------
    /// Reset some body data to zero
    ///
    /// Body data of fields in \a f, which are allocated, are set to zero
    /// \param[in] f body data fields to be reset
    void reset_data(fieldset f) const falcON_THROWING
    { for(const block*p=FIRST; p; p=p->next()) p->reset_data(f); }
    //--------------------------------------------------------------------------
  protected:
    /// swap the bytes (little-endian <-> big-endian)
    ///
    /// \param[in] b body data field to swap bytes for
    void swap_bytes(fieldbit b) falcON_THROWING;
  public:
    //==========================================================================
    // \name methods using another set of bodies
#if(0)
    ///
    /// Copy: replace our data (if any) with those of other bodies (\b not \b 
    /// yet \b implemented)
    ///
    /// We only copy data specified by 2nd argument and copy only those bodies
    /// whose flag matches 3rd argument. Equivalent to (but more efficient than)
    /// destruction followed by constructor 2.
    ///
    /// \note \b NOT \b YET \b IMPLEMENTED
    ///
    /// \param[in] B  bodies to be copied
    /// \param[in] Bd body data to be copied
    /// \param[in] f  flag for bodies to be copied
    void copy(bodies const&B, fieldset Bd=fieldset::all, int f=0)
      falcON_THROWING;
    //--------------------------------------------------------------------------
    /// Add: add other bodies (\b not \b yet \b implemented)
    /// 
    /// We only copy data specified by 2nd argument and copy only those bodies
    /// whose flag matches 3rd argument. Previously existing bodies are
    /// preserved.
    ///
    /// \note \b NOT \b YET \b IMPLEMENTED
    ///
    /// \param[in] B  bodies to be added
    /// \param[in] Bd body data to be copied
    /// \param[in] f  flag for bodies to be copied
    void add(bodies const&B,  fieldset Bd=fieldset::all, int f=0)
      falcON_THROWING;
#endif
    //--------------------------------------------------------------------------
    /// Merge: merge two sets of bodies
    ///
    /// The bodies merged in will be emptied. Body data not previously supported
    /// by this will be deleted before merging. Body data previously supported
    /// by this but not by \a B will be allocated but remain unitialized.
    ///
    /// \param[in] B bodies to merge with, will be empty after call
    void merge(bodies&B) falcON_THROWING;
    //==========================================================================
    // \name removal and creation of bodies
    /// Remove bodies which have been flagged for removal, eg by
    /// iterator::flag_for_removal().
    ///
    /// We will fill in the gaps left by the removed bodies by bodies taken
    /// from the end. This implies that the order of bodies is altered (however,
    /// bodies of the same bodytype are kept together) and that the number of
    /// bodies used is less than that allocated.
    /// \note alters the body order, affecting iterator::bodyindex()
    void remove();
    //--------------------------------------------------------------------------
    /// Remove bodies of type t which have been flagged for removal, eg by
    /// iterator::flag_for_removal().
    ///
    /// We will fill in the gaps left by the removed bodies by bodies taken
    /// from the end. This implies that the order of bodies is altered and that
    /// the number of bodies used is less than that allocated. 
    /// \note alters to body order.
    void remove(bodytype t);
    //--------------------------------------------------------------------------
    /// activate a body of type \a t.
    /// If no body is available for activation, we allocate at least \a Na
    /// bodies (though only one will be activated).
    ///
    /// \return a valid bodies::iterator to a body of type \a t
    /// \param[in] t   bodytype of body requested
    /// \param[in] Na  if no body available, allocate at least this many
    iterator new_body(bodytype t, unsigned Na=0) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// activate \a Nb bodies of type \a t.
    /// If not enough bodies are available for activation, we allocate at
    /// least \a Na bodies (though only \a Nb will be activated).
    ///
    /// \return a valid bodies::iterator to the first of \a Nb new bodies
    /// \param[in] Nb  # bodies to activate
    /// \param[in] t   bodytype of body requested
    /// \param[in] Na  # bodies to allocate (but at least Nb)
    iterator new_bodies(unsigned Nb, bodytype t, unsigned Na=0) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// \brief # bodies of type \a t created by new_body() or new_bodies()
    /// since last call of reset_Nnew(), if any.
    unsigned const&N_new(bodytype t) const {
      return NNEW[t];
    }
    //--------------------------------------------------------------------------
    /// \brief # bodies of types \a t created by new_body() or new_bodies()
    /// since last call of reset_Nnew(), if any.
    unsigned N_new(bodytypes t) const {
      unsigned N(0u);
      for(bodytype b; b; ++b)
	if(t.contain(b))
	  N += NNEW[int(t)];
      return N;
    }
    //--------------------------------------------------------------------------
    /// \brief # bodies of all types created by new_body() or new_bodies()
    /// since last call of reset_Nnew(), if any.
    unsigned N_new() const
    { return sum<bodytype::NUM>(NNEW); }
    //--------------------------------------------------------------------------
    /// resets the counters of bodies created by new_body()
    void reset_Nnew() {
      for(bodytype t; t; ++t)
	NNEW[t] = 0u;
    }
    //--------------------------------------------------------------------------
    /// \brief returns the number of bodies of type \a t removed by remove() 
    /// since last call of reset_Ndel(), if any.
    unsigned const&N_del(bodytype t) const {
      return NDEL[t];
    }
    //--------------------------------------------------------------------------
    /// \brief returns the number of bodies of types \a t created by new_body()
    /// since last call of reset_Nnew(), if any.
    unsigned N_del(bodytypes t) const {
      unsigned N(0u);
      for(bodytype b; b; ++b)
	if(t.contain(b))
	  N += NDEL[int(t)];
      return N;
    }
    //--------------------------------------------------------------------------
    /// \brief returns the number of bodies of any type removed by remove()
    /// since last call of reset_Ndel(), if any.
    unsigned N_del() const
    { return sum<bodytype::NUM>(NDEL); }
    //--------------------------------------------------------------------------
    /// resets the counters of bodies removed by remove()
    void reset_Ndel() {
      for(bodytype t; t; ++t)
	NDEL[t] = 0u;
    }
    //--------------------------------------------------------------------------
    /// \brief ensure bodies of type \a t are contiguous, i.e. free bodies, if
    /// any, are beyond the last activated body.
    ///
    /// This is done without memory allocation by copying bodies across blocks.
    void joinup(bodytype t) falcON_THROWING; 
    //--------------------------------------------------------------------------
    /// shrink allocation non-agressively.
    /// After joining up bodies of all types, we de-allocate blocks with no
    /// active bodies in them. Thus, at most one block per type will have free
    /// space.
    void shrink() falcON_THROWING {
      for(bodytype t; t; ++t) joinup(t);
      remove_empty_blocks();
    }
    //==========================================================================
    // \name I/O to file in various formats
#ifdef falcON_NEMO
    //--------------------------------------------------------------------------
    /// Read from a single NEMO snapshot
    ///
    /// We start reading into the body \a start. There must be enough space
    /// allocated to accomodate input (we do not reset the number of bodies).
    /// Body data fields are added, if required, by calls to
    /// bodies::add_field().
    ///
    /// \note this routine is used from snapshot::read_nemo()
    ///
    /// \return           body data fields that have been read
    /// \param[in] Si     snapshot input to read from
    /// \param[in] Bd     body data fields to be read
    /// \param[in] start  start position for reading
    /// \param[in] Nread  read this many (default: as many as possible)
    /// \param[in] warn   issue falcON::warning() for missing data
  protected:
    fieldset read_snapshot(snap_in  const&Si,
			   fieldset       Bd,
			   iterator const&start,
			   unsigned       Nread = 0,
			   bool           warn  = 1) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// Write to a single NEMO snapshot
    ///
    /// We start writing from the body \a start. The snapshot output must
    /// accomodate the required number of bodies to write.
    ///
    /// \note this routine is used from snapshot::write_nemo()
    ///
    /// \param[in] So     snapshot output to write to
    /// \param[in] Bd     body data fields to be written
    /// \param[in] start  start position for writing
    /// \param[in] Nwrite write this many (default: all in output)
    void write_snapshot(snap_out const&So,
			fieldset       Bd,
			iterator const&start,
			unsigned       Nwrite= 0) const falcON_THROWING;
#endif // falcON_NEMO
  public:
    //--------------------------------------------------------------------------
    /// Read simple ascii formatted input                                       
    ///
    /// Read data from an ASCII file. Any lines starting with '\#' are ignored.
    /// the number of lines must match (or exceed) \a N. The order and number
    /// of entries in each line must match \a Bd and \a Nd, respectively. Data
    /// are assumed to have \a N[bodytype::sink] lines for sink bodies, followed
    /// by \a N[bodytype::gas] lines for SPH bodies, and \a N[bodytype::std]
    /// lines of standard bodies.
    ///
    /// \param[in] In (input) input stream to read from
    /// \param[in] Bd (input) array of body data fields to read in given order
    /// \param[in] Nd (input) number of entries in that array
    /// \param[in] N  (input) number of lines to read per body type
    void read_simple_ascii(std::istream  &In,
			   const fieldbit*Bd,			   
			   unsigned       Nd,
			   const unsigned N[bodytype::NUM]);
#ifdef falcON_REAL_IS_FLOAT
    //--------------------------------------------------------------------------
    /// Reads a single snapshot from file(s) written in gadget2 data format 1
    ///
    /// Gadget allows a single snapshot to be distributed over more than one
    /// file. In this case the file names are assumed to be derived from the one
    /// given as fname.0 fname.1 etc.
    ///
    /// \return simulation time of snapshot
    /// \param[in]  fname basename of data file(s).
    /// \param[in]  read  data to read (maximum is mxvkURH)
    /// \param[out] got   data which actually have been read
    /// \param[in]  rec   size of FORTRAN record header; must be 4 or 8
    double read_gadget(const char*fname, fieldset read, fieldset&got,
		       unsigned rec=4) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// writes snapshot in gadget format to a single data file
    ///
    /// gadget data files MUST contain at least m,x,v,k,U. In addition, they
    /// may contain R,H or R,H,p,a, partly dependent on whether the file is a
    /// full dump or a initial-conditions file only. Here we allow for all of
    /// these three options.\n
    ///
    /// Data not existent in the snapshot but obliged to be written out will be
    /// reset to zero and a warning will issued, unless explicitly suppressed.
    ///
    /// \param[in] out   falcON::output to write to
    /// \param[in] time  time to write to snapshot
    /// \param[in] write data to write (minimum: mxvkU, maximum mxvkURHpa
    /// \param[in] warn  warn about missing data?
    /// \param[in] rec   size of FORTRAN record header; must be 4 or 8
    void write_gadget(output&out, double time, fieldset write, bool warn=0,
		      unsigned rec=4)
      const falcON_THROWING;
#endif
    //==========================================================================
    // \name miscellaneous
    /// set body keys to bodyindeces
    /// \note the bodyindex is defined via that of the first body in each
    /// block. In a MPI parallel situation, this is not identical to the LOCAL
    /// running number (available as localindex()), but reflects the global
    /// running number!
    void reset_keys() {
      if(BITS.contain(fieldbit::k))
	for(iterator b = begin_all_bodies(); b; ++b) 
	  b.key() = bodyindex(b);
    }
    //--------------------------------------------------------------------------
    /// reset body flags to hold only body type information
    void reset_flags() const {
      if(BITS.contain(fieldbit::f))
	for(const block* p=FIRST; p; p=p->next())
	  p->reset_flags();
    }
    //--------------------------------------------------------------------------
    void mark_sph_data_read    () const { SPHC = 0; }
    void mark_sph_data_changed () const { SPHC = 1; }
    //--------------------------------------------------------------------------
    void mark_srce_data_read   () const { SRCC = 0; }
    void mark_srce_data_changed() const { SRCC = 1; }
    //--------------------------------------------------------------------------
    /// return the total mass in bodies of given type
    real TotalMass(bodytype) const;
    /// return to total mass in bodies of all types
    real TotalMass() const {
      real M(zero);
      for(bodytype t; t; ++t)
	M += TotalMass(t);
      return M;
    }
    //==========================================================================
    // \name creating sorted index tables
    //
    /// \brief Create an index table sorted in \a func(body) for all bodies 
    /// flagged not to be ignored (in_subset()).
    ///
    /// \param[out] T    table of bodies::index, sorted
    /// \param[in]  func function for property to be sorted
    void sorted(Array<index>&T, real(*func)(iterator const&))
      const falcON_THROWING;
    /// \brief Create an index table sorted in \a func(body) for all bodies
    /// flagged not to be ignored (in_subset()) and also generate a sorted table
    /// of quantities
    ///
    /// \param[out] T    table of bodies::index, sorted
    /// \param[out] Q    table of quantity, sorted
    /// \param[in]  func function for property to be sorted
    void sorted(Array<index>&T, Array<real>&Q,
		real(*func)(iterator const&)) const falcON_THROWING;
    /// \brief Create an index table of the K bodies in_subset() which are
    /// nearest to a given body, 
    /// \param[in]  B  body to find neighbours for
    /// \param[in]  K  number of neighbours (including body B itself) to find
    /// \param[out] I  indices of neighbours, sorted in ascending distance to B
    /// \return  actual number of neighbours found in subset (could be < K).
    /// \note Involves a (single) loop over all bodies, hence not very
    ///       efficient. To be used occasionally only.
    /// \note Used by manipulator bound_centre
    unsigned findNeighbours(const iterator&B, unsigned K, Array<index>&I) const
      falcON_THROWING;
    //--------------------------------------------------------------------------
    bool CheckData(fieldset s, const char*f, int l) const
    {
      if(debug(6) && !have_all(s)) {
	DebugInfoN(" [%s:%d]: bodies data required but not present: \"%s\"\n",
		  f,l,word(all_data().missing(s)));
	return false;
      }
      return true;
    }
  protected:
    //==========================================================================
    /// method used by ebodies
    bodies(char, const unsigned[bodytype::NUM]) falcON_THROWING;
    /// method used by ebodies
    void reset(char, fieldbit, void*) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// set block::FIRST
    /// The block's FIRST data are set such that the bodyindices of bodies of
    /// type t on start at F[t].
    /// \param F array with first bodyindex per body type.
    void reset_firsts(unsigned F[bodytype::NUM]);
  private:
    //==========================================================================
    //                                                                          
    // private member methods                                                   
    //                                                                          
    //==========================================================================
    // set up blocks to hold N[t] bodies of type t                              
    void set_data(const unsigned[bodytype::NUM]) falcON_THROWING;
    //--------------------------------------------------------------------------
    // set block::FIRST and NALL, NBOD & NTOT
    void set_firsts();
    //--------------------------------------------------------------------------
    // delete all blocks and reset related data                                 
    void del_data() falcON_THROWING;
    //--------------------------------------------------------------------------
    // add a block into our linking & reset block::FIRST for all blocks
    /// \note we do NOT check for block number overflow
    void add_block(block*B);
    //--------------------------------------------------------------------------
    /// \brief find or create a block whose first free body is also the first
    /// of a contiguous set of N free bodies.
    ///
    /// If no chunk of \a N contiguous free bodies of type \a t can be found,
    /// we allocate new bodies such that a total of at least \a Na
    /// bodies of type \a t are free.\n
    ///
    /// \param[in] N   # bodies to guarantee avtivatablity of
    /// \param[in] t   type of bodies
    /// \param[in] Na  minimum # bodies to allocate should we need to allocate
    block*ensure_contiguous(unsigned N, bodytype t, unsigned Na=0);
    //--------------------------------------------------------------------------
    // erase a block from our linking & reset block::FIRST for all blocks
    void erase_block(block*B);
    //--------------------------------------------------------------------------
    /// create new block and link it in
    /// \return       new block
    /// \param[in] t  bodytype for new block
    /// \param[in] Na # bodies to allocate
    /// \param[in] Nb # bodies to activate (Nb <= Na)
    /// \param[in] f  data fields to support (data will be allocated)
    block*new_block(bodytype t, unsigned Na=0, unsigned Nb=0,
		    fieldset f=fieldset::empty) falcON_THROWING;
    //--------------------------------------------------------------------------
    void remove_empty_blocks(bool=1) falcON_THROWING;
    //==========================================================================
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // inline definitions of friends of class bodies and its sub-types and      //
  // related functions                                                        //
  // also serve to inject these functions into namespace falcON               //
  //                                                                          //
  // ///////////////////////////////////////////////////////////////////////////
  /// \relates falcON::bodies::iterator
  /// return const datum
  template<int BIT>
  typename field_traits<BIT>::type const&const_datum(bodies::iterator const&i) {
    return i. template const_dat<BIT>();
  }
  // ///////////////////////////////////////////////////////////////////////////
  inline std::ostream& operator<<(std::ostream&o, const bodies::index&i) {
    return o << i.no() << ':' << i.in();
  }
  // ///////////////////////////////////////////////////////////////////////////
  inline const bodies*const&bodies_of(bodies::iterator const&i) {
    return i.my_bodies();
  }
  inline const bodies::block*const&block_of(bodies::iterator const&i) {
    return i.my_block();
  }
  inline unsigned const&block_No(bodies::iterator const&i) {
    return i.my_block_No();
  }
  inline unsigned const&subindex(bodies::iterator const&i) {
    return i.K;
  }
  inline unsigned bodyindex(bodies::iterator const&i) {
    return i.my_index();
  }
  inline unsigned localindex(bodies::iterator const&i) {
    return i.my_localindex();
  }
  inline bool has_field(bodies::iterator const&i, fieldbit f) {
    return i.B && i.B->has_field(f);
  }
#define HasDatum(Bit,Name)				\
  inline bool has_##Name(bodies::iterator const&i) {	\
    return falcON::has_field(i,Bit);			\
  }
  DEF_NAMED(HasDatum)
#undef HasDatum
#define ConstAccess(Bit,Name)						\
  inline field_traits<Bit>::type const&Name(bodies::iterator const&i) {	\
    return i.B-> const_datum<Bit>(i.K);					\
  }
  DEF_NAMED(ConstAccess)
#undef ConstAccess
  inline bodies::iterator& bodies::iterator::next_in_subset() {
    do ++(*this);
    while(is_valid() && falcON::has_flag(*this) && !in_subset());
    return*this;
  }
  inline bodies::iterator& bodies::iterator::next_active() {
    do ++(*this);
    while(is_valid() && !is_active());
    return*this;
  }
  inline vect angmom(bodies::iterator const&i) {
    return falcON::pos(i) ^ falcON::vel(i);
  }
  inline bodytype const&type(bodies::iterator const&i) {
    return i.B->type();
  }
  inline double const&tau(bodies::iterator const&i) {
    return i.my_bodies()->have_level()? 
      i.my_bodies()->tau (level(i)) : i.my_bodies()->tau(0);
  }
  inline double const&tauq(bodies::iterator const&i) {
    return i.my_bodies()->have_level()? 
      i.my_bodies()->tauq(level(i)) : i.my_bodies()->tauq(0);
  }
  inline double const&tauh(bodies::iterator const&i) {
    return i.my_bodies()->have_level()? 
      i.my_bodies()->tauh(level(i)) : i.my_bodies()->tauh(0);
  }
  inline bool bodies::iterator::is_source() const {
    return falcON::mass(*this) != zero;
  }
  inline bool is_active     (const bodies::iterator&i) { return i.is_active(); }
  inline bool to_remove     (const bodies::iterator&i) { return i.to_remove(); }
  inline bool is_sph        (const bodies::iterator&i) { return i.is_sph(); }
  inline bool is_std        (const bodies::iterator&i) { return i.is_std(); }
  inline bool is_sink       (const bodies::iterator&i) { return i.is_sink(); }
  inline bool is_sticky     (const bodies::iterator&i) { return i.is_sticky(); }
  inline bool is_new        (const bodies::iterator&i) { return i.is_new(); }
  inline bool is_marked     (const bodies::iterator&i) { return i.is_marked(); }
  inline bool in_subset     (const bodies::iterator&i) { return i.in_subset(); }
  inline bool is_source     (const bodies::iterator&i) { return i.is_source(); }
  inline bool may_go_longer (const bodies::iterator&i) {
    return i.may_go_longer();
  }
  inline bool may_go_shorter(const bodies::iterator&i) {
    return i.may_go_shorter();
  }
  inline std::ostream& operator<<(std::ostream&o, const bodies::iterator&i) {
    return i? o << i.my_index() : o << "nil";
  }
  // ///////////////////////////////////////////////////////////////////////////
  inline double const&bodies::TimeSteps::tau(bodies::iterator const&i) const {
    return tau (falcON::level(i));
  }
  inline double const&bodies::TimeSteps::tauq(bodies::iterator const&i) const {
    return tauq(falcON::level(i));
  }
  inline double const&bodies::TimeSteps::tauh(bodies::iterator const&i) const {
    return tauh(falcON::level(i));
  }
#define CheckMissingBodyData(B,F) (B)->CheckData((F),__FILE__,__LINE__)
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // falcON::body                                                               
  //                                                                            
  /// a body is simply synonymous to bodies::iterator                           
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  typedef bodies::iterator body;
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  //  class falcON::snapshot                                                    
  //                                                                            
  /// Holds and manages snapshot data (= body data + time)                      
  ///                                                                           
  /// Derived from class falcON::bodies, inheriting all public members.         
  /// In addition, methods for snapshot I/O to file are provided.               
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class snapshot : public bodies
  {
    bool             INIT;
    double           TINI;
    mutable double   TIME;
    void            *PBNK;
    ParallelSnapshot*PARA;                         // parent if MPI parallel    
    //--------------------------------------------------------------------------
  protected:
    void set_parallel(ParallelSnapshot*P) { PARA = P; }
  public:
    ParallelSnapshot      *parallel()       { return PARA; }
    const ParallelSnapshot*parallel() const { return PARA; }
    //==========================================================================
    // \name initial-time information and manipulation                         
    bool const&has_initial_time() const { return INIT; }  ///< has initial time?
    double const&initial_time() const { return TINI; }    ///< get initial time
    void init_time(double t) {                            ///< set initial time
      INIT = true;
      TINI = t;
    }
    //==========================================================================
    // \name time information and manipulation                                 
    double const&time() const { return TIME; }            ///< get time
    void advance_time_by(double d) const { TIME += d; }   ///< advance time
    void set_time(double t) const { TIME  = t; }          ///< set time
    void reset_time() const { TIME  = 0.; }               ///< set time to zero
  private:
    void __add_pointer(const void*, const char*, size_t, const char*)
      const falcON_THROWING;
    void __set_pointer(const void*, const char*, size_t, const char*)
      const falcON_THROWING;
    void __rep_pointer(const void*, const char*, size_t, const char*)
      const falcON_THROWING;
    const void*__get_pointer(const char*, size_t, const char*) const
      falcON_THROWING;
  public:
    //==========================================================================
    // \name pointer bank: support for communication between manipulators      
    //
    /// add a pointer to pointer bank (see also set_pointer()).  The idea is
    /// that two or more manipulators may communicate data, for instance one
    /// manipulator determines the centre with respect to which another
    /// manipulator analyses the snapshot. For this purpose, the first
    /// manipulator registers a pointer to the centre position with the pointer
    /// bank under a handle, say \c "xcen", using members add_pointer() or
    /// set_pointer(). Any other manipulator may access that pointer with the
    /// same handle using the member pointer() or get_pointer().\n
    ///
    /// The sizeof(T) and nameof(T) are also remembered in the pointer bank.
    /// If the handle is already known to the bank, an error is thrown.
    /// A NULL pointer will not be added to the bank.
    /// \param T        (template parameter) type of object pointed to
    /// \param[in] pter pointer to be stored in bank
    /// \param[in] hand handle for the pointer
    template<typename T>
    void add_pointer(const T*pter, const char*hand) const falcON_THROWING {
      __add_pointer(pter, hand, sizeof(T), nameof(T));
    }
    /// set a pointer to pointer bank (see also add_pointer()).
    ///
    /// If the handle is already known to the bank, the old pointer is
    /// overridden. nameof(T) and sizeof(T) must match with old entry (if
    /// present). In case of a NULL pointer, any existing pointer will be
    /// deleted.
    /// \param T        (template parameter) type of object pointed to
    /// \param[in] pter pointer to be stored in bank
    /// \param[in] hand handle for the pointer
    template<typename T>
    void set_pointer(const T*pter, const char*hand) const falcON_THROWING {
      __set_pointer(pter, hand, sizeof(T), nameof(T));
    }
    /// return pointer to a given handle (see also add_pointer()).
    ///
    /// If the handle is unknown in the pointer bank, a zero pointer is
    /// returned. If the handle is known, but either the sizeof(T) or nameof(T)
    /// don't match, an error is thrown.
    /// \return         pointer referred to by handle (or null)
    /// \param[in] hand handle for pointer wanted
    ///
    /// Example code: \code
    ///   const snapshot*shot;
    ///   const vect    *x0 = shot->pointer<vect>("xcen");
    /// \endcode
    /// \note A template specification "<vect>" is required, because the
    /// compiler cannot determine the return type from the arguments.
    template<typename T>
    const T*pointer(const char*hand) const falcON_THROWING { 
      return static_cast<const T*>(__get_pointer(hand, sizeof(T), nameof(T)));
    }
    /// get a pointer to a given handle, routine version of pointer()
    ///
    /// If the handle is unknown in the pointer bank, a zero pointer is
    /// returned. If the handle is known, but either the sizeof(T) or nameof(T)
    /// don't match, an error is thrown.
    /// \param[out] pter pointer referred to by handle (or null)
    /// \param[in]  hand handle for pointer wanted
    ///
    /// Example code: \code
    ///   const snapshot*shot;
    ///   const vect    *x0;
    ///   shot->get_pointer(x0,"xcen");
    /// \endcode
    /// \note In contrast to pointer(), no template specification "<vect>" is
    ///       required.
    template<typename T>
    void get_pointer(const T* &pter, const char*hand) const falcON_THROWING {
      pter = static_cast<const T*>(__get_pointer(hand, sizeof(T), nameof(T)));
    }
    /// delete an entry from the pointer bank
    ///
    /// if key not found in pointer bank, no action is taken
    /// \param[in] hand handle for pointer to be deleted from bank
    void del_pointer(const char*hand) const;
    //@}
    //==========================================================================
    // \name construction & related
    //--------------------------------------------------------------------------
    /// Constructor 0: just give the fields to be supported,
    /// used in NBodyCode::NBodyCode() of file nbody.h
    explicit
    snapshot(fieldset Bd= fieldset(DefaultBits)) falcON_THROWING :
    bodies(Bd), INIT(false), TINI(0.), TIME(0.), PBNK(0), PARA(0) {}
    //--------------------------------------------------------------------------
    /// Constructor 1 (new version)
    ///
    /// \param[in] t   time of snapshot
    /// \param[in] N   number of bodies per bodytype
    /// \param[in] Bd  body data fields to be allocated
    explicit 
    snapshot(double         t,
	     const unsigned N[bodytype::NUM],
	     fieldset       Bd= fieldset(DefaultBits)) falcON_THROWING :
    bodies(N,Bd),
    INIT(true), TINI(t), TIME(t), PBNK(0), PARA(0) {}
    //--------------------------------------------------------------------------
    /// copy constructor from bodies
    ///
    /// make a partial copy, only copying data specified by 3rd arg
    ///
    /// \param[in] t   time of snapshot
    /// \param[in] B   set of bodies
    /// \param[in] Bd  only copy these body data fields
    /// \param[in] F   only copy bodies matching this flag
    /// \param[in] T   only copy bodies of type contained in T
    snapshot(double       t,
	     bodies const&B,
	     fieldset     Bd=fieldset::all,
	     flags        F =flags::empty,
	     bodytypes    T =bodytypes::all) falcON_THROWING :
    bodies(B,Bd,F,T), INIT(true), TINI(t), TIME(t), PBNK(0), PARA(0) {}
    //--------------------------------------------------------------------------
    /// copy constructor from snapshot
    ///
    /// make a partial copy, only copying data specified by 3rd arg
    ///
    /// \param[in] S   set of bodies
    /// \param[in] Bd  only copy these body data fields
    /// \param[in] F   only copy bodies matching this flag
    /// \param[in] T   only copy bodies of type contained in T
    explicit
    snapshot(snapshot const&S,
	     fieldset       Bd=fieldset::all,
	     flags          F =flags::empty,
	     bodytypes      T =bodytypes::all) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// Copy: replace our data (if any) with those of other bodies (\b not 
    /// \b yet \b implemented)
    ///
    /// We only copy data specified by 2nd argument and copy only those bodies
    /// whose flag matches 3rd argument. Equivalent to (but more efficient than)
    /// destruction followed by constructor 2.
    ///
    /// \note \b NOT \b YET \b IMPLEMENTED
    ///
    /// \param[in] S   snapshot to be copied
    /// \param[in] Bd  body data to be copied
    /// \param[in] F   flag for bodies to be copied
    void copy(snapshot const&S,
	      fieldset       Bd=fieldset::all,
	      flags          F =flags::empty) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// destruction
    ~snapshot();
    //==========================================================================
#ifdef falcON_NEMO
    //==========================================================================
    // \name NEMO data I/O
    //==========================================================================
    /// Read a single NEMO snapshot if time in input is in desired time range;
    /// on return this mirrors the snapshot on disk (if it has been read)
    ///
    /// If time in input is in time range (4th arg), we read all data indicated
    /// by 3rd argument for all bodies in input (re-allocating if necessary).
    /// Previous data are lost.
    ///
    /// If the initial time was previously unset, we set it now.
    ///
    /// \return       true if time in input was in desired time range
    /// \param[in]  Is     NEMO input stream containing snapshot
    /// \param[out] Read   body data fields which have been read
    /// \param[in]  Get    body data fields to be read
    /// \param[in]  range  time range in NEMO range format
    ///                    (default: read next snapshot in stream)
    /// \param[in]  warn   issue falcON::warning()s about missing data
    ///                    (default: issue warning)
    bool read_nemo(nemo_in const&Is,
		   fieldset     &Read,
		   fieldset      Get,
		   const char   *range=0,
		   bool          warn=1) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// Read a part of a NEMO snapshot into a given body range
    ///
    /// Here, we read data for \e Nread (4th argument) bodies from the snapshot
    /// input (1st arg) into the bodies starting at \e start (3rd argument).
    /// Enough bodies must be have been allocated previously.
    ///
    /// If the initial time was previously unset, we set it now.
    ///
    /// \return      body data fields that have been read
    /// \param[in] Si     snapshot input
    /// \param[in] Bd     body data fields to be read
    /// \param[in] start  first body to be read into
    /// \param[in] warn   issue falcON::warning()s about missing data
    ///                   (default: issue warning)
    /// \param[in] Nread  read this many (default: all in input)
    fieldset read_part(snap_in const &Si,
		       fieldset       Bd,
		       iterator const&start,
		       bool           warn=1,
		       unsigned       Nread=0) falcON_THROWING;
    //--------------------------------------------------------------------------
    /// Generate a NEMO snapshot on disk from all bodies
    ///
    /// \param[in] Os  nemo output stream
    /// \param[in] Bd  body data fields to be written
    void write_nemo(nemo_out const&Os,
		    fieldset       Bd = fieldset::nemo) const falcON_THROWING;
    //--------------------------------------------------------------------------
    /// Generate a NEMO snapshot on disk from a subset of all bodies
    ///
    /// We write body data fields specified by \e Bd (2nd argument) for \e 
    /// Nwrite (4th argument) bodies starting at \e start (3rd argument) in
    /// NEMO snapshot format.
    ///
    /// \param[in] Os     nemo output stream
    /// \param[in] Bd     body data fields to be written
    /// \param[in] start  first body to write
    /// \param[in] Nwrite write this many (default: all)
    void write_nemo(nemo_out const&Os,
		    fieldset       Bd,
		    iterator const&start,
		    unsigned       Nwrite=0) const falcON_THROWING;
#endif // falcON_NEMO
#ifdef falcON_REAL_IS_FLOAT
    //--------------------------------------------------------------------------
    /// Reads a single snapshot from file(s) written in gadget2 data format 1
    ///
    /// Gadget allows a single snapshot to be distributed over more than one
    /// file. In this case the file names are assumed to be derived from the one
    /// given as fname.0 fname.1 etc.
    ///
    /// \return     data which actually have been read
    /// \param[in]  fname basename of data file(s).
    /// \param[in]  read  data to read (maximum is mxvkURH)
    /// \param[in]  rec   size of FORTRAN record header; must be 4 or 8
    fieldset read_gadget(const char*fname, fieldset read, unsigned rec=4)
      falcON_THROWING
    {
      fieldset got;
      double t = bodies::read_gadget(fname, read, got, rec);
      set_time(t);
      if(!has_initial_time()) init_time(t);
      return got;
    }
    //--------------------------------------------------------------------------
    /// writes snapshot in gadget format to a single data file
    ///
    /// gadget data files MUST contain at least m,x,v,k,U. In addition, they
    /// may contain R,H or R,H,p,a, partly dependent on whether the file is a
    /// full dump or only a initial-conditions file. Here we allow for three
    /// options.\n
    /// Data not existent in the snapshot but obliged to be written out will be
    /// reset to zero and, if \a warn=true, a warning will issued.
    ///
    /// \param[in]  out   falcON::output to write to
    /// \param[in]  write data to write (minimum: mxvkU, maximum mxvkURHpa
    /// \param[in]  warn  warn  about missing data?
    /// \param[in]  rec   size of FORTRAN record header; must be 4 or 8
    void write_gadget(output &out, fieldset write, bool warn=0, unsigned rec=4)
      const falcON_THROWING
    {
      bodies::write_gadget(out, time(), write, warn, rec);
    }
#endif // falcON_REAL_IS_FLOAT
  };// class snapshot
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
falcON_TRAITS(falcON::bodies::index,"bodies::index");
falcON_TRAITS(falcON::bodies,"bodies");
falcON_TRAITS(falcON::snapshot,"snapshot");
// /////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \name macros for looping over bodies                                        
//@{                                                                            
////////////////////////////////////////////////////////////////////////////////
#ifndef LoopAllBodies           /* loop all bodies                           */
/// This macro provides an easy way to loop over all bodies.
///
/// A typical usage would look like this \code
///   snapshot *S;
///   LoopAllBodies(S,b) {
///     b.pos() += dt*vel(b);
///     b.vel() += dt*acc(b);
///   } \endcode
///
/// \param PTER  valid pointer to falcON::bodies (or falcON::snapshot)
/// \param NAME  name given to loop variable (of type falcON::body)
#define LoopAllBodies(PTER,NAME)			\
  for(falcON::body NAME=(PTER)->begin_all_bodies();	\
                   NAME; ++NAME)
#endif
//------------------------------------------------------------------------------
#ifndef LoopTypedBodies         /* loop all bodies of a given type           */
/// This macro provides an easy way to loop over all bodies of a given 
/// falcON::bodytype
///
/// A typical usage would look like this \code
///   snapshot *S;
///   bodytype  t;
///   LoopTypedBodies(S,b,t) {
///     b.pos() += dt*vel(b);
///     b.vel() += dt*acc(b);
///   } \endcode
///
/// \param PTER  valid pointer to falcON::bodies (or falcON::snapshot)
/// \param NAME  name given to loop variable (of type falcON::body)
/// \param TYPE  the bodytype of the bodies to be looped over
#define LoopTypedBodies(PTER,NAME,TYPE)					\
  for(falcON::body NAME= (PTER)->begin_typed_bodies(TYPE);		\
                   NAME  != (PTER)->end_typed_bodies(TYPE); ++NAME)
#endif
//------------------------------------------------------------------------------
#ifndef LoopSPHBodies           /* loop all SPH bodies                       */
/// This macro provides an easy way to loop over all SPH bodies
///
/// A typical usage would look like this \code
///   snapshot *S;
///   LoopSPHBodies(S,b,t) {
///     b.pos() += dt*vel(b);
///     b.vel() += dt*acc(b);
///     b.uin() += dt*udin(b);
///   } \endcode
///
/// \param PTER  valid pointer to falcON::bodies (or falcON::snapshot)
/// \param NAME  name given to loop variable (of type falcON::body)
#define LoopSPHBodies(PTER,NAME)		\
  LoopTypedBodies(PTER,NAME,bodytype::gas)
#endif
//------------------------------------------------------------------------------
#ifndef LoopSTDBodies           /* loop all STD bodies                   */
/// This macro provides an easy way to loop over all STD bodies
///
/// A typical usage would look like this \code
///   snapshot *S;
///   LoopSTDBodies(S,b,t) {
///     b.pos() += dt*vel(b);
///     b.vel() += dt*acc(b);
///   } \endcode
///
/// \param PTER  valid pointer to falcON::bodies (or falcON::snapshot)
/// \param NAME  name given to loop variable (of type falcON::body)
#define LoopSTDBodies(PTER,NAME)		\
  LoopTypedBodies(PTER,NAME,bodytype::std)
#endif
//------------------------------------------------------------------------------
#ifndef LoopSinkBodies           /* loop all STD bodies                   */
/// This macro provides an easy way to loop over all SINK bodies
///
/// A typical usage would look like this \code
///   snapshot *S;
///   LoopSinkBodies(S,b,t) {
///     b.pos() += dt*vel(b);
///     b.vel() += dt*acc(b);
///   } \endcode
///
/// \param PTER  valid pointer to falcON::bodies (or falcON::snapshot)
/// \param NAME  name given to loop variable (of type falcON::body)
#define LoopSinkBodies(PTER,NAME)		\
  LoopTypedBodies(PTER,NAME,bodytype::sink)
#endif
//------------------------------------------------------------------------------
#ifndef LoopBodyPairs           /* loop all bodies after a given body        */
/// This macro provides an easy way to loop over all body pairs
///
/// A typical usage would look like this \code
///   snapshot *S;
///   LoopAllBodies(S,bi)
///     LoopBodyPairs(bi,bj) {
///       real Rij_squared = dist_sq(pos(bi),pos(bj));
///     } \endcode
/// \param NAME1  name of the body in an outer loop
/// \param NAME2  name given to inner-loop variable (of type falcON::body)
#define LoopBodyPairs(NAME1,NAME2)		\
  for(falcON::body NAME2(NAME1,1); NAME2; ++NAME2)
#endif
//------------------------------------------------------------------------------
#ifndef LoopActiveBodies        /* loop all active bodies */
/// This macro provides an easy way to loop over all active bodies
///
/// A typical usage would look like this \code
///   snapshot *S;
///   LoopActiveBodies(S,b) {
///     b.pos() += dt*vel(b);
///     b.vel() += dt*acc(b);
///   } \endcode
///
/// \param PTER  valid pointer to falcON::bodies (or falcON::snapshot)
/// \param NAME  name given to loop variable (of type falcON::body)
#define LoopActiveBodies(PTER,NAME)					\
  for(falcON::body NAME=(PTER)->begin_active_bodies();			\
      NAME; NAME.next_active())
#endif
//------------------------------------------------------------------------------
#ifndef LoopSubsetBodies        /* loop all bodies which are in_subset() */
/// This macro provides an easy way to loop over all bodies in a subset
///
/// A typical usage would look like this \code
///   snapshot *S;
///   LoopSubsetBodies(S,b) {
///     b.pos() += dt*vel(b);
///     b.vel() += dt*acc(b);
///   } \endcode
///
/// \param PTER  valid pointer to falcON::bodies (or falcON::snapshot)
/// \param NAME  name given to loop variable (of type falcON::body)
/// \version 04-nov-08 WD debugged
#define LoopSubsetBodies(PTER,NAME)					\
  for(falcON::body NAME=(PTER)->begin_subset_bodies();			\
      NAME; NAME.next_in_subset())
#endif
//------------------------------------------------------------------------------
//@}
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_body_h
