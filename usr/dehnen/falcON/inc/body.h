// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// body.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_nbod_h
#define falcON_included_nbod_h 1

#ifndef falcON_included_auxx_h
#  include <public/auxx.h>
#endif
#ifndef falcON_included_flag_h
#  include <public/flag.h>
#endif
#ifndef falcON_included_nbio_h
#  include <public/nbio.h>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::bodies_data<>                                                //
  //                                                                          //
  // holds the body data as a structure of arrays (SoA) but provides an       //
  // iterator that mimicks an array of structures (AoS) implementation.       //
  //                                                                          //
  // design april 2004, using bits & pieces (iterator) from older designs     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<int BTYPE> class bodies_data {
    bodies_data& operator=(bodies_data const&);    // not provided              
    //--------------------------------------------------------------------------
    // static data                                                              
    //--------------------------------------------------------------------------
  public:
    static const int bodies_tag = BTYPE;           // used in templates         
    //--------------------------------------------------------------------------
    // data                                                                     
    //                                                                          
    // the interpretation of DATA[] depends on the template argument            
    // If N=0, the BITS may indicate DATA[] support, even if DATA[b]=0.         
    //--------------------------------------------------------------------------
  private:
    uint         NBOD;                             // # bodies                  
    uint         NSPH;                             // # SPH particles           
    void        *DATA[IO_MAX];                     // pointers to body data     
    io           BITS;                             // which data do we support? 
    mutable bool SRCC;                             // source data changed?      
    mutable bool SPHC;                             // SPH data changed?         
    mutable bool CUSE;                             // flag for tree usage change
    //--------------------------------------------------------------------------
    // construction: only set N, no data yet                                    
  protected:
    bodies_data(uint nb = 0u,                      //[I: # bodies]              
		uint ns = 0u,                      //[I: # SPH particles]       
		bool cu = 0) :                     //[I: CUSE flag]             
      NBOD(nb), NSPH(ns), BITS(io::o), SRCC(1), SPHC(1), CUSE(cu) {
      for(int i=0; i!=IO_MAX; ++i) DATA[i] = 0;
    }
    //--------------------------------------------------------------------------
    // reset N (the user has to ensure that the data arrays reflect this)       
    void setN(uint nb, uint ns = 0u) {
      NBOD = nb;
      NSPH = ns;
    }
    //--------------------------------------------------------------------------
    void add_data_void(int bit, void*D) {
      DATA[bit] = D;
      BITS |= 1<<bit;
    }
    //--------------------------------------------------------------------------
    void del_data_void(int bit) {
      DATA[bit] = 0;
      BITS &=~(1<<bit);
    }
    //--------------------------------------------------------------------------
    // set data & bits                                                          
    void set_data_void(int bit, void*D) {
      DATA[bit] = D;
      if(D) BITS |=  1<<bit;
      else  BITS &=~(1<<bit);
    }
    //--------------------------------------------------------------------------
  public:
    // are all the fields in i supported?                                       
    bool has(                                      // R: fields suppored?       
	     io i) const {                         // I: bits of fields         
      return BITS.contains(i);                     // fields contained in BITS  
    }
    //--------------------------------------------------------------------------
    // give the bits of the data currently supported                            
    io  const &all_bits() const {
      return BITS;
    }
    //--------------------------------------------------------------------------
    // give the bits of the data currently NOT supported                        
    io  has_not(io const&from) const {
      return io((~BITS) & from);
    }
    //--------------------------------------------------------------------------
    void after_tree_growth     () const { CUSE = 0; }
    void mark_tree_usage_change() const { CUSE = 1; }
    //--------------------------------------------------------------------------
    void mark_sph_data_read    () const { SPHC = 0; }
    void mark_sph_data_changed () const { SPHC = 1; }
    //--------------------------------------------------------------------------
    void mark_srce_data_read   () const { SRCC = 0; }
    void mark_srce_data_changed() const { SRCC = 1; }
    //==========================================================================
    //                                                                          
    // 2. Data Access via Members Methods                                       
    //                                                                          
    // Note, we do not sanity or bound checks whatoever. The user has to ensure 
    // (using has()) that a given data field is non null.                       
    //                                                                          
    //--------------------------------------------------------------------------
    bool const&changes_in_tree_usage_flags () const { return CUSE; }
    uint const&N_bodies                    () const { return NBOD; }
    uint const&N_sph                       () const { return NSPH; }
    bool const&srce_data_changed           () const { return SRCC; }
    bool const&sph_data_changed            () const { return SPHC; }
    //--------------------------------------------------------------------------
  protected:
    void*const&data_void(int bit) const {
      return DATA[bit];
    }
  public:
    //--------------------------------------------------------------------------
    // 2.1 access to data arrays via templates                                  
    // 2.1.1 via templates taking the bit# (0 <= ... < IO_NQUANT)               
    template<int BIT>
    typename field_bit<BIT,BTYPE>::d_type* data_bit() const {
      return
	static_cast<typename field_bit<BIT,BTYPE>::d_type*>(DATA[BIT]);
   }
    //--------------------------------------------------------------------------
    // 2.1.2 via templates taking io = 1<<bit                                   
    template<int IO>
    typename field_io<IO,BTYPE>::d_type* data_io() const {
      return data_bit<field_io<IO,BTYPE>::bit>();
    }
    //--------------------------------------------------------------------------
    // 2.2 const & non-const access to data elements via templates              
    // 2.2.1 via templates taking the bit# (0 <= ... < IO_NQUANT)               
    template<int BIT>
    typename field_bit<BIT,BTYPE>::r_type datum_bit(int i) const {
      return field_bit<BIT,BTYPE>::element(data_void(BIT),i);
    }
    template<int BIT>
    typename field_bit<BIT,BTYPE>::cr_type const_datum_bit(int i) const {
      return field_bit<BIT,BTYPE>::c_element(data_void(BIT),i);
    }
    //--------------------------------------------------------------------------
    // 2.2.2 via templates taking io = 1<<bit                                   
    //--------------------------------------------------------------------------
    template<int IO>
    typename field_io<IO,BTYPE>::r_type datum_io(int i) const {
      return field_io<IO,BTYPE>::element(data_void(field_io<IO,BTYPE>::bit),i);
    }
    template<int IO>
    typename field_io<IO,BTYPE>::cr_type const_datum_io(int i) const {
      return field_io<IO,BTYPE>::c_element(data_void(field_io<IO,BTYPE>::bit),
					   i);
    }
    //--------------------------------------------------------------------------
#define ACCESS(NAME,BIT)						\
    typename field_bit<BIT,BTYPE>::d_type*     NAME##_s()      const {	\
      return data_bit<BIT>(); }						\
    typename field_bit<BIT,BTYPE>::r_type      NAME    (int i) const {	\
      return datum_bit<BIT>(i); }					\
    typename field_bit<BIT,BTYPE>::cr_type c_##NAME(int i) const {	\
      return const_datum_bit<BIT>(i); }
    //--------------------------------------------------------------------------
    // 2.3 access to data arrays and elments via individual named member methods
    ACCESS(mass,0); 
    ACCESS(pos,1);
    ACCESS(vel,2);
    ACCESS(eps,3);
    ACCESS(flg,4);
    ACCESS(key,5);
    ACCESS(pot,6);
    ACCESS(pex,7);
    ACCESS(acc,8);
    ACCESS(rho,9);
    ACCESS(aux,10);
    ACCESS(level,11);
    ACCESS(num,12);
    ACCESS(size,13);
#ifdef falcON_SPH
    ACCESS(snum,14);
    ACCESS(uin,15);
    ACCESS(uprd,16);
    ACCESS(udex,17);
    ACCESS(udin,18);
    ACCESS(ent,19);
    ACCESS(srho,20);
    ACCESS(drho,21);
    ACCESS(vprd,22);
    ACCESS(sigq,23);
    ACCESS(temp,24);
    ACCESS(hdot,25);
#endif
#undef ACCESS
    amom angmon (int i) const { return c_pos(i) ^ c_vel(i); }
    //==========================================================================
    //                                                                          
    // 4. Miscellaneous                                                         
    //                                                                          
    //--------------------------------------------------------------------------
    uint const&Number(int b) const {
      return b < IO_NOSPH? NBOD : NSPH;
    }
    //--------------------------------------------------------------------------
    void reset_flags() const {
      if(flg_s()) {
	for(register int i=0; i!=NBOD; ++i) flg(i).reset();
	after_tree_growth();
      } else
	falcON_WarningF("flags not supported","bodies::reset_flags()");
    }
    //--------------------------------------------------------------------------
    void flag_all_as_active() const {              // flag all bodies as active 
      if(flg_s())
	for(register int i=0; i!=NBOD; ++i) flg(i).add(flag::ACTIVE);
      else falcON_ErrorF("flags not supported","bodies::flag_all_as_active()");
    }
    //--------------------------------------------------------------------------
    void flag_as_sph() const {                     // flag the first NSPH as sph
      if(flg_s())
	for(register int i=0; i!=NSPH; ++i) flg(i).add(flag::SPH);
      else falcON_ErrorF("flags not supported","bodies::flag_as_sph()");
    }
    //--------------------------------------------------------------------------
    void allow_tree_usage(int i) const {
      if(!is_in_tree(flg(i))) {                    // IF body not in tree       
	flg(i).un_set(flag::NOT_IN_TREE);          //   unflag 'not in tree'    
	mark_tree_usage_change();                  //   set flag for change     
      }                                            // ENDIF                     
    }
    //--------------------------------------------------------------------------
    void forbid_tree_usage(int i) const {
      if(is_in_tree(flg(i))) {                     // IF body was in tree       
	flg(i).add(flag::NOT_IN_TREE);             //   flag 'not in tree'      
	mark_tree_usage_change();                  //   set flag for change     
      }                                            // ENDIF                     
    }
  };
  //////////////////////////////////////////////////////////////////////////////
#if defined(falcON_NEMO) && !defined(falcON_included_nmio_h)
  class nemo_in;                                   // forward declaration       
  class nemo_out;                                  // forward declaration       
#endif
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::bodies                                                       //
  //                                                                          //
  // holds the body data as a structure of arrays (SoA) but provides an       //
  // iterator that mimicks an array of structures (AoS) implementation.       //
  //                                                                          //
  // design april 2004, using bits & pieces (iterator) from older designs     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class bodies : public bodies_data<0> {
    bodies& operator=(bodies const&);              // not provided              
    //==========================================================================
    //                                                                          
    // 1. Data Management                                                       
    //                                                                          
    // On construction data arrays specified by the bits argument will be       
    // allocated. The member method add_fields() allows to allocate data arrays 
    // that are not yet supported, whereas del_fields() deallocates data arrays 
    // that were supported. The member method resetN(), if called with          
    // deviating N, deletes all old data arrays and then re-allocates them      
    // with the new sizes. changeN() is similar, but copies the old data into   
    // the new arrays as far as possible.                                       
    //                                                                          
    //--------------------------------------------------------------------------
  public:
    enum {
      NbyMaxBits = io::source | io::sink,          // gravity etc: max fields   
      NbyDefBits = io::gravity,                    // gravity etc: default      
      NbyNIOBits = NbyMaxBits & io::NEMO,          // gravity etc: nemo I/O     
      SPHMaxBits = io::sphmax,                     // SPH etc: max fields       
      SPHDefBits = io::sphmin,                     // SPH etc: default          
      SPHNIOBits = SPHMaxBits & io::NEMO,          // SPH etc: emo I/O          
      MaxBits    = NbyMaxBits | SPHMaxBits,        // maximum fields            
      DefBits    = NbyDefBits | SPHDefBits,        // default fields            
      NIOBits    = NbyNIOBits | SPHNIOBits         // nemo I/O fields           
    };
    //--------------------------------------------------------------------------
    // construction: set bits, N, allocate DATA[] if N!=0                       
    bodies(uint = 0,                               //[I: # bodies]              
	   io   = DefBits,                         //[I: data bits]             
	   uint = 0,                               //[I: # SPH particles]       
	   bool = false);                          //[I: CUSE flag]             
    //--------------------------------------------------------------------------
    // construction: make a partial copy, only copying data specified by 2nd arg
    bodies(bodies const&,                          // I: bodies                 
	   io   = io::all);                        //[I: which data to copy]    
    //--------------------------------------------------------------------------
    // copy body data indicated. If NBOD, NSPH deviating, adjust first          
    void copy(bodies const&,                       // I: bodies                 
	      io   = io::all);                     //[I: which data to copy]    
    //--------------------------------------------------------------------------
    // destruction: delete all supported data fields                            
    ~bodies();
    //--------------------------------------------------------------------------
    // adds unsupported DATA[] fields indicated by arg                          
    // supported DATA[] fields also contained in the arg are not affected       
    void add_fields(io);                           // I: which fields to add    
    //--------------------------------------------------------------------------
    // adds unsupported DATA[] field as indicated by arg                        
    // supported DATA[] fields are not affected                                 
    void add_field(int);                           // I: which field to add     
    //--------------------------------------------------------------------------
    // removes supported DATA[] fields inicated by arg                          
    void del_fields(io);                           // I: which fields to delete 
    //--------------------------------------------------------------------------
    // removes supported DATA[] field as indicated by arg                       
    void del_field(int);                           // I: which fields to delete 
    //--------------------------------------------------------------------------
    // resets DATA[] fields to arg                                              
    // DATA[] fields already supported AND contained in arg are not affected    
    void reset_fields(io);                         // I: which fields to support
    //--------------------------------------------------------------------------
    // reset NBODIES & NSPH. If different from old, previous data are deleted   
    void resetN(uint,                              // I: # bodies               
		uint = 0);                         //[I: # SPH particles]       
    //--------------------------------------------------------------------------
    // change NBODIES & NSPH. Old data are preserved for n < min(N_new,N_old)   
    void changeN(uint,                             // I: # bodies               
		 uint = 0);                        //[I: # SPH particles]       
    //==========================================================================
    //                                                                          
    // 2. Access to bodies DATA[] via sub-type iterator                         
    //                                                                          
    // Note that we cannot put this into the base class bodies_data, because    
    // this is a template which confuses the definition of friends of iterator. 
    //                                                                          
    //--------------------------------------------------------------------------
    class iterator {
      //........................................................................
      // basic stuff (construction, iteration, etc)                             
      //........................................................................
      friend class bodies;
    private:
      iterator();                                  // not implemented           
      //........................................................................
      // data members                                                           
      //........................................................................
      const bodies *B;                             // pointer to OBJECTS        
      unsigned      K;                             // index                     
      //........................................................................
      // private types                                                          
      //........................................................................
      enum {
	MARK_BODY   = 1<<9,                        // use  9th bit of flag      
	NOT_LONGER  = 1<<10,                       // use 10th bit of flag      
	NOT_SHORTER = 1<<11                        // use 11th bit of flag      
      };
    protected:
      //........................................................................
      // construction                                                           
      //........................................................................
      iterator(const bodies*b, const uint k) : B(b), K(k) {}
    public:
      iterator(const iterator&I)             : B(I.B), K(I.K) {}
      //........................................................................
      // forward iteration                                                      
      //........................................................................
      iterator       next ()            const { return iterator(B,K+1); }
      iterator& operator++()                  { ++K; return *this; }
      iterator  operator++(int)               { return iterator(B,K++); }
      iterator& operator+=(const int&k)       { K+=k; return *this; }
      iterator  operator+ (const int&k) const { return iterator(B,K+k); }
      //........................................................................
      // boolean methods                                                        
      //........................................................................
      bool  operator== (const iterator&I) const { return K == I.K; }
      bool  operator!= (const iterator&I) const { return K != I.K; }
      bool  operator<  (const iterator&I) const { return K <  I.K; }
      bool  operator<= (const iterator&I) const { return K <= I.K; }
      bool  operator>  (const iterator&I) const { return K >  I.K; }
      bool  operator>= (const iterator&I) const { return K >= I.K; }
      //........................................................................
      // access to bodies                                                       
      //........................................................................
      const bodies*const&mybodies() const { return B; }
      //........................................................................
      // data access:                                                           
      //   non-const:  via member methods                                       
      //   const    :  via friends                                              
      //........................................................................
      friend uint const&index(const iterator&I) { return I.K; }
      //........................................................................
      // non-const data access via template taking bit#                         
      template<int BIT>
      typename field_bit<BIT>::r_type data_bit() {
	return B-> template datum_bit<BIT>(K);
      }
      //........................................................................
      // non-const data access via template taking io                           
      template<int IO>
      typename field_io<IO>::r_type data_io() {
	return B-> template datum_io<IO>(K);
      }
      //........................................................................
      // const data access via template friend taking bit#                      
      template<int BIT> friend
      typename field_bit<BIT>::cr_type data_bit(iterator const&I) {
	return I.B-> template const_datum_bit<BIT>(I.K);
      }
      //........................................................................
      // const data access via template friend taking io                        
      template<int IO> friend
      typename field_io<IO>::cr_type data_io(iterator const&I) {
	return I.B-> template const_datum_io<IO>(I.K);
      }
      //........................................................................
#define ACCESS(NAME,BIT)					\
             field_bit<BIT>::r_type  NAME() {			\
	return iterator::data_bit<BIT>(); }			\
      friend field_bit<BIT>::cr_type NAME(iterator const&I) {	\
	return I.B-> const_datum_bit<BIT>(I.K); }
      //........................................................................
      // data access via named member methods (non-const) and friends (const)   
      ACCESS(mass,0);
      ACCESS(pos,1);
      ACCESS(vel,2);
      ACCESS(eps,3);
      ACCESS(flg,4);
      ACCESS(key,5);
      ACCESS(pot,6);
      ACCESS(pex,7);
      ACCESS(acc,8);
      ACCESS(rho,9);
      ACCESS(aux,10);
      ACCESS(level,11);
      ACCESS(num,12);
      ACCESS(size,13);
#ifdef falcON_SPH
      ACCESS(snum,14);
      ACCESS(uin,15);
      ACCESS(uprd,16);
      ACCESS(udex,17);
      ACCESS(udin,18);
      ACCESS(ent,19);
      ACCESS(srho,20);
      ACCESS(drho,21);
      ACCESS(vprd,22);
      ACCESS(sigq,23);
      ACCESS(temp,24);
      ACCESS(hdot,25);
#endif
#undef ACCESS
      friend amom angmon (iterator const&B) {
	return nbdy::pos(B) ^ nbdy::vel(B);
      }
      //........................................................................
      // other const methods                                                    
      //........................................................................
      friend amom angmom(iterator const&I) {
	return nbdy::pos(I) ^ nbdy::vel(I); }
      //........................................................................
      // flag manipulations                                                     
      //........................................................................
      bool flag_is_set       (int const&F) const {
	return nbdy::flg(*this).is_set(F); }
      operator const flag&   ()            const {
	return nbdy::flg(*this); }
      void flag_as_active    () { flg().add    (flag::ACTIVE); }
      void unflag_active     () { flg().un_set (flag::ACTIVE); }
      void flag_as_sticky    () { flg().add    (flag::STICKY); }
      void unflag_sticky     () { flg().un_set (flag::STICKY); }
      void flag_as_sph       () { flg().add    (flag::SPH); }
      void unflag_sph        () { flg().un_set (flag::SPH); }
      void allow_tree_usage  () { B->allow_tree_usage(K); }
      void forbid_tree_usage () { B->forbid_tree_usage(K); }
      void allow_longer      () { flg().un_set (NOT_LONGER); }
      void forbid_longer     () { flg().add    (NOT_LONGER); }
      void allow_shorter     () { flg().un_set (NOT_SHORTER); }
      void forbid_shorter    () { flg().add    (NOT_SHORTER); }
      void mark              () { flg().add    (MARK_BODY); }
      void unmark            () { flg().un_set (MARK_BODY); }
      //........................................................................
      // const boolean informations via mambers                                 
      //........................................................................
      bool is_source     () const { return nbdy::mass(*this) != zero; }
      bool may_go_longer () const { return !flag_is_set(NOT_LONGER);}
      bool may_go_shorter() const { return !flag_is_set(NOT_SHORTER); }
      bool is_marked     () const { return  flag_is_set(MARK_BODY); }
      //........................................................................
      // const boolean informations via friends                                 
      //........................................................................
      friend bool is_source     (const iterator&I) { return I.is_source(); }
      friend bool may_go_longer (const iterator&I) { return I.may_go_longer(); }
      friend bool may_go_shorter(const iterator&I) { return I.may_go_shorter();}
      friend bool is_marked     (const iterator&I) { return I.is_marked(); }
      //........................................................................
      // output: give just the index (an iterator acts like a pointer)          
      //........................................................................
      friend std::ostream& operator<<(std::ostream&o, const iterator&I) {
	return o<<nbdy::index(I);
      }
    };
    //--------------------------------------------------------------------------
    // access to iterators                                                      
    iterator body_no     (uint i) const { return iterator(this,i); }
    iterator bodyNo      (uint i) const { return iterator(this,i); }
    iterator begin_bodies()       const { return iterator(this,0); }
    iterator end_bodies  ()       const { return iterator(this,N_bodies()); }
    iterator last_bodies ()       const { return iterator(this,N_bodies()-1); }
    iterator back_bodies ()       const { return iterator(this,N_bodies()-1); }
    iterator begin_SPH   ()       const { return iterator(this,0); }
    iterator end_SPH     ()       const { return iterator(this,N_sph()); }
    iterator back_SPH    ()       const { return iterator(this,N_sph()-1); }
#ifdef falcON_NEMO
    //==========================================================================
    //                                                                          
    // 3. NEMO Data I/O                                                         
    //                                                                          
    // Input is done as follows:                                                
    //  - if the time is not in the range, we return but set time;              
    //  - we resetN(), i.e. old data are deleted if NBOD, NSPH differ;          
    //  - for fields with NEMO I/O which are both wanted AND available:         
    //    - add_fields()                                                        
    //    - read fields                                                         
    //                                                                          
    // Note that fields wanted but not available are NOT allocated.             
    // Differently from an earlier implementation, old data which are not       
    // superseeded will not be deleted, unless N differs.                       
    //                                                                          
    //--------------------------------------------------------------------------
    // read a single nemo particle set into bodies begin to begin+N             
  private:
    io read_nemo(                                  // R: data read              
		 nemo_in const&,                   // I: nemo input             
		 io            ,                   // I: what to read           
		 size_t       = 0);                //[I: first body to get]     
    //--------------------------------------------------------------------------
    // read a single nemo snapshot, supposed to be open                         
  public:
    bool read_nemo_particles(                      // R: was time in range?     
			     nemo_in const&,       // I: nemo input             
			     io           &,       // O: what has been read     
			     double       *,       // O: time                   
			     io            ,       // I: what to read           
			     char*        = 0,     //[I: time range]            
			     bool         = 1);    //[I: warn: missing data]    
    //--------------------------------------------------------------------------
    // read a single particle set combining several nemo input streams          
    bool read_nemo_particles(                      // R: was time in range?     
			     const nemo_in*,       // I: nemo inputs            
			     int   const  &,       // I: # nemo inputs          
			     io           *,       // O: what has been read     
			     uint         *,       // O: how many have ---      
			     double       * = 0,   //[O: time]                  
			     const io     * = 0,   //[I: what to read]          
			     char         * = 0,   //[I: time range]            
			     const bool     = 1);  //[I: warn: missing data]    
    //--------------------------------------------------------------------------
    // read a single nemo snapshop ignoring diagnose etc.                       
    bool read_nemo_snapshot (                      // R: was time in range?     
			     nemo_in const&,       // I: nemo input             
			     io           &,       // O: what has been read     
			     double       * = 0,   //[O: time]                  
			     io             = io::mxv, //[I: what to read]      
			     char         * = 0,   //[I: time range]            
			     bool           = 1);  //[I: warn: missing data]    
    //--------------------------------------------------------------------------
    // read a single snapshot combining several nemo input streams              
    bool read_nemo_snapshots(                      // R: was time in range?     
			     const nemo_in*,       // I: nemo inputs            
			     int           ,       // I: # nemo inputs          
			     io           *,       // O: what has been read     
			     uint         *,       // O: how many have ---      
			     double       * = 0,   //[O: time]                  
			     const io     * = 0,   //[I: what to read]          
			     char         * = 0,   //[I: time range]            
			     bool           = 1);  //[I: warn: missing data]    
    //--------------------------------------------------------------------------
    // write a single snapshot                                                  
    void write_nemo_snapshot(                      // write snapshot            
			     nemo_out const&,      // I: nemo output            
			     const double  * = 0,  //[I: write time]            
			     io              = io::mxv, //[I: what to write]    
			     uint            = 0,  //[I: only write K]          
			     uint            = 0)  //[I: begin with this]       
      const;
    //--------------------------------------------------------------------------
    // write a single particle set                                              
    void write_nemo_particles(                     // write bodies to output    
			      nemo_out const&,     // I: nemo output            
			      const double  * = 0, //[I: write time]            
			      io              = io::mxv, //[I: what to write]   
			      uint            = 0, //[I: only write K]          
			      uint            = 0) //[I: begin with this]       
      const;
#endif                                             // falcON_NEMO               
    //--------------------------------------------------------------------------
    // read simple ascii formatted input                                        
    void read_simple_ascii(                        // read simple ascii file    
			   std::istream  &,        // I: input stream           
			   const io*const&,        // I: array: data items      
			   uint           ,        // I: # total lines          
			   uint          = 0);     // I: # lines with SPH       
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::ebodies                                                      //
  //                                                                          //
  // holds N-body data as a structure of arrays (SoA), very similar to class  //
  // bodies above, except that data management (and I/O) is not done. Class   //
  // ebodies merely serves as an interface for falcON of N-body data hold in  //
  // arrays in external programs, eg in a FORTRAN application callin falcON.  //
  //                                                                          //
  // Two other important differences to bodies are as follows.                //
  // 1. Real valued scalars are hold in a type called "areal", to allow for   //
  //    double when falcON uses float and vice versa.                         //
  // 2. Vector valued quantities (position, velocity, acceleration) are hold  //
  //    in Ndim arrays (one for each dimension).                              //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class ebodies : public bodies_data<1> {
    ebodies& operator=(ebodies const&);            // not provided              
  public:
    //--------------------------------------------------------------------------
    // construction: only set N, no data yet                                    
    //--------------------------------------------------------------------------
    ebodies(uint nb = 0u, uint ns = 0u) : bodies_data<1>(nb,ns) {}
    //--------------------------------------------------------------------------
    // set data                                                                
    //--------------------------------------------------------------------------
    bodies_data<1>::setN;
    void set_mass (areal *D) { set_data_void( 0,static_cast<void*>(D)); }
    void set_pos  (areal**D) { set_data_void( 1,static_cast<void*>(D)); }
    void set_vel  (areal**D) { set_data_void( 2,static_cast<void*>(D)); }
    void set_eps  (areal *D) { set_data_void( 3,static_cast<void*>(D)); }
    void set_flg  (int   *D) { set_data_void( 4,static_cast<void*>(D)); }
    void set_key  (int   *D) { set_data_void( 5,static_cast<void*>(D)); }
    void set_pot  (areal *D) { set_data_void( 6,static_cast<void*>(D)); }
    void set_pex  (areal *D) { set_data_void( 7,static_cast<void*>(D)); }
    void set_acc  (areal**D) { set_data_void( 8,static_cast<void*>(D)); }
    void set_rho  (areal *D) { set_data_void( 9,static_cast<void*>(D)); }
    void set_aux  (areal *D) { set_data_void(10,static_cast<void*>(D)); }
    void set_level(short *D) { set_data_void(11,static_cast<void*>(D)); }
    void set_num  (int   *D) { set_data_void(12,static_cast<void*>(D)); }
    void set_size (areal *D) { set_data_void(13,static_cast<void*>(D)); }
#ifdef falcON_SPH
    void set_snum (int   *D) { set_data_void(14,static_cast<void*>(D)); }
    void set_uin  (areal *D) { set_data_void(15,static_cast<void*>(D)); }
    void set_uprd (areal *D) { set_data_void(16,static_cast<void*>(D)); }
    void set_udex (areal *D) { set_data_void(17,static_cast<void*>(D)); }
    void set_udin (areal *D) { set_data_void(18,static_cast<void*>(D)); }
    void set_end  (areal *D) { set_data_void(19,static_cast<void*>(D)); }
    void set_srho (areal *D) { set_data_void(20,static_cast<void*>(D)); }
    void set_drho (areal *D) { set_data_void(21,static_cast<void*>(D)); }
    void set_vprd (areal**D) { set_data_void(22,static_cast<void*>(D)); }
    void set_sigq (areal *D) { set_data_void(23,static_cast<void*>(D)); }
    void set_temp (areal *D) { set_data_void(24,static_cast<void*>(D)); }
    void set_hdot (areal *D) { set_data_void(25,static_cast<void*>(D)); }
#endif
    //==========================================================================
    //                                                                          
    // 2. Access to bodies DATA[] via sub-type iterator                         
    //                                                                          
    // Note that we cannot put this into the base class bodies_data, because   
    // this is a template which confuses the definition of friends of iterator.
    //                                                                          
    //--------------------------------------------------------------------------
    class iterator {
      //........................................................................
      // basic stuff (construction, iteration, etc)                             
      //........................................................................
      friend class ebodies;
    private:
      iterator();                                  // not implemented           
      //........................................................................
      // data members                                                           
      //........................................................................
      const ebodies *B;                            // pointer to OBJECTS        
      unsigned       K;                            // index                     
      //........................................................................
      // private types                                                          
      //........................................................................
      enum {
	MARK_BODY   = 1<<9,                        // use  9th bit of flag      
	NOT_LONGER  = 1<<10,                       // use 10th bit of flag      
	NOT_SHORTER = 1<<11                        // use 11th bit of flag      
      };
    protected:
      //........................................................................
      // construction                                                           
      //........................................................................
      iterator(const ebodies*b, const uint k) : B(b), K(k) {}
    public:
      iterator(const iterator&I)              : B(I.B), K(I.K) {}
      //........................................................................
      // forward iteration                                                      
      //........................................................................
      iterator       next ()            const { return iterator(B,K+1); }
      iterator& operator++()                  { ++K; return *this; }
      iterator  operator++(int)               { return iterator(B,K++); }
      iterator& operator+=(const int&k)       { K+=k; return *this; }
      iterator  operator+ (const int&k) const { return iterator(B,K+k); }
      //........................................................................
      // boolean methods                                                        
      //........................................................................
      bool  operator== (const iterator&I) const { return K == I.K; }
      bool  operator!= (const iterator&I) const { return K != I.K; }
      bool  operator<  (const iterator&I) const { return K <  I.K; }
      bool  operator<= (const iterator&I) const { return K <= I.K; }
      bool  operator>  (const iterator&I) const { return K >  I.K; }
      bool  operator>= (const iterator&I) const { return K >= I.K; }
      //........................................................................
      // access to bodies                                                       
      //........................................................................
      const ebodies*const&mybodies() const { return B; }
      //........................................................................
      // data access:                                                           
      //   non-const:  via member methods                                       
      //   const    :  via friends                                              
      //........................................................................
      friend uint const&index(const iterator&I) { return I.K; }
      //........................................................................
      // non-const data access via template taking bit#                         
      template<int BIT>
      typename field_bit<BIT,1>::r_type data_bit() {
	return B-> template datum_bit<BIT>(K);
      }
      //........................................................................
      // non-const data access via template taking io                           
      template<int IO>
      typename field_io<IO,1>::r_type data_io() {
	return B-> template datum_io<IO>(K);
      }
      //........................................................................
      // const data access via template friend taking bit#                      
      template<int BIT> friend
      typename field_bit<BIT,1>::cr_type data_bit(iterator const&I) {
	return I.B-> template const_datum_bit<BIT>(I.K);
      }
      //........................................................................
      // const data access via template friend taking io                        
      template<int IO> friend
      typename field_io<IO,1>::cr_type data_io(iterator const&I) {
	return I.B-> template const_datum_io<IO>(I.K);
      }
      //........................................................................
#define ACCESS(NAME,BIT)					\
             field_bit<BIT,1>::r_type  NAME() {			\
	return iterator::data_bit<BIT>(); }			\
      friend field_bit<BIT,1>::cr_type NAME(iterator const&I) {	\
	return I.B-> const_datum_bit<BIT>(I.K); }
      //........................................................................
      // data access via named member methods (non-const) and friends (const)   
      ACCESS(mass,0);
      ACCESS(pos,1);
      ACCESS(vel,2);
      ACCESS(eps,3);
      ACCESS(flg,4);
      ACCESS(key,5);
      ACCESS(pot,6);
      ACCESS(pex,7);
      ACCESS(acc,8);
      ACCESS(rho,9);
      ACCESS(aux,10);
      ACCESS(level,11);
      ACCESS(num,12);
      ACCESS(size,13);
#ifdef falcON_SPH
      ACCESS(snum,14);
      ACCESS(uin,15);
      ACCESS(uprd,16);
      ACCESS(udex,17);
      ACCESS(udin,18);
      ACCESS(ent,19);
      ACCESS(srho,20);
      ACCESS(drho,21);
      ACCESS(vprd,22);
      ACCESS(sigq,23);
      ACCESS(temp,24);
      ACCESS(hdot,25);
#endif
#undef ACCESS
      friend amom angmon (iterator const&B) {
	return nbdy::pos(B) ^ nbdy::vel(B);
      }
      //........................................................................
      // other const methods                                                    
      //........................................................................
      friend amom angmom(iterator const&I) {
	return nbdy::pos(I) ^ nbdy::vel(I);
      }
      //........................................................................
      // flag manipulations                                                     
      //........................................................................
      bool flag_is_set       (int const&F) const {
	return nbdy::flg(*this).is_set(F); }
      operator const flag&   ()            const {
	return nbdy::flg(*this); }
      void flag_as_active    () { flg().add    (flag::ACTIVE); }
      void unflag_active     () { flg().un_set (flag::ACTIVE); }
      void flag_as_sticky    () { flg().add    (flag::STICKY); }
      void unflag_sticky     () { flg().un_set (flag::STICKY); }
      void flag_as_sph       () { flg().add    (flag::SPH); }
      void unflag_sph        () { flg().un_set (flag::SPH); }
      void allow_tree_usage  () { B->allow_tree_usage(K); }
      void forbid_tree_usage () { B->forbid_tree_usage(K); }
      void allow_longer      () { flg().un_set (NOT_LONGER); }
      void forbid_longer     () { flg().add    (NOT_LONGER); }
      void allow_shorter     () { flg().un_set (NOT_SHORTER); }
      void forbid_shorter    () { flg().add    (NOT_SHORTER); }
      void mark              () { flg().add    (MARK_BODY); }
      void unmark            () { flg().un_set (MARK_BODY); }
      //........................................................................
      // const boolean informations via mambers                                 
      //........................................................................
      bool is_source     () const { return nbdy::mass(*this) != zero; }
      bool may_go_longer () const { return !flag_is_set(NOT_LONGER);}
      bool may_go_shorter() const { return !flag_is_set(NOT_SHORTER); }
      bool is_marked     () const { return  flag_is_set(MARK_BODY); }
      //........................................................................
      // const boolean informations via friends                                 
      //........................................................................
      friend bool is_source     (const iterator&I) { return I.is_source(); }
      friend bool may_go_longer (const iterator&I) { return I.may_go_longer(); }
      friend bool may_go_shorter(const iterator&I) { return I.may_go_shorter();}
      friend bool is_marked     (const iterator&I) { return I.is_marked(); }
      //........................................................................
      // output: give just the index (an iterator acts like a pointer)          
      //........................................................................
      friend std::ostream& operator<<(std::ostream&o, const iterator&I) {
	return o<<nbdy::index(I);
      }
    };
    //--------------------------------------------------------------------------
    // access to iterators                                                      
    iterator body_no     (uint i) const { return iterator(this,i); }
    iterator bodyNo      (uint i) const { return iterator(this,i); }
    iterator begin_bodies()       const { return iterator(this,0); }
    iterator end_bodies  ()       const { return iterator(this,N_bodies()); }
    iterator last_bodies ()       const { return iterator(this,N_bodies()-1); }
    iterator back_bodies ()       const { return iterator(this,N_bodies()-1); }
    iterator begin_SPH   ()       const { return iterator(this,0); }
    iterator end_SPH     ()       const { return iterator(this,N_sph()); }
    iterator back_SPH    ()       const { return iterator(this,N_sph()-1); }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::LoopIO<>                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<template<int IO> class K, int BIT=0, int END=IO_NQUANT>
  struct LoopIO {
    static void loop(bodies::iterator&B) {
      K<1 << BIT>::act_on_body(B);
      LoopIO<K,BIT+1>::loop(B);
    }
    static void loop(bodies*const&B) {
      K<1 << BIT>::act_on_bodies(B);
      LoopIO<K,BIT+1>::loop(B);
    }
    static void const_loop(const bodies*const&B) {
      K<1 << BIT>::act_on_bodies(B);
      LoopIO<K,BIT+1>::const_loop(B);
    }
  };
  template<template<int IO> class K, int BIT> struct LoopIO<K,BIT,BIT> {
    static void loop(bodies::iterator  &) {}
    static void loop(bodies*const&) {}
    static void const_loop(const bodies*const&) {}
  };
  //////////////////////////////////////////////////////////////////////////////
}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// macros for looping bodies                                                    
//                                                                              
////////////////////////////////////////////////////////////////////////////////
#ifndef LoopBodies
#define LoopBodies(BODIES_TYPE,                    /* type of bodies         */\
		   BODIES_PTER,                    /* pointer to bodies      */\
		   NAME)                           /* name for body          */\
  for(register BODIES_TYPE::iterator               /* type of body           */\
      NAME  =(BODIES_PTER)->begin_bodies();        /* from first body        */\
      NAME !=(BODIES_PTER)->end_bodies();          /* until beyond last body */\
    ++NAME)                                        /* get next body          */
#endif
//------------------------------------------------------------------------------
#ifndef LoopBodiesRange
#define LoopBodiesRange(BODIES_TYPE,               /* type of bodies         */\
		        BODIES_PTER,               /* pointer to bodies      */\
		        NAME,                      /* name for body          */\
                        FIRST,                     /* first body             */\
                        END)                       /* end body               */\
  for(register BODIES_TYPE::iterator               /* type of body           */\
      NAME  =(BODIES_PTER)->bodyNo(FIRST);         /* from first body        */\
      NAME !=(END? (BODIES_PTER)->bodyNo(END)  :   /* until end OR           */\
	           (BODIES_PTER)->end_bodies());   /* until beyond last body */\
    ++NAME)                                        /* get next body          */
#endif
//------------------------------------------------------------------------------
#ifndef LoopBodyPairs
#define LoopBodyPairs(BODIES_TYPE,                 /* type of bodies         */\
		      BODIES_PTER,                 /* pointer to bodies      */\
		      NAME1,                       /* name of  1st body      */\
		      NAME2)                       /* name for 2nd body      */\
  for(register BODIES_TYPE::iterator               /* type of body           */\
      NAME2  = NAME1 + 1;                          /* from first body        */\
      NAME2 !=(BODIES_PTER)->end_bodies();         /* until beyond last body */\
    ++NAME2)                                       /* get next body          */
#endif
//------------------------------------------------------------------------------
#ifndef LoopSPHBodies
#define LoopSPHBodies(BODIES_TYPE,                 /* type of bodies         */\
		      BODIES_PTER,                 /* pointer to bodies      */\
		      NAME)                        /* name for body          */\
  for(register BODIES_TYPE::iterator               /* type of body           */\
      NAME  =(BODIES_PTER)->begin_SPH();           /* from first SPH body    */\
      NAME !=(BODIES_PTER)->end_SPH();             /* until last SPH body    */\
    ++NAME)                                        /* get next body          */
#endif
//------------------------------------------------------------------------------
#ifndef LoopNonSPHBodies
#define LoopNonSPHBodies(BODIES_TYPE,              /* type of bodies         */\
		         BODIES_PTER,              /* pointer to bodies      */\
		         NAME)                     /* name for body          */\
  for(register BODIES_TYPE::iterator               /* type of body           */\
      NAME  =(BODIES_PTER)->end_SPH();             /* from first non-SPH body*/\
      NAME !=(BODIES_PTER)->end_bodies();          /* until beyond last body */\
    ++NAME)                                        /* get next body          */
#endif
//------------------------------------------------------------------------------
#endif                                             // falcON_included_body_h    
