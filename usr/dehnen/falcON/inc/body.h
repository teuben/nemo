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
//                                                                             |
// defines                                                                     |
//                                                                             |
// namespace nbdy {                                                            |
//   class sbodies;                                                            |
//   class sbodies::iterator;                                                  |
//   typedef sbodies::iterator body;                                           |
//   template struct nbody_io<>;                                               |
//   template struct LoopIO<>;                                                 |
//   class abodies;                                                            |
// }                                                                           |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_body_h
#define falcON_included_body_h 1

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
#if defined(falcON_NEMO) && !defined(falcON_included_nmio_h)
  class nemo_in;                                   // forward declaration       
  class nemo_out;                                  // forward declaration       
#endif
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::BodyComm                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class BodyComm {
    //--------------------------------------------------------------------------
  private:
    uint  NBODIES, NSPH;
    io    ALLBITS;
    //--------------------------------------------------------------------------
    // construction/destruction and such                                        
    //--------------------------------------------------------------------------
  protected:
    BodyComm(uint const&nb,                        // I: N_bodies               
	     uint const&ns,                        // I: N_sph                  
	     io   const&b ) :                      // I: all data bits          
      NBODIES( nb ), NSPH( min(ns,nb) ), ALLBITS( b ) {}
    //--------------------------------------------------------------------------
    BodyComm(BodyComm const&B) : 
      NBODIES(B.NBODIES), ALLBITS(B.ALLBITS) {}
    //--------------------------------------------------------------------------
    BodyComm& operator= (BodyComm const&B) 
    { NBODIES = B.NBODIES; ALLBITS=B.ALLBITS; return *this; }
    //--------------------------------------------------------------------------
    void reset(uint const&nb,                      // I: N_bodies               
	       uint const&ns,                      // I: N_sph                  
	       io   const&b )                      // I: all data bits          
    { NBODIES = nb; NSPH = ns; ALLBITS = b; }
    //--------------------------------------------------------------------------
    void add_all_bits(io const&b) { ALLBITS |= b; }
    void change_bits (io const&mask, io const&b) {
      ALLBITS = (b & mask) | (ALLBITS & ~mask);
    }
    //--------------------------------------------------------------------------
    // const methods                                                            
    //--------------------------------------------------------------------------
  public:
    const uint&N_bodies    ()           const { return NBODIES; }
    const uint&N_sph       ()           const { return NSPH; }
    const io  &all_bits    ()           const { return ALLBITS; }
    bool  all_bits_contain (io const&b) const { return b == (b & ALLBITS); }
    bool  has              (io const&b) const { return b == (b & ALLBITS); }
    io    has_not          (io const&b) const { return b & ~ALLBITS; }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // structs defining the data holdable in the various data classes           //
  //                                                                          //
  // NOTE due to a compiler bug (gcc 3.2.2), we must define the bits below as //
  //      static functions rather than static const. Otherwise, the compiler  //
  //      spits references to them out, which the linker cannot resolve.      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  struct BodySrceBits {
    static const int MINBITS_  = io::mx;
    static const int MAXBITS_  = io::mxv|io::f|io::e|io::k;
    static const int DEFBITS_  = io::mxv|io::f;
    static const int NEMOBITS_ = MAXBITS_;
    static int MAXBITS () { return MAXBITS_; }
    static int MINBITS () { return MINBITS_; }
    static int DEFBITS () { return DEFBITS_; }
    static int NEMOBITS() { return NEMOBITS_; }
    static int data_bits(io const&b) { return MINBITS() | (b & MAXBITS()); }
    static int part_bits(io const&b) { return b & MAXBITS(); }
  };
  //////////////////////////////////////////////////////////////////////////////
  struct BodySinkBits {
    static const int MINBITS_  = io::a|io::p;
    static const int MAXBITS_  = io::a|io::p|io::q|io::r|io::y|io::n|io::l;
    static const int DEFBITS_  = MINBITS_;
    static const int NEMOBITS_ = MINBITS_ |io::r|io::y|io::l|io::n;
    static int MAXBITS () { return MAXBITS_; }
    static int MINBITS () { return MINBITS_; }
    static int DEFBITS () { return DEFBITS_; }
    static int NEMOBITS() { return NEMOBITS_; }
    static int data_bits(io const&b) { return MINBITS() | (b & MAXBITS()); }
    static int part_bits(io const&b) { return b & MAXBITS(); }
  };
  //////////////////////////////////////////////////////////////////////////////
  struct BodyPsphBits {
    static const int MINBITS_  = io::sphmin;
    static const int MAXBITS_  = io::sphmax;
    static const int DEFBITS_  = io::sphmin;
    static const int NEMOBITS_ = io::sphnemo;
    static int MAXBITS () { return MAXBITS_; }
    static int MINBITS () { return MINBITS_; }
    static int DEFBITS () { return DEFBITS_; }
    static int NEMOBITS() { return NEMOBITS_; }
    static int data_bits(io const&b) { return MINBITS() | (b & MAXBITS()); }
    static int part_bits(io const&b) { return b & MAXBITS(); }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::BodyDataBase                                                 //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename BodyBits>
  class BodyDataBase {
    BodyDataBase           (BodyDataBase const&);  // not implemented           
    BodyDataBase& operator=(BodyDataBase const&);  // not implemented           
    //--------------------------------------------------------------------------
    // data members                                                             
    //--------------------------------------------------------------------------
  protected:
    //--------------------------------------------------------------------------
    io           DATABITS;                         // which data do we hold?    
    uint         BBYTES;                           // bytes per body            
    uint         NB;                               // # bodies                  
    char        *ALLOC;                            // data                      
    //--------------------------------------------------------------------------
    // protected member methods                                                 
    //--------------------------------------------------------------------------
    void allocate(uint n = 0) {                    //[I: N_bodies]              
      NB    = n;
      ALLOC = (NB && BBYTES) ? new char[NB*BBYTES] : 0;
    }
    //--------------------------------------------------------------------------
    void set_data_bits(io const&b) {
      DATABITS = BodyBits::data_bits(b);
      BBYTES   = DATABITS.bytes();
    }
    //--------------------------------------------------------------------------
    BodyDataBase(uint const&n,                     // I: N_particles            
		 io   const&b) {                   // I: data bits              
      set_data_bits(b);                            // order important, 'cause   
      allocate     (n);                            // BBYTES needed for allocate
    }
    //--------------------------------------------------------------------------
    ~BodyDataBase() {                              // destructor                
      if(ALLOC) delete[] ALLOC;                    //   delete memory           
    }
    //--------------------------------------------------------------------------
    bool reset  (uint const&n,                     // I: N_bodies               
		 io   const&b) {                   // I: data bits              
      if(n != NB || BodyBits::data_bits(b) != DATABITS) {
	if(ALLOC) delete[] ALLOC;
	set_data_bits(b);                          // order important, 'cause   
	allocate     (n);                          // BBYTES needed for allocate
	return true;
      } else
	return false;
    }
    //--------------------------------------------------------------------------
    // public member methods                                                    
    //--------------------------------------------------------------------------
  public:
    const io  &my_databits () const { return DATABITS; }
    //--------------------------------------------------------------------------
    bool is_supported(const io&b) const
    {
      return DATABITS.contains(BodyBits::part_bits(b));
    }
  };
  //////////////////////////////////////////////////////////////////////////////
#define SETPTER(TYPE,NAME,BIT)						       \
    if(DATABITS & BIT) { NAME = static_cast<TYPE*>(static_cast<void*>(M));     \
                         M   += sizeof(TYPE) * NB; }			       \
    else 		 NAME = 0;
  //////////////////////////////////////////////////////////////////////////////
#define ACCESS(NAME,TYPE,ARRAY)						\
    TYPE      &NAME(const uint&i) const { return ARRAY[i]; }		\
    TYPE*const&NAME##_s() const { return ARRAY; } 
#define SUBACCESS(NAME,TYPE,ARRAY)					\
    TYPE &NAME(const uint&i, const int&d) const { return ARRAY[i][d]; }
#ifdef DEBUG
#  define ACCESS_CHECK(NAME,TYPE,ARRAY)					\
    TYPE &NAME(const uint&i) const {					\
      if(ARRAY) return ARRAY[i];					\
      else falcON_Error("data wanted not allocated");			\
    }									\
    TYPE*const&NAME##_s() const { return ARRAY; } 
#  define SUBACCESS_CHECK(NAME,TYPE,ARRAY)				\
    TYPE &NAME(const uint&i, const int&d) const { 			\
      if(ARRAY) return ARRAY[i][d];  					\
      else falcON_Error("data wanted not allocated");			\
    }
#else
#  define ACCESS_CHECK(NAME,TYPE,ARRAY)    ACCESS(NAME,TYPE,ARRAY)
#  define SUBACCESS_CHECK(NAME,TYPE,ARRAY) SUBACCESS(NAME,TYPE,ARRAY)
#endif
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::BodySrce                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // - holds those data that every gravity source particle needs              //
  // - these are mass, position, velocity [, flag, eps, key]                  //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class BodySrce :
    public BodySrceBits,
    public BodyComm,
    public BodyDataBase<BodySrceBits> {
    //--------------------------------------------------------------------------
    BodySrce           (const BodySrce&);          // not implemented           
    BodySrce& operator=(const BodySrce&);          // not implemented           
    //--------------------------------------------------------------------------
    mutable bool CUSE;                             // flag for tree usage change
    mutable uint NA;                               // # bodies in tree          
  protected:
    real  *MAS;                                    // pointer to masses         
    vect  *POS;                                    // pointer to positions      
    vect  *VEL;                                    // pointer to velocities     
    real  *EPS;                                    // pointer to eps_i          
    flag  *FLG;                                    // pointer to flags          
    int   *KEY;                                    // pointer to keyes          
  private:
    mutable bool   DATA_CHANGED;
    //--------------------------------------------------------------------------
    void set_pointers() {
      if(ALLOC) {
	register char* M = ALLOC;
	SETPTER(real,MAS,io::m)
	SETPTER(vect,POS,io::x)
	SETPTER(vect,VEL,io::v)
	SETPTER(real,EPS,io::e)
	SETPTER(flag,FLG,io::f)
	SETPTER(int, KEY,io::k)
      } else {
	MAS  = 0;
	POS  = 0;
	VEL  = 0;
	EPS  = 0;
	FLG  = 0;
	KEY  = 0;
      }
    }
    //--------------------------------------------------------------------------
  public:
    BodySrce(const uint nb,                        // I: # bodies               
	     const io   b ,                        // I: data bits              
	     const uint ns = 0,                    //[I: # SPH particles]       
	     const uint na = 0,                    //[I: # bodies in tree]      
	     const bool cu = false) :              //[I: CUSE flag]             
      BodyComm                   ( nb,ns,b ),
      BodyDataBase<BodySrceBits> ( nb,b ),
      CUSE                       ( cu ),
      NA                         ( na > 0? na : NB ),
      DATA_CHANGED               ( 1 )
    {
      set_pointers();                              //   set pointers to data    
    }
    //--------------------------------------------------------------------------
    void reset(const uint&nb,                      // I: # bodies               
	       const io  &b,                       // I: data bits              
	       const uint ns = 0) {                //[I: # SPH particles]       
      BodyComm::reset(nb,ns,b);
      if(BodyDataBase<BodySrceBits>::reset(nb,b)) {// IF(memory reset)          
	set_pointers();                            //   set pointers to data    
	if(FLG) reset_flags ();                    //   reset flags             
      }                                            // ENDIF                     
      DATA_CHANGED = 1;
    }
    //--------------------------------------------------------------------------
    // further member methods                                                   
    //--------------------------------------------------------------------------
    const uint &N_intree  ()             const { return NA; }
    bool        has_flg   ()             const { return FLG!=0; }
    bool        has_flag  ()             const { return FLG!=0; }
    bool        has_eps   ()             const { return EPS!=0; }
    bool        has_key   ()             const { return KEY!=0; }
    ACCESS      (pos,vect,POS) SUBACCESS(pos,real,POS)
    ACCESS      (vel,vect,VEL) SUBACCESS(vel,real,VEL)
    ACCESS      (mass,real,MAS)
    ACCESS_CHECK(flg,flag,FLG)
    ACCESS_CHECK(eps,real,EPS)
    ACCESS_CHECK(key,int, KEY)
    amom        angmom    (const uint&i) const { return POS[i] ^ VEL[i]; }
#ifdef falcON_NEMO
    //--------------------------------------------------------------------------
    // NEMO Input                                                               
    // to be used by bodies; does not reset(); reset() needs to be called before
    //--------------------------------------------------------------------------
    io read_nemo(                                  // R: data read              
		 nemo_in const&,                   // I: nemo input             
		 io      const& =MINBITS(),        //[I: what to read]          
		 bool    const& =true,             //[I: warn upon missing data]
		 size_t  const& =0) const;         //[I: first body to get]     
#endif
    //--------------------------------------------------------------------------
    // flag manipulations                                                       
    //--------------------------------------------------------------------------
  public:
    void reset_flags() const {
      if(FLG) {
	for(register uint i=0; i!=NB; ++i) FLG[i].reset();
	NA   = NB;
	CUSE = false;
      } else
	falcON_WarningF("flags not supported","sbodies::reset_flags()");
    }
    //--------------------------------------------------------------------------
    void flag_all_as_active() const {              // flag all bodies as active 
      if(FLG) for(register uint i=0; i!=NB; ++i) FLG[i].add(flag::ACTIVE);
      else falcON_ErrorF("flags not supported","sbodies::flag_all_as_active()");
    }
    //--------------------------------------------------------------------------
    void flag_as_sph(uint const&ns) const {        // flag the first ns as sph  
      const uint Ns = min(ns,NB);
      if(FLG) for(register uint i=0; i!=Ns; ++i) FLG[i].add(flag::SPH);
      else falcON_ErrorF("flags not supported","sbodies::flag_as_sph()");
    }
    //--------------------------------------------------------------------------
    void allow_tree_usage  (const uint&i) const {
#ifdef DEBUG
      if(FLG) {
#endif
      if(! is_in_tree(FLG[i])) {                   // IF(body not in tree)     >
	FLG[i].un_set(flag::NOT_IN_TREE);          //   unflag 'not in tree'    
	NA++;                                      //   increment N_intree      
	CUSE = true;                               //   set flag for change     
      }                                            // <                         
#ifdef DEBUG
      } else falcON_ErrorF("flags not supported","sbodies::allow_tree_usage()");
#endif
    }
    //--------------------------------------------------------------------------
    void forbid_tree_usage (const uint&i) const {
#ifdef DEBUG
      if(FLG) {
#endif
      if(is_in_tree(FLG[i])) {                     // IF(body was in tree)     >
	FLG[i].add(flag::NOT_IN_TREE);             //   flag 'not in tree'      
	NA--;                                      //   decrement N_intree      
	CUSE = true;                               //   set flag for change     
      }                                            // <                         
#ifdef DEBUG
      } else falcON_ErrorF("flags not supported",
			   "sbodies::forbid_tree_usage()");
#endif
    }
    //--------------------------------------------------------------------------
    void set_mempter(char ** const&MEM) const
    {
      register char* M=ALLOC;
      for(register int b=1<<io::P_srce(), i=0; i!=io::N_srce(); ++i, b<<=1) {
	if(DATABITS & b) {
	  MEM[i] = M;
	  M     += io::Z_quant(i);
	} else
	  MEM[i] = 0;
      }
    }
    //--------------------------------------------------------------------------
    bool changes_in_tree_usage_flags() const { return CUSE; }
    void after_tree_growth()           const { CUSE = false;  }
    //--------------------------------------------------------------------------
    void mark_srce_data_read   () const { DATA_CHANGED = 0; }
    void mark_srce_data_changed() const { DATA_CHANGED = 1; }
    bool srce_data_changed     () const { return DATA_CHANGED; }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class nbdy::BodySink                                                       
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // - holds those data that need not to be initialized                         
  // - these are acceleration,potential,[external potential,density,aux,level]  
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class BodySink : 
    public BodySinkBits,
    public BodyDataBase<BodySinkBits>
  {
    BodySink           (const BodySink&);          // not implemented           
    BodySink& operator=(const BodySink&);          // not implemented           
    //--------------------------------------------------------------------------
    // data members                                                             
    //--------------------------------------------------------------------------
  protected:
    real        *POT;                              // pointer to N-body pots    
    real        *PEX;                              // pointer to external pots  
    vect        *ACC;                              // pointer to acceleration   
    real        *RHO;                              // pointer to densities      
    real        *AUX;                              // pointer to aux data       
    indx        *LEV;                              // pointer to levels         
    uint        *NUM;                              // pointer to # neighbours   
  private:
    //--------------------------------------------------------------------------
    void set_pointers() {
      if(ALLOC) {
	register char* M = ALLOC;
	SETPTER(real,POT,io::p)
	SETPTER(real,PEX,io::q)
	SETPTER(vect,ACC,io::a)
	SETPTER(real,RHO,io::r)
	SETPTER(real,AUX,io::y)
	SETPTER(indx,LEV,io::l)
	SETPTER(uint,NUM,io::n)
      } else {
	POT  = 0;
	PEX  = 0;
	ACC  = 0;
	RHO  = 0;
	AUX  = 0;
	LEV  = 0;
	NUM  = 0;
      }
    }
    //--------------------------------------------------------------------------
  protected:
    BodySink(const uint n = 0,                     //[I: N_bodies]              
	     const io   b = DEFBITS()) :           //[I: data bits]             
      BodyDataBase<BodySinkBits> ( n,b )
    {
      set_pointers();                              // set pointers to data      
    }
    //--------------------------------------------------------------------------
    void reset  (const uint n,                     // I: N_bodies               
		 const io   b)                     // I: data bits              
    {
      if(BodyDataBase<BodySinkBits>::reset(n,b))   // IF(memory reset)          
	set_pointers();                            //   set pointers to data    
    }
#ifdef falcON_NEMO
    //--------------------------------------------------------------------------
    // NEMO Input                                                               
    // to be used by bodies; does not reset(); reset() needs to be called before
    //--------------------------------------------------------------------------
  protected:
    io read_nemo(                                  // R: data read              
		 nemo_in const&,                   // I: nemo input             
		 io      const& =MINBITS(),        //[I: what to read]          
		 bool    const& =true,             //[I: warn upon missing data]
		 size_t  const& =0) const;         //[I: first body to get]     
#endif
    //--------------------------------------------------------------------------
    // const member methods                                                     
    //--------------------------------------------------------------------------
  public:
    ACCESS      (acc,vect,ACC)
    ACCESS      (pot,real,POT)
    ACCESS_CHECK(pex,real,PEX)
    ACCESS_CHECK(aux,real,AUX)
    ACCESS_CHECK(rho,real,RHO)
    ACCESS_CHECK(num,uint,NUM)
    ACCESS_CHECK(level,indx,LEV)
    //--------------------------------------------------------------------------
    void set_mempter(char ** const&MEM) const
    {
      register char* M=ALLOC;
      for(register int b=1<<io::P_sink(), i=0; i!=io::N_sink(); ++i, b<<=1) {
	if(DATABITS & b) {
	  MEM[i] = M;
	  M     += io::Z_quant(i);
	} else
	  MEM[i] = 0;
      }
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class nbdy::BodyPsph                                                       
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // - holds those data that belong to SPH particles only                       
  // - these are size [,temperature,inner energy]                               
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class BodyPsph : 
    public BodyPsphBits,
    public BodyDataBase<BodyPsphBits> {
    BodyPsph           (const BodyPsph&);          // not implemented           
    BodyPsph& operator=(const BodyPsph&);          // not implemented           
    //--------------------------------------------------------------------------
    // data members                                                             
    //--------------------------------------------------------------------------
  protected:
    real        *SIZ;                              // pointer to SPH sizes      
#ifdef falcON_SPH
    uint        *NSP;                              // pointer to # SPH partners 
    real        *UIN;                              // pointer to SPH U          
    real        *UPR;                              // pointer to SPH U_predicted
    real        *UDI;                              // pointer to SPH dU/dt_in   
    real        *UDE;                              // pointer to SPH dU/dt_ex   
    real        *ENT;                              // pointer to SPH entropies  
    real        *SRH;                              // pointer to SPH densities  
    real        *DRH;                              // pointer to SPH drho/dt    
    vect        *VPR;                              // pointer to SPH v_predicted
    real        *SIQ;                              // pointer to SPH sigma^2    
    real        *TEM;                              // pointer to SPH Temperature
    real        *DHT;                              // pointer to SPH dsize/dt   
#endif
  private:
    mutable bool DATA_CHANGED;
    //--------------------------------------------------------------------------
    void set_pointers() {
      if(ALLOC) {
	register char* M = ALLOC;
	SETPTER(real,SIZ,io::H)
#ifdef falcON_SPH
	SETPTER(uint,NSP,io::N)
	SETPTER(real,UIN,io::U)
	SETPTER(real,UPR,io::Y)
	SETPTER(real,UDI,io::I)
	SETPTER(real,UDE,io::E)
	SETPTER(real,ENT,io::S)
	SETPTER(real,SRH,io::R)
	SETPTER(real,DRH,io::D)
	SETPTER(vect,VPR,io::V)
	SETPTER(real,SIQ,io::Q)
	SETPTER(real,TEM,io::T)
	SETPTER(real,DHT,io::J)
#endif
      } else {
	SIZ   = 0;
#ifdef falcON_SPH
	NSP   = 0;
	UIN   = 0;
	UPR   = 0;
	UDI   = 0;
	UDE   = 0;
	ENT   = 0;
	SRH   = 0;
	DRH   = 0;
	VPR   = 0;
	SIQ   = 0;
	TEM   = 0;
	DHT   = 0;
#endif
      }
    }
    //--------------------------------------------------------------------------
  protected:
    BodyPsph(const uint n = 0,                     //[I: N_sph]                 
	     const io   b = DEFBITS()) :           //[I: data bits]             
      BodyDataBase<BodyPsphBits> ( n,b ),
      DATA_CHANGED               ( 1 )
    {
      set_pointers();                              // set pointers to data      
    }
    //--------------------------------------------------------------------------
    void reset  (const uint n,                     // I: N_sph                  
		 const io   b)                     // I: data bits              
    {
      if(BodyDataBase<BodyPsphBits>::reset(n,b))   // IF(memory reset)          
	set_pointers();                            //   set pointers to data    
      DATA_CHANGED = 1;
    }
#ifdef falcON_NEMO
    //--------------------------------------------------------------------------
    // NEMO Input                                                               
    // to be used by bodies; does not reset(); reset() needs to be called before
    //--------------------------------------------------------------------------
  protected:
    io read_nemo(                                  // R: data read              
		 nemo_in const&,                   // I: nemo input             
		 io      const& =MINBITS(),        //[I: what to read]          
		 bool    const& =true,             //[I: warn upon missing data]
		 size_t  const& =0) const;         //[I: first body to get]     
#endif
    //--------------------------------------------------------------------------
    // const member methods                                                     
    //--------------------------------------------------------------------------
  public:
    ACCESS      (size,real,SIZ)
#ifdef falcON_SPH
    ACCESS      (srho,real,SRH)
    ACCESS_CHECK(snum,uint,NSP)
    ACCESS_CHECK(uin, real,UIN)
    ACCESS_CHECK(uprd,real,UPR)
    ACCESS_CHECK(udin,real,UDI)
    ACCESS_CHECK(udex,real,UDE)
    ACCESS_CHECK(ent, real,ENT)
    ACCESS_CHECK(drho,real,DRH)
    ACCESS_CHECK(vprd,vect,VPR)
    ACCESS_CHECK(sigq,real,SIQ)
    ACCESS_CHECK(temp,real,TEM)
    ACCESS_CHECK(hdot,real,DHT)
#endif
    //--------------------------------------------------------------------------
    void set_mempter(char ** const&MEM) const
    {
      register char* M=ALLOC;
      for(register int b=1<<io::P_psph(), i=0; i!=io::N_psph(); ++i, b<<=1) {
	if(DATABITS & b) {
	  MEM[i] = M;
	  M     += io::Z_quant(i);
	} else
	  MEM[i] = 0;
      }
    }
    //--------------------------------------------------------------------------
    void mark_sph_data_read   () const { DATA_CHANGED = 0; }
    void mark_sph_data_changed() const { DATA_CHANGED = 1; }
    bool sph_data_changed     () const { return DATA_CHANGED; }
  };
  //////////////////////////////////////////////////////////////////////////////
#undef ACCESS
#undef ACCESS_CHECK
#undef SUBACCESS
#undef SUBACCESS_CHECK
#undef SETPTER
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class nbdy::sbodies                                                        
  //                                                                            
  // bodies for serial code                                                     
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // - holds all body data                                                      
  // - provides iterator 'body' for sequential access to all bodies             
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class sbodies :
    public BodySrce,                               // bits,N,x,m,v,[f,e,k]      
    public BodySink,                               // p,a,[P,r,y,l]             
    public BodyPsph {                              // [s]                       
  private:
    sbodies           (const sbodies&);            // not implemented           
    sbodies& operator=(const sbodies&);            // not implemented           
    //--------------------------------------------------------------------------
    // data member: pointers to memory associated with data arrays              
    //--------------------------------------------------------------------------
    char * MEMORY[io::N_quant];
    //--------------------------------------------------------------------------
    // static constants                                                         
    //                                                                          
    // NOTE: see note above before definition of BodySrceBits.                  
    //--------------------------------------------------------------------------
  public:
    static int MINBITS() { return 
			     BodySrceBits::MINBITS() | 
			     BodySinkBits::MINBITS() | 
			     BodyPsphBits::MINBITS(); }
    //--------------------------------------------------------------------------
    static int MAXBITS() { return 
			     BodySrceBits::MAXBITS() | 
			     BodySinkBits::MAXBITS() | 
			     BodyPsphBits::MAXBITS(); }
    //--------------------------------------------------------------------------
    static int DEFBITS() { return 
			     BodySrceBits::DEFBITS() | 
			     BodySinkBits::DEFBITS() | 
			     BodyPsphBits::DEFBITS(); }
    //--------------------------------------------------------------------------
    static int NEMOBITS() { return 
			     BodySrceBits::NEMOBITS() | 
			     BodySinkBits::NEMOBITS() | 
			     BodyPsphBits::NEMOBITS(); }
    //--------------------------------------------------------------------------
    // construction                                                             
    //--------------------------------------------------------------------------
    // The second argument specifies the data that shall be supported. This can 
    // subsequently only be reduced by a call to reset() [or reset_sink() and   
    // reset_sph()], but not by any of the I/O. If only the first argument is   
    // given, then it is ensured that subsequent I/O will allocate the specified
    // data (and, possibly, additional ones that are to be read in).            
    // For example, after                                                       
    //                                                                          
    //   sbodies BB(0,io::mxv | io::y);                                         
    //   nemo_in input("file.snp");                                             
    //   io      read;                                                          
    //   real    time;                                                          
    //   BB.read_nemo_snapshot(input,read,&time,io::mxve);                      
    //                                                                          
    // the number N of bodies as well as N masses, positions, velocities and    
    // eps have been read from file "file.snp". Moreover, N aux data have been  
    // allocated (but not initialized), as indicated in the constructor.        
    //                                                                          
    // NOTE, however, that the call to read_nemo_snapshot() may well have       
    // deleted previous data, even if these were not read in, like the aux data 
    // in the above example. This is because the number of bodies has changed   
    // and all data arrays have had to be re-allocated.                         
    //--------------------------------------------------------------------------
  public:
    sbodies(const uint nb = 0,                     //[I: # bodies]              
	    const io   b  = DEFBITS(),             //[I: data bits]             
	    const uint ns = 0,                     //[I: # SPH particles]       
	    const uint na = 0,                     //[I: # bodies in tree]      
	    const bool cu = false) :               //[I: CUSE flag]             
      BodySrce(nb,ns? b|io::f:b,ns,na,cu),
      BodySink(N_bodies(),b),
      BodyPsph(N_sph(),b)
    {
      BodySrce::set_mempter(MEMORY+io::P_srce());
      BodySink::set_mempter(MEMORY+io::P_sink());
      BodyPsph::set_mempter(MEMORY+io::P_psph());
      if(N_sph()) flag_as_sph(N_sph());
    }
    //--------------------------------------------------------------------------
    // enum nbdy::sbodies::part                                                 
    //--------------------------------------------------------------------------
    enum part { Srce=1, Sink=2, Psph=4 };
    //--------------------------------------------------------------------------
    // reset partial data bits, not N                                           
    //--------------------------------------------------------------------------
    void reset(const part p,                       // I: what to reset          
	       const io  &b = DEFBITS())           // I: data bits              
    {
      switch(p) {
      case Srce:
	BodySrce::reset(N_bodies(),b); 
	change_bits(BodySrce::MAXBITS(),b);
	break;
      case Sink:
	BodySink::reset(N_bodies(),b);
	change_bits(BodySink::MAXBITS(),b);
	break;
      case Psph:
	BodyPsph::reset(N_sph   (),b);
	change_bits(BodyPsph::MAXBITS(),b);
	break;
      }
      BodySrce::set_mempter(MEMORY+io::P_srce());
      BodySink::set_mempter(MEMORY+io::P_sink());
      BodyPsph::set_mempter(MEMORY+io::P_psph());
    }
    //--------------------------------------------------------------------------
    // reset N_bodies, data bits, [and N_sph]  deletes all data                 
    //--------------------------------------------------------------------------
    void reset(const uint&nb,                      // I: # bodies               
	       const io  &b,                       // I: data bits              
	       const uint ns = 0)                  //[I: # SPH particles]       
    {
      BodySrce::reset(nb,b,ns);
      BodySink::reset(nb,b);
      BodyPsph::reset(ns,b);
      BodySrce::set_mempter(MEMORY+io::P_srce());
      BodySink::set_mempter(MEMORY+io::P_sink());
      BodyPsph::set_mempter(MEMORY+io::P_psph());
    }
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    // copy a body n times                                                      
    //--------------------------------------------------------------------------
    void copy_body(sbodies const& B,
		   int     const& fr,
		   int     const& to,
		   int     const& n = 0)
    {
      if(n) {
	if(B.has(io::x) && has(io::x))
	  for(register int i=0; i!=n; ++i)  pos (to+i) = B.pos (fr);
	if(B.has(io::m) && has(io::m))
	  for(register int i=0; i!=n; ++i)  mass(to+i) = B.mass(fr);
	if(B.has(io::v) && has(io::v))
	  for(register int i=0; i!=n; ++i)  vel (to+i) = B.vel (fr);
	if(B.has(io::e) && has(io::e))
	  for(register int i=0; i!=n; ++i)  eps (to+i) = B.eps (fr);
	if(B.has(io::k) && has(io::k))
	  for(register int i=0; i!=n; ++i)  key (to+i) = B.key (fr);
	if(B.has(io::p) && has(io::p))
	  for(register int i=0; i!=n; ++i)  pot (to+i) = B.pot (fr);
	if(B.has(io::q) && has(io::q))
	  for(register int i=0; i!=n; ++i)  pex (to+i) = B.pex (fr);
	if(B.has(io::a) && has(io::a))
	  for(register int i=0; i!=n; ++i)  acc (to+i) = B.acc (fr);
	if(B.has(io::f) && has(io::f))
	  for(register int i=0; i!=n; ++i)  flg (to+i) = B.flg (fr);
	if(B.has(io::l) && has(io::l))
	  for(register int i=0; i!=n; ++i)  level(to+i)= B.level(fr);
	if(B.has(io::r) && has(io::r))
	  for(register int i=0; i!=n; ++i)  rho (to+i) = B.rho (fr);
	if(B.has(io::y) && has(io::y))
	  for(register int i=0; i!=n; ++i)  aux (to+i) = B.aux (fr);
	if(B.has(io::n) && has(io::n))
	  for(register int i=0; i!=n; ++i)  num (to+i) = B.num (fr);
	if(B.has(io::H) && has(io::H) && to+n<N_sph())
	  for(register int i=0; i!=n; ++i)  size(to+i) = B.size(fr);
#ifdef falcON_SPH
	if(B.has(io::N) && has(io::N) && to+n<N_sph()) 
	  for(register int i=0; i!=n; ++i)  snum(to+i) = B.snum(fr);
	if(B.has(io::U) && has(io::U) && to+n<N_sph()) 
	  for(register int i=0; i!=n; ++i)  uin (to+i) = B.uin (fr);
	if(B.has(io::Y) && has(io::Y) && to+n<N_sph()) 
	  for(register int i=0; i!=n; ++i)  uprd(to+i) = B.uprd(fr);
	if(B.has(io::I) && has(io::I) && to+n<N_sph()) 
	  for(register int i=0; i!=n; ++i)  udin(to+i) = B.udin(fr);
	if(B.has(io::E) && has(io::E) && to+n<N_sph()) 
	  for(register int i=0; i!=n; ++i)  udex(to+i) = B.udex(fr);
	if(B.has(io::S) && has(io::S) && to+n<N_sph()) 
	  for(register int i=0; i!=n; ++i)  ent (to+i) = B.ent (fr);
	if(B.has(io::R) && has(io::R) && to+n<N_sph()) 
	  for(register int i=0; i!=n; ++i)  srho(to+i) = B.srho(fr);
	if(B.has(io::D) && has(io::D) && to+n<N_sph()) 
	  for(register int i=0; i!=n; ++i)  drho(to+i) = B.drho(fr);
	if(B.has(io::V) && has(io::V) && to+n<N_sph()) 
	  for(register int i=0; i!=n; ++i)  vprd(to+i) = B.vprd(fr);
	if(B.has(io::Q) && has(io::Q) && to+n<N_sph()) 
	  for(register int i=0; i!=n; ++i)  sigq(to+i) = B.sigq(fr);
	if(B.has(io::T) && has(io::T) && to+n<N_sph()) 
	  for(register int i=0; i!=n; ++i)  temp(to+i) = B.temp(fr);
	if(B.has(io::J) && has(io::J) && to+n<N_sph()) 
	  for(register int i=0; i!=n; ++i)  hdot(to+i) = B.hdot(fr);
#endif
      } else {
	if(B.has(io::x) && has(io::x)) pos (to) = B.pos (fr);
	if(B.has(io::m) && has(io::m)) mass(to) = B.mass(fr);
	if(B.has(io::v) && has(io::v)) vel (to) = B.vel (fr);
	if(B.has(io::e) && has(io::e)) eps (to) = B.eps (fr);
	if(B.has(io::k) && has(io::k)) key (to) = B.key (fr);
	if(B.has(io::p) && has(io::p)) pot (to) = B.pot (fr);
	if(B.has(io::q) && has(io::q)) pex (to) = B.pex (fr);
	if(B.has(io::a) && has(io::a)) acc (to) = B.acc (fr);
	if(B.has(io::f) && has(io::f)) flg (to) = B.flg (fr);
	if(B.has(io::l) && has(io::l)) level(to)= B.level(fr);
	if(B.has(io::r) && has(io::r)) rho (to) = B.rho (fr);
	if(B.has(io::y) && has(io::y)) aux (to) = B.aux (fr);
	if(B.has(io::n) && has(io::n)) num (to) = B.num (fr);
	if(B.has(io::H) && has(io::H) && to<N_sph()) size(to) = B.size(fr);
#ifdef falcON_SPH
	if(B.has(io::N) && has(io::N) && to<N_sph()) snum(to) = B.snum(fr);
	if(B.has(io::U) && has(io::U) && to<N_sph()) uin (to) = B.uin (fr);
	if(B.has(io::Y) && has(io::Y) && to<N_sph()) uprd(to) = B.uprd(fr);
	if(B.has(io::I) && has(io::I) && to<N_sph()) udin(to) = B.udin(fr);
	if(B.has(io::E) && has(io::E) && to<N_sph()) udex(to) = B.udex(fr);
	if(B.has(io::S) && has(io::S) && to<N_sph()) ent (to) = B.ent (fr);
	if(B.has(io::R) && has(io::R) && to<N_sph()) srho(to) = B.srho(fr);
	if(B.has(io::D) && has(io::D) && to<N_sph()) drho(to) = B.drho(fr);
	if(B.has(io::V) && has(io::V) && to<N_sph()) vprd(to) = B.vprd(fr);
	if(B.has(io::Q) && has(io::Q) && to<N_sph()) sigq(to) = B.sigq(fr);
	if(B.has(io::T) && has(io::T) && to<N_sph()) temp(to) = B.temp(fr);
	if(B.has(io::J) && has(io::J) && to<N_sph()) hdot(to) = B.hdot(fr);
#endif
      }
    }
    //--------------------------------------------------------------------------
    // I/O  both nemo and yanc format                                           
    //--------------------------------------------------------------------------
#ifdef falcON_NEMO
    //--------------------------------------------------------------------------
    // INPUT:                                                                   
    // The general way of action is as follows:                                 
    // 1. if some of the data wanted are not supported or if N differs from old 
    //    value, de-allocate and re-allocate the corresponding BodyData, whereby
    //    retaining old support.                                                
    // 2. trying to read the data wanted.                                       
    // That implies that old data are preserved by an input only if             
    // -  there is no change in N   and                                         
    // -  there is no change in the partial databits due to input               
    // -  they are not overridden by freshly read data.                         
    // For instance, if before input mxv are supported and N_new = N_old, but   
    // a is wanted to be read, then the old mxv will not be deleted, for they   
    // are in BodySrce, while a is in BodySink, which only must be reset()ed.   
    //--------------------------------------------------------------------------
    bool read_nemo_snapshot (                          // R: was time in range? 
			     nemo_in const&,           // I: nemo input         
			     io           &,           // O: what has been read 
			     real*        =0,          //[O: time]              
			     const io     =io::mxv,    //[I: what to read]      
			     char*        =0,          //[I: time range]        
			     const bool   =true);      //[I: warn: missing data]
    //--------------------------------------------------------------------------
    bool read_nemo_particles(                          // R: was time in range? 
			     nemo_in const&,           // I: nemo input         
			     io           &,           // O: what has been read 
			     real*        =0,          //[O: time]              
			     const io     =io::mxv,    //[I: what to read]      
			     char*        =0,          //[I: time range]        
			     const bool   =true);      //[I: warn: missing data]
    //--------------------------------------------------------------------------
    bool read_nemo_snapshots(                          // R: was time in range? 
			     const nemo_in*,           // I: nemo inputs        
			     int   const  &,           // I: # nemo inputs      
			     io           *,           // O: what has been read 
			     uint         *,           // O: how many have ---  
			     real*        =0,          //[O: time]              
			     const io*    =0,          //[I: what to read]      
			     char*        =0,          //[I: time range]        
			     const bool   =true);      //[I: warn: missing data]
    //--------------------------------------------------------------------------
    bool read_nemo_particles(                          // R: was time in range? 
			     const nemo_in*,           // I: nemo inputs        
			     int   const  &,           // I: # nemo inputs      
			     io           *,           // O: what has been read 
			     uint         *,           // O: how many have ---  
			     real*        =0,          //[O: time]              
			     const io*    =0,          //[I: what to read]      
			     char*        =0,          //[I: time range]        
			     const bool   =true);      //[I: warn: missing data]
    //--------------------------------------------------------------------------
    void write_nemo_snapshot(                          // write snapshot        
			     nemo_out const&,          // I: nemo output        
			     const real    * =0,       //[I: write time]        
			     io    const   & =io::mxv, //[I: what to write]     
			     uint  const   & =0,       //[I: only write K]      
			     uint  const   & =0)const; //[I: begin with this]   
    //--------------------------------------------------------------------------
    void write_nemo_particles(                         // write bodies to output
			      nemo_out const&,         // I: nemo output        
			      const real    * =0,      //[I: write time]        
			      io    const   & =io::mxv,//[I: what to write]     
			      uint  const   & =0,      //[I: only write K]      
			      uint  const   & =0)const;//[I: begin with this]   
#endif
    //--------------------------------------------------------------------------
    void write_yanc_ascii(                             // write ascii to ostream
			  std::ostream&,               // I: output stream      
			  const io=io::mxv) const;     //[I: what to write out ]
    //--------------------------------------------------------------------------
    io   read_yanc_ascii(                              // R: what has been read 
			 std::istream&,                // I: input stream       
			 const io=io::mxv);            //[I: what to read in]   
    //--------------------------------------------------------------------------
    void write_yanc_binary(                            // write bin to ostream  
			   std::ostream&,              // I: output stream      
			   const io =io::mxv) const;   //[I: what to write]     
    //--------------------------------------------------------------------------
    void read_yanc_binary (                            // read bin from istream 
			   std::istream&,              // I: input stream       
			   const io =io::mxv);         //[I: what to read]      
    //--------------------------------------------------------------------------
    void read_simple_ascii(                            // read simple ascii file
			   std::istream  &,            // I: input stream       
			   const io*const&,            // I: array: data items  
			   uint     const&,            // I: # total lines      
			   uint     const& = 0);       // I: # lines with SPH   
    //--------------------------------------------------------------------------
    // iterator                                                                 
    //--------------------------------------------------------------------------
  public:
    class iterator {
      //........................................................................
      // basic stuff (construction, iteration, etc)                             
      //........................................................................
      friend class sbodies;
    private:
      iterator();                                  // not implemented           
      //........................................................................
      // data members                                                           
      //........................................................................
      const sbodies *B;                            // pointer to OBJECTS        
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
      iterator(const sbodies*b, const uint k) : B(b), K(k) {}
    public:
      iterator(const iterator&I)              : B(I.B), K(I.K) {}
      //........................................................................
      // forward iteration                                                      
      //........................................................................
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
      const sbodies*const&mybodies() const { return B; }
      //........................................................................
      // data access                                                            
      //........................................................................
#define DATA_ACCESS_REF(TYPE,NAME)					\
             TYPE      &NAME()                 { return B->NAME(K); }	\
             TYPE const&NAME()           const { return B->NAME(K); }	\
      friend TYPE const&NAME(iterator const&I) { return I.NAME(); }
      //........................................................................
  public:
      iterator next   () const { return iterator(B,K+1); }
      uint         const&index() const       { return K; }
      friend uint  const&index(const iterator&I) { return I.index(); }
      DATA_ACCESS_REF(real,mass)
      DATA_ACCESS_REF(vect,pos)
      real&pos  (const int i) { return pos()[i]; }
      DATA_ACCESS_REF(vect,vel)
      real&vel  (const int i) { return vel()[i]; }
      DATA_ACCESS_REF(vect,acc)
      real&acc  (const int i) { return acc()[i]; }
      DATA_ACCESS_REF(real,eps)
      DATA_ACCESS_REF(flag,flg)
      DATA_ACCESS_REF(int, key)
      DATA_ACCESS_REF(real,pot)
      DATA_ACCESS_REF(real,pex)
      DATA_ACCESS_REF(real,rho)
      DATA_ACCESS_REF(real,aux)
      DATA_ACCESS_REF(indx,level)
      DATA_ACCESS_REF(uint,num)
      DATA_ACCESS_REF(real,size)
#ifdef falcON_SPH
      DATA_ACCESS_REF(uint,snum)
      DATA_ACCESS_REF(real,uin)
      DATA_ACCESS_REF(real,uprd)
      DATA_ACCESS_REF(real,udin)
      DATA_ACCESS_REF(real,udex)
      DATA_ACCESS_REF(real,ent)
      DATA_ACCESS_REF(real,srho)
      DATA_ACCESS_REF(real,drho)
      DATA_ACCESS_REF(vect,vprd)
      DATA_ACCESS_REF(real,sigq)
      DATA_ACCESS_REF(real,temp)
      DATA_ACCESS_REF(real,hdot)
#endif
#undef DATA_ACCESS_REF
      //........................................................................
      // other const methods                                                    
      //........................................................................
      amom        angmom()           const { return pos() ^ vel(); }
      friend amom angmom(iterator const&I) { return I.angmom(); }
      //........................................................................
      // flag manipulations                                                     
      //........................................................................
      bool flag_is_set       (int const&F) const { return flg().is_set(F); }
      operator const flag&   ()            const { return flg(); }
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
      bool is_source     () const { return mass() != zero; }
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
	return o<<I.index();
      }
    };
    //--------------------------------------------------------------------------
    // types                                                                    
    //--------------------------------------------------------------------------
    typedef class iterator body;                   // type of body              
    typedef vect     const&const_vect_type;
    typedef vect          &vect_type;
    //--------------------------------------------------------------------------
    // initializing iterators                                                   
    //--------------------------------------------------------------------------
    iterator body_no(uint const&i) const { return iterator(this,i); }
    iterator bodyNo (uint const&i) const { return iterator(this,i); }
    iterator begin_bodies()        const { return iterator(this,0); }
    iterator end_bodies  ()        const { return iterator(this,N_bodies()); }
    iterator last_bodies ()        const { return iterator(this,N_bodies()-1); }
    iterator back_bodies ()        const { return iterator(this,N_bodies()-1); }
    iterator begin_SPH   ()        const { return iterator(this,0); }
    iterator end_SPH     ()        const { return iterator(this,N_sph()); }
    iterator back_SPH    ()        const { return iterator(this,N_sph()-1); }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::nbody_io<>                                                   //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<int, typename bodies_type=sbodies> struct nbody_io;
  //----------------------------------------------------------------------------
#define DEF_NBODY_IO(I,TYPE,NAME)					\
  template<typename bodies_type> struct nbody_io< 1<<I, bodies_type > {	\
    static const bool is_sph= (I >=IO_NOSPH);				\
    typedef TYPE                       io_type;				\
    typedef typename bodies_type::body body;				\
    static io_type const& const_access(body const&B) {			\
      return NAME(B); }							\
    static io_type      &       access(body      &B) {			\
      return B.NAME(); }						\
    static io_type*const& array_access(const bodies_type*const&B) {	\
      return B->NAME##_s(); }						\
    static io_type      &       access(const bodies_type*const&B,	\
				       uint const&i) {			\
      return B->NAME(i); }						\
    static uint    const& number      (const bodies_type*const&B) {	\
      return is_sph? B->N_sph() : B->N_bodies(); }			\
  };
  //----------------------------------------------------------------------------
  DEF_NBODY_IO(0,real,mass)
  DEF_NBODY_IO(1,vect,pos)
  DEF_NBODY_IO(2,vect,vel)
  DEF_NBODY_IO(3,real,eps)
  DEF_NBODY_IO(4,flag,flg)
  DEF_NBODY_IO(5,int,key)
  DEF_NBODY_IO(6,real,pot)
  DEF_NBODY_IO(7,real,pex)
  DEF_NBODY_IO(8,vect,acc)
  DEF_NBODY_IO(9,real,rho)
  DEF_NBODY_IO(10,real,aux)
  DEF_NBODY_IO(11,indx,level)
  DEF_NBODY_IO(12,uint,num)
  DEF_NBODY_IO(13,real,size)
#ifdef falcON_SPH
  DEF_NBODY_IO(14,uint,snum)
  DEF_NBODY_IO(15,real,uin)
  DEF_NBODY_IO(16,real,uprd)
  DEF_NBODY_IO(17,real,udin)
  DEF_NBODY_IO(18,real,udex)
  DEF_NBODY_IO(19,real,ent)
  DEF_NBODY_IO(20,real,srho)
  DEF_NBODY_IO(21,real,drho)
  DEF_NBODY_IO(22,vect,vprd)
  DEF_NBODY_IO(23,real,sigq)
  DEF_NBODY_IO(24,real,temp)
  DEF_NBODY_IO(25,real,hdot)
#endif
#undef DEF_NBODY_IO
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::LoopIO<>                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<template<int IO> class K, int BIT=0, int END=IO_NQUANT>
  struct LoopIO {
    static void loop(sbodies::body &Bi) {
      K< 1<<BIT >::act_on_body(Bi);
      LoopIO<K,BIT+1>::loop(Bi);
    }
    static void loop(const sbodies*const&BB) {
      K< 1<<BIT >::act_on_bodies(BB);
      LoopIO<K,BIT+1>::loop(BB);
    }
  };
  template<template<int IO> class K, int BIT> struct LoopIO<K,BIT,BIT> {
    static void loop(sbodies::body const&Bi) {}
    static void loop(const sbodies*const&BB) {}
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // struct nbdy::abodies                                                       
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // - holds externally allocated arrays of type areal with the body data       
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class abodies {
    uint         N, NSPH;
    int          BITS;
    const int   *FLG;
    areal       *POS[Ndim], *VEL[Ndim], *MAS;
#ifdef falcON_INDI
    areal       *EPS;
#endif
    uint        *NUM;
    areal       *ACC[Ndim], *POT, *RHO, *SIZ;
#ifdef falcON_SPH
    uint        *NSP;
    areal       *UIN, *UPR, *UDI, *UDE, *ENT, *SRH, *DRH, *VPR[Ndim],
                *SIQ, *TEM, *DHT;
#endif
    mutable bool SRC_DATA_CHANGED, SPH_DATA_CHANGED;
  public:
    //--------------------------------------------------------------------------
    typedef const_pseudo_tupel<Ndim,areal> cvect;
    typedef pseudo_tupel<Ndim,areal>       pvect;
    typedef const_pseudo_tupel<Ndim,areal> const_vect_type;
    typedef pseudo_tupel<Ndim,areal>       vect_type;
    //--------------------------------------------------------------------------
    abodies() : N(0), BITS(0), FLG(0), MAS(0), SIZ(0), NUM(0), POT(0), RHO(0),
#ifdef falcON_SPH
		NSP(0), UIN(0), UPR(0), UDI(0), UDE(0), ENT(0), SRH(0), DRH(0),
		SIQ(0), TEM(0), DHT(0),
#endif
		SRC_DATA_CHANGED(1), SPH_DATA_CHANGED(1)
    {
      for(register int d=0; d!=Ndim; ++d) {
	POS[d] = 0;
	VEL[d] = 0;
	ACC[d] = 0;
#ifdef falcON_SPH
	VPR[d] = 0;
#endif
      }
    }
    //--------------------------------------------------------------------------
    void mark_sph_data_read    () const { SPH_DATA_CHANGED = 0; }
    void mark_sph_data_changed () const { SPH_DATA_CHANGED = 1; }
    bool sph_data_changed      () const { return SPH_DATA_CHANGED; }
    //--------------------------------------------------------------------------
    void mark_srce_data_read   () const { SRC_DATA_CHANGED = 0; }
    void mark_srce_data_changed() const { SRC_DATA_CHANGED = 1; }
    bool srce_data_changed     () const { return SRC_DATA_CHANGED; }
    //--------------------------------------------------------------------------
    void set_N(uint const&n) { N=n; }
    //--------------------------------------------------------------------------
    void set_Nsph(uint const&n) { NSPH=n; }
    //--------------------------------------------------------------------------
    void set_flag(const int*const&x)
    { FLG=x; if(x) BITS |= io::f; else BITS &= ~io::f; }
    //--------------------------------------------------------------------------
    void set_num(uint*const&x)
    { NUM=x; if(x) BITS |= io::n; else BITS &= ~io::n; }
    //--------------------------------------------------------------------------
    void set_mass(areal*const&x)
    { MAS=x; if(x) BITS |= io::m; else BITS &= ~io::m; }
    //--------------------------------------------------------------------------
    void set_size(areal*const&x)
    { SIZ=x; if(x) BITS |= io::H; else BITS &= ~io::H; }
    //--------------------------------------------------------------------------
    void set_pot(areal*const&x)
    { POT=x; if(x) BITS |= io::p; else BITS &= ~io::p; }
    //--------------------------------------------------------------------------
    void set_rho(areal*const&x)
    { RHO=x; if(x) BITS |= io::r; else BITS &= ~io::r; }
#ifdef falcON_INDI
    //--------------------------------------------------------------------------
    void set_eps(areal*const&x)
    { EPS=x; if(x) BITS |= io::e; else BITS &= ~io::e; }
#endif
    //--------------------------------------------------------------------------
#if falcON_NDIM == 3
    void set_pos(areal*const&x, areal*const&y, areal*const&z)
    { 
      POS[0] = x; POS[1] = y; POS[2] = z;
      if(x) BITS |= io::x; else BITS &= ~io::x;
    }
    //--------------------------------------------------------------------------
    void set_vel(areal*const&x, areal*const&y, areal*const&z)
    { 
      VEL[0] = x; VEL[1] = y; VEL[2] = z;
      if(x) BITS |= io::v; else BITS &= ~io::v;
    }
    //--------------------------------------------------------------------------
    void set_acc(areal*const&x, areal*const&y, areal*const&z)
    { 
      ACC[0] = x; ACC[1] = y; ACC[2] = z;
      if(x) BITS |= io::a; else BITS &= ~io::a;
    }
#else
    //--------------------------------------------------------------------------
    void set_pos(const areal*const&x, const areal*const&y)
    { 
      POS[0] = x; POS[1] = y;
      if(x) BITS |= io::x; else BITS &= ~io::x;
    }
    //--------------------------------------------------------------------------
    void set_vel(const areal*const&x, const areal*const&y)
    { 
      VEL[0] = x; VEL[1] = y;
      if(x) BITS |= io::v; else BITS &= ~io::v;
    }
    //--------------------------------------------------------------------------
    void set_acc(areal*const&x, areal*const&y)
    { 
      ACC[0] = x; ACC[1] = y;
      if(x) BITS |= io::a; else BITS &= ~io::a;
    }
#endif
#ifdef falcON_SPH
    //--------------------------------------------------------------------------
    void set_nsp(uint*const&x)
    { NSP=x; if(x) BITS |= io::N; else BITS &= ~io::N; }
    //--------------------------------------------------------------------------
    void set_uin(areal*const&x)
    { UIN=x; if(x) BITS |= io::U; else BITS &= ~io::U; }
    //--------------------------------------------------------------------------
    void set_upr(areal*const&x)
    { UPR=x; if(x) BITS |= io::Y; else BITS &= ~io::Y; }
    //--------------------------------------------------------------------------
    void set_udin(areal*const&x)
    { UDI=x; if(x) BITS |= io::I; else BITS &= ~io::I; }
    //--------------------------------------------------------------------------
    void set_udex(areal*const&x)
    { UDE=x; if(x) BITS |= io::E; else BITS &= ~io::E; }
    //--------------------------------------------------------------------------
    void set_ent(areal*const&x)
    { ENT=x; if(x) BITS |= io::S; else BITS &= ~io::S; }
    //--------------------------------------------------------------------------
    void set_srh(areal*const&x)
    { SRH=x; if(x) BITS |= io::R; else BITS &= ~io::R; }
    //--------------------------------------------------------------------------
    void set_dvv(areal*const&x)
    { DRH=x; if(x) BITS |= io::D; else BITS &= ~io::D; }
    //--------------------------------------------------------------------------
    void set_siq(areal*const&x)
    { SIQ=x; if(x) BITS |= io::Q; else BITS &= ~io::Q; }
    //--------------------------------------------------------------------------
    void set_tem(areal*const&x)
    { TEM=x; if(x) BITS |= io::T; else BITS &= ~io::T; }
    //--------------------------------------------------------------------------
    void set_hdot(areal*const&x)
    { DHT=x; if(x) BITS |= io::J; else BITS &= ~io::J; }
    //--------------------------------------------------------------------------
# if falcON_NDIM == 3
    void set_vpr(areal*const&x, areal*const&y, areal*const&z)
    { 
      VPR[0] = x; VPR[1] = y; VPR[2] = z;
      if(x) BITS |= io::V; else BITS &= ~io::V;
    }
# else
    void set_vpr(areal*const&x, areal*const&y)
    { 
      VPR[0] = x; VPR[1] = y;
      if(x) BITS |= io::V; else BITS &= ~io::V;
    }
# endif
#endif
    //--------------------------------------------------------------------------
    bool has    (io const&b) const { return b == (b & BITS); }
    io   has_not(io const&b) const { return b & ~BITS; }
    //--------------------------------------------------------------------------
    const uint &N_bodies()   const { return N; }
    const uint &N_sph   ()   const { return NSPH; }
    const flag  flg   (int const&i) const { return flag(FLG[i]); }
          pvect pos   (int const&i) const { return pvect(POS,i); }
          pvect vel   (int const&i) const { return pvect(VEL,i); }
          areal&mass  (int const&i) const { return MAS[i]; }
#ifdef falcON_INDI
          areal&eps   (int const&i) const { return EPS[i]; }
#endif
          pvect acc   (int const&i) const { return pvect(ACC,i); }
          areal&pot   (int const&i) const { return POT[i]; }
          areal&rho   (int const&i) const { return RHO[i]; }
          uint &num   (int const&i) const { return NUM[i]; }
          areal&size  (int const&i) const { return SIZ[i]; }
#ifdef falcON_SPH
          areal&srho  (int const&i) const { return SRH[i]; }
          uint &snum  (int const&i) const { return NSP[i]; }
          areal&uin   (int const&i) const { return UIN[i]; }
          areal&uprd  (int const&i) const { return UPR[i]; }
          areal&udin  (int const&i) const { return UDI[i]; }
          areal&udex  (int const&i) const { return UDE[i]; }
          areal&ent   (int const&i) const { return ENT[i]; }
          areal&drho  (int const&i) const { return DRH[i]; }
          pvect vprd  (int const&i) const { return pvect(VPR,i); }
          areal&sigq  (int const&i) const { return SIQ[i]; }
          areal&temp  (int const&i) const { return TEM[i]; }
          areal&hdot  (int const&i) const { return DHT[i]; }
#endif
  };
}                                                  // namespace nbdy            
////////////////////////////////////////////////////////////////////////////////
// useful macros                                                                
////////////////////////////////////////////////////////////////////////////////
#ifndef LoopBodies
#define LoopBodies(BODIES_TYPE,                    /* type of bodies         */\
		   BODIES_PTER,                    /* pointer to bodies      */\
		   NAME)                           /* name for body          */\
  for(register BODIES_TYPE::body                   /* type of body           */\
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
  for(register BODIES_TYPE::body                   /* type of body           */\
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
  for(register BODIES_TYPE::body                   /* type of body           */\
      NAME2  = NAME1 + 1;                          /* from first body        */\
      NAME2 !=(BODIES_PTER)->end_bodies();         /* until beyond last body */\
    ++NAME2)                                       /* get next body          */
#endif
//------------------------------------------------------------------------------
#ifndef LoopSPHBodies
#define LoopSPHBodies(BODIES_TYPE,                 /* type of bodies         */\
		      BODIES_PTER,                 /* pointer to bodies      */\
		      NAME)                        /* name for body          */\
  for(register BODIES_TYPE::body                   /* type of body           */\
      NAME  =(BODIES_PTER)->begin_SPH();           /* from first SPH body    */\
      NAME !=(BODIES_PTER)->end_SPH();             /* until last SPH body    */\
    ++NAME)                                        /* get next body          */
#endif
//------------------------------------------------------------------------------
#ifndef LoopNonSPHBodies
#define LoopNonSPHBodies(BODIES_TYPE,              /* type of bodies         */\
		         BODIES_PTER,              /* pointer to bodies      */\
		         NAME)                     /* name for body          */\
  for(register BODIES_TYPE::body                   /* type of body           */\
      NAME  =(BODIES_PTER)->end_SPH();             /* from first non-SPH body*/\
      NAME !=(BODIES_PTER)->end_bodies();          /* until beyond last body */\
    ++NAME)                                        /* get next body          */
#endif
//------------------------------------------------------------------------------
#endif                                             // falcON_included_body_h    
