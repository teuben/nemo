// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// body.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2003                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// defines                                                                     |
//                                                                             |
// sbodies                                                                     |
// sbodies::iterator -> body                                                   |
// typedef sbodies bodies;                                                     |
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
#ifdef falcON_NEMO
#ifndef falcON_included_nmio_h
namespace nbdy { class nemo_in;                    // forward declaration       
                 class nemo_out; }                 // forward declaration       
#endif
#endif
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class nbdy::BodyComm                                                       
  //                                                                            
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
  //                                                                            
  // class nbdy::BodyDataBase                                                   
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  template<typename BodyData>
  class BodyDataBase {
    BodyDataBase           (BodyDataBase const&);  // not implemented           
    BodyDataBase& operator=(BodyDataBase const&);  // not implemented           
    //--------------------------------------------------------------------------
    // data members                                                             
    //--------------------------------------------------------------------------
  protected:
    static int data_bits(io const&b) {
      return BodyData::MINBITS() | (b & BodyData::MAXBITS());
    }
    //--------------------------------------------------------------------------
    uint         NB;                               // # bodies                  
    io           DATABITS;                         // which data do we hold?    
    uint         BBYTES;                           // bytes per body            
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
      DATABITS = data_bits(b);
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
      if(n != NB || data_bits(b) != DATABITS) {
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
    const bool is_supported(const io&b) const
    {
      return DATABITS.contains(BodyData::MAXBITS() & b);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
#define SETPTER(TYPE,NAME,BIT)						       \
    if(DATABITS & BIT) { NAME = static_cast<TYPE*>(static_cast<void*>(M));     \
                         M   += sizeof(TYPE) * NB; }			       \
    else 		 NAME = 0;
  //////////////////////////////////////////////////////////////////////////////
#define ACCESS(NAME,TYPE,ARRAY)						\
    TYPE &NAME(const uint&i) const { return ARRAY[i]; }
#define SUBACCESS(NAME,TYPE,ARRAY)					\
    TYPE &NAME(const uint&i, const int&d) const { return ARRAY[i][d]; }
#ifdef DEBUG
#  define ACCESS_CHECK(NAME,TYPE,ARRAY)					\
    TYPE &NAME(const uint&i) const { 					\
      if(ARRAY) return ARRAY[i];					\
      else falcON_Error("data wanted not allocated");			\
    }
#  define SUBACCESS_CHECK(NAME,TYPE,ARRAY)				\
    TYPE &NAME(const uint&i, const int&d) const { 			\
      if(ARRAY) return ARRAY[i][d]; }					\
      else falcON_Error("data wanted not allocated");			\
    }
#else
#  define ACCESS_CHECK(NAME,TYPE,ARRAY)    ACCESS(NAME,TYPE,ARRAY)
#  define SUBACCESS_CHECK(NAME,TYPE,ARRAY) SUBACCESS(NAME,TYPE,ARRAY)
#endif
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class nbdy::BodySrce                                                       
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // - holds those data that every gravity source particle needs                
  // - these are mass, position, velocity [, flag, eps, key]                    
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class BodySrce :
    public BodyComm,
    public BodyDataBase<BodySrce> {
    //--------------------------------------------------------------------------
    BodySrce           (const BodySrce&);          // not implemented           
    BodySrce& operator=(const BodySrce&);          // not implemented           
    //--------------------------------------------------------------------------
  public:
    static const io MINBITS () { return io::mxv; }
    static const io MAXBITS () { return io::mxv|io::f|io::e|io::k; }
    static const io DEFBITS () { return io::mxv|io::f; }
    static const io NEMOBITS() { return MAXBITS(); }
  private:
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
      BodyComm(nb,ns,b),
      BodyDataBase<BodySrce>(nb,b),
      CUSE ( cu ),
      NA   ( na > 0? na : NB )
    {
      set_pointers();                              //   set pointers to data    
    }
    //--------------------------------------------------------------------------
    void reset(const uint&nb,                      // I: # bodies               
	       const io  &b,                       // I: data bits              
	       const uint ns = 0) {                //[I: # SPH particles]       
      BodyComm::reset(nb,ns,b);
      if(BodyDataBase<BodySrce>::reset(nb,b)) {    // IF(memory reset)         >
	set_pointers();                            //   set pointers to data    
	reset_flags ();                            //   reset flags             
      }                                            // <                         
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
    public BodyDataBase<BodySink> {
    BodySink           (const BodySink&);          // not implemented           
    BodySink& operator=(const BodySink&);          // not implemented           
    //--------------------------------------------------------------------------
    // data members                                                             
    //--------------------------------------------------------------------------
  public:
    static const io MINBITS () { return io::a|io::p; }
    static const io MAXBITS () { return io::a|io::p|io::P|
	        		        io::r|io::y|io::n|io::l; }
    static const io DEFBITS () { return MINBITS(); }
    static const io NEMOBITS() { return MINBITS() |io::r|io::y|io::l; }
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
	SETPTER(real,PEX,io::P)
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
      BodyDataBase<BodySink>(n,b) 
    {
      set_pointers();                              // set pointers to data      
    }
    //--------------------------------------------------------------------------
    void reset  (const uint n,                     // I: N_bodies               
		 const io   b)                     // I: data bits              
    {
      if(BodyDataBase<BodySink>::reset(n,b))       // IF(memory reset)         >
	set_pointers();                            //   set pointers to data   <
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
    public BodyDataBase<BodyPsph> {
    BodyPsph           (const BodyPsph&);          // not implemented           
    BodyPsph& operator=(const BodyPsph&);          // not implemented           
    //--------------------------------------------------------------------------
    // data members                                                             
    //--------------------------------------------------------------------------
  public:
    static const io MINBITS () { return io::sphmin; }
    static const io MAXBITS () { return io::sphmax; }
    static const io DEFBITS () { return io::sphmin; }
    static const io NEMOBITS() { return io::o; }
  protected:
    real        *SIZ;                              // pointer to SPH sizes      
#ifdef falcON_SPH
    real        *TEM;                              // pointer to SPH Temperature
    real        *EIN;                              // pointer to SPH in. energe 
    real        *ENT;                              // pointer to SPH entropies  
    real        *SRH;                              // pointer to SPH densities  
    real        *DVV;                              // pointer to SPH div(v)     
    vect        *RTV;                              // pointer to SPH rot(v)     
#endif
  private:
    //--------------------------------------------------------------------------
    void set_pointers() {
      if(ALLOC) {
	register char* M = ALLOC;
	SETPTER(real,SIZ,io::s)
#ifdef falcON_SPH
	SETPTER(real,TEM,io::T)
	SETPTER(real,EIN,io::u)
	SETPTER(real,ENT,io::S)
	SETPTER(real,SRH,io::R)
	SETPTER(real,DVV,io::D)
	SETPTER(vect,RTV,io::V)
#endif
      } else {
	SIZ   = 0;
#ifdef falcON_SPH
	TEM   = 0;
	EIN   = 0;
	ENT   = 0;
	SRH   = 0;
	DVV   = 0;
	RTV   = 0;
#endif
      }
    }
    //--------------------------------------------------------------------------
  protected:
    BodyPsph(const uint n = 0,                     //[I: N_sph]                 
	     const io   b = DEFBITS()) :           //[I: data bits]             
      BodyDataBase<BodyPsph>(n,b) 
    {
      set_pointers();                              // set pointers to data      
    }
    //--------------------------------------------------------------------------
    void reset  (const uint n,                     // I: N_sph                  
		 const io   b)                     // I: data bits              
    {
      if(BodyDataBase<BodyPsph>::reset(n,b))       // IF(memory reset)         >
	set_pointers();                            //   set pointers to data   <
    }
    //--------------------------------------------------------------------------
    // const member methods                                                     
    //--------------------------------------------------------------------------
  public:
    ACCESS      (size,real,SIZ)
#ifdef falcON_SPH
    ACCESS      (srho,real,SRH)
    ACCESS_CHECK(temp,real,TEM)
    ACCESS_CHECK(ein, real,EIN)
    ACCESS_CHECK(ent, real,ENT)
    ACCESS_CHECK(divv,real,DVV)
    ACCESS_CHECK(rotv,vect,RTV)
    bool        has_temp  ()             const { return TEM!=0; }
    bool        has_ent   ()             const { return ENT!=0; }
    bool        has_divv  ()             const { return DVV!=0; }
    bool        has_rotv  ()             const { return RTV!=0; }
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
    // data member: pointer to memory                                           
    //--------------------------------------------------------------------------
    char * MEMORY[IO_NQUANT];
    //--------------------------------------------------------------------------
    // static constants                                                         
    //--------------------------------------------------------------------------
  public:
    static const io MINBITS () {
      return BodySrce::MINBITS() | BodySink::MINBITS() | BodyPsph::MINBITS(); }
    static const io MAXBITS () {
      return BodySrce::MAXBITS() | BodySink::MAXBITS() | BodyPsph::MAXBITS(); }
    static const io DEFBITS () {
      return BodySrce::DEFBITS() | BodySink::DEFBITS() | BodyPsph::DEFBITS(); }
    static const io NEMOBITS() {
      return BodySrce::NEMOBITS() |BodySink::NEMOBITS() |BodyPsph::NEMOBITS(); }
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
    // in the above example.                                                    
    //--------------------------------------------------------------------------
  public:
    sbodies(const uint nb = 0,                     //[I: # bodies]              
	    const io   b  = DEFBITS(),             //[I: data bits]             
	    const uint ns = 0,                     //[I: # SPH particles]       
	    const uint na = 0,                     //[I: # bodies in tree]      
	    const bool cu = false) :               //[I: CUSE flag]             
      BodySrce(nb,b,ns,na,cu),
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
	if(B.has(io::P) && has(io::P))
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
	if(B.has(io::s) && has(io::s) && to+n<N_sph())
	  for(register int i=0; i!=n; ++i)  size(to+i) = B.size(fr);
#ifdef falcON_SPH
	if(B.has(io::T) && has(io::T) && to+n<N_sph()) 
	  for(register int i=0; i!=n; ++i)  temp(to+i) = B.temp(fr);
	if(B.has(io::u) && has(io::u) && to+n<N_sph()) 
	  for(register int i=0; i!=n; ++i)  ein (to+i) = B.ein (fr);
	if(B.has(io::S) && has(io::S) && to+n<N_sph()) 
	  for(register int i=0; i!=n; ++i)  ent (to+i) = B.ent (fr);
	if(B.has(io::R) && has(io::R) && to+n<N_sph()) 
	  for(register int i=0; i!=n; ++i)  srho(to+i) = B.srho(fr);
	if(B.has(io::D) && has(io::D) && to+n<N_sph()) 
	  for(register int i=0; i!=n; ++i)  divv(to+i) = B.divv(fr);
	if(B.has(io::V) && has(io::V) && to+n<N_sph()) 
	  for(register int i=0; i!=n; ++i)  rotv(to+i) = B.rotv(fr);
#endif
      } else {
	if(B.has(io::x) && has(io::x)) pos (to) = B.pos (fr);
	if(B.has(io::m) && has(io::m)) mass(to) = B.mass(fr);
	if(B.has(io::v) && has(io::v)) vel (to) = B.vel (fr);
	if(B.has(io::e) && has(io::e)) eps (to) = B.eps (fr);
	if(B.has(io::k) && has(io::k)) key (to) = B.key (fr);
	if(B.has(io::p) && has(io::p)) pot (to) = B.pot (fr);
	if(B.has(io::P) && has(io::P)) pex (to) = B.pex (fr);
	if(B.has(io::a) && has(io::a)) acc (to) = B.acc (fr);
	if(B.has(io::f) && has(io::f)) flg (to) = B.flg (fr);
	if(B.has(io::l) && has(io::l)) level(to)= B.level(fr);
	if(B.has(io::r) && has(io::r)) rho (to) = B.rho (fr);
	if(B.has(io::y) && has(io::y)) aux (to) = B.aux (fr);
	if(B.has(io::n) && has(io::n)) num (to) = B.num (fr);
	if(B.has(io::s) && has(io::s) && to<N_sph()) size(to) = B.size(fr);
#ifdef falcON_SPH
	if(B.has(io::T) && has(io::T) && to<N_sph()) temp(to) = B.temp(fr);
	if(B.has(io::u) && has(io::u) && to<N_sph()) ein (to) = B.ein (fr);
	if(B.has(io::S) && has(io::S) && to<N_sph()) ent (to) = B.ent (fr);
	if(B.has(io::R) && has(io::R) && to<N_sph()) srho(to) = B.srho(fr);
	if(B.has(io::D) && has(io::D) && to<N_sph()) divv(to) = B.divv(fr);
	if(B.has(io::V) && has(io::V) && to<N_sph()) rotv(to) = B.rotv(fr);
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
			     nemo_out     &,           // I: nemo output        
			     const real* =0,           //[I: write time]        
			     const io    =io::mxv,     //[I: what to write]     
			     const uint  =0,           //[I: only write K]      
			     const uint  =0) const;    //[I: begin with this]   
    //--------------------------------------------------------------------------
    void write_nemo_particles(                         // write bodies to output
			      nemo_out    &,           // I: nemo output        
			      const real* =0,          //[I: write time]        
			      const io    =io::mxv,    //[I: what to write]     
			      const uint  =0,          //[I: only write K]      
			      const uint  =0) const;   //[I: begin with this]   
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
			   std::istream&,              // I: input stream       
			   uint   const&,              // I: # lines to read    
			   const io    *);             // I: array: data items  
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
      DATA_ACCESS_REF(real,pot)
      DATA_ACCESS_REF(real,pex)
      DATA_ACCESS_REF(real,aux)
      DATA_ACCESS_REF(int, key)
      DATA_ACCESS_REF(real,rho)
      DATA_ACCESS_REF(uint,num)
      DATA_ACCESS_REF(indx,level)
      DATA_ACCESS_REF(flag,flg)
      DATA_ACCESS_REF(real,size)
#ifdef falcON_SPH
      DATA_ACCESS_REF(real,srho)
      DATA_ACCESS_REF(real,temp)
      DATA_ACCESS_REF(real,ein)
      DATA_ACCESS_REF(real,ent)
      DATA_ACCESS_REF(real,divv)
      DATA_ACCESS_REF(vect,rotv)
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
    typedef class iterator     body_type;          // type of body              
    //--------------------------------------------------------------------------
    // initializing iterators                                                   
    //--------------------------------------------------------------------------
    iterator body_no(uint const&i) const { return iterator(this,i); }
    iterator begin_bodies()        const { return iterator(this,0); }
    iterator end_bodies  ()        const { return iterator(this,N_bodies()); }
    iterator last_bodies ()        const { return iterator(this,N_bodies()-1); }
    iterator back_bodies ()        const { return iterator(this,N_bodies()-1); }
    iterator end_SPH     ()        const { return iterator(this,N_sph()); }
    iterator back_SPH    ()        const { return iterator(this,N_sph()-1); }
  };
  //////////////////////////////////////////////////////////////////////////////
  // typedefs                                                                   
  //////////////////////////////////////////////////////////////////////////////
  typedef sbodies            bodies;               // define nbdy::bodies       
  typedef sbodies::body_type body;                 // define nbdy::body         
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // struct nbdy::abodies                                                       
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // - holds externally allocated arrays of type areal with the body data       
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class barrays {
    uint         N;
    int          BITS;
    const int   *FLG;
    const areal *POS[Ndim], *VEL[Ndim], *MAS;
#ifdef falcON_INDI
    areal       *EPS;
#endif
    uint        *NUM;
    areal       *ACC[Ndim], *POT, *RHO, *SIZ;
  public:
    //--------------------------------------------------------------------------
    barrays() : N(0), BITS(0), FLG(0), MAS(0), SIZ(0), NUM(0), POT(0), RHO(0)
    {
      for(register int d=0; d!=Ndim; ++d) POS[d] = VEL[d] = ACC[d] = 0;
    }
    //--------------------------------------------------------------------------
    void set_N(uint const&n) { N=n; }
    //--------------------------------------------------------------------------
    void set_flag(const int*const&x)
    { FLG=x; if(x) BITS |= io::f; else BITS &= ~io::f; }
    //--------------------------------------------------------------------------
    void set_num(uint*const&x)
    { NUM=x; if(x) BITS |= io::n; else BITS &= ~io::n; }
    //--------------------------------------------------------------------------
    void set_mass(const areal*const&x)
    { MAS=x; if(x) BITS |= io::m; else BITS &= ~io::m; }
    //--------------------------------------------------------------------------
    void set_size(areal*const&x)
    { SIZ=x; if(x) BITS |= io::s; else BITS &= ~io::s; }
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
    void set_pos(const areal*const&x, const areal*const&y, const areal*const&z)
    { 
      POS[0] = x; POS[1] = y; POS[2] = z;
      if(x) BITS |= io::x; else BITS &= ~io::x;
    }
    //--------------------------------------------------------------------------
    void set_vel(const areal*const&x, const areal*const&y, const areal*const&z)
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
    //--------------------------------------------------------------------------
    bool has    (io const&b) const { return b == (b & BITS); }
    io   has_not(io const&b) const { return b & ~BITS; }
    //--------------------------------------------------------------------------
    const uint &N_bodies()   const { return N; }
    const flag  flg  (int const&i) const { return flag(FLG[i]); }
    const areal&pos_x(int const&i) const { return POS[0][i]; }
    const areal&pos_y(int const&i) const { return POS[1][i]; }
    const areal&pos_z(int const&i) const { return POS[2][i]; }
    const areal*pos  (int const&d) const { return POS[d]; }
    const areal&vel_x(int const&i) const { return VEL[0][i]; }
    const areal&vel_y(int const&i) const { return VEL[1][i]; }
    const areal&vel_z(int const&i) const { return VEL[2][i]; }
    const areal&size (int const&i) const { return SIZ[i]; }
    const areal&mass (int const&i) const { return MAS[i]; }
#ifdef falcON_INDI
          areal&eps  (int const&i) const { return EPS[i]; }
#endif
          areal&acc_x(int const&i) const { return ACC[0][i]; }
          areal&acc_y(int const&i) const { return ACC[1][i]; }
          areal&acc_z(int const&i) const { return ACC[2][i]; }
          areal&pot  (int const&i) const { return POT[i]; }
          areal&rho  (int const&i) const { return RHO[i]; }
          uint &num  (int const&i) const { return NUM[i]; }
  };
}                                                  // namespace nbdy            
////////////////////////////////////////////////////////////////////////////////
// useful macros                                                                
////////////////////////////////////////////////////////////////////////////////
#ifndef LoopBodies
#define LoopBodies(BODIES_TYPE,                    /* type of bodies         */\
		   BODIES_PTER,                    /* pointer to bodies      */\
		   NAME)                           /* name for body          */\
  for(register BODIES_TYPE::body_type              /* type of body           */\
      NAME  = BODIES_PTER->begin_bodies();         /* from first body        */\
      NAME != BODIES_PTER->end_bodies();           /* until beyond last body */\
    ++NAME)                                        /* get next body          */
#endif
//------------------------------------------------------------------------------
#ifndef LoopSBodies
#define LoopSBodies(SBODIES_PTER,                  /* pointer to sbodies     */\
		    NAME)                          /* name for body          */\
  for(register sbodies::body_type                  /* type of body           */\
      NAME  = SBODIES_PTER->begin_bodies();        /* from first body        */\
      NAME != SBODIES_PTER->end_bodies();          /* until beyond last body */\
    ++NAME)                                        /* get next body          */
#endif
//------------------------------------------------------------------------------
#ifndef LoopSPHBodies
#define LoopSPHBodies(BODIES_TYPE,                 /* type of bodies         */\
		      BODIES_PTER,                 /* pointer to bodies      */\
		      NAME)                        /* name for body          */\
  for(register BODIES_TYPE::body_type              /* type of body           */\
      NAME  = BODIES_PTER->begin_SPH();            /* from first body        */\
      NAME != BODIES_PTER->end_SPH();              /* until beyond last body */\
    ++NAME)                                        /* get next body          */
#endif
//------------------------------------------------------------------------------
#endif                                             // falcON_included_body_h    
