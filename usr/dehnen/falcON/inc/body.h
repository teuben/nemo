// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// body.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2002                                          |
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
#ifndef included_body_h
#define included_body_h 1

#ifndef included_auxx_h
#  include <public/auxx.h>
#endif
#ifndef included_flag_h
#  include <public/flag.h>
#endif
#ifndef included_nbio_h
#  include <public/nbio.h>
#endif
////////////////////////////////////////////////////////////////////////////////
#ifdef ALLOW_NEMO
namespace nbdy { class nemo_in;                    // forward declaration       
                 class nemo_out; }                 // forward declaration       
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
      NBODIES( nb ), NSPH( ns ), ALLBITS( b ) {}
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
    bool  all_bits_contain (io const&b) const { return b == b & ALLBITS; }
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
#define SETPTER(TYPE,NAME,BIT)							\
    if(DATABITS & BIT) { NAME = static_cast<TYPE*>(static_cast<void*>(M));	\
                         M   += sizeof(TYPE) * NB; }				\
    else 		 NAME = 0;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class nbdy::BodySrce                                                       
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
    static const int MINBITS () { return io::mxv; }
    static const int MAXBITS () { return io::mxv|io::f|io::e|io::k; }
    static const int DEFBITS () { return io::mxv|io::f; }
    static const int NEMOBITS() { return MAXBITS(); }
  private:
    mutable bool CUSE;                             // flag for tree usage change
    mutable uint NA;                               // # active bodies           
  protected:
    vect        *POS;                              // pointer to positions      
    flag        *FLG;                              // pointer to flags          
    real        *MAS;                              // pointer to masses         
    vect        *VEL;                              // pointer to velocities     
    real        *EPS;                              // pointer to eps_i          
    int         *KEY;                              // pointer to keyes          
  private:
    //--------------------------------------------------------------------------
    void set_pointers() {
      if(ALLOC) {
	register char* M = ALLOC;
	SETPTER(vect,POS,io::x)
	SETPTER(flag,FLG,io::f)
	SETPTER(real,MAS,io::m)
	SETPTER(vect,VEL,io::v)
	SETPTER(real,EPS,io::e)
	SETPTER(int, KEY,io::k)
      } else {
	POS  = 0;
	FLG  = 0;
	MAS  = 0;
	VEL  = 0;
	EPS  = 0;
	KEY  = 0;
      }
    }
    //--------------------------------------------------------------------------
  public:
    BodySrce(const uint nb,                        // I: # bodies               
	     const io   b ,                        // I: data bits              
	     const uint ns = 0,                    //[I: # SPH particles]       
	     const uint na = 0,                    //[I: # active bodies]       
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
    const uint &N_active  ()             const { return NA; }
    bool        has_flg   ()             const { return FLG != 0; }
    bool        has_eps   ()             const { return EPS!=0; }
    bool        has_key   ()             const { return KEY!=0; }
    flag       &flg       (const uint&i) const { return FLG[i]; }
    vect       &pos       (const uint&i) const { return POS[i]; }
    real       &pos       (const uint&i,
			   const  int&d) const { return POS[i][d]; }
    real       &mas       (const uint&i) const { return MAS[i]; }
    real       &mass      (const uint&i) const { return MAS[i]; }
    vect       &vel       (const uint&i) const { return VEL[i]; }
    real       &vel       (const uint&i,
			   const  int&d) const { return VEL[i][d]; }
    real       &eps       (const uint&i) const { return EPS[i]; }
    int        &key       (const uint&i) const { return KEY[i]; }
#ifdef ALLOW_NEMO
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
      }
    }
    //--------------------------------------------------------------------------
    void flag_all_as_sink() const {                // flag all bodies as sinks  
      for(register uint i=0; i!=NB; ++i) FLG[i].add(flag::SINK);
    }
    //--------------------------------------------------------------------------
    void allow_tree_usage  (const uint&i) const {
      if(! is_in_tree(FLG[i])) {                   // IF(body not in tree)     >
	FLG[i].un_set(flag::NOT_IN_TREE);          //   unflag 'not in tree'    
	NA++;                                      //   increment N_active      
	CUSE = true;                               //   set flag for change     
      }                                            // <                         
    }
    //--------------------------------------------------------------------------
    void forbid_tree_usage (const uint&i) const {
      if(is_in_tree(FLG[i])) {                     // IF(body was in tree)     >
	FLG[i].add(flag::NOT_IN_TREE);             //   flag 'not in tree'      
	NA--;                                      //   decrement N_active      
	CUSE = true;                               //   set flag for change     
      }                                            // <                         
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
  // - whether or not key and eps are contained is must be told                 
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
    static const int MINBITS () { return io::a|io::p; }
    static const int MAXBITS () { return io::a|io::p|io::r|
				    io::y|io::n|io::l; }
    static const int DEFBITS () { return MINBITS(); }
    static const int NEMOBITS() { return MINBITS() |io::r|io::y|io::l; }
  protected:
    real        *POT;                              // pointer to N-body pots    
    vect        *ACC;                              // pointer to acceleration   
    real        *RHO;                              // pointer to densities      
    real        *AUX;                              // pointer to aux data       
    uint        *NUM;                              // pointer to # neighbours   
    indx        *LEV;                              // pointer to levels         
  private:
    //--------------------------------------------------------------------------
    void set_pointers() {
      if(ALLOC) {
	register char* M = ALLOC;
	SETPTER(real,POT,io::p)
	SETPTER(vect,ACC,io::a)
	SETPTER(real,RHO,io::r)
	SETPTER(real,AUX,io::y)
	SETPTER(uint,NUM,io::n)
	SETPTER(indx,LEV,io::l)
      } else {
	POT  = 0;
	ACC  = 0;
	RHO  = 0;
	AUX  = 0;
	NUM  = 0;
	LEV  = 0;
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
#ifdef ALLOW_NEMO
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
    bool        has_rho   ()             const { return RHO!=0; }
    bool        has_aux   ()             const { return AUX!=0; }
    bool        has_num   ()             const { return NUM!=0; }
    bool        has_lev   ()             const { return LEV!=0; }
    bool        has_level ()             const { return LEV!=0; }
    real       &pot       (const uint&i) const { return POT[i]; }
    vect       &acc       (const uint&i) const { return ACC[i]; }
    real       &rho       (const uint&i) const { return RHO[i]; }
    real       &aux       (const uint&i) const { return AUX[i]; }
    uint       &num       (const uint&i) const { return NUM[i]; }
    indx       &lev       (const uint&i) const { return LEV[i]; }
    indx       &level     (const uint&i) const { return LEV[i]; }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class nbdy::BodyPsph                                                       
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // - holds those data that belong to SPH particles only                       
  // - these are size (more to be added later)                                  
  // - whether or not key and eps are contained is must be told                 
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
    static const int MINBITS () { return io::o; }
    static const int MAXBITS () { return io::s; }
    static const int DEFBITS () { return io::s; }
    static const int NEMOBITS() { return io::o; }
  protected:
    real        *SIZ;                              // pointer to SPH sizes      
  private:
    //--------------------------------------------------------------------------
    void set_pointers() {
      if(ALLOC) {
	register char* M = ALLOC;
	SETPTER(real,SIZ,io::s)
      } else {
	SIZ   = 0;
      }
    }
#undef SETPTER
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
    bool        has_siz   ()             const { return SIZ!=0; }
    bool        has_size  ()             const { return SIZ!=0; }
    real       &siz       (const uint&i) const { return SIZ[i]; }
    real       &size      (const uint&i) const { return SIZ[i]; }
  };
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
    // static constants                                                         
    //--------------------------------------------------------------------------
  public:
    static const int MINBITS () {
      return BodySrce::MINBITS() | BodySink::MINBITS() | BodyPsph::MINBITS(); }
    static const int MAXBITS () {
      return BodySrce::MAXBITS() | BodySink::MAXBITS() | BodyPsph::MAXBITS(); }
    static const int DEFBITS () {
      return BodySrce::DEFBITS() | BodySink::DEFBITS() | BodyPsph::DEFBITS(); }
    static const int NEMOBITS() {
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
    //   sbodies BB(io::mxv | io::y);                                           
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
	    const uint na = 0,                     //[I: # active bodies]       
	    const bool cu = false) :               //[I: CUSE flag]             
      BodySrce(nb,b,ns,na,cu),
      BodySink(nb,b),
      BodyPsph(ns,b) {}
    //--------------------------------------------------------------------------
    // enum nbdy::sbodies::part                                                 
    //--------------------------------------------------------------------------
    enum part { Srce=1, Sink=2, Psph=4 };
    //--------------------------------------------------------------------------
    // reset partial data bits, not N                                           
    //--------------------------------------------------------------------------
    void reset(const part p,                       // I: what to reset          
	       const io  &b = DEFBITS())           // I: data bits              
    { switch(p) {
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
    } }
    //--------------------------------------------------------------------------
    // reset N_bodies, data bits, [and N_sph]  deletes all data                 
    //--------------------------------------------------------------------------
    void reset(const uint&nb,                      // I: # bodies               
	       const io  &b,                       // I: data bits              
	       const uint ns = 0) {                //[I: # SPH particles]       
      BodySrce::reset(nb,b,ns);
      BodySink::reset(nb,b);
      BodyPsph::reset(ns,b);
    }
    //--------------------------------------------------------------------------
    // copy a body                                                              
    //--------------------------------------------------------------------------
    void copy_body(sbodies const& B,
		   int     const& fr,
		   int     const& to,
		   int     const& n = 0)
    {
      if(n) {
	if(B.all_bits() & io::x && all_bits() & io::x) 
	  for(register int i=0; i!=n; ++i)  pos(to+i) = B.pos(fr);
	if(B.all_bits() & io::m && all_bits() & io::m) 
	  for(register int i=0; i!=n; ++i)  mas(to+i) = B.mas(fr);
	if(B.all_bits() & io::v && all_bits() & io::v) 
	  for(register int i=0; i!=n; ++i)  vel(to+i) = B.vel(fr);
	if(B.all_bits() & io::e && all_bits() & io::e) 
	  for(register int i=0; i!=n; ++i)  eps(to+i) = B.eps(fr);
	if(B.all_bits() & io::k && all_bits() & io::k) 
	  for(register int i=0; i!=n; ++i)  key(to+i) = B.key(fr);
	if(B.all_bits() & io::p && all_bits() & io::p) 
	  for(register int i=0; i!=n; ++i)  pot(to+i) = B.pot(fr);
	if(B.all_bits() & io::a && all_bits() & io::a) 
	  for(register int i=0; i!=n; ++i)  acc(to+i) = B.acc(fr);
	if(B.all_bits() & io::f && all_bits() & io::f) 
	  for(register int i=0; i!=n; ++i)  flg(to+i) = B.flg(fr);
	if(B.all_bits() & io::l && all_bits() & io::l) 
	  for(register int i=0; i!=n; ++i)  lev(to+i) = B.lev(fr);
	if(B.all_bits() & io::r && all_bits() & io::r) 
	  for(register int i=0; i!=n; ++i)  rho(to+i) = B.rho(fr);
	if(B.all_bits() & io::y && all_bits() & io::y) 
	  for(register int i=0; i!=n; ++i)  aux(to+i) = B.aux(fr);
	if(B.all_bits() & io::n && all_bits() & io::n) 
	  for(register int i=0; i!=n; ++i)  num(to+i) = B.num(fr);
	if(B.all_bits() & io::s && all_bits() & io::s) 
	  for(register int i=0; i!=n; ++i)  siz(to+i) = B.siz(fr);
      } else {
	if(B.all_bits() & io::x && all_bits() & io::x) pos(to) = B.pos(fr);
	if(B.all_bits() & io::m && all_bits() & io::m) mas(to) = B.mas(fr);
	if(B.all_bits() & io::v && all_bits() & io::v) vel(to) = B.vel(fr);
	if(B.all_bits() & io::e && all_bits() & io::e) eps(to) = B.eps(fr);
	if(B.all_bits() & io::k && all_bits() & io::k) key(to) = B.key(fr);
	if(B.all_bits() & io::p && all_bits() & io::p) pot(to) = B.pot(fr);
	if(B.all_bits() & io::a && all_bits() & io::a) acc(to) = B.acc(fr);
	if(B.all_bits() & io::f && all_bits() & io::f) flg(to) = B.flg(fr);
	if(B.all_bits() & io::l && all_bits() & io::l) lev(to) = B.lev(fr);
	if(B.all_bits() & io::r && all_bits() & io::r) rho(to) = B.rho(fr);
	if(B.all_bits() & io::y && all_bits() & io::y) aux(to) = B.aux(fr);
	if(B.all_bits() & io::n && all_bits() & io::n) num(to) = B.num(fr);
	if(B.all_bits() & io::s && all_bits() & io::s) siz(to) = B.siz(fr);
      }
    }
    //--------------------------------------------------------------------------
    // I/O  both nemo and yanc format                                           
    //--------------------------------------------------------------------------
#ifdef ALLOW_NEMO
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
    // constant methods                                                         
    //--------------------------------------------------------------------------
    bool support_timesteps () const { return has_lev(); }
    bool support_SPH       () const { return has_siz(); }
    bool support_aux       () const { return has_aux(); }
    bool support_key       () const { return has_key(); }
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
      DATA_ACCESS_REF(real,mas)
      DATA_ACCESS_REF(real,mass)
      DATA_ACCESS_REF(vect,pos)
      real&pos  (const int i) { return pos()[i]; }
      DATA_ACCESS_REF(vect,vel)
      real&vel  (const int i) { return vel()[i]; }
      DATA_ACCESS_REF(vect,acc)
      real&acc  (const int i) { return acc()[i]; }
      DATA_ACCESS_REF(real,eps)
      DATA_ACCESS_REF(real,pot)
      DATA_ACCESS_REF(real,aux)
      DATA_ACCESS_REF(int, key)
      DATA_ACCESS_REF(real,rho)
      DATA_ACCESS_REF(uint,num)
      DATA_ACCESS_REF(indx,lev)
      DATA_ACCESS_REF(indx,level)
      DATA_ACCESS_REF(real,size)
      DATA_ACCESS_REF(flag,flg)
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
      void flag_as_sink      () { flg().add    (flag::SINK); }
      void unflag_sink       () { flg().un_set (flag::SINK); }
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
      friend bool may_go_shorter(const iterator&I) { return I.may_go_shorter(); }
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
    iterator back_bodies ()        const { return iterator(this,N_bodies()-1); }
    iterator end_SPH     ()        const { return iterator(this,N_sph()); }
    iterator back_SPH    ()        const { return iterator(this,N_sph()-1); }
  };
  //////////////////////////////////////////////////////////////////////////////
  // typedefs                                                                   
  //////////////////////////////////////////////////////////////////////////////
  typedef sbodies            bodies;               // define nbdy::bodies       
  typedef sbodies::body_type body;                 // define nbdy::body         
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
#endif                                             // included_body_h           
