// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// nmio.h                                                                      |
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
//-----------------------------------------------------------------------------+
#ifndef included_nmio_h
#define included_nmio_h

#ifndef included_auxx_h
#  include <public/auxx.h>
#endif

#ifdef  ALLOW_NEMO
//------------------------------------------------------------------------------

#ifdef TWODIMENSIONAL
#  define NDM 2
#  define TDM 4
#else
#  define NDM 3
#  define TDM 6
#endif

//------------------------------------------------------------------------------
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  // class nbdy::nemo_io                                                        
  //////////////////////////////////////////////////////////////////////////////
  class nemo_io {
  public:
    enum Set          {snap=0,param=1,
		       bodies=2,diags=3};          // type of nemo set          
    enum CoSys        {cart,spher,scat};           // type of coordinate system 
    enum SingleScalar {time,cputime};              // 1               real      
    enum SingleVector {energy};                    // 1 * Ndim        reals     
    enum SinglePhases {cofm};                      // 1 * 2 * Nsim    reals     
    enum SingleMatrix {KinT,PotT,AmT};             // 1 * Ndim * Ndim reals     
    enum BodiesScalar {mass,pot,rho,aux,eps};      // N               reals     
    enum BodiesVector {acc,pos,vel};               // N * Ndim        reals     
    enum BodiesPhases {posvel};                    // N * 2 * Ndim    reals     
    enum BodiesInteger{key,flag};                  // N               integers  
    enum BodiesShort  {level};                     // N               shorts    
  private:
    nemo_io(const nemo_io&);
    mutable bool   OPEN[4];
  protected:
    void          *STREAM;
    mutable int    N,CS;
    mutable float  SINGLEVECTOR[NDM], SINGLEMATRIX[NDM*NDM], SINGLEPHASES[TDM];
    mutable float *BODIESSCALAR, *BODIESARRAYS;
    mutable int   *BODIESINTEGER;
    mutable short *BODIESSHORT;
    //-------------------------------------------------------------------------+
    // Note on the memory access                                               |
    //                                                                         |
    // The idea is to first write the body data into the arrays and then use   |
    // the write() methods below to write them to the nemo file. For input,    |
    // just do the reverse: first read() them into the arrays and then read    |
    // them out in your own memory.                                            |
    // The arrays are allocated upon calling read_N() or write_N() and deleted |
    // upon calling reset() or the destructor.                                 |
    //                                                                         |
    // IMPORTANT NOTE                                                          |
    // There is only one array for N numbers and one for N*2*Ndim numbers.     |
    // That means you have to use the above scheme for each item in turn,      |
    // because otherwise read() or write() will overwrite the contents of      |
    // these arrays. However, we can hold simultaneously N masses and N        |
    // phase-space positions, the latter either as phases or as positions &    |
    // velocities, but not both at the same time                               |
    //-------------------------------------------------------------------------+
    nemo_io         () : STREAM(0),                // nemo_io w/o stream        
			 BODIESSCALAR (0),
    			 BODIESARRAYS (0),
    			 BODIESINTEGER(0),
    			 BODIESSHORT  (0) {}
    nemo_io         (const char*, const char*);    // nemo_io to file           
    ~nemo_io        ();                            // close file & free memory  
    void  open      (const char*, const char*);    // open new file (close old) 
    //--------------------------------------------------------------------------
    void  open_set  (const Set, const bool) const; // open a nemo set           
    void  close_set (const Set, const bool) const; // close a nemo set          
    //--------------------------------------------------------------------------
    bool  is_present(const Set)          const;    // can we open: nemo set?    
    bool  is_present(const SingleScalar) const;    // can we read: 1 scalar?    
    bool  is_present(const SingleVector) const;    // can we read: 1 vector?    
    bool  is_present(const SinglePhases) const;    // can we read: 1 phases?    
    bool  is_present(const SingleMatrix) const;    // can we read: 1 matrix?    
    bool  is_present(const BodiesScalar) const;    // can we read: many scalars?
    bool  is_present(const BodiesVector) const;    // can we read: many vectors?
    bool  is_present(const BodiesPhases) const;    // can we read: many phases? 
    bool  is_present(const BodiesInteger)const;    // can we read: many ints?   
    bool  is_present(const BodiesShort)  const;    // can we read: many shorts? 
    //--------------------------------------------------------------------------
    int   read_N    () const;                      // de-alloc, read N, set N   
    float read      (const SingleScalar,           // read 1          -> return 
		     float* =0)          const;    //[O: scalar read]           
    void  read      (const SingleVector,           // read Ndim       -> storage
		     float* =0)          const;    //[O: OR           -> here ] 
    void  read      (const SinglePhases,           // read 2*Ndim     -> storage
		     float* =0)          const;    //[O: OR           -> here ] 
    void  read      (const SingleMatrix,           // read Ndim*Ndim  -> storage
		     float* =0)          const;    //[O: OR           -> here ] 
    void  read      (const BodiesScalar,           // read N          -> array  
		     float* =0)          const;    //[O: OR           -> here ] 
    void  read      (const BodiesVector,           // read N*Ndim     -> array  
		     float* =0)          const;    //[O: OR           -> here ] 
    void  read      (const BodiesPhases,           // read N*2*Ndim   -> array  
		     float* =0)          const;    //[O: OR           -> here ] 
    void  read      (const BodiesInteger,          // read N          -> array  
		     int*   =0)          const;    //[O: OR           -> here ] 
    void  read      (const BodiesShort,            // read N          -> array  
		     short* =0)          const;    //[O: OR           -> here ] 
    void  read_history () const;                   // read nemo history         
    //--------------------------------------------------------------------------
    void  write_N   (const int)          const;    // de-alloc, write N, set N  
    void  write     (const CoSys)        const;    // write type: coord. system 
    void  write     (const SingleScalar,           // write 2nd arg -> 1        
		     const float)        const;    //                           
    void  write     (const SingleVector,           // write storage -> Ndim     
		     float* =0)          const;    //[O:    here    -> out    ] 
    void  write     (const SinglePhases,           // write storage -> 2*Ndim   
		     float* =0)          const;    //[O:    here    -> out    ] 
    void  write     (const SingleMatrix,           // write storage -> Ndim*Ndim
		     float* =0)          const;    //[O:    here    -> out    ] 
    void  write     (const BodiesScalar,           // write array   -> N        
		     float* =0)          const;    //[O:    here    -> out    ] 
    void  write     (const BodiesVector,           // write array   -> N*Ndim   
		     float* =0)          const;    //[O:    here    -> out    ] 
    void  write     (const BodiesPhases,           // write array   -> N*2*Ndim 
		     float* =0)          const;    //[O:    here    -> out    ] 
    void  write     (const BodiesInteger,          // write array   -> N        
		     int*   =0)          const;    //[O:    here    -> out    ] 
    void  write     (const BodiesShort,            // write array   -> N        
		     short* =0)          const;    //[O:    here    -> out    ] 
    void  write_history() const;                   // write nemo history        
  public:
    //--------------------------------------------------------------------------
    void  close     ();                            // close open file           
    //--------------------------------------------------------------------------
    bool  is_open   () const { return STREAM!=0; } // are we ready for output ? 
    //--------------------------------------------------------------------------
    void  allocscalar() const {                    // allocate BODYSCALAR       
      if(!BODIESSCALAR) { MemoryCheck(BODIESSCALAR = new float[N]); }
    }
    //--------------------------------------------------------------------------
    void  allocarrays() const {                    // allocate BODYARRAYS       
      if(!BODIESARRAYS) { MemoryCheck(BODIESARRAYS = new float[N*TDM]); }
    }
    //--------------------------------------------------------------------------
    void  allocinteger() const {                   // allocate BODYINTEGER      
      if(!BODIESINTEGER) { MemoryCheck(BODIESINTEGER = new int[N]); }
    }
    //--------------------------------------------------------------------------
    void  allocshort() const {                     // allocate BODYSHORT        
      if(!BODIESSHORT) { MemoryCheck(BODIESSHORT = new short[N]); }
    }
    //==========================================================================
    void  de_allocscalar() const {                 // de-allocate BODYSCALAR    
      if(BODIESSCALAR) { delete[] BODIESSCALAR; BODIESSCALAR=0; }
    }
    //--------------------------------------------------------------------------
    void  de_allocarrays() const {                 //de-allocate BODYARRAYS     
      if(BODIESARRAYS) { delete[] BODIESARRAYS; BODIESARRAYS=0; }
    }
    //--------------------------------------------------------------------------
    void  de_allocinteger() const {                // de-allocate BODYINTEGER   
      if(BODIESINTEGER) { delete[] BODIESINTEGER; BODIESINTEGER=0; }
    }
    //--------------------------------------------------------------------------
    void  de_allocshort() const {                  // allocate BODYSHORT        
      if(BODIESSHORT) { delete[] BODIESSHORT; BODIESSHORT=0; }
    }
    //==========================================================================
    void  reset      () const {                    // de-allocate arrays        
      de_allocscalar();
      de_allocarrays();
      de_allocinteger();
      de_allocshort();
    }
    //--------------------------------------------------------------------------
    bool  is_open_set(const Set S) const { return OPEN[S]; }
    //--------------------------------------------------------------------------
    const int& Number() const { return N; }
  };
  //////////////////////////////////////////////////////////////////////////////
  // class nbdy::nemo_in                                                        
  //////////////////////////////////////////////////////////////////////////////
#define c_i const int i
#define c_j const int j
  class nemo_in : public nemo_io {
  private:
    nemo_in(const nemo_in&);
  public:
    nemo_in() {}
    nemo_in(const char* file, const char* mode="r") : nemo_io(file,mode) {}
    //--------------------------------------------------------------------------
    nemo_io::is_present;
    nemo_io::read;
    nemo_io::read_N;
    nemo_io::read_history;
    //--------------------------------------------------------------------------
    void open (const char* file) { nemo_io::open(file,"r"); }
    //--------------------------------------------------------------------------
    void open_set (const Set S)   const { return nemo_io::open_set (S,true); }
    void close_set(const Set S)   const { return nemo_io::close_set(S,true); }
    //--------------------------------------------------------------------------
    const float& single_vec(c_i)     const { return SINGLEVECTOR[i]; }
    const float& single_mat(c_i,c_j) const { return SINGLEMATRIX[i*NDM+j]; }
    const float& single_phs(c_i,c_j) const { return SINGLEPHASES[i*NDM+j]; }
    const float& bodies_scl(c_i)     const { return BODIESSCALAR[i]; }
    const float* bodies_vec(c_i)     const { return BODIESARRAYS+NDM*i; }
    const float* bodies_vel(c_i)     const { return BODIESARRAYS+NDM*(N+i); }
    const float* bodies_phs(c_i)     const { return BODIESARRAYS+TDM*i; }
    const int  & bodies_int(c_i)     const { return BODIESINTEGER[i]; }
    const short& bodies_sht(c_i)     const { return BODIESSHORT[i]; }
  };
  //////////////////////////////////////////////////////////////////////////////
  // class nbdy::nemo_out                                                       
  //////////////////////////////////////////////////////////////////////////////
  class nemo_out : public nemo_io {
  private:
    nemo_out(const nemo_out&);
  public:
    nemo_out() {}
    nemo_out(const char* file, const char* mode="w") : nemo_io(file,mode) {}
    //--------------------------------------------------------------------------
    nemo_io::write;
    nemo_io::write_N;
    nemo_io::write_history;
    //--------------------------------------------------------------------------
    void open           (const char* file) { nemo_io::open(file,"w"); }
    void open_to_append (const char* file) { nemo_io::open(file,"a"); }
    //--------------------------------------------------------------------------
    void open_set (const Set S)   const { return nemo_io::open_set (S,false); }
    void close_set(const Set S)   const { return nemo_io::close_set(S,false); }
    //--------------------------------------------------------------------------
    float& single_vec(c_i)     { return SINGLEVECTOR[i]; }
    float& single_mat(c_i,c_j) { return SINGLEMATRIX[i*NDM+j]; }
    float& single_phs(c_i,c_j) { return SINGLEPHASES[i*NDM+j]; }
    float& bodies_scl(c_i)     { return BODIESSCALAR[i]; }
    float* bodies_vec(c_i)     { return BODIESARRAYS+NDM*i; }
    float* bodies_vel(c_i)     { return BODIESARRAYS+NDM*(N+i); }
    float* bodies_phs(c_i)     { return BODIESARRAYS+TDM*i; }
    int  & bodies_int(c_i)     { return BODIESINTEGER[i]; }
    short& bodies_sht(c_i)     { return BODIESSHORT[i]; }
  };
  //----------------------------------------------------------------------------
  bool time_in_range(const real&, const char*);
}
////////////////////////////////////////////////////////////////////////////////
#undef c_i
#undef c_j
#undef NDM
#undef TDM
////////////////////////////////////////////////////////////////////////////////
#endif // ALLOW_NEMO
#endif // included_nmio_h
