// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// nmio.h                                                                      |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_nmio_h
#define falcON_included_nmio_h

#ifndef falcON_included_auxx_h
#  include <public/auxx.h>
#endif

#ifdef  falcON_NEMO
//------------------------------------------------------------------------------
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::nemo_io                                                      //
  //                                                                          //
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
    enum SPHScalar    {uin,udin,udex,entr,srho,h}; // NS              reals     
    enum BodiesVector {acc,pos,vel,auxv};          // N * Ndim        reals     
    enum BodiesPhases {posvel};                    // N * 2 * Ndim    reals     
    enum BodiesInteger{key,flag,numb,numbSPH};     // N               integers  
    enum BodiesShort  {level};                     // N               shorts    
  private:
    nemo_io(const nemo_io&);
    mutable bool   OPEN[4];
  protected:
    static const int NDM = Ndim, TDM=2*NDM;
    void          *STREAM;
    mutable int    N,NS,CS;
    mutable real   SINGLEVECTOR[NDM], SINGLEMATRIX[NDM*NDM], SINGLEPHASES[TDM];
    mutable real  *BODIESSCALAR, *BODIESARRAYS;
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
    // Alternatively, if your data are stored in an array of reals, you may    |
    // read and write them directly, avoiding the copying. To this end, the    |
    // second argumend of the read() and write() functions is required.        |
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
    bool  is_present(const SPHScalar)    const;    // can we read: many scalars?
    bool  is_present(const BodiesVector) const;    // can we read: many vectors?
    bool  is_present(const BodiesPhases) const;    // can we read: many phases? 
    bool  is_present(const BodiesInteger)const;    // can we read: many ints?   
    bool  is_present(const BodiesShort)  const;    // can we read: many shorts? 
    //--------------------------------------------------------------------------
    void  read_N()                 const;          // de-alloc, read N & NS     
    real  read(const SingleScalar) const;          // read 1          -> return 
    void  read(const SingleVector) const;          // read Ndim       -> storage
    void  read(const SinglePhases) const;          // read 2*Ndim     -> storage
    void  read(const SingleMatrix) const;          // read Ndim*Ndim  -> storage
    void  read(const BodiesScalar) const;          // read N          -> array  
    void  read(const SPHScalar)    const;          // read NS         -> array  
    void  read(const BodiesVector) const;          // read N*Ndim     -> array  
    void  read(const BodiesPhases) const;          // read N*2*Ndim   -> array  
    void  read(const BodiesInteger)const;          // read N          -> array  
    void  read(const BodiesShort)  const;          // read N          -> array  
    //--------------------------------------------------------------------------
    void  read(const SingleScalar, real*) const;   // read 1          -> return 
    void  read(const SingleVector, real*) const;   // read Ndim       -> pter   
    void  read(const SinglePhases, real*) const;   // read 2*Ndim     -> pter   
    void  read(const SingleMatrix, real*) const;   // read Ndim*Ndim  -> pter   
    void  read(const BodiesScalar, real*) const;   // read N          -> pter   
    void  read(const SPHScalar,    real*) const;   // read NS         -> pter   
    void  read(const BodiesVector, real*) const;   // read N*Ndim     -> pter   
    void  read(const BodiesPhases, real*) const;   // read N*2*Ndim   -> pter   
    void  read(const BodiesInteger,int *) const;   // read N          -> pter   
    void  read(const BodiesShort, short*) const;   // read N          -> pter   
    //==========================================================================
    void  write_N (int const&,                     // de-alloc, write N, set N  
		   int const& = 0)   const;        //[I: N_sph]                 
    void  write(const CoSys)         const;        // write type: coord. system 
    void  write(const SingleScalar,
		real const&)         const;        // write 2nd arg -> 1        
    //--------------------------------------------------------------------------
    void  write(const SingleVector)  const;        // write storage -> Ndim     
    void  write(const SinglePhases)  const;        // write storage -> 2*Ndim   
    void  write(const SingleMatrix)  const;        // write storage -> Ndim*Ndim
    void  write(const BodiesScalar)  const;        // write array   -> N        
    void  write(const SPHScalar)     const;        // write array   -> NS       
    void  write(const BodiesVector)  const;        // write array   -> N*Ndim   
    void  write(const BodiesPhases)  const;        // write array   -> N*2*Ndim 
    void  write(const BodiesInteger) const;        // write array   -> N        
    void  write(const BodiesShort)   const;        // write array   -> N        
    //--------------------------------------------------------------------------
    void  write(const SingleVector, real*)  const; // write pter    -> Ndim     
    void  write(const SinglePhases, real*)  const; // write pter    -> 2*Ndim   
    void  write(const SingleMatrix, real*)  const; // write pter    -> Ndim*Ndim
    void  write(const BodiesScalar, real*)  const; // write pter    -> N        
    void  write(const SPHScalar,    real*)  const; // write pter    -> NS       
    void  write(const BodiesVector, real*)  const; // write pter    -> N*Ndim   
    void  write(const BodiesPhases, real*)  const; // write pter    -> N*2*Ndim 
    void  write(const BodiesInteger,int *)  const; // write pter    -> N        
    void  write(const BodiesShort, short*)  const; // write pter    -> N        
  public:
    //--------------------------------------------------------------------------
    void  close     ();                            // close open file           
    //--------------------------------------------------------------------------
    bool  is_open   () const { return STREAM!=0; } // are we ready for output ? 
    //--------------------------------------------------------------------------
    void  allocscalar() const {                    // allocate BODYSCALAR       
      if(!BODIESSCALAR) { BODIESSCALAR = falcON_New(real,max(N,NS)); }
    }
    //--------------------------------------------------------------------------
    void  allocarrays() const {                    // allocate BODYARRAYS       
      if(!BODIESARRAYS) { BODIESARRAYS = falcON_New(real,max(N,NS)*TDM); }
    }
    //--------------------------------------------------------------------------
    void  allocinteger() const {                   // allocate BODYINTEGER      
      if(!BODIESINTEGER) { BODIESINTEGER = falcON_New(int,max(N,NS)); }
    }
    //--------------------------------------------------------------------------
    void  allocshort() const {                     // allocate BODYSHORT        
      if(!BODIESSHORT) { BODIESSHORT = falcON_New(short,max(N,NS)); }
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
    int const&Number   () const { return N; }
    int const&NumberSPH() const { return NS; }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::nemo_in                                                      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
#define c_i int const&i
#define c_j int const&j
  class nemo_in : public nemo_io {
  private:
    nemo_in(nemo_in const&);
  public:
    nemo_in() {}
    nemo_in(const char* file, const char* mode="r") : nemo_io(file,mode) {}
    //--------------------------------------------------------------------------
    nemo_io::is_present;
    nemo_io::read;
    nemo_io::read_N;
    //--------------------------------------------------------------------------
    void open (const char* file) { nemo_io::open(file,"r"); }
    //--------------------------------------------------------------------------
    void open_set (const Set S)   const { return nemo_io::open_set (S,true); }
    void close_set(const Set S)   const { return nemo_io::close_set(S,true); }
    //--------------------------------------------------------------------------
    real const &single_vec(c_i)     const { return SINGLEVECTOR[i]; }
    real const &single_mat(c_i,c_j) const { return SINGLEMATRIX[i*NDM+j]; }
    real const &single_phs(c_i,c_j) const { return SINGLEPHASES[i*NDM+j]; }
    real const &bodies_scl(c_i)     const { return BODIESSCALAR[i]; }
    const real *bodies_vec(c_i)     const { return BODIESARRAYS+NDM*i; }
    const real *bodies_vel(c_i)     const { return BODIESARRAYS+NDM*(N+i); }
    const real *bodies_phs(c_i)     const { return BODIESARRAYS+TDM*i; }
    int   const&bodies_int(c_i)     const { return BODIESINTEGER[i]; }
    short const&bodies_sht(c_i)     const { return BODIESSHORT[i]; }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::nemo_out                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class nemo_out : public nemo_io {
  private:
    nemo_out(nemo_out const&);
  public:
    nemo_out() {}
    nemo_out(const char* file, const char* mode="w") : nemo_io(file,mode) {}
    //--------------------------------------------------------------------------
    nemo_io::write;
    nemo_io::write_N;
    //--------------------------------------------------------------------------
    void open           (const char* file) { nemo_io::open(file,"w"); }
    void open_to_append (const char* file) { nemo_io::open(file,"a"); }
    //--------------------------------------------------------------------------
    void open_set (const Set S)   const { return nemo_io::open_set (S,false); }
    void close_set(const Set S)   const { return nemo_io::close_set(S,false); }
    //--------------------------------------------------------------------------
    real &single_vec(c_i)     const { return SINGLEVECTOR[i]; }
    real &single_mat(c_i,c_j) const { return SINGLEMATRIX[i*NDM+j]; }
    real &single_phs(c_i,c_j) const { return SINGLEPHASES[i*NDM+j]; }
    real &bodies_scl(c_i)     const { return BODIESSCALAR[i]; }
    real *bodies_vec(c_i)     const { return BODIESARRAYS+NDM*i; }
    real *bodies_vel(c_i)     const { return BODIESARRAYS+NDM*(N+i); }
    real *bodies_phs(c_i)     const { return BODIESARRAYS+TDM*i; }
    int  &bodies_int(c_i)     const { return BODIESINTEGER[i]; }
    short&bodies_sht(c_i)     const { return BODIESSHORT[i]; }
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
#endif // falcON_NEMO
#endif // falcON_included_nmio_h
