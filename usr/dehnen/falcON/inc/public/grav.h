// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// grav.h                                                                      |
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
// class grav_mac                                                              |
// class grav_leaf                                                             |
// class grav_cell                                                             |
// class grav_estimator                                                        |
// class grav_stat                                                             |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_grav_h
#define falcON_included_grav_h

#ifndef falcON_included_deft_h
#  include <public/deft.h>
#endif
#ifndef falcON_included_tree_h
#  include <public/tree.h>
#endif
#ifndef falcON_included_ionl_h
#  include <public/ionl.h>
#endif
#ifndef falcON_included_memo_h
#  include <public/memo.h>
#endif
#ifndef falcON_included_tset_h
#  include <public/tset.h>
#endif
////////////////////////////////////////////////////////////////////////////////
#define ENHANCED_IACT_STATS
#undef  ENHANCED_IACT_STATS
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // namespace nbdy::grav                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  namespace grav {
    static const int ORDER = falcON_ORDER;         // expansion order is fixed  
    static const int D_DIM = ORDER+1;              // # terms in D[]            
    static const int P_ORD = ORDER-1;              // order of highest multipole
    typedef symset3D<ORDER,real>      Cset;        // set of Taylor coeffs      
    typedef poles3D <P_ORD,real>      Mset;        // set of multipoles         
    static const int NCOEF = Cset::NDAT;           // # reals in taylor coeffs  
  }
  //////////////////////////////////////////////////////////////////////////////
  class InvertZ;                                   // forward declaration       
  class grav_estimator;                            // forward declaration       
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::grav_mac                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class grav_mac {
    //--------------------------------------------------------------------------
    // data:                                                                    
    //--------------------------------------------------------------------------
  private:
    MAC_type          MAC;                         // type of MAC               
    uint              P;                           // expansion order           
    real              TH0, iTH0;                   // params for theta(Y)       
    InvertZ          *IZ;                          // inversion method          
    //--------------------------------------------------------------------------
    // private methods:                                                         
    //--------------------------------------------------------------------------
  private:
    real inv_theta_of_M(                           // R: 1/theta(M)             
			const real) const;         // I: log(M/M_tot)           
    //--------------------------------------------------------------------------
    // public methods:                                                          
    //--------------------------------------------------------------------------
  public:
    grav_mac  (                                    // constructor               
	       const MAC_type,                     // I: type of MAC            
	       const real,                         // I: parameter: theta_0     
	       const uint);                        // I: expansion order        
    //--------------------------------------------------------------------------
    void reset(const MAC_type,                     // I: type of MAC            
	       const real,                         // I: parameter: theta_0     
	       const uint);                        // I: expansion order        
    //--------------------------------------------------------------------------
    void reset_theta(const real);                  // I: parameter: theta_0     
    //--------------------------------------------------------------------------
    ~grav_mac  ();                                 // destructor                
    //--------------------------------------------------------------------------
    // const method                                                             
    //--------------------------------------------------------------------------
    void set_rcrit(const grav_estimator*) const;   // set rcrit for all cells   
    //--------------------------------------------------------------------------
    // const inlined methods                                                    
    //--------------------------------------------------------------------------
    real        theta_min() const { return TH0; }  // theta_0                   
    const MAC_type&method() const { return MAC; }  // MAC                       
    const char* describe_method() const            // describes MAC             
    { return describe(MAC); }
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::grav_leaf                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class grav_leaf : public basic_leaf {
    grav_leaf           (const grav_leaf&);        // not implemented           
    grav_leaf& operator=(const grav_leaf&);        // not implemented           
    //--------------------------------------------------------------------------
    // friends                                                                  
    //--------------------------------------------------------------------------
    friend class grav_estimator;
    //--------------------------------------------------------------------------
    // data of class grav_leaf (only static)                                    
    //--------------------------------------------------------------------------
    struct sink_data : public symset3D<1,real> {
      void reset() {        symset3D<1,real>::operator=(zero); }
      real&pot  () { return symset3D<1,real>::tensor<0>(); }
      vect&acc  () { return symset3D<1,real>::tensor<1>(); }
      uint&num  () { return *(static_cast<uint*>
			      (static_cast<void*>(&pot()))); }
    };
    //--------------------------------------------------------------------------
    // private data access                                                      
    //--------------------------------------------------------------------------
#define __sink  static_cast<sink_data*>(PROP)
    real      &mass ()       { return SCAL; }
    real const&mass () const { return SCAL; }
    vect const&acc  () const { return __sink->acc(); }
    vect const&cofm () const { return pos(); }
    real const&sizeq() const { return acc()[0]; }
    real const&pot  () const { return __sink->pot(); }
    real const&rho  () const { return pot(); }
    uint const&num  () const { return __sink->num(); }
#ifdef falcON_INDI
    real const&eph  () const { return AUXR; }
#endif
    //--------------------------------------------------------------------------
    // non-const data access via members                                        
    //--------------------------------------------------------------------------
  public:
    vect&acc  () { return __sink->acc(); }
    vect&cofm () { return pos(); }
    real&sizeq() { return acc()[0]; }
    real&pot  () { return __sink->pot(); }
    real&rho  () { return pot(); }
    uint&num  () { return __sink->num(); }
#ifdef falcON_INDI
    real&eph  () { return AUXR; }
#endif
    void inc  () { ++num(); }
    symset3D<1,real>& Coeffs() { return *__sink; }
    //--------------------------------------------------------------------------
    // const data access via friends                                            
    //--------------------------------------------------------------------------
    friend vect const&cofm  (const grav_leaf*const&L) {return L->pos(); }
    friend vect const&acc   (const grav_leaf*const&L) {return L->acc(); }
    friend real const&pot   (const grav_leaf*const&L) {return L->pot(); }
    friend real const&rho   (const grav_leaf*const&L) {return L->rho(); }
    friend uint const&num   (const grav_leaf*const&L) {return L->num(); }
#ifdef falcON_INDI
    friend real const&eph   (const grav_leaf*const&L) {return L->eph(); }
    friend real       size  (const grav_leaf*const&L) {return twice(L->eph()); }
    friend real const&sizeq (const grav_leaf*const&L) {return L->sizeq(); }
#endif
    friend real const&mass  (const grav_leaf*const&L) {return L->mass(); }
    //--------------------------------------------------------------------------
    // stuff needed for MPI code                                                
    //--------------------------------------------------------------------------
#ifdef falcON_MPI
    void copy_prune(const grav_leaf*L) {
      basic_leaf::copy_prune(L);
      mass() = nbdy::mass(L);
#ifdef falcON_INDI
      eph () = nbdy::eph (L);
#endif
    }
#endif
    //--------------------------------------------------------------------------
    // copy data from body to leaf                                              
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void copy_from_bodies_mass(const bodies_type*const&B) {
      mass() = B->mass(mybody());
    }
    //--------------------------------------------------------------------------
#ifdef falcON_INDI
    template<typename bodies_type>
    void copy_from_bodies_eph(const bodies_type*const&B) {
      eph() = half*B->eps(mybody());
    }
#endif
    //--------------------------------------------------------------------------
    // copy data to body from leaf                                              
    //--------------------------------------------------------------------------
#ifdef falcON_ADAP
    template<typename bodies_type>
    void copy_to_bodies_eps(const bodies_type*const&B) {
      B->eps(mybody()) = twice(eph());
    }
#endif
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void copy_to_bodies_rho(const bodies_type*const&B) const {
      B->rho(mybody()) = rho();
    }
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void copy_to_bodies_acc(const bodies_type*const&B) const {
      B->acc(mybody()) = acc();
    }
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void copy_to_bodies_pot(const bodies_type*const&B) const {
      B->pot(mybody()) = pot();
    }
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void copy_to_bodies_grav(const bodies_type*const&B) const {
      copy_to_bodies_pot(B);
      copy_to_bodies_acc(B);
    }
    //--------------------------------------------------------------------------
    // with non-unity constant G of gravity                                     
    template<typename bodies_type>
    void copy_to_bodies_acc(const bodies_type*const&B,
			    real              const&G) const {
      B->acc(mybody()).ass_times(acc(),G);
    }
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void copy_to_bodies_pot(const bodies_type*const&B,
			    real              const&G) const {
      B->pot(mybody()) = G * pot();
    }
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void copy_to_bodies_grav(const bodies_type*const&B,
			     real              const&G) const {
      copy_to_bodies_pot(B,G);
      copy_to_bodies_acc(B,G);
    }
    //--------------------------------------------------------------------------
    // reset body gravity data (needed if G=0)                                  
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void reset_bodies_acc(const bodies_type*const&B) const {
      B->acc(mybody()) = zero;
    }
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void reset_bodies_pot(const bodies_type*const&B) const {
      B->pot(mybody()) = zero;
    }
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void reset_bodies_grav(const bodies_type*const&B) const {
      reset_bodies_pot(B);
      reset_bodies_acc(B);
    }
    //--------------------------------------------------------------------------
    // simple manipulations                                                     
    //--------------------------------------------------------------------------
    void set_sink  (sink_data*const&sink) { PROP = static_cast<void*>(sink); }
    void reset_sink()                     { __sink->reset(); }
    //--------------------------------------------------------------------------
    void normalize_grav () {                       // acc,pot     /= mass       
      if(mass()>zero) {
	register real im = one/mass();
	pot() *= im;
	acc() *= im;
      }
    }
    //--------------------------------------------------------------------------
    // boolean information                                                      
    //--------------------------------------------------------------------------
    friend bool is_source(const grav_leaf*const&S) { return S->mass() != zero; }
    //--------------------------------------------------------------------------
    // dump leaf data                                                           
    //--------------------------------------------------------------------------
    static void dump_head(std::ostream& o) {
      basic_leaf::dump_head(o);
      o<<"              mass";
    }
    //--------------------------------------------------------------------------
    void dump(std::ostream&o) const {
      basic_leaf::dump(o);
      o<<' '<<setw(8)<<mass();
    }
    //--------------------------------------------------------------------------
  };
#undef __sink
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // nbdy::class grav_cell                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class grav_cell : public basic_cell {
    grav_cell           (const grav_cell&);        // not implemented           
    grav_cell& operator=(const grav_cell&);        // not implemented           
    //--------------------------------------------------------------------------
    // friendships                                                              
    //--------------------------------------------------------------------------
    friend class grav_estimator;                   // for alloc of srce_data    
    //--------------------------------------------------------------------------
    // types and static data                                                    
    //--------------------------------------------------------------------------
  public:
    typedef grav_leaf leaf_type;                   // type of associated leafs  
    typedef grav::Mset Mset;
    typedef grav::Cset Cset;
    //--------------------------------------------------------------------------
    static const int N_SINK = Cset::NDAT;
    //--------------------------------------------------------------------------
  private:
    struct srce_data {
      real MASS;
      vect COFM;
#ifdef falcON_INDI
      real EPH;
#endif
      Mset POLS;
      void normalize_poles()         { POLS.normalize(MASS); }
    };
    //--------------------------------------------------------------------------
    // private data access                                                      
    //--------------------------------------------------------------------------
#define __srce  static_cast<srce_data*>(AUX1.PTER)
    real const&mass  () const { return __srce->MASS; }
    vect const&cofm  () const { return __srce->COFM; }
    real const&rmax  () const { return AUX3.SCAL; }
    real const&size  () const { return AUX3.SCAL; }
    real const&rcrit () const { return AUX3.SCAL; }
    real const rcrit2() const { return square(rcrit()); }
#ifdef falcON_INDI
    real const&eph   () const { return __srce->EPH; }
#endif
          Cset&Coeffs() const { return *static_cast<Cset*>(AUX2.PTER); }
    const Mset&poles () const { return __srce->POLS; }
    //--------------------------------------------------------------------------
    // simple manipulations                                                     
    //--------------------------------------------------------------------------
  public:
    void set_rcrit   (real const&it) { AUX3.SCAL *= it; }
    //--------------------------------------------------------------------------
    void set_srce    (srce_data*const&srce)
    {
      AUX1.PTER = static_cast<void*>(srce);
    }
    //--------------------------------------------------------------------------
    void setCoeffs   (Cset     *const&sink)
    {
      AUX2.PTER = static_cast<void*>(sink); 
    }
    //--------------------------------------------------------------------------
    void*returnCoeffs()                     { return AUX2.PTER; }
    void resetCoeffs ()                     { AUX2.PTER=0; }
    bool hasCoeffs   () const               { return AUX2.PTER != 0; }
    //--------------------------------------------------------------------------
    // non-const data access via members                                        
    //--------------------------------------------------------------------------
    real &mass  () { return __srce->MASS; }
    vect &cofm  () { return __srce->COFM; }
    real &rmax  () { return AUX3.SCAL; }
    real &size  () { return AUX3.SCAL; }
#ifdef falcON_INDI
    real &eph   () { return __srce->EPH; }
#endif
    Cset &Coeffs() { return *static_cast<Cset*>(AUX2.PTER); }
    Mset &poles () { return __srce->POLS; }
#undef __srce
    //--------------------------------------------------------------------------
    // const data access via friends                                            
    //--------------------------------------------------------------------------
    friend real const&mass  (const grav_cell*const&C) { return C->mass(); }
    friend vect const&cofm  (const grav_cell*const&C) { return C->cofm(); }
    friend vect const&pos   (const grav_cell*const&C) { return C->cofm(); }
    friend real const&rmax  (const grav_cell*const&C) { return C->rmax(); }
    friend real const&rcrit (const grav_cell*const&C) { return C->rcrit(); }
    friend real const&size  (const grav_cell*const&C) { return C->size(); }
    friend real const rcrit2(const grav_cell*const&C) { return C->rcrit2(); }
#ifdef falcON_INDI
    friend real const&eph   (const grav_cell*const&C) { return C->eph(); }
#endif
    friend Cset      &Coeffs(const grav_cell*const&C) { return C->Coeffs(); }
    friend Mset const&poles (const grav_cell*const&C) { return C->poles();}
    friend bool hasCoeffs   (const grav_cell*const&C) { return C->hasCoeffs();}
    //--------------------------------------------------------------------------
    // boolean information via friends                                          
    //--------------------------------------------------------------------------
    friend bool is_source   (const grav_cell*const&C) { 
      return nbdy::mass(C)!=zero; }
    //--------------------------------------------------------------------------
    // other const methods and friends                                          
    //--------------------------------------------------------------------------
    friend real xmin(const grav_cell*const&C) {
      return nbdy::cofm(C).min() - nbdy::rmax(C); }
    friend real xmax(const grav_cell*const&C) {
      return nbdy::cofm(C).max() + nbdy::rmax(C); }
    //--------------------------------------------------------------------------
    // dump cell data                                                           
    //--------------------------------------------------------------------------
    static void dump_head(std::ostream&o) {
      basic_cell::dump_head(o);
      o<<
	"              mass"
	"              cofm         "
	"         rmax"
	"        rcrit";
    }
    //--------------------------------------------------------------------------
    void dump(std::ostream&o) const {
      basic_cell::dump(o);
      o<<' '<<setw(8)<<mass();
      for(register int i=0; i<Ndim; i++)
	o<<' '<<setw(8)<<setprecision(4)<<cofm()[i];
      o<<' '<<setw(12)<<rmax()
       <<' '<<setw(12)<<rcrit();
    }
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::grav_estimator                                               //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class grav_stat;                                 // forward declaration       
  class grav_estimator {
    //--------------------------------------------------------------------------
    grav_estimator            (grav_estimator const&);
    grav_estimator& operator= (grav_estimator const&);
    //--------------------------------------------------------------------------
    // data:                                                                    
    //--------------------------------------------------------------------------
    const oct_tree       *TREE;                    // the tree to be used       
    bool                  CELLS_UPTODATE;          // are cell srces up to date?
    bool                  LEAFS_UPTODATE;          // are leaf srces up to date?
#ifdef falcON_INDI
    const bool            INDI_SOFT;               // use individual eps_i      
#endif
    const int             DIR[4];                  // direct loop control       
    kern_type             KERNEL;                  // softening kernel          
    grav_stat            *STATS;                   // interaction statistics    
    real                  EPS;                     // global softening length   
    real                  GRAV;                    // Newton's G (can be 0)     
    uint                  Ncoeffs,Nchunks,Ncsize;  // # coeffs, chunks, bytes/c 
    grav_cell::srce_data *CELL_SRCE;               // memory for cell srce      
#ifdef falcON_MPI
    grav::Cset           *CELL_COEF;               // memory for cell coeffs    
#endif
    grav_leaf::sink_data *LEAF_SINK;               // memory for leafs          
    uint                  NCT, NCA, NLA;           // # allocation of these     
    uint                  NLA_needed;              // # active leafs            
    //--------------------------------------------------------------------------
    // private methods                                                          
    //--------------------------------------------------------------------------
#ifdef falcON_MPI
    inline void   set_cell_sink(bool const& =0);   // for parallel code only    
#endif
    //--------------------------------------------------------------------------
    // update leafs' source fields (requires no allocation here)                
    // - mass, leaf flags[, eps/2]                                              
    // - count active leafs                                                     
    void update_leafs();
    //--------------------------------------------------------------------------
    // adjust leafs' eph_i                                                      
#ifdef falcON_ADAP
    void adjust_eph(bool const&,                   // I: for all or active only 
		    real const&,                   // I: Nsoft: adjust actives  
		    real const&,                   // I: emin:  adjust actives  
		    real const&,                   // I: emax:  adjust actives  
		    uint const&,                   // I: Nref:  adjust actives  
		    real const&);                  // I: max change in eps_i    
#endif
    //--------------------------------------------------------------------------
    // - passes up the tree: flag, mass, cofm, rmax[, eph], multipoles          
    // - set the cells r_crit & r_crit^2                                        
    uint pass_up(                                  // R: # active cells         
		 const grav_mac*const&,            // I: MAC                    
		 bool           const&);           // I: reused old tree?       
    //--------------------------------------------------------------------------
    // prepare for interactions                                                 
    //  - allocate memory for leaf sink properties for active leafs             
    //  - allocate memory for cell source properties & reset their coeff pter   
    //  - pass source properties up the tree, count active cells                
    // [- allocate memory for cell sink properties, (falcON_MPI only)]          
    bool prepare(                                  // R: ALL || all are active  
		 const grav_mac*const&,            // I: MAC                    
		 bool           const&,            // I: all or active only?    
		 bool           const&);           // I: allocate cell coeffs?  
    //--------------------------------------------------------------------------
    // tree stuff to be superseeded                                             
    //--------------------------------------------------------------------------
  public:
    typedef grav_cell                       cell_type;
    typedef grav_leaf                       leaf_type;
    typedef oct_tree::CellIter<grav_cell>   cell_iterator;
    typedef leaf_type*                      leaf_iterator;
    //--------------------------------------------------------------------------
    const oct_tree*const&my_tree() const { return TREE; }
    cell_iterator root          () const {
      return cell_iterator(TREE,static_cast<grav_cell*>(TREE->FstCell())); }
    //--------------------------------------------------------------------------
    // dump cell and leaf data                                                  
    //--------------------------------------------------------------------------
    void dump_cells(std::ostream&) const;
    void dump_leafs(std::ostream&) const;
    //--------------------------------------------------------------------------
    // construction: allocate memory for leafs' & cells' source properties      
    //--------------------------------------------------------------------------
    grav_estimator(const oct_tree*const&T,         // I: tree to be used        
		   kern_type      const&k,         // I: kernel to be used      
		   grav_stat*     const&st,        // I: statistics             
		   real           const&e,         // I: global/max eps         
		   real           const&g = one,   //[I: Newton's G]            
		   bool           const&s = 0,     //[I: use individual eps?]   
		   const int d[4]=Default::direct)://[I: N_direct for gravity]  
      TREE           ( T ),
      LEAFS_UPTODATE ( 0 ),
      CELLS_UPTODATE ( 0 ),
#ifdef falcON_INDI
      INDI_SOFT      ( s ),
#endif
      KERNEL         ( k ),
      STATS          ( st ),
      EPS            ( e ),
      GRAV           ( g ),
      Ncoeffs        ( 0u ),
      Nchunks        ( 0u ),
      Ncsize         ( 0u ),
      CELL_SRCE      ( 0 ),
#ifdef falcON_MPI
      CELL_COEF      ( 0 ), 
#endif
      LEAF_SINK      ( 0 ),
      NCT            ( 0u ),
      NCA            ( 0u ),
      NLA            ( 0u ),
      NLA_needed     ( 0u )
    {
      const_cast<int*>(DIR)[0]=d[0];
      const_cast<int*>(DIR)[1]=d[1];
      const_cast<int*>(DIR)[2]=d[2];
      const_cast<int*>(DIR)[3]=d[3];
    }
    //--------------------------------------------------------------------------
    // reset allocations                                                        
    void reset() {
#ifdef falcON_MPI
      if(CELL_COEF) { delete[] CELL_COEF; CELL_COEF=0; }
#endif
      LEAFS_UPTODATE = 0;
      CELLS_UPTODATE = 0;
    }
    //--------------------------------------------------------------------------
    // change tree entry, set flags                                             
    void new_tree(const oct_tree*const&T) {
      TREE = T;
      reset();
    }
    //--------------------------------------------------------------------------
    void reset_softening(const real      e,
			 const kern_type k) {
      EPS    = e;
      KERNEL = k;
    }
    //--------------------------------------------------------------------------
    void reset_NewtonsG(const real g) {
      GRAV = g;
    }
    //--------------------------------------------------------------------------
    // destruction                                                              
    //--------------------------------------------------------------------------
    ~grav_estimator() {
      if(CELL_SRCE) delete[] CELL_SRCE;
      if(LEAF_SINK) delete[] LEAF_SINK;
#ifdef falcON_MPI
      if(CELL_COEF) delete[] CELL_COEF;
#endif
    }
    //--------------------------------------------------------------------------
    // compute gravity by direct summation                                      
    //--------------------------------------------------------------------------
    void exact(bool      const& =false             //[I: for all or active only]
#ifdef falcON_ADAP
	      ,real      const& =zero,             //[I: Nsoft: adjust eps_i]   
	       uint      const& =0u,               //[I: Nref:  adjust eps_i]   
	       real      const& =zero,             //[I: eps_min]               
	       real      const& =zero              //[I: max change of eps]     
#endif
	       );
    //--------------------------------------------------------------------------
    // compute gravity by the approximate method of Dehnen (2002)               
    //--------------------------------------------------------------------------
    void approx(const grav_mac*const&,             // I: MAC                    
		bool           const& =false,      //[I: for all or active only]
		bool           const& =true        //[I: combine phases]        
#ifdef falcON_ADAP
		,
		real           const& =zero,       //[I: Nsoft: adjust eps_i]   
		uint           const& =0u,         //[I: Nref:  adjust eps_i]   
		real           const& =zero,       //[I: eps_min]               
		real           const& =zero        //[I: max change of eps]     
#endif
		);
    //--------------------------------------------------------------------------
    // density estimation (number-, surface-, mass-density)                     
    //--------------------------------------------------------------------------
    void estimate_nd(bool const&,                  // I: all or active only?    
		     uint const&) const;           // I: critical cell size     
    void estimate_sd(bool const&,                  // I: all or active only?    
		     uint const&);                 // I: critical cell size     
    void estimate_md(bool const&,                  // I: all or active only?    
		     uint const&);                 // I: critical cell size     
    //--------------------------------------------------------------------------
    // public const data access                                                 
    //--------------------------------------------------------------------------
    const oct_tree*const&tree            () const { return TREE; }
    kern_type      const&kernel          () const { return KERNEL; }
    grav_stat     *const&stats           () const { return STATS; }
    real           const&NewtonsG        () const { return GRAV; }
    bool                 ZeroGravity     () const { return GRAV==zero; }
    real           const&softening_length() const { return EPS; }
    uint           const&N_active_cells  () const { return NCA; }
    uint           const&N_active_leafs  () const { return NLA; }
    uint           const&N_coeffs        () const { return Ncoeffs; }
    uint           const&N_chunks        () const { return Nchunks; }
    uint           const&N_elems_in_chunk() const { return Ncsize; }
#ifdef falcON_INDI
    bool           const&use_indiv_eps   () const { return INDI_SOFT; }
#endif
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // namespace nbdy::grav                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  namespace grav {
    typedef grav_cell                      cell_type;
    typedef grav_leaf                      leaf_type;
    typedef grav_cell                     *cell_pter;
    typedef grav_leaf                     *leaf_pter;
    typedef grav_estimator::cell_iterator  cell_iter;
    typedef grav_estimator::leaf_iterator  leaf_iter;
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class grav_stat                                                          //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class grav_stat {
  private:
    typedef grav::cell_iter cell_iter;
    typedef grav::leaf_iter leaf_iter;
    uint D_BB, D_CB, D_CC, D_CX;                   // # direct interactions     
    uint A_CB, A_CC;                               // # approximate interactions
#ifdef ENHANCED_IACT_STATS
    uint P_CB, P_CC, P_CX;                         // # BB pairs in direct      
#  define ADD_SS(COUNTER,NUMBER)  COUNTER += NUMBER;
    static const char* trick_SS(uint const&n, int&w)
    {
      if(n < 10) { w=1; return "         ("; }
      if(n < 100) { w=2; return "        ("; }
      if(n < 1000) { w=3; return "       ("; }
      if(n < 10000) { w=4; return "      ("; }
      if(n < 100000) { w=5; return "     ("; }
      if(n < 1000000) { w=6; return "    ("; }
      if(n < 10000000) { w=7; return "   ("; }
      if(n < 100000000) { w=8; return "  ("; }
      if(n < 1000000000) { w=9; return " ("; }
      w=10; return "(";
    }
#else
#  define ADD_SS(COUNTER,NUMBER) 
#endif
    static int trick(real const&x, int const&w)
    {
      if(x<0.001)  return max(1,w-5);
      if(x<0.01)   return max(1,w-4);
      if(x<0.1)    return max(1,w-3);
      if(x<one)    return max(1,w-2);
      if(x<ten)    return max(1,w-1);
      return w;
    }
  public:
    // 1 recording                                                              
    void reset    () {
      D_BB=0, D_CB=0, D_CC=0, D_CX=0;
      A_CB=0, A_CC=0;
#ifdef ENHANCED_IACT_STATS
      P_CB=0, P_CC=0, P_CX=0;
#endif
    }
    void record_BB() { ++D_BB; }
    void record_approx_CB(cell_iter const&A, leaf_iter const&B)
    { ++A_CB; }
    void record_approx_CC(cell_iter const&A, cell_iter const&B)
    { ++A_CC; }
    void record_direct_CB(cell_iter const&A, leaf_iter const&B)
    { ++D_CB; ADD_SS(P_CB, number(A)) }
    void record_direct_CC(cell_iter const&A, cell_iter const&B)
    { ++D_CC; ADD_SS(P_CC, number(A)*number(B)) }
    void record_direct_CX(cell_iter const&A)
    { ++D_CX; ADD_SS(P_CX, (number(A)*(number(A)-1))>>1 ) }
#ifdef ENHANCED_IACT_STATS
#  undef ADD_SS
#endif
    // 2 reporting                                                              
    uint const&BB_direct_iacts   () const { return D_BB; }
    uint const&CB_direct_iacts   () const { return D_CB; }
    uint const&CC_direct_iacts   () const { return D_CC; }
    uint const&CX_direct_iacts   () const { return D_CX; }
    uint       total_direct_iacts() const { return
					      BB_direct_iacts() +
					      CB_direct_iacts() +
					      CC_direct_iacts() +
					      CX_direct_iacts();
    }
#ifdef ENHANCED_IACT_STATS
    uint const&BB_direct_pairs   () const { return D_BB; }
    uint const&CB_direct_pairs   () const { return P_CB; }
    uint const&CC_direct_pairs   () const { return P_CC; }
    uint const&CX_direct_pairs   () const { return P_CX; }
    uint       total_direct_pairs() const { return
					      BB_direct_pairs() +
					      CB_direct_pairs() +
					      CC_direct_pairs() +
					      CX_direct_pairs();
    }
#endif
    uint       BB_approx_iacts   () const { return 0u; }
    uint const&CB_approx_iacts   () const { return A_CB; }
    uint const&CC_approx_iacts   () const { return A_CC; }
    uint       CX_approx_iacts   () const { return 0u; }
    uint       total_approx_iacts() const { return
					      BB_approx_iacts() +
					      CB_approx_iacts() +
					      CC_approx_iacts() +
					      CX_approx_iacts();
    }
    // 3 writing stats to ostream                                               
    void write(std::ostream&out) const {
      register uint part, total=total_approx_iacts()+total_direct_iacts();
#ifdef ENHANCED_IACT_STATS
      register int wSS;
#endif
      register real percent;
      out <<
	" interaction statitics:\n"
	"     type          approx   direct"
#ifdef ENHANCED_IACT_STATS
	"      (pairs)"
#endif
	"      total\n"
	" # body-body :          - ";
      part    = BB_direct_iacts();
      percent = 100.*part/real(total);
      out<<setw( 8)<<BB_direct_iacts()<<' '
#ifdef ENHANCED_IACT_STATS
	 <<trick_SS(BB_direct_pairs(),wSS)<<setw(wSS)<<BB_direct_pairs()<<") "
#endif
	 <<setw(10)<<part<<" = "
	 <<setprecision(trick(percent,5))<<setw(8)<<percent<<"%\n"
	 <<" # cell-body : ";
      part    = CB_approx_iacts()+CB_direct_iacts();
      percent = 100.*part/real(total);
      out<<setw(10)<<CB_approx_iacts()<<" "
	 <<setw( 8)<<CB_direct_iacts()<<" "
#ifdef ENHANCED_IACT_STATS
	 <<trick_SS(CB_direct_pairs(),wSS)<<setw(wSS)<<CB_direct_pairs()<<") "
#endif
	 <<setw(10)<<part<<" = "
	 <<setprecision(trick(percent,5))<<setw(8)<<percent<<"%\n"
	 <<" # cell-cell : ";
      part    = CC_approx_iacts() + CC_direct_iacts();
      percent = 100.*part/real(total);
      out<<setw(10)<<CC_approx_iacts()<<" "
	 <<setw( 8)<<CC_direct_iacts()<<" "
#ifdef ENHANCED_IACT_STATS
	 <<trick_SS(CC_direct_pairs(),wSS)<<setw(wSS)<<CC_direct_pairs()<<") "
#endif
	 <<setw(10)<<part<<" = "
	 <<setprecision(trick(percent,5))<<setw(8)<<percent<<"%\n"
	 <<" # cell-self :          - ";
      part    = CX_direct_iacts();
      percent = 100.*part/real(total);
      out<<setw( 8)<<CX_direct_iacts()<<" "
#ifdef ENHANCED_IACT_STATS
	 <<trick_SS(CX_direct_pairs(),wSS)<<setw(wSS)<<CX_direct_pairs()<<") "
#endif
	 <<setw(10)<<part<<" = "
	 <<setprecision(trick(percent,5))<<setw(8)<<percent<<"%\n"
	 <<" # total     : ";
      out<<setw(10)<<total_approx_iacts()<<" "
	 <<setw( 8)<<total_direct_iacts()<<" "
#ifdef ENHANCED_IACT_STATS
	 <<trick_SS(total_direct_pairs(),wSS)
	 <<setw(wSS)<<total_direct_pairs()<<") "
#endif
	 <<setw(10)<<total<<" =  100.000%\n";
    }
  };
}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
#undef  ENHANCED_IACT_STATS
#endif                                             // falcON_included_grav_h    
