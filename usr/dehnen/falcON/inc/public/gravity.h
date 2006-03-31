// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// gravity.h                                                                   |
//                                                                             |
// Copyright (C) 2000-2006  Walter Dehnen                                      |
//                                                                             |
// This program is free software; you can redistribute it and/or modify        |
// it under the terms of the GNU General Public License as published by        |
// the Free Software Foundation; either version 2 of the License, or (at       |
// your option) any later version.                                             |
//                                                                             |
// This program is distributed in the hope that it will be useful, but         |
// WITHOUT ANY WARRANTY; without even the implied warranty of                  |
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU           |
// General Public License for more details.                                    |
//                                                                             |
// You should have received a copy of the GNU General Public License           |
// along with this program; if not, write to the Free Software                 |
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                   |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// defines                                                                     |
//                                                                             |
// class GravMAC                                                               |
// class GravEstimator                                                         |
// class GravStats                                                             |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_gravity_h
#define falcON_included_gravity_h

#ifndef falcON_included_tree_h
#  include <public/tree.h>
#endif
#ifndef falcON_included_tensor_set_h
#  include <public/tensor_set.h>
#endif
////////////////////////////////////////////////////////////////////////////////
#define ENHANCED_IACT_STATS
#undef  ENHANCED_IACT_STATS
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // namespace falcON::grav                                                   //
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
  class GravEstimator;                             // forward declaration       
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::GravMAC                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class GravMAC {
    //--------------------------------------------------------------------------
    // data:                                                                    
    //--------------------------------------------------------------------------
  private:
    MAC_type          MAC;                         // type of MAC               
    unsigned          P;                           // expansion order           
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
    GravMAC   (                                    // constructor               
	       MAC_type,                           // I: type of MAC            
	       real,                               // I: parameter: theta_0     
	       unsigned);                          // I: expansion order        
    //--------------------------------------------------------------------------
    void reset(MAC_type,                           // I: type of MAC            
	       real,                               // I: parameter: theta_0     
	       unsigned);                          // I: expansion order        
    //--------------------------------------------------------------------------
    void reset_theta(real);                        // I: parameter: theta_0     
    //--------------------------------------------------------------------------
    ~GravMAC   ();                                 // destructor                
    //--------------------------------------------------------------------------
    // const method                                                             
    //--------------------------------------------------------------------------
    void set_rcrit(const GravEstimator*) const;    // set rcrit for all cells   
    //--------------------------------------------------------------------------
    // const inlined methods                                                    
    //--------------------------------------------------------------------------
    real const &theta_min() const { return TH0; }  // theta_0                   
    MAC_type const&method() const { return MAC; }  // MAC                       
    const char* describe_method() const            // describes MAC             
    { return describe(MAC); }
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::GravEstimator                                              //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class GravStats;                                 // forward declaration       
  class GravEstimator {
    //--------------------------------------------------------------------------
    GravEstimator            (GravEstimator const&);
    GravEstimator& operator= (GravEstimator const&);
    //--------------------------------------------------------------------------
    //                                                                          
    // sub-type Leaf                                                            
    //                                                                          
    //--------------------------------------------------------------------------
  public:
    class Leaf : public OctTree::Leaf {
      Leaf           (const Leaf&);                // not implemented           
      Leaf& operator=(const Leaf&);                // not implemented           
      //------------------------------------------------------------------------
      // data of class Leaf (only static)                                       
      //------------------------------------------------------------------------
#if defined(__GNUC__) && (__GNUC__ < 3 || __GNUC__ == 3 && __GNUC_MINOR__ < 4)
      // patch to fix a bug with gcc version < 3.4
    public:
#endif
      struct sink_data : public symset3D<1,real> {
	void     reset() {        symset3D<1,real>::operator=(zero); }
	real    &pot  () { return symset3D<1,real>::tensor<0>(); }
	vect    &acc  () { return symset3D<1,real>::tensor<1>(); }
	unsigned&num  () { return *(static_cast<unsigned*>
				   (static_cast<void*>(&pot()))); }
      };
      //------------------------------------------------------------------------
      // friends                                                                
      //------------------------------------------------------------------------
      friend class GravEstimator;
      friend class falcON::traits<sink_data>;
      //------------------------------------------------------------------------
      // private const data access                                              
      //------------------------------------------------------------------------
      sink_data*       sink () const { return static_cast<sink_data*>(PROP); }
      real            &mass ()       { return SCAL; }
      real       const&mass () const { return SCAL; }
      vect       const&acc  () const { return sink()->acc(); }
      vect       const&cofm () const { return pos(); }
      real       const&sizeq() const { return acc()[0]; }
      real       const&pot  () const { return sink()->pot(); }
      real       const&rho  () const { return pot(); }
      unsigned   const&num  () const { return sink()->num(); }
#ifdef falcON_INDI
      real       const&eph  () const { return AUXR; }
#endif
      //------------------------------------------------------------------------
      // non-const data access via members                                      
      //------------------------------------------------------------------------
    public:
      vect    &acc  () { return sink()->acc(); }
      vect    &cofm () { return pos(); }
      real    &sizeq() { return acc()[0]; }
      real    &pot  () { return sink()->pot(); }
      real    &rho  () { return pot(); }
      unsigned&num  () { return sink()->num(); }
#ifdef falcON_INDI
      real    &eph  () { return AUXR; }
#endif
      void     inc  () { ++num(); }
      symset3D<1,real>& Coeffs() { return *(sink()); }
      //------------------------------------------------------------------------
      // const data access via friends                                          
      //------------------------------------------------------------------------
      friend vect     const&cofm  (const Leaf*);
      friend vect     const&acc   (const Leaf*);
      friend real     const&pot   (const Leaf*);
      friend real     const&rho   (const Leaf*);
      friend unsigned const&num   (const Leaf*);
#ifdef falcON_INDI
      friend real     const&eph   (const Leaf*);
      friend real           size  (const Leaf*);
      friend real     const&sizeq (const Leaf*);
#endif
      friend real     const&mass  (const Leaf*);
      //------------------------------------------------------------------------
      // copy data from body to leaf                                            
      //------------------------------------------------------------------------
      template<typename bodies_type>
      void copy_from_bodies_mass(const bodies_type*const&B) {
	mass() = B->mass(mybody());
      }
      //------------------------------------------------------------------------
#ifdef falcON_INDI
      template<typename bodies_type>
      void copy_from_bodies_eph(const bodies_type*const&B) {
	eph() = half*B->eps(mybody());
      }
#endif
      //------------------------------------------------------------------------
      // copy data to body from leaf                                            
      //------------------------------------------------------------------------
#ifdef falcON_ADAP
      template<typename bodies_type>
      void copy_to_bodies_eps(const bodies_type*const&B) {
	B->eps(mybody()) = twice(eph());
      }
#endif
      //------------------------------------------------------------------------
      template<typename bodies_type>
      void copy_to_bodies_rho(const bodies_type*const&B) const {
	B->rho(mybody()) = rho();
      }
      //------------------------------------------------------------------------
      template<typename bodies_type>
      void copy_to_bodies_acc(const bodies_type*const&B) const {
	B->acc(mybody()) = acc();
      }
      //------------------------------------------------------------------------
      template<typename bodies_type>
      void copy_to_bodies_pot(const bodies_type*const&B) const {
	B->pot(mybody()) = pot();
      }
      //------------------------------------------------------------------------
      template<typename bodies_type>
      void copy_to_bodies_grav(const bodies_type*const&B) const {
	copy_to_bodies_pot(B);
	copy_to_bodies_acc(B);
      }
      //------------------------------------------------------------------------
      // with non-unity constant G of gravity                                   
      template<typename bodies_type>
      void copy_to_bodies_acc(const bodies_type*const&B,
			      real              const&G) const {
	B->acc(mybody()) = G * acc();
      }
      //------------------------------------------------------------------------
      template<typename bodies_type>
      void copy_to_bodies_pot(const bodies_type*const&B,
			      real              const&G) const {
	B->pot(mybody()) = G * pot();
      }
      //------------------------------------------------------------------------
      template<typename bodies_type>
      void copy_to_bodies_grav(const bodies_type*const&B,
			       real              const&G) const {
	copy_to_bodies_pot(B,G);
	copy_to_bodies_acc(B,G);
      }
      //------------------------------------------------------------------------
      // reset body gravity data (needed if G=0)                                
      //------------------------------------------------------------------------
      template<typename bodies_type>
      void reset_bodies_acc(const bodies_type*const&B) const {
	B->acc(mybody()) = zero;
      }
      //------------------------------------------------------------------------
      template<typename bodies_type>
      void reset_bodies_pot(const bodies_type*const&B) const {
	B->pot(mybody()) = zero;
      }
      //------------------------------------------------------------------------
      template<typename bodies_type>
      void reset_bodies_grav(const bodies_type*const&B) const {
	reset_bodies_pot(B);
	reset_bodies_acc(B);
      }
      //------------------------------------------------------------------------
      // simple manipulations                                                   
      //------------------------------------------------------------------------
      void set_sink  (sink_data*const&sink) { PROP = static_cast<void*>(sink); }
      void reset_sink()                     { sink()->reset(); }
      //------------------------------------------------------------------------
      void normalize_grav () {                     // acc,pot     /= mass       
	if(mass()>zero) {
	  register real im = one/mass();
	  pot() *= im;
	  acc() *= im;
	}
      }
      //------------------------------------------------------------------------
      // dump leaf data                                                         
      //------------------------------------------------------------------------
      static void dump_head(std::ostream& o) {
	OctTree::Leaf::dump_head(o);
	o<<"              mass";
      }
      //------------------------------------------------------------------------
      void dump(std::ostream&o) const {
	OctTree::Leaf::dump(o);
	o<<' '<<setw(8)<<mass();
      }
      //------------------------------------------------------------------------
    };// class Leaf {
    //--------------------------------------------------------------------------
    //                                                                          
    // sub-type Cell                                                            
    //                                                                          
    //--------------------------------------------------------------------------
    //                                                                          
    // On the data usage for data beyond those used by OctTree etc              
    //                                                                          
    // variable            |  datum                                             
    // --------------------+-------------------------------------------------   
    // OctTree::Cell::POS  |  centre of mass                                    
    // OctTree::Cell::RAD  |  rmax, rcrit, size                                 
    // OctTree::Cell::AUX1 |  pointer to srce_data                              
    // OctTree::Cell::AUX2 |  pointer to Coeffs (sink data)                     
    // OctTree::Cell::AUX3 |  eps/2 (only used with indiv softening lengths)    
    // srce_data::MASS     |  mass                                              
    // srce_data::POLS     |  specific multipole moments                        
    //                                                                          
    //--------------------------------------------------------------------------
    class Cell : public OctTree::Cell {
      Cell           (const Cell&);                // not implemented           
      Cell& operator=(const Cell&);                // not implemented           
      //------------------------------------------------------------------------
      // types and static data                                                  
      //------------------------------------------------------------------------
    public:
      typedef Leaf leaf_type;                      // type of associated leafs  
      typedef grav::Mset Mset;
      typedef grav::Cset Cset;
      //------------------------------------------------------------------------
      static const int N_SINK = Cset::NDAT;
      //------------------------------------------------------------------------
    private:
#if defined(__GNUC__) && (__GNUC__ < 3 || __GNUC__ == 3 && __GNUC_MINOR__ < 4)
      // patch for fix a bug with gcc version < 3.4
    public:
#endif
      struct srce_data {
	real MASS;
	Mset POLS;
	void normalize_poles()         { POLS.normalize(MASS); }
      };
      //------------------------------------------------------------------------
      // friendships                                                            
      //------------------------------------------------------------------------
      friend class GravEstimator;                  // for alloc of srce_data    
      friend class falcON::traits<srce_data>;
      //------------------------------------------------------------------------
      // private data access                                                    
      //------------------------------------------------------------------------
#define __srce  static_cast<srce_data*>(AUX1.PTER)
      real const&mass  () const { return __srce->MASS; }
      vect const&cofm  () const { return POS; }
      real const&rmax  () const { return RAD; }
      real const&size  () const { return RAD; }
      real const&rcrit () const { return RAD; }
      real       rcrit2() const { return square(rcrit()); }
      real const&eph   () const { return AUX3.SCAL; }
      Cset&Coeffs() const { return *static_cast<Cset*>(AUX2.PTER); }
      const Mset&poles () const { return __srce->POLS; }
      //------------------------------------------------------------------------
      // simple manipulations                                                   
      //------------------------------------------------------------------------
    public:
      void set_rcrit(real const&it) { RAD *= it; }
      //------------------------------------------------------------------------
      void set_srce(srce_data*const&srce) {
	AUX1.PTER = static_cast<void*>(srce);
      }
      //------------------------------------------------------------------------
      void setCoeffs(Cset*const&sink) {
	AUX2.PTER = static_cast<void*>(sink); 
      }
      //------------------------------------------------------------------------
      void*returnCoeffs()       { return AUX2.PTER; }
      void resetCoeffs ()       { AUX2.PTER=0; }
      bool hasCoeffs   () const { return AUX2.PTER != 0; }
      //------------------------------------------------------------------------
      // non-const data access via members                                      
      //------------------------------------------------------------------------
      real &mass  () { return __srce->MASS; }
      vect &cofm  () { return POS; }
      real &rmax  () { return RAD; }
      real &size  () { return RAD; }
      real &eph   () { return AUX3.SCAL; }
      Cset &Coeffs() { return *static_cast<Cset*>(AUX2.PTER); }
      Mset &poles () { return __srce->POLS; }
#undef __srce
      //------------------------------------------------------------------------
      // const data access via friends                                          
      //------------------------------------------------------------------------
      friend real const&mass  (const Cell*);
      friend vect const&cofm  (const Cell*);
      friend vect const&pos   (const Cell*);
      friend real const&rmax  (const Cell*);
      friend real const&rcrit (const Cell*);
      friend real const&size  (const Cell*);
      friend real       rcrit2(const Cell*);
      friend real const&eph   (const Cell*);
      friend Cset      &Coeffs(const Cell*);
      friend Mset const&poles (const Cell*);
      friend bool hasCoeffs   (const Cell*);
      //------------------------------------------------------------------------
      // boolean information via friends                                        
      //------------------------------------------------------------------------
      friend bool is_source (const Cell*C);
      //------------------------------------------------------------------------
      // other const methods and friends                                        
      //------------------------------------------------------------------------
      friend real xmin(const Cell*);
      friend real xmax(const Cell*);
      //------------------------------------------------------------------------
      // dump cell data                                                         
      //------------------------------------------------------------------------
      static void dump_head(std::ostream&o) {
	OctTree::Cell::dump_head(o);
	o<<
	  "              mass"
	  "              cofm         "
	  "         rmax"
	  "        rcrit";
      }
      //------------------------------------------------------------------------
      void dump(std::ostream&o) const {
	OctTree::Cell::dump(o);
	o<<' '<<setw(8)<<mass();
	for(register int i=0; i<Ndim; i++)
	  o<<' '<<setw(8)<<setprecision(4)<<cofm()[i];
	o<<' '<<setw(12)<<rmax()
	 <<' '<<setw(12)<<rcrit();
      }
      //------------------------------------------------------------------------
    };// class Cell {
  private:
    //--------------------------------------------------------------------------
    // data:                                                                    
    //--------------------------------------------------------------------------
    const OctTree        *TREE;                    // the tree to be used       
    bool                  CELLS_UPTODATE;          // are cell srces up to date?
    bool                  LEAFS_UPTODATE;          // are leaf srces up to date?
#ifdef falcON_INDI
    const bool            INDI_SOFT;               // use individual eps_i      
#endif
    const int             DIR[4];                  // direct loop control       
    kern_type             KERNEL;                  // softening kernel          
    GravStats            *STATS;                   // interaction statistics    
    real                  EPS;                     // global softening length   
    real                  GRAV;                    // Newton's G (can be 0)     
    unsigned              Ncoeffs,Nchunks,Ncsize;  // # coeffs, chunks, bytes/c 
    Cell::srce_data      *CELL_SRCE;               // memory for cell srce      
    Leaf::sink_data      *LEAF_SINK;               // memory for leafs          
    unsigned              NCT, NCA, NLA;           // # allocation of these     
    unsigned              NLA_needed;              // # active leafs            
    //--------------------------------------------------------------------------
    // private methods                                                          
    //--------------------------------------------------------------------------
    // update leafs' source fields (requires no allocation here)                
    // - mass, leaf flags[, eps/2]                                              
    // - count active leafs                                                     
    void update_leafs();
    //--------------------------------------------------------------------------
    // adjust leafs' eph_i                                                      
#ifdef falcON_ADAP
    void adjust_eph(bool     const&,               // I: for all or active only 
		    real     const&,               // I: Nsoft: adjust actives  
		    real     const&,               // I: emin:  adjust actives  
		    real     const&,               // I: emax:  adjust actives  
		    unsigned const&,               // I: Nref:  adjust actives  
		    real     const&);              // I: max change in eps_i    
#endif
    //--------------------------------------------------------------------------
    // - passes up the tree: flag, mass, cofm, rmax[, eph], multipoles          
    // - set the cells r_crit & r_crit^2                                        
    unsigned pass_up(                              // R: # active cells         
		     const GravMAC*const&,         // I: MAC                    
		     bool          const&);        // I: reused old tree?       
    //--------------------------------------------------------------------------
    // prepare for interactions                                                 
    //  - allocate memory for leaf sink properties for active leafs             
    //  - allocate memory for cell source properties & reset their coeff pter   
    //  - pass source properties up the tree, count active cells                
    bool prepare(                                  // R: ALL || all are active  
		 const GravMAC*const&,             // I: MAC                    
		 bool          const&,             // I: all or active only?    
		 bool          const&);            // I: allocate cell coeffs?  
    //--------------------------------------------------------------------------
    // tree stuff to be superseeded                                             
    //--------------------------------------------------------------------------
  public:
    typedef Cell                         cell_type;
    typedef Leaf                         leaf_type;
    typedef OctTree::CellIter<cell_type> cell_iterator;
    typedef leaf_type*                   leaf_iterator;
    //--------------------------------------------------------------------------
    const OctTree*const&my_tree() const { return TREE; }
    cell_iterator root         () const {
      return cell_iterator(TREE,static_cast<Cell*>(TREE->FstCell())); }
    //--------------------------------------------------------------------------
    // dump cell and leaf data                                                  
    //--------------------------------------------------------------------------
    void dump_cells(std::ostream&) const;
    void dump_leafs(std::ostream&) const;
    //--------------------------------------------------------------------------
    // construction: allocate memory for leafs' & cells' source properties      
    //--------------------------------------------------------------------------
    GravEstimator(const OctTree* const&T,          // I: tree to be used        
		  kern_type      const&k,          // I: kernel to be used      
		  GravStats*     const&st,         // I: statistics             
		  real           const&e,          // I: global/max eps         
		  real           const&g = one,    //[I: Newton's G]            
		  bool           const&s = 0,      //[I: use individual eps?]   
		  const int d[4]=Default::direct): //[I: N_direct for gravity]  
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
      LEAF_SINK      ( 0 ),
      NCT            ( 0u ),
      NCA            ( 0u ),
      NLA            ( 0u ),
      NLA_needed     ( 0u ),
      DIR            ()
    {
      const_cast<int*>(DIR)[0]=d[0];
      const_cast<int*>(DIR)[1]=d[1];
      const_cast<int*>(DIR)[2]=d[2];
      const_cast<int*>(DIR)[3]=d[3];
    }
    //--------------------------------------------------------------------------
    // reset allocations                                                        
    void reset() {
      LEAFS_UPTODATE = 0;
      CELLS_UPTODATE = 0;
    }
    //--------------------------------------------------------------------------
    // change tree entry, set flags                                             
    void new_tree(const OctTree*const&T) {
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
    ~GravEstimator() {
      if(CELL_SRCE) delete[] CELL_SRCE;
      if(LEAF_SINK) delete[] LEAF_SINK;
    }
    //--------------------------------------------------------------------------
    // compute gravity by direct summation                                      
    //--------------------------------------------------------------------------
    void exact(bool      const& =false             //[I: for all or active only]
#ifdef falcON_ADAP
	      ,real      const& =zero,             //[I: Nsoft: adjust eps_i]   
	       unsigned  const& =0u,               //[I: Nref:  adjust eps_i]   
	       real      const& =zero,             //[I: eps_min]               
	       real      const& =zero              //[I: max change of eps]     
#endif
	       );
    //--------------------------------------------------------------------------
    // compute gravity by the approximate method of Dehnen (2002)               
    //   if enabled, also adapt the individual softening lengths                
    //--------------------------------------------------------------------------
    void approx(const GravMAC*const&,              // I: MAC                    
		bool          const& =false,       //[I: for all or active only]
		bool          const& =true         //[I: combine phases]        
#ifdef falcON_ADAP
	       ,real          const& =zero,        //[I: Nsoft: adjust eps_i]   
		unsigned      const& =0u,          //[I: Nref:  adjust eps_i]   
		real          const& =zero,        //[I: eps_min]               
		real          const& =zero         //[I: max change of eps]     
#endif
		);
    //--------------------------------------------------------------------------
    // density estimation (number-, surface-, mass-density)                     
    //--------------------------------------------------------------------------
    void estimate_nd(bool     const&,              // I: all or active only?    
		     unsigned const&) const;       // I: critical cell size     
    void estimate_sd(bool     const&,              // I: all or active only?    
		     unsigned const&);             // I: critical cell size     
    void estimate_md(bool     const&,              // I: all or active only?    
		     unsigned const&);             // I: critical cell size     
    //--------------------------------------------------------------------------
    // public const data access                                                 
    //--------------------------------------------------------------------------
    const OctTree *const&tree            () const { return TREE; }
    kern_type      const&kernel          () const { return KERNEL; }
    GravStats     *const&stats           () const { return STATS; }
    real           const&NewtonsG        () const { return GRAV; }
    bool                 ZeroGravity     () const { return GRAV==zero; }
    real           const&softening_length() const { return EPS; }
    unsigned       const&N_active_cells  () const { return NCA; }
    unsigned       const&N_active_leafs  () const { return NLA; }
    unsigned       const&N_coeffs        () const { return Ncoeffs; }
    unsigned       const&N_chunks        () const { return Nchunks; }
    unsigned       const&N_elems_in_chunk() const { return Ncsize; }
#ifdef falcON_INDI
    bool           const&use_indiv_eps   () const { return INDI_SOFT; }
#endif
  };// class GravEstimator {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // inline definitions of friends of class GravEstimator::Leaf               //
  // also serve to inject these functions into namespace falcON               //
  //                                                                          //
  // ///////////////////////////////////////////////////////////////////////////
  inline vect const&cofm(const GravEstimator::Leaf*L) { return L->pos(); }
  inline vect const&acc(const GravEstimator::Leaf*L) { return L->acc(); }
  inline real const&pot(const GravEstimator::Leaf*L) { return L->pot(); }
  inline real const&rho(const GravEstimator::Leaf*L) { return L->rho(); }
  inline unsigned const&num(const GravEstimator::Leaf*L) { return L->num(); }
#ifdef falcON_INDI
  inline real const&eph(const GravEstimator::Leaf*L) { return L->eph(); }
  inline real size(const GravEstimator::Leaf*L) {return twice(L->eph()); }
  inline real const&sizeq(const GravEstimator::Leaf*L) { return L->sizeq(); }
#endif
  inline real const&mass(const GravEstimator::Leaf*L) { return L->mass(); }
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // inline definitions of friends of class GravEstimator::Cell               //
  // also serve to inject these functions into namespace falcON               //
  //                                                                          //
  // ///////////////////////////////////////////////////////////////////////////
  inline real const&mass(const GravEstimator::Cell*C) {
    return C->mass();
  }
  inline vect const&cofm(const GravEstimator::Cell*C) {
    return C->cofm();
  }
  inline vect const&pos(const GravEstimator::Cell*C) {
    return C->cofm();
  }
  inline real const&rmax(const GravEstimator::Cell*C) {
    return C->rmax();
  }
  inline real const&rcrit(const GravEstimator::Cell*C) {
    return C->rcrit();
  }
  inline real const&size(const GravEstimator::Cell*C) {
    return C->size();
  }
  inline real rcrit2(const GravEstimator::Cell*C) {
    return C->rcrit2();
  }
  inline real const&eph(const GravEstimator::Cell*C) {
    return C->eph();
  }
  inline grav::Cset&Coeffs(const GravEstimator::Cell*C) {
    return C->Coeffs();
  }
  inline grav::Mset const&poles(const GravEstimator::Cell*C) {
    return C->poles();
  }
  inline bool hasCoeffs(const GravEstimator::Cell*C) {
    return C->hasCoeffs();
  }
  inline bool is_source(const GravEstimator::Cell*C) { 
    return mass(C)!=zero;
  }
  inline real xmin(const GravEstimator::Cell*C) {
    return cofm(C).min() - rmax(C);
  }
  inline real xmax(const GravEstimator::Cell*C) {
    return cofm(C).max() + rmax(C);
  }
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // namespace falcON::grav                                                   //
  //                                                                          //
  // ///////////////////////////////////////////////////////////////////////////
  namespace grav {
    typedef GravEstimator::Cell           cell;
    typedef GravEstimator::Leaf           leaf;
    typedef GravEstimator::Cell          *cell_pter;
    typedef GravEstimator::Leaf          *leaf_pter;
    typedef GravEstimator::cell_iterator  cell_iter;
    typedef GravEstimator::leaf_iterator  leaf_iter;
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class GravStats                                                          //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class GravStats {
  private:
    typedef grav::cell_iter cell_iter;
    typedef grav::leaf_iter leaf_iter;
    unsigned D_BB, D_CB, D_CC, D_CX;               // # direct interactions     
    unsigned A_CB, A_CC;                           // # approximate interactions
#ifdef ENHANCED_IACT_STATS
    unsigned P_CB, P_CC, P_CX;                     // # BB pairs in direct      
#  define ADD_SS(COUNTER,NUMBER)  COUNTER += NUMBER;
    static const char* trick_SS(unsigned const&n, int&w)
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
    unsigned const&BB_direct_iacts   () const { return D_BB; }
    unsigned const&CB_direct_iacts   () const { return D_CB; }
    unsigned const&CC_direct_iacts   () const { return D_CC; }
    unsigned const&CX_direct_iacts   () const { return D_CX; }
    unsigned       total_direct_iacts() const {
      return
	BB_direct_iacts() +
	CB_direct_iacts() +
	CC_direct_iacts() +
	CX_direct_iacts();
    }
#ifdef ENHANCED_IACT_STATS
    unsigned const&BB_direct_pairs   () const { return D_BB; }
    unsigned const&CB_direct_pairs   () const { return P_CB; }
    unsigned const&CC_direct_pairs   () const { return P_CC; }
    unsigned const&CX_direct_pairs   () const { return P_CX; }
    unsigned       total_direct_pairs() const {
      return
	BB_direct_pairs() +
	CB_direct_pairs() +
	CC_direct_pairs() +
	CX_direct_pairs();
    }
#endif
    unsigned       BB_approx_iacts   () const { return 0u; }
    unsigned const&CB_approx_iacts   () const { return A_CB; }
    unsigned const&CC_approx_iacts   () const { return A_CC; }
    unsigned       CX_approx_iacts   () const { return 0u; }
    unsigned       total_approx_iacts() const {
      return
	BB_approx_iacts() +
	CB_approx_iacts() +
	CC_approx_iacts() +
	CX_approx_iacts();
    }
    // 3 writing stats to ostream                                               
    void write(std::ostream&out) const {
      unsigned part, total=total_approx_iacts()+total_direct_iacts();
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
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
falcON_TRAITS(falcON::grav::Cset,"grav::Cset","grav::Csets");
falcON_TRAITS(falcON::grav::Mset,"grav::Mset","grav::Msets");
falcON_TRAITS(falcON::GravEstimator,"GravEstimator","GravEstimators");
falcON_TRAITS(falcON::GravEstimator::Cell,
	      "GravEstimator::Cell","GravEstimator::Cells");
falcON_TRAITS(falcON::GravEstimator::Leaf,
	      "GravEstimator::Leaf","GravEstimator::Leafs");
falcON_TRAITS(falcON::GravEstimator::Cell::srce_data,
	      "GravEstimator::Cell::srce_data",
	      "GravEstimator::Cell::srce_data");
falcON_TRAITS(falcON::GravEstimator::Leaf::sink_data,
	      "GravEstimator::Leaf::sink_data",
	      "GravEstimator::Leaf::sink_data");
////////////////////////////////////////////////////////////////////////////////
#undef  ENHANCED_IACT_STATS
#endif // falcON_included_gravity_h
    
