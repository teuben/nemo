// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// grat.h                                                                      |
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
// class grav_mac                                                              |
// class grav_soul                                                             |
// class grav_cell                                                             |
// class grav_tree                                                             |
// class grav_stat                                                             |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_grat_h
#define falcON_included_grat_h

#ifndef falcON_included_tree_h
#  include <public/tree.h>
#endif
#ifndef falcON_included_deft_h
#  include <public/deft.h>
#endif
#ifndef falcON_included_ionl_h
#  include <public/ionl.h>
#endif
#ifndef falcON_included_memo_h
#  include <public/memo.h>
#endif

////////////////////////////////////////////////////////////////////////////////
#define ENHANCED_IACT_STATS
#undef  ENHANCED_IACT_STATS
////////////////////////////////////////////////////////////////////////////////
#define falcON_ORDER 3                             // expansion order is fixed !
////////////////////////////////////////////////////////////////////////////////
namespace nbdy {
  //////////////////////////////////////////////////////////////////////////////
  class InvertZ;                                   // forward declaration       
  class grav_tree;                                 // forward declaration       
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
	       const uint = 3);                    //[I: expansion order]       
    //--------------------------------------------------------------------------
    void reset(const MAC_type,                     // I: type of MAC            
	       const real,                         // I: parameter: theta_0     
	       const uint = 3);                    //[I: expansion order]       
    //--------------------------------------------------------------------------
    void reset_theta(const real);                  // I: parameter: theta_0     
    //--------------------------------------------------------------------------
    ~grav_mac  ();                                 // destructor                
    //--------------------------------------------------------------------------
    // const method                                                             
    //--------------------------------------------------------------------------
    void set_rcrit(const grav_tree*) const;        // set rcrit for all cells   
    //--------------------------------------------------------------------------
    // const inlined methods                                                    
    //--------------------------------------------------------------------------
    const real  theta_min() const { return TH0; }  // theta_0                   
    const MAC_type&method() const { return MAC; }  // MAC                       
    const char* describe_method() const            // describes MAC             
    { return describe(MAC); }
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::grav_soul                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class grav_soul : public basic_soul {
    grav_soul           (const grav_soul&);        // not implemented           
    grav_soul& operator=(const grav_soul&);        // not implemented           
    //--------------------------------------------------------------------------
    // friends                                                                  
    //--------------------------------------------------------------------------
    friend class grav_tree;
    //--------------------------------------------------------------------------
    // data of class grav_soul                                                  
    //--------------------------------------------------------------------------
  private:
    real    MASS;                                  // mass of body              
    real   *SINKPT;                                // sink part of soul         
#ifdef falcON_INDI
    real    EPH;                                   // epsi/2                    
#endif
  public:
    //--------------------------------------------------------------------------
    // stuff needed for MPI code                                                
    //--------------------------------------------------------------------------
#ifdef falcON_MPI
    void copy_prune(const grav_soul*S) {
      basic_soul::copy_prune(S);
      MASS = nbdy::mass(S);
#ifdef falcON_INDI
      EPH  = nbdy::eph (S);
#endif
    }
#endif
    //--------------------------------------------------------------------------
    // non-const data access via members                                        
    //--------------------------------------------------------------------------
    ten1 acc  ()            { return ten1(SINKPT+1); }
    real&acc  (int const&i) { return SINKPT[i+1]; }
    vect&cofm ()            { return pos(); }
    real&cofm (int const&i) { return pos(i); }
    real&sizeq()            { return SINKPT[1]; }
    real&pot  ()            { return SINKPT[0]; }
    real&rho  ()            { return SINKPT[0]; }
    uint&num  ()            { return static_cast<uint*>(
				     static_cast<void*>(SINKPT))[0]; }
#ifdef falcON_INDI
    real&eph  ()            { return EPH; }
#endif
    void inc  ()            { ++num(); }
    //--------------------------------------------------------------------------
    // const data access via friends                                            
    //--------------------------------------------------------------------------
#define CGCS const grav_soul*const&S
    friend vect const&cofm  (CGCS) { return S->pos(); }
    friend vect const&pos   (CGCS) { return S->pos(); }
    friend uint const&mybody(CGCS) { return S->mybody(); }
    friend ten1 const acc   (CGCS) { return ten1(S->SINKPT+1); }
    friend real const&acc   (CGCS,
			     int i){ return S->SINKPT[i+1]; }
    friend real const&pot   (CGCS) { return S->SINKPT[0]; }
    friend real const&rho   (CGCS) { return S->SINKPT[0]; }
    friend uint const&num   (CGCS) { return static_cast<uint*>(
				 	    static_cast<void*>(S->SINKPT))[0]; }
#ifdef falcON_INDI
    friend real const&eph   (CGCS) { return S->EPH; }
    friend real const size  (CGCS) { return twice(S->EPH); }
    friend real const&sizeq (CGCS) { return S->SINKPT[1]; }
#endif
    friend real const&mass  (CGCS) { return S->MASS; }
#undef CGCS
    //--------------------------------------------------------------------------
    // simple manipulations                                                     
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void set_mass(const bodies_type*const&B) {
      MASS = B->mass(mybody());
    }
    //--------------------------------------------------------------------------
    void set_mass_and_pos(const barrays*const&B) {
      MASS      = B->mass (mybody());
      cofm()[0] = B->pos_x(mybody());
      cofm()[1] = B->pos_y(mybody());
#if falcON_NDIM==3
      cofm()[2] = B->pos_z(mybody());
#endif
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    template<typename bodies_type>
    void set_mass_and_pos(const bodies_type*const&B) {
      MASS   = B->mass(mybody());
      cofm() = B->pos (mybody());
    }
#ifdef falcON_INDI
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void copy_eph(const bodies_type*const&B) {
      eph() = half*B->eps(mybody());
    }
#endif
    //--------------------------------------------------------------------------
    void reset_srce() {
      SINKPT[0] = zero;
      SINKPT[1] = zero;
      SINKPT[2] = zero;
#if falcON_NDIM > 2
      SINKPT[3] = zero;
#endif
    }
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void update_dens(const bodies_type*const&B) const {
      B->rho(mybody()) = nbdy::rho(this);
    }
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void update_grav(const bodies_type*const&B) const {
      B->pot(mybody()) = nbdy::pot(this);
      B->acc(mybody()) = nbdy::acc(this);
    }
    //--------------------------------------------------------------------------
    void update_grav(const barrays*const&B) const {
      if(B->has(io::p)) B->pot(mybody()) = nbdy::pot(this);
      B->acc_x(mybody()) = nbdy::acc(this,0);
      B->acc_y(mybody()) = nbdy::acc(this,1);
#if falcON_NDIM>2
      B->acc_z(mybody()) = nbdy::acc(this,2);
#endif
    }
    //--------------------------------------------------------------------------
#ifdef falcON_INDI
    template<typename bodies_type>
    void update_eps(const bodies_type*const&B) const {
      B->eps(mybody()) = twice(nbdy::eph(this));
    }
#endif
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void update_num(const bodies_type*const&B) const {
      B->num(mybody()) = nbdy::num(this);
    }
    //--------------------------------------------------------------------------
    void normalize_grav () {                       // acc,pot     /= mass       
      if(MASS>zero) {
	register real im = one/MASS;
	pot() *= im;
	acc() *= im;
      }
    }
    //--------------------------------------------------------------------------
    // boolean information                                                      
    //--------------------------------------------------------------------------
    friend bool is_source(const grav_soul*const&S) { return mass(S) != zero; }
    //--------------------------------------------------------------------------
    // dump soul data                                                           
    //--------------------------------------------------------------------------
    static void dump_head(std::ostream& o) {
      o<<"     #";
      basic_soul::dump_head(o);
      o<<"              mass";
    }
    //--------------------------------------------------------------------------
    static void dump(std::ostream&o, const grav_soul*const&S, const int&index) {
      o<<' '<<setw(5)<<index;
      basic_soul::dump(o,S);
      o<<' '<<setw(8)<<nbdy::mass(S);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // nbdy::class grav_cell                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class grav_cell : public basic_cell {
    grav_cell           (const grav_cell&);        // not implemented           
    grav_cell& operator=(const grav_cell&);        // not implemented           
    //--------------------------------------------------------------------------
    // types                                                                    
    //--------------------------------------------------------------------------
  public:
    typedef grav_soul soul_type;                   // type of associated souls  
    //--------------------------------------------------------------------------
    // friends                                                                  
    //--------------------------------------------------------------------------
    friend class grav_tree;
    friend class grap_cell;
    //--------------------------------------------------------------------------
    // static data and methods                                                  
    //--------------------------------------------------------------------------
  private:
    static const int
      N_MASS  = 0,
      N_IMAS  = 1,
      N_RMAX  = 2,
      N_RCRT  = 3,
      N_RCR2  = 4,
      N_COFM  = 5,
      N_P2    = N_COFM + ten1::NDAT,
#if   falcON_ORDER > 3
      N_P3    = N_P2   + ten2::NDAT,
# if  falcON_ORDER > 4
      N_P4    = N_P3   + ten3::NDAT,
#  if falcON_ORDER > 5
#   error "expansion order > 5 not supported in public/grav.h"
#  endif
      N_EPH   = N_P4   + ten4::NDAT,
# else
      N_EPH   = N_P3   + ten3::NDAT,
# endif
#else
      N_EPH   = N_P2   + ten2::NDAT,
#endif
#ifdef falcON_INDI
      N_TOTAL = N_EPH  + 1;
    static const int N_eph() { return N_EPH; }
#else
      N_TOTAL = N_EPH;
#endif
    static const int N_tot() { return N_TOTAL; }
    //--------------------------------------------------------------------------
    static const int N_COEFF = 1 + ten1::NDAT + ten2::NDAT + ten3::NDAT
#if  falcON_ORDER > 3
    + ten4::NDAT
#if  falcON_ORDER > 4
    + ten5::NDAT
#if  falcON_ORDER > 5
    + ten6::NDAT
#endif
#endif
#endif
    ;
    static int const&N_coeff() { return N_COEFF; }
    //--------------------------------------------------------------------------
    // data of class grav_cell                                                  
    //--------------------------------------------------------------------------
    real   *SOURCE;                                // source data for cell      
    real   *COEFFS;                                // sink data for cell        
  public:
    //--------------------------------------------------------------------------
    // simple manipulations                                                     
    //--------------------------------------------------------------------------
    void set_rcrit(real const&ith) {
      SOURCE[N_RCRT] = ith * SOURCE[N_RMAX];
      SOURCE[N_RCR2] = square(SOURCE[N_RCRT]);
    }
    //--------------------------------------------------------------------------
    // non-const data access via members                                        
    //--------------------------------------------------------------------------
    ten1   cofm  () { return ten1(SOURCE+N_COFM); }
    real  &mass  () { return SOURCE[N_MASS]; }
    real  &imass () { return SOURCE[N_IMAS]; }
    real  &rmax  () { return SOURCE[N_RMAX]; }
    real  &size  () { return SOURCE[N_RCRT]; }
#ifdef falcON_INDI
    real  &eph   () { return SOURCE[N_EPH]; }
#endif
    real* &coeffs() { return COEFFS; }
    ten2   quad  () { return ten2(SOURCE+N_P2); }
#if   falcON_ORDER > 3
    ten3   octo  () { return ten3(SOURCE+N_P3); }
# if  falcON_ORDER > 4
    ten4   hexa  () { return ten4(SOURCE+N_P4); }
# endif
#endif
    //--------------------------------------------------------------------------
    // const data access via friends                                            
    //--------------------------------------------------------------------------
#define CCCC const grav_cell*const&C
    friend real const&radius(CCCC) { return C->radius(); }
    friend vect const&center(CCCC) { return C->center(); }
    friend int  const&number(CCCC) { return C->number(); }
    friend ten1 const cofm  (CCCC) { return ten1(C->SOURCE+N_COFM); }
    friend ten1 const pos   (CCCC) { return ten1(C->SOURCE+N_COFM); }
    friend real const&mass  (CCCC) { return C->SOURCE[N_MASS]; }
    friend real const&imass (CCCC) { return C->SOURCE[N_IMAS]; }
    friend real const&rmax  (CCCC) { return C->SOURCE[N_RMAX]; }
    friend real const&rcrit (CCCC) { return C->SOURCE[N_RCRT]; }
    friend real const&size  (CCCC) { return C->SOURCE[N_RCRT]; }
    friend real const&rcrit2(CCCC) { return C->SOURCE[N_RCR2]; }
#ifdef falcON_INDI
    friend real const&eph   (CCCC) { return C->SOURCE[N_EPH]; }
#endif
    friend real*const&coeffs(CCCC) { return C->COEFFS; }
    friend ten2       quad  (CCCC) { return ten2(C->SOURCE+N_P2); }
#if   falcON_ORDER > 3
    friend ten3       octo  (CCCC) { return ten3(C->SOURCE+N_P3); }
# if  falcON_ORDER > 4
    friend ten4       hexa  (CCCC) { return ten4(C->SOURCE+N_P4); }
# endif
#endif
    //--------------------------------------------------------------------------
    // boolean information via friends                                          
    //--------------------------------------------------------------------------
    friend bool is_source    (CCCC) { return nbdy::mass(C)!=zero; }
    friend bool has_cell_kids(CCCC) { return C->has_cell_kids(); }
    friend bool has_soul_kids(CCCC) { return C->has_soul_kids(); }
    friend bool is_twig      (CCCC) { return C->is_twig(); }
    friend bool is_branch    (CCCC) { return C->is_branch(); }
    //--------------------------------------------------------------------------
    // other const methods and friends                                          
    //--------------------------------------------------------------------------
    friend real xmin(CCCC) { return nbdy::cofm(C).min() - nbdy::rmax(C); }
    friend real xmax(CCCC) { return nbdy::cofm(C).max() + nbdy::rmax(C); }
#undef CCCC
    //--------------------------------------------------------------------------
    // flag manipulation                                                        
    //--------------------------------------------------------------------------
    void add_active_flag(const grav_soul* const&S) {
      add_active_flag_from_soul(S);
    }
    void add_active_flag(const grav_cell* const&C) {
      add_active_flag_from_cell(C);
    }
    //--------------------------------------------------------------------------
    // dump cell data                                                           
    //--------------------------------------------------------------------------
    static void dump_head(std::ostream&o) {
      o<<"     #";
      basic_cell::dump_head(o);
      o<<"              mass              cofm                  rmax";
    }
    //--------------------------------------------------------------------------
    static void dump(std::ostream&o, const grav_cell*const&C, const int&index) {
      o<<' '<<setw(5)<<index;
      basic_cell::dump(o,C);
      o<<' '<<setw(8)<<nbdy::mass(C);
      for(register int i=0; i<Ndim; i++)
	o<<" "<<setw(8)<<setprecision(4)<<nbdy::cofm(C)[i];
      o<<" "<<setw(12)<<nbdy::rmax(C);
    }
    //--------------------------------------------------------------------------
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class nbdy::grav_tree                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class grav_stat;                                 // forward declaration       
  class grav_tree : public basic_tree<grav_tree,grav_cell> {
    //--------------------------------------------------------------------------
    friend class tree_transport<grav_tree>;        // for parallel code         
    //--------------------------------------------------------------------------
    grav_tree            (grav_tree const&);       // not implemented           
    grav_tree& operator= (grav_tree const&);       // not implemented           
    //--------------------------------------------------------------------------
    // data:                                                                    
    //--------------------------------------------------------------------------
#ifdef falcON_INDI
    const soft_type   SOFT;                        // global / individual       
#endif
    const int         nCS;                         // source:    memory / cell  
    static const int  nSS = Ndim+1;                // sink part: memory / soul  
    static const int  nCC = grav_cell::N_COEFF;    // coeffs:    memory / cell  
    int               Ncs, Nss;                    // # active cell/soul        
    uint              Ncoeffs;                     // # coeffs used             
    real             *CELL_SOURCE;                 // memory for cell source    
    real             *CELL_COEFFS;                 // memory for cell sink part 
    real             *SOUL_SINKPT;                 // memory for soul sink part 
    //--------------------------------------------------------------------------
    inline void reset_soul_sinkpt(bool const& =0); // set grav_soul::SINKPT     
    inline void reset_cell_source();               // set grav_cell::SOURCE     
    inline void reset_cell_coeffs(bool const& =0); // for parallel code only    
    inline void   set_soul_sinkpt(bool const& =0); // set grav_soul::SINKPT     
    inline void   set_cell_source();               // set grav_cell::SOURCE     
    inline void   set_cell_coeffs(bool const& =0); // for parallel code only    
    //--------------------------------------------------------------------------
    // prepare for neighbour counting:                                          
    // - update souls' and cells' active flags                                  
    // - give sinkpt memory to active souls                                     
    // - update souls' eph_i                       (individual softening only)  
    // - souls: set size^2=eph^2, num=0                                         
    // - cell actives: update size (pass up tree)                               
    // - update cell's eph_i (pass up the tree)    (individual softening only)  
    void prepare_count_neighbours(                 //                           
				  real const& =0,  //[I: global body size]      
				  bool const& =0); //[I: re-use memory?]        
    //--------------------------------------------------------------------------
    // prepare for exact gravity computation:                                   
    // - update souls' and cells' active flags                                  
    // - update souls' eph_i                       (individual softening only)  
    // - optionally: adjust active's eph_i         (individual softening only)  
    // - update cell's eph_i (pass up the tree)    (individual softening only)  
    // - give sinkpt memory to active souls                                     
    // - reset souls' pot & acc                                                 
    void prepare_grav_exact (                      // prepare for direct sums   
			     bool const&,          // I: for all or active only 
#ifdef falcON_INDI
			     real const&,          // I: Nsoft: adjust actives  
			     real const&,          // I: emin:  adjust actives  
			     real const&,          // I: emax:  adjust actives  
			     uint const&,          // I: Nref:  adjust actives  
			     real const&,          // I: max change in eps_i    
#endif
			     bool const&);         // I: re-use memory?         
    //--------------------------------------------------------------------------
    // prepare for gravity approximation:                                       
    // - update souls' and cells' active flags                                  
    // - update souls' eph_i                       (individual softening only)  
    // - optionally: adjust active's eph_i         (individual softening only)  
    // - update cell's eph_i (pass up the tree)    (individual softening only)  
    // - give sinkpt memory to active souls                                     
    // - reset souls' pot & acc                                                 
    // - optionally give coeffs memory to active cells                          
    // - recursively compute the multipoles                                     
    // - set the cells r_crit & r_crit^2                                        
    void prepare_grav_approx(const grav_mac*const&,// I: MAC                    
			     bool const&,          // I: for all or active only 
			     bool const&,          // I: give cell coeffs?      
#ifdef falcON_INDI
			     real const&,          // I: Nsoft: adjust actives  
			     real const&,          // I: emin:  adjust actives  
			     real const&,          // I: emax:  adjust actives  
			     uint const&,          // I: Nref:  adjust actives  
			     real const&,          // I: max change in eps_i    
#endif
			     bool const&);         // I: re-use memory?         
    //--------------------------------------------------------------------------
    // update gravity [& eps] of bodies from souls                              
    void update_grav_eps(bool const&,              //[I: update eps, too?]      
			 bool const&);             //[I: for all or active only]
  public:
    //--------------------------------------------------------------------------
    // construction:                                                            
    // - sets up the oct tree structure from scratch                            
    // - gives source memory to cells (souls have their own)                    
    // - recursively computes for every cell: mass, cofm, rmax                  
    // - does NOT give sink memory to cells nor to souls                        
    //--------------------------------------------------------------------------
    //   construction from list of bodies                                       
    grav_tree(const sbodies* const&,               // I: sbodies                
#ifdef falcON_INDI
	      const soft_type  = Default::soften,  //[I: global/individual]     
#endif
	      const int        = Default::Ncrit);  //[I: N_crit]                
    //--------------------------------------------------------------------------
#ifdef falcON_MPI
    grav_tree(const pbodies* const&,               // I: pbodies                
#ifdef falcON_INDI
	      const soft_type  = Default::soften,  //[I: global/individual]     
#endif
	      const int        = Default::Ncrit);  //[I: N_crit]                
#endif
    //--------------------------------------------------------------------------
    //   construction from list of bodies and boxing position x_min, x_max      
    grav_tree(const sbodies* const&,               // I: sbodies                
	      vect           const&,               // I: x_min                  
	      vect           const&,               // I: x_max                  
#ifdef falcON_INDI
	      const soft_type  = Default::soften,  //[I: global/individual]     
#endif
	      const int        = Default::Ncrit);  //[I: N_crit]                
    //--------------------------------------------------------------------------
#ifdef falcON_MPI
    grav_tree(const pbodies* const&,               // I: pbodies                
	      vect           const&,               // I: x_min                  
	      vect           const&,               // I: x_max                  
#ifdef falcON_INDI
	      const soft_type  = Default::soften,  //[I: global/individual]     
#endif
	      const int        = Default::Ncrit);  //[I: N_crit]                
#endif
    //--------------------------------------------------------------------------
    //   construction from barrays                                              
    grav_tree(const barrays* const&,               // I: body arrays            
#ifdef falcON_INDI
	      const soft_type  = Default::soften,  //[I: global/individual]     
#endif
	      const int        = Default::Ncrit);  //[I: N_crit]                
    //--------------------------------------------------------------------------
    // destruction                                                              
    ~grav_tree() {
      if(CELL_SOURCE) delete[] CELL_SOURCE;
      if(SOUL_SINKPT) delete[] SOUL_SINKPT;
      if(CELL_COEFFS) delete[] CELL_COEFFS;
    }
    //--------------------------------------------------------------------------
    // re_build:                                                                
    // - deletes old cell source memory                                         
    // - deletes old soul sinkpt memory                                         
    // - sets up the oct tree structure, whereby trying to use old structure    
    // - gives source memory to cells (souls have their own)                    
    // - recursively computes for every cell: mass, cofm, rmax                  
    // - does NOT give sink memory to cells nor to souls                        
    void rebuild(const int = Default::Ncrit,       //[I: N_crit]                
		 const int = 0);                   //[I: N_cut; 0: rebuild]     
    //--------------------------------------------------------------------------
    // re_use:                                                                  
    // - re-use old tree structure, do not establish new tree                   
    // - deletes old soul sinkpt memory                                         
    // - recursively computes for every cell: mass, cofm, rmax                  
    void reuse();
    //--------------------------------------------------------------------------
    // prepare for density estimation:                                          
    // - update souls active flags                                              
    // - pass active flags up the tree                                          
    // - give sinkpt memory to active souls                                     
    void prepare_density    (bool = false);        //[I: re-use memory?]        
    //--------------------------------------------------------------------------
    // compute gravity by direct summation                                      
    //--------------------------------------------------------------------------
    void exact_gravity(kern_type const&,           // I: kernel to be used      
		       grav_stat*const&,           // I: statistics             
		       real      const&,           // I: global/max eps         
		       bool      const& =false     //[I: for all or active only]
#ifdef falcON_INDI
		      ,real      const& =zero,     //[I: Nsoft: adjust eps_i]   
		       uint      const& =0u,       //[I: Nref:  adjust eps_i]   
		       real      const& =zero,     //[I: eps_min]               
		       real      const& =zero      //[I: max change of eps]     
#endif
		       );
    //--------------------------------------------------------------------------
    // compute gravity by the approximate method of Dehnen (2002)               
    //--------------------------------------------------------------------------
    void approx_gravity(const grav_mac*const&,     // I: MAC                    
			kern_type const&,          // I: kernel to be used      
			grav_stat*const&,          // I: statistics             
			real      const&,          // I: global/max eps         
			bool      const& =false,   //[I: for all or active only]
			bool      const& =true,    //[I: combine phases]        
#ifdef falcON_INDI
			real      const& =zero,    //[I: Nsoft: adjust eps_i]   
			uint      const& =0u,      //[I: Nref:  adjust eps_i]   
			real      const& =zero,    //[I: eps_min]               
			real      const& =zero,    //[I: max change of eps]     
#endif
			const int[4]=Default::direct); //[I: direct sum control]
    //--------------------------------------------------------------------------
    // count neighbours                                                         
    //--------------------------------------------------------------------------
    void count_neighbours(                         // count neighbours          
			  real const& = zero);     //[I: eps if global softning]
    //--------------------------------------------------------------------------
    // informations                                                             
    //--------------------------------------------------------------------------
    int         const& N_active_cells() const { return Ncs; }
    int         const& N_active_souls() const { return Nss; }
    uint        const& N_coeffs      () const { return Ncoeffs; }
  };
#if defined(__GNUC__) && (__GNUC__ < 3)
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // const data access to grav_cell from cell_iterator                        //
  //                                                                          //
  // these are not really necessary, since they are already defined for taking//
  // pointer to cell and cell_iterator has a type-conversion operator defined,//
  // but some older version of the GNU compiler gcc (gcc 2.95.3) get confused //
  // if we omit them.                                                         //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
# define CCCC grav_tree::cell_iterator const&I
  inline ten1 const cofm     (CCCC) { return cofm     (I.c_pter()); }
  inline ten1 const pos      (CCCC) { return pos      (I.c_pter()); }
  inline real const&mass     (CCCC) { return mass     (I.c_pter()); }
  inline real const&imass    (CCCC) { return imass    (I.c_pter()); }
  inline real const&rmax     (CCCC) { return rmax     (I.c_pter()); }
  inline real const&rcrit    (CCCC) { return rcrit    (I.c_pter()); }
  inline real const&size     (CCCC) { return size     (I.c_pter()); }
  inline real const&rcrit2   (CCCC) { return rcrit2   (I.c_pter()); }
#ifdef falcON_INDI
  inline real const&eph      (CCCC) { return eph      (I.c_pter()); }
#endif
  inline real*const&coeffs   (CCCC) { return coeffs   (I.c_pter()); }
  inline ten2       quad     (CCCC) { return quad     (I.c_pter()); }
# if   falcON_ORDER > 3
  inline ten3       octo     (CCCC) { return octo     (I.c_pter()); }
#  if  falcON_ORDER > 4
  inline ten4       hexa     (CCCC) { return hexa     (I.c_pter()); }
#  endif
# endif
  inline bool       is_source(CCCC) { return is_source(I.c_pter()); }
  inline real       xmin     (CCCC) { return xmin     (I.c_pter()); }
  inline real       xmax     (CCCC) { return xmax     (I.c_pter()); }
# undef CCCC
#endif
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // namespace nbdy::grav                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  namespace grav {
    //--------------------------------------------------------------------------
    // static data                                                              
    //--------------------------------------------------------------------------
    static const int
      N_C1    = 1,
      N_C2    = N_C1 + ten1::NDAT,
      N_C3    = N_C2 + ten2::NDAT,
#if   falcON_ORDER > 3
      N_C4    = N_C3 + ten3::NDAT,
# if  falcON_ORDER > 4
      N_C5    = N_C4 + ten4::NDAT,
#  if falcON_ORDER > 5
      N_C6    = N_C5 + ten5::NDAT,
      N_COEFF = N_C6 + ten6::NDAT;
#  else
      N_COEFF = N_C5 + ten5::NDAT;
#  endif
# else
      N_COEFF = N_C4 + ten4::NDAT;
# endif
#else
      N_COEFF = N_C3 + ten3::NDAT;
#endif
    //--------------------------------------------------------------------------
    // types                                                                    
    //--------------------------------------------------------------------------
    typedef grav_tree::cell_iterator  cell_iter;   // cell iterator             
    typedef grav_tree::soul_iterator  soul_iter;   // soul iterator             
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class grav_stat                                                          //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class grav_stat {
  private:
    typedef grav_tree::cell_iterator cell_iter;
    typedef grav_tree::soul_iterator soul_iter;
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
    void record_approx_CB(cell_iter const&A, soul_iter const&B)
    { ++A_CB; }
    void record_approx_CC(cell_iter const&A, cell_iter const&B)
    { ++A_CC; }
    void record_direct_CB(cell_iter const&A, soul_iter const&B)
    { ++D_CB; ADD_SS(P_CB, number(A)) }
    void record_direct_CC(cell_iter const&A, cell_iter const&B)
    { ++D_CC; ADD_SS(P_CC, number(A)*number(B)) }
    void record_direct_CX(cell_iter const&A)
    { ++D_CX; ADD_SS(P_CX, (number(A)*(number(A)-1))>>1 ) }
#ifdef ENHANCED_IACT_STATS
#  undef ADD_SS
#endif
    // 2 reporting                                                             
    const uint& BB_direct_iacts() const { return D_BB; }
    const uint& CB_direct_iacts() const { return D_CB; }
    const uint& CC_direct_iacts() const { return D_CC; }
    const uint& CX_direct_iacts() const { return D_CX; }
    const uint  total_direct_iacts() const { return
					       BB_direct_iacts() +
					       CB_direct_iacts() +
					       CC_direct_iacts() +
					       CX_direct_iacts();
    }
#ifdef ENHANCED_IACT_STATS
    const uint& BB_direct_pairs() const { return D_BB; }
    const uint& CB_direct_pairs() const { return P_CB; }
    const uint& CC_direct_pairs() const { return P_CC; }
    const uint& CX_direct_pairs() const { return P_CX; }
    const uint  total_direct_pairs() const { return
					       BB_direct_pairs() +
					       CB_direct_pairs() +
					       CC_direct_pairs() +
					       CX_direct_pairs();
    }
#endif
    const uint  BB_approx_iacts() const { return 0u; }
    const uint& CB_approx_iacts() const { return A_CB; }
    const uint& CC_approx_iacts() const { return A_CC; }
    const uint  CX_approx_iacts() const { return 0u; }
    const uint  total_approx_iacts() const { return
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
#endif                                             // falcON_included_grat_h    
