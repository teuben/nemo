// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// grat.h                                                                      |
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
// class grav_mac                                                              |
// class grav_soul                                                             |
// class grav_cell                                                             |
// class grav_tree                                                             |
// class grav_stat                                                             |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef included_grat_h
#define included_grat_h

#ifndef included_body_h
#  include <body.h>
#endif
#ifndef included_tree_h
#  include <public/tree.h>
#endif
#ifndef included_deft_h
#  include <public/deft.h>
#endif
#ifndef included_ionl_h
#  include <public/ionl.h>
#endif
#ifndef included_memo_h
#  include <public/memo.h>
#endif

////////////////////////////////////////////////////////////////////////////////
#define ENHANCED_IACT_STATS
#undef  ENHANCED_IACT_STATS
////////////////////////////////////////////////////////////////////////////////
#define P_ORDER 3                                  // expansion order is fixed !
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
#ifdef ALLOW_INDI
    real    EPH;                                   // epsi/2                    
#endif
    //--------------------------------------------------------------------------
    // constructor                                                              
    //--------------------------------------------------------------------------
  public:
    grav_soul() {}
    //--------------------------------------------------------------------------
    // non-const data access via members                                        
    //--------------------------------------------------------------------------
    ten1 acc  ()            { return ten1(SINKPT+1); }
    real&acc  (int const&i) { return SINKPT[i+1]; }
    real&sizeq()            { return SINKPT[1]; }
    real&pot  ()            { return SINKPT[0]; }
    real&rho  ()            { return SINKPT[0]; }
    uint&num  ()            { return static_cast<uint*>(
				     static_cast<void*>(SINKPT))[0]; }
#ifdef ALLOW_INDI
    real&eph  ()            { return EPH; }
#endif
    void inc  ()            { ++num(); }
    //--------------------------------------------------------------------------
    // const data access via friends                                            
    //--------------------------------------------------------------------------
#define CGCS const grav_soul*const&S
    friend vect const&cofm  (CGCS) { return S->cofm(); }
    friend vect const&pos   (CGCS) { return S->cofm(); }
    friend uint const&mybody(CGCS) { return S->mybody(); }
    friend ten1 const acc   (CGCS) { return ten1(S->SINKPT+1); }
    friend real const&acc   (CGCS,
			     int i){ return S->SINKPT[i+1]; }
    friend real const&pot   (CGCS) { return S->SINKPT[0]; }
    friend real const&rho   (CGCS) { return S->SINKPT[0]; }
    friend uint const&num   (CGCS) { return static_cast<uint*>(
				 	    static_cast<void*>(S->SINKPT))[0]; }
#ifdef ALLOW_INDI
    friend real const&eph   (CGCS) { return S->EPH; }
    friend real const size  (CGCS) { return twice(S->EPH); }
    friend real const&sizeq (CGCS) { return S->SINKPT[1]; }
#endif
    friend real const&mass  (CGCS) { return S->MASS; }
#undef CGCS
    //--------------------------------------------------------------------------
    // simple manipulations                                                     
    //--------------------------------------------------------------------------
    void set_mass(const areal*const&m) {
      MASS = m[mybody()];
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    template<typename bodies_type>
    void set_mass(const bodies_type*const&B) {
      MASS = B->mas(mybody());
    }
    //--------------------------------------------------------------------------
    void set_mass_and_pos(const areal*const&m, const areal*x[NDIM]) {
      MASS      = m   [mybody()];
      cofm()[0] = x[0][mybody()];
      cofm()[1] = x[1][mybody()];
#if NDIM==3
      cofm()[2] = x[2][mybody()];
#endif
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    template<typename bodies_type>
    void set_mass_and_pos(const bodies_type*const&B) {
      MASS   = B->mas(mybody());
      cofm() = B->pos(mybody());
    }
#ifdef ALLOW_INDI
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void copy_eph(const bodies_type*const&B) {
      eph() = half*B->eps(mybody());
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void copy_eph(const areal*const&E) { 
      eph() = half*E[mybody()];
    }
#endif
    //--------------------------------------------------------------------------
    void reset_srce() {
      SINKPT[0] = zero;
      SINKPT[1] = zero;
      SINKPT[2] = zero;
#if NDIM > 2
      SINKPT[3] = zero;
#endif
    }
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void update_dens(const bodies_type*const&B) const {
      B->rho(mybody()) = nbdy::rho(this);
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void update_dens(areal*const&rh) const {
      rh[mybody()] = nbdy::rho(this);
    }
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void update_grav(const bodies_type*const&B) const {
      B->pot(mybody()) = nbdy::pot(this);
      B->acc(mybody()) = nbdy::acc(this);
    }
    //--------------------------------------------------------------------------
    void update_grav(areal*a[NDIM], areal*const&p) const {
      if(p) p[mybody()] = nbdy::pot(this);
      a[0][mybody()] = nbdy::acc(this,0);
      a[1][mybody()] = nbdy::acc(this,1);
#if NDIM>2
      a[2][mybody()] = nbdy::acc(this,2);
#endif
    }
#ifdef ALLOW_INDI
    //--------------------------------------------------------------------------
    void update_eps(const sbodies*const&B) const {
     B->eps(mybody()) = twice(nbdy::eph(this));
    }
#ifdef ALLOW_MPI
    //--------------------------------------------------------------------------
    void update_eps(const pbodies*const&B) const {
     B->eps(mybody()) = twice(nbdy::eph(this));
    }
#endif
    //--------------------------------------------------------------------------
    void update_eps(areal*const&ep) const {
      ep[mybody()] = twice(nbdy::eph(this));
    }
#endif
    //--------------------------------------------------------------------------
    template<typename bodies_type>
    void update_num(const bodies_type*const&B) const {
      B->num(mybody()) = nbdy::num(this);
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void update_num(int*const&num) const {
      num[mybody()] = nbdy::num(this);
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
#if   P_ORDER > 3
      N_P3    = N_P2   + ten2::NDAT,
# if  P_ORDER > 4
      N_P4    = N_P3   + ten3::NDAT,
#  if P_ORDER > 5
#   error "expansion order > 5 not supported in public/grav.h"
#  endif
      N_EPH   = N_P4   + ten4::NDAT,
# else
      N_EPH   = N_P3   + ten3::NDAT,
# endif
#else
      N_EPH   = N_P2   + ten2::NDAT,
#endif
#ifdef ALLOW_INDI
      N_TOTAL = N_EPH  + 1;
    static const int N_eph() { return N_EPH; }
#else
      N_TOTAL = N_EPH;
#endif
    static const int N_tot() { return N_TOTAL; }
    //--------------------------------------------------------------------------
    static const int N_COEFF = 1 + ten1::NDAT + ten2::NDAT + ten3::NDAT
#if  P_ORDER > 3
    + ten4::NDAT
#if  P_ORDER > 4
    + ten5::NDAT
#if  P_ORDER > 5
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
    //--------------------------------------------------------------------------
    // constructor                                                              
    //--------------------------------------------------------------------------
  public:
    grav_cell() : COEFFS(0) {}
    //--------------------------------------------------------------------------
    // simple manipulations                                                     
    //--------------------------------------------------------------------------
  public:
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
#ifdef ALLOW_INDI
    real  &eph   () { return SOURCE[N_EPH]; }
#endif
    real* &coeffs() { return COEFFS; }
    ten2   quad  () { return ten2(SOURCE+N_P2); }
#if   P_ORDER > 3
    ten3   octo  () { return ten3(SOURCE+N_P3); }
# if  P_ORDER > 4
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
#ifdef ALLOW_INDI
    friend real const&eph   (CCCC) { return C->SOURCE[N_EPH]; }
#endif
    friend real*const&coeffs(CCCC) { return C->COEFFS; }
    friend ten2       quad  (CCCC) { return ten2(C->SOURCE+N_P2); }
#if   P_ORDER > 3
    friend ten3       octo  (CCCC) { return ten3(C->SOURCE+N_P3); }
# if  P_ORDER > 4
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
    void add_sink_flag(const grav_soul* const&S) { add_sink_flag_from_soul(S); }
    void add_sink_flag(const grav_cell* const&C) { add_sink_flag_from_cell(C); }
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
      for(register int i=0; i<NDIM; i++)
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
  class grav_tree : public basic_tree<grav_tree,grav_cell> {
    //--------------------------------------------------------------------------
    friend class tree_transport<grav_tree>;        // for parallel code         
    //--------------------------------------------------------------------------
    grav_tree            (grav_tree const&);       // not implemented           
    grav_tree& operator= (grav_tree const&);       // not implemented           
    //--------------------------------------------------------------------------
    // data:                                                                    
    //--------------------------------------------------------------------------
#ifdef ALLOW_INDI
    const soft_type   SOFT;                        // global / individual       
#endif
    const areal      *M;                           // array with masses         
#ifdef ALLOW_INDI
    const areal      *EP;                          // array with eps_i          
#endif
    const int         nCS;                         // source:    memory / cell  
    static const int  nSS = NDIM+1;                // sink part: memory / soul  
    static const int  nCC = grav_cell::N_COEFF;    // coeffs:    memory / cell  
    int               Ncs, Nss;                    // # cell/soul sinks         
    real             *CELL_SOURCE;                 // memory for cell source    
    real             *CELL_COEFFS;                 // memory for cell sink part 
    real             *SOUL_SINKPT;                 // memory for soul sink part 
    //--------------------------------------------------------------------------
    inline void reset_soul_sinkpt();               // set grav_soul::SINKPT     
    inline void reset_cell_source();               // set grav_cell::SOURCE     
    inline void reset_cell_coeffs();               // for parallel code only    
    inline void   set_soul_sinkpt();               // set grav_soul::SINKPT     
    inline void   set_cell_source();               // set grav_cell::SOURCE     
    inline void   set_cell_coeffs();               // for parallel code only    
    //--------------------------------------------------------------------------
    // construction:                                                            
    // - sets up the oct tree structure from scratch                            
    // - gives source memory to cells (souls have their own)                    
    // - recursively computes for every cell: mass, cofm, rmax                  
    // - does NOT give sink memory to cells nor to souls                        
    //--------------------------------------------------------------------------
  public:
    //--------------------------------------------------------------------------
    //   construction from list of bodies                                       
    grav_tree(const sbodies* const&,               // I: sbodies                
#ifdef ALLOW_INDI
	      const soft_type  = Default::soften,  //[I: global/individual]     
#endif
	      const int        = Default::Ncrit);  //[I: N_crit]                
#ifdef ALLOW_MPI
    //--------------------------------------------------------------------------
    grav_tree(const pbodies* const&,               // I: pbodies                
#ifdef ALLOW_INDI
	      const soft_type  = Default::soften,  //[I: global/individual]     
#endif
	      const int        = Default::Ncrit);  //[I: N_crit]                
    //--------------------------------------------------------------------------
    grav_tree(const pbodies* const&,               // I: pbodies                
	      vect           const&,               // I: x_min                  
	      vect           const&,               // I: x_max                  
#ifdef ALLOW_INDI
	      const soft_type  = Default::soften,  //[I: global/individual]     
#endif
	      const int        = Default::Ncrit);  //[I: N_crit]                
#endif
    //--------------------------------------------------------------------------
    //   construction from list of bodies and boxing position x_min, x_max      
    grav_tree(const sbodies* const&,               // I: sbodies                
	      vect           const&,               // I: x_min                  
	      vect           const&,               // I: x_max                  
#ifdef ALLOW_INDI
	      const soft_type  = Default::soften,  //[I: global/individual]     
#endif
	      const int        = Default::Ncrit);  //[I: N_crit]                
    //--------------------------------------------------------------------------
    //   construction from arrays of flags and positions                        
    grav_tree(const int        *,                  // I: array with flags       
	      const areal      *[NDIM],            // I: arrays with x,y,z      
	      const areal      *,                  // I: array  with m_i        
#ifdef ALLOW_INDI
	      const areal      *,                  // I: array  with eps_i      
#endif
	      const uint        ,                  // I: size of arrays         
#ifdef ALLOW_INDI
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
    // - update souls sink flags                                                
    // - pass sink flags up the tree                                            
    // - give sinkpt memory to sink souls                                       
    void prepare_density    (bool = false);        //[I: re-use memory?]        
    //--------------------------------------------------------------------------
    // prepare for exact gravity computation:                                   
    // - update souls' and cells' sink flags                                    
    // - update souls' eph_i                       (individual softening only)  
    // - optionally: adjust sink's eph_i           (individual softening only)  
    // - update cell's eph_i (pass up the tree)    (individual softening only)  
    // - give sinkpt memory to sink souls                                       
    // - reset souls' pot & acc                                                 
    void prepare_grav_exact (                      // prepare for direct sums   
#ifdef ALLOW_INDI
			     real = zero,          //[I: Nsoft: adjust sinks]   
			     real = one,           //[I: emax:  adjust sinks]   
			     uint = 0u,            //[I: Nref:  adjust sinks]   
			     real = two,           //[I: max change in eps_i]   
#endif
			     bool = false);        //[I: re-use memory?]        
    //--------------------------------------------------------------------------
    // prepare for gravity approximation:                                       
    // - update souls' and cells' sink flags                                    
    // - update souls' eph_i                       (individual softening only)  
    // - optionally: adjust sink's eph_i           (individual softening only)  
    // - update cell's eph_i (pass up the tree)    (individual softening only)  
    // - give sinkpt memory to sink souls                                       
    // - reset souls' pot & acc                                                 
    // - optionally give coeffs memory to sink cells                            
    // - recursively compute the multipoles                                     
    // - set the cells r_crit & r_crit^2                                        
    void prepare_grav_approx(const grav_mac*,      // I: MAC                    
			     bool = false,         //[I: give cell coeffs?]     
#ifdef ALLOW_INDI
			     real = zero,          //[I: Nsoft: adjust sinks]   
			     real = one,           //[I: emax:  adjust sinks]   
			     uint = 0u,            //[I: Nref:  adjust sinks]   
			     real = two,           //[I: max change in eps_i]   
#endif
			     bool = false);        //[I: re-use memory?]        
    //--------------------------------------------------------------------------
    // prepare for neighbour counting:                                          
    // - update souls' and cells' sink flags                                    
    // - give sinkpt memory to sink souls                                       
    // - update souls' eph_i                       (individual softening only)  
    // - souls: set size^2=eph^2, num=0                                         
    // - cell sinks: update size (pass up tree)                                 
    // - update cell's eph_i (pass up the tree)    (individual softening only)  
    void prepare_neighbour_counting(               //                           
				    const real* =0,//[I: global body size]      
				    bool = false); //[I: re-use memory?]        
    //--------------------------------------------------------------------------
    int         const& N_cell_sinks() const { return Ncs; }
    int         const& N_soul_sinks() const { return Nss; }
    const areal*const& my_masses   () const { return M; }
#ifdef ALLOW_INDI
    const areal*const& my_eps      () const { return EP; }
#endif
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
#ifdef ALLOW_INDI
  inline real const&eph      (CCCC) { return eph      (I.c_pter()); }
#endif
  inline real*const&coeffs   (CCCC) { return coeffs   (I.c_pter()); }
  inline ten2       quad     (CCCC) { return quad     (I.c_pter()); }
# if   P_ORDER > 3
  inline ten3       octo     (CCCC) { return octo     (I.c_pter()); }
#  if  P_ORDER > 4
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
  // class grav_stat                                                          //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class grav_stat {
  private:
    typedef grav_tree::cell_iterator cell_iter;
    typedef grav_tree::soul_iterator soul_iter;
    uint N_BB, N_CB, M_CB, N_CC, M_CC, M_CS;
#ifdef ENHANCED_IACT_STATS
    double D_C, D_CC;
#endif
    static int trick(const real x, const int w)
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
    void reset    () { N_BB=0; N_CB=0; M_CB=0; N_CC=0; M_CC=0; M_CS=0;
#ifdef ENHANCED_IACT_STATS
                       D_CC=0.;
#endif
    }
    void record_BB() { N_BB++; }
    void record_taylor_CB(cell_iter const&A, soul_iter const&B) { N_CB++; }
    void record_taylor_CC(cell_iter const&A, cell_iter const&B) { N_CC++;
#ifdef ENHANCED_IACT_STATS
      D_CC += abs(A.index()-B.index());
#endif
    }
    void record_direct_CB(cell_iter const&A, soul_iter const&B) { M_CB++; }
    void record_direct_CC(cell_iter const&A, cell_iter const&B) { M_CC++; }
    void record_direct_CS(cell_iter const&A)                    { M_CS++; }
    // 2 reporting                                                             
    const uint& BB_iacts()        const { return N_BB; }
    const uint& CB_taylor_iacts() const { return N_CB; }
    const uint& CC_taylor_iacts() const { return N_CC; }
    const uint& CB_direct_iacts() const { return M_CB; }
    const uint& CC_direct_iacts() const { return M_CC; }
    const uint& CS_direct_iacts() const { return M_CS; }
#ifdef ENHANCED_IACT_STATS
    const int   CC_taylor_dist() const { return int(D_CC / N_CC); }
#endif
    // 3 writing stats to ostream                                               
    void write(std::ostream&out) const {
      register uint total = BB_iacts()
	+ CB_taylor_iacts()
	+ CB_direct_iacts()
	+ CC_taylor_iacts()
	+ CC_direct_iacts()
	+ CS_direct_iacts(), part;
      register real percent;
      out<<" interaction statitics:\n"
	 <<"     type          approx   direct      total\n"
	 <<" # body-body :          - ";
      part    = BB_iacts();
      percent = 100.*part/real(total);
      out<<setw( 8)<<BB_iacts()<<" "
	 <<setw(10)<<part<<" = "
	 <<setprecision(trick(percent,5))<<setw(8)<<percent<<"%\n"
	 <<" # cell-body : ";
      part    = CB_taylor_iacts()
	+CB_direct_iacts();
      percent = 100.*part/real(total);
      out<<setw(10)<<CB_taylor_iacts()<<" "
	 <<setw( 8)<<CB_direct_iacts()<<" "
	 <<setw(10)<<part<<" = "
	 <<setprecision(trick(percent,5))<<setw(8)<<percent<<"%\n"
	 <<" # cell-cell : ";
      part    = CC_taylor_iacts() + CC_direct_iacts();
      percent = 100.*part/real(total);
      out<<setw(10)<<CC_taylor_iacts()<<" "
	 <<setw( 8)<<CC_direct_iacts()<<" "
	 <<setw(10)<<part<<" = "
	 <<setprecision(trick(percent,5))<<setw(8)<<percent<<"%\n"
	 <<" # cell-self :          - ";
      part    = CS_direct_iacts();
      percent = 100.*part/real(total);
      out<<setw( 8)<<CS_direct_iacts()<<" "
	 <<setw(10)<<part<<" = "
	 <<setprecision(trick(percent,5))<<setw(8)<<percent<<"%\n"
	 <<" # total     : ";
      part=BB_iacts()
	+CB_direct_iacts()
	+CC_direct_iacts()
	+CS_direct_iacts();
      out<<setw(10)<<CB_taylor_iacts()+CC_taylor_iacts()<<" "
	 <<setw( 8)<<part<<" "
	 <<setw(10)<<total<<" =  100.000%\n";
#ifdef ENHANCED_IACT_STATS
      out<<" mean memory distance between interacting cells: "
	 <<CC_taylor_dist()<<"\n";
#endif
    }
  };
}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
#endif                                             // included_grat_h           
