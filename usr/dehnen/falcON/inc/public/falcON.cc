//-----------------------------------------------------------------------------+
//                                                                             |
// falcON.cc                                                                   |
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
// defines all function members of class falcON                                |
//                                                                             |
// Contents                                                                    |
//                                                                             |
// 1  Constructor, destructor and such                                         |
// 2  tree growth and re-use                                                   |
// 3  information and statistics                                               |
// 4  estimation of density and such                                           |
// 5  approximation of gravity                                                 |
// 6  SPH support                                                              |
// 7  neighbour finding and counting                                           |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef falcON_included_falcON_h
#  include <falcON.h>                      // class falcON                      
#endif
#ifndef falcON_included_iomanip_h
#  include <iomanip>                       // C++ file I/O                      
#  define falcON_included_iomanip_h
#endif
#ifndef falcON_included_fstream_h
#  include <fstream>                       // C++ file I/O                      
#  define falcON_included_fstream_h
#endif
#ifndef falcON_included_grav_h
#  include <public/grav.h>                 // gravity stuff                     
#endif
#ifndef falcON_included_stic_h
#  include <public/stic.h>                 // collision partner search          
#endif
#ifdef falcON_SPH
#  ifndef falcON_included_spht_h
#    include <sph/spht.h>                  // full-fledged SPH support          
#  endif
#endif

namespace nbdy {
  //---------------------------------------------------------------------------+
  //                                                                           |
  // CONSTRUCTORS, DESTRUCTOR AND SUCH                                         |
  //                                                                           |
  //---------------------------------------------------------------------------+
#ifdef falcON_INDI
#  define I_SOFT i_soft
#else
#  define I_SOFT 0
#endif
  inline falcON::falcON(const sbodies*  b,
			const real      e,
			const real      th,
			const kern_type k,
#ifdef falcON_INDI
			const bool      i_soft,
#endif
			const real      g,
			const MAC_type  mac,
			const int       gd[4]
#ifdef falcON_SPH
		       ,const int       sd[3]
#endif
			) :
    STATS   ( falcON_Memory(new grav_stat()) ),
    BODIES  ( b ),
    ARRAYS  ( 0 ),
    Ncrit   ( 0 ),
    TREE    ( 0 ),
    GMAC    ( new grav_mac(mac, abs(th), falcON_ORDER) ),
    GRAV    ( new grav_estimator(TREE,k,STATS,e,g,I_SOFT,gd) ),
    STSP    ( new stsp_estimator(TREE) )
#ifdef falcON_SPH
					     ,
    SPHT    ( new sph_estimator(TREE,sd) )
#endif
  {
#ifdef falcON_INDI
    if( BODIES->has(io::e) && !i_soft)
      warning("eps_i given but global softening chosen");
    if(!BODIES->has(io::e) &&  i_soft)
      error("no eps_i given, but individual softening chosen");
#endif
  }
  //----------------------------------------------------------------------------
  inline falcON::falcON(const abodies  *b,
			const real      e,
			const real      th,
			const kern_type k,
#ifdef falcON_INDI
			const bool      i_soft,
#endif
			const real      g,
			const MAC_type  mac,
			const int       gd[4]
#ifdef falcON_SPH
		       ,const int       sd[3]
#endif
			) :
    STATS   ( falcON_Memory(new grav_stat()) ),
    BODIES  ( 0 ),
    ARRAYS  ( b ),
    Ncrit   ( 0 ),
    TREE    ( 0 ),
    GMAC    ( new grav_mac(mac, abs(th), falcON_ORDER) ),
    GRAV    ( new grav_estimator(TREE,k,STATS,e,g,I_SOFT,gd) ),
    STSP    ( new stsp_estimator(TREE) )
#ifdef falcON_SPH
					     ,
    SPHT    ( new sph_estimator(TREE,sd) )
#endif
  {
#ifdef falcON_INDI
    if( ARRAYS->has(io::e) && !i_soft)
      warning("eps_i given but global softening chosen");
    if(!ARRAYS->has(io::e) &&  i_soft)
      error("no eps_i given, but individual softening chosen");
#endif
  }
#undef I_SOFT
  //----------------------------------------------------------------------------
  inline falcON::~falcON()
  {
    if(TREE) delete TREE;
    delete GMAC;
    delete STATS;
    delete GRAV;
    delete STSP;
#ifdef falcON_SPH
    delete SPHT;
#endif
  }
  //----------------------------------------------------------------------------
#ifdef falcON_INDI
  inline bool const &falcON::use_individual_eps() const
  {
    return GRAV->use_indiv_eps();
  }
#endif
  //----------------------------------------------------------------------------
  inline void falcON::reset_softening (const real      e,
				       const kern_type k) const
  {
    GRAV->reset_softening(e,k);
  }
  //----------------------------------------------------------------------------
  inline void falcON::reset_NewtonsG (const real g) const
  {
    GRAV->reset_NewtonsG(g);
  }
  //----------------------------------------------------------------------------
  inline void falcON::reset_opening(const real     th,
				    const MAC_type mac) const
  {
    GMAC->reset(mac, abs(th), falcON_ORDER);
  }
  //---------------------------------------------------------------------------+
  //                                                                           |
  // TREE GROWTH AND RE-USE                                                    |
  //                                                                           |
  //---------------------------------------------------------------------------+
  inline void falcON::grow(int        const&Ncr,
			   const vect*const&x0)
    // This routines grows a fresh tree or re-grows an existing tree            
  {
    SET_I
      Ncrit = max(Ncr,1);
    if(TREE) {
      TREE->build(Ncrit,x0);
      GRAV->reset();
      STSP->reset();
#ifdef falcON_SPH
      SPHT->reset();
#endif
      SET_T(" time: oct_tree::build():               ");
    } else {
      TREE = BODIES? new oct_tree(BODIES,Ncrit,x0) : 
	new oct_tree(ARRAYS,Ncrit,x0) ;
      GRAV->new_tree(TREE);
      STSP->new_tree(TREE);
#ifdef falcON_SPH
      SPHT->new_tree(TREE);
#endif
      SET_T(" time: oct_tree::oct_tree():            ");
    }
  }
  //----------------------------------------------------------------------------
#if (0)
  inline void falcON::re_grow(int        const&Ncut,
			      int        const&Ncr,
			      const vect*const&x0)
    // This routines re-builds an existing tree                                 
  {
    SET_I
      Ncrit = max(Ncr,1);
    if(TREE) {
      TREE->rebuild(Ncut,Ncrit);
      GRAV->reset();
      STSP->reset();
#ifdef falcON_SPH
      SPHT->reset();
#endif
      SET_T(" time: oct_tree::rebuild():             ");
    } else {
      TREE = BODIES? new oct_tree(BODIES,Ncrit,x0) :
	new oct_tree(ARRAYS,Ncrit,x0) ;
      GRAV->new_tree(TREE);
      STSP->new_tree(TREE);
#ifdef falcON_SPH
      SPHT->new_tree(TREE);
#endif
      SET_T(" time: oct_tree::oct_tree():            ");
    }
  }
#endif
  //----------------------------------------------------------------------------
  inline void falcON::reuse()
  {
    // This routine merely updates the leaf's positions                         
    SET_I
      if(TREE == 0) {
	warning("no old tree to be re-used","falcON::reuse()");
	return grow();
      }
    TREE->reuse();
    GRAV->reset();
    STSP->reset();
#ifdef falcON_SPH
    SPHT->reset();
#endif
  }
  //---------------------------------------------------------------------------+
  //                                                                           |
  // INFORMATION AND STATISTICS                                                |
  //                                                                           |
  //---------------------------------------------------------------------------+
  inline kern_type const &falcON::kernel() const
  {
    return GRAV->kernel();
  }
  //----------------------------------------------------------------------------
  inline const char*  falcON::describe_kernel () const
  {
    return describe(GRAV->kernel());
  }
  //----------------------------------------------------------------------------
  inline real const &falcON::eps() const
  {
    return GRAV->softening_length();
  }
  //----------------------------------------------------------------------------
  inline real const &falcON::softening_length() const
  {
    return GRAV->softening_length();
  }
  //----------------------------------------------------------------------------
  inline real const &falcON::NewtonsG() const
  {
    return GRAV->NewtonsG();
  }
  //----------------------------------------------------------------------------
  inline const MAC_type &falcON::MAC() const
  {
    return GMAC->method();
  }
  //----------------------------------------------------------------------------
  inline const char* falcON::describe_MAC() const
  {
    return GMAC->describe_method();
  }
  //----------------------------------------------------------------------------
  inline const uint &falcON::No_bodies() const
  {
    return BODIES->N_bodies();
  }
  //----------------------------------------------------------------------------
  inline uint falcON::No_bodies_used() const
  {
    return TREE->N_leafs();
  }
  //----------------------------------------------------------------------------
  inline uint falcON::No_cells_used() const
  {
    return TREE->N_cells();
  }
  //----------------------------------------------------------------------------
  inline void falcON::dump_nodes(const char* fcells, const char* fleafs) const
  {
    if(fcells) {
      std::ofstream o(fcells);
      GRAV->dump_cells(o);
    }
    if(fleafs) {
      std::ofstream o(fleafs);
      GRAV->dump_leafs(o);
    }
  }
  //----------------------------------------------------------------------------
  inline void falcON::stats(std::ostream& out) const
  {
    out<<"\n state:                ";
    if(TREE) {
      if     (TREE->is_fresh()   ) out<<" tree new\n";
      else if(TREE->is_re_grown()) out<<" tree re-grown\n"; 
      else if(TREE->is_re_used() ) out<<" tree re-used\n";
      out<<" root center:           "<<TREE->root_center()    <<'\n'
	 <<" root radius:           "<<radius(TREE->root())   <<'\n'
	 <<" bodies loaded:         "<<number(TREE->root())   <<'\n';
      if(GRAV->N_coeffs())
	out<<" total mass:            "<<mass  (GRAV->root())<<'\n';
      out<<" N_crit:                "<<Ncrit                  <<'\n'
	 <<" cells used:            "<<TREE->N_cells()        <<'\n';
      if(GRAV->N_coeffs())
	out<<" of which were active   "<<GRAV->N_active_cells()<<'\n';
      out<<" maximum depth:         "<<TREE->depth()          <<'\n'
	 <<" current theta:         "<<GMAC->theta_min()      <<'\n'
	 <<" current MAC:           "<<GMAC->describe_method()<<'\n';
#ifdef falcON_INDI
      if(GRAV->use_indiv_eps())
	out<<" softening:             individual\n";
      else
	out<<" softening:             global\n"
	   <<" softening length:      "<<GRAV->softening_length()<<'\n';
#else
      out<<" softening length:      "<<GRAV->softening_length()<<'\n';
#endif
      out<<" softening kernel:      "<<describe_kernel()      <<'\n';
      if(TREE->is_used_for_grav()) {               // IF just had gravity       
	out
	  <<" Taylor coeffs used:    "<<GRAV->N_coeffs()
	  <<" in "<<GRAV->N_chunks()
	  <<" chunks of "<<GRAV->N_elems_in_chunk()<<'\n';
	STATS->write(out);
      }
    } else out<<" no tree\n";
  }
  //----------------------------------------------------------------------------
  inline uint const&falcON::No_coeffs_used() const
  {
    return GRAV->N_coeffs();
  }
  //----------------------------------------------------------------------------
  inline uint falcON::BB_interactions() const
  {
    return STATS->BB_direct_iacts();
  }
  //----------------------------------------------------------------------------
  inline uint falcON::MB_interactions() const
  {
    return
      STATS->CB_direct_iacts()+
      STATS->CC_direct_iacts()+
      STATS->CX_direct_iacts();
  }
  //----------------------------------------------------------------------------
  inline uint falcON::CB_interactions() const
  {
    return STATS->CB_approx_iacts();
  }
  //----------------------------------------------------------------------------
  inline uint falcON::CC_interactions() const
  {
    return STATS->CC_approx_iacts();
  }
  //----------------------------------------------------------------------------
  inline uint falcON::total_interactions () const
  {
    return 
      BB_interactions()+
      MB_interactions()+
      CB_interactions()+
      CC_interactions();
  }
  //----------------------------------------------------------------------------
  inline const vect& falcON::root_center() const
  {
    return TREE->root_center(); }
  //----------------------------------------------------------------------------
  inline const real& falcON::root_radius() const
  {
    return radius(TREE->root());
  }
  //----------------------------------------------------------------------------
  inline const uint& falcON::root_depth () const
  {
    return TREE->depth();
  }
  //----------------------------------------------------------------------------
  inline const int&        falcON::root_number() const
  {
    return number(TREE->root());
  }
  //----------------------------------------------------------------------------
  inline const real& falcON::root_mass  () const
  {
    return mass(GRAV->root());
  }
  //---------------------------------------------------------------------------+
  //                                                                           |
  // ESTIMATION OF DENSITIES                                                   |
  //                                                                           |
  //---------------------------------------------------------------------------+
  inline void falcON::estimate_n(uint const&Nx, bool const&all)
  {
    if(BODIES && !BODIES->has(io::r) || ARRAYS && !ARRAYS->has(io::r))
      error("[falcON::estimate_nd()]: nobody has memory for rho");
    GRAV->estimate_nd(all,Nx);
  }
  //----------------------------------------------------------------------------
  inline void falcON::estimate_sd(uint const&Nx, bool const&all)
  {
    if(BODIES && !BODIES->has(io::r) || ARRAYS && !ARRAYS->has(io::r))
      error("[falcON::estimate_sd()]: nobody has memory for rho");
    GRAV->estimate_sd(all,Nx);
  }
  //----------------------------------------------------------------------------
  inline void falcON::estimate_rho(uint const&Nx, bool const&all)
  {
    if(BODIES && !BODIES->has(io::r) || ARRAYS && !ARRAYS->has(io::r))
      error("[falcON::estimate_md()]: nobody has memory for rho");
    GRAV->estimate_md(all,Nx);
  }
  //---------------------------------------------------------------------------+
  //                                                                           |
  // APPROXIMATION OF GRAVITY                                                  |
  //                                                                           |
  //---------------------------------------------------------------------------+
  inline void falcON::approximate_gravity(bool const&split,
					  bool const&all
#ifdef falcON_ADAP
					  ,
					  real const&Nsoft,
					  uint const&Nref,
					  real const&emin,
					  real const&efac
#  define ADAP_ARGS ,Nsoft,Nref,emin,efac
#else
#  define ADAP_ARGS
#endif
					  )
  {
    GRAV->approx(GMAC,all,split ADAP_ARGS);
  }
  //----------------------------------------------------------------------------
  inline void falcON::exact_gravity(bool const&all
#ifdef falcON_ADAP
				    ,real const&Nsoft,
				    uint const&Nref,
				    real const&emin,
				    real const&efac
#endif
				    )
  {
    GRAV->exact(all ADAP_ARGS);
  }
  //---------------------------------------------------------------------------+
  //                                                                           |
  // SPH SUPPORT                                                               |
  //                                                                           |
  //---------------------------------------------------------------------------+
#ifdef falcON_SPH
  inline void falcON::count_SPH_partners(bool const&al,
					 bool const&up)
  {
    SPHT->count_partners(al,up);
  }
  //----------------------------------------------------------------------------
  inline void falcON::adjust_SPH_sizes(uint const&nh,
				       real const&hm,
				       bool const&al,
				       uint       ni,
				       uint       nx,
				       uint const&ix)
  {
    SPHT->adjust_sizes(nh,hm,al,ni,nx,ix);
  }
  //----------------------------------------------------------------------------
  inline void falcON::SPH_sweep_one(const eq_of_state*const&es,
				    uint              const&nh,
				    bool              const&al,
				    bool              const&up)
  
  {
    SPHT->sweep_one(es,nh,al,up);
  }
  //----------------------------------------------------------------------------
  inline void falcON::SPH_sweep_two()
  
  {
    SPHT->sweep_two();
  }
#endif
  //---------------------------------------------------------------------------+
  //                                                                           |
  // PARTNER FINDING AND COUNTING                                              |
  //                                                                           |
  //---------------------------------------------------------------------------+
  inline void falcON::make_iaction_list(elem_pair *bl,
					uint const&nl,
					uint      &na,
					bool const&Mx,
					real const&tau)
  {
    if(tau<zero) STSP->make_sph_list   (bl,nl,na,Mx);
    else         STSP->make_sticky_list(bl,nl,na,tau);
  }
  //---------------------------------------------------------------------------+
  inline void falcON::count_sph_partners()
  {
    STSP->count_sph_partners();
  }
  //////////////////////////////////////////////////////////////////////////////
}                                                  // END: namespace nbdy       
////////////////////////////////////////////////////////////////////////////////
