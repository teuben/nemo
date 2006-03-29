//-----------------------------------------------------------------------------+
//                                                                             |
// FAlCON.cc                                                                   |
//                                                                             |
// Copyright (C) 2000-2005  Walter Dehnen                                      |
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
// defines all function members of class FAlCON                                |
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
#ifndef falcON_included_FAlCON_cc
#  define falcON_included_FAlCON_cc

#ifndef falcON_included_FAlCON_h
#  include <FAlCON.h>                      // class FAlCON                      
#endif
#ifndef falcON_included_iomanip_h
#  include <iomanip>                       // C++ I/O formatting                
#  define falcON_included_iomanip_h
#endif
#ifndef falcON_included_fstream_h
#  include <fstream>                       // C++ file I/O                      
#  define falcON_included_fstream_h
#endif
#ifndef falcON_included_body_h
#  include <body.h>                        // bodies                            
#endif
#ifndef falcON_included_gravity_h
#  include <public/gravity.h>              // gravity stuff                     
#endif
#ifndef falcON_included_partner_h
#  include <public/partner.h>              // partner stuff                     
#endif
#ifdef falcON_SPH
#  ifndef falcON_included_sph_h
#    include <sph/sph.h>                   // full-fledged SPH support          
#  endif
#endif

namespace falcON {
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
  inline FAlCON::FAlCON(const bodies*   b,
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
    STATS   ( new GravStats() ),
    BODIES  ( b ),
    Ncrit   ( 0 ),
    TREE    ( 0 ),
    GMAC    ( new GravMAC(mac, abs(th), falcON_ORDER) ),
    GRAV    ( new GravEstimator(TREE,k,STATS,e,g,I_SOFT,gd) ),
    PAES    ( 0 )
#ifdef falcON_SPH
			  ,
    SPHT    ( new SphEstimator(TREE,sd) )
#endif
  {}
  //----------------------------------------------------------------------------
  inline FAlCON::~FAlCON()
  {
    if(TREE) delete TREE;
    delete GMAC;
    delete STATS;
    delete GRAV;
    if(PAES) delete PAES;
#ifdef falcON_SPH
    delete SPHT;
#endif
  }
  //----------------------------------------------------------------------------
#ifdef falcON_INDI
  inline bool const &FAlCON::use_individual_eps() const
  {
    return GRAV->use_indiv_eps();
  }
#endif
  //----------------------------------------------------------------------------
  inline void FAlCON::reset_softening (const real      e,
				       const kern_type k) const
  {
    GRAV->reset_softening(e,k);
  }
  //----------------------------------------------------------------------------
  inline void FAlCON::reset_NewtonsG (const real g) const
  {
    GRAV->reset_NewtonsG(g);
  }
  //----------------------------------------------------------------------------
  inline void FAlCON::reset_opening(const real     th,
				    const MAC_type mac) const
  {
    GMAC->reset(mac, abs(th), falcON_ORDER);
  }
  //---------------------------------------------------------------------------+
  //                                                                           |
  // TREE GROWTH AND RE-USE                                                    |
  //                                                                           |
  //---------------------------------------------------------------------------+
  inline void FAlCON::grow(int        const&Ncr,
			   const vect*const&x0)
    // This routines grows a fresh tree or re-grows an existing tree            
  {
    SET_I
      Ncrit = max(Ncr,1);
    if(TREE) {
      TREE->build(Ncrit,x0);
      GRAV->reset();
#ifdef falcON_SPH
      SPHT->reset();
#endif
      SET_T(" time: OctTree::build():               ");
    } else {
      TREE = new OctTree(BODIES,Ncrit,x0);
      GRAV->new_tree(TREE);
#ifdef falcON_SPH
      SPHT->new_tree(TREE);
#endif
      SET_T(" time: OctTree::OctTree():            ");
    }
  }
  //----------------------------------------------------------------------------
  inline void FAlCON::reuse()
  {
    // This routine merely updates the leaf's positions                         
    SET_I
      if(TREE == 0) {
	warning("no old tree to be re-used","FAlCON::reuse()");
	return grow();
      }
    TREE->reuse();
    GRAV->reset();
#ifdef falcON_SPH
    SPHT->reset();
#endif
  }
  //---------------------------------------------------------------------------+
  //                                                                           |
  // INFORMATION AND STATISTICS                                                |
  //                                                                           |
  //---------------------------------------------------------------------------+
  inline kern_type const &FAlCON::kernel() const
  {
    return GRAV->kernel();
  }
  //----------------------------------------------------------------------------
  inline const char*  FAlCON::describe_kernel () const
  {
    return describe(GRAV->kernel());
  }
  //----------------------------------------------------------------------------
  inline real const &FAlCON::eps() const
  {
    return GRAV->softening_length();
  }
  //----------------------------------------------------------------------------
  inline real const &FAlCON::softening_length() const
  {
    return GRAV->softening_length();
  }
  //----------------------------------------------------------------------------
  inline real const &FAlCON::NewtonsG() const
  {
    return GRAV->NewtonsG();
  }
  //----------------------------------------------------------------------------
  inline const MAC_type &FAlCON::MAC() const
  {
    return GMAC->method();
  }
  //----------------------------------------------------------------------------
  inline const char* FAlCON::describe_MAC() const
  {
    return GMAC->describe_method();
  }
  //----------------------------------------------------------------------------
  inline unsigned FAlCON::No_bodies() const
  {
    return BODIES->N_bodies();
  }
  //----------------------------------------------------------------------------
  inline unsigned FAlCON::No_bodies_used() const
  {
    return TREE->N_leafs();
  }
  //----------------------------------------------------------------------------
  inline unsigned FAlCON::No_cells_used() const
  {
    return TREE->N_cells();
  }
  //----------------------------------------------------------------------------
  inline void FAlCON::dump_nodes(const char* fcells, const char* fleafs) const
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
  inline void FAlCON::stats(std::ostream& out) const
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
  inline unsigned const&FAlCON::No_coeffs_used() const
  {
    return GRAV->N_coeffs();
  }
  //----------------------------------------------------------------------------
  inline unsigned FAlCON::BB_interactions() const
  {
    return STATS->BB_direct_iacts();
  }
  //----------------------------------------------------------------------------
  inline unsigned FAlCON::MB_interactions() const
  {
    return
      STATS->CB_direct_iacts()+
      STATS->CC_direct_iacts()+
      STATS->CX_direct_iacts();
  }
  //----------------------------------------------------------------------------
  inline unsigned FAlCON::CB_interactions() const
  {
    return STATS->CB_approx_iacts();
  }
  //----------------------------------------------------------------------------
  inline unsigned FAlCON::CC_interactions() const
  {
    return STATS->CC_approx_iacts();
  }
  //----------------------------------------------------------------------------
  inline unsigned FAlCON::total_interactions () const
  {
    return 
      BB_interactions()+
      MB_interactions()+
      CB_interactions()+
      CC_interactions();
  }
  //----------------------------------------------------------------------------
  inline vect const&FAlCON::root_center() const
  {
    return TREE->root_center();
  }
  //----------------------------------------------------------------------------
  inline real const&FAlCON::root_radius() const
  {
    return radius(TREE->root());
  }
  //----------------------------------------------------------------------------
  inline unsigned const&FAlCON::root_depth () const
  {
    return TREE->depth();
  }
  //----------------------------------------------------------------------------
  inline int const&FAlCON::root_number() const
  {
    return number(TREE->root());
  }
  //----------------------------------------------------------------------------
  inline real const&FAlCON::root_mass  () const
  {
    return mass(GRAV->root());
  }
  //---------------------------------------------------------------------------+
  //                                                                           |
  // ESTIMATION OF DENSITIES                                                   |
  //                                                                           |
  //---------------------------------------------------------------------------+
  inline void FAlCON::estimate_n(unsigned const&Nx, bool const&all)
  {
    if(BODIES && !BODIES->have(fieldbit::r))
      error("[FAlCON::estimate_nd()]: nobody has memory for rho");
    GRAV->estimate_nd(all,Nx);
  }
  //----------------------------------------------------------------------------
  inline void FAlCON::estimate_sd(unsigned const&Nx, bool const&all)
  {
    if(BODIES && !BODIES->have(fieldbit::r))
      error("[FAlCON::estimate_sd()]: nobody has memory for rho");
    GRAV->estimate_sd(all,Nx);
  }
  //----------------------------------------------------------------------------
  inline void FAlCON::estimate_rho(unsigned const&Nx, bool const&all)
  {
    if(BODIES && !BODIES->have(fieldbit::r))
      error("[FAlCON::estimate_md()]: nobody has memory for rho");
    GRAV->estimate_md(all,Nx);
  }
  //---------------------------------------------------------------------------+
  //                                                                           |
  // APPROXIMATION OF GRAVITY                                                  |
  //                                                                           |
  //---------------------------------------------------------------------------+
  inline void FAlCON::approximate_gravity(bool     const&split,
					  bool     const&all
#ifdef falcON_ADAP
					 ,real     const&Nsoft,
					  unsigned const&Nref,
					  real     const&emin,
					  real     const&efac
#  define ADAP_ARGS ,Nsoft,Nref,emin,efac
#else
#  define ADAP_ARGS
#endif
					  )
  {
    GRAV->approx(GMAC,all,split ADAP_ARGS);
  }
  //----------------------------------------------------------------------------
  inline void FAlCON::exact_gravity(bool     const&all
#ifdef falcON_ADAP
				   ,real     const&Nsoft,
				    unsigned const&Nref,
				    real     const&emin,
				    real     const&efac
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
  inline void FAlCON::adjust_SPH_sizes(real     const&mu,
				       real     const&hm,
				       real     const&dm,
				       bool     const&al,
				       unsigned const&ix)
  {
    SPHT->adjust_sizes(mu,hm,dm,al,ix);
  }
  //----------------------------------------------------------------------------
  inline int FAlCON::SPH_sweep_one(real const&mu,
				   real const&dm,
				   real const&hm,
				   real const&wf,
				   bool const&al)
  
  {
    return SPHT->sweep_one(mu,dm,hm,wf,al);
  }
  //----------------------------------------------------------------------------
  inline void FAlCON::SPH_sweep_two(const EquationOfState*const&e,
				    real                  const&a)
  
  {
    SPHT->sweep_two(e,a);
  }
  //----------------------------------------------------------------------------
  inline unsigned const&FAlCON::N_MuSmall() const {
    return SPHT->N_MuSmall();
  }
  //----------------------------------------------------------------------------
  inline unsigned const&FAlCON::N_MuLarge() const {
    return SPHT->N_MuLarge();
  }
  //----------------------------------------------------------------------------
  inline unsigned const&FAlCON::N_HatMax () const {
    return SPHT->N_HatMax ();
  }
  //----------------------------------------------------------------------------
  inline unsigned FAlCON::N_SPH_active(bool all) const {
    if(BODIES->N_bodies(bodytype::gas) == 0) 
      return 0;
    if(all) 
      return BODIES->N_bodies(bodytype::gas);
    unsigned n=0;
    LoopSPHBodies(BODIES,B) if(is_active(B)) ++n;
    return n;
  }
#endif
  //---------------------------------------------------------------------------+
  //                                                                           |
  // NEIGHBOUR FINDING AND COUNTING                                            |
  //                                                                           |
  //---------------------------------------------------------------------------+
  inline void FAlCON::make_iaction_list(indx_pair *bl,
                                        uint       nl,
                                        uint      &na,
                                        bool       Max,
                                        real       tau,
					bool       count)
  {
    if(PAES==0) PAES = new PartnerEstimator(TREE);
    if(tau < zero) PAES->make_sph_list   (bl,nl,na,Max,count);
    else           PAES->make_sticky_list(bl,nl,na,tau,count);
  }
  //---------------------------------------------------------------------------+
  inline void FAlCON::count_sph_partners(bool Max)
  {
    if(PAES==0) PAES = new PartnerEstimator(TREE);
    PAES->count_sph_partners(Max);
  }
  //////////////////////////////////////////////////////////////////////////////
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_FAlCON_cc
