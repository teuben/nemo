//-----------------------------------------------------------------------------+
//                                                                             |
// falcON.cc                                                                   |
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
// defines all function members of class falcON                                |
//                                                                             |
// Contents                                                                    |
//                                                                             |
// 1  Constructor, destructor and such                                         |
// 2  tree growth and re-use                                                   |
// 3  information and statistics                                               |
// 4  estimation of density and such                                           |
// 5  approximation of gravity                                                 |
// 6  neighbour finding and counting                                           |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <falcON.h>                        // class falcON                      
#include <iomanip>                         // C++ file I/O                      
#include <fstream>                         // C++ file I/O                      
#include <public/grav.h>                   // gravity stuff                     
#include <public/stic.h>                   // support for SPH & sticky particles

using namespace nbdy;                      // make all nbdy stuff available     
namespace nbdy {
  typedef grav_tree                tree_type;
  typedef tree_type::cell_iterator cell_iter;
  typedef tree_type::soul_iterator soul_iter;
}
//=============================================================================+
//                                                                             |
// 1 CONSTRUCTORS, DESTRUCTOR AND SUCH                                         |
// ===================================                                         |
//                                                                             |
//=============================================================================+
falcON::falcON(const bodies*   b,
	       const real      e,
	       const real      th,
	       const kern_type k,
#ifdef falcON_INDI
	       const soft_type s,
#endif
	       const MAC_type  mac) :
  STATE     ( unset ),
  BODIES    ( b ),
  ARRAYS    ( 0 ),
  TREE      ( 0 ),
  Ncrit     ( 0 ),
  EPS       ( e ),
  GMAC      ( falcON_Memory(new grav_mac(mac, abs(th))) ),
  STATS     ( falcON_Memory(new grav_stat()) ),
#ifdef falcON_INDI
  SOFTENING ( s ),
#endif
  KERNEL    ( k )
{
#ifdef falcON_INDI
  if(BODIES->has(io::e) && SOFTENING==global)
    warning("eps_i given but global softening chosen");
  if(!BODIES->has(io::e) && SOFTENING==individual)
    error("no eps_i given, but individual softening chosen");
#endif
}
//------------------------------------------------------------------------------
falcON::falcON(const barrays  *b,
	       const real      e,
	       const real      th,
	       const kern_type k,
#ifdef falcON_INDI
	       const soft_type s,
#endif
	       const MAC_type  mac) :
  STATE     ( unset ),
  BODIES    ( 0 ),
  ARRAYS    ( b ),
  TREE      ( 0 ),
  Ncrit     ( 0 ),
  EPS       ( e ),
  GMAC      ( falcON_Memory(new grav_mac(mac, abs(th))) ),
  STATS     ( falcON_Memory(new grav_stat()) ),
#ifdef falcON_INDI
  SOFTENING ( s ),
#endif
  KERNEL    ( k )
{
#ifdef falcON_INDI
  if(ARRAYS->has(io::e) && SOFTENING==global)
    warning("eps_i given but global softening chosen");
  if(!ARRAYS->has(io::e) && SOFTENING==individual)
    error("no eps_i given, but individual softening chosen");
#endif
}
//------------------------------------------------------------------------------
falcON::~falcON() {
  if(TREE) delete TREE;
  delete GMAC;
  delete STATS;
}
//------------------------------------------------------------------------------
void falcON::reset_softening (const real      e,
			      const kern_type t) const
{
  EPS    = e;
  KERNEL = t;
}
//------------------------------------------------------------------------------
void falcON::reset_opening(const real     th,
			   const MAC_type mac) const
{
  GMAC->reset(mac, abs(th));
}
//=============================================================================+
//                                                                             |
// 2 TREE GROWTH AND RE-USE                                                    |
// ========================                                                    |
//                                                                             |
//=============================================================================+
void falcON::grow(int const&Ncr)
  // This routines grows a fresh tree & evaluates the cells' basic source       
  // properties.                                                                
{
  Ncrit = max(Ncr,1);
  if(TREE)
    TREE->rebuild(Ncrit);
  else
    TREE = BODIES?
#ifdef falcON_INDI
      falcON_Memory(new tree_type(BODIES,SOFTENING,Ncrit)) :
      falcON_Memory(new tree_type(ARRAYS,SOFTENING,Ncrit)) ;
#else
      falcON_Memory(new tree_type(BODIES,Ncrit)) :
      falcON_Memory(new tree_type(ARRAYS,Ncrit)) ;
#endif
  STATE=built;                                     // set tree state            
}
//------------------------------------------------------------------------------
void falcON::re_grow(int const&Ncut, int const&Ncr)
  // This routines grows a fresh tree & evaluates the cells' basic source       
  // properties.                                                                
{
  Ncrit = max(Ncr,1);
  if(TREE)
    TREE->rebuild(Ncrit,Ncut);
  else
    TREE = BODIES?
#ifdef falcON_INDI
      falcON_Memory(new tree_type(BODIES,SOFTENING,Ncrit)) :
      falcON_Memory(new tree_type(ARRAYS,SOFTENING,Ncrit)) ;
#else
      falcON_Memory(new tree_type(BODIES,Ncrit)) :
      falcON_Memory(new tree_type(ARRAYS,Ncrit)) ;
#endif
  STATE=built;                                     // set tree state            
}
//------------------------------------------------------------------------------
void falcON::reuse()
{
  // This routine uses the existing tree and                                    
  // -  re-sets the souls' positions, masses, and source flags                  
  // -  re-evaluates the basic source properties of the cells                   
  // and pre-compute basic source properties of cells                           
  if(TREE == 0) {
    warning("no old tree to be re-used","falcON::reuse()");
    return grow();
  }
  TREE->reuse();
  STATE=built;                                     // set tree state            
}
//=============================================================================+
//                                                                             |
// 3 INFORMATION AND STATISTICS                                                |
// ============================                                                |
//                                                                             |
//=============================================================================+
const MAC_type &falcON::MAC() const {
  return GMAC->method();
}
//------------------------------------------------------------------------------
const char* falcON::describe_MAC() const {
  return GMAC->describe_method();
}
//------------------------------------------------------------------------------
const nbdy::uint &falcON::No_bodies() const {
  return BODIES->N_bodies();
}
//------------------------------------------------------------------------------
nbdy::uint falcON::No_bodies_used() const {
  return TREE->N_souls();
}
//------------------------------------------------------------------------------
nbdy::uint falcON::No_cells_used() const {
  return TREE->N_cells();
}
//------------------------------------------------------------------------------
nbdy::uint falcON::No_zombie_bodies() const {
  register uint  n=0;
  if(BODIES) {
    LoopSouls(tree_type,TREE,Si) if(!is_in_tree(BODIES->flg(mybody(Si)))) n++;
  } else {
    LoopSouls(tree_type,TREE,Si) if(!is_in_tree(ARRAYS->flg(mybody(Si)))) n++;
  }
  return n;
}
//------------------------------------------------------------------------------
namespace nbdy {
  inline int trick(const real x, const int w)
  {
    if(x<0.001)  return max(1,w-5);
    if(x<0.01)   return max(1,w-4);
    if(x<0.1)    return max(1,w-3);
    if(x<one)    return max(1,w-2);
    if(x<ten)    return max(1,w-1);
    return w;
  }
}
//------------------------------------------------------------------------------
void falcON::dump_nodes(const char* fcells, const char* fleafs) const
{
  std::ofstream *oc=0, *os=0;
  if(fcells) oc = falcON_Memory(new std::ofstream(fcells));
  if(fleafs) os = falcON_Memory(new std::ofstream(fleafs));
  TREE->dump_nodes(oc,os);
  if(oc) delete oc;
  if(os) delete os;
}
//------------------------------------------------------------------------------
void falcON::stats(std::ostream& out) const
{
  out<<"\n state:                ";
  if(STATE&unset)  out<<" no tree";
  if(STATE&built)  out<<" tree built";
  if(STATE&reused) out<<" tree reused";
  out<<"\n";
  if(STATE!=unset) {
    out<<" root center:           "<<center(TREE->root())   <<"\n"
       <<" root radius:           "<<radius(TREE->root())   <<"\n"
       <<" bodies loaded:         "<<number(TREE->root())   <<"\n"
       <<" total mass:            "<<mass  (TREE->root())   <<"\n"
       <<" N_crit:                "<<Ncrit                  <<"\n"
       <<" cells used:            "<<TREE->N_cells()        <<"\n"
       <<" maximum depth:         "<<TREE->depth()          <<"\n"
       <<" current theta:         "<<GMAC->theta_min()      <<"\n"
       <<" current MAC:           "<<GMAC->describe_method()<<"\n";
#ifdef falcON_INDI
    if(SOFTENING == global) 
      out<<" softening:             global\n"
	 <<" softening length:      "<< EPS                 <<"\n";
    else
      out<<" softening:             individual\n";
#else
    out<<" softening length:      "<< EPS                   <<"\n";
#endif
    out<<" softening kernel:      "<<describe_kernel()      <<"\n"
       <<" Taylor coeffs used:    "<<TREE->N_coeffs()       <<"\n";
  }
  STATS->write(out);
}
//-----------------------------------------------------------------------------+
nbdy::uint falcON::No_coeffs_used()   const {
  return TREE? TREE->N_coeffs() : 0;
}
//------------------------------------------------------------------------------
nbdy::uint falcON::BB_interactions() const {
  return STATS->BB_direct_iacts();
}
//------------------------------------------------------------------------------
nbdy::uint falcON::MB_interactions() const {
  return STATS->CB_direct_iacts() +STATS->CC_direct_iacts()
        +STATS->CX_direct_iacts();
}
//------------------------------------------------------------------------------
nbdy::uint falcON::CB_interactions() const {
  return STATS->CB_approx_iacts();
}
//------------------------------------------------------------------------------
nbdy::uint falcON::CC_interactions() const {
  return STATS->CC_approx_iacts();
}
//------------------------------------------------------------------------------
nbdy::uint falcON::total_interactions () const {
  return BB_interactions()+MB_interactions()
        +CB_interactions()+CC_interactions();
}
//------------------------------------------------------------------------------
const nbdy::vect& falcON::root_center() const { return center(TREE->root()); }
const nbdy::real& falcON::root_radius() const { return radius(TREE->root()); }
const nbdy::uint& falcON::root_depth () const { return TREE->depth(); }
const int&        falcON::root_number() const { return number(TREE->root()); }
const nbdy::real& falcON::root_mass  () const { return mass  (TREE->root()); }
//-----------------------------------------------------------------------------+
//                                                                             |
// 4 ESTIMATION OF DENSITY AND SUCH                                            |
// ================================                                            |
//                                                                             |
//-----------------------------------------------------------------------------+
namespace nbdy {
  inline void guess_N(cell_iter const&Ci, real n, int const&Nx)
  {
#if falcON_NDIM==3
    if(number(Ci)>Nx || n==zero) n = number(Ci)/(Nsub*cube(radius(Ci)));
#else  // falcON_NDIM==3
    if(number(Ci)>Nx || n==zero) n = number(Ci)/(Nsub*square(radius(Ci)));
#endif // falcON_NDIM==3
    LoopCellKids(cell_iter,Ci,c) if(is_active(c)) guess_N(c,n,Nx);
    LoopSoulKids(cell_iter,Ci,s) if(is_active(s)) s->rho() = n;
  }
  //----------------------------------------------------------------------------
  inline void guess_sd(cell_iter const&Ci, real sd, int const&Nx)
  {
    if(number(Ci)>Nx || sd==zero) sd = mass(Ci) / (4*square(radius(Ci)));
    LoopCellKids(cell_iter,Ci,c) if(is_active(c)) guess_sd(c,sd,Nx);
    LoopSoulKids(cell_iter,Ci,s) if(is_active(s)) s->rho() = sd;
  }
  //----------------------------------------------------------------------------
  inline void guess_rho(cell_iter const&Ci, real rh, int const&Nx)
  {
#if falcON_NDIM==3
    if(number(Ci)>Nx || rh==zero) rh = mass(Ci)/(Nsub*cube(radius(Ci)));
#else  // falcON_NDIM==3
    if(number(Ci)>Nx || rh==zero) rh = mass(Ci)/(Nsub*square(radius(Ci)));
#endif // falcON_NDIM==3
    LoopCellKids(cell_iter,Ci,c) if(is_active(c)) guess_rho(c,rh,Nx);
    LoopSoulKids(cell_iter,Ci,s) if(is_active(s)) s->rho() = rh;
  }
}
//------------------------------------------------------------------------------
void falcON::estimate_n(int const&Nx)
{
  if(BODIES && !BODIES->has(io::r) || ARRAYS && !ARRAYS->has(io::r))
    error("[falcON::estimate_n()]: nobody has memory for rho");
  TREE->prepare_density();                         // update flags, set memory  
  if(TREE->N_active_souls() == 0) {                // IF(no active souls) THEN  
    warning("[falcON::estimate_n()]: nobody wants updating");
    return;                                        //   we are done already     
  }
  guess_N(TREE->root(),zero,Nx);                   // guess n_i                 
  if(BODIES) {
    LoopSouls(tree_type,TREE,Si)
      if(is_active(Si)) Si->update_dens(BODIES);   // copy n_i into body        
  } else {
    LoopSouls(tree_type,TREE,Si)
      if(is_active(Si)) Si->update_dens(ARRAYS);   // copy n_i into array RH    
  }
}
//------------------------------------------------------------------------------
void falcON::estimate_sd(int const&Nx)
{
  if(BODIES && !BODIES->has(io::r) || ARRAYS && !ARRAYS->has(io::r))
    error("[falcON::estimate_sd()]: nobody has memory for rho");
  TREE->prepare_density();                         // update flags, set memory  
  if(TREE->N_active_souls() == 0) {                // IF(no active souls) THEN  
    warning("[falcON::estimate_sd()]: nobody wants updating");
    return;                                        //   we are done already     
  }
  guess_sd(TREE->root(),zero,Nx);                  // guess SD_i                
  if(BODIES) {
    LoopSouls(tree_type,TREE,Si)
      if(is_active(Si)) Si->update_dens(BODIES);   // copy n_i into body        
  } else {
    LoopSouls(tree_type,TREE,Si)
      if(is_active(Si)) Si->update_dens(ARRAYS);   // copy n_i into array RH    
  }
}
//------------------------------------------------------------------------------
void falcON::estimate_rho(int const&Nx)
{
  if(BODIES && !BODIES->has(io::r) || ARRAYS && !ARRAYS->has(io::r))
    error("[falcON::estimate_rho()]: nobody has memory for rho");
  TREE->prepare_density();                         // update flags, set memory  
  if(TREE->N_active_souls() == 0) {                // IF(no active souls) THEN  
    warning("[falcON::estimate_rho()]: no body wants updating");
    return;                                        //   we are done already     
  }
  guess_rho(TREE->root(),zero,Nx);                 // guess n_i                 
  if(BODIES) {
    LoopSouls(tree_type,TREE,Si)
      if(is_active(Si)) Si->update_dens(BODIES);   // copy n_i into body        
  } else {
    LoopSouls(tree_type,TREE,Si)
      if(is_active(Si)) Si->update_dens(ARRAYS);   // copy n_i into array RH    
  }
}
//-----------------------------------------------------------------------------+
//                                                                             |
// 5 APPROXIMATION OF GRAVITY                                                  |
// ==========================                                                  |
//                                                                             |
//-----------------------------------------------------------------------------+
void falcON::approximate_gravity(bool const&split,
				 bool const&all,
#ifdef falcON_INDI
				 real const&Nsoft,
				 uint const&Nref,
				 real const&emin,
				 real const&efac,
#endif
				 const int direct[4])
{
  TREE->approx_gravity(GMAC,KERNEL,STATS,EPS,all,split,
#ifdef falcON_INDI
		       Nsoft,Nref,emin,efac,
#endif
		       direct);
}
//------------------------------------------------------------------------------
void falcON::exact_gravity(bool const&all
#ifdef falcON_INDI
			  ,real const&Nsoft,
			   uint const&Nref,
			   real const&emin,
			   real const&efac
#endif
			   )
{
  TREE->exact_gravity(KERNEL,STATS,EPS,all
#ifdef falcON_INDI
		      ,Nsoft,Nref,emin,efac
#endif
		      );
}
//-----------------------------------------------------------------------------+
//                                                                             |
// 6 NEIGHBOUR FINDING AND COUNTING                                            |
// ================================                                            |
//                                                                             |
//-----------------------------------------------------------------------------+
void falcON::make_iaction_list(elem_pair *bl,
			       uint const&nl,
			       uint      &na,
			       real const&tau) const
{
  if(BODIES) LoopSouls(tree_type,TREE,Si) Si->add_body_flag(BODIES);
  else       LoopSouls(tree_type,TREE,Si) Si->add_body_flag(ARRAYS);
  sticky_tree      stree(TREE, tau<zero? flag::SPH : flag::STICKY);
  stree.make_iaction_list(bl,nl,na,tau);
}
//-----------------------------------------------------------------------------+
void falcON::count_sph_neighbours()
{
  if(BODIES && !BODIES->has(io::n) || ARRAYS && !ARRAYS->has(io::n))
    falcON_ErrorF("nobody has memory for num","falcON::count_sph_neighbours()");
  if(BODIES) LoopSouls(tree_type,TREE,Si) Si->add_body_flag(BODIES);
  else       LoopSouls(tree_type,TREE,Si) Si->add_body_flag(ARRAYS);
  sticky_tree stree (TREE, flag::SPH);
  stree.count_neighbours();
}
//-----------------------------------------------------------------------------+
void falcON::count_neighbours()
{
#ifdef falcON_INDI
  if(SOFTENING == individual)
    return TREE->count_neighbours();
  else
#endif
    return TREE->count_neighbours(EPS);
}
////////////////////////////////////////////////////////////////////////////////
