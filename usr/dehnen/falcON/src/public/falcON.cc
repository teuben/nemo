//-----------------------------------------------------------------------------+
//                                                                             |
// falcON.cc                                                                   |
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
#include <body.h>                          // bodies                            
#include <public/grav.h>                   // gravity stuff                     
#include <public/iact.h>                   // generic interaction algorithm     
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
#ifdef ALLOW_INDI
	       const soft_type s,
#endif
	       const MAC_type  mac) :
  STATE     ( unset ),
  BODIES    ( b ),
  TREE      ( 0 ),
  Nb        ( BODIES->N_bodies() ),
  Ncrit     ( 0 ),
  Ncoeffs   ( 0 ),
  FL        ( 0 ),
  M         ( 0 ),
#ifdef ALLOW_INDI
  EP        ( 0 ),
#endif
  X         ( 0 ),
  PH        ( 0 ),
  RH        ( 0 ),
  A         ( 0 ),
  EPS       ( e ),
  GMAC      ( new grav_mac(mac, abs(th)) ),
  STATS     ( new grav_stat() ),
#ifdef ALLOW_INDI
  SOFTENING ( s ),
#endif
  KERNEL    ( k )
{
  MemoryCheck(GMAC);
  MemoryCheck(STATS);
#ifdef ALLOW_INDI
  if(BODIES->has_eps() && SOFTENING==global)
    warning("eps_i given but global softening chosen");
  if(!BODIES->has_eps() && SOFTENING==individual)
    error("no eps_i given, but individual softening chosen");
#endif
}
//------------------------------------------------------------------------------
falcON::falcON(const int      *fl,
	       const areal    *m,
	       const areal    *x[NDIM],
#ifdef ALLOW_INDI
	             areal    *ep,
#endif
	             areal    *a[NDIM],
	             areal    *ph,
	             areal    *rh,
	       const unsigned  n,
	       const real      e,
	       const real      th,
	       const kern_type k,
#ifdef ALLOW_INDI
	       const soft_type s,
#endif
	       const MAC_type  mac) :
  STATE     ( unset ),
  BODIES    ( 0 ),
  TREE      ( 0 ),
  Nb        ( n ),
  Ncrit     ( 0 ),
  Ncoeffs   ( 0 ),
  FL        ( fl ),
  M         ( m ),
#ifdef ALLOW_INDI
  EP        ( ep ),
#endif
  X         ( x ),
  PH        ( ph ),
  RH        ( rh ),
  A         ( a ),
  EPS       ( e ),
  GMAC      ( new grav_mac(mac, abs(th)) ),
  STATS     ( new grav_stat() ),
#ifdef ALLOW_INDI
  SOFTENING ( s ),
#endif
  KERNEL    ( k )
{
  MemoryCheck(GMAC);
  MemoryCheck(STATS);
#ifdef ALLOW_INDI
  if(EP!=0 && SOFTENING==global)
    warning("eps_i given but global softening chosen");
  if(EP==0 && SOFTENING==individual)
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
void falcON::grow(const int Ncr)
  // This routines grows a fresh tree & evaluates the cells' basic source       
  // properties.                                                                
{
  Ncrit = max(Ncr,1);
  if(TREE)
    TREE->rebuild(Ncrit);
  else {
    TREE = BODIES?
#ifdef ALLOW_INDI
      new tree_type(BODIES,      SOFTENING,Ncrit) :
      new tree_type(FL,X,M,EP,Nb,SOFTENING,Ncrit) ;
#else
      new tree_type(BODIES,   Ncrit) :
      new tree_type(FL,X,M,Nb,Ncrit) ;
#endif
    MemoryCheck(TREE);
  }
  STATE=built;                                     // set tree state            
}
//------------------------------------------------------------------------------
void falcON::re_grow(const int Ncut, const int Ncr)
  // This routines grows a fresh tree & evaluates the cells' basic source       
  // properties.                                                                
{
  Ncrit = max(Ncr,1);
  if(TREE)
    TREE->rebuild(Ncrit,Ncut);
  else {
    TREE = BODIES?
#ifdef ALLOW_INDI
      new tree_type(BODIES,      SOFTENING,Ncrit) :
      new tree_type(FL,X,M,EP,Nb,SOFTENING,Ncrit) ;
#else
      new tree_type(BODIES,   Ncrit) :
      new tree_type(FL,X,M,Nb,Ncrit) ;
#endif
    MemoryCheck(TREE);
  }
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
    LoopSouls(tree_type,TREE,Si) if(!is_in_tree(FL[mybody(Si)])) n++;
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
  std::ofstream *oc = fcells? new std::ofstream(fcells) : 0; MemoryCheck(oc);
  std::ofstream *os = fleafs? new std::ofstream(fleafs) : 0; MemoryCheck(os);
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
#ifdef ALLOW_INDI
    if(SOFTENING == global) 
      out<<" softening:             global\n"
	 <<" softening length:      "<< EPS                 <<"\n";
    else
      out<<" softening:             individual\n";
#else
    out<<" softening length:      "<< EPS                   <<"\n";
#endif
    out<<" softening kernel:      "<<describe_kernel()      <<"\n"
       <<" Taylor coeffs used:    "<<Ncoeffs                <<"\n";
  }
  STATS->write(out);
}
//------------------------------------------------------------------------------
nbdy::uint falcON::BB_interactions() const {
  return STATS->BB_iacts();
}
//------------------------------------------------------------------------------
nbdy::uint falcON::MB_interactions() const {
  return STATS->CB_direct_iacts() +STATS->CC_direct_iacts()
        +STATS->CS_direct_iacts();
}
//------------------------------------------------------------------------------
nbdy::uint falcON::CB_interactions() const {
  return STATS->CB_taylor_iacts();
}
//------------------------------------------------------------------------------
nbdy::uint falcON::CC_interactions() const {
  return STATS->CC_taylor_iacts();
}
//------------------------------------------------------------------------------
nbdy::uint falcON::total_interactions () const {
  return BB_interactions()+MB_interactions()
        +CB_interactions()+CC_interactions();
}
//------------------------------------------------------------------------------
const nbdy::vect& falcON::root_center() const { return center(TREE->root()); }
const nbdy::real& falcON::root_radius() const { return radius(TREE->root()); }
const int&        falcON::root_depth () const { return TREE->depth(); }
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
#if NDIM==3
    if(number(Ci)>Nx || n==zero) n = number(Ci)/(NSUB*cube(radius(Ci)));
#else  // NDIM==3
    if(number(Ci)>Nx || n==zero) n = number(Ci)/(NSUB*square(radius(Ci)));
#endif // NDIM==3
    LoopCellKids(cell_iter,Ci,c) if(is_sink(c)) guess_N(c,n,Nx);
    LoopSoulKids(cell_iter,Ci,s) if(is_sink(s)) s->rho() = n;
  }
  //----------------------------------------------------------------------------
  inline void guess_sd(cell_iter const&Ci, real sd, int const&Nx)
  {
    if(number(Ci)>Nx || sd==zero) sd = mass(Ci) / (4*square(radius(Ci)));
    LoopCellKids(cell_iter,Ci,c) if(is_sink(c)) guess_sd(c,sd,Nx);
    LoopSoulKids(cell_iter,Ci,s) if(is_sink(s)) s->rho() = sd;
  }
  //----------------------------------------------------------------------------
  inline void guess_rho(cell_iter const&Ci, real rh, int const&Nx)
  {
#if NDIM==3
    if(number(Ci)>Nx || rh==zero) rh = mass(Ci)/(NSUB*cube(radius(Ci)));
#else  // NDIM==3
    if(number(Ci)>Nx || rh==zero) rh = mass(Ci)/(NSUB*square(radius(Ci)));
#endif // NDIM==3
    LoopCellKids(cell_iter,Ci,c) if(is_sink(c)) guess_rho(c,rh,Nx);
    LoopSoulKids(cell_iter,Ci,s) if(is_sink(s)) s->rho() = rh;
  }
}
//------------------------------------------------------------------------------
void falcON::estimate_n(const int Nx)
{
  if(BODIES && !BODIES->has_rho())
    error("[falcON::estimate_n()]: nobody has memory for rho");
  if(!BODIES && RH==0)
    error("[falcON::estimate_n()]: no array for density given");
  TREE->prepare_density();                         // update flags, set memory  
  if(TREE->N_soul_sinks() == 0) {                  // IF(no sinks)             >
    warning("[falcON::estimate_n()]: nobody wants updating");
    return;                                        //   we are done already     
  }
  guess_N(TREE->root(),zero,Nx);                   // guess n_i                 
  if(BODIES) {
    LoopSouls(tree_type,TREE,Si)
      if(is_sink(Si)) Si->update_dens(BODIES);     // copy n_i into body        
  } else {
    LoopSouls(tree_type,TREE,Si)
      if(is_sink(Si)) Si->update_dens(RH);         // copy n_i into array RH    
  }
}
//------------------------------------------------------------------------------
void falcON::estimate_sd(const int Nx)
{
  if(BODIES && !BODIES->has_rho())
    error("[falcON::estimate_sd()]: nobody has memory for rho");
  if(!BODIES && RH==0)
    error("[falcON::estimate_sd()]: no array for density given");
  TREE->prepare_density();                         // update flags, set memory  
  if(TREE->N_soul_sinks() == 0) {                  // IF(no sinks)             >
    warning("[falcON::estimate_sd()]: nobody wants updating");
    return;                                        //   we are done already     
  }
  guess_sd(TREE->root(),zero,Nx);                  // guess SD_i                
  if(BODIES) {
    LoopSouls(tree_type,TREE,Si)
      if(is_sink(Si)) Si->update_dens(BODIES);     // copy n_i into body        
  } else {
    LoopSouls(tree_type,TREE,Si)
      if(is_sink(Si)) Si->update_dens(RH);         // copy n_i into array RH    
  }
}
//------------------------------------------------------------------------------
void falcON::estimate_rho(const int Nx)
{
  if(BODIES && !BODIES->has_rho())
    error("[falcON::estimate_rho()]: nobody has memory for rho");
  if(!BODIES && RH==0)
    error("[falcON::estimate_rho()]: no array for density given");
  TREE->prepare_density();                         // update flags, set memory  
  if(TREE->N_soul_sinks() == 0) {                  // IF(no sinks)             >
    warning("[falcON::estimate_rho()]: no body wants updating");
    return;                                        //   we are done already     
  }
  guess_rho(TREE->root(),zero,Nx);                 // guess n_i                 
  if(BODIES) {
    LoopSouls(tree_type,TREE,Si)
      if(is_sink(Si)) Si->update_dens(BODIES);     // copy n_i into body        
  } else {
    LoopSouls(tree_type,TREE,Si)
      if(is_sink(Si)) Si->update_dens(RH);         // copy n_i into array RH    
  }
}
//-----------------------------------------------------------------------------+
//                                                                             |
// 5 APPROXIMATION OF GRAVITY                                                  |
// ==========================                                                  |
//                                                                             |
//-----------------------------------------------------------------------------+
#define ROOT TREE->root()
void falcON::approximate_gravity(bool      split,
#ifdef ALLOW_INDI
				 real      Nsoft,
				 uint      Nref,
				 real      efac,
#endif
				 const int direct[4])
{
  if(split) {                                      // IF(interweaving)         >
#ifdef ALLOW_INDI
    TREE->prepare_grav_approx(GMAC,0,Nsoft,EPS,Nref,efac); // prepare for grav  
#else
    TREE->prepare_grav_approx(GMAC,0);             // prepare for grav          
#endif
    if(TREE->N_cell_sinks()==0) {                  // IF(no sink cells)        >
      warning("[falcON::approximate_gravity()]: nobody wants updating");
      return;                                      //   we are done already     
    }                                              // <                         
    STATS->reset();                                // reset iaction statistics  
    register uint NP = split? 4+TREE->N_cell_sinks()/8 : TREE->N_cell_sinks();
                                                   // initial size of C_i pool  
#ifdef ALLOW_INDI
    grav_iact_s GK(STATS,EPS,NP,KERNEL,SOFTENING,direct,false);
#else
    grav_iact_s GK(STATS,EPS,NP,KERNEL,direct,false);
#endif
	                                           //   init gravity kernel     
    MutualInteractor<grav_iact_s> MI(&GK,TREE->depth()-1);
                                                   //   init mutual interactor  
    LoopCellKids(cell_iter,ROOT,c1) {              //   loop cell kids c1      >
      MI.cell_self(c1);                            //     self-interaction c1   
      LoopCellSecd(cell_iter,ROOT,c1+1,c2)         //     loop cell kids c2 > c1
	MI.cell_cell(c1,c2);                       //       interaction c1,c2   
      LoopSoulKids(cell_iter,ROOT,s2)              //     loop soul kids s      
	MI.cell_soul(c1,s2);                       //       interaction c1,s    
      GK.evaluate(c1);                             //     evaluation phase      
    }                                              //   <                       
    LoopSoulKids(cell_iter,ROOT,s1) {              //   loop soul kids s1      >
      LoopSoulSecd(cell_iter,ROOT,s1+1,s2)         //     loop soul kids s2 > s1
	GK.interact(s1,s2);                        //       interaction phase   
      s1->normalize_grav();                        //     evaluation phase      
    }                                              //   <                       
    Ncoeffs = GK.coeffs_used();                    //   remember # coeffs used  
  // case 2.2: non-interweaving, as in the JCP paper                            
  } else {                                         // < ELSE(not splitting)    >
#ifdef _OPENMP
#ifdef ALLOW_INDI
    TREE->prepare_grav_approx(GMAC,1,Nsoft,EPS,Nref,efac); // prepare for grav  
#else
    TREE->prepare_grav_approx(GMAC,1,);            // prepare for grav          
#endif
    if(TREE->N_cell_sinks()==0) {                  // IF(no sink cells)        >
      warning("[falcON::approximate_gravity()]: nobody wants updating");
      return;                                      //   we are done already     
    }                                              // <                         
    STATS->reset();                                // reset iaction statistics  
#ifdef ALLOW_INDI
    grav_iact GK(STATS,EPS,KERNEL,SOFTENING,direct,false);
#else
    grav_iact GK(STATS,EPS,KERNEL,direct,false);   //   init gravity kernel     
#endif
    SelfInteractorP<grav_iact> MI(&GK,ROOT,TREE->depth());
    MI.interact_parallel();
    GK.evaluate(ROOT);                             //   evaluation phase        
    Ncoeffs = TREE->N_cell_sinks();                //   # coeffs used           
#else
#ifdef ALLOW_INDI
    TREE->prepare_grav_approx(GMAC,0,Nsoft,EPS,Nref,efac); // prepare for grav  
#else
    TREE->prepare_grav_approx(GMAC,0);             // prepare for grav approx   
#endif
    if(TREE->N_cell_sinks()==0) {                  // IF(no sink cells)        >
      warning("[falcON::approximate_gravity()]: nobody wants updating");
      return;                                      //   we are done already     
    }                                              // <                         
    STATS->reset();                                // reset iaction statistics  
    register uint NP = split? 4+TREE->N_cell_sinks()/8 : TREE->N_cell_sinks();
                                                   // initial size of C_i pool  
#ifdef ALLOW_INDI
    grav_iact_s GK(STATS,EPS,NP,KERNEL,SOFTENING,direct,false);
#else
    grav_iact_s GK(STATS,EPS,NP,KERNEL,direct,false); //init gravity kernel     
#endif
    MutualInteractor<grav_iact_s> MI(&GK,TREE->depth());
                                                   //   init mutual interactor  
    MI.cell_self(ROOT);                            //   interaction phase       
    GK.evaluate(ROOT);                             //   evaluation phase        
    Ncoeffs = GK.coeffs_used();                    //   remember # coeffs used  
#endif
  }                                                // <                         
  // 3. copy gravity and eps_i from the souls into the bodies                   
  if(BODIES) {
#ifdef ALLOW_INDI
    if(Nsoft) {
      LoopSouls(tree_type,TREE,Si) if(is_sink(Si)) {
	Si->update_grav(BODIES);
	Si->update_eps (BODIES);
      }
    } else
#endif
      LoopSouls(tree_type,TREE,Si) if(is_sink(Si)) Si->update_grav(BODIES);
  } else {
#ifdef ALLOW_INDI
    if(Nsoft) {
      LoopSouls(tree_type,TREE,Si) if(is_sink(Si)) {
	Si->update_grav(A,PH);
	Si->update_eps (EP);
      }
    } else
#endif
      LoopSouls(tree_type,TREE,Si) if(is_sink(Si)) Si->update_grav(A,PH);
  }
}
#undef ROOT
//------------------------------------------------------------------------------
void falcON::exact_gravity(
#ifdef ALLOW_INDI
			   real Nsoft, uint Nref, real efac
#endif
			   )
{
#ifdef ALLOW_INDI
  TREE->prepare_grav_exact(Nsoft,EPS,Nref,efac);
#else
  TREE->prepare_grav_exact();
#endif
  if(TREE->N_cell_sinks()==0) {
    warning("[falcON::exact_gravity()]: nobody wants updating");
    return;
  }
  STATS->reset();
#ifdef ALLOW_INDI
  grav_iact K(STATS,EPS,KERNEL,SOFTENING);
#else
  grav_iact K(STATS,EPS,KERNEL);
#endif
  K.direct_summation(TREE->root());
  if(BODIES) {
#ifdef ALLOW_INDI
    if(Nsoft) {
      LoopSouls(tree_type,TREE,Si) if(is_sink(Si)) {
	Si->normalize_grav();
	Si->update_grav(BODIES);
	Si->update_eps (BODIES);
      }
    } else 
#endif
      LoopSouls(tree_type,TREE,Si) if(is_sink(Si)) {
	Si->normalize_grav();
	Si->update_grav(BODIES);
      }
  } else {
#ifdef ALLOW_INDI
    if(Nsoft) {
      LoopSouls(tree_type,TREE,Si) if(is_sink(Si)) {
	Si->normalize_grav();
	Si->update_grav(A,PH);
	Si->update_eps (EP);
      }
    } else 
#endif
      LoopSouls(tree_type,TREE,Si) if(is_sink(Si)) {
	Si->normalize_grav();
	Si->update_grav(A,PH);
      }
  }
}
//-----------------------------------------------------------------------------+
//                                                                             |
// 6 NEIGHBOUR FINDING AND COUNTING                                            |
// ================================                                            |
//                                                                             |
//-----------------------------------------------------------------------------+
void falcON::make_iaction_list(      elem_pair*bl,
			       const uint      nl,
			             uint     &na,
			       const real      tau)
{
  if(BODIES==0) error("[falcON::make_iaction_list()]: bodies/arrays mismatch");
  if(tau < zero) {             // SPH
    LoopSouls(tree_type,TREE,Si) if(is_sph(BODIES->flg(mybody(Si)))) {
      flag_for_subtree(Si);
      Si->add_body_flag(BODIES);
    }
    sticky_tree      stree(TREE);
    neighbour_finder sfind(nl,bl);
    stree.prepare_sph();
    MutualInteractor<neighbour_finder> MI(&sfind,stree.depth());
    MI.cell_self(stree.root());
    na = sfind.actual_size_of_list();
  } else {                     // sticky particles
    LoopSouls(tree_type,TREE,Si) if(is_sticky(BODIES->flg(mybody(Si)))) {
      flag_for_subtree(Si);
      Si->add_body_flag(BODIES);
    }
    sticky_tree   stree(TREE);
    sticky_finder sfind(tau,nl,bl);
    stree.prepare_sticky();
    MutualInteractor<sticky_finder> MI(&sfind,stree.depth());
    MI.cell_self(stree.root());
    na = sfind.actual_size_of_list();
  }
}
//-----------------------------------------------------------------------------+
void falcON::count_sph_neighbours()
{
  if(BODIES==0) error("[falcON::count_sph_neighbours()]: bodies/arrays mismatch");
  LoopSouls(tree_type,TREE,Si) if(is_sph(BODIES->flg(mybody(Si)))) {
    flag_for_subtree(Si);
    Si->add_body_flag(BODIES);
  }
  sticky_tree       stree (TREE);
  neighbour_counter<sticky_tree,individual> scount;
  stree.prepare_sph();
  MutualInteractor<neighbour_counter<sticky_tree,individual> > 
    MI(&scount,stree.depth());
  MI.cell_self(stree.root());
  stree.set_num(TREE->begin_souls());
}
//-----------------------------------------------------------------------------+
void falcON::make_iaction_list(      elem_pair*bl,
			       const uint      nl,
			             uint     &na,
			       const areal    *s,
			       const areal    *v[NDIM],
			       const real      tau)
{
  if(BODIES) error("[falcON::make_iaction_list()]: bodies/arrays mismatch");
  if(tau < zero) {             // SPH
    LoopSouls(tree_type,TREE,Si) if(is_sph(FL[mybody(Si)])) {
      flag_for_subtree(Si);
      Si->add_body_flag(FL);
    }
    sticky_tree      stree(TREE,s,v);
    neighbour_finder sfind(nl,bl);
    stree.prepare_sph();
    MutualInteractor<neighbour_finder> MI(&sfind,stree.depth());
    MI.cell_self(stree.root());
    na = sfind.actual_size_of_list();
  } else {                     // sticky particles
    LoopSouls(tree_type,TREE,Si) if(is_sticky(FL[mybody(Si)])) {
      flag_for_subtree(Si);
      Si->add_body_flag(FL);
    }
    sticky_tree   stree(TREE,s,v);
    sticky_finder sfind(tau,nl,bl);
    stree.prepare_sticky();
    MutualInteractor<sticky_finder> MI(&sfind,stree.depth());
    MI.cell_self(stree.root());
    na = sfind.actual_size_of_list();
  }
}
//-----------------------------------------------------------------------------+
void falcON::count_sph_neighbours(const areal* s, int* NUM)
{
  if(BODIES) error("[falcON::sph_neighbours()]: bodies/arrays mismatch");
  LoopSouls(tree_type,TREE,Si) if(is_sph(FL[mybody(Si)])) {
    flag_for_subtree(Si);
    Si->add_body_flag(FL);
  }
  sticky_tree       stree(TREE,s,0);
  neighbour_counter<sticky_tree,individual> scount;
  stree.prepare_sph();
  MutualInteractor<neighbour_counter<sticky_tree,individual> >
    MI(&scount,stree.depth());
  MI.cell_self(stree.root());
  stree.set_num(NUM);
}
//-----------------------------------------------------------------------------+
void falcON::count_neighbours() {
  if(!BODIES) error("[falcON::count_neighbours()]: bodies/arrays mismatch");
#ifdef ALLOW_INDI
  switch(SOFTENING) {
  case individual: {
    TREE->prepare_neighbour_counting();
    if(TREE->N_soul_sinks()==0) {
      warning("[falcON::count_neighbours()]: nobody wants updating");
      return;
    }
    neighbour_counter<grav_tree,individual> count;
    MutualInteractor<neighbour_counter<grav_tree,individual> > 
    MI(&count,TREE->depth());
    MI.cell_self(TREE->root());
  } break;
  case global: {
#endif
    TREE->prepare_neighbour_counting(&EPS);
    if(TREE->N_soul_sinks()==0) {
      warning("[falcON::count_neighbours()]: nobody wants updating");
      return;
    }
    neighbour_counter<grav_tree,global> count(EPS);
    MutualInteractor<neighbour_counter<grav_tree,global> > 
    MI(&count,TREE->depth());
    MI.cell_self(TREE->root());
#ifdef ALLOW_INDI
  } break;
  }
#endif
  LoopSouls(tree_type,TREE,Si) if(is_sink(Si)) Si->update_num(BODIES);
}
//-----------------------------------------------------------------------------+
void falcON::count_neighbours(int* NUM) {
  if(BODIES) error("[falcON::count_neighbours()]: bodies/arrays mismatch");
#ifdef ALLOW_INDI
  switch(SOFTENING) {
  case individual: {
    TREE->prepare_neighbour_counting();
    if(TREE->N_soul_sinks()==0) {
      warning("[falcON::count_neighbours()]: nobody wants updating");
      return;
    }
    neighbour_counter<grav_tree,individual> count;
    MutualInteractor<neighbour_counter<grav_tree,individual> > 
    MI(&count,TREE->depth());
    MI.cell_self(TREE->root());
  } break;
  case global: {
#endif
    TREE->prepare_neighbour_counting(&EPS);
    if(TREE->N_soul_sinks()==0) {
      warning("[falcON::count_neighbours()]: nobody wants updating");
      return;
    }
    neighbour_counter<grav_tree,global> count(EPS);
    MutualInteractor<neighbour_counter<grav_tree,global> > 
    MI(&count,TREE->depth());
    MI.cell_self(TREE->root());
#ifdef ALLOW_INDI
  } break;
  }
#endif
  LoopSouls(tree_type,TREE,Si) if(is_sink(Si)) Si->update_num(NUM);
}
////////////////////////////////////////////////////////////////////////////////
