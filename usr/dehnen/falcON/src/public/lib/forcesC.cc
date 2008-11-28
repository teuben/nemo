//-----------------------------------------------------------------------------+
//                                                                             |
// forcesC.cc                                                                  |
//                                                                             |
// Copyright (C) 2000-2008 Walter Dehnen                                       |
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
#define forcesC_cc
#include <forces_C.h>
#include <forces.h>
#include <body.h>
using namespace falcON;
//=============================================================================#
// auxiliary data                                                              |
//=============================================================================#
namespace {
  struct ebodies : public bodies {
    ebodies(const unsigned n[BT_NUM]) : bodies('e', n) {}
    void reset(falcON::fieldbit f, void *D) { bodies::reset('e',f,D); }
  }     *BODIES = 0;
  forces*FALCON = 0;
  bool   BUILT  = 0;
  //===========================================================================#
  // auxiliary routines                                                        |
  //===========================================================================#
  inline falcON::kern_type ktype(const int K) {
    switch(K%10) {
    case  0: return falcON::p0;
    case  1: return falcON::p1;
    case  2: return falcON::p2;
    case  3: return falcON::p3;
    case  9: return falcON::newton;
    default: falcON_Warning("unknown kernel, "
			    "using Newtonian greens function\n");
      return falcON::newton;
    }
  }
  //----------------------------------------------------------------------------
  inline int ktype(const falcON::kern_type K) {
    switch(K) {
    case falcON::p0:     return 0;
    case falcON::p1:     return 1;
    case falcON::p2:     return 2;
    case falcON::p3:     return 3;
    case falcON::newton: return 9;
    default: falcON_Warning(" unknown kernel type, defaulting to newton\n");
      return 9;
    }
  }
  //----------------------------------------------------------------------------
  inline void __falcON_error(const char* name) {
    if(FALCON==0)
      falcON_Error("%s() called before falcON_initialize()\n",name);
  }
  //----------------------------------------------------------------------------
  inline int __falcON_warning(const char* name) {
    if(FALCON==0) {
      falcON_Warning("%s() called before falcON_initialize()\n",name);
      return 1;
    }
    return 0;
  }
  //----------------------------------------------------------------------------
  inline void __falcON_grown(const char* name) {
    if(!BUILT) {
      falcON_Warning("%s() called before a tree has been grown\n"
		     "      I will grow the tree (via falcON_grow()) first\n",
	      name);
      FALCON->grow();
      BUILT = true;
    }
  }
  //===========================================================================#
  inline void __falcON_initialize(int *F,
				  real*M,
				  real*X,
				  real*E,
				  real*A, 
				  real*P,
				  real*R,
				  int  Ntot,
				  int  Nsph,
				  real EPS,
				  real TH,
				  int  K,
				  real G)
  {
    if(BODIES) falcON_DEL_O(BODIES);
    if(FALCON) falcON_DEL_O(FALCON);
    if(Nsph > Ntot)
      falcON_Error("falcON_initialize(): Ntot (%d) < Nsph (%d)\n", Ntot, Nsph);
    unsigned Nbod[BT_NUM] = {Nsph, Ntot-Nsph};
    BODIES = new ebodies(Nbod);
    BODIES->reset(fieldbit::f,F);
    BODIES->reset(fieldbit::m,M);
    BODIES->reset(fieldbit::x,X);
    BODIES->reset(fieldbit::e,E);
    BODIES->reset(fieldbit::a,A);
    BODIES->reset(fieldbit::p,P);
    BODIES->reset(fieldbit::r,R);
    FALCON = new forces(BODIES,
			abs(real(EPS)),
			abs(real(TH)),
			ktype(K),
			E != 0,
			real(G),
			TH<0? const_theta : theta_of_M);
    BUILT   = false;
  }
  //===========================================================================#
  typedef unsigned elem_pair[2];

#define INPAIR(ELPAIR)							\
  (static_cast<falcON::forces::indx_pair*>(static_cast<void*>(ELPAIR)))

  inline void __falcON_iactionlist(elem_pair *IL,
				   unsigned   NL,
				   int       *NS,
				   const real*H,
				   bool       Max,
				   real       tau,
				   const real*V,
				   int       *N,
				   unsigned   shift)
  {
    __falcON_error("falcon_iactionlist");
    __falcON_grown("falcon_iactionlist");
    if(BODIES->N_sph() == 0)
      falcON_Error("falcON_iactionlist(): no SPH particles registered with "
		   "falcON_initialize(): no action taken\n");
    if(tau>=zero && V==0)
      falcON_Error("falcON_iactionlist(): tau<0 but no velocities given\n");
    BODIES->reset(fieldbit::v,const_cast<real*>(V));
    BODIES->reset(fieldbit::H,const_cast<real*>(H));
    BODIES->reset(fieldbit::N,N);
    unsigned NA;
    FALCON->make_iaction_list(INPAIR(IL),NL,NA,Max,tau, N!=0);
    if(NL < NA) NA = NL;
    for(unsigned i=0; i!=NA; ++i) {
      IL[i][0] = shift + BODIES->bodyindex(INPAIR(IL)[i][0]);
      IL[i][1] = shift + BODIES->bodyindex(INPAIR(IL)[i][1]);
    }
    *NS = NA;
  }
} // namespace {
//=============================================================================#
// routines defined in falcON_C.h and falcon.f                                 |
//=============================================================================#

extern "C" {

  //===========================================================================#
  real falcON_default_theta() {
    return falcON::Default::theta;
  }
  //---------------------------------------------------------------------------+
  int falcON_default_Ncrit() {
    return falcON::Default::Ncrit;
  }
  //---------------------------------------------------------------------------+
  int falcON_default_kernel() {
    return ktype(falcON::Default::kernel);
  }
  //===========================================================================#
  void falcON_initialize(const int *F,
			 const real*M,
			 const real*X,
			 const real*E,
			 real      *A,
			 real      *P,
			 real      *R,
			 int        Ntot,
			 int        Nsph,
			 real       EPS,
			 real       TH,
			 int        K,
			 real       G)
  {
    __falcON_initialize(const_cast<int *>(F),
			const_cast<real*>(M),
			const_cast<real*>(X),
			const_cast<real*>(E),
			A,P,R,Ntot,Nsph,EPS,TH,K,G);
  }
  //----------------------------------------------------------------------------
  void falcon_initialize_ (int *F, real*M, real*X, real*E,
			   real*A, real*P, real*R,
			   int *Nt, int*Ns, real*EPS,
			   real*TH, int*K, real*G)
  {
    __falcON_initialize(F,M,X,E,A,P,R,*Nt,*Ns,*EPS,*TH,*K,*G);
  }
  //----------------------------------------------------------------------------
  void falcon_initialize__(int *F, real*M, real*X, real*E,
			   real*A, real*P, real*R,
			   int *Nt, int*Ns, real*EPS,
			   real*TH, int*K, real*G)
  {
    __falcON_initialize(F,M,X,E,A,P,R,*Nt,*Ns,*EPS,*TH,*K,*G);
  }
  //===========================================================================#
  void falcON_resetsoftening(real EPS, int K)
  {
    if(__falcON_warning("falcON_resetsoftening")) return;
    FALCON->reset_softening(real(EPS),ktype(K));
  }
  //----------------------------------------------------------------------------
  void falcon_resetsoftening_(real*EPS, int*K)
  {
    if(__falcON_warning("falcon_resetsoftening")) return;
    FALCON->reset_softening(real(*EPS),ktype(*K));
  }
  //----------------------------------------------------------------------------
  void falcon_resetsoftening__(real*EPS, int*K)
  {
    if(__falcON_warning("falcon_resetsoftening")) return;
    FALCON->reset_softening(real(*EPS),ktype(*K));
  }
  //===========================================================================#
  void falcON_resetopening(real TH)
  {
    if(__falcON_warning("falcON_resetopening")) return;
    FALCON->reset_opening(real(TH),falcON::theta_of_M);
  }  
  //----------------------------------------------------------------------------
  void falcon_resetopening_(real *TH)
  {
    if(__falcON_warning("falcon_resetopening")) return;
    FALCON->reset_opening(real(*TH),falcON::theta_of_M);
  }  
  //----------------------------------------------------------------------------
  void falcon_resetopening__(real *TH)
  {
    if(__falcON_warning("falcon_resetopening")) return;
    FALCON->reset_opening(real(*TH),falcON::theta_of_M);
  }  
  //===========================================================================#
  void falcON_clearup()
  {
    if(FALCON) falcON_DEL_O(FALCON);
    FALCON = 0;
    if(BODIES) falcON_DEL_O(BODIES);
    BODIES = 0;
    BUILT  = 0;
  }
  //----------------------------------------------------------------------------
  void falcon_clearup_()
  {
    if(FALCON) falcON_DEL_O(FALCON);
    FALCON = 0;
    if(BODIES) falcON_DEL_O(BODIES);
    BODIES = 0;
    BUILT  = 0;
  }
  //----------------------------------------------------------------------------
  void falcon_clearup__()
  {
    if(FALCON) falcON_DEL_O(FALCON);
    FALCON = 0;
    if(BODIES) falcON_DEL_O(BODIES);
    BODIES = 0;
    BUILT  = 0;
  }
  //===========================================================================#
  void falcON_grow(int Nc)
  {
    __falcON_error("falcON_grow");
    FALCON->grow(Nc);
    BUILT = true;
  }
  //----------------------------------------------------------------------------
  void falcon_grow_(int *Nc)
  {
    __falcON_error("falcon_grow");
    FALCON->grow(*Nc);
    BUILT = true;
  }
  //----------------------------------------------------------------------------
  void falcon_grow__(int *Nc)
  {
    __falcON_error("falcon_grow");
    FALCON->grow(*Nc);
    BUILT = true;
  }
  //===========================================================================#
  void falcON_reuse()
  {
    __falcON_error("falcON_reuse");
    if(!BUILT) {
      falcON_Warning(" faclON WARNING: falcON_reuse()"
		     " called before a tree has been grown\n"
		     "   I will grow the tree (via falcON_grow()) instead\n");
      FALCON->grow();
      BUILT = true;
    } else
      FALCON->reuse();
  }
  //----------------------------------------------------------------------------
  void falcon_reuse_()
  {
    __falcON_error("falcon_reuse");
    if(!BUILT) {
      falcON_Warning(" faclON WARNING: falcON_reuse()"
		     " called before a tree has been grown\n"
		     "   I will grow the tree (via falcON_grow()) instead\n");
      FALCON->grow();
      BUILT = true;
    } else
      FALCON->reuse();
  }
  //----------------------------------------------------------------------------
  void falcon_reuse__()
  {
    __falcON_error("falcon_reuse");
    if(!BUILT) {
      falcON_Warning(" faclON WARNING: falcON_reuse()"
		     " called before a tree has been grown\n"
		     "   I will grow the tree (via falcON_grow()) instead\n");
      FALCON->grow();
      BUILT = true;
    } else
      FALCON->reuse();
  }
  //===========================================================================#
  void falcON_approx_grav()
  {
    __falcON_error("falcON_approx_gravity");
    __falcON_grown("falcON_approx_gravity");
    FALCON->approximate_gravity();
  }  
  //----------------------------------------------------------------------------
  void falcon_approx_grav_()
  {
    __falcON_error("falcon_approx_gravity");
    __falcON_grown("falcon_approx_gravity");
    FALCON->approximate_gravity();
  }  
  //----------------------------------------------------------------------------
  void falcon_approx_grav__()
  {
    __falcON_error("falcon_approx_gravity");
    __falcON_grown("falcon_approx_gravity");
    FALCON->approximate_gravity();
  }  
#ifdef falcON_ADAP
  //===========================================================================#
  void falcON_adjust_epsi_and_approx_grav(real Nsoft, int Nref, real fac)
  {
    __falcON_error("falcON_adjust_epsi_and_approx_gravity");
    __falcON_grown("falcON_adjust_epsi_and_approx_gravity");
    FALCON->approximate_gravity(true,false,real(Nsoft),Nref,real(fac));
  }
  //----------------------------------------------------------------------------
  void falcon_adjust_epsi_and_approx_grav_(real*Nsoft, int*Nref, real*fac)
  {
    __falcON_error("falcon_adjust_epsi_and_approx_gravity");
    __falcON_grown("falcon_adjust_epsi_and_approx_gravity");
    FALCON->approximate_gravity(true,false,real(*Nsoft),*Nref,real(*fac));
  }
  //----------------------------------------------------------------------------
  void falcon_adjust_epsi_and_approx_grav__(real*Nsoft, int*Nref, real*fac)
  {
    __falcON_error("falcon_adjust_epsi_and_approx_gravity");
    __falcON_grown("falcon_adjust_epsi_and_approx_gravity");
    FALCON->approximate_gravity(true,false,real(*Nsoft),*Nref,real(*fac));
  }
  //===========================================================================#
#endif
  void falcON_estimate_rho(int nx)
  {
    __falcON_error("falcON_estimate_rho");
    __falcON_grown("falcON_estimate_rho");
    FALCON->estimate_rho(nx);
  }  
  //----------------------------------------------------------------------------
  void falcon_estimate_rho_(int*nx)
  {
    __falcON_error("falcon_estimate_rho");
    __falcON_grown("falcon_estimate_rho");
    FALCON->estimate_rho(*nx);
  }  
  //----------------------------------------------------------------------------
  void falcon_estimate_rho__(int*nx)
  {
    __falcON_error("falcon_estimate_rho");
    __falcON_grown("falcon_estimate_rho");
    FALCON->estimate_rho(*nx);
  }  
  //===========================================================================#
  void falcON_estimate_n(int nx)
  {
    __falcON_error("falcON_estimate_n");
    __falcON_grown("falcON_estimate_n");
    FALCON->estimate_n(nx);
  }  
  //----------------------------------------------------------------------------
  void falcon_estimate_n_(int*nx)
  {
    __falcON_error("falcon_estimate_n");
    __falcON_grown("falcon_estimate_n");
    FALCON->estimate_n(*nx);
  }  
  //----------------------------------------------------------------------------
  void falcon_estimate_n__(int*nx)
  {
    __falcON_error("falcon_estimate_n");
    __falcON_grown("falcon_estimate_n");
    FALCON->estimate_n(*nx);
  }  
  //===========================================================================#
  void falcON_iactionlist(int*il, int nl, int*na,
			  const real*H, int Max, real tau, const real*V,
			  int*N)
  {
    __falcON_iactionlist(static_cast<elem_pair*>(static_cast<void*>(il)),
			 nl, na, H, Max!=0, tau, V, N, 0u);
  }
  //----------------------------------------------------------------------------
  void falcON_ialist_ (int*il, int*nl, int*na,
		       real*H, int*Max, real*tau, real*V,
		       int*N)
  {
    __falcON_iactionlist(static_cast<elem_pair*>(static_cast<void*>(il)),
			 *nl, na, H, *Max!=0, *tau, V, N, 1u);
  }
  //----------------------------------------------------------------------------
  void falcON_ialist__(int*il, int*nl, int*na,
		       real*H, int*Max, real*tau, real*V,
		       int*N)
  {
    __falcON_iactionlist(static_cast<elem_pair*>(static_cast<void*>(il)),
			 *nl, na, H, *Max!=0, *tau, V, N, 1u);
  }
  //===========================================================================#
  void falcON_sph_count(const real*H, int Max, int*N)
  {
    __falcON_error("falcon_sph_count");
    __falcON_grown("falcon_sph_count");
    if(BODIES->N_sph() == 0)
      falcON_Error("falcON_sph_count(): no SPH particles registered with "
		   "falcON_initialize(): no action taken\n");
    BODIES->reset(fieldbit::H,const_cast<real*>(H));
    BODIES->reset(fieldbit::N,N);
    FALCON->count_sph_partners(Max != 0);
  }
  //----------------------------------------------------------------------------
  void falcON_sph_count_(real*H, int*Max, int*N)
  {
    __falcON_error("falcon_sph_count");
    __falcON_grown("falcon_sph_count");
    if(BODIES->N_sph() == 0)
      falcON_Error("falcON_sph_count(): no SPH particles registered with "
		   "falcON_initialize(): no action taken\n");
    BODIES->reset(fieldbit::H,H);
    BODIES->reset(fieldbit::N,N);
    FALCON->count_sph_partners( *Max!=0 );
  }
  //----------------------------------------------------------------------------
  void falcON_sph_count__(real*H, int*Max, int*N)
  {
    __falcON_error("falcon_sph_count");
    __falcON_grown("falcon_sph_count");
    if(BODIES->N_sph() == 0)
      falcON_Error("falcON_sph_count(): no SPH particles registered with "
		   "falcON_initialize(): no action taken\n");
    BODIES->reset(fieldbit::H,H);
    BODIES->reset(fieldbit::N,N);
    FALCON->count_sph_partners( *Max!=0 );
  }
  //===========================================================================#
  real falcON_root_center(int i)
  {
    if(__falcON_warning("falcON_root_center")) return 0.;
    return real((FALCON->root_center())[i]);
  }
  //----------------------------------------------------------------------------
  float falcon_root_center_(int*i)
  {
    if(__falcON_warning("falcon_root_center")) return 0.;
    return float((FALCON->root_center())[*i]);
  }
  //----------------------------------------------------------------------------
  float falcon_root_center__(int*i)
  {
    if(__falcON_warning("falcon_root_center")) return 0.;
    return float((FALCON->root_center())[*i]);
  }
  //===========================================================================#
  real falcON_root_radius()
  {
    if(__falcON_warning("falcON_root_radius")) return 0.;
    return real((FALCON->root_radius()));
  }
  //----------------------------------------------------------------------------
  float falcon_root_radius_()
  {
    if(__falcON_warning("falcon_root_radius")) return 0.;
    return float((FALCON->root_radius()));
  }
  //----------------------------------------------------------------------------
  float falcon_root_radius__()
  {
    if(__falcON_warning("falcon_root_radius")) return 0.;
    return float((FALCON->root_radius()));
  }
  //===========================================================================#
  real falcON_current_eps()
  {
    if(__falcON_warning("falcON_current_eps")) return 0.;
    return real(FALCON->softening_length());
  }
  //----------------------------------------------------------------------------
  float falcon_current_eps_()
  {
    if(__falcON_warning("falcon_current_eps")) return 0.;
    return float(FALCON->softening_length());
  }
  //----------------------------------------------------------------------------
  float falcon_current_eps__()
  {
    if(__falcON_warning("falcon_current_eps")) return 0.;
    return float(FALCON->softening_length());
  }
  //===========================================================================#
  int falcON_current_kernel()
  {
    if(__falcON_warning("falcON_current_kernel")) return 0;
    return ktype(FALCON->kernel());
  }
  //----------------------------------------------------------------------------
  int falcon_current_kernel_()
  {
    if(__falcON_warning("falcon_current_kernel")) return 0;
    return ktype(FALCON->kernel());
  }
  //----------------------------------------------------------------------------
  int falcon_current_kernel__()
  {
    if(__falcON_warning("falcon_current_kernel")) return 0;
    return ktype(FALCON->kernel());
  }
  //===========================================================================#
  int falcON_softening()
  {
    if(__falcON_warning("falcON_softening")) return 0;
    return FALCON->use_individual_eps();
  }
  //----------------------------------------------------------------------------
  int falcon_softening_()
  {
    if(__falcON_warning("falcon_softening")) return 0;
    return FALCON->use_individual_eps();
  }
  //----------------------------------------------------------------------------
  int falcon_softening__()
  {
    if(__falcON_warning("falcon_softening")) return 0;
    return FALCON->use_individual_eps();
  }
  //===========================================================================#
  int falcON_No_cells()
  {
    if(__falcON_warning("falcON_No_cells")) return 0;
    return FALCON->No_cells_used();
  }
  //----------------------------------------------------------------------------
  int falcon_No_cells_()
  {
    if(__falcON_warning("falcon_No_cells")) return 0;
    return FALCON->No_cells_used();
  }
  //----------------------------------------------------------------------------
  int falcon_No_cells__()
  {
    if(__falcON_warning("falcon_No_cells")) return 0;
    return FALCON->No_cells_used();
  }
  //===========================================================================#
  int falcON_depth()
  {
    if(__falcON_warning("falcON_depth")) return 0;
    return FALCON->root_depth();
  }
  //----------------------------------------------------------------------------
  int falcon_depth_()
  {
    if(__falcON_warning("falcon_depth")) return 0;
    return FALCON->root_depth();
  }
  //----------------------------------------------------------------------------
  int falcon_depth__()
  {
    if(__falcON_warning("falcon_depth")) return 0;
    return FALCON->root_depth();
  }
  //===========================================================================#
  void falcON_stats()
  {
    if(__falcON_warning("falcON_stats")) return;
    FALCON->stats(std::cout);
  }
  //----------------------------------------------------------------------------
  void falcon_stats_()
  {
    if(__falcON_warning("falcon_stats")) return;
    FALCON->stats(std::cout);
  }
  //----------------------------------------------------------------------------
  void falcon_stats__()
  {
    if(__falcON_warning("falcon_stats")) return;
    FALCON->stats(std::cout);
  }
  //===========================================================================#
  void falcON_set_debug_level(int d)
  {
    falcON::RunInfo::set_debug_level(d);
  }
  //----------------------------------------------------------------------------
  void falcon_set_debug_level_(int*d)
  {
    falcON::RunInfo::set_debug_level(*d);
  }
  //----------------------------------------------------------------------------
  void falcon_set_debug_level__(int*d)
  {
    falcON::RunInfo::set_debug_level(*d);
  }
  //===========================================================================#

} // extern "C"

#include <time.h>
#include <stdlib.h>

extern "C" {

  int    clock_()                { return clock(); }
  int    clock__()               { return clock(); }
  double drand48_()              { return drand48(); }  
  double drand48__()             { return drand48(); }
  double clocks_per_second_()    { return CLOCKS_PER_SEC; }
  double clocks_per_second__()   { return CLOCKS_PER_SEC; }
  void   srand48_(int*S)         { srand48(*S); }
  void   srand48__(int*S)        { srand48(*S); }

} // extern "C"

//=============================================================================#
