//-----------------------------------------------------------------------------+
//                                                                             |
// falcONC.cc                                                                  |
//                                                                             |
// copyright Walter Dehnen, 2000-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#define falcONC_cc
#include <falcON_C.h>
#include <falcON.h>
#include <body.h>
#if falcON_NDIM != 3
#  error falcON_NDIM != 3 in falcONC.cc
#endif
using namespace nbdy;
//=============================================================================#
// auxiliary data                                                              |
//=============================================================================#
namespace {

  static ebodies      ARRAYS;                 // body data arrays               
  static falcON      *FALCON    = 0;          // pointer  to 'our' falcON       
  static bool         BUILT     = false;      // is our tree build or not?      
  //===========================================================================#
  // auxiliary routines                                                        |
  //===========================================================================#
  inline void __falcON_error(const char* name)
  {
    if(FALCON==0) {
      std::cerr<<" falcON ERROR: "<<name
	       <<"() called before falcON_initialize()\n";
      nbdy::exit(1);
    }
  }
  //----------------------------------------------------------------------------
  inline int __falcON_warning(const char* name)
  {
    if(FALCON==0) {
      std::cerr<<" falcON WARNING: "<<name
	       <<"() called before falcON_initialize()\n";
      return 1;
    }
    return 0;
  }
  //----------------------------------------------------------------------------
  inline void __falcON_grown(const char* name)
  {
    if(!BUILT) {
      std::cerr<<" faclON WARNING: "<<name<<"()"
	       <<" called before a tree has been grown\n"
	       <<"        I will grow the tree (via falcON_grow()) first\n";
      FALCON->grow();
      BUILT=true;
    }
  }
  //----------------------------------------------------------------------------
  inline nbdy::kern_type ktype(const int K)
  {
    switch(K%10) {
    case  0: return nbdy::p0;
    case  1: return nbdy::p1;
    case  2: return nbdy::p2;
    case  3: return nbdy::p3;
    case  9: return nbdy::newton;
    default: std::cerr<<" unknown kernel, using Newtonian greens function\n";
      return nbdy::newton;
    }
  }
  //----------------------------------------------------------------------------
  inline int ktype(const nbdy::kern_type K)
  {
    switch(K) {
    case nbdy::p0:     return 0;
    case nbdy::p1:     return 1;
    case nbdy::p2:     return 2;
    case nbdy::p3:     return 3;
    case nbdy::newton: return 9;
    default:
      std::cerr<<" unknown kernel type, defaulting to newton\n";
      return 9;
    }
  }
  //===========================================================================#
  void __falcON_initialize(const int  * F,
			   const areal* M,
			   const areal**X,
#ifdef falcON_INDI
			   const areal* E,
#endif
			         areal**A, 
			         areal* P,
			         areal* R,
			   const int    N,
			   const areal  EPS,
			   const areal  TH,
			   const int    K,
			   const areal  G)
  {
    if(FALCON) delete FALCON;
    ARRAYS.setN    (abs(N));
    ARRAYS.set_mass(const_cast<areal *>(M));
    ARRAYS.set_pos (const_cast<areal**>(X));
#ifdef falcON_INDI
    ARRAYS.set_eps (const_cast<areal *>(E));
#endif
    ARRAYS.set_flg (const_cast<int   *>(F));
    ARRAYS.set_pot (P);
    ARRAYS.set_acc (A);
    ARRAYS.set_rho (R);
    FALCON = new falcON(&ARRAYS,
			abs(real(EPS)),
			abs(real(TH)),
			ktype(K),
#ifdef falcON_INDI
			E != 0,
#endif
			real(G),
			TH<0? const_theta : theta_of_M);
    BUILT  = false;
  }
  //----------------------------------------------------------------------------
  void __falcON_ialist(      int  * I1,
		             int  * I2,
		       const int    Ni,
		             int  * Nt,
		       const areal* S,
		       const bool   Mx,
		       const areal  tau,
		       const             areal**V,
		       const int    shift)
  {
    __falcON_error("falcON_ialist");
    __falcON_grown("falcON_ialist");
    nbdy::uint na;
    falcON::elem_pair *bl = new falcON::elem_pair[Ni];
    ARRAYS.set_vel (const_cast<areal**>(V));
    ARRAYS.set_size(const_cast<areal *>(S));
    FALCON->make_iaction_list(bl,Ni,na,Mx,tau);
    const nbdy::uint iend = min(na,nbdy::uint(Ni));
    for(register nbdy::uint i=0; i<iend; i++) {
      I1[i] = shift+bl[i][0];           // FORTRAN indices = 1 + C indices
      I2[i] = shift+bl[i][1];
    }
    delete[] bl;
    *Nt=na;
  }
}                                                 // END: unnamed namespace    
//=============================================================================#
// routines defined in falcON_C.h and falcon.f                                 |
//=============================================================================#

extern "C" {

//=============================================================================#
void falcON_initialize(const int   *F,
		       const areal *M,
		       const areal**X,
#ifdef falcON_INDI
		       const areal *E,
#endif
		             areal**A,
		             areal *P,
		             areal *R,
		       const int    N,
		       const areal  EPS,
		       const areal  TH,
		       const int    K,
		       const areal  G)
{
  __falcON_initialize(F,M,X,
#ifdef falcON_INDI
		      E,
#endif
		      A,P,R,N,EPS,TH,K,G);
}
//------------------------------------------------------------------------------
void falcon_initialize_(int *F, areal*M, areal**X,
#ifdef falcON_INDI
			areal*E,
#endif
			areal**A,
			areal*P, areal*R, int*N, areal*EPS,
			areal*TH, int*K, areal*G)
{
  __falcON_initialize(F,M,const_cast<const areal**>(X),
#ifdef falcON_INDI
		      E,
#endif
		      A,P,R,*N,*EPS,*TH,*K,*G);
}
//------------------------------------------------------------------------------
void falcon_initialize__(int *F, areal*M, areal**X,
#ifdef falcON_INDI
			 areal*E,
#endif
			 areal**A,
			 areal*P, areal*R, int*N, areal*EPS,
			 areal*TH, int*K, areal*G)
{
  __falcON_initialize(F,M,const_cast<const areal**>(X),
#ifdef falcON_INDI
		      E,
#endif
		      A,P,R,*N,*EPS,*TH,*K,*G);
}
//=============================================================================#
void falcON_resetsoftening(const areal EPS, const int K)
{
  if(__falcON_warning("falcON_resetsoftening")) return;
  FALCON->reset_softening(real(EPS),ktype(K));
}
//------------------------------------------------------------------------------
void falcon_resetsoftening_(areal*EPS, int*K)
{
  if(__falcON_warning("falcon_resetsoftening")) return;
  FALCON->reset_softening(real(*EPS),ktype(*K));
}
//------------------------------------------------------------------------------
void falcon_resetsoftening__(areal*EPS, int*K)
{
  if(__falcON_warning("falcon_resetsoftening")) return;
  FALCON->reset_softening(real(*EPS),ktype(*K));
}
//=============================================================================#
void falcON_resetopening(const areal TH)
{
  if(__falcON_warning("falcON_resetopening")) return;
  FALCON->reset_opening(real(TH),nbdy::theta_of_M);
}  
//------------------------------------------------------------------------------
void falcon_resetopening_(areal *TH)
{
  if(__falcON_warning("falcon_resetopening")) return;
  FALCON->reset_opening(real(*TH),nbdy::theta_of_M);
}  
//------------------------------------------------------------------------------
void falcon_resetopening__(areal *TH)
{
  if(__falcON_warning("falcon_resetopening")) return;
  FALCON->reset_opening(real(*TH),nbdy::theta_of_M);
}  
//=============================================================================#
void falcON_clearup()
{
  if(FALCON) delete FALCON;
  FALCON = 0;
  BUILT  = false;
}
//------------------------------------------------------------------------------
void falcon_clearup_()
{
  if(FALCON) delete FALCON;
  FALCON = 0;
  BUILT  = false;
}
//------------------------------------------------------------------------------
void falcon_clearup__()
{
  if(FALCON) delete FALCON;
  FALCON = 0;
  BUILT  = false;
}
//=============================================================================#
void falcON_grow(const int Nc)
{
  __falcON_error("falcON_grow");
  FALCON->grow(Nc);
  BUILT = true;
}
//------------------------------------------------------------------------------
void falcon_grow_(int *Nc)
{
  __falcON_error("falcon_grow");
  FALCON->grow(*Nc);
  BUILT = true;
}
//------------------------------------------------------------------------------
void falcon_grow__(int *Nc)
{
  __falcON_error("falcon_grow");
  FALCON->grow(*Nc);
  BUILT = true;
}
//=============================================================================#
void falcON_reuse()
{
  __falcON_error("falcON_reuse");
  if(!BUILT) {
    std::cerr<<" faclON WARNING: falcON_reuse()"
	     <<" called before a tree has been grown\n"
	     <<"        I will grow the tree (via falcON_grow()) instead\n";
    FALCON->grow();
    BUILT = true;
  } else
    FALCON->reuse();
}
//------------------------------------------------------------------------------
void falcon_reuse_()
{
  __falcON_error("falcon_reuse");
  if(!BUILT) {
    std::cerr<<" faclON WARNING: falcon_reuse()"
	     <<" called before a tree has been grown\n"
	     <<"        I will grow the tree (via falcon_grow()) instead\n";
    FALCON->grow();
    BUILT = true;
  } else
    FALCON->reuse();
}
//------------------------------------------------------------------------------
void falcon_reuse__()
{
  __falcON_error("falcon_reuse");
  if(!BUILT) {
    std::cerr<<" faclON WARNING: falcon_reuse()"
	     <<" called before a tree has been grown\n"
	     <<"        I will grow the tree (via falcon_grow()) instead\n";
    FALCON->grow();
    BUILT = true;
  } else
    FALCON->reuse();
}
//=============================================================================#
void falcON_approx_grav()
{
  __falcON_error("falcON_approx_gravity");
  __falcON_grown("falcON_approx_gravity");
  FALCON->approximate_gravity();
}  
//------------------------------------------------------------------------------
void falcon_approx_grav_()
{
  __falcON_error("falcon_approx_gravity");
  __falcON_grown("falcon_approx_gravity");
  FALCON->approximate_gravity();
}  
//------------------------------------------------------------------------------
void falcon_approx_grav__()
{
  __falcON_error("falcon_approx_gravity");
  __falcON_grown("falcon_approx_gravity");
  FALCON->approximate_gravity();
}  
#ifdef falcON_ADAP
//=============================================================================#
void falcON_adjust_epsi_and_approx_grav(areal Nsoft, int Nref, areal fac)
{
  __falcON_error("falcON_adjust_epsi_and_approx_gravity");
  __falcON_grown("falcON_adjust_epsi_and_approx_gravity");
  FALCON->approximate_gravity(true,false,real(Nsoft),Nref,real(fac));
}
//------------------------------------------------------------------------------
void falcon_adjust_epsi_and_approx_grav_(areal*Nsoft, int*Nref, areal*fac)
{
  __falcON_error("falcon_adjust_epsi_and_approx_gravity");
  __falcON_grown("falcon_adjust_epsi_and_approx_gravity");
  FALCON->approximate_gravity(true,false,real(*Nsoft),*Nref,real(*fac));
}
//------------------------------------------------------------------------------
void falcon_adjust_epsi_and_approx_grav__(areal*Nsoft, int*Nref, areal*fac)
{
  __falcON_error("falcon_adjust_epsi_and_approx_gravity");
  __falcON_grown("falcon_adjust_epsi_and_approx_gravity");
  FALCON->approximate_gravity(true,false,real(*Nsoft),*Nref,real(*fac));
}
//=============================================================================#
#endif
void falcON_estimate_rho(const int nx)
{
  __falcON_error("falcON_estimate_rho");
  __falcON_grown("falcON_estimate_rho");
  FALCON->estimate_rho(nx);
}  
//------------------------------------------------------------------------------
void falcon_estimate_rho_(int*nx)
{
  __falcON_error("falcon_estimate_rho");
  __falcON_grown("falcon_estimate_rho");
  FALCON->estimate_rho(*nx);
}  
//------------------------------------------------------------------------------
void falcon_estimate_rho__(int*nx)
{
  __falcON_error("falcon_estimate_rho");
  __falcON_grown("falcon_estimate_rho");
  FALCON->estimate_rho(*nx);
}  
//=============================================================================#
void falcON_estimate_n(const int nx)
{
  __falcON_error("falcON_estimate_n");
  __falcON_grown("falcON_estimate_n");
  FALCON->estimate_n(nx);
}  
//------------------------------------------------------------------------------
void falcon_estimate_n_(int*nx)
{
  __falcON_error("falcon_estimate_n");
  __falcON_grown("falcon_estimate_n");
  FALCON->estimate_n(*nx);
}  
//------------------------------------------------------------------------------
void falcon_estimate_n__(int*nx)
{
  __falcON_error("falcon_estimate_n");
  __falcON_grown("falcon_estimate_n");
  FALCON->estimate_n(*nx);
}  
//=============================================================================#
void falcON_iactionlist(int*I1, int*I2, const int Ni, int* Nt, const areal* S,
			const bool Mx, const areal tau, const areal**V)
{
  __falcON_ialist(I1,I2,Ni,Nt,S,Mx,tau,V,0);
}
//=============================================================================#
void falcon_sticky_ (areal*t,int*i1,int*i2,int*ni,int*na,
		     areal**V,areal*s) {
  __falcON_ialist(i1,i2,*ni,na,s,1,*t,const_cast<const areal**>(V),1); }
void falcon_sticky__(areal*t,int*i1,int*i2,int*ni,int*na,
		     areal**V,areal*s) {
  __falcON_ialist(i1,i2,*ni,na,s,1,*t,const_cast<const areal**>(V),1); }
//=============================================================================#
void falcon_sph_ (int*i1,int*i2,int*ni,int*na,int*mx,areal*s) {
  __falcON_ialist(i1,i2,*ni,na,s,*mx,-1,0,1); }
void falcon_sph__(int*i1,int*i2,int*ni,int*na,int*mx,areal*s) {
  __falcON_ialist(i1,i2,*ni,na,s,*mx,-1,0,1); }

//=============================================================================#
areal falcON_root_center(const int i)
{
  if(__falcON_warning("falcON_root_center")) return 0.;
  return areal((FALCON->root_center())[i]);
}
//------------------------------------------------------------------------------
float falcon_root_center_(int*i)
{
  if(__falcON_warning("falcon_root_center")) return 0.;
  return float((FALCON->root_center())[*i]);
}
//------------------------------------------------------------------------------
float falcon_root_center__(int*i)
{
  if(__falcON_warning("falcon_root_center")) return 0.;
  return float((FALCON->root_center())[*i]);
}
//=============================================================================#
areal falcON_root_radius()
{
  if(__falcON_warning("falcON_root_radius")) return 0.;
  return areal((FALCON->root_radius()));
}
//------------------------------------------------------------------------------
float falcon_root_radius_()
{
  if(__falcON_warning("falcon_root_radius")) return 0.;
  return float((FALCON->root_radius()));
}
//------------------------------------------------------------------------------
float falcon_root_radius__()
{
  if(__falcON_warning("falcon_root_radius")) return 0.;
  return float((FALCON->root_radius()));
}
//=============================================================================#
areal falcON_current_eps()
{
  if(__falcON_warning("falcON_current_eps")) return 0.;
  return areal(FALCON->softening_length());
}
//------------------------------------------------------------------------------
float falcon_current_eps_()
{
  if(__falcON_warning("falcon_current_eps")) return 0.;
  return float(FALCON->softening_length());
}
//------------------------------------------------------------------------------
float falcon_current_eps__()
{
  if(__falcON_warning("falcon_current_eps")) return 0.;
  return float(FALCON->softening_length());
}
//=============================================================================#
int falcON_current_kernel()
{
  if(__falcON_warning("falcON_current_kernel")) return 0;
  return ktype(FALCON->kernel());
}
//------------------------------------------------------------------------------
int falcon_current_kernel_()
{
  if(__falcON_warning("falcon_current_kernel")) return 0;
  return ktype(FALCON->kernel());
}
//------------------------------------------------------------------------------
int falcon_current_kernel__()
{
  if(__falcON_warning("falcon_current_kernel")) return 0;
  return ktype(FALCON->kernel());
}
//=============================================================================#
int falcON_softening()
{
  if(__falcON_warning("falcON_softening")) return 0;
#ifdef falcON_INDI
  return FALCON->use_individual_eps();
#else
  return 0;
#endif
}
//------------------------------------------------------------------------------
int falcon_softening_()
{
  if(__falcON_warning("falcon_softening")) return 0;
#ifdef falcON_INDI
  return FALCON->use_individual_eps();
#else
  return 0;
#endif
}
//------------------------------------------------------------------------------
int falcon_softening__()
{
  if(__falcON_warning("falcon_softening")) return 0;
#ifdef falcON_INDI
  return FALCON->use_individual_eps();
#else
  return 0;
#endif
}
//=============================================================================#
int falcON_No_cells()
{
  if(__falcON_warning("falcON_No_cells")) return 0;
  return FALCON->No_cells_used();
}
//------------------------------------------------------------------------------
int falcon_No_cells_()
{
  if(__falcON_warning("falcon_No_cells")) return 0;
  return FALCON->No_cells_used();
}
//------------------------------------------------------------------------------
int falcon_No_cells__()
{
  if(__falcON_warning("falcon_No_cells")) return 0;
  return FALCON->No_cells_used();
}
//=============================================================================#
int falcON_depth()
{
  if(__falcON_warning("falcON_depth")) return 0;
  return FALCON->root_depth();
}
//------------------------------------------------------------------------------
int falcon_depth_()
{
  if(__falcON_warning("falcon_depth")) return 0;
  return FALCON->root_depth();
}
//------------------------------------------------------------------------------
int falcon_depth__()
{
  if(__falcON_warning("falcon_depth")) return 0;
  return FALCON->root_depth();
}
//=============================================================================#
void falcON_stats()
{
  if(__falcON_warning("falcON_stats")) return;
  FALCON->stats(std::cout);
}
//------------------------------------------------------------------------------
void falcon_stats_()
{
  if(__falcON_warning("falcon_stats")) return;
  FALCON->stats(std::cout);
}
//------------------------------------------------------------------------------
void falcon_stats__()
{
  if(__falcON_warning("falcon_stats")) return;
  FALCON->stats(std::cout);
}
//=============================================================================#

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
