//-----------------------------------------------------------------------------+
//                                                                             |
// falcONC.cc                                                                  |
//                                                                             |
// copyright Walter Dehnen, 2000-2002                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
#define falcONC_cc
#include <falcON_C.h>
#include <falcON.h>
using namespace nbdy;
//==============================================================================
#if NDIM==3
#  define POSIS   const areal*X,  const areal*Y,  const areal*Z
#  define POSS    X,Y,Z
#  define VELIS   const areal*VX, const areal*VY, const areal*VZ
#  define VELS    VX,VY,VZ
#  define ACCIS         areal*AX,       areal*AY,       areal*AZ
#  define ACCS    AX,AY,AZ
#  define ZEROS   0,0,0
#else
#  define POSIS   const areal*X,  const areal*Y
#  define POSS    X,Y
#  define VELIS   const areal*VX, const areal*VY
#  define VELS    VX,VY
#  define ACCIS         areal*AX,       areal*AY
#  define ACCS    AX,AY
#  define ZEROS   0,0
#endif
//=============================================================================#
// auxiliary data                                                              |
//=============================================================================#
static falcON      *FALCON    = 0;            // pointer  to 'our' falcON       
static const areal *POS[NDIM] = {0};          // pointers to position arrays    
static const areal *VEL[NDIM] = {0};          // pointers to velocity arrays    
static areal       *ACC[NDIM] = {0};          // pointers to acceleration arrays
static bool         BUILT     = false;        // is our tree build or not?      
//=============================================================================#
// auxiliary routines                                                          |
//=============================================================================#
inline void falcON_error(const char* name)
{
  if(FALCON==0) {
    std::cerr<<" falcON ERROR: "<<name
	     <<"() called before falcON_initialize()\n";
    nbdy::exit(1);
  }
}
//------------------------------------------------------------------------------
inline int falcON_warning(const char* name)
{
  if(FALCON==0) {
    std::cerr<<" falcON WARNING: "<<name
	     <<"() called before falcON_initialize()\n";
    return 1;
  }
  return 0;
}
//------------------------------------------------------------------------------
inline void falcON_grown(const char* name)
{
  if(!BUILT) {
    std::cerr<<" faclON WARNING: "<<name<<"()"
	     <<" called before a tree has been grown\n"
	     <<"        I will grow the tree (via falcON_grow()) first\n";
    FALCON->grow();
    BUILT=true;
  }
}
//------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------
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
//=============================================================================#
inline
void __falcON_initialize(const int *F, const areal*M, POSIS, areal*E, ACCIS,
			 areal*PH, areal*RH, const int N, const areal EPS,
			 const areal TH, const int K)
{
  if(FALCON) delete FALCON;
  POS[0] = X; ACC[0] = AX;
  POS[1] = Y; ACC[1] = AY;
#if NDIM==3
  POS[2] = Z; ACC[2] = AZ;
#endif
#ifdef ALLOW_INDI
  FALCON = new falcON(F,M,POS,E,ACC,PH,RH,
		      abs(N),
		      abs(real(EPS)),
		      abs(real(TH)),
		      ktype(K),
		      E==0? global      : individual,
		      TH<0? const_theta : theta_of_M);
#else
  if(E) warning("[falcON_initialize()]: individual eps_i given, but not enabled"
		"in this version (see make.defs).\n"
		"                We will use a global eps = %g",abs(real(EPS)));
  FALCON = new falcON(F,M,POS,ACC,PH,RH,
		      abs(N),
		      abs(real(EPS)),
		      abs(real(TH)),
		      ktype(K),
		      TH<0? const_theta : theta_of_M);
#endif
  BUILT  = false;
}
//------------------------------------------------------------------------------
inline 
void __falcON_ialist(int*I1, int*I2, const int Ni, int* Nt, const areal* S,
		     const areal tau, VELIS)
{
  falcON_error("falcON_ialist");
  falcON_grown("falcON_ialist");
  nbdy::uint na;
  falcON::elem_pair *bl = new falcON::elem_pair[Ni];
  VEL[0] = VX;
  VEL[1] = VY;
#if NDIM==3
  VEL[2] = VZ;
#endif
  FALCON->make_iaction_list(bl,Ni,na,S,VEL,tau);
  const nbdy::uint iend = min(na,nbdy::uint(Ni));
  for(register nbdy::uint i=0; i<iend; i++) {
    I1[i] = 1+bl[i][0];           // FORTRAN indices = 1 + C indices
    I2[i] = 1+bl[i][1];
  }
  delete[] bl;
  *Nt=na;
}
//=============================================================================#
// routines defined in falcON_C.h and falcon.f                                 |
//=============================================================================#

extern "C" {

//=============================================================================#
void falcON_initialize(const int *F, const areal*M, POSIS, areal*E, ACCIS,
		       areal*PH, areal*RH, const int N, const areal EPS,
		       const areal TH, const int K)
{
  __falcON_initialize(F,M,POSS,E,ACCS,PH,RH,N,EPS,TH,K);
}
//------------------------------------------------------------------------------
void falcon_init_     (int *F, areal*M, POSIS, areal*E, ACCIS,
		       areal*PH, areal*RH, int*N, areal*EPS,
		       areal*TH, int*K)
{
  __falcON_initialize(F,M,POSS,E,ACCS,PH,RH,*N,*EPS,*TH,*K);
}
//------------------------------------------------------------------------------
void falcon_init__    (int *F, areal*M, POSIS, areal*E, ACCIS,
		       areal*PH, areal*RH, int*N, areal*EPS,
		       areal*TH, int*K)
{
  __falcON_initialize(F,M,POSS,E,ACCS,PH,RH,*N,*EPS,*TH,*K);
}
//=============================================================================#
void falcON_resetsoftening(const areal EPS, const int K)
{
  if(falcON_warning("falcON_resetsoftening")) return;
  FALCON->reset_softening(real(EPS),ktype(K));
}
//------------------------------------------------------------------------------
void falcon_resetsoftening_(areal*EPS, int*K)
{
  if(falcON_warning("falcon_resetsoftening")) return;
  FALCON->reset_softening(real(*EPS),ktype(*K));
}
//------------------------------------------------------------------------------
void falcon_resetsoftening__(areal*EPS, int*K)
{
  if(falcON_warning("falcon_resetsoftening")) return;
  FALCON->reset_softening(real(*EPS),ktype(*K));
}
//=============================================================================#
void falcON_resetopening(const areal TH)
{
  if(falcON_warning("falcON_resetopening")) return;
  FALCON->reset_opening(real(TH),nbdy::theta_of_M);
}  
//------------------------------------------------------------------------------
void falcon_resetopening_(areal *TH)
{
  if(falcON_warning("falcon_resetopening")) return;
  FALCON->reset_opening(real(*TH),nbdy::theta_of_M);
}  
//------------------------------------------------------------------------------
void falcon_resetopening__(areal *TH)
{
  if(falcON_warning("falcon_resetopening")) return;
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
  falcON_error("falcON_grow");
  FALCON->grow(Nc);
  BUILT = true;
}
//------------------------------------------------------------------------------
void falcon_grow_(int *Nc)
{
  falcON_error("falcon_grow");
  FALCON->grow(*Nc);
  BUILT = true;
}
//------------------------------------------------------------------------------
void falcon_grow__(int *Nc)
{
  falcON_error("falcon_grow");
  FALCON->grow(*Nc);
  BUILT = true;
}
//=============================================================================#
void falcON_reuse()
{
  falcON_error("falcON_reuse");
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
  falcON_error("falcon_reuse");
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
  falcON_error("falcon_reuse");
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
  falcON_error("falcON_approx_gravity");
  falcON_grown("falcON_approx_gravity");
  FALCON->approximate_gravity();
}  
//------------------------------------------------------------------------------
void falcon_approx_grav_()
{
  falcON_error("falcon_approx_gravity");
  falcON_grown("falcon_approx_gravity");
  FALCON->approximate_gravity();
}  
//------------------------------------------------------------------------------
void falcon_approx_grav__()
{
  falcON_error("falcon_approx_gravity");
  falcON_grown("falcon_approx_gravity");
  FALCON->approximate_gravity();
}  
#ifdef ALLOW_INDI
//=============================================================================#
void falcON_adjust_epsi_and_approx_grav(areal Nsoft, int Nref, areal fac)
{
  falcON_error("falcON_adjust_epsi_and_approx_gravity");
  falcON_grown("falcON_adjust_epsi_and_approx_gravity");
  FALCON->approximate_gravity(true,real(Nsoft),Nref,real(fac));
}
//------------------------------------------------------------------------------
void falcon_adjust_epsi_and_approx_grav_(areal*Nsoft, int*Nref, areal*fac)
{
  falcON_error("falcon_adjust_epsi_and_approx_gravity");
  falcON_grown("falcon_adjust_epsi_and_approx_gravity");
  FALCON->approximate_gravity(true,real(*Nsoft),*Nref,real(*fac));
}
//------------------------------------------------------------------------------
void falcon_adjust_epsi_and_approx_grav__(areal*Nsoft, int*Nref, areal*fac)
{
  falcON_error("falcon_adjust_epsi_and_approx_gravity");
  falcON_grown("falcon_adjust_epsi_and_approx_gravity");
  FALCON->approximate_gravity(true,real(*Nsoft),*Nref,real(*fac));
}
//=============================================================================#
#endif
void falcON_estimate_rho(const int nx)
{
  falcON_error("falcON_estimate_rho");
  falcON_grown("falcON_estimate_rho");
  FALCON->estimate_rho(nx);
}  
//------------------------------------------------------------------------------
void falcon_estimate_rho_(int*nx)
{
  falcON_error("falcon_estimate_rho");
  falcON_grown("falcon_estimate_rho");
  FALCON->estimate_rho(*nx);
}  
//------------------------------------------------------------------------------
void falcon_estimate_rho__(int*nx)
{
  falcON_error("falcon_estimate_rho");
  falcON_grown("falcon_estimate_rho");
  FALCON->estimate_rho(*nx);
}  
//=============================================================================#
void falcON_estimate_n(const int nx)
{
  falcON_error("falcON_estimate_n");
  falcON_grown("falcON_estimate_n");
  FALCON->estimate_n(nx);
}  
//------------------------------------------------------------------------------
void falcon_estimate_n_(int*nx)
{
  falcON_error("falcon_estimate_n");
  falcON_grown("falcon_estimate_n");
  FALCON->estimate_n(*nx);
}  
//------------------------------------------------------------------------------
void falcon_estimate_n__(int*nx)
{
  falcON_error("falcon_estimate_n");
  falcON_grown("falcon_estimate_n");
  FALCON->estimate_n(*nx);
}  
//=============================================================================#
void falcON_iactionlist(int*I1, int*I2, const int Ni, int* Nt, const areal* S,
			const areal tau, VELIS)
{
  falcON_error("falcON_iactionlist");
  falcON_grown("falcON_iactionlist");
  nbdy::uint na;
  falcON::elem_pair *bl = new falcON::elem_pair[Ni];
  VEL[0] = VX;
  VEL[1] = VY;
#if NDIM==3
  VEL[2] = VZ;
#endif
  FALCON->make_iaction_list(bl,Ni,na,S,VEL,tau);
  const nbdy::uint iend = min(na,nbdy::uint(Ni));
  for(register nbdy::uint i=0; i<iend; i++) {
    I1[i] = bl[i][0];
    I2[i] = bl[i][1];
  }
  delete[] bl;
  *Nt=na;
}
//=============================================================================#
void falcon_sticky_ (areal*t,int*i1,int*i2,int*ni,int*na,
		     VELIS,areal*s) {
  __falcON_ialist(i1,i2,*ni,na,s,*t,VELS); }
void falcon_sticky__(areal*t,int*i1,int*i2,int*ni,int*na,
		     VELIS,areal*s) {
  __falcON_ialist(i1,i2,*ni,na,s,*t,VELS); }
//=============================================================================#
void falcon_sph_ (int*i1,int*i2,int*ni,int*na,areal*s) {
  __falcON_ialist(i1,i2,*ni,na,s,-1,ZEROS); }
void falcon_sph__(int*i1,int*i2,int*ni,int*na,areal*s) {
  __falcON_ialist(i1,i2,*ni,na,s,-1,ZEROS); }

//=============================================================================#
areal falcON_root_center(const int i)
{
  if(falcON_warning("falcON_root_center")) return 0.;
  return areal((FALCON->root_center())[i]);
}
//------------------------------------------------------------------------------
float falcon_root_center_(int*i)
{
  if(falcON_warning("falcon_root_center")) return 0.;
  return float((FALCON->root_center())[*i]);
}
//------------------------------------------------------------------------------
float falcon_root_center__(int*i)
{
  if(falcON_warning("falcon_root_center")) return 0.;
  return float((FALCON->root_center())[*i]);
}
//=============================================================================#
areal falcON_root_radius()
{
  if(falcON_warning("falcON_root_radius")) return 0.;
  return areal((FALCON->root_radius()));
}
//------------------------------------------------------------------------------
float falcon_root_radius_()
{
  if(falcON_warning("falcon_root_radius")) return 0.;
  return float((FALCON->root_radius()));
}
//------------------------------------------------------------------------------
float falcon_root_radius__()
{
  if(falcON_warning("falcon_root_radius")) return 0.;
  return float((FALCON->root_radius()));
}
//=============================================================================#
areal falcON_current_eps()
{
  if(falcON_warning("falcON_current_eps")) return 0.;
  return areal(FALCON->softening_length());
}
//------------------------------------------------------------------------------
float falcon_current_eps_()
{
  if(falcON_warning("falcon_current_eps")) return 0.;
  return float(FALCON->softening_length());
}
//------------------------------------------------------------------------------
float falcon_current_eps__()
{
  if(falcON_warning("falcon_current_eps")) return 0.;
  return float(FALCON->softening_length());
}
//=============================================================================#
int falcON_current_kernel()
{
  if(falcON_warning("falcON_current_kernel")) return 0;
  return ktype(FALCON->kernel());
}
//------------------------------------------------------------------------------
int falcon_current_kernel_()
{
  if(falcON_warning("falcon_current_kernel")) return 0;
  return ktype(FALCON->kernel());
}
//------------------------------------------------------------------------------
int falcon_current_kernel__()
{
  if(falcON_warning("falcon_current_kernel")) return 0;
  return ktype(FALCON->kernel());
}
//=============================================================================#
int falcON_softening()
{
  if(falcON_warning("falcON_softening")) return 0;
#ifdef ALLOW_INDI
  switch(FALCON->softening()) {
  case nbdy::global:              return 0;
  case nbdy::individual: default: return 1;
  }
#else
  return 0;
#endif
}
//------------------------------------------------------------------------------
int falcon_softening_()
{
  if(falcON_warning("falcon_softening")) return 0;
#ifdef ALLOW_INDI
  switch(FALCON->softening()) {
  case nbdy::global:              return 0;
  case nbdy::individual: default: return 1;
  }
#else
  return 0;
#endif
}
//------------------------------------------------------------------------------
int falcon_softening__()
{
  if(falcON_warning("falcon_softening")) return 0;
#ifdef ALLOW_INDI
  switch(FALCON->softening()) {
  case nbdy::global:              return 0;
  case nbdy::individual: default: return 1;
  }
#else
  return 0;
#endif
}
//=============================================================================#
int falcON_No_cells()
{
  if(falcON_warning("falcON_No_cells")) return 0;
  return FALCON->No_cells_used();
}
//------------------------------------------------------------------------------
int falcon_No_cells_()
{
  if(falcON_warning("falcon_No_cells")) return 0;
  return FALCON->No_cells_used();
}
//------------------------------------------------------------------------------
int falcon_No_cells__()
{
  if(falcON_warning("falcon_No_cells")) return 0;
  return FALCON->No_cells_used();
}
//=============================================================================#
int falcON_depth()
{
  if(falcON_warning("falcON_depth")) return 0;
  return FALCON->root_depth();
}
//------------------------------------------------------------------------------
int falcon_depth_()
{
  if(falcON_warning("falcon_depth")) return 0;
  return FALCON->root_depth();
}
//------------------------------------------------------------------------------
int falcon_depth__()
{
  if(falcON_warning("falcon_depth")) return 0;
  return FALCON->root_depth();
}
//=============================================================================#
void falcON_stats()
{
  if(falcON_warning("falcON_stats")) return;
  FALCON->stats(std::cout);
}
//------------------------------------------------------------------------------
void falcon_stats_()
{
  if(falcON_warning("falcon_stats")) return;
  FALCON->stats(std::cout);
}
//------------------------------------------------------------------------------
void falcon_stats__()
{
  if(falcON_warning("falcon_stats")) return;
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
