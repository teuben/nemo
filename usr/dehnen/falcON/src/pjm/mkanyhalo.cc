// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// mkanyhalo.cc                                                                |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2005                                          |
//           Paul McMillan, 2004-2005                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
//           paul.mcmillan@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// creates N-body initial conditions for any halo with an external potential   |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 1.0    21/03/2005  PJM bastardized mkedcm to create this                  |
// v 1.1    29/09/2005  PJM added Ossipkov-Merritt type anisotropy radius      |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "2.3"
#define falcON_VERSION_D "17-may-2004 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile mkanyhalo
#endif
#define falcON_RepAction 0                         // no action reporting       
//-----------------------------------------------------------------------------+
#include <public/sample.h>                         // my N-body sampler         
#include <body.h>                                  // N-body bodies             
#include <public/io.h>                             // my I/O utilities          
#include <pjm/anyhalo.h>                           // my ya know, thing   
#include <fstream>                                 // C++ file I/O              
#include <iomanip>                                 // C++ I/O formatting        
#include <main.h>                                  // main & NEMO stuff         
////////////////////////////////////////////////////////////////////////////////
namespace {
  using namespace falcON;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class AnyHaloSampler                                                    //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class AnyHaloSampler :
    private AnyHalo,
    public  SphericalSampler
  {
  public:
    AnyHaloSampler( double const&__MM,             // I: Mdisc/Mhalo         
		    double const&__inn,            // I: inner density exponent         
		    double const&__out,            // I: outer density exponent        
		    double const&__r_s,            // I: scale radius        
		    double const&__r_t,            // I: trunc. radius       
		    double const&__beta,           // I: anisotropy             
		    double const&__r_a,            // I: Ossipkov-Merritt radius
		    const double*__r,              // I: mass adaption: radii   
		    int    const&__n,              // I: mass adaption: # --    
		    double const&__f,              // I: mass adaption: factor  
		    bool   const&__p,              // I: mass adaption: R_-/Re
		    bool   const&__care,            // I: non-monotonic DF?
		    const acceleration *aex) :      // I: external potential  
      AnyHalo         ( __MM,__inn,__out,__r_s,__r_t, __beta,__r_a,aex ),
      SphericalSampler( total_mass(),__r_a, __beta,__r,__n,__f,__p,__care ) {}
    //--------------------------------------------------------------------------
    AnyHalo const&TCDM() const { return *this; }
    //--------------------------------------------------------------------------
    double DF(double q)           const { return AnyHalo::g_E(q); }
    double Ps(double r)           const { return AnyHalo::psi(r); }
    double rM(double m)           const { return AnyHalo::rM_h(m); }
    double Re(double e)           const { return AnyHalo::RcE(e); } 
    double Rp(double e, double l) const { return AnyHalo::Rp(e,l*l); }
  };
}                                                  // END: unnamed namespace    
////////////////////////////////////////////////////////////////////////////////
string defv[] = {
  "out=???\n          output file                                        ",
  "nbody=???\n        number of bodies                                   ",
  "MonM=1\n           Mexternal/Mhalo                                    ",
  "inner=1\n          inner density exponent                             ",
  "outer=4\n          outer density exponent                             ",
  "r_s=1\n            Scale radius                                       ",
  "r_t=0\n            truncation radius (if 0, no truncation)            ",
  "beta=0\n           anisotropy (if r_a =0); -1.5 <= beta <= g/2        ",
  "r_a=\n             Ossipkov-Merritt anisotropy radius                 ",
  "seed=0\n           seed for the randum number generator               ",
  "q-ran=f\n          use quasi- instead of pseudo-random numbers        ",
  "time=0\n           simulation time of snapshot                        ",
  "f_pos=0.5\n        fraction of bodies with positive sense of rotation ",
  "Rp=\n              for mass adaption: list of R in increasing order   ",
  "fac=1.2\n          for mass adaption: factor between mass bins        ",
  "peri=f\n           for mass adaption: R_peri(E,L) rather than R_c(E)  ",
  "epar=\n            if given, set eps_i = epar * sqrt(m_i/M_tot)       ",
  "WD_units=f\n       input:  kpc, km/s\n"
  "                   output: kpc, kpc/Gyr, G=1 (-> mass unit)           ",
  "outputs=f\n        give some global quantities to stderr              ",
  "tabfile=\n         write table with r, pot, vc, sigma to file         ",
  "potname=\n         name of external acceleration field                ",
  "potpars=\n         parameters of external acceleration field          ",
  "potfile=\n         file required by external acceleration field       ",
  "careful=f\n        possibly non-monotonic DF?                         ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string usage =
"mkanyhalo -- construct a truncated CDM model with density profile\n"
"\n"
"                               sech(r/r_t)\n"
"           rho(r) = C ------------------------------\n"
"                       r^inner (r+r_s)^(outer-inner)\n" 
"\n"
"\n"
"         Model has constant anisotropy beta = 1-(sigma_theta/sigma_r)^2\n"
"    if r_a=0, else it follows Cuddeford's (1991) generalised models\n"
"\n"
"         N.B. requires external monopole field\n"
"\n";
//------------------------------------------------------------------------------
void falcON::main()
{
  const double vf = 0.977775320024919;             // km/s  in WD_units         
  const double mf = 2.2228847e5;                   // M_sun in WD_units         
  //----------------------------------------------------------------------------
  // 1. set some parameters                                                     
  //----------------------------------------------------------------------------
  const bool   WD (getbparam("WD_units"));         // using WD_units?           
  const bool care (getbparam("careful"));          // non-monotonic DF?        
  const Random Ran(getparam("seed"),6);
  const fieldset data( 
    (hasvalue("epar")  ? fieldset::e : fieldset::o) |
    fieldset::basic);

  const nemo_acc*aex = hasvalue("potname")?           // IF(potname given) THEN 
    new nemo_acc(getparam  ("potname"),               //   initialize external  
		 getparam_z("potpars"),               //   accelerations        
		 getparam_z("potfile")) : 0;          // ELSE: no potential     
  const int    nbmax(100);
  double       Rad[nbmax];
  int          nb=0;
  if(hasvalue("Rp")) {
    nb=nemoinpd(getparam("Rp"),Rad,nbmax);
    if(nb+1>nbmax)
      ::error("exceeding expected number of radii in mass adaption");
    Rad[nb] = Rad[nb-1] * 1.e20;
    nb++;
  }
  //----------------------------------------------------------------------------
  // 2. create initial conditions from a ecdm model using mass adaption       
  //----------------------------------------------------------------------------
  AnyHaloSampler TS( getdparam("MonM"),
		     getdparam("inner"),
		     getdparam("outer"),
		     getdparam("r_s"),
		     getdparam("r_t"),
		     getdparam("beta"),
		     getdparam("r_a"),
		     Rad,nb,
		     getdparam("fac"),
		     getbparam("peri"),
		     care, aex);
  snapshot shot(getdparam("time"), getiparam("nbody"), data);
  TS.sample(shot,getbparam("q-ran"),Ran,getdparam("f_pos"),getdparam("epar"),false);
  //----------------------------------------------------------------------------
  // 3. output of snapshot                                                      
  //----------------------------------------------------------------------------
  nemo_out out(getparam("out"));
  shot.write_nemo(out,data);

};
////////////////////////////////////////////////////////////////////////////////
