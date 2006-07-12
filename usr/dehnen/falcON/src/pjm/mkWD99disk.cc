// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// mkWD99disk.cc                                                                |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2004                                          |
//           Paul McMillan, 2005                                               |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
//     creates N-body initial conditions for a galaxy disk based on the        |
//       distribution function f_new from Dehnen 1999 (ApJ 118, 1201)          |
//                                                                             |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v1.0    24/06/2005   PJM started the process of making this                 |
// v1.1       09/2005   PJM added iteration of the density profile             |
// v1.2       01/2006   PJM added iteration of the velocity distribution       |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "2.3"
#define falcON_VERSION_D "24-jun-2005 Paul McMillan                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile mkWD99disk
#endif
//-----------------------------------------------------------------------------+
#include <body.h>                                  // N-body bodies             
#include <public/io.h>                             // my I/O utilities          
#include <pjm/WD99disk.h>                          // my disk models   
#include <fstream>                                 // C++ file I/O              
#include <iomanip>                                 // C++ I/O formatting        
#include <main.h>                                  // main & NEMO stuff         
////////////////////////////////////////////////////////////////////////////////
string defv[] = {
  "out=???\n          output file                                        ",
  "nbody=???\n        number of bodies                                   ",
  "nbpero=100\n       number of bodies per orbit                         ",
  "r_d=1\n            disk scale radius: Surf dens=Sig_0*(e^(-R/r_d))    ",
  "Sig_0=???\n        central surface density, see above                 ",
  "r_sig=0\n          vel. disp. scale radius, sigr propto e^(-R/r_sig)  ",
  "Qmin=???\n         Toomre's Q, const if r_sig=0, else Q at rsig       ",
  "z_0=0.1\n          vertical scale height                              ",
  "eps=0\n            particle smoothing length (undefined if 0)         ",
  "seed=0\n           seed for the randum number generator               ",
  "q-ran=f\n          use quasi- instead of pseudo-random numbers        ",
  "time=0\n           simulation time of snapshot                        ",
  "ni=3\n             no. iterations of disk surf dens & vel disp (min 1)",
  "WD_units=f\n       input:  kpc, km/s\n"
  "                   output: kpc, kpc/Gyr, G=1 (-> mass unit)           ",
  "outputs=f\n        give some global quantities to stderr              ",
  "tabfile=\n         write table with r, pot, vc, sigma to file         ",
  "potname=\n         name of external acceleration field (required)     ",
  "potpars=\n         parameters of external acceleration field          ",
  "potfile=\n         file required by external acceleration field       ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string usage =
"mkWD99disk: construct a disk according to Dehnen (1999)\n"
"            given external NEMO potential field corresponding to the\n"
"            potential of the entire system (including desired disk)\n"
"\n"
"    Surface density:             Sigma(R) = Sigma_0 * exp(R/r_d)\n"
"\n"
"    Radial velocity dispersion: defined by constant Q (if r_sig=0)\n"
"\n"
"                                else sigma(R) = sigma_0 * exp(R/r_sig)\n"
"                                defined such that Q is Qmin at r_sig\n"
"\n"
"              N.B. possible confusion between sigma/Sigma;\n"
"          and that Qmin isn't always a minimum since I changed the code\n"
"\n"
"                           sigma(R) * kappa(R)\n"
"                     Q =   -------------------\n"
"                           3.36 * G * Sigma(R)\n"
"\n"
"                with kappa(R) the epicycle frequency\n"
"\n"
"  z-component is an isothemal sheet with rhoz propto sech^2(z/z_0)\n"
"\n"
" BEWARE: known \"feature\". On occasion there is an error report from\n"
" \"find\". Can be avoided by reducing ni, or increasing nbpero \n";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  const double vf = 0.977775320024919;             // km/s  in WD_units         
  const double mf = 2.2228847e5;                   // M_sun in WD_units    
  //----------------------------------------------------------------------------
  // 1. set some parameters                                                     
  //----------------------------------------------------------------------------
  if (getiparam("ni") < 1) error("Code requires at least one iteration");

  const bool   WD (getbparam("WD_units"));         // using WD_units?           
  const Random Ran(getparam("seed"),8);

  const fieldset data((getdparam("eps")? fieldset::e : fieldset::o) |
		      fieldset::basic);
  const nemo_acc*aex = hasvalue("potname")?           // IF(potname given) THEN 
    new nemo_acc(getparam  ("potname"),               //   initialize external  
		 getparam_z("potpars"),               //   accelerations        
		 getparam_z("potfile")) : 0;          // ELSE: no potential     
  WD99disk DM(getiparam("nbpero"),
		getdparam("r_d"),
		getdparam("Sig_0"),
		getdparam("r_sig"),
		getdparam("Qmin"),
		getdparam("z_0"),
		getdparam("eps"),aex);
  snapshot shot(getdparam("time"),getiparam("nbody"), data);
  if(getdparam("Qmin"))    
    DM.sample(shot,getiparam("ni"),getbparam("q-ran"),Ran);
  else
    DM.coldsample(shot,getbparam("q-ran"),Ran);

  //----------------------------------------------------------------------------
  // 3. output of snapshot                                                      
  //----------------------------------------------------------------------------
  nemo_out out(getparam("out"));
  shot.write_nemo(out,data);
}
