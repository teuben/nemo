// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// GalPot.cc                                                                   |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002                                               |
// e-mail:   dehnen@mpia.de                                                    |
// address:  Max-Planck Institut fuer Astronomie,                              |
//           Koenigstuhl 17, D-69117 Heidelberg, Germany                       |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// defines Dehnen & Binney's (1998) Galaxy potential as NEMO potential         |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// The units differ from those used in class GalPot: we have                   |
//                                                                             |
//  quantity              | here          | GalPot        | unit               |
// -----------------------+---------------+---------------+--------------      |
// unit of length         | 1             | 1             | kpc                |
// unit of time           | 1 E+9         | 1 E+6         | yr                 |
// unit of mass           | 2.2228847 E+5 | 1             | solar masses       |
// unit of velocity       | 0.9777753     | 9.777753  E+2 | km/s               |
// unit of potential      | 0.95604457    | 9.5604457 E+5 | (km/s)^2           |
// unit of acceleration   | 0.9777753 E-9 | 0.977753  E-6 | km/s/yr            |
//                                                                             |
// The implied size of Newton's constant G are                                 |
// G = 1                 here                                                  |
// G = 4.49865897 E-12   in GalPot units                                       |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// Versions                                                                    |
// 0.0   06-jun-2002    created                                           WD   |
// 0.1   28-dec-2002    replaced NEMO's <stdinc.h> for gcc3              PJT   |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <iostream>
#include <fstream>
using namespace std;		
//=============================================================================#
void dummy()
{
  std::cerr<<" this is a dummy\n";
}
//=============================================================================#
// declare externally linkable C routines                                      |
//=============================================================================#
extern "C" {
  //  #include <stdinc.h>     this code really only need warning/error
  //  gcc3 using string, which NEMO also has (but different)
  extern void warning(char *, ...);
  extern void error(char *, ...);
  extern void nemo_dprintf(int, char *, ...);
  void inipotential    (int*, double*, char *);
  void potential_double(int*, double*, double*, double*, double*);
  void potential_float (int*, float *, float *, float *, float *);
}
//=============================================================================#
// define C++ implementation of Galaxy potential                               |
//=============================================================================#
#include "inc/GalPot.cc"
//=============================================================================#
// now define externally linkable C routines                                   |
//=============================================================================#
static GalaxyPotential *POT = 0;
//------------------------------------------------------------------------------
void inipotential(int *npar, double *par, char *file) {
  double omega = (*npar>0)? par[0] : 0.;
  if (*npar>1) warning("Skipped potential parameters for GalPot beyond 1");
  if (file==0) 
    error("Need potfile to initialize GalPot, please consult $NEMODAT/GalPot");
  ifstream from(file);
  if (!from)   error  ("Cannot open potfile in initialization of GalPot");
  if(POT) delete POT;
  POT = new GalaxyPotential(from);
  from.close();
  nemo_dprintf (1,"INI_POTENTIAL Dehnen & Binney (1998) Galaxy potential\n");
  nemo_dprintf (1,"  Parameters read from file %s\n",file);
}
//------------------------------------------------------------------------------
void potential_double (int    *NDIM,
		       double *X,
		       double *F,
		       double *P,
		       double *T)
{
  if(POT) {
    register double fR,fz,R=hypot(X[0],X[1]);
    *P   = 1.e6 * (*POT)(R,(*NDIM>2)? X[2]:0., fR, fz);
    fR  /= R;
    F[0] =-1.e6 * fR * X[0];
    F[1] =-1.e6 * fR * X[1];
    if(*NDIM>2) 
      F[2] =-1.e6 * fz;
  }
  else
    error("potential GalPot not initialized");
//     cerr<<"potential GalPot not initialized\n";
}
//------------------------------------------------------------------------------
void potential_float  (int    *NDIM,
		       float  *X,
		       float  *F,
		       float  *P,
		       float  *T)
{
  if(POT) {
    register double fR,fz,R=hypot(X[0],X[1]);
    *P   = 1.e6 * (*POT)(R,(*NDIM>2)? X[2]:0., fR, fz);
    fR  /= R;
    F[0] =-1.e6 * fR * X[0];
    F[1] =-1.e6 * fR * X[1];
    if(*NDIM>2) 
      F[2] =-1.e6 * fz;
  }
  else
    error("potential GalPot not initialized");
//     cerr<<"potential GalPot not initialized\n";
}
//------------------------------------------------------------------------------
