// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// GalPot.cc                                                                   |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
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
// 0.2   23-aug-2004    added acceleration support using defacc.h         WD   |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <iostream>
#include <fstream>
using namespace std;
#define POT_DEF
#include <defacc.h>
////////////////////////////////////////////////////////////////////////////////
#include "inc/GalPot.cc"
////////////////////////////////////////////////////////////////////////////////
namespace {
  using GalPot::GalaxyPotential;

  class GalaxyFile {
  protected:
    ifstream from;
    GalaxyFile(const char*file)
    {
      if(file==0 || file[0]==0) 
	error("Need data file to initialize GalPot");
      from.open(file);
      if(!from.is_open())
	error("GalPot: cannot open file \"%s\"",file);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  class GalaxyPot : 
    private GalaxyFile,
    private GalaxyPotential
  {
  public:
    //--------------------------------------------------------------------------
    static const char* name() { return "GalPot"; }
    bool NeedMass() const { return false; }
    bool NeedVels() const { return false; }
    //--------------------------------------------------------------------------
    GalaxyPot(const double*pars,
	      int          npar,
	      const char  *file)
      : GalaxyFile     ( file ),
	GalaxyPotential( from )
    {
      double omega = (npar>0)? pars[0] : 0.;
      if (npar>1) warning("Skipped potential parameters for GalPot beyond 1");
      from.close();
    }
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar>
    void set_time(double       ,
		  int          ,
		  const scalar*,
		  const scalar*,
		  const scalar*) const {}
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar>
    void acc(const scalar*,
	     const scalar*X,
	     const scalar*,
	     scalar      &P,
	     scalar      *A) const
    {
      register double fR,fz,R=hypot(X[0],X[1]);
      if(NDIM > 2) {
	P    = 1.e6 * GalaxyPotential::operator()(R, X[2], fR, fz);
	fR  *=-1.e6/R;
	A[0] = fR * X[0];
	A[1] = fR * X[1];
	A[2] =-1.e6 * fz;
      } else {
	P    = 1.e6 * GalaxyPotential::operator()(R, 0.,   fR, fz);
	fR  *=-1.e6/R;
	A[0] = fR * X[0];
	A[1] = fR * X[1];
      }
    }
  };
} // namespace {
//------------------------------------------------------------------------------

__DEF__ACC(GalaxyPot)
__DEF__POT(GalaxyPot)

//------------------------------------------------------------------------------
