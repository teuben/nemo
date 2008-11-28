// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// DiscPot.cc                                                                  |
//                                                                             |
// Copyright (C) 2007 Walter Dehnen                                            |
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
// The gravitational potential of a disc, using GalPot.                        |
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
// 0.0   18-sep-2007    created                                             WD |
// 0.1   19-sep-2007    small change in debug info                          WD |
// 0.2   12-jun-2008    "using namespace std" behind "#include <defacc.h>"  WD |
//                                                                             |
//-----------------------------------------------------------------------------+
#include <iostream>
#include <fstream>
#define POT_DEF
#include <defacc.h>
////////////////////////////////////////////////////////////////////////////////
using namespace std;
#include "acc/GalPot.cc"
////////////////////////////////////////////////////////////////////////////////
namespace {
  using GalPot::GalaxyPotential;
  //////////////////////////////////////////////////////////////////////////////
  class DiscPot : 
    private GalPot::DiskPar,
    private GalPot::GalaxyPotential
  {
  public:
    //--------------------------------------------------------------------------
    static const char* name() { return "DiscPot"; }
    bool NeedMass() const { return false; }
    bool NeedVels() const { return false; }
    //--------------------------------------------------------------------------
    DiscPot(const double*pars,
	    int          npar,
	    const char  *file)
      : GalPot::DiskPar         ( npar>1? pars[1]:1.,   // Sigma_0
				  npar>2? pars[2]:1.,   // R_d
				  npar>3? pars[3]:0.1,  // z_d
				  npar>4? pars[4]:0.,   // R_0
				  npar>5? pars[5]:0.),  // eps
	GalPot::GalaxyPotential ( 1, this, 0, 0 )
    {
      if((npar<5 && nemo_debug(1)) || nemo_debug(2) )
	std::cerr<<
	  "### nemo debug info:\n"
	  " external potential \"DiscPot\" recognizes 5 parameters:\n"
	  "   omega   pattern speed (ignored)\n"
	  "   Sig_0   central surface density   [1]\n"
	  "   R_d     disc scale length         [1]\n"
	  "   z_d     disc scale height         [0.1]\n"
	  "   R_0     radius of central hole    [0]\n"
	  "   eps     cosine modulation term    [0]\n"
	  " The disc density is given by rho(R,z)=Sigma(R)*h(z) with\n\n"
	  "   Sigma(R) = Sig_0 exp(-R_0/R-R/R_d+eps*cos[R/R_d])\n\n"
	  "          exp(-z/z_d) / (2 z_d)          if z_d > 0\n"
	  "   h(z) = delta(z)                       if z_d = 0\n"
	  "          sech^2(z/2|z_d|) / (4|z_d|)    if z_d < 0\n\n";
      if (npar>6)
	warning("Skipped potential parameters for DiscPot beyond %d", npar);
      if(nemo_debug(2))
	std::cerr<<
	  "### nemo debug info:\n"
	  " external potential \"DiscPot\" initialized with:\n"
	  "   Sig_0 = "<<(*this)[0]<<"\n"
	  "   R_d   = "<<(*this)[1]<<"\n"
	  "   z_d   = "<<(*this)[2]<<"\n"
	  "   R_0   = "<<(*this)[3]<<"\n"
	  "   eps   = "<<(*this)[4]<<"\n";
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
	if(R) {
	  fR  *=-1.e6/R;
	  A[0] = fR * X[0];
	  A[1] = fR * X[1];
	} else {
	  A[0] = 0.;
	  A[1] = 0.;
	}
	A[2] =-1.e6 * fz;
      } else {
	P    = 1.e6 * GalaxyPotential::operator()(R, 0.,   fR, fz);
	if(R) {
	  fR  *=-1.e6/R;
	  A[0] = fR * X[0];
	  A[1] = fR * X[1];
	} else {
	  A[0] = 0.;
	  A[1] = 0.;
	}
      }
    }
  };
} // namespace {
//------------------------------------------------------------------------------

__DEF__ACC(DiscPot)
__DEF__POT(DiscPot)

//------------------------------------------------------------------------------
