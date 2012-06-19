// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// LogPot.cc                                                                   |
//                                                                             |
// Copyright (C) 2011,2012 Walter Dehnen                                       |
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
// Versions                                                                    |
// 0.0   23-jun-2011    created                                             WD |
// 0.1   19-jun-2012    avoid warning about uninitialised vars              WD |
//-----------------------------------------------------------------------------+
#include <iostream>
#include <fstream>
#include <cmath>
#define POT_DEF
#include <defacc.h>
////////////////////////////////////////////////////////////////////////////////
using namespace std;
////////////////////////////////////////////////////////////////////////////////
namespace {
  //////////////////////////////////////////////////////////////////////////////
  class LogPot
  {
    double VQ;    ///< Vcirc^2
    double RQ;    ///< Rcore^2
    double iQ[3]; ///< inverse axis ratios squared
  public:
    //--------------------------------------------------------------------------
    static const char* name() { return "LogPot"; }
    bool NeedMass() const { return false; }
    bool NeedVels() const { return false; }
    //--------------------------------------------------------------------------
    LogPot(const double*pars, int npar, const char*)
      : VQ ( npar>1? pars[1] : 1),
      	RQ ( npar>2? pars[2] : 0)
    {
      if(VQ<=0) error("LogPot: circular speed <=0\n");
      if(RQ< 0) error("LogPot: core radius <0\n");
      double CA = npar>3? pars[3] : 1;
      double BA = npar>4? pars[4] : 1;
      if(CA<=0 || CA> 1) error("LogPot: major axis ratio = %g",CA);
      if(BA<=0) error("LogPot: minor axis ratio = %g",BA);
      if(BA>CA) error("LogPot: minor axis ratio = %g > %g = major axis ratio",
		      BA,CA);
      iQ[0] = std::pow(CA*BA,2./3.);
      iQ[1] = iQ[0] / (BA*BA);
      iQ[2] = iQ[0] / (CA*CA);
      if((npar<3 && nemo_debug(1)) || nemo_debug(2) )
	std::cerr<<
	"### LogPot: external logarithmic potential\n\n"
	"                  2      2  x^2   y^2   z^2\n"
	"      Phi = 0.5 vc log(rc + --- + --- + ---)\n"
	"                            a^2   b^2   c^2\n\n"
	"      par[0] ignored\n"
	"      par[1] vc : circular speed (default: 1)\n"
	"      par[2] rc : core radius (default: 0)\n"
	"      par[3] c/a: minor to major axis ratio (default: 1)\n"
	"      par[4] b/a: intermediate to major axis ratio (default: 1)\n";
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
      if(NDIM == 3) {
	double sQ = RQ;
	A[0] = iQ[0]*X[0]; sQ += A[0]*X[0];
	A[1] = iQ[1]*X[1]; sQ += A[1]*X[1];
	A[2] = iQ[2]*X[2]; sQ += A[2]*X[2];
	P    = 0.5 * VQ * log(sQ);
	sQ   =-VQ/sQ;
	A[0]*= sQ;
	A[1]*= sQ;
	A[2]*= sQ;
      } else
	P = A[0] = A[1] = A[2] = 0;
	error("LogPot: wrong number (%d) of dimensions, only allow 3D\n",NDIM);
    }
  };
} // namespace {
//------------------------------------------------------------------------------
__DEF__ACC(LogPot)
__DEF__POT(LogPot)
//------------------------------------------------------------------------------
