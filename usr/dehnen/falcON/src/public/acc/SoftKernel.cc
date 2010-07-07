//------------------------------------------------------------------------------
//
// SoftKernel.cc
//
// Copyright (C) 2010 Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 675
// Mass Ave, Cambridge, MA 02139, USA.
//
//------------------------------------------------------------------------------
// Versions
// 0.1    07/07/2010  WD  created.
//------------------------------------------------------------------------------
#include <cmath>             // for std::sqrt()
#include <public/default.h>  // $FALCON/inc/public/default.h
#define POT_DEF
#include <defacc.h>          // $NEMOINC/defacc.h
//
namespace {
  using namespace falcON;
  class SoftKernel {
    const double    GM,EQ,HQ,QQ;
    const kern_type K;
    //
  public:
    static const char* name() { return "SoftKernel"; }
    //
    SoftKernel(const double*pars,
	       int          npar,
	       const char  *file)
      : GM (     - (npar>1? pars[1] : 1) ),
	EQ ( square(npar>2? pars[2] : 1) ),
	HQ ( 0.5*EQ ),
	QQ ( 0.5*HQ ),
	K  ( npar>3? kern(indx(pars[3])) : Default::kernel )
    {
      if((npar<4 && nemo_debug(1)) || nemo_debug(2) )
	DebugInfo("%s: recognising 4 parameters:\n"
		  " omega    pattern speed (ignored)   [0]\n"
		  " G*M      mass;                     [1]\n"
		  " eps      scale radius              [1]\n"
		  " kern     kernel type (0,1,2,3)     [" 
		  falcON_KERNEL_TEXT "]\n\n"
		  "the potential is given by that of a softenen massive\n"
		  "particle with mass M and softening length eps using the\n"
		  "softening kernel kern.\n",name());
      if(file && file[0])
	falcON_WarningN("%s: file \"%s\" ignored",name(),file);
    }
    //
    template<typename scalar>
    void potacc(scalar const&Rq, scalar&P, scalar&T) const
    {
      double x  = 1/(Rq+EQ);
      double D0 = GM*std::sqrt(x);
      double D1 = D0*x;
      switch(K) {
      case p1: {
	double D2 = 3*x*D1;
	D0 += HQ*D1;
	D1 += HQ*D2;
      } break;
      case p2: {
	double D2 = 3*x*D1;
	double D3 = 5*x*D2;
	D0 += HQ*(D1+HQ*D2);
	D1 += HQ*(D2+HQ*D3);
      } break;
      case p3: {
	double D2 = 3*x*D1;
	double D3 = 5*x*D2;
	double D4 = 7*x*D3;
	D0 += HQ*(D1+QQ*(D2+HQ*D3));
	D1 += HQ*(D2+QQ*(D3+HQ*D4));
      } break;
      default: break;
      }
      P = D0;
      T = D1;
    }
  };
}
//------------------------------------------------------------------------------
    
__DEF__ACC(SphericalPot<SoftKernel>)
__DEF__POT(SphericalPot<SoftKernel>)

//------------------------------------------------------------------------------

