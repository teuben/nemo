//-----------------------------------------------------------------------------+
//                                                                             |
// NFW.cc                                                                      |
//                                                                             |
// Copyright (C) 2004 Walter Dehnen                                            |
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
// defines as NEMO potential and acceleration field                            |
//                                                                             |
//             M0                           ln(1+r/a)                          |
// rho(r) = ---------    Phi(r) = - 4 Pi M0 ---------                          |
//          r (r+a)^2                           r                              |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// Versions                                                                    |
// 0.0   14-aug-2002    created                                           WD   |
// 0.1   18-nov-2002    converted from C++ to C                           WD   |
// 1.0   23-aug-2004    converted back to C++, implementing acceleration  WD   |
//                                                                             |
//-----------------------------------------------------------------------------+
#define POT_DEF
#include <cmath>
#include <acc/defacc.h>
////////////////////////////////////////////////////////////////////////////////
namespace {
  class NFWPot {
    double A,iA,Fac;
  public:
    static const char* name() { return "NFW"; }
    //--------------------------------------------------------------------------
    NFWPot(const double*pars,
	   int          npar,
	   const char  *file)
    {
      if(npar < 3)
	warning("%s: recognizing 3 parameters:\n"
		" omega        pattern speed (ignored)           [0]\n"
		" a            scale radius;                     [1]\n"
		" vc_max       maximum circular speeed           [1]\n"
		"the density is proportional to\n\n"
		"        1    \n"
                "    -------- \n" 
		"    r (r+a)^2\n\n",name());
      if(file)
	warning("%s: file \"%s\" ignored",name(),file);
      double
	o   = npar>0? pars[0] : 0.;
      A     = npar>1? pars[1] : 1.;
      iA    = 1./A;
      double
	v   = npar>2? pars[2] : 1.;
      Fac   = A * v*v / 0.2162165954;
      if(npar>3) warning("%s: skipped parameters beyond 3",name());
      nemo_dprintf (1,"initializing %s\n",name());
      nemo_dprintf (1," parameters : pattern speed = %f (ignored)\n",o);
      nemo_dprintf (1,"              scale radius  = %f\n",A);
      nemo_dprintf (1,"              vc_max        = %f\n",v);
    }
    //--------------------------------------------------------------------------
    template<typename scalar>
    void potacc(scalar const&Rq,
		scalar      &P,
		scalar      &fR) const
    {
      register scalar
	R = std::sqrt(Rq),
	iR= 1/R;
      fR  = std::log(1+iA*R) * iR;
      P   =-Fac*fR;
      fR *= iR;
      fR -= iR/(R+A);
      fR *=-Fac*iR;
    }
  }; // class NFW {
} // namespace {
//------------------------------------------------------------------------------

__DEF__ACC(SphericalPot<NFWPot>)
__DEF__POT(SphericalPot<NFWPot>)

//------------------------------------------------------------------------------
