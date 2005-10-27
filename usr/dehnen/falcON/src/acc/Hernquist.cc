//-----------------------------------------------------------------------------+
//                                                                             |
// Hernquist.cc                                                                |
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
#define POT_DEF
#include <cmath>
#include <defacc.h>
//=============================================================================#
// define C++ implementation of potential                                      |
//=============================================================================#
namespace {
  class Hernquist {
    double Rs, Eq, GM;
  public:
    static const char* name() { return "Hernquist"; }
    Hernquist(const double*pars,
	      int          npar,
	      const char  *file)
    {
      if(npar < 3)
	warning("%s: recognizing 4 parameters:\n"
		"     omega -- pattern speed (ignored)\n"
		"     G*M   -- mass; defaults to 1\n"
		"     a     -- scale radius; defaults to 1\n"
		"     e     -- softening or core radius; default to 0\n"
		" the potential is given by\n\n"
		"                    G M\n"
		"    Phi = - ----------------- .\n"
		"            sqrt(r^2+e^2) + a\n\n",name());
      if(file)
	warning("%s: file \"%s\" ignored",name(),file);
      double o,e,m;
      o  = npar>0? pars[0] : 0.;
      m  = npar>1? pars[1] : 1.;
      Rs = npar>2? pars[2] : 1.;
      e  = npar>3? pars[3] : 0.;
      GM =-m;
      Eq = e*e;
      if(npar>4) warning("%s: skipped parameters beyond 4",name());
      nemo_dprintf (1,
		    "initializing %s\n"
		    " parameters : pattern speed = %f (ignored)\n"
		    "              mass          = %f\n"
		    "              scalelength   = %f\n"
		    "              core-radius   = %f\n",
		    name(),o,m,Rs,e);
    }    
    template<typename scalar>
    void potacc(scalar const&Rq,
		scalar      &P,
		scalar      &T) const
    {
      register scalar R = std::sqrt(Eq+Rq);
      T  = 1/(R+Rs);
      P  = GM * T;
      T *= P/R;
    }
  };
}
//------------------------------------------------------------------------------

__DEF__ACC(SphericalPot<Hernquist>)
__DEF__POT(SphericalPot<Hernquist>)

//------------------------------------------------------------------------------

  
