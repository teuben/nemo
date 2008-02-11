//-----------------------------------------------------------------------------+
//                                                                             |
// Plummer.cc                                                                  |
//                                                                             |
// Copyright (C) 2004-2006 Walter Dehnen                                       |
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
// 0.1    09/08/2006  WD use $NEMOINC/defacc.h                                 |
//-----------------------------------------------------------------------------+
#define POT_DEF
#include <cmath>
#include <defacc.h> // $NEMOINC/defacc.h
////////////////////////////////////////////////////////////////////////////////
namespace {
  class Plummer {
    double GM, Aq;
  public:
    static const char* name() { return "Plummer"; }
    Plummer(const double*pars,
	    int          npar,
	    const char  *file)
    {
      if(npar < 3)
	warning("%s: recognizing 3 parameters:\n"
		" omega        pattern speed (ignored)           [0]\n"
		" G*M          mass;                             [1]\n"
		" a            scale radius;                     [1]\n"
		"\nthe potential is given by\n\n"
		"                 G*M\n"
                "    Phi = - ---------------\n"
		"            sqrt(a^2 + r^2)\n\n",name());
      if(file && file[0])
	warning("%s: file \"%s\" ignored",name(),file);
      double
	o   = npar>0? pars[0] : 0.,
	m   = npar>1? pars[1] : 1.,
	a   = npar>2? pars[2] : 1.;
      Aq    = a*a;
      GM    = -m;
      if(npar>3) warning("%s: skipped parameters beyond 3",name());
      nemo_dprintf (1,
		    "initializing %s\n"
		    " parameters : pattern speed = %f (ignored)\n"
		    "              mass          = %f\n"
		    "              scalelength   = %f\n",
		    name(),o,m,a);
    }
    template<typename scalar>
    void potacc(scalar const&Rq,
		scalar      &P,
		scalar      &T) const
    {
      T = 1/(Aq+Rq);
      P = GM * sqrt(T);
      T*= P;
    }
  };
}
//------------------------------------------------------------------------------

__DEF__ACC(SphericalPot<Plummer>)
__DEF__POT(SphericalPot<Plummer>)

//------------------------------------------------------------------------------
