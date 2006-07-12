//-----------------------------------------------------------------------------+
//                                                                             |
// Dehnen.cc                                                                   |
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
////////////////////////////////////////////////////////////////////////////////
namespace {
  class Dehnen {
    enum { hernquist, jaffe, general } model;
    double g, g2, Pf, GM, Rs, Eq;
  public:
    static const char* name() { return "Dehnen"; }
    Dehnen(const double*pars,
	   int          npar,
	   const char  *file)
    {
      if(npar < 4)
	warning("%s: recognizing 5 parameters:\n"
		" omega        pattern speed (ignored)           [0]\n"
		" gamma        inner density slope               [1]\n"
		" G*M          mass;                             [1]\n"
		" a            scale radius;                     [1]\n"
		" e            softening or core radius;         [0]\n"
		"the potential is given by\n\n"
		"             G*M    1           r  2-g\n"
                "    Phi = - ----- ----- [ 1 - (---)   ]\n"
		"              a    2-g         r+a\n\n"
		"with\n\n"
		"    r   = sqrt(x^2+e^2)\n\n",name());
      if(file)
	warning("%s: file \"%s\" ignored",name(),file);
      double o,e,m;
      o     = npar>0? pars[0] : 0.;
      g     = npar>1? pars[1] : 1.;
      m     = npar>2? pars[2] : 1.;
      Rs    = npar>3? pars[3] : 1.;
      e     = npar>4? pars[4] : 0.;
      Eq    = e*e;
      g2    = 2.-g;
      GM    = -m;
      Pf    = g==1? GM        : g==2? m/Rs  : -m/Rs/g2;
      model = g==1? hernquist : g==2? jaffe : general;
      if(npar>5) warning("%s: skipped parameters beyond 5",name());
      nemo_dprintf (1,
		    "initializing %s\n"
		    " parameters : pattern speed = %f (ignored)\n"
		    "              gamma         = %f\n"
		    "              mass          = %f\n"
		    "              scalelength   = %f\n"
		    "              core-radius   = %f\n",
		    name(),o,g,m,Rs,e);
    }
    template<typename scalar>
    void potacc(scalar const&Rq,
		scalar      &P,
		scalar      &T) const
    {
      register scalar R=std::sqrt(Eq+Rq);
      T = 1/(R+Rs);
      switch(model) {
      case hernquist: {
	// gamma=1: softened Hernquist model
	P  = Pf * T;
	T *= P/R;
      } break;
      case jaffe: {
	// gamma==2: softened Jaffe model
	P  = Pf * std::log(R*T);
	T *= GM/(R*R);
      } break;
      default: {
	// general gamma!=2
	register scalar Q = std::pow(double(R*T),g2);
	P  = Pf*(1-Q);
	T *= GM*Q/(R*R);
      } break;
      }
    }
  };
}
//------------------------------------------------------------------------------

__DEF__ACC(SphericalPot<Dehnen>)
__DEF__POT(SphericalPot<Dehnen>)

//------------------------------------------------------------------------------
