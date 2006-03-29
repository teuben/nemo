//-----------------------------------------------------------------------------+
//                                                                             |
// Point.cc                                                                    |
//                                                                             |
// Copyright (C) 2006 Walter Dehnen                                            |
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
  class Point {
    static double gm(double x) { return x<0.? x : -x; }
    static double aq(double x) { return x*x; }
    const double GM, Aq, Hq, Qq;
    const int    Kl;
  public:
    static const char* name() { return "Point"; }
    Point(const double*pars,
	  int          npar,
	  const char  *file) :
      GM ( gm (npar>1? pars[1] : 1.) ),
      Aq ( aq (npar>2? pars[2] : 1.) ),
      Hq ( 0.5 * Aq ),
      Qq ( 0.5 * Hq ),
      Kl ( npar>3? int(pars[3]) : 0 )
    {
      if(npar < 3)
	warning("%s: recognizing 4 parameters:\n"
		" omega        pattern speed (ignored)           [0]\n"
		" G*M          mass;                             [1]\n"
		" a            scale radius;                     [1]\n"
		" kernel       use falcON's P_n kernel           [0]\n",
		name());
      if(file)
	warning("%s: file \"%s\" ignored",name(),file);
      if(npar>4) warning("%s: skipped parameters beyond 4",name());
      nemo_dprintf (1,
		    "initializing %s\n"
		    " parameters : mass          = %f\n"
		    "              scalelength   = %f\n",
		    "              kernel        = %d\n",
		    name(),GM,Aq,Kl);
    }
    template<typename scalar>
    void potacc(scalar const&Rq,
		scalar      &D0,
		scalar      &D1) const
    {
      scalar x = 1/(Aq+Rq);
      D0 = GM * sqrt(x);
      D1 = D0 * x;
      switch(Kl) {
      case 1:
	D0 += Hq*D1;
	D1 += Hq*D1*3*x;
	break;
      case 2: {
	scalar D2=3*x*D1, D3=5*x*D2;
	D0 += Hq*(D1+Hq*D2);
	D1 += Hq*(D2+Hq*D3);
      } break;
      case 3: {
	scalar D2=3*x*D1, D3=5*x*D2, D4=7*x*D3;;
	D0 += Hq*(D1+Qq*(D2+Hq*D3));
	D1 += Hq*(D2+Qq*(D3+Hq*D4));
      } break;
      case 0: default: break;
      }
    }
  };
}
//------------------------------------------------------------------------------

__DEF__ACC(SphericalPot<Point>)
__DEF__POT(SphericalPot<Point>)

//------------------------------------------------------------------------------
