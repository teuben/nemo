// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// randomize_azimuth.cc                                                        |
//                                                                             |
// Copyright (C) 2004, 2005 Walter Dehnen                                      |
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
#include <public/defman.h>
#include <public/random.h>
#include <ctime>
#include <cmath>

namespace {
  using namespace falcON;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class randomize_azimuth                                                  //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class randomize_azimuth : public manipulator {
  private:
    const   Random3 Ran;
    const   int     K, N;
    //--------------------------------------------------------------------------
    static void rotate(vect&x, real const&c, real const&s) {
      register real 
	t0 = c*x[0] + s*x[1];
      x[1] = c*x[1] - s*x[0];
      x[0] = t0;
    }
    //--------------------------------------------------------------------------
  public:
    const char* name    () const { return "randomize_azimuth"; }
    const char* describe() const {
      if(N>0) return message("randomize azimuth for bodies [%d:%d]",K,N);
      if(K>0) return message("randomize azimuth for bodies [i>=%d]",K);
      else    return "randomize azimuth for all bodies";
    }
    //--------------------------------------------------------------------------
    fieldset need   () const { return fieldset::o; }
    fieldset provide() const { return fieldset::o; }
    fieldset change () const { return fieldset::vectors; }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*) const;
    //--------------------------------------------------------------------------
    randomize_azimuth(const double*pars,
		      int          npar,
		      const char  *file) :
      Ran (npar>0? long(pars[0]) : long(time(0)) ),
      K   (npar>1? int (pars[1]) : 0),
      N   (npar>2? int (pars[2]) : 0)
    {
      if(npar<3 && nemo_debug(1) || nemo_debug(2))
	std::cerr<<
	  " Manipulator \"randomize_azimuth\":\n"
	  " randomizes (using seed=par[0], default: secs since 1970)"
	  " azimuth in xvaw"
	  " for bodies with k in (par[1],par[2]), default: all\n";
      if(K>N)
	error(" Manipulator \"randomize_azimuth\":"
	      " par[0]=%d >= par[1]=%d",K,N);
      if(npar>3 && nemo_debug(1))
	::warning(" Manipulator \"randomize_azimuth\":"
		  " skipping parameters beyond 3\n");
    }
    //--------------------------------------------------------------------------
    ~randomize_azimuth() {}
  };
  //////////////////////////////////////////////////////////////////////////////
  bool randomize_azimuth::manipulate(const snapshot*S) const {
    if(S->have_not(fieldset::vectors)) return false;
    const body bMin = S->bodyNo(K>0? K : 0);
    const body bMax = S->bodyNo(N<=0?            S->N_bodies() : 
				N<S->N_bodies()? N : S->N_bodies() );
    const bool 
      do_x = S->have(fieldbit::x),
      do_v = S->have(fieldbit::v),
      do_a = S->have(fieldbit::a),
      do_j = S->have(fieldbit::j),
      do_w = S->have(fieldbit::w);
    for(body b(bMin); b!=bMax; ++b) {
      double phi = TPi * Ran();
      register real
	c = std::cos(phi),
	s = std::sin(phi);
      if(do_x) rotate(b.pos (), c,s);
      if(do_v) rotate(b.vel (), c,s);
      if(do_a) rotate(b.acc (), c,s);
      if(do_j) rotate(b.jerk(), c,s);
      if(do_w) rotate(b.vprd(), c,s);
    }
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} // namespace {

__DEF__MAN(randomize_azimuth)
