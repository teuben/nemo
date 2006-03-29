// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// LMXBkick.cc                                                                 |
//                                                                             |
// Copyright (C) 2006 Jalpesh Sachania & Walter Dehnen                         |
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
  // class kick_binary                                                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class kick_binary : public manipulator {
  private:
    // constant data members
    const   Random3 Ran;
    const   real    vmin, vmax;
    const   double  rate, irate;
    const   int     Nmax;
    // non-const mutable data members
    mutable bool    first;
    mutable body    Bi, Bn;
    mutable double  tkick;
    //--------------------------------------------------------------------------
    // applies a kick with vkick in [vmin,vmax] and with random orientation
    void kick_velocity(vect &vel) const {
      real phi = Ran(zero,TPi);
      real cth = Ran(-one,one);
      real sth = sqrt(one-cth*cth);
      real vk  = Ran(vmin,vmax);
      vel[0]  += vk * sth * cos(phi);
      vel[1]  += vk * sth * sin(phi);
      vel[2]  += vk * cth;
    }
    //--------------------------------------------------------------------------
  public:
    const char* name    () const { return "kick_binary"; }
    const char* describe() const {
      return message("apply a randomly oriented velocity kick with |vkick|"
		     "uniform in [%f,%f]\n to bodies 0 to %d at a rate of %f"
		     "per time unit",vmin,vmax,Nmax-1,rate);
    }
    //--------------------------------------------------------------------------
    fieldset need   () const { return fieldset::v; }
    fieldset provide() const { return fieldset::o; }
    fieldset change () const { return fieldset::v; }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*) const;
    //--------------------------------------------------------------------------
    kick_binary(const double*pars,
		int          npar,
		const char  *file) :
      Ran   ( npar>0? long(pars[0]) : long(time(0)) ),
      vmin  ( npar>1? pars[1] : 0.),
      vmax  ( npar>2? pars[2] : 100.),
      rate  ( npar>3? pars[3] : 1.),
      irate ( one/rate ),
      Nmax  ( npar>4? int (pars[4]) : 0),
      first ( true ),
      Bi    ( bodies::bodyNil() ),
      Bn    ( bodies::bodyNil() )
    {
      if(npar<5 && nemo_debug(1) || nemo_debug(2))
	std::cerr<<
	  " Manipulator \"kick_binary\":\n"
	  " apply a randomly oriented velocity kick with |vkick| in"
	  " uniform [vmin,vmax] to bodies 0 to Nmax at a given rate\n"
	  " par[0] : seed for random number generator [time]\n"
	  " par[1] : vmin [0]\n"
	  " par[2] : vmax [100]\n"
	  " par[3] : Nmax [0]\n"
	  " par[4] : rate: bodies/time unit [1]\n";
      if(Nmax <= 0)
	warning("Manipulator \"kick_binary\": nobody to be kicked");
      if(rate <= zero)
	error("Manipulator \"kick_binary\": rate=%f <= 0",rate);
      if(vmin < zero)
	error("Manipulator \"kick_binary\": vmin=%f < 0",vmin);
      if(vmin > vmax)
	error("Manipulator \"kick_binary\": vmin=%f > vmax=%f",vmin,vmax);
      if(npar>5 && nemo_debug(1))
	::warning(" Manipulator \"kick_binary\":"
		  " skipping parameters beyond 5\n");
    }
    //--------------------------------------------------------------------------
    ~kick_binary() {}
  };
  //////////////////////////////////////////////////////////////////////////////
  bool kick_binary::manipulate(const snapshot*S) const {
    // if this is the first ever call, initialize some things
    if(first) {
      first = false;
      Bi    = S->begin_all_bodies();
      Bn    = body(Bi,Nmax);
      tkick = S->time();
    }
    // if some time has passed since last kick, we may want to kick again
    if(S->time() > tkick && Bi != Bn) {
      int k = int(rate * (S->time()-tkick));
      for(int i=0; i!=k && Bi!=Bn; ++i, ++Bi) {
	kick_velocity(Bi.vel());
	tkick += irate;
      }
    }
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} // namespace {

__DEF__MAN(kick_binary)
