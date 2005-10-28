// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// randomize_azimuth.cc                                                        |
//                                                                             |
// Copyright (C) 2005 Matthew Hanson                                           |
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
#include <public/Pi.h>
#include <ctime>
#include <cmath>

namespace {
  using namespace falcON;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class removesph                                                          //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class removesph : public manipulator {
  private:
    const double rin, rout, k;
    //--------------------------------------------------------------------------
  public:
    const char* name    () const { return "remove sph"; }
    const char* describe() const {
      return message("remove SPH bodies inside max{rin,k*h} or outside rout");
    }
    //--------------------------------------------------------------------------
    fieldset need   () const { return fieldset::o; }
    fieldset provide() const { return fieldset::o; }
    fieldset change () const { return fieldset::o; }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*) const;
    //--------------------------------------------------------------------------
    removesph(const double*pars,
	      int          npar,
	      const char  *file) :
      rin    (npar>0? double(pars[0]) : 0.05 ),
      rout   (npar>1? double(pars[1]) : 100),
      k      (npar>2? double(pars[2]) : 3.00)

    {
      if(npar>3 && nemo_debug(1))
	::warning(" Manipulator \"removesph\":"
		  " skipping parameters beyond 3\n"); 
    }
    //--------------------------------------------------------------------------
    ~removesph() {}
  };
  //////////////////////////////////////////////////////////////////////////////
  bool removesph::manipulate(const snapshot*S) const falcON_THROWING {
    const real rq_in = rin*rin, kq=k*k, rq_out=rout*rout;
    bool remove = false;
    LoopSPHBodies(S,B) {
      real rq = norm(pos(B));
      if(rq < rq_in || rq < kq * square(size(B)) || rq > rq_out) {
	B.flag_for_removal();
	remove = true;
      }
    }
    if(remove)
      const_cast<snapshot*>(S)->remove();
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} // namespace {

__DEF__MAN(removesph)
