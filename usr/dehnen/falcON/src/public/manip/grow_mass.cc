// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   src/public/manip/grow_mass.cc
///
/// \author Walter Dehnen
///
/// \date   2011
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Walter Dehnen
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
////////////////////////////////////////////////////////////////////////////////
//
// history:
//
// v 0.0    15/06/2011  WD created
// v 0.0.0  16/06/2011  WD debugged
////////////////////////////////////////////////////////////////////////////////
#include <public/defman.h>
#include <acc/timer.h>

namespace falcON {
  namespace Manipulate {
    // /////////////////////////////////////////////////////////////////////////
    //
    // class grow_mass
    //
    /// manipulator: grow mass of all bodies @a in_subset()
    ///
    /// This manipulator grows the mass of all bodies @a in_subset() according
    /// to @f$ m(t) = m(t_0) * [1 + (F-1) * A(t)] @f$ with @f$ A(t) @f$ the
    /// 'adiabatic' timing function of acc/timer.h.
    ///
    /// Meaning of the parameters:\n
    /// par[0] : t0 : simulation time: begin of mass reduction (required)\n
    /// par[1] : tau: time scale of mass reduction (required)\n
    /// par[2] : F  : factor by which to grow masses (required)\n
    /// file   : not used\n
    ///
    /// Usage of pointers: none\n
    /// Usage of flags:    uses in_subset()\n
    ///
    // /////////////////////////////////////////////////////////////////////////
    class grow_mass : public manipulator, private timer
    {
      const real   F1;                    ///< F-1: mass growth factor
      bool         FS;                    ///< first ever manipulation?
      mutable real MF;                    ///< previous mass-reduction factor
    public:
      const char*name() const
      { return "grow_mass"; }
      const char*describe() const
      { return "grows mass of bodies in subset"; }
      fieldset need() const 
      { return fieldset::f | fieldset::m; }
      fieldset provide() const
      { return fieldset::empty; }
      fieldset change() const
      { return fieldset::m; }
      /// ctor
      grow_mass(const double*pars, int npar, const char*file)
	: timer(timer::adiabatic,
		npar>0? pars[0] : 0.0,
		npar>1? pars[1] : 1.0),
	  F1   (npar>2? pars[2]-1 : 0.0),
	  FS   (1),
	  MF   (one)
      {
	if(npar<3)
	  falcON_THROW("Manipulator \"%s\": "
		       "three parameters (t0, tau, F) are obligatory\n",name());
	if(F1<=0)
	  falcON_THROW("Manipulator \"%s\": m_min (par[3]) <= 1\n",name());
	if(file && file[0])
	  falcON_Warning("Manipulator \"%s\": file '%s' given but not used\n",
			 name(),file);
	if(npar>3)
	  falcON_Warning("Manipulator \"%s\": ignoring parameters beyond 3\n",
			 name());
      }
      /// manipulation
      bool manipulate(const snapshot*S) const
      {
	// adiabatic factor
	real A = timer::operator()(S->time());
	if(FS) {
	  // first ever call: perform sanity checks
	  const_cast<bool&>(FS)     = 0;
	  if(T0() < S->time())
	    falcON_Warning("Manipulator \"%s\": t0=%g < t_ini=%g\n",
			   name(),T0(),S->time());
	  if(A>zero) // this should not happen
	    falcON_Warning("Manipulator \"%s\": A=%g > 0 at t_ini=%g\n",
			   name(),A,S->time());
	}
	if(S->N_subset() && A>zero) {
	  real     Mf = one+F1*A;     // mass reduction from t=t0
	  real     fm = Mf/MF;        // mass reduction form last manipulation
	  MF          = Mf;           // remember for next manipulation
	  LoopSubsetBodies(S,b)
	    b.mass() *= fm;           // adjust mass
	}
	return false;
      }
    };
} }
////////////////////////////////////////////////////////////////////////////////
__DEF__MAN(falcON::Manipulate::grow_mass);
