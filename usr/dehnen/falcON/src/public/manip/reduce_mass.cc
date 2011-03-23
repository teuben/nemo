// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   src/public/manip/reduce_mass.cc
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
// v 0.0    23/03/2011  WD created                                              
////////////////////////////////////////////////////////////////////////////////
#include <public/defman.h>
#include <acc/timer.h>

namespace falcON {
  namespace Manipulate {
    // /////////////////////////////////////////////////////////////////////////
    //
    // class reduce_mass
    //
    /// manipulator: reduce mass of all bodies @a in_subset()
    ///
    /// This manipulator reduces the mass of all bodies @a in_subset()
    /// according to @f$ m(t) = m(t_0) [1-A(t)] @f$ with @f$ A(t) @f$ the
    /// 'adiabatic' timing function of acc/timer.h.
    ///
    /// \note The falcON force solver cannot tolerate particles with zero
    ///       mass. Therefore, we must either reduce the mass only to a
    ///       non-zero level, or remove the body when their mass hits zero.
    ///
    /// \warning The user must ensure that the subset is defined appropriately
    ///       so as not to change meaning when a particle is removed. For
    ///       example, using manipulator set_subset with filter "i<1" is
    ///       flawed: after particle i=0 is removed the next particle will be
    ///       removed and so on, until all particles are removed N
    ///       manipulations after t=t0+tau. In such a case, the number of
    ///       bodies in the subset after the removal of the first body is
    ///       still one, while we expect none, and a warning will be issued.
    ///
    /// Meaning of the parameters:\n
    /// par[0] : t0 : simulation time: begin of mass reduction (required)\n
    /// par[1] : tau: time scale of mass reduction (required)\n
    /// par[2] : m0 : don't reduce masses below this (default: 0)\n
    /// file   : not used\n
    ///
    /// Usage of pointers: none\n
    /// Usage of flags:    uses in_subset()\n
    ///
    // /////////////////////////////////////////////////////////////////////////
    class reduce_mass : public manipulator, private timer
    {
      const real   M0;                    ///< minimum mass for bodies
      bool         F1;                    ///< first ever manipulation?
      unsigned     NS;                    ///< # particles in subset
      mutable real MF;                    ///< previous mass-reduction factor
    public:
      const char*name() const
      { return "reduce_mass"; }
      const char*describe() const
      { return "reduces mass of bodies in subset"; }
      fieldset need() const 
      { return fieldset::f | fieldset::m; }
      fieldset provide() const
      { return fieldset::empty; }
      fieldset change() const
      { return fieldset::m; }
      /// ctor
      reduce_mass(const double*pars, int npar, const char*file)
	: timer(timer::adiabatic,
		npar>0? pars[0] : 0.0,
		npar>1? pars[1] : 1.0),
	  M0   (npar>2? pars[2] : zero),
	  F1   (1),
	  NS   (0),
	  MF   (one)
      {
	if(npar<2)
	  falcON_THROW("Manipulator \"%s\": "
		       "two parameters (t0, tau) are obligatory\n",name());
	if(M0<0)
	  falcON_THROW("Manipulator \"%s\": m_min (par[3]) < 0\n",name());
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
	if(F1) {
	  // first ever call: remember N_subset and perform sanity checks
	  const_cast<unsigned&>(NS) = S->N_subset();
	  const_cast<bool&>(F1)     = 0;
	  if(T0() < S->time())
	    falcON_Warning("Manipulator \"%s\": t0=%g < t_ini=%g\n",
			   name(),T0(),S->time());
	  if(A>zero) // this should never happen
	    falcON_Warning("Manipulator \"%s\": A=%g > 0 at t_ini=%g\n",
			   name(),A,S->time());
	    
	}
	// check for unexpected change of N_subset
	if(NS != S->N_subset()) {
	  falcON_Warning("Manipulator \"%s\": # bodies in subset changed "
			 "unexpectedly from %u to %u; "
			 "this suggests a problem with the subset definition\n",
			 name(),NS,S->N_subset());
	  const_cast<unsigned&>(NS) = S->N_subset();
	}
	if(S->N_subset() && A>zero) {
	  real     Mf = one-A;        // mass reduction from t=t0
	  real     fm = Mf/MF;        // mass reduction form last manipulation
	  MF          = Mf;           // remember for next manipulation
	  unsigned rem= 0;            // counter: bodies to be removed
	  LoopSubsetBodies(S,b) {
	    b.mass() *= fm;                      // adjust mass
	    if(M0 > zero && mass(b)<M0)          // but not below M0
	      b.mass() = M0;
	    if(mass(b)<= zero) {                 // zero mass?
	      b.flag_for_removal();              //   prepare for removal
	      ++rem;                             //   count
	    }
	  }
	  if(rem) {                              // remove zero-mass bodies
	    const_cast<snapshot*>(S)->remove();
	    const_cast<unsigned&>(NS) -= rem;    //   keep count of subset
	  }
	}
	return false;
      }
    };
} }
////////////////////////////////////////////////////////////////////////////////
__DEF__MAN(falcON::Manipulate::reduce_mass);
