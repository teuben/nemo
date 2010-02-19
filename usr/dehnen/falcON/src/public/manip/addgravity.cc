// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   src/public/manip/addgravity.cc
///
/// \author Walter Dehnen
/// \date   2010
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2010 Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
////////////////////////////////////////////////////////////////////////////////
//
// history:                                                                     
//
// v 0.0    19/02/2010  WD created
////////////////////////////////////////////////////////////////////////////////
#include <public/defman.h>
#include <forces.h>

namespace falcON {
namespace Manipulate {

  // ///////////////////////////////////////////////////////////////////////////
  //
  // class addgravity
  //
  /// manipulator: adds self-gravity (acceleration and potential)
  ///
  /// Using the falcON force solver, this manipulator computes for forces on
  /// all bodies in subset (default: all) from all particles.
  ///
  /// Meaning of the parameters:\n
  /// par[0]: eps    gravitational softening length (<0: use individual e_i)\n
  /// par[1]: theta  tolerance parameter at M=M_tot (default: 0.6)\n
  /// par[2]: kernel softening kernel (default: 1->P1)\n
  /// par[3]: Grav   Newton's constant of gravity (default: 1)\n
#ifdef falcON_PROPER
  /// par[4]: fsink theta_sink/theta (default: 1)\n
#endif
  ///
  /// Usage of pointers: none\n
  /// Usage of flags:    uses in_subset()\n
  ///
  // ///////////////////////////////////////////////////////////////////////////
  class addgravity : public manipulator {
  private:
#ifdef falcON_PROPER
    const static int MaxNumPar = 5;
#else
    const static int MaxNumPar = 4;
#endif
    real            Eps, The;
    kern_type       Kern;
    real            Grav, fsnk;
    mutable forces *FALC;
  public:
    const char* name    () const { return "addgravity"; }
    const char* describe() const
    { return message("adds self-gravity using "
		     "eps=%g, theta=%g, kernel=%s, G=%g",
		     Eps,The,falcON::describe(Kern),Grav);
    }
    //
    fieldset need   () const
    { return fieldset::m | fieldset::x | (Eps<0?fieldset::e:fieldset::empty); }
    fieldset provide() const { return fieldset::a | fieldset::p; }
    fieldset change () const { return fieldset::empty; }
    //
    bool manipulate(const snapshot*) const;
    //
    addgravity(const double*pars, int npar, const char*) falcON_THROWING
    : Eps ( npar>0? pars[0] : 0.05 ),
      The ( npar>1? pars[1] : Default::theta ),
      Kern( npar>2? kern_type(int(pars[2])) : Default::kernel ),
      Grav( npar>3? pars[3] : one ),
      fsnk( 
#ifdef falcON_PROPER
	    npar>4? pars[4] :
#endif
	    one ),
      FALC( 0 )
    {
      if((npar==0 && debug(1)) || debug(2))
	std::cerr<<
	  " Manipulator \""<<name()<<"\":\n"
	  " adds N-body gravity for bodies in_subset() (default: all)"
	  " par[0]: # nearest neighbours used in density estimate (def: 16)\n"
	  " par[1]: delta time between estimation (which takes long; def: 0)\n"
	  "par[0]: gravitational softening length (<0: use individual e_i)\n"
	  "par[1]: tolerance parameter at M=M_tot (default: 0.6)\n"
	  "par[2]: softening kernel (default: 1->P1)\n"
	  "par[3]: Newton's constant of gravity (default: 1)\n"
#ifdef falcON_PROPER
	  "par[4]: theta_sink/theta (default: 1)\n"
#endif
	  ;
      if(npar > MaxNumPar)
	falcON_WarningN("Manipulator \"%s\": "
			"skipping parameters beyond %d\n",name(),MaxNumPar);

    }
    //
    ~addgravity()
    { if(FALC) falcON_DEL_O(FALC); FALC=0; }
  };
  //
  bool addgravity::manipulate(const snapshot*S) const
  {
    if(FALC && FALC != S->Forces()) { falcON_DEL_O(FALC); FALC=0; }
    if(S->Forces() == 0)
      FALC = new forces(S,Eps,The,Kern,Eps<0,Grav,theta_of_M,fsnk);
    forces*Forces = const_cast<forces*>(S->Forces());
    Forces->grow();
    bool all = true;
    if(S->have(fieldbit::f)) {
      unsigned Nsub=0u;
      LoopAllBodies(S,b)
	if(in_subset(b)) {
	  ++Nsub;
	  b.flag_as_active();
	} else
	  b.unflag_active();
      all = Nsub == S->N_bodies();
    }
    Forces->approximate_gravity(true,all);
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} }

__DEF__MAN(falcON::Manipulate::addgravity)
    
