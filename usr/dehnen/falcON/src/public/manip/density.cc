// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   src/public/manip/density.cc
///
/// \author Walter Dehnen
/// \date   2006-2010
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006-2008 Walter Dehnen
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
// v 0.0    27/04/2006  WD created
// v 0.1    07/07/2006  WD replacing bodyset with flags::ignore
// v 0.2    27/07/2006  WD made public
// v 0.3    04/09/2007  WD new neighbours.h
// v 0.4    06/11/2007  WD order of Ferrers kernel from param, K=32 default
// v 0.4.1  20/05/2008  WD renamed routine from neighbour.h
// v 0.4.2  11/06/2008  WD new DebugInfo and falcON_Warning
// v 0.4.3  25/09/2008  WD debugged (debug output only)
// v 0.5    05/11/2008  WD register time of manipulation under trho
// v 0.5.1  13/04/2010  WD removed use of initial_time()
// v 0.6    25/06/2010  WD manipulation if time=integer*step (or step==0)
////////////////////////////////////////////////////////////////////////////////
#include <public/defman.h>
#include <public/neighbours.h>
#include <ctime>

namespace falcON {
namespace Manipulate {
  //////////////////////////////////////////////////////////////////////////////
  //
  // auxiliary data and function
  //
  //////////////////////////////////////////////////////////////////////////////
  real FAC; ///< normalisation factor for kernel
  int  NF;  ///< order of Ferrers kernel
  void prepare(int n) {
    NF  = n;
    FAC = 0.75/Pi;
    for(n=1; n<=NF; ++n)
      FAC *= double(n+n+3)/double(n+n);
  }
  void SetDensity(const bodies       *B,
		  const OctTree::Leaf*L,
		  const Neighbour    *N,
		  int                 K)
  {
    real iHq = one/N[K-1].Q;
    real rho = zero;
    for(int k=0; k!=K-1; ++k)
      rho += scalar(N[k].L) * std::pow(one-iHq*N[k].Q,NF);
    rho *= FAC * std::pow(sqrt(iHq),3);
    B->rho(mybody(L)) = rho;
  }
  // ///////////////////////////////////////////////////////////////////////////
  //
  // class density
  //
  /// manipulator: estimates density using nearest neighbour count
  ///
  /// This manipulator finds for each body in_subset() (default: all, see
  /// set_subset) the Kth nearest neighbours (from the same set of bodies) and
  /// estimates the mass density, which is written into field 'r'.\n
  ///
  /// As this procedure is computationally expensive (more than computing
  /// gravity) it done only every T time units. The simulation time of the
  /// density estimation is put in 'trho' for the benefit of subsequent
  /// manipulators, such as densprof.\n
  ///
  /// Meaning of the parameters:\n
  /// par[0]: K: # nearest neighbours used in density estimate (def: 32)\n
  /// par[1]: N: order of Ferrers kernel: W = (1-r^2/h^2)^N (def: 1)\n
  /// par[2]: T: time between estimation (which takes long; def: 0)\n
  ///
  /// Usage of pointers: 'trho'\n
  /// Usage of flags:    uses in_subset()\n
  ///
  // ///////////////////////////////////////////////////////////////////////////
  class density : public manipulator {
  private:
    int             K;        ///< # neighbours
    int             N;        ///< # order of Ferrers kernel
    double          STEP;     ///< delta time between manipulations
    mutable double  TRHO;     ///< time of actual density.
    //--------------------------------------------------------------------------
  public:
    const char* name    () const { return "density"; }
    const char* describe() const {
      return message("estimates density using %dth nearest neighbour",K);
    }
    //--------------------------------------------------------------------------
    fieldset need   () const { return fieldset::m | fieldset::x; }
    fieldset provide() const { return fieldset::r; }
    fieldset change () const { return fieldset::empty; }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*) const;
    //--------------------------------------------------------------------------
    density(const double*pars,
	    int          npar,
	    const char  *file) falcON_THROWING
    : K    ( npar>0?     int(pars[0])    : 16 ),
      N    ( npar>1?     int(pars[1])    : 1  ),
      STEP ( npar>2?         pars[2]     : 0. ),
      TRHO ( 0.0 )
    {
      if((npar==0 && debug(1)) || debug(2))
	std::cerr<<
	  " Manipulator \""<<name()<<"\":\n"
	  " estimates density of bodies in_subset() (default: all)"
	  " and writes it to field r\n parameters:\n"
	  " par[0]: # nearest neighbours used in density estimate (def: 16)\n"
	  " par[1]: delta time between estimation (which takes long; def: 0)\n";
      if(K<=0)
	falcON_THROW("Manipulator \"%s\": "
		     "# neighbours (%d) must be positive",name(),K);
      if(file && file[0])
	falcON_WarningN("Manipulator \"%s\": "
			"file given but not used\n",name());
      if(npar>3)
	falcON_WarningN("Manipulator \"%s\": "
			"skipping parameters beyond 3\n",name());
    }
    //--------------------------------------------------------------------------
    ~density() {}
  };
  //////////////////////////////////////////////////////////////////////////////
  bool density::manipulate(const snapshot*S) const
  {
    // 0. preliminaries
    S->set_pointer(&TRHO,"trho");
    if(STEP && S->time() != STEP*int(S->time()/STEP))
      return false;
    TRHO = S->time();
    // 1. establish tree
    clock_t CPU0 = clock();
    flags FLAG   = flags::empty;
    // 1.1 if limited range of bodies, mark them
    if(S->N_bodies() != S->N_subset()) {
      FLAG = flags::marked;
      if(!S->have(fieldbit::f))
	const_cast<snapshot*>(S)->add_field(fieldbit::f);
      LoopAllBodies(S,b)
	if(in_subset(b)) b.mark(); else b.unmark();
    }
    // 1.2 build tree
    OctTree TREE(S, K+1, 0, Default::MaxDepth, FLAG);
    if(falcON::debug(1)) {
      clock_t CPU1 = clock();
      DebugInfo("density::manipulate(): %f sec needed for tree build\n",
		(CPU1 - CPU0)/real(CLOCKS_PER_SEC));
      CPU0 = CPU1;
    }
    // 2. find Kth nearest neighbours and estimate density
    if(!S->have(fieldbit::r))
      const_cast<snapshot*>(S)->add_field(fieldbit::r);
    prepare(N);
    unsigned NIAC;
    ProcessNearestNeighbours(&TREE,K,&SetDensity,NIAC,true);
    if(falcON::debug(1)) {
      clock_t CPU1 = clock();
      DebugInfo("density::manipulate(): %f sec needed for density estimation;"
		" %d neighbour updates\n",
		(CPU1 - CPU0)/real(CLOCKS_PER_SEC),NIAC);
    }
    //
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} }

__DEF__MAN(falcON::Manipulate::density)
