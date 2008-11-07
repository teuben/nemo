// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   src/public/manip/bound_centre.cc
///
/// \author Walter Dehnen
/// \date   2006-2008
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
// v 0.0    11/07/2006  WD created
// v 0.1    06/11/2007  WD deBUGged
// v 0.2    29/10/2008  WD append to existing output files, output format
// v 1.0    06/11/2008  WD new algorithm for centre
////////////////////////////////////////////////////////////////////////////////
#include <utils/inline_io.h>
#include <public/defman.h>
#include <public/io.h>
#include <utils/numerics.h>

namespace {
  using namespace falcON;
  const unsigned K_default = 256;
  const double   A_default = 3;
  inline real ppp(body const&b) {
    real p = zero;
    if(has_pot(b)) p += pot(b);
    if(has_pex(b)) p += pex(b);
    return p;
  }
  inline real ppp(const snapshot*S, bodies::index i) {
    real p = zero;
    if(S->have_pot()) p += S->pot(i);
    if(S->have_pex()) p += S->pex(i);
    return p;
  }
}
namespace falcON {
namespace Manipulate {
  // ///////////////////////////////////////////////////////////////////////////
  //
  // class bound_centre
  //
  /// manipulator: sets the centre related to the most bound region.
  ///                                                                           
  /// This manipulator finds the most bound body from those in_subset()
  /// (default: all) and finds its \f$K\f$ nearest neighbours (including
  /// itself). From the \f$K/8\f$ (at least 4) most bound of these, it
  /// computes the centre weighted by \f$ (E_{\mathrm{max}}-E)^\alpha\f$.
  ///
  /// The centre position and velocity are put in 'xcen' and 'vcen' to be used
  /// by subsequent manipulators, such as sphereprof.
  ///
  /// Meaning of the parameters:\n
  /// pars[0]: number \f$K\f$ of particles to consider (default: 256)\n
  /// pars[1]: power \f$\alpha\f$ in energy weighting (default: 3)\n
  /// file: write centre position to file.
  ///
  /// Usage of pointers: sets 'xcen' and 'vcen'\n
  /// Usage of flags:    uses in_subset()\n
  ///
  // ///////////////////////////////////////////////////////////////////////////
  class bound_centre : public manipulator {
  private:
    unsigned       K;
    double         A;
    mutable output OUT;
    mutable vect   XCEN,VCEN;
    mutable bool   FIRST;
    mutable Array<bodies::index> Nb;
    mutable Array<double>        En;
    mutable Array<int>           In;
    //--------------------------------------------------------------------------
    static std::ostream&print_line(std::ostream&to) {
      return to << 
	"#-------------"
	"----------------"
	"----------------"
	"----------------"
	"----------------"
	"----------------"
	"----------------\n";
    }
    static std::ostream&print_head(std::ostream&to) {
      return to << 
	"#        time "
	"              x "
	"              y "
	"              z "
	"            v_x "
	"            v_y "
	"            v_z\n";
    }
    //--------------------------------------------------------------------------
  public:
    const char*name    () const { return "bound_centre"; }
    const char*describe() const {
      return 
	"provides 'xcen' and 'vcen' as position and velocity of the"
        "most bound body in_subset() (default: all)";
    }
    //--------------------------------------------------------------------------
    fieldset need   () const { return fieldset::p | fieldset::basic; }
    fieldset provide() const { return fieldset::o; }
    fieldset change () const { return fieldset::o; }
    //--------------------------------------------------------------------------
    bound_centre(const double*pars,
		 int          npar,
		 const char  *file) falcON_THROWING
    : K    ( npar>0? max(1,int(pars[0])) : K_default ),
      A    ( npar>1? pars[1]  : A_default ),
      OUT  ( file, true ),
      XCEN ( vect(zero) ),
      VCEN ( vect(zero) ),
      FIRST( true ),
      Nb   ( K ),
      En   ( K ),
      In   ( K )
    {
      if(debug(2) || npar > 2)
	std::cerr <<
	  " Manipulator \"bound_centre\" centre:\n"
	  " find centre of most bound region and put it in 'xcen' and 'vcen'\n"
	  " Of the K nearest neighbours of the most bound particle, we \n"
	  " compute the |Phi|^A weighted centre.\n"
	  " pars[0] K (default "<<K_default<<")\n"
	  " pars[1] A (default "<<A_default<<")\n";
      if(A < 0)
	falcON_THROW("Manipulator \"%s\": "
		     "A=%f < 0",name(),A);
      if(npar>2)
	falcON_WarningN("Manipulator \"%s\": "
			"skipping parameters beyond 2\n",name());
    }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*) const;
    //--------------------------------------------------------------------------
    ~bound_centre() {}
  };// class bound_centre
  //////////////////////////////////////////////////////////////////////////////
  bool bound_centre::manipulate(const snapshot*S) const
  {
    // 0  first call: write header to output file, if any
    if(FIRST && OUT) {
      FIRST = false;
      print_line(OUT);
      OUT  << "#\n# output from Manipulator \"" << name() << "\"\n#\n";
      RunInfo::header(OUT);
      print_head(OUT);
      print_line(OUT);
    }
    // 1  find most bound body
    real Emin = 1.e10;
    body B;
    LoopSubsetBodies(S,b) {
      real E = half*norm(vel(b)) + ppp(b);
      if(E < Emin) {
	Emin = E;
	B    = b;
      }
    }
    if(K>1) {
      // 2  find K nearest neighbours of most bound body
      unsigned n = S->findNeighbours(B,K,Nb);
      // 3   get |Phi|^A weigted centre of K/2 most bound of those
      // 3.1 sort neighbours according to their energy
      for(int i=0; i!=n; ++i)
	En[i] = half*norm(S->vel(Nb[i])) + ppp(S,Nb[i]);
      HeapIndex(En.array(),n,In.array());
      double Emax = En[In[n-1]];
      // 3.2 loop K/8 most bound and find weighted centre
      double W(0.);
      vect_d X(0.), V(0.);
      for(int i=0; i!=min(n,max(4u,K/8)); ++i) {
	double w = pow(Emax-En[In[i]],A);
	W += w;
	X += w * S->pos(Nb[In[i]]);
	V += w * S->vel(Nb[In[i]]);
      }
      X   /= W;
      V   /= W;
      XCEN = X;
      VCEN = V;
    } else {
      XCEN = pos(B);
      VCEN = vel(B);
    }
    // 4  output
    if(OUT)
      OUT << ' '
	  << print(S->time(),12,8) << ' '
	  << print(XCEN,15,8) << ' '
	  << print(VCEN,15,8) << std::endl;
    // 5  put centre under 'xcen' and 'vcen'
    S->set_pointer(&XCEN,"xcen");
    S->set_pointer(&VCEN,"vcen");
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} }

__DEF__MAN(falcON::Manipulate::bound_centre)

////////////////////////////////////////////////////////////////////////////////
