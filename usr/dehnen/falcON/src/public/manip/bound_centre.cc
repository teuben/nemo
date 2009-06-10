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
// v 1.1    21/11/2008  WD subtract sinks' contribution to gravity
////////////////////////////////////////////////////////////////////////////////
#include <public/defman.h>
#include <public/default.h>
#include <public/kernel.h>
#include <utils/numerics.h>

namespace {
  using namespace falcON;
  const unsigned K_default = 256;
  const double   A_default = 3;
  inline real&poq(body&b) {
    return has_pot(b)? b.pot() : b.pex();
  }
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
  /// If sink particles are present and are not in_subset() and have
  /// individual softening lengths, we subtract their contribution to the
  /// potential off before the above process.
  ///
  /// The centre position and velocity are put in 'xcen' and 'vcen' to be used
  /// by subsequent manipulators, such as sphereprof.
  ///
  /// Meaning of the parameters:\n
  /// pars[0]: number \f$K\f$ of particles to consider (default: 256)\n
  /// pars[1]: power \f$\alpha\f$ in energy weighting (default: 3)\n
  /// pars[2]: (0,1,2,3): kernel for subtracting off satellite
  /// file: write centre position to file.
  ///
  /// Usage of pointers: sets 'xcen' and 'vcen'\n
  /// Usage of flags:    uses in_subset()\n
  ///
  // ///////////////////////////////////////////////////////////////////////////
  class bound_centre : public manipulator {
  private:
    unsigned        K;
    double          A;
    const kern_type KERN;
    mutable output  OUT;
    mutable vect    XCEN,VCEN;
    mutable bool    FIRST;
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
    static kern_type kernel(int par) {
      switch(par) {
      case 0: return falcON::p0;
      case 1: return falcON::p1;
      case 2: return falcON::p2;
      case 3: return falcON::p3;
      default: 
	falcON_THROW("Manipulator \"bound_centre\": "
		     "unknown kern=%d\n",par);
	return Default::kernel;
      }
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
      KERN ( npar>2? kernel(int(pars[2])) : Default::kernel ),
      OUT  ( file, true ),
      XCEN ( vect(zero) ),
      VCEN ( vect(zero) ),
      FIRST( true ),
      Nb   ( K ),
      En   ( K ),
      In   ( K )
    {
      if(debug(2) || npar > 3)
	std::cerr <<
	  " Manipulator \"bound_centre\" centre:\n"
	  " find centre of most bound region and put it in 'xcen' and 'vcen'\n"
	  " Of the K nearest neighbours of the most bound particle, we \n"
	  " compute the |Phi|^A weighted centre.\n"
	  " If sink particles not in_subset() are present and have eps_i,\n"
	  " we subtract their contribution to the energy, assuming a kernel\n"
	  " as indicated by par[2]\n"
	  " pars[0] K (default "<<K_default<<")\n"
	  " pars[1] A (default "<<A_default<<")\n"
	  " pars[2] k (default 1)\n";
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
    // 1   subtract sink particles' contribution to potential
    // 1.1 are there sini particles not in_subset() with eps_i ?
    bool correct_sink = true;
    if(!S->have_eps()) correct_sink = false;
    if(correct_sink) {
      unsigned Nss = 0;  // sinks not in subset
      LoopSinkBodies(S,b)
	if(!in_subset(b)) ++Nss;
      if(Nss == 0) correct_sink = false;
    }
    // 1.2 correct for potential of sink particles
    Array<real,1> Ps(S->N_subset());
    unsigned i=0;
    LoopSubsetBodies(S,b) {
      Ps[i]  = zero;
      LoopSinkBodies(S,s) if(!in_subset(b)) {
	Ps[i] -= mass(s) * GravKernBase::Psi(KERN,dist_sq(pos(b),pos(s)),
					     square(half*(eps(b)+eps(s))));
      }
      poq(b) -= Ps[i++];
    }
    // 2   find most bound body
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
      // 3   find K nearest neighbours of most bound body
      unsigned n = S->findNeighbours(B,K,Nb);
      // 4   get |Phi|^A weigted centre of K/2 most bound of those
      // 4.1 sort neighbours according to their energy
      for(i=0; i!=n; ++i)
	En[i] = half*norm(S->vel(Nb[i])) + ppp(S,Nb[i]);
      HeapIndex(En.array(),n,In.array());
      double Emax = En[In[n-1]];
      // 4.2 loop K/8 most bound and find weighted centre
      double W(0.);
      vect_d X(0.), V(0.);
      for(i=0; i!=min(n,max(4u,K/8)); ++i) {
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
    // 5   add sink particles' contribution to potential back
    i=0;
    LoopSubsetBodies(S,b)
      poq(b) += Ps[i++];
    // 6   output
    if(OUT)
      OUT << ' '
	  << print(S->time(),12,8) << ' '
	  << print(XCEN,15,8) << ' '
	  << print(VCEN,15,8) << std::endl;
    // 7   put centre under 'xcen' and 'vcen'
    S->set_pointer(&XCEN,"xcen");
    S->set_pointer(&VCEN,"vcen");
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} }

__DEF__MAN(falcON::Manipulate::bound_centre)

////////////////////////////////////////////////////////////////////////////////
