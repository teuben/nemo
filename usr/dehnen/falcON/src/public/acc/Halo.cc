// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// Halo.cc                                                                     |
//                                                                             |
// Copyright (C) 2007 Walter Dehnen                                            |
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
#include <defacc.h>
#include <ctime>
#include "public/lib/halo.cc"
////////////////////////////////////////////////////////////////////////////////
namespace {
  using namespace falcON;

  class Halo : 
    private ModifiedDoublePowerLawHalo,
    private HaloPotential
  {
    typedef DoublePowerLawHalo DPLH;
    static double inner(const char*model, const double*par)
    {
      DPLH::Model M = DPLH::model(model);
      double def = DPLH::inner_default(M);
      if(def == DPLH::null_value())
	return par?  *par : 7./9.;
      else if(par && *par != def)
	falcON_WarningN("external potential \"Halo\": "
			"par[3]=g_i=%g ignored, using g_i=%g for model '%s'",
			*par, def, model);
      return def;
    }
    static double outer(const char*model, const double*par)
    {
      DPLH::Model M = DPLH::model(model);
      double def = DPLH::outer_default(M);
      if(def == DPLH::null_value())
	return par?  *par : 31./9.;
      else if(par && *par != def)
	falcON_WarningN("external potential \"Halo\": "
			"par[4]=g_o=%g ignored, using g_o=%g for model '%s'",
			*par, def, model);
      return def;
    }
    static double trans(const char*model, const double*par)
    {
      DPLH::Model M = DPLH::model(model);
      double def = DPLH::trans_default(M);
      if(def == DPLH::null_value())
	return par?  *par : 4./9.;
      else if(par && *par != def)
	falcON_WarningN("external potential \"Halo\": "
			"par[5]=eta=%g ignored, using eta=%g for model '%s'",
			*par, def, model);
      return def;
    }
  public:
    //
    static const char* name() { return "Halo"; }
    //
    Halo(const double*pars, int npar, const char*file) :
      ModifiedDoublePowerLawHalo(npar>1? pars[1] : 1.,
				 npar>7? pars[7] : 0.,
				 npar>6? pars[6] : 0.,
				 npar>2? pars[2] : 1.,
				 inner(file, npar>3? pars+3:0),
				 outer(file, npar>4? pars+4:0),
				 trans(file, npar>5? pars+5:0)),
      HaloPotential(*this,0)
    {
      if((npar<8 && nemo_debug(1)) || nemo_debug(2) ) {
	DebugInfo("external potential \"Halo\" recognizes 8 parameters:\n");
	std::cerr<<
	  "   omega   pattern speed (ignored)   [0]\n"
	  "   r_s     scale radius              [1]\n"
	  "   m_t     total mass                [1]\n"
	  "   g_i     inner power-law slope     [7/9]\n"
	  "   g_o     outer power-law slope     [31/9]\n"
	  "   eta     transition steepness      [4/9]\n"
	  "   r_t     truncation radius (0->oo) [0]\n"
	  "   r_c     core radius               [0]\n"
	  " The halo density proportional to\n"
	  "\n"
	  "      Model(x) * sech(r/r_t)\n"
	  "\n"
	  " with x=sqrt(r^2+r_c^2)/r_s,\n"
	  " and\n"
	  "                  -gi   eta    [gi-go]/eta\n"
	  "      Model(x) = x    (x    + 1).\n"
	  "\n"
	  " If file is given, it's interpreted as follows:\n"
	  "   Plummer:   gi=0, go=5, eta=2\n"
	  "   Jaffe:     gi=2, go=4, eta=1\n"
	  "   Hernquist: gi=1, go=4, eta=1\n"
	  "   Dehnen:          go=4, eta=1\n"
	  "   NFW:       gi=1, go=4, eta=1\n"
	  "   Moore:     gi=3/2, go=3, eta=3/2\n"
	  "   DM:        gi=7/9, go=31/9, eta=4/9\n"
	  " and differing values are ignored.\n";
      }
      if(nemo_debug(2)) {
	DebugInfo("external potential \"Halo\" initialized with:\n");
	std::cerr<<
	  "   r_s = "<<scale_radius()<<"\n"
	  "   m_t = "<<total_mass()  <<"\n"
	  "   g_i = "<<inner_gamma() <<"\n"
	  "   g_o = "<<outer_gamma() <<"\n"
	  "   eta = "<<transition()  <<"\n"
	  "   r_t = "<<trunc_radius()<<"\n"
	  "   r_c = "<<core_radius() <<'\n';
      }
      if(npar>8)
	falcON_Warning("external potential \"Halo\": "
		       "skipping parameters beyond 8\n");
    }
    //--------------------------------------------------------------------------
    template<typename scalar>
    void potacc(scalar const&Rq, scalar&P, scalar&T) const;
  };
  //----------------------------------------------------------------------------
  template<> inline
  void Halo::potacc<double>(double const&Rq, double&P, double&T) const
  {
    P = PotAcc(Rq,T);
  }
  //----------------------------------------------------------------------------
  template<> inline
  void Halo::potacc<float>(float const&Rq, float&P, float&T) const
  {
    double Td;
    P = PotAcc(double(Rq),Td);
    T = Td;
  }
}
//------------------------------------------------------------------------------

__DEF__ACC(SphericalPot<Halo>)
__DEF__POT(SphericalPot<Halo>)

//------------------------------------------------------------------------------
