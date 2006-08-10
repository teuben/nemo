//-----------------------------------------------------------------------------+
//                                                                             |
// DehnenMcLaughlin.cc                                                         |
//                                                                             |
// Copyright (C) 2004-2006 Walter Dehnen                                       |
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
//                                                                             |
// Versions                                                                    |
// 0.1    09/08/2006  WD use $NEMOINC/defacc.h                                 |
//-----------------------------------------------------------------------------+
#undef POT_DEF
#include <iostream>
#include <utils/WDMath.h>
#include <defacc.h> // $NEMOINC/defacc.h
////////////////////////////////////////////////////////////////////////////////
namespace {
  using namespace WDutils;
  //////////////////////////////////////////////////////////////////////////////
  class DehnenMcLaughlin {
    const double   GM,a,b0,e;            // defininng parameters
    const double   ia,nM,eta,kk,Pf;      // derived parameters
    const BetaFunc BPs;                  // Beta_u( 1/eta, (1-beta0)/eta + 1/2 )
  public:
    static const char* name() { return "DehnenMcLaughlin"; }
    DehnenMcLaughlin(const double*pars,
		     int          npar,
		     const char  *file)
      : GM  ( npar>1? -std::abs(pars[1]) : -1. ),
	a   ( npar>2? pars[2] : 1. ),
	b0  ( npar>3? pars[3] : 0. ),
	e   ( npar>4? pars[4] : 3. ),
	ia  ( 1./a ),
	nM  ( (e+2)/(e-2) ),
	eta ( 2*(e-2)*(2-b0)/(6+e) ),
	kk  ( eta*nM-3 ),
	Pf  ( GM*ia/eta ),
	BPs ( 1./eta, (1-b0)/eta+0.5 )
    {
      if(npar>5 || debug(2))
	std::cerr<<
	  "falcON Debug Info: acceleration \"DehnenMcLaughlin\""
	  "recognizing 5 parameters:\n"
	  "   omega      pattern speed (ignored)            [0]\n"
	  "   G*M        mass;                              [1]\n"
	  "   a          scale radius;                      [1]\n"
	  "   b0         inner halo anisotropy              [0]\n"
	  "   e          assume rho/sigma_r^e is power law  [3]\n"
	  "  the potential is given by the density\n\n"
	  "              4+eta-2b0  M   -g0     eta -2e/(2-e)\n"
	  "      rho = - --------- --- x    (1+x   )\n"
	  "                  8 Pi   a\n\n"
	  "  with x=r/a and\n\n"
	  "    eta= 2*(e-2)*(2-b0)/(6+e)  [4/9]\n"
	  "    g0 = 1-eta/2+b0            [7/9]\n\n";
      if(file)
	warning("acceleration \"%s\": file \"%s\" ignored",name(),file);
      if(a<=0.)
	error("acceleration \"%s\": a=%f <= 0\n",name(),a);
      if(e<=2.)
	error("acceleration \"%s\": e=%f <= 2\n",name(),e);
      if(b0>1.)
	error("acceleration \"%s\": b0=%f > 1\n",name(),b0);
      if(npar>5) warning("acceleration \"%s\":"
			 " skipped parameters beyond 5",name());
    }
    ///                                                                         
    /// routine used by class SphericalPot<DehnenMcLaughlin> of defacc.h        
    ///                                                                         
    /// Note that the potential for a Dehnen & McLaughlin (2005, hereafter D&M) 
    /// model is given by (in units of GM/a)                                    
    ///                                                                         
    ///    Psi = Psi_0 - Beta_y ( (1-beta0)/eta + 1/2, 1/eta ) / eta            
    ///                                                                         
    /// with                                                                    
    ///                                                                         
    ///    Psi_0 = Beta( (1-beta0)/eta + 1/2, 1/eta ) / eta                     
    ///                                                                         
    /// and                                                                     
    ///                                                                         
    ///    y = x^eta / (1 + x^eta).                                             
    ///                                                                         
    /// Here, Beta_u(p,q) the incomplete beta function, see eq (40h) of D&M.    
    /// At small radii (note that there is a typo in eq 40j of D&M)             
    ///                                                                         
    /// Psi = Psi_0 - 2/(2+eta-2*beta0) x^(eta/2+1-beta0)                       
    ///                                                                         
    template<typename scalar>
    void potacc(scalar const&rq,
		scalar      &P,
		scalar      &T) const
    {
      if(rq == 0.) {
	P = Pf * BPs(1.);
	T = 0;
      } else {
	double r  = std::sqrt(rq);
	double u  = std::pow(ia*r,eta);
	double u1 = 1./(1+u);
	P = Pf * BPs(u1);
	T = GM * std::pow(u*u1,nM) / (rq*r);
      }
    }
  };
}
//------------------------------------------------------------------------------

__DEF__ACC(SphericalPot<DehnenMcLaughlin>)

//------------------------------------------------------------------------------
