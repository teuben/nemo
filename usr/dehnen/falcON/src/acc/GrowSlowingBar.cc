//-----------------------------------------------------------------------------+
//                                                                             |
// GrowSlowingBar.cc                                                           |
//                                                                             |
// Copyright (C) 2004 Walter Dehnen                                            |
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
// defines a rotating quadrupole potential as used by Katz & Weinberg and      |
// Sellwood (note, however an error in the latter paper)                       |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// Versions                                                                    |
// 0.1    22/07/2004  WD     adapted from GrowRigidBar.cc                      |
//-----------------------------------------------------------------------------+
#include <cmath>
#include <fstream>
#include <timer.h>
#include <center.h>
#include <defacc.h>
//=============================================================================#
// define C++ implementation of potential                                      |
//=============================================================================#
namespace {
  //////////////////////////////////////////////////////////////////////////////
  class slowing_bar {
    mutable double time_om;                        // time of omega             
    mutable double omega;                          // pattern speed             
    double Ib;                                     // moment of intertia        
    double GM_asq_alpha;                           // GM * alpha * a^2          
    double gamma,five_over_gamma,gamma_half;       // gamma, 5/gamma, gamma/2   
    bool   gamma_is_one;                           // is gamma == 1             
    bool   gamma_is_two;                           // is gamma == 2             
    double betaa_to_gamma;                         // (beta * a)^gamma          
    mutable double cos_omt, sin_omt;               // current cos/sin Omega*time
    //--------------------------------------------------------------------------
    template<typename scalar>
    static scalar twice (scalar const&x) { return x+x; }
    //--------------------------------------------------------------------------
  protected:
    slowing_bar() {}
    void init(double const&o,
	      double const&gm,
	      double const&alpha,
	      double const&beta,
	      double const&a,
	      double const&g)
    {
      time_om         = -1.e15;
      omega           = o;
      Ib              = ???????????????;
      GM_asq_alpha    = gm*alpha*a*a;
      gamma           = g;
      five_over_gamma = 5./g;
      gamma_half      = 0.5*g;
      gamma_is_one    = g == 1.;
      gamma_is_two    = g == 2.;
      betaa_to_gamma  = pow(beta*a,gamma);
      cos_omt         = 1.;
      sin_omt         = 0.; 
    }
    //--------------------------------------------------------------------------
    void set_time(double t, double deltaL) const
    {
      if(time_om == -1.e15) { // first call ever: don't change omega
	time_om = t;
	cos_omt = cos(omega * t);
	sin_omt = sin(omega * t);
      } else {                // subsequent calls:
	double dOm  = -deltaL / Ib;
	double dphi = (omega + 0.5*dOm) * (t-time_om);
	omega   += dOm;
	time_om  = t;
	double c = cos(dphi), s=sin(dphi), cot=cos_omt, sot=sin_omt;
	cos_omt  = cot * c - sot * s;
	sin_omt  = sot * c + cot * s;
      }
    }
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar>
    void force(const scalar*x,
	       scalar      *f,
	       scalar      &p) const
    {
      scalar 
	X    = x[1] * sin_omt + x[0] * cos_omt,
	Y    = x[1] * cos_omt - x[0] * sin_omt,
	tmp1 = X * X,
	tmp2 = Y * Y,
	rq   = NDIM==2? tmp1+tmp2 : tmp1+tmp2+x[2]*x[2];
      scalar rg = gamma_is_two? rq : 
	gamma_is_one? sqrt(rq) :
	pow(double(rq),gamma_half);
      tmp2   = tmp1 - tmp2;
      tmp1   = 1 / (rg + betaa_to_gamma);
      scalar
	pot  = gamma_is_two? GM_asq_alpha*sqrt(tmp1)*tmp1*tmp1 :
	gamma_is_one? GM_asq_alpha*pow(tmp1,5)          :
	GM_asq_alpha*pow(double(tmp1),five_over_gamma),
	dpdx = twice((X * cos_omt + Y * sin_omt) * pot),
	dpdy = twice((X * sin_omt - Y * cos_omt) * pot);
      pot   *= tmp2;
      p      =-pot;
      if(rq==0. && gamma < 2.0)
	v_set<NDIM>(f,scalar(0));
      else {
	tmp1 *= gamma_is_two? -5 * pot : -5 * pot * rg/rq;
	v_asstimes<NDIM>(f,x,tmp1);
	f[0] += dpdx;
	f[1] += dpdy;
      }
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  class DeltaAngMom {
    mutable bool   first;
    mutable double L;
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar>
    void SetL(int          N,
	      const scalar*M,
	      const scalar*X,
	      const scalar*V)const {
      L = 0;
      for(int n=0,nn=0; n!=N; ++n,nn+=NDIM)
	L += M[n] * (X[nn]*V[nn+1] - X[nn+1]*V[nn]);
    }
    //--------------------------------------------------------------------------
  public:
    DeltaAngMom() : first(true) {}
    //--------------------------------------------------------------------------
    double const&AngMom() const { return L; }
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar>
    double deltaL(int          N,
		  const scalar*M,
		  const scalar*X,
		  const scalar*V) const {
      if(first) {
	SetL<NDIM>(N,M,X,V);
	first = false;
	return 0.;
      } else {
	double Lold = L;
	SetL<NDIM>(N,M,X,V);
	return L-Lold;
      }
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  class GrowSlowingBar :
    private slowing_bar,
    private timer,
    private centering,
    private DeltaAngMom
  {
    mutable double         Ampl;
    bool                   centring;
    mutable std::ofstream  Cout;
  public:
    //--------------------------------------------------------------------------
    static const char* name() {
      return "GrowSlowingBar";
    }
    //--------------------------------------------------------------------------
    GrowSlowingBar(const double*pars,
		 int          npar,
		 const char  *file) {
      if(npar < 10) 
	warning("GrowSlowingBar potential recognizes 10 parameters:\n"
		" par[0] = pattern speed                                [1]\n"
		" par[1] = G*M                                          [1]\n"
		" par[2] = size of bar                                  [1]\n"
		" par[3] = parameter alpha                              [0.1404]\n"
		" par[4] = dimensionless scale radius beta              [0.4372]\n"
		" par[5] = gamma: controlling shape of quadrupole       [2]\n"
		" par[6] = Ncen: if>0, the origin of the bar quadrupole\n"
		"                is taken to be the centre of the snap-\n"
		"                shot obtained as density maximum with\n"
		"                Ncen bodies (see density_centre)       [0]\n"
		" par[7] = controlling growth factor, see below         [0]\n"
		" par[8] = t0: start time for growth                    [0]\n"
		" par[9] = tau: time scale for growth                   [1]\n"
		" with par[7]=0: %s\n"
		"      par[7]=1: %s\n"
		"      par[7]=2: %s\n"
		"      par[7]=3: %s\n"
		" if accfile is given and centering is done, the centre\n"
		" position will be written to accfile at every time,\n"
		" appending to an existing file, if applicable\n",
		timer::describe(timer::adiabatic),
		timer::describe(timer::saturate),
		timer::describe(timer::quasi_linear),
		timer::describe(timer::linear));
	    
      double
	omega = npar>0? pars[0] : 1.,
	GM    = npar>1? pars[1] : 1.,
	a     = npar>2? pars[2] : 1.,
	alpha = npar>3? pars[3] : 0.1404,
	beta  = npar>4? pars[4] : 0.4372,
	gamma = npar>5? pars[5] : 2.;
      int
	Nmin  = npar>6? int(pars[6]) : 0;
      timer::index
	timin = (timer::index)(npar>7? int(pars[7]) : 0);
      double
	t0    = npar>8? pars[8] : 0.,
	tau   = npar>9? pars[9] : 1.;

      slowing_bar::init(omega, GM, alpha, beta, a, gamma);
      timer::init(timin,t0,tau);
      centering::init(Nmin);
      Ampl = timer::operator() (0.);
      centring = Nmin > 0;
      bool appending = true;
      if(file) {
	Cout.open(file,std::ios::out | std::ios::app);
	if(!Cout.is_open()) {
	  appending = false;
	  Cout.open(file,std::ios::out);
	}
      }
      if(Cout)
	if(appending)
	  Cout  <<"# appending: acceleration GrowSlowingBar"<<std::endl;
	else {
	  Cout  <<"#\n"
		<<"# file \""<<file<<"\"\n"
		<<"# generated by acceleration GrowSlowingBar\n"
		<<"#\n"
		<<"# time, L_halo ";
	  if(centring)
	    Cout<<" centre position, centering algorithm converged?";
	  Cout  <<"\n#"<<std::endl;
	}
      if(npar>10) warning("Skipped potential parameters beyond 10");
    }
    //--------------------------------------------------------------------------
    bool NeedMass() const { return centring; }
    bool NeedVels() const { return true; }
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar>
    void set_time(double       t,
		  int          N,
		  const scalar*M,
		  const scalar*X,
		  const scalar*V) const
    {
      slowing_bar::set_time(t, deltaL<NDIM>(N,M,X,V) );
      Ampl = timer::operator() (t);
      bool okay = true;
      if(Ampl!=0. && centring) {
	okay = centering:: template update<NDIM,scalar>(M,X,N);
	if(!okay)
	  warning("problems centering snapshot for acceleration at time %f",t);
      }
      if(Cout) {
	Cout  <<t<<"  "<<AngMom();
	if(centring)
	  if(Ampl!=0)
	    Cout<<"  "
		<<centre()[0]<<' '
		<<centre()[1]<<' '
		<<centre()[2]<<"  "
		<< (okay? "yes" : "NO");
	  else
	    Cout<<"  0 0 0 ---";
	Cout  <<std::endl;
      }
    }
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar>
    void acc(const scalar*,
	     const scalar*X,
	     const scalar*,
	     scalar      &P,
	     scalar      *A) const
    {
      if(Ampl) {
	if(centring) {
	  scalar x[NDIM];
	  v_assdif<NDIM>(x,X,centre());
	  slowing_bar::force<NDIM>(x,A,P);
	} else
	  slowing_bar::force<NDIM>(X,A,P);
	if(Ampl != 1.) {
	  P *= Ampl;
	  v_mul<NDIM>(A,Ampl);
	}
      } else {
	P = scalar(0);
	v_set<NDIM>(A,scalar(0));
      }
    }
  };
}                                                  // END: unnamed namespace    
//------------------------------------------------------------------------------

__DEF__ACC(GrowSlowingBar)

//------------------------------------------------------------------------------
