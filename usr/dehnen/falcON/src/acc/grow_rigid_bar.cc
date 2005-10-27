//-----------------------------------------------------------------------------+
//                                                                             |
// grow_rigid_bar.cc                                                           |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// defines a rotating quadrupole potential as used by Katz & Weinberg and      |
// Sellwood.                                                                   |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// Versions                                                                    |
// 0.0    20/01/2004  WD     created from rigid_bar.c                          |
// 1.0    16/06/2004  WD/PJM corrected for a bug in Sellwood's paper           |
// 1.1    16/06/2004  WD/PJM added parameter gamma                             |
//-----------------------------------------------------------------------------+
#define POT_DEF
#include <cmath>
#include <timer.h>
#include <defacc.h>
//=============================================================================#
// define C++ implementation of potential                                      |
//=============================================================================#
namespace {
  using namespace growth;
  //////////////////////////////////////////////////////////////////////////////
  template<typename scalar>
  scalar twice(scalar const&x) { return x+x; }
  //////////////////////////////////////////////////////////////////////////////
  class rigid_bar {
    double omega;                                  // pattern speed             
    double GM_asq_alpha;                           // GM * alpha * a^2          
    double gamma,five_over_gamma,gamma_half;       // gamma, 5/gamma, gamma/2   
    bool   gamma_is_one;                           // is gamma == 1             
    bool   gamma_is_two;                           // is gamma == 2             
    double betaa_to_gamma;                         // (beta * a)^gamma          
    double cos_omt;                                // current cos Omega*time    
    double sin_omt;                                // current sin Omega*time    
    //--------------------------------------------------------------------------
  protected:
    rigid_bar() {}
    void init(double const&o,
	      double const&gm,
	      double const&alpha,
	      double const&beta,
	      double const&a,
	      double const&g) :
    {
      omega           = o;
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
    void set_time(double t)
    {
      cos_omt = cos(omega * time);
      sin_omt = sin(omega * time);
    }
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar>
    void force(const scalar*x,
	       scalar      *f,
	       scalar      &p)
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
  class GrowRigidBar :
    private rigid_bar,
    private timer {
    double Ampl
  public:
    static const char* name()
    {
      return "GrowRigidBar";
    }
    GrowRigidBar(const double*pars,
		 int          npar,
		 const char  *)
    {
      if(*npar < 6) 
	warning("GrowRigidBar potential recognizes 8 parameters:\n"
		" par[0] = pattern speed                             [1]\n"
		" par[1] = G*M                                       [1]\n"
		" par[2] = size of bar                               [1]\n"
		" par[3] = parameter alpha                           [0.1404]\n"
		" par[4] = dimensionless scale radius beta           [0.4372]\n"
		" par[5] = controlling growth factor, see below      [0]\n"
		" par[6] = t0: start time for growth                 [0]\n"
		" par[7] = tau: time scale for growth                [1]\n"
		" par[8] = gamma: controlling shape of quadrupole    [2]\n"
		" with par[5]=0: %s\n"
		"      par[5]=1: %s\n"
		"      par[5]=2: %s\n"
		"      par[5]=3: %s\n",
		npar,
		description[0],
		description[1],
		description[2],
		description[3]);
	    
      double
	omega = npar>0? pars[0] : 1.,
	GM    = npar>1? pars[1] : 1.,
	a     = npar>2? pars[2] : 1.,
	alpha = npar>3? pars[3] : 0.1404,
	beta  = npar>4? pars[4] : 0.4372;
      index
	timin = (index)(npar>5? int(pars[5]) : 0);
      double
	t0    = npar>6? pars[6] : 0.,
	tau   = npar>7? pars[7] : 1.,
	gamma = npar>8? pars[8] : 2.;

      rigid_bar::init(omega, GM, alpha, beta, a, gamma);
      timer::init(timin,t0,tau);
      Ampl = timer::operator() (0.);

      if(npar>9) warning("Skipped potential parameters beyond 9");
    }
    bool NeedMass() const { return false; }
    bool NeedVels() const { return false; }
    void set_time(double t) const
    {
      rigid_bar::set_time(t);
      Ampl = timer::operator() (t);
    }
    template<int NDIM, typename scalar>
    void acc(const scalar*,
	     const scalar*X,
	     const scalar*,
	     scalar      &P,
	     scalar      *A) const
    {
      if(Ampl) {
	rigid_bar::force<NDIM>(X,A,P);
	if(Ampl != 1.) {
	  P *= Ampl;
	  v_mul(A,Ampl);
	}
      } else {
	P = scalar(0);
	v_set(A,scalar(0));
      }
    }
  };
}                                                  // END: unnamed namespace    
//------------------------------------------------------------------------------

__DEF__ACC(GrowRigidBar)
__DEF__POT

//------------------------------------------------------------------------------
