//-----------------------------------------------------------------------------+
//                                                                             |
// truncNFW.cc                                                                 |
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
// defines as NEMO potential and acceleration field                            |
//                                                                             |
//          M0 sech(r/r_t)                                                     |
// rho(r) = ---------                                                          |
//          r^inner (r+r_s)^outer                                              |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// Versions                                                                    |
// 1.0   21-nov-2005    created from NFW.cc                               PJM  |
//                                                                             |
//-----------------------------------------------------------------------------+
#define POT_DEF

#include <iostream>
#include <cmath>
#include <defacc.h>
#include <utils/spline.h>
#include <utils/numerics.h>
#include <utils/Pi.h>

using std::exp;
using std::log;
using std::pow;
using std::sqrt;
////////////////////////////////////////////////////////////////////////////////
namespace {
using namespace WDutils;
// Halo Model -------------------------------------------------------------
//
class Halo {
const double b,g,rs,irs,rt;
  double rsont;
public:
  Halo(const double _g,
       const double _b,
       const double _rs,
       const double _rt):
    b   (_b),                 // Outer density exponent
    g   (_g),                 // Inner density exponent
    rs  (_rs),                // Scale radius
    irs (1./rs),
    rt  (_rt)                 // Truncation radius (if required)
  {
    rsont=(rt)? rs/rt : 0;
  }

double operator() (const double r) const
    {
      register double x=r*irs,
	x1=1.+ x;
      if(rt){return pow(x,-g)*pow(x1,(g-b))/(exp(x*rsont)+exp(-x*rsont));}
      else  return pow(x,-g)*pow(x1,(g-b)); // put density here
    }

  ~Halo();

};
//--------------------------------------------------------------------------------

  Halo *RHO;
  spline<double> *SPLINEM,*SPLINEP;
  inline double dM(double const&lr, double const&M)// dM/dln r                  
{
  register double r = exp(lr);                     // ln r as independent var   
  return r*r*r*(*RHO)(r);
}
  inline double dP(double const&l, double const&F)
  {
    const double x = exp(l);
    return (*SPLINEM)(l) / x;                      // dP       = M/r * dln r    
  }

  class truncNFWPot {
  public:
    int n;
    double rmin,rmax,M;
    double *lr,*mh,*phi,*RMIN,*RMAX;

    static const char* name() { return "truncNFW"; }
    //--------------------------------------------------------------------------
    truncNFWPot(const double*pars,
	   int          npar,
	   const char  *file)
    {
      if(npar < 6)
	warning("%s: recognizing 6 parameters:\n"
		" omega        pattern speed (ignored)           [0 ]\n"
		" M            Halo mass;                        [1 ]\n"
		" r_s          scale radius;                     [1 ]\n"
		" r_t          truncation radius;                [10]\n"
		" inner        inner density exponent;           [1]\n"
		" outer        outer density exponent;           [3]\n"
		"the density is proportional to\n\n"
		"            M sech(r/r_t)     \n"
                "        --------------------  \n" 
		"       r^inner (r+r_s)^outer\n\n",name());
      if(file)
	warning("%s: file \"%s\" ignored",name(),file);
      double
	o     = npar>0? pars[0] : 0.,
	M     = npar>1? pars[1] : 1.,
	r_s   = npar>2? pars[2] : 1.,
	r_t   = npar>3? pars[3] : 10.,
	inner = npar>4? pars[4] : 1.,
	outer = npar>5? pars[5] : 3.;
      RHO = new Halo(inner,outer,r_s,r_t);
      const double rmin=0.00000001*r_s, rmax=100000000*r_s;
      int n1 = int(200*log10(rmax/rmin));
      int  n  = 1+n1;
      double Mtmp;
      const double dlr= log(rmax/rmin)/double(n1);
      lr   = new double[n];
      mh   = new double[n];
      phi  = new double[n];
      RMIN = new double;
      RMAX = new double;
      *RMIN = rmin;
      *RMAX = rmax;
      Mtmp=rmin*rmin*rmin*pow(rmin/r_s,-inner)/(3-inner)
	+ rmin*rmin*rmin*(inner-outer)*pow(rmin/r_s,1.-inner)/(4-inner);
      mh[0]=Mtmp;
      lr[0]=log(rmin);
      for(int i=1; i!=n; ++i) {
	lr[i] = lr[0] + i * dlr;
	Mtmp  = rk4(Mtmp,lr[i-1],dlr,dM); //Integrate mass outwards
	mh[i]  = FPi * Mtmp;
      }

      Mtmp=M/mh[n1];
      for(int i=0; i!=n; ++i) {mh[i] *= Mtmp;}
      double
	s0h     = FPi * Mtmp * pow(rmin ,3) * (*RHO)(rmin),
	sNh     = FPi * Mtmp * pow(rmax ,3) * (*RHO)(rmax);
      SPLINEM = new spline<double>(n,lr,mh,&s0h,&sNh);   // static spline: m_halo(ln r)
      double phitmp = -mh[n1]/exp(lr[n1]);
      phi[n1]=phitmp;
      for(int i=n-2; i!=-1; --i) {
	phitmp     = rk4(phitmp,lr[i+1],-dlr,dP);         // integrate phi inwards
      	phi[i] = phitmp;
      }
      double
	s0p =-mh[0] /exp(lr[0]) ,
	sNp =-mh[n1]/exp(lr[n1]);
      SPLINEP = new spline<double>(n,lr,phi,&s0p,&sNp);   // static spline: phi_halo(ln r)
      if(npar>6) warning("%s: skipped parameters beyond 6",name());
      nemo_dprintf (1,"initializing %s\n",name());
      nemo_dprintf (1," parameters : pattern speed  = %f (ignored)\n",o);
      nemo_dprintf (1,"              Halo mass      = %f\n",M);
      nemo_dprintf (1,"              scale radius   = %f\n",r_s);
      nemo_dprintf (1,"              trunc radius   = %f\n",r_t);
      nemo_dprintf (1,"              inner dens exp = %f\n",inner);
      nemo_dprintf (1,"              outer dens exp = %f\n",outer);
    }
    //--------------------------------------------------------------------------

    template<typename scalar>
    void potacc(scalar const&Rq,
		scalar      &P,
		scalar      &fR) const
    {
      register scalar
	R = sqrt(Rq),
	iR= 1/R,
	lR= log(R);
      fR  = (R>(*RMAX))? (*SPLINEM)(double(log(*RMAX))) :
	    (R<(*RMIN))? (*SPLINEM)(double(log(*RMIN))) 
	    * pow(R,3) * pow((*RMIN),-3)                :
	    (*SPLINEM)(double(lR));
      fR *= -iR*iR*iR;
      P   = (R>(*RMAX))? (*SPLINEM)(double(log(*RMAX)))*iR :
	    (R<(*RMIN))? (*SPLINEP)(double(log(*RMIN)))    :
	    (*SPLINEP)(double(lR));
    }

    ~truncNFWPot() {
      if (RHO)  delete   RHO;
      if (lr)   delete[] lr;
      if (mh)   delete[] mh;
      if (phi)  delete[] phi;
      if (SPLINEM) delete SPLINEM;
      if (SPLINEP) delete SPLINEP;
    }

  }; // class NFW {
} // namespace {
//------------------------------------------------------------------------------

__DEF__ACC(SphericalPot<truncNFWPot>)
__DEF__POT(SphericalPot<truncNFWPot>)

//------------------------------------------------------------------------------
