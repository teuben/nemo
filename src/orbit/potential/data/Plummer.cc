//-----------------------------------------------------------------------------+
//                                                                             |
// Plummer.cc                                                                  |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2004                                               |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
#define POT_DEF
#include <cmath>
#include <defacc.h>
////////////////////////////////////////////////////////////////////////////////
namespace {
  class Plummer {
    double GM, Aq;
  public:
    static const char* name() { return "Plummer"; }
    Plummer(const double*pars,
	    int          npar,
	    const char  *file)
    {
      if(npar < 3)
	warning("%s: recognizing 3 parameters:\n"
		" omega        pattern speed (ignored)           [0]\n"
		" G*M          mass;                             [1]\n"
		" a            scale radius;                     [1]\n"
		"\nthe potential is given by\n\n"
		"                 G*M\n"
                "    Phi = - ---------------\n"
		"            sqrt(a^2 + r^2)\n\n",name());
      if(file)
	warning("%s: file \"%s\" ignored",name(),file);
      double
	o   = npar>0? pars[0] : 0.,
	m   = npar>2? pars[1] : 1.,
	a   = npar>3? pars[2] : 1.;
      Aq    = a*a;
      GM    = -m;
      if(npar>3) warning("%s: skipped parameters beyond 3",name());
      nemo_dprintf (1,"initializing %s\n",name());
      nemo_dprintf (1," parameters : pattern speed = %f (ignored)\n",o);
      nemo_dprintf (1,"              mass          = %f\n",m);
      nemo_dprintf (1,"              scalelength   = %f\n",a);
    }
    template<typename scalar>
    void potacc(scalar const&Rq,
		scalar      &P,
		scalar      &T) const
    {
      T = 1/(Aq+Rq);
      P = GM * sqrt(T);
      T*= P;
    }
  };
}
//------------------------------------------------------------------------------

__DEF__ACC(SphericalPot<Plummer>)
__DEF__POT(SphericalPot<Plummer>)

//------------------------------------------------------------------------------
