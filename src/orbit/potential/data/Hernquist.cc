//-----------------------------------------------------------------------------+
//                                                                             |
// Hernquist.cc                                                                |
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
//=============================================================================#
// define C++ implementation of potential                                      |
//=============================================================================#
namespace {
  class Hernquist {
    double Rs, Eq, GM;
  public:
    static const char* name() { return "Hernquist"; }
    Hernquist(const double*pars,
	      int          npar,
	      const char  *file)
    {
      if(npar < 4)
	warning("%s: recognizing 4 parameters:\n"
		"     omega -- pattern speed (ignored)\n"
		"     G*M   -- mass; defaults to 1\n"
		"     a     -- scale radius; defaults to 1\n"
		"     e     -- softening or core radius; default to 0\n"
		" the potential is given by\n\n"
		"                    G M\n"
		"    Phi = - ----------------- .\n"
		"            sqrt(r^2+e^2) + a\n\n",name());
      if(file)
	warning("%s: file \"%s\" ignored",name(),file);
      double o,e,m;
      o  = npar>0? pars[0] : 0.;
      m  = npar>1? pars[1] : 1.;
      Rs = npar>2? pars[2] : 1.;
      e  = npar>3? pars[3] : 0.;
      GM =-m;
      Eq = e*e;
      if(npar>4) warning("%s: skipped parameters beyond 4",name());
      nemo_dprintf (1,"initializing %s\n",name());
      nemo_dprintf (1," parameters : pattern speed = %f (ignored)\n",o);
      nemo_dprintf (1,"              mass          = %f\n",m);
      nemo_dprintf (1,"              scalelength   = %f\n",Rs);
      nemo_dprintf (1,"              core-radius   = %f\n",e);
    }    
    template<typename scalar>
    void potacc(scalar const&Rq,
		scalar      &P,
		scalar      &T) const
    {
      register scalar R = std::sqrt(Eq+Rq);
      T  = 1/(R+Rs);
      P  = GM * T;
      T *= P/R;
    }
  };
}
//------------------------------------------------------------------------------

__DEF__ACC(SphericalPot<Hernquist>)
__DEF__POT(SphericalPot<Hernquist>)

//------------------------------------------------------------------------------

  
