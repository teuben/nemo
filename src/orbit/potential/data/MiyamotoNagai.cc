//-----------------------------------------------------------------------------+
//                                                                             |
// MiyamotoNagai.cc                                                            |
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

  struct MiyamotoNagai {
    double A,Bq,ABq,GM;
    //--------------------------------------------------------------------------
  public:
    static const char* name() { return "MiyamotoNagai"; }
    bool NeedMass() const { return false; }
    bool NeedVels() const { return false; }
    template<int NDIM, typename scalar>
    void set_time(double,int,const scalar*,const scalar*,const scalar*) const {}
    //--------------------------------------------------------------------------
    MiyamotoNagai(const double*pars,
		  int          npar,
		  const char  *file)
    {
      if(npar < 4)
	warning("%s: recognizing 4 parameters:\n"
		"     omega -- pattern speed (ignored)\n"
		"     G*M   -- mass; defaults to 1\n"
		"     a     -- scale radius; defaults to 1\n"
		"     b     -- scale height; defaults to 0.1\n"
		" the potential is given by\n\n"
		"                          G*M\n"
		"    Phi = - --------------------------------- .\n"
		"            sqrt(x^2+y^2+(a+sqrt(z^2+b^2))^2)\n\n",name());
      if(file)
	warning("%s: file \"%s\" ignored",name(),file);
      double
	o = npar>0? pars[0] : 0.,
	m = npar>1? pars[1] : 1.;
      A   = npar>2? pars[2] : 1.;
      double
	b = npar>3? pars[3] : 0.1;
      GM  =-m;
      Bq  = b*b;
      ABq = (A+b)*(A+b);
      if(npar > 4) warning("%s: skipped parameters beyond 4",name());
      nemo_dprintf (1,"initializing %s\n",name());
      nemo_dprintf (1," parameters : pattern speed = %f (ignored)\n",o);
      nemo_dprintf (1,"              mass          = %f\n",m);
      nemo_dprintf (1,"              scale radius  = %f\n",A);
      nemo_dprintf (1,"              scale height  = %f\n",b);
    }
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar>
    inline void acc(const scalar*,
		    const scalar*,
		    const scalar*,
		    scalar      &,
		    scalar      *) const;
  }; // class MiyamotoNagai {
  //////////////////////////////////////////////////////////////////////////////
  template<int NDIM> struct __helper {
    template<typename scalar> static
    void __acc(MiyamotoNagai const&POT,
	       const scalar       *pos,
	       scalar             &pot,
	       scalar             *acc)
    {
      register scalar
	ZB   = std::sqrt(pos[2]*pos[2]+POT.Bq),
	AZ   = POT.A+ZB,
	F    = 1/(pos[0]*pos[0]+pos[1]*pos[1]+AZ*AZ);
      pot    = POT.GM * std::sqrt(F);
      F     *= pot;
      acc[0] = pos[0] * F;
      acc[1] = pos[1] * F;
      acc[2] = ZB? pos[2]*F*AZ/ZB : 0;
    }
  }; // struct __helper
  //////////////////////////////////////////////////////////////////////////////
  template<> struct __helper<2> {
    template<typename scalar> static
    void __acc(MiyamotoNagai const&POT,
	       const scalar       *pos,
	       scalar             &pot,
	       scalar             *acc)
    {
      register scalar
	F    = 1/(pos[0]*pos[0]+pos[1]*pos[1]+POT.ABq);
      pot    = POT.GM * sqrt(F);
      F     *= pot;
      acc[0] = pos[0] * F;
      acc[1] = pos[1] * F;
    }
  }; // struct __helper<2>
  //////////////////////////////////////////////////////////////////////////////
  template<int NDIM, typename scalar>
  inline void MiyamotoNagai::acc(const scalar*,
				 const scalar*pos,
				 const scalar*,
				 scalar      &pot,
				 scalar      *acc) const
  { 
    __helper<NDIM>::__acc(*this,pos,pot,acc);
  }
} // namespace {
//------------------------------------------------------------------------------

__DEF__ACC(MiyamotoNagai)
__DEF__POT(MiyamotoNagai)

//------------------------------------------------------------------------------
    
      
