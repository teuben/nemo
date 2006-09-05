// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    utils/inc/random.h                                                 
///                                                                             
/// \author  Walter Dehnen                                                      
/// \author  Paul McMillan                                                      
///                                                                             
/// \date    1994-2006                                                          
///                                                                             
/// \todo    add doxygen documentation                                          
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 1994-2006  Walter Dehnen                                       
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
// class RandomNumberGenerator base class for random number generators          
// class Random3               pseudo-random number generator (Press et al.)    
// class PseudoRandom          Random3 OR NEMO's random (if using NEMO)         
// class Sobol                 quasi-random number generator (Press et al.)     
// class Random                PseudoRandom OR Sobol                            
//                                                                              
// class RandomDeviate         base class for a random distribution             
// class Uniform               P(x) = 1/ b-a  for x in (a,b)                    
// class Gaussian              P(x) = Exp[-x^2/(2 sigma^2)]; x in [-oo,oo]      
// class Gaussian1D            P(x) =     Exp[-x^2/2]; x in [-oo,oo]            
// class Gaussian2D            P(x) = x * Exp[-x^2/2]; x in [0  ,oo]            
// class Exponential           P(x) = Exp[-x]; x in [0,oo)                      
// class ExpDisk               P(x) = x Exp[-x]; x in [0,oo)                    
//                                                                              
////////////////////////////////////////////////////////////////////////////////
#ifndef WDutils_included_random_h
#define WDutils_included_random_h

#ifndef WDutils_included_iostream
# define WDutils_included_iostream
# include <iostream>
#endif
#ifndef WDutils_included_cmath
# define WDutils_included_cmath
# include <cmath>
#endif
#ifndef WDutils_included_inline_h
#  include <inline.h>
#endif
#ifndef WDutils_included_Pi_h
#  include <Pi.h>
#endif
#ifndef WDutils_included_traits_h
#  include <traits.h>
#endif

namespace WDutils {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class WDutils::RandomNumberGenerator                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class RandomNumberGenerator {
  public:
    virtual double RandomDouble()  const = 0;
    void   RandomDouble(double& x) const { x = RandomDouble(); }
    double operator()  ()          const { return RandomDouble(); }
    void   operator()  (double& x) const { RandomDouble(x); }
    double operator()  (double const&a,
		        double const&b) const {
      return a<b? a + (b-a) * RandomDouble() :
	          b + (a-b) * RandomDouble() ;
    }
  }; 
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class WDutils::RandomDeviate                                             //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class RandomDeviate {
  public:
    virtual double operator() () const = 0;
    virtual double value      (const double) const = 0;
  }; 
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class WDutils::Random3                                                   //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class Random3 : public RandomNumberGenerator {
    // pseudo-random number generator
  private:
    mutable int  inext, inextp;
    mutable long ma[56];
  public:
    typedef int  seed_type;
    explicit Random3(long const&);
    double RandomDouble() const;
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class WDutils::Sobol                                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class Sobol : public RandomNumberGenerator {
    // quasi-random number generator
  private:
    mutable int           in;
    mutable unsigned long ix;
    int                   actl;
    unsigned long         bits, *v;
    double                fac;
  public:
    int const& actual() const { return actl; }
    explicit Sobol(const int=-1, const int=0);
    virtual ~Sobol();
    double RandomDouble() const;
    void   Reset       ()   { ix = in = 0; }
    int    actual      ()   { return actl; }
  };

  WDutils_TRAITS(Sobol,"Sobol");
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class WDutils::Uniform                                                   //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class Uniform : public RandomDeviate {
    // gives x uniformly in [a,b]
  private:
    RandomNumberGenerator &r;
    double   a,b,ba;
  public:
    explicit
    Uniform(RandomNumberGenerator* R,                 // random number generator
	    const double A=0., const double B=1.)     // a,b                    
      : r(*R), a(min(A,B)), b(max(A,B)), ba(b-a) {}
    double lower_bound() { return a; }
    double upper_bound() { return b; }
    double operator() ()         const { return a+ba*r(); }
    double value(const double x) const { return (a<=x && x<=b)? 1. : 0.; }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class WDutils::Gaussian                                                  //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class Gaussian : public RandomDeviate {
    // gives x in [-oo,oo] with probability proportional to exp[-x^2/(2sigma^2)]
  private:
    mutable int    iset;
    mutable double gset;
    double         sig, norm;
    const RandomNumberGenerator *R1, *R2;
  public:
    Gaussian(const RandomNumberGenerator*,            // 1st RNG                
	     const RandomNumberGenerator*,            // 2nd RNG                
	     const double=1.);                        // sigma                  
    double operator() () const;
    double sigma   () { return sig; }
    double value(const double) const;
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class WDutils::Gaussian1D                                                //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class Gaussian1D : public RandomDeviate {
    // gives x in [0,oo] with probability proportional to exp[-x^2/2]           
  private:
    mutable bool   iset;
    mutable double gset;
    const RandomNumberGenerator *R1;
    const RandomNumberGenerator *R2;
  public:
    Gaussian1D(const RandomNumberGenerator*r1,
	       const RandomNumberGenerator*r2) : iset(0), R1(r1), R2(r2) {}
    double operator() () const {
      if(iset) {
	iset = 0;
	return gset;
      } else {
	register double v1,v2,rsq,fac;
	do {
	  v1  = 2 * R1->RandomDouble() - 1.;
	  v2  = 2 * R2->RandomDouble() - 1.;
	  rsq = v1*v1 + v2*v2;
	} while (rsq>=1. || rsq <=0. );
	fac  = sqrt(-2.*log(rsq)/rsq);
	gset = v1*fac;
	iset = 1;
	return v2*fac;
      }
    }
    double value(const double x) const {
      return std::exp(-0.5*square(x)) / WDutils::STPi;
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class WDutils::Gaussian2D                                                //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class Gaussian2D : public RandomDeviate {
    // gives x in [0,oo] with probability proportional to x*exp[-x^2/2]         
  private:
    const RandomNumberGenerator *R;
  public:
    explicit 
    Gaussian2D(const RandomNumberGenerator*r) : R(r) {}
    double operator() () const {
      register double x=R->RandomDouble();
      while(x>1. || x<=0.) x=R->RandomDouble();
      return std::sqrt(-2*std::log(x));
    }
    double value(const double x) const {
      return x * std::exp(-0.5*square(x));
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class WDutils::Exponential                                               //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class Exponential : public RandomDeviate {
    // gives x in [0,oo] with probability proportional to Exp[-x/alpha]         
  private:
    double alf;
    RandomNumberGenerator *Rn;
  public:
    explicit 
    Exponential(RandomNumberGenerator* R,             // RNG                    
                const double a=1.)                    // alpha                  
      : alf(a), Rn(R) {}
    double operator() () const;
    double value(const double) const;
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class WDutils::ExpDisk                                                   //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class ExpDisk: public RandomDeviate {
    // gives x in [0,oo] with probability proportional to  x*Exp[-x/h]
  private:
    const int N,N1;
    RandomNumberGenerator *R;
    double    h, hi, hqi, *Y, *P;
    double    ranvar() const;
  public:
    explicit 
    ExpDisk(RandomNumberGenerator*,                   // RNG                    
	    const double=1.);                         // h                      
    virtual ~ExpDisk();
    double value      (const double) const;
    double radius     (const double) const;
    double operator() ()          const { return ranvar(); }
    void   operator() (double& x) const { x = ranvar(); }
  };
  //////////////////////////////////////////////////////////////////////////////
} // namespace WDutils {
//-----------------------------------------------------------------------------+
#endif // WDutils_included_random_h
