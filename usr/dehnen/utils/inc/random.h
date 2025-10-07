// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    utils/inc/random.h                                                 
///                                                                             
/// \author  Walter Dehnen                                                      
/// \author  Paul McMillan                                                      
///                                                                             
/// \date    1994-2011                                                          
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 1994-2011  Walter Dehnen                                       
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
#ifndef WDutils_included_random_h
#define WDutils_included_random_h

#ifndef WDutils_included_iostream
#  define WDutils_included_iostream
#  include <iostream>
#endif
#ifndef WDutils_included_cmath
#  define WDutils_included_cmath
#  include <cmath>
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

#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wweak-vtables"
#endif
namespace WDutils {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class WDutils::RandomNumberGenerator                                       
  //                                                                            
  /// abstract base class for a random number generator                         
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class RandomNumberGenerator {
  public:
    /// generate a random number
    /// \return a random number in interval [0,1]
    virtual double RandomDouble()  const = 0;
    /// generate a random number
    /// \param[out] x a random number in interval [0,1]
    void   RandomDouble(double& x) const { x = RandomDouble(); }
    /// generate a random number
    /// \return a random number in interval [0,1]
    double operator()  ()          const { return RandomDouble(); }
    /// generate a random number
    /// \param[out] x a random number in interval [0,1]
    void   operator()  (double& x) const { RandomDouble(x); }
    /// generate a random number in given interval
    /// \return    x a random number in interval [a,b]
    /// \param[in] a lower interval limit
    /// \param[in] b upper interval limit
    double operator() (double a, double b) const {
      return a<b? a + (b-a) * RandomDouble() :
	          b + (a-b) * RandomDouble() ;
    }
    /// is this a pseudo- (or quasi-) random number generator?
    virtual bool is_pseudo() const { return true; }
    /// noop dtor
    virtual~RandomNumberGenerator() {}
  }; 
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class WDutils::RandomDeviate                                               
  //                                                                            
  /// abstract base class for a random distribution                             
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class RandomDeviate {
  public:
    /// generate a random number from a distribution
    virtual double operator() () const = 0;
    /// compute value of parent distribution
    /// \param[in] x value at which to compute \f$ p(x)\f$
    /// \return value of \f$ p(x)\f$
    virtual double value(double x) const = 0;
    /// noop dtor
    virtual~RandomDeviate() {}
  }; 
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class WDutils::Random3                                                     
  //                                                                            
  /// instance of a (pseudo-) random number generator.                          
  ///                                                                           
  /// This is essentially an algorithm from Numercial Recipes                   
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class Random3 : public RandomNumberGenerator {
  private:
    mutable int  inext, inextp;
    mutable long ma[56];
  public:
    typedef long seed_type;
    /// construction random from seed
    /// \param[in] seed value to seed deterministic pseudo-random sequence
    explicit Random3(seed_type seed);
    /// generate a random number in [0,1]
    double RandomDouble() const;
    /// noop dtor
    virtual~Random3() {}
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class WDutils::Sobol                                                       
  //                                                                            
  /// instance of a (quasi-) random number generator.                           
  ///                                                                           
  /// Quasi random numbers are not meant to be random, but to fill in the gaps  
  /// left by previous numbers of the sequence. This is very useful, if Poisson 
  /// noise is to be avoided, see Numerical Recipes for details.                
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class Sobol : public RandomNumberGenerator {
  private:
    mutable unsigned      in;
    mutable unsigned long ix;
    unsigned              actl,bits;
    unsigned long         *v;
    unsigned long         *v_allocated;
    double                fac;
  public:
    /// which sequence is ours?
    /// \return number of sequence
    unsigned const& actual() const { return actl; }
    /// construction
    /// \param[in] _actl (optional) frequence to implement
    /// \note By default, we implement the next available sequence for which no
    /// instance exists (yet).
    /// \note With this implementation, actl is restricted to [0,51]
    /// \param[in] _bits (optional) number of bits used
    explicit Sobol(int _actl=-1, unsigned _bits=30);
    /// destruction: frees our sequence for next construction
    virtual ~Sobol();
    /// generate random number in [0,1]
    double RandomDouble() const;
    /// reset the sequence to start again
    void Reset() { ix = in = 0; }
    /// is this a pseudo- (or quasi-) random number generator?
    bool is_pseudo() const { return false; }
    /// sets the default value for number of bits used, if 2nd arg to ctor is 0
    static void set_bits(const unsigned);
  };
  WDutils_TRAITS(Sobol,"Sobol");
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class WDutils::Uniform                                                     
  //                                                                            
  /// random distribution: gives x uniformly in [a,b]                           
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class Uniform : public RandomDeviate {
  private:
    RandomNumberGenerator &r;
    double   a,b,ba;
  public:
    /// construction
    /// \param[in] R (pointer to a) random number generator
    /// \param[in] _a (optional) lower limit of interval
    /// \param[in] _b (optional) upper limit of interval
    explicit Uniform(RandomNumberGenerator* R, double _a=0., double _b=1.)
      : r(*R), a(min(_a,_b)), b(max(_a,_b)), ba(b-a) {}
    /// give lower bound
    double lower_bound() { return a; }
    /// give upper bound
    double upper_bound() { return b; }
    /// generate random number
    /// \return a randomm number in [a,b]
    double operator() () const { return a+ba*r(); }
    /// compute parent distribution
    /// \return distribution \f$ p(x) \f$
    /// \param[in] x potential value for random number
    double value(double x) const { return (a<=x && x<=b)? 1/ba : 0.; }
  };
  WDutils_TRAITS(Uniform,"Uniform");
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class WDutils::Normal                                                      
  //                                                                            
  /// random distribution: the normal distribution                              
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class Normal : public RandomDeviate {
  private:
    mutable int    iset;
    mutable double gset;
    const RandomNumberGenerator *R1, *R2;
  public:
    /// construction
    /// \param[in] r1 (pointer to a) random number generator
    /// \param[in] r2 (optional) (pointer to a) random number generator,
    ///               default is to use r1
    /// \note The algorithm generates two random numbers at a time (and re-uses
    /// the second one at the next call) by sampling a 2D normal distribution.
    /// One random number generator is used for the polar angle, the other for
    /// the radius. If quasi-random numbers are wanted, these to random number
    /// generator \b must not be encoding identical sequences. However, with 
    /// ordinary (pseudo-) random numbers, r1 and r2 may well be the same. The
    /// code shall issue a fatal error if r1==r2 and is not a pseudo-RNG.
    Normal(const RandomNumberGenerator*r1,
	   const RandomNumberGenerator*r2 = 0);
    /// generate random number
    /// \return x normally distributed in [-oo,oo]
    double operator() () const;
    /// give parent distribution
    /// \return distribution \f$ p(x)\f$
    /// \param[in] x potential value for random number
    double value(double x) const
    {
      return std::exp(-0.5*square(x)) / WDutils::STPi;
    }
  };
  WDutils_TRAITS(Normal,"Normal")
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class WDutils::Gaussian2D                                                  
  //                                                                            
  /// random distribution: radius of a 2D normal distribution                   
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class Gaussian2D : public RandomDeviate {
  private:
    const RandomNumberGenerator *R;
  public:
    /// construction
    /// \param[in] r random number generator
    explicit Gaussian2D(const RandomNumberGenerator*r) : R(r) {}
    /// generate radius of 2D normal distribution
    /// \note The return value is in [0,oo]
    double operator() () const {
      double x;
      do { x=R->RandomDouble(); } while(x>1. || x<=0.);
      return std::sqrt(-2*std::log(x));
    }
    /// give parent distribution
    /// \return distribution \f$ p(x)\f$
    /// \param[in] x potential value for random number
    double value(const double x) const {
      return x * std::exp(-0.5*square(x));
    }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class WDutils::Exponential                                                 
  //                                                                            
  /// random distribution: \f$ p(x) = 1/a \exp(-x/a) \f$                        
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class Exponential : public RandomDeviate {
  private:
    double alf;
    const RandomNumberGenerator *Rn;
  public:
    /// construction
    /// \param[in] R random number generator
    /// \param[in] a (optional) decay rate: \f$ p(x) = 1/a \exp(-x/a) \f$
    explicit  Exponential(const RandomNumberGenerator* R, double a=1.)
      : alf(a), Rn(R) {}
    /// generate random number in [0,oo]
    double operator() () const { return -alf * log( Rn->RandomDouble() ); }
    /// give parent distribution
    double value(double x) const { return exp(-x/alf)/alf; }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class WDutils::SechSquared                                                 
  //                                                                            
  /// random distribution: \f$ p(x) = sech^2(x) / 2 \f$                        
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class SechSquared : public RandomDeviate {
  private:
    const RandomNumberGenerator *Rn;
  public:
    /// construction
    /// \param[in] R random number generator
    explicit SechSquared(const RandomNumberGenerator* R)
      : Rn(R) {}
    /// generate random number in [-oo,oo]
    double operator() () const
    {
#if(0)
      return std::atanh( 2*Rn->RandomDouble()-1 );
#else
      double p,x;
      do {
	p=Rn->RandomDouble();        // p in [0,1]
	x=0.5*std::log(p/(1-p));
      } while(isinf(x) || isnan(x));
      return x;
#endif
    }
    /// give parent distribution f(x)=(1/2)*sech^2(x)
    double value(double x) const
    {
      double ex=x<0? std::exp(x) : std::exp(-x);
      return 2*square(ex/(1+ex*ex));
    }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class WDutils::ExpDisk                                                     
  //                                                                            
  /// random distribution: \f$ p(r) \propto r \exp(-r/h) \f$                    
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class ExpDisk: public RandomDeviate {
    // gives x in [0,oo] with probability proportional to  x*Exp[-x/h]
  private:
    static const unsigned N=256, N1=N+1;
    const RandomNumberGenerator*const R;
    const double h,hi,hqi;
    double Y[N1],P[N1];
    double ranvar() const;
  public:
    /// construction
    /// \param[in] rng          random number generator
    /// \param[in] scale_height scale length
    explicit ExpDisk(const RandomNumberGenerator*rng, double scale_height=1.);
    virtual ~ExpDisk() {};
    /// give parent distribution
    /// \return p(r)
    /// \param[in] r potential random value
    double value (double r) const;
    /// inverse of p(r)
    /// \return r
    /// \param[in] p probability density p(r)
    double radius(double p) const;
    /// generate random number in [0,oo]
    double operator() () const { return ranvar(); }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //
  // class WDutils::PowerLawDist
  // 
  /// random distribution: \f$ p(r) \propto r^alpha \f$
  //
  // ///////////////////////////////////////////////////////////////////////////
  class PowerLawDist: public RandomDeviate
  {
    const RandomNumberGenerator*R;
    const double xmin,xmax;
    const double p,p1,ip1;
    const bool   islog;
    const double ranfc,pnorm;
    double ranvar() const;
  public:
    /// ctor
    /// \param[in] rng   random number generator
    /// \param[in] alpha power
    /// \param[in] xmin  mininum value for x
    /// \param[in] xmax  maximum valud for x
    PowerLawDist(const RandomNumberGenerator*rng, double alpha,
		 double xmin, double xmax);
    /// dtor
    virtual ~PowerLawDist() {}
    /// give parent distribution
    /// \return p(r)
    /// \param[in] x potential random value
    double value (double x) const
    { return pnorm*std::pow(x,p); }
    /// generate random number in [xmin,xmax]
    double operator() () const
    { return ranvar(); }
  };// PowerLawDist
  //////////////////////////////////////////////////////////////////////////////
} // namespace WDutils {
#ifdef __clang__
#  pragma clang diagnostic pop
#endif
//-----------------------------------------------------------------------------+
#endif // WDutils_included_random_h
