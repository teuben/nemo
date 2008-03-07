// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    utils/inc/random.h                                                 
///                                                                             
/// \author  Walter Dehnen                                                      
/// \author  Paul McMillan                                                      
///                                                                             
/// \date    1994-2008                                                          
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 1994-2008  Walter Dehnen                                       
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
    typedef int  seed_type;
    /// construction random from seed
    /// \param[in] seed value to seed deterministic pseudo-random sequence
    explicit Random3(long seed);
    /// generate a random number in [0,1]
    double RandomDouble() const;
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
    mutable int           in;
    mutable unsigned long ix;
    int                   actl;
    unsigned long         bits, *v;
    double                fac;
  public:
    /// which sequence is ours?
    /// \return number of sequence
    int const& actual() const { return actl; }
    /// construction
    /// \param[in] actl (optional) frequence to implement
    /// \note By default, we implement the next available sequence for which no
    /// instance exists (yet).
    /// \note With this implementation, actl is restricted to [0,51]
    /// \param[in] bits (optional) number of bits used
    explicit Sobol(int actl=-1, int=30);
    /// destruction: frees our sequence for next construction
    virtual ~Sobol();
    /// generate random number in [0,1]
    double RandomDouble() const;
    /// reset the sequence to start again
    void Reset() { ix = in = 0; }
    /// is this a pseudo- (or quasi-) random number generator?
    bool is_pseudo() const { return false; }
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
    /// \return \f$ p(x) \f$
    /// \param[in] x potential value for random number
    double value(double x) const { return (a<=x && x<=b)? 1. : 0.; }
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
    /// \construction
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
	   const RandomNumberGenerator*r2 = 0) WDutils_THROWING;
    /// generate random number
    /// \return x normally distributed in [-oo,oo]
    double operator() () const;
    /// give parent distribution
    /// \return \f$ p(x)\f$
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
    /// \param[in] R random number generator
    explicit Gaussian2D(const RandomNumberGenerator*r) : R(r) {}
    /// generate radius of 2D normal distribution
    /// \note The return value is in [0,oo]
    double operator() () const {
      double x;
      do { x=R->RandomDouble(); } while(x>1. || x<=0.);
      return std::sqrt(-2*std::log(x));
    }
    /// give parent distribution
    /// \return \f$ p(x)\f$
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
    RandomNumberGenerator *Rn;
  public:
    /// construction
    /// \param[in] R random number generator
    /// \param[in] a (optional) decay rate: \f$ p(x) = 1/a \exp(-x/a) \f$
    explicit  Exponential(RandomNumberGenerator* R, double a=1.)
      : alf(a), Rn(R) {}
    /// generate random number in [0,oo]
    double operator() () const { return -alf * log( Rn->RandomDouble() ); }
    /// give parent distribution
    double value(double x) const { return exp(-x/alf)/alf; }
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
    const int N,N1;
    RandomNumberGenerator *R;
    double    h, hi, hqi, *Y, *P;
    double    ranvar() const;
  public:
    /// construction
    /// \param[in] R random number generator
    /// \param[in] h scale length
    explicit ExpDisk(RandomNumberGenerator*R, double h=1.);
    virtual ~ExpDisk();
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
//     void   operator() (double& x) const { x = ranvar(); }
  };
  //////////////////////////////////////////////////////////////////////////////
} // namespace WDutils {
//-----------------------------------------------------------------------------+
#endif // WDutils_included_random_h
