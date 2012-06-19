// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    inc/public/random.h                                                
///                                                                             
/// \author  Walter Dehnen                                                      
/// \author  Paul McMillan                                                      
///                                                                             
/// \date    1994-2007                                                          
///                                                                             
/// \todo    add doxygen documentation                                          
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 1994-2007  Walter Dehnen                                       
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
// class PseudoRandom          Random3 OR NEMO's random (if using NEMO)         
// class Random                PseudoRandom OR Sobol                            
//                                                                              
////////////////////////////////////////////////////////////////////////////////
#ifndef falcON_included_random_h
#define falcON_included_random_h

#ifndef falcON_included_basic_h
#  include <public/basic.h>
#endif
#ifndef WDutils_included_random_h
#  include <random.h>
#endif

#ifdef falcON_NEMO
extern "C" {
  int set_xrandom(int);
  int init_xrandom(char*);
  double xrandom(double,double);
}
#endif

namespace falcON {
  using namespace WDutils;
  //
  // class falcON::PseudoRandom
  //
  /// supplies Random3 or NEMO's xrandom (depending on falcON_NEMO)
  //
#ifdef falcON_NEMO

  class PseudoRandom : public RandomNumberGenerator {
  public:
    typedef long  seed_type;
    double RandomDouble() const { return xrandom(0.,1.); }
    explicit PseudoRandom(seed_type  s) { set_xrandom(s); }
    explicit PseudoRandom(const char*s) { init_xrandom(const_cast<char*>(s)); }
  };

#else

  typedef Random3 PseudoRandom;

#endif

  //
  // class falcON::Random
  //
  /// supplies falcON::PseudoRandom AND falcON::Sobol[]
  //
  class Random : public PseudoRandom {
  private:
    const unsigned N;
    const Sobol   *S;
  public:
    /// construction
    /// \param[in] s random seed for pseudo RNG
    /// \param[in] n number of Sobol RNGs to be generated
    Random(seed_type s, unsigned  n)
      : PseudoRandom(s), N(n), S(falcON_NEW(Sobol,n)) {}
#ifdef falcON_NEMO
    /// construction
    /// \param[in] s random seed encoded as character string
    /// \param[in] n number of Sobol RNGs to be generated
    Random(const char*s, unsigned n)
      : PseudoRandom(s), N(n), S(falcON_NEW(Sobol,n)) {}
#endif
    /// destruction
    ~Random() { falcON_DEL_A(S); }
    /// pseudo random numbers
    using PseudoRandom::operator();
    /// number of Sobol RNGs
    unsigned const& Nsob() const { return N; }
    /// quasi random number in [0,1]
    /// \param[in] i  index for Sobol RNG
    double operator()(int i) const {
      return (S+i)->RandomDouble();
    }
    /// quasi random number in [a,b]
    /// \param[in] i  index for Sobol RNG
    /// \param[in] a  lower interval limit
    /// \param[in] b  upper interval limit
    double operator()(int i, double a, double b) const
    {
      return a<b? a + (b-a) * (S+i)->RandomDouble() :
	          b + (a-b) * (S+i)->RandomDouble() ;
    }
    /// give RNG
    /// \param[in] i index of Sobol RNG
    /// \param[in] q return Sobol or PseudoRanom?
    /// \note if q=true, i must be in range [0,Nsob()-1]
    //--------------------------------------------------------------------------
    const RandomNumberGenerator*rng(int i, bool q) const
    {
      return q?
	static_cast<const RandomNumberGenerator*>(S+i) :
	static_cast<const RandomNumberGenerator*>(this);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
} // namespace falcON {
//-----------------------------------------------------------------------------+
#endif // falcON_included_inline_h
