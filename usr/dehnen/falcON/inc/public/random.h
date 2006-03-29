// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    inc/random.h                                                       
///                                                                             
/// \author  Walter Dehnen                                                      
/// \author  Paul McMillan                                                      
///                                                                             
/// \date    1994-2005                                                          
///                                                                             
/// \todo    add doxygen documentation                                          
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 1994-2005  Walter Dehnen                                       
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
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::PseudoRandom                                               //
  //                                                                          //
  // supplies Random3 or NEMO's xrandom (depending on falcON_NEMO)            //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
#ifdef falcON_NEMO

  class PseudoRandom : public RandomNumberGenerator {
  public:
    typedef long  seed_type;
    double RandomDouble() const { return xrandom(0.,1.); }
    explicit PseudoRandom(seed_type const&s) { set_xrandom(s); }
    explicit PseudoRandom(const char*s) { init_xrandom(const_cast<char*>(s)); }
  };

#else

  typedef Random3 PseudoRandom;

#endif
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::Random                                                     //
  //                                                                          //
  // supplies falcON::PseudoRandom AND falcON::Sobol[]                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class Random : public PseudoRandom {
  private:
    const unsigned  N;
    const Sobol    *S;
    //--------------------------------------------------------------------------
    // construction:                                                            
    //--------------------------------------------------------------------------
  public:
    Random(seed_type const&s,                      // I: seed for pseudo RNG    
	   unsigned  const&n) :                    // I: # Sobols               
//       PseudoRandom(s), N(n), S( new Sobol[n] ) {}
      PseudoRandom(s), N(n), S(falcON_NEW(Sobol,n)) {}
#ifdef falcON_NEMO
    Random(const     char *s,                      // I: seed for NEMO::xrandom 
	   unsigned  const&n) :                    // I: # Sobols               
//       PseudoRandom(s), N(n), S( new Sobol[n] ) {}
      PseudoRandom(s), N(n), S(falcON_NEW(Sobol,n)) {}
#endif
    //--------------------------------------------------------------------------
    ~Random() { delete[] S; }
    //--------------------------------------------------------------------------
    // pseudo random numbers                                                    
    //--------------------------------------------------------------------------
    PseudoRandom::operator();
    //--------------------------------------------------------------------------
    // quasi random numbers: must give No of Sobol to be used                   
    //--------------------------------------------------------------------------
    unsigned const& Nsob() const { return N; }
    //--------------------------------------------------------------------------
    double operator()(                             // R: ith RNG in (0,1)       
		      int    const&i) const {      // I: i                      
      return (S+i)->RandomDouble();
    }
    //--------------------------------------------------------------------------
    double operator()(                             // R: ith RNG in (a,b)       
		      int    const&i,              // I: i                      
		      double const&a,              // I: a = lower limit        
		      double const&b) const {      // I: b = upper limit        
      return a<b? a + (b-a) * (S+i)->RandomDouble() :
	          b + (a-b) * (S+i)->RandomDouble() ;
    }
    //--------------------------------------------------------------------------
    // miscellaneous                                                            
    //--------------------------------------------------------------------------
    const RandomNumberGenerator*rng(               // R: pter to RNG            
			            int  const&i,  // I: No of RNG              
			            bool const&q)  // I: quasi or pseudo?       
      const {
      if(q) return S+i;
      else  return this;
    }
  };
  //////////////////////////////////////////////////////////////////////////////
} // namespace falcON {
//-----------------------------------------------------------------------------+
#endif // falcON_included_inline_h
