// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   inc/public/types.h                                                  
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2000-2005                                                           
///                                                                             
/// \brief  contains declaration of some basic types, such as falcON::real      
///         and falcON::vect, and constants.                                    
///                                                                             
////////////////////////////////////////////////////////////////////////////////
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
#ifndef falcON_included_types_h
#define falcON_included_types_h

#ifndef falcON_included_iostream
#  include <iostream>
#  define falcON_included_iostream
#endif
#ifndef falcON_included_cmath
#  include <cmath>
#  define falcON_included_cmath
#endif
#ifndef falcON_included_cstdlib
#  include <cstdlib>
#  define falcON_included_cstdlib
#endif
#ifndef falcON_included_basic_h
#  include <public/basic.h>
#endif

////////////////////////////////////////////////////////////////////////////////
//                                                                              
// general falcON macros                                                        
//                                                                              
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
// SPH code?                                                                    
//------------------------------------------------------------------------------
#if defined(falcON_SPH) && !defined(falcON_INDI)
#  warning 'falcON_SPH' #defined but not 'falcON_INDI'
#  warning will #define 'falcON_INDI' now
#  define falcON_INDI
#endif

//------------------------------------------------------------------------------
// adaptive individual softening?                                               
//------------------------------------------------------------------------------
#if defined(falcON_ADAP) && !defined(falcON_INDI)
#  warning
#  warning 'falcON_ADAP' #defined but not 'falcON_INDI' 
#  warning we will #undef 'falcON_ADAP'
#  warning
#endif
//------------------------------------------------------------------------------
// expansion order used in gravity                                              
//------------------------------------------------------------------------------
#if defined(Walter) && !defined(falcON_included_ordr_h)
#  include <walter/ordr.h>
#endif

#ifndef falcON_ORDER
#  define falcON_ORDER 3
#endif
//------------------------------------------------------------------------------
// falcON_REAL_IS_FLOAT  and type falcON::real                                  
//------------------------------------------------------------------------------
namespace falcON {
#if defined(falcON_DOUBLE)
#  undef falcON_REAL_IS_FLOAT
#  define falcON_REAL double
#  define falcON_NOTREAL float
#else
#  define falcON_REAL_IS_FLOAT
#  define falcON_REAL float
#  define falcON_NOTREAL double
#endif
  /// floating point type used by falcON for body and internal data;
  /// determined at compile time by macro falcON_DOUBLE or falcON_FLOAT.
  typedef falcON_REAL real;
  /// the "other" type (double if real is float, otherwise float)
  typedef falcON_NOTREAL notreal;
#undef falcON_REAL
#undef falcON_NOTREAL
}
//------------------------------------------------------------------------------
// SSE code?                                                                    
//------------------------------------------------------------------------------
#if defined(falcON_REAL_IS_FLOAT) && defined(falcON_SSE)
#  define  falcON_SSE_CODE
#else
#  undef   falcON_SSE_CODE
#endif
//------------------------------------------------------------------------------
// Dimensionality                                                               
//------------------------------------------------------------------------------
#ifdef falcON_NDIM
#  if falcON_NDIM != 3
#    error 
#    error falcON only works in 3D
#    error 
#  endif
#  warning
#  warning falcON_NDIM deprecated, falcON is ALWAYS 3D
#  warning
#  undef falcON_NDIM
#endif
namespace falcON {
  const int Ndim = 3;                            ///< falcON is 3D only         
  const int Nsub = 1<<Ndim;                      ///< 2^(# dimensions)          
}
//------------------------------------------------------------------------------
// macros used to test timings                                                  
//------------------------------------------------------------------------------
#ifdef TEST_TIMING
#  include <ctime>
#  define SET_I        std::clock_t __C0_TIMING = std::clock();
#  define SET_T(TEXT)  std::cerr<<TEXT					\
                                <<(std::clock()-__C0_TIMING)/		\
                                   double(CLOCKS_PER_SEC)<<'\n';	\
                       __C0_TIMING = std::clock();
#else
#  define SET_I        {}
#  define SET_T(TEXT)  {}
#endif

////////////////////////////////////////////////////////////////////////////////
//                                                                              
// elementary falcON types and constants                                        
//                                                                              
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
//  vectors and tensors                                                         
//------------------------------------------------------------------------------
#ifndef falcON_included_inline_h
#  include <public/inline.h>
#endif
#ifndef falcON_included_tupel_h
#  include <public/tupel.h>
#endif
#ifndef falcON_included_tensor_h
#  include <public/tensor.h>
#endif

namespace falcON {
  typedef tupel<Ndim,real  > vect;               ///< a vector of 3 reals       
  typedef tupel<Ndim,float>  vect_f;             ///< a vector of 3 floats      
  typedef tupel<Ndim,double> vect_d;             ///< a vector of 3 doubles     
  //----------------------------------------------------------------------------
  falcON_TRAITS(vect,"vect","vects");
#ifdef falcON_REAL_IS_FLOAT
  falcON_TRAITS(vect_d,"vect_d","vect_d");
#else
  falcON_TRAITS(vect_f,"vect_f","vect_f");
#endif
  //----------------------------------------------------------------------------
}
//------------------------------------------------------------------------------
// integer types of given size                                                  
//------------------------------------------------------------------------------
namespace falcON {
  namespace meta {
    template<typename I, typename U> struct __ISIZE {
      typedef I integer_s;
      typedef U integer_u;
      static const int size_s = sizeof(integer_s);
      static const int size_u = sizeof(integer_u);
    };

    template<int I> struct __ITRAITS {};
    template<> struct __ITRAITS<0> : 
      public __ISIZE<char, unsigned char> {};
    template<> struct __ITRAITS<1> :
      public __ISIZE<short, unsigned short> {};
    template<> struct __ITRAITS<2> :
      public __ISIZE<int, unsigned int> {};
    template<> struct __ITRAITS<3> :
      public __ISIZE<long, unsigned long> {};
    template<> struct __ITRAITS<4> :
      public __ISIZE<long long, unsigned long long> {};

    template<int WORDS> struct __IWORDS {
      static const int STYPE = 
    __ITRAITS<0>::size_s == WORDS ? 0 :
    __ITRAITS<1>::size_s == WORDS ? 1 :
    __ITRAITS<2>::size_s == WORDS ? 2 :
    __ITRAITS<3>::size_s == WORDS ? 3 :
    __ITRAITS<4>::size_s == WORDS ? 4 : 5;
      typedef typename __ITRAITS<STYPE>::integer_s integer_s;
      static const int UTYPE = 
    __ITRAITS<0>::size_u == WORDS ? 0 :
    __ITRAITS<1>::size_u == WORDS ? 1 :
    __ITRAITS<2>::size_u == WORDS ? 2 :
    __ITRAITS<3>::size_u == WORDS ? 3 :
    __ITRAITS<4>::size_u == WORDS ? 4 : 5;
      typedef typename __ITRAITS<UTYPE>::integer_u integer_u;
    };
  } // namespace meta {

  typedef meta::__IWORDS<1>::integer_s int8;   ///<   signed integer of  8 bytes
  typedef meta::__IWORDS<2>::integer_s int16;  ///<   signed integer of 16 bytes
  typedef meta::__IWORDS<4>::integer_s int32;  ///<   signed integer of 32 bytes
  typedef meta::__IWORDS<8>::integer_s int64;  ///<   signed integer of 64 bytes
  typedef meta::__IWORDS<1>::integer_u uint8;  ///< unsigned integer of  8 bytes
  typedef meta::__IWORDS<2>::integer_u uint16; ///< unsigned integer of 16 bytes
  typedef meta::__IWORDS<4>::integer_u uint32; ///< unsigned integer of 32 bytes
  typedef meta::__IWORDS<8>::integer_u uint64; ///< unsigned integer of 64 bytes
} // namespace falcON {
//------------------------------------------------------------------------------
// useful constants and typedefs                                                
//------------------------------------------------------------------------------
namespace falcON {
  typedef uint16 indx;                         ///< unsigned integer of 16 bytes
  typedef uint64 peanokey;                     ///< type of Peano-Hilbert key   
  /// \name some general floating-point constants                               
  //@{
  //----------------------------------------------------------------------------
  const real zero       = 0.,                         ///< real: zero           
             sixth      = 0.166666666666666666666667, ///< real: 1/6            
             fifth      = 0.2,                        ///< real: 1/5            
             quarter    = 0.25,                       ///< real: 1/4            
             third      = 0.333333333333333333333333, ///< real: 1/3            
             half       = 0.5,                        ///< real: 1/2            
             one        = 1.,                         ///< real: 1              
             threehalfs = 1.5,                        ///< real: 3/2            
             two        = 2.,                         ///< real: 2              
             three      = 3.,                         ///< real: 3              
             four       = 4.,                         ///< real: 4              
             six        = 6.,                         ///< real: 6              
             eight      = 8.,                         ///< real: 8              
             ten        = 10.,                        ///< real: 10             
             twelve     = 12.;                        ///< real: 12             
  //@}
  //----------------------------------------------------------------------------
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_types_h
