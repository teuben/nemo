// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/traits.h
///
/// \author Walter Dehnen
///
/// \date   2005-2007, 2010
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005-2007, 2010 Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but          
// WITHOUT ANY WARRANTY; without even the implied warranty of                   
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
// General Public License for more details.                                     
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 675
// Mass Ave, Cambridge, MA 02139, USA.
//
////////////////////////////////////////////////////////////////////////////////
#ifndef WDutils_included_traits_h
#define WDutils_included_traits_h

#ifndef WDutils_included_typeinfo_h
#  include <typeinfo>
#  define WDutils_included_typeinfo_h
#endif
#ifndef WDutils_included_exception_h
#  include <exception.h>
#endif
////////////////////////////////////////////////////////////////////////////////
#ifdef __INTEL_COMPILER
#pragma warning (disable:424) /* extra ";" ignored */
#endif
////////////////////////////////////////////////////////////////////////////////
namespace WDutils {
  // ///////////////////////////////////////////////////////////////////////////
  //
  // integer types of given size
  //
  // ///////////////////////////////////////////////////////////////////////////
  namespace meta {
    template<typename I, typename U> struct __ISIZE {
      typedef I integer_s;
      typedef U integer_u;
      static const int size_s = sizeof(integer_s);
      static const int size_u = sizeof(integer_u);
    };

    template<int> struct __ITRAITS;
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

      WDutilsStaticAssert(sizeof(integer_s) == WORDS &&
			  sizeof(integer_u) == WORDS   );

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
  // ///////////////////////////////////////////////////////////////////////////
  //
  // class WDutils::traits<type>
  //
  /// provides very basic information about the name and size of any type
  //
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T> struct traits {
    static const char*name () { return typeid(T).name(); }
  };
  //////////////////////////////////////////////////////////////////////////////
  template<typename T> struct traits<T*> {
    static const char*name () {
      return message("%s*",traits<T>::name()); }
  };
  //////////////////////////////////////////////////////////////////////////////
  template<typename T> struct traits<const T*> {
    static const char*name () {
      return message("const %s*",traits<T>::name()); }
  };
  //////////////////////////////////////////////////////////////////////////////
  /// macro  returning the name of a given type.
#define nameof(TYPE) (static_cast<const char*>(WDutils::traits<TYPE>::name()))
  //////////////////////////////////////////////////////////////////////////////
#define WDutils_TRAITS(TYPE,NAME)			\
  template<> struct traits<TYPE> {			\
    static const char  *name () { return NAME; }	\
  };
  WDutils_TRAITS(void,"void");
  WDutils_TRAITS(bool,"bool");
  WDutils_TRAITS(char,"char");
  WDutils_TRAITS(unsigned char,"unsigned char");
  WDutils_TRAITS(int,"int");
  WDutils_TRAITS(unsigned,"unsigned");
  WDutils_TRAITS(short,"short");
  WDutils_TRAITS(unsigned short,"unsigned short");
  WDutils_TRAITS(long,"long");
  WDutils_TRAITS(unsigned long,"unsigned long");
  WDutils_TRAITS(long long,"long long");
  WDutils_TRAITS(unsigned long long,"unsigned long long");
  WDutils_TRAITS(float,"float");
  WDutils_TRAITS(double,"double");
  WDutils_TRAITS(long double,"long double");
  //////////////////////////////////////////////////////////////////////////////
} // namespace WDutils {
////////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_traits_h
