// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/traits.h
///
/// \author Walter Dehnen
///
/// \date   2005-2007, 2010-2011
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005-2007, 2010-2011 Walter Dehnen
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
      WDutilsStaticAssert(sizeof(integer_s) == sizeof(integer_u));
      static const int  size = sizeof(integer_s);
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

    template<int WORDS> struct __IWORDS
    {
    private:
      static const int TYPE = 
    __ITRAITS<0>::size == WORDS ? 0 :
    __ITRAITS<1>::size == WORDS ? 1 :
    __ITRAITS<2>::size == WORDS ? 2 :
    __ITRAITS<3>::size == WORDS ? 3 :
    __ITRAITS<4>::size == WORDS ? 4 : 5;
    public:
      typedef typename __ITRAITS<TYPE>::integer_s integer_s;
      typedef typename __ITRAITS<TYPE>::integer_u integer_u;
    private:
      WDutilsStaticAssert(sizeof(integer_s) == WORDS);
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
  //
  // class WDutils::traits<type>
  //
  /// provides the name and suitable printf-style format for any type
  //
  template<typename T> struct traits {
    static const char*name () { return typeid(T).name(); }
    static const char*fmt  () { return 0; }
  };
  //
  template<typename T> struct traits<T*> {
    static const char*name () { return message("%s*",traits<T>::name()); }
    static const char*fmt  () { return "%p"; }
  };
  //
  template<typename T> struct traits<const T*> {
    static const char*name () { return message("const %s*",traits<T>::name()); }
    static const char*fmt  () { return "%p"; }
  };
  //
  template<> struct traits<char*> {
    static const char*name () { return "char*"; }
    static const char*fmt  () { return "%s"; }
  };
  //
  template<> struct traits<const char*> {
    static const char*name () { return "const char*"; }
    static const char*fmt  () { return "%s"; }
  };
  //
  /// macro  returning the name of a given type.
#define nameof(TYPE) (static_cast<const char*>(WDutils::traits<TYPE>::name()))
  /// macro  returning the suitable printf-style format for a given type
#define fmtof(TYPE) (static_cast<const char*>(WDutils::traits<TYPE>::fmt()))
  //
#define WDutils_TRAITS_FORMAT(TYPE,NAME,FMT)	\
  template<> struct traits<TYPE> {		\
    static const char*name () { return NAME; }	\
    static const char*fmt  () { return FMT; }	\
  };
  WDutils_TRAITS_FORMAT(void,"void",0);
  WDutils_TRAITS_FORMAT(bool,"bool","%d");
  WDutils_TRAITS_FORMAT(char,"char","%c");
  WDutils_TRAITS_FORMAT(unsigned char,"unsigned char","%uc");
  WDutils_TRAITS_FORMAT(int,"int","%d");
  WDutils_TRAITS_FORMAT(unsigned,"unsigned","%ud");
  WDutils_TRAITS_FORMAT(short,"short","%sd");
  WDutils_TRAITS_FORMAT(unsigned short,"unsigned short","%sud");
  WDutils_TRAITS_FORMAT(long,"long","%ld");
  WDutils_TRAITS_FORMAT(unsigned long,"unsigned long","%lud");
  WDutils_TRAITS_FORMAT(long long,"long long","%lld");
  WDutils_TRAITS_FORMAT(unsigned long long,"unsigned long long","%llud");
  WDutils_TRAITS_FORMAT(float,"float","%g");
  WDutils_TRAITS_FORMAT(double,"double","%g");
  WDutils_TRAITS_FORMAT(long double,"long double","%g");
#define WDutils_TRAITS(TYPE,NAME)		\
  template<> struct traits<TYPE> {		\
    static const char*name () { return NAME; }	\
    static const char*fmt  () { return 0; }	\
  };
  //////////////////////////////////////////////////////////////////////////////
} // namespace WDutils {
////////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_traits_h
