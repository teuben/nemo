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
// Copyright (C) 2005-2007, 2010-2012 Walter Dehnen
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

//
#ifndef WDutils_included_typeinfo
#  include <typeinfo>
#  define WDutils_included_typeinfo
#endif
#ifndef WDutils_included_exception_h
#  include <exception.h>
#endif

#if __cplusplus >= 201103L
# ifndef WDutils_included_cstdint
#  include <cstdint>
#  define WDutils_included_cstdint
# endif
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
#if __cplusplus >= 201103L
  using std::int8_t;
  using std::int16_t;
  using std::int32_t;
  using std::int64_t;
  using std::uint8_t;
  using std::uint16_t;
  using std::uint32_t;
  using std::uint64_t;
  namespace meta {
    template<int WORDS> struct IntTypeWords;
    template<> struct IntTypeWords<1>
    { typedef  int8_t  signed_integer; typedef uint8_t  unsigned_integer; };
    template<> struct IntTypeWords<2>
    { typedef  int16_t signed_integer; typedef uint16_t unsigned_integer; };
    template<> struct IntTypeWords<4>
    { typedef  int32_t signed_integer; typedef uint32_t unsigned_integer; };
    template<> struct IntTypeWords<8>
    { typedef  int64_t signed_integer; typedef uint64_t unsigned_integer; };
  }
#else
  namespace meta {
    template<typename I, typename U> struct _ISIZE {
      typedef I signed_integer;
      typedef U unsigned_integer;
      WDutilsStaticAssert(sizeof(signed_integer) == sizeof(unsigned_integer));
      static const int  size = sizeof(signed_integer);
    };

    template<int> struct _ITRAITS;
    template<> struct _ITRAITS<0> : 
      public _ISIZE<char, unsigned char> {};
    template<> struct _ITRAITS<1> :
      public _ISIZE<short, unsigned short> {};
    template<> struct _ITRAITS<2> :
      public _ISIZE<int, unsigned int> {};
    template<> struct _ITRAITS<3> :
      public _ISIZE<long, unsigned long> {};
    template<> struct _ITRAITS<4> :
      public _ISIZE<long long, unsigned long long> {};

    template<int WORDS> struct IntTypeWords
    {
    private:
      static const int TYPE = 
    _ITRAITS<0>::size == WORDS ? 0 :
    _ITRAITS<1>::size == WORDS ? 1 :
    _ITRAITS<2>::size == WORDS ? 2 :
    _ITRAITS<3>::size == WORDS ? 3 :
    _ITRAITS<4>::size == WORDS ? 4 : 5;
    public:
      typedef typename _ITRAITS<TYPE>::signed_integer signed_integer;
      typedef typename _ITRAITS<TYPE>::unsigned_integer unsigned_integer;
    private:
      WDutilsStaticAssert(sizeof(signed_integer) == WORDS);
    };
  }// namespace meta {
  typedef meta::IntTypeWords<1>::signed_integer int8_t;
  typedef meta::IntTypeWords<2>::signed_integer int16_t;
  typedef meta::IntTypeWords<4>::signed_integer int32_t;
  typedef meta::IntTypeWords<8>::signed_integer int64_t;
  typedef meta::IntTypeWords<1>::unsigned_integer uint8_t;
  typedef meta::IntTypeWords<2>::unsigned_integer uint16_t;
  typedef meta::IntTypeWords<4>::unsigned_integer uint32_t;
  typedef meta::IntTypeWords<8>::unsigned_integer uint64_t;
#endif
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
  template<int K, typename T> struct traits<T[K]> {
    static const char*name ()
    { return message("%s[%d]",traits<T>::name(),K); }
    static const char*fmt  () { return "%p"; }
  };
  //
  template<int K, typename T> struct traits<const T[K]> {
    static const char*name ()
    { return message("const %s[%d]",traits<T>::name(),K); }
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
