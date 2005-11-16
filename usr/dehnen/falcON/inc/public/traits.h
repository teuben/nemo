// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// traits.h                                                                    |
//                                                                             |
// Copyright (C) 2005  Walter Dehnen                                           |
//                                                                             |
// This program is free software; you can redistribute it and/or modify        |
// it under the terms of the GNU General Public License as published by        |
// the Free Software Foundation; either version 2 of the License, or (at       |
// your option) any later version.                                             |
//                                                                             |
// This program is distributed in the hope that it will be useful, but         |
// WITHOUT ANY WARRANTY; without even the implied warranty of                  |
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU           |
// General Public License for more details.                                    |
//                                                                             |
// You should have received a copy of the GNU General Public License           |
// along with this program; if not, write to the Free Software                 |
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                   |
//                                                                             |
//------------------------------------------------------------------------------
#ifndef falcON_included_traits_h
#define falcON_included_traits_h

////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class falcON::traits<type>                                               //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  template<typename T> struct traits {
    static const char  *name () { return "unknown type"; }
    static const char  *names() { return "unknown types"; }
    static const size_t size = sizeof(T);
  };
  //////////////////////////////////////////////////////////////////////////////
#define falcON_TRAITS(TYPE,NAME,NAMES)			\
  template<> struct traits<TYPE> {			\
    static const char  *name () { return NAME; }	\
    static const char  *names() { return NAMES; }	\
    static const size_t size = sizeof(TYPE);		\
  }
  falcON_TRAITS(char,"char","chars");
  falcON_TRAITS(unsigned char,"unsigned char","unsigned chars");
  falcON_TRAITS(int,"int","ints");
  falcON_TRAITS(unsigned,"unsigned","unsigneds");
  falcON_TRAITS(short,"short","shorts");
  falcON_TRAITS(unsigned short,"unsigned short","unsigned shorts");
  falcON_TRAITS(long,"long","longs");
  falcON_TRAITS(unsigned long,"unsigned long","unsigned longs");
  falcON_TRAITS(long long,"long long","long longs");
  falcON_TRAITS(unsigned long long,"unsigned long long","unsigned long longs");
  falcON_TRAITS(float,"float","floats");
  falcON_TRAITS(double,"double","doubles");
  falcON_TRAITS(long double,"long double","long doubles");
  //////////////////////////////////////////////////////////////////////////////
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_traits_h
