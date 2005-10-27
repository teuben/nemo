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
    static const char  *name() { return "unknown type"; }
    static const size_t size = sizeof(T);
  };
  //////////////////////////////////////////////////////////////////////////////
#define falcON_TRAITS(TYPE,NAME)		\
  template<> struct traits<TYPE> {		\
    static const char  *name() { return NAME; }	\
    static const size_t size = sizeof(TYPE);	\
  }
  falcON_TRAITS(char,"char");
  falcON_TRAITS(unsigned char,"unsigned char");
  falcON_TRAITS(int,"int");
  falcON_TRAITS(unsigned,"unsigned");
  falcON_TRAITS(short,"short");
  falcON_TRAITS(unsigned short,"unsigned short");
  falcON_TRAITS(long,"long");
  falcON_TRAITS(unsigned long,"unsigned long");
  falcON_TRAITS(long long,"long long");
  falcON_TRAITS(unsigned long long,"unsigned long long");
  falcON_TRAITS(float,"float");
  falcON_TRAITS(double,"double");
  falcON_TRAITS(long double,"long double");
  //////////////////////////////////////////////////////////////////////////////
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_traits_h
