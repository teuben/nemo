// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    inc/exception.h                                                    
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    2005                                                               
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2005  Walter Dehnen                                            
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
#ifndef WDutils_included_traits_h
#define WDutils_included_traits_h

////////////////////////////////////////////////////////////////////////////////
namespace WDutils {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class WDutils::traits<type>                                              //
  //                                                                          //
  /// provides very basic information about the name and size of any type     //
  //                                                                          //
  // ///////////////////////////////////////////////////////////////////////////
  template<typename T> struct traits {
    static const char  *name () { return "unknown type"; }
    static const char  *names() { return "unknown types"; }
    static const unsigned size = sizeof(T);
  };
  //////////////////////////////////////////////////////////////////////////////
#define WDutils_TRAITS(TYPE,NAME,NAMES)			\
  template<> struct traits<TYPE> {			\
    static const char  *name () { return NAME; }	\
    static const char  *names() { return NAMES; }	\
    static const unsigned size = sizeof(TYPE);		\
  };
  WDutils_TRAITS(char,"char","chars");
  WDutils_TRAITS(unsigned char,"unsigned char","unsigned chars");
  WDutils_TRAITS(int,"int","ints");
  WDutils_TRAITS(unsigned,"unsigned","unsigneds");
  WDutils_TRAITS(short,"short","shorts");
  WDutils_TRAITS(unsigned short,"unsigned short","unsigned shorts");
  WDutils_TRAITS(long,"long","longs");
  WDutils_TRAITS(unsigned long,"unsigned long","unsigned longs");
  WDutils_TRAITS(long long,"long long","long longs");
  WDutils_TRAITS(unsigned long long,"unsigned long long","unsigned long longs");
  WDutils_TRAITS(float,"float","floats");
  WDutils_TRAITS(double,"double","doubles");
  WDutils_TRAITS(long double,"long double","long doubles");
  //////////////////////////////////////////////////////////////////////////////
} // namespace WDutils {
////////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_traits_h
