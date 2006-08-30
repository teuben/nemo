// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    utils/inc/traits.h                                                 
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    2005-2006                                                          
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2005-2006  Walter Dehnen                                       
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

#ifndef WDutils_included_basic_h
# include <exception.h>
#endif
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
