// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    inc/public/utils.h                                                 
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    2005                                                               
///                                                                             
/// \brief   includes utilities from WDutils into namespace falcON              
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
#ifndef falcON_included_utils_h
#define falcON_included_utils_h 1

// include WDutils
#ifndef WDutils_included_traits_h
#  include <traits.h>
#endif
#ifndef WDutils_included_Pi_h
#  include <Pi.h>
#endif
#ifndef WDutils_included_inline_h
#  include <inline.h>
#endif
#ifndef WDutils_included_inline_io_h
#  include <inline_io.h>
#endif
#ifndef WDutils_included_tupel_h
#  include <tupel.h>
#endif
#ifndef WDutils_included_exception_h
#  include <exception.h>
#endif
#ifndef WDutils_included_memory_h
#  include <memory.h>
#endif
// declare WDutils in namespace falcON
namespace falcON {
  using namespace WDutils;
  using WDutils::pow;
  using WDutils::exception;
  using WDutils::snprintf;
  namespace meta {
    using namespace WDutils::meta;
  }
}
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_utils_h
