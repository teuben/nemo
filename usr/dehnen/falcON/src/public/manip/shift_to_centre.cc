// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/public/manip/shift_to_centre.cc                                 
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2007                                                                
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2007 Walter Dehnen                                             
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
#include <public/defman.h>

namespace falcON { namespace Manipulate {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class shift_to_centre                                                      
  //                                                                            
  /// manipulator: shifts snapshot to centre 'xcen' and 'vcen'                  
  ///                                                                           
  /// No parameters or file used.                                               
  ///                                                                           
  /// Usage of pointers: uses 'xcen' and 'vcen'\n                               
  /// Usage of flags:    none\n                                                 
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class shift_to_centre : public manipulator {
  public:
    const char*name    () const { return "shift_to_centre"; }
    const char*describe() const {
      return 
	"shifts snapshot to 'xcen' and 'vcen' (default: origin)";
    }
    //--------------------------------------------------------------------------
    fieldset need   () const { return fieldset::empty; }
    fieldset provide() const { return fieldset::empty; }
    fieldset change () const { return fieldset::empty; }
    //--------------------------------------------------------------------------
    shift_to_centre(const double*pars,
		    int          npar,
		    const char  *file) falcON_THROWING
    {
      if(debug(2) || npar  || file)
	std::cerr<<" Manipulator \"shift_to_centre\":\n"
		 <<" shifts snapshot to 'xcen' and 'vcen' (default: origin)\n";
    }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*) const;
    //--------------------------------------------------------------------------
    ~shift_to_centre() {}
    //--------------------------------------------------------------------------
  };  
  //////////////////////////////////////////////////////////////////////////////
  bool shift_to_centre::manipulate(const snapshot*S) const
  {
    const vect*C = S->pointer<vect>("xcen");
    if(C && S->have_pos())
      LoopAllBodies(S,B)
	B.pos() -= *C;
    
    C = S->pointer<vect>("vcen");
    if(C && S->have_vel())
      LoopAllBodies(S,B)
	B.vel() -= *C;
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} }

__DEF__MAN(falcON::Manipulate::shift_to_centre)

////////////////////////////////////////////////////////////////////////////////
