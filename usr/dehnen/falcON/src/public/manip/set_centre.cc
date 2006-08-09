// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/public/manip/set_centre.cc                                      
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2006                                                                
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2006 Walter Dehnen                                             
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
//                                                                              
// history:                                                                     
//                                                                              
// v 0.0    30/05/2006  WD created                                              
// v 0.1    31/05/2006  WD velocity centre added                                
// v 1.0    04/07/2006  WD reset (rather than set to origin) if not given       
////////////////////////////////////////////////////////////////////////////////
#include <public/defman.h>
#include <public/basic.h>
#include <public/io.h>

namespace falcON { namespace Manipulate {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class set_centre                                                           
  //                                                                            
  /// manipulator: sets 'xcen' and 'vcen' to parameters 0-2 and 3-5             
  ///                                                                           
  /// This manipulator simply sets the position and velocity centre pointers    
  /// 'xcen' and 'vcen' to refer to the vectors given by the parameters 0-2 and 
  /// 3-5, respectively. If no parameters are given, 'xcen' and 'vcen' are      
  /// reset to NULL pointers (no centre known, which for most subsquent         
  /// manipulators implies using the origin). The number of parameters must     
  /// either be 0 ('xcen' and 'vcen' will be reset), or 3 ('xcen' will be set   
  /// and 'vcen' reset) or 6 (both 'xcen' and 'vcen' will be set).              
  ///                                                                           
  /// Usage of pointers: sets 'xcen' and 'vcen'.                                
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // class set_centre                                                         //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  class set_centre : public manipulator {
  private:
    vect XC,VC,*X0,*V0;
    //--------------------------------------------------------------------------
  public:
    const char* name    () const { return "set_centre"; }
    const char* describe() const {
      if(X0)
	if(V0)
	  return message("sets 'xcen' to (%f,%f,%f)"
			 " and 'vcen' to (%f,%f,%f)",
			 XC[0],XC[1],XC[2], VC[0],VC[1],VC[2]);
	else
	  return message("sets 'xcen' to (%f,%f,%f)"
			 " and resets 'vcen' to null\n",
			 XC[0],XC[1],XC[2]);
      else
	return "resets 'xcen' and 'vcen' to null";
    }
    //--------------------------------------------------------------------------
    fieldset need   () const { return fieldset::o; }
    fieldset provide() const { return fieldset::o; }
    fieldset change () const { return fieldset::o; }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*) const;
    //--------------------------------------------------------------------------
    set_centre(const double*pars,
	       int          npar,
	       const char  *file) :
      X0(0), V0(0)
    {
      if(npar == 0 || debug(2)) {
	std::cerr<<
	  " Manipulator \"set_centre\":\n"
	  " provides 'xcen'=(par[0],par[1],par[2]) (default: none)\n"
	  " and      'vcen'=(par[3],par[4],par[5]) (default: none)\n";
      }
      if(npar!=0 && npar!=3 && npar<6)
	error("Manipulator \"set_centre\": #pars must be 0,3, or 6\n");
      if(npar >= 3) {
	XC = vect(pars);
	X0 = &XC;
      }
      if(npar >= 6)  {
	VC = vect(pars+3);
	V0 = &VC;
      }
    }
    //--------------------------------------------------------------------------
    ~set_centre() {}
  };
  //////////////////////////////////////////////////////////////////////////////
  bool set_centre::manipulate(const snapshot*S) const {
    S->set_pointer(X0,"xcen");
    S->set_pointer(V0,"vcen");
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} }

__DEF__MAN(falcON::Manipulate::set_centre)

////////////////////////////////////////////////////////////////////////////////
    
