// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/public/manip/set_filter.cc                                      
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
// v 0.0    06/07/2006  WD created                                              
////////////////////////////////////////////////////////////////////////////////
#include <public/defman.h>
#include <public/bodyfunc.h>

namespace falcON { namespace Manipulate {
  // this macro controls the behaviour in case of a parse or syntax error       
#undef WARNING_ON_ERROR
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class set_filter                                                           
  //                                                                            
  /// manipulator: creates a BodyFilter and registers it under 'filter'.        
  ///                                                                           
  /// This manipulator creates a BodyFilter from a boolean bodyfunc expression  
  /// (see man pages (5bodyfunc), given in manipfile, with parameters from      
  /// manippars and registers it under 'filter' for usage by subsequent         
  /// manipulators. \n                                                          
  /// A simple usage is to restrict the range of bodies. For instance the       
  /// bodyfunc expression \n                                                    
  ///      "#0<=i && i<#1" \n                                                   
  /// filters bodies with index between parameters #0 and #1 (note that the     
  /// the indices of bodies may not be preserved; in this case use the key      
  /// instead, that is 'k' instead of 'i' in the expression).\n                 
  /// An empty expression makes the open filter: all bodies are accepted.       
  ///                                                                           
  /// Usage of pointers: sets 'filter'                                          
  /// Usage of flags:    none                                                   
  ///                                                                           
  /// \note If the expression cannot be parsed to a boolean filter function,    
  ///       a fatal error will result.                                          
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class set_filter : public manipulator {
  private:
    BodyFilter *BF;
    //--------------------------------------------------------------------------
  public:
    const char*name    () const { return "set_filter"; }
    const char*describe() const {
      if(BF && *BF) 
	return message("filters bodies according to filter \"%s\" "
		       "with parameters %s",
		       BF->expression(), BF->parameters());
      else return "open filter: no effect\n";
    }
    //--------------------------------------------------------------------------
    fieldset need   () const { 
      return (BF && *BF)? BF->need() - fieldset::k : fieldset::o;
    }
    fieldset provide() const { return fieldset::o; }
    fieldset change () const { return fieldset::o; }
    //--------------------------------------------------------------------------
    set_filter(const char*params, const char*filter) falcON_THROWING
    : BF(0) {
      try {
	BF = new BodyFilter(filter, params);
      } catch (falcON::exception E) {
#ifdef WARNING_ON_ERROR
	falcON::warning("Manipulator \"%s\": "
			"could not generate a filter from expression \"%s\" "
			"(error: \"%s\"); "
			"will make an open filter (accepting all bodies)\n",
			name(),filter,text(E));
	BF = 0;
#else
	falcON_THROW("Manipulator \"%s\": "
		     "could not generate a filter from expression \"%s\" "
		     "(error: \"%s\")\n",name(),filter,text(E));
#endif
      }
    }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*S) const
    {
      if(BF) {
	// if body keys are needed, make sure they are supported
	if(BF->need(fieldbit::k) && !S->have(fieldbit::k))
	  const_cast<snapshot*>(S)->add_field(fieldbit::f);
	// set time (filter may depend on time)
	BF->set_time(S->time());
      }
      S->set_pointer(BF,"filter");
      return false;
    }
    //--------------------------------------------------------------------------
    ~set_filter() {
      if(BF) falcON_DEL_O(BF);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
} }
////////////////////////////////////////////////////////////////////////////////
__DEF__MAN__ALT(falcON::Manipulate::set_filter);
