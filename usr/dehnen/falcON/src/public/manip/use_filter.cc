// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/public/manip/use_filter.cc                                      
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
// v 0.0    07/07/2006  WD created                                              
////////////////////////////////////////////////////////////////////////////////
#include <public/defman.h>
#include <public/bodyfunc.h>
#include <sstream>

namespace falcON { namespace Manipulate {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class use_filter                                                           
  //                                                                            
  /// manipulator: using a BodyFilter in 'filter' to define subset of all bodies
  ///                                                                           
  /// This manipulator filters the bodies with a BodyFilter in 'filter'. Bodies 
  /// passing the filter will be in_subset(). \n                                
  /// In order to set 'filter', use manipulator set_filter.\n                   
  ///                                                                           
  /// Usage of pointers: uses 'filter'                                          
  /// Usage of flags:    sets subset flag (in fact flags::ignore)               
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class use_filter : public manipulator {
    mutable std::string DS;
  public:
    const char*name    () const { return "use_filter"; }
    const char*describe() const {
      return "chooses subset of bodies according to 'filter'";
    }
    //--------------------------------------------------------------------------
    fieldset need   () const { return fieldset::o; }
    fieldset provide() const { return fieldset::f; }
    fieldset change () const { return fieldset::f; }
    //--------------------------------------------------------------------------
    use_filter(const char*p, const char*f) falcON_THROWING  {
      if(p) warning("Manipulator \"%s\": parameters \"%s\" ignored",name(),p);
      if(f) warning("Manipulator \"%s\": file \"%s\" ignored",name(),f);
    }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*S) const {
      const BodyFilter*BF = S->pointer<BodyFilter>("filter");
      if(BF && *BF) {
	// make sure flags are supported
	if(!S->have(fieldbit::f))
	  const_cast<snapshot*>(S)->add_field(fieldbit::f);
	// check if all data needed are supported
	if(!S->have_all(BF->need()))
	  falcON_THROW("set_subset::manipulate(): "
		       "filter needs '%s' but snapshot has only '%s'\n",
		       word(BF->need()), word(S->all_data()));
	// loop bodies and run them through the filter
	unsigned sub(0);
	LoopAllBodies(S,b)
	  if((*BF)(b)) {
	    b.into_subset();
	    ++sub;
	  } else
	    b.outof_subset();
	// create description of subset
	std::ostringstream ost;
	ost  << sub << " bodies with \""
	     << BF->expression() << '\"';
	if(BF->npar()) {
	  ost<< " where";
	  for(int n=0; n!=BF->npar(); ++n)
	    ost<<" #"<<n<<'='<<BF->param(n);
	}
	DS = ost.str();
	// put pointer to description into pointer bank
	S->set_pointer(&DS,"subset_description");
      }
      return false;
    }
  };
} }
////////////////////////////////////////////////////////////////////////////////
__DEF__MAN__ALT(falcON::Manipulate::use_filter);
