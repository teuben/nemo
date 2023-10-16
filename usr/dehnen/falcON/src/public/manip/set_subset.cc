// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/public/manip/set_subset.cc                                      
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
// v 0.0    26/06/2006  WD created                                              
// v 0.1    04/07/2006  WD renamed to set_subset (previously: subset)           
// v 1.0    04/07/2006  WD BodyFilter define subset; flags::ignore to mark it   
////////////////////////////////////////////////////////////////////////////////
#include <public/defman.h>
#include <public/bodyfunc.h>
#include <sstream>

namespace falcON {
// /////////////////////////////////////////////////////////////////////////////
//                                                                              
/// \brief                                                                      
/// namespace falcON::Manipulate contains implementations of                    
/// falcON::manipulator prepared for loading at run time via                    
/// falcON::Manipulator                                                         
//                                                                              
// /////////////////////////////////////////////////////////////////////////////
namespace Manipulate {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class set_subset                                                           
  //                                                                            
  /// manipulator: using a BodyFilter to define subset of all bodies            
  ///                                                                           
  /// This manipulator creates a BodyFilter from a boolean bodyfunc expression  
  /// (see man pages (5bodyfunc), given in manipfile, with parameters from      
  /// manippars. At each manipulation, bodies not passing the filter are flagged
  /// to be ignored, while those passing the filter will be flagged not to be   
  /// ignored and hence are in_subset().\n                                      
  /// This effects subsequent analysing manipulators which only use bodies      
  /// in_subset().\n                                                            
  /// Essentially, set_subset is like set_filter followed by use_filter.\n      
  /// A simple usage is to restrict the range of bodies. For instance the       
  /// bodyfunc expression                                                       
  ///      "#0<=i&&i<#1"                                                      
  /// filters bodies with index between parameters #0 and #1 (note that the     
  /// the indices of bodies may not be preserved; in this case use the key      
  /// instead, that is 'k' instead of 'i' in the expression).\n                 
  /// An empty expression makes the open filter: all bodies are accepted.       
  ///                                                                           
  /// Meaning of the parameters:\n                                              
  /// par[] : parameters used in bodyfunc expression (if any)\n                 
  /// file  : bodyfunc expression\n                                             
  ///                                                                           
  /// Usage of pointers: none\n                                                 
  /// Usage of flags:    sets subset flag (in fact flags::ignore)\n             
  ///                                                                           
  /// \note If the expression cannot be parsed to a boolean filter function,    
  ///       a fatal error will result.                                          
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class set_subset : public manipulator {
  private:
    BodyFilter*BF;
    mutable std::string DS;
    mutable char DESC[1024];
    //--------------------------------------------------------------------------
  public:
    const char*name    () const { return "set_subset"; }
    const char*describe() const {
      if(BF && *BF) {
	if(DESC[0]==0)
	  sprintf(DESC,"chooses subset of bodies according to filter \"%s\" "
		  "with parameters %s",BF->expression(), BF->parameters());
	return DESC;
      } else
	return "chooses subset of bodies: all bodies";
    }
    //--------------------------------------------------------------------------
    fieldset need() const { 
      return (BF && *BF)?  BF->need() - fieldset::k : fieldset::empty;
    }
    fieldset provide() const { return fieldset::f; }
    fieldset change () const { return fieldset::f; }
    //--------------------------------------------------------------------------
    /// construction from boolean bodyfunc expression
    /// \param filter boolean bodyfunc expression
    /// \param params parameters (if any) required by filter
    set_subset(const char*params, const char*filter) falcON_THROWING 
    : BF(0) {
      DESC[0]=0;
      try {
	BF = new BodyFilter(filter, params);
      } catch (falcON::exception& E) {
#ifdef WARNING_ON_ERROR
	falcON::warning("Manipulator \"%s\": "
			"could not generate a filter from expression \"%s\" "
			"(error: \"%s\"); "
			"will take an open filter (accepting all bodies)\n",
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
    /// manipulate snapshot: flag filtered bodies to be "in_subset()"
    /// \param S snapshot
    /// \return always false (to continue simulation)
    bool manipulate(const snapshot*S) const
    {
      if(BF && *BF) {
	// make sure flags are supported
	if(!S->have(fieldbit::f))
	  const_cast<snapshot*>(S)->add_field(fieldbit::f);
	// if body keys are needed, make sure they are supported
	if(BF->need(fieldbit::k) && !S->have(fieldbit::k))
	  const_cast<snapshot*>(S)->add_field(fieldbit::f);
	// check if all data needed are supported
	if(!S->have_all(BF->need()))
	  falcON_THROW("set_subset::manipulate(): "
		       "filter needs '%s' but snapshot has only '%s'\n",
		       word(BF->need()), word(S->all_data()));
	// set time (filter may depend on time)
	BF->set_time(S->time());
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
    //--------------------------------------------------------------------------
    ~set_subset() {
      if(BF) falcON_DEL_O(BF);
    }
  };
} }
////////////////////////////////////////////////////////////////////////////////
__DEF__MAN__ALT(falcON::Manipulate::set_subset);
