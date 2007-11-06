// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/public/manip/bound_centre.cc                                    
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
// v 0.0    11/07/2006  WD created                                              
// v 0.1    06/11/2007  WD deBUGged                                             
////////////////////////////////////////////////////////////////////////////////
#include <public/defman.h>
#include <public/io.h>

namespace falcON { namespace Manipulate {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class bound_centre                                                         
  //                                                                            
  /// manipulator: sets the centre position & velocity: that of most bound body 
  ///                                                                           
  /// This manipulator finds the most bound body from those in_subset()         
  /// (default: all) and puts its position and velocity in 'xcen' and 'vcen'    
  /// to be used by subsequent manipulators, such as sphereprof.                
  ///                                                                           
  /// Meaning of the parameters:\n                                              
  /// file: write centre position to file.                                      
  ///                                                                           
  /// Usage of pointers: sets 'xcen' and 'vcen'\n                               
  /// Usage of flags:    uses in_subset()\n                                     
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class bound_centre : public manipulator {
  private:
    mutable output OUT;
    mutable vect   XCEN,VCEN;
    mutable bool   FIRST;
    //--------------------------------------------------------------------------
  public:
    const char*name    () const { return "bound_centre"; }
    const char*describe() const {
      return 
	"provides 'xcen' and 'vcen' as position and velocity of the"
        "most bound body in_subset() (default: all)";
    }
    //--------------------------------------------------------------------------
    fieldset need   () const { return fieldset::p | fieldset::basic; }
    fieldset provide() const { return fieldset::o; }
    fieldset change () const { return fieldset::o; }
    //--------------------------------------------------------------------------
    bound_centre(const double*pars,
		 int          npar,
		 const char  *file) falcON_THROWING
    : OUT  ( file ),
      XCEN ( vect(zero) ),
      VCEN ( vect(zero) ),
      FIRST( true )
    {
      if(debug(2) || npar)
	std::cerr<<" Manipulator \"bound_centre\" centre:\n"
		 <<" find most bound body in_subset() (default: all)\n";
    }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*) const;
    //--------------------------------------------------------------------------
    ~bound_centre() {}
  };// class bound_centre
  //////////////////////////////////////////////////////////////////////////////
  bool bound_centre::manipulate(const snapshot*S) const
  {
    if(FIRST && OUT) {
      FIRST = false;
      OUT  << "#\n"
	   << "# output from Manipulator \"bound_centre\"\n#\n";
      if(RunInfo::cmd_known ()) OUT<<"# command: \""<<RunInfo::cmd ()<<"\"\n";
      OUT  << "# run at "<<RunInfo::time()<<'\n';
      if(RunInfo::user_known())
	OUT<< "#     by \""<<RunInfo::user()<<"\"\n";
      if(RunInfo::host_known())
	OUT<< "#     on \""<<RunInfo::host()<<"\"\n";
      if(RunInfo::pid_known())
	OUT<<"#     pid "<<RunInfo::pid()<<'\n';
      OUT  << "#\n# "
	   << "           time  "
	   << "              x               y               z  "
	   << "             vx              vy              vz  "
	   << "\n# ---------------"
	   << "-------------------------------------------------"
	   << "-------------------------------------------------\n";
    }
    real Emin = 1.e10;
    if(S->have(fieldbit::q)) {
      LoopSubsetBodies(S,b) {
	real E = half*norm(vel(b)) + pot(b) + pex(b);
	if(E < Emin) {
	  Emin = E;
	  XCEN = pos(b);
	  VCEN = vel(b);
	}
      }
    } else {
      LoopSubsetBodies(S,b) {
	real E = half*norm(vel(b)) + pot(b);
	if(E < Emin) {
	  Emin = E;
	  XCEN = pos(b);
	  VCEN = vel(b);
	}
      }
    }
    if(OUT)
      OUT <<"  "
	  << std::setw(15) << std::setprecision(8) << S->time() << "  "
	  << std::setw(15) << std::setprecision(8) << XCEN[0]   << ' '
	  << std::setw(15) << std::setprecision(8) << XCEN[1]   << ' '
	  << std::setw(15) << std::setprecision(8) << XCEN[2]   << "  "
	  << std::setw(15) << std::setprecision(8) << VCEN[0]   << ' '
	  << std::setw(15) << std::setprecision(8) << VCEN[1]   << ' '
	  << std::setw(15) << std::setprecision(8) << VCEN[2]   << std::endl;
    S->set_pointer(&XCEN,"xcen");
    S->set_pointer(&VCEN,"vcen");
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} }

__DEF__MAN(falcON::Manipulate::bound_centre)

////////////////////////////////////////////////////////////////////////////////
