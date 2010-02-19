// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/public/manip/dens_centre.cc                                     
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2006,2008                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2006,2008 Walter Dehnen                                        
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
// v 0.0    21/06/2006  WD created                                              
// v 0.1    26/06/2006  WD set_pointer() & replace pointer()                    
// v 0.2    27/06/2006  WD using 'subset'                                       
// v 1.0    04/07/2006  WD made public (along with the tools used)              
// v 1.1    06/07/2006  WD using 'filter' rather than 'subset'                  
// v 1.2    06/07/2006  WD using useful bodies, rather than filter function     
// v 2.0    11/07/2006  WD warning and no output out if no converging           
// v 2.0.1  11/06/2008  WD new DebugInfo and falcON_Warning                     
////////////////////////////////////////////////////////////////////////////////
#include <public/defman.h>
#include <public/tools.h>

namespace falcON { namespace Manipulate {
  //////////////////////////////////////////////////////////////////////////////
  const int W_default = 500;
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class dens_centre                                                          
  //                                                                            
  /// manipulator: determines the density centre of all bodies in_subset()      
  ///                                                                           
  /// This manipulator computes the position at which the density, estimated in 
  /// an SPH-like manner, is maximal. Only bodies in_subset() (default: all) are
  /// considered (see set_subset and use_filter). The centre position and       
  /// velocity are put into 'xcen' and 'vcen' to be used by subsequent          
  /// manipulators, such as sphereprof.                                         
  ///                                                                           
  /// \note The algorithm determining the centre may not always converge, in    
  /// particular if the particle distribution has a flat core. In this case,    
  /// a warning is written to stderr, no output is made to file and no centre   
  /// position and velocity are put into 'xcen' and 'vcen'.                     
  ///                                                                           
  /// Meaning of the parameters:\n                                              
  /// par[0]: number of neighbours in SPH-like estimation (def: 500)\n          
  /// file: write centre position to file.                                      
  ///                                                                           
  /// Usage of pointers: sets 'xcen' and 'vcen'\n                               
  /// Usage of flags:    uses in_subset()\n                                     
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class dens_centre : public manipulator {
  private:
    const   int    W;
    mutable output OUT;
    mutable real   RAD;
    mutable vect   XCEN,VCEN;
    mutable bool   FIRST;
    //--------------------------------------------------------------------------
  public:
    const char*name    () const { return "dens_centre"; }
    const char*describe() const {
      return 
	"provides 'xcen' and 'vcen' as density centre "
	"for bodies in_subset() (default: all)";
    }
    //--------------------------------------------------------------------------
    fieldset need   () const { return fieldset::basic; }
    fieldset provide() const { return fieldset::empty; }
    fieldset change () const { return fieldset::empty; }
    //--------------------------------------------------------------------------
    dens_centre(const double*pars,
		int          npar,
		const char  *file) falcON_THROWING
    : W    ( npar>0? int(pars[0]) : W_default ),
      OUT  ( file ),
      RAD  ( one ),
      XCEN ( vect(zero) ),
      VCEN ( vect(zero) ),
      FIRST( true )
    {
      if(debug(2) || npar>1)
	std::cerr<<" Manipulator \"dens_centre\" centre:\n"
		 <<" find density maximum for bodies in_subset()"
		 <<" (default: all)\n"
		 <<" parameter:\n"
		 <<" par[0]: # bodies in centre (def: "<<W_default<<")\n"
		 <<" file  : write centre time pos, vel & density\n";
      if(W <0) falcON_THROW("Manipulator \"dens_centre\": W = %d < 0\n",W);
    }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*) const;
    //--------------------------------------------------------------------------
    ~dens_centre() {}
  };// class dens_centre
  //////////////////////////////////////////////////////////////////////////////
  bool dens_centre::manipulate(const snapshot*S) const
  {
    if(FIRST && debug(1))
      DebugInfo(2,"dens_centre::manipulate(): "
		"first call: obtain initial guess using tree:\n");
    if(FIRST) {
      // get initial guess for centre position
      // (after first manipulation, we will use the previous centre)
      clock_t CPU0 = clock();
      flags F = flags::empty;
      if(S->N_bodies() != S->N_subset()) {
	DebugInfo(2,"dens_centre::manipulate(): subset < all, so"
		  "must first mark subset bodies\n");
	F = flags::marked;
	if(!S->have(fieldbit::f))
	  const_cast<snapshot*>(S)->add_field(fieldbit::f);
	LoopAllBodies(S,b)
	  if(in_subset(b)) b.mark(); else b.unmark();
      }
      DebugInfo(2,"dens_centre::manipulate(): "
		"   ... then build a tree ... \n");
      OctTree TREE(S,W/4,0,Default::MaxDepth,F);
      if(debug(1)) {
	clock_t CPU1 = clock();
	DebugInfo(1,"dens_centre::manipulate(): "
		  "   tree build took %f sec\n",
		   (CPU1 - CPU0)/real(CLOCKS_PER_SEC));
	CPU0 = CPU1;
      }
      DebugInfo(2,"dens_centre::manipulate():    ... and finally"
		" estimate position of density peak ...\n");
      estimate_density_peak(&TREE,0u,W,XCEN,RAD);
      RAD *= three;
      FIRST = false;
      if(OUT) {
	OUT  << "#\n"
	     << "# output from Manipulator \"dens_centre\"\n#\n";
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
	     << "        density"
	     << "\n# ---------------"
	     << "-------------------------------------------------"
	     << "-------------------------------------------------"
	     << "------------------\n";
      }
    }
    // any call: refine centre estimate
    RAD *= three;
    real RHO;
    DebugInfo(2,"dens_centre::manipulate(): initial guess: (%f,%f,%f)\n",
	      XCEN[0],XCEN[1],XCEN[2]);
    if(find_density_centre(S,W,XCEN,RAD,&VCEN,&RHO)) {
      DebugInfo(2,"dens_centre::manipulate():    found density "
		"centre: x0=(%f,%f,%f), v0=(%f,%f,%f), radius=%f, rho=%f\n",
		XCEN[0],XCEN[1],XCEN[2],VCEN[0],VCEN[1],VCEN[2],RAD,RHO);
      if(OUT)
	OUT <<"  "
	    << std::setw(15) << std::setprecision(8) << S->time() << "  "
	    << std::setw(15) << std::setprecision(8) << XCEN[0]   << ' '
	    << std::setw(15) << std::setprecision(8) << XCEN[1]   << ' '
	    << std::setw(15) << std::setprecision(8) << XCEN[2]   << "  "
	    << std::setw(15) << std::setprecision(8) << VCEN[0]   << ' '
	    << std::setw(15) << std::setprecision(8) << VCEN[1]   << ' '
	    << std::setw(15) << std::setprecision(8) << VCEN[2]   << "  "
	    << std::setw(15) << std::setprecision(8) << RHO       << std::endl;
      // register pointers to centre
      S->set_pointer(&XCEN,"xcen");
      S->set_pointer(&VCEN,"vcen");
    } else {
      falcON_WarningN("dens_centre::manipulate(): no convergence at time %f\n",
		      S->time());
      if(OUT)
	OUT << "# WARNING: could not find centre at time " 
	    << S->time() << std::endl;
      S->del_pointer("xcen");
      S->del_pointer("vcen");
    }
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} }

__DEF__MAN(falcON::Manipulate::dens_centre)

////////////////////////////////////////////////////////////////////////////////
