// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/public/manip/lagrange.cc                                        
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2004-2006                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2004-2006 Walter Dehnen                                        
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
// v 1.0    27/06/2006  WD using bodyset in 'subset'                            
// v 1.1    27/06/2006  WD using filter in 'filter' instead of 'subset'         
// v 1.2    07/07/2006  WD using flags::ignore instead of filter                
////////////////////////////////////////////////////////////////////////////////
#include <public/defman.h>
#include <public/basic.h>
#include <public/tools.h>
#include <public/io.h>
#include <ctime>
#include <cmath>

namespace falcON { namespace Manipulate {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class lagrange                                                             
  //                                                                            
  /// manipulator: computes Lagrange radii for subset, writes them to file      
  ///                                                                           
  /// This manipulator computes the Lagrange radii w.r.t. centre 'xcen' for     
  /// all bodies in_subset() (default: all, see set_subset), and writes them    
  /// to file.                                                                  
  ///                                                                           
  /// Meaning of the parameters:\n                                              
  /// par[0-n]: mass fractions for which Lagrange radii shall be computed     \n
  /// file:     file name for output of table with Lagrange radii             \n
  ///                                                                           
  /// Usage of pointers: uses 'xcen'\n                                          
  /// Usage of flags:    uses in_subset()\n                                     
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class lagrange : public manipulator {
  private:
    int                 N;
    double             *M;
    mutable double     *R;
    mutable output      OUT;
    mutable bool        FST;
    //--------------------------------------------------------------------------
  public:
    const char* name    () const { return "lagrange"; }
    const char* describe() const {
      return
	"computes lagrange radii of bodies passing 'filter' (default: all) "
	"w.r.t. \"xcen\" (default: origin)";
    }
    //--------------------------------------------------------------------------
    fieldset need   () const { return fieldset::m | fieldset::x; }
    fieldset provide() const { return fieldset::o; }
    fieldset change () const { return fieldset::o; }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*) const;
    //--------------------------------------------------------------------------
    lagrange(const double*pars,
	     int          npar,
	     const char  *file) :
      N   ( npar ),
      M   ( npar>0? falcON_NEW(double,npar) : 0 ),
      R   ( npar>0? falcON_NEW(double,npar) : 0 ),
      OUT ( file? file : "." ),
      FST ( true )
    {
      for(int i=0; i!=N; ++i) M[i] = pars[i];
      if(file == 0 || npar == 0) {
	std::cerr<<
	  " Manipulator \""<<name()<<"\":\n"
	  " computes lagrange radii w.r.t. 'xcen' (default: origin)\n"
	  " at given mass fractions (parameters) for bodies passing\n"
	  " 'filter' (default: all) and writes them given data file.\n";
      }
      if(file == 0)
	error("Manipulator \"%s\": no output file given\n",name());
      if(!OUT.is_open())
	error("Manipulator \"%s\": couldn't open output\n",name());
      if(npar == 0)
	error("Manipulator \"%s\": no mass fractions given\n",name());
    }
    //--------------------------------------------------------------------------
    ~lagrange() {
      if(M) falcON_DEL_A(M); 
      if(R) falcON_DEL_A(R);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  bool lagrange::manipulate(const snapshot*S) const {
    if(S->have_not(need())) {
      warning("Manipulator \"%s\" insufficient data\n",name());
      return false;
    }
    if(FST) {
      FST = false;
      OUT  << "#\n"
	   << "# output from Manipulator \"lagrange\"\n";
      if(RunInfo::cmd_known())
	OUT<< "#\n# command: \""<<RunInfo::cmd() <<"\"\n";
      OUT  << "# run at "<<RunInfo::time()<<'\n';
      if(RunInfo::user_known())
	OUT<< "#     by \""<<RunInfo::user()<<"\"\n";
      if(RunInfo::host_known())
	OUT<< "#     on \""<<RunInfo::host()<<"\"\n";
      if(RunInfo::pid_known())
	OUT<< "#     pid "<<RunInfo::pid()<<'\n';
      OUT  << "#\n# time      ";
      for(int i=0; i!=N; ++i)
	OUT << " r[" << std::setw(4) << 100*M[i] << "%]";
      OUT  << "\n#-----------";
      for(int i=0; i!=N; ++i)
	OUT << "---------";
      OUT << '\n';
      OUT.stream().setf(std::ios::left, std::ios::adjustfield);
    }
    const vect*X0= S->pointer<vect>("xcen");
    find_lagrange_rad(S,N,M,R,X0);
    OUT  << std::setw(12) << S->time();
    for(int i=0; i!=N; ++i)
      OUT<< ' ' << std::setw(8) << R[i];
    OUT  << std::endl;
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} }

__DEF__MAN(falcON::Manipulate::lagrange)
