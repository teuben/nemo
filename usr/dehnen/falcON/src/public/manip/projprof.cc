// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/public/manip/projprof.cc                                        
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2006,2010                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2006,2010 Walter Dehnen                                        
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
// v 0.0    30/06/2006  WD created                                              
// v 0.1    03/07/2006  JS changed default behavious                            
// v 0.2    04/07/2006  WD made public (along with the tools used)              
// v 0.3    07/07/2006  WD using flags::ignore (in_subset()) instead of subset  
// v 0.3.1  11/06/2008  WD new DebugInfo                                        
// v 0.3.2  13/04/2010  WD set initial file index to non-existing file
////////////////////////////////////////////////////////////////////////////////
#include <public/defman.h>
#include <public/profile.h>

namespace falcON { namespace Manipulate {
  //////////////////////////////////////////////////////////////////////////////
  const int    W_default = 500;
  const double L_default = 0.1;
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class projprof                                                             
  //                                                                            
  /// manipulator: estimated projected radial profiles for bodies in_subset()   
  ///                                                                           
  /// This manipulator estimates the projected radial profiles of surface       
  /// density; mean, rotational, and dispersion of the line-of-sight velocity;  
  /// the axis ratio; position and rotation angle for all bodies in_subset()    
  /// (default: all, see set_subset) w.r.t. 'xcen' and 'vcen'.\n                
  /// See falcON::projected_profile for details.                                
  ///                                                                           
  /// Meaning of the parameters:\n                                              
  /// par[0-2]: position P. line of sight has direction P - 'xcen'              
  /// par[3]: minimum # bodies in radial bin (def: 500)                       \n
  /// par[4]: minimum bin size in log(r)     (def: 0.1)                       \n
  /// file:   a C-type format string to form file name for writing profiles     
  ///                                                                           
  /// Usage of pointers: uses 'xcen' and 'vcen'\n                               
  /// Usage of flags:    uses in_subset()\n                                     
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class projprof : public manipulator {
  private:
    const unsigned   W;
    const double     L;
    vect             P;
    char*  const     FILE;
    mutable int      I;
    mutable output   OUT;
    mutable bool     FRST;
    void print_line(bool) const;
    //--------------------------------------------------------------------------
  public:
    const char*name    () const { return "projprof"; }
    const char*describe() const {
      return
	"measures projected radial profile w.r.t. 'xcen' and 'vcen' "
	"(default: origin) for 'subset' (default: all)";
    }
    //--------------------------------------------------------------------------
    fieldset need   () const { return fieldset::basic; }
    fieldset provide() const { return fieldset::empty; }
    fieldset change () const { return fieldset::empty; }
    //--------------------------------------------------------------------------
    projprof(const double*pars,
	     int          npar,
	     const char  *file) falcON_THROWING
    : 
      W    ( npar>3?  int(pars[3]): W_default ),
      L    ( npar>4?      pars[4] : L_default ),
      P    ( npar>2? vect(pars)   : vect(zero)),
      FILE ( (file && file[0])? falcON_NEW(char,strlen(file)+1) : 0 ),
      I    ( 0 ),
      FRST ( true )
    {
      if(debug(2) || file==0 || file[0]==0 || npar>5)
	std::cerr<<
	  " Manipulator \"projprof\" measures projected radial profile w.r.t.\n"
	  " 'xcen' and 'vcen' (default: origin) for bodies in 'subset'"
	  " (default: all):\n"
	  " sufrace density; mean, rotation & disp velocity; axis ratios are"
	  " computed as function of projected radius and written to file.\n"
	  " parameters:\n"
	  " (par[0],par[1],par[2]): position from which to project (to oo)\n"
	  " par[3]: minimum # bodies in radial bin (def: "<<W_default<<")\n"
	  " par[4]: minimum bin size in log(r)     (def: "<<L_default<<")\n"
	  " file:   format string to build file name\n";
      if(FILE) strcpy(FILE,file);
      if(file==0 || file[0]==0)
	falcON_THROW("Manipulator \"projprof\": no file given\n");
      if(npar>3 && pars[3]<0)
        falcON_THROW("Manipulator \"projprof\": W = %d < 0\n",W);
      if(L<0.)
	falcON_THROW("Manipulator \"projprof\": L = %f < 0\n",L);
    }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*) const;
    //--------------------------------------------------------------------------
    ~projprof() {
      if(FILE) falcON_DEL_A(FILE);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  inline void projprof::print_line(bool V) const {
    OUT  <<"#--------------------------------------------------";
    if(V)
      OUT<<"----------------------------------------------------";
    OUT  <<'\n';
  }
  //////////////////////////////////////////////////////////////////////////////
  bool projprof::manipulate(const snapshot*S) const
  {
    DebugInfo(2,"projprof::manipulate(): start\n");
    // make sure I refers to a non-existing file
    if(FRST) {
      if(std::strchr(FILE,'%'))
	while(output::file_exists(FILE,I)) ++I;
      FRST = false;
    }
    const vect*X0= S->pointer<vect>("xcen");
    const vect*V0= S->pointer<vect>("vcen");
    vect proj = P; if(X0) proj -= *X0;
    projected_profile PP(S,proj,W,L,X0,V0);
    if(OUT.reopen(FILE,I++,1)) {
      print_line(PP.has_vels());
      OUT  << "#\n"
	   << "# output from Manipulator \"projprof\"\n#\n";
      if(RunInfo::cmd_known ()) OUT<<"# command: \""<<RunInfo::cmd ()<<"\"\n";
      OUT  <<"# run at "<<RunInfo::time()<<'\n';
      if(RunInfo::user_known())
	OUT<<"#     by \""<<RunInfo::user()<<"\"\n";
      if(RunInfo::host_known())
	OUT<<"#     on \""<<RunInfo::host()<<"\"\n";
      if(RunInfo::pid_known())
	OUT<<"#     pid "<<RunInfo::pid()<<'\n';
      OUT  <<"#\n";
    }
    print_line(PP.has_vels());
    OUT  <<"#\n"
	 <<"# time    = "<<S->time()<<'\n';
    const std::string *SUB = S->pointer<std::string>("subset_description");
    if(SUB)
      OUT<<"# subset: "<< (*SUB) <<'\n';
    OUT  <<"# project = "<<proj;
    if(X0) OUT<<" = ("<<P<<") - xcen\n";
    else   OUT<<'\n';
    if(X0) OUT<<"# xcen    = "<<(*X0)<<'\n';
    if(V0) OUT<<"# vcen    = "<<(*V0)<<'\n';
    OUT  <<"#\n"
	 <<"#     radius        Sigma";
    if(PP.has_vels())
      OUT<<"        <v_l>      <v_rot>      sigma_l";
    OUT  <<"          b/a           PA";
    if(PP.has_vels())
      OUT<<"           RA";
    OUT  <<'\n';
    print_line(PP.has_vels());
    for(int i=0; i!=PP.N(); ++i) {
      OUT  << std::setw(12) << PP.rad(i) <<' '
	   << std::setw(12) << PP.Sig(i) <<' ';
      if(PP.has_vels())
	OUT<< std::setw(12) << PP.vlos(i)<<' '
	   << std::setw(12) << PP.vrot(i)<<' '
	   << std::setw(12) << PP.sigl(i)<<' ';
      OUT  << std::setw(12) << PP.bova(i)<<' '
	   << std::setw(12) << PP.posa(i);
      if(PP.has_vels())
	OUT<< ' '
	   << std::setw(12) << PP.rota(i);
      OUT  <<'\n';
    }
    OUT.flush();
    DebugInfo(2,"projprof::manipulate(): finished\n");
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} }

__DEF__MAN(falcON::Manipulate::projprof)
