// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/public/manip/sphereprof.cc                                      
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
// v 0.0    31/05/2006  WD created                                              
// v 1.0    27/06/2006  WD using bodyset in 'subset' instead of pars[1,2]       
// v 1.0.1  03/07/2006  WD renamed: radprof -> sphereprof                       
// v 1.1    04/07/2006  WD doxygen documented and made public                   
// v 1.2    07/07/2006  WD using flags::ignore (in_subset()) instead of subset  
// v 1.3    09/08/2006  WD warn if non-spherical                                
// v 1.3.1  06/11/2006  WD change in profile.h                                  
// v 1.3.2  13/05/2008  WD debugged error with minor & major axes               
// v 1.3.2  11/06/2008  WD new DebugInfo                                        
////////////////////////////////////////////////////////////////////////////////
#include <public/defman.h>
#include <public/profile.h>
#include <public/io.h>

namespace falcON { namespace Manipulate {
  //////////////////////////////////////////////////////////////////////////////
  const char*Line ="----------------------------------------";
  //////////////////////////////////////////////////////////////////////////////
  class PrintSmall {
    int N,P;
  public:
    PrintSmall(int n) : N(n), P(1) {
      for(int i=0; i!=n; ++i) P *= 10;
    }
    template<typename X>
    std::ostream&print_pos(std::ostream&o, X const&x) const {
      int I = int(P*x+0.5);
      if(I >= P) return o << "1." << std::setw(N-1) << std::setfill('0') << 0
			  << std::setfill(' ');
      else       return o << '.' << std::setw(N) << std::setfill('0') << I
			  << std::setfill(' ');
    }
    template<typename X>
    std::ostream&print(std::ostream&o, X const&x) const {
      if(x < 0) { o << '-'; return print_pos(o,-x); }
      else      { o << ' '; return print_pos(o, x); }
    }
    template<typename X>
    std::ostream&print_dir(std::ostream&o, tupel<3,X> const&x) const {
      return print(print(print(o,x[0])<<' ',x[1])<<' ',x[2]);
    }
    template<typename X>
    std::ostream&print_dir(std::ostream&o, const X x[3]) const {
      return print(print(print(o,x[0])<<' ',x[1])<<' ',x[2]);
    }
    const char*line_pos() const { return Line+(39-N); }
    const char*line    () const { return Line+(38-N); }
    const char*line_dir() const { return Line+(38-3*(N+2)); }
  };
  //////////////////////////////////////////////////////////////////////////////
  const int    W_default = 500;
  const double L_default = 0.1;
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class sphereprof                                                           
  //                                                                            
  /// manipulator: estimated radial spherical profiles for bodies in_subset()   
  ///                                                                           
  /// This manipulator estimates the radial profiles of density, velocity,      
  /// velocity dispersion, axis ratios, and orientations from spherical binning 
  /// of all bodies in_subset() (default: all, see set_subset) w.r.t. position  
  /// 'xcen' and velocity 'vcen'.                                             \n
  /// See falcON::spherical_profile for details.                                
  /// \note                                                                     
  /// Don't use this manipulator for analysis of non-spherical components as    
  /// then the shape information will be biased and the density profile may     
  /// also be affected. A warning will be issued if c/a<0.8 at any radius.      
  ///                                                                           
  /// Meaning of the parameters:\n                                              
  /// par[0]: minimum # bodies in radial bin (def: 500)                       \n
  /// par[1]: minimum bin size in log(r)     (def: 0.1)                       \n
  /// file:   a C-type format string to form file name for writing profiles     
  ///                                                                           
  /// Usage of pointers: uses 'xcen' and 'vcen'\n                               
  /// Usage of flags:    uses in_subset()\n                                     
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class sphereprof : public manipulator {
  private:
    const int        W;
    const double     L;
    char*  const     FILE;
    mutable int      I;
    mutable output   OUT;
    const PrintSmall PS;
    //--------------------------------------------------------------------------
    void print_line(bool) const;
  public:
    const char*name    () const { return "sphereprof"; }
    const char*describe() const {
      return
	"measures, by spherical averaging, the radial profile w.r.t. 'xcen' "
	"and 'vcen' (default: origin) for bodies in_subset() (default: all)";
    }
    //--------------------------------------------------------------------------
    fieldset need   () const { return fieldset::basic; }
    fieldset provide() const { return fieldset::o; }
    fieldset change () const { return fieldset::o; }
    //--------------------------------------------------------------------------
    sphereprof(const double*pars,
	       int          npar,
	       const char  *file) falcON_THROWING
    : W    ( npar>0?     int(pars[0])    : W_default ),
      L    ( npar>1?         pars[1]     : L_default ),
      I    ( 0 ),
      FILE ( (file && file[0])? falcON_NEW(char,strlen(file)+1) : 0 ),
      PS   ( 3 )
    {
      if(debug(2) || file==0 || file[0]==0 || npar>2)
	std::cerr<<
	  " Manipulator \""<<name()<<"\" measures radial profile w.r.t."
	  " 'xcen' and\n"
	  " 'vcen' (default: origin) for bodies in 'subset' (default: all):\n"
	  " density; circ, mean & disp velocity; axis ratios are computed as\n"
	  " function of radius and written to file.\n"
	  " parameters:\n"
	  " par[0]: minimum # bodies in radial bin (def: "<<W_default<<")\n"
	  " par[1]: minimum bin size in log(r)     (def: "<<L_default<<")\n"
	  " file:   format string to build file name\n";
      if(FILE) strcpy(FILE,file);
      if(file==0 || file[0]==0)
	falcON_THROW("Manipulator \"%s\": no file given\n",name());
      if(W <0) falcON_THROW("Manipulator \"%s\": W = %d < 0\n",name(),W);
      if(L<0.) falcON_THROW("Manipulator \"%s\": L = %f < 0\n",name(),L);
    }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*) const;
    //--------------------------------------------------------------------------
    ~sphereprof() {
      if(FILE) falcON_DEL_A(FILE);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  inline void sphereprof::print_line(bool vels) const
  {
    OUT  <<"#-------------------------------------";
    if(vels)
      OUT<<"--------------------------"
	 <<"---------------------------------------";
    OUT  <<"----------------------------" << PS.line_dir() << PS.line_dir();
    if(vels)
      OUT<< "--" << PS.line_dir();
    OUT  <<'\n';
  }
  //////////////////////////////////////////////////////////////////////////////
  bool sphereprof::manipulate(const snapshot*S) const
  {
    DebugInfo(2,"sphereprof::manipulate(): start\n");
    const vect*X0= S->pointer<vect>("xcen");
    const vect*V0= S->pointer<vect>("vcen");
    spherical_profile SP(S,W,L,X0,V0);
    if(OUT.reopen(FILE,I++,1)) {
      print_line(SP.has_vels());
      OUT << "#\n"
	  << "# output from Manipulator \"sphereprof\"\n#\n";
      if(RunInfo::cmd_known ())
	OUT<<"# command: \""<<RunInfo::cmd ()<<"\"\n";
      OUT  <<"# run at "<<RunInfo::time()<<'\n';
      if(RunInfo::user_known())
	OUT<<"#     by \""<<RunInfo::user()<<"\"\n";
      if(RunInfo::host_known())
	OUT<<"#     on \""<<RunInfo::host()<<"\"\n";
      if(RunInfo::pid_known())
	OUT<<"#     pid "<<RunInfo::pid()<<'\n';
      OUT <<"#\n";
    }
    print_line(SP.has_vels());
    OUT  <<"#\n"
	 <<"# time   = "<<S->time()<<'\n';
    const std::string *SUB = S->pointer<std::string>("subset_description");
    if(SUB)
      OUT<<"# subset: "<< (*SUB) <<'\n';
    if(X0) OUT<<"# xcen   = "<<(*X0)<<'\n';
    if(V0) OUT<<"# vcen   = "<<(*V0)<<'\n';
    OUT  <<"#\n"
	 <<"#   radius        rho      vcirc";
    if(SP.has_vels())
      OUT<<"      <v_r>    <v_phi>"
	 <<"    sigma_r   sigma_th  sigma_phi";
    OUT  <<"     c/a      b/a"
	 <<"        major axis         minor axis";
    if(SP.has_vels())
      OUT<<"      rotation axis";
    OUT  <<'\n';
    print_line(SP.has_vels());
    bool nonspherical = false;
    for(int i=0; i!=SP.N(); ++i) {
      OUT  << print(SP.rad(i),10,4) << ' '
	   << print(SP.rho(i),10,4) << ' '
	   << print(sqrt(SP.vcq(i)),10,4) << ' ';
      if(SP.has_vels())
	OUT<< print(SP.vrad(i),10,4) << ' '
	   << print(SP.vphi(i),10,4) << ' '
	   << print(SP.sigr(i),10,4) << ' '
	   << print(SP.sigt(i),10,4) << ' '
	   << print(SP.sigp(i),10,4) << ' ';
      OUT  << print(SP.cova(i), 8,3) << ' '
	   << print(SP.bova(i), 8,3) << ' ';
      PS.print_dir(OUT, SP.major_axis(i)) << "  ";
      PS.print_dir(OUT, SP.minor_axis(i));
      if(SP.has_vels()) {
	OUT<< "  ";
	PS.print_dir(OUT, SP.drot(i));
      }
      OUT<<'\n';
      if(SP.cova(i) < 0.8) nonspherical = true;
    }
    OUT.flush();
    if(nonspherical)
      warning("Manipulator sphereprof: "
	      "shape seems significantly non-spherical at t=%f\n",S->time());
    DebugInfo(2,"sphereprof::manipulate(): finished\n");
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} }
////////////////////////////////////////////////////////////////////////////////
__DEF__MAN(falcON::Manipulate::sphereprof)
