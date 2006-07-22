// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   src/manip/diskanalysis.cc                                           
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2004-2006                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2004-2006 Walter Dehnen                                        
//                                                                              
// This is a non-public part of the code.                                       
// It is property of its author and not to be made public without his written   
// consent.                                                                     
//                                                                              
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// history:                                                                     
//                                                                              
// v 0.0    29/06/2006  WD created drawing from DiskAnalysis.cc                 
// v 0.1    07/07/2006  WD replacing bodyset with flags::ignore                 
////////////////////////////////////////////////////////////////////////////////
#include <public/defman.h>
#include <proper/diskanalysis.h>
#include <public/io.h>

namespace falcON { namespace Manipulate {
  //////////////////////////////////////////////////////////////////////////////
  const int NRAD_default = 20;
  const int MMAX_default = 8;
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class diskanalysis                                                         
  //                                                                            
  /// manipulator: performs azimuthal mode analysis of bodies in_subset()       
  ///                                                                           
  /// This manipulator computes the amplitudes and phases of the m=0,...,m_max  
  /// azimuthal modes in n_rad cylindrical annuli for all bodies in_subset      
  /// (default: all, see set_subset) and writes the result to file in XDR       
  /// format (which can be read by, say, PrintAnalysis). The centre is taken    
  /// from 'xcen' and 'vcen' (default: origin).                                 
  ///                                                                           
  /// Meaning of the parameters:\n                                              
  /// par[0]: number of radial bins (default: 20)\n                             
  /// par[1]: maximal azimuthal wave number m (default: 8)\n                    
  /// file  : name for output file                                              
  ///                                                                           
  /// Usage of pointers: uses 'xcen'\n                                          
  /// Usage of flags:    uses in_subset()\n                                     
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class diskanalysis : public manipulator {
  private:
    const   unsigned       NR, MM;
    mutable DiskComponents DC;
    mutable ComponentsO    CO;
    mutable DiskAnalysis   DA;
    //--------------------------------------------------------------------------
  public:
    const char*name    () const { return "diskanalysis"; }
    const char*describe() const {
      return message("performs modal disk analysis up to m=%d "
		     "using %d radial bins "
		     "w.r.t. 'xcen' and 'vcen' (default: origin) "
		     "for 'subset' (default: all) of bodies", MM,NR);
    }
    //--------------------------------------------------------------------------
    fieldset need   () const { return fieldset::basic; }
    fieldset provide() const { return fieldset::o; }
    fieldset change () const { return fieldset::o; }
    //--------------------------------------------------------------------------
    diskanalysis(const double*pars,
		 int          npar,
		 const char  *file) falcON_THROWING
    : NR ( npar>0? int(pars[0]) : NRAD_default ),
      MM ( npar>1? int(pars[1]) : MMAX_default ),
      DC ( NR, MM ),
      CO ( file, NR, MM )
    { if(debug(2) || file==0 || file[0]==0 || npar>2)
	std::cerr<<
	  " Manipulator \""<<name()<<"\" measures disk profile w.r.t. 'xcen'"
	  " and 'vcen' (default: origin) for bodies in_subset() (def: all);\n"
	  " results are written in XDR format to file.\n"
	  " parameters:\n"
	  " par[0]: number of radial bins (def: "<<NRAD_default<<")\n"
	  " par[1]: maximal azimuthal wave number m (def: "<<MMAX_default<<")\n"
	  " file:   file name\n";
      if(!CO.is_open())
	falcON_THROW("Manipulator \"diskanalysis\": cannot open output file\n");
    }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*) const;
    //--------------------------------------------------------------------------
    ~diskanalysis() {}
  };
  //////////////////////////////////////////////////////////////////////////////
  bool diskanalysis::manipulate(const snapshot*S) const
  {
    const vect *X0 = S->pointer<vect>("xcen");
    const vect *V0 = S->pointer<vect>("vcen");
    DA.analyse(S,DC,X0,V0);
    CO.write(DC);
    return false;
  }
  //////////////////////////////////////////////////////////////////////////////
} }
__DEF__MAN(falcON::Manipulate::diskanalysis)
