// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// myradial_profile.cc                                                         |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2005                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// This is a non-public part of the code.                                      |
// It is property of its author and not to be made public without his written  |
// consent.                                                                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 1.0    27/08/2002 WD created                                              |
// v 1.1    28/08/2002 WD improved center finder                               |
// v 1.2    30/08/2002 WD adapted this file for usage of MPI otherwise         |
// v 1.3    13/12/2002 WD abandoned centering, bound/unbount                   |
// v 1.4    20/03/2003 WD                                                      |
// v 1.5    23/05/2003 WD automated NEMO history                               |
// v 1.6    07/11/2003 WD automated NEMO version & compile time                |
// v 1.7    27/05/2004 WD significant improvements in smod.h/cc                |
// v 1.8    02/06/2004 WD this file largely re-written, new options out & give |
// v 1.9    04/06/2004 WD renamed (previously: dens_prof), new smod.h          |
// v 2.0    24/06/2004 WD new falcON                                           |
// v 2.1    13/06/2005 WD changes in fieldset                                  |
// v 3.0      /11/2005 PJM added skew and kurtosis, created new version        |
//                                                                             |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "2.1"
#define falcON_VERSION_D "13-jul-2005 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                   // this is a NEMO program 
#  error You need NEMO to compile "radial_profile"
#endif
#define falcON_RepAction 0                            // no action reporting    
//-----------------------------------------------------------------------------+
#include <iostream>                                   // C++ I/O                
#include <fstream>                                    // C++ file I/O           
#include <iomanip>
#include <body.h>                                     // the bodies             
#include <pjm/myprofile.h>                            // profiling utilities    
#include <public/io.h>                                // my NEMO I/O            
#include <utils/inline_io.h>                          // my I/O utilities       
#include <main.h>                                     // main & NEMO stuff      
using namespace falcON;
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n           input file                                         ",
  "tabfile=-\n        file for output (may contain format string)        ",
  "times=all\n        times to process                                   ",
  "Lmax=0.1\n         width in log(r) of smoothing windows               ",
  "Nmin=100\n         mininum # bodies in smoothing window               ",
  "out=\n             output file [default: no output]                   ",
  "step=0\n           step between outputs                               ",
  "give=\n            output only these; if not given, output all we got ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string
usage = "myradial_profile -- writes table with spherically averaged quantities\n";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  // set up parameters, check for need to do outputs                            
  const bool    DO_SNP(file_for_output("out"));
  const bool    DO_TAB(file_for_output("tabfile"));
  if(!DO_SNP && !DO_TAB) {
    ::warning("no outputs wanted: nothing to be done; I'll finish\n");
    return;
  }
  const nemo_in  IN (getparam("in"));
  const fieldset GIVE(DO_SNP? getioparam_a("give") : fieldset(fieldset::o));
  const fieldset WANT(DO_TAB? fieldset(fieldset::basic | GIVE) : GIVE);
  nemo_out       SNP;
  output         TAB;
  fieldset       READ;
  snapshot       SHOT;
  double         TSNP;
  const double   STEP(getdparam("step"));
  bool           FIRST=true;
  int            index=0;
  // loop snapshots, process them and make outputs                              
  while(IN.has_snapshot()) {
    // read time, read snapshot if in times; take first simulation time         
    const bool IN_TIMES = SHOT.read_nemo(IN,READ,WANT,getparam("times"),0);
    if(FIRST) {
      TSNP  = SHOT.time(); - 1.e-10*STEP;
      FIRST = false;
    }
    if(!IN_TIMES)  continue;
    // process to find radial profiles, write them to tabfile                   
    if(DO_TAB) {
      if(TAB.reopen(getparam("tabfile"),index++,1)) {
	if(!TAB) ::error("cannot open tabfile\n");
	TAB << "#\n"
	    << "# \""<< (*(ask_history())) <<"\"\n"
	    << "#\n";
      }
      check_sufficient(READ,fieldset(fieldset::m|fieldset::x));
      spherical_profile SP(&SHOT,
			   getiparam("Nmin"),
			   getdparam("Lmax"),
			   READ.contain(fieldset::v));
      TAB << "#\n"
	  << "# time = "<<SHOT.time()<<": "<<SHOT.N_bodies()<<" bodies\n"
	  << "#     radius          rho        vcirc";
      if(SP.has_vels())
	TAB<<"        <v_r>    <v_theta>      <v_phi>"
	   <<"      sigma_r     sigma_th    sigma_phi"
	   <<"       skew_r      skew_th     skew_phi"
	   <<"       kurt_r      kurt_th     kurt_phi";
      TAB  <<"\n#\n";
      for(int i=0; i!=SP.N(); ++i) {
	TAB  << std::setw(12) << SP.rad(i) <<' '
	     << std::setw(12) << SP.rho(i) <<' '
	     << std::setw(12) << sqrt(SP.vcq(i));
	if(SP.has_vels())
	  TAB<< ' '
	     << std::setw(12) << SP.vrad(i) <<' '
	     << std::setw(12) << SP.vthe(i) <<' '
	     << std::setw(12) << SP.vphi(i) <<' '
	     << std::setw(12) << SP.sigr(i) <<' '
	     << std::setw(12) << SP.sigt(i) <<' '
	     << std::setw(12) << SP.sigp(i) <<' '
	     << std::setw(12) << SP.sker(i) <<' '
	     << std::setw(12) << SP.sket(i) <<' '
	     << std::setw(12) << SP.skep(i) <<' '
	     << std::setw(12) << SP.kurr(i) <<' '
	     << std::setw(12) << SP.kurt(i) <<' '
	     << std::setw(12) << SP.kurp(i);
	TAB  << '\n';
      }
      TAB.flush();
    }
    // output of snapshot                                                       
    if(DO_SNP && SHOT.time() >= TSNP) {
      if(!SNP.is_open()) {
	SNP.open(getparam_z("out"));
	if(!SNP) falcON_THROW("cannot open out\n");
      }
      SHOT.write_nemo(SNP,READ&GIVE);
      TSNP += STEP;
    }
  }
}
//------------------------------------------------------------------------------
