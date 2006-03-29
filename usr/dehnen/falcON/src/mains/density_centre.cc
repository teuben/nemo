// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// density_centre.cc                                                           |
//                                                                             |
// Copyright (C) 2002-2005 Walter Dehnen                                       |
//                                                                             |
// This program is free software; you can redistribute it and/or modify        |
// it under the terms of the GNU General Public License as published by        |
// the Free Software Foundation; either version 2 of the License, or (at       |
// your option) any later version.                                             |
//                                                                             |
// This program is distributed in the hope that it will be useful, but         |
// WITHOUT ANY WARRANTY; without even the implied warranty of                  |
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU           |
// General Public License for more details.                                    |
//                                                                             |
// You should have received a copy of the GNU General Public License           |
// along with this program; if not, write to the Free Software                 |
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                   |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 1.0   18/11/2002  WD created                                              |
// v 1.1   20/03/2003  WD                                                      |
// v 1.2   23/05/2003  WD automated NEMO history                               |
// v 1.2.1 23/10/2003  WD automated version, compiler in version etc           |
// v 1.3   12/02/2004  WD debugged, added central density to output            |
// v 1.4   13/02/2004  WD added option give                                    |
// v 1.5   11/03/2004  WD changed default fac=0.95 (was 0.8)                   |
// v 2.0   30/03/2004  WD new tool.h, new center.cc, new option XXX            |
// v 2.1   20/04/2004  WD abolish XXX, pre-center using tree                   |
// v 2.2   21/04/2004  WD abandon alpha & fac, use SPH type density estimate   |
// v 2.3   12/05/2004  WD made PUBLIC; renamed (original: "center")            |
// v 2.4   19/05/2004  WD change of give: if not given, we write all we got    |
//                        added option step to reduce snapshot outputs         |
// v 2.5   20/05/2004  WD improved iterative algorithm (cg)                    |
// v 2.5.1 26/05/2004  WD slight changes in iterative algorithm                |
// v 2.5.2 01/06/2004  WD new output; automated stdout control                 |
// v 2.5.3 02/06/2004  WD create output files with first output, not earlier   |
// v 2.5.4 20/05/2005  WD several minor updates                                |
// v 3.0   14/06/2005  WD new falcON                                           |
// v 3.1   13/06/2005  WD changes in fieldset                                  |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "3.1"
#define falcON_VERSION_D "13-jul-2005 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile "density_centre"
#endif
#define falcON_RepAction 0                         // no action reporting       
#define falcON_NOT_USING_PROPER                    // not using PROPRIETARY code
//-----------------------------------------------------------------------------+
#include <iostream>                                // C++ I/O                   
#include <fstream>                                 // C++ file I/O              
#include <iomanip>                                 // C++ I/O formatting        
#include <body.h>                                  // the bodies                
#include <public/io.h>                             // my NEMO I/O               
#include <public/tools.h>                          // my tools for bodies       
#include <main.h>                                  // main & NEMO stuff         
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n           input file                                         ",
  "out=\n             output file for snapshots                          ",
  "give=\n            output only these; if not given, output all we got ",
  "step=0\n           step between snapshot outputs                      ",
  "times=all\n        times to process                                   ",
  "Ncen=200\n         # bodies in center                                 ",
  "centrefile=-\n     file to write center position & velocity to        ",
  "centerfile=\n      superseeds centrefile                              ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string
usage = "density_centre -- find the global density maximum; center snapshots\n";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  // open I/O                                                                   
  const bool     DO_OUT   (file_for_output("out"));
  const bool     DO_CENOUT(hasvalue("centerfile")?
			   file_for_output("centerfile") :
			   file_for_output("centrefile") );
  if(!DO_OUT && !DO_CENOUT) {
    ::warning("no outputs wanted: nothing to be done; I'll finish\n");
    return;
  }
  const nemo_in  IN   (getparam("in"));
  const fieldset GIVE (DO_OUT   ? getioparam_a("give") : fieldset::o);
  const fieldset WANT (DO_CENOUT? fieldset(fieldset::basic | GIVE) : GIVE);
  nemo_out       OUT;
  output         CENOUT;
  // loop snapshots & process them                                              
  snapshot       SHOT;
  int            NOLD = 0u;
  OctTree*       TREE = 0;
  fieldset       READ;
  real           RHO,RAD(zero);
  vect           XCEN(zero),VCEN(zero);
  bool           FIRST = true;
  double         TOUT, STEP(getdparam("step"));
  const unsigned NCEN  (getiparam("Ncen"));
  while(IN.has_snapshot()) {
    // read time, read snapshot if in times                                     
    const bool IN_TIMES = SHOT.read_nemo(IN,READ,WANT,getparam("times"),0);
    if(FIRST) {
      TOUT  = SHOT.time() - 1.e-10*STEP;
      FIRST = false;
    }
    if(!IN_TIMES) continue;
    // build tree & find first estimate for density centre                      
    check_sufficient(READ,fieldset(fieldset::m|fieldset::x));
    if(TREE==0 || NOLD != SHOT.N_bodies()) {
      if(TREE) delete TREE;
      TREE = new OctTree(&SHOT,NCEN/4);
    } else {
      TREE->build(NCEN/4);
    }
    estimate_density_peak(TREE,0u,NCEN,XCEN,RAD);
    RAD *= 3;
    // iterate to improve density centre                                        
    if(READ.contain(fieldset::v))
      find_centre(&SHOT,NCEN,XCEN,RAD,&VCEN,&RHO);
    else {
      find_centre(&SHOT,NCEN,XCEN,RAD,0,&RHO);
      VCEN = zero;
    }
    // output of density centre: position & velocity                            
    if(DO_CENOUT) {
      if(!CENOUT.is_open()) {
	CENOUT.open(hasvalue("centerfile")?
		    getparam("centerfile") :
		    getparam("centrefile"),1);
	if(!CENOUT)
	  ::error("cannot open %s\n",
		  hasvalue("centerfile")? "centerfile" : "centrefile");
	CENOUT << "#\n"
	       << "# \""<< (*(ask_history())) <<"\"\n"
	       << "#\n# "
	       << "           time  "
	       << "              x               y               z  "
	       << "             vx              vy              vz  "
	       << "        density\n";
      }
      CENOUT <<"  "
	     << std::setw(15) << std::setprecision(8) << SHOT.time() << "  "
	     << std::setw(15) << std::setprecision(8) << XCEN[0]     << ' '
	     << std::setw(15) << std::setprecision(8) << XCEN[1]     << ' '
	     << std::setw(15) << std::setprecision(8) << XCEN[2]     << "  "
	     << std::setw(15) << std::setprecision(8) << VCEN[0]     << ' '
	     << std::setw(15) << std::setprecision(8) << VCEN[1]     << ' '
	     << std::setw(15) << std::setprecision(8) << VCEN[2]     << "  "
	     << std::setw(15) << std::setprecision(8) << RHO     << std::endl;
    }
    // output of snapshot                                                       
    if(DO_OUT && SHOT.time() >= TOUT) {
      if(!OUT.is_open()) {
	OUT.open(getparam_z("out"));
	if(!OUT) ::error("cannot open out\n");
      }
      LoopAllBodies(&SHOT,Bi) {
	Bi.pos() -= XCEN;
	Bi.vel() -= VCEN;
      }
      SHOT.write_nemo(OUT,READ&GIVE);
      TOUT += STEP;
    }
  }
  if(TREE) delete TREE;
}
