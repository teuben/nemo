// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// lagrange_radii.cc                                                           |
//                                                                             |
// Copyright (C) 2002-2006 Walter Dehnen                                       |
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
// v 1.0    21/06/2002  WD created                                             |
// v 1.1    12/07/2002  WD added option masses                                 |
// v 1.2    20/08/2002  WD added options centering, alpha, and times           |
// v 1.3    29/08/2002  WD improved find_center()                              |
// v 1.4    30/08/2002  WD adapted this file for usage of MPI otherwise        |
// v 1.5    06/12/2002  WD debugged                                            |
// v 1.6    08/01/2003  WD abandoned centering, use 'center' instead           |
// v 2.0    13/03/2003  WD added stop mechanism; interpolate in R^2,           |
//                          => the actual radii may be slightly different      |
// v 2.1    20/03/2003  WD action reporting                                    |
// v 2.2    23/05/2003  WD automated NEMO history                              |
// v 2.3    29/10/2003  WD automated version etc; changed give: no default     |
// v 2.4    13/02/2004  WD debugged error with give=                           |
// v 2.5    24/03/2004  WD added 3% to the default masses list                 |
// v 3.0    07/05/2004  WD new sorting, 10 times faster than old version       |
// v 3.1    10/05/2004  WD fixed two bugs (allow for >1 @ same r, r=0)         |
// v 3.2    11/05/2004  WD made PUBLIC; changed name (origina: "lagrange_rad") |
// v 3.2.1  19/05/2004  WD change of give: if not given, we write all we got   |
// v 3.3    19/05/2004  WD re-written completely; sanity check for stdouts     |
// v 3.3.1  01/06/2004  WD new I/O; no tabfile output before 1st snapshot      |
// v 3.3.2  02/06/2004  WD create output files with first output, not earlier  |
//                         returns immediately if nothing will be done anyway  |
// v 3.3.3  20/05/2005  WD several minor updates                               |
// v 4.0    16/06/2005  WD new falcON                                          |
// v 4.1    13/06/2005  WD changes in fieldset                                 |
// v 4.1.1  27/06/2006  WD usage of bodyset                                    |
// v 4.1.2  27/06/2006  WD back to before bodyset                              |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "4.1.2"
#define falcON_VERSION_D "07-jul-2006 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile "lagrange_radii"
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
  "masses=0.01,0.03,0.1,0.3,0.5,0.7,0.9,0.97,0.99\n"
  "                   table of masses; will be parsed by nemoinp         ",
  "times=all\n        times to process                                   ",
  "tabfile=-\n        file to append table to                            ",
  "out=\n             output file [default: no output]                   ",
  "give=\n            output only these; if not given, output all we got ",
  "step=0\n           step between outputs                               ",
  "stopfile=\n        create stop file                                   ",
  "stopindex=0\n      index for Lagrange radius used in stop, first=0    ",
  "stopvalue=\n       stop if R_i < value ( or > |value| )               ",
  "stoprelative=f\n   use R_i / R_i[t=t_initial] in stop criterion       ",
  "stopafter=\n       don't stop before this time (default: t_initial)   ",
  "stopdelay=0\n      delay stopping after condition is satisfied        ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string
usage = "lagrange_radii -- find the Lagrange radii of a stellar system\n";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  // set parameters for Lagrange radii                                          
  const int      NMAX    = 100;
  double         M[NMAX] = {0.}, R[NMAX];
  int            N       = nemoinpd(getparam("masses"),M,NMAX);
  if(N>NMAX) falcON_THROW("too many mass shells specified");
  // set stop parameters                                                        
  bool           STOPPING     ( hasvalue   ("stopfile") );
  const int      STOPINDEX    ( getiparam  ("stopindex") );
  const real     STOPVALUE    ( getdparam_z("stopvalue") );
  const bool     STOPRELATIVE ( getbparam  ("stoprelative") );
  const double   STOPAFTER    ( getdparam_z("stopafter") );
  const double   STOPDELAY    ( getdparam  ("stopdelay") );
  if(STOPPING && (STOPINDEX < 0 || STOPINDEX >= N)) {
    warning("will not use stop condition, since stopindex not in [0,%d]",N-1);
    STOPPING = false;
  }
  if(STOPPING && STOPVALUE == zero) {
    warning("will not use stop condition, since stopvalue=0");
    STOPPING = false;
  }
  if(STOPPING && 
     STOPRELATIVE && (STOPVALUE>one || (STOPVALUE<zero && STOPVALUE>-one)))
    warning("stop condition already true at initial time");
  // deal with I/O                                                              
  const bool     DO_OUT(file_for_output("out"));
  const bool     DO_TAB(file_for_output("tabfile"));
  if(!DO_TAB && !DO_OUT && !STOPPING) {
    warning("no outputs wanted: nothing to be done; I'll finish\n");
    return;
  }
  const nemo_in  IN  (getparam("in"));
  const fieldset GIVE(DO_OUT? getioparam_a("give") : fieldset::o);
  const fieldset srcedata(fieldset::m|fieldset::x);
  const fieldset WANT(DO_TAB? (srcedata | GIVE) : GIVE);
  nemo_out       OUT;
  output         TAB;
  // loop snapshots & process them                                              
  snapshot       SHOT;
  fieldset       READ;
  bool           FIRST=true, STOPPED=false, HAS_RI0=false;
  double         TOUT, STEP(getdparam("step")), STOPTIME, LASTTIME, RI, RI0;
  while(IN.has_snapshot()) {
    // read time, read snapshot if in times; take first simulation time         
    const bool IN_TIMES = SHOT.read_nemo(IN,READ,WANT,getparam("times"),0);
    if(FIRST) {
      TOUT  = SHOT.time() - 1.e-10*STEP;
      FIRST = false;
    }
    if(!IN_TIMES)  continue;
    // find & write Lagrange radii                                              
    check_sufficient(READ,srcedata);
    find_lagrange_rad(&SHOT,N,M,R);
    // write table with Lagrange radii                                          
    if(DO_TAB) {
      if(!TAB.is_open()) {
	TAB.open(getparam("tabfile"),1);
	if(!TAB) falcON_THROW("cannot open tabfile\n");
	TAB << "#\n"
	    << "# \""<< (*(ask_history())) <<"\"\n"
	    << "#\n# time       ";
	for(int j=0; j!=N; ++j)
	  TAB << "r[" << std::setw(4) << 100*M[j] << "%] ";
	TAB << std::endl;
	TAB.stream().setf(std::ios::left, std::ios::adjustfield);
      }
      TAB << std::setw(12) << SHOT.time() <<' ';
      for(int i=0; i!=N; ++i)
	TAB << std::setw(8) << R[i] << ' ';
      TAB << std::endl;
    }
    // output of snapshot                                                       
    if(DO_OUT && SHOT.time() >= TOUT) {
      if(!OUT.is_open()) {
	OUT.open(getparam_z("out"));
	if(!OUT) falcON_THROW("cannot open out\n");
      }
      SHOT.write_nemo(OUT,READ&GIVE);
      TOUT += STEP;
    }
    // deal with STOPPING of simulation                                         
    if(STOPPING) {
      RI = R[STOPINDEX];
      if(!HAS_RI0) { RI0 = RI; HAS_RI0 = true; }
      if(!STOPPED) {
	register real STOPX = STOPRELATIVE? RI/RI0 : RI;
	if((STOPVALUE > zero && STOPX <     STOPVALUE  ) ||
	   (STOPVALUE < zero && STOPX > abs(STOPVALUE) ) ) {
	  STOPPED  = true;
	  STOPTIME = SHOT.time();
	  LASTTIME = max(SHOT.time()+STOPDELAY,STOPAFTER);
	}
      }
      if(STOPPED && SHOT.time() >= LASTTIME) {
	std::ofstream STOP(getparam("stopfile"));
	STOP<<" stopfile \""<<getparam("stopfile")<<"\"\n"
	    <<" generated by \""<< (*(ask_history())) <<"\"\n"
	    <<" on simulation time "<<SHOT.time()
	    <<", because the stop condition\n"
	    <<"\n   R(M="<<M[STOPINDEX]<<"*M_tot)";
	if(STOPRELATIVE) STOP<<"/R(M="<<M[STOPINDEX]<<"*M_tot, t=t_initial)";
	STOP<< ((STOPVALUE<zero)? " > " : " < ")
	    << abs(STOPVALUE) <<"\n\n"
	    <<" was satisfied at simulation time "<<STOPTIME<<"\n";
      }
    }
  }
}
