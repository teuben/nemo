// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// SerialMBodyCode.cc                                                          |
//                                                                             |
// Copyright Walter Dehnen, 2006                                               |
//                                                                             |
// This is a non-public part of the code.                                      |
// It is property of its author and not to be made public without his written  |
// consent.                                                                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 1.0    17/03/2006 WD  created, modelled on gyrfalcON and FewBodyCode      |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "1.0"
#define falcON_VERSION_D "17-mar-2006 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need "NEMO" to compile FewBodyCode
#endif
#include <proper/mbody.h>                          // potential expansion       
#include <public/manip.h>                          // run-time manipulator      
#include <main.h>                                  // main & NEMO stuff         
#include <iomanip>                                 // C++ formating             
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n           input file                                         ",
  "out=???\n          snapshout output file                              ",
  "tstop=\n           final integration time [default: never]            ",
  "step=1\n           time between primary outputs; 0 -> every step      ",
  "logfile=-\n        file for log output                                ",
  "logstep=1\n        # blocksteps between log outputs                   ",
  "hmax=4\n           tau_max = (1/2)^hmax                               ",
  "fac=???\n          time step control parameter                        ",
  "give=mxv\n         list of output specifications. Recognizing:\n"
  "                    m: mass                              (default)\n"
  "                    x: position                          (default)\n"
  "                    v: velocity                          (default)\n"
  "                    a: acceleration\n"
  "                    p: N-body potential                               ",
  "Grav=4*Pi**2\n     Newton's constant of gravity (must be >0)          ",
  "manipname=\n       name of run-time manipulator                       ",
  "manippars=\n       parameters for manipulator                         ",
  "manipfile=\n       data file required by manipulator                  ",
  "manippath=\n       path to search for manipulator                     ",
  "manipinit=t\n      manipulate initial snapshot?                       ",
  "flevel=.\n         file for ascii output of N_level                   ",
//   "fout=.\n           file for ascii output every time step              ",
//   "prec=8\n           precision used                                     ",
  "startout=t\n       primary output for t=tstart?                       ",
  "lastout=t\n        primary output for t=tstop?                        ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string usage = "SerialMBodyCode -- direct multi-body code";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  // 1. set some parameters                                                     
  const bool
    never_ending = !hasvalue("tstop"),                // integrate forever?     
    lastout      = getbparam("lastout");              // write out last snapshot
  const double
    t_end   = getdparam_z("tstop"),                   // integrate until t=t_end
    dt_out  = getdparam("step");                      // primary output interval
  const int
    logstep = getiparam("logstep");                   // # blocksteps/logoutput 
  const fieldset
    write   = getioparam("give");                     // what to output 0?      
  const Manipulator MANIP(getparam_z("manipname"),    // IF(manipname given)    
			  getparam_z("manippars"),    //   THEN                 
			  getparam_z("manipfile"),    //   initialize N-body    
			  getparam_z("manippath"));   //   manipulator          
  output LEVOUT(getparam("flevel"));                  // output for N_level     
  // 2. initialize N-body integrator etc                                        
  SerialMBodyCode MBDY(getparam("in"),
		       getiparam("hmax"),
		       getrparam("fac"),
		       getrparam("Grav"));
  if(!never_ending && t_end < MBDY.initial_time()) {  // IF(t_end < t_start)    
    warning("tstop < t_ini: nothing to be done\n");   //   THEN we are done     
    return;                                           //   and can stop here    
  }                                                   // ENDIF                  
  if(MANIP) {                                         // IF manipulating        
    const_cast<snapshot*>(MBDY.my_snapshot())->       //   supported required   
      add_fields(MANIP.provide());                    //   data                 
    if(getbparam("manipinit"))                        //   IF(manipulating 1st) 
      MANIP(MBDY.my_snapshot());                      //     do it now          
  }                                                   // ENDIF                  
  // 3. open output streams & make initial outputs                              
  output LOGOUT(getparam("logfile"));                 // log output stream      
  nemo_out OUT(hasvalue("out")?                       // open NEMO output       
	       getparam("out") : getparam("in"),      //   if not out given,    
	       !hasvalue("out"));                     //   append to input      
  bool written=false;                                 // current snapshot wrtten
  if(getbparam("startout")) {                         // IF(start output)       
    MBDY.write(OUT,write);                            //   write snapshot       
    written = true;                                   //   record writing       
  }                                                   // ENDIF                  
  if(LOGOUT) {                                        // IF any log output      
    MBDY.describe  (LOGOUT);                          //   put history to logout
    MBDY.stats_head(LOGOUT);                          //   header for logout    
    if(!LOGOUT.is_appending())                        //   Unless appending     
      MBDY.stats   (LOGOUT);                          //     statistics ->logout
  }                                                   // ENDIF                  
  if(LEVOUT) MBDY.plot_Nlev(LEVOUT);                  // output N_level         
  // 4. time integration & outputs                                              
  double t_out = MBDY.initial_time()+0.999999*dt_out; // time for next output 0 
  for(int steps=1;                                    // blockstep counter      
      (never_ending || MBDY.time() < t_end);          // WHILE t < t_end        
      ++steps) {                                      //   increment counter    
    MBDY.full_step();                                 //   make full block step 
    if(LEVOUT) MBDY.plot_Nlev(LEVOUT);                //   output N_level       
    if(LOGOUT && steps%logstep ==0)                   //   IF(time for logout)  
      MBDY.stats(LOGOUT);                             //     statistics output  
    if(MANIP)                                         //   IF(manipulating)     
      if(MANIP(MBDY.my_snapshot()))                   //     IF(manip says so)  
	break;                                        //       STOP simulation  
    if(OUT && MBDY.time() >= t_out) {                 //   IF(t >= t_out)       
      MBDY.write(OUT,write);                          //     primary output     
      t_out  += dt_out;                               //     increment t_out    
      written = true;                                 //     written out        
    } else                                            //   ELSE                 
      written = false;                                //     not written        
  }                                                   // END: WHILE             
  if(OUT && !written && lastout)                      // IF not yet done:       
    MBDY.write(OUT,write);                            //   write last snapshot  
}
////////////////////////////////////////////////////////////////////////////////
