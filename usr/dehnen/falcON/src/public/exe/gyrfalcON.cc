// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// gyrfalcON.cc                                                                |
//                                                                             |
// Copyright (C) 2000-2008 Walter Dehnen                                       |
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
// v 1.0.0  28/05/2001  WD created with the help of PJT, based on YANC         |
// v 1.0.1  01/06/2001  WD bug in this file removed, added history output      |
// v 1.0.2  13/06/2001  WD Ncrit option added (before Ncrit==1)                |
// v 1.0.3  20/06/2001  WD tree::make() back to adding-leafs algorithm         |
//                         Coordinatesystem given with snapshot                |
// v 1.0.4  25/06/2001  WD bug (builder::link_tree()) removed (thanks to PJT)  |
// v 1.0.5  28/06/2001  WD changes in tree::grow(). Avoid box overflow.        |
//                         bug (in builder) removed                            |
// v 1.0.6  01/08/2001  WD substantial code upgrade in the tree (falcON)       |
// v 1.0.7  02/08/2001  WD Nreuse option added (before Nreuse==0)              |
// v 1.0.8  08/08/2001  WD new data lay-out for bodies, enabling special       |
//                         treatment for leap-frog, saves memory & time (?)    |
// v 1.0.9  17/09/2001  WD added support for dynamically linked external       |
//                         potential (potential.h)                             |
// v 1.0.10 20/09.2001  WD added option ("resume") for resuming an old or      |
//                         interrupted simulation. Appends to input file       |
// v 1.0.11 21/09/2001  WD improved handling of nemo I/O; can now also read    |
//                         PosTag & VelTag instead of PhaseSpaceTag            |
// v 1.0.12 21/09/2001  WD added parameter give_acc                            |
// v 1.0.13 11/10/2001  WD changed give_pot option to allow only N-body pot    |
//                         added option startout                               |
// v 1.0.14 18/10/2001  WD interweaving interaction & evaluation phase of      |
//                         gravity approximation: saves about 20b/body         |
// v 1.0.15 23/10/2001  WD some changes in I/O handling (body, yanc)           |
// v 1.0.16 23/11/2001  WD added option logout (to allow output pipe)          |
// v 1.0.17 08/01/2002  WD changed body lay-out (gives 2% speed-up)            |
// v 1.1.0  21/01/2002  WD minor changes                                       |
// v 1.1.1  25/01/2002  WD minor changes                                       |
// v 1.1.2  29/01/2002  WD minor changes (1.5% speed-up)                       |
// v 1.1.3  26/02/2002  WD time-symmetric stepping criterion tau=fac/acc       |
// v 1.1.4  27/02/2002  WD added options: out2, step2                          |
// v 1.1.5  01/03/2002  WD added time stepping criterion tau=fp/pot            |
// v 1.1.6  05/03/2002  WD added time stepping criterion tau=fc*sqrt(pot)/acc  |
// v 1.2    07/06/2002  WD replaced giveacc,givepot,giverho with give, give2   |
// v 1.2.1  11/06/2002  WD support for kernels Fn and Kn withdrawn             |
// v 1.2.2  14/06/2002  WD added output option for flag & level                |
// v 1.2.3  17/06/2002  WD allow for individual adaptive softening lengths     |
//                         changed option Nreuse -> hgrow                      |
// v 1.2.4  26/08/2002  WD bug in MAC removed,                                 |
//                         tree re-build (saves 40% CPU time on build)         |
// v 1.2.5  28/08/2002  WD bug with hgrow removed                              |
// v 1.2.6  29/08/2002  WD improved initialization of external potential       |
// v 1.2.7  30/08/2002  WD adapted this file for usage of MPI otherwise        |
// v 1.2.8  09/09/2002  WD further adaption to MPI usage                       |
// v 1.2.9  05/11/2002  WD added option Grav (for comparison with GADGET)      |
// v 1.3.0  15/11/2002  WD various updates (SSE code, eps_i treatment)         |
// v 1.4    20/11/2002  WD splitted between public and proprietary code        |
// v 1.4.1  26/11/2002  WD debugged some features                              |
// v 1.5    04/12/2002  WD re-named "gyrfalcON" (previously "YancNemo")        |
//                         several changes in file layout for public version   |
// v 1.5.1  09/01/2003  WD saver C-macros, default parameters                  |
// v 1.5.2  13/01/2003  WD never ending; logstep; stopfile                     |
// v 1.5.3  24/01/2003  WD option lastout added.                               |
// v 1.5.4  03/03/2003  WD debugged: handling of external pot in total energy  |
// v 1.5.5  13/03/2003  WD added check for over-using of stdout/pipe           |
// v 1.5.6  17/03/2003  WD options emin, fea, limits on |eps_new/eps_old|      |
// v 1.5.7  20/03/2003  WD changes in gravity, action reporting (proper only)  |
// v 1.6    02/06/2003  WD allow for individual but fixed eps_i by eps<0       |
// v 1.6.1  28/07/2003  WD happy gcc 3.3 (about 6% faster than gcc 3.2)        |
// v 1.6.2  08/08/2003  WD version provides compiler info, automated           |
// v 1.6.3  13/08/2003  WD fixed bug with individual_fixed; thanks to J.Bailin |
// v 1.7.0  05/09/2003  WD individual eps made public; changes in tensors      |
// v 1.7.1  17/09/2003  WD changes in tree: avoiding template specs            |
// v 1.7.2  07/10/2003  WD changes in gravity: using common basic_tree         |
// v 1.7.3  23/10/2003  WD changes in design of gravity, tree, kernel, falcON  |
// v 1.8    05/11/2003  WD changes in grav; changed io::P to io::q             |
// v 1.8.1  10/02/2004  WD minor change in this file (logfile)                 |
// v 1.8.2  11/02/2004  WD allow for Grav=0 (Grav now handled in gravity.h)    |
// v 1.8.3  18/02/2004  WD bug with Grav=0 fixed; avoid tree building if G=0   |
// v 1.9    19/02/2004  WD added option root_center                            |
// v 1.9.1  23/02/2004  WD improved diagnose output (new T, V_in, [V_ex,] W)   |
// v 1.9.2  27/02/2004  WD use nemo::error() instead of falcON::error()        |
// v 2.0    11/03/2004  WD elimated yanc.h & yanc.cc                           |
// v 2.0.1  31/03/2004  WD log format changed slightly; change in ext pot      |
// v 2.1    30/04/2004  WD happy icc 8.0; new body.h;                          |
// v 2.1.1  02/06/2004  WD falcON::output: moved stdout tests to io.cc         |
// v 2.1.2  28/06/2004  WD file_exists() in main.h                             |
// v 2.2    30/06/2004  WD changed options potname, etc to accname, etc        |
// v 2.2.1  08/07/2004  WD if resume, output ->out if given, else append ->in  |
// v 2.2.2  15/07/2004  WD if resume, append log output, no logout at t=t_ini  |
// v 2.2.3  16/07/2004  WD deBUGged (initialization of body flags)             |
// v 2.2.4  21/07/2004  WD velocity prediction IF external acc needs velocities|
// v 2.3    13/09/2004  WD new time integrator layout (nbody.h)                |
// v 2.4    17/09/2004  WD added manipulator                                   |
// v 2.4.1  15/10/2004  WD minor bug (!never_ending && tstop < t_ini)          |
// v 2.5    10/11/2004  WD new keyword time                                    |
// v 2.5.1  10/11/2004  WD new keyword manipinit                               |
// v 2.6    21/01/2005  WD manip may stop simulation                           |
// v 2.6.1  01/04/2005  WD adapt for change of scope of soft_type              |
// v 3.0.0  10/06/2005  WD new falcON: new bodies, new nemo_io ...             |
// v 3.0.1  23/06/2005  WD changes in warning/error                            |
// v 3.0.2  25/06/2005  WD option manippath added                              |
// v 3.0.3  20/07/2005  WD output file name appending with '!' or '@'          |
// v 3.0.4  20/07/2005  WD omit list of possible output specs                  |
// v 3.0.5  22/03/2006  WD minor bug in this file (logstep) fixed              |
// v 3.0.6  11/08/2006  WD made hmin obligatory (no default)                   |
// v 3.0.7  07/12/2006  WD renamed hmin to kmin (to prepare for SPH)           |
// v 3.0.8  28/02/2007  WD replaced kmin by kmax=kmin-Nlev+1                   |
// v 3.0.9  28/02/2007  WD added diagnostic output: log_2(root radius), depth  |
// v 3.1    15/10/2007  WD added keyword fsink                                 |
// v 3.1.1  16/10/2007  WD added keyword ksink                                 |
// v 3.2    21/11/2007  WD abandoned secondary output (use a manipulator)      |
//                         fixed minor problem with output when resume=t       |
// v 3.2.1  05/02/2008  WD read data needed by manipulator from input file     |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "3.2.1"
#define falcON_VERSION_D "05-feb-2008 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need "NEMO" to compile gyrfalcON
#endif
#define falcON_RepAction 1                         // do action reporting       
#include <public/nbody.h>                          // the N-body code           
#include <public/manip.h>                          // N-body manipulators       
#include <main.h>                                  // main & NEMO stuff         
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n           input file                                         ",
  "out=\n             file for primary output; required, unless resume=t ",
  "tstop=\n           final integration time [default: never]            ",
  "step=1\n           time between primary outputs; 0 -> every step      ",
  "logfile=-\n        file for log output                                ",
  "stopfile=\n        stop simulation as soon as file exists             ",
  "logstep=1\n        # blocksteps between log outputs                   ",
  "theta="falcON_THETA_TEXT
  "\n                 tolerance parameter at M=M_tot                     ",
#ifdef falcON_PROPER
  "fsink=1\n          theta_sink / theta (only values <=1 are used)      ",
#endif
  "hgrow=0\n          grow fresh tree every 2^hgrow smallest steps       ",
  "Ncrit="falcON_NCRIT_TEXT
  "\n                 max # bodies in un-split cells                     ",
#ifdef falcON_INDI
#  ifdef falcON_ADAP
  "eps=0.05\n         >=0: softening length OR maximum softening length\n"
  "                   < 0: use individual but FIXED softening lengths    ",
#  else
  "eps=0.05\n         >=0: softening length\n"
  "                   < 0: use individual fixed softening lengths        ",
#  endif
#else
  "eps=0.05\n         softening length                                   ",
#endif
  "kernel="falcON_KERNEL_TEXT
  "\n                 softening kernel of family P_n (P_0=Plummer)       ",
#ifdef falcON_ADAP
  "Nsoft=0\n          if >0: use individual adaptive eps_i with\n"
  "                   approx Nsoft bodies in eps spheres                 ",
  "Nref=16\n          if using eps_i: size of cell for estimating n      ",
  "emin=0\n           if using eps_i: lower limit for eps_i              ",
#endif
  "kmax=???\n         tau_max = (1/2)^kmax  MUST be given                ",
  "Nlev=1\n           # time-step levels                                 ",
#ifdef falcON_PROPER
  "ksink=\n           min tau for sinks=(1/2)^ksink (default: kmax)      ",
#endif
  "fac=\n             tau = fac / acc           \\   If more than one of  ",
  "fph=\n             tau = fph / pot            |  these is non-zero,   ",
  "fpa=\n             tau = fpa * sqrt(pot)/acc  |  we use the minimum   ",
  "fea=\n             tau = fea * sqrt(eps/acc) /   tau.                 ",
  "time=\n            time of input snapshot (default: first)            ",
  "resume=f\n         resume old simulation?  that implies:\n"
  "                   - read last snapshot from input file\n"
  "                   - append primary output to input (unless out given)",
  "give=mxv\n         list of output specifications.                     ",
//   "give2=mxv\n        list of specifications for secondary output        ",
  "Grav=1\n           Newton's constant of gravity (0-> no self-gravity) ",
  "root_center=\n     if given (3 numbers), forces tree-root centering   ",
  "accname=\n         name of external acceleration field                ",
  "accpars=\n         parameters of external acceleration field          ",
  "accfile=\n         file required by external acceleration field       ",
  "manipname=\n       name of run-time manipulator                       ",
  "manippars=\n       parameters for manipulator                         ",
  "manipfile=\n       data file required by manipulator                  ",
  "manippath=\n       path to search for manipulator                     ",
  "manipinit=f\n      manipulate initial snapshot?                       ",
  "startout=t\n       primary output for t=tstart?                       ",
  "lastout=t\n        primary output for t=tstop?                        ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string usage = "gyrfalcON -- a superb N-body code";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  // 1. set some parameters                                                     
  const bool never_ending = !hasvalue("tstop");       // integrate forever?     
  const bool resume = getbparam("resume");            // resume old simulation ?
  const bool stopfile = hasvalue ("stopfile");        // use stopfile?          
  const bool lastout = getbparam("lastout");          // write out last snapshot
  const double t_end = getdparam_z("tstop");          // integrate until t=t_end
  const double dt_out = getdparam("step");            // output interval        
  const int Nlev = getiparam("Nlev");                 // # time step levels     
  const int logstep = getiparam("logstep");           // # blocksteps/logoutput 
  const fieldset write = getioparam("give");          // what to output?        
  const nemo_acc*aex = hasvalue("accname")?           // IF(potname given) THEN 
    new nemo_acc(getparam  ("accname"),               //   initialize external  
		 getparam_z("accpars"),               //   accelerations        
		 getparam_z("accfile")) : 0;          // ELSE: no potential     
  const Manipulator MANIP(getparam_z("manipname"),    // IF(manipname given)    
			  getparam_z("manippars"),    //   THEN                 
			  getparam_z("manipfile"),    //   initialize N-body    
			  getparam_z("manippath"));   //   manipulator          
  const fieldset need(MANIP? MANIP.need() :           // additional data to be  
		      fieldset(fieldset::empty));     //   read from data input 
  if(Nlev>1 &&                                        // IF(more than one level)
     ! (hasvalue("fac") || hasvalue("fph") ||         //   we need fac or fph   
	hasvalue("fpa") || hasvalue("fea") ))         //        or fpa or fea   
     falcON_THROW("fac, fph, fap, or fea required if Nlev>1"); //error otherwise
  // 2. initialize N-body integrator                                            
#ifdef falcON_ADAP
  if(getrparam  ("eps")   <  0 &&                     // IF(indiv fixed eps_i   
     getrparam  ("Nsoft") != 0)                       // AND Nsoft != 0         
    falcON_THROW("eps<0 && Nsoft!=0: combination not sensible"); // issue error 
#endif
  vect X0;                                            // potential root center  
  FalcONCode NBDY(getparam   ("in"),                  //   snapshot input       
		  resume,                             //   resume old simul?    
		  getiparam  ("kmax"),                //   tau_max=2^(-kmax)    
		  getiparam  ("Nlev"),                //   # time steps         
#ifdef falcON_PROPER
		  hasvalue   ("ksink")?               //   tau_sink=(1/2)^ksink 
		  getiparam  ("ksink") :
#endif
		  getiparam  ("kmax"),
 		  getrparam_z("fac"),                 //   fac in adapting steps
		  getrparam_z("fph"),                 //   fph in adapting steps
		  getrparam_z("fpa"),                 //   fpa in adapting steps
		  getrparam_z("fea"),                 //   fea in adapting steps
		  getiparam  ("Ncrit"),               //   min# bodies/leaf cell
		  getiparam  ("hgrow"),               //   growing fresh tree   
		  getvparam_z("root_center",X0),      //   root centering ?     
		  getrparam  ("eps"),                 //   softening length     
		  kern(getiparam  ("kernel")),        //   softening kernel     
		  aex,                                //   external potential   
		  getrparam  ("theta"),               //   tolerance parameter  
		  getrparam  ("Grav"),                //   Newton's constant    
#ifdef falcON_PROPER
		  getrparam  ("fsink"),               //   theta_sink/theta     
#else
		  one,
#endif
#ifdef falcON_INDI
#  ifdef falcON_ADAP
		  getrparam  ("Nsoft"),               //   #in eps sphere       
		  getiparam  ("Nref"),                //   #bodies in n-estimate
		  getiparam  ("emin"),                //   lower limit for eps_i
#  endif
		  getrparam ("eps")   <zero?          //   type of softening    
		  individual_fixed    :
#ifdef falcON_ADAP
		  getrparam ("Nsoft") >zero?
		  individual_adaptive :
#endif
		  global_fixed,
#endif
		  getparam_z("time"),                 //   initial time         
		  need);                              //   read these data too  
  if(!never_ending && t_end < NBDY.initial_time()) {  // IF(t_end < t_start)    
    warning("tstop < t_ini: nothing to be done\n");   //   THEN we are done     
    return;                                           //   and can stop here    
  }                                                   // ENDIF                  
  if(MANIP) {                                         // IF manipulating        
    const_cast<snapshot*>(NBDY.my_snapshot())->       //   supported required   
      add_fields(MANIP.provide());                    //   data                 
    if(getbparam("manipinit"))                        //   IF(manipulating 1st) 
      MANIP(NBDY.my_snapshot());                      //     do it now          
  }                                                   // ENDIF                  
  // 3. open output streams & make initial outputs                              
  output LOGOUT(getparam("logfile"),resume);          // log output stream      
  if(!resume && !hasvalue("out"))                     // check for out=         
    falcON_THROW("you must provide an output file");  //   error if not given   
  nemo_out OUT(hasvalue("out")?                       // open NEMO output       
	       getparam("out") : getparam("in"),      //   if not out given,    
	       !hasvalue("out"));                     //   append to input      
  bool written=false;                                 // current snapshot wrtten
  if(!resume && getbparam("startout")) {              // IF(not resuming)       
    NBDY.write(OUT,write);                            //   write snapshot       
    written = true;                                   //   record writing       
  }                                                   // ENDIF                  
  if(LOGOUT) {                                        // IF any log output      
    NBDY.describe  (LOGOUT);                          //   put history to logout
    NBDY.stats_head(LOGOUT);                          //   header for logout    
    if(!LOGOUT.is_appending())                        //   Unless appending     
      NBDY.stats   (LOGOUT);                          //     statistics ->logout
  }                                                   // ENDIF                  
  // 4. time integration & outputs                                              
  double t_out = NBDY.time()+0.999999*dt_out;         // time for next output   
  for(int steps=1;                                    // blockstep counter      
      (never_ending || NBDY.time() < t_end) &&        // WHILE t < t_end        
      (!stopfile ||                                   //   AND                  
       !file_exists(getparam("stopfile")));           //    no stopfile         
      ++steps) {                                      //   increment counter    
    NBDY.full_step();                                 //   make full block step 
    if(LOGOUT && steps%logstep ==0)                   //   IF(time for logout)  
      NBDY.stats(LOGOUT);                             //     statistics output  
    if(MANIP)                                         //   IF(manipulating)     
      if(MANIP(NBDY.my_snapshot()))                   //     IF(manip says so)  
	break;                                        //       STOP simulation  
    if(OUT && NBDY.time() >= t_out) {                 //   IF(t >= t_out)       
      NBDY.write(OUT,write);                          //     primary output     
      t_out += dt_out;                                //     increment t_out    
      written = true;                                 //     written out        
    } else                                            //   ELSE                 
      written = false;                                //     not written        
  }                                                   // END: WHILE             
  if(OUT && !written && lastout)                      // IF not yet done:       
    NBDY.write(OUT,write);                            //   write last snapshot  
  if(LOGOUT && stopfile && file_exists(getparam("stopfile")))
    LOGOUT <<"# simulation STOPPED because file \""
	   << getparam("stopfile") << "\" found to exist\n";
  // 5. cleaning up (including implicit call of destructors)                    
  if(aex) falcON_DEL_O(aex);                          // delete external accs   
}
//---------------------end-of-gyrfalcON.cc------that's-it-!---------------------
