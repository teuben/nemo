// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// nemo/gyrfalcON.cc                                                           |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2002                                          |
// e-mail:   wdehnen@aip.de                                                    |
// address:  Astrophysikalisches Institut Potsdam,                             |
//           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 1.0.0   28/05/2001  WD created with the help of NEMO caption PJT          |
// v 1.0.1   01/06/2001  WD bug in this file removed, added history output     |
// v 1.0.2   13/06/2001  WD Ncrit option added (before Ncrit==1)               |
// v 1.0.3   20/06/2001  WD tree::make() back to adding-leafs algorithm        |
//                          Coordinatesystem given with snapshot               |
// v 1.0.4   25/06/2001  WD bug (builder::link_tree()) removed (thanks to PJT) |
// v 1.0.5   28/06/2001  WD changes in tree::grow(). Avoid box overflow.       |
//                          bug (in builder) removed                           |
// v 1.0.6   01/08/2001  WD substantial code upgrade in the tree (falcON)      |
// v 1.0.7   02/08/2001  WD Nreuse option added (before Nreuse==0)             |
// v 1.0.8   08/08/2001  WD new data lay-out for bodies, enabling special      |
//                          treatment for leap-frog, saves memory & time (?)   |
// v 1.0.9   17/09/2001  WD added support for dynamically linked external      |
//                          potential (potential.h)                            |
// v 1.0.10  20/09.2001  WD added option ("resume") for resuming an old or     |
//                          interrupted simulation. Appends to input file      |
// v 1.0.11  21/09/2001  WD improved handling of nemo I/O; can now also read   |
//                          PosTag & VelTag instead of PhaseSpaceTag           |
// v 1.0.12  21/09/2001  WD added parameter give_acc                           |
// v 1.0.13  11/10/2001  WD changed give_pot option to allow only N-body pot   |
//                          added option startout                              |
// v 1.0.14  18/10/2001  WD interweaving interaction & evaluation phase of     |
//                          gravity approximation: saves about 20b/body        |
// v 1.0.15  23/10/2001  WD some changes in I/O handling (body, yanc)          |
// v 1.0.16  23/11/2001  WD added option logout (to allow output pipe)         |
// v 1.0.17  08/01/2002  WD changed body lay-out (gives 2% speed-up)           |
// v 1.1.0   21/01/2002  WD minor changes                                      |
// v 1.1.1   25/01/2002  WD minor changes                                      |
// v 1.1.2   29/01/2002  WD minor changes (1.5% speed-up)                      |
// v 1.1.3   26/02/2002  WD time-symmetric stepping criterion tau=fac/acc      |
// v 1.1.4   27/02/2002  WD added options: out2, step2                         |
// v 1.1.5   01/03/2002  WD added time stepping criterion tau=fp/pot           |
// v 1.1.6   05/03/2002  WD added time stepping criterion tau=fc*sqrt(pot)/acc |
// v 1.2     07/06/2002  WD replaced giveacc,givepot,giverho with give, give2  |
// v 1.2.1   11/06/2002  WD support for kernels Fn and Kn withdrawn            |
// v 1.2.2   14/06/2002  WD added output option for flag & level               |
// v 1.2.3   17/06/2002  WD allow for individual adaptive softening lengths    |
//                          changed option Nreuse -> hgrow                     |
// v 1.2.4   26/08/2002  WD bug in MAC removed,                                |
//                          tree re-build (saves 40% CPU time on build)        |
// v 1.2.5   28/08/2002  WD bug with hgrow removed                             |
// v 1.2.6   29/08/2002  WD improved initialization of external potential      |
// v 1.2.7   30/08/2002  WD adapted this file for usage of MPI otherwise       |
// v 1.2.8   09/09/2002  WD further adaption to MPI usage                      |
// v 1.2.9   05/11/2002  WD added option Grav (for comparison with GADGET)     |
// v 1.3.0   15/11/2002  WD various updates (SSE code, eps_i treatment)        |
// v 1.4     20/11/2002  WD splitted between public and proprietary code       |
// v 1.4.1   26/11/2002  WD debugged some features                             |
// v 1.5     04/12/2002  WD re-named "gyrfalcON" (previously "YancNemo")       |
//                          several changes in file layout for public version  |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef ALLOW_NEMO                                    // this is a NEMO program 
#  error You need "ALLOW_NEMO" to compile __FILE__
#endif
#include <yanc.h>                                     // the N-body code        
#include <public/pext.h>                              // external potential     
#include <iostream>                                   // C++ I/O                
#include <fstream>                                    // C++ file I/O           
#include <nemomain.h>                                 // NEMO main              
#include <nemo.h>                                     // NEMO basics & main     

using namespace nbdy;
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n            input file                                         ",
  "out=\n              file for primary output; required, unless resume=t ",
  "tstop=10\n          final integration time                             ",
  "step=1\n            time between primary outputs; 0 -> every step      ",
  "logfile=-\n         file for log output                                ",
  "out2=\n             file for secondary output stream                   ",
  "step2=0\n           time between secondary outputs; 0 -> every step    ",
  "theta=0.6\n         tolerance parameter at M=M_tot                     ",
  "hgrow=0\n           grow fresh tree every 2^hgrow smallest steps       ",
#ifndef DEFAULT_NCRIT
# error "unknown default for Ncrit"
#elif DEFAULT_NCRIT ==  1
  "Ncrit=1\n           max # bodies in un-split cells                     ",
#elif DEFAULT_NCRIT ==  6
  "Ncrit=6\n           max # bodies in un-split cells                     ",
#elif DEFAULT_NCRIT ==  8
  "Ncrit=8\n           max # bodies in un-split cells                     ",
#elif DEFAULT_NCRIT == 12
  "Ncrit=12\n          max # bodies in un-split cells                     ",
#elif DEFAULT_NCRIT == 16
  "Ncrit=16\n          max # bodies in un-split cells                     ",
#else 
# error "unknown default for Ncrit"
#endif
  "eps=0.05\n          softening length OR maximum softening length       ",
  "kernel=1\n          softening kernel of family P_n (P_0=Plummer)       ",
#ifdef ALLOW_INDI
  "Nsoft=0\n           if 0:  use global softening length eps\n"
  "                   else:  use eps_i with ~Nsoft bodies in eps spheres ",
  "Nref=16\n           if using eps_i: size of cell for estimating n      ",
#endif
  "hmin=6\n            tau_min = (1/2)^hmin                               ",
  "Nlev=1\n            # time-step levels                                 ",
  "fac=\n              tau = fac/acc           \\   If more than one of    ",
  "fph=\n              tau = fph/pot            |  these is non-zero, we  ",
  "fpa=\n              tau = fpa*sqrt(pot)/acc /   use the minimum tau.   ",
  "resume=f\n          resume old simulation?  that implies:\n"
  "                   - read last snapshot from input file\n"
  "                   - append primary output to input (unless out given)",
  "give=mxv\n          list of output specifications. Recognizing:\n"
  "                     m: mass                              (default)\n"
  "                     x: position                          (default)\n"
  "                     v: velocity                          (default)\n"
  "                     a: acceleration\n"
#ifdef ALLOW_INDI
  "                     e: individual eps_i (if they exist)\n"
#endif
  "                     p: N-body potential\n"
  "                     l: time-step level (if they exist)\n"
#if(0)
  "                     k: body key (if given with input)\n"
#endif
#ifdef ALLOW_INDI
  "                     f: body flag\n"
  "                     r: density estimate                              ",
#else
  "                     f: body flag                                     ",
#endif
  "give2=mxv\n         list of specifications for secondary output        ",
  "Grav=1\n            Newton's constant of gravity                       ",
  "potname=\n          name of external potential                         ",
  "potpars=\n          parameters of external potential                   ",
  "potfile=\n          file required by external potential                ",
  "startout=t\n        primary output for t=tstart?                       ",
  "VERSION=1.5\n       04/December/2002  WD\n"
  "                   compiled  " __DATE__ ", " __TIME__ "                    ",
  NULL};
//------------------------------------------------------------------------------
string usage = "gyrfalcON -- a superberb N-body code";
//------------------------------------------------------------------------------
void nemo::main()
{
  // 1. set some parameters                                                     
  register double t_out0,t_out1,                      // t_out0, t_out1         
    t_end   = getdparam("tstop"),                     // integrate until t=t_end
    dt_out0 = getdparam("step"),                      // primary output interval
    dt_out1 = getdparam("step2");                     // secondary  --------    
  register int Nlev = getiparam("Nlev");              // # time step levels     
  if(Nlev>1 &&                                        // IF(more than one level)
     ! (hasvalue("fac") || hasvalue("fph")))          //   we need fac or fph   
    error("fac or fph required if Nlev > 1");         //   error otherwise      
  nemo_pot *pot = hasvalue("potname")?                // IF(potname given) THEN 
    new nemo_pot(getparam  ("potname"),               //   initialize external  
		 getparam_z("potpars"),               //   potential            
		 getparam_z("potfile"))               //                        
    : 0;                                              // ELSE: no potential     
  bool resume = getbparam("resume");                  // resume old simulation ?
  io wr0(getparam("give")), wr1(getparam("give2"));   // what to output?        
  // 2. initialize N-body integrator                                            
  yanc YANC(getparam   ("in"),                        //   snapshot input       
	    true,                                     //   do nemo I/O          
	    getdparam  ("theta"),                     //   tolerance parameter  
	    getiparam  ("hgrow"),                     //   growing fresh tree   
	    getiparam  ("Ncrit"),                     //   min# bodies/leaf cell
	    getdparam  ("eps"),                       //   softening length     
	    getiparam  ("kernel"),                    //   softening kernel     
	    getiparam  ("hmin"),                      //   -log_2(time step)    
	    Nlev,                                     //   # time-step levels   
	    getdparam_z("fac"),                       //   fac in adapting steps
	    getdparam_z("fph"),                       //   fph in adapting steps
	    getdparam_z("fpa"),                       //   fpa in adapting steps
#ifdef ALLOW_INDI
	    getdparam  ("Nsoft"),                     //   #in eps sphere       
	    getiparam  ("Nref"),                      //   #bodies in n-estimate
	    getdparam  ("Nsoft")!=0 ? 2 : 0,          //   softening type       
#endif
	    getdparam  ("Grav"),                      //   Newton's constant    
	    resume,                                   //   resume old simul?    
	    pot);                                     //   external potential   
  if(!YANC.okay()) error("initialization error");     // Initialization okay?   
  if(t_end < YANC.initial_time()) {                   // IF(t_end < t_start)    
    warning("tstop < t_ini: nothing to be done\n");   //   THEN we are done     
    return;                                           //   and can stop here    
  }                                                   // ENDIF                  
  // 3. initialize output streams                                               
  // 3.1 primary output                                                         
  if(hasvalue("out"))                                 // IF(out given) THEN     
    YANC.open_nemo(0,getparam("out"),resume);         //   open NEMO output     
  else if(resume)                                     // ELSE IF(resuming) THEN 
    YANC.open_nemo(0,0,resume);                       //   append to NEMO input 
  else                                                // ELSE(no out, no resume)
    error("out required if resume=f");                //   error out            
  if(!YANC.nemo_is_open(0))                           // IF(no nemo output)     
    error("opening for primary output failed");       //   error out            
  if(!resume && getbparam("startout"))                // IF(not resuming)       
    YANC.write_nemo(wr0,0);                           //   write snapshot       
  // 3.2 secondary output                                                       
  if(hasvalue("out2")) {                              // IF(out2 given) THEN    
    YANC.open_nemo(1,getparam("out2"));               //   try to open output   
      if(!YANC.nemo_is_open(1))                       //   IF(no nemo output)   
	error("opening for secondary output failed"); //     error out          
    YANC.write_nemo(wr1,1);                           //   write snapshot       
  }                                                   // ENDIF                  
  // 3.3 log output                                                             
  std::ostream *logout;                               // pointer to log ostream 
  if(0 == strcmp(getparam("logfile"),"."))
    warning("option \"logfile=.\" supresses any log output");
  if(strcmp(getparam("logfile"),"-"))                 // IF desired             
    logout = new std::ofstream(getparam("logfile"));  //   log output to file   
  else                                                // ELSE                   
    logout = &std::clog;                              //   log output to stdlog 
  if(logout) {
    YANC.describe_nemo(*logout,*(ask_history()));     // put history to logout  
    YANC.stats_head   (*logout);                      // header for logout      
    YANC.stats  (*logout);                            // statistics -> logout   
  }
  // 4. time integration & outputs                                              
  t_out0 = YANC.initial_time()+0.999999*dt_out0;      // time for next output   
  t_out1 = YANC.initial_time()+0.999999*dt_out1;      // time for next output   
  while(YANC.time() < t_end) {                        // WHILE( t < t_end )     
    YANC.full_step();                                 //   make full block step 
    if(logout) YANC.stats(*logout);                   //   statistics output    
    if(YANC.time() >= t_out0) {                       //   IF(t >= t_out0)      
      YANC.write_nemo(wr0,0);                         //     primary output     
      t_out0 += dt_out0;                              //     increment t_out0   
    }                                                 //   ENDIF                
    if(YANC.nemo_is_open(1) &&                        //   IF(secondary output  
       YANC.time() >= t_out1) {                       //   AND t >= t_out1)     
      YANC.write_nemo(wr1,1);                         //     secondary output   
      t_out1 += dt_out1;                              //     increment t_out1   
    }                                                 //   ENDIF                
  }                                                   // END: WHILE             
  // 5. cleaning up (including implicit call of destructors)                    
  if(logout && strcmp(getparam("logfile"),"-"))       // IF log output to file  
    static_cast<std::ofstream*>(logout)->close();     //   close log output file
  if(pot) delete pot;                                 // delete external pot    
}
//---------------------end-of-YancNemo.cc------that's-it-!----------------------
