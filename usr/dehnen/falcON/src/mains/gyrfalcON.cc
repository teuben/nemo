// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// gyrfalcON.cc                                                                |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2000-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 1.0.0   28/05/2001  WD created with the help of PJT, based on YANC        |
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
// v 1.5.1   09/01/2003  WD saver C-macros, default parameters                 |
// v 1.5.2   13/01/2003  WD never ending; logstep; stopfile                    |
// v 1.5.3   24/01/2003  WD option lastout added.                              |
// v 1.5.4   03/03/2003  WD debugged: handling of external pot in total energy |
// v 1.5.5   13/03/2003  WD added check for over-using of stdout/pipe          |
// v 1.5.6   17/03/2003  WD options emin, fea, limits on |eps_new/eps_old|     |
// v 1.5.7   20/03/2003  WD changes in gravity, action reporting (proper only) |
// v 1.6     02/06/2003  WD allow for individual but fixed eps_i by eps<0      |
// v 1.6.1   28/07/2003  WD happy gcc 3.3 (about 6% faster than gcc 3.2)       |
// v 1.6.2   08/08/2003  WD version provides compiler info, automated          |
// v 1.6.3   13/08/2003  WD fixed bug with individual_fixed; thanks to J.Bailin|
// v 1.7.0   05/09/2003  WD individual eps made public; changes in tensors     |
// v 1.7.1   17/09/2003  WD changes in tree: avoiding template specs           |
// v 1.7.2   07/10/2003  WD changes in gravity: using common basic_tree        |
// v 1.7.3   23/10/2003  WD changes in design of gravity, tree, kernel, falcON |
// v 1.8     05/11/2003  WD changes in grav; changed io::P to io::q            |
// v 1.8.1   10/02/2004  WD minor change in this file (logfile)                |
// v 1.8.2   11/02/2004  WD allow for Grav=0 (Grav now handled in grav.h)      |
// v 1.8.3   18/02/2004  WD bug with Grav=0 fixed; avoid tree building if G=0  |
// v 1.9     19/02/2004  WD added option root_center                           |
// v 1.9.1   23/02/2004  WD improved diagnose output (new T, V_in, [V_ex,] W)  |
// v 1.9.2   27/02/2004  WD use nemo::error() instead of nbdy::error()         |
// v 2.0     11/03/2004  WD elimated yanc.h & yanc.cc                          |
// v 2.0.1   31/03/2004  WD log format changed slightly; change in ext pot     |
// v 2.1     30/04/2004  WD happy icc 8.0; new body.h;                         |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "2.1"
#define falcON_VERSION_D "30-apr-2004 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need "NEMO" to compile gyrfalcON
#endif
#include <nbdy.h>                                  // the N-body code           
#include <iostream>                                // C++ I/O                   
#include <fstream>                                 // C++ file I/O              
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
  "out2=\n            file for secondary output stream                   ",
  "step2=0\n          time between secondary outputs; 0 -> every step    ",
  "theta="falcON_THETA_TEXT
  "\n                 tolerance parameter at M=M_tot                     ",
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
  "hmin=6\n           tau_min = (1/2)^hmin                               ",
  "Nlev=1\n           # time-step levels                                 ",
  "fac=\n             tau = fac / acc           \\   If more than one of  ",
  "fph=\n             tau = fph / pot            |  these is non-zero,   ",
  "fpa=\n             tau = fpa * sqrt(pot)/acc  |  we use the minimum   ",
  "fea=\n             tau = fea * sqrt(eps/acc) /   tau.                 ",
  "resume=f\n         resume old simulation?  that implies:\n"
  "                   - read last snapshot from input file\n"
  "                   - append primary output to input (unless out given)",
  "give=mxv\n         list of output specifications. Recognizing:\n"
  "                    m: mass                              (default)\n"
  "                    x: position                          (default)\n"
  "                    v: velocity                          (default)\n"
  "                    a: acceleration\n"
  "                    p: N-body potential\n"
  "                    q: add external potential before output\n"
#ifdef falcON_INDI
  "                    e: individual eps_i (if they exist)\n"
#endif
  "                    l: time-step level (if they exist)\n"
#if(0)
  "                    k: body key (if given with input)\n"
#endif
#ifdef falcON_ADAP
  "                    r: density estimate\n"
#endif
  "                    f: body flag                                      ",
  "give2=mxv\n        list of specifications for secondary output        ",
  "Grav=1\n           Newton's constant of gravity (0-> no self-gravity) ",
  "root_center=\n     if given (3 numbers), forces tree-root centering   ",
  "potname=\n         name of external potential                         ",
  "potpars=\n         parameters of external potential                   ",
  "potfile=\n         file required by external potential                ",
  "startout=t\n       primary output for t=tstart?                       ",
  "lastout=t\n        primary output for t=tstop?                        ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string usage = "gyrfalcON -- a superberb N-body code";
//------------------------------------------------------------------------------
inline bool exists(const char* stop_file)
{
  static std::ifstream IN;
  IN.open(stop_file);
  if(IN.is_open()) {
    IN.close();
    return true;
  } else
    return false;
}
//------------------------------------------------------------------------------
void nbdy::main()
{
  // 1. set some parameters                                                     
  const bool
    never_ending = !hasvalue("tstop"),                // integrate forever?     
    resume       = getbparam("resume"),               // resume old simulation ?
    stopfile     = hasvalue ("stopfile"),             // use stopfile?          
    lastout      = getbparam("lastout");              // write out last snapshot
  const double
    t_end   = getdparam_z("tstop"),                   // integrate until t=t_end
    dt_out0 = getdparam("step"),                      // primary output interval
    dt_out1 = getdparam("step2");                     // secondary  --------    
  const int
    Nlev    = getiparam("Nlev"),                      // # time step levels     
    logstep = getiparam("logstep");                   // # blocksteps/logoutput 
  const io
    wr0 = getioparam("give"),                         // what to output 0?      
    wr1 = getioparam("give2");                        // what to output 1?      
  const nemo_grav*pex = hasvalue("potname")?          // IF(potname given) THEN 
    new nemo_grav(getparam  ("potname"),              //   initialize external  
		  getparam_z("potpars"),              //   gravity              
		  getparam_z("potfile")) : 0;         // ELSE: no potential     
  if(Nlev>1 &&                                        // IF(more than one level)
     ! (hasvalue("fac") || hasvalue("fph") ||         //   we need fac or fph   
	hasvalue("fpa") || hasvalue("fea") ))         //        or fpa or fea   
    ::error("fac, fph, fap, or fea required if Nlev>1"); //error otherwise      
  // 2. initialize N-body integrator                                            
#ifdef falcON_ADAP
  if(getdparam  ("eps")   <  0 &&                     // IF(indiv fixed eps_i   
     getdparam  ("Nsoft") != 0)                       // AND Nsoft != 0         
    ::error("eps<0 && Nsoft!=0: combination not sensible"); // issue error      
#endif
  vect X0;                                            // potential root center  
  NbodyCode NBDY(getparam   ("in"),                   //   snapshot input       
		 resume,                              //   resume old simul?    
		 getiparam  ("hmin"),                 //   -log_2(time step)    
		 getiparam  ("Nlev"),                 //   # time steps         
		 getdparam_z("fac"),                  //   fac in adapting steps
		 getdparam_z("fph"),                  //   fph in adapting steps
		 getdparam_z("fpa"),                  //   fpa in adapting steps
		 getdparam_z("fea"),                  //   fea in adapting steps
		 getiparam  ("Ncrit"),                //   min# bodies/leaf cell
		 getiparam  ("hgrow"),                //   growing fresh tree   
		 getvparam_z("root_center",X0),       //   root centering ?     
		 getdparam  ("eps"),                  //   softening length     
		 kern(getiparam  ("kernel")),         //   softening kernel     
		 pex,                                 //   external potential   
		 getdparam  ("theta"),                //   tolerance parameter  
		 getdparam  ("Grav")                  //   Newton's constant    
#ifdef falcON_INDI
		 ,
#  ifdef falcON_ADAP
		 getdparam  ("Nsoft"),                //   #in eps sphere       
		 getiparam  ("Nref"),                 //   #bodies in n-estimate
		 getiparam  ("emin"),                 //   lower limit for eps_i
#  endif
		 getdparam ("eps")   <zero?
		 basic_nbody::individual_fixed    :
#ifdef falcON_ADAP
		 getdparam ("Nsoft") >zero?
		 basic_nbody::individual_adaptive :
#endif
		 basic_nbody::global_fixed
#endif
		 );
  if(!NBDY.okay()) ::error("initialization error");   // Initialization okay?   
  if(t_end < NBDY.initial_time()) {                   // IF(t_end < t_start)    
    ::warning("tstop < t_ini: nothing to be done\n"); //   THEN we are done     
    return;                                           //   and can stop here    
  }                                                   // ENDIF                  
  // 3. initialize output streams                                               
  // 3.0 check that at most one output is to stdout                             
  int nstdout = 0;                                    // counter for stdout     
  if(0==strcmp(getparam("logfile"),"-") ) ++nstdout;  // is logfile=- ?         
  if(0==strcmp(getparam("out")    ,"-") ) ++nstdout;  // is out=- ?             
  if(0==strcmp(getparam("out2")   ,"-") ) ++nstdout;  // is out2=- ?            
  if(nstdout > 1)                                     // IF # stdout > 1 ERROR  
    ::error("more than one of \"logfile\", \"out\", \"out2\" equals \"-\"");
  // 3.1 primary output                                                         
  if(hasvalue("out"))                                 // IF(out given) THEN     
    NBDY.open_nemo(0,getparam("out"),resume);         //   open NEMO output     
  else if(resume)                                     // ELSE IF(resuming) THEN 
    NBDY.open_nemo(0,0,resume);                       //   append to NEMO input 
  else                                                // ELSE(no out, no resume)
    ::error("out required if resume=f");              //   error out            
  if(!NBDY.nemo_is_open(0))                           // IF(no nemo output)     
    ::error("opening for primary output failed");     //   error out            
  bool written=false;                                 // current snapshot wrtten
  if(!resume && getbparam("startout")) {              // IF(not resuming)       
    NBDY.write_nemo(wr0,0);                           //   write snapshot       
    written = true;                                   //   record writing       
  }                                                   // ENDIF                  
  // 3.2 secondary output                                                       
  if(hasvalue("out2")) {                              // IF(out2 given) THEN    
    NBDY.open_nemo(1,getparam("out2"));               //   try to open output   
      if(!NBDY.nemo_is_open(1))                       //   IF(no nemo output)   
	::error("opening for secondary output failed");
    NBDY.write_nemo(wr1,1);                           //   write snapshot       
  }                                                   // ENDIF                  
  // 3.3 log output                                                             
  std::ostream *logout = 0;                           // pointer to log ostream 
  std::ofstream logfile;                              // file for log output    
  if(0 == strcmp(getparam("logfile"),"."))            // IF no log output       
    ::warning("option \"logfile=.\" supresses any log output");
  else if(strcmp(getparam("logfile"),"-")) {          // ELIF desired:          
    logfile.open(getparam("logfile"));                //   open logfile         
    logout = &logfile;                                //   log output to file   
  } else                                              // ELSE                   
    logout = &std::clog;                              //   log output to stdlog 
  if(logout) {                                        // IF any log output      
    NBDY.describe_nemo(*logout,*(ask_history()));     //   put history to logout
    NBDY.stats_head   (*logout);                      //   header for logout    
    NBDY.stats        (*logout);                      //   statistics -> logout 
  }                                                   // ENDIF                  
  // 4. time integration & outputs                                              
  double
    t_out0 = NBDY.initial_time()+0.999999*dt_out0,    // time for next output 0 
    t_out1 = NBDY.initial_time()+0.999999*dt_out1;    // time for next output 1 
  for(int steps=0;                                    // blockstep counter      
      (never_ending || NBDY.time() < t_end) &&        // WHILE t < t_end        
      (!stopfile || !exists(getparam("stopfile")));   //   AND no stopfile      
      ++steps) {                                      //   increment counter    
    NBDY.full_step();                                 //   make full block step 
    if(logout && steps%logstep ==0)                   //   IF(time for logout)  
      NBDY.stats(*logout);                            //     statistics output  
    if(NBDY.time() >= t_out0) {                       //   IF(t >= t_out0)      
      NBDY.write_nemo(wr0,0);                         //     primary output     
      t_out0 += dt_out0;                              //     increment t_out0   
      written = true;                                 //     written out        
    } else                                            //   ELSE                 
      written = false;                                //     not written        
    if(NBDY.nemo_is_open(1) &&                        //   IF(secondary output  
       NBDY.time() >= t_out1) {                       //   AND t >= t_out1)     
      NBDY.write_nemo(wr1,1);                         //     secondary output   
      t_out1 += dt_out1;                              //     increment t_out1   
    }                                                 //   ENDIF                
  }                                                   // END: WHILE             
  if(!written && lastout) NBDY.write_nemo(wr0,0);     // write last snapshot    
  if(logout && stopfile && exists(getparam("stopfile")))
    (*logout) <<"# simulation STOPPED because file \""
	      << getparam("stopfile") << "\" found to exist\n";
  // 5. cleaning up (including implicit call of destructors)                    
  if(pex) delete pex;                                 // delete external pot    
}
//---------------------end-of-gyrfalcON.cc------that's-it-!---------------------
