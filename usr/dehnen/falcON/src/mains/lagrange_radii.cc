// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// lagrange_radii.cc                                                           |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 1.0   21/06/2002  WD created                                              |
// v 1.1   12/07/2002  WD added option masses                                  |
// v 1.2   20/08/2002  WD added options centering, alpha, and times            |
// v 1.3   29/08/2002  WD improved find_center()                               |
// v 1.4   30/08/2002  WD adapted this file for usage of MPI otherwise         |
// v 1.5   06/12/2002  WD debugged                                             |
// v 1.6   08/01/2003  WD abandoned centering, use 'center' instead            |
// v 2.0   13/03/2003  WD added stop mechanism; interpolate in R^2,            |
//                         => the actual radii may be slightly different       |
// v 2.1   20/03/2003  WD action reporting                                     |
// v 2.2   23/05/2003  WD automated NEMO history                               |
// v 2.3   29/10/2003  WD automated version etc; changed give: no default      |
// v 2.4   13/02/2004  WD debugged error with give=                            |
// v 2.5   24/03/2004  WD added 3% to the default masses list                  |
// v 3.0   07/05/2004  WD new sorting, 10 times faster than old version        |
// v 3.1   10/05/2004  WD fixed two bugs (allow for >1 @ same r, r=0)          |
// v 3.2   11/05/2004  WD made PUBLIC; changed name (origina: "lagrange_rad")  |
// v 3.2.1 19/05/2004  WD change of give: if not given, we write all we got    |
// v 3.3   19/05/2004  WD re-written completely; sanity check for stdouts      |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "3.3"
#define falcON_VERSION_D "19-may-2004 Walter Dehnen                          "
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
#include <public/nmio.h>                           // my NEMO I/O               
#include <public/tool.h>                           // my tools for bodies       
#include <public/ionl.h>                           // my I/O utilities          
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
namespace {
  using namespace nbdy;
  inline std::ostream& open_ofstream(const char* file)
  {
    std::ofstream *out = new std::ofstream();
    if(! open_to_append(*out,file) ) nbdy::exit(1);
    return *(static_cast<std::ostream*>(out));
  }
}
//------------------------------------------------------------------------------
void nbdy::main()
{
  // sanity check for number of outputs to stdout                               
  if( ( 0 == strcmp(getparam("tabfile"), "-") ) &&
      ( hasvalue("out") && 0==strcmp(getparam("out"),"-") ) )
    ::error("\"out=-\" and \"tabfile=-\": cannot write both to stdout\n");
  // set snapshot I/O                                                           
  const bool     SOUT(hasvalue("out") && strcmp(getparam("out"),"."));
  const nemo_in  IN  (getparam("in"));
  const nemo_out OUT (SOUT? getparam("out") : ".");
  const io       GIVE(hasvalue("give")? getioparam("give") : io::all);
  const io       WANT(SOUT? GIVE : io::mx);
  // set parameters for Lagrange radii                                          
  const int NMAX    = 100;
  double    M[NMAX] = {0.}, R[NMAX];
  int       N       = nemoinpd(getparam("masses"),M,NMAX);
  if(N>NMAX) ::error("too many mass shells specified");
  // open output for Lagrange radii                                             
  std::ostream &TABOUT = strcmp(getparam("tabfile"), "-")?
    open_ofstream(getparam("tabfile")) : std::cout;
  TABOUT << "# lagrange_radii in="<<getparam("in")
	 << " tabfile="<<getparam("tabfile")
	 << ":\n# time       ";
  for(int j=0; j!=N; ++j)
    TABOUT << "r[" << std::setw(4) << 100*M[j] << "%] ";
  TABOUT << std::endl;
  TABOUT.setf(std::ios::left, std::ios::adjustfield);
  // set stop parameters                                                        
  bool         STOPPING       = hasvalue   ("stopfile");
  const int    STOPINDEX      = getiparam  ("stopindex");
  const real   STOPVALUE      = getdparam_z("stopvalue");
  const bool   STOPRELATIVE   = getbparam  ("stoprelative");
  double       STOPAFTER      = getdparam_z("stopafter");
  const double STOPDELAY      = getdparam  ("stopdelay");
  if(STOPPING && (STOPINDEX < 0 || STOPINDEX >= N)) {
    ::warning("will not use stop condition, since stopindex not in [0,%d]",N-1);
    STOPPING = false;
  }
  if(STOPPING && STOPVALUE == zero) {
    ::warning("will not use stop condition, since stopvalue=0");
    STOPPING = false;
  }
  if(STOPPING && 
     STOPRELATIVE && (STOPVALUE>one || (STOPVALUE<zero && STOPVALUE>-one)))
    ::warning("stop condition already true at initial time");
  bool STOPPED = false, HAS_RI0 =false;
  register double RI, RI0;
  real STOPTIME, LASTTIME;
  // loop snapshots & process them                                              
  bodies   BODIES;
  char     WORD[30];
  io       READ;
  double   TIME;
  bool     FIRST = true;
  double   TOUT, STEP(getdparam("step"));
  while(IN.is_present(nemo_io::snap)) {
    // read time, read snapshot if in times                                     
    bool IN_TIMES = BODIES.read_nemo_snapshot(IN,READ,&TIME,WANT,
					      getparam("times"),0);
    if(FIRST) {
      TOUT  = TIME - 1.e-10*STEP;
      FIRST = false;
    }
    if(! IN_TIMES ) continue;
    // output of snapshot                                                       
    if(SOUT && TIME >= TOUT) {
      BODIES.write_nemo_snapshot(OUT,&TIME,READ&GIVE);
      TOUT += STEP;
    }
    // find & write Lagrange radii                                              
    if(!READ.contains(io::mx)) {                   //   too few data: abort     

      ::error("insufficient data: need mx, got %s", READ.make_word(WORD));
    }
    TABOUT << std::setw(12) << TIME <<' ';
    find_lagrange_rad(&BODIES,N,M,R);
    for(int i=0; i!=N; ++i)
      TABOUT << std::setw(8) << R[i] << ' ';
    TABOUT << std::endl;
    // deal with STOPPING of simulation                                         
    if(STOPPING) {
      RI = R[STOPINDEX];
      if(!HAS_RI0) { RI0 = RI; HAS_RI0 = true; }
      if(!STOPPED) {
	register real STOPX = STOPRELATIVE? RI/RI0 : RI;
	if((STOPVALUE > zero && STOPX <     STOPVALUE  ) ||
	   (STOPVALUE < zero && STOPX > abs(STOPVALUE) ) ) {
	  STOPPED  = true;
	  STOPTIME = TIME;
	  LASTTIME = max(TIME+STOPDELAY,STOPAFTER);
	}
      }
      if(STOPPED && TIME >= LASTTIME) {
	std::ofstream STOP(getparam("stopfile"));
	STOP<<" stopfile \""<<getparam("stopfile")<<"\"\n"
	    <<" generated by program \"lagrange_radii\"\n"
	    <<" on simulation time "<<TIME<<", because the stop condition\n"
	    <<"\n   R(M="<<M[STOPINDEX]<<"*M_tot)";
	if(STOPRELATIVE) STOP<<"/R(M="<<M[STOPINDEX]<<"*M_tot, t=t_initial)";
	STOP<< ((STOPVALUE<zero)? " > " : " < ")
	    << abs(STOPVALUE) <<"\n\n"
	    <<" was satisfied at simulation time "<<STOPTIME<<"\n";
      }
    }
  }
}
