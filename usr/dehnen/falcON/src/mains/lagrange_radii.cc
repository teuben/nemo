// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// lagrange_radii.cc                                                           |
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
// v 1.0    21/06/2002  WD created                                             |
// v 1.1    12/07/2002  WD added option masses                                 |
// v 1.2    20/08/2002  WD added options centering, alpha, and times           |
// v 1.3    29/08/2002  WD improved find_center()                              |
// v 1.4    30/08/2002  WD adapted this file for usage of MPI otherwise        |
// v 1.5    06/12/2002  WD debugged                                            |
// v 1.6    08/01/2003  WD abandoned centering, use 'center' instead           |
// v 2.0    13/03/2003  WD added stop mechanism; interpolate in R^2,           |
//                         => the actual radii may be slightly different       |
// v 2.1    20/03/2003  WD action reporting                                    |
// v 2.2    23/05/2003  WD automated NEMO history                              |
// v 2.3    29/10/2003  WD automated version etc; changed give: no default     |
// v 2.4    13/02/2004  WD debugged error with give=                           |
// v 2.5    24/03/2004  WD added 3% to the default masses list                 |
// v 3.0    07/05/2004  WD new sorting, 10 times faster than old version       |
// v 3.1    10/05/2004  WD fixed two bugs (allow for >1 @ same r, r=0)         |
// v 3.2    11/05/2004  WD made PUBLIC; changed name (origina: "lagrange_rad") |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "3.2"
#define falcON_VERSION_D "11-may-2004 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#error You need NEMO to compile "src/mains/lagrange_radii.cc"
#endif
#define falcON_RepAction 0                         // no action reporting       
#include <iostream>                                // C++ I/O                   
#include <fstream>                                 // C++ file I/O              
#include <iomanip>                                 // C++ I/O formatting        
#include <body.h>                                  // the bodies                
#include <public/nmio.h>                           // my NEMO I/O               
#include <public/ionl.h>                           // my I/O utilities          
#include <public/tool.h>                           // my N-body tools           
#include <main.h>                                  // main & NEMO stuff         
using namespace nbdy;
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n           input file                                         ",
  "masses=0.01,0.03,0.1,0.3,0.5,0.7,0.9,0.97,0.99\n"
  "                   table of masses; will be parsed by nemoinp         ",
  "times=all\n        times to process                                   ",
  "tabfile=\n         file to append table to; [stdout]                  ",
  "out=\n             output file [default: no output]                   ",
  "step=0\n           step between outputs                               ",
  "give=\n            what data to put out                               ",
  "stopfile=\n        create stop file                                   ",
  "stopindex=0\n      index for lagrange rad used in stop, first=0       ",
  "stopvalue=\n       stop if R_i < value ( or > |value| )               ",
  "stoprelative=f\n   use R_i / R_i[t=t_initial] in stop criterion       ",
  "stopafter=\n       don't stop before this time (default: t_initial)   ",
  "stopdelay=0\n      delay stopping after condition is satisfied        ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string
usage = "lagrange_radii -- writes table with Lagrange radii\n";
//------------------------------------------------------------------------------
inline std::ostream& open_ofstream(const char* file)
{
  std::ofstream *out = new std::ofstream();
  if(! open_to_append(*out,file) ) nbdy::exit(1);
  return *(static_cast<std::ostream*>(out));
}
//------------------------------------------------------------------------------
void nbdy::main()
{
  // 1. set some general parameters                                             
  const int Nmax = 100;
  double    M[Nmax] = {0.}, R[Nmax];
  double    T =-1., dT=getdparam("step");
  int       N = nemoinpd(getparam("masses"),M,Nmax);
  if(N>Nmax) error("too many mass shells specified");
  io        need=io::mx, read;
  const io  give = getioparam_z("give");
  bodies    BB;
  nemo_in in(getparam("in"));                      // open nemo input stream    
  nemo_out*snp  = hasvalue("out")?                 // open nemo output stream   
    new nemo_out(getparam("out")) : 0;             //   for primary output      
  std::ostream &out = hasvalue("tabfile")?
    open_ofstream(getparam("tabfile")) : std::cout;
  out<<"# lagrange_rad in="<<getparam("in");
  if(hasvalue("tabfile"))    out<<" tabfile="<<getparam("tabfile");
  out<<":\n# time       ";
  for(int j=0; j!=N; ++j) out<<"r["<<std::setw(4)<<100*M[j]<<"%] ";
  out<<std::endl;
  out.setf(std::ios::left, std::ios::adjustfield);
  // 2. set stop parameters                                                     
  bool       stopping = hasvalue   ("stopfile");
  const int  stop_i   = getiparam  ("stopindex");
  const real stop_val = getdparam_z("stopvalue");
  const bool stop_rel = getbparam  ("stoprelative");
  real       stop_aft = getdparam_z("stopafter");
  bool have_stop_after= hasvalue   ("stopafter");
  const real stop_del = getdparam  ("stopdelay");
  if(hasvalue ("stopfile") && (stop_i < 0 || stop_i >= N)) {
    warning("will not use stop condition, since stopindex not in [0,%d]",N-1);
    stopping = false;
  }
  if(hasvalue ("stopfile") && stop_val == zero) {
    warning("will not use stop condition, since stopvalue=0");
    stopping = false;
  }
  if(hasvalue ("stopfile") && 
     stop_rel && (stop_val>one || (stop_val<zero && stop_val>-one)))
    warning("stop condition already true at initial time");
  bool stopped = false, has_Ri0=false;
  register double Ri,Ri0;
  real stop_time, last_time;
  // 3. loop snapshots                                                          
  bool   in_time_interval;
  real   Mt;
  double t;
  while(in.is_present(nemo_io::snap)) {            // WHILE(snapshots)    >     
    // 3.1 read data, set stop_after if needed, skip if not in time interval    
    in_time_interval = BB.read_nemo_snapshot(in,read,&t,need|give,
					     getparam("times"),0);
    if(stopping && !have_stop_after) {
      stop_aft        = t;
      have_stop_after = true;
    }
    if(!in_time_interval) continue;
    // 3.2 check for snapshot output                                            
    if(T==-1.) T=t-1.e-6*dT;
    if(snp && t >= T) {                            //   IF(time for output) >   
      BB.write_nemo_snapshot(*snp,&t,read);        //     output snapshot       
      T += dT;                                     //     time for next O       
    }                                              //   <                       
    if(!read.contains(need)) {                     //   too few data: abort     
      char w_read[30], w_need[30];
      read.make_word(w_read);
      need.make_word(w_need);
      error("insufficient data: got only %s, but need %s",w_read,w_need);
    }
    // 3.3 compute lagrange radii and write them out                            
    out<<std::setw(12)<<t<<" ";
    find_lagrange_rad(&BB,N,M,R);
    for(int i=0; i!=N; ++i) out << std::setw(8) << R[i] << " ";
    out<<std::endl;
    Ri = R[stop_i];
    // 3.4 deal with stopping of simulation                                     
    if(stopping) {
      if(!has_Ri0) { Ri0 = Ri; has_Ri0 = true; }
      if(!stopped) {
	register real stop_x = stop_rel? Ri/Ri0 : Ri;
	if((stop_val > zero && stop_x <     stop_val  ) ||
	   (stop_val < zero && stop_x > abs(stop_val) ) ) {
	  stopped   = true;
	  stop_time = t;
	  last_time = max(t+getdparam("stopdelay"),getdparam("stopafter"));
	}
      }
      if(stopped && t >= last_time) {
	std::ofstream out(getparam("stopfile"));
	out<<" stopfile \""<<getparam("stopfile")<<"\"\n"
	   <<" generated by program \"lagrange_rad\"\n"
	   <<" on simulation time "<<t<<", because the stop condition\n"
	   <<"\n   R("<<M[stop_i]<<"M_tot)";
	if(stop_rel) out<<"/R("<<M[stop_i]<<"M_tot, t=t_initial)";
	out<< ((stop_val<zero)? " > " : " < ")
	   << abs(stop_val) <<"\n\n"
	   <<" was satisfied at simulation time "<<stop_time<<"\n";
      }
    }
  }
}
//------------------------------------------------------------------------------
