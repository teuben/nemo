// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// density_centre.cc                                                           |
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
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "2.3"
#define falcON_VERSION_D "12-may-2004 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile "src/mains/density_centre.cc"
#endif
#define falcON_RepAction 0                         // no action reporting       
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
  "out=\n             output file                                        ",
  "times=all\n        times to process                                   ",
  "Ncen=200\n         # bodies in center                                 ",
  "centrefile=-\n     file to write center position & velocity to        ",
  "centerfile=\n      same as centrefile                                 ",
  "give=\n            what to give to out                                ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string
usage = "density_centre -- determine density centre & centers snapshots\n"
        "The center is defined as the position of the density maximum\n"
        "If that is not well-defined, the routine may fail\n";
//------------------------------------------------------------------------------
void nbdy::main()
{
  // get parameters                                                             
  const bool     SOUT(hasvalue("out") && strcmp(getparam("out"),"."));
  const unsigned NCEN(getiparam("Ncen"));
  const nemo_in  IN  (getparam("in"));
  const nemo_out OUT (SOUT? getparam("out") : ".");
  const io       GIVE(getioparam_z("give"));
  const io       NEED(io::mxv | GIVE);
  // open output for center position                                            
  std::ostream  *COUT(0);
  std::ofstream  CFILE;
  int stdout=0;
  if(hasvalue("out") && 0==strcmp(getparam("out"),"-")) stdout+= 1;
  if(0 == strcmp(hasvalue("centerfile")? 
		 getparam("centerfile") : 
		 getparam("centrefile"), "-")) {
    COUT   = &(std::cout);
    stdout+= 1;
  } else if( 0 != strcmp(hasvalue("centerfile")? 
			 getparam("centerfile") : 
			 getparam("centrefile"), ".")) {
    open_to_append(CFILE,hasvalue("centerfile")?
		   getparam("centerfile") : getparam("centrefile"));
    COUT  = &CFILE;
  }
  if(stdout > 1) error("more than one output to \"-\", ie. stdout\n");
  if(COUT)
    (*COUT) << "#\n"
	    << "# \""<< (*(ask_history())) <<"\"\n"
	    << "#\n"
	    << "# t x y z vx vy vz rho"<<std::endl;
  // loop snapshots & process them                                              
  bodies   BODIES;
  int      NOLD = 0u;
  oct_tree*TREE = 0;
  io       READ;
  real     RHO,RAD(zero);
  double   TIME;
  vect     XCEN(zero),VCEN(zero);
  while(IN.is_present(nemo_io::snap)) {
    if(! BODIES.read_nemo_snapshot(IN,READ,&TIME,NEED,getparam("times"),0))
      continue;
    if(! READ.contains(NEED)) error("insufficient data");
    if(TREE==0 || NOLD != BODIES.N_bodies()) {
      if(TREE) delete TREE;
      TREE = new oct_tree(&BODIES,NCEN/4);
    } else {
      TREE->build(NCEN/4);
    }
    estimate_density_peak(TREE,0u,NCEN,XCEN,RAD);
    RAD *= 3;
    find_centre(&BODIES,NCEN,XCEN,RAD,&VCEN,&RHO);
    if(COUT)
      (*COUT)<<"  "
	     <<std::setw(15)<<std::setprecision(8)<<TIME   <<"  "
	     <<std::setw(15)<<std::setprecision(8)<<XCEN[0]<<" "
	     <<std::setw(15)<<std::setprecision(8)<<XCEN[1]<<" "
	     <<std::setw(15)<<std::setprecision(8)<<XCEN[2]<<"  "
	     <<std::setw(15)<<std::setprecision(8)<<VCEN[0]<<" "
	     <<std::setw(15)<<std::setprecision(8)<<VCEN[1]<<" "
	     <<std::setw(15)<<std::setprecision(8)<<VCEN[2]<<"  "
	     <<std::setw(15)<<std::setprecision(8)<<RHO    <<std::endl;
    if(SOUT) {
      LoopBodies(bodies,&BODIES,Bi) {
	Bi.pos() -= XCEN;
	Bi.vel() -= VCEN;
      }
      BODIES.write_nemo_snapshot(OUT,&TIME,READ);
    }
  }
  if(TREE) delete TREE;
}
