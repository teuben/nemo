// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// splitstream.cc                                                              |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2002-2003                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 1.0     09/10/2003  WD created                                            |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "1.0"
#define falcON_VERSION_D "10-oct-2003 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need "NEMO" to compile src/mains/gyrfalcON.cc
#endif

#include <cstdio>                                  // C I/O                     
#include <body.h>                                  // bodies                    
#include <public/nmio.h>                           // nemo I/O                  
#include <main.h>                                  // main & NEMO stuff         
using namespace nbdy;
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n           input file, not \"-\" or \".\"                    ",
  "out=out%03d.snp\n  format for output file name        ",
  "i0=0\n             index of first output file name                   ",
  "times=all\n        times to process                              ",
  "copy=mxv\n         list of items to copy. Recognizing:\n"
  "                    m: mass                                (default)\n"
  "                    x: position                            (default)\n"
  "                    v: velocity                            (default)\n"
  "                    a: acceleration\n"
  "                    e: individual eps_i\n"
  "                    p: N-body potential\n"
  "                    l: time-step level\n"
  "                    f: body flag\n"
  "                    r: density estimate                              ",
  version, compiled, NULL };
//------------------------------------------------------------------------------
string 
usage = "streamsplit -- splits a snapshot file into one file per snapshot";
//------------------------------------------------------------------------------
void nbdy::main()
{
  nemo_in   input(getparam("in"));
  nemo_out  output;
  real      time;
  io        copy = io(getparam("copy")), got;
  bodies    bs;
  bool      is_in_times;
  char      out[200];
  int       index = getiparam("i0");
  while(input.is_present(nemo_io::snap)) 
    if(bs.read_nemo_snapshot(input,got,&time,copy,getparam("times"),0)) {
      if(!got.contains(copy)) 
	warning("couldn't read body data: %s",word(got.missing(copy)));
      sprintf(out,getparam("out"),index++);
      output.open(out);
      bs.write_nemo_snapshot(output,&time,copy&got);
      output.close();
    }
}
//------------------------------------------------------------------------------
