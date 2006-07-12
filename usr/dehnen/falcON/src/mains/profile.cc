// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// profile.cc                                                                  |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 2006                                               |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// This is a non-public part of the code.                                      |
// It is property of its author and not to be made public without his written  |
// consent.                                                                    |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 0.0    26/04/2006 WD created                                              |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "0.0"
#define falcON_VERSION_D "26-apr-2006 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                   // this is a NEMO program 
#  error You need NEMO to compile "profile"
#endif
#define falcON_RepAction 0                            // no action reporting    
//-----------------------------------------------------------------------------+
#include <fstream>                                    // C++ file I/O           
#include <iomanip>                                    // C++ I/O manipulators   
#include <body.h>                                     // the bodies             
#include <public/io.h>                                // my NEMO I/O            
#include <main.h>                                     // main & NEMO stuff      
using namespace falcON;
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n           input file                                         ",
  "file=-\n           file for output (may contain format string)        ",
  "times=all\n        times to process                                   ",
  "K=16\n             # bodies in nearest-neighbour density estimate     ",
  "Nwin=1000\n        # bodies per bin                                   ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string
usage = "profile -- writes table with averaged quantities\n";
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  if(! file_for_output("file") )
    return warning("no output to be made\n");
  const nemo_in  IN(getparam("in"));
  const fieldset WANT(fieldset(fieldset::basic));
  output         TAB;
  fieldset       READ;
  snapshot       SHOT;
  int            I=0;
  // loop snapshots
  while(IN.has_snapshot()) {
    if(! SHOT.read_nemo(IN,READ,WANT,getparam("times"),0) ) continue;
    if(  warn_insufficient(READ,WANT) ) continue;
    // 1. estimate density

    // 2. create table with bodies sorted in descending density

    // 3. open output file

    // 4. loops bins in descending density

    // 4.1  compute centre and mean density
    
    // 4.2  compute moment of intertia,
    //              angular momentum,
    //              mean rotation and mean radial motion

    // 4.3  compute velocity dispersions (r,theta,phi), correlations

  }
}
