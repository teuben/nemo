// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// TestNeigh.cc                                                                |
//                                                                             |
// Copyright (C) 2006 Walter Dehnen                                            |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// history:                                                                    |
//                                                                             |
// v 0.0    04/04/2006  WD created as TestNeigh                                |
// v 1.0    05/04/2006  WD implemented density estimation                      |
//-----------------------------------------------------------------------------+
#define falcON_VERSION   "0.1"
#define falcON_VERSION_D "04-apr-2006 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile "TestNeigh"
#endif
#define falcON_RepAction 0                         // no action reporting       
#define falcON_NOT_USING_PROPER                    // not using PROPRIETARY code
//-----------------------------------------------------------------------------+
#include <iostream>                                // C++ I/O                   
#include <fstream>                                 // C++ file I/O              
#include <iomanip>                                 // C++ I/O formatting        
#include <ctime>
#include <body.h>                                  // the bodies                
#include <public/io.h>                             // my NEMO I/O               
#include <proper/neighbours.h>                     // getting dist to neighbours
#include <main.h>                                  // main & NEMO stuff         
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n           input file                                         ",
  "out=???\n          output file                                        ",
  "K=10\n             number of neighbours                               ",
  "times=all\n        times to process                                   ",
  "give=\n            output only these, but at least rho (r)            ",
  "verbose=f\n        write some blurb about CPU time etc                ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string
usage = "estimates density based on distance to Kth neighbour\n";
//------------------------------------------------------------------------------
namespace {
  using namespace falcON;
  snapshot SHOT;
  real     FAC;
  void dens(const NeighbourLister::Leaf*     L,
	    const NeighbourLister::Neighbour*N,
	    int K)
  {
    real m(zero);
    const NeighbourLister::Neighbour*NK = N+K;
    for(NeighbourLister::Neighbour*n=N; n!=NK; ++n)
      m += SHOT.mass(mybody(n->L));
    SHOT.rho(mybody(L)) = Fac*m/cube(max_dist(L));
  }
}
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  // deal with I/O                                                              
  nemo_in        IN  (getparam("in"));
  nemo_out       OUT (getparam("out"));
  const fieldset GIVE(getioparam_z("give") | fieldset::r);
  const fieldset SRCE(fieldset::m|fieldset::x);
  const fieldset WANT(SRCE | getioparam_z("give"));
  const unsigned K   (getiparam("K"));
  const bool     VERB(getbparam("verbatim"));
  clock_t        CPU0, CPU1;
  // loop snapshots & process them                                              
  fieldset       READ;
  FAC          = (K-0.5)/(K*FPit);
  while(IN.has_snapshot()) {
    // read time, read snapshot if in times; take first simulation time         
    if(!SHOT.read_nemo(IN,READ,WANT,getparam("times"),0))  continue;
    // find Kth nearest neighbours
    check_sufficient(READ,SRCE);
    SHOT.add_field(fieldbit::r);
    CPU0 = clock();
    OctTree         TREE(&SHOT, K+1);
    if(VERB) {
      CPU1 = clock();
      std::cerr<<setprecision(6)
	       <<" time needed for tree build:         "
	       <<(CPU1 - CPU0)/real(CLOCKS_PER_SEC)<<std::endl;
      CPU0 = CPU1;
    }
    NeighbourLister NELI(&TREE,K);
    NELI.Estimate(&dens,true);
    if(VERB) {
      CPU1 = clock();
      std::cerr<<setprecision(6)
	       <<" time needed for density estimation: "
	       <<(CPU1 - CPU0)/real(CLOCKS_PER_SEC)<<'\n'
	       <<" number of neighbour updates:        "
	       <<NELI.N_interact()<<std::endl;
    }
    // write output
    if(OUT) SHOT.write_nemo(OUT,GIVE);
  }
}
