// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    src/mains/density.cc                                               
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    2006                                                               
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2006  Walter Dehnen                                            
//                                                                              
// This program is free software; you can redistribute it and/or modify         
// it under the terms of the GNU General Public License as published by         
// the Free Software Foundation; either version 2 of the License, or (at        
// your option) any later version.                                              
//                                                                              
// This program is distributed in the hope that it will be useful, but          
// WITHOUT ANY WARRANTY; without even the implied warranty of                   
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
// General Public License for more details.                                     
//                                                                              
// You should have received a copy of the GNU General Public License            
// along with this program; if not, write to the Free Software                  
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                    
//                                                                              
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// history:                                                                     
//                                                                              
// v 0.0    04/04/2006  WD created as TestNeigh                                 
// v 0.1    05/04/2006  WD implemented density estimation                       
// v 0.2    06/04/2006  WD added choose & params                                
// v 1.0    27/04/2006  WD default K=16                                         
// v 1.1    02/05/2006  WD output is opened only when needed                    
// v 1.2    26/07/2006  WD BodyFilter                                           
// v 2.0    27/07/2006  WD made public                                          
// v 2.1    28/07/2006  WD keyword findmax added                                
////////////////////////////////////////////////////////////////////////////////
#define falcON_VERSION   "2.0"
#define falcON_VERSION_D "27-jul-2006 Walter Dehnen                          "
//-----------------------------------------------------------------------------+
#ifndef falcON_NEMO                                // this is a NEMO program    
#  error You need NEMO to compile "TestNeigh"
#endif
#define falcON_RepAction 0                         // no action reporting       
//-----------------------------------------------------------------------------+
#include <iostream>                                // C++ I/O                   
#include <fstream>                                 // C++ file I/O              
#include <iomanip>                                 // C++ I/O formatting        
#include <ctime>
#include <body.h>                                  // the bodies                
#include <public/io.h>                             // my NEMO I/O               
#include <public/neighbours.h>                     // getting dist to neighbours
#include <public/bodyfunc.h>                       // body functions            
#include <main.h>                                  // main & NEMO stuff         
//------------------------------------------------------------------------------
string defv[] = {
  "in=???\n           input file                                         ",
  "out=???\n          output file                                        ",
  "K=16\n             number of neighbours                               ",
  "times=all\n        times to process                                   ",
  "choose=\n          a boolean bodyfunc expression to choose bodies     ",
  "params=\n          any parameters required by choose                  ",
  "give=\n            output only these, but at least rho (r)            ",
  "findmax=f\n        write density maximum to stderr                    ",
  "verbose=f\n        write some blurb about CPU time etc                ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
string
usage = "estimate mass density based on distance to Kth neighbour";
//------------------------------------------------------------------------------
namespace {
  using namespace falcON;
  falcON::real  FAC;        ///< factor used in estimation of density
  falcON::real  RHO = zero; ///< maximum density 
  bodies::index BRH;        ///< body with maximum density

  void dens(const bodies*                    B,
	    const NeighbourLister::Leaf*     L,
	    const NeighbourLister::Neighbour*N,
	    int                              K)
  {
    falcON::real m(zero);
    const NeighbourLister::Neighbour*NK = N+K;
    for(const NeighbourLister::Neighbour*n=N; n!=NK; ++n)
      m += B->mass(mybody(n->L));
    m *= FAC/cube(max_dist(L));
    B->rho(mybody(L)) = m;
    if(m > RHO) {
      RHO = m;
      BRH = mybody(L);
    }
  }
}
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING
{
  // deal with I/O                                                              
  nemo_in        IN  (getparam("in"));
  nemo_out       OUT;
  BodyFilter     BF  (getparam_z("choose"), getparam_z("params"));
  const fieldset SRCE(fieldset::m | fieldset::x);
  const fieldset GIVE(getioparam_z("give") | fieldset::r);
  const fieldset WANT(getioparam_z("give") | SRCE |
		      (BF? BF.need() : fieldset(fieldset::o)));
  const unsigned K   (getiparam("K"));
  const bool     VERB(getbparam("verbose"));
  // loop snapshots & process them                                              
  fieldset       READ;
  snapshot       SHOT;
  FAC          = (K-0.5)/(K*FPit);
  while(IN.has_snapshot()) {
    // read time, read snapshot
    if(!SHOT.read_nemo(IN,READ,WANT,getparam("times"),0)) continue;
    // check data availability
    check_sufficient(READ,WANT);
    SHOT.add_field(fieldbit::r);
    BF.set_time(SHOT.time());
    // build tree
    clock_t CPU0 = clock();
    flags   FLAG = flags::empty;
    if(BF) {
      FLAG = flags::marked;
      int Nmark(0);
      SHOT.add_field(fieldbit::f);
      LoopAllBodies(&SHOT,B)
	if(BF(B)) {
	  B.mark();
	  ++Nmark;
	} else
	  B.unmark();
      if(VERB)
	std::cerr<<' '<<Nmark<<" bodies satisfy criterion"<<std::endl;
      if(Nmark == 0) continue;
    }
    OctTree TREE(&SHOT, K+1, 0, Default::MaxDepth, FLAG);
    if(VERB) {
      clock_t CPU1 = clock();
      std::cerr<<setprecision(6)
	       <<" time needed for tree build:         "
	       <<(CPU1 - CPU0)/real(CLOCKS_PER_SEC)<<std::endl;
      CPU0 = CPU1;
    }
    // find Kth nearest neighbours and estimate density
    NeighbourLister NELI(&TREE,K);
    NELI.Estimate(&dens,true);
    // for non-chosen bodies set density to zero
    if(BF)
      LoopAllBodies(&SHOT,B)
	if(!is_marked(B)) B.rho() = zero;
    if(VERB) {
      clock_t CPU1 = clock();
      std::cerr<<setprecision(6)
	       <<" time needed for density estimation: "
	       <<(CPU1 - CPU0)/real(CLOCKS_PER_SEC)<<'\n'
	       <<" number of neighbour updates:        "
	       <<NELI.N_interact()<<std::endl;
    }
    // write output
    if(OUT || OUT.open(getparam("out"))) SHOT.write_nemo(OUT,GIVE);
    // find maximum density if wanted
    if(getbparam("findmax"))
      std::cerr<<" maximum density "<< RHO
	       <<" found for body #"<< SHOT.bodyindex(BRH)
	       <<" at x = "<< SHOT.pos(BRH)
	       <<"\n";
  }
}
