// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// tools.h                                                                     |
//                                                                             |
// Copyright (C) 2002-2005  Walter Dehnen                                      |
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
#ifndef falcON_included_tools_h
#define falcON_included_tools_h 1

#ifndef falcON_included_body_h
#  include <body.h>
#endif
#ifndef falcON_included_tree_h
#  include <public/tree.h>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  //----------------------------------------------------------------------------
  // compute centre of mass of bodies                                           
  //----------------------------------------------------------------------------
  vect centre_of_mass(const bodies*const&);        // I: bodies                 
  //----------------------------------------------------------------------------
  inline                                           // american spelling         
  vect center_of_mass(const bodies*const&B)
  {
    return centre_of_mass(B);
  }
  //----------------------------------------------------------------------------
  // compute iteratively the |Phi|^alpha weighted centre                        
  //                                                                            
  //   The centre is found iteratively as the m*|phi|^alpha weighted position   
  //   of all bodies that are within r_cut from the last centre position. At    
  //   each iteration r_cut is reduce by the factor f. The iteration is stopped 
  //   if the number of bodies in the final sphere is the smallest greater than 
  //   Nmin.                                                                    
  //   If the bool argument is true, we use the input values for centre position
  //   and radius as starting point for the iteration. Otherwise, we use the    
  //   mean position and mean-square radius of the whole distribution.          
  //----------------------------------------------------------------------------
  void find_centre(const bodies*const&,            // I  : bodies               
		   real         const&,            // I  : reduction factor f   
		   unsigned     const&,            // I  : alpha                
		   unsigned     const&,            // I  : Nmin                 
		   vect              &,            // I/O: centre position      
		   real              &,            // I/O: centre radius        
		   bool                = 0,        //[I  : use input as initial 
		   vect              * = 0,        //[O  : centre velocity]     
		   real              * = 0,        //[O  : centre density]      
		   unsigned            = 0,        //[I  : first body]          
		   unsigned            = 0)        //[I  : # bodies]            
    falcON_THROWING;
  //----------------------------------------------------------------------------
  // estimate position of density maximum weighted by w = mass * |pot|^alpha    
  //     given an OctTree which holds the bodies we                             
  //     - find the cell with the maximum weight density                        
  //     - compute its weighted center                                          
  //     - estimate h (center radius) to yield Nmin bodies                      
  //----------------------------------------------------------------------------
  void estimate_density_peak(OctTree *const&,      // I  : OctTree              
			     unsigned const&,      // I  : alpha                
			     unsigned const&,      // I  : Nmin                 
			     vect          &,      // O  : center position      
			     real          &);     // O  : center radius        
  //----------------------------------------------------------------------------
  // find the lagrange radii around the origin for a given set of cumulative    
  // masses (relative to total).                                                
  //----------------------------------------------------------------------------
  void find_lagrange_rad(body  const &,            // I: first body to take     
			       int    ,            // I: size of tables below   
			 const double*,            // I: table: masses          
			       double*,            // O: table: lagrange radii  
			 const vect  * =0,         //[I: centre offset]         
			 unsigned      =0);        //[I: # bodies, def: all]    
#ifdef falcON_PROPER
  //----------------------------------------------------------------------------
  // compute iteratively the density centre                                     
  // defined as the position X0 where                                           
  //     Sum_i Mi * K(|X0-Xi|/R)                                                
  // becomes maximal, where R is the centre radius, determined such that        
  //     Ncen = bodies with |X0-Xi| < R                                         
  //                                                                            
  //----------------------------------------------------------------------------
  void find_centre(const bodies*const&,            // I  : bodies               
		   unsigned     const&,            // I  : Ncen                 
		   vect              &,            // I/O: centre position      
		   real              &,            // I/O: centre radius        
		   vect              * = 0,        //[O  : centre velocity]     
		   real              * = 0);       //[O  : centre density]      
  //----------------------------------------------------------------------------
  // estimate phase-space density using an algorithm very similar to that of    
  // Ascasibar & Binney (2005)                                                  
  //----------------------------------------------------------------------------
  void EstimateVol(const snapshot*,                // I: snapshot               
		   real*    const&,                // O: phase-space volumes    
		   int, int, int,                  // I: me,wx,wv               
		   real* = 0,                      //[I: alpha,theta,phi]       
		   bool  = false,                  //[I: add (or assign) vol?]  
		   bool  = true    );              //[I: boundary correction?]  
  //----------------------------------------------------------------------------
  inline
  void EstimatePhD(const snapshot*S,               // I: snapshot               
		   real*    const&F,               // O: phase-space densities  
		   int me, int wx, int wv,
		   bool           C=true) {        // I: boundary correction?   
    EstimateVol(S,F,me,wx,wv,0,0,C);
    LoopAllBodies(S,b)
      F[bodyindex(b)] = mass(b) / F[bodyindex(b)];
  }
#endif
  //----------------------------------------------------------------------------
} // namespace {
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_tool_h
