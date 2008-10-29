// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   inc/public/tools.h                                                  
///                                                                             
/// \author Walter Dehnen                                                       
/// \date   2002-2006                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2002-2006 Walter Dehnen                                        
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
  //============================================================================
  //                                                                            
  /// \name some tools for computing properties of subset of bodies             
  //@{                                                                          
  //----------------------------------------------------------------------------
  /// compute centre of mass for set of all bodies in_subset()                  
  /// \param B  (input) a set of bodies                                         
  /// \return centre of mass of all bodies in set                               
  vect centre_of_mass(const bodies*B);
  //----------------------------------------------------------------------------
  /// compute iteratively the |Phi|^alpha weighted centre of all bodies         
  /// in_subset()                                                               
  ///                                                                           
  /// \param B      (input) bodies to consider                                  
  /// \param f      (input) reduction factor                                    
  /// \param alpha  (input) power for weighting                                 
  /// \param Nmin   (input) stop iteration before # bodies drops below this     
  /// \param xc     (in/output) centre position (initial guess/improved value)  
  /// \param rc     (in/output) centre radius (initial guess/improved value)    
  /// \param ini    (input, optional): use x0,r0 on input as initial guess?     
  /// \param vc     (output, optional): return centre velocity if non-null      
  /// \param rhc    (output, optional): return centre density if non-null       
  ///                                                                           
  /// The centre is found iteratively as the m*|phi|^alpha weighted position    
  /// of all bodies that are within r_cut from the last centre position. At     
  /// each iteration r_cut is reduced by the factor \a f. The iteration is      
  /// stopped if the number of bodies in the final sphere is the smallest       
  /// greater than Nmin. If \a use is true, we use the input values for centre  
  /// position and radius as starting point for the iteration. Otherwise, we    
  /// use the mean position and mean-square radius of the whole distribution.   
  void find_centre_alpha(const bodies*B,
			 real         f,
			 unsigned     alpha,
			 unsigned     Nmin,
			 vect        &xc,
			 real        &rc,
			 bool         ini = 0,
			 vect        *vc  = 0,
			 real        *rhc = 0) falcON_THROWING;
  //----------------------------------------------------------------------------
  /// estimate position of density maximum weighted by w = mass * |pot|^alpha   
  ///                                                                           
  /// \param Tree  (input)  OctTree holding the relevant bodies                 
  /// \param alpha (input)  power for weighting                                 
  /// \param Nmin  (input)  min # bodies in centre                              
  /// \param xc    (output) estimated density centre position                   
  /// \param rc    (output) estimated radius holding \a Nmin bodies             
  ///                                                                           
  /// The centre \a xc is taken to be the mass * |pot|^alpha weighted position  
  /// of all bodies inside the tree cell with the maximum average weight        
  /// density.                                                                  
  void estimate_density_peak(OctTree *Tree,
			     unsigned alpha,
			     unsigned Nmin,
			     vect    &xc,
			     real    &rc);
  //----------------------------------------------------------------------------
  /// find the lagrange radii around the origin for a given set of cumulative
  /// masses (relative to total) for all bodies in_subset().
  ///                                                                           
  /// \param B[in]      set of bodies to consider
  /// \param n[in]      size of tables for masses and Lagrange radii
  /// \param M[in]      table with relative masses
  /// \param R[out]     table with Lagrange radi
  /// \param off[in]    (optional) offset of centre from origin
  void find_lagrange_rad(const bodies*B,
			 unsigned     n,
			 const double*M,
			 double      *R,
			 const vect  *off = 0);
  //----------------------------------------------------------------------------
  /// compute iteratively the density centre of all bodies in_subset().         
  /// \return       true if density centre has been found (convergence)         
  /// \param B      (input) set of bodies to consider                           
  /// \param Ncen   (input) # bodies in centre                                  
  /// \param xc     (in/output) centre position                                 
  /// \param rc     (in/output) centre radius                                   
  /// \param vc     (output, optional) centre velocity                          
  /// \param rhc    (putput, optional) centre density                           
  ///                                                                           
  /// The density centre is found as the position X0 where                      
  ///     Sum_i Mi * K(|X0-Xi|/R)                                               
  /// becomes maximal, where R is the centre radius, determined such that       
  ///     Ncen = bodies with |X0-Xi| < R                                        
  /// We find this position iteratively starting from the initial values for xc 
  /// and rhc.                                                                  
  bool find_density_centre(const bodies*B,
			   unsigned     Ncen,
			   vect        &xc,
			   real        &rc,
			   vect        *vc  = 0,
			   real        *rhc = 0);
#ifdef falcON_PROPER
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
  //@}
  //----------------------------------------------------------------------------
} // namespace {
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_tool_h
