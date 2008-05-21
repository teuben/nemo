// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    inc/public/neighbours.h                                            
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    2008                                                               
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2008  Walter Dehnen                                            
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
// \version 31/08/2007 WD initial hack (based on code from code/dens)           
// \version 04/09/2007 WD debugged and working version                          
// \version 04/09/2007 WD moved to neighbours.h (which is superseeded)          
// \version 20/05/2008 WD added FindLeaf(), ProcessNeighbours()                 
//                                                                              
////////////////////////////////////////////////////////////////////////////////
#ifndef falcON_included_neighbours_h
#define falcON_included_neighbours_h

#ifndef falcON_included_tree_h
#  include <public/tree.h>
#endif
// /////////////////////////////////////////////////////////////////////////////
namespace falcON {
  //----------------------------------------------------------------------------
  /// represents entry in a neighbour list
  struct Neighbour {
    real Q;                  ///< distance^2 to leaf below
    const OctTree::Leaf*L;   ///< neighbouring leaf
  };
  //----------------------------------------------------------------------------
  /// For each active body's leaf: find list of neighbours and process it.
  ///
  /// We first copy the body mass (to OctTree::Leaf::scalar()), and, if
  /// all=false, also the body flags to leafs.\n
  /// We then find for each active leaf the K nearest neighbours in ascending
  /// order of distance.\n
  /// This list is then processed by a user-supplied function.\n
  /// Used in density estimation.
  /// \note The user function's 4th argument is the size K of the list.
  /// \note The user function can assume scalar(Leaf) to hold the mass.
  /// \note The tree must not be have been re-used, but grown.
  /// \param T tree to use
  /// \param K number of neighbours to find
  /// \param f function for processing neighbour list
  /// \param Ni (output) number of leaf interactions
  /// \param all (optional) consider all bodies to be active?
  void ProcessNearestNeighbours(const OctTree*T, int K,
				void(*f)(const bodies*, const OctTree::Leaf*,
					 const Neighbour*, int),
				unsigned&Ni, bool all=0) falcON_THROWING;
  //----------------------------------------------------------------------------
  /// given a body, find smallest surrounding cell
  /// \return (pter to) smallest cell containing body
  /// \param[in]  T tree to search
  /// \param[in]  B body to find leaf for
  /// \note The tree must not be have been re-used, but grown.
  const OctTree::Cell*FindCell(const OctTree*T, body const&B);
  //----------------------------------------------------------------------------
  /// Process all leafs within a certain distance from body B
  /// \param[in] T tree to use
  /// \param[in] B body to process neighbours for
  /// \param[in] Q radius^2 of sphere to look for neighbours in
  /// \param[in] f function to call for each neighbour found
  /// \note The leaf, if any, corresponding to B is not processed
  /// \note The tree must not be have been re-used, but grown.
  /// \note The 3rd argument of f() is the distance to the leaf
  void ProcessNeighbours(const OctTree*T, body const&B, real Q,
			 void(*f)(const bodies*, const OctTree::Leaf*, real q))
    falcON_THROWING;
} // namespace falcON
////////////////////////////////////////////////////////////////////////////////
falcON_TRAITS(falcON::Neighbour,"Neighbour");
// /////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_neighbours_h
