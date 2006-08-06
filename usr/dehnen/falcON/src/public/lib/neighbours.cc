// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    src/public/neighbours.cc                                           
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
// \version 27/07/2006 WD made public                                           
//                                                                              
////////////////////////////////////////////////////////////////////////////////
#include <public/neighbours.h>
#include <public/interact.h>

using namespace falcON;
namespace {
  NeighbourLister::cell_iterator 
  parent(NeighbourLister::cell_iterator const&C) {
    return C.CellNo(fcparent(C));
  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class NeighbourFinderBase                                                  
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  class NeighbourFinderBase {
  public:
    typedef NeighbourLister::leaf_iterator leaf_iter;
    typedef NeighbourLister::cell_iterator cell_iter;
    //--------------------------------------------------------------------------
  protected:
    const unsigned  K;   ///< number of neighbours listed
    mutable uint64  I;   ///< counter: number of neighbour updates
    //--------------------------------------------------------------------------
    NeighbourFinderBase(int k) : K(k), I(0u) {}
    //--------------------------------------------------------------------------
    void many_YA(leaf_iter const&A,
		 leaf_iter const&B0,
		 leaf_iter const&BN,
		 unsigned       &iA,
		 unsigned       &iB) const {
      const vect xA(pos(A));
      for(register leaf_iter B=B0; B!=BN; ++B) {
	real Rq = dist_sq(xA,pos(B));
	if(Rq < max_dist_sq(A))
	  if(A->replace_furthest(B,Rq,K)) ++iA;
	if(Rq < max_dist_sq(B))
	  if(B->replace_furthest(A,Rq,K)) ++iB;
      }
    }
    //--------------------------------------------------------------------------
    void many_YS(leaf_iter const&A,
		 leaf_iter const&B0,
		 leaf_iter const&BN,
		 unsigned       &iA,
		 unsigned       &iB) const {
      const vect xA(pos(A));
      for(register leaf_iter B=B0; B!=BN; ++B) {
	real Rq = dist_sq(xA,pos(B));
	if(Rq < max_dist_sq(A))
	  if(A->replace_furthest(B,Rq,K)) ++iA;
	if(is_active(B) && Rq < max_dist_sq(B))
	  if(B->replace_furthest(A,Rq,K)) ++iB;
      }
    }
    //--------------------------------------------------------------------------
    void many_YN(leaf_iter const&A,
		 leaf_iter const&B0,
		 leaf_iter const&BN,
		 unsigned       &iA,
		 unsigned       &iB) const {
      const vect xA(pos(A));
      for(register leaf_iter B=B0; B!=BN; ++B) {
	real Rq = dist_sq(xA,pos(B));
	if(Rq < max_dist_sq(A))
	  if(A->replace_furthest(B,Rq,K)) ++iA;
      }
    }
    //--------------------------------------------------------------------------
    void many_NA(leaf_iter const&A,
		 leaf_iter const&B0,
		 leaf_iter const&BN,
		 unsigned       &iA,
		 unsigned       &iB) const {
      const vect xA(pos(A));
      for(register leaf_iter B=B0; B!=BN; ++B) {
	real Rq = dist_sq(xA,pos(B));
	if(Rq < max_dist_sq(B))
	  if(B->replace_furthest(A,Rq,K)) ++iB;
      }
    }
    //--------------------------------------------------------------------------
    void many_NS(leaf_iter const&A,
		 leaf_iter const&B0,
		 leaf_iter const&BN,
		 unsigned       &iA,
		 unsigned       &iB) const {
      const vect xA(pos(A));
      for(register leaf_iter B=B0; B!=BN; ++B) if(is_active(B)) {
	real Rq = dist_sq(xA,pos(B));
	if(Rq < max_dist_sq(B))
	  if(B->replace_furthest(A,Rq,K)) ++iB;
      }
    }
    //--------------------------------------------------------------------------
  public:
    //--------------------------------------------------------------------------
    uint64 const&N_interact() const { return I; }
    //--------------------------------------------------------------------------
    bool split_first(cell_iter const&A, cell_iter const&B) const {
      if(is_twig(A)) return false;
      if(is_twig(B)) return true;
      return number(A) > number(B);
    }
    //--------------------------------------------------------------------------
    static bool take(cell_iter const&) { return true; }
    static bool take(leaf_iter const&) { return true; }
  };
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class NeighbourFinder<AL_ACTIVE>                                           
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  template<bool> class NeighbourFinder;
  template<> class NeighbourFinder<false> : public NeighbourFinderBase {
    void many(cell_iter const&A,
	      cell_iter const&B,
	      unsigned       &iA,
	      unsigned       &iB) const {
      if(al_active(A)) {
	if(al_active(B)) {
	  LoopAllLeafs(cell_iter,A,Ai) {
	    real Rq = dist_sq(pos(Ai),pos(B));
	    bool uA = Rq < square(max_dist(A)+size(B));
	    bool uB = Rq < square(max_dist(B));
	    if(uA || uB)
	      many_YA(Ai, B.begin_leafs(), B.end_leaf_desc(), iA,iB);
	  }
	} else if(is_active(B)) {
	  LoopAllLeafs(cell_iter,A,Ai) {
	    real Rq = dist_sq(pos(Ai),pos(B));
	    bool uA = Rq < square(max_dist(A)+size(B));
	    bool uB = Rq < square(max_dist(B));
	    if(uA || uB)
	      many_YS(Ai, B.begin_leafs(), B.end_leaf_desc(), iA,iB);
	  }
	} else {
	  LoopAllLeafs(cell_iter,A,Ai) {
	    real Rq = dist_sq(pos(Ai),pos(B));
	    bool uA = Rq < square(max_dist(A)+size(B));
	    if(uA)
	      many_YN(Ai, B.begin_leafs(), B.end_leaf_desc(), iA,iB);
	  }
	}
      } else if(is_active(A)) {
	if(al_active(B)) {
	  LoopAllLeafs(cell_iter,A,Ai) 
	    if(is_active(Ai)) {
	      real Rq = dist_sq(pos(Ai),pos(B));
	      bool uA = Rq < square(max_dist(A)+size(B));
	      bool uB = Rq < square(max_dist(B));
	      if(uA || uB)
		many_YA(Ai, B.begin_leafs(), B.end_leaf_desc(), iA,iB);
	    } else {
	      real Rq = dist_sq(pos(Ai),pos(B));
	      bool uB = Rq < square(max_dist(B));
	      if(uB)
		many_NA(Ai, B.begin_leafs(), B.end_leaf_desc(), iA,iB);
	    }
	} else if(is_active(B)) {
	  LoopAllLeafs(cell_iter,A,Ai) 
	    if(is_active(Ai)) {
	      real Rq = dist_sq(pos(Ai),pos(B));
	      bool uA = Rq < square(max_dist(A)+size(B));
	      bool uB = Rq < square(max_dist(B));
	      if(uA || uB)
		many_YS(Ai, B.begin_leafs(), B.end_leaf_desc(), iA,iB);
	    } else {
	      real Rq = dist_sq(pos(Ai),pos(B));
	      bool uB = Rq < square(max_dist(B));
	      if(uB)
		many_NS(Ai, B.begin_leafs(), B.end_leaf_desc(), iA,iB);
	    }
	} else {
	  LoopAllLeafs(cell_iter,A,Ai) if(is_active(A)) {
	    real Rq = dist_sq(pos(Ai),pos(B));
	    bool uA = Rq < square(max_dist(A)+size(B));
	    if(uA)
	      many_YN(Ai, B.begin_leafs(), B.end_leaf_desc(), iA,iB);
	  }
	}
      } else {
	if(al_active(B)) {
	  LoopAllLeafs(cell_iter,A,Ai) {
	    real Rq = dist_sq(pos(Ai),pos(B));
	    bool uB = Rq < square(max_dist(B));
	    if(uB)
	      many_NA(Ai, B.begin_leafs(), B.end_leaf_desc(), iA,iB);
	  }
	} else if(is_active(B)) {
	  LoopAllLeafs(cell_iter,A,Ai) {
	    real Rq = dist_sq(pos(Ai),pos(B));
	    bool uB = Rq < square(max_dist(B));
	    if(uB)
	      many_NS(Ai, B.begin_leafs(), B.end_leaf_desc(), iA,iB);
	  }
	}
      }
    }
    //--------------------------------------------------------------------------
    void many(cell_iter const&A,
	      unsigned       &iA) const {
      if(al_active(A)) {
	LoopLstLeafs(cell_iter,A,Ai)
	  many_YA(Ai,Ai+1,A.end_leaf_desc(),iA,iA);
      } else if(is_active(A)) {
	LoopLstLeafs(cell_iter,A,Ai)
	  if(is_active(Ai))
	    many_YS(Ai,Ai+1,A.end_leaf_desc(),iA,iA);
	  else
	    many_NS(Ai,Ai+1,A.end_leaf_desc(),iA,iA);
      }
    }
    //--------------------------------------------------------------------------
    void many(leaf_iter const&A,
	      cell_iter const&B,
	      unsigned       &iA,
	      unsigned       &iB) const {
      if(is_active(A))
	return al_active(B)? many_YA(A,B.begin_leafs(),B.end_leaf_desc(),iA,iB)
	  :    is_active(B)? many_YS(A,B.begin_leafs(),B.end_leaf_desc(),iA,iB)
	  :                  many_YN(A,B.begin_leafs(),B.end_leaf_desc(),iA,iB);
      else
	return al_active(B)? many_NA(A,B.begin_leafs(),B.end_leaf_desc(),iA,iB)
	  :    is_active(B)? many_NS(A,B.begin_leafs(),B.end_leaf_desc(),iA,iB)
	  :                  void(0);
    }
    //--------------------------------------------------------------------------
    static bool update_max_dist(cell_iter const&C) {
      real d(zero);
      vect x(pos(C));
      LoopLeafKids(cell_iter,C,l)
	if(is_active(l)) update_max(d, dist(x,pos(l))+max_dist(l));
      LoopCellKids(cell_iter,C,c)
	if(is_active(c)) update_max(d, dist(x,pos(c))+max_dist(c));
      if(d < max_dist(C)) {
	C->max_dist() = d;
	return true;
      }
      return false;
    }
    //--------------------------------------------------------------------------
    static bool update_and_parent(cell_iter&P) {
      if(update_max_dist(P)) {
	P = parent(P);
	return true;
      } else
	return false;
    }
    //--------------------------------------------------------------------------
    static void update(cell_iter A) {
      while(update_max_dist(A) && level(A)>0) A = parent(A);
    }
    //--------------------------------------------------------------------------
    static void update(cell_iter A, cell_iter B) {
      bool uA(true), uB(true);
      // bring to same level
      while(uA && level(A) > level(B))
	uA = update_and_parent(A);
      while(uB && level(B) > level(A))
	uB = update_and_parent(B);
      // advance in parallel
      while(A != B && (uA || uB) && level(A) > 0 && level(B) > 0) {
	if(uA) uA = update_and_parent(A);
	if(uB) uB = update_and_parent(B);
      }
      // advance singly
      if(uA && level(A)>0) update(A);
    }
    //--------------------------------------------------------------------------
  public:
    //--------------------------------------------------------------------------
    NeighbourFinder<false>(int k) : NeighbourFinderBase(k) {}
    //--------------------------------------------------------------------------
    bool interact(cell_iter const&A) const
    {
      // if no one is active, nothing to be done
      if(!is_active(A)) return true;
      // if not twig cell, split
      if(!is_twig(A)) return false;
      // otherwise, we must loop over body pairs and check for neighbours
      unsigned iA(0);
      many(A,iA);
      // finally update the max_dist entries in the cells
      if(iA) update(A);
      I += iA;
      return true;
    }
    //--------------------------------------------------------------------------
    bool interact(cell_iter const&A,
		  cell_iter const&B) const
    {
      // if no one is active, nothing to be done
      if(!is_active(A) && !is_active(B)) return true;
      // are there possibly any neighbouring pairs of bodies?
      real Rq = dist_sq(pos(A),pos(B));
      bool uA = is_active(A) && Rq < square(max_dist(A) + size(B));
      bool uB = is_active(B) && Rq < square(max_dist(B) + size(A));
      if(!uA && !uB) return true;
      // if not both are twig cells, split
      if(!is_twig(A) || !is_twig(B)) return false;
      // otherwise, we must loop over body pairs and check for neighbours
      unsigned iA(0), iB(0);
      if(number(A) < number(B))
	many(A,B,iA,iB);
      else
	many(B,A,iB,iA);
      // finally update the max_dist entries in the cells
      if(iA &&iB) update(A,B);
      else if(iA) update(A);
      else if(iB) update(B);
      I += iA+iB;
      return true;
    }
    //--------------------------------------------------------------------------
    bool interact(cell_iter const&A,
		  leaf_iter const&B) const
    {
      // if no one is active, nothing to be done
      if(!is_active(A) && !is_active(B)) return true;
      // are there possibly any neighbouring pairs of bodies?
      real Rq = dist_sq(pos(A),pos(B));
      bool uA = is_active(A) && Rq < square(max_dist(A));
      bool uB = is_active(B) && Rq < square(max_dist(B) + size(A));
      if(!uA && !uB) return true;
      // if A not twig cell, split
      if(!is_twig(A)) return false;
      // otherwise, we must loop over body pairs and check for neighbours
      unsigned iA(0), iB(0);
      many(B,A,iB,iA);
      // finally update the max_dist entries in the cell
      if(iA) update(A);
      I += iA+iB;
      return true;
    }
    //--------------------------------------------------------------------------
    bool interact(leaf_iter const&A,
		  cell_iter const&B) const {
      return interact(B,A);
    }
    //--------------------------------------------------------------------------
    void interact(leaf_iter const&A,
		  leaf_iter const&B) const {
      register real Rq = dist_sq(pos(A),pos(B));
      if(is_active(A) && Rq > max_dist_sq(A)) {
	A->replace_furthest(B,Rq,K);
	++I;
      }
      if(is_active(B) && Rq > max_dist_sq(B)) {
	B->replace_furthest(A,Rq,K);
	++I;
      }
    }
  };// class NeighbourFinder<false>
  //////////////////////////////////////////////////////////////////////////////
  template<> class NeighbourFinder<true> : public NeighbourFinderBase {
    void many(cell_iter const&A,
	      cell_iter const&B,
	      unsigned       &iA,
	      unsigned       &iB) const {
      LoopAllLeafs(cell_iter,A,Ai) {
	real Rq = dist_sq(pos(Ai),pos(B));
	bool uA = Rq < square(max_dist(A)+size(B));
	bool uB = Rq < square(max_dist(B));
	if(uA || uB)
	  many_YA(Ai, B.begin_leafs(), B.end_leaf_desc(), iA,iB);
      }
    }
    //--------------------------------------------------------------------------
    void many(cell_iter const&A,
	      unsigned       &iA) const {
      LoopLstLeafs(cell_iter,A,Ai)
	many_YA(Ai,Ai+1,A.end_leaf_desc(),iA,iA);
    }
    //--------------------------------------------------------------------------
    void many(leaf_iter const&A,
	      cell_iter const&B,
	      unsigned       &iA,
	      unsigned       &iB) const {
      many_YA(A,B.begin_leafs(),B.end_leaf_desc(),iA,iB);
    }
    //--------------------------------------------------------------------------
    static bool update_max_dist(cell_iter const&C) {
      real d(zero);
      vect x(pos(C));
      LoopLeafKids(cell_iter,C,l)
	update_max(d, dist(x,pos(l))+max_dist(l));
      LoopCellKids(cell_iter,C,c)
	update_max(d, dist(x,pos(c))+max_dist(c));
      if(d < max_dist(C)) {
	C->max_dist() = d;
	return true;
      }
      return false;
    }
    //--------------------------------------------------------------------------
    static bool update_and_parent(cell_iter&P) {
      if(update_max_dist(P)) {
	P = parent(P);
	return true;
      } else
	return false;
    }
    //--------------------------------------------------------------------------
    static void update(cell_iter A) {
      while(update_max_dist(A) && level(A)>0) A = parent(A);
    }
    //--------------------------------------------------------------------------
    static void update(cell_iter A, cell_iter B) {
      bool uA(true), uB(true);
      // bring to same level
      while(uA && level(A) > level(B))
	uA = update_and_parent(A);
      while(uB && level(B) > level(A))
	uB = update_and_parent(B);
      // advance in parallel
      while(A != B && (uA || uB) && level(A) > 0 && level(B) > 0) {
	if(uA) uA = update_and_parent(A);
	if(uB) uB = update_and_parent(B);
      }
      // advance singly
      if(uA && level(A)>0) update(A);
    }
  public:
    //--------------------------------------------------------------------------
    NeighbourFinder<true>(int k) : NeighbourFinderBase(k) {}
    //--------------------------------------------------------------------------
    bool interact(cell_iter const&A) const
    {
#ifdef TESTING
      std::cerr<<" CS: "<<A.index()<<": ";
      if(!is_twig(A))
	std::cerr<<" not twig ==> split\n";
#endif
      // if not twig cell, split
      if(!is_twig(A)) return false;
      // otherwise, we must loop over body pairs and check for neighbours
#ifdef TESTING
      std::cerr<<" calling many()... ";
#endif
      unsigned iA(0);
      many(A,iA);
#ifdef TESTING
      std::cerr<<" interactions: "<<iA<<'\n';
#endif
      // finally update the max_dist entries in the cells
      if(iA) update(A);
      I += iA;
      return true;
    }
    //--------------------------------------------------------------------------
    bool interact(cell_iter const&A,
		  cell_iter const&B) const
    {
#ifdef TESTING
      std::cerr<<" CC: A:"<<A.index()<<" s="<<size(A)<<" d="<<max_dist(A)<<'\n'
	       <<"     B:"<<B.index()<<" s="<<size(B)<<" d="<<max_dist(B)<<'\n';
#endif
      // are there possibly any neighbouring pairs of bodies?
      real Rq = dist_sq(pos(A),pos(B));
      bool uA = Rq < square(max_dist(A) + size(B));
      bool uB = Rq < square(max_dist(B) + size(A));
#ifdef TESTING
      std::cerr<<"     R^2="<<Rq<<" ==> uA="<<uA<<" uB="<<uB;
      if(!uA && !uB) 
	std::cerr<<" ==> no updates required\n";
      else if(!is_twig(A))
	std::cerr<<" A is not twig ==> split\n";
      else if(!is_twig(B))
	std::cerr<<" B is not twig ==> split\n";
#endif
      if(!uA && !uB) return true;
      // if not both are twig cells, split
      if(!is_twig(A) || !is_twig(B)) return false;
      // otherwise, we must loop over body pairs and check for neighbours
#ifdef TESTING
      std::cerr<<" calling many()... ";
#endif
      unsigned iA(0), iB(0);
      if(number(A) < number(B))
	many(A,B,iA,iB);
      else
	many(B,A,iB,iA);
#ifdef TESTING
      std::cerr<<" interactions: "<<iA<<':'<<iB<<'\n';
#endif
      // finally update the max_dist entries in the cells
      if(iA &&iB) update(A,B);
      else if(iA) update(A);
      else if(iB) update(B);
      I += iA+iB;
      return true;
    }
    //--------------------------------------------------------------------------
    bool interact(cell_iter const&A,
		  leaf_iter const&B) const
    {
#ifdef TESTING
      std::cerr<<" CL: A:"<<A.index()<<" s="<<size(A)<<" d="<<max_dist(A)<<'\n'
	       <<"     B:"<<A.my_tree()->index(B)
	       <<" d="<<max_dist(B)<<'\n';
#endif
      // are there possibly any neighbouring pairs of bodies?
      real Rq = dist_sq(pos(A),pos(B));
      bool uA = Rq < square(max_dist(A));
      bool uB = Rq < square(max_dist(B) + size(A));
#ifdef TESTING
      std::cerr<<"     R^2="<<Rq<<" ==> uA="<<uA<<" uB="<<uB;
      if(!uA && !uB) 
	std::cerr<<" ==> no updates required\n";
      else if(!is_twig(A))
	std::cerr<<" A is not twig ==> split\n";
#endif
      if(!uA && !uB) return true;
      // if A not twig cell, split
      if(!is_twig(A)) return false;
      // otherwise, we must loop over body pairs and check for neighbours
#ifdef TESTING
      std::cerr<<" calling many()... ";
#endif
      unsigned iA(0), iB(0);
      many(B,A,iB,iA);
#ifdef TESTING
      std::cerr<<" interactions: "<<iA<<':'<<iB<<'\n';
#endif
      // finally update the max_dist entries in the cell
      if(iA) update(A);
      I += iA+iB;
      return true;
    }
    //--------------------------------------------------------------------------
    bool interact(leaf_iter const&A,
		  cell_iter const&B) const
    {
      return interact(B,A);
    }
    //--------------------------------------------------------------------------
    void interact(leaf_iter const&A,
		  leaf_iter const&B) const
    {
#ifdef TESTING
      std::cerr<<" LL: A: d="<<max_dist(A)<<'\n'
	       <<"     B: d="<<max_dist(A)<<'\n';
#endif
      register real Rq = dist_sq(pos(A),pos(B));
#ifdef TESTING
      std::cerr<<"     R^2="<<Rq<<" ==> uA="<< (Rq > max_dist_sq(A))
	       <<" uB="<< (Rq > max_dist_sq(B)) <<'\n';
#endif
      if(Rq > max_dist_sq(A)) {
	A->replace_furthest(B,Rq,K);
	++I;
      }
      if(Rq > max_dist_sq(B)) {
	B->replace_furthest(A,Rq,K);
	++I;
      }
    }
  };// class NeighbourFinder<true>
} // namespace {
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class falcON::NeighbourLister                                                
//                                                                              
////////////////////////////////////////////////////////////////////////////////
void NeighbourLister::init_neigh(leaf_iterator const&L,
				 cell_iterator       C) const
{
  // 1. make sure cell contains enough leafs
  while(number(C) < K1) C = parent(C);
  // 2. assign neighbours
  leaf_iterator lo = L - hK;
  if(lo      < C.begin_leafs  ()) lo = C.begin_leafs();
  if(lo + K1 > C.end_leaf_desc()) lo = C.end_leaf_desc() - K1;
  unsigned k = 0;
  for(leaf_iterator l=lo; k!=K; ++l)
    if(l!=L) L->reset_neighbour(k++,l);
  // 3. update furthest
  L->update_list(K);
}
////////////////////////////////////////////////////////////////////////////////
void NeighbourLister::prepare(bool all) const
{
  unsigned NA = TREE->N_leafs();
  // - if not all, copy activity flags from bodies
  // - provide the leafs with memory to hold the neighbour list
  if(NEIB) falcON_DEL_A(NEIB);
  NEIB = falcON_NEW(Leaf::Neighbour,K*NA);
  Leaf::Neighbour*N = NEIB;
  if(all) {
    LoopLeafs(Leaf,TREE,Li) {
      Li->set_list(N);
      N += K;
    }
  } else {
    NA = 0u;
    LoopLeafs(Leaf,TREE,Li) {
      Li->copy_from_bodies_flag(TREE->my_bodies());
      if(is_active(Li)) {
	++NA;
	Li->set_list(N);
	N += K;
      }
    }
  }
  // - provide the cells with a datum to find their parent cell
  LoopCellsDown(cell_iterator,TREE,C) {
    unsigned p = TREE->NoCell(C);
    LoopCellKids(cell_iterator,C,c)
      c->set_parent(p);
  }
  root()->set_parent(0);
  // - set an initial neighbour list and find the furthest one
  // - set the cells' position, size, and distance to furthest
  // - pass activity flags up the tree (only if not all active)
  if(all) {
    LoopCellsUp(cell_iterator,TREE,C) {
      vect x(zero);
      LoopLeafKids(cell_iterator,C,l) {
	init_neigh(l,C);
	x += pos(l);
      }
      LoopCellKids(cell_iterator,C,c)
	x += number(c) * pos(c);
      x /= real(number(C));
      C->pos() = x;
      real s=zero, m=zero;
      LoopLeafKids(cell_iterator,C,l) {
	real d = dist(x,pos(l));
	update_max(s,d);
	update_max(m,d+max_dist(l));
      }
      LoopCellKids(cell_iterator,C,c) {
	real d = dist(x,pos(c));
	update_max(s,d+size(c));
	update_max(m,d+max_dist(c));
      }
      C->size    () = s;
      C->max_dist() = m;
    }
  } else {
    LoopCellsUp(cell_iterator,TREE,C) {
      vect x(zero);
      C->reset_active_flag();
      LoopLeafKids(cell_iterator,C,l) {
	if(is_active(l)) init_neigh(l,C);
	C->add_active_flag(l);
	x += pos(l);
      }
      LoopCellKids(cell_iterator,C,c) {
	C->add_active_flag(c);
	x += number(c) * pos(c);
      }
      x /= real(number(C));
      C->pos() = x;
      real s(zero), m(zero);
      LoopLeafKids(cell_iterator,C,l) {
	real d = dist(x,pos(l));
	update_max(s,d);
	if(is_active(l)) update_max(m,d+max_dist(l));
      }
      LoopCellKids(cell_iterator,C,c) {
	real d = dist(x,pos(c));
	update_max(s,d+size(c));
	if(is_active(c)) update_max(m,d+max_dist(c));
      }
      C->size    () = s;
      C->max_dist() = m;
    }
  }
#ifdef TESTING
  std::ofstream cells("cells.dat");
  TREE->dump_cells<Cell>(cells);
  std::ofstream leafs("leafs.dat");
  TREE->dump_leafs<Leaf>(leafs);
#endif
}
////////////////////////////////////////////////////////////////////////////////
void NeighbourLister::copy_dist_to_bodies(bool all) const
{
  if(all)
    LoopLeafs(Leaf,TREE,Li)
      Li->copy_to_bodies_dist(TREE->my_bodies());
  else
    LoopLeafs(Leaf,TREE,Li) if(is_active(Li))
      Li->copy_to_bodies_dist(TREE->my_bodies());
}
////////////////////////////////////////////////////////////////////////////////
void NeighbourLister::EstimateDist(bool all) const
{
  prepare(all);
  if(all) {
    NeighbourFinder<true> NF(K);
    MutualInteractor< NeighbourFinder<true> > MI(&NF,TREE->depth());
    MI.cell_self(root());
  } else {
    NeighbourFinder<false> NF(K);
    MutualInteractor< NeighbourFinder<false> > MI(&NF,TREE->depth());
    MI.cell_self(root());
  }
  copy_dist_to_bodies(all);
}
////////////////////////////////////////////////////////////////////////////////
void NeighbourLister::
Estimate(void(*f)(const bodies*, const Leaf*, const Neighbour*, int),
	 bool all) const
{
  prepare(all);
  if(all) {
    NeighbourFinder<true> NF(K);
    MutualInteractor< NeighbourFinder<true> > MI(&NF,TREE->depth());
    MI.cell_self(root());
    I = NF.N_interact();
    LoopLeafs(Leaf,TREE,Li)
      f(TREE->my_bodies(),Li,first(Li),K);
  } else {
    NeighbourFinder<false> NF(K);
    MutualInteractor< NeighbourFinder<false> > MI(&NF,TREE->depth());
    MI.cell_self(root());
    I = NF.N_interact();
    LoopLeafs(Leaf,TREE,Li) if(is_active(Li))
      f(TREE->my_bodies(),Li,first(Li),K);
  }
  
}
////////////////////////////////////////////////////////////////////////////////
