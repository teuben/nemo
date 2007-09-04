// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    src/public/neighbours.cc                                           
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    2007                                                               
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2007  Walter Dehnen                                            
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
// \version 04/09/2007 WD moved to neighbours.cc (which is superseeded)         
//                                                                              
////////////////////////////////////////////////////////////////////////////////
#include <public/neighbours.h>
#include <utils/heap.h>

using namespace falcON;

namespace falcON {
  inline bool operator<(Neighbour const&a, Neighbour const&b) {return a.Q<b.Q;}
  inline bool operator>(Neighbour const&a, Neighbour const&b) {return a.Q>b.Q;}
  inline bool operator<(real q, Neighbour const&b) {return q<b.Q;}
  inline bool operator>(real q, Neighbour const&b) {return q>b.Q;}
  inline bool operator<(Neighbour const&a, real q) {return a.Q<q;}
  inline bool operator>(Neighbour const&a, real q) {return a.Q>q;}
}
namespace {

  typedef OctTree::Leaf leaf;
  typedef OctTree::Cell cell;

  /// class providing functionality for neighbour finding
  class NeighbourSearch {
    struct CellQ {
      real Q;
      const cell*C;
      void set(real q, const cell*c) { Q=q; C=c; }
      bool operator<(CellQ const&x) const { return Q < x.Q; }
      bool operator>(CellQ const&x) const { return Q > x.Q; }
      bool operator<(real q) const { return Q < q; }
      bool operator>(real q) const { return Q > q; }
    };
    //--------------------------------------------------------------------------
    /// \name data
    //@{
    const OctTree *TREE;          ///< the tree to be used 
    const int      NDIR;          ///< direct-loop control
    const real     BIGQ;          ///< larger than largest possible Q
    unsigned       NIAC;          ///< interaction counter;
    int            M;             ///< (K - N_iac) for current search
    vect           X;             ///< position to find neighbour around
    int            K;             ///< size of list
    Neighbour     *LIST;          ///< neighbour list
    const leaf    *L;             ///< leaf for which list is made (if any)
    const cell    *C;             ///< cell already done
    //@}
    //--------------------------------------------------------------------------
    real const&radius (const cell*c) const { return TREE->rad(level(c)); }
    const cell*parent (const cell*c) const {
      return pacell(c)<0? 0 : TREE->CellNo(pacell(c));
    }
    //--------------------------------------------------------------------------
#undef LoopAllLeafs
#undef LoopLeafKids
#undef LoopCellKids
#define LoopAllLeafs(TREE,CELL,NAME)				\
    for(const leaf*NAME =TREE->LeafNo(fcleaf(CELL));		\
                   NAME!=TREE->LeafNo(ncleaf(CELL)); ++NAME)
#define LoopLeafKids(TREE,CELL,NAME)				\
    for(const leaf*NAME =TREE->LeafNo(fcleaf(CELL));		\
                   NAME!=TREE->LeafNo(ecleaf(CELL)); ++NAME)
#define LoopCellKids(TREE,CELL,NAME)				\
    for(const cell*NAME =TREE->CellNo(fccell(CELL));		\
                   NAME!=TREE->CellNo(eccell(CELL)); ++NAME)
    //--------------------------------------------------------------------------
    /// prepare tree: copy body flags & masses
    void copy_masses() const {
      for(leaf*l=TREE->begin_leafs(); l!=TREE->end_leafs(); ++l) {
	l->scalar()=TREE->my_bodies()->mass(mybody(l));
	l->copy_from_bodies_flag(TREE->my_bodies());
      }
    }
    //--------------------------------------------------------------------------
    /// distance^2 from centre of search sphere to the nearest point on a cube
    /// \param Z centre of cube
    /// \param R radius (half side) of cube
    real outside_dist_sq(vect const&Z, real R) const {
      real q=0,D;
      D=abs(Z[0]-X[0]); if(D>R) q+=square(D-R);
      D=abs(Z[1]-X[1]); if(D>R) q+=square(D-R);
      D=abs(Z[2]-X[2]); if(D>R) q+=square(D-R);
      return q;
    }
    //--------------------------------------------------------------------------
    /// distance^2 from centre of search sphere to a cell
    real outside_dist_sq(const cell*c) const {
      return outside_dist_sq(centre(c),radius(c));
    }
    //--------------------------------------------------------------------------
    /// is a sphere centred on X outside of a cubic box?
    /// \return sphere is outside of cube
    /// \param Q radius^2 of sphere
    /// \param Z centre of cube
    /// \param R radius (half side) of cube
    /// \note equivalent to, but faster than,
    /// \code Q < outside_dist_sq(C,R) \endcode
    bool outside(real Q, vect const&Z, real R) const {
      real q=0,D;
      D=abs(Z[0]-X[0]); if(D>R && Q<(q+=square(D-R))) return true;
      D=abs(Z[1]-X[1]); if(D>R && Q<(q+=square(D-R))) return true;
      D=abs(Z[2]-X[2]); if(D>R && Q<(q+=square(D-R))) return true;
      return false;
    }
    //--------------------------------------------------------------------------
    /// is search sphere outside of a cell?
    bool outside(const cell*c) const {
      return outside(LIST->Q,centre(c),radius(c));
    }
    //--------------------------------------------------------------------------
    /// is a sphere centred on X inside of a cubic box?
    /// \return sphere is inside of cube
    /// \param Q radius^2 of sphere
    /// \param Z centre of cube
    /// \param R radius (half side) of cube
    bool inside(real Q, vect const&Z, real R) const {
      real q=0,D;
      D=abs(Z[0]-X[0]); if(D>R || Q>square(R-D)) return false;
      D=abs(Z[1]-X[1]); if(D>R || Q>square(R-D)) return false;
      D=abs(Z[2]-X[2]); if(D>R || Q>square(R-D)) return false;
      return true;
    }
    //--------------------------------------------------------------------------
    /// is search sphere inside of a cell?
    bool inside(const cell*c) const {
      return inside(LIST->Q,centre(c),radius(c));
    }
    //--------------------------------------------------------------------------
    /// updates the list w.r.t. a leaf
    void add_leaf(const leaf*l) {
      real q = dist_sq(X,pos(l));
      if(LIST->Q > q) {
	LIST->Q = q;
	LIST->L = l;
	MaxHeap::after_top_replace(LIST,K);
	--M;
	++NIAC;
      }
    }
    //--------------------------------------------------------------------------
    /// updates the list w.r.t. a cell, recursive
    void add_cell(const cell*Ci, int cL=0, int cC=0) {
      if(cC==0 && number(Ci) <= max(M,NDIR)) {
	if(cL) { LoopAllLeafs(TREE,Ci,l) if(l!=L) add_leaf(l); }
	else   { LoopAllLeafs(TREE,Ci,l) add_leaf(l); }
      } else {
	if(nleafs(Ci)) {
	  if(cL) { LoopLeafKids(TREE,Ci,l) if(l!=L) add_leaf(l); }
	  else   { LoopLeafKids(TREE,Ci,l) add_leaf(l); }
	}
	if(ncells(Ci)>cC+1) {
	  CellQ Z[Nsub];
	  int   J(0);
	  LoopCellKids(TREE,Ci,c)
	    if(c!=C) Z[J++].set(outside_dist_sq(c),c);
	  MinHeap::build(Z,J);
	  while(J && LIST->Q > Z->Q) {
	    add_cell(Z->C);
	    Z[0] = Z[J-1];
	    MinHeap::after_top_replace(Z,--J);
	  }
	} else if(ncells(Ci)>cC) {
	  if(cC) { LoopCellKids(TREE,Ci,c) if(c!=C && !outside(c)) add_cell(c);}
	  else   { LoopCellKids(TREE,Ci,c) if(!outside(c)) add_cell(c); }
	}
      }
    }
    //--------------------------------------------------------------------------
    /// make list of {leaf, dist^2}, given leaf and surrouding cell
    /// \param l leaf to make list for
    /// \param P cell surrounding leaf (this is not checked, but required)
  public:
    void make_list(const leaf*l, const cell*P, Neighbour*List, int k)
    {
      K = k;
      LIST = List;
      L = l;
      X = pos(L);
      C = P;
      M = K;
      for(int k=0; k!=K; ++k)
	LIST[k].Q = BIGQ;
      for(; P && !inside(C); C=P,P=parent(C))
	add_cell(P,C==P,C!=P);
      MaxHeap::sort(LIST,K);
    }
    //--------------------------------------------------------------------------
    /// ctor
    /// \param t the tree to be used 
    /// \param k number of neighbours
    /// \param n direct-loop control; default: k/4
    NeighbourSearch(const OctTree*t, int n)
      : TREE(t), NDIR(max(1,n)),
	BIGQ(12*square(TREE->root_radius())), LIST(0), NIAC(0)
    {
      copy_masses();
    }
    //--------------------------------------------------------------------------
    unsigned const&N_iact() const { return NIAC; }
  };
}
// /////////////////////////////////////////////////////////////////////////////
void falcON::ProcessNeighbourList(const OctTree*T, int K,
				  void(*f)(const bodies*, const OctTree::Leaf*,
					   const Neighbour*, int),
				  unsigned&Ni, bool all) falcON_THROWING
{
  if(T->is_re_used())
    falcON_THROW("ProcessNeighbourList(): tree has been re-used\n");
  NeighbourSearch NS(T,K/4);
  Array<Neighbour> E(K);
  LoopCellsUp(OctTree::CellIter<OctTree::Cell>,T,C) {
    LoopLeafKids(T,C,L) if(all || is_active(L)) {
      NS.make_list(L,C,E.array(),K);
      f(T->my_bodies(),L,E.array(),K);
    }
  }
  Ni = NS.N_iact();
}
// /////////////////////////////////////////////////////////////////////////////
