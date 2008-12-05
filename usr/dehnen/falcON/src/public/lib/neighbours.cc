// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    src/public/neighbours.cc                                           
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
// \version 04/09/2007 WD moved to neighbours.cc (which is superseeded)         
// \version 06/11/2007 WD copy flags only if required.                          
// \version 20/05/2008 WD added FindLeaf(), ProcessNeighbours()                 
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

  // ///////////////////////////////////////////////////////////////////////////
  // we modify these macros to work on pointers (originals work on CellIter)
#undef LoopAllLeafs
#undef LoopLeafKids
#undef LoopCellKids
#define LoopAllLeafs(TREE,CELL,NAME)			\
    for(const leaf*NAME =TREE->LeafNo(fcleaf(CELL));	\
	NAME!=TREE->LeafNo(ncleaf(CELL)); ++NAME)
#define LoopLeafKids(TREE,CELL,NAME)			\
    for(const leaf*NAME =TREE->LeafNo(fcleaf(CELL));	\
	NAME!=TREE->LeafNo(ecleaf(CELL)); ++NAME)
#define LoopCellKids(TREE,CELL,NAME)			\
    for(const cell*NAME =TREE->CellNo(fccell(CELL));	\
	NAME!=TREE->CellNo(eccell(CELL)); ++NAME)
  // ///////////////////////////////////////////////////////////////////////////
  /// class providing basic functionality for neighbour finding
  class NeighbourSearchBase {
  protected:
    const OctTree *TREE;    ///< the tree to be used 
    vect           X;       ///< position to find neighbour around
  public:
    /// ctor
    /// \param t the tree to be used 
    NeighbourSearchBase(const OctTree*t) : TREE(t)
    {
      if(TREE->is_re_used())
	falcON_THROW("NeighbourSearchBase: cannot work with re-used tree\n");
    }
    /// ctor
    /// \param t the tree to be used 
    NeighbourSearchBase(const OctTree*t, vect const&x) : TREE(t), X(x)
    {
      if(TREE->is_re_used())
	falcON_THROW("NeighbourSearchBase: cannot work with re-used tree\n");
    }
    /// radius of given cell in TREE
    real const&radius (const cell*c) const
    {
      return TREE->rad(level(c));
    }
    /// parent cell, if any, of given cell in TREE
    const cell*parent (const cell*c) const
    {
      return pacell(c) == cell::INVALID? 0 : TREE->CellNo(pacell(c));
    }
    /// does a cell contain point X?
    /// \return point X is within cell C
    /// \param X point
    /// \param C cell
    bool contains(const cell*C) const
    {
      return abs(centre(C)[0] - X[0]) <= radius(C)
	&&   abs(centre(C)[1] - X[1]) <= radius(C)
	&&   abs(centre(C)[2] - X[2]) <= radius(C);
    }
    /// distance^2 from X to the nearest point on a cube
    /// \param Z centre of cube
    /// \param R radius (half side) of cube
    real outside_dist_sq(vect const&Z, real R) const
    {
      real q=0,D;
      D=abs(Z[0]-X[0]); if(D>R) q+=square(D-R);
      D=abs(Z[1]-X[1]); if(D>R) q+=square(D-R);
      D=abs(Z[2]-X[2]); if(D>R) q+=square(D-R);
      return q;
    }
    /// is a sphere centred on X outside of a cubic box?
    /// \return sphere is outside of cube
    /// \param Q radius^2 of sphere
    /// \param Z centre of cube
    /// \param R radius (half side) of cube
    /// \note equivalent to, but faster than,
    /// \code Q < outside_dist_sq(X,C,R) \endcode
    bool outside(real Q, vect const&Z, real R) const
    {
      real q=0,D;
      D=abs(Z[0]-X[0]); if(D>R && Q<(q+=square(D-R))) return true;
      D=abs(Z[1]-X[1]); if(D>R && Q<(q+=square(D-R))) return true;
      D=abs(Z[2]-X[2]); if(D>R && Q<(q+=square(D-R))) return true;
      return false;
    }
    /// is a sphere centred on X inside of a cubic box?
    /// \return sphere is inside of cube
    /// \param Q radius^2 of sphere
    /// \param Z centre of cube
    /// \param R radius (half side) of cube
    bool inside(real Q, vect const&Z, real R) const
    {
      real D;
      D=abs(Z[0]-X[0]); if(D>R || Q>square(R-D)) return false;
      D=abs(Z[1]-X[1]); if(D>R || Q>square(R-D)) return false;
      D=abs(Z[2]-X[2]); if(D>R || Q>square(R-D)) return false;
      return true;
    }
    /// distance^2 from X to the nearest point of cell
    /// \param C cell
    real outside_dist_sq(const cell*C) const
    {
      return outside_dist_sq(centre(C),radius(C));
    }
    /// is a sphere centred on X outside of a cell
    /// \return sphere is outside of cell
    /// \param Q radius^2 of sphere
    /// \param C cell
    bool outside(real Q, const cell*C) const
    {
      return outside(Q,centre(C),radius(C));
    }
    /// is a sphere centred on X inside of a cell
    /// \return sphere is inside of cell
    /// \param Q radius^2 of sphere
    /// \param C cell
    bool inside(real Q, const cell*C) const
    {
      return inside(Q,centre(C),radius(C));
    }
    /// given a body, find smallest surrounding cell, or root cell
    /// \note if the root cell does not contain the position return root
    const cell* findcell() const
    {
      const cell*C = TREE->FstCell();
      if(! contains(C) ) return C;
      for(;;) {
	bool incellkid=false;
	LoopCellKids(TREE,C,D) {
	  incellkid=contains(D);
	  if(incellkid) { C=D; break; }
	}
	if(!incellkid) return C;
      }
    }
  };// class NeighbourSearchBase
} // namespace P
// /////////////////////////////////////////////////////////////////////////////
const cell* falcON::FindCell(const OctTree*T, const body&B) falcON_THROWING
{
  NeighbourSearchBase N(T,pos(B));
  return N.findcell();
}
// /////////////////////////////////////////////////////////////////////////////
namespace {
  // ///////////////////////////////////////////////////////////////////////////
  /// class providing functionality for neighbour search
  class NeighbourSearch : public NeighbourSearchBase {
    /// \name data
    //@{
    const unsigned NDIR;          ///< direct-loop control
    const real  Q;                ///< radius^2 of search sphere
    const leaf *L;                ///< leaf for which list is made (if any)
    const cell *C;                ///< cell already done
    /// function to call, if any
    void      (*F)(const bodies*, const leaf*, real);
    //@}
    //--------------------------------------------------------------------------
    /// is search sphere outside of a cell?
    bool outside(const cell*c) const {
      return NeighbourSearchBase::outside(Q,c);
    }
    //--------------------------------------------------------------------------
    /// is search sphere inside of a cell?
    bool inside(const cell*c) const {
      return NeighbourSearchBase::inside(Q,c);
    }
    //--------------------------------------------------------------------------
    /// process a leaf
    void add_leaf(const leaf*l) {
      real q = dist_sq(X,pos(l));
      if(q < Q) F(TREE->my_bodies(),l,q);
    }
    //--------------------------------------------------------------------------
    /// process a cell, recursive
    /// \param Ci cell to be processed
    /// \param cL does Ci contain L?
    /// \param cC does Ci contain C?
    void add_cell(const cell*Ci, int cL=0, int cC=0) {
      if(cC==0 && number(Ci) <= NDIR) {
	if(cL) { LoopAllLeafs(TREE,Ci,l) if(l!=L) add_leaf(l); }
	else   { LoopAllLeafs(TREE,Ci,l) add_leaf(l); }
      } else {
	if(nleafs(Ci)) {
	  if(cL) { LoopLeafKids(TREE,Ci,l) if(l!=L) add_leaf(l); }
	  else   { LoopLeafKids(TREE,Ci,l) add_leaf(l); }
	}
	if(ncells(Ci)>cC) {
	  if(cC) { LoopCellKids(TREE,Ci,c) if(c!=C && !outside(c)) add_cell(c);}
	  else   { LoopCellKids(TREE,Ci,c) if(!outside(c)) add_cell(c); }
	}
      }
    }
  public:
    //--------------------------------------------------------------------------
    /// ctor
    /// \param t the tree to be used 
    /// \param n direct-loop control; default: k/4
    NeighbourSearch(const OctTree*t, vect const&x, real q,
		    void(*f)(const bodies*, const leaf*, real), unsigned n=1)
      : NeighbourSearchBase(t,x), NDIR(n), Q(q), F(f) {}
    //--------------------------------------------------------------------------
    /// process neighbours
    /// \param B body to process neighbours for
    /// \param q radius^2 of search sphere
    void process(bodies::index i)
    {
      C = findcell();
      L = 0;
      LoopLeafKids(TREE,C,l)
	if(i == mybody(l))  { L=l; break; }
      for(const cell*P=C; P && !inside(C); C=P,P=parent(C))
	add_cell(P, L && C==P, C!=P);
    }
  };// class NeighbourSearch
} // namespace {
// /////////////////////////////////////////////////////////////////////////////
void falcON::ProcessNeighbours(const OctTree*T, const body&B, real Q, 
			       void(*F)(const bodies*, const leaf*,real))
  falcON_THROWING
{
  NeighbourSearch NS(T,pos(B),Q,F);
  NS.process(static_cast<bodies::index>(B));
}
// /////////////////////////////////////////////////////////////////////////////
namespace {
  // ///////////////////////////////////////////////////////////////////////////
  /// class providing functionality for nearest neighbour search
  class NearestNeighbourSearch : public NeighbourSearchBase {
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
    const unsigned NDIR;             ///< direct-loop control
    const real     BIGQ;             ///< larger than largest possible Q
    unsigned       NIAC;             ///< interaction counter;
    int            M;                ///< (K - N_iac) for current search
    int            K;                ///< size of list
    Neighbour     *LIST;             ///< neighbour list
    const leaf    *L;                ///< leaf for which list is made (if any)
    const cell    *C;                ///< cell already done
    //@}
    //--------------------------------------------------------------------------
    /// is search sphere outside of a cell?
    bool outside(const cell*c) const {
      return NeighbourSearchBase::outside(LIST->Q,c);
    }
    //--------------------------------------------------------------------------
    /// is search sphere inside of a cell?
    bool inside(const cell*c) const {
      return NeighbourSearchBase::inside(LIST->Q,c);
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
    /// \param Ci cell to be processed
    /// \param cL does Ci contain L?
    /// \param cC does Ci contain C?
    void add_cell(const cell*Ci, int cL=0, int cC=0) {
      if(cC==0 &&
	 number(Ci) <= (M<=0? NDIR:max(static_cast<unsigned>(M),NDIR))) {
	// direct loop
	if(cL) { LoopAllLeafs(TREE,Ci,l) if(l!=L) add_leaf(l); }
	else   { LoopAllLeafs(TREE,Ci,l) add_leaf(l); }
      } else {
	// process leaf kids
      	if(nleafs(Ci)) {
	  if(cL) { LoopLeafKids(TREE,Ci,l) if(l!=L) add_leaf(l); }
	  else   { LoopLeafKids(TREE,Ci,l) add_leaf(l); }
	}
	// process cell kids
	if(ncells(Ci)>cC+1) {
	  // more than one sub-cell to process: sort in distance, then process
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
	  // only 1 sub-cell to process
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
    /// \param copy_flags if true, body flags will be copied (if present)
    NearestNeighbourSearch(const OctTree*t, unsigned n, bool copy_flags=0)
      : NeighbourSearchBase(t), NDIR(max(1u,n)),
	BIGQ(12*square(TREE->root_radius())), NIAC(0), LIST(0) 
    {
      if(copy_flags && TREE->my_bodies()->have_flag())
	for(leaf*l=TREE->begin_leafs(); l!=TREE->end_leafs(); ++l) {
	  l->scalar()=TREE->my_bodies()->mass(mybody(l));
	  l->copy_from_bodies_flag(TREE->my_bodies());
	}
      else
	for(leaf*l=TREE->begin_leafs(); l!=TREE->end_leafs(); ++l)
	  l->scalar()=TREE->my_bodies()->mass(mybody(l));
    }
    //--------------------------------------------------------------------------
    unsigned const&N_iact() const { return NIAC; }
  };// class NearestNeighbourSearch
}
// /////////////////////////////////////////////////////////////////////////////
void falcON::ProcessNearestNeighbours(const OctTree*T, int K,
				      void(*f)(const bodies*, const leaf*,
					       const Neighbour*, int),
				      unsigned&Ni, bool all) falcON_THROWING
{
  NearestNeighbourSearch NNS(T,K/4,!all);
  Array<Neighbour> E(K);
  LoopCellsUp(OctTree::CellIter<OctTree::Cell>,T,C) {
    LoopLeafKids(T,C,L) if(all || is_active(L)) {
      NNS.make_list(L,C,E.array(),K);
      f(T->my_bodies(),L,E.array(),K);
    }
  }
  Ni = NNS.N_iact();
}
// /////////////////////////////////////////////////////////////////////////////
