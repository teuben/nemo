// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    inc/public/neighbours.h                                            
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
#ifndef falcON_included_neighbours_h
#define falcON_included_neighbours_h

#ifndef falcON_included_tree_h
#  include <public/tree.h>
#endif

#define TESTING_FIRST 1
#undef  TESTING_FIRST

#ifdef  TESTING_FIRST
#define TESTING 1
#undef  TESTING
#endif
////////////////////////////////////////////////////////////////////////////////
namespace falcON {
  // ///////////////////////////////////////////////////////////////////////////
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::NeighbourLister                                              
  //                                                                            
  /// Estimator class for listing the nearest K neighbours.                     
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  // ///////////////////////////////////////////////////////////////////////////
  class NeighbourLister {
    NeighbourLister           (const NeighbourLister&);
    NeighbourLister& operator=(const NeighbourLister&);
    //--------------------------------------------------------------------------
    //                                                                          
    // sub-type NeighbourLister::Leaf                                           
    //                                                                          
    /// acts as Leaf type, holds relevant interaction data                      
    //                                                                          
    //--------------------------------------------------------------------------
#define USE_POINTER
// #undef  USE_POINTER
  public:
    class Leaf : public OctTree::Leaf {
      Leaf           (const Leaf&);
      Leaf& operator=(const Leaf&);
    public:
      /// holds information about neighbours
      struct Neighbour {
	const Leaf*L;  ///< pointer to neighbour
	real       Q;  ///< squared distance to neighbour
	/// reset data
	void reset(const Leaf*l, real q) { L=l; Q=q; }
      };
    private:
      /// const access to pointer to first neighbour
      const Neighbour*first() const { 
	return static_cast<const Neighbour*>(PROP);
      }
      /// const access to pointer to kth neighbour
      const Neighbour*neigh(int k) const { 
	return static_cast<const Neighbour*>(PROP) + k;
      }
      /// const access to pointer to furthest neighbour
      const Neighbour*furthest() const {
#ifdef USE_POINTER
 	return static_cast<const Neighbour*>(AUXP);
#else
	return neigh(AUXI);
#endif
      }
      /// const access to distance to furthest neighbour
      real      const&max_dist() const {
	return SCAL;
      }
      /// const access to distance^2 to furthest neighbour
      real      const&max_dist_sq() const {
	return furthest()->Q;
      }
    public:
      /// set Neighbour list
      void set_list(Neighbour*N) {
	PROP = static_cast<void*>(N);
      }
      /// reset Neighbour in list
      void reset_neighbour(int i, const Leaf*l) {
	real q = dist_sq(falcON::pos(this),falcON::pos(l));
	(static_cast<Neighbour*>(PROP)+i)->reset(l,q);
      }
      /// const access to distance to most distant neighbour
      friend real const&max_dist(const Leaf*);
      /// const access to distance^2 to most distant neighbour
      friend real const&max_dist_sq(const Leaf*);
      /// const access to pointer to first neighbour
      friend const Neighbour*first(const Leaf*);
      /// const access to pointer to most distant neighbour
      friend const Neighbour*furthest(const Leaf*);
      /// find most distant neighbour and update data
      /// \param K  number of neighbours in list
      void update_list(int K) {
#ifdef USE_POINTER
 	const Neighbour*const END = first()+K;
 	const Neighbour*      F   = first();
 	for(const Neighbour*N=F+1; N!=END; ++N) 
 	  if(N->Q > F->Q) F = N;
 	AUXP = const_cast<Neighbour*>(F);
 	SCAL = sqrt(F->Q);
#else
	int f=0;
	for(int k=1; k!=K; ++k)
	  if(neigh(k)->Q > neigh(f)->Q) f = k;
	AUXI = f;
 	SCAL = sqrt(neigh(f)->Q);
#endif
#ifdef TESTING_FIRST
	if(mybody().in() == 0)
	  std::cerr<<" --> update_list body 0: max_dist = "<<SCAL<<'\n';
#endif
      }
      /// if the candidate is not in the list already: replace most distant 
      /// neighbour and find new most distant
      /// \param l  leaf to replace most distant neighbour with
      /// \param d  distance^2 to that leaf
      /// \param K  number of neighbours in list
      /// \return   true if we had to replace (l was not in list already)
      bool replace_furthest(const Leaf*l, real q, int K) {
#ifdef USE_POINTER
 	// scan list for l
 	const Neighbour*const END = first()+K;
 	for(const Neighbour*N=first(); N!=END; ++N) 
 	  if(N->L == l) return false;
 	// replace furthest with new datum
 	const_cast<Neighbour*>(furthest())->reset(l,q);
 	// find new furthest
 	const Neighbour*      F   = first();
 	for(const Neighbour*N=F+1; N!=END; ++N) 
 	  if(N->Q > F->Q) F = N;
 	AUXP = const_cast<Neighbour*>(F);
 	SCAL = sqrt(F->Q);
# ifdef TESTING_FIRST
 	if(mybody().in() == 0)
 	  std::cerr<<" --> update_list body 0: max_dist = "<<SCAL<<'\n';
# endif
 	return true;
#else
 	// scan list for l
	for(int k=0; k!=K; ++k)
	  if(neigh(k)->L == l) return false;
	// replace furthest with new datum
	const_cast<Neighbour*>(furthest())->reset(l,q);
 	// find new furthest
	update_list(K);
	return true;
#endif
      }
      /// copy distance to furthest neighbour to bodies::aux()
      void copy_to_bodies_dist(const bodies*B) const {
	B->aux(mybody()) = max_dist();
      }
      //------------------------------------------------------------------------
      // dump leaf data                                                         
      //------------------------------------------------------------------------
      static void dump_head(std::ostream& o) {
	OctTree::Leaf::dump_head(o);
	o<<"        max_dist";
      }
      //------------------------------------------------------------------------
      void dump(std::ostream&o) const {
	OctTree::Leaf::dump(o);
	o<<' '<<setw(8)<<max_dist();
      }
#undef  USE_POINTER
    }; // class NeighbourLister::Leaf
    //--------------------------------------------------------------------------
    typedef Leaf::Neighbour Neighbour;
    //--------------------------------------------------------------------------
    //                                                                          
    // sub-type NeighbourLister::Cell                                           
    //                                                                          
    /// acts as Leaf type, holds relevant interaction data                      
    //                                                                          
    //--------------------------------------------------------------------------
  public:
    class Cell : public OctTree::Cell {
      Cell           (const Cell&);
      Cell& operator=(const Cell&);
      //------------------------------------------------------------------------
    public:
      typedef Leaf leaf_type;
      //------------------------------------------------------------------------
    private:
      vect        const&pos     () const { return POS; }
      real        const&size    () const { return RAD; }
      real        const&max_dist() const { return AUX1.SCAL; }
      unsigned    const&fcparent() const { return AUX2.NUMB; }
      //------------------------------------------------------------------------
    public:
      vect                 &pos     () { return POS; }
      real                 &size    () { return RAD; }
      real                 &max_dist() { return AUX1.SCAL; }
      friend vect     const&pos       (const Cell*);
      friend real     const&size      (const Cell*);
      friend real     const&max_dist  (const Cell*);
      friend unsigned const&fcparent  (const Cell*);
      //------------------------------------------------------------------------
      void mark() { flags::add(flags::marked); }
      void set_parent(unsigned p) { AUX2.NUMB = p; }
      //------------------------------------------------------------------------
      static void dump_head(std::ostream&o) {
	OctTree::Cell::dump_head(o);
	o<<"                           pos                   size   max_dist";
      }
      //------------------------------------------------------------------------
      void dump(std::ostream &o) const
      {
	OctTree::Cell::dump(o);
	o<<' '<<setw(10)<<setprecision(5)<<pos()
	 <<' '<<setw(10)<<setprecision(6)<<size()
	 <<' '<<setw(10)<<setprecision(6)<<max_dist();
      }
    };// class NeighbourLister::Cell
    //--------------------------------------------------------------------------
    // tree stuff to be superseeded                                             
    //--------------------------------------------------------------------------
    typedef Cell                         cell_type;
    typedef Leaf                         leaf_type;
    typedef OctTree::CellIter<Cell>      cell_iterator;
    typedef leaf_type*                   leaf_iterator;
  private:
    //--------------------------------------------------------------------------
    // data:                                                                    
    //--------------------------------------------------------------------------
    const OctTree           *TREE;               ///< the tree to be used       
    const unsigned           K,K1,hK;            ///< number of neighbours      
    mutable Leaf::Neighbour *NEIB;               ///< memory for neighbour lists
    mutable unsigned         I;                  ///< number of interactions    
    //--------------------------------------------------------------------------
    /// \name private methods                                                   
    //@{ -----------------------------------------------------------------------
    /// set initial guess for the K nearest neighbours from leafs in same cell
    void init_neigh(leaf_iterator const&, cell_iterator) const;
    /// prepare the tree for the interactions:
    /// - if not all, copy activity flags from bodies
    /// - provide the leafs with memory to hold the neighbour list
    /// - provide the cells with a datum to find their parent cell
    /// - set an initial neighbour list and find the furthest one
    /// - set the cells' position, size, and distance to furthers neighbour
    /// - pass activity flags up the tree
    void prepare(bool) const;
    /// copies distance to Kth neighbour to bodies
    void copy_dist_to_bodies(bool) const;
    //@} -----------------------------------------------------------------------
  public:
    //--------------------------------------------------------------------------
    const OctTree*const&my_tree() const { return TREE; }
    cell_iterator root         () const {
      return cell_iterator(TREE,static_cast<Cell*>(TREE->FstCell())); }
    //--------------------------------------------------------------------------
    // dump cell and leaf data                                                  
    //--------------------------------------------------------------------------
    void dump_cells(std::ostream&) const;
    void dump_leafs(std::ostream&) const;
    //--------------------------------------------------------------------------
    // public methods                                                           
    //--------------------------------------------------------------------------
    NeighbourLister(const OctTree*t,
		    unsigned      k) :
      TREE(t), K(k), K1(K+1), hK(K/2), NEIB(0), I(0u) {}
    //--------------------------------------------------------------------------
    ~NeighbourLister() {
      if(NEIB) {
	falcON_DEL_A(NEIB);
	NEIB = 0;
      }
    }
    //--------------------------------------------------------------------------
    /// report number of interactions
    unsigned const&N_interact() const { return I; }
    //--------------------------------------------------------------------------
    /// for active bodies: writes distance to Kth neighbour into aux()
    void EstimateDist(bool) const;
    //--------------------------------------------------------------------------
    /// apply user provided function to neighbour list
    /// \param f    function called for all active leafs
    /// \param all  are all bodies supposed to be active?
    void Estimate(void(*f)(const bodies*, const Leaf*, const Neighbour*, int),
		  bool all) const;
  };
  //////////////////////////////////////////////////////////////////////////////
  inline real const&
  max_dist(const NeighbourLister::Leaf*L) { return L->max_dist(); }
  inline real const&
  max_dist_sq(const NeighbourLister::Leaf*L) { return L->max_dist_sq(); }
  inline const NeighbourLister::Leaf::Neighbour*
  first(const NeighbourLister::Leaf*L) { return L->first(); }
  inline const NeighbourLister::Leaf::Neighbour*
  furthest(const NeighbourLister::Leaf*L) { return L->furthest(); }
  //////////////////////////////////////////////////////////////////////////////
  inline vect     const&pos     (const NeighbourLister::Cell*C) {
    return C->POS; }
  inline real     const&size    (const NeighbourLister::Cell*C) {
    return C->RAD; }
  inline real     const&max_dist(const NeighbourLister::Cell*C) {
    return C->AUX1.SCAL; }
  inline unsigned const&fcparent(const NeighbourLister::Cell*C) {
    return C->AUX2.NUMB; }
} // namespace falcON {
////////////////////////////////////////////////////////////////////////////////
falcON_TRAITS(falcON::NeighbourLister::Leaf,"NeighbourLister::Leaf");
falcON_TRAITS(falcON::NeighbourLister::Cell,"NeighbourLister::Cell");
////////////////////////////////////////////////////////////////////////////////
#endif // falcON_included_neighbours_h  
