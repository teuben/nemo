// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/octtree.h
///
/// \brief  methods for building and walking an octtree in 2D or 3D as well as
///         interaction and neighbour-search algorithms using this tree
///
/// \author Walter Dehnen
///
/// \date   2009,2010
///
/// \version May-2009 WD  first tested version
/// \version Oct-2009 WD  new design using indices for leafs and cells
/// \version Nov-2009 WD  removed redundant template parameter
/// \version Jan-2010 WD  renamed methods in TreeAccess; added leaf's parent
/// \version Jan-2010 WD  neighbour search methods
/// \version Feb-2010 WD  new initialisation: removed need for OctalTree::Dot
/// \version Mar-2010 WD  class FastNeighbourFinder
/// \version Apr-2010 WD  class TreeWalkAlgorithm, tree pruning
/// \version Apr-2010 WD  faster and memory-leaner tree-building algorithm
/// \version Jun-2010 WD  16-byte alignement, using geometry.h
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009,2010 Walter Dehnen
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
#ifndef WDutils_included_octtree_h
#define WDutils_included_octtree_h

#ifndef WDutils_included_iostream
#  include <iostream>
#  define WDutils_included_iostream
#endif
#ifndef WDutils_included_iomanip
#  include <iomanip>
#  define WDutils_included_iomanip
#endif
#ifndef WDutils_included_memory_h
#  include <memory.h>
#endif
#ifndef WDutils_included_geometry_h
#  include <geometry.h>
#endif

namespace {
  template<int, typename> struct BoxDotTree;
}
namespace WDutils {
  template<typename> struct TreeAccess;
  ///
  /// A spatial tree of square (2D) or cubic (3D) cells
  ///
  /// Cells with more than @a nmax (argument to constructor and member
  /// rebuild) are split and octants (or quarters for Dim=2) with more than
  /// @a nmin particles are assigned a new cell. After tree construction, the
  /// octants are dissolved (though each cell still knows its octant in its
  /// parent cell) and any leafs from single-leaf octants become direct
  /// leaf-children of their cell. Thus, a non-final cell (cell with daughter
  /// cells) can also contain leafs which are in none of their daughter cells;
  /// these leafs are placed first in the array of a cell's leafs and are
  /// referred to as 'leaf kids' as opposed to 'leaf descendants', which
  /// includes all leafs contained within a cell.
  ///
  /// \note Access to leafs and cells via struct TreeAccess
  /// \note Implementations for @a __X = float,double and @a __D = 2,3
  template<int __D, typename __X>
  class OctalTree {
    friend struct BoxDotTree<__D,__X>;
    friend struct TreeAccess<OctalTree>;
    /// ensure that the only valid instantinations are those in octtree.cc
    WDutilsStaticAssert
    (( (  __D == 2 || __D == 3 )              &&
       meta::TypeInfo<__X>::is_floating_point    ));
    // disable default and copy ctor
    OctalTree           (const OctalTree&);  // not implemented
    OctalTree& operator=(const OctalTree&);  // not implemented
  public:
    /// \name public constants and types
    //@{
    typedef __X             real;             ///< floating point type: position
    typedef tupel<__D,real> point;            ///< type: positions
    typedef uint32          particle_key;     ///< type: indexing particles
    typedef uint32          count_type;       ///< type: indexing leafs & cells
    typedef uint8           depth_type;       ///< type: tree depth & level
    typedef uint8           octant_type;      ///< type: octant and # cell kids
    typedef uint16          local_count;      ///< type: # leaf kids
    typedef Geometry::cube<__D,real> cube;    ///< type: cubic box
    typedef SSE::Extend16<point> point16;
    typedef SSE::Extend16<cube>  cube16;
    static const depth_type Dim     = __D;    ///< number of dimensions
    static const depth_type Nsub    = 1<<Dim; ///< number of octants per cell
    static const depth_type MAXNMAX = 250;
    /// the maximum allowed tree depth
    /// \note Set to the number of digits in our floating point type.
    static const depth_type MaximumDepth=std::numeric_limits<real>::digits;
    /// Interface for initialising particles at tree building
    ///
    /// The interface here is intended to allow for a general particle data
    /// layout in the user application and at the same time ensure that at
    /// OctalTree::rebuild() we keep the particle order (as much as possible),
    /// enabling significant speed-up of tree re-building.
    class Initialiser {
      friend struct BoxDotTree<__D,__X>;
    protected:
      // virtual dtor: only needed with older compiler versions
      virtual ~Initialiser() {}
      /// initialise key and position for one particle.
      /// \param[out] I  key for particle
      /// \param[out] X  position for particle
      /// \note Will be called @a N times in the first constructor of class
      ///       OctalTree. The implementation (the non-abstract version in any
      ///       derived class) must ensure that each call initialises another
      ///       particle.
      virtual void Initialise(particle_key&I, point&X) const = 0;
      /// re-initialise position for valid key only
      /// \param[in]  I  original particle key, may have become invalid
      /// \param[out] X  if key @a I is valid: position for associated particle
      /// \return        was key @a I valid and position initialised?
      /// \note Will be called min(@a Nnew, @a Nold) times in
      ///       OctalTree::rebuild() in an attempt to re-initialise all
      ///       particles in the particle order of the old tree. The option to
      ///       return false allows for the possibility that a particle key has
      ///       become invalid by some data re-arrangement. See also the
      ///       documentation for ReInitialiseInvalid() below.
      virtual bool ReInitialiseValid(particle_key I, point&X) const = 0;
      /// re-initialise key and position.
      /// \param[out] I  valid key for particle
      /// \param[out] X  position for particle
      /// \note Will be called in OctalTree::rebuild() @b after calling
      ///       ReInitialiseValid() min(@a Nnew, @a Nold) times. This is to
      ///       re-initialise (1) particles for which the original key from
      ///       the tree prior to rebuilding has become invalid and (2) any
      ///       surplus particles (if @a Nnew > @a Nold).
      /// \note The implementation (the non-abstract version in any derived
      ///       class) must ensure that each call initialises another particle.
      virtual void ReInitialiseInvalid(particle_key&I, point&X) const = 0;
      /// pick particles for building a pruned tree
      /// \param[in] I  particle key of leaf in parent tree
      /// \return       shall we include this leaf/particle in pruned tree?
      virtual bool Pick(particle_key I) const = 0;
    };
    //@}
    /// \name general data for class OctalTree
    //@{
  private:
    char*             ALLOC;          ///< actually allocated memory
    size_t            NALLOC;         ///< # bytes allocated at ALLOC
    const depth_type  NMAX;           ///< N_max
    const depth_type  NMIN;           ///< N_min
    const bool        AVSPC;          ///< avoid single parent cells?
    depth_type        DEPTH;          ///< tree depth
    //@}
    /// \name leaf data (access via struct TreeAccess )
    //@{
    count_type    NLEAF;              ///< total number of leafs
    point16      *XL;                 ///< leaf positions
    particle_key *PL;                 ///< index of associated particle
    count_type   *PC;                 ///< index of parent cell
    //@}
    /// \name cell data (access via struct TreeAccess )
    //@{
    count_type    NCELL;              ///< total number of cells
    depth_type   *LE;                 ///< cells' tree level
    octant_type  *OC;                 ///< cells' octant in parent cell
    cube16       *XC;                 ///< cells' centre and radius (of cube)
    count_type   *L0;                 ///< cells' first leaf
    local_count  *NL;                 ///< number of cells' leaf kids
    count_type   *NM;                 ///< number of cells' leaf descendants
    count_type   *CF;                 ///< first daughter cell
    octant_type  *NC;                 ///< number of cells' daughter cells
    count_type   *PA;                 ///< parent cell
    //@}
    void allocate();
    void build(char, count_type, const Initialiser*, const OctalTree*)
      WDutils_THROWING;
  public:
#ifdef __INTEL_COMPILER
#pragma warning (disable:2259) /* non-pointer conversion from A to B may lose significant bits */
#endif
    /// \name tree building
    //@{
    /// build tree from scratch.
    ///
    /// The tree is build in three stages. First, Initialiser::Initialise() is
    /// called @a N times to obtain key and position for all @a N particles.
    /// \n
    /// Second, a 'box-dot' tree is built using an algorithm which adds one
    /// particle at a time and splits octants in excess of @a nmax particles.
    /// \n
    /// Third, the cell-leaf tree is established by mapping boxes and octants
    /// with at least @a nmin particles to cells and dots to leafs in such a
    /// way that any cell's leaf and cell descendants are contiguous in memory.
    ///
    /// \param[in] N      number of particles
    /// \param[in] init   Initialiser for particle keys and positions
    /// \param[in] nmax   maximum number of particles in unsplit octants
    /// \param[in] nmin   minimum number of particles in cell
    /// \note We require @a nmin,nmax<=250 but map @a nmin=0 to @a
    ///       nmin=min(2,nmax).
    ///       \n
    ///       For @a nmin=2 (and @a nmax>1) the cell-leaf tree reflects the
    ///       depth of the original box-dot tree. However, for @a nmin > 2,
    ///       while the cell-leaf tree is not as deep, the tree order of the
    ///       deeper box-dot tree is preserved in the tree order of the leafs.
    ///       \n
    ///       For @a nmax=nmin=1 the tree is build to maximum depth when each
    ///       leaf has its own final cell (and tree building is most CPU time
    ///       consuming), which is unlikely to be required by any application.
    ///
    /// \param[in] avspc  avoid single-parent cells
    /// \note Single-parent cells occur if all particles are in just one
    ///       octant. If @a avspc=true (the default) such cells are eliminated
    ///       in favour of their only daughter cell. In this case, the parent
    ///       and daughter cell may be more than one tree level apart and the
    ///       tree depth less than the maximum tree level of any cell.
    OctalTree(count_type N, const Initialiser*init, depth_type nmax,
	      depth_type nmin=0, bool avspc=true) WDutils_THROWING
    : ALLOC ( 0 ),
      NALLOC( 0 ),
      NMAX  ( nmax<1? 1 : (nmax>MAXNMAX? MAXNMAX:nmax) ),
      NMIN  ( nmin? (nmin<NMAX? nmin:NMAX) : (NMAX<2?NMAX:2) ),
      AVSPC ( avspc )
    { 
      if(N == 0)
	WDutils_THROWN  ("OctalTree<%d,%s>: N=0\n",Dim,nameof(real));
      if(init == 0)
	WDutils_THROWN  ("OctalTree<%d,%s>: init=0\n",Dim,nameof(real));
      if(nmax == 0)
	WDutils_WarningN("OctalTree<%d,%s>: "
			 "nmax=%d; will use nmax=%d instead\n",
			 Dim,nameof(real), nmax,int(NMAX));
      if(nmax > 250)
	WDutils_WarningN("OctalTree<%d,%s>: "
			 "nmax=%d exceeds 250; will use nmax=%d instead\n",
			 Dim,nameof(real),nmax, int(NMAX));
      if(nmin > 250)
	WDutils_WarningN("OctalTree<%d,%s>: "
			 "nmin=%d exceeds 250; will use nmin=%d instead\n",
			 Dim,nameof(real),nmin, int(NMIN));
      if(nmin > nmax)
	WDutils_WarningN("OctalTree<%d,%s>: "
			 "nmin=%d exceeds nmax=%d; will use nmin=%d\n",
			 Dim,nameof(real), nmin, nmax, int(NMIN));
      build('n', N, init, 0);
    }
    /// re-build the tree after particles have changed (position or number).
    ///
    /// The tree is build exactly in the same way as with the constructor, only
    /// the order in which the particles are added to the 'box-dot' tree is the
    /// leaf order of the original tree (rather than ascending particle index).
    /// This change alone results in a significant speed-up for the whole
    /// process, because it avoids cache misses.
    ///
    /// \param[in] init    Initialiser required to re-initialise particle data
    /// \param[in] Nnew    new number of particles (if @a Nnew=0 we assume the
    ///                    number has not changed)
    /// \param[in] nmax    maximum number of particles in unsplit octants
    /// \param[in] nmin    minimum number of particles in cell
    /// \note If @a nmax=0 we take the old values for both @a nmin and @a nmax.
    ///       Otherwise, we require @a nmin,nmax<=250 and map @a nmin=0 to @a
    ///       nmin=min(2,nmax).
    ///
    /// \note First, Initialiser::ReInitialiseValid() is called min(@a Nnew, @a
    ///       Nold) times in an attempt to re-initialise all particles in the
    ///       original particle order of the existing tree. Then,
    ///       Initialiser::ReInitialiseInvalid() is called for all particles
    ///       whose keys have become invalid as indicated by the return value
    ///       of Initialiser::ReInitialiseValid(). Finally,
    ///       Initialiser::ReInitialiseInvalid() is called to initialise any
    ///       additional particles.
    void rebuild(const Initialiser*init, count_type Nnew=0,
		 depth_type nmax=0, depth_type nmin=0) WDutils_THROWING
    {
      if(0==init)
	WDutils_THROW("OctalTree<%d,%s>::rebuild(): init=0\n",Dim,nameof(real));
      if(nmax!=0) {
	const_cast<depth_type&>(NMAX) = nmax>MAXNMAX? MAXNMAX:nmax;
	const_cast<depth_type&>(NMIN) = 
	  nmin ? (nmin<NMAX? nmin:NMAX) : (NMAX<2?NMAX:2);
	if(nmax > 250)
	  WDutils_WarningN("OctalTree<%d,%s>::rebuild(): "
			   "nmax=%d exceeds 250; will use nmax=%d instead\n",
			   Dim,nameof(real), nmax,int(NMAX));
	if(nmin > 250)
	  WDutils_WarningN("OctalTree<%d,%s>::rebuild(): "
			   "nmin=%d exceeds 250; will use nmin=%d instead\n",
			   Dim,nameof(real), nmin,int(NMIN));
	if(nmin > nmax)
	  WDutils_WarningN("OctalTree<%d,%s>::rebuild(): "
			   "nmin=%d exceeds nmax=%d; will use nmin=%d\n",
			   Dim,nameof(real), nmin,nmax,int(NMIN));
      }
      build('r', Nnew?Nnew:NLEAF, init, this);
    }
    /// make a pruned version of another octtree
    ///
    /// The new tree contains all leafs of the parent tree for which
    /// Initialiser::Pick(l) returns true.
    ///
    /// \note If all parent-tree leafs are picked, a warning is issued.
    ///       Conversely, if none is picked, an error is thrown.
    ///
    /// \param[in] parent  parent tree to prune
    /// \param[in] init    Initialiser, used to pick leafs
    /// \param[in] nsub    number of particles in pruned tree
    /// \note If @a nsub==0, the number is actually counted (and
    ///       Initialiser::Pick() called twice for each leaf of the parent
    ///       tree). Otherwise, if @a nsub>0, it is expected that at most @a
    ///       nsub particles are in the pruned tree (an error is thrown if more
    ///       are found).
    ///
    /// \param[in] nmax    maximum number of leafs in unsplit octants
    /// \param[in] nmin    minimum number of particles in cell
    /// \note If @a nmax=0 we take the values for both @a nmin and @a nmax from
    ///       the parent tree. Otherwise, we require @a nmin,nmax<=250 but map
    ///       @a nmin=0 to @a nmin=min(2,nmax).
    ///
    /// \param[in] avspc   avoid single-parent cells
    OctalTree(const OctalTree*parent, const Initialiser*init,
	      count_type nsub=0, depth_type nmax=0, depth_type nmin=0,
	      bool avspc=true) WDutils_THROWING
    : ALLOC ( 0 ),
      NALLOC( 0 ),
      NMAX  ( nmax==0? parent->Nmax() : (nmax>MAXNMAX? MAXNMAX:nmax) ),
      NMIN  ( nmax==0? parent->Nmin() : 
	      nmin? (nmin<NMAX? nmin:NMAX) : (NMAX<2?NMAX:2) ),
      AVSPC ( avspc )
    {
      if(0==init)
	WDutils_THROW   ("OctalTree<%d,%s>: init=0\n",Dim,nameof(real));
      if(nmax > 250)
	WDutils_WarningN("OctalTree<%d,%s>: "
			 "nmax=%d exceeds 250; will use nmax=%d instead\n",
			 Dim,nameof(real),nmax, int(NMAX));
      if(nmin > 250)
	WDutils_WarningN("OctalTree<%d,%s>: "
			 "nmin=%d exceeds 250; will use nmin=%d instead\n",
			 Dim,nameof(real),nmin, int(NMIN));
      if(nmin > nmax)
	WDutils_WarningN("OctalTree<%d,%s>: "
			 "nmin=%d exceeds nmax=%d; will use nmin=%d\n",
			 Dim,nameof(real), nmin, nmax, int(NMIN));
      build('p', nsub, init, parent);
    }
#ifdef __INTEL_COMPILER
#pragma warning (restore:2259) /* non-pointer conversion from A to B may lose significant bits */
#endif
    /// establish as pruned version of an existing octtree
    ///
    /// This is equivalent to destruction followed by construction as pruned
    /// version of another octtree: the tree will contain all leafs of the
    /// parent tree for which Initialiser::Pick() returns true.
    ///
    /// \note If all parent-tree leafs are picked, a warning is issued.
    ///       Conversely, if none is picked, an error is thrown.
    ///
    /// \param[in] parent  parent tree to prune
    /// \param[in] init    Initialiser, used to pick leafs
    /// \param[in] nsub    number of particles in pruned tree
    /// \note If @a nsub==0, the number is actually counted (and
    ///       Initialiser::Pick() called twice for each leaf of the parent
    ///       tree). Otherwise, if @a nsub>0, it is expected that at most @a
    ///       nsub particles are in the pruned tree (an error is thrown if more
    ///       are found).
    ///
    /// \param[in] nmax    maximum number of leafs in unsplit octants
    /// \param[in] nmin    minimum number of particles in cell
    /// \note If @a nmax=0 we take the values for both @a nmin and @a nmax from
    ///       the parent tree. Otherwise, we require @a nmin,nmax<=250 but map
    ///       @a nmin=0 to @a nmin=min(2,nmax).
    void reprune(const OctalTree*parent, const Initialiser*init,
		 count_type nsub=0, depth_type nmax=0, depth_type nmin=0)
      WDutils_THROWING
    {
      if(0==init)
	WDutils_THROW("OctalTree<%d,%s>::reprune(): init=0\n",Dim,nameof(real));
      if(nmax==0) {
	const_cast<depth_type&>(NMAX) = parent->Nmax();
	const_cast<depth_type&>(NMIN) = parent->Nmin();
      } else {
	const_cast<depth_type&>(NMAX) = nmax>MAXNMAX? MAXNMAX:nmax;
	const_cast<depth_type&>(NMIN) = 
	  nmin ? (nmin<NMAX? nmin:NMAX) : (NMAX<2?NMAX:2);
	if(nmax > 250)
	  WDutils_WarningN("OctalTree<%d,%s>::reprune(): "
			   "nmax=%d exceeds 250; will use nmax=%d instead\n",
			   Dim,nameof(real), nmax, int(NMAX));
	if(nmin > 250)
	  WDutils_WarningN("OctalTree<%d,%s>::reprune(): "
			   "nmin=%d exceeds 250; will use nmin=%d instead\n",
			   Dim,nameof(real), nmin, int(NMIN));
	if(nmin > nmax)
	  WDutils_WarningN("OctalTree<%d,%s>::reprune(): "
			   "nmin=%d exceeds nmax=%d; will use nmin=%d\n",
			   Dim,nameof(real), nmin, nmax, int(NMIN));
      }
      build('p', nsub, init, parent);
    }
    //@}
    /// dtor
    ~OctalTree()
    {
      if(ALLOC) WDutils_DEL16(ALLOC);
      ALLOC = 0;
      NALLOC= 0;
      NLEAF = 0;
      NCELL = 0;
      DEPTH = 0;
    }
    /// was @a avspc set during tree building?
    bool const&AvoidedSingleParentCells() const
    { return AVSPC; }
    /// tree depth
    depth_type const&Depth() const
    { return DEPTH; }
    /// root radius
    real const&RootRadius() const
    { return XC->H; }
    /// N_max
    depth_type const&Nmax() const
    { return NMAX; }
    /// N_min
    depth_type const&Nmin() const
    { return NMIN; }
    /// total number of leafs
    count_type const&NLeaf() const
    { return NLEAF; }
    /// total number of cells
    count_type const&NCell() const
    { return NCELL; }
    /// centre suitable for root
    /// \return position suitable for root centre
    /// \param[in] xave  some average or mean position within the root box
    static point RootCentre(point const&xave);
    /// radius suitable for root
    /// \return radius suitable for root
    /// \param[in] xcen  root centre
    /// \param[in] xmin  most extreme position with smallest coordinates
    /// \param[in] xmax  most extreme position with largest coordinates
    static real RootRadius(point const&xcen, point const&xmin,
			   point const&xmax);
  };// class OctalTree<>

  template<typename Tree> struct TreeAccess; // forward declaration
  ///
  /// support for using an OctalTree.
  ///
  /// Holds just a pointer to an OctalTree and provides access to leaf and
  /// cell data, as well as member methods for walking the tree.
  ///
  /// Essential as base class for tree-walking algorithms.
  ///
  /// \relates WDutils::OctalTree
  template<int __D, typename __X>
  struct TreeAccess< OctalTree<__D,__X> >
  {
    /// tree
    typedef OctalTree<__D,__X> Tree;
    /// floating-point type
    typedef typename Tree::real real;
    /// type for positions
    typedef typename Tree::point point;
    /// type for cubic boxes
    typedef typename Tree::cube cube;
    /// type for particle index
    typedef typename Tree::particle_key particle_key;
    /// type for indexing leafs and cells
    typedef typename Tree::count_type count_type;
    /// type for tree depth & level
    typedef typename Tree::depth_type depth_type;
    /// type for cell octants
    typedef typename Tree::octant_type octant_type;
    /// type for number of leaf kids and Nmax
    typedef typename Tree::local_count local_count;
    /// number of dimensions
    static const depth_type Dim  = Tree::Dim;
    /// number of octants per cell
    static const depth_type Nsub = Tree::Nsub;
    /// pointer to tree
    const Tree*const TREE;
    /// ctor
    TreeAccess(const Tree*t) : TREE(t) {}
    /// copy ctor
    TreeAccess(const TreeAccess&t) : TREE(t.TREE) {}
    // virtual dtor required by some old compilers
    virtual ~TreeAccess() {}
    /// iterator used for leafs.
    /// A simple wrapper around an index, which, being a separate type, avoids
    /// confusion with other indices or variables of type count_type.
    //
    //  NOTE A conversion to count_type is not a good idea, as it allows an
    //       implicit conversion to bool, which almost certainly results in
    //       behaviour that is not intended, i.e. instead of IsInvalid()
    struct Leaf {
      count_type I;
      /// default ctor
      Leaf()
      {}
      /// ctor from index
      explicit Leaf(count_type i)
	: I(i) {}
      /// prefix increment
      Leaf& operator++()
      { ++I; return*this; }
      /// comparison <
      bool operator < (Leaf l) const
      { return I< l.I; }
      /// comparison <=
      bool operator <=(Leaf l) const
      { return I<=l.I; }
      /// comparison ==
      bool operator ==(Leaf l) const
      { return I==l.I; }
      /// comparison !=
      bool operator !=(Leaf l) const
      { return I!=l.I; }
    };
    /// type for compact indexing of leafs
    /// \note Useful for templates over tree type (OctalTree or ParallelTree)
    typedef Leaf LeafIdentifier;
    //
    static LeafIdentifier identifier(Leaf l)
    { return l; }
    //
    static Leaf leaf(LeafIdentifier l)
    { return l; }
    /// \name leaf data access
    //@{
    /// leaf position, 16-byte aligned
    point const&position(Leaf l) const
    { return TREE->XL[l.I]; }
    /// index of associated particle
    particle_key const&particle(Leaf l) const
    { return TREE->PL[l.I]; }
    //@}
    /// iterator used for cells.
    /// A simple wrapper around an index, which, being a separate type, avoids
    /// confusion with other indices or variables of type count_type.
    //
    //  NOTE A conversion to count_type is not a good idea, as it allows an
    //       implicit conversion to bool, which almost certainly results in
    //       behaviour that is not intended, i.e. instead of IsInvalid()
    struct Cell {
      count_type I;
      /// default ctor
      Cell()
      {}
      /// ctor from index
      explicit Cell(count_type i)
	: I(i) {}
      /// prefix increment
      Cell& operator++()
      { ++I; return*this; }
      /// prefix decrement
      Cell& operator--()
      { --I; return*this; }
      /// comparison <
      bool operator < (Cell c) const
      { return I< c.I; }
      /// comparison <=
      bool operator <=(Cell c) const
      { return I<=c.I; }
      /// comparison >
      bool operator > (Cell c) const
      { return I> c.I; }
      /// comparison >=
      bool operator >=(Cell c) const
      { return I>=c.I; }
      /// comparison ==
      bool operator ==(Cell c) const
      { return I==c.I; }
      /// comparison !=
      bool operator !=(Cell c) const
      { return I!=c.I; }
    };
    /// \name cell data access
    //@{
    /// tree level of cell
    depth_type const&level(Cell c) const
    { return TREE->LE[c.I]; }
    /// octant of cell in parent
    octant_type const&octant(Cell c) const
    { return TREE->OC[c.I]; }
    /// cell's cubic box, 16-byte aligned
    cube const&box(Cell c) const
    { return TREE->XC[c.I]; }
    /// cell's geometric centre (of cubic box)
    point const&centre(Cell c) const
    { return box(c).X; }
    /// radius (half-side-length of box) of cell
    real const&radius(Cell c) const
    { return box(c).H; }
    /// number of leaf kids
    local_count const&NLeafKids(Cell c) const
    { return TREE->NL[c.I]; }
    /// total number of leafs
    count_type const&Number(Cell c) const
    { return TREE->NM[c.I]; }
    /// number of daughter cells
    octant_type const&NCell(Cell c) const
    { return TREE->NC[c.I]; }
    /// index of cell's first leaf
    count_type const&firstleafindex(Cell c) const
    { return TREE->L0[c.I]; }
    /// index of cell's first daughter cell, if any
    count_type const&firstcellindex(Cell c) const
    { return TREE->CF[c.I]; }
    //@}
    /// \name tree walking and related
    //@{
    /// root radius
    real const&RootRadius() const
    { return TREE->RootRadius(); }
    /// was @a avspc set during tree building?
    bool const&AvoidedSingleParentCells() const
    { return TREE->AvoidedSingleParentCells(); }
    /// N_max
    depth_type const&Nmax() const
    { return TREE->Nmax(); }
    /// N_min
    depth_type const&Nmin() const
    { return TREE->Nmin(); }
    /// tree depth
    depth_type const&Depth() const
    { return TREE->Depth(); }
    /// # leafs
    count_type const&NLeaf() const
    { return TREE->NLeaf(); }
    /// next leaf
    static Leaf next(Leaf l)
    { return Leaf(l.I+1); }
    /// first leaf
    static Leaf BeginLeaf()
    { return Leaf(0u); }
    /// end of leafs (beyond last leaf)
    Leaf EndLeaf() const
    { return Leaf(TREE->NLEAF); }
    /// an invalid leaf
    Leaf InvalidLeaf() const
    { return Leaf(TREE->NLEAF); }
    /// is @a l a valid leaf?
    bool IsValid(Leaf l) const
    { return l.I < TREE->NLEAF; }
    /// # cells
    count_type const&NCell() const
    { return TREE->NCell(); }
    /// next cell
    static Cell next(Cell c)
    { return Cell(c.I+1); }
    /// root cell
    static Cell Root()
    { return Cell(0u); }
    /// first cell
    static Cell BeginCell()
    { return Cell(0u); }
    /// end of cells (beyond last cell)
    Cell EndCell() const
    { return Cell(TREE->NCELL); }
    /// first cell in reversed order: last cell
    Cell RBeginCell() const
    { return Cell(TREE->NCELL-1); }
    /// end cell in reversed order: invalid Cell
    static Cell REndCell()
    { return --(Cell(0)); }
    /// an invalid cell
    Cell InvalidCell() const
    { return Cell(TREE->NCELL); }
    /// is cell @a c valid?
    bool IsValid(Cell c) const
    { return c.I < TREE->NCELL; }
    /// first of cell @a c's leafs
    Leaf BeginLeaf(Cell c) const
    { return Leaf(firstleafindex(c)); }
    /// end of cell @a c's leaf children 
    Leaf EndLeafKids(Cell c) const
    { return Leaf(firstleafindex(c)+NLeafKids(c)); }
    /// end of all of cell @a c's leafs
    Leaf EndLeafDesc(Cell c) const
    { return Leaf(firstleafindex(c)+Number(c)); }
    /// last of all of cell @a c's leafs
    Leaf LastLeafDesc(Cell c) const
    { return Leaf(firstleafindex(c)+Number(c)-1); }
    /// first of cell @a c's daughter cells
    /// \note if @a c has no daughter cells, this returns @a c itself
    Cell BeginCell(Cell c) const
    { return Cell(firstcellindex(c)); }
    /// end of cell @a c's daughter cells
    /// \note if @a c has no daughter cells, this returns @a c itself
    Cell EndCell(Cell c) const
    { return Cell(firstcellindex(c)+NCell(c)); }
    /// last of cell @a c's daughter cells
    Cell RBeginCell(Cell c) const
    { return --(EndCell(c)); }
    /// before first of cell @a c's daughter cells
    Cell REndCell(Cell c) const
    { return --(BeginCell(c)); }
    /// cell @a c's parent cell
    Cell Parent(Cell c) const
    { return Cell(c.I? TREE->PA[c.I] : TREE->NCELL); }
    /// leaf @a l's parent cell
    Cell Parent(Leaf l) const
    { return Cell(TREE->PC[l.I]); }
    /// does cell @a c contain leaf @a l ?
    bool Contains(Cell c, Leaf l) const
    { return BeginLeaf(c) <= l && l < EndLeafDesc(c); }
    /// is either cell ancestor of the other?
    bool IsAncestor(Cell a, Cell b) const
    { return maxnorm(centre(a)-centre(b)) < max(radius(a),radius(b)); }
    /// find smallest cell containing a given position
    /// \param[in] x   position to find cell for
    /// \return        smallest tree cell containing @a x
    /// \note If @a x is outside the root cell, an invalid cell is returned.
    /// \note For @a x == position(Leaf), this is equivalent to, but slower
    ///       than, Parent(Leaf).
    Cell SmallestContainingCell(point const&x) const;
    //@}
    /// \name macros for tree walking from within a TreeAccess<OctalTree>
    //@{
    /// loop cells down: root first
    /// \relates WDutils::TreeAccess
    /// \note useful for a down-ward pass
#ifndef LoopCellsDown
# define LoopCellsDown(NAME)			\
    for(Cell NAME = this->BeginCell();		\
	NAME != this->EndCell(); ++NAME)
#endif
    /// loop cells up: root last
    /// \relates WDutils::TreeAccess
    /// \note useful for an up-ward pass
#ifndef LoopCellsUp
# define LoopCellsUp(NAME)			\
    for(Cell NAME = this->RBeginCell();	\
	NAME != this->REndCell(); --NAME)
#endif
    /// loop leafs
    /// \relates WDutils::TreeAccess
#ifndef LoopLeafs
# define LoopLeafs(NAME)			\
    for(Leaf NAME = this->BeginLeaf();		\
	NAME != this->EndLeaf(); ++NAME)
#endif
    /// loop cell kids of a given cell
    /// \relates WDutils::TreeAccess
#ifndef LoopCellKids
# define LoopCellKids(CELL,NAME)		\
    for(Cell NAME = this->BeginCell(CELL);	\
	NAME != this->EndCell(CELL); ++NAME)
#endif
    /// loop cell kids of a given cell in reverse order
    /// \relates WDutils::TreeAccess
#ifndef LoopCellKidsReverse
# define LoopCellKidsReverse(CELL,NAME)		\
    for(Cell NAME = this->RBeginCell(CELL);	\
	NAME != this->REndCell(CELL); --NAME)
#endif
    /// loop cell kids of a given cell, starting somewhere
    /// \relates WDutils::TreeAccess
#ifndef LoopCellSecd
# define LoopCellSecd(CELL,START,NAME)		\
    for(Cell NAME = START;			\
	NAME != this->EndCell(CELL); ++NAME)
#endif
    /// loop leaf kids of a given cell
    /// \relates WDutils::TreeAccess
#ifndef LoopLeafKids
# define LoopLeafKids(CELL,NAME)			\
    for(Leaf NAME = this->BeginLeaf(CELL);		\
	NAME != this->EndLeafKids(CELL); ++NAME)
#endif
    /// loop leaf kids of a given cell, starting somewhere
    /// \relates WDutils::TreeAccess
#ifndef LoopLeafSecd
# define LoopLeafSecd(CELL,START,NAME)			\
    for(Leaf NAME = START;				\
	NAME != this->EndLeafKids(CELL); ++NAME)
#endif
    /// loop leaf descendants of a given cell
    /// \relates WDutils::TreeAccess
#ifndef LoopAllLeafs
# define LoopAllLeafs(CELL,NAME)			\
    for(Leaf NAME = this->BeginLeaf(CELL);		\
	NAME != this->EndLeafDesc(CELL); ++NAME)
#endif
    /// loop leaf descendants of a given cell, starting somewhere
    /// \relates WDutils::TreeAccess
#ifndef LoopSecLeafs
# define LoopSecLeafs(CELL,START,NAME)			\
    for(Leaf NAME = START;				\
	NAME != this->EndLeafDesc(CELL); ++NAME)
#endif
    /// loop all except the last leaf descendants of a given cell
    /// \relates WDutils::TreeAccess
#ifndef LoopLstLeafs
# define LoopLstLeafs(CELL,NAME)			\
    for(Leaf NAME = this->BeginLeaf(CELL);		\
	NAME != this->LastLeafDesc(CELL); ++NAME)
#endif
    //@}
    /// \name dumping leaf and cell data (e.g. for debugging purposes)
    //@{
  protected:
    /// header for leaf dump
    /// \note virtual: may be overridden/extended in derived class
    virtual std::ostream&Head(Leaf, std::ostream&out) const
    {
      return out << "   Leaf      I                     X            "
		 << " up     " ;
    }
    /// dump leaf data
    /// \note virtual: may be overridden/extended in derived class
    virtual std::ostream&Data(Leaf l, std::ostream&out) const
    {
      return out << 'L'  << std::setfill('0') << std::setw(6) << l.I
		 << ' '  << std::setfill(' ') << std::setw(6) << particle(l)
		 << ' '                       << std::setw(10)<< position(l)
		 << " C" << std::setfill('0') << std::setw(6) << TREE->PC[l.I];
    }
    /// header for cell dump
    /// \note virtual: may be overridden/extended in derived class
    virtual std::ostream&Head(Cell, std::ostream&out) const
    {
      return out << "Cell    le up     oc Cf     Nc Lf      Nl      N "
		 << "       R                 X          ";
    }
    /// dump cell data
    /// \note virtual: may be overridden/extended in derived class
    virtual std::ostream&Data(Cell c, std::ostream&out) const
    {
      out  << 'C' << std::setfill('0') << std::setw(6) << c.I <<' '
	   << std::setw(2) << std::setfill(' ') << int(level(c)) <<' ';
      if(c.I > 0u)
	out<< 'C' << std::setfill('0') << std::setw(6) << TREE->PA[c.I] <<' ';
      else
	out<< "nil     ";
      out  << int(octant(c)) << ' ';
      if(NCell(c))
	out<< 'C' << std::setfill('0') << std::setw(6)
	   << firstcellindex(c) << ' '
	   << std::setfill(' ') << std::setw(1) << int(NCell(c)) << ' ';
      else
	out<< "nil     0 ";
      return 
	out<< 'L' << std::setfill('0') << std::setw(6)
	   << firstleafindex(c) << ' ' 
	   << std::setfill(' ')
	   << std::setw(2) << NLeafKids(c) << ' '
	   << std::setw(6) << Number   (c) << ' '
	   << std::setw(8) << radius   (c) << ' '
	   << std::setw(8) << centre   (c);
    }
  public:
    /// dump leaf data
    /// \param[in] out  ostream to write to
    void DumpLeafs(std::ostream&out) const
    {
      Head(Leaf(0),out) << '\n';
      LoopLeafs(L) Data(L,out) << '\n';
      out.flush();
    }
    /// dump cell data
    /// \param[in] out  ostream to write to
    void DumpCells(std::ostream&out) const
    {
      Head(Cell(0),out) << '\n';
      LoopCellsDown(C) Data(C,out) << '\n';
      out.flush();
    }
    //@}
  };// struct TreeAccess<>

  ///
  /// The standard tree-walking algorithm.
  ///
  /// To use it, create a class derived from this and provide the abstract
  /// functions.
  template<int __D, typename __X>
  struct TreeWalkAlgorithm : public TreeAccess<OctalTree<__D,__X> > {
    typedef TreeAccess<OctalTree<__D,__X> > Base;
    typedef typename Base::Tree Tree;
    typedef typename Base::Leaf Leaf;
    typedef typename Base::Cell Cell;
    /// \name interaction interface 
    //@{
    /// perform an interaction with a single leaf
    /// \param[in] L  leaf to interact with
    virtual void interact(Leaf L) = 0;
    /// try to perform an interaction with a cell
    /// \param[in] C  cell to interact with
    /// \return  whether interaction was possible, otherwise cell will be split
    virtual bool interact(Cell C) = 0;
    /// perform all leaf interactions with a set of leafs
    /// \param[in] Li  first leaf to interact with
    /// \param[in] Ln  beyond last leaf to interact with
    /// \note Not abstract, but virtual; default: calls single leaf interact
    ///       on each leaf (perhaps somewhat inefficient).
    virtual void interact_many(Leaf Li, Leaf Ln)
    { for(; Li<Ln; ++Li) interact(Li); }
    //@}
    /// ctor: allocate memory for stack
    TreeWalkAlgorithm(const Tree*t)
      : Base(t), S(Tree::Nsub*Base::Depth()) {}
    /// virtual dtor
    virtual ~TreeWalkAlgorithm() {}
    /// perform a tree walk.
    void walk()
    {
      if(!interact(Cell(0)))
	S.push(Cell(0));
      while(!S.is_empty()) {
	Cell C = S.pop();
	if(NLeafKids(C))
	  interact_many(BeginLeaf(C),EndLeafKids(C));
	if(NCell(C))
	  LoopCellKidsReverse(C,Ci)
	    if(!interact(C))
	      S.push(C);
      }
    }
    //
  private:
    Stack<Cell> S;                    ///< stack of cells to interact with
  };

  ///
  /// The mutual tree-walking interaction algorithm of Dehnen (2002, JCP,179,27)
  ///
  /// We implement an "early-testing" mutual tree walk, which means that we try
  /// to perform any interaction as soon as it is generated and only stack it
  /// if it needs splitting. Consequently, any interaction taken from stack is
  /// splitted without further ado. This is usually faster than "late testing".
  ///
  /// \note fully inline
  template<int __D, typename __X>
  struct MutualInteractionAlgorithm : public TreeAccess<OctalTree<__D,__X> > {
  public:
    typedef TreeAccess<OctalTree<__D,__X> > Base;
    typedef typename Base::Tree Tree;
    typedef typename Base::Leaf Leaf;
    typedef typename Base::Cell Cell;
    /// \name interaction interface 
    //@{
    /// perform single leaf-leaf interaction
    virtual void interact(Leaf, Leaf) = 0;
    /// try to perform cell-leaf interaction, return true if success
    virtual bool interact(Cell, Leaf) = 0;
    /// try to perform cell-cell interaction, return true if success
    virtual bool interact(Cell, Cell) = 0;
    /// try to perform cell self-interaction, return true if success
    virtual bool interact(Cell) = 0;
    /// which of two cells of an interaction to split?
    /// \note Not abstract, but virtual; default: split cell with more leafs
    virtual bool split_left(Cell A, Cell B) const
    { return Number(A) > Number(B); }
    /// perform all leaf-leaf interactions between a set of leafs
    /// \note Not abstract, but virtual; default: calls individual leaf-leaf
    ///       interact N*(N-1)/2 times (probably very inefficient).
    virtual void interact_many(Leaf li, Leaf ln)
    { for(; li<ln; ++li) for(Leaf lj=next(li); lj<ln; ++lj) interact(li,lj); }
    /// perform all leaf-leaf interactions between one left and some right
    /// \note Not abstract, but virtual; default: calls individual leaf-leaf
    ///       interact N times (probably not too efficient)
    virtual void interact_many(Leaf ll, Leaf lr, Leaf ln)
    { for(; lr<ln; ++lr) interact(ll,lr); }
    //@}
    /// ctor: allocate memory for stacks
    /// \param[in] t  (pter to) tree
    MutualInteractionAlgorithm(const Tree*t)
      : Base ( t ),
	CL   ( Tree::Nsub*Base::Depth()),
	CC   ( 2*(Tree::Nsub-1)*(Base::Depth()+1)+1 )
    {}
    //  virtual dtor (required by some old compilers)
    virtual ~MutualInteractionAlgorithm() {}
    /// perform a mutual tree walk.
    void walk()
    {
      CC.reset();
      CL.reset();
      perform(Cell(0));
      clear_CC_stack();
    }
  private:
    //
    typedef std::pair<Cell,Leaf> pCL;   ///< represents a cell-leaf interaction
    typedef std::pair<Cell,Cell> pCC;   ///< represents a cell-cell interaction
    //
    Stack<pCL> CL;                      ///< stack of cell-leaf interactions
    Stack<pCC> CC;                      ///< stack of cell-cell interactions
    /// perform a cell-leaf interaction: try to interact, then stack it
    void perform(Cell A, Leaf B)
    { if(!interact(A,B)) CL.push(pCL(A,B)); }
    /// perform a cell-cell interaction: try to interact, then stack it
    void perform(Cell A, Cell B)
    { if(!interact(A,B)) CC.push(pCC(A,B)); }
    /// perform a cell self-interaction: try to interact, then stack it
    void perform(Cell A)
    { if(!interact(A)) CC.push(pCC(A,Cell(0u))); }
    /// clear the stack of cell-leaf interactions
    void clear_CL_stack()
    {
      while(!CL.is_empty()) {
	pCL p = CL.pop();
	if(NLeafKids(p.first))
	  interact_many(p.second,BeginLeaf(p.first),EndLeafKids(p.first));
	if(NCell(p.first))
	  LoopCellKids(p.first,Ci)
	    perform(Ci,p.second);
      }
    }
    /// split a mutual cell-cell interaction
    /// \param[in] A cell to be split
    /// \param[in] B cell to be kept
    void split(Cell A, Cell B)
    {
      if(NLeafKids(A))
	LoopLeafKids(A,Li)
	  perform(B,Li);
      if(NCell(A))
	LoopCellKids(A,Ci)
	  perform(Ci,B);
    }
    /// split a cell self-interaction
    /// \param[in] A cell to be split
    void split(Cell A)
    {
      // leaf-leaf sub-interactions
      if(NLeafKids(A) > 1u)
	interact_many(BeginLeaf(A),EndLeafKids(A));
      // self-interactions between sub-cells
      if(NCell(A)) {
	LoopCellKids(A,Ci)
	  perform(Ci);
      // mutual interactions between sub-cells and sub-leafs
	LoopCellKids(A,Ci) {
	  LoopCellSecd(A,next(Ci),Cj)
	    perform(Ci,Cj);
	  LoopLeafKids(A,Li)
	    perform(Ci,Li);
	}
      }
    }
    /// clear the stack of cell-cell interactions, keep cell-leaf stack clear
    void clear_CC_stack()
    {
      while(!CC.is_empty()) {
	pCC p = CC.pop();
	if     (Cell(0) == p.second)          split(p.first);
	else if(split_left(p.first,p.second)) split(p.first,p.second);
	else                                  split(p.second,p.first);
	clear_CL_stack();
      }
    }
  };
} // namespace WDutils
#endif // WDutils_included_octtree_h
