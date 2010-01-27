// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/octtree.h
///
/// \brief  contains definition of classes @a OctalTree<D,X>,
///         @a TreeWalker<OctTree>, and @a MutualInteractionAlgorithm<OctTree>
///
/// \author Walter Dehnen
///
/// \date   2009,2010
///
/// \version May-2009 WD  first tested version
/// \version Oct-2009 WD  new design using indices for leafs and cells
/// \version Nov-2009 WD  removed redundant template parameter
/// \version Jan-2010 WD  renamed methods in TreeWalker; added leaf's parent
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

#ifndef WDutils_included_fstream
#  include <fstream>
#  define WDutils_included_fstream
#endif
#ifndef WDutils_included_iomanip
#  include <iomanip>
#  define WDutils_included_iomanip
#endif
#ifndef WDutils_included_memory_h
#  include <memory.h>
#endif
#ifndef WDutils_included_tupel_h
#  include <tupel.h>
#endif

namespace { template<int, typename> struct BoxDotTree; }
namespace WDutils {
  template<typename> struct TreeWalker;
  template<typename> struct DumpTreeData;
  //
  /// A spatial tree of square (2D) or cubic (3D) cells
  //
  /// Cells with more than \a nmax (argument to constructor and member
  /// rebuild) are split and octants (or quarters for Dim=2) with more than
  /// one particle are assigned a new cell. After tree construction, the
  /// octants are dissolved (though each cell still knows its octant in its
  /// parent cell) and any leafs from single-leaf octants become direct
  /// leaf-children of their cell. Thus, a non-final cell (cell with daughter
  /// cells) can also contain leafs which are in none of their daughter cells;
  /// these leafs are placed first in the array of a cell's leafs and are
  /// referred to as 'leaf kids' as opposed to 'leaf descendants', which
  /// includes all leafs contained within a cell.
  ///
  /// \note Access to leafs and cells via struct TreeWalker<OctalTree> below
  /// \note Implementations for @a __X = float,double and @a __D = 2,3
  template<int __D, typename __X>
  class OctalTree {
    friend struct BoxDotTree<__D,__X>;
    friend struct TreeWalker<OctalTree>;
    /// ensure that the only valid instantinations are those implemented in
    /// octtree.cc
    WDutilsStaticAssert
    (( (  __D == 2 || __D == 3 )              &&
       meta::TypeInfo<__X>::is_floating_point    ));
    // disable default and copy ctor
    OctalTree           (const OctalTree&);  // not implemented
    OctalTree& operator=(const OctalTree&);  // not implemented
  public:
    /// \name public constants and types
    //@{
    const static int Dim = __D;              ///< number of dimensions
    const static int Nsub= 1<<Dim;           ///< number of octants per cell
    typedef __X             Real;            ///< floating point type: position
    typedef tupel<Dim,Real> Point;           ///< type: positions
    typedef uint32          particle_index;  ///< type: indexing particles
    typedef uint32          node_index;      ///< type: indexing leafs & cells
    typedef uint32          depth_type;      ///< type: tree depth & related
    /// holds the particle data required to build a tree
    //
    //  Not to be messed with! this data type MUST match DotBase<__D,__X> in
    //  octtree.cc for the code to work.
    struct Dot {
      Point           X;             ///< position
      particle_index  I;             ///< identifier of associated particle
    };
    /// used to initialize and/or re-initialize Dots in tree building
    class Initialiser {
    public:
      // virtual dtor: only needed with older compiler versions
      virtual ~Initialiser() {}
      /// initializes @a Dot::I and @a Dot::X
      /// \note used in tree construction (and possibly in OctalTree::rebuild())
      /// \param[in] D  Dot to be initialized
      virtual void Init(Dot*D) const = 0;
      /// re-initializes @a Dot::X, may also change @a Dot::I
      /// \note will be called in OctalTree::rebuild()
      /// \param[in] D  Dot to be re-initialized
      virtual void ReInit(Dot*D) const = 0;
    };
    //@}
    /// \name general data for class OctalTree
    //@{
  protected:
    const Initialiser*INIT;           ///< initialising
  private:
    char*             ALLOC;          ///< actually allocated memory
    unsigned          NALLOC;         ///< # bytes allocated at ALLOC
    const depth_type  MAXD;           ///< maximum tree depth
    const bool        AVSPC;          ///< avoid single-parent cells?
    depth_type        NMAX;           ///< N_max
    depth_type        DEPTH;          ///< tree depth
    //@}
    /// \name leaf data (access via TreeWalker<OctalTree>)
    //@{
    node_index        NLEAF;          ///< total number of leafs
    Point            *XL;             ///< leaf positions
    particle_index   *PL;             ///< index of associated particle
    node_index       *PC;             ///< index of parent cell
    //@}
    /// \name cell data (access via TreeWalker<OctalTree>)
    //@{
    node_index        NCELL;          ///< total number of cells
    uint8            *LE;             ///< cells' tree level
    uint8            *OC;             ///< cells' octant in parent cell
    Point            *XC;             ///< cells' centre (of cube)
    node_index       *L0;             ///< cells' first leaf
    uint16           *NL;             ///< number of cells' leaf kids
    node_index       *NM;             ///< number of cells' leaf descendants
    node_index       *CF;             ///< first daughter cell
    uint8            *NC;             ///< number of cells' daughter cells
    node_index       *PA;             ///< parent cell
    Real             *RAD;            ///< table: radius[level]
    //@}
    void Allocate();                  ///< allocates memory
  public:
    /// ctor: build octtree from scratch.
    ///
    /// The tree is build in two stages. First, a `box-dot' tree is built
    /// using an algorithm which adds one dot (representing a particle) at a
    /// time. Second, this tree is linked to leafs and cells such that 
    /// leaf and cell descendants of any cell are contiguous in memory.
    ///
    /// \param[in] n    number of particles
    /// \param[in] init Initialiser to set initial particle position and index
    /// \note Calls Initialiser::Init(), to set Dot::I and Dot::X for all @a n
    ///       particles.
    /// \param[in] nmax maximum number of leafs in unsplit cells
    /// \param[in] avsc avoid single-parent cells?
    /// \note Single-parent cells occur if all leafs of a cell live in just
    ///       one octant. When @a avsc is true, such a cell is eliminated in
    ///       favour of its only daughter. With this option on, mother and
    ///       daughter cells may be more than one level apart and the tree
    ///       depth may be less than the highest cell level.
    /// \param[in] maxd maximum tree depth
    OctalTree(node_index n, const Initialiser*init, depth_type nmax,
	      bool avsc=true, depth_type maxd=100) WDutils_THROWING;
    /// dtor
    ~OctalTree();
    /// build the tree again, re-initialising the dots
    ///
    /// The tree is build exactly in the same way as with the constructor,
    /// only the order in which the dots are added to the `box-dot' tree is
    /// that of the original tree rather than increasing index. This change
    /// alone results in a speed-up by about a factor 2 for the whole process
    /// (including linking the final tree), because it avoids cache misses.
    ///
    /// \param[in] n    new number of particles
    /// \param[in] nmax maximum number of leafs in unsplit cell
    ///
    /// \note If any of the arguments equals 0, we take the old value instead
    /// \note Calls Initialiser::ReInit(), to set Dot::X, but possibly also
    ///       change Dot::I. If @a n is larger than previously,
    ///       Initialiser::Init() is called on the extra particles to set
    ///       Dot::I as well as Dot::X.
    void rebuild(node_index n=0, depth_type nmax=0) WDutils_THROWING;
    /// tree depth
    depth_type const&Depth() const
    { return DEPTH; }
    /// are single-parent cells avoided?
    bool const&AvoidedSingleParentCells() const
    { return AVSPC; }
    /// N_max
    depth_type const&Nmax() const
    { return NMAX; }
    /// total number of leafs
    node_index const&Nleafs() const
    { return NLEAF; }
    /// total number of cells
    node_index const&Ncells() const
    { return NCELL; }
    /// initialiser
    const Initialiser*Init() const
    { return INIT; }
  };// class OctalTree<>

  ///
  /// support for walking an OctalTree.
  ///
  /// Holds just a pointer to an OctalTree and provides access to leaf and
  /// cell data, as well as member methods for walking the tree.\n
  ///
  /// Useful as base class for tree-walking algorithms
  template<typename OctTree>
  struct TreeWalker {
    typedef typename OctTree::Real Real;
    typedef typename OctTree::Point Point;
    typedef typename OctTree::particle_index particle_index;
    typedef typename OctTree::node_index node_index;
    typedef typename OctTree::depth_type depth_type;
    /// pointer to tree
    const OctTree*const TREE;
    /// ctor
    TreeWalker(const OctTree*t) : TREE(t) {}
    /// copy ctor
    TreeWalker(const TreeWalker&t) : TREE(t.TREE) {}
    /// virtual dtor (to make gcc version 4.1.0 happy)
    virtual ~TreeWalker() {}
    /// \name leaf and leaf data access
    //@{
    /// iterator used for leafs.
    /// A simple wrapper around an index, which, being a separate type, avoids
    /// confusion with other indices or variables of node_index.
    struct Leaf {
      node_index I;
      /// default ctor
      Leaf() {}
      /// ctor from index
      explicit Leaf(node_index i)
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
      /// conversion to leaf's index within tree
      operator node_index() const
      { return I; }
    };
    /// leaf position
    Point const&position(Leaf l) const
    { return TREE->XL[l.I]; }
    /// index of associated particle
    particle_index const&particle(Leaf l) const
    { return TREE->PL[l.I]; }
    /// index of parent cell
    particle_index const&parentcellindex(Leaf l) const
    { return TREE->PC[l.I]; }
    //@}
    /// \name cell and cell data access
    //@{
    /// iterator used for cells.
    /// A simple wrapper around an index, which, being a separate type, avoids
    /// confusion with other indices or variables of node_index.
    struct Cell {
      node_index I;
      /// default ctor
      Cell()
      {}
      /// ctor from index
      explicit Cell(node_index i)
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
      /// conversion to cell's index within tree
      operator node_index() const
      { return I; }
    };
    /// tree level of cell
    uint8 const&level(Cell c) const
    { return TREE->LE[c.I]; }
    /// octant of cell in parent
    uint8 const&octant(Cell c) const
    { return TREE->OC[c.I]; }
    /// cell's geometric centre (of cubic box)
    Point const&centre(Cell c) const
    { return TREE->XC[c.I]; }
    /// radius (half-side-length of box) of cell
    Real const&radius(Cell c) const
    { return TREE->RAD[level(c)]; }
    /// number of leaf kids
    uint16 const&Nleafkids(Cell c) const
    { return TREE->NL[c.I]; }
    /// total number of leafs
    node_index const&Nleafs(Cell c) const
    { return TREE->NM[c.I]; }
    /// number of daughter cells
    uint8  const&Ncells(Cell c) const
    { return TREE->NC[c.I]; }
    /// index of cell's first leaf
    /// \note included for completeness, not usually needed in tree walks
    node_index const&firstleafindex(Cell c) const
    { return TREE->L0[c.I]; }
    /// index of cell's first daughter cell, if any
    /// \note included for completeness, not usually needed in tree walks
    node_index const&firstcellindex(Cell c) const
    { return TREE->CF[c.I]; }
    /// index of parent cell, if any
    /// \note included for completeness, not usually needed in tree walks
    node_index const&parentcellindex(Cell c) const
    { return TREE->PA[c.I]; }
    //@}
    /// \name tree walking and related
    //@{
    /// tree depth
    depth_type const&Depth() const
    { return TREE->DEPTH; }
    /// # leafs
    node_index const&Nleafs() const
    { return TREE->NLEAF; }
    /// next leaf
    static Leaf next(Leaf l)
    { return Leaf(l.I+1); }
    /// first leaf
    static Leaf BeginLeafs()
    { return Leaf(0u); }
    /// end of leafs (beyond last leaf)
    Leaf EndLeafs() const
    { return Leaf(TREE->NLEAF); }
    /// is this a valid leaf?
    bool IsValid(Leaf l) const
    { return l.I < TREE->NLEAF; }
    /// # cells
    node_index const&Ncells() const
    { return TREE->NCELL; }
    /// next cell
    static Cell next(Cell c)
    { return Cell(c.I+1); }
    /// first cell
    static Cell BeginCells()
    { return Cell(0u); }
    /// end of cells (beyond last cell)
    Cell EndCells() const
    { return Cell(TREE->NCELL); }
    /// first cell in reversed order: last cell
    Cell RBeginCells() const
    { return Cell(TREE->NCELL-1); }
    /// end cell in reversed order: invalid Cell
    static Cell REndCells()
    { return --(Cell(0)); }
    /// is a cell index valid (refers to an actual cell)?
    bool IsValid(Cell c) const
    { return c.I < TREE->NCELL; }
    /// first of cell's leafs
    Leaf BeginLeafs(Cell c) const
    { return Leaf(firstleafindex(c)); }
    /// end of cell's leaf children 
    Leaf EndLeafKids(Cell c) const
    { return Leaf(firstleafindex(c)+Nleafkids(c)); }
    /// end of all of cell's leafs
    Leaf EndLeafDesc(Cell c) const
    { return Leaf(firstleafindex(c)+Nleafs(c)); }
    /// first of cell's daughter cells
    /// \note if there are no daughter cells, this returns c
    Cell BeginCells(Cell c) const
    { return Cell(firstcellindex(c)); }
    /// end of cell's daughter cells
    /// \note if there are no daughter cells, this returns c
    Cell EndCells(Cell c) const
    { return Cell(firstcellindex(c)+Ncells(c)); }
    /// cell's parent cell
    Cell Parent(Cell c) const
    { return Cell(c.I? parentcellindex(c):TREE->NCELL); }
    /// leaf's parent cell
    Cell Parent(Leaf l) const
    { return Cell(parentcellindex(l)); }
    /// does this cell contain a certain leaf
    bool Contains(Cell c, Leaf l) const
    { return BeginLeafs(c) <= l && l < EndLeafDesc(c); }
    /// is either cell ancestor of the other?
    bool IsAncestor(Cell a, Cell b) const
    { return maxnorm(centre(a)-centre(b)) < max(radius(a),radius(b)); }
    //@}
    /// \name macros for tree walking from within a TreeWalker
    //@{
    /// loop cells down: root first
    /// \note useful for a down-ward pass
#ifndef LoopCellsDown
# define LoopCellsDown(NAME)			\
    for(Cell NAME = this->BeginCells();		\
	NAME != this->EndCells(); ++NAME)
#endif
    /// loop cells up: root last
    /// \note useful for an up-ward pass
#ifndef LoopCellsUp
# define LoopCellsUp(NAME)			\
    for(Cell NAME = this->RBeginCells();	\
	NAME != this->REndCells(); --NAME)
#endif
    /// loop leafs
#ifndef LoopLeafs
# define LoopLeafs(NAME)			\
    for(Leaf NAME = this->BeginLeafs();		\
	NAME != this->EndLeafs(); ++NAME)
#endif
    /// loop cell kids of a given cell
#ifndef LoopCellKids
# define LoopCellKids(CELL,NAME)		\
    for(Cell NAME = this->BeginCells(CELL);	\
    NAME != this->EndCells(CELL); ++NAME)
#endif
    /// loop cell kids of a given cell, starting somewhere
#ifndef LoopCellSecd
# define LoopCellSecd(CELL,START,NAME)		\
    for(Cell NAME = START;			\
	NAME != this->EndCells(CELL); ++NAME)
#endif
    /// loop leaf kids of a given cell
#ifndef LoopLeafKids
# define LoopLeafKids(CELL,NAME)			\
    for(Leaf NAME = this->BeginLeafs(CELL);		\
	NAME != this->EndLeafKids(CELL); ++NAME)
#endif
    /// loop leaf kids of a given cell, starting somewhere
#ifndef LoopLeafSecd
# define LoopLeafSecd(CELL,START,NAME)			\
    for(Leaf NAME = START;				\
	NAME != this->EndLeafKids(CELL); ++NAME)
#endif
    /// loop leaf descendants of a given cell
#ifndef LoopAllLeafs
# define LoopAllLeafs(CELL,NAME)			\
    for(Leaf NAME = this->BeginLeafs(CELL);		\
	NAME != this->EndLeafDesc(CELL); ++NAME)
#endif
    /// loop leaf descendants of a given cell, starting somewhere
#ifndef LoopSecLeafs
# define LoopSecLeafs(CELL,START,NAME)			\
    for(Leaf NAME = START;				\
	NAME != this->EndLeafDesc(CELL); ++NAME)
#endif
    /// loop all except the last leaf descendants of a given cell
#ifndef LoopLstLeafs
# define LoopLstLeafs(CELL,NAME)			\
    for(Leaf NAME = this->BeginLeafs(CELL);		\
	NAME != this->LastLeafDesc(CELL); ++NAME)
#endif
    //@}
    /// \name dumping leaf and cell data (for debugging purposes)
    //@{
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
      return out << 'L' << std::setfill('0') << std::setw(6) << l.I
		 << ' ' << std::setfill(' ') << std::setw(6) << particle(l)
		 << ' ' << std::setw(10) << position(l)
		 << " C" << std::setfill('0') << std::setw(6)
		 << parentcellindex(l) ;
    }
    /// header for cell dump
    /// \note virtual: may be overridden/extended in derived class
    virtual std::ostream&Head(Cell, std::ostream&out) const
    {
      return out << "Cell    le up      Cf     Nc Lf      Nl      N "
		 << "       R                 X          ";
    }
    /// dump cell data
    /// \note virtual: may be overridden/extended in derived class
    virtual std::ostream&Data(Cell c, std::ostream&out) const
    {
      out  << 'C' << std::setfill('0') << std::setw(6) << c.I <<' '
	   << std::setw(2) << std::setfill(' ') << int(level(c)) <<' ';
      if(c.I > 0u)
	out<< 'C' << std::setfill('0') << std::setw(6)
	   << parentcellindex(c) <<' ';
      else
	out<< "nil     ";
      if(Ncells(c))
	out<< 'C' << std::setfill('0') << std::setw(6)
	   << firstcellindex(c) << ' '
	   << std::setfill(' ') << std::setw(1) << int(Ncells(c)) << ' ';
      else
	out<< "nil     0 ";
      return 
	out<< 'L' << std::setfill('0') << std::setw(6)
	   << firstleafindex(c) << ' ' 
	   << std::setfill(' ')
	   << std::setw(2) << Nleafkids(c) << ' '
	   << std::setw(6) << Nleafs   (c) << ' '
	   << std::setw(8) << radius   (c) << ' '
	   << std::setw(8) << centre   (c);
    }
    /// dump leaf data
    /// \param[in] out   ostream to write to
    void DumpLeafs(std::ostream&out) const
    {
      Head(Leaf(0),out) << '\n';
      LoopLeafs(L) Data(L,out) << '\n';
      out.flush();
    }
    /// dump leaf data
    /// \param[in] file name of file to write to
    void DumpLeafs(const char*file) const
    {
      std::ofstream out(file);
      DumpLeafs(out);
    }
    /// dump cell data
    /// \param[in] out   ostream to write to
    void DumpCells(std::ostream&out) const
    {
      Head(Cell(0),out) << '\n';
      LoopCellsDown(C) Data(C,out) << '\n';
      out.flush();
    }
    /// dump cell data
    /// \param[in] file name of file to write to
    void DumpCells(const char*file) const
    {
      std::ofstream out(file);
      DumpCells(out);
    }
    //@}
  };// struct TreeWalker<>

  ///
  /// The mutual tree-walking interaction algorithm of Dehnen (2002)
  ///
  /// We implement an "Early-testing" mutual tree walk, which means that we try
  /// to perform any interaction as soon as it is generated and only stack it
  /// if it needs splitting. Consequently, any interaction taken from stack is
  /// splitted without further ado. This is faster than "late testing".
  ///
  /// \note fully inline
  template<typename OctTree>
  class MutualInteractionAlgorithm : private TreeWalker<OctTree> {
  public:
    typedef typename TreeWalker<OctTree>::Leaf Leaf;
    typedef typename TreeWalker<OctTree>::Cell Cell;
    /// specifies how to walk the tree and what kind of interactions to do.
    class Interactor : public TreeWalker<OctTree> {
    public:
      /// ctor
      Interactor(const OctTree*t) : TreeWalker<OctTree>(t) {}
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
      { return N(A) > N(B); }
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
    };
  private:
    //
    typedef std::pair<Cell,Leaf> pCL;   ///< represents a cell-leaf interaction
    typedef std::pair<Cell,Cell> pCC;   ///< represents a cell-cell interaction
    //
    Interactor*const IA;                ///< pter to interactor
    Stack<pCL>       CL;                ///< stack of cell-leaf interactions
    Stack<pCC>       CC;                ///< stack of cell-cell interactions
    //
    void perform(Cell A, Leaf B)
    { if(!IA->interact(A,B)) CL.push(pCL(A,B)); }
    void perform(Cell A, Cell B)
    { if(!IA->interact(A,B)) CC.push(pCC(A,B)); }
    void perform(Cell A)
    { if(!IA->interact(A)) CC.push(pCC(A,Cell(0u))); }
    /// clear the stack of cell-leaf interactions
    void clear_CL_stack()
    {
      while(!CL.is_empty()) {
	pCL p = CL.pop();
	if(Nl(p.first))
	  IA->interact_many(p.second,BeginLeafs(p.first),EndLeafKids(p.first));
	if(Nc(p.first))
	  LoopCellKids(p.first,Ci)
	    perform(Ci,p.second);
      }
    }
    /// split a mutual cell-cell interaction
    /// \param[in] A cell to be split
    /// \param[in] B cell to be kept
    void split(Cell A, Cell B)
    {
      if(Nl(A))
	LoopLeafKids(A,Li)
	  perform(B,Li);
      if(Nc(A))
	LoopCellKids(A,Ci)
	  perform(Ci,B);
    }
    /// split a cell self-interaction
    /// \param[in] A cell to be split
    void split(Cell A)
    {
      // leaf-leaf sub-interactions
      if(Nl(A) > 1u)
	IA->interact_many(BeginLeafs(A),EndLeafKids(A));
      // self-interactions between sub-cells
      if(Nc(A)) {
	LoopCellKids(A,Ci)
	  perform(Ci);
      // mutual interaction between sub-cells and sub-leafs
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
	if     (Cell(0) == p.second)              split(p.first);
	else if(IA->split_left(p.first,p.second)) split(p.first,p.second);
	else                                      split(p.second,p.first);
	clear_CL_stack();
      }
    }
  public:
    /// construction for interaction within octal tree
    /// \param[in] i  pter to interactor
    MutualInteractionAlgorithm(Interactor*i)
      : TreeWalker<OctTree>(*i),
	IA   ( i ), 
	CL   ( OctTree::Nsub*this->Depth()),
	CC   ( 2*(OctTree::Nsub-1)*(this->Depth()+1)+1 )
    {}
    /// perform a mutual tree walk.
    void walk()
    {
      CC.reset();
      CL.reset();
      perform(Cell(0));
      clear_CC_stack();
    }
  };
  /// walk an oct-tree: simple wrapper around MutualInteractionAlgorithm.
  template<typename OctTree>
  void MutualInteraction(typename MutualInteractionAlgorithm<OctTree>::Interactor*I)
  {
    MutualInteractionAlgorithm<OctTree> MIA(I);
    MIA.walk();
  }
} // namespace WDutils

#endif // WDutils_included_octtree_h
