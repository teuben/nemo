// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/octtree.h
///
/// \brief  contains definition of class \a OctalTree<D,X> and 
///         of class \s MutualOctTreeWalker<TREE> for mutual tree walking
///
/// \author Walter Dehnen
///
/// \date   2009
///
/// \version        2009 WD  based on falcON's tree.h and interact.h
/// \version    May-2009 WD  has been tested and seems okay (May 2009)
/// \version 14-Oct-2009 WD  new design avoiding classes for leafs and cells
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Walter Dehnen
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

#ifndef WDutils_included_iomanip
#  include <iomanip>
#  define WDutils_included_iomanip
#endif
#ifndef WDutils_included_memory_h
#  include <memory.h>
#endif
#ifndef WDutils_included_traits_h
#  include <traits.h>
#endif
#ifndef WDutils_included_tupel_h
#  include <tupel.h>
#endif

namespace { template<int, typename, typename> struct BoxDotTree; }
namespace WDutils {
  //
  /// A spatial tree of square (2D) or cubic (3D) cells
  //
  /// Cells and Leafs are only represented by indices, wrapped in member
  /// classes Cell and Leaf (access to the leaf and cell data via member
  /// methods of class OctalTree).\n
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
  /// \note type __X for positions, type __F for other foating point variables
  /// \note implementations for __X,__F = float,double and Dim=2,3
  template<int __D, typename __X=float, typename __F=__X>
  class OctalTree {
    friend struct BoxDotTree<__D,__X,__F>;
    /// ensure that the only valid instantinations are those implemented in
    /// octtree.cc
    WDutilsStaticAssert
    (( (  __D == 2 || __D == 3 )              &&
       meta::TypeInfo<__X>::is_floating_point &&
       meta::TypeInfo<__F>::is_floating_point    ));
    // disable default and copy ctor
    OctalTree           (const OctalTree&);  // not implemented
    OctalTree& operator=(const OctalTree&);  // not implemented
  public:
    /// \name public static constants some types
    //@{
    const static int        Dim = __D;       ///< number of dimensions
    const static uint32     Nsub= 1<<Dim;    ///< number of octants per cell
    typedef __X             Real;            ///< floating point type: position
    typedef __F             Float;           ///< floating point type: other
    typedef tupel<Dim,Real> point;           ///< type for positions
    typedef uint32          index_type;      ///< type used to index particles
    /// holds the particle data required to build a tree
    struct Dot {
      point       X; ///< position
      index_type  I; ///< index to identify the associated particle
    };
    /// used to initialize and/or re-initialize Dots in tree building
    class Initialiser {
    public:
      /// virtual dtor: only needed with older compiler versions
      virtual ~Initialiser() {}
      /// initializes @a Dot::I and @a Dot::X
      /// used in tree construction (and possibly in rebuild())
      /// \param[in] D Dot to be initialized
      virtual void Init(Dot*D) const = 0;
      /// re-initializes @a Dot::X, may also change @a Dot::I
      /// will be called in rebuild()
      /// \param[in] D Dot to be re-initialized
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
    void Allocate();                  ///< allocates memory
    const uint32      MAXD;           ///< maximum tree depth
    const bool        AVSPC;          ///< avoid single-parent cells?
    uint32            NMAX;           ///< N_max
    uint32            DEPTH;          ///< tree depth
    //@}
    /// \name leafs and related: types, data, and methods
    //@{
    uint32            NLEAF;          ///< # leafs
    point            *XL;             ///< leaf positions
    index_type       *IL;             ///< index of associated particle
  public:
    /// iterator used for leafs
    struct Leaf {
      uint32 I;
      /// ctor
      explicit Leaf(uint32 i) : I(i) {}
      /// increment
      Leaf& operator++() { ++I; return*this; }
      /// comparison <
      bool operator < (Leaf l) const { return I< l.I; }
      /// comparison <=
      bool operator <=(Leaf l) const { return I<=l.I; }
      /// comparison ==
      bool operator ==(Leaf l) const { return I==l.I; }
      /// comparison !=
      bool operator !=(Leaf l) const { return I!=l.I; }
      /// conversion to leaf's index within tree
      operator uint32() const { return I; }
    };
    /// # leafs
    uint32 const&Nleafs() const
    { return NLEAF; }
    /// first leaf
    Leaf BeginLeafs() const
    { return Leaf(0u); }
    /// end of leafs (beyond last leaf)
    Leaf EndLeafs() const
    { return Leaf(NLEAF); }
    /// is this a valid leaf?
    bool IsValid(Leaf l) const
    { return l.I < NLEAF; }
    /// const access to leaf position
    point const&X(Leaf l) const
    { return XL[l]; }
    /// const access to index of associated particle
    index_type const&I(Leaf l) const
    { return IL[l]; }
    //@}
    /// \name cells and related: data and methods
    //@{
  private:
    uint32            NCELL;          ///< # leafs
    uint8            *LE;             ///< cells' tree level
    uint8            *OC;             ///< cells' octant in parent cell
    point            *XC;             ///< cells' centre (of cube)
    uint32           *L0;             ///< cells' first leaf
    uint16           *NL;             ///< number of cells' daughter leafs
    uint32           *NM;             ///< number of cells' leafs
    uint32           *CF;             ///< difference to first daughter cell
    uint8            *NC;             ///< number of cells' daughter cells
    uint32           *PA;             ///< difference to parent cell
    Real             *RAD;            ///< table: radius[level]
  public:
    /// iterator used for cells
    struct Cell {
      uint32 I;
      /// ctor
      explicit Cell(uint32 i) : I(i) {}
      /// increment
      Cell& operator++() { ++I; return*this; }
      /// decrement
      Cell& operator--() { --I; return*this; }
      /// comparison <
      bool operator < (Cell c) const { return I< c.I; }
      /// comparison <=
      bool operator <=(Cell c) const { return I<=c.I; }
      /// comparison >
      bool operator > (Cell c) const { return I> c.I; }
      /// comparison >=
      bool operator >=(Cell c) const { return I>=c.I; }
      /// comparison ==
      bool operator ==(Cell c) const { return I==c.I; }
      /// comparison !=
      bool operator !=(Cell c) const { return I!=c.I; }
      /// conversion to cell's index within tree
      operator uint32() const { return I; }
    };
    /// # cells
    uint32 const&Ncells() const
    { return NCELL; }
    /// root cell
    Cell Root() const
    { return Cell(0u); }
    /// first cell
    Cell BeginCells() const
    { return Cell(0u); }
    /// end of cells (beyond last cell)
    Cell EndCells() const
    { return Cell(NCELL); }
    /// first cell in reversed order: last cell
    Cell RBeginCells() const
    { return Cell(NCELL-1); }
    /// end cell in reversed order: invalid Cell
    Cell REndCells() const
    { return --(Cell(0)); }
    /// is a cell index valid (refers to an actual cell)?
    bool IsValid(Cell c) const
    { return c.I < NCELL; }
    /// tree level of cell
    uint8 const&Le(Cell c) const
    { return LE[c.I]; }
    /// radius (half-side-length of box) of cell
    Real const&Rad(Cell c) const
    { return RAD[LE[c.I]]; }
    /// octant of cell in parent
    uint8 const&Oc(Cell c) const
    { return OC[c.I]; }
    /// cell's geometric centre (of cubic box)
    point const&X(Cell c) const
    { return XC[c.I]; }
    /// first of cell's leafs
    Leaf BeginLeafs(Cell c) const
    { return Leaf(L0[c.I]); }
    /// end of cell's leaf children 
    Leaf EndLeafKids(Cell c) const
    { return Leaf(L0[c.I]+NL[c.I]); }
    /// end of all of cell's leafs
    Leaf EndLeafDesc(Cell c) const
    { return Leaf(L0[c.I]+NM[c.I]); }
    /// number of cell's leaf children
    uint16 const&NumLeafKids(Cell c) const
    { return NL[c.I]; }
    /// total number of cell's leafs
    uint32 const&Number(Cell c) const
    { return NM[c.I]; }
    /// has this cell any daughter leafs?
    bool HasLeafKids(Cell c) const
    { return NL[c.I] > 0u; }
    /// number of daughter cells
    uint8 const&NumCells(Cell c) const
    { return NC[c.I]; }
    /// first of cell's daughter cells
    Cell BeginCells(Cell c) const
    { return Cell(CF[c.I]? c.I+CF[c.I]:NCELL); }
    /// end of cell's daughter cells
    Cell EndCells(Cell c) const
    { return Cell(CF[c.I]? c.I+CF[c.I]+NC[c.I]:NCELL); }
    /// cell's parent cell
    Cell Parent(Cell c) const
    { return Cell(PA[c.I]? c.I-PA[c.I]:NCELL); }
    /// does this cell contain a certain leaf
    bool Contains(Cell c, Leaf l) const 
    { return l >= BeginLeafs(c) && l < EndLeafDesc(c); }
    //@}
    /// ctor: build octtree from scratch.
    ///
    /// The tree is build in two stages. First, a `box-dot' tree is built
    /// using an algorithm which adds one dot (representing a particle) at a
    /// time. Second, this tree is linked to leafs and cells such that 
    /// leaf and cell descendants of any cell are contiguous in memory.
    ///
    /// \param[in] n    number of particles
    /// \param[in] init Initialiser to set initial particle position and index
    /// \note Calls Initialiser::Init(), to set Dot::I and Dot::X for all \a n
    ///       particles.
    /// \param[in] nmax maximum number of leafs in unsplit cells
    /// \param[in] avsc avoid single-parent cells?
    /// \note Single-parent cells occur if all leafs of a cell live in just
    ///       one octant. When \a avsc is set to true, such a cell is
    ///       eliminated in favour of its only daughter. With this option on,
    ///       mother and daughter cells may be more than one level apart and
    ///       the tree depth may be less than the highest cell level.
    /// \param[in] maxd maximum tree depth
    OctalTree(uint32 n, const Initialiser*init, uint32 nmax,
	      bool avsc=true, uint32 maxd=100) WDutils_THROWING;
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
    /// \note Calls Initialiser::ReInit(), to set Dot::X. If \a n is larger
    ///       than previously, Initialiser::Init() is called on the extra
    ///       particles to set Dot::I as well as Dot::X.
    void rebuild(uint32 n=0, uint32 nmax=0) WDutils_THROWING;
    /// tree depth
    uint32 const&Depth() const { return DEPTH; }
    /// are single-parent cells avoided?
    bool const&AvoidedSingleParentCells() const { return AVSPC; }
    /// N_max
    uint32 const&Nmax() const { return NMAX; }
    /// \name dump leaf or cell data to output (for debugging purposes)
    //@{
    /// class whose member specify the data dumped
    struct DumpTreeData {
      DumpTreeData() {}
      virtual ~DumpTreeData() {} // make gcc version 4.1.0 happy
      virtual
      void Head(Leaf, std::ostream&out) const {
	out << " Leaf    I                     X           ";
      }
      virtual
      void Data(Leaf L, const OctalTree*T, std::ostream&out) const {
	out << 'L' << std::setfill('0') << std::setw(4) << L.I
	    << ' ' << std::setfill(' ') << std::setw(4) << T->I(L) << ' '
	    << std::setw(10) << T->X(L);
      }
      virtual
      void Head(Cell, std::ostream&out) const {
	out << "Cell  le   up    Cf Nc  Lf   Nl     N "
	    << "       R                 X         ";
      }
      virtual
      void Data(Cell C, const OctalTree*T, std::ostream&out) const {
	out  << 'C' << std::setfill('0') << std::setw(4) << C.I <<' '
	     << std::setw(2) << std::setfill(' ') << int(T->Le(C)) <<' ';
	if(C > Cell(0))
	  out<< 'C' << std::setfill('0') << std::setw(4) << T->Parent(C) <<' ';
	else
	  out<< "  nil ";
	if(T->NumCells(C))
	  out<< 'C' << std::setfill('0') << std::setw(4)
	     << static_cast<uint32>(T->BeginCells(C)) <<' '
	     << std::setfill(' ') << std::setw(1) << int(T->NumCells(C)) << ' ';
	else
	  out<< " nil 0 ";
	out  << 'L' << std::setfill('0') << std::setw(4)
	     << static_cast<uint32>(T->BeginLeafs(C))
	     << ' ' << std::setfill(' ')
	     << std::setw(2) << T->NumLeafKids(C) << ' '
	     << std::setw(5) << T->Number(C) << ' '
	     << std::setw(8) << T->Rad(C) << ' '
	     << std::setw(8) << T->X(C);
      }
    } DUMP;  ///< no memory required, purely abstract data member
    /// dump leaf data
    /// \param[in] out   ostream to write to
    /// \param[in] dump  pointer to DumpTreeData or derived
    void DumpLeafs(std::ostream&out, const DumpTreeData*dump=0) const {
      if(dump==0) dump=&DUMP;
      dump->Head(BeginLeafs(),out);
      out <<'\n';
      for(Leaf L=BeginLeafs(); L!=EndLeafs(); ++L) {
	dump->Data(L,this,out);
	out <<'\n';
      }
      out.flush();
    }
    /// dump cell data
    /// \param[in] out   ostream to write to
    /// \param[in] dump  pointer to DumpTreeData or derived
    void DumpCells(std::ostream&out, const DumpTreeData*dump=0) const
    {
      if(dump==0) dump=&DUMP;
      dump->Head(BeginCells(),out);
      out <<'\n';
      for(Cell C=BeginCells(); C!=EndCells(); ++C) {
	dump->Data(C,this,out);
	out <<'\n';
      }
      out.flush();
    }
    //@}
  };
  /// running index of a tree leaf
  template<int __D, typename __X, typename __F>
  inline uint32 No(typename OctalTree<__D,__X,__F>::Leaf l) { return l.I; }
  /// running index of a tree cell
  template<int __D, typename __X, typename __F>
  inline uint32 No(typename OctalTree<__D,__X,__F>::Cell c) { return c.I; }
  /// \name macros for looping leafs and cells in an OctalTree    
  //@{
  /// loop cells down: root first
#ifndef LoopCellsDown
# define LoopCellsDown(TREE,NAME)		\
  for(Cell *NAME = (TREE)->BeginCells();	\
      NAME != (TREE)->EndCells(); ++NAME)
#endif
  /// loop cells up: root last
  /// \note useful for an up-ward pass, e.g. computation of centre of mass
#ifndef LoopCellsUp
# define LoopCellsUp(TREE,NAME)			\
  for(Cell *NAME = (TREE)->RBeginCells();	\
      NAME != (TREE)->REndCells(); --NAME)
#endif
  /// loop leafs
#ifndef LoopLeafs
# define LoopLeafs(TREE,NAME)			\
  for(Leaf *NAME = (TREE)->BeginLeafs();	\
      NAME != (TREE)->EndLeafs(); ++NAME)
#endif
  /// loop cell kids of a given cell
#ifndef LoopCellKids
# define LoopCellKids(TREE,CELL,NAME)		\
  for(Cell *NAME = (TREE)->BeginCells(CELL);	\
      NAME != (TREE)->EndCells(CELL); ++NAME)
#endif
  /// loop cell kids of a given cell, starting somewhere
#ifndef LoopCellSecd
# define LoopCellSecd(TREE,CELL,START,NAME)	\
  for(Cell *NAME = START;			\
      NAME != (TREE)->EndCells(CELL); ++NAME)
#endif
  /// loop leaf kids of a given cell
#ifndef LoopLeafKids
# define LoopLeafKids(TREE,CELL,NAME)		\
  for(Leaf *NAME = (TREE)->BeginLeafs(CELL);	\
      NAME != (TREE)->EndLeafKids(CELL); ++NAME)
#endif
  /// loop leaf kids of a given cell, starting somewhere
#ifndef LoopLeafSecd
# define LoopLeafSecd(TREE,CELL,START,NAME)	\
  for(Leaf *NAME = START;			\
      NAME != (TREE)->EndLeafKids(CELL); ++NAME)
#endif
  /// loop leaf descendants of a given cell
#ifndef LoopAllLeafs
# define LoopAllLeafs(TREE,CELL,NAME)		\
  for(Leaf *NAME = (TREE)->BeginLeafs(CELL);	\
      NAME != (TREE)->EndLeafDesc(CELL); ++NAME)
#endif
  /// loop leaf descendants of a given cell, starting somewhere
#ifndef LoopSecLeafs
# define LoopSecLeafs(TREE,CELL,START,NAME)	\
  for(Leaf *NAME = START;			\
      NAME != (TREE)->EndLeafDesc(CELL); ++NAME)
#endif
  /// loop all except the last leaf descendants of a given cell
#ifndef LoopLstLeafs
# define LoopLstLeafs(TREE,CELL,NAME)		\
  for(Leaf *NAME = (TREE)->BeginLeafs(CELL);	\
      NAME != (TREE)->LastLeafDesc(CELL); ++NAME)
#endif
  //@}

  ///
  /// A mutual walk of an OctalTree
  ///
  /// We implement an "Early-testing" mutual tree walk, which means that we try
  /// to perform any interaction as soon as it is generated and only stack it
  /// if it needs splitting. Consequently, any interaction taken from stack is
  /// splitted without further ado. This is faster than "late testing".
  ///
  /// \note fully inline
  template<typename OctTree>
  class MutualOctTreeWalker {
  public:
    typedef typename OctTree::Leaf Leaf;
    typedef typename OctTree::Cell Cell;
    /// specifies how to walk the tree and what kind of interactions to do.
    class Interactor {
      friend class MutualOctTreeWalker;
    protected:
      const OctTree*TREE;
      Interactor(const OctTree*t) : TREE(t) {}
    public:
      /// perform single leaf-leaf interaction
      virtual void interact(Leaf, Leaf) = 0;
      /// try to perform cell-leaf interaction, return true if success
      virtual bool interact(Cell, Leaf) = 0;
      /// try to perform cell-cell interaction, return true if success
      virtual bool interact(Cell, Cell) = 0;
      /// try to perform cell self-interaction, return true if success
      virtual bool interact(Cell) = 0;
      /// which of two cells of an interaction to split?
      virtual bool split_left(Cell A, Cell B) const {
	return TREE->Number(A) > TREE->Number(B);
      }
      /// perform all leaf-leaf interactions between a set of leafs
      /// \note Not abstract, but virtual: default calls individual leaf-leaf
      ///       interact N*(N-1)/2 times (probably very inefficient).
      virtual void interact_many(Leaf L0, Leaf LN) {
	for(Leaf Li=L0; Li<LN; ++Li)
	  for(Leaf Lj=Li+1; Lj<LN; ++Lj)
	    interact(Li,Lj);
      }
    };
  private:
    //
    typedef std::pair<Cell,Leaf> pCL;   ///< represents a cell-leaf interaction
    typedef std::pair<Cell,Cell> pCC;   ///< represents a cell-cell interaction
    //
    const OctTree   *TREE;              ///< tree to walk
    Interactor*const IA;                ///< pter to interactor
    Stack<pCL>       CL;                ///< stack of cell-leaf interactions
    Stack<pCC>       CC;                ///< stack of cell-cell interactions
    //
    void perform(Cell A, Leaf B) { if(!IA->interact(A,B)) CL.push(pCL(A,B)); }
    void perform(Cell A, Cell B) { if(!IA->interact(A,B)) CC.push(pCC(A,B)); }
    void perform(Cell A)         { if(!IA->interact(A))   CC.push(pCC(A,0)); }
    /// clear the stack of cell-leaf interactions
    void clear_CL_stack()
    {
      while(!CL.is_empty()) {
	pCL P = CL.pop();
	if(TREE->HasLeafKids(P.first))
	  LoopLeafKids(TREE,P.first,Li)
	    IA->interact(Li,P.second);
	if(TREE->NumCells(P.first))
	  LoopCellKids(TREE,P.first,Ci)
	    perform(Ci,P.second);
      }
    }
    /// split a mutual cell-cell interaction
    /// \param[in] A cell to be split
    /// \param[in] B cell to be kept
    void split(Cell*A, Cell*B)
    {
      if(TREE->HasLeafKids(A))
	LoopLeafKids(TREE,A,Li)
	  perform(B,Li);
      if(TREE->NumCells(A))
	LoopCellKids(TREE,A,Ci)
	  perform(Ci,B);
    }
    /// split a cell self-interaction
    /// \param[in] A cell to be split
    void split(Cell*A)
    {
      // leaf-leaf sub-interactions
      if(TREE->NumLeafKids(A) > 1)
	IA->interact_many(TREE->BeginLeafs(A),TREE->EndLeafKids(A));
      // self-interactions between sub-cells
      if(TREE->NumCells(A)) {
	LoopCellKids(TREE,A,Ci)
	  perform(Ci);
      // mutual interaction between sub-cells and sub-leafs
	LoopCellKids(TREE,A,Ci) {
	  LoopCellSecd(TREE,A,Ci+1,Cj)
	    perform(Ci,Cj);
	  LoopLeafKids(TREE,A,Li)
	    perform(Ci,Li);
	}
      }
    }
    /// clear the stack of cell-cell interactions, keep cell-leaf stack clear
    void clear_CC_stack()
    {
      while(!CC.is_empty()) {
	pCC P = CC.pop();
	if     (0 == P.second)                    split(P.first);
	else if(IA->split_left(P.first,P.second)) split(P.first,P.second);
	else                                      split(P.second,P.first);
	clear_CL_stack();
      }
    }
  public:
    /// construction for interaction within one octal tree
    /// \param[in] i  pter to interactor
    /// \param[in] d  depth of tree
    MutualOctTreeWalker(Interactor*i)
      : TREE ( i.TREE ),
	IA   ( i ), 
	CL   ( OctTree::Nsub*TREE->Depth()),
	CC   ( 2*(OctTree::Nsub-1)*(TREE->Depth()+1)+1 )
    {}
    /// perform a mutual tree walk.
    void walk()
    {
      CC.reset();
      CL.reset();
      perform(TREE->Root());
      clear_CC_stack();
    }
  };
  /// walk an oct-tree: simple wrapper around MutualTreeWalker.
  template<typename OctTree>
  void MutualOctTreeWalk(typename MutualOctTreeWalker<OctTree>::Interactor*I)
  {
    MutualOctTreeWalker<OctTree> MTW(I);
    MTW.walk();
  }
} // namespace WDutils
////////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_octtree_h
