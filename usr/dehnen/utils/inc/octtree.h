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
/// \note   based on falcON's tree.h and interact.h
/// \note   has been tested and seems okay (May 2009)
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

namespace WDutils {
  //
  /// A spatial tree of square (2D) or cubic (3D) cells
  //
  /// Cells and particles are represented by member types Cell and Leaf,
  /// respectively. Cells with more than \a nmax (argument to constructor and
  /// member rebuild) are split and octants (or quarters for Dim=2) with more
  /// than one particle are assigned a new cell. After tree construction, the
  /// octants are dissolved (though each cell still knows its octant in its
  /// parent cell) and any leafs from single-leaf octants become direct
  /// leaf-children of their cell. Thus, a non-final cell (cell with daughter
  /// cells) can also contain leafs which are in none of their daughter
  /// cells; these leafs are placed first in the array of a cell's leafs and
  /// are referred to as 'leaf kids' as opposed to 'leaf descendants', which
  /// includes all leafs contained within a cell.
  ///
  /// \note type __X for positions, type __F for other foating point variables
  /// \note implementations for __X,__F = float,double and Dim=2,3
  template<int __D, typename __X=float, typename __F=__X>
  class OctalTree {
    /// ensure that the only valid instantinations are those implemented in
    /// octtree.cc
    WDutilsStaticAssert
    (( (  __D == 2 || __D == 3 )              &&
       meta::TypeInfo<__X>::is_floating_point &&
       meta::TypeInfo<__F>::is_floating_point    ));
    // disable default and copy ctor
    OctalTree           (const OctalTree&); // not implemented
    OctalTree& operator=(const OctalTree&); // not implemented
  public:
    /// \name public static constants and types
    //@{
    const static int        Dim = __D;    ///< number of dimensions
    const static unsigned   Nsub= 1<<Dim; ///< number of octants per cell
    typedef __X             Real;         ///< floating point type for position
    typedef __F             Float;        ///< floating point type for other
    typedef tupel<Dim,Real> point;        ///< type for positions
    typedef unsigned        index_type;   ///< type used to index particles
    /// represents the basic relevant particle data
    struct Dot {
      point       X; ///< position
      index_type  I; ///< index to identify the associated particle
    };
    /// represents a particle in the tree
    /// \note not memory critical during tree build.
    /// \note Only @a Dot data @a X and @a I will be set at tree build.
    /// \note On 64bit machines and in 3D, the sizeof(Leaf) is 32 and 48 for
    ///       float and double.
    struct Leaf : public Dot 
    {
      // additional auxiliary data
      union {
	int      aI; ///< auxiliary int
	unsigned aU; ///< auxiliary unsigned
	Real     aR; ///< auxiliary scalar
	void    *aP; ///< auxiliary pointer to more data
      };
      int         F; ///< any int, e.g. bitfield of flags
      Float       M; ///< any scalar, e.g. mass
    };
    /// represents a cubic (square for 2D) cell 
    /// \note not memory critical during tree build.
    /// \note auxiliary data will not be set at tree build.
    struct Cell {
      // data set at tree build (20+sizeof(*void)+Dim*sizeof(Float) bytes)
      uint8       Le; ///< level in tree
      uint8       Oc; ///< octant in parent Cell
      uint16      Rk; ///< rank of mother domain (currently not used)
      uint16      Nc; ///< # cell children
      uint16      Nl; ///< # leaf children (coming first amongst the desc)
      uint32      N;  ///< # leaf descendants
      Leaf       *L0; ///< pointer to first leaf
      uint32      Cf; ///< pointer difference to first cell daughter, if any
      uint32      Pa; ///< pointer difference to parent cell
      point       X;  ///< center position of cube
      // additional auxiliary data
      void       *P;  ///< auxiliary pointer to more data
      point       Z;  ///< auxiliary position, e.g. centre of mass
      Float       R;  ///< auxiliary scalar, e.g. radius
      Float       M;  ///< auxiliary scalar, e.g. mass
      int         F;  ///< auxiliary integer, e.g. bitfield
      union {
	int      aI;  ///< auxiliary int
	unsigned aU;  ///< auxiliary unsigned
	Float    aR;  ///< auxiliary scalar
	void    *aP;  ///< auxiliary pointer to more data
      };
      /// is *this a final cell (has no daughter cells)?
      bool IsFinal() const { return Cf==0u; }
      /// number of cell kids
      uint16 const&NumCells() const { return Nc; }
      /// has this cell leaf kids (direct leaf descendants)?
      bool HasLeafKids() const { return Nl!=0u; }
      /// number of leaf kids (direct leaf descendants)
      uint16 const&NumLeafKids() const { return Nl; }
      /// number of leaf descendants
      uint32 const&Number() const { return N; }
      /// pointer to first leaf descendant (always exists)
      Leaf*BeginLeafs() const { return L0; }
      /// pointer to end of leaf kids
      Leaf*EndLeafKids() const { return L0+Nl; }
      /// pointer to end of leaf descendants
      Leaf*EndLeafDesc() const { return L0+N; }
      /// pointer to last leaf descendants
      Leaf*LastLeafDesc() const { return EndLeafDesc()-1; }
      /// pointer to first daughter cell, if any
      Cell*BeginCells() const { return Cf? const_cast<Cell*>(this)+Cf : 0; }
      /// pointer to last daugher cell, if any
      Cell*EndCells() const { return Cf? const_cast<Cell*>(this)+(Cf+Nc) : 0; }
      /// pointer to parent cell (null pointer if this is the root cell)
      Cell*Parent() const { return Pa? const_cast<Cell*>(this)-Pa : 0; }
      /// does this cell contain a certain leaf?
      bool Contains(const Leaf*L) const {
	return L >= BeginLeafs() && L < EndLeafDesc(); 
      }
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
    /// \name data for class OctalTree
    //@{
  protected:
    const Initialiser*INIT;           ///< initialising
  private:
    const unsigned    MAXD;           ///< maximum tree depth
    const bool        AVSPC;          ///< avoid single-parent cells?
    unsigned          NMAX;           ///< N_max
    unsigned          NLEAF;          ///< # leafs
    unsigned          NCELL;          ///< # cells
    unsigned          DEPTH;          ///< tree depth
    Leaf*             LEAFS;          ///< array of leafs
    Cell*             CELLS;          ///< array of cells
    Real*             RAD;            ///< table: radius[level]
    char*             ALLOC;          ///< actually allocated memory
    unsigned          NALLOC;         ///< # bytes allocated at ALLOC
    //@}
    /// ensures LEAFS, CELLS, RAD point to sufficient memory
    void Allocate();
  public:
    /// \name construction and related
    //@{
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
    OctalTree(unsigned n, const Initialiser*init, unsigned nmax,
	      bool avsc=true, unsigned maxd=100) WDutils_THROWING;
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
    void rebuild(unsigned n=0, unsigned nmax=0) WDutils_THROWING;
    //@}
    /// \name functionality of a build tree
    //@{
    /// root cell
    Cell*const&Root() const { return CELLS; }
    /// first cell
    Cell*const&BeginCells() const { return CELLS; }
    /// end cell
    Cell*EndCells() const { return CELLS+NCELL; }
    /// first cell in reverse order: last cell
    Cell*RBeginCells() const { return CELLS+NCELL-1; }
    /// end cell in reverse order: before first cell
    Cell*REndCells() const { return CELLS-1; }
    /// running number of a given cell
    unsigned NoCell(const Cell*C) const { return C-CELLS; }
    /// running number of a given leaf
    unsigned NoLeaf(const Leaf*L) const { return L-LEAFS; }
    /// first leaf
    Leaf*const&BeginLeafs() const { return LEAFS; }
    /// end of leafs
    Leaf*EndLeafs() const { return LEAFS+NLEAF; }
    /// leaf of given index
    Leaf*LeafNo(int i) const { return LEAFS+i; }
    /// tree depth
    unsigned const&Depth() const { return DEPTH; }
    /// are single-parent cells avoided?
    bool const&AvoidedSingleParentCells() const { return AVSPC; }
    /// N_max
    unsigned const&Nmax() const { return NMAX; }
    /// number of leafs
    unsigned const&Nleafs() const { return NLEAF; }
    /// number of cells
    unsigned const&Ncells() const { return NCELL; }
    /// radius (half-side-length) of cell
    Real Radius(const Cell*C) const { return RAD[C->Le]; }
    /// is A ancestor of B or vice versa
    bool IsAncestor(const Cell*a, const Cell*b) const {
      return maxnorm(a->X - b->X) < max(Radius(a),Radius(b));
    }
    //@}
    /// \name dump leaf or cell data to output (for debugging purposes)
    //@{
    /// class whose member specify the data dumped
    struct DumpTreeData {
      DumpTreeData() {}
      virtual
      void Head(const Leaf*, std::ostream&out) const {
	out << " Leaf    I                     X           ";
      }
      virtual
      void Data(const Leaf*L, const OctalTree*T, std::ostream&out) const {
	out << 'L' << std::setfill('0') << std::setw(4) << T->NoLeaf(L)
	    << ' ' << std::setfill(' ') << std::setw(4) << L->I << ' '
	    << std::setw(10) << L->X;
      }
      virtual
      void Head(const Cell*, std::ostream&out) const {
	out << "Cell le  up   Cf Nc  Lf   Nl    N "
	    << "       R                 X         ";
      }
      virtual
      void Data(const Cell*C, const OctalTree*T, std::ostream&out) const {
	unsigned no = T->NoCell(C);
	out  << 'C' << std::setfill('0') << std::setw(3) << no <<' '
	     << std::setw(2) << std::setfill(' ') << int(C->Le) <<' ';
	if(C->Pa)
	  out<< 'C' << std::setfill('0') << std::setw(3) << no-C->Pa <<' ';
	else
	  out<< " nil ";
	if(C->Cf)
	  out<< 'C' << std::setfill('0') << std::setw(3) << no+C->Cf <<' ';
	else
	  out<< " nil ";
	out  << std::setfill(' ') << std::setw(1) << C->Nc << ' '
	     << 'L' << std::setfill('0') << std::setw(4) << T->NoLeaf(C->L0)
	     << ' ' << std::setfill(' ')
	     << std::setw(2) << C->Nl << ' '
	     << std::setw(4) << C->N  << ' '
	     << std::setw(8) << T->Radius(C) << ' '
	     << std::setw(8) << C->X;
      }
    } DUMP;  ///< no memory required, purely abstract data member
    /// dump leaf data
    /// \param[in] out   ostream to write to
    /// \param[in] dump  pointer to DumpTreeData or derived
    void DumpLeafs(std::ostream&out, const DumpTreeData*dump=0) const {
      if(dump==0) dump=&DUMP;
      dump->Head(BeginLeafs(),out);
      out <<'\n';
      for(const Leaf*L=BeginLeafs(); L!=EndLeafs(); ++L) {
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
      for(const Cell*C=BeginCells(); C!=EndCells(); ++C) {
	dump->Data(C,this,out);
	out <<'\n';
      }
      out.flush();
    }
    //@}
  };
#if(0)
  //
  template<int __D, typename __X, typename __F>
  struct traits< OctalTree<__D,__X,__F> > {
    static const char *name () {
      return message("OctalTree<%d,%s,%s>",
		     __D,traits<__X>::name(),traits<__F>::name());
    }
  };
  //
  template<int __D, typename __X, typename __F>
  struct traits< typename OctalTree<__D,__X,__F>::Dot> {
    static const char *name () {
      return message("OctalTree<%d,%s,%s>::Dot",
		     __D,traits<__X>::name(),traits<__F>::name());
    }
  };
  //
  template<int __D, typename __X, typename __F>
  struct traits< typename OctalTree<__D,__X,__F>::Leaf> {
    static const char *name () {
      return message("OctalTree<%d,%s,%s>::Leaf",
		     __D,traits<__X>::name(),traits<__F>::name());
    }
  };
  //
  template<int __D, typename __X, typename __F>
  struct traits< typename OctalTree<__D,__X,__F>::Cell> {
    static const char *name () {
      return message("OctalTree<%d,%s,%s>::Cell",
		     __D,traits<__X>::name(),traits<__F>::name());
    }
  };
#endif
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
# define LoopCellKids(CELL,NAME)		\
  for(Cell *NAME = (CELL)->BeginCells();	\
      NAME != (CELL)->EndCells(); ++NAME)
#endif
  /// loop cell kids of a given cell, starting somewhere
#ifndef LoopCellSecd
# define LoopCellSecd(CELL,START,NAME)		\
  for(Cell *NAME = START;			\
      NAME != (CELL)->EndCells(); ++NAME)
#endif
  /// loop leaf kids of a given cell
#ifndef LoopLeafKids
# define LoopLeafKids(CELL,NAME)		\
  for(Leaf *NAME = (CELL)->BeginLeafs();	\
      NAME != (CELL)->EndLeafKids(); ++NAME)
#endif
  /// loop leaf kids of a given cell, starting somewhere
#ifndef LoopLeafSecd
# define LoopLeafSecd(CELL,START,NAME)		\
  for(Leaf *NAME = START;			\
      NAME != (CELL)->EndLeafKids(); ++NAME)
#endif
  /// loop leaf descendants of a given cell
#ifndef LoopAllLeafs
# define LoopAllLeafs(CELL,NAME)		\
  for(Leaf *NAME = (CELL)->BeginLeafs();	\
      NAME != (CELL)->EndLeafDesc(); ++NAME)
#endif
  /// loop leaf descendants of a given cell, starting somewhere
#ifndef LoopSecLeafs
# define LoopSecLeafs(CELL,START,NAME)		\
  for(Leaf *NAME = START;			\
      NAME != (CELL)->EndLeafDesc(); ++NAME)
#endif
  /// loop all except the last leaf descendants of a given cell
#ifndef LoopLstLeafs
# define LoopLstLeafs(CELL,NAME)		\
  for(Leaf *NAME = (CELL)->BeginLeafs();	\
      NAME != (CELL)->LastLeafDesc(); ++NAME)
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
    public:
      /// perform single leaf-leaf interaction
      virtual void interact(Leaf*, Leaf*) = 0;
      /// try to perform cell-leaf interaction, return true if success
      virtual bool interact(Cell*, Leaf*) = 0;
      /// try to perform cell-cell interaction, return true if success
      virtual bool interact(Cell*, Cell*) = 0;
      /// try to perform cell self-interaction, return true if success
      virtual bool interact(Cell*) = 0;
      /// which of two cells of an interaction to split?
      virtual bool split_left(const Cell*A, const Cell*B) const {
	return A->N > B->N;
      }
      /// perform all leaf-leaf interactions between a set of leafs
      /// \note Not abstract, but virtual: default calls individual leaf-leaf
      ///       interact N*(N-1)/2 times.
      virtual void interact_many(Leaf*L0, Leaf*LN) {
	for(Leaf*Li=L0; Li<LN; ++Li)
	  for(Leaf*Lj=Li+1; Lj<LN; ++Lj)
	    interact(Li,Lj);
      }
    };
  private:
    //
    typedef std::pair<Cell*,Leaf*> pCL; ///< represents a cell-leaf interaction
    typedef std::pair<Cell*,Cell*> pCC; ///< represents a cell-cell interaction
    //
    Interactor*const IA;                ///< pter to interactor
    Stack<pCL>       CL;                ///< stack of cell-leaf interactions
    Stack<pCC>       CC;                ///< stack of cell-cell interactions
    //
    void perform(Cell*A, Leaf*B) { if(!IA->interact(A,B)) CL.push(pCL(A,B)); }
    void perform(Cell*A, Cell*B) { if(!IA->interact(A,B)) CC.push(pCC(A,B)); }
    void perform(Cell*A)         { if(!IA->interact(A))   CC.push(pCC(A,0)); }
    /// clear the stack of cell-leaf interactions
    void clear_CL_stack()
    {
      while(!CL.is_empty()) {
	pCL P = CL.pop();
	LoopLeafKids(P.first,Li)
	  IA->interact(Li,P.second);
	LoopCellKids(P.first,Ci)
	  perform(Ci,P.second);
      }
    }
    /// split a mutual cell-cell interaction
    /// \param[in] A cell to be split
    /// \param[in] B cell to be kept
    void split(Cell*A, Cell*B)
    {
      LoopLeafKids(A,Li)
	perform(B,Li);
      LoopCellKids(A,Ci)
	perform(Ci,B);
    }
    /// split a cell self-interaction
    /// \param[in] A cell to be split
    void split(Cell*A)
    {
      // leaf-leaf sub-interactions
      if(A->NumLeafKids() > 1)
	IA->interact_many(A->BeginLeafs(),A->EndLeafKids());
      // self-interactions between sub-cells
      LoopCellKids(A,Ci)
	perform(Ci);
      // mutual interaction between sub-cells and sub-leafs
      LoopCellKids(A,Ci) {
	LoopCellSecd(A,Ci+1,Cj)
	  perform(Ci,Cj);
	LoopLeafKids(A,Li)
	  perform(Ci,Li);
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
    MutualOctTreeWalker(Interactor*i, size_t d) :
      IA(i), CL(OctTree::Nsub*d), CC(2*(OctTree::Nsub-1)*(d+1)+1) {}
    /// construction for interaction between two octal trees
    /// \param[in] i  pter to interactor
    /// \param[in] d1 depth of tree of sources
    /// \param[in] d2 depth of tree of sinks
    MutualOctTreeWalker(Interactor*i, size_t d1, size_t d2) :
      IA(i), CL(OctTree::Nsub*max(d1,d2)), CC((OctTree::Nsub-1)*(d1+d2-1)+1) {}
    /// perform a mutual tree walk.
    /// If two distinct cells are given, the mutual interaction between these is
    /// performed, including splitting them until no unresolved interaction
    /// remains.\n
    /// If only one cell is given (or if both are the same), the mutual
    /// self-interaction of this cell is performed, including splitting it until
    /// no unresolved interaction remains.
    /// \param[in] A interaction cell
    /// \param[in] B (optional) second interacting cell
    void walk(Cell*A, Cell*B=0)
    {
      CC.reset();
      CL.reset();
      if(B && B!=A) perform(A,B);
      else          perform(A);
      clear_CC_stack();
    }
  };
  /// walk an oct-tree: simple wrapper around MutualTreeWalker.
  template<typename OctTree>
  void MutualOctTreeWalk(const OctTree*T,
			 typename MutualOctTreeWalker<OctTree>::Interactor*I)
  {
    MutualOctTreeWalker<OctTree> MTW(I,T->Depth());
    MTW.walk(T->Root());
  }
} // namespace WDutils
////////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_octtree_h
