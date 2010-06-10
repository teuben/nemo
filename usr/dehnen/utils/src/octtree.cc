// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file    utils/src/octtree.cc
/// 
/// \brief   implements utils/inc/octtree.h
///
/// \author  Walter Dehnen
///
/// \date    2009,2010
///
/// \note    originally based on falcON's tree.cc (by the same author)
/// 
/// \version 08-may-2009 WD  real test: debugged error in linking
/// \version 13-may-2009 WD  abolished Peano-Hilbert support
/// \version 25-sep-2009 WD  new version using indices for Leaf & Cell
/// \version 14-oct-2009 WD  new version tested against old, old abolished.
/// \version 27-jan-2010 WD  added leaf's parent cell, changed Node magic
/// \version 28-jan-2010 WD  TreeAccess::SmallestContainingCell tested
/// \version 29-jan-2010 WD  NeighbourFinder, NearestNeighbourFinder tested
/// \version 26-feb-2010 WD  new initialiser interface
/// \version 22-mar-2010 WD  FastNeighbourFinder
/// \version 15-apr-2010 WD  tree pruning, OctTree::build()
/// \version 22-apr-2010 WD  parameter nmin, changes in tree building
/// \version 24-apr-2010 WD  changes in tree building code.
/// \version 09-jun-2010 WD  16-byte alignement, using geometry.h
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
#include <octtree.h>
#include <heap.h>
#include <cstring>

#define TESTING
#undef  TESTING

#ifdef TESTING
# include <fstream>
# warning compilation for TESTING purposes only
#endif

#ifndef WD_HOT
# if defined(__GNUC__) && (__GNUC__ > 3) && (__GNUC_MINOR__ > 1)
#  define WD_HOT __attribute__((hot))
# else
#  define WD_HOT
# endif
#endif

//
// Wdutils::OctalTree<Dim,real>
//
namespace {
  using std::setw;
  using std::setfill;
  using namespace WDutils;
  ///
  template<int D> struct TreeHelper;
  template<> struct TreeHelper<2> {
    template<typename real>
    static tupel<2,real> Integer(tupel<2,real> const&x)
    {
      tupel<2,real> c;
      c[0]=int(x[0]+real(0.5));
      c[1]=int(x[1]+real(0.5));
      return c;
    }
    template<typename real>
    static real RootRadius(tupel<2,real> const&X,
			   tupel<2,real>const&Xmin,
			   tupel<2,real> const&Xmax)
    {
      real D  = max(Xmax[0]-X[0], X[0]-Xmin[0]);
      real R1 = max(Xmax[1]-X[1], X[1]-Xmin[1]);
      if(R1>D) D=R1;
      return pow(real(2), int(1+std::log(D)/M_LN2));
    }
  };
  template<> struct TreeHelper<3> {
    template<typename real>
    static tupel<3,real> Integer(tupel<3,real> const&x)
    {
      tupel<3,real> c;
      c[0]=int(x[0]+real(0.5));
      c[1]=int(x[1]+real(0.5));
      c[2]=int(x[2]+real(0.5));
      return c;
    }
    template<typename real>
    static real RootRadius(tupel<3,real> const&X,
			   tupel<3,real>const&Xmin,
			   tupel<3,real> const&Xmax)
    {
      real D  = max(Xmax[0]-X[0], X[0]-Xmin[0]);
      real R1 = max(Xmax[1]-X[1], X[1]-Xmin[1]); if(R1>D) D=R1;
      R1 = max(Xmax[2]-X[2], X[2]-Xmin[2]); if(R1>D) D=R1;
      return pow(real(2), int(1+std::log(D)/M_LN2));
    }
  };
  /// type to estimate the number of tree boxes needed.
  class EstimateNalloc
  {
    const size_t Ndots;
    const size_t Nsofar;
  public:
    EstimateNalloc(size_t a, size_t b) : Ndots(a), Nsofar(b) {}
    size_t operator() (size_t Nused) const
    {
      double x = Nused*(double(Ndots)/double(Nsofar)-1);
      return size_t(x+4*std::sqrt(x)+16);
    }
  };
  /// estimate number of boxes needed for BoxDotTree
  /// \param[in] Ndot  # dots
  /// \param[in] Nmax  max # dots/octant
  template<int Dim>
  inline size_t NBoxes(size_t Ndot, uint8 Nmax)
  {
    double
      x = double(400*Nmax+Ndot)/double(100*Nmax+Ndot);
    x  *= Ndot;
    x  /= Nmax+Nmax;
    x  += 4*std::sqrt(x);
    if(x<100) x=100;
    return size_t(Dim==2? 1.5*x : x);
  }
  /// tree of boxes and dots.
  ///
  /// An OctalTree is build by first making a BoxDotTree, then mapping it to
  /// an OctalTree via BoxDotTree::Link()
  template<int Dim, typename real>
  struct BoxDotTree
  {
    const static int Nsub = 1<<Dim; ///< number of octants per cell
    //
    typedef OctalTree<Dim,real>            OctTree;
    typedef Geometry::Algorithms<1>        GeoAlg;
    typedef typename OctTree::Initialiser  Initialiser;
    typedef typename OctTree::node_index   node_index;
    typedef typename OctTree::depth_type   depth_type;
    typedef typename OctTree::local_count  local_count;
    typedef typename OctTree::particle_key particle_key;
    typedef typename OctTree::point        point;
    typedef typename OctTree::cube         cube;
    //
    const static depth_type MAXD = OctTree::MaximumDepth;
    /// represents a particle in the BoxDotTree
    struct __Dot {
      point             X;             ///< position
      particle_key      I;             ///< identifier of associated particle
      mutable void     *Next;          ///< next dot in a linked list
    };
    typedef SSE::Extend16<__Dot> Dot;
    /// represents a cubic cell in the BoxDotTree
    struct __Box {
      typedef typename meta::__IWORDS<Nsub>::integer_u ndl_type;
      cube              CUB;           ///< box cubus
      void             *OCT[Nsub];     ///< octants
      // if 0 == NDL[i]          the octant is empty
      // if 0 <  NDL[i] <= NMAX  the octant holds a list of NDL[] leafs
      // if      NDL[i]  > NMAX  the octant holds a box
      union {
	ndl_type        NDl;
	uint8           NDL[Nsub];     ///< # dots/octant
      };
      node_index        NUM;           ///< total # dots in box
      uint8             NBX;           ///< number of daughter boxes  <= Nsub
      uint8             NOC;           ///< number of octant occupied <= NBX
      uint8             LEV;           ///< tree level of box
      uint8             PEA;           ///< Peano key, not currently used
    };
    typedef SSE::Extend16<__Box> Box;
    //
#define  pDOT(d) static_cast<Dot*>((d))
#define cpDOT(d) static_cast<const Dot*>((d))
#define  pBOX(d) static_cast<Box*>((d))
#define cpBOX(d) static_cast<const Box*>((d))
    /// \name data
    //@{
    node_index          NDOT;          ///< number of dots to load
    Dot          *const D0;            ///< begin of dots
    const uint8         NMAX;          ///< maximum # particles/octant
    const uint8         NMIN;          ///< minimum # particles/cell
    const bool          AVSPC;         ///< avoid single-parent cells
    const node_index    NMAX1;         ///< NMAX + 1
    depth_type          DEPTH;         ///< depth of linked tree
    block_alloc<Box>    BM;            ///< allocator for boxes
    Box                *P0;            ///< root box
    /// after building: # cells required; after linking: # cells actually used
    node_index          NCELL;
    mutable node_index  CF;            ///< free cells during linking
    mutable node_index  LF;            ///< free leafs during linking
    mutable size_t      ND;            ///< # dots added sofar
    const OctTree      *TREE;          ///< tree to be linked
    //@}
    /// replace (pter to) box by (pter to) its non-single-parent descendant
    /// \param[in,out] P  (pter to) given box
    void EnsureNonSingleParent(const Box*&P) const
    {
      while(P->NOC==1 && P->NBX==1)
	for(int i=0; i!=Nsub; ++i)
	  if(P->NDL[i]) {
	    P = pBOX(P->OCT[i]);
	    break;
	  }
    }
    /// \name methods used in construction
    //@{
    /// provides a new empty (daughter) box in the i th octant of B
    /// \return new box
    /// \param[in] B  parent box
    /// \param[in] i  parent box's octant
    Box*MakeSubBox(const Box*B, int i) WDutils_THROWING
    {
      Box *P = BM.new_element(EstimateNalloc(NDOT,ND));
      P->NDl = 0;
      P->NBX = 0;
      P->NOC = 0;
      P->LEV = B->LEV + 1;
      if(P->LEV >= MAXD)
	WDutils_THROW("exceeding maximum tree depth of %d\n         "
		      "(perhaps more than Nmax=%du positions are identical "
		      "within floating-point precision)\n", MAXD, NMAX);
      GeoAlg::copy(B->CUB,P->CUB);
      GeoAlg::shrink(P->CUB,i);
      return P;
    }
    /// adds dot @a Di to octant @a b of box @a P
    void AddDotToOctant(Box*P, Dot*Di, int b)
    {
      if(P->NDL[b] == 0) {
	P->NOC ++;
	Di->Next = 0;
      } else
	Di->Next = P->OCT[b];
      P->OCT[b]  = Di;
      P->NDL[b] ++;
      if(P->NDL[b] == NMIN) ++NCELL;
    }
    /// adds a dot to a box.
    ///
    /// If the dot's octant @a b contains another daughter box as indicated by
    /// Box::NDL[b] > NMAX, the process is repeated on that box.\n
    ///
    /// Otherwise, the dot is added to the linked list of dots in octant @a b.
    /// If this makes the number of dots to exceed @a NMAX, a new daughter box
    /// is created to hold the dots in the linked list. This latter process is
    /// possibly repeated (if all dots are in just one octant of the daughter
    /// box).
    ///
    /// \param[in] P    box to add to
    /// \param[in] Di   dot to add
    void AddDot(Box*P, Dot*Di) WDutils_THROWING
    {
      for(;;) {
	// loop boxes down the tree
	P->NUM++;
	int b = GeoAlg::octant(P->CUB,Di->X);
	if(P->NDL[b] > NMAX)
	  // octant is a box
	  P = pBOX(P->OCT[b]);
	else {
	  // octant either empty or has a list of dots
	  AddDotToOctant(P,Di,b);              // add dot to octant in box
	  while(P->NDL[b] > NMAX) {            // WHILE dot list too long
	    P->NBX++;                          //   increment box counter
	    Di        = pDOT(P->OCT[b]);       //   get first dot from octant
	    P->OCT[b] = MakeSubBox(P,b);       //   new empty box in octant
	    P         = pBOX(P->OCT[b]);       //   P = daughter box
	    P->NUM    = NMAX1;                 //   set P->NUM
	    for(Dot*Dn; Di; Di=Dn) {           //   sort dots into octants
	      Dn= pDOT(Di->Next);
	      b = GeoAlg::octant(P->CUB,Di->X);
	      AddDotToOctant(P,Di,b);
	    }
	  }
	  return;
	}
      }
    }
    /// ctor: build a BoxDotTree
    /// If old tree given, we add the dots in the order of its leafs
    /// \param[in] build type of tree build to perform
    /// \param[in] Ndot  number of positions
    /// \param[in] Init  initializer to re-initiliase particle data
    /// \param[in] nmax  max # dots / octant in building
    /// \param[in] nmin  min leaf / cell in linking
    /// \param[in] tree  old OctalTree, needed if build='r' or 'p'
    BoxDotTree(char build, node_index Ndot, const Initialiser*Init,
	       depth_type nmax, depth_type nmin, bool avspc,
	       const OctTree*tree=0) WDutils_THROWING WD_HOT;
    /// dtor
    ~BoxDotTree();
    //@}
    /// \name methods used in linking to OctalTree
    //@{
    /// copy Dot data to Leaf
    /// \param[in] L  leaf index in tree
    /// \param[in] D  dot to link leaf to
    /// \param[in] P  parent cell index for leaf
    void LinkLeaf(node_index L, const Dot*D, node_index P) const
    {
      static_cast<point&>(TREE->XL[L]) = D->X;
      TREE->PL[L] = D->I;
      TREE->PC[L] = P;
    }
    /// tree linking: make a cell from an octant of a parent box
    /// \param[in] P  parent box
    /// \param[in] C  index of current cell to be linked
    /// \param[in] o  octant of cell in @a P
    void LinkOctantCell(const Box*P, node_index C, int i) const
    {
      TREE->L0[C] = LF;
      TREE->OC[C] = i;
      TREE->LE[C] = P->LEV + 1;
      GeoAlg::copy(P->CUB,static_cast<cube&>(TREE->XC[C]));
      GeoAlg::shrink(TREE->XC[C],i);
      TREE->NM[C] = P->NDL[i];
      TREE->NL[C] = P->NDL[i];
      TREE->NC[C] = 0;
      TREE->CF[C] = C;
      for(const Dot*Di=cpDOT(P->OCT[i]); Di; Di=pDOT(Di->Next))
	LinkLeaf(LF++,Di,C);
    }
    /// tree linking: leaf & cell descendants are continuous in memory
    /// \param[in] P  box to link with C
    /// \param[in] C  index of current cell to be linked
    /// \param[in] o  octant of box @a P in parent
    /// \return       tree depth of cell C
    /// \note recursive.
    /// \node uses data CF and LF
    depth_type LinkCell(const Box*P, node_index C, int o) const WD_HOT;
    /// frontend for tree linking
    /// \param[in] tree  tree to be linked
    /// \param[in] nmin  only make cells from boxes with at least @a nmin dots
    void Link(const OctTree*tree) const
    {
      const_cast<const OctTree*&>(TREE) = tree;
      CF = 1;
      LF = 0;
      TREE->PA[0] = 0;
      const_cast<depth_type&>(DEPTH) = LinkCell(P0,0,0);
      const_cast<node_index&>(NCELL) = CF;
    }
    //
    node_index const&NCell() const { return NCELL; }
    //@}
#ifdef TESTING
    /// header for dot data dump
    void DumpHeadDot(std::ostream&out)
    {
      out<<" Dot    Next       I                     X\n";
    }
    /// dump dot data
    void Dump(const Dot*D, std::ostream&out)
    {
      out<<" D"<<setfill('0')<<setw(5)<<int(D-D0);
      if(D->Next)
	out<<" D"<<setfill('0')<<setw(5)<<int(pDOT(D->Next)-D0);
      else
	out<<" nil   ";
      out<<' '<<setfill(' ')<<setw(5)<<D->I<<' '<<setw(10)<<D->X<<'\n';
    }
    /// header for box data dump
    void DumpHeadBox(std::ostream&out)
    {
      out<<" Box        N No Nb"
	 <<(Dim==2? "    NDL[]   " : "          NDL[]         ");
      for(int i=0; i!=Nsub; ++i) out<<" OCT["<<i<<']';
      out<<" S NSP         Rad                 C\n";
    }
    /// Find first non-single parent descendant
    const Box* NonSingleParent(const Box*B) const
    {
      const Box*P = B;
      EnsureNonSingleParent(P);
      return P;
    }
    /// dump box data
    /// \param[out] ns  counter for single-parent boxes
    void Dump(const Box*B, std::ostream&out, node_index&ns)
    {
      out<<" B"<<setfill('0')<<setw(5)<<BM.number_of_element(B)
	 <<' ' <<setfill(' ')<<setw(5)<<B->NUM
	 <<' ' <<setw(2)<<int(B->NOC)
	 <<' ' <<setw(2)<<int(B->NBX);
      for(int i=0; i!=Nsub; ++i)
	out<<' '<<setw(2)<<int(B->NDL[i]);
      for(int i=0; i!=Nsub; ++i) {
	if(B->NDL[i]==0)
	  out<<" nil   ";
	else if(B->NDL[i] <= NMAX)
	  out<<" D"<<setfill('0')<<setw(5)
	     <<int(pDOT(B->OCT[i])-D0);
	else
	  out<<" B"<<setfill('0')<<setw(5)
	     <<BM.number_of_element(pBOX(B->OCT[i]));
      }
      if(B->NOC==1 && B->NBX==1) {
	++ns;
	out<<" S";
      } else
	out<<" M";
      out<<" B"<<setfill('0')<<setw(5)
	 <<BM.number_of_element(NonSingleParent(const_cast<Box*>(B)));
      out<<' '<<setfill(' ')<<setw(8)<<B->CUB.H
	 <<' '<<setw(8)<<B->CUB.X<<'\n';
    }
    /// dump tree, recursive
    void Dump(const Box*B, std::ostream&outd, std::ostream&outb, node_index&ns)
    {
      Dump(B,outb,ns);
      for(int i=0; i!=Nsub; ++i)
	if     (B->NDL[i] >  NMAX) Dump(pBOX(B->OCT[i]),outd,outb,ns);
	else if(B->NDL[i] != 0   ) Dump(pDOT(B->OCT[i]),outd);
    }
    /// dump tree
    /// \param[in] outd ostream for dumping Dot data
    /// \param[in] outb ostream for dumping Box data
    void Dump(std::ostream&outd, std::ostream&outb) {
      node_index ns=0;
      DumpHeadDot(outd);
      DumpHeadBox(outb);
      Dump(P0,outd,outb,ns);
      std::cerr<<" # single parent boxes: "<<ns<<'\n';
    }
#endif // TESTING
  };
}
//
namespace WDutils {
  //
  template<> struct traits< ::BoxDotTree<2,float>::Dot >
  { static const char*name() { return "BoxDotTree<2,float>::Dot"; } };
  template<> struct traits< ::BoxDotTree<3,float>::Dot >
  { static const char*name() { return "BoxDotTree<3,float>::Dot"; } };
  template<> struct traits< ::BoxDotTree<2,double>::Dot >
  { static const char*name() { return "BoxDotTree<2,double>::Dot"; } };
  template<> struct traits< ::BoxDotTree<3,double>::Dot >
  { static const char*name() { return "BoxDotTree<3,double>::Dot"; } };
  //
  template<> struct traits< ::BoxDotTree<2,float>::Box >
  { static const char*name() { return "BoxDotTree<2,float>::Box"; } };
  template<> struct traits< ::BoxDotTree<3,float>::Box >
  { static const char*name() { return "BoxDotTree<3,float>::Box"; } };
  template<> struct traits< ::BoxDotTree<2,double>::Box >
  { static const char*name() { return "BoxDotTree<2,double>::Box"; } };
  template<> struct traits< ::BoxDotTree<3,double>::Box >
  { static const char*name() { return "BoxDotTree<3,double>::Box"; } };
  //
  template<> struct traits<block_alloc< ::BoxDotTree<2,float>::Box > >
  { static const char*name()
    { return "block_alloc<BoxDotTree<2,float>::Box>"; }
  };
  template<> struct traits<block_alloc< ::BoxDotTree<3,float>::Box > >
  { static const char*name()
    { return "block_alloc<BoxDotTree<3,float>::Box>"; }
  };
  template<> struct traits<block_alloc< ::BoxDotTree<2,double>::Box > >
  { static const char*name()
    { return "block_alloc<BoxDotTree<2,double>::Box>"; }
  };
  template<> struct traits<block_alloc< ::BoxDotTree<3,double>::Box > >
  { static const char*name()
    { return "block_alloc<BoxDotTree<3,double>::Box>"; }
  };
  //
  template<> struct traits<block_alloc< ::BoxDotTree<2,float>::Box >::block >
  { static const char*name()
    { return "block_alloc<BoxDotTree<2,float>::Box>::block"; }
  };
  template<> struct traits<block_alloc< ::BoxDotTree<3,float>::Box >::block >
  { static const char*name()
    { return "block_alloc<BoxDotTree<3,float>::Box>::block"; }
  };
  template<> struct traits<block_alloc< ::BoxDotTree<2,double>::Box >::block >
  { static const char*name()
    { return "block_alloc<BoxDotTree<2,double>::Box>::block"; }
  };
  template<> struct traits<block_alloc< ::BoxDotTree<3,double>::Box >::block >
  { static const char*name()
    { return "block_alloc<BoxDotTree<3,double>::Box>::block"; }
  };
}
//
namespace {
  //
  template<int Dim, typename real>
  BoxDotTree<Dim,real>::BoxDotTree(char build, node_index Ndot,
				   const Initialiser*Init,
				   depth_type nmax, depth_type nmin,
				   bool avspc, const OctTree*Tree)
    WDutils_THROWING : 
    NDOT  ( Ndot ),
    D0    ( NDOT? new16<Dot>(NDOT) : 0 ),
    NMAX  ( nmax ),
    NMIN  ( nmin ),
    AVSPC ( avspc ),
    NMAX1 ( nmax + 1 ),
    BM    ( NBoxes<Dim>(Ndot,nmax) ),
    P0    ( BM.new_element() ),
    NCELL ( 1 )
  {
    Dot *const DN=D0+NDOT;
    // 1 initialise dots
    switch(build) {
    case 'n':
      // 1.1    build from scratch: initialise all dots
      for(Dot*Di=D0; Di!=DN; ++Di)
	Init->Initialise(Di->I,Di->X);
      break;
    case 'r': {
      // 1.2    re-building using the old tree order
      Dot*Di=D0, *List=0, *DN1=D0+min(NDOT,Tree->Nleafs());
      // 1.2.1  try to re-initialise dot positions from particle key in old tree
      //        put uninitialised dots in linked list
      for(node_index i=0; Di!=DN1; ++Di,++i) {
	Di->I = Tree->PL[i];
	if(! Init->ReInitialiseValid(Di->I,Di->X) ) { Di->Next=List; List=Di; }
      }
      // 1.2.2  initialise particle keys and positions for dots in linked list
      for(Di=List; Di; Di=pDOT(Di->Next))
	Init->ReInitialiseInvalid(Di->I,Di->X);
      // 1.2.3  initialise particle keys and positions for any remaining dots
      for(Di=DN1; Di!=DN; ++Di)
	Init->ReInitialiseInvalid(Di->I,Di->X);
    } break;
    case 'p': {
      // 1.3     pruning a given (parent) tree
      if(0==NDOT) {
	// 1.3A  unknown number of particles in pruned tree
	// 1.3A.1  count number of particles in pruned tree
	for(node_index i=0; i!=Tree->Nleafs(); ++i)
	  if(Init->Pick(Tree->PL[i])) ++NDOT;
	if(0==NDOT)
	  WDutils_THROW("OctalTree<%d,%s>::build(): empty tree\n",
			Dim,nameof(real));
	// 1.3A.2  allocate dots
	const_cast<Dot*&>(D0) = WDutils_NEW(Dot,NDOT);
	const_cast<Dot*&>(DN) = D0+NDOT;
	// 1.3A.3  initialise dots
	Dot*Di=D0;
	for(node_index i=0; i!=Tree->Nleafs(); ++i)
	  if(Init->Pick(Tree->PL[i])) {
	    Di->I = Tree->PL[i];
	    Di->X = static_cast<point const&>(Tree->XL[i]);
	    ++Di;
	  }
      } else {
	// 1.3B  assume no more than Ndot leafs of parent tree are picked
	Dot*Di=D0;
	for(node_index i=0; i!=Tree->Nleafs(); ++i)
	  if(Init->Pick(Tree->PL[i])) {
	    if(Di==DN)
	      WDutils_THROW("OctalTree<%d,%s>::build(): "
			    "more leafs in pruned tree than expected (%d)\n",
			    Dim,nameof(real),NDOT);
	    Di->I = Tree->PL[i];
	    Di->X = static_cast<point const&>(Tree->XL[i]);
	    ++Di;
	  }
	if(Di < DN) NDOT = Di-D0;
      }
    } break;
    default: WDutils_THROW("OctalTree<%d,%s>::build(): unknown build '%c'\n",
			   Dim,nameof(real),build);
    }
    // 2  find min, max and average position
    Dot*Di=D0;
    point Xmin(D0->X), Xmax(D0->X), Xave(D0->X);
    for(++Di; Di!=DN; ++Di) {
      Di->X.up_min_max(Xmin,Xmax);
      Xave += Di->X;
    }
    if(isnan(Xave) || isinf(Xave)) {
      for(Di=D0; Di!=DN; ++Di)
	if(isnan(Di->X) || isinf(Di->X))
	  Dim==2?
	    WDutils_THROW("OctalTree<%d,%s>: particle %d: "
			  "invalid X=%g,%g\n",Dim,nameof(real),
			  Di->I, Di->X[0],Di->X[1])
	    :
	    WDutils_THROW("OctalTree<%d,%s>: particle %d: "
			  "invalid X=%g,%g,%g\n",Dim,nameof(real),
			  Di->I, Di->X[0],Di->X[1],Di->X[2]);
      WDutils_THROW("OctalTree<%d,%s>: unidentified invalid position(s)\n",
		    Dim,nameof(real));
    }
    Xave   /= real(NDOT);
    // 3  set (empty) root box, RA[]
    P0->NDl   = 0;
    P0->NBX   = 0;
    P0->NOC   = 0;
    P0->CUB.X = TreeHelper<Dim>::Integer(Xave);
    P0->CUB.H = TreeHelper<Dim>::RootRadius(P0->CUB.X,Xmin,Xmax);
    P0->LEV   = 0;
    P0->NUM   = 0;
    // 4  add dots
    ND = 0;
    for(Di=D0; Di!=DN; ++Di,++ND)
      AddDot(P0,Di);
#ifdef TESTING
    std::ofstream dumpD("dots.dat"), dumpB("boxs.dat");
    Dump(dumpD, dumpB);
#endif
  }
  //
  template<int D, typename real>
  inline BoxDotTree<D,real>::~BoxDotTree()
  { 
    if(D0) delete16(D0);
  }
  //
  template<int D, typename real>
  typename BoxDotTree<D,real>::depth_type
  BoxDotTree<D,real>::LinkCell(const Box*P, node_index C, int o) const
  {
    // 1 if single-parent replace by non-single-parent descendant
    if(AVSPC) EnsureNonSingleParent(P);
    // 2 copy some data, set octant
    TREE->L0[C] = LF;
    TREE->OC[C] = o;
    TREE->LE[C] = P->LEV;
    GeoAlg::copy(P->CUB,static_cast<cube&>(TREE->XC[C]));
    TREE->NM[C] = P->NUM;
    // 3 loop octants: link leaf kids (from octants with < NMIN), count cells
    depth_type tmp = 0;
    for(int i=0; i!=Nsub; ++i)
      if(P->NDL[i] >= NMIN)
	++tmp;
      else if(P->NDL[i])
	for(const Dot*Di=cpDOT(P->OCT[i]); Di; Di=pDOT(Di->Next))
	  LinkLeaf(LF++,Di,C);
    // 4 set number of leaf and cell kids, first daughter cell, if any
    TREE->NL[C] = LF - TREE->L0[C];
    TREE->NC[C] = tmp;
    TREE->CF[C] = tmp? CF : C;
    if(tmp==0) return 1;
    // 5 link daughter cells
    node_index Ci = CF;
    CF += tmp;
    tmp = 1;
    for(int i=0; i!=Nsub; ++i) {
      if(P->NDL[i] > NMAX) {
	TREE->PA[Ci] = C;
	update_max(tmp, LinkCell(cpBOX(P->OCT[i]), Ci++, i));
      } else if(P->NDL[i] >= NMIN) {
	TREE->PA[Ci] = C;
	LinkOctantCell(P, Ci++, i);
      }
    }
    return ++tmp;
  }
  //
#undef   pDOT
#undef  cpDOT
#undef   pBOX
#undef  cpBOX
  /// the next multiple of 16 to n*sizeof(T)
  template<typename T>
  inline size_t Next16(size_t n)
  {
    return WDutils::next_aligned16(n*sizeof(T));
  }
}
//
namespace WDutils {
  template<int D, typename real> 
  void OctalTree<D,real>::allocate()
  {
    // all arrays are to be 16-byte aligned
    unsigned need =
      Next16<point16>       (NLEAF) +  // XL
      Next16<particle_key>  (NLEAF) +  // PL
      Next16<node_index>    (NLEAF) +  // PC
      Next16<depth_type>    (NCELL) +  // LE
      2*Next16<octant_type> (NCELL) +  // OC,NC
      Next16<cube16>        (NCELL) +  // XC
      Next16<local_count>   (NCELL) +  // NL
      4*Next16<node_index>  (NCELL);   // L0,NM,CF,PA
    if((need > NALLOC) || (3*need < 2*NALLOC)) {
      if(ALLOC) delete16(ALLOC);
      ALLOC  = new16<char>(need);
      NALLOC = need;
    }
    char* A = ALLOC;

#define SET_POINTER(Name,Number,Type)				\
    Name=reinterpret_cast<Type*>(A);  A+=Next16<Type>(Number);
    
    SET_POINTER(XL,NLEAF,point16);
    SET_POINTER(PL,NLEAF,particle_key);
    SET_POINTER(PC,NLEAF,node_index);
    SET_POINTER(LE,NCELL,depth_type);
    SET_POINTER(OC,NCELL,octant_type);
    SET_POINTER(XC,NCELL,cube16);
    SET_POINTER(L0,NCELL,node_index);
    SET_POINTER(NL,NCELL,local_count);
    SET_POINTER(NM,NCELL,node_index);
    SET_POINTER(CF,NCELL,node_index);
    SET_POINTER(NC,NCELL,octant_type);
    SET_POINTER(PA,NCELL,node_index);
#undef SET_POINTER
  }
  //
  template<int D, typename real>
  void OctalTree<D,real>::build(char building, node_index n,
				const Initialiser*init, const OctalTree*tree)
    WDutils_THROWING
  {
#undef GIVE_TIMING
#ifdef GIVE_TIMING
    clock_t cpu0, cpu1;
    cpu0 = clock();
#endif
    BoxDotTree<D,real> BDT(building,n,init,NMAX,NMIN,AVSPC,tree);
#ifdef GIVE_TIMING
    cpu1 = clock();
    std::cerr<<" OctalTree::build(): BoxDotTree::BoxDotTree took "
	     <<double(cpu1 - cpu0)/double(CLOCKS_PER_SEC)<<" sec\n";
#endif
    NLEAF = BDT.NDOT;
    NCELL = BDT.NCell();
    allocate();
#ifdef GIVE_TIMING
    cpu0 = clock();
#endif
    BDT.Link(this);
#ifdef GIVE_TIMING
    cpu1 = clock();
    std::cerr<<" OctalTree::build(): BoxDotTree::Link()     took "
	     <<double(cpu1 - cpu0)/double(CLOCKS_PER_SEC)<<" sec\n";
#endif
    DEPTH = BDT.DEPTH;
    NCELL = BDT.NCELL;
//     std::memcpy(RAD,BDT.RA,(MaximumDepth+1)*sizeof(real));
  }
}
//
// Wdutils::TreeAccess<OctTree>
//
namespace WDutils {
#ifdef __SSE__
#define  PF(__X) static_cast<float*>(__X)
#define cPF(__X) static_cast<const float*>(__X)
  template<> 
  TreeAccess<OctalTree<2,float> >::Cell
  TreeAccess<OctalTree<2,float> >::SmallestContainingCell(point const&x) const
  {
    // start with root cell
    Cell c = Root();
    __m128 X = _mm_loadu_ps(cPF(x));
    __m128 C = _mm_load_ps(cPF(box(c).X));
    __m128 H = _mm_shuffle_ps(C,C,_MM_SHUFFLE(2,2,2,2));
    // if x not in root cell, return invalid cell
    if(3 != (3&_mm_movemask_ps(_mm_and_ps(_mm_cmple_ps(_mm_sub_ps(C,H),X),
					  _mm_cmpgt_ps(_mm_add_ps(C,H),X)))))
      return InvalidCell();
    // descend down the tree
    for(;;) {
      // if cell has no sub-cells: we are done
      if(Ncells(c)==0)
	return c;
      // find sub-cell in same octant as x
      uint8 o = 3&_mm_movemask_ps(_mm_cmplt_ps(C,X));
      Cell cc = BeginCells(c), ce=EndCells(c);
      while(cc!=ce && o!=octant(cc)) ++cc;
      // non found: return c
      if(cc == ce)
	return c;
      C = _mm_load_ps(cPF(box(cc).X));
      // ensure sub-cell does contain x (otherwise return c)
      // [this may not be so if the sub-cell does not cover the full of the
      // octant, which can occur when avoiding single-parent cells]
      if(AvoidedSingleParentCells() && level(cc) > level(c)+1 &&
	 3 != (3&_mm_movemask_ps(_mm_and_ps(_mm_cmple_ps(_mm_sub_ps(C,H),X),
					    _mm_cmpgt_ps(_mm_add_ps(C,H),X)))))
	return c;
      // continue with sub-cell
      c = cc;
    }
  }
  template<> 
  TreeAccess<OctalTree<3,float> >::Cell
  TreeAccess<OctalTree<3,float> >::SmallestContainingCell(point const&x) const
  {
    Cell c = Root();
    __m128 X = _mm_loadu_ps(cPF(x));
    __m128 C = _mm_load_ps(cPF(box(c).X));
    __m128 H = _mm_shuffle_ps(C,C,_MM_SHUFFLE(3,3,3,3));
    if(7 != (7&_mm_movemask_ps(_mm_and_ps(_mm_cmple_ps(_mm_sub_ps(C,H),X),
					  _mm_cmpgt_ps(_mm_add_ps(C,H),X)))))
      return InvalidCell();
    for(;;) {
      if(Ncells(c)==0) return c;
      uint8 o = 7&_mm_movemask_ps(_mm_cmplt_ps(C,X));
      Cell cc = BeginCells(c), ce=EndCells(c);
      while(cc!=ce && o!=octant(cc)) ++cc;
      if(cc==ce)
	return c;
      C = _mm_load_ps(cPF(box(cc).X));
      if(AvoidedSingleParentCells() && level(cc) > level(c)+1 &&
	 7 != (7&_mm_movemask_ps(_mm_and_ps(_mm_cmple_ps(_mm_sub_ps(C,H),X),
					    _mm_cmpgt_ps(_mm_add_ps(C,H),X)))))
	return c;
      c = cc;
    }
  }
#undef  PF
#undef cPF
#ifdef __SSE2__
#define  PD(__X) static_cast<double*>(__X)
#define cPD(__X) static_cast<const double*>(__X)
#define  PD2(__X) static_cast<double*>(__X)+2
#define cPD2(__X) static_cast<const double*>(__X)+2
  template<> 
  TreeAccess<OctalTree<2,double> >::Cell
  TreeAccess<OctalTree<2,double> >::SmallestContainingCell(point const&x) const
  {
    Cell c = Root();
    __m128d X = _mm_loadu_pd(cPD(x));
    __m128d C = _mm_load_pd(cPD(box(c).X));
    __m128d H = _mm_set1_pd(box(c).H);
    if(3 != _mm_movemask_pd(_mm_and_pd(_mm_cmple_pd(_mm_sub_pd(C,H),X),
				       _mm_cmpgt_pd(_mm_add_pd(C,H),X))))
      return InvalidCell();
    for(;;) {
      if(Ncells(c)==0) return c;
      uint8 o = _mm_movemask_pd(_mm_cmplt_pd(C,X));
      Cell cc = BeginCells(c), ce=EndCells(c);
      while(cc!=ce && o!=octant(cc)) ++cc;
      if(cc==ce) return c;
      C = _mm_load_pd(cPD(box(cc).X));
      if(AvoidedSingleParentCells() && level(cc) > level(c)+1 &&
	 3 != _mm_movemask_pd(_mm_and_pd(_mm_cmple_pd(_mm_sub_pd(C,H),X),
					 _mm_cmpgt_pd(_mm_add_pd(C,H),X))))
	return c;
      c = cc;
    }
  }
  template<> 
  TreeAccess<OctalTree<3,double> >::Cell
  TreeAccess<OctalTree<3,double> >::SmallestContainingCell(point const&x) const
  {
    Cell c = Root();
    __m128d X0 = _mm_loadu_pd(cPD(x));
    __m128d C0 = _mm_load_pd(cPD(box(c).X));
    __m128d H  = _mm_set1_pd(box(c).H);
    if(3 != _mm_movemask_pd(_mm_and_pd(_mm_cmple_pd(_mm_sub_pd(C0,H),X0),
				       _mm_cmpgt_pd(_mm_add_pd(C0,H),X0))))
      return InvalidCell();
    __m128d C1 = _mm_load_pd(cPD2(box(c).X));
    __m128d X1 = _mm_loadu_pd(cPD2(x));
    if(! (1&_mm_movemask_pd(_mm_and_pd(_mm_cmple_pd(_mm_sub_pd(C1,H),X1),
				       _mm_cmpgt_pd(_mm_add_pd(C1,H),X1)))))
      return InvalidCell();
    for(;;) {
      if(Ncells(c)==0) return c;
      uint8 o = _mm_movemask_pd(_mm_cmplt_pd(C0,X0)) |
	(    (1&_mm_movemask_pd(_mm_cmplt_pd(C1,X1)))<<2) ;
      Cell cc = BeginCells(c), ce=EndCells(c);
      while(cc!=ce && o!=octant(cc)) ++cc;
      if(cc==ce) return c;
      C0 = _mm_load_pd(cPD (box(cc).X));
      C1 = _mm_load_pd(cPD2(box(cc).X));
      if(AvoidedSingleParentCells() && level(cc) > level(c)+1 &&
	 ( 3 != _mm_movemask_pd(_mm_and_pd(_mm_cmple_pd(_mm_sub_pd(C0,H),X0),
					   _mm_cmpgt_pd(_mm_add_pd(C0,H),X0)))
	   ||
	   ! (1&_mm_movemask_pd(_mm_and_pd(_mm_cmple_pd(_mm_sub_pd(C1,H),X1),
					   _mm_cmpgt_pd(_mm_add_pd(C1,H),X1))))
	   )
	 )
	return c;
      c  = cc;
    }
  }
#undef  PD
#undef cPD
#undef  PD2
#undef cPD2
#else  // no __SSE2__ but __SSE__
  template<int D>
  typename TreeAccess<OctalTree<Dim,double> >::Cell
  TreeAccess<OctalTree<Dim,double> >::SmallestContainingCell(point const&x)
    const
  {
    Cell c=Root();
    if(! Geometry::Algorithms<0>::contains(box(c),x))
      return InvalidCell();
    for(;;) {
      if(Ncells(c)==0)
	return c;
      uint8 o=::octant(centre(c),x);
      Cell cc=BeginCells(c), ce=EndCells(c);
      while(cc!=ce && o!=octant(cc)) ++cc;
      if(cc==ce || 
	 AvoidedSingleParentCells() &&
	 level(cc) > level(c)+1     &&
	 ! Geometry::Algorithms<0>::contains(box(cc),x))
	return c;
      c=cc;
    }
  }  
#endif // __SSE2__
#else  // no __SSE__
  template<typename OctTree>
  typename TreeAccess<OctTree>::Cell
  TreeAccess<OctTree>::SmallestContainingCell(point const&x) const
  {
    Cell c=Root();
    if(! Geometry::Algorithms<0>::contains(box(c),x))
      return InvalidCell();
    for(;;) {
      if(Ncells(c)==0)
	return c;
      uint8 o=::octant(centre(c),x);
      Cell cc=BeginCells(c), ce=EndCells(c);
      while(cc!=ce && o!=octant(cc)) ++cc;
      if(cc==ce || 
	 AvoidedSingleParentCells() &&
	 level(cc) > level(c)+1     &&
	 ! Geometry::Algorithms<0>::contains(box(cc),x))
	return c;
      c=cc;
    }
  }
#endif // __SSE__
}
//
// Wdutils::NeighbourLoop<OctTree>
//
namespace WDutils {
  template<typename OctTree>
  void NeighbourLoop<OctTree>::ProcessCell(Cell Ci, node_index cC) const
  {
    if(cC==0 && Number(Ci) <= NDIR)
      ProcessLeafs(BeginLeafs(Ci),EndLeafDesc(Ci));
    else {
      if(Nleafkids(Ci))
	ProcessLeafs(BeginLeafs(Ci),EndLeafKids(Ci));
      if(Ncells(Ci)>cC) {
	if(cC) { LoopCellKids(Ci,c) if(c!=C && !Outside(c)) ProcessCell(c);}
	else   { LoopCellKids(Ci,c) if(!Outside(c)) ProcessCell(c); }
      }
    }
  }
  //
  template<typename OctTree> inline
  void NeighbourLoop<OctTree>::Process()
  {
    if(IsValid(C))
      for(Cell P=C; IsValid(P) && (C==P || !Inside(C)); C=P,P=Parent(C))
	ProcessCell(P, C!=P);
    else
      ProcessCell(Base::Root(), 0);
  }
}
//
// Wdutils::NeighbourFinder<OctTree>
//
namespace {
  /// Processor: add neighbours to neighbour list
  template<typename OctTree>
  struct Lister : public NeighbourFinder<OctTree>::Processor {
    typedef typename NeighbourFinder<OctTree>::Leaf Leaf;
    typedef typename NeighbourFinder<OctTree>::real real;
    typedef typename NeighbourFinder<OctTree>::node_index node_index;
    //
    Neighbour<OctTree> *LIST;    ///< neighbour list
    const   node_index  K;       ///< size of list
    mutable node_index  I;       ///< index of current element
    /// ctor: take data
    Lister(Neighbour<OctTree>*list, node_index size)
      : LIST(list), K(size), I(0) {}
    /// process: add neighbour to list
    void process(Leaf l, real q) const
    {
      if(I<K) {
	LIST[I].Q = q;
	LIST[I].L = l;
      }
      ++I;
    }
  };
}
//
namespace WDutils {
  template<typename OctTree> inline
  void NeighbourFinder<OctTree>::ProcessLeafs(Leaf b, Leaf e) const
  {
    for(Leaf l=b; l!=e; ++l) {
      real q = dist_sq(S.X,position(l));
      if(q<S.Q) PROC->process(l,q);
    }
  }
  //
  template<typename OctTree>
  typename NeighbourFinder<OctTree>::node_index
  NeighbourFinder<OctTree>::Find(Leaf l, real q, Neighbour<OctTree>*nb,
				 node_index m)
  {
    Lister<OctTree> LL(nb,m);
    PROC =&LL;
    S.Q  = q;
    S.X  = position(l);
    C    = Parent(l);
    Base::Process();
    return LL.I;
  }
  //
  template<typename OctTree>
  typename NeighbourFinder<OctTree>::node_index
  NeighbourFinder<OctTree>::Find(point const&x, real q, Neighbour<OctTree>*nb,
				 node_index m)
  {
    Lister<OctTree> LL(nb,m);
    PROC =&LL;
    S.Q  = q;
    S.X  = x;
    C    = SmallestContainingCell(S.X);
    Base::Process();
    return LL.I;
  }
  //
  template<typename OctTree>
  void NeighbourFinder<OctTree>::Process(Leaf l, real q, const Processor*p)
    WDutils_THROWING
  {
    if(0==p) WDutils_THROW("NeighbourFinder::Process(): p=0\n");
    PROC = p;
    S.Q  = q;
    S.X  = position(l);
    C    = Parent(l);
    Base::Process();
  }
  //
  template<typename OctTree>
  void NeighbourFinder<OctTree>::Process(point const&x, real q,
					 const Processor*p) WDutils_THROWING
  {
    if(0==p) WDutils_THROW("NeighbourFinder::Process(): p=0\n");
    PROC = p;
    S.Q  = q;
    S.X  = x;
    C    = SmallestContainingCell(S.X);
    Base::Process();
  }
}
#ifdef __SSE__
//
// Wdutils::PositionsSSE<OctTree>
//
namespace WDutils {
  template<typename OctTree>
  void PositionsSSE<OctTree>::Allocate(typename Access::node_index nl)
  {
    const_cast<unsigned&>(N16) = (nl+L) & nL;
    size_t nx  = N16 * sizeof(real);
    size_t na  = N16 + OctTree::Dim * nx;
    if(na > NALLOC || 2*na < 3*NALLOC) {
      if(ALLOC) free16(ALLOC);
      const_cast<size_t&>(NALLOC) = na;
      const_cast<char* &>(ALLOC)  = NALLOC? new16<char>(na) : 0;
    }
    char* A = ALLOC;
    const_cast<real*&>(XX)=reinterpret_cast<real*>(A);  A+=nx;
    const_cast<real*&>(YY)=reinterpret_cast<real*>(A);  A+=nx;
    if(OctTree::Dim == 3) {
      const_cast<real*&>(ZZ)=reinterpret_cast<real*>(A);  A+=nx;
    }
    const_cast<bool*&>(PP)=reinterpret_cast<bool*>(A);
  }
  //
  template<typename OctTree>
  void PositionsSSE<OctTree>::Update(Access const*T)
  {
    if(T->TREE) {
      Allocate(T->Nleafs());
      unsigned i=0;
      if(OctTree::Dim==2)
	for(Leaf l=T->BeginLeafs(); l!=T->EndLeafs(); ++l) {
	  PP[i] = 0;
	  XX[i] = T->position(l)[0];
	  YY[i] = T->position(l)[1];
	  ++i;
	}
      else
	for(Leaf l=T->BeginLeafs(); l!=T->EndLeafs(); ++l) {
	  PP[i] = 0;
	  XX[i] = T->position(l)[0];
	  YY[i] = T->position(l)[1];
	  ZZ[i] = T->position(l)[2];
	  ++i;
	}
    }
  }
  //
  template<typename OctTree>
  PositionsSSE<OctTree>::PositionsSSE(Access const*T)
    : N16(0), PP(0), XX(0), YY(0), ZZ(0), ALLOC(0), NALLOC(0)
  {
    Update(T);
  }
  //
  template<typename OctTree>
  PositionsSSE<OctTree>::~PositionsSSE()
  {
    if(ALLOC) free16(ALLOC);
    const_cast<char* &>(ALLOC) =0;
    const_cast<size_t&>(NALLOC)=0;
  }
}
//
// Wdutils::FastNeighbourFinder<OctTree>
//
namespace WDutils {
  //
  template<typename OctTree>
  FastNeighbourFinder<OctTree>::FastNeighbourFinder(OctTree const*tree,
						    node_index ndir)
    : NLoop  ( tree, ndir ),
      PosSSE ( this ),
      C0     ( new16<chunk>(PosSSE::N16/K) )
  {}
  //
  template<typename OctTree>
  FastNeighbourFinder<OctTree>::FastNeighbourFinder(Access const&tree,
						    node_index ndir)
    : NLoop  ( tree, ndir ),
      PosSSE ( this ),
      C0     ( new16<chunk>(PosSSE::N16/K) )
  {}
  //
  template<typename OctTree>
  void FastNeighbourFinder<OctTree>::UpdatePositions()
  {
    unsigned n16old = PosSSE::N16;
    PosSSE::Update(this);
    if(n16old < PosSSE::N16 || 2*n16old > 3*PosSSE::N16) {
      if(C0) free16(C0);
      const_cast<chunk*&>(C0) = new16<chunk>(PosSSE::N16/K);
    }
  }
  //
  template<typename OctTree>
  FastNeighbourFinder<OctTree>::~FastNeighbourFinder()
  {
    if(C0) free16(C0); const_cast<chunk*&>(C0)=0;
  }
  //
  template<>
  unsigned FastNeighbourFinder<OctalTree<2,float> >::
  ProcessBlocks(qandi*List, unsigned Nl) const
  {
    __m128 HQ = _mm_set1_ps(S.Q);
    __m128 X0 = _mm_set1_ps(S.X[0]);
    __m128 Y0 = _mm_set1_ps(S.X[1]);
    SSE::Traits<real>::vector QQ;
    unsigned n=0;
    ++CL;
    for(chunk*Ch=C0; Ch!=CL; ++Ch) {
      for(unsigned iK=Ch->I0; iK!=Ch->IN; iK+=K) {
	PP[iK] = 0;
	__m128
	  __D = _mm_sub_ps(X0,_mm_load_ps(XX+iK)),
	  __Q = _mm_mul_ps(__D,__D);
	__D   = _mm_sub_ps(Y0,_mm_load_ps(YY+iK));
	__Q   = _mm_add_ps(__Q,_mm_mul_ps(__D,__D));
	int s = _mm_movemask_ps(_mm_cmplt_ps(__Q,HQ));
	if(s) {
	  _mm_store_ps(QQ,__Q);
	  if(s&1) { if(n<Nl) { List[n].Q=QQ[0]; List[n].I=iK  ; } ++n; }
	  if(s&2) { if(n<Nl) { List[n].Q=QQ[1]; List[n].I=iK+1; } ++n; }
	  if(s&4) { if(n<Nl) { List[n].Q=QQ[2]; List[n].I=iK+2; } ++n; }
	  if(s&8) { if(n<Nl) { List[n].Q=QQ[3]; List[n].I=iK+3; } ++n; }
	}
      }
    }
    return n;
  }
  //
  template<>
  unsigned FastNeighbourFinder<OctalTree<3,float> >::
  ProcessBlocks(qandi*List, unsigned Nl) const
  {
    __m128 HQ = _mm_set1_ps(S.Q);
    __m128 X0 = _mm_set1_ps(S.X[0]);
    __m128 Y0 = _mm_set1_ps(S.X[1]);
    __m128 Z0 = _mm_set1_ps(S.X[2]);
    SSE::Traits<real>::vector QQ;
    unsigned n=0;
    ++CL;
    for(chunk*Ch=C0; Ch!=CL; ++Ch) {
      for(unsigned iK=Ch->I0; iK!=Ch->IN; iK+=K) {
	PP[iK] = 0;
	__m128
	  __D = _mm_sub_ps(X0,_mm_load_ps(XX+iK)),
	  __Q = _mm_mul_ps(__D,__D);
	__D   = _mm_sub_ps(Y0,_mm_load_ps(YY+iK));
	__Q   = _mm_add_ps(__Q,_mm_mul_ps(__D,__D));
	__D   = _mm_sub_ps(Z0,_mm_load_ps(ZZ+iK));
	__Q   = _mm_add_ps(__Q,_mm_mul_ps(__D,__D));
	int s = _mm_movemask_ps(_mm_cmplt_ps(__Q,HQ));
	if(s) {
	  _mm_store_ps(QQ,__Q);
	  if(s&1) { if(n<Nl) { List[n].Q=QQ[0]; List[n].I=iK  ; } ++n; }
	  if(s&2) { if(n<Nl) { List[n].Q=QQ[1]; List[n].I=iK+1; } ++n; }
	  if(s&4) { if(n<Nl) { List[n].Q=QQ[2]; List[n].I=iK+2; } ++n; }
	  if(s&8) { if(n<Nl) { List[n].Q=QQ[3]; List[n].I=iK+3; } ++n; }
	}
      }
    }
    return n;
  }
  //
#ifdef __SSE2__
  template<>
  unsigned FastNeighbourFinder<OctalTree<2,double> >::
  ProcessBlocks(qandi*List, unsigned Nl) const
  {
    __m128d HQ = _mm_set1_pd(S.Q);
    __m128d X0 = _mm_set1_pd(S.X[0]);
    __m128d Y0 = _mm_set1_pd(S.X[1]);
    SSE::Traits<real>::vector QQ;
    unsigned n=0;
    ++CL;
    for(chunk*Ch=C0; Ch!=CL; ++Ch) {
      for(unsigned iK=Ch->I0; iK!=Ch->IN; iK+=K) {
	PP[iK] = 0;
	__m128d
	  __D = _mm_sub_pd(X0,_mm_load_pd(XX+iK)),
	  __Q = _mm_mul_pd(__D,__D);
	__D   = _mm_sub_pd(Y0,_mm_load_pd(YY+iK));
	__Q   = _mm_add_pd(__Q,_mm_mul_pd(__D,__D));
	int s = _mm_movemask_pd(_mm_cmplt_pd(__Q,HQ));
	if(s) {
	  _mm_store_pd(QQ,__Q);
	  if(s&1) { if(n<Nl) { List[n].Q=QQ[0]; List[n].I=iK  ; } ++n; }
	  if(s&2) { if(n<Nl) { List[n].Q=QQ[1]; List[n].I=iK+1; } ++n; }
	}
      }
    }
    return n;
  }
  //
  template<>
  unsigned FastNeighbourFinder<OctalTree<3,double> >::
  ProcessBlocks(qandi*List, unsigned Nl) const
  {
    __m128d HQ = _mm_set1_pd(S.Q);
    __m128d X0 = _mm_set1_pd(S.X[0]);
    __m128d Y0 = _mm_set1_pd(S.X[1]);
    __m128d Z0 = _mm_set1_pd(S.X[2]);
    SSE::Traits<real>::vector QQ;
    unsigned n=0;
    ++CL;
    for(chunk*Ch=C0; Ch!=CL; ++Ch) {
      for(unsigned iK=Ch->I0; iK!=Ch->IN; iK+=K) {
	PP[iK] = 0;
	__m128d
	  __D = _mm_sub_pd(X0,_mm_load_pd(XX+iK)),
	  __Q = _mm_mul_pd(__D,__D);
	__D   = _mm_sub_pd(Y0,_mm_load_pd(YY+iK));
	__Q   = _mm_add_pd(__Q,_mm_mul_pd(__D,__D));
	__D   = _mm_sub_pd(Z0,_mm_load_pd(ZZ+iK));
	__Q   = _mm_add_pd(__Q,_mm_mul_pd(__D,__D));
	int s = _mm_movemask_pd(_mm_cmplt_pd(__Q,HQ));
	if(s) {
	  _mm_store_pd(QQ,__Q);
	  if(s&1) { if(n<Nl) { List[n].Q=QQ[0]; List[n].I=iK  ; } ++n; }
	  if(s&2) { if(n<Nl) { List[n].Q=QQ[1]; List[n].I=iK+1; } ++n; }
	}
      }
    }
    return n;
  }
#endif // __SSE2__
  //
  template<typename OctTree> inline
  void FastNeighbourFinder<OctTree>::AddBlocks(unsigned i0, unsigned iN) const
  {
    if     (CL->IN == i0)  // try to add at end of current block
      CL->IN = iN;
    else if(CL->I0 == iN)  // try to add at begin of current block
      CL->I0 = i0;
    else {                 // must open a new block
      ++CL;
      CL->I0 = i0;
      CL->IN = iN;
    }
  }
  //
  template<typename OctTree> inline
  void FastNeighbourFinder<OctTree>::ProcessLeafs(Leaf L0, Leaf LN) const
  {
    unsigned i0K = L0.I    & nL;     // first block
    unsigned iEK =(LN.I-1) & nL;     // last block
    if(i0K==iEK) {
      // 1   range within a single block
      if(!PP[i0K]) {
	PP[i0K] = 1;
	AddBlocks(i0K,iEK+=K);
      }
    } else {
      // 2   range in more than one block 
      // 2.1 deal with incomplete first block: ignore if already in list
      if(L0.I&L) {
	if(PP[i0K]) i0K    += K;
	else        PP[i0K] = 1;
      }
      // 2.2 deal with incomplete last block: ignore if already in list
      if(LN.I&L) {
	if(PP[iEK]) iEK    -= K;
	else        PP[iEK] = 1;
      }
      // 2.3 add blocks to chunk list of blocks
      if(i0K<=iEK)
	AddBlocks(i0K,iEK+=K);
    }
  }
  //
  template<typename OctTree>
  typename FastNeighbourFinder<OctTree>::node_index
  FastNeighbourFinder<OctTree>::
  Find(Leaf l, real q, Neighbour<OctTree>*nb, node_index m)
  {
    S.Q  = q;
    S.X  = position(l);
    C    = Parent(l);
    while(C != NLoop::Root() && IsValid(C) && Number(Parent(C)) < NDIR )
      C = Parent(C);
    CL     = C0;
    CL->I0 = 0;
    CL->IN = 0;
    NLoop::Process();
    return ProcessBlocks(reinterpret_cast<qandi*>(nb),m);
  }
  //
  template<typename OctTree>
  typename FastNeighbourFinder<OctTree>::node_index
  FastNeighbourFinder<OctTree>::
  Find(point const&x, real q, Neighbour<OctTree>*nb, node_index m)
  {
    S.Q  = q;
    S.X  = x;
    C    = SmallestContainingCell(S.X);
    while(C != NLoop::Root() && IsValid(C) && Number(Parent(C)) < NDIR )
      C = Parent(C);
    CL     = C0;
    CL->I0 = 0;
    CL->IN = 0;
    NLoop::Process();
    return ProcessBlocks(reinterpret_cast<qandi*>(nb),m);
  }
}
#endif
//
// Wdutils::NearestNeighbourFinder<OctTree>
//
namespace {
  /// type used for sorting daughter cells in NearestNeighbourFinder::AddCell
  template<typename OctTree>
  struct CellQ {
    typename TreeAccess<OctTree>::real Q;
    typename TreeAccess<OctTree>::Cell C;
    void set(typename TreeAccess<OctTree>::real q,
	     typename TreeAccess<OctTree>::Cell c)
    { Q=q; C=c; }
    bool operator>(CellQ const&x) const 
    { return Q > x.Q; }
  };
}
//
namespace WDutils {
  template<typename OctTree> inline
  typename NearestNeighbourFinder<OctTree>::node_index
  NearestNeighbourFinder<OctTree>::Ndir() const
  {
    return M<=0? NDIR : max(static_cast<node_index>(M),NDIR);
  }
  //
  template<typename OctTree> inline
  void NearestNeighbourFinder<OctTree>::AddLeaf(Leaf l) const
  {
    // testing after adding the contributions to q from each dimension does
    // make the code run more slowly
    real q = dist_sq(S.X,position(l));
    if(S.Q > q) {
      LIST->Q = q;
      LIST->L = l;
      MaxHeap::after_top_replace(LIST,K);
      S.Q = LIST->Q;
      --M;
      ++NIAC;
    }
  }
  //
  template<typename OctTree>
  void NearestNeighbourFinder<OctTree>::AddCell(Cell Ci, node_index cC) const
  {
    if(cC==0 && Number(Ci) <= Ndir())
      // direct loop
      LoopAllLeafs(Ci,l) AddLeaf(l);
    else {
      // process leaf kids
      if(Nleafkids(Ci))
	LoopLeafKids(Ci,l) AddLeaf(l);
      // process cell kids
      if(Ncells(Ci)>cC+1) {
	// more than one sub-cell: process in order of increasing distance
	CellQ<OctTree> Z[Base::Nsub];
	int J(0);
	LoopCellKids(Ci,c)
	  if(c!=C) Z[J++].set(GeoAlgos::dist_sq(box(c),S.X),c);
	MinHeap::build(Z,J);
	while(J && S.Q > Z->Q) {
	  AddCell(Z->C);
	  Z[0] = Z[J-1];
	  MinHeap::after_top_replace(Z,--J);
	}
      } else if(Ncells(Ci)>cC) {
	// only 1 sub-cell to process
	LoopCellKids(Ci,c)
	  if((!cC || c!=C) && !GeoAlgos::outside(box(c),S))
	    AddCell(c);
      }
    }
  }
  /// comparison of Neighbours on distance, needed below
  template<typename OctTree> inline
  bool operator< (Neighbour<OctTree> const&x, Neighbour<OctTree> const&y)
  { return x.Q < y.Q; }
  //
  template<typename OctTree>
  void NearestNeighbourFinder<OctTree>::FillList()
  {
    NIAC = 0;
    M    = K;
    S.Q  = 12*square(Base::RootRadius());
    for(node_index k=0; k!=K; ++k) {
      LIST[k].L = Base::InvalidLeaf();
      LIST[k].Q = S.Q;
    }
    for(Cell P=C; IsValid(P) && !GeoAlgos::inside(box(C),S); C=P,P=Parent(C))
      AddCell(P, C!=P);
    MaxHeap::sort(LIST,K);
  }
}
//
// instantinations
//
namespace WDutils {
#define  INST(DIM,TYPE)						\
  template class  OctalTree <DIM,TYPE>;				\
  template struct TreeAccess<OctalTree<DIM,TYPE> >;		\
  template struct NeighbourFinder<OctalTree<DIM,TYPE> >;	\
  template struct NearestNeighbourFinder<OctalTree<DIM,TYPE> >;

  INST(2,float)
  INST(2,double)
  INST(3,float)
  INST(3,double)

#if defined(__GNUC__) && defined(__SSE__)
  template struct PositionsSSE<OctalTree<2,float> >;
  template struct PositionsSSE<OctalTree<3,float> >;
  template struct FastNeighbourFinder<OctalTree<2,float> >;
  template struct FastNeighbourFinder<OctalTree<3,float> >;
#ifdef __SSE2__
  template struct PositionsSSE<OctalTree<2,double> >;
  template struct PositionsSSE<OctalTree<3,double> >;
  template struct FastNeighbourFinder<OctalTree<2,double> >;
  template struct FastNeighbourFinder<OctalTree<3,double> >;
#endif
#endif
}
//
