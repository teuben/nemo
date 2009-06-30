// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file    utils/src/octtree.cc
///                                                                             
/// \brief   implements utils/inc/octtree.h
///                                                                             
/// \author  Walter Dehnen                                                      
///                                                                             
/// \date    2009                                                          
///                                                                             
/// \note    based on falcON's tree.cc
///                                                                             
/// \version 08-may-2009 WD  real test: debugged error in linking
/// \version 13-may-2009 WD  abolished Peano-Hilbert support
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
#include <octtree.h>
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

namespace {
  using std::setw;
  using std::setfill;
  using namespace WDutils;
  //////////////////////////////////////////////////////////////////////////////
  // these macros are mainly to avoid confusion for emacs with '<' and '>'.
#define  pDOT(d) static_cast<Dot*>((d))
#define cpDOT(d) static_cast<const Dot*>((d))
#define  pBOX(d) static_cast<Box*>((d))
#define cpBOX(d) static_cast<const Box*>((d))
  //////////////////////////////////////////////////////////////////////////////
  /// base class for Dot and Box
  template<int Dim, typename Real>
  struct Node
  {
    tupel<Dim,Real> X;    ///< position
  };
  /// position and index: same as first two data members of OctalTree::Leaf
  template<int Dim, typename Real>
  struct DotBase : public Node<Dim,Real>
  {
    unsigned I;           ///< index
  };
  /// represents leafs
  template<int Dim, typename Real>
  struct Dot : public DotBase<Dim,Real>
  {
    mutable Dot *Next;    ///< next dot in a linked lisst
    /// add this to a linked list of dots
    void AddToList(Dot* &List, unsigned&Counter)
    {
      Next = List;
      List = this;
      ++Counter;
    }
    /// add this to a linked list of nodes
    void AddToList(Node<Dim,Real>*&List, unsigned&Counter)
    {
      Next = pDOT(List);
      List = this;
      ++Counter;
    }
    /// copy Dot data to Leaf
    template<typename Float>
    void SetLeaf(typename OctalTree<Dim,Real,Float>::Leaf*L) const
    {
      std::memcpy(L,this,sizeof(DotBase<Dim,Real>));
    }
  };
  /// represents cells
  template<int Dim, typename Real>
  struct Box : public Node<Dim,Real>
  {
    typedef ::Node<Dim,Real> Node;
    const static unsigned Nsub = 1<<Dim;
    typedef ::Dot<Dim,Real> Dot;
    uint16 TYP;          ///< bitfield: 1=cell, 0=dot
    uint8  LEV;          ///< tree level of box
    uint8  PEA;          ///< Peano-Hilbert map (currently not used)
    Node  *OCT[Nsub];    ///< octants
    uint32 NUM;          ///< # dots
    Dot   *DOT;          ///< linked list of dots, if any.
    /// is octant i a box?
    bool MarkedAsBox(unsigned i) const { return TYP & (1<<i); }
    /// is octant i a dot?
    bool MarkedAsDot(unsigned i) const { return !MarkedAsBox(i); }
    /// octant of dot within box (not checked)
    inline unsigned octant(const Dot*D) const;
#ifdef TESTING
    /// Is box a single parent, i.e. has only one sub-box
    bool IsSingleParent() const
    {
      if(DOT) return false;
      unsigned n=0;
      for(Node*const*B=OCT; B!=OCT+Nsub; ++B)
	if(*B && ++n>1) return false;
      return true;
    }
#endif
    /// mark octant i as being a box
    void MarkAsBox(unsigned i) { TYP |= (1<<i); }
    /// reset octants
    Box&ResetOctants()
    {
      for(Node**B=OCT; B!=OCT+Nsub; ++B) *B=0;
      return*this;
    }
    /// reset data
    Box&Reset()
    {
      TYP = 0;
      NUM = 0;
      DOT = 0;
      return ResetOctants();
    }
    /// add Dot to linked list
    Box*AddDotToList(Dot*D)
    {
      D->AddToList(DOT,NUM);
      return this;
    }
  };
  //////////////////////////////////////////////////////////////////////////////
  namespace meta {
    template<int Dim, typename Real> struct Helper;
    template<typename Real> struct Helper<2,Real>
    {
      typedef tupel<2,Real> point;
      typedef ::Box<2,Real> Box;
      static unsigned octant(point const&cen, point const&pos)
      {
	unsigned oct(0);
	if(pos[0] > cen[0]) oct |= 1;
	if(pos[1] > cen[1]) oct |= 2;
	return oct;
      }
      static tupel<2,Real> Integer(point const&x)
      {
	tupel<2,Real> c;
	c[0]=int(x[0]+Real(0.5));
	c[1]=int(x[1]+Real(0.5));
	return c;
      }
      static bool ShrinkToOctant(Box*B, unsigned i, uint8 m, const Real*ra)
      {
	uint8 l = ++(B->LEV);
	if(l > m) return false;
	Real rad=ra[l];
	if(i&1) B->X[0] += rad;  else  B->X[0] -= rad;
	if(i&2) B->X[1] += rad;  else  B->X[1] -= rad;
	return true;
      }      
      static Real RootRadius(point const&X, point const&Xmin, point const&Xmax)
      {
	Real D  = max(Xmax[0]-X[0], X[0]-Xmin[0]);
	Real R1 = max(Xmax[1]-X[1], X[1]-Xmin[1]);
	if(R1>D) D=R1;
	return pow(Real(2), int(1+std::log(D)/M_LN2));
      }
    };
    template<typename Real> struct Helper<3,Real>
    {
      typedef tupel<3,Real> point;
      typedef ::Box<3,Real> Box;
      static unsigned octant(point const&cen, point const&pos)
      {
	unsigned oct(0);
	if(pos[0] > cen[0]) oct |= 1;
	if(pos[1] > cen[1]) oct |= 2;
	if(pos[2] > cen[2]) oct |= 4;
	return oct;
      }
      static tupel<3,Real> Integer(point const&x)
      {
	tupel<3,Real> c;
	c[0]=int(x[0]+Real(0.5));
	c[1]=int(x[1]+Real(0.5));
	c[2]=int(x[2]+Real(0.5));
	return c;
      }
      static bool ShrinkToOctant(Box*B, unsigned i, uint8 md, const Real*ra)
      {
	uint8 l = ++(B->LEV);
	if(l > md) return false;
	Real rad=ra[l];
	if(i&1) B->X[0] += rad;  else  B->X[0] -= rad;
	if(i&2) B->X[1] += rad;  else  B->X[1] -= rad;
	if(i&4) B->X[2] += rad;  else  B->X[2] -= rad;
	return true;
      }      
      static Real RootRadius(point const&X, point const&Xmin, point const&Xmax)
      {
	Real D  = max(Xmax[0]-X[0], X[0]-Xmin[0]);
	Real R1 = max(Xmax[1]-X[1], X[1]-Xmin[1]); if(R1>D) D=R1;
	R1 = max(Xmax[2]-X[2], X[2]-Xmin[2]); if(R1>D) D=R1;
	return pow(Real(2), int(1+std::log(D)/M_LN2));
      }
    };
  } // namespace meta
  //////////////////////////////////////////////////////////////////////////////
  template<int Dim, typename Real> inline
  unsigned Box<Dim,Real>::octant(const Dot*D) const
  {
    return meta::Helper<Dim,Real>::octant(this->X,D->X);
  }
  //////////////////////////////////////////////////////////////////////////////
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
  //////////////////////////////////////////////////////////////////////////////
  //
  /// auxiliary class for OctalTree.
  ///
  /// An OctalTree is build by first making a BoxDotTree, then mapping it to
  /// an OctalTree via BoxDotTree::Link()
  template<int Dim, typename Real, typename Float>
  struct BoxDotTree
  {
    const static unsigned Nsub = 1<<Dim; ///< number of octants per cell
    //
    typedef OctalTree<Dim,Real,Float> OctTree;
    typedef typename OctTree::point   point;
    typedef typename OctTree::Leaf    Leaf;
    typedef typename OctTree::Cell    Cell;
    typedef ::Node<Dim,Real>          Node;
    typedef ::Dot <Dim,Real>          Dot;
    typedef ::Box <Dim,Real>          Box;
    /// \name data
    //@{
    unsigned          NMAX;          ///< maximum number dots/box
    uint8             MAXD;          ///< maximum tree depth
    unsigned          NDOT;          ///< number of dots to load
    Dot        *const D0;            ///< begin of dots
    Dot        *const DN;            ///< end of dots
    block_alloc<Box>  BM;            ///< allocator for boxes
    Real             *RA;            ///< array with radius(level)  
    Box              *P0;            ///< root box
    mutable Cell     *CF;            ///< free cells during linking
    mutable Leaf     *LF;            ///< free leafs during linking
    mutable unsigned  NCELL;         ///< # cells linked
    mutable unsigned  DEPTH;         ///< depth of linked tree
    //@}
#ifdef TESTING
    /// Find first non-single parent descendant
    /// \parent P     given box
    /// \return first non-single parent descendant
    Box*NonSingleParent(Box*P) const
    {
      for(;;) {
	if(P->DOT) return P;
	Node**N=P->OCT,**B;
	for(unsigned n=0; N!=P->OCT+Nsub; ++N)
	  if(*N) {
	    if(++n>1) return P;
	    B = N;
	  }
	P = pBOX(*B);
      }
    }
#endif
    /// shrink box to its octant i
    /// \param[in,out] B Box to shrink
    /// \param[in]     i octant
    bool ShrinkToOctant(Box*B, unsigned i)
    {
      return meta::Helper<Dim,Real>::ShrinkToOctant(B,i,MAXD,RA);
    }
    /// get a new box, resetted.
    Box* NewBox(size_t nl)
    {
      return &(BM.new_element(EstimateNalloc(NDOT,nl))->Reset());
    }
    /// provides a new empty (daughter) box in the i th octant of B
    /// \return new box
    /// \param[in] B  parent box
    /// \param[in] i  parent box's octant
    /// \param[in] nl # dots added sofar
    Box*MakeSubBox(const Box*B, unsigned i, size_t nl) WDutils_THROWING
    {
      Box*P = NewBox(nl);
      P->LEV  = B->LEV;
      P->X    = B->X;
      if(!ShrinkToOctant(P,i))
	WDutils_THROW("exceeding maximum tree depth of %d\n         "
		      "(perhaps more than Nmax=%du positions are identical "
		      "within floating-point precision)\n", MAXD, NMAX);
#ifdef OctalTreeSupportPeano
      P->PEA = B->PEA;
      P->PEA.shift_to_kid(i);
#endif
      return P;                                    // return new box            
    }
    /// provides a new (daughter) box in octant i of B containing dot D
    Box*MakeSubBoxN(const Box*B, unsigned i, Dot*D, size_t nl) WDutils_THROWING
    {
      return MakeSubBox(B,i,nl)->AddDotToList(D);
    }
    /// splits a box.
    ///
    /// The Dots in the Box's linked list are sorted into octants. Octants
    /// with one Dot will just hold that Dot, octants with many Dots will be
    /// Boxes with the Dots in the linked list. If all Dots happen to be in
    /// just one octant, the process is repeated on the Box of this octant.
    ///
    /// \param[in,out] P Box to be splitted
    /// \param[in] nl Dots added so far
    void SplitBox(Box*P, size_t nl) WDutils_THROWING
    {
      unsigned NUM[Nsub];
      unsigned b,ne;
      Box*S=0;
      Dot*Di,*Dn;
      do {
	// split linked list into one per octant
	for(b=0; b!=Nsub; ++b)
	  NUM[b] = 0;
	for(Di=P->DOT; Di; Di=Dn) {
	  Dn = Di->Next;
	  b  = P->octant(Di);
	  Di->AddToList(P->OCT[b], NUM[b]);
	}
	P->DOT = 0;
	// loop octants: make box if octant contains > 1 Dot
	for(ne=b=0; b!=Nsub; ++b) if(NUM[b]) {
	  ne++;
	  if(NUM[b]>1) {
	    S      = MakeSubBox(P,b,nl);
	    S->DOT = pDOT(P->OCT[b]);
	    S->NUM = NUM[b];
	    P->OCT[b] = S;
	    P->MarkAsBox(b);
	  }
	}
	P = S;
      } while(ne==1);
    }
    /// box-adding tree building algorithm
    /// \param[in] base base box to add
    /// \param[in] Di dot to add
    /// \param[in] nl # dots added so far
    void AddDotN(Box*base, Dot*Di, size_t nl) WDutils_THROWING
    {
      // loop boxes, starting with root
      for(Box*P=base;;) {
	if(P->DOT) {
	  // box is final: add dot, split box if N > NMAX  -->  done
	  P->AddDotToList(Di);
	  if(P->NUM > NMAX) SplitBox(P,nl);
	  return;
	} else {
	  // non-final box: find octant, increment N
	  unsigned b = P->octant(Di);
	  Node  **oc = P->OCT+b;
	  P->NUM++;
	  if((*oc)==0) {
	    // octant empty: put dot into octant  -->  done
	    *oc = Di;
	    return;
	  } else if(P->MarkedAsDot(b)) {
	    // octant contains another dot: make new box, P=new box
	    Dot*Do = pDOT(*oc);
	    P->MarkAsBox(b);
	    P = MakeSubBoxN(P,b,Do,nl);
	    *oc = P;
	  } else
	    // octant is box: P = box
	    P = pBOX(*oc);
	}
      }
    }
    /// tree linking: leaf & cell descendants are continuous in memory
    /// \param[in,out] C  current cell to be linked
    /// \param[in]     B  current box to link with C
    /// \param[in]     o  octant of box B in parent
    /// \param[in]     k  local peano key
    /// \return           tree depth of cell C
    /// \note recursive.
    /// \node uses data CF and LF
    unsigned LinkCellsS(Cell*C, const Box*B, unsigned o) const WDutils_THROWING
      WD_HOT;
    /// same as LinkCellsS(), but for avoiding single-parent cells
    unsigned LinkCellsA(Cell*C, const Box*B, unsigned o) const WDutils_THROWING
      WD_HOT;
    /// link a final box to a cell
    /// \param[in,out] C  current cell to be linked
    /// \param[in]     B  current box to link with C
    /// \param[in]     o  octant of box B in parent
    /// \param[in]     k  local peano key
    /// \return           tree depth of cell C, always 1 for final boxes
    unsigned LinkFinal(Cell*C, const Box*P, unsigned o) const
    {
      C->Le = P->LEV;
      C->Oc = o;
      C->X  = P->X;
      C->N  = P->NUM;
      C->L0 = LF;
      C->Cf = 0;
      C->Nc = 0;
      C->Nl = P->NUM;
      for(Dot*Di = P->DOT; Di; Di=Di->Next)
	Di-> template SetLeaf<Float>(LF++);
      return 1;
    }
    /// frontend for tree linking
    /// \return      tree depth
    /// \param[in] C array of NBox() cells
    /// \param[in] L array of NDOT leafs
    void Link(bool avspc, Cell*C, Leaf*L) const WDutils_THROWING
    {
      CF    = C+1;
      LF    = L;
      C->Pa = 0;
      DEPTH = avspc? LinkCellsA(C,P0,0) : LinkCellsS(C,P0,0);
      NCELL = unsigned(CF-C);
    }
    /// ctor: build a BoxDotTree
    /// If old tree given, we add the dots in the order of its leafs
    /// \param[in] Ndot number of positions
    /// \param[in] Init initializer to re-initiliase particle data
    /// \param[in] Nmax max positions / cell
    /// \param[in] maxD max tree depth
    /// \param[in] Tree old OctalTree
    BoxDotTree(unsigned Ndot, const typename OctTree::Initialiser*Init,
	       unsigned Nmax, unsigned maxD, const OctTree*Tree=0)
    WDutils_THROWING WD_HOT;
    /// dtor: de-allocate data
    ~BoxDotTree()
    {
      if(D0) WDutils_DEL_A(D0);
      if(RA) WDutils_DEL_A(RA);
    }
    /// number of boxes
    size_t NBox() const { return BM.N_used(); }
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
	out<<" D"<<setfill('0')<<setw(5)<<int(D->Next-D0);
      else
	out<<" nil   ";
      out<<' '<<setfill(' ')<<setw(5)<<D->I<<' '<<setw(10)<<D->X<<'\n';
    }
    /// header for box data dump
    void DumpHeadBox(std::ostream&out)
    {
      out<<" Box        N D     ";
      for(unsigned i=0; i!=Nsub; ++i) out<<" OCT["<<i<<']';
      out<<" S NSP         Rad                 C\n";
    }
    /// dump box data
    void Dump(const Box*B, std::ostream&out, unsigned&ns)
    {
      out<<" B"<<setfill('0')<<setw(5)<<BM.number_of_element(B)
	 <<' ' <<setfill(' ')<<setw(5)<<B->NUM;
      if(B->DOT)
	out<<" D"<<setfill('0')<<setw(5)<<int(B->DOT-D0);
      else
	out<<" nil   ";
      for(unsigned i=0; i!=Nsub; ++i) {
	if(B->OCT[i] == 0)
	  out<<" nil   ";
	else if(B->MarkedAsBox(i))
	  out<<" B"<<setfill('0')<<setw(5)
	     <<BM.number_of_element(pBOX(B->OCT[i]));
	else
	  out<<" D"<<setfill('0')<<setw(5)
	     <<int(pDOT(B->OCT[i])-D0);
      }
      if(B->IsSingleParent()) {
	++ns;
	out<<" S";
      } else
	out<<" M";
      out<<" B"<<setfill('0')<<setw(5)
	 <<BM.number_of_element(NonSingleParent(const_cast<Box*>(B)));
      out<<' '<<setfill(' ')<<setw(8)<<RA[B->LEV]
	 <<' '<<setw(8)<<B->X<<'\n';
    }
    /// dump tree, recursive
    void Dump(const Box*B, std::ostream&outd, std::ostream&outb, unsigned&ns)
    {
      Dump(B,outb,ns);
      if(B->DOT)
	for(Dot*Di=B->DOT; Di; Di=Di->Next)
	  Dump(Di,outd);
      else {
	for(unsigned i=0; i!=Nsub; ++i)
	  if(B->OCT[i] && !B->MarkedAsBox(i))
	    Dump(pDOT(B->OCT[i]),outd);
	for(unsigned i=0; i!=Nsub; ++i)
	  if(B->OCT[i] && B->MarkedAsBox(i))
	    Dump(pBOX(B->OCT[i]),outd,outb,ns);
      }
    }
    /// dump tree
    /// \param[in] outd ostream for dumping Dot data
    /// \param[in] outb ostream for dumping Box data
    void Dump(std::ostream&outd, std::ostream&outb) {
      unsigned ns=0;
      DumpHeadDot(outd);
      DumpHeadBox(outb);
      Dump(P0,outd,outb,ns);
      std::cerr<<" # single parent boxes: "<<ns<<'\n';
    }
#endif // TESTING
  };
  //
  template<int D, typename Real, typename Float>
  unsigned BoxDotTree<D,Real,Float>::LinkCellsS(Cell*C, const Box*P, unsigned o)
    const WDutils_THROWING
  {
    // final box: use LinkFinal()
    if(P->DOT)
      return LinkFinal(C,P,o);
    // set some data
    C->L0 = LF;
    C->Oc = o;
    // loop octants: count dots & boxes
    unsigned i,nbox,ndot;
    const Node*const*N;
    nbox=0,ndot=0;
    for(i=0,N=P->OCT; i!=Nsub; ++i,++N)
      if(*N) {
	if(P->MarkedAsBox(i)) {
	  ++nbox;
	} else {
	  ++ndot;
	  cpDOT(*N)-> template SetLeaf<Float>(LF++); 
	}
      }
    // copy more data
    C->Le = P->LEV;
    C->X  = P->X;
    C->N  = P->NUM;
    C->Nl = ndot;
    C->Nc = nbox;
    if(nbox) {
      unsigned dep = 1;
      C->Cf = CF-C;
      Cell*Ci = CF;
      CF += nbox;
      for(i=0,N=P->OCT; i!=Nsub; ++i,++N)
	if(*N && P->MarkedAsBox(i)) {
	  Ci->Pa = Ci-C;
	  unsigned de = LinkCellsS(Ci++,cpBOX(*N),i);
          if(de>dep) dep=de;
        }
      return ++dep;
    } else {
      C->Cf = 0;
      if(NMAX > unsigned(Nsub))
	WDutils_THROW("LinkCells: found non-final box without daughters\n");
      return 1;
    }
  }
  //
  template<int D, typename Real, typename Float>
  unsigned BoxDotTree<D,Real,Float>::LinkCellsA(Cell*C, const Box*P, unsigned o)
    const WDutils_THROWING
  {
    // final box: use LinkFinal()
    if(P->DOT)
      return LinkFinal(C,P,o);
    // set some data
    C->L0 = LF;
    C->Oc = o;
    // loop octants: count dots & boxes, replace P by single-parent descendant
    unsigned i,nbox,ndot;
    const Node*const*N;
    for(const Box *B=0;;) {
      nbox=0,ndot=0;
      for(i=0,N=P->OCT; i!=Nsub; ++i,++N)
	if(*N) {
	  if(P->MarkedAsBox(i)) {
	    ++nbox;
	    B = cpBOX(*N);
	  } else {
	    ++ndot;
	    cpDOT(*N)-> template SetLeaf<Float>(LF++); 
	  }
	}
      if(ndot || nbox>1) break;
      P = B;
      if(P->DOT) return LinkFinal(C,P,o);
    }
    // copy some simple data
    C->Le = P->LEV;
    C->X  = P->X;
    C->N  = P->NUM;
    C->Nl = ndot;
    C->Nc = nbox;
    if(nbox) {
      unsigned dep = 1;
      C->Cf = CF-C;
      Cell*Ci = CF;
      CF += nbox;
      for(i=0,N=P->OCT; i!=Nsub; ++i,++N)
	if(*N && P->MarkedAsBox(i)) {
	  Ci->Pa = Ci-C;
	  unsigned de = LinkCellsA(Ci++,cpBOX(*N),i);
          if(de>dep) dep=de;
        }
      return ++dep;
    } else {
      C->Cf = 0;
      if(NMAX > unsigned(Nsub))
	WDutils_THROW("LinkCells: found non-final box without daughters\n");
      return 1;
    }
  }
  //
  template<int D, typename Real, typename Float>
  BoxDotTree<D,Real,Float>::BoxDotTree(unsigned Ndot,
				       const typename OctTree::Initialiser*Init,
				       unsigned Nmax, unsigned maxD,
				       const OctTree*Tree)
    WDutils_THROWING
    : NMAX(Nmax), MAXD(maxD),
      NDOT(Ndot), D0(WDutils_NEW(Dot,NDOT)), DN(D0+NDOT),
      BM  (Tree? Tree->Ncells() : 1+NDOT/4),
      RA  (WDutils_NEW(Real,MAXD+1))
  {
    if(NDOT == 0) WDutils_THROW("OctalTree: N=0\n");
    point Xmin, Xmax, Xave;
    // 1  set dots
    if(Tree && Tree->Nleafs()) {
      // 1.1 in old tree order
      const Leaf*Li=Tree->BeginLeafs(), *LN=Tree->EndLeafs();
      D0->I = (Li++)->I;
      Init->ReInit(reinterpret_cast<typename OctTree::Dot*>(D0));
      if(isnan(D0->X)) WDutils_THROW("OctalTree: position contains NaN\n");
      Xmin = Xmax = Xave = D0->X;
      Dot*Di=D0+1;
      for(; Di!=DN && Li!=LN; ++Di) {
	Di->I = (Li++)->I;
	Init->ReInit(reinterpret_cast<typename OctTree::Dot*>(Di));
	if(isnan(Di->X)) WDutils_THROW("OctalTree: position contains NaN\n");
	Di->X.up_min_max(Xmin,Xmax);
	Xave += Di->X;
      }
      for(; Di!=DN; ++Di) {
	Init->Init(reinterpret_cast<typename OctTree::Dot*>(Di));
	if(isnan(Di->X)) WDutils_THROW("OctalTree: position contains NaN\n");
	Di->X.up_min_max(Xmin,Xmax);
	Xave += Di->X;
      }
    } else {
      // 1.2 from scratch: loop dots: initialise, find Xmin, Xmax
      Init->Init(reinterpret_cast<typename OctTree::Dot*>(D0));
      if(isnan(D0->X)) WDutils_THROW("OctalTree: position contains NaN\n");
      Xmin = Xmax = Xave = D0->X;
      for(Dot*Di=D0+1; Di!=DN; ++Di) {
	Init->Init(reinterpret_cast<typename OctTree::Dot*>(Di));
	if(isnan(Di->X)) WDutils_THROW("OctalTree: position contains NaN\n");
	Di->X.up_min_max(Xmin,Xmax);
	Xave += Di->X;
      }
    }
    // 2 set root box, RA[]
    Xave   /= Real(NDOT);
    P0      = NewBox(1);
    P0->X   = meta::Helper<D,Real>::Integer(Xave);
    P0->LEV = 0;
    RA[0]   = meta::Helper<D,Real>::RootRadius(P0->X,Xmin,Xmax);
    for(unsigned l=0; l!=MAXD; ++l)
      RA[l+1] = Real(0.5)*RA[l];
#ifdef OctalTreeSupportPeano
    P0->PEA.set_root();
#endif
    // add dots
    size_t nl = 0;
    for(Dot*Di=D0; Di!=DN; ++Di,++nl)
      AddDotN(P0,Di,nl);
#ifdef TESTING
    std::ofstream dumpD("dots.dat"), dumpB("boxs.dat");
    Dump(dumpD, dumpB);
#endif
  }
}
//
namespace WDutils {
  template<int D, typename Real, typename Float>
  void OctalTree<D,Real,Float>::Allocate()
  {
    unsigned need =
      NLEAF*sizeof(Leaf) +
      NCELL*sizeof(Cell) +
      (MAXD+1)*sizeof(Real);
    if((need > NALLOC) || (need+need < NALLOC)) {
      if(ALLOC) delete16(ALLOC);
      ALLOC  = new16<char>(need);
      NALLOC = need;
    }
    char* A = ALLOC;
    LEAFS=reinterpret_cast<Leaf*>(A); A+=NLEAF*sizeof(Leaf);
    CELLS=reinterpret_cast<Cell*>(A); A+=NCELL*sizeof(Cell);
    RAD  =reinterpret_cast<Real*>(A);
  }
  //
  template<int D, typename Real, typename Float>
  OctalTree<D,Real,Float>::OctalTree(unsigned n, const Initialiser*init,
				     unsigned nmax, bool avspc, unsigned maxd)
    WDutils_THROWING
    : INIT(init), MAXD(maxd), AVSPC(avspc), NMAX(nmax), ALLOC(0), NALLOC(0)
  {
    if(NMAX<2)
      WDutils_THROW("OctalTree<%d,%s>: nmax=%du < 2\n",D,nameof(Real),nmax);
    BoxDotTree<D,Real,Float> BDT(n,INIT,NMAX,MAXD);
    NLEAF = BDT.NDOT;
    NCELL = BDT.NBox();
    Allocate();
    BDT.Link(AVSPC,CELLS,LEAFS);
    DEPTH = BDT.DEPTH;
    NCELL = BDT.NCELL;
    std::memcpy(RAD,BDT.RA,(MAXD+1)*sizeof(Real));
//     for(unsigned d=0; d<=MAXD; ++d)
//       RAD[d] = BDT.RA[d];
  }
  //
  template<int D, typename Real, typename Float>
  void OctalTree<D,Real,Float>::rebuild(unsigned n, unsigned nmax)
    WDutils_THROWING
  {
    if(nmax) NMAX = nmax;
    if(NMAX<2)
      WDutils_THROW("OctalTree<%d,%s>: nmax=%du < 2\n",D,nameof(Real),nmax);
    BoxDotTree<D,Real,Float> BDT(n? n:NLEAF,INIT,NMAX,MAXD,this);
    NLEAF = BDT.NDOT;
    NCELL = BDT.NBox();
    Allocate();
    BDT.Link(AVSPC,CELLS,LEAFS);
    DEPTH = BDT.DEPTH;
    NCELL = BDT.NCELL;
    std::memcpy(RAD,BDT.RA,(MAXD+1)*sizeof(Real));
//     for(unsigned d=0; d<=MAXD; ++d)
//       RAD[d] = BDT.RA[d];
  }
  //
  template<int D, typename Real, typename Float>
  OctalTree<D,Real,Float>::~OctalTree()
  {
    if(ALLOC) delete16(ALLOC);
    ALLOC = 0;
    NALLOC= 0;
    CELLS = 0;
    LEAFS = 0;
    RAD   = 0;
  }
}
//
namespace WDutils {
  template class OctalTree<2,float,float>;
  template class OctalTree<2,double,float>;
  template class OctalTree<2,double,double>;
  template class OctalTree<3,float,float>;
  template class OctalTree<3,double,float>;
  template class OctalTree<3,double,double>;
}
//
