// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file    utils/src/octtree.cc
/// 
/// \brief   implements utils/inc/octtree.h
///
/// \author  Walter Dehnen
///
/// \date    2009-2012
///
/// \note    originally based on falcON's tree.cc (by the same author)
/// 
/// \version 08-may-2009 WD  real test: debugged error in linking
/// \version 13-may-2009 WD  abolished Peano-Hilbert support
/// \version 25-sep-2009 WD  new version using indices for Leaf & Cell
/// \version 14-oct-2009 WD  new version tested against old, old abolished.
/// \version 27-jan-2010 WD  added leaf's parent cell, changed Node magic
/// \version 28-jan-2010 WD  TreeAccess::SmallestContainingCell tested
/// \version 26-feb-2010 WD  new initialiser interface
/// \version 15-apr-2010 WD  tree pruning, OctalTree::build()
/// \version 22-apr-2010 WD  parameter nmin, changes in tree building
/// \version 24-apr-2010 WD  changes in tree building code.
/// \version 09-jun-2010 WD  16-byte alignement, using geometry.h
/// \version 15-jul-2010 WD  maximum tree depth = numeric_limits<real>::digits
/// \version 02-jul-2012 WD  added periodic boundary condition
/// \version 10-aug-2012 WD  completely new, OMP parallel in non-public part
///                          requires C++11
/// \version 10-aug-2012 WD  removed OMP stuff to non-public part of library
/// \version 30-aug-2012 WD  replaced typedef with using
/// \version 19-sep-2012 WD  alternative memory layout using std::vector
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009-2012 Walter Dehnen
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
#ifndef WD_HOT
# if defined(__GNUC__) && (__GNUC__ > 3) && (__GNUC_MINOR__ > 1)
#  define WD_HOT __attribute__((hot))
# else
#  define WD_HOT
# endif
#endif
//
#include <octtree.h>
#include <peano.h>
#include <timer.h>
#define WDutils_octtree_cc
//
#define TESTING
#undef  TESTING

#ifdef TESTING
# include <fstream>
# warning compilation for TESTING purposes
#endif

#ifdef OCTALTREE_USE_OPENMP
# define WALLCLOCK OMP::WallClock
#else
# define WALLCLOCK RunInfo::WallClock
#endif
//
namespace {
  using std::setw;
  using std::setfill;
  using namespace WDutils;
  using namespace WDutils::octtree;
  /// some general auxiliaries
  template<int D, typename _Tp> struct TreeHelper;
  template<typename _Tp> struct TreeHelper<2,_Tp>
  {
    using point = typename Geometry::cube<2,_Tp>::point;
    static point RootCentre(point const&x)
    {
      point c;
      c[0]=int(x[0]+_Tp(0.5));
      c[1]=int(x[1]+_Tp(0.5));
      return c;
    }
    static _Tp RootRadius(point const&X,
			  point const&Xmin,
			  point const&Xmax)
    {
      _Tp D  = max(Xmax[0]-X[0], X[0]-Xmin[0]);
      _Tp R1 = max(Xmax[1]-X[1], X[1]-Xmin[1]);
      if(R1>D) D=R1;
      return pow(_Tp(2), int(1+std::log(D)/M_LN2));
    }
  };
  //
  template<typename _Tp> struct TreeHelper<3,_Tp>
  {
    using point = typename Geometry::cube<3,_Tp>::point;
    static point RootCentre(point const&x)
    {
      point c;
      c[0]=int(x[0]+_Tp(0.5));
      c[1]=int(x[1]+_Tp(0.5));
      c[2]=int(x[2]+_Tp(0.5));
      return c;
    }
    static _Tp RootRadius(point const&X,
			  point const&Xmin,
			  point const&Xmax)
    {
      _Tp D  = max(Xmax[0]-X[0], X[0]-Xmin[0]);
      _Tp R1 = max(Xmax[1]-X[1], X[1]-Xmin[1]); if(R1>D) D=R1;
      R1 = max(Xmax[2]-X[2], X[2]-Xmin[2]); if(R1>D) D=R1;
      return pow(_Tp(2), int(1+std::log(D)/M_LN2));
    }
  };
  /// estimate number of boxes needed for box-dot tree
  /// \param[in] Ndot  # dots
  /// \param[in] Nmax  max # dots/octant
  template<int Dim>
  inline size_t NBoxes(size_t Ndot, uint8_t Nmax)
  {
    double
      x = double(400*Nmax+Ndot)/double(100*Nmax+Ndot);
    x  *= Ndot;
    x  /= Nmax+Nmax;
    x  += 4*std::sqrt(x);
    if(x<100) x=100;
    return size_t(Dim==2? 1.5*x : x);
  }
  /// the next multiple of 16 to n*sizeof(T)
  template<typename T> inline size_t Next16(size_t n)
  { return WDutils::next_aligned16(n*sizeof(T)); }
  //
#ifdef OCTALTREE_DATA_IN_ONE_BLOCK
  template<typename _Tp>
  static _Tp*p_data(_Tp*X) { return X; }
#else
  template<typename _Tp>
  static _Tp*p_data(octtree::storage<_Tp> const&X)
  { return const_cast<_Tp*>(X.data()); }
  template<typename _Tp>
  static _Tp*p_data(octtree::store16<_Tp> const&X)
  { return const_cast<_Tp*>(X.data()); }
#endif
  //
  /// tree of boxes and dots, to be used inside an openMP parallel region
  //
  template<int __D, typename __X> 
  struct TreeImplementer< OctalTree<__D,__X> >
  {
    //
    using OTree            = OctalTree<__D,__X>;
    using Helper           = TreeHelper<__D,__X>;
    using Dot              = typename OTree::Dot;
    using Initialiser      = typename OTree::Initialiser;
    using GeoAlg           = Geometry::Algorithms<1>;
    using pos_type         = typename OTree::pos_type;
    using cube             = typename OTree::cube;
    using cube16           = typename OTree::cube16;
    using point            = typename Helper::point;
    using point16          = typename OTree::point16;
    using PerBoundary      = typename OTree::PerBoundary;
    using BuilderInterface = typename OTree::BuilderInterface;
    using peano_map        = typename Peano<__D>::Map;
    using peano_dig        = typename Peano<__D>::Digit;
    using peano_oct        = typename Peano<__D>::Octant;
    //
    static const peano_oct  Dim          = __D;
    static const peano_oct  Nsub         = 1<<Dim;
    static const depth_type MaximumDepth = OTree::MaximumDepth;
    /// \name types used in tree construction
    //@{
    /// represents a cubic cell in the box-dot tree
    /// \note We cannot rely on a default constructor being called when
    ///       allocating boxes for the tree (see documentation for block_alloc).
    /// \note This is carefully designed to have no more than 64, 80, 96, and
    ///       112 bytes for (Dim,pos_type)=(2,float), (2,double), (3,float),
    ///       and (3,double), respectively.
    struct __Box {
      using ndl_type = typename meta::IntTypeWords<Nsub>::integer_u;
      cube        CUB;                 ///< box cubus
      void       *OCT[Nsub];           ///< octants
      // if 0 == NDL[i]          the octant is empty
      // if 0 <  NDL[i] <= NMAX  the octant holds a list of NDL[] leaves
      // if      NDL[i]  > NMAX  the octant holds a box
      union {
	ndl_type  NDl;
	uint8_t   NDL[Nsub];           ///< # dots/octant
      };
      count_type  NUM;                 ///< total # dots in box
      uint8_t     NOC;                 ///< # octants occupied
      uint8_t     NBX;                 ///< # sub-boxes
      uint8_t     LEV;                 ///< tree level of box
      peano_map   MAP;                 ///< peano map octant <-> digit
    };
    using Box = SSE::Extend16<__Box>;
    //
#define  pDOT(d) static_cast<Dot*>((d))
#define cpDOT(d) static_cast<const Dot*>((d))
#define  pBOX(d) static_cast<Box*>((d))
#define cpBOX(d) static_cast<const Box*>((d))
    //
#ifdef OCTALTREE_USE_OPENMP
# define OCTALTREE_IMPLEMENTER_ADDITIONAL_MEMBER_TYPES
# include <../devel/octtree_omp.cc>
# undef  OCTALTREE_IMPLEMENTER_ADDITIONAL_MEMBER_TYPES
#endif
    //@}
    /// \name data set by constructor
    //@{
    mutable double        TIME;         ///< wallclock time at start
    const uint8_t         NMAX;         ///< maximum # particles/octant
    const uint8_t         NMIN;         ///< maximum # particles/octant
    const bool            ASCC;         ///< avoid single-child cells
#ifdef OCTALTREE_USE_OPENMP
    depth_type            IT,NT;        ///< this thread, # threads
#endif// OCTALTREE_USE_OPENMP
    const count_type      NMX1;         ///< NMAX + 1
    //@}
    /// \name data set by SetDots()
    //@{
    Dot                  *DOTS;         ///< begin: locally loaded dots
    count_type            NDOT;         ///< # locally loaded dots
    //@}
    /// \name data set by Build()
    //@{
    block_alloc<Box>     *BXAL;         ///< allocator for boxes
    Box                  *ROOT;         ///< root box
    //@}
#ifdef OCTALTREE_USE_OPENMP
# define OCTALTREE_IMPLEMENTER_ADDITIONAL_DATA_MEMBERS
# include <../devel/octtree_omp.cc>
# undef  OCTALTREE_IMPLEMENTER_ADDITIONAL_DATA_MEMBERS
#endif
    OTree                *TREE;         ///< tree to be linked
    mutable count_type    CF;           ///< cell counter during linking
    mutable count_type    LF;           ///< leaf counter during linking
#ifdef OCTALTREE_USE_OPENMP
    mutable count_type    BF;           ///< branch counter during linking
#endif// OCTALTREE_USE_OPENMP
    //@}
    /// set dots; cumulate Xmin, Xmax, Xsum; called from Build()
    /// \param[in]  init  initialiser for particle data
    /// \param[in]  tree  domain to take order from
    /// \param[out] Xmin  in each dimension: minimum of positions loaded
    /// \param[out] Xmin  in each dimension: maximum of positions loaded
    /// \param[out] Xsum  sum of all positions loaded
    /// \note sets @a DOTS and @a NDOT
    void SetDots(const Initialiser*init, const OTree*tree,
		 point&Xmin, point&Xmax, point&Xave) WD_HOT;
    /// build tree, called from ctor
    void Build(const Initialiser*init, const OTree*tree) WD_HOT;
    /// replace (pter to) box by (pter to) its non-single-child descendant
    /// \param[in,out] P  (pter to) given box
    void EnsureNonSingleChildBox(const Box*&P) const
    {
      while(P->NOC==1 && P->NBX==1)
	for(int i=0; i!=Nsub; ++i)
	  if(P->NDL[i]) {
	    P = pBOX(P->OCT[i]);
	    break;
	  }
    }
    /// provides a new empty (daughter) box in the i th octant of B
    /// \return new box
    /// \param[in] B  parent box
    /// \param[in] i  parent box's octant
    Box*MakeSubBox(const Box*B, int i)
    {
      Box *P = BXAL->new_element();
      P->NDl = 0;
      P->NOC = 0;
      P->NBX = 0;
      P->MAP = B->MAP.daughter(i);
      P->LEV = B->LEV + 1;
      GeoAlg::copy(B->CUB,P->CUB);
      GeoAlg::shrink(P->CUB,i);
      return P;
    }
    /// adds dot @a Di to octant @a b of box @a P
    void AddDotToOctantList(Box*P, Dot*Di, int b)
    {
      if(P->NDL[b] == 0) {
	P->NOC ++;
	Di->Next = 0;
      } else
	Di->Next = P->OCT[b];
      P->OCT[b]  = Di;
      P->NDL[b] ++;
    }
    /// adds a dot to given octant of given box, deepens box if necessary
    /// \param[in] P   box to add to
    /// \param[in] Di  dot to add
    /// \parma[in] b   octant to add to
    /// \note The octant must not contain a box already (P->NDL[b] <= NMAX)
    void AddDotToOctant(Box*P, Dot*Di, int b)
    {
      AddDotToOctantList(P,Di,b);          // add dot to octant in box
      while(P->NDL[b] > NMAX) {            // WHILE dot list too long
	Di        = pDOT(P->OCT[b]);       //   get first dot from octant
	P->NBX++;                          //   increment box counter
	P->OCT[b] = MakeSubBox(P,b);       //   new empty box in octant
	P         = pBOX(P->OCT[b]);       //   P = daughter box
	if(P->LEV>= MaximumDepth) {
#ifdef TESTING
	  DumpNodes();
#endif
	  std::ostringstream out;
	  for(;Di; Di=pDOT(Di->Next))
	    out<<"         particle #"<<std::setw(6)<<Di->I
	       <<" @ "<<std::setw(10)<<Di->X<<'\n';
	  WDutils_ErrorN("exceeding maximum tree depth of %d\n         "
			 "perhaps more than Nmax=%d positions are identical:\n"
			 "%s", MaximumDepth,NMAX,out.str().c_str());
	}
	P->NUM    = NMX1;                  //   set P->NUM
	for(Dot*Dn; Di; Di=Dn) {           //   sort dots into octants
	  Dn= pDOT(Di->Next);
	  b = GeoAlg::octant(P->CUB,Di->X);
	  AddDotToOctantList(P,Di,b);
	}
      }
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
    void AddDotToBox(Box*P, Dot*Di)
    {
      for(;;) {
	P->NUM++;
	int b = GeoAlg::octant(P->CUB,Di->X);
	if(P->NDL[b] > NMAX) P = pBOX(P->OCT[b]);
	else return AddDotToOctant(P,Di,b);
      }
    }
    /// add list of dots to box
    void AddDotListToBox(Box*P, Dot*Di)
    {
      for(Dot*Dn; Di; Di=Dn) {
	Dn= pDOT(Di->Next);
	AddDotToBox(P,Di);
      }
    }
    /// add list of dots to octant in box
    /// \note We must allow for the octant to be(come) a box.
    void AddDotListToOctant(Box*P, Dot*Di, int b)
    {
      for(Dot*Dn; Di; Di=Dn) {
	Dn= pDOT(Di->Next);
	if(P->NDL[b] > NMAX) AddDotToBox(pBOX(P->OCT[b]),Di);
	else AddDotToOctant(P,Di,b);
      }
    }
    // 
#ifdef OCTALTREE_USE_OPENMP
# define OCTALTREE_IMPLEMENTER_ADDITIONAL_MEMBER_FUNC_DECL
# include <../devel/octtree_omp.cc>
# undef  OCTALTREE_IMPLEMENTER_ADDITIONAL_MEMBER_FUNC_DECL
#endif 
    /// copy Dot data to Leaf
    /// \param[in] L  leaf index in tree
    /// \param[in] D  dot to link leaf to
    /// \param[in] P  parent cell index for leaf
    void LinkLeaf(const count_type L, const Dot*const D, const count_type P)
      const
    {
      static_cast<point&>(TREE->_XL[0][L]) = D->X;
      TREE->_PL[0][L] = D->I;
      TREE->_PC[L] = P;
    }
    /// tree linking: make a cell from an octant of a parent box
    /// \param[in] P    parent box
    /// \param[in] C    index of current cell to be linked
    /// \param[in] oct  octant of cell in @a P
    void LinkOctantCell(const Box*const P, const count_type C, octant_type oct,
			const count_type 
#ifdef TESTING
			Cn
#endif
			) const
    {
#ifdef TESTING
      if(C > Cn)
	WDutils_ErrorN("error in cell linking: C=%d > NC=%d\n",C,Cn);
#endif
      TREE->_LE[C] = P->LEV + 1;
      GeoAlg::copy(P->CUB,static_cast<cube&>(TREE->_XC[C]));
      GeoAlg::shrink(TREE->_XC[C],oct);
      TREE->_L0[C] = LF;
      TREE->_NL[C] = P->NDL[oct];
      TREE->_NM[C] = P->NDL[oct];
      TREE->_C0[C] = C;
      TREE->_NC[C] = 0;
#ifdef OCTALTREE_HAVE_D0
      TREE->_D0[C] = IT;
#endif
      TREE->_DP[C] = 0;
      for(const Dot*Di=cpDOT(P->OCT[oct]); Di; Di=pDOT(Di->Next))
	LinkLeaf(LF++,Di,C);
    }
    //
    // PassUpNactive
    //
    static void PassUpNactive(OTree*const tree,
			      bool all_active, count_type C0, count_type CN)
    {
      if(all_active)
	for(count_type c=C0; c!=CN; ++c)
	  tree->_NA[c] = tree->_NM[c];
      else
	for(count_type c=CN-1; c!=(C0-1); --c) {
	  tree->_NA[c] = 0;
	  if(tree->_NL[c])
	    for(count_type ll=tree->_L0[c]; ll!=tree->_L0[c]+tree->_NL[c]; ++ll)
	      if(tree->_FL[0][ll]) tree->_NA[c]++;
	  if(tree->_NC[c])
	    for(count_type cc=tree->_C0[c]; cc!=tree->_C0[c]+tree->_NC[c]; ++cc)
	      tree->_NA[c] += tree->_NA[cc];
	}
    }
    /// link ordinary cell
    /// \param[in] P   box to link
    /// \param[in] C   index for domain cell
    /// \param[in] cN  maximum allowed value for @a C
    depth_type LinkCell(const Box*B, const count_type C, const count_type cN)
      const WD_HOT;
    /// new ctor: build a TreeImplementer
    /// \param[in] wt    wallclock time
    /// \param[in] init  Initialiser for dots
#ifdef OCTALTREE_USE_OPENMP
    /// \param[in] iT    rank of this thread
    /// \param[in] nT    # parallel threads
#endif
    /// \param[in] tree  if non-null: tree to take dot order from
    /// \param[in] nmax  max # dots / unsplit octant
    /// \param[in] nmin  min # leaf / cell
    /// \param[in] ascc  avoid single-child cell
    TreeImplementer(double wt, const Initialiser*init, 
#ifdef OCTALTREE_USE_OPENMP
	       count_type iT, count_type nT,
#endif
	       const OTree*tree, depth_type nmax, depth_type nmin, bool ascc)
    : TIME  ( wt )
    , NMAX  ( nmax )
    , NMIN  ( nmin )
    , ASCC  ( ascc )
#ifdef OCTALTREE_USE_OPENMP
    , IT    ( iT )
    , NT    ( nT )
#endif
    , NMX1  ( nmax + 1 )
    , DOTS  ( 0 )
    , BXAL  ( 0 )
    , ROOT  ( 0 )
#ifdef OCTALTREE_USE_OPENMP
    , TBAL  ( 0 )
#endif
#ifdef TESTING
    , TRES  ( 0 )
#endif
    {
      WDutilsAssert(NMAX >= NMIN);
      Build(init,tree);
    }
    /// count boxes and dots in box-dot tree
    void CountNodesSerial(count_type&nbox, count_type&ndot) const;
#ifdef OCTALTREE_USE_OPENMP
    /// count boxes and dots in 'our' domain of the tree
    /// \note branch boxes are not counted here
    void CountNodesDomain(count_type&nbox, count_type&ndot) const;
#endif
    /// link serial
    /// \param[in,out] tree  tree to link
    void LinkSerial(OTree*tree) const;
    /// dtor
    ~TreeImplementer();
    /// allocate data for an OctalTree
    static void AllocateTree(OTree*tree, count_type nc, count_type nl,
			     count_type ne, bool fr) WDutils_THROWING;
    /// serial build of an OctalTree
    static void BuildTreeSerial(OTree*tree,bool fresh,
				const BuilderInterface*builder);
    /// build an OctalTree
    static void BuildTree(OTree*tree,bool fresh, count_type 
#ifdef OCTALTREE_USE_OPENMP
			  tol
#endif
			  , const BuilderInterface*builder=0)
    {
#ifdef OCTALTREE_USE_OPENMP
      if(tree->NDOM>1)
	BuildTreeParallel(tree,fresh,tol,builder);
      else 
#endif
	BuildTreeSerial  (tree,fresh,builder);
    }
#ifdef TESTING
    /// header for dot data dump
    void DumpHeadDot(std::ostream&out)
    { out<<"Dot          Next        I                     X\n"; }
    /// print dot
    std::ostream&Print(const Dot*D, std::ostream&out)
    {
      if(D==0)
	return out <<"Dunknown";
#ifdef OCTALTREE_USE_OPENMP
      if(TRES!=0) {
	depth_type i=0;
	for(; i!=NT; ++i)
	  if(D >= TRES[i]->DOTS && TRES[i]->DOTS+TRES[i]->ROOT->NUM > D) {
	    return out<<'D'<<setfill('0')
		      <<setw(2)<<int(i)<<'.'
		      <<setw(4)<<int(D-TRES[i]->DOTS)
		      <<setfill(' ');
	  }
	return out<<"D(X="<<setfill(' ')<<setw(8)<<D->X<<')';
      } else 
#endif// OCTALTREE_USE_OPENMP
	return out<<'D'<<setfill('0')<<setw(7)<<int(D-DOTS)
		  <<setfill(' ');
    }
    /// dump dot data
    void Dump(const Dot*D, std::ostream&out)
    {
      Print(D,out) <<' ';
      if(D->Next) Print(pDOT(D->Next),out) << ' ';
      else out<<" nil    ";
      out<<setfill(' ')<<setw(5)<<D->I<<' '<<setw(10)<<D->X<<'\n';
      if(D->Next)
	Dump(pDOT(D->Next),out);
    }
    /// header for box data dump
    void DumpHeadBox(std::ostream&out)
    {
      out<<" Box          N No Nb Mp Le"
	 <<(Dim==2? "  NDL[]   " : "        NDL[]         ");
      for(int i=0; i!=Nsub; ++i) out<<"   OCT["<<i<<']';
      out<<"   S NSP           Rad                 C\n";
    }
    /// Find first non-single parent descendant
    const Box* NonSingleChildBox(const Box*B) const
    {
      const Box*P = B;
      EnsureNonSingleChildBox(P);
      return P;
    }
    /// running number of local box
    int BoxLocalNo(const Box*B) const
    { return BXAL->number_of_element(B); }
    /// print box
    std::ostream&Print(const Box*B, std::ostream&out)
    {
      int num = BoxLocalNo(B);
#ifdef OCTALTREE_USE_OPENMP
      depth_type thr=IT;
      if(TRES && num < 0)
	for(thr=0; thr != NT; ++thr) 
	  if(thr!=IT) {
	    num=TRES[thr]->BoxLocalNo(B);
	    if(num >= 0) break;
	  }
#endif// OCTALTREE_USE_OPENMP
      return num<0?
	out<<"Bunknown" :
	out<<'B'<<setfill('0')
#ifdef OCTALTREE_USE_OPENMP
	   <<setw(2)<<int(thr)<<'.'<<setw(4)<<num
#else // OCTALTREE_USE_OPENMP
	   <<setw(7)<<num
#endif// OCTALTREE_USE_OPENMP
	   <<setfill(' ');
    }
    /// dump box data
    /// \param[out] ns  counter for single-child boxes
    void Dump(const Box*B, std::ostream&out, count_type&ns, bool root=false)
    {
      Print(B, out<<' ')
	<<' ' << setfill(' ')<<setw(5)<<B->NUM
	<<' ' << setw(2)<<int(B->NOC)
	<<' ' << setw(2)<<int(B->NBX)
	<<' ' << setw(2)<<int(B->MAP.Rot())
	<<' ' << setw(2)<<int(B->LEV);
      for(int i=0; i!=Nsub; ++i)
	out<<' ' << setw(2)<<int(B->NDL[i]);
      for(int i=0; i!=Nsub; ++i) {
	if(B->NDL[i]==0)
	  out<<" nil     ";
	else if(B->NDL[i] <= NMAX)
	  Print(pDOT(B->OCT[i]),out<<' ');
	else
	  Print(pBOX(B->OCT[i]),out<<' ');
      }
      if(B->NOC==1 && B->NBX==1) {
	++ns;
	out<<" S";
      } else
	out<<" M";
      Print(root? B:NonSingleChildBox(const_cast<Box*>(B)),out<<' ')
	<<' ' << setfill(' ')<<setw(8)<<B->CUB.H
	<<' ' << setw(8)<<B->CUB.X<<'\n';
    }
    /// dump tree, recursive
    void Dump(const Box*B, std::ostream&outd, std::ostream&outb, count_type&ns,
	      bool root=false)
    {
      Dump(B,outb,ns,root);
      for(int i=0; i!=Nsub; ++i)
	if     (B->NDL[i] >  NMAX) Dump(pBOX(B->OCT[i]),outd,outb,ns);
	else if(B->NDL[i] != 0   ) Dump(pDOT(B->OCT[i]),outd);
    }
    /// dump tree
    /// \note to be used after construction
    void DumpNodes() {
      count_type ns=0;
      output dumpD, dumpB;
      dumpD.reopen("fdot%d.dat",IT);
      dumpB.reopen("fbox%d.dat",IT);
      DumpHeadDot(dumpD);
      DumpHeadBox(dumpB);
      Dump(ROOT,dumpD,dumpB,ns,true);
      DebugInfoN(5," %d single-child boxes\n",ns);
    }
#endif// TESTING
  };// class TreeImplementer<>
} // namespace {
//
namespace WDutils {
#define IMPL(DIM,TYPE)							\
  template<>								\
  struct traits< ::TreeImplementer<OctalTree<DIM,TYPE> > >		\
  { static const char*name() {						\
      return "TreeImplementer<OctalTree<" __STRING(DIM) ","		\
	__STRING(TYPE) ">>";						\
  } };									\
  template<>								\
  struct traits<typename ::TreeImplementer<OctalTree<DIM,TYPE> >::Box >	\
  { static const char*name() {						\
      return "TreeImplementer<OctalTree<" __STRING(DIM) ","		\
	__STRING(TYPE) ">>::Box";					\
  } };
  IMPL(2,float);
  IMPL(3,float);
  IMPL(2,double);
  IMPL(3,double);
#undef IMPL
#ifdef OCTALTREE_USE_OPENMP
#define IMPL(DIM,TYPE)							\
  template<>								\
  struct traits<typename ::TreeImplementer<OctalTree<DIM,TYPE>>::TopBox > \
  { static const char*name() {						\
      return "TreeImplementer<OctalTree<" __STRING(DIM) ","		\
	__STRING(TYPE) ">>::TopBox";					\
  } };
  IMPL(2,float);
  IMPL(3,float);
  IMPL(2,double);
  IMPL(3,double);
#undef IMPL
#endif
}
//
namespace {
  //
  // set dots for local partition from old domain
  //
  template<int __D, typename __X>
  void TreeImplementer<OctalTree<__D,__X> >::
  SetDots(const Initialiser*init, const OTree*tree,
	  point&Xmin, point&Xmax, point&Xsum)
  {
    // 1 allocate & initialise dots
    DOTS = init->InitIntern(NDOT,
#ifdef OCTALTREE_USE_OPENMP
			    IT,NT,
#endif// OCTALTREE_USE_OPENMP
			    tree? tree->template particle_keys<0>() : 0,
			    tree? tree->NLEAF : 0);
    if(DOTS==0) WDutils_ErrorN("0 returned by Initialiser::InitIntern()\n");
    if(NDOT==0) WDutils_ErrorN("# Dots = 0 from Initialiser::InitIntern()\n");
    // 2 cumulate Xmin,Xmax,Xsum
    Dot*Di=DOTS, *const DN=DOTS+NDOT;
    Xmin = Xmax = Xsum = Di->X;
    for(++Di; Di!=DN; ++Di) {
      Di->X.up_min_max(Xmin,Xmax);
      Xsum += Di->X;
    }
    // 3 check for nan and inf
    if(isnan(Xsum) || isinf(Xsum)) {
      for(Di=DOTS; Di!=DN; ++Di)
	if(isnan(Di->X) || isinf(Di->X)) {
	  std::ostringstream out;
	  out<<"OctalTree: x="<<Di->X<<" for particle "<<Di->I;
	  WDutils_ErrorN("%s\n",out.str().c_str());
	}
      WDutils_ErrorN("OctalTree: found non-numeric position sum\n");
    }
  }
  //
  // build box-dot tree from all 'our' dots
  //
  template<int __D, typename __X>
  void TreeImplementer< OctalTree<__D,__X> >::
  Build(const Initialiser*init, const OTree*tree)
  {
    // 1 load data into dots
    point Xmin,Xmax,Xave;
    SetDots(init,tree,Xmin,Xmax,Xave);
    // 2 establish global domain
    unsigned Ntot = NDOT;
#ifdef OCTALTREE_USE_OPENMP
    if(NT>1) {
      OMP::AllReduce<Parallel::Sum>(NDOT,Ntot);
      OMP::AllReduce<Parallel::Min>(static_cast<pos_type*>(Xmin),Dim);
      OMP::AllReduce<Parallel::Max>(static_cast<pos_type*>(Xmax),Dim);
      OMP::AllReduce<Parallel::Sum>(static_cast<pos_type*>(Xave),Dim);
    }
#endif// OCTALTREE_USE_OPENMP
    Xave /= pos_type(Ntot);
    // 3 set empty box: local root
    BXAL = new block_alloc<Box>(NBoxes<Dim>(NDOT,NMAX));
    ROOT = BXAL->new_element();
    ROOT->NDl   = 0;
    ROOT->NOC   = 0;
    ROOT->NBX   = 0;
    ROOT->CUB.X = Helper::RootCentre(Xave);
    ROOT->CUB.H = Helper::RootRadius(ROOT->CUB.X,Xmin,Xmax);
    ROOT->LEV   = 0;
    ROOT->NUM   = 0;
    ROOT->MAP   = peano_map::Root();
    if(debug(3)) {
      double time=WALLCLOCK();
      DebugInfoN("OctalTree: dots initialised in   %f sec\n", time-TIME);
      TIME = time;
    }
    // 4 add all local dots to local root box
    for(Dot*Di=DOTS,*DN=DOTS+NDOT; Di!=DN; ++Di)
      AddDotToBox(ROOT,Di);
    if(debug(3)) {
      double time=WALLCLOCK();
      DebugInfoN("OctalTree: box-dot tree grown in %f sec (NDOT=%d)\n",
		time-TIME,NDOT);
      TIME = time;
    }
  }
  //
  // count boxes (including root) and dots in the tree
  //
  template<int __D, typename __X>
  void TreeImplementer< OctalTree<__D,__X> >::
  CountNodesSerial(count_type&NB, count_type&ND) const
  {
    NB=1, ND=NDOT;
    Stack<const Box*> BS(Nsub*MaximumDepth);
    BS.push(ROOT);
    while(! BS.is_empty() ) {
      const Box*P=BS.pop();
      if(ASCC) EnsureNonSingleChildBox(P);
      for(int i=0; i!=Nsub; ++i)
	if(P->NDL[i] >= NMIN) {
	  ++NB;
	  if(P->NDL[i] > NMAX)
	    BS.push(cpBOX(P->OCT[i]));
	} 
    }
  }
  //
  // link a cell (recursive)
  //
  template<int __D, typename __X>
  depth_type TreeImplementer< OctalTree<__D,__X> >::
  LinkCell(const Box*P, const count_type C, const count_type cN) const
  {
#ifdef TESTING
    if(C >= cN)
      WDutils_ErrorN("error in cell linking: C=%d >= cN=%d\n",C,cN);
#endif
    // 1 if single-child box replace by non-single-child box descendant
    if(ASCC) EnsureNonSingleChildBox(P);
    // 2 copy some data, set octant
    TREE->_LE[C] = P->LEV;
    GeoAlg::copy(P->CUB,static_cast<cube&>(TREE->_XC[C]));
    TREE->_L0[C] = LF;
    TREE->_NM[C] = P->NUM;
#ifdef OCTALTREE_HAVE_D0
    TREE->_D0[C] = IT;
#endif
    // 3 loop octants: link leaf kids (from octants with < NMIN), count cells
    depth_type tmp = 0;
    for(peano_dig d=0; d!=Nsub; ++d) {
      peano_oct i = P->MAP.octant(d);
      if(P->NDL[i] >= NMIN)
	++tmp;
      else if(P->NDL[i])
	for(const Dot*Di=cpDOT(P->OCT[i]); Di; Di=pDOT(Di->Next))
	  LinkLeaf(LF++,Di,C);
    }
    // 4 set number of leaves kids and daughter cells
    TREE->_NL[C] = LF - TREE->_L0[C];
    TREE->_NC[C] = tmp;
    // 5 link daughter cells, if any
    if(tmp) {
      TREE->_C0[C] = CF;
      count_type Ci = CF;
      CF += tmp;
      tmp = 0;
      for(peano_dig d=0; d!=Nsub; ++d) {
	peano_oct oc = P->MAP.octant(d);
	if(P->NDL[oc] > NMAX) {
	  TREE->_PA[Ci] = C;
	  TREE->_OC[Ci] = oc;
	  update_max(tmp, LinkCell(cpBOX(P->OCT[oc]),Ci++,cN));
	} else if(P->NDL[oc] >= NMIN) {
	  TREE->_PA[Ci] = C;
	  TREE->_OC[Ci] = oc;
	  LinkOctantCell(P,Ci++,oc,cN);
	}
      }
      ++tmp;
    } else 
      TREE->_C0[C] = C;
    // 6 set and return tree depth of cell
    return TREE->_DP[C] = tmp;
  }
  //
  // link serial
  //
  template<int __D, typename __X>
  void TreeImplementer< OctalTree<__D,__X> >::
  LinkSerial(OTree*tree) const
  {
    const_cast<OTree*&>(TREE) = tree;
    CF = 1;
    LF = 0;
#ifdef OCTALTREE_USE_OPENMP
    DomainData*DOM=const_cast<DomainData*>(tree->DOM);
    tree->NTOPC = 1;
    DOM->BRANCH.resize(1);
    DOM->BRANCH[0] = 0;
    DOM->C0= CF;
    DOM->L0= LF;
#endif// OCTALTREE_USE_OPENMP
    TREE->_PA[0] = 0;
    TREE->_OC[0] = 0;
    LinkCell(ROOT,0,TREE->NCELL);
    tree->NCELL = CF;
#ifdef OCTALTREE_USE_OPENMP
    DOM->CN= CF; DOM->NC= CF-DOM->C0;
    DOM->LN= LF; DOM->NL= LF-DOM->L0;
    DOM->DEPTH= TREE->_DP[0];
#endif// OCTALTREE_USE_OPENMP
    if(debug(3)) {
      double time=WALLCLOCK();
      DebugInfoN("OctalTree: tree linked in       %f sec\n", time-TIME);
      const_cast<double&>(TIME) = time;
    }
  }
  //
  // destructor
  //
  template<int __D, typename __X>
  TreeImplementer< OctalTree<__D,__X> >::~TreeImplementer()
  {
    if(DOTS) WDutils_DEL16(DOTS); DOTS=0;
    if(BXAL) WDutils_DEL_O(BXAL); BXAL=0;
#ifdef OCTALTREE_USE_OPENMP
    if(TBAL) WDutils_DEL_O(TBAL); TBAL=0;
#endif// OCTALTREE_USE_OPENMP
    if(debug(3)) {
      double time=WALLCLOCK();
      DebugInfoN("OctalTree: finished after        %f sec\n", time-TIME);
    }
  }
  //
  // allocate OctalTree data
  //
  template<int __D, typename __X>
  void TreeImplementer< OctalTree<__D,__X> >::
  AllocateTree(OTree*const tree, const count_type nc, const count_type nl,
	       const count_type ne, const bool flag_and_rung) WDutils_THROWING
  {
    tree->NLEAFEXT = ne;
    tree->NLEAF    = nl;
    tree->NCELL    = nc;
#ifdef OCTALTREE_DATA_IN_ONE_BLOCK
    size_t need = Next16<point16>      (ne)        // _XL[1]
      +           Next16<particle_key> (ne)        // _PL[1]
      +           Next16<point16>      (nl)        // _XL[0]
      +           Next16<particle_key> (nl)        // _PL[0]
      +           Next16<count_type>   (nl)        // _PC
      +           Next16<depth_type>   (nc)        // _LE
      +           Next16<octant_type>  (nc)        // _OC
      +           Next16<cube16>       (nc)        // _XC
      +           Next16<count_type>   (nc)        // _L0
      +           Next16<local_count>  (nc)        // _NL
      +           Next16<count_type>   (nc)        // _NM
      +           Next16<count_type>   (nc)        // _C0
      +           Next16<octant_type>  (nc)        // _NC
      +           Next16<count_type>   (nc)        // _PA
      +           Next16<depth_type>   (nc)        // _DP
# ifdef OCTALTREE_HAVE_D0
      +           Next16<domain_id>    (nc)        // _D0
# endif// OCTALTREE_USE_OPENMP
      ;
    if(flag_and_rung)
      need     += Next16<float>        (ne)        // _RG[1]
	+         Next16<uint8_t>      (ne)        // _FL[1]
	+         Next16<float>        (nl)        // _RG[0]
	+         Next16<uint8_t>      (nl)        // _FL[0]
	+         Next16<count_type>   (nc);       // _NA
    tree->_MEM. template reset_conditional<3,2>(need);
# define SET_DATA(Name,Number,Type)				\
    if(Number) {						\
      tree->Name = reinterpret_cast<Type*>(A);			\
      A += Next16<Type>(Number);				\
    } else							\
      tree->Name = 0;						\
    DebugInfoN(6,"OctalTree: set %s = %p for %d %s\n",		\
	       __STRING(Name),tree->Name,Number,__STRING(Type));
    char*A = tree->_MEM.begin();
    SET_DATA(_XL[1],ne,point16);
    SET_DATA(_PL[1],ne,particle_key);
    SET_DATA(_RG[1],flag_and_rung?ne:0,float);
    SET_DATA(_FL[1],flag_and_rung?ne:0,uint8_t);
    SET_DATA(_XL[0],nl,point16);
    SET_DATA(_PL[0],nl,particle_key);
    SET_DATA(_PC   ,nl,count_type);
    SET_DATA(_RG[0],flag_and_rung?nl:0,float);
    SET_DATA(_FL[0],flag_and_rung?nl:0,uint8_t);
    SET_DATA(_LE,nc,depth_type);
    SET_DATA(_OC,nc,octant_type);
    SET_DATA(_XC,nc,cube16);
    SET_DATA(_L0,nc,count_type);
    SET_DATA(_NL,nc,local_count);
    SET_DATA(_NM,nc,count_type);
    SET_DATA(_C0,nc,count_type);
    SET_DATA(_NC,nc,octant_type);
    SET_DATA(_PA,nc,count_type);
    SET_DATA(_NA,flag_and_rung?nc:0,count_type);
    SET_DATA(_DP,nc,depth_type);
# ifdef OCTALTREE_HAVE_D0
    SET_DATA(_D0,nc,domain_id);
# endif// OCTALTREE_USE_OPENMP
    WDutilsAssert(A == tree->_MEM.begin()+need && A <= tree->_MEM.end());
#else // OCTALTREE_DATA_IN_ONE_BLOCK
    const count_type nc4=(nc+3)&~3;              // multiple of 4 >= nc
    const count_type ne4=(ne+3)&~3;              // multiple of 4 >= ne
    const count_type nl4=(nl+3)&~3;              // multiple of 4 >= nl
# define SET_DATA(Name,Number)					\
    if(Number) {						\
      tree->Name.resize(Number);				\
    } else {							\
      tree->Name.clear();					\
      tree->Name.shrink_to_fit();				\
    }								\
    DebugInfoN(6,"OctalTree: set %s = %p for %d %s\n",		\
	       __STRING(Name),tree->Name.data(),Number,		\
	       nameof(typename std::remove_reference<		\
		      decltype(tree->Name[0])>::type));
    SET_DATA(_XL[1],ne4);
    SET_DATA(_PL[1],ne);
    SET_DATA(_RG[1],flag_and_rung?ne4:0);
    SET_DATA(_FL[1],flag_and_rung?ne4:0);
    SET_DATA(_XL[0],nl4);
    SET_DATA(_PL[0],nl);
    SET_DATA(_RG[0],flag_and_rung?nl4:0);
    SET_DATA(_FL[0],flag_and_rung?nl4:0);
    SET_DATA(_PC,nl);
    SET_DATA(_LE,nc);
    SET_DATA(_OC,nc);
    SET_DATA(_XC,nc4);
    SET_DATA(_L0,nc);
    SET_DATA(_NL,nc);
    SET_DATA(_NM,nc);
    SET_DATA(_NA,flag_and_rung?nc:0);
    SET_DATA(_C0,nc);
    SET_DATA(_NC,nc);
    SET_DATA(_PA,nc);
    SET_DATA(_DP,nc);
# ifdef OCTALTREE_HAVE_D0
    SET_DATA(_D0,nc);
# endif// OCTALTREE_USE_OPENMP
#endif
#undef SET_DATA
  }
  //
  // BuildTreeSerial
  // 
  template<int __D, typename __X>
  void TreeImplementer< OctalTree<__D,__X> >::
  BuildTreeSerial(OTree*tree,bool fresh, const BuilderInterface*builder=0)
  {
    double    wtime = WALLCLOCK();
    count_type Next = tree->INIT->Nexternal();
    bool load_flags = tree->INIT->LoadRungFlag();
    {
      // 1 build box-dot-tree
      TreeImplementer TI(wtime, tree->INIT, 
#ifdef OCTALTREE_USE_OPENMP
			 0, 1, 
#endif
			 fresh?0:tree, tree->NMAX, tree->NMIN, tree->ASCC);
#ifdef TESTING
      if(debug(4)) TI.DumpNodes();
#endif
      // 2 allocate top domain and link cells and leaves
      TI.CountNodesSerial(tree->NCELL,tree->NLEAF);
      AllocateTree(tree,tree->NCELL,tree->NLEAF,Next,load_flags);
      TI.LinkSerial(tree);
      // 3 load rungs & flags and pass _NA[] up the tree
      if(load_flags) {
	bool all_active = tree->INIT->
	  InitRungFlag(tree->template particle_keys<0>(),
		       p_data(tree->_RG[0]),
		       p_data(tree->_FL[0]),
		       tree->NLEAF);
	PassUpNactive(tree,all_active,0,tree->NCELL);
      }
      // 4 destruct TreeImplementer
    }
    // 5 load external particles
    if(Next) {
      for(count_type e=0; e!=Next; ++e)
	tree->INIT->InitExtern(e, static_cast<point&>(tree->_XL[1][e]),
			       tree->_PL[1][e]);
      if(load_flags)
	tree->INIT->
	  InitRungFlag(tree->template particle_keys<1>(),
		       p_data(tree->_RG[1]),
		       p_data(tree->_FL[1]),
		       Next);
    }
    // 6 set derived data
    if(builder)
      builder->Serial(fresh,tree);
    DebugInfoN(2,"  OctalTree: build finished\n");
  }
  //
#ifdef OCTALTREE_USE_OPENMP
# define OCTALTREE_IMPLEMENTER_ADDITIONAL_MEMBER_FUNC_DEFS
# include <../devel/octtree_omp.cc>
# undef  OCTALTREE_IMPLEMENTER_ADDITIONAL_MEMBER_FUNC_DEFS
#endif
} // namespace {
//
namespace WDutils {
  //
  // OctalTree::~OctalTree()
  //
  template<int __D, typename __X>
  OctalTree<__D,__X>::~OctalTree()
  {
#ifdef OCTALTREE_USE_OPENMP
    WDutilsAssertE(!OMP::IsParallel());
    WDutils_DEL_AN(DOM,NDOM);
#endif // OCTALTREE_USE_OPENMP
  }  
  //
  // OctalTree::OctalTree()
  //
  template<int __D, typename __X>
  OctalTree<__D,__X>::OctalTree(const Initialiser*ini, Data const&data)
    throw(exception)
    : INIT  ( ini )
    , PERB  ( INIT->Boundary() )
#ifdef OCTALTREE_USE_OPENMP
    , NDOM  ( RunInfo::omp_threads() )
    , DOM   ( WDutils_NEW(DomainData,NDOM) )
#endif// OCTALTREE_USE_OPENMP
    , NBUILD( 0 )
    , NMAX  ( data.NMAX<1? 1 : (data.NMAX>MAXNMAX? MAXNMAX:data.NMAX) )
    , NMIN  ( data.NMIN? (data.NMIN<NMAX? data.NMIN:NMAX) : (NMAX<2?NMAX:2) )
    , ASCC  ( data.ASCC )
  {
    // sanity checks
    if(0==INIT)
      WDutils_THROWF("null pter to Initialiser provided\n");
    if(data.NMAX == 0)
      WDutils_WarningF("nmax=0; will use nmax=%d instead\n",int(NMAX));
    if(data.NMAX > MAXNMAX)
      WDutils_WarningF("nmax=%d exceeds %d; will use nmax=%d instead\n",
		       int(data.NMAX),int(MAXNMAX),int(NMAX));
    if(data.NMIN > MAXNMAX)
      WDutils_WarningF("nmin=%d exceeds %d; will use nmin=%d instead\n",
		       int(data.NMIN),int(MAXNMAX),int(NMIN));
    else if(data.NMIN > data.NMAX)
      WDutils_WarningF("nmin=%d exceeds nmax=%d; will use nmin=%d\n",
		       int(data.NMIN),int(data.NMAX),int(NMIN));
#ifdef OCTALTREE_HAVE_DATA
    // reset DATA[]
    for(int i=0; i!=NDAT; ++i) DATA[0][i] = DATA[1][i] = 0;
#endif
    // build tree
    TreeImplementer<OctalTree>::BuildTree(this,1,data.TOL);
    ++NBUILD;
  }
  //
  // OctalTree::OctalTree()
  //
  template<int __D, typename __X>
  OctalTree<__D,__X>::OctalTree(const Initialiser*ini, Data const&data,
				BuilderInterface const&builder)
    throw(exception)
    : INIT  ( ini )
    , PERB  ( INIT->Boundary()  )
#ifdef OCTALTREE_USE_OPENMP
    , NDOM  ( RunInfo::omp_threads() )
    , DOM   ( WDutils_NEW(DomainData,NDOM) )
#endif// OCTALTREE_USE_OPENMP
    , NBUILD( 0 )
    , NMAX  ( data.NMAX<1? 1 : (data.NMAX>MAXNMAX? MAXNMAX:data.NMAX) )
    , NMIN  ( data.NMIN? (data.NMIN<NMAX? data.NMIN:NMAX) : (NMAX<2?NMAX:2) )
    , ASCC  ( data.ASCC )
  {
    // sanity checks
    if(0==INIT)
      WDutils_THROWF("null pter to Initialiser provided\n");
    if(data.NMAX == 0)
      WDutils_WarningF("nmax=0; will use nmax=%d instead\n",int(NMAX));
    if(data.NMAX > MAXNMAX)
      WDutils_WarningF("nmax=%d exceeds %d; will use nmax=%d instead\n",
		       int(data.NMAX),int(MAXNMAX),int(NMAX));
    if(data.NMIN > MAXNMAX)
      WDutils_WarningF("nmin=%d exceeds %d; will use nmin=%d instead\n",
		       int(data.NMIN),int(MAXNMAX),int(NMIN));
    else if(data.NMIN > data.NMAX)
      WDutils_WarningF("nmin=%d exceeds nmax=%d; will use nmin=%d\n",
		       int(data.NMIN),int(data.NMAX),int(NMIN));
#ifdef OCTALTREE_HAVE_DATA
    // reset DATA[]
    for(int i=0; i!=NDAT; ++i) DATA[0][i] = DATA[1][i] = 0;
#endif
    // build tree
    TreeImplementer<OctalTree>::BuildTree(this,1,data.TOL,&builder);
    ++NBUILD;
  }
  //
  // OctalTree::rebuild
  //
  template<int __D, typename __X>
  void OctalTree<__D,__X>::rebuild(Data const&data) throw(exception)
  {
    // sanity checks
    if(data.NMAX!=0) {
      NMAX = data.NMAX>MAXNMAX? MAXNMAX:data.NMAX;
      NMIN =(data.NMIN? (data.NMIN<NMAX?data.NMIN:NMAX) : (NMAX<2?NMAX:2) );
      if(data.NMAX > MAXNMAX)
	WDutils_WarningF("nmax=%d exceeds %d; will use nmax=%d instead\n",
			 int(data.NMAX),int(MAXNMAX),int(NMAX));
      if(data.NMIN > MAXNMAX)
	WDutils_WarningF("nmin=%d exceeds %d; will use nmin=%d instead\n",
			 int(data.NMIN),int(MAXNMAX),int(NMIN));
      else if(data.NMIN > data.NMAX)
	WDutils_WarningF("nmin=%d exceeds nmax=%d; will use nmin=%d\n",
			 int(data.NMIN),int(data.NMAX),int(NMIN));
    }
#ifdef OCTALTREE_HAVE_DATA
    // reset DATA[]
    for(int i=0; i!=NDAT; ++i) DATA[0][i] = DATA[1][i] = 0;
#endif
    // rebuild tree
    TreeImplementer<OctalTree>::BuildTree(this,0,data.TOL);
    ++NBUILD;
  }
  //
  // OctalTree::rebuild
  //
  template<int __D, typename __X>
  void OctalTree<__D,__X>::rebuild(Data const&data,
				   BuilderInterface const&builder)
    throw(exception)
  {
    // sanity checks
    if(data.NMAX!=0) {
      NMAX = data.NMAX>MAXNMAX? MAXNMAX:data.NMAX;
      NMIN =(data.NMIN? (data.NMIN<NMAX?data.NMIN:NMAX) : (NMAX<2?NMAX:2) );
      if(data.NMAX > MAXNMAX)
	WDutils_WarningF("nmax=%d exceeds %d; will use nmax=%d instead\n",
			 int(data.NMAX),int(MAXNMAX),int(NMAX));
      if(data.NMIN > MAXNMAX)
	WDutils_WarningF("nmin=%d exceeds %d; will use nmin=%d instead\n",
			 int(data.NMIN),int(MAXNMAX),int(NMIN));
      else if(data.NMIN > data.NMAX)
	WDutils_WarningF("nmin=%d exceeds nmax=%d; will use nmin=%d\n",
			 int(data.NMIN),int(data.NMAX),int(NMIN));
    }
#ifdef OCTALTREE_HAVE_DATA
    // reset DATA[]
    for(int i=0; i!=NDAT; ++i) DATA[0][i] = DATA[1][i] = 0;
#endif
    // build tree and derived
    TreeImplementer<OctalTree>::BuildTree(this,0,data.TOL,&builder);
    ++NBUILD;
  }
  //
  // OctalTree::scc()
  //
#ifdef __SSE__
#  define cPF(__X) static_cast<const float*>(__X)
  //
  // SSE version for 2D, single precision
  template<>
  count_type OctalTree<2,float>::scc(point const&x) const noexcept
  {
    count_type c=0;
    __m128 X = _mm_loadu_ps(cPF(x));
    __m128 C = _mm_load_ps(cPF(_XC[c].X));
    __m128 H = _mm_shuffle_ps(C,C,_MM_SHUFFLE(2,2,2,2));
#define NotInCell							\
    3 != (3&_mm_movemask_ps(_mm_and_ps(_mm_cmple_ps(_mm_sub_ps(C,H),X),	\
				       _mm_cmpgt_ps(_mm_add_ps(C,H),X))))
    if(NotInCell) return c;
    while(_NC[c]) {
      uint8_t o = 3&_mm_movemask_ps(_mm_cmplt_ps(C,X));
      count_type cc=_C0[c], ce=cc+_NC[c];
      while(cc!=ce && o!=_OC[cc]) ++cc;
      if(cc==ce)
	return c;
      C = _mm_load_ps(cPF(_XC[cc].X));
      if(ASCC && _LE[cc] > _LE[c]+1 && NotInCell)
	return c;
      c = cc;
    }
    return c;
#undef NotInCell
  }
  // SSE version for 3D, single precision
  template<>
  count_type OctalTree<3,float>::scc(point const&x) const noexcept
  {
    count_type c=0;
    __m128 X = _mm_loadu_ps(cPF(x));
    __m128 C = _mm_load_ps(cPF(_XC[c].X));
    __m128 H = _mm_shuffle_ps(C,C,_MM_SHUFFLE(3,3,3,3));
#define NotInCell							\
    7 != (7&_mm_movemask_ps(_mm_and_ps(_mm_cmple_ps(_mm_sub_ps(C,H),X),	\
				       _mm_cmpgt_ps(_mm_add_ps(C,H),X))))
    if(NotInCell) return c;
    while(_NC[c]) {
      uint8_t o = 7&_mm_movemask_ps(_mm_cmplt_ps(C,X));
      count_type cc=_C0[c], ce=cc+_NC[c];
      while(cc!=ce && o!=_OC[cc]) ++cc;
      if(cc==ce)
	return c;
      C = _mm_load_ps(cPF(_XC[cc].X));
      if(ASCC && _LE[cc] > _LE[c]+1 && NotInCell)
	return c;
      c = cc;
    }
    return c;
#undef NotInCell
  }
#ifdef __SSE2__
#  define cPD0(__X) static_cast<const double*>(__X)
#  define cPD1(__X) static_cast<const double*>(__X)+2
  // SSE version for 2D, double precision
  template<>
  count_type OctalTree<2,double>::scc(point const&x) const noexcept
  {
    count_type c=0;
    __m128d X = _mm_loadu_pd(cPD0(x));
    __m128d C = _mm_load_pd(cPD0(_XC[c].X));
    __m128d H = _mm_set1_pd(_XC[c].H);
#define NotInCell							\
    3 != _mm_movemask_pd(_mm_and_pd(_mm_cmple_pd(_mm_sub_pd(C,H),X),	\
				    _mm_cmpgt_pd(_mm_add_pd(C,H),X)))
    if(NotInCell) return c;
    while(_NC[c]) {
      uint8_t o = _mm_movemask_pd(_mm_cmplt_pd(C,X));
      count_type cc=_C0[c], ce=cc+_NC[c];
      while(cc!=ce && o!=_OC[cc]) ++cc;
      if(cc==ce)
	return c;
      C = _mm_load_pd(cPD0(_XC[cc].X));
      if(ASCC && _LE[cc] > _LE[c]+1 && NotInCell)
	return c;
      c = cc;
    }
    return c;
#undef NotInCell
  }
  // SSE version for 3D, double precision
  template<>
  count_type OctalTree<3,double>::scc(point const&x) const noexcept
  {
    count_type c=0;
    __m128d Xa = _mm_loadu_pd(cPD0(x));
    __m128d Ca = _mm_load_pd(cPD0(_XC[c].X));
    __m128d H  = _mm_set1_pd(_XC[c].H);
#define NotInCell_a							\
    3 != _mm_movemask_pd(_mm_and_pd(_mm_cmple_pd(_mm_sub_pd(Ca,H),Xa),	\
				    _mm_cmpgt_pd(_mm_add_pd(Ca,H),Xa)))
    if(NotInCell_a)
      return c;
    __m128d Xb = _mm_loadu_pd(cPD1(x));
    __m128d Cb = _mm_load_pd(cPD1(_XC[c].X));
#define NotInCell_b							\
    ! (1&_mm_movemask_pd(_mm_and_pd(_mm_cmple_pd(_mm_sub_pd(Cb,H),Xb),	\
				    _mm_cmpgt_pd(_mm_add_pd(Cb,H),Xb))))
    if(NotInCell_b)
      return c;
    while(_NC[c]) {
      uint8_t o = _mm_movemask_pd(_mm_cmplt_pd(Ca,Xa)) |
	(      (1&_mm_movemask_pd(_mm_cmplt_pd(Cb,Xb)))<<2) ;
      count_type cc=_C0[c], ce=cc+_NC[c];
      while(cc!=ce && o!=_OC[cc]) ++cc;
      if(cc==ce)
	return c;
      Ca = _mm_load_pd(cPD0(_XC[cc].X));
      Cb = _mm_load_pd(cPD1(_XC[cc].X));
      if(ASCC && _LE[cc] > _LE[c]+1 && (NotInCell_a || NotInCell_b))
	return c;
      c = cc;
    }
    return c;
#undef NotInCell_a
#undef NotInCell_b
  }
#endif // __SSE2__
#endif // __SSE__
  // generic non-SSE version
#ifndef __SSE2__
#ifdef __SSE__
  //   only for double
  template<int __D>
  count_type OctalTree<__D,double>::scc(point const&x) const noexcept
#else // __SSE__
  //   for float and double
  template<int __D, typename __X>
  count_type OctalTree<__D,__X>::scc(point const&x) const noexcept
#endif// __SSE__
  {
    count_type c=0;
    if(! Geometry::Algorithms<0>::contains(_XC[c],x))
      return c;
    while(_NC[c]) {
      uint8_t o=Geometry::Algorithms<0>::octant(_XC[c],x);
      count_type cc=_C0[c], ce=cc+_NC[c];
      while(cc!=ce && o!=_OC[cc]) ++cc;
      if(cc==ce || ( ASCC && _LE[cc] > _LE[c]+1 &&
		     ! Geometry::Algorithms<0>::contains(_XC[cc],x)))
	return c;
      c = cc;
    }
    return c;
  }
#endif// no __SSE2__
  //
  // OctalTree::LeafHeader()
  //
  template<int __Dim, typename __X> template<int __D>
  typename std::enable_if<__D==2,std::ostream&>::type
  OctalTree<__Dim,__X>::__LeafHeader(std::ostream&s) const
  {
    s  << "Leaf     ";
    if(have_flag())
      s<< " f rg";
    s  << "      k up       ";
    s  << "             X            ";
    return s;
  }
  //
  template<int __Dim, typename __X> template<int __D>
  typename std::enable_if<__D==3,std::ostream&>::type
  OctalTree<__Dim,__X>::__LeafHeader(std::ostream&s) const
  {
    s  << "Leaf     ";
    if(have_flag())
      s<< " f rg";
    s  << "       k up       ";
    s  << "                         X             ";
    return s;
  }
  //
  // OctalTree::Dump(Leaf)
  //
  template<int __Dim, typename __X>
  std::ostream&OctalTree<__Dim,__X>::Dump(const Leaf l, std::ostream&s)
    const
  {
    s  << l;
    if(have_flag())
      s<< ' ' << std::setw(1) << (flag(l)!=0)
       << ' ' << std::setw(2) << int(rung(l));
    s  << ' ' << std::setw(7) << _PL[0][l.I]
       << ' ' << std::setfill(' ') << Cell(_PC[l.I])
       << ' ' << std::setw(12) << _XL[0][l.I];
    return s;
  }
  //
  // OctalTree::Dump(ExtLeaf)
  //
  template<int __Dim, typename __X>
  std::ostream&OctalTree<__Dim,__X>::Dump(const ExtLeaf l, std::ostream&s)
    const
  {
    s  << l << ' ';
    if(have_flag())
      s<< ' ' << std::setw(1) << (flag(l)!=0)
       << ' ' << std::setw(2) << int(rung(l));
    s  << std::setw(7) << _PL[1][l.I] << ' '
       << std::setfill(' ') 
       << "nil        "
       << std::setw(12) << _XL[1][l.I];
    return s;
  }
  //
  // OctalTree::CellHeader()
  //
  template<int __Dim, typename __X> template<int __D>
  typename std::enable_if<__D==2,std::ostream&>::type
  OctalTree<__Dim,__X>::__CellHeader(std::ostream&s) const
  {
    s  << "Cell     ";
    s  << " le dp";
#ifdef OCTALTREE_HAVE_D0
    s  << " d0";
#endif
    s  << " up       ";
    s  << "oc C0       ";
    s  << "Nc L0        ";
    s  << " Nl ";
    s  << "      N ";
    if(have_flag())
      s<< "     Na ";
    s  << "         R ";
    s  << "          X             ";
    return s;
  }
  //
  template<int __Dim, typename __X> template<int __D>
  typename std::enable_if<__D==3,std::ostream&>::type
  OctalTree<__Dim,__X>::__CellHeader(std::ostream&s) const
  {
    s  << "Cell     ";
    s  << " le dp";
#ifdef OCTALTREE_HAVE_D0
    s  << " d0";
#endif
    s  << " up       ";
    s  << "oc C0       ";
    s  << "Nc L0        ";
    s  << " Nl ";
    s  << "      N ";
    if(have_flag())
      s<< "     Na ";
    s  << "         R ";
    s  << "                      X             ";
    return s;
  }
  //
  // OctalTree::Dump(Cell)
  //
  template<int __Dim, typename __X>
  std::ostream&OctalTree<__Dim,__X>::Dump(const Cell c, std::ostream&s)
    const
  {
    // cell id
    s   << c << ' ';
    // level
    s   << std::setw(2) << int (_LE[c.I]) << ' '
	<< std::setw(2) << int (_DP[c.I]) << ' ';
#ifdef OCTALTREE_HAVE_D0
    // domain and # domains
    s   << std::setw(2) << int (_D0[c.I]) << ' ';
#endif
    // parent
    if(c.I)
      s << Cell(_PA[c.I]) << ' ';
    else 
      s << "nil       ";
    // octant
    s   << int (_OC[c.I]) << ' ';
    // cells
    if(_NC[c.I])
      s << Cell(_C0[c.I]) << ' '
	<< int (_NC[c.I]) << ' ';
    else
      s << "nil       0 ";
    // leaves
    s   << Leaf(_L0[c.I]) << ' '
	<< std::setw(3)  << int(_NL[c.I])  << ' '
	<< std::setw(7)  << _NM[c.I] << ' ';
    if(have_flag())
      s << std::setw(7)  << _NA[c.I] << ' ';
    // radius, centre
    s   << std::setw(10) << _XC[c.I].H << ' '
	<< std::setw(11) << _XC[c.I].X;
    return s;
  }
#ifdef OCTALTREE_USE_OPENMP
# define OCTALTREE_MEMBER_FUNC_DEFS
# include <../devel/octtree_omp.cc>
# undef  OCTALTREE_MEMBER_FUNC_DEFS
#endif
  //
  template class OctalTree<2,float>;
  template class OctalTree<3,float>;
  template class OctalTree<2,double>;
  template class OctalTree<3,double>;
  template std::ostream&OctalTree<2,double>::__CellHeader<2>(std::ostream&)
    const;
  template std::ostream&OctalTree<3,double>::__CellHeader<3>(std::ostream&)
    const;
  template std::ostream&OctalTree<2,float>::__CellHeader<2>(std::ostream&)
    const;
  template std::ostream&OctalTree<3,float>::__CellHeader<3>(std::ostream&)
    const;
  //
} // namespace WDutils
////////////////////////////////////////////////////////////////////////////////
namespace {
#ifdef __SSE__

#  define __TemplateDecl						\
  template<int __Dim, typename __PosType, typename __PropType, bool __Anc>
#  define __TemplateList __Dim,__PosType,__PropType,__Anc

#else  // __SSE__

#  define __TemplateDecl					\
  template<int __Dim, typename __PosType, typename __PropType>
#  define __TemplateList __Dim,__PosType,__PropType

#endif
  ///
  /// auxiliary for InteractionTree<>::Builder
  ///
  __TemplateDecl struct TreeImplementer
  < octtree::InteractionTreeData<__TemplateList> >
  {
    using ITree = InteractionTree<__TemplateList>;
    static const int  Dim       = __Dim;
#ifdef __SSE__
    static const bool Anchoring = ITree::Anchoring;
#endif
    using IData        = typename ITree::DataBase;
    using OTree        = typename ITree::OTree;
    using Leaf         = typename ITree::Leaf;
    using LeafRange    = typename ITree::LeafRange;
    using ExtLeaf      = typename ITree::ExtLeaf;
    using prop_type    = typename ITree::prop_type;
    using Initialiser  = typename ITree::PropInitialiser;
#ifdef OCTALTREE_USE_OPENMP
    using cp_domain    = const typename ITree::DomainData*;
#endif

#ifdef __SSE__
    using Anchor       = typename ITree::Anchor;
    using sse_pos_type = typename ITree::sse_pos_type;
    /// set SSE position without anchoring
    template<bool __A>
    static typename std::enable_if<!__A>::type
    __set_sse(const OTree*tree, IData*data,
# ifdef OCTALTREE_USE_OPENMP
	      count_type dom,
# endif
	      float*, float)
    {
      DebugInfoN(2,"InteractionTree: setting sse positions w/o anchoring\n");
      tree->
# ifdef OCTALTREE_USE_OPENMP
	domain(dom)->
# endif
	loop_leaves([=](const Leaf l) noexcept
		   {
		     for(int d=0; d!=Dim; ++d)
		       data->_X16[d][l.I] = tree->position(l)[d];
		   } );
# ifdef OCTALTREE_USE_OPENMP
      if(dom+1 == tree->n_dom())
# endif
	for(auto l=tree->end_leaf(); l.I<tree->n_leaf4(); ++l)
	  for(int d=0; d!=Dim; ++d)
	    data->_X16[d][l.I] = 0;
    }
    /// set SSE positions and anchors
    template<bool __A>
    static typename std::enable_if< __A>::type
    __set_sse(const OTree*tree, IData*data,
# ifdef OCTALTREE_USE_OPENMP
	      count_type dom,
# endif
	      float*QI, float del)
    {
# ifdef OCTALTREE_USE_OPENMP
      cp_domain DOM = tree->domain(dom);
# else
#  define DOM tree
# endif
      DebugInfoN(2,"InteractionTree: setting anchors and sse positions ...\n");
      // 1 up-ward pass: set QI[cell] = min(Q)
      DebugInfoN(4,"InteractionTree: passing qmin up the domain ...\n");
      DOM -> loop_cells_up([&](const Cell c) noexcept {
	  float qmin;
	  tree->loop_leaf_kids(c,
			       [&](const Leaf l0) noexcept
			       { qmin = data->_HSQ[0][l0.I]; },
			       [&](const Leaf ll) noexcept
			       { update_min(qmin,data->_HSQ[0][ll.I]); }
			       );
	  tree->loop_daughters(c,[&](const Cell cc) noexcept
			       { update_min(qmin,QI[cc.I]); }
			       );
	  QI[c.I] = qmin;
	});
      DebugInfoN(5,"InteractionTree: DONE (passing qmin up the domain)\n");
      // 2 count anchors needed
      DebugInfoN(4,"InteractionTree: counting anchors needed ...\n");
      const double fac = square(double(del)/
				std::numeric_limits<float>::epsilon());
      Stack<Cell> CS(OTree::Nsub * DOM->depth());
      count_type ianc=0,ntot=0;
# ifdef OCTALTREE_USE_OPENMP
      DOM->loop_branches([&](const Cell b) {
	  CS.push(b);
# else
	  CS.push(tree->root());
# endif
	  while(!CS.is_empty()) {
	    auto c = CS.pop();
	    WDutilsAssertE(QI[c.I]>0);
	    if(square(tree->radius(c)) < fac * QI[c.I])
	      ++ntot;
	    else
	      tree->loop_daughters(c,[&](Cell cc) noexcept { CS.push(cc); });
	  }
# ifdef OCTALTREE_USE_OPENMP
	});
# endif
      DebugInfoN(6,"InteractionTree: DONE (counting anchors needed: %d)\n",
		 ntot);
# ifdef OCTALTREE_USE_OPENMP
      if(tree->n_dom()>1) {
	OMP::Cumulate<Parallel::Sum>(ianc=ntot,ntot);
	ianc = ntot-ianc;
      }
      // 3 set anchors and SSE positions
      if(dom+1 == tree->n_dom())
# endif
	// 3.1 last thread allocates anchors
	{
	  DebugInfoN(5,"InteractionTree: %d anchors needed in total\n",ntot);
	  data->_ANCH.resize(ntot);
	}
      DebugInfoN(4,"InteractionTree: setting anchors and SSE pos ...\n");
# ifdef OCTALTREE_USE_OPENMP
# pragma omp barrier
      DOM->loop_branches([&](const Cell b) {
	  CS.push(b);
# else
	  CS.push(tree->root());
# endif
	  while(!CS.is_empty()) {
	    Cell c = CS.pop();
	    if(square(tree->radius(c)) < fac * QI[c.I]) {
	      Anchor*An = &(data->_ANCH[ianc++]);
	      An->Z = tree->centre(c);
	      An->B = tree->begin_leaf(c);
	      An->E = tree->end_leaf_desc(c);
	      tree->loop_all_leaves(c,[&](const Leaf l) noexcept {
		  data->_ANC[l.I] = An;
		  for(int d=0; d!=Dim; ++d)
		    data->_X16[d][l.I] = float(tree->position(l)[d] - An->Z[d]);
		} );
	    } else {
	      tree->loop_leaf_kids(c,[&](const Leaf l) noexcept {
		  data->_ANC[l.I] = 0;
		  for(int d=0; d!=Dim; ++d)
		    data->_X16[d][l.I] = float(tree->position(l)[d] -
					      tree->centre  (c)[d] );
		} );
	      tree->loop_daughters(c,[&](const Cell cc) noexcept {
		  CS.push(cc);
		});
	    }
	  }
# ifdef OCTALTREE_USE_OPENMP
	});
# endif
      DebugInfoN(5,"InteractionTree: DONE (setting anchors and SSE pos)\n");
# ifdef OCTALTREE_USE_OPENMP
      // 4 pad unused final positions
      if(dom+1 == tree->n_dom())
# endif
	for(auto l=tree->end_leaf(); l.I<tree->n_leaf4(); ++l)
	  for(int d=0; d!=Dim; ++d)
	    data->_X16[d][l.I] = 0;
# undef DOM
      DebugInfoN(3,"InteractionTree: DONE "
		 "(setting anchors and SSE positions)\n");
    }
    /// set SSE positions and anchors (if required)
    static void
    set_sse(const OTree*tree, IData*data, 
# ifdef OCTALTREE_USE_OPENMP
	    domain_id dom,
# endif
	    float*QI, float del)
    { return __set_sse<__Anc>(tree,data,
# ifdef OCTALTREE_USE_OPENMP
			      dom,
# endif
			      QI,del); }
    /// (not) dump anchor
    template<bool __A>
    static typename std::enable_if<!__A>::type
    dumpanchor(const ITree*, const Leaf, std::ostream&) {}
    /// (not) dump anchor
    template<bool __A>
    static typename std::enable_if<!__A>::type
    dumpanchor(const ITree*, const ExtLeaf, std::ostream&) {}
    /// dump anchor
    template<bool __A>
    static typename std::enable_if<__A>::type
    dumpanchor(const ITree*tree, const Leaf l, std::ostream&s)
    {
      if(tree->have_anchor()) {
	if(tree->anchor(l)) s << ' ' << std::setw(8) << tree->anchor(l)->Z;
	else                s << "               ---          ";
      }
    }
    /// dump anchor
    template<bool __A>
    static typename std::enable_if<__A>::type
    dumpanchor(const ITree*tree, const ExtLeaf, std::ostream&s)
    {
      if(tree->have_anchor())
	s << "               ---          ";
    }
    /// no anchor head
    template<bool __A>
    static typename std::enable_if<!__A>::type
    headanchor(const ITree*, std::ostream&) {}
    /// anchor head
    template<bool __A>
    static typename std::enable_if<__A>::type
    headanchor(const ITree*tree, std::ostream&s)
    {
      if(tree->have_anchor())
	s << "            anchor        ";
    }
#endif // __SSE__
    /// load masses for internal leaves
    static void load_mass(const OTree*tree, const IData*data, 
			  const Initialiser*init, count_type
#ifdef OCTALTREE_USE_OPENMP
			  dom
#endif
			  )
    {
      DebugInfoN(2,"InteractionTree: loading masses ...\n");
      // 1 load mass for domain leaves
#ifdef OCTALTREE_USE_OPENMP
      cp_domain DOM = tree->domain(dom);
      init->InitMass(tree->template particle_keys<0>()+DOM->L0,
		     p_data(data->_MASS[0])           +DOM->L0,
		     DOM ->n_leaf());
#else
      init->InitMass(tree->template particle_keys<0>(),
		     p_data(data->_MASS[0]),
		     tree->n_leaf());
#endif
      // 2 pad last few
#ifdef OCTALTREE_USE_OPENMP
      if(dom+1 == tree->n_dom())
#endif
	for(count_type l=tree->n_leaf(); l<tree->n_leaf4(); ++l)
	  p_data(data->_MASS[0])[l] = 0;
      DebugInfoN(3,"InteractionTree: DONE (loading masses)\n");
    }
    /// load masses for external leaves
    static void load_mass_ext(const OTree*tree, const IData*data, 
			      const Initialiser*init)
    {
      DebugInfoN(2,"InteractionTree: loading masses for ext leaves ...\n");
      // 1 load mass for external leaves
      init->InitMass(tree->template particle_keys<1>(),
		     p_data(data->_MASS[1]),
		     tree->n_extleaf());
      // 2 pad last few
      for(count_type l=tree->n_extleaf(); l<tree->n_extleaf4(); ++l)
	p_data(data->_MASS[1])[l] = 0;
      DebugInfoN(3,"InteractionTree: DONE (loading masses for ext leaves)\n");
    }
    /// load size^2 for domain leaves
    static void load_size(const OTree*tree, const IData*data, 
			  const Initialiser*init, count_type 
#ifdef OCTALTREE_USE_OPENMP
			  dom
#endif
			  )
    {
      DebugInfoN(2,"InteractionTree: loading sizes^2 ...\n");
      // 1 load size^2 for domain leaves
#ifdef OCTALTREE_USE_OPENMP
      cp_domain DOM = tree->domain(dom);
      init->InitSizeQ(tree->template particle_keys<0>()+DOM->L0,
		      p_data(data->_HSQ[0])            +DOM->L0,
		      DOM ->n_leaf());
#else
      init->InitSizeQ(tree->template particle_keys<0>(),
		      p_data(data->_HSQ[0]),
		      tree->n_leaf());
#endif
      // 2 pad last few rungs
#ifdef OCTALTREE_USE_OPENMP
      if(dom+1 == tree->n_dom())
#endif
      for(count_type l=tree->n_leaf(); l<tree->n_leaf4(); ++l)
	p_data(data->_HSQ[0])[l] = 0;
      DebugInfoN(3,"InteractionTree: DONE (loading sizes^2)\n");
    }
    /// load size^2 for external leaves
    static void load_size_ext(const OTree*tree, const IData*data, 
			      const Initialiser*init)
    {
      DebugInfoN(2,"InteractionTree: loading sizes^2 for ext leaves ...\n");
      // 1 load size^2 for external leaves
      init->InitSizeQ(tree->template particle_keys<1>(),
		      p_data(data->_HSQ[1]),
		      tree->n_extleaf());
      // 2 pad last few rungs
      for(count_type l=tree->n_extleaf(); l<tree->n_extleaf4(); ++l)
	p_data(data->_HSQ[1])[l] = 0;
      DebugInfoN(3,"InteractionTree: DONE (loading sizes^2 for ext leaves)\n");
    }
    /// allocate
    static void allocate(IData*DATA, 
			 int load, count_type NLeaf, count_type NExtLeaf)
    { DATA->allocate(load,NLeaf,NExtLeaf); }
  };// TreeImplenter<InteractionTree>
} // namespace {
////////////////////////////////////////////////////////////////////////////////
namespace WDutils {
  namespace octtree {
    ///
    /// InteractionTreeData::allocate()
    /// 
    template<int __Dim, typename __PosType, typename __PropType
# ifdef __SSE__
	     , bool __Anchoring
# endif
	     >
    void InteractionTreeData<__Dim,__PosType,__PropType
# ifdef __SSE__
			     ,__Anchoring
# endif
			     >::
    allocate(const int load, const count_type nl, const count_type ne)
    {
      static_assert(!__Anchoring,
		    " InteractionTreeData::allocate(): not anchoring");
      DebugInfoN(2,"InteractionTree<>: allocating additional data ...\n");
#ifdef OCTALTREE_DATA_IN_ONE_BLOCK
      count_type need = 0;
      if(load & p_mass)
	need +=
	  Next16<prop_type> (nl) +                 // _MASS[0]
	  Next16<prop_type> (ne) ;                 // _MASS[1]
      if(load & p_size)
	need +=
	  Next16<prop_type> (nl) +                 // _HSQ [0]
	  Next16<prop_type> (ne) ;                 // _HSQ [1] 
# ifdef __SSE__
      if(load & p_ssex)
	need +=
	  Dim * Next16<sse_pos_type>(nl) ;         // _X16 [DIM]
# endif
      _MEM. template reset_conditional<3,2>(need);
# define SET_DATA(Name,Number,Type)					\
      {									\
	if(Number) {							\
	  Name = reinterpret_cast<Type*>(A);				\
	  A += Next16<Type>(Number);					\
	  DebugInfoN(6,"InteractionTree: set %s = %p for %d %s "	\
		     "(%lu bytes)\n",					\
		     __STRING(Name),Name,				\
		     Number,__STRING(Type),Number*sizeof(Type));	\
	} else {							\
	  Name = 0;							\
	  DebugInfoN(6,"set %s = 0\n",__STRING(Name));			\
	}								\
      }
      char*A=_MEM.begin();
      SET_DATA(_MASS[0],(load&p_mass)?nl:0,prop_type);
      SET_DATA(_MASS[1],(load&p_mass)?ne:0,prop_type);
      SET_DATA(_HSQ [0],(load&p_size)?nl:0,prop_type);
      SET_DATA(_HSQ [1],(load&p_size)?ne:0,prop_type);
# ifdef __SSE__
      for(int d=0; d!=Dim; ++d)
	SET_DATA(_X16[d],(load&p_ssex)?nl:0,sse_pos_type);
# endif
      WDutilsAssert(A == _MEM.begin()+need && A <= _MEM.end());
#else // OCTALTREE_DATA_IN_ONE_BLOCK
      const count_type ne4=(ne+3)&~3;              // multiple of 4 >= ne
      const count_type nl4=(nl+3)&~3;              // multiple of 4 >= nl
# define SET_DATA(Name,Number)					\
      if(Number) {						\
	Name.resize(Number);					\
      } else {							\
	Name.clear();						\
	Name.shrink_to_fit();					\
      }								\
      DebugInfoN(6,"InteractionTree: set %s = %p for %d %s\n",	\
		 __STRING(Name),Name.data(),Number,		\
		 nameof(typename std::remove_reference<		\
			decltype(Name[0])>::type));
      SET_DATA(_MASS[0],(load&p_mass)?nl4:0);
      SET_DATA(_MASS[1],(load&p_mass)?ne4:0);
      SET_DATA( _HSQ[0],(load&p_size)?nl4:0);
      SET_DATA( _HSQ[1],(load&p_size)?ne4:0);
# ifdef __SSE__
      SET_DATA(_x16    ,(load&p_ssex)?Dim*nl4:0);
      if(load&p_ssex) {
	_X16[0] = _x16.data();
	for(int d=1; d<Dim; ++d) _X16[d] = _X16[d-1] + nl4;
      } else
	for(int d=0; d<Dim; ++d) _X16[d] = 0;
# endif
#endif// OCTALTREE_DATA_IN_ONE_BLOCK
      DebugInfoN(3,"InteractionTree<>: DONE "
		   "(allocating additional data)\n");
    }
#ifdef __SSE__
    //
    // NOTE: typically, our sse_pos_type is different from base!
    //
    template<int __Dim>
    void InteractionTreeData<__Dim,double,float,true>::
    allocate(const int load, const count_type nl, const count_type ne)
    {
      DebugInfoN(2,"InteractionTree<>: allocating additional data ...\n");


#ifdef OCTALTREE_DATA_IN_ONE_BLOCK
      count_type need = 0;
      if(load & p_mass)
	need +=
	  Next16<prop_type> (nl) +                 // _MASS[0]
	  Next16<prop_type> (ne) ;                 // _MASS[1]
      if(load & p_size)
	need +=
	  Next16<prop_type> (nl) +                 // _HSQ [0]
	  Next16<prop_type> (ne) ;                 // _HSQ [1] 
      if(load & p_ssex)
	need +=
	  Dim * Next16<sse_pos_type>(nl) +         // _X16 [DIM]
	  Next16<Anchor*>   (nl) ;                 // _ANC
      _MEM. template reset_conditional<3,2>(need);
      char*A=_MEM.begin();
      SET_DATA(_MASS[0],(load&p_mass)?nl:0,prop_type);
      SET_DATA(_MASS[1],(load&p_mass)?ne:0,prop_type);
      SET_DATA(_HSQ [0],(load&p_size)?nl:0,prop_type);
      SET_DATA(_HSQ [1],(load&p_size)?ne:0,prop_type);
      for(int d=0; d!=Dim; ++d)
	SET_DATA(_X16[d],(load&p_ssex)?nl:0,sse_pos_type);
      SET_DATA(_ANC,(load&p_ssex)?nl:0,Anchor*);
      WDutilsAssert(A == _MEM.begin()+need && A <= _MEM.end());
#else // OCTALTREE_DATA_IN_ONE_BLOCK
      const count_type ne4=(ne+3)&~3;              // multiple of 4 >= ne
      const count_type nl4=(nl+3)&~3;              // multiple of 4 >= nl
      SET_DATA(_MASS[0],(load&p_mass)?nl4:0);
      SET_DATA(_MASS[1],(load&p_mass)?ne4:0);
      SET_DATA( _HSQ[0],(load&p_size)?nl4:0);
      SET_DATA( _HSQ[1],(load&p_size)?ne4:0);
      SET_DATA(_x16    ,(load&p_ssex)?Dim*nl4:0);
      if(load&p_ssex) {
	_X16[0] = _x16.data();
	for(int d=1; d<Dim; ++d) _X16[d] = _X16[d-1] + nl4;
      } else
	for(int d=0; d<Dim; ++d) _X16[d] = 0;
      SET_DATA(_ANC    ,(load&p_ssex)?nl:0);
#endif// OCTALTREE_DATA_IN_ONE_BLOCK
      DebugInfoN(3,"InteractionTree<>: DONE "
		 "(allocating additional data)\n");
    }
#endif// __SSE__
#undef SET_DATA
  } // namespace WDutils::octtree
  ///
  /// InteractionTree::Builder::Builder
  ///
  __TemplateDecl InteractionTree<__TemplateList>::
  Builder::Builder(DataBase*t, const Initialiser*i, Data const&d)
    : INIT(i)
    , DATA(t)
    , LOAD(d.LOAD)
    , DEL (d.DEL)
#ifdef __SSE__
    , QI  (0)
#endif
  {
    if(LOAD == 0)
      WDutils_WarningN("InteractionTree: load=0\n");
#ifdef __SSE__
    if(LOAD & octtree::p_ssex) {
      if(!(LOAD & octtree::p_size))
	WDutils_THROW("InteractionTree: SSE positions require sizes\n");
      if(DEL >= 1.f)
	WDutils_THROWN("InteractionTree: del=%f >= 1\n",DEL);
      if(DEL < std::numeric_limits<float>::epsilon())
	WDutils_THROWN("InteractionTree: "
		       "del=%f < single precision (%f)\n",
		       DEL,std::numeric_limits<float>::epsilon());
      if(DEL > 0.1f)
	WDutils_WarningN("InteractionTree: "
			 "del=%f close to 1\n",DEL);
      if(DEL < 10*std::numeric_limits<float>::epsilon())
	WDutils_WarningN("InteractionTree: "
			 "del=%f close to %f (single precision)\n",
			 DEL,std::numeric_limits<float>::epsilon());
    }
#endif
  }
  ///
  /// InteractionTree::Builder::Before
  ///
  __TemplateDecl void InteractionTree<__TemplateList>::
  Builder::Before(bool, count_type NCell,
		  count_type NLeaf, count_type NExtLeaf) const
  {
    DebugInfoN(4,"InteractionTree<>::Builder::Before() ...\n");
    using Impl = TreeImplementer<typename InteractionTree::DataBase>;
    Impl::allocate(DATA,LOAD,NLeaf,NExtLeaf);
#ifdef __SSE__
    if(LOAD & p_ssex) {
      DebugInfoN(6,"InteractionTree<>::Builder::Before(): "
		 "resetting QI (NCELL=%d)\n",NCell);
      const_cast<float*&>(QI) = WDutils_NEW(float,NCell);
    }
#endif
    DebugInfoN(5,"InteractionTree<>::Builder::Before(): done\n");
  }
  ///
  /// InteractionTree::Builder::SetSub
  ///
  __TemplateDecl void InteractionTree<__TemplateList>::
  Builder::SetSub(const OTree*tree, count_type dom) const
  {
    using Impl = TreeImplementer<typename InteractionTree::DataBase>;
    // 1 load masses
    if(LOAD & p_mass)
      Impl::load_mass(tree,DATA,INIT,dom);
    // 2 load size^2
    if(LOAD & p_size)
      Impl::load_size(tree,DATA,INIT,dom);
#ifdef __SSE__
    // 3 set SSE positions
    if(LOAD & p_ssex)
      Impl::set_sse(tree,DATA,
# ifdef OCTALTREE_USE_OPENMP
		    dom,
# endif
		    QI,DEL);
#endif
  }
  ///
  /// InteractionTree::Builder::SetTop
  ///
  __TemplateDecl void InteractionTree<__TemplateList>::
  Builder::SetTop(const OTree*tree) const
  {
    using Impl = TreeImplementer<typename InteractionTree::DataBase>;
    if(tree->n_extleaf()) {
      // 1 load masses
      if(LOAD & p_mass)
	Impl::load_mass_ext(tree,DATA,INIT);
      // 2 load size^2
      if(LOAD & p_size)
	Impl::load_size_ext(tree,DATA,INIT);
    }
  }
#undef itree
  ///
  /// InteractionTree::LeafHeader
  ///
  __TemplateDecl std::ostream&InteractionTree<__TemplateList>::
  LeafHeader(std::ostream&s) const
  {
    OTree::LeafHeader(s);
    if(this->have_mass())
      s << "       mass";
    if(this->have_Hsq())
      s << "     size^2";
#ifdef __SSE__
    using Impl = TreeImplementer<typename InteractionTree::DataBase>;
    Impl:: template headanchor<Anchoring>(this,s);
#endif
    return s;
  }
  ///
  /// InteractionTree::Dump(internal leaf)
  ///
  __TemplateDecl std::ostream&InteractionTree<__TemplateList>::
  Dump(const Leaf l, std::ostream&s) const
  {
    OTree::Dump(l,s);
    if(this->have_mass())
      s << ' ' << print(this->mass(l),10,3);
    if(this->have_Hsq())
      s << ' ' << print(this->Hsq(l),10,4);
#ifdef __SSE__
    using Impl = TreeImplementer<typename InteractionTree::DataBase>;
    Impl:: template dumpanchor<Anchoring>(this,l,s);
#endif
    return s;
  }
  ///
  /// InteractionTree::Dump(external leaf)
  ///
  __TemplateDecl std::ostream&InteractionTree<__TemplateList>::
  Dump(const ExtLeaf l, std::ostream&s) const
  {
    OTree::Dump(l,s);
    if(this->have_mass())
      s << ' ' << print(this->mass(l),10,3);
    if(this->have_Hsq())
      s << ' ' << print(this->Hsq(l),10,4);
#ifdef __SSE__
    using Impl = TreeImplementer<typename InteractionTree::DataBase>;
    Impl:: template dumpanchor<Anchoring>(this,l,s);
#endif
    return s;
  }
  //////////////////////////////////////////////////////////////////////////////
#ifdef __SSE__
  template class InteractionTree<2,float ,float ,0>;
  template class InteractionTree<3,float ,float ,0>;
  template class InteractionTree<2,double,float ,0>;
  template class InteractionTree<3,double,float ,0>;
  template class InteractionTree<2,double,double,0>;
  template class InteractionTree<3,double,double,0>;
  template class InteractionTree<2,double,float ,1>;
  template class InteractionTree<3,double,float ,1>;
#else
  template class InteractionTree<2,float ,float >;
  template class InteractionTree<3,float ,float >;
  template class InteractionTree<2,double,float >;
  template class InteractionTree<3,double,float >;
  template class InteractionTree<2,double,double>;
  template class InteractionTree<3,double,double>;
#endif
} // namespace WDutils  
////////////////////////////////////////////////////////////////////////////////
