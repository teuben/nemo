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
/// \version 09-jun-2010 WD  16-byte alignment, using geometry.h
/// \version 15-jul-2010 WD  maximum tree depth = numeric_limits<real>::digits
/// \version 02-jul-2012 WD  added periodic boundary condition
/// \version 10-aug-2012 WD  completely new, OMP parallel in non-public part
///                          requires C++11
/// \version 10-aug-2012 WD  removed OMP stuff to non-public part of library
/// \version 19-sep-2012 WD  alternative memory layout using std::vector
/// \version 30-oct-2012 WD  memory alignment for AVX
/// \version 08-nov-2012 WD  data types pointXX, cubeXX, Dot, Box, TopBox
///                            correctly aligned for Geometry::Algorithms<>
///                          using SSE::packed<> to implement OctalTree::scc()
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
  //
#ifdef OCTALTREE_DATA_IN_ONE_BLOCK
  template<typename _Tp>
  static _Tp*p_data(_Tp*X) { return X; }
#else
  template<typename _Tp>
  static _Tp*p_data(octtree::storage<_Tp> const&X)
  { return const_cast<_Tp*>(X.data()); }
  template<typename _Tp>
  static _Tp*p_data(octtree::aligned_store<_Tp> const&X)
  { return const_cast<_Tp*>(X.data()); }
#endif
  //
  /// tree of boxes and dots, to be used inside an openMP parallel region
  //
  template<int _D, typename _X> 
  struct TreeImplementer< OctalTree<_D,_X> >
  {
    //
    using OTree              = OctalTree<_D,_X>;
    using Helper             = TreeHelper<_D,_X>;
    using Dot                = typename OTree::Dot;
    using Initialiser        = typename OTree::Initialiser;
    using GeoAlg             = Geometry::Algorithms<1>;
    using pos_type           = typename OTree::pos_type;
    using cube               = typename OTree::cube;
    using cubeXX             = typename OTree::cubeXX;
    using point              = typename Helper::point;
    using pointXX            = typename OTree::pointXX;
    using PerBoundary        = typename OTree::PerBoundary;
    using TreeBuilder        = typename OTree::TreeBuilderInterface;
    using peano_map          = typename Peano<_D>::Map;
    using peano_dig          = typename Peano<_D>::Digit;
    using peano_oct          = typename Peano<_D>::Octant;
    //
    static const peano_oct  Dim            = _D;
    static const peano_oct  Nsub           = 1<<Dim;
    static const depth_type MaximumDepth   = OTree::MaximumDepth;
    static const count_type PosAlignment   = OTree::PosAlignment;
    /// \name types used in tree construction
    //@{
    /// represents a cubic cell in the box-dot tree
    /// \note We cannot rely on a default constructor being called when
    ///       allocating boxes for the tree (see documentation for block_alloc).
    /// \note This is carefully designed to have no more than 64, 80, 96, and
    ///       112 bytes for (Dim,pos_type)=(2,float), (2,double), (3,float),
    ///       and (3,double), respectively.
    struct _Box {
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
    using Box = SSE::Extend<_Box,PosAlignment>;
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
    /// \name data
    // data set by constructor
    mutable double                   TIME; ///< wallclock time at start
    const uint8_t                    NMAX; ///< maximum # particles/octant
    const uint8_t                    NMIN; ///< maximum # particles/octant
    const bool                       ASCC; ///< avoid single-child cells
#ifdef OCTALTREE_USE_OPENMP
    depth_type                       IT,NT;///< this thread, # threads
#endif// OCTALTREE_USE_OPENMP
    const count_type                 NMX1; ///< NMAX + 1
    // data set by SetDots()
    octtree::aligned_vector<Dot>     vDOT; ///< locally loaded dots
    Dot                             *DOTS; ///< pter to locally loaded dots
    // data set by Build()
    block_alloc<Box,PosAlignment>   *BXAL; ///< allocator for boxes
    Box                             *ROOT; ///< root box
#ifdef OCTALTREE_USE_OPENMP
# define OCTALTREE_IMPLEMENTER_ADDITIONAL_DATA_MEMBERS
# include <../devel/octtree_omp.cc>
# undef  OCTALTREE_IMPLEMENTER_ADDITIONAL_DATA_MEMBERS
#endif
    // data used in tree linking
    OTree                *TREE;         ///< tree to be linked
    mutable count_type    CF;           ///< cell counter during linking
    mutable count_type    LF;           ///< leaf counter during linking
    //@}
    /// set dots; cumulate Xmin, Xmax, Xsum; called from Build()
    /// \param[in]  init  initialiser for particle data
    /// \param[in]  tree  domain to take order from
    /// \param[out] Xmin  in each dimension: minimum of positions loaded
    /// \param[out] Xmin  in each dimension: maximum of positions loaded
    /// \param[out] Xsum  sum of all positions loaded
    /// \note sets @a vDOT & @a DOTS
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
      static_cast<point&>(TREE->XL[0][L]) = D->X;
      TREE->PL[0][L] = D->I;
      TREE->PC[L] = P;
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
	WDutils_Error("OctalTree: error in cell linking: C=%d > NC=%d\n",C,Cn);
#endif
      TREE->LE[C] = P->LEV + 1;
      GeoAlg::copy(P->CUB,static_cast<cube&>(TREE->XC[C]));
      GeoAlg::shrink(TREE->XC[C],oct);
      TREE->L0[C] = LF;
      TREE->NL[C] = P->NDL[oct];
      TREE->NM[C] = P->NDL[oct];
      TREE->C0[C] = C;
      TREE->NC[C] = 0;
#ifdef OCTALTREE_HAVE_D0
      TREE->D0[C] = IT;
#endif
      TREE->DP[C] = 0;
      for(const Dot*Di=cpDOT(P->OCT[oct]); Di; Di=pDOT(Di->Next))
	LinkLeaf(LF++,Di,C);
    }
    /// load RG[] and set FL[] for domain
    /// \return # active leafs
    /// \param[in] tree   tree to load rungs and flags for
    /// \param[in] ract   set flag = rung>=ract? 0xff:0
    /// \param[in] L0     begin of leaves
    /// \param[in] NL     # leaves
    template<int EXT>
    static count_type LoadRungFlag(OTree*const tree, float ract,
				   count_type L0, count_type NL)
    {
      static_assert(EXT==0 || EXT==1,"unexpected value for EXT");
      tree->INIT->InitRung(tree->template particle_keys<EXT>() + L0,
			   p_data(tree->RG[EXT])+ L0, NL);
      count_type Nact=0;
      for(count_type l=L0; l!=L0+NL; ++l) {
	tree->FL[EXT][l] = tree->RG[EXT][l] >= ract? 0xff : 0;
	if(tree->FL[EXT][l]) ++Nact;
      }
      return Nact;
    }
    /// pass up NA[] for sub-domain
    /// \param[in] tree         tree to pass NA[] up for
    /// \param[in] none_active  are none of the leafs active?
    /// \param[in] all_active   are all of the leafs active?
    /// \param[in] C0           begin of cells to pass NA[] up for
    /// \param[in] CN           end of cells to pass NA[] up for
    static void PassUpNactive(OTree*const tree,
			      bool none_active, bool all_active,
			      count_type C0, count_type CN)
    {
      if(none_active)
	for(count_type c=C0; c!=CN; ++c)
	  tree->NA[c] = 0;
      else if(all_active)
	for(count_type c=C0; c!=CN; ++c)
	  tree->NA[c] = tree->NM[c];
      else
	for(count_type c=CN-1; c!=(C0-1); --c) {
	  tree->NA[c] = 0;
	  if(tree->NL[c])
	    for(count_type ll=tree->L0[c]; ll!=tree->L0[c]+tree->NL[c]; ++ll)
	      if(tree->FL[0][ll]) tree->NA[c]++;
	  if(tree->NC[c])
	    for(count_type cc=tree->C0[c]; cc!=tree->C0[c]+tree->NC[c]; ++cc)
	      tree->NA[c] += tree->NA[cc];
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
			     count_type ne, bool fr,
			     const TreeBuilder*builder=0) WDutils_THROWING;
    /// serial build of an OctalTree
    static void BuildTreeSerial(OTree*tree, bool fresh, bool load_rungs,
				const TreeBuilder*builder);
    /// build an OctalTree
    static count_type BuildTree(OTree*tree, bool fresh, bool load_rungs,
				count_type 
#ifdef OCTALTREE_USE_OPENMP
				tol
#endif
				, const TreeBuilder*builder=0)
    {
#ifdef OCTALTREE_USE_OPENMP
      if(tree->NDOM>1)
	return BuildTreeParallel(tree,fresh,load_rungs,tol,builder);
      else 
#endif
	BuildTreeSerial  (tree,fresh,load_rungs,builder);
      return 0;
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
      DebugInfo(5," %d single-child boxes\n",ns);
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
  template<int _D, typename _X>
  void TreeImplementer<OctalTree<_D,_X> >::
  SetDots(const Initialiser*init, const OTree*tree,
	  point&Xmin, point&Xmax, point&Xsum)
  {
    // 1 allocate & initialise dots
    vDOT = init->InitIntern(
#ifdef OCTALTREE_USE_OPENMP
			    IT,NT,
#endif// OCTALTREE_USE_OPENMP
			    tree? tree->template particle_keys<0>() : 0,
			    tree? tree->NLEAF : 0);
    if(vDOT.empty())
      WDutils_ErrorN("OctalTree: no Dots from by Initialiser::InitIntern()\n");
    DOTS = vDOT.data();
    if(!is_aligned<PosAlignment>(DOTS))
      WDutils_ErrorN("OctalTree: Dots returned by Initialiser::InitIntern() "
		     "not aligned to %d bytes\n",PosAlignment);
    // 2 cumulate Xmin,Xmax,Xsum
    Xmin = Xmax = Xsum = vDOT.cbegin()->X;
    for(auto d:vDOT) {
      d.X.up_min_max(Xmin,Xmax);
      Xsum += d.X;
    }
    // 3 check for nan and inf
    if(isnan(Xsum) || isinf(Xsum)) {
      for(auto d:vDOT)
	if(isnan(d.X) || isinf(d.X)) {
	  std::ostringstream out;
	  out<<"OctalTree: x="<<d.X<<" for particle "<<d.I;
	  WDutils_ErrorN("%s\n",out.str().c_str());
	}
      WDutils_ErrorN("OctalTree: found non-numeric position sum\n");
    }
  }
  //
  // build box-dot tree from all 'our' dots
  //
  template<int _D, typename _X>
  void TreeImplementer< OctalTree<_D,_X> >::
  Build(const Initialiser*init, const OTree*tree)
  {
    // 1 load data into dots
    point Xmin,Xmax,Xave;
    SetDots(init,tree,Xmin,Xmax,Xave);
    // 2 establish global domain
    unsigned Ntot = vDOT.size();
#ifdef OCTALTREE_USE_OPENMP
    if(Ntot <= NMAX)
      WDutils_ErrorN("OctTree: found #particles=%d < nmax=%d\n",Ntot,NMAX);
    if(NT>1) {
      OMP::AllReduce<Parallel::Sum>(Ntot);
      OMP::AllReduce<Parallel::Min>(static_cast<pos_type*>(Xmin),Dim);
      OMP::AllReduce<Parallel::Max>(static_cast<pos_type*>(Xmax),Dim);
      OMP::AllReduce<Parallel::Sum>(static_cast<pos_type*>(Xave),Dim);
    }
#endif// OCTALTREE_USE_OPENMP
    Xave /= pos_type(Ntot);
    // 3 set empty box: local root
    BXAL = new block_alloc<Box,PosAlignment>(NBoxes<Dim>(vDOT.size(),NMAX));
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
    for(Dot*Di=DOTS, *DN=Di+vDOT.size(); Di!=DN; ++Di)
      AddDotToBox(ROOT,Di);
    if(debug(3)) {
      double time=WALLCLOCK();
      DebugInfoN("OctalTree: box-dot tree grown in %f sec (#Dots=%d)\n",
		 time-TIME,int(vDOT.size()));
      TIME = time;
    }
  }
  //
  // count boxes (including root) and dots in the tree
  //
  template<int _D, typename _X>
  void TreeImplementer< OctalTree<_D,_X> >::
  CountNodesSerial(count_type&NB, count_type&ND) const
  {
    NB=1, ND=vDOT.size();
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
  template<int _D, typename _X>
  depth_type TreeImplementer< OctalTree<_D,_X> >::
  LinkCell(const Box*P, const count_type C, const count_type cN) const
  {
#ifdef TESTING
    if(C >= cN)
      WDutils_ErrorN("error in cell linking: C=%d >= cN=%d\n",C,cN);
#endif
    // 1 if single-child box replace by non-single-child box descendant
    if(ASCC) EnsureNonSingleChildBox(P);
    // 2 copy some data, set octant
    TREE->LE[C] = P->LEV;
    GeoAlg::copy(P->CUB,static_cast<cube&>(TREE->XC[C]));
    TREE->L0[C] = LF;
    TREE->NM[C] = P->NUM;
#ifdef OCTALTREE_HAVE_D0
    TREE->D0[C] = IT;
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
    TREE->NL[C] = LF - TREE->L0[C];
    TREE->NC[C] = tmp;
    // 5 link daughter cells, if any
    if(tmp) {
      TREE->C0[C] = CF;
      count_type Ci = CF;
      CF += tmp;
      tmp = 0;
      for(peano_dig d=0; d!=Nsub; ++d) {
	peano_oct oc = P->MAP.octant(d);
	if(P->NDL[oc] > NMAX) {
	  TREE->PA[Ci] = C;
	  TREE->OC[Ci] = oc;
	  update_max(tmp, LinkCell(cpBOX(P->OCT[oc]),Ci++,cN));
	} else if(P->NDL[oc] >= NMIN) {
	  TREE->PA[Ci] = C;
	  TREE->OC[Ci] = oc;
	  LinkOctantCell(P,Ci++,oc,cN);
	}
      }
      ++tmp;
    } else 
      TREE->C0[C] = C;
    // 6 set and return tree depth of cell
    return TREE->DP[C] = tmp;
  }
  //
  // link serial
  //
  template<int _D, typename _X>
  void TreeImplementer< OctalTree<_D,_X> >::
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
    TREE->PA[0] = 0;
    TREE->OC[0] = 0;
    LinkCell(ROOT,0,TREE->NCELL);
    tree->NCELL = CF;
#ifdef OCTALTREE_USE_OPENMP
    DOM->CN= CF; DOM->NC= CF-DOM->C0;
    DOM->LN= LF; DOM->NL= LF-DOM->L0;
    DOM->DEPTH= TREE->DP[0];
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
  template<int _D, typename _X>
  TreeImplementer< OctalTree<_D,_X> >::~TreeImplementer()
  {
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
  template<int _D, typename _X>
  void TreeImplementer< OctalTree<_D,_X> >::
  AllocateTree(OTree*const tree, count_type nc, count_type nl,
	       count_type ne, const bool flag_and_rung,
	       const TreeBuilder*builder) WDutils_THROWING
  {
    tree->NLEAFEXT = ne;
    tree->NLEAF    = nl;
    tree->NCELL    = nc;
    // ensure that ne,nl,nc are multiples of LeafDataBlockSize
    ne= OctalTree<_D,_X>::data_blocked(ne);
    nl= OctalTree<_D,_X>::data_blocked(nl);
#ifdef OCTALTREE_DATA_IN_ONE_BLOCK
    size_t need = NextAligned<pointXX>      (ne)        // XL[1]
      +           NextAligned<particle_key> (ne)        // PL[1]
      +           NextAligned<pointXX>      (nl)        // XL[0]
      +           NextAligned<particle_key> (nl)        // PL[0]
      +           NextAligned<count_type>   (nl)        // PC
      +           NextAligned<depth_type>   (nc)        // LE
      +           NextAligned<octant_type>  (nc)        // OC
      +           NextAligned<cubeXX>       (nc)        // XC
      +           NextAligned<count_type>   (nc)        // L0
      +           NextAligned<local_count>  (nc)        // NL
      +           NextAligned<count_type>   (nc)        // NM
      +           NextAligned<count_type>   (nc)        // C0
      +           NextAligned<octant_type>  (nc)        // NC
      +           NextAligned<count_type>   (nc)        // PA
      +           NextAligned<depth_type>   (nc)        // DP
# ifdef OCTALTREE_HAVE_D0
      +           NextAligned<domain_id>    (nc)        // D0
# endif// OCTALTREE_USE_OPENMP
      ;
    if(flag_and_rung)
      need     += NextAligned<float>        (ne)        // RG[1]
	+         NextAligned<uint8_t>      (ne)        // FL[1]
	+         NextAligned<float>        (nl)        // RG[0]
	+         NextAligned<uint8_t>      (nl)        // FL[0]
	+         NextAligned<count_type>   (nc);       // NA
    // account for derived data
    size_t need_derived_leaves = builder? builder->bytes_for_leaves(nl,ne) : 0;
    size_t need_derived_cells  = builder? builder->bytes_for_cells(nc) : 0;
    need += need_derived_leaves + need_derived_cells;
    // ensure enough memory
    tree->BUF. template reset_conditional<3,2>(need);
# define SET_DATA(Name,Number,Type)				\
    if(Number) {						\
      tree->Name = reinterpret_cast<Type*>(A);			\
      A += NextAligned<Type>(Number);				\
    } else							\
      tree->Name = 0;						\
    DebugInfoN(6,"OctalTree: set %s = %p for %d %s\n",		\
	       __STRING(Name),tree->Name,Number,__STRING(Type));
    // set pointers to memory
    char*A = tree->BUF.begin();
    SET_DATA(XL[1],ne,pointXX);
    SET_DATA(PL[1],ne,particle_key);
    SET_DATA(RG[1],flag_and_rung?ne:0,float);
    SET_DATA(FL[1],flag_and_rung?ne:0,uint8_t);
    SET_DATA(XL[0],nl,pointXX);
    SET_DATA(PL[0],nl,particle_key);
    SET_DATA(PC   ,nl,count_type);
    SET_DATA(RG[0],flag_and_rung?nl:0,float);
    SET_DATA(FL[0],flag_and_rung?nl:0,uint8_t);
    if(builder && need_derived_leaves) {
      builder->set_leaf_memory(A,nl,ne);
      A += need_derived_leaves;
    }
    SET_DATA(LE,nc,depth_type);
    SET_DATA(OC,nc,octant_type);
    SET_DATA(XC,nc,cubeXX);
    SET_DATA(L0,nc,count_type);
    SET_DATA(NL,nc,local_count);
    SET_DATA(NM,nc,count_type);
    SET_DATA(C0,nc,count_type);
    SET_DATA(NC,nc,octant_type);
    SET_DATA(PA,nc,count_type);
    SET_DATA(NA,flag_and_rung?nc:0,count_type);
    SET_DATA(DP,nc,depth_type);
# ifdef OCTALTREE_HAVE_D0
    SET_DATA(D0,nc,domain_id);
# endif// OCTALTREE_USE_OPENMP
    if(builder && need_derived_cells) {
      builder->set_cell_memory(A,nc);
      A += need_derived_cells;
    }
    WDutilsAssert(A == tree->BUF.begin()+need && A <= tree->BUF.end());
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
    SET_DATA(XL[1],ne4);
    SET_DATA(PL[1],ne);
    SET_DATA(RG[1],flag_and_rung?ne4:0);
    SET_DATA(FL[1],flag_and_rung?ne4:0);
    SET_DATA(XL[0],nl4);
    SET_DATA(PL[0],nl);
    SET_DATA(RG[0],flag_and_rung?nl4:0);
    SET_DATA(FL[0],flag_and_rung?nl4:0);
    SET_DATA(PC,nl);
    SET_DATA(LE,nc);
    SET_DATA(OC,nc);
    SET_DATA(XC,nc4);
    SET_DATA(L0,nc);
    SET_DATA(NL,nc);
    SET_DATA(NM,nc);
    SET_DATA(NA,flag_and_rung?nc:0);
    SET_DATA(C0,nc);
    SET_DATA(NC,nc);
    SET_DATA(PA,nc);
    SET_DATA(DP,nc);
# ifdef OCTALTREE_HAVE_D0
    SET_DATA(D0,nc);
# endif// OCTALTREE_USE_OPENMP
#endif
#undef SET_DATA
  }
  //
  // BuildTreeSerial
  // 
  template<int _D, typename _X>
  void TreeImplementer< OctalTree<_D,_X> >::
  BuildTreeSerial(OTree*tree, bool fresh, bool load_rungs,
		  const TreeBuilder*builder=0)
  {
    double    wtime = WALLCLOCK();
    count_type Next = tree->INIT->Nexternal();
    rung_type  Ract = tree->RACT;
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
      if(builder)
	builder->Before(fresh,tree->NCELL,tree->NLEAF,Next);
      AllocateTree(tree,tree->NCELL,tree->NLEAF,Next,load_rungs,builder);
      TI.LinkSerial(tree);
      // 3 load rungs & flags and pass NA[] up the tree
      if(load_rungs) {
	count_type Nact = LoadRungFlag<0>(tree,float(Ract),0,tree->NLEAF);
	if(Nact==0)
	  WDutils_WarningN("OctalTree: no particle has rung >= RACT=%d\n",Ract);
	PassUpNactive(tree,Nact==0,Nact==tree->NLEAF,0,tree->NCELL);
      }
      // 4 destruct TreeImplementer
    }
    // 5 load external particles
    if(Next) {
      for(count_type e=0; e!=Next; ++e)
	tree->INIT->InitExtern(e, static_cast<point&>(tree->XL[1][e]),
			       tree->PL[1][e]);
      if(load_rungs)
	LoadRungFlag<1>(tree,float(Ract),0,Next);
    }
    // 6 set derived data
    if(builder) {
      builder->SetSub(tree,0);
      builder->SetTop(tree);
    }
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
  template<int _D, typename _X>
  OctalTree<_D,_X>::~OctalTree()
  {
#ifdef OCTALTREE_USE_OPENMP
    WDutilsAssertE(!OMP::IsParallel());
    WDutils_DEL_AN(DOM,NDOM);
#endif // OCTALTREE_USE_OPENMP
  }  
  //
  // OctalTree::OctalTree()
  //
  template<int _D, typename _X>
  OctalTree<_D,_X>::OctalTree(const Initialiser*ini, Data const&data)
    throw(exception)
    : INIT  ( ini )
    , PERB  ( INIT->Boundary() )
#ifdef OCTALTREE_USE_OPENMP
    , NDOM  ( RunInfo::omp_threads() )
    , DOM   ( WDutils_NEW(DomainData,NDOM) )
#endif// OCTALTREE_USE_OPENMP
    , NBUILD( 0 )
    , NMAX  ( std::min( MAXNMAX, std::max(depth_type(1),data.NMAX) ) )
    , NMIN  ( std::min( NMAX   , data.NMIN? data.NMIN : depth_type(2) ) )
    , RACT  ( data.RACT )
    , ASCC  ( data.ASCC )
  {
    // sanity checks
    if(0==INIT)
      WDutils_THROW("OctalTree: null pter to Initialiser provided\n");
    if(data.NMAX == 0)
      WDutils_WarningN("OctalTree: nmax=0; will use nmax=%d instead\n",
		      int(NMAX));
    if(data.NMAX > MAXNMAX)
      WDutils_WarningN("OctalTree: nmax=%d exceeds %d; "
		      "will use nmax=%d instead\n",
		      int(data.NMAX),int(MAXNMAX),int(NMAX));
    if(data.NMIN > MAXNMAX)
      WDutils_WarningN("OctalTree: nmin=%d exceeds %d; "
		      "will use nmin=%d instead\n",
		      int(data.NMIN),int(MAXNMAX),int(NMIN));
    else if(data.NMIN > data.NMAX)
      WDutils_WarningN("OctalTree: nmin=%d exceeds nmax=%d; will use nmin=%d\n",
		      int(data.NMIN),int(data.NMAX),int(NMIN));
    // loading rungs?
    bool load_rungs = data.LOAD & octtree::rung_tag;
    if(RACT>0 && !load_rungs)
      WDutils_WarningN("OctalTree: ract=%d > 0 but no rungs to be loaded\n",
		       RACT);
#ifdef OCTALTREE_HAVE_DATA
    // reset DATA[]
    for(int i=0; i!=NDAT; ++i) DATA[0][i] = DATA[1][i] = 0;
#endif
    // build tree
#ifdef OCTALTREE_USE_OPENMP
    TOL = 
#endif// OCTALTREE_USE_OPENMP
      TreeImplementer<OctalTree>::BuildTree(this,1,load_rungs,data.TOL);
    ++NBUILD;
  }
  //
  // OctalTree::OctalTree()
  //
  template<int _D, typename _X>
  OctalTree<_D,_X>::OctalTree(const Initialiser*ini, Data const&data,
				TreeBuilderInterface const&builder)
    throw(exception)
    : INIT  ( ini )
    , PERB  ( INIT->Boundary()  )
#ifdef OCTALTREE_USE_OPENMP
    , NDOM  ( RunInfo::omp_threads() )
    , DOM   ( WDutils_NEW(DomainData,NDOM) )
#endif// OCTALTREE_USE_OPENMP
    , NBUILD( 0 )
    , NMAX  ( std::min( MAXNMAX, std::max(depth_type(1),data.NMAX) ) )
    , NMIN  ( std::min( NMAX   , data.NMIN? data.NMIN : depth_type(2) ) )
//     , NMAX  ( data.NMAX<1? 1 : (data.NMAX>MAXNMAX? MAXNMAX:data.NMAX) )
//     , NMIN  ( data.NMIN? (data.NMIN<NMAX? data.NMIN:NMAX) : (NMAX<2?NMAX:2) )
    , RACT  ( data.RACT )
    , ASCC  ( data.ASCC )
  {
    // sanity checks
    if(0==INIT)
      WDutils_THROW("OctalTree: null pter to Initialiser provided\n");
    if(data.NMAX == 0)
      WDutils_WarningN("OctalTree: nmax=0; will use nmax=%d instead\n",
		      int(NMAX));
    if(data.NMAX > MAXNMAX)
      WDutils_WarningN("OctalTree: nmax=%d > %d; will use nmax=%d instead\n",
		      int(data.NMAX),int(MAXNMAX),int(NMAX));
    if(data.NMIN > MAXNMAX)
      WDutils_WarningN("OctalTree: nmin=%d > %d; will use nmin=%d instead\n",
		      int(data.NMIN),int(MAXNMAX),int(NMIN));
    else if(data.NMIN > data.NMAX)
      WDutils_WarningN("OctalTree: nmin=%d > nmax=%d; will use nmin=%d\n",
		      int(data.NMIN),int(data.NMAX),int(NMIN));
    // loading rungs?
    bool load_rungs = data.LOAD & octtree::rung_tag;
    if(RACT>0 && !load_rungs)
      WDutils_WarningN("OctalTree: ract=%d > 0 but no rungs to be loaded\n",
		       RACT);
#ifdef OCTALTREE_HAVE_DATA
    // reset DATA[]
    for(int i=0; i!=NDAT; ++i) DATA[0][i] = DATA[1][i] = 0;
#endif
    // build tree
#ifdef OCTALTREE_USE_OPENMP
    TOL = 
#endif// OCTALTREE_USE_OPENMP
      TreeImplementer<OctalTree>::BuildTree(this,1,load_rungs,data.TOL,
					    &builder);
    ++NBUILD;
  }
  //
  // OctalTree::rebuild
  //
  template<int _D, typename _X>
  void OctalTree<_D,_X>::rebuild(Data const&data) throw(exception)
  {
    // reset NMAX, NMIN
    if(data.NMAX!=0) {
      NMAX = std::min( MAXNMAX, std::max(depth_type(1),data.NMAX) );
      NMIN = std::min( NMAX   , data.NMIN? data.NMIN : depth_type(2) );
//       NMAX = data.NMAX>MAXNMAX? MAXNMAX:data.NMAX;
//       NMIN =(data.NMIN? (data.NMIN<NMAX?data.NMIN:NMAX) : (NMAX<2?NMAX:2) );
      if(data.NMAX > MAXNMAX)
	WDutils_WarningN("OctalTree::rebuild(): nmax=%d exceeds %d; "
			"will use nmax=%d instead\n",
			int(data.NMAX),int(MAXNMAX),int(NMAX));
      if(data.NMIN > MAXNMAX)
	WDutils_WarningN("OctalTree::rebuild(): nmin=%d exceeds %d; "
			"will use nmin=%d instead\n",
			int(data.NMIN),int(MAXNMAX),int(NMIN));
      else if(data.NMIN > data.NMAX)
	WDutils_WarningN("OctalTree::rebuild(): nmin=%d exceeds nmax=%d; "
			"will use nmin=%d\n",
			int(data.NMIN),int(data.NMAX),int(NMIN));
    }
    // check ASCC
    if(ASCC != data.ASCC)
      WDutils_WarningN("OctalTree: cannot change data.ASCC (from %d to %d) "
		       "in rebuild()\n",ASCC,data.ASCC);
    // reset RACT
    RACT = data.RACT;
    // loading rungs?
    bool load_rungs = data.LOAD & octtree::rung_tag;
    if(RACT>0 && !load_rungs)
      WDutils_WarningN("OctalTree: ract=%d > 0 but no rungs to be loaded\n",
		       RACT);
#ifdef OCTALTREE_HAVE_DATA
    // reset DATA[]
    for(int i=0; i!=NDAT; ++i) DATA[0][i] = DATA[1][i] = 0;
#endif
    // rebuild tree
#ifdef OCTALTREE_USE_OPENMP
    TOL = 
#endif// OCTALTREE_USE_OPENMP
      TreeImplementer<OctalTree>::BuildTree(this,0,load_rungs,data.TOL);
    ++NBUILD;
  }
  //
  // OctalTree::rebuild
  //
  template<int _D, typename _X>
  void OctalTree<_D,_X>::rebuild(Data const&data,
				 TreeBuilderInterface const&builder)
    throw(exception)
  {
    // reset NMAX, NMIN
    if(data.NMAX!=0) {
      NMAX = std::min( MAXNMAX, std::max(depth_type(1),data.NMAX) );
      NMIN = std::min( NMAX   , data.NMIN? data.NMIN : depth_type(2) );
      if(data.NMAX > MAXNMAX)
	WDutils_WarningN("OctalTree::rebuild(): nmax=%d exceeds %d; "
			"will use nmax=%d instead\n",
			int(data.NMAX),int(MAXNMAX),int(NMAX));
      if(data.NMIN > MAXNMAX)
	WDutils_WarningN("OctalTree::rebuild(): nmin=%d exceeds %d; "
			"will use nmin=%d instead\n",
			int(data.NMIN),int(MAXNMAX),int(NMIN));
      else if(data.NMIN > data.NMAX)
	WDutils_WarningN("OctalTree::rebuild(): nmin=%d exceeds nmax=%d; "
			"will use nmin=%d\n",
			int(data.NMIN),int(data.NMAX),int(NMIN));
    }
    // check ASCC
    if(ASCC != data.ASCC)
      WDutils_WarningN("OctalTree: cannot change data.ASCC (from %d to %d) "
		       "in rebuild()\n",ASCC,data.ASCC);
    // reset RACT
    RACT = data.RACT;
    // loading rungs?
    bool load_rungs = data.LOAD & octtree::rung_tag;
    if(RACT>0 && !load_rungs)
      WDutils_WarningN("OctalTree: ract=%d > 0 but no rungs to be loaded\n",
		       RACT);
#ifdef OCTALTREE_HAVE_DATA
    // reset DATA[]
    for(int i=0; i!=NDAT; ++i) DATA[0][i] = DATA[1][i] = 0;
#endif
    // build tree and derived
#ifdef OCTALTREE_USE_OPENMP
    TOL = 
#endif// OCTALTREE_USE_OPENMP
    TreeImplementer<OctalTree>::BuildTree(this,0,load_rungs,data.TOL,&builder);
    ++NBUILD;
  }
  //
  // OctalTree::scc()
  //
  template<int _D, typename _X>
  count_type OctalTree<_D,_X>::scc(point const&x) const noexcept
  {
    using geo = Geometry::Algorithms<1>;
    count_type c=0;
    if(!geo::contains(XC[c],x))
      return c;
    while(NC[c]) {
      octant_type oc=geo::octant(XC[c],x);
      count_type  cc=C0[c], ce=cc+NC[c];
      while(cc!=ce && oc!=OC[cc]) ++cc;
      if(cc==ce || (ASCC && LE[cc] > LE[c]+1 && !geo::contains(XC[cc],x)))
	return c;
      c = cc;
    }
    return c;
  }
  //
  // OctalTree::LeafHeader()
  //
  template<int _Dim, typename _X> template<int _D>
  typename std::enable_if<_D==2,std::ostream&>::type
  OctalTree<_Dim,_X>::_LeafHeader(std::ostream&s) const
  {
    s  << "Leaf     ";
    if(have_flag())
      s<< " f rg";
    s  << "       k up       ";
    s  << "             X            ";
    return s;
  }
  //
  template<int _Dim, typename _X> template<int _D>
  typename std::enable_if<_D==3,std::ostream&>::type
  OctalTree<_Dim,_X>::_LeafHeader(std::ostream&s) const
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
  template<int _Dim, typename _X>
  std::ostream&OctalTree<_Dim,_X>::Dump(const Leaf l, std::ostream&s)
    const
  {
    s  << l;
    if(have_flag())
      s<< ' ' << std::setw(1) << (flag(l)!=0)
       << ' ' << std::setw(2) << int(rung(l));
    s  << ' ' << std::setw(7) << PL[0][l.I]
       << ' ' << std::setfill(' ') << Cell(PC[l.I])
       << ' ' << std::setw(12) << XL[0][l.I];
    return s;
  }
  //
  // OctalTree::Dump(ExtLeaf)
  //
  template<int _Dim, typename _X>
  std::ostream&OctalTree<_Dim,_X>::Dump(const ExtLeaf l, std::ostream&s)
    const
  {
    s  << l << ' ';
    if(have_flag())
      s<< ' ' << std::setw(1) << (flag(l)!=0)
       << ' ' << std::setw(2) << int(rung(l));
    s  << std::setw(7) << PL[1][l.I] << ' '
       << std::setfill(' ') 
       << "nil        "
       << std::setw(12) << XL[1][l.I];
    return s;
  }
  //
  // OctalTree::CellHeader()
  //
  template<int _Dim, typename _X> template<int _D>
  typename std::enable_if<_D==2,std::ostream&>::type
  OctalTree<_Dim,_X>::_CellHeader(std::ostream&s) const
  {
    s  << "Cell     ";
    s  << " le dp";
#ifdef OCTALTREE_USE_OPENMP
    s  << " d0 nd";
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
  template<int _Dim, typename _X> template<int _D>
  typename std::enable_if<_D==3,std::ostream&>::type
  OctalTree<_Dim,_X>::_CellHeader(std::ostream&s) const
  {
    s  << "Cell     ";
    s  << " le dp";
#ifdef OCTALTREE_USE_OPENMP
    s  << " d0 nd";
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
  template<int _Dim, typename _X>
  std::ostream&OctalTree<_Dim,_X>::Dump(const Cell c, std::ostream&s)
    const
  {
    // cell id
    s   << c << ' ';
    // level
    s   << std::setw(2) << int (LE[c.I]) << ' '
	<< std::setw(2) << int (DP[c.I]) << ' ';
#ifdef OCTALTREE_USE_OPENMP
    // domain and # domains
    s   << std::setw(2) << int (first_domain(c)) << ' '
	<< std::setw(2) << int (n_domain(c)) << ' ';
#endif
    // parent
    if(c.I)
      s << Cell(PA[c.I]) << ' ';
    else 
      s << "nil       ";
    // octant
    s   << int (OC[c.I]) << ' ';
    // cells
    if(NC[c.I])
      s << Cell(C0[c.I]) << ' '
	<< int (NC[c.I]) << ' ';
    else
      s << "nil       0 ";
    // leaves
    s   << Leaf(L0[c.I]) << ' '
	<< std::setw(3)  << int(NL[c.I])  << ' '
	<< std::setw(7)  << NM[c.I] << ' ';
    if(have_flag())
      s << std::setw(7)  << NA[c.I] << ' ';
    // radius, centre
    s   << std::setw(10) << XC[c.I].H << ' '
	<< std::setw(11) << XC[c.I].X;
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
  template std::ostream&OctalTree<2,double>::_CellHeader<2>(std::ostream&)
    const;
  template std::ostream&OctalTree<3,double>::_CellHeader<3>(std::ostream&)
    const;
  template std::ostream&OctalTree<2,float>::_CellHeader<2>(std::ostream&)
    const;
  template std::ostream&OctalTree<3,float>::_CellHeader<3>(std::ostream&)
    const;
  //
} // namespace WDutils
////////////////////////////////////////////////////////////////////////////////
namespace {
#ifdef __SSE__

#  define _TemplateDecl						\
  template<int _Dim, typename _PosType, typename _PropType, bool _Anc>
#  define _TemplateList _Dim,_PosType,_PropType,_Anc

#else  // __SSE__

#  define _TemplateDecl					\
  template<int _Dim, typename _PosType, typename _PropType>
#  define _TemplateList _Dim,_PosType,_PropType

#endif
  ///
  /// auxiliary for InteractionTree<>::Builder
  ///
  _TemplateDecl struct TreeImplementer
  < octtree::InteractionTreeData<_TemplateList> >
  {
    using ITree = InteractionTree<_TemplateList>;
    static const int  Dim       = _Dim;
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
    using packed_pos_type= typename ITree::packed_pos_type;
    using pos_comp_block = typename ITree::pos_comp_block;
    using pos_vect_block = typename ITree::pos_vect_block;
    /// set SSE position without anchoring
    template<bool _A>
    static typename std::enable_if<!_A>::type
    _set_packed_pos(const OTree*tree, IData*data,
# ifdef OCTALTREE_USE_OPENMP
	     count_type dom,
# endif
	     float*, float)
    {
      DebugInfoN(2,"InteractionTree: setting sse positions ...\n");
# ifdef OCTALTREE_USE_OPENMP
      auto sub=tree->domain(dom);
# else
      auto sub=tree;
# endif
      for(auto l=sub->begin_leaf(); l!=sub->end_leaf(); ++l) {
	count_type j = packed_pos_type::block_index(l.I);
	count_type k = packed_pos_type::sub_index(l.I);
	for(int d=0; d!=Dim; ++d)
	  data->PPOS[j][d][k] = tree->position(l)[d];
      }
# ifdef OCTALTREE_USE_OPENMP
      if(dom+1 == tree->n_dom())
# endif
	for(auto l=tree->end_leaf(); l.I<tree->leaf_capacity(); ++l) {
	  count_type j = packed_pos_type::block_index(l.I);
	  count_type k = packed_pos_type::sub_index(l.I);
	  for(int d=0; d!=Dim; ++d)
	    data->PPOS[j][d][k] = 0;
	}
      DebugInfoN(3,"InteractionTree: DONE (sse positions)\n");
    }
    //
    template<typename FuncOfCell>
    static void loop_cells_up(const OTree*tree, FuncOfCell f)
      noexcept(noexcept(f))
    { for(auto c=tree->rbegin_cell(); c!=tree->rend_cell(); --c) f(c); }
# ifdef OCTALTREE_USE_OPENMP
    //
    template<typename FuncOfCell>
    static void loop_cells_up(cp_domain dom, FuncOfCell f) noexcept(noexcept(f))
    { dom->loop_cells_up(f); }
# endif
    /// set SSE positions and anchors
    template<bool _A>
    static typename std::enable_if< _A>::type
    _set_packed_pos(const OTree*tree, IData*data,
# ifdef OCTALTREE_USE_OPENMP
		    count_type dom,
# endif
		    float*QI, float del)
    {
      DebugInfoN(2,"InteractionTree: setting anchors and sse positions ...\n");
# ifdef OCTALTREE_USE_OPENMP
      auto sub=tree->domain(dom);
# else
      auto sub=tree;
# endif
      // 1 up-ward pass: set QI[cell] = min(Q)
      DebugInfoN(4,"InteractionTree: passing qmin up the domain ...\n");
      loop_cells_up(sub,[&](const Cell c) noexcept
		    {
		      auto ll=tree->begin_leaf(c);
		      float qmin = data->SQ[0][ll.I];
		      for(++ll; ll<tree->end_leaf_kids(c); ++ll)
			update_min(qmin,data->SQ[0][ll.I]);
		      for(auto cc=tree->begin_cell(c); cc!=tree->end_cell(c);
			  ++cc)
			update_min(qmin,QI[cc.I]);
		      QI[c.I] = qmin;
		    });
      DebugInfoN(5,"InteractionTree: DONE (passing qmin up the domain)\n");
      // 2 count anchors needed
      DebugInfoN(4,"InteractionTree: counting anchors needed ...\n");
      const double fac = square(double(del)/
				std::numeric_limits<float>::epsilon());
      Stack<Cell> CS(OTree::Nsub * sub->depth());
      count_type ianc=0,ntot=0;
# ifdef OCTALTREE_USE_OPENMP
      sub->loop_branches([&](const Cell b) {
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
	      for(auto cc=tree->begin_cell(c); cc!=tree->end_cell(c); ++cc)
		CS.push(cc);
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
	  data->ANCH.resize(ntot);
	}
      DebugInfoN(4,"InteractionTree: setting anchors and SSE pos ...\n");
# ifdef OCTALTREE_USE_OPENMP
# pragma omp barrier
      sub->loop_branches([&](const Cell b) {
	  CS.push(b);
# else
	  CS.push(tree->root());
# endif
	  while(!CS.is_empty()) {
	    auto c = CS.pop();
	    if(square(tree->radius(c)) < fac * QI[c.I]) {
	      Anchor*An = &(data->ANCH[ianc++]);
	      An->Z = tree->centre(c);
	      An->B = tree->begin_leaf(c);
	      An->E = tree->end_leaf_desc(c);
	      for(auto l=tree->begin_leaf(c); l!=tree->end_leaf_desc(c); ++l) {
		data->ANC[l.I] = An;
		count_type j = packed_pos_type::block_index(l.I);
		count_type k = packed_pos_type::sub_index(l.I);
		for(int d=0; d!=Dim; ++d)
		  data->PPOS[j][d][k] = float(tree->position(l)[d] - An->Z[d]);
	      }
	    } else {
	      for(auto l=tree->begin_leaf(c); l<tree->end_leaf_kids(c); ++l) {
		data->ANC[l.I] = 0;
		count_type j = packed_pos_type::block_index(l.I);
		count_type k = packed_pos_type::sub_index(l.I);
		for(int d=0; d!=Dim; ++d)
		  data->PPOS[j][d][k] = float(tree->position(l)[d] -
					      tree->centre  (c)[d] );
	      }
	      for(auto cc=tree->begin_cell(c); cc!=tree->end_cell(c); ++cc)
		CS.push(cc);
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
	for(auto l=tree->end_leaf(); l.I<tree->leaf_capacity(); ++l) {
	  count_type b = packed_pos_type::block_index(l.I);
	  count_type k = packed_pos_type::sub_index(l.I);
	  for(int d=0; d!=Dim; ++d)
	    data->PPOS[b][d][k] = 0;
	}
      DebugInfoN(3,"InteractionTree: DONE "
		 "(setting anchors and SSE positions)\n");
    }
    /// set SSE positions and anchors (if required)
    static void
    set_packed_pos(const OTree*tree, IData*data, 
# ifdef OCTALTREE_USE_OPENMP
		   domain_id dom,
# endif
		   float*QI, float del)
    { return _set_packed_pos<_Anc>(tree,data,
# ifdef OCTALTREE_USE_OPENMP
				   dom,
# endif
				   QI,del); }
    /// (not) dump anchor
    template<bool _A>
    static typename std::enable_if<!_A>::type
    dumpanchor(const ITree*, const Leaf, std::ostream&) {}
    /// (not) dump anchor
    template<bool _A>
    static typename std::enable_if<!_A>::type
    dumpanchor(const ITree*, const ExtLeaf, std::ostream&) {}
    /// dump anchor
    template<bool _A>
    static typename std::enable_if<_A>::type
    dumpanchor(const ITree*tree, const Leaf l, std::ostream&s)
    {
      if(tree->have_anchor()) {
	if(tree->anchor(l)) s << ' ' << std::setw(8) << tree->anchor(l)->Z;
	else                s << "               ---          ";
      }
    }
    /// dump anchor
    template<bool _A>
    static typename std::enable_if<_A>::type
    dumpanchor(const ITree*tree, const ExtLeaf, std::ostream&s)
    {
      if(tree->have_anchor())
	s << "               ---          ";
    }
    /// no anchor head
    template<bool _A>
    static typename std::enable_if<!_A>::type
    headanchor(const ITree*, std::ostream&) {}
    /// anchor head
    template<bool _A>
    static typename std::enable_if<_A>::type
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
      auto sub = tree->domain(dom);
      init->InitMass(tree->template particle_keys<0>()+sub->L0,
		     p_data(data->MS[0])           +sub->L0,
		     sub ->n_leaf());
#else
      init->InitMass(tree->template particle_keys<0>(),
		     p_data(data->MS[0]),
		     tree->n_leaf());
#endif
      // 2 pad last few
#ifdef OCTALTREE_USE_OPENMP
      if(dom+1 == tree->n_dom())
#endif
	for(count_type l=tree->n_leaf(); l<tree->leaf_capacity(); ++l)
	  p_data(data->MS[0])[l] = 0;
      DebugInfoN(3,"InteractionTree: DONE (loading masses)\n");
    }
    /// load masses for external leaves
    static void load_mass_ext(const OTree*tree, const IData*data, 
			      const Initialiser*init)
    {
      DebugInfoN(2,"InteractionTree: loading masses for ext leaves ...\n");
      // 1 load mass for external leaves
      init->InitMass(tree->template particle_keys<1>(),
		     p_data(data->MS[1]),
		     tree->n_extleaf());
      // 2 pad last few
      for(count_type l=tree->n_extleaf(); l<tree->extleaf_capacity(); ++l)
	p_data(data->MS[1])[l] = 0;
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
      auto sub = tree->domain(dom);
      init->InitSizeQ(tree->template particle_keys<0>()+sub->L0,
		      p_data(data->SQ[0])              +sub->L0,
		      sub ->n_leaf());
#else
      init->InitSizeQ(tree->template particle_keys<0>(),
		      p_data(data->SQ[0]),
		      tree->n_leaf());
#endif
      // 2 pad last few rungs
#ifdef OCTALTREE_USE_OPENMP
      if(dom+1 == tree->n_dom())
#endif
      for(count_type l=tree->n_leaf(); l<tree->leaf_capacity(); ++l)
	p_data(data->SQ[0])[l] = 0;
      DebugInfoN(3,"InteractionTree: DONE (loading sizes^2)\n");
    }
    /// load size^2 for external leaves
    static void load_size_ext(const OTree*tree, const IData*data, 
			      const Initialiser*init)
    {
      DebugInfoN(2,"InteractionTree: loading sizes^2 for ext leaves ...\n");
      // 1 load size^2 for external leaves
      init->InitSizeQ(tree->template particle_keys<1>(),
		      p_data(data->SQ[1]),
		      tree->n_extleaf());
      // 2 pad last few rungs
      for(count_type l=tree->n_extleaf(); l<tree->extleaf_capacity(); ++l)
	p_data(data->SQ[1])[l] = 0;
      DebugInfoN(3,"InteractionTree: DONE (loading sizes^2 for ext leaves)\n");
    }
  };// TreeImplenter<InteractionTree>
} // namespace {
////////////////////////////////////////////////////////////////////////////////
namespace WDutils {
  namespace octtree {
#ifdef OCTALTREE_DATA_IN_ONE_BLOCK
    //
    // InteractionTreeData<not anchoring>::bytes_needed()
    //    NOTE: typically, our sse_pos_type is different from base!
    //
    _TemplateDecl size_t InteractionTreeData<_TemplateList>::
    bytes_needed(int load, count_type nl, count_type ne) const
    {
      static_assert(!_Anc," InteractionTreeData::allocate(): not anchoring");
      size_t need = 0;
      if(load & mass_tag)
	need +=
	  NextAligned<prop_type> (nl) +                 // MS[0]
	  NextAligned<prop_type> (ne) ;                 // MS[1]
      if(load & size_tag)
	need +=
	  NextAligned<prop_type> (nl) +                 // SQ[0]
	  NextAligned<prop_type> (ne) ;                 // SQ[1] 
# ifdef __SSE__
      if(load & ppos_tag)
	need += NextAligned<pos_vect_block>
	  (packed_pos_type::num_blocks(nl));            // PPOS
# endif
      DebugInfoN(2,"InteractionTree<>: need %ld bytes of data\n",need);
      return need;
    }
# ifdef __SSE__
    //
    // InteractionTreeData<anchoring>::bytes_needed()
    //
    template<int _Dim>
    size_t InteractionTreeData<_Dim,double,float,true>::
    bytes_needed(int load, count_type nl, count_type ne) const
    {
      size_t need = 0;
      if(load & mass_tag)
	need +=
	  NextAligned<prop_type> (nl) +                 // MS[0]
	  NextAligned<prop_type> (ne) ;                 // MS[1]
      if(load & size_tag)
	need +=
	  NextAligned<prop_type> (nl) +                 // SQ[0]
	  NextAligned<prop_type> (ne) ;                 // SQ[1] 
      if(load & ppos_tag)
	need += NextAligned<pos_vect_block>
	  (packed_pos_type::num_blocks(nl)) +           // PPOS
	  NextAligned<Anchor*>   (nl) ;                 // ANC
      DebugInfoN(4,"InteractionTree<>: need %ld bytes of data\n",need);
      return need;
    }
# endif// __SSE__
    //
    // InteractionTreeData<not anchoring>::set_memory()
    //
    _TemplateDecl void InteractionTreeData<_TemplateList>::
    set_memory(int load, char*buf, count_type nl, count_type ne)
    {
      static_assert(!_Anc," InteractionTreeData::allocate(): not anchoring");
# define SET_DATA(Name,Number,Type)					\
      {									\
	if(Number) {							\
	  Name = reinterpret_cast<Type*>(A);				\
	  A += NextAligned<Type>(Number);				\
	  DebugInfoN(6,"InteractionTree: set %s = %p for %d %s "	\
		     "(%lu bytes)\n",					\
		     __STRING(Name),Name,				\
		     Number,__STRING(Type),Number*sizeof(Type));	\
	} else {							\
	  Name = 0;							\
	  DebugInfoN(6,"InteractionTree: set %s = 0\n",			\
		     __STRING(Name));					\
	}								\
      }
      char*A=buf;
      SET_DATA(MS[0],((load&mass_tag)?nl:0),prop_type);
      SET_DATA(MS[1],((load&mass_tag)?ne:0),prop_type);
      SET_DATA(SQ[0],((load&size_tag)?nl:0),prop_type);
      SET_DATA(SQ[1],((load&size_tag)?ne:0),prop_type);
# ifdef __SSE__
      count_type nb = packed_pos_type::num_blocks(nl);
      SET_DATA(PPOS ,((load&ppos_tag)?nb:0),pos_vect_block);
# endif
      WDutilsAssert(A==buf+bytes_needed(load,nl,ne));
      DebugInfoN(5,"InteractionTree<>: setting data: DONE\n");
    }
# ifdef __SSE__
    //
    // InteractionTreeData<anchoring>::set_memory()
    //    NOTE: typically, our sse_pos_type is different from base!
    //
    template<int _Dim>
    void InteractionTreeData<_Dim,double,float,true>::
    set_memory(int load, char*buf, count_type nl, count_type ne)
    {
      char*A=buf;
      count_type nb = packed_pos_type::num_blocks(nl);
      SET_DATA(MS[0],((load&mass_tag)?nl:0),prop_type);
      SET_DATA(MS[1],((load&mass_tag)?ne:0),prop_type);
      SET_DATA(SQ[0],((load&size_tag)?nl:0),prop_type);
      SET_DATA(SQ[1],((load&size_tag)?ne:0),prop_type);
      SET_DATA(PPOS ,((load&ppos_tag)?nb:0),pos_vect_block);
      SET_DATA(ANC  ,((load&ppos_tag)?nl:0),Anchor*);
      WDutilsAssert(A==buf+bytes_needed(load,nl,ne));
      DebugInfoN(3,"InteractionTree<>: setting data: DONE\n");
    }
# endif// __SSE__
#endif // OCTALTREE_DATA_IN_ONE_BLOCK
    //
    // InteractionTreeData<>::allocate()
    //
#ifndef INTERACTIONTREE_USES_BASE_BLOCK
# ifdef OCTALTREE_DATA_IN_ONE_BLOCK
    // allocate() using own block
    _TemplateDecl void InteractionTreeData<_TemplateList>::
    allocate(int load, const count_type nl, const count_type ne)
    {
      DebugInfoN(5,"InteractionTreeDatat<non-anchoring>::allocate() ...\n");
      size_t need = bytes_needed(load,nl,ne);
      BUF. template reset_conditional<3,2>(need);
      set_memory(load,BUF.begin(),nl,ne);
    }
    template<int _Dim>
    void InteractionTreeData<_Dim,double,float,true>::
    allocate(const int load, const count_type nl, const count_type ne)
    {
      DebugInfoN(5,"InteractionTreeDatat<anchoring>::allocate() ...\n");
      size_t need = bytes_needed(load,nl,ne);
      BUF. template reset_conditional<3,2>(need);
      set_memory(load,BUF.begin(),nl,ne);
    }
# else// OCTALTREE_DATA_IN_ONE_BLOCK
    // allocating using vector<>, non-anchoring
    _TemplateDecl void InteractionTreeData<_TemplateList>::
    allocate(const int load, const count_type nl, const count_type ne)
    {
      static_assert(!_Anc," InteractionTreeData::allocate(): not anchoring");
      DebugInfoN(2,"InteractionTree<>: allocating additional data ...\n");
      const count_type neb=packed_prop_type::blocked_num(ne);
      const count_type nlb=packed_prop_type::blocked_num(nl);
      const count_type nbb=packed_prop_type::num_blocks(nl);
      #  define SET_DATA(Name,Number)  				\
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
      SET_DATA(MS[0],((load&mass_tag)?nlb:0));
      SET_DATA(MS[1],((load&mass_tag)?neb:0));
      SET_DATA(SQ[0],((load&size_tag)?nlb:0));
      SET_DATA(SQ[1],((load&size_tag)?neb:0));
#  ifdef __SSE__
      SET_DATA(PPOS ,((load&ppos_tag)?nbb:0));
#  endif
      DebugInfoN(3,"InteractionTree<>: DONE "
		   "(allocating additional data)\n");
    }
#  ifdef __SSE__
    // allocating using vector<>, anchoring
    //   NOTE: typically, our sse_pos_type is different from base!
    template<int _Dim>
    void InteractionTreeData<_Dim,double,float,true>::
    allocate(const int load, const count_type nl, const count_type ne)
    {
      DebugInfoN(2,"InteractionTree<>: allocating additional data ...\n");
      const count_type neb=packed_prop_type::blocked_num(ne);
      const count_type nlb=packed_prop_type::blocked_num(nl);
      const count_type nbb=packed_prop_type::num_blocks(nl);
      SET_DATA(MS[0],((load&mass_tag)?nlb:0));
      SET_DATA(MS[1],((load&mass_tag)?neb:0));
      SET_DATA(SQ[0],((load&size_tag)?nlb:0));
      SET_DATA(SQ[1],((load&size_tag)?neb:0));
      SET_DATA(PPOS ,((load&ppos_tag)?nbb:0));
      SET_DATA(ANC  ,((load&ppos_tag)?nl:0));
      DebugInfoN(3,"InteractionTree<>: DONE "
		 "(allocating additional data)\n");
    }
#  endif// __SSE__
# endif// OCTALTREE_DATA_IN_ONE_BLOCK
#endif// INTERACTIONTREE_USES_BASE_BLOCK
#undef SET_DATA
  } // namespace WDutils::octtree
  //
  // InteractionTree::Builder::Builder()
  //
  _TemplateDecl InteractionTree<_TemplateList>::
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
    if(LOAD & octtree::ppos_tag) {
      if(!(LOAD & octtree::size_tag))
	WDutils_THROW("InteractionTree: SSE positions require sizes\n");
      if(DEL >= 1.f)
	WDutils_THROW("InteractionTree: del=%f >= 1\n",DEL);
      if(DEL < std::numeric_limits<float>::epsilon())
	WDutils_THROW("InteractionTree: "
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
    DebugInfoN(6,"InteractionTree<>::Builder::Builder(): DATA=%p\n",DATA);
  }
  //
  // InteractionTree::Builder::Before()
  //
  _TemplateDecl void InteractionTree<_TemplateList>::
  Builder::Before(bool, count_type NCell,
		  count_type NLeaf, count_type NExtLeaf) const
  {
    DebugInfoN(4,"InteractionTree<>::Builder::Before() ...\n");
    WDutilsAssert(NLeaf>0);
    WDutilsAssert(NCell>0);
#ifndef INTERACTIONTREE_USES_BASE_BLOCK
    DATA->allocate(LOAD,NLeaf,NExtLeaf); 
#endif
#ifdef __SSE__
    if(LOAD & ppos_tag) {
      DebugInfoN(6,"InteractionTree<>::Builder::Before(): "
		 "resetting QI (NCELL=%d)\n",NCell);
      const_cast<float*&>(QI) = WDutils_NEW(float,NCell);
    }
#endif
    DebugInfoN(5,"InteractionTree<>::Builder::Before(): done\n");
  }
  //
  // InteractionTree::Builder::SetSub()
  //
  _TemplateDecl void InteractionTree<_TemplateList>::
  Builder::SetSub(const OTree*tree, count_type dom) const
  {
    using Impl = TreeImplementer<typename InteractionTree::DataBase>;
    // 1 load masses
    if(LOAD & mass_tag)
      Impl::load_mass(tree,DATA,INIT,dom);
    // 2 load size^2
    if(LOAD & size_tag)
      Impl::load_size(tree,DATA,INIT,dom);
#ifdef __SSE__
    // 3 set SSE positions
    if(LOAD & ppos_tag)
      Impl::set_packed_pos(tree,DATA,
# ifdef OCTALTREE_USE_OPENMP
			   dom,
# endif
			   QI,DEL);
#endif
  }
  //
  // InteractionTree::Builder::SetTop()
  //
  _TemplateDecl void InteractionTree<_TemplateList>::
  Builder::SetTop(const OTree*tree) const
  {
    using Impl = TreeImplementer<typename InteractionTree::DataBase>;
    if(tree->n_extleaf()) {
      // 1 load masses
      if(LOAD & mass_tag)
	Impl::load_mass_ext(tree,DATA,INIT);
      // 2 load size^2
      if(LOAD & size_tag)
	Impl::load_size_ext(tree,DATA,INIT);
    }
  }
#undef itree
  //
  // InteractionTree::LeafHeader()
  //
  _TemplateDecl std::ostream&InteractionTree<_TemplateList>::
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
  //
  // InteractionTree::Dump(internal leaf)
  //
  _TemplateDecl std::ostream&InteractionTree<_TemplateList>::
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
  //
  // InteractionTree::Dump(external leaf)
  //
  _TemplateDecl std::ostream&InteractionTree<_TemplateList>::
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
