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
/// \version 20-apr-2010 WD  parameter nmin, changed in tree building
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
// Wdutils::OctalTree<Dim,Real>
//
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
  /// dummy base class for Dot and Box
  /// \note must be empty so that a cast from Dot to OctTree::Dot is valid!
  struct Node {};
  /// represents leafs
  template<int Dim, typename Real>
  struct Dot : public Node
  {
    typedef typename OctalTree<Dim,Real>::node_index node_index;
    typedef typename OctalTree<Dim,Real>::particle_key particle_key;
    typedef typename OctalTree<Dim,Real>::point point;
    //
    point        X;               ///< position
    particle_key I;               ///< identifier of associated particle
    mutable Dot *Next;            ///< next dot in a linked list
    /// add this to a linked list of dots and increment counter
    void AddToList(Dot*&List, node_index&Counter)
    {
      Next = List;
      List = this;
      ++Counter;
    }
    /// add this to a linked list of nodes
    void AddToList(Node*&List)
    {
      Next = pDOT(List);
      List = this;
    }
  };
  /// type used for octants
  typedef unsigned octant_type;
  /// represents cells
  template<int Dim, typename Real>
  struct Box : public Node
  {
    const static octant_type Nsub = 1<<Dim;
    typedef ::Dot<Dim,Real> Dot;
    typedef typename OctalTree<Dim,Real>::node_index node_index;
    //
    tupel<Dim,Real> CEN;          ///< centre position
    uint8           TYP;          ///< bitfield: 1=cell, 0=dot
    uint8           NBX;          ///< number of daughter boxes
    uint8           NOC;          ///< number of octant occupied
    uint8           LEV;          ///< tree level of box
    Node*           OCT[Nsub];    ///< octants
    node_index      NUM;          ///< # dots
    Dot            *DOT;          ///< linked list of dots, if any.
    /// is octant i a box?
    bool MarkedAsBox(octant_type i) const { return TYP & (1<<i); }
    /// is octant i a dot?
    bool MarkedAsDot(octant_type i) const { return !MarkedAsBox(i); }
    /// octant of dot within box (not checked)
    inline octant_type octant(const Dot*D) const;
    /// Is box a single parent, i.e. has only one sub-box
    bool IsSingleParent() const
    { return DOT==0 && NOC==1; }
    /// mark octant i as being a box and increment daughter box counter
    void MarkAsBox(octant_type i)
    {
      TYP |= (1<<i);
      ++NBX;
    }
    /// reset octants
    Box&ResetOctants()
    {
      for(Node**B=OCT; B!=OCT+Nsub; ++B) *B=0;
      return*this;
    }
    /// reset data
    Box&Reset()
    {
//       reinterpret_cast<uint32&>(TYP) = 0;
      TYP = 0;
      NBX = 0;
      NOC = 0;
      NUM = 0;
      DOT = 0;
      return ResetOctants();
    }
  };
  /// \name some geometrical methods, only for D=2,3
  //@{
  /// is pos in the interval [cen-rad, cen+rad) ?
  template<typename Real> inline
  bool ininterval(Real cen, Real rad, Real pos)
  // necessary to trick emacs, which otherwise confused "<" for a bracket
#define LessThan <
  { return pos LessThan cen?  cen <= pos+rad : pos < cen+rad; }
#undef  LessThan
  /// contains geometrical methods in Dim dimensions, only D=2,3
  template<int Dim, typename Real> struct Helper;
  //
  template<> struct Helper<2,float>
  {
    typedef float         Real;
    typedef tupel<2,Real> point;
    typedef ::Box<2,Real> Box;
    static octant_type octant(point const&cen, point const&pos)
    {
#ifdef __SSE__
      return 3&_mm_movemask_ps(_mm_cmplt_ps(_mm_loadu_ps(cen),
					    _mm_loadu_ps(pos)));
#else
      octant_type oct(0);
      if(pos[0] > cen[0]) oct |= 1;
      if(pos[1] > cen[1]) oct |= 2;
      return oct;
#endif
    }
    static bool contains(point const&cen, Real rad, point const&pos)
    {
      return ininterval(cen[0],rad,pos[0])
	&&   ininterval(cen[1],rad,pos[1]);
    }
    static Real outside_dist_sq(point const&cen, Real rad, point const&pos)
    {
      Real q(0),D;
      D = abs(cen[0]-pos[0]); if(D>rad) q+=square(D-rad);
      D = abs(cen[1]-pos[1]); if(D>rad) q+=square(D-rad);
      return q;
    }
    static bool outside(point const&cen, Real rad, point const&pos, Real Q)
    {
      Real q(0),D;
      D=abs(cen[0]-pos[0]); if(D>rad && Q<(q+=square(D-rad))) return true;
      D=abs(cen[1]-pos[1]); if(D>rad && Q<(q+=square(D-rad))) return true;
      return false;
    }
    static bool inside(point const&cen, Real rad, point const&pos, Real Q)
    {
      Real D;
      D=abs(cen[0]-pos[0]); if(D>rad || Q>square(D-rad)) return false;
      D=abs(cen[1]-pos[1]); if(D>rad || Q>square(D-rad)) return false;
      return true;
    }
    static tupel<2,Real> Integer(point const&x)
    {
      tupel<2,Real> c;
      c[0]=int(x[0]+Real(0.5));
      c[1]=int(x[1]+Real(0.5));
      return c;
    }
    static bool ShrinkToOctant(Box*B, octant_type i, uint8 m, const Real*ra)
    {
      uint8 l = ++(B->LEV);
      if(l > m) return false;
      Real rad=ra[l];
      if(i&1) B->CEN[0] += rad;  else  B->CEN[0] -= rad;
      if(i&2) B->CEN[1] += rad;  else  B->CEN[1] -= rad;
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
  //
  template<> struct Helper<2,double>
  {
    typedef double        Real;
    typedef tupel<2,Real> point;
    typedef ::Box<2,Real> Box;
    static octant_type octant(point const&cen, point const&pos)
    {
      octant_type oct(0);
      if(pos[0] > cen[0]) oct |= 1;
      if(pos[1] > cen[1]) oct |= 2;
      return oct;
    }
    static bool contains(point const&cen, Real rad, point const&pos)
    {
      return ininterval(cen[0],rad,pos[0])
	&&   ininterval(cen[1],rad,pos[1]);
    }
    static Real outside_dist_sq(point const&cen, Real rad, point const&pos)
    {
      Real q(0),D;
      D = abs(cen[0]-pos[0]); if(D>rad) q+=square(D-rad);
      D = abs(cen[1]-pos[1]); if(D>rad) q+=square(D-rad);
      return q;
    }
    static bool outside(point const&cen, Real rad, point const&pos, Real Q)
    {
      Real q(0),D;
      D=abs(cen[0]-pos[0]); if(D>rad && Q<(q+=square(D-rad))) return true;
      D=abs(cen[1]-pos[1]); if(D>rad && Q<(q+=square(D-rad))) return true;
      return false;
    }
    static bool inside(point const&cen, Real rad, point const&pos, Real Q)
    {
      Real D;
      D=abs(cen[0]-pos[0]); if(D>rad || Q>square(D-rad)) return false;
      D=abs(cen[1]-pos[1]); if(D>rad || Q>square(D-rad)) return false;
      return true;
    }
    static tupel<2,Real> Integer(point const&x)
    {
      tupel<2,Real> c;
      c[0]=int(x[0]+Real(0.5));
      c[1]=int(x[1]+Real(0.5));
      return c;
    }
    static bool ShrinkToOctant(Box*B, octant_type i, uint8 m, const Real*ra)
    {
      uint8 l = ++(B->LEV);
      if(l > m) return false;
      Real rad=ra[l];
      if(i&1) B->CEN[0] += rad;  else  B->CEN[0] -= rad;
      if(i&2) B->CEN[1] += rad;  else  B->CEN[1] -= rad;
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
  //
  template<> struct Helper<3,float>
  {
    typedef float         Real;
    typedef tupel<3,Real> point;
    typedef ::Box<3,Real> Box;
    static octant_type octant(point const&cen, point const&pos)
    {
#ifdef __SSE__
      return 7 & _mm_movemask_ps(_mm_cmplt_ps(_mm_loadu_ps(cen),
					      _mm_loadu_ps(pos)));
#else
      octant_type oct(0);
      if(pos[0] > cen[0]) oct |= 1;
      if(pos[1] > cen[1]) oct |= 2;
      if(pos[2] > cen[2]) oct |= 4;
      return oct;
#endif
    }
    static bool contains(point const&cen, Real rad, point const&pos)
    {
#ifdef __SSE__
      __m128 CC = _mm_loadu_ps(cen);
      __m128 RR = _mm_set1_ps (rad);
      __m128 XX = _mm_loadu_ps(pos);
      return 7 == (7 & 
      _mm_movemask_ps(_mm_and_ps(_mm_cmple_ps(CC,_mm_add_ps(XX,RR)),
				 _mm_cmplt_ps(XX,_mm_add_ps(CC,RR)))));
#else
      return ininterval(cen[0],rad,pos[0])
	&&   ininterval(cen[1],rad,pos[1])
	&&   ininterval(cen[2],rad,pos[2]);
#endif
    }
    static Real outside_dist_sq(point const&cen, Real rad, point const&pos)
    {
#ifdef __SSE__
      SSE::Traits<float>::vector Q;
      __m128 R = _mm_set1_ps (rad);
      __m128 C = _mm_loadu_ps(cen);
      __m128 X = _mm_loadu_ps(pos);
      __m128 D = _mm_sub_ps(_mm_max_ps(C,X),_mm_min_ps(C,X));
      __m128 DR= _mm_sub_ps(D,R);
      _mm_store_ps(Q,_mm_and_ps(_mm_cmplt_ps(R,D),_mm_mul_ps(DR,DR)));
      return Q[0]+Q[1]+Q[2];
#else
      Real q(0),D;
      D=abs(cen[0]-pos[0]); if(D>rad) q+=square(D-rad);
      D=abs(cen[1]-pos[1]); if(D>rad) q+=square(D-rad);
      D=abs(cen[2]-pos[2]); if(D>rad) q+=square(D-rad);
      return q;
#endif
    }
    static bool outside(point const&cen, Real rad, point const&pos, Real Q)
    {
#ifdef __SSE__
      return  outside_dist_sq(cen,rad,pos) > Q;
#else
      Real q(0),D;
      D=abs(cen[0]-pos[0]); if(D>rad && Q<(q+=square(D-rad))) return true;
      D=abs(cen[1]-pos[1]); if(D>rad && Q<(q+=square(D-rad))) return true;
      D=abs(cen[2]-pos[2]); if(D>rad && Q<(q+=square(D-rad))) return true;
      return false;
#endif
    }
    static bool inside(point const&cen, Real rad, point const&pos, Real Q)
    {
      Real D;
      D=abs(cen[0]-pos[0]); if(D>rad || Q>square(D-rad)) return false;
      D=abs(cen[1]-pos[1]); if(D>rad || Q>square(D-rad)) return false;
      D=abs(cen[2]-pos[2]); if(D>rad || Q>square(D-rad)) return false;
      return true;
    }
    static tupel<3,Real> Integer(point const&x)
    {
      tupel<3,Real> c;
      c[0]=int(x[0]+Real(0.5));
      c[1]=int(x[1]+Real(0.5));
      c[2]=int(x[2]+Real(0.5));
      return c;
    }
    static bool ShrinkToOctant(Box*B, octant_type i, uint8 md, const Real*ra)
    {
      uint8 l = ++(B->LEV);
      if(l > md) return false;
      Real rad=ra[l];
      if(i&1) B->CEN[0] += rad;  else  B->CEN[0] -= rad;
      if(i&2) B->CEN[1] += rad;  else  B->CEN[1] -= rad;
      if(i&4) B->CEN[2] += rad;  else  B->CEN[2] -= rad;
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
  //
  template<> struct Helper<3,double>
  {
    typedef double        Real;
    typedef tupel<3,Real> point;
    typedef ::Box<3,Real> Box;
    static octant_type octant(point const&cen, point const&pos)
    {
      octant_type oct(0);
      if(pos[0] > cen[0]) oct |= 1;
      if(pos[1] > cen[1]) oct |= 2;
      if(pos[2] > cen[2]) oct |= 4;
      return oct;
    }
    static bool contains(point const&cen, Real rad, point const&pos)
    {
      return ininterval(cen[0],rad,pos[0])
	&&   ininterval(cen[1],rad,pos[1])
	&&   ininterval(cen[2],rad,pos[2]);
    }
    static Real outside_dist_sq(point const&cen, Real rad, point const&pos)
    {
      Real q(0),D;
      D = abs(cen[0]-pos[0]); if(D>rad) q+=square(D-rad);
      D = abs(cen[1]-pos[1]); if(D>rad) q+=square(D-rad);
      D = abs(cen[2]-pos[2]); if(D>rad) q+=square(D-rad);
      return q;
    }
    static bool outside(point const&cen, Real rad, point const&pos, Real Q)
    {
      Real q(0),D;
      D=abs(cen[0]-pos[0]); if(D>rad && Q<(q+=square(D-rad))) return true;
      D=abs(cen[1]-pos[1]); if(D>rad && Q<(q+=square(D-rad))) return true;
      D=abs(cen[2]-pos[2]); if(D>rad && Q<(q+=square(D-rad))) return true;
      return false;
    }
    static bool inside(point const&cen, Real rad, point const&pos, Real Q)
    {
      Real D;
      D=abs(cen[0]-pos[0]); if(D>rad || Q>square(D-rad)) return false;
      D=abs(cen[1]-pos[1]); if(D>rad || Q>square(D-rad)) return false;
      D=abs(cen[2]-pos[2]); if(D>rad || Q>square(D-rad)) return false;
      return true;
    }
    static tupel<3,Real> Integer(point const&x)
    {
      tupel<3,Real> c;
      c[0]=int(x[0]+Real(0.5));
      c[1]=int(x[1]+Real(0.5));
      c[2]=int(x[2]+Real(0.5));
      return c;
    }
    static bool ShrinkToOctant(Box*B, octant_type i, uint8 md, const Real*ra)
    {
      uint8 l = ++(B->LEV);
      if(l > md) return false;
      Real rad=ra[l];
      if(i&1) B->CEN[0] += rad;  else  B->CEN[0] -= rad;
      if(i&2) B->CEN[1] += rad;  else  B->CEN[1] -= rad;
      if(i&4) B->CEN[2] += rad;  else  B->CEN[2] -= rad;
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
  /// octant of pos with respect to cen.
  /// \note if pos[i]>=cen[i], the ith bit of octant is set to 1, otherwise 0
  template<int D, typename Real> inline
  octant_type octant(tupel<D,Real> const&cen,
		     tupel<D,Real> const&pos)
  { return Helper<D,Real>::octant(cen,pos); }
  /// does a cubic box contain a given position
  /// \param[in] cen  geometric centre of cube
  /// \param[in] rad  radius = half side length of cube
  /// \param[in] pos  position to test for containment
  /// \return         is pos[i] in [cen[i]-rad, cen[i]+rad) for i=0...D-1 ?
  /// \note This definition of containment matches the way positions are
  ///       sorted into the OctalTree.
  template<int D, typename Real> inline
  bool contains(tupel<D,Real> const&cen, Real rad, tupel<D,Real> const&pos)
  { return Helper<D,Real>::contains(cen,rad,pos); }
  /// distance^2 from given position to the nearest point on a cube
  /// \param[in] cen  geometric centre of cube
  /// \param[in] rad  radius = half side length of cube
  /// \param[in] pos  position to compute distance^2 for
  /// \return    squared distance of @a pos to cube
  /// \note If pos is inside the cube, zero is returned
  template<int D, typename Real> inline
  Real outside_dist_sq(tupel<D,Real> const&cen, Real rad,
		       tupel<D,Real> const&pos)
  { return Helper<D,Real>::outside_dist_sq(cen,rad,pos); }
  /// is a sphere outside of a cubic box?
  /// \param[in] cen  geometric centre of cube
  /// \param[in] rad  radius = half side length of cube
  /// \param[in] pos  centre of sphere
  /// \param[in] q    radius^2 of sphere
  /// \return is sphere outside cube?
  /// \note Equivalent to, but on average faster than, 
  ///       \code q < outside_dist_sq(cen,rad,pos) \endcode
  template<int D, typename Real> inline
  bool outside(tupel<D,Real> const&cen, Real rad,
	       tupel<D,Real> const&pos, Real q)
  { return Helper<D,Real>::outside(cen,rad,pos,q); }
  /// is a sphere inside of a cubic box?
  /// \param[in] cen  geometric centre of cube
  /// \param[in] rad  radius = half side length of cube
  /// \param[in] pos  centre of sphere
  /// \param[in] q    radius^2 of sphere
  template<int D, typename Real> inline
  bool inside(tupel<D,Real> const&cen, Real rad,
	      tupel<D,Real> const&pos, Real q)
  { return Helper<D,Real>::inside(cen,rad,pos,q); }
  //
  template<int Dim, typename Real> inline
  octant_type Box<Dim,Real>::octant(const Dot*D) const
  { return Helper<Dim,Real>::octant(this->CEN,D->X); }
  /// base for BoxDotTree
  template<int Dim, typename Real>
  struct DotInitialiser
  {
    typedef OctalTree<Dim,Real>           OctTree;
    typedef typename OctTree::Initialiser Initialiser;
    typedef typename OctTree::node_index  node_index;
    typedef ::Dot<Dim,Real>               Dot;
    node_index  NDOT;          ///< number of dots to load
    Dot  *const D0;            ///< begin of dots
    Dot  *const DN;            ///< end of dots
    /// ctor: initialise dots
    DotInitialiser(char build, node_index Ndot,
		   const Initialiser*Init,const OctTree*Tree)
    WDutils_THROWING WD_HOT;
    /// dtor: de-allocate data
    ~DotInitialiser()
    { if(D0) WDutils_DEL_A(D0); }
  };
  template<int D, typename Real>
  DotInitialiser<D,Real>::DotInitialiser(char build, node_index Ndot,
					 const Initialiser*Init,
					 const OctTree*Tree) WDutils_THROWING
    : NDOT(Ndot), D0(NDOT? WDutils_NEW(Dot,NDOT) : 0), DN(D0+NDOT)
  {
    switch(build) {
    case 'n':
      // build from scratch: initialise all dots
      for(Dot*Di=D0; Di!=DN; ++Di)
	Init->Initialise(Di->I,Di->X);
      break;
    case 'r': {
      // re-building using the old tree order
      Dot*Di=D0, *List=0, *DN1=D0+min(NDOT,Tree->Nleafs());
      // 1  try to re-initialise dot positions from particle key in old tree
      //    put uninitialised dots in linked list
      for(node_index i=0; Di!=DN1; ++Di,++i) {
	Di->I = Tree->PL[i];
	if(! Init->ReInitialiseValid(Di->I,Di->X) ) { Di->Next=List; List=Di; }
      }
      // 2  initialise particle keys and positions for dots in linked list
      for(Di=List; Di; Di=Di->Next)
	Init->ReInitialiseInvalid(Di->I,Di->X);
      // 3  initialise particle keys and positions for any remaining dots
      for(Di=DN1; Di!=DN; ++Di)
	Init->ReInitialiseInvalid(Di->I,Di->X);
    } break;
    case 'p': {
      // pruning a given (parent) tree
      if(0==NDOT) {
	// A  unknown number of particles in pruned tree
	// 1  count number of particles in pruned tree
	for(node_index i=0; i!=Tree->Nleafs(); ++i)
	  if(Init->Pick(Tree->PL[i])) ++NDOT;
	if(0==NDOT)
	  WDutils_THROW("OctalTree<%d,%s>::build(): empty tree\n",
			D,nameof(Real));
	// 2  allocate dots
	const_cast<Dot*&>(D0) = WDutils_NEW(Dot,NDOT);
	const_cast<Dot*&>(DN) = D0+NDOT;
	// 3  initialise dots
	Dot*Di=D0;
	for(node_index i=0; i!=Tree->Nleafs(); ++i)
	  if(Init->Pick(Tree->PL[i])) {
	    Di->I = Tree->PL[i];
	    Di->X = Tree->XL[i];
	    ++Di;
	  }
      } else {
	// B  assume no more than Ndot leafs of parent tree are in pruned tree
	Dot*Di=D0;
	for(node_index i=0; i!=Tree->Nleafs(); ++i)
	  if(Init->Pick(Tree->PL[i])) {
	    if(Di==DN)
	      WDutils_THROW("OctalTree<%d,%s>::build(): "
			    "more leafs in pruned tree than expected (%d)\n",
			    D,nameof(Real),NDOT);
	    Di->I = Tree->PL[i];
	    Di->X = Tree->XL[i];
	    ++Di;
	  }
	if(Di < DN) NDOT = Di-D0;
      }
    } break;
    default: WDutils_THROW("OctalTree<%d,%s>::build(): unknown build '%c'\n",
			   D,nameof(Real),build);
    }
  }
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
  //
  /// auxiliary class for OctalTree.
  ///
  /// An OctalTree is build by first making a BoxDotTree, then mapping it to
  /// an OctalTree via BoxDotTree::Link()
  template<int Dim, typename Real>
  struct BoxDotTree : public DotInitialiser<Dim,Real>
  {
    const static octant_type Nsub = 1<<Dim; ///< number of octants per cell
    //
    typedef DotInitialiser<Dim,Real>      Base;
    typedef OctalTree<Dim,Real>           OctTree;
    typedef typename OctTree::Initialiser Initialiser;
    typedef typename OctTree::node_index  node_index;
    typedef typename OctTree::depth_type  depth_type;
    typedef typename OctTree::local_count local_count;
    typedef typename OctTree::point       point;
    typedef ::Dot<Dim,Real>               Dot;
    typedef ::Box<Dim,Real>               Box;
    const static depth_type MAXD = OctTree::MaximumDepth;
    /// \name data
    //@{
    Base::NDOT;
    Base::D0;
    Base::DN;
    node_index          NCRIT;         ///< critical number dots/box
    /// \note in constructor: NCRIT=Nmax in linking: NCRIT=Nmin
    block_alloc<Box>    BM;            ///< allocator for boxes
    Real                RA[MAXD+1];    ///< array with radius(level)  
    Box                *P0;            ///< root box
    node_index          NCELL;         ///< # cells linked
    depth_type          DEPTH;         ///< depth of linked tree
    mutable node_index  CF;            ///< free cells during linking
    mutable node_index  LF;            ///< free leafs during linking
    mutable size_t      ND;            ///< # dots added sofar
    const OctTree      *TREE;          ///< tree to be linked
    //@}
    /// replace (pter to) box by (pter to) its non-single-parent descendant
    /// \param[in,out] P  (pter to) given box
    void EnsureNonSingleParent(const Box*&P) const
    {
      while(P->IsSingleParent())
	for(octant_type i=0; i!=Nsub; ++i)
	  if(P->OCT[i]) {
	    P = pBOX(P->OCT[i]);
	    break;
	  }
    }
#ifdef TESTING
    /// Find first non-single parent descendant
    /// \param[in] P   given box
    /// \return    P or its first non-single-parent descendant
    const Box* NonSingleParent(const Box*B) const
    {
      const Box*P = B;
      EnsureNonSingleParent(P);
      return P;
    }
#endif
    /// \name method used in construction
    //@{
    /// shrink box to its octant i
    /// \param[in,out] B Box to shrink
    /// \param[in]     i octant
    bool ShrinkToOctant(Box*B, octant_type i)
    { return Helper<Dim,Real>::ShrinkToOctant(B,i,MAXD,RA); }
    /// get a new box, resetted.
    Box* NewBox()
    { return &(BM.new_element(EstimateNalloc(NDOT,ND))->Reset()); }
    /// provides a new empty (daughter) box in the i th octant of B
    /// \return new box
    /// \param[in] B  parent box
    /// \param[in] i  parent box's octant
    Box*MakeSubBox(const Box*B, octant_type i) WDutils_THROWING
    {
      Box *P = NewBox();
      P->LEV = B->LEV;
      P->CEN = B->CEN;
      if(!ShrinkToOctant(P,i))
	WDutils_THROW("exceeding maximum tree depth of %d\n         "
		      "(perhaps more than Nmax=%du positions are identical "
		      "within floating-point precision)\n", MAXD, NCRIT);
      return P;                                    // return new box            
    }
    /// splits a box.
    ///
    /// The Dots in the Box's linked list are sorted into octants. Octants
    /// with one Dot will just hold that Dot, octants with many Dots will be
    /// Boxes with the Dots in the linked list. If all Dots happen to be in
    /// just one octant, the process is repeated on the Box of this octant.
    ///
    /// \param[in,out] P  box to be splitted
    void SplitBox(Box*P) WDutils_THROWING
    {
      for(;;) {
	Box*S=0;
	for(Dot*Dn,*Di=P->DOT; Di; Di=Dn) {
	  Dn=Di->Next;
	  octant_type b = P->octant(Di);
	  if(P->OCT[b]==0) {
	    // octant empty: add dot and increment box's octant counter
	    Di->AddToList(P->OCT[b]);
	    P->NOC++;
	  } else if(P->MarkedAsDot(b)) {
	    // octant contains another dot: add dot to list, then make new box
	    Di->AddToList(P->OCT[b]);
	    S         = MakeSubBox(P,b);
	    S->DOT    = pDOT(P->OCT[b]);
	    S->NUM    = 2;
	    P->OCT[b] = S;
	    P->MarkAsBox(b);
	  } else {
	    // octant holds box: add dot to box
	    S = pBOX(P->OCT[b]);
	    Di->AddToList(S->DOT, S->NUM);
	  }
	}
	P->DOT=0;
	if(P->NOC > 1) break;
	P = S;
      }
    }
    /// adds a dot to a box.
    ///
    /// If the dot's octant is empty, the dot is put into that octant and we
    /// are done. If there is already another dot, a new daughter box is
    /// created hold just this dot and the process then repeated on the daughter
    /// box. Finally, if the octant contains already a box, the process is
    /// repeated on that box.
    /// \param[in] P     box to add to
    /// \param[in] Di    dot to add
    void AddDot(Box*P, Dot*Di) WDutils_THROWING
    {
      // loop boxes, starting with root
      for(;;) {
	if(P->DOT) {
	  // box is final: add dot, split box if N > NCRIT  -->  done
	  Di->AddToList(P->DOT,P->NUM);
	  if(P->NUM > NCRIT) SplitBox(P);
	  return;
	} else {
	  // non-final box: find octant, increment N
	  octant_type b = P->octant(Di);
	  Node **oc = P->OCT+b;
	  P->NUM++;
	  if((*oc)==0) {
	    // octant empty: put dot into octant, increment oct counter, return
	    *oc = Di;
	    P->NOC++;
	    return;
	  } else if(P->MarkedAsDot(b)) {
	    // octant contains another dot: make new box, then set P=new box
	    P->MarkAsBox(b);
	    P= MakeSubBox(P,b);
	    pDOT(*oc)->AddToList(P->DOT,P->NUM);
	    *oc = P;
	  } else
	    // octant is box: P = box
	    P = pBOX(*oc);
	}
      }
    }
    /// ctor: build a BoxDotTree
    /// If old tree given, we add the dots in the order of its leafs
    /// \param[in] build type of tree build to perform
    /// \param[in] Ndot  number of positions
    /// \param[in] Init  initializer to re-initiliase particle data
    /// \param[in] Nmax  max positions / cell
    /// \param[in] Tree  old OctalTree
    BoxDotTree(char build, node_index Ndot, const Initialiser*Init,
	       local_count Nmax, const OctTree*Tree=0) WDutils_THROWING WD_HOT;
    //@}
    /// \name method used in linking to OctalTree
    //@{
    /// tree linking: link leafs in box not to be split
    /// \param[in] P  box to link dots of
    /// \param[in] C  index of parent cell for box
    /// \note recursive.
    /// \note assumes the P->NUM < NCRIT
    /// \node changes LF
    void LinkLeafs(const Box*B, node_index C) const WD_HOT;
    /// tree linking: leaf & cell descendants are continuous in memory
    /// \param[in] P  box to link with C
    /// \param[in] C  index of current cell to be linked
    /// \param[in] o  octant of box B in parent
    /// \return       tree depth of cell C
    /// \note recursive.
    /// \note assumes the P->NUM >= NCRIT
    /// \node uses data CF and LF
    depth_type LinkCell(const Box*P, node_index C, octant_type o) const WD_HOT;
    /// copy Dot data to Leaf
    /// \param[in] L  leaf index in tree
    /// \param[in] D  dot to link leaf to
    /// \param[in] P  parent cell index for leaf
    void LinkLeaf(node_index L, const Dot*D, node_index P) const
    {
      TREE->XL[L] = D->X;
      TREE->PL[L] = D->I;
      TREE->PC[L] = P;
    }
    /// frontend for tree linking
    /// \param[in] tree  tree to be linked
    /// \param[in] nmin  only make cells from boxes with at least @a nmin dots
    void Link(const OctTree*tree, node_index nmin) const
    {
      const_cast<const OctTree*&>(TREE) = tree;
      CF = 1;
      LF = 0;
      TREE->PA[0] = 0;
      const_cast<node_index&>(NCRIT) = nmin;
      const_cast<depth_type&>(DEPTH) = LinkCell(P0,0,0);
      const_cast<node_index&>(NCELL) = CF;
    }
    //@}
    /// minimum number of cells needed to link
    size_t NCells(node_index nmin) const
    {
      if(nmin<=2)
	return BM.N_used();
      size_t Nb=0;
      for(typename block_alloc<Box>::iterator i=BM.begin(); i; ++i)
	if(i->NUM >= nmin) ++Nb;
      return Nb? Nb:1;
    }
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
      for(octant_type i=0; i!=Nsub; ++i) out<<" OCT["<<i<<']';
      out<<" S NSP         Rad                 C\n";
    }
    /// dump box data
    /// \param[out] ns  counter for single-parent boxes
    void Dump(const Box*B, std::ostream&out, node_index&ns)
    {
      out<<" B"<<setfill('0')<<setw(5)<<BM.number_of_element(B)
	 <<' ' <<setfill(' ')<<setw(5)<<B->NUM;
      if(B->DOT)
	out<<" D"<<setfill('0')<<setw(5)<<int(B->DOT-D0);
      else
	out<<" nil   ";
      for(octant_type i=0; i!=Nsub; ++i) {
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
	 <<' '<<setw(8)<<B->CEN<<'\n';
    }
    /// dump tree, recursive
    void Dump(const Box*B, std::ostream&outd, std::ostream&outb, node_index&ns)
    {
      Dump(B,outb,ns);
      if(B->DOT)
	for(Dot*Di=B->DOT; Di; Di=Di->Next)
	  Dump(Di,outd);
      else {
	for(octant_type i=0; i!=Nsub; ++i)
	  if(B->OCT[i] && !B->MarkedAsBox(i))
	    Dump(pDOT(B->OCT[i]),outd);
	for(octant_type i=0; i!=Nsub; ++i)
	  if(B->OCT[i] && B->MarkedAsBox(i))
	    Dump(pBOX(B->OCT[i]),outd,outb,ns);
      }
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
  //
  template<int D, typename Real>
  void BoxDotTree<D,Real>::LinkLeafs(const Box*P, node_index C) const
  {
    if(P->DOT)
      for(Dot*Di = P->DOT; Di; Di=Di->Next)
	LinkLeaf(LF++, Di, C);
    else
      for(octant_type i=0; i!=Nsub; ++i) if(P->OCT[i]) {
	  if(P->MarkedAsBox(i)) LinkLeafs(      cpBOX(P->OCT[i]), C);
	  else                  LinkLeaf (LF++, cpDOT(P->OCT[i]), C);
	}
  }
  //
  template<int D, typename Real>
  typename BoxDotTree<D,Real>::depth_type
  BoxDotTree<D,Real>::LinkCell(const Box*P, node_index C, octant_type o) const
  {
    // 1 if single-parent replace by non-single-parent descendant
    EnsureNonSingleParent(P);
    // 2 copy some data, set octant
    TREE->L0[C] = LF;
    TREE->OC[C] = o;
    TREE->LE[C] = P->LEV;
    TREE->XC[C] = P->CEN;
    TREE->NM[C] = P->NUM;
    // 3 final box: link dots and return
    if(P->DOT) {
      TREE->NL[C] = P->NUM;
      TREE->CF[C] = 0;
      TREE->NC[C] = 0;
      for(Dot*Di = P->DOT; Di; Di=Di->Next)
	LinkLeaf(LF++, Di, C);
      return 1;
    }
    // 4 loop octants: link leaf kids, count cells
    node_index ncell=0;
    for(octant_type i=0; i!=Nsub; ++i)
      if(P->OCT[i]) {
	if(P->MarkedAsDot(i))
	  LinkLeaf (LF++, cpDOT(P->OCT[i]), C);
	else if(cpBOX(P->OCT[i])->NUM < NCRIT)
	  LinkLeafs(      cpBOX(P->OCT[i]), C);
	else
	  ++ncell;
      }
    // 5 set number of leaf and cell kids, first daughter cell, if any
    TREE->NL[C] = LF - TREE->L0[C];
    TREE->NC[C] = ncell;
    TREE->CF[C] = ncell? CF : 0;
    if(ncell==0) return 1;
    // 6 link daughter cells
    depth_type dp = 1;
    node_index Ci = CF;
    CF += ncell;
    for(octant_type i=0; i!=Nsub; ++i)
      if(P->OCT[i] && P->MarkedAsBox(i) && cpBOX(P->OCT[i])->NUM >= NCRIT) {
	TREE->PA[Ci] = C;
	update_max(dp, LinkCell(cpBOX(P->OCT[i]), Ci++, i));
      }
    return ++dp;
  }
  //
  template<int Dim, typename Real>
  BoxDotTree<Dim,Real>::BoxDotTree(char build, node_index Ndot,
				   const Initialiser*Init, local_count Nmax,
				   const OctTree*Tree) WDutils_THROWING
    : Base  ( build,Ndot,Init,Tree ),
      NCRIT ( Nmax ),
      BM    ( build=='r'? Tree->Ncells() : 1+NDOT/4 )
  {
    // 1 find min, max and average position
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
	    WDutils_THROW("OctalTree: Dot %d (particle %d): "
			  "invalid X=%g,%g\n",
			  int(Di-D0), Di->I, Di->X[0],Di->X[1])
	    :
	    WDutils_THROW("OctalTree: Dot %d (particle %d): "
			  "invalid X=%g,%g,%g\n",
			  int(Di-D0), Di->I, Di->X[0],Di->X[1],Di->X[2]);
      WDutils_THROW("OctalTree: unidentified invalid position(s)\n");
    }
    // 2 set root box, RA[]
    Xave   /= Real(NDOT);
    P0      = NewBox();
    P0->CEN = Helper<Dim,Real>::Integer(Xave);
    P0->LEV = 0;
    RA[0]   = Helper<Dim,Real>::RootRadius(P0->CEN,Xmin,Xmax);
    for(unsigned l=0; l!=MAXD; ++l)
      RA[l+1] = Real(0.5)*RA[l];
    // 3 add dots
    ND = 0;
    for(Di=D0; Di!=DN; ++Di,++ND)
      AddDot(P0,Di);
#ifdef TESTING
    std::ofstream dumpD("dots.dat"), dumpB("boxs.dat");
    Dump(dumpD, dumpB);
#endif
  }
  //
#undef   pDOT
#undef  cpDOT
#undef   pBOX
#undef  cpBOX
}
//
namespace WDutils {
  template<int D, typename Real>
  void OctalTree<D,Real>::Allocate()
  {
    unsigned need =
      NLEAF * (sizeof(point) + sizeof(particle_key) + sizeof(node_index)) +
      NCELL * (3*sizeof(uint8) + sizeof(uint16) + 4*sizeof(node_index) +
	       sizeof(point));
    if((need > NALLOC) || (need+need < NALLOC)) {
      if(ALLOC) delete16(ALLOC);
      ALLOC  = new16<char>(need);
      NALLOC = need;
    }
    char* A = ALLOC;
    XL = reinterpret_cast<point*>       (A); A += NLEAF * sizeof(point);
    PL = reinterpret_cast<particle_key*>(A); A += NLEAF * sizeof(particle_key);
    PC = reinterpret_cast<node_index *> (A); A += NLEAF * sizeof(node_index);
    LE = reinterpret_cast<uint8*>       (A); A += NCELL * sizeof(uint8);
    OC = reinterpret_cast<uint8*>       (A); A += NCELL * sizeof(uint8);
    XC = reinterpret_cast<point*>       (A); A += NCELL * sizeof(point);
    L0 = reinterpret_cast<node_index*>  (A); A += NCELL * sizeof(node_index);
    NL = reinterpret_cast<uint16*>      (A); A += NCELL * sizeof(uint16);
    NM = reinterpret_cast<node_index*>  (A); A += NCELL * sizeof(node_index);
    CF = reinterpret_cast<node_index*>  (A); A += NCELL * sizeof(node_index);
    NC = reinterpret_cast<uint8*>       (A); A += NCELL * sizeof(uint8);
    PA = reinterpret_cast<node_index*>  (A); A += NCELL * sizeof(node_index);
  }
  //
  template<int D, typename Real>
  void OctalTree<D,Real>::build(char building, node_index n,
				const Initialiser*init, const OctalTree*tree)
    WDutils_THROWING
  {
    if(0==init)
      WDutils_THROW("OctalTree<%d,%s>::build(): init=0\n",D,nameof(Real));
    BoxDotTree<D,Real> BDT(building,n,init,NMAX,tree);
    NLEAF = BDT.NDOT;
    NCELL = BDT.NCells(NMIN);
    Allocate();
    BDT.Link(this,NMIN);
    DEPTH = BDT.DEPTH;
    NCELL = BDT.NCELL;
    std::memcpy(RAD,BDT.RA,(MaximumDepth+1)*sizeof(Real));
  }
}
//
// Wdutils::TreeAccess<OctTree>
//
namespace WDutils {
  template<typename OctTree>
  typename TreeAccess<OctTree>::Cell
  TreeAccess<OctTree>::SmallestContainingCell(point const&x) const
  {
    Cell c=Root();
    if(!::contains(centre(c),radius(c),x) )
      return InvalidCell();
    for(;;) {
      if(Ncells(c)==0) return c;
      uint8 o=::octant(centre(c),x);
      Cell cc=BeginCells(c), ce=EndCells(c);
      while(cc!=ce && o!=octant(cc)) ++cc;
      if(cc==ce) return c;
      c=cc;
    }
  }
}
//
// Wdutils::NeighbourLoop<OctTree>
//
namespace WDutils {
  template<typename OctTree> inline
  bool NeighbourLoop<OctTree>::Outside(Cell c) const
  { return ::outside(centre(c),radius(c),X,Q); }
  //
  template<typename OctTree> inline
  bool NeighbourLoop<OctTree>::Inside(Cell c) const
  { return ::inside(centre(c),radius(c),X,Q); }
  //
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
    typedef typename NeighbourFinder<OctTree>::Real Real;
    typedef typename NeighbourFinder<OctTree>::node_index node_index;
    //
    Neighbour<OctTree> *LIST;    ///< neighbour list
    const   node_index  K;       ///< size of list
    mutable node_index  I;       ///< index of current element
    /// ctor: take data
    Lister(Neighbour<OctTree>*list, node_index size)
      : LIST(list), K(size), I(0) {}
    /// process: add neighbour to list
    void process(Leaf l, Real q) const
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
      Real q = dist_sq(X,position(l));
      if(q<Q) PROC->process(l,q);
    }
  }
  //
  template<typename OctTree>
  typename NeighbourFinder<OctTree>::node_index
  NeighbourFinder<OctTree>::Find(Leaf l, Real q, Neighbour<OctTree>*nb,
				 node_index m)
  {
    Lister<OctTree> LL(nb,m);
    PROC =&LL;
    Q    = q;
    X    = position(l);
    C    = Parent(l);
    Base::Process();
    return LL.I;
  }
  //
  template<typename OctTree>
  typename NeighbourFinder<OctTree>::node_index
  NeighbourFinder<OctTree>::Find(point const&x, Real q, Neighbour<OctTree>*nb,
				 node_index m)
  {
    Lister<OctTree> LL(nb,m);
    PROC =&LL;
    Q    = q;
    X    = x;
    C    = SmallestContainingCell(X);
    Base::Process();
    return LL.I;
  }
  //
  template<typename OctTree>
  void NeighbourFinder<OctTree>::Process(Leaf l, Real q, const Processor*p)
    WDutils_THROWING
  {
    if(0==p) WDutils_THROW("NeighbourFinder::Process(): p=0\n");
    PROC = p;
    Q    = q;
    X    = position(l);
    C    = Parent(l);
    Base::Process();
  }
  //
  template<typename OctTree>
  void NeighbourFinder<OctTree>::Process(point const&x, Real q,
					 const Processor*p) WDutils_THROWING
  {
    if(0==p) WDutils_THROW("NeighbourFinder::Process(): p=0\n");
    PROC = p;
    Q    = q;
    X    = x;
    C    = SmallestContainingCell(X);
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
    size_t nx  = N16 * sizeof(Real);
    size_t na  = N16 + OctTree::Dim * nx;
    if(na > NALLOC || 2*na < 3*NALLOC) {
      if(ALLOC) free16(ALLOC);
      const_cast<size_t&>(NALLOC) = na;
      const_cast<char* &>(ALLOC)  = NALLOC? new16<char>(na) : 0;
    }
    char* A = ALLOC;
    const_cast<Real*&>(XX)=reinterpret_cast<Real*>(A);  A+=nx;
    const_cast<Real*&>(YY)=reinterpret_cast<Real*>(A);  A+=nx;
    if(OctTree::Dim == 3) {
      const_cast<Real*&>(ZZ)=reinterpret_cast<Real*>(A);  A+=nx;
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
    __m128 HQ = _mm_set1_ps(Q);
    __m128 X0 = _mm_set1_ps(X[0]);
    __m128 Y0 = _mm_set1_ps(X[1]);
    SSE::Traits<Real>::vector QQ;
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
    __m128 HQ = _mm_set1_ps(Q);
    __m128 X0 = _mm_set1_ps(X[0]);
    __m128 Y0 = _mm_set1_ps(X[1]);
    __m128 Z0 = _mm_set1_ps(X[2]);
    SSE::Traits<Real>::vector QQ;
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
    __m128d HQ = _mm_set1_pd(Q);
    __m128d X0 = _mm_set1_pd(X[0]);
    __m128d Y0 = _mm_set1_pd(X[1]);
    SSE::Traits<Real>::vector QQ;
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
    __m128d HQ = _mm_set1_pd(Q);
    __m128d X0 = _mm_set1_pd(X[0]);
    __m128d Y0 = _mm_set1_pd(X[1]);
    __m128d Z0 = _mm_set1_pd(X[2]);
    SSE::Traits<Real>::vector QQ;
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
  Find(Leaf l, Real q, Neighbour<OctTree>*nb, node_index m)
  {
    Q    = q;
    X    = position(l);
    C    = Parent(l);
    while(C != NLoop::Root() && IsValid(C) && Number(Parent(C)) < NLoop::NDIR )
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
  Find(point const&x, Real q, Neighbour<OctTree>*nb, node_index m)
  {
    Q    = q;
    X    = x;
    C    = SmallestContainingCell(NLoop::X);
    while(C != NLoop::Root() && IsValid(C) && Number(Parent(C)) < NLoop::NDIR )
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
    typename TreeAccess<OctTree>::Real Q;
    typename TreeAccess<OctTree>::Cell C;
    void set(typename TreeAccess<OctTree>::Real q,
	     typename TreeAccess<OctTree>::Cell c)
    { Q=q; C=c; }
//     bool operator<(CellQ const&x) const
//     { return Q < x.Q; }
    bool operator>(CellQ const&x) const 
    { return Q > x.Q; }
//     bool operator<(typename TreeAccess<OctTree>::Real q) const
//     { return Q < q; }
//     bool operator>(typename TreeAccess<OctTree>::Real q) const
//     { return Q > q; }
  };
}
//
namespace WDutils {
  template<typename OctTree> inline
  typename NearestNeighbourFinder<OctTree>::Real
  NearestNeighbourFinder<OctTree>::OutsideDistSq(Cell c) const
  { return outside_dist_sq(centre(c),radius(c),X); }
  //
  template<typename OctTree> inline
  bool NearestNeighbourFinder<OctTree>::Outside(Cell c) const
  { return outside(centre(c),radius(c),X,LIST->Q); }
  //
  template<typename OctTree> inline
  bool NearestNeighbourFinder<OctTree>::Inside(Cell c) const
  { return inside(centre(c),radius(c),X,LIST->Q); }
  //
  template<typename OctTree> inline
  typename NearestNeighbourFinder<OctTree>::node_index
  NearestNeighbourFinder<OctTree>::Ndir() const
  { return M<=0? NDIR : max(static_cast<node_index>(M),NDIR); }
  //
  template<typename OctTree> inline
  void NearestNeighbourFinder<OctTree>::AddLeaf(Leaf l) const
  {
    // testing after adding the contributions to q from each dimension does
    // make the code run more slowly
    Real q = dist_sq(X,position(l));
    if(LIST->Q > q) {
      LIST->Q = q;
      LIST->L = l;
      MaxHeap::after_top_replace(LIST,K);
      --M;
      ++NIAC;
    }
  }
  //
  template<typename OctTree>
  void NearestNeighbourFinder<OctTree>::AddCell(Cell Ci, node_index cC) const
  {
    if(cC==0 && Number(Ci) <= Ndir()) {
      // direct loop
      LoopAllLeafs(Ci,l) AddLeaf(l);
    } else {
      // process leaf kids
      if(Nleafkids(Ci)) {
	LoopLeafKids(Ci,l) AddLeaf(l);
      }
      // process cell kids
      if(Ncells(Ci)>cC+1) {
	// more than one sub-cell: process in order of increasing distance
	CellQ<OctTree> Z[Base::Nsub];
	int J(0);
	LoopCellKids(Ci,c)
	  if(c!=C) Z[J++].set(OutsideDistSq(c),c);
	MinHeap::build(Z,J);
	while(J && LIST->Q > Z->Q) {
	  AddCell(Z->C);
	  Z[0] = Z[J-1];
	  MinHeap::after_top_replace(Z,--J);
	}
      } else if(Ncells(Ci)>cC) {
	// only 1 sub-cell to process
	if(cC) { LoopCellKids(Ci,c) if(c!=C && !Outside(c)) AddCell(c);}
	else   { LoopCellKids(Ci,c) if(!Outside(c)) AddCell(c); }
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
//     // TEST
//     {
//       Real  r;
//       point c,x;
//       for(;;) {
// 	std::cerr<<" c="; std::cin>>c;
// 	std::cerr<<" r="; std::cin>>r;
// 	std::cerr<<" x="; std::cin>>x;
// 	Real q(0),D;
// 	D = abs(c[0]-x[0]); if(D>r) q+=square(D-r);
// 	D = abs(c[1]-x[1]); if(D>r) q+=square(D-r);
// 	D = abs(c[2]-x[2]); if(D>r) q+=square(D-r);
// 	std::cerr<<" outside_dist_sq(c,r,x) = "<<q<<'\n';
// 	std::cerr<<" SSE version:             "<<outside_dist_sq(c,r,x)<<'\n';
//       }
//     }
//     // TSET
    NIAC = 0;
    M    = K;
    Real Q = 12*square(Base::RootRadius());
    for(node_index k=0; k!=K; ++k)
      LIST[k].Q = Q;
    for(Cell P=C; IsValid(P) && !Inside(C); C=P,P=Parent(C))
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
