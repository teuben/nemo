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
// Note. I tried to accelerate the tree building by making Dot and Box 16-byte
//       aligned (and having size a multiple of 16) plus using aligned memory
//       loading (into SSE) to find the octant. This trial was not successful.
//
namespace {
  using std::setw;
  using std::setfill;
  using namespace WDutils;
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
    static int octant(point const&cen, point const&pos)
    {
#ifdef __SSE__
      return 3&_mm_movemask_ps(_mm_cmplt_ps(_mm_loadu_ps(cen),
					    _mm_loadu_ps(pos)));
#else
      int oct(0);
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
    static void ShrinkToOctant(point&cen, int i, Real rad)
    {
      if(i&1) cen[0] += rad;  else  cen[0] -= rad;
      if(i&2) cen[1] += rad;  else  cen[1] -= rad;
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
    static int octant(point const&cen, point const&pos)
    {
      int oct(0);
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
    static void ShrinkToOctant(point&cen, int i, Real rad)
    {
      if(i&1) cen[0] += rad;  else  cen[0] -= rad;
      if(i&2) cen[1] += rad;  else  cen[1] -= rad;
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
    static int octant(point const&cen, point const&pos)
    {
#ifdef __SSE__
      return 7 & _mm_movemask_ps(_mm_cmplt_ps(_mm_loadu_ps(cen),
					      _mm_loadu_ps(pos)));
#else
      int oct(0);
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
    static void ShrinkToOctant(point&cen, int i, Real rad)
    {
      if(i&1) cen[0] += rad;  else  cen[0] -= rad;
      if(i&2) cen[1] += rad;  else  cen[1] -= rad;
      if(i&4) cen[2] += rad;  else  cen[2] -= rad;
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
    static int octant(point const&cen, point const&pos)
    {
      int oct(0);
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
    static void ShrinkToOctant(point&cen, int i, Real rad)
    {
      if(i&1) cen[0] += rad;  else  cen[0] -= rad;
      if(i&2) cen[1] += rad;  else  cen[1] -= rad;
      if(i&4) cen[2] += rad;  else  cen[2] -= rad;
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
  int octant(tupel<D,Real> const&cen, tupel<D,Real> const&pos)
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
  /// move @a cen by @a rad in direction of octant @a i
  /// \note used for setting centre of daughter box from that of parent box
  template<int D, typename Real> inline
  void ShrinkToOctant(tupel<D,Real>&cen, int i, Real rad)
  { return Helper<D,Real>::ShrinkToOctant(cen,i,rad); }
  /// represents a particle in the BoxDotTree
  template<int Dim, typename Real>
  struct Dot
  {
    typedef typename OctalTree<Dim,Real>::particle_key particle_key;
    typedef typename OctalTree<Dim,Real>::point point;
    point        X;               ///< position
    particle_key I;               ///< identifier of associated particle
    mutable Dot *Next;            ///< next dot in a linked list
  };
  /// represents a cubic cell in the BoxDotTree
  template<int Dim, typename Real>
  struct Box
  {
    const static int Nsub = 1<<Dim;
    typedef typename OctalTree<Dim,Real>::node_index node_index;
    typedef typename meta::__IWORDS<Nsub>::integer_u ndl_type;
    typedef typename OctalTree<Dim,Real>::point point;
    //
    point           CEN;          ///< centre position
    void*           OCT[Nsub];    ///< octants
    union {
      ndl_type      NDl;
      /// if 0 == NDL[i]          the octant is empty                        \n
      /// if 0 <  NDL[i] <= NMAX  the octant holds a list of NDL[] leafs     \n
      /// if      NDL[i]  > NMAX  the octant holds a box
      uint8         NDL[Nsub];
    };
    uint8           NBX;          ///< number of daughter boxes  <= Nsub
    uint8           NOC;          ///< number of octant occupied <= NBX
    uint8           LEV;          ///< tree level of box
    uint8           PEA;          ///< Peano key, not currently used
    node_index      NUM;          ///< total # dots in box
  };
}
//
namespace WDutils {
  template<int D, typename R> struct traits< ::Dot<D,R> >
  { static const char*name () { return message("Dot<%d,%s>",D,nameof(R)); } };
  template<int D, typename R> struct traits< ::Box<D,R> >
  { static const char*name () { return message("Box<%d,%s>",D,nameof(R)); } };
}
//
namespace {
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
    { if(D0) delete16(D0); }
  };
  //
#define  pDOT(d) static_cast<Dot*>((d))
#define cpDOT(d) static_cast<const Dot*>((d))
#define  pBOX(d) static_cast<Box*>((d))
#define cpBOX(d) static_cast<const Box*>((d))
  //
  template<int D, typename Real>
  DotInitialiser<D,Real>::DotInitialiser(char build, node_index Ndot,
					 const Initialiser*Init,
					 const OctTree*Tree) WDutils_THROWING
    : NDOT(Ndot), D0(NDOT? new16<Dot>(NDOT) : 0), DN(D0+NDOT)
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
      for(Di=List; Di; Di=pDOT(Di->Next))
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
  template<int Dim, typename Real>
  struct BoxDotTree : public DotInitialiser<Dim,Real>
  {
    const static int Nsub = 1<<Dim; ///< number of octants per cell
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
    const uint8         NMAX;          ///< maximum # particles/octant
    const uint8         NMIN;          ///< minimum # particles/cell
    const bool          AVSPC;         ///< avoid single-parent cells
    const node_index    NMAX1;         ///< NMAX + 1
    depth_type          DEPTH;         ///< depth of linked tree
    block_alloc<Box>    BM;            ///< allocator for boxes
    Real                RA[MAXD+1];    ///< array with radius(level)  
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
      P->CEN = B->CEN;
      ShrinkToOctant(P->CEN,i,RA[P->LEV]);
      return P;
    }
    /// adds dot @a Di to octant @a b of box @a P
    void AddDotToOctant(Box*P, Dot*Di, int b)
    {
      if(P->NDL[b] == 0) {
	P->NOC ++;
	Di->Next = 0;
      } else
	Di->Next = pDOT(P->OCT[b]);
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
	int b = octant(P->CEN,Di->X);
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
	      Dn= Di->Next;
	      b = octant(P->CEN,Di->X);
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
    /// \param[in] Tree  old OctalTree
    BoxDotTree(char build, node_index Ndot, const Initialiser*Init,
	       depth_type nmax, depth_type nmin, bool avspc,
	       const OctTree*Tree=0)
    WDutils_THROWING WD_HOT;
    //@}
    /// \name methods used in linking to OctalTree
    //@{
    /// tree linking: make a cell from an octant of a parent box
    /// \param[in] P  parent box
    /// \param[in] C  index of current cell to be linked
    /// \param[in] o  octant of cell in @a P
    void LinkOctantCell(const Box*P, node_index C, int i) const
    {
      TREE->L0[C] = LF;
      TREE->OC[C] = i;
      TREE->LE[C] = P->LEV + 1;
      TREE->XC[C] = P->CEN;
      ShrinkToOctant(TREE->XC[C], i, RA[TREE->LE[C]]);
      TREE->NM[C] = P->NDL[i];
      TREE->NL[C] = P->NDL[i];
      TREE->NC[C] = 0;
      TREE->CF[C] = 0;
      for(const Dot*Di=cpDOT(P->OCT[i]); Di; Di=Di->Next)
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
	out<<" D"<<setfill('0')<<setw(5)<<int(D->Next-D0);
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
      out<<' '<<setfill(' ')<<setw(8)<<RA[B->LEV]
	 <<' '<<setw(8)<<B->CEN<<'\n';
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
  //
  template<int D, typename Real>
  typename BoxDotTree<D,Real>::depth_type
  BoxDotTree<D,Real>::LinkCell(const Box*P, node_index C, int o) const
  {
    // 1 if single-parent replace by non-single-parent descendant
    if(AVSPC) EnsureNonSingleParent(P);
    // 2 copy some data, set octant
    TREE->L0[C] = LF;
    TREE->OC[C] = o;
    TREE->LE[C] = P->LEV;
    TREE->XC[C] = P->CEN;
    TREE->NM[C] = P->NUM;
    // 3 loop octants: link leaf kids (from octants with < NMIN), count cells
    depth_type tmp = 0;
    for(int i=0; i!=Nsub; ++i)
      if(P->NDL[i] >= NMIN)
	++tmp;
      else if(P->NDL[i])
	for(const Dot*Di=cpDOT(P->OCT[i]); Di; Di=Di->Next)
	  LinkLeaf(LF++,Di,C);
    // 4 set number of leaf and cell kids, first daughter cell, if any
    TREE->NL[C] = LF - TREE->L0[C];
    TREE->NC[C] = tmp;
    TREE->CF[C] = tmp? CF : 0;
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
  template<int Dim, typename Real>
  BoxDotTree<Dim,Real>::BoxDotTree(char build, node_index Ndot,
				   const Initialiser*Init,
				   depth_type nmax, depth_type nmin,
				   bool avspc, const OctTree*Tree)
    WDutils_THROWING
    : Base  ( build,Ndot,Init,Tree ),
      NMAX  ( nmax ),
      NMIN  ( nmin ),
      AVSPC ( avspc ),
      NMAX1 ( nmax + 1 ),
      BM    ( NBoxes<Dim>(Ndot,nmax) ),
      P0    ( BM.new_element() ),
      NCELL ( 1 )
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
	    WDutils_THROW("OctalTree<%d,%s>: particle %d: "
			  "invalid X=%g,%g\n",Dim,nameof(Real),
			  Di->I, Di->X[0],Di->X[1])
	    :
	    WDutils_THROW("OctalTree<%d,%s>: particle %d: "
			  "invalid X=%g,%g,%g\n",Dim,nameof(Real),
			  Di->I, Di->X[0],Di->X[1],Di->X[2]);
      WDutils_THROW("OctalTree<%d,%s>: unidentified invalid position(s)\n",
		    Dim,nameof(Real));
    }
    Xave   /= Real(NDOT);
    // 2 set (empty) root box, RA[]
    P0->NDl = 0;
    P0->NBX = 0;
    P0->NOC = 0;
    P0->CEN = Helper<Dim,Real>::Integer(Xave);
    P0->LEV = 0;
    P0->NUM = 0;
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
      NCELL * (sizeof(depth_type) + 2*sizeof(int) +
	       sizeof(local_count) + 4*sizeof(node_index) + sizeof(point));
    if((need > NALLOC) || (3*need < 2*NALLOC)) {
      if(ALLOC) delete16(ALLOC);
      ALLOC  = new16<char>(need);
      NALLOC = need;
    }
    char* A = ALLOC;
    XL = reinterpret_cast<point*>       (A); A += NLEAF * sizeof(point);
    PL = reinterpret_cast<particle_key*>(A); A += NLEAF * sizeof(particle_key);
    PC = reinterpret_cast<node_index*>  (A); A += NLEAF * sizeof(node_index);
    LE = reinterpret_cast<depth_type*>  (A); A += NCELL * sizeof(depth_type);
    OC = reinterpret_cast<octant_type*> (A); A += NCELL * sizeof(octant_type);
    XC = reinterpret_cast<point*>       (A); A += NCELL * sizeof(point);
    L0 = reinterpret_cast<node_index*>  (A); A += NCELL * sizeof(node_index);
    NL = reinterpret_cast<local_count*> (A); A += NCELL * sizeof(local_count);
    NM = reinterpret_cast<node_index*>  (A); A += NCELL * sizeof(node_index);
    CF = reinterpret_cast<node_index*>  (A); A += NCELL * sizeof(node_index);
    NC = reinterpret_cast<octant_type*> (A); A += NCELL * sizeof(octant_type);
    PA = reinterpret_cast<node_index*>  (A); A += NCELL * sizeof(node_index);
  }
  //
  template<int D, typename Real>
  void OctalTree<D,Real>::build(char building, node_index n,
				const Initialiser*init, const OctalTree*tree)
    WDutils_THROWING
  {
#undef GIVE_TIMING
#ifdef GIVE_TIMING
    clock_t cpu0, cpu1;
    cpu0 = clock();
#endif
    BoxDotTree<D,Real> BDT(building,n,init,NMAX,NMIN,AVSPC,tree);
#ifdef GIVE_TIMING
    cpu1 = clock();
    std::cerr<<" OctalTree::build(): BoxDotTree::BoxDotTree took "
	     <<double(cpu1 - cpu0)/double(CLOCKS_PER_SEC)<<" sec\n";
#endif
    NLEAF = BDT.NDOT;
    NCELL = BDT.NCell();
    Allocate();
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
