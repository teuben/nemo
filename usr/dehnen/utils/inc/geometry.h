// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/geometry.h
///
/// \brief  simple geometrical algorithms for 2D and 3D.
///
/// \author Walter Dehnen
///
/// \date   2010,2013
///
/// \version 07-06-2010 WD  tested: octant,contains,dist_sq,outside,inside
/// \version 08-06-2010 WD  new template layout, struct Algorithms<al,sse>
/// \version 10-06-2010 WD  struct SearchSphere<Dim,real>
/// \version 11-06-2010 WD  cuboid<> supported in Algorithms<>, SearchSphere<>
/// \version 11-01-2012 WD  contains(cuboid,cuboid)
/// \version 18-06-2012 WD  adapted for tupel.h -> vector.h (using C++11)
/// \version 06-11-2012 WD  adapted for AVX usage; using noexcept decls
/// \version 09-04-2013 WD  added contains_open() 
/// \version 10-04-2013 WD  octant() & contains() consistent; always_inline
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2010,2013 Walter Dehnen
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
#ifndef WDutils_included_geometry_h
#define WDutils_included_geometry_h

#ifndef WDutils_included_sse_h
#  include <sse.h>
#endif
#if __cplusplus >= 201103L
# ifndef WDutils_included_vector_h
#  include <vector.h>
# endif
#else
# ifndef WDutils_included_tupel_h
#  include <tupel.h>
# endif
# define noexcept
#endif

namespace WDutils {
#if defined(__clang__) || (defined(__GNUC__) && !defined(__INTEL_COMPILER))
#  define always_inline inline __attribute__((__always_inline__))
#else
#  define always_inline inline
#endif
  ///
  /// contains support for geometry algorithms in 2D and 3D
  ///
  namespace Geometry {
#if __cplusplus < 201103L
#  define  GeoVec tupel
#else
#  define  GeoVec vector
#endif
    ///
    /// a cubic box
    ///
    template<int Dim, typename real>
    struct cube {
      typedef GeoVec<Dim,real> point;
      point  X;                            ///< centre of cube
      real   H;                            ///< half side length of cube
    };
    ///
    /// a sphere
    ///
    template<int Dim, typename real>
    struct sphere {
      typedef GeoVec<Dim,real> point;
      point  X;                            ///< centre of cube
      real   Q;                            ///< radius-squared of sphere
    };
    //
    template<int Dim, typename real>
    struct alignment {
      static_assert(Dim==2 || Dim==3,
		    "only supporting Dim=2,3");
      static_assert(std::is_same<real,float >::value ||
		    std::is_same<real,double>::value,
		    "only supporting float and double");
      static const size_t value = 
#ifdef __AVX__
	Dim==3? 4*sizeof(real) :
#endif
	16;
    };
    ///
    /// alignment-padded types
    ///
    template<int Dim, typename real>
    struct aligned 
    {
      static const size_t align = alignment<Dim,real>::value;
      /// a point padded to multiple of value
      typedef SSE::Extend< GeoVec<Dim,real>, align > point;
      /// a cube padded to multiple of value
      typedef SSE::Extend< Geometry::cube<Dim,real>, align> cube;
      /// a sphere padded to multiple of value
      typedef SSE::Extend< Geometry::sphere<Dim,real>, align> sphere;
      /// a pair of points padded to multiple of value
      struct point_pair { point X,Y; };
    };
    //  special case Dim=2, real=float
    template<>
    struct aligned<2,float>
    {
      typedef float real;
      static const size_t align = alignment<2,real>::value;
      typedef SSE::Extend< GeoVec<2,real>, align > point;
      typedef SSE::Extend< Geometry::cube<2,real>, align> cube;
      typedef SSE::Extend< Geometry::sphere<2,real>, align> sphere;
      struct point_pair { GeoVec<2,real> X,Y; };
    };
    //
    template<int Dim, typename real>
    struct PointPair : aligned<Dim,real>::point_pair {};
    ///
    /// cuboid: a rectangular box
    ///
    /// essentially a pair of vectors interpreted as geometrical centre C
    /// and offset H in each dimension to the side.
    ///
    /// \note Alternatively, one could code the cuboid as the positions of
    ///       the lower left and upper right corners C-H and C+H. However,
    ///       the geometric algorithms are simpler to code with (C,H).
    ///
    template<int Dim, typename real>
    struct cuboid : PointPair<Dim,real>
    {
      typedef PointPair<Dim,real> base;
      using base::X;
      using base::Y;
      typedef GeoVec<Dim,real> point;
      /// centre of cuboid
      always_inline point const&centre() const { return X; }
      /// centre of cuboid
      always_inline point      &centre()       { return X; }
      /// size of cuboid
      always_inline point const&size  () const { return Y; }
      /// size of cuboid
      always_inline point      &size  ()       { return Y; }
    };
    ///
    /// simple geometry algorithms implemented using explicit SSE
    ///
    /// \note If template parameter @a aligned is true, we assume that all
    ///       data are aligned to alignment<dim,type>::value bytes, which
    ///       often results in considerable speed-up. If data are unaligned
    ///       but both @a aligned and @a use_sse are true, a run-time error
    ///       will result.
    ///
    /// \note The SSE algorithms are at least as fast as the non-SSE ones, and
    ///       usually faster in particular if data are aligned. It is
    ///       recommended to use the default setting for @a use_sse (if @a
    ///       use_sse is set, but SSE is unavailable, the non-SSE instructions
    ///       are used automatically). The explicit non-SSE versions are used
    ///       mainly for validation of the SSE versions.
    ///
    /// \note Unless otherwise stated, the algorithms (static methods below) are
    ///       implemented only for 2D and 3D and for single (float) and double
    ///       precision. Using other parameters causes a compile-time error.
    ///
    template<bool aligned, bool use_sse = true>
    struct Algorithms
    {
      ///
      /// is a certain datum appropriately aligned for usage?
      ///
      template<int Dim, typename real> static always_inline 
      bool is_aligned(cube<Dim,real> const&c)
      {
	return aligned
	  && WDutils::is_aligned<alignment<Dim,real>::value>
	  (static_cast<const real*>(c.X));
      }
      //
      template<int Dim, typename real> static always_inline
      bool is_aligned(cuboid<Dim,real> const&c)
      {
	return aligned
	  && WDutils::is_aligned<alignment<Dim,real>::value>
	  (static_cast<const real*>(c.X));
      }
      //
      template<int Dim, typename real> static always_inline
      bool is_aligned(sphere<Dim,real> const&s)
      {
	return aligned
	  && WDutils::is_aligned<alignment<Dim,real>::value>
	  (static_cast<const real*>(s.X));
      }
      //
      template<int Dim, typename real> static always_inline
      bool is_aligned(GeoVec<Dim,real> const&x)
      {
	return aligned
	  && WDutils::is_aligned<alignment<Dim,real>::value>
	  (static_cast<const real*>(x));
      }
      ///
      /// copy a cube
      ///
      /// \param[in]  in  cube to copy
      /// \param[out] out cube copied
      template<int Dim, typename real> static always_inline
      void copy(cube<Dim,real> const&in, cube<Dim,real> &out) noexcept;
      ///
      /// copy a sphere
      ///
      /// \param[in]  in  sphere to copy
      /// \param[out] out sphere copied
      template<int Dim, typename real> static always_inline
      void copy(sphere<Dim,real> const&in, sphere<Dim,real> &out) noexcept;
      ///
      /// move centre position to octant
      ///
      /// \param[in,out] c cube
      /// \param[in]     i octant
      /// \param[in]     r amount to move by (= radius of shrunk cube)
      template<int Dim, typename real> static always_inline
      void move_to_octant(GeoVec<Dim,real>&c, int i, real r) noexcept;
      ///
      /// shrink cube to its octant
      ///
      /// \param[in,out] c cube
      /// \param[in]     i octant
      template<int Dim, typename real> static always_inline
      void shrink(cube<Dim,real>&c, int i) noexcept
      { c.H *= real(0.5); move_to_octant(c.X,i,c.H); }
      ///
      /// octant of a point w.r.t. a centre.
      ///
      /// \param[in] c  centre
      /// \param[in] x  point
      /// \return integer: ith bit equals x[i]>c[i]; bits @a D and beyond are 0.
      template<int Dim, typename real> static always_inline
      int octant(GeoVec<Dim,real> const&c, GeoVec<Dim,real> const&x) noexcept;
      ///
      /// octant of a point w.r.t. a cube's centre
      ///
      /// \param[in] c  cube
      /// \param[in] x  point
      /// \return integer: ith bit equals x[i]>=c[i]
      template<int Dim, typename real> static always_inline
      int octant(cube<Dim,real> const&c, GeoVec<Dim,real> const&x) noexcept
      { return octant(c.X,x); }
      ///
      /// distance^2 between two points
      ///
      /// \param[in] x    point
      /// \param[in] y    point
      /// \return |x-y|^2
      template<int Dim, typename real> static always_inline
      real dist_sq(GeoVec<Dim,real> const&x, GeoVec<Dim,real> const&y) noexcept;
      /// does a cubic box contain a given position (half-open intervals)
      /// \param[in] c    cube
      /// \param[in] x    position
      /// \return is @a x[i] in [c.X[i]-c.R[i], c.X[i]+c.R[i]) for i=0..D-1?
      /// \note consistent with @a octant() and @a shrink() in the sense that
      ///       if contains(c,x)==true, then contains(cc,x)==true where cc=c;
      ///       shrink(cc,octant(c,x)); (assuming exact cube representation,
      ///       e.g. c.X=0, c.H=1), even if x=c.X.
      template<int Dim, typename real> static always_inline
      bool contains(cube<Dim,real> const&c, GeoVec<Dim,real> const&x) noexcept;
      /// does a cubic box contain a given position (open intervals)
      /// \param[in] c    cube
      /// \param[in] x    position
      /// \return is @a x[i] in [c.X[i]-c.R[i], c.X[i]+c.R[i]] for i=0..D-1?
      template<int Dim, typename real> static always_inline
      bool contains_open(cube<Dim,real> const&c, GeoVec<Dim,real> const&x)
	  noexcept;
      ///
      /// distance^2 from given point to the nearest point on a cube
      ///
      /// \param[in] c    cube
      /// \param[in] x    point
      /// \return squared distance of @a x to @a c; zero if @a x is inside @a c.
      template<int Dim, typename real> static always_inline
      real dist_sq(cube<Dim,real> const&c, GeoVec<Dim,real> const&x) noexcept;
      ///
      /// is a sphere completely outside of a cube?
      ///
      /// \param[in] c    cube
      /// \param[in] s    sphere
      /// \return is sphere outside cube?
      /// \note Equivalent to, but on average faster than, 
      ///       \code s.Q < dist_sq(c,s.X) \endcode
      template<int Dim, typename real> static always_inline
      bool outside(cube<Dim,real> const&c, sphere<Dim,real> const&s) noexcept;
      ///
      /// is a sphere completely inside of a cube?
      ///
      /// \param[in] c    cube
      /// \param[in] s    sphere
      /// \return is sphere completely inside cube?
      template<int Dim, typename real> static always_inline
      bool inside(cube<Dim,real> const&c, sphere<Dim,real> const&s) noexcept;
      ///
      /// does a cuboid contain a given position (half-open intervals)
      ///
      /// \param[in] c    cuboid
      /// \param[in] x    position
      /// \return is @a x[i] in [c.X[i]-c.R[i], c.X[i]+c.R[i]) for i=0..D-1?
      template<int Dim, typename real> static always_inline
      bool contains(cuboid<Dim,real> const&c,
		    GeoVec<Dim,real> const&x) noexcept;
      ///
      /// does a cuboid contain a given position (open intervals)
      ///
      /// \param[in] c    cuboid
      /// \param[in] x    position
      /// \return is @a x[i] in [c.X[i]-c.R[i], c.X[i]+c.R[i]] for i=0..D-1?
      template<int Dim, typename real> static always_inline
      bool contains_open(cuboid<Dim,real> const&c,
			 GeoVec<Dim,real> const&x) noexcept;
      ///
      /// does a cuboid contain another cuboid
      ///
      /// \param[in] c    cuboid
      /// \param[in] b    inner cuboid
      /// \return is @a b inside @a c (common boundary accepted)?
      template<int Dim, typename real> static always_inline
      bool contains(cuboid<Dim,real> const&c,
		    cuboid<Dim,real> const&b) noexcept;
      ///
      /// distance^2 from given point to the nearest point on a cuboid
      ///
      /// \param[in] c    cuboid
      /// \param[in] x    point
      /// \return squared distance of @a x to @a c; zero if @a x is inside @a c.
      template<int Dim, typename real> static always_inline
      real dist_sq(cuboid<Dim,real>  const&c,
		   GeoVec<Dim,real> const&x) noexcept;
      ///
      /// is a sphere completely outside of a cuboid?
      ///
      /// \param[in] c    cuboid
      /// \param[in] s    sphere
      /// \return is sphere outside cube?
      /// \note Equivalent to, but on average faster than, 
      ///       \code s.Q < dist_sq(c,s.X) \endcode
      template<int Dim, typename real> static always_inline
      bool outside(cuboid<Dim,real> const&c,
		   sphere<Dim,real> const&s) noexcept;
      ///
      /// is a sphere completely inside of a cuboid?
      ///
      /// \param[in] c    cuboid
      /// \param[in] s    sphere
      /// \return is sphere completely inside cube?
      template<int Dim, typename real> static always_inline
      bool inside(cuboid<Dim,real> const&c,
		  sphere<Dim,real> const&s) noexcept;
      ///
      /// converting (Xmin,Xmax) to (centre,size)
      ///
      /// \note centre,size = 0.5 (Xmax +/- Xmin)
      /// \param[in,out] p  input: (Xmin,Xmax) output: cuboid=(centre,size)
      template<int Dim, typename real> static always_inline
      void convert2cuboid(PointPair<Dim,real> &p) noexcept;
    };// struct WDutils::Geometry::Algorithms<aligned,sse>

    ///
    /// apropriately aligned search-sphere data
    ///
    template<int _tD, typename _tX> struct SearchSphereData
    {
    protected: WDutils__align16 sphere<_tD,_tX> S;
    };
#ifdef __AVX__
    template<> struct SearchSphereData<3,double>
    {
    protected: WDutils__align32 sphere<3,double> S;
    };
#endif
    ///
    /// a search sphere
    ///
    /// \note If SSE instructions for @a _tX are available, we implement in
    ///       geometry_inl.h specialisation which exploit them and, to this
    ///       end, have a different data representation. Therefore, nothing
    ///       must be assumed about the private member layout.
    ///
    /// \note In case of SSE support, 16-byte aligned data allow for faster
    ///       implementation. Therefore, two versions are provided for most
    ///       methods: one assuming 16-byte alignment of data provided and
    ///       another not assuming any aligment. To distinguish them, these
    ///       second versions take an additional integer argument, which is
    ///       not used. That is, the default is to assume 16-byte alignement.
    ///
    ///       Providing data not aligned to 16 bytes to methods assuming
    ///       16-byte alignement will cause a run-time error.
    ///
    /// \note Only @a _tD = 2 or 3 and @a _tX = float or double are possible.
    template<int _tD, typename _tX> struct SearchSphere
    : private AlignedDatum<
#ifdef __AVX__
      _tD==3? 4*sizeof(_tX) :
#endif
      16, sphere<_tD,_tX> >
    {
      static const int         Dim = _tD; ///< number of spatial dimensions
      typedef _tX              real;      ///< floating point type
#if __cplusplus < 201103L
      typedef tupel<Dim,real>  point;     ///< position type
#else
      typedef vector<Dim,real> point;     ///< position type
#endif
      ///
      /// default ctor
      ///
#if __cplusplus < 201103L
      SearchSphere() {}
#else
      SearchSphere() = default;
#endif
      ///
      /// ctor from sphere
      ///
      /// \param[in] s  sphere, need not be aligned
      always_inline
      SearchSphere(sphere<Dim,real> const&s, int) noexcept
      { reset(s,0); }
      ///
      /// ctor from sphere, assuming 16-byte alignment
      ///
      /// \param[in] s  sphere, 16-byte aligned
      explicit always_inline  
      SearchSphere(sphere<Dim,real> const&s) noexcept
      { reset(s); }
      ///
      /// ctor from centre position and radius^2
      ///
      /// \param[in] x  centre of sphere
      /// \param[in] q  radius^2 of sphere
      always_inline
      SearchSphere(point const&x, real q, int) noexcept
      { reset(x,q,0); }
      ///
      /// ctor from centre position and radius^2, assuming 16-byte alignment
      ///
      /// \param[in] x  centre of sphere, 16-byte aligned
      /// \param[in] q  radius^2 of sphere
      always_inline
      SearchSphere(point const&x, real q) noexcept
      { reset(x,q); }
      ///
      /// reset centre and radius^2
      ///
      /// \param[in] x  new centre of sphere
      /// \param[in] q  new radius^2 of sphere
      always_inline
      void reset(point const&x, real q, int) noexcept
      { Datum.X=x; Datum.Q=q; }
      ///
      /// reset centre and radius^2, assuming 16-byte alignment
      ///
      /// \param[in] x  new centre of sphere, 16-byte aligned
      /// \param[in] q  new radius^2 of sphere
      always_inline
      void reset(point const&x, real q) noexcept
      { Datum.X=x; Datum.Q=q; }
      ///
      /// reset centre and radius^2
      ///
      /// \param[in] s  new sphere
      always_inline
      void reset(sphere<Dim,real> const&s, int) noexcept
      { Datum.X=s.X; Datum.Q=s.Q; }
      ///
      /// reset centre and radius^2, assuming 16-byte alignment
      ///
      /// \param[in] s  new sphere, 16-byte aligned
      always_inline
      void reset(sphere<Dim,real> const&s) noexcept
      { Datum.X=s.X; Datum.Q=s.Q; }
      ///
      /// reset just the radius^2
      ///
      /// \param[in] q  new radius^2 of sphere
      always_inline
      void reset(real q) noexcept
      { Datum.Q = q; }
      ///
      /// reset just the centre
      ///
      /// \param[in] x  new centre of sphere
      always_inline
      void reset(point const&x, int) noexcept
      { Datum.X=x; }
      ///
      /// reset just the centre, assuming 16-byte alignment
      ///
      /// \param[in] x  new centre of sphere
      always_inline
      void reset(point const&x) noexcept
      { Datum.X=x; }
      ///
      /// centre of search sphere
      ///
      always_inline
      point const&Centre() const noexcept
      { return Datum.X; }
      ///
      /// centre of search sphere
      ///
      always_inline
      point const&centre() const noexcept
      { return Datum.X; }
      ///
      /// ith co-ordinate of centre of search sphere
      ///
      always_inline
      real const&Centre(int i) const noexcept
      { return Datum.X[i]; }
      ///
      /// radius^2 of search sphere
      ///
      always_inline
      real const&RadSq() const noexcept
      { return Datum.Q; }
      ///
      /// radius^2 of search sphere
      ///
      always_inline
      real const&radius_squared() const noexcept
      { return Datum.Q; }
      ///
      /// \name geometric relations with position
      ///
      //@{
      ///
      /// distance^2 from centre of sphere to some point
      ///
      /// \param[in]  x  position
      always_inline
      real dist_sq(point const&x, int) const noexcept
      { return Algorithms<0>::dist_sq(Datum.X,x); }
      ///
      /// distance^2 from centre of sphere to some point
      ///
      /// \param[in]  x  position, 16-byte aligned
      always_inline
      real dist_sq(point const&x) const noexcept
      { return Algorithms<1>::dist_sq(Datum.X,x); }
      ///
      /// is a position contained within the search sphere?
      ///
      /// \param[in]  x  position
      always_inline
      bool contains(point const&x, int) const noexcept
      { return dist_sq(x,0) < Datum.Q; }
      ///
      /// is a position contained within the search sphere?
      ///
      /// \param[in]  x  position, 16-byte aligned
      always_inline
      bool contains(point const&x) const noexcept
      { return dist_sq(x) < Datum.Q; }
      //@}
      ///
      /// \name geometric relations with cubic box
      ///
      //@{
      ///
      /// distance^2 from cube to centre of sphere (zero if inside cube)
      ///
      /// \param[in]  c  cubic box
      always_inline
      real dist_sq(cube<Dim,real> const&c, int) const noexcept
      { return Algorithms<0>::dist_sq(c,Datum.X); }
      ///
      /// distance^2 from cube to centre of sphere (zero if inside cube)
      ///
      /// \param[in]  c  cubic box, 16-byte aligned
      always_inline
      real dist_sq(cube<Dim,real> const&c) const noexcept
      { return Algorithms<1>::dist_sq(c,Datum.X); }
      ///
      /// is search sphere outside of a cubic box?
      ///
      /// \param[in]  c  cubic box
      always_inline
      bool outside(cube<Dim,real> const&c, int) const noexcept
      { return Algorithms<0>::outside(c,Datum); }
      ///
      /// is search sphere outside of a cubic box?
      ///
      /// \param[in]  c  cubic box, 16-byte aligned
      always_inline
      bool outside(cube<Dim,real> const&c) const noexcept
      { return Algorithms<1>::outside(c,Datum); }
      ///
      /// is search sphere inside of a cubic box?
      ///
      /// \param[in]  c  cubic box
      always_inline
      bool inside(cube<Dim,real> const&c, int) const noexcept
      { return Algorithms<0>::inside(c,Datum); }
      ///
      /// is search sphere inside of a cubic box?
      ///
      /// \param[in]  c  cubic box, 16-byte aligned
      always_inline
      bool inside(cube<Dim,real> const&c) const noexcept
      { return Algorithms<1>::inside(c,Datum); }
      //@}
      ///
      /// \name geometric relations with rectangular box
      ///
      //@{
      ///
      /// distance^2 from cuboid to centre of sphere (zero if inside cube)
      ///
      /// \param[in]  c  rectangular box
      always_inline
      real dist_sq(cuboid<Dim,real> const&c, int) const noexcept
      { return Algorithms<0>::dist_sq(c,Datum.X); }
      ///
      /// distance^2 from cuboid to centre of sphere (zero if inside cuboid)
      ///
      /// \param[in]  c  rectangular box, 16-byte aligned
      always_inline
      real dist_sq(cuboid<Dim,real> const&c) const noexcept
      { return Algorithms<1>::dist_sq(c,Datum.X); }
      ///
      /// is search sphere outside of a cubic box?
      ///
      /// \param[in]  c  rectangular box
      always_inline
      bool outside(cuboid<Dim,real> const&c, int) const noexcept
      { return Algorithms<0>::outside(c,Datum); }
      ///
      /// is search sphere outside of a cubic box?
      ///
      /// \param[in]  c  rectangular box, 16-byte aligned
      always_inline
      bool outside(cuboid<Dim,real> const&c) const noexcept
      { return Algorithms<1>::outside(c,Datum); }
      ///
      /// is search sphere inside of a cubic box?
      ///
      /// \param[in]  c  rectangular box
      always_inline
      bool inside(cuboid<Dim,real> const&c, int) const noexcept
      { return Algorithms<0>::inside(c,Datum); }
      ///
      /// is search sphere inside of a cubic box?
      ///
      /// \param[in]  c  rectangular box, 16-byte aligned
      always_inline
      bool inside(cuboid<Dim,real> const&c) const noexcept
      { return Algorithms<1>::inside(c,Datum); }
      //@}
    private:
      WDutilsStaticAssert( ( _tD == 2 || _tD == 3 )               &&
			   is_floating_point<_tX>::value   );
      using AlignedDatum<
#ifdef __AVX__
	_tD==3? 4*sizeof(_tX) :
#endif
	16, sphere<_tD,_tX> >::Datum;
    };// struct SearchSphere
  } // namespace WDutils::Geometry
#define Geometry_TRAITS(TYPE,NAME)					\
  template<> struct traits<TYPE<2,float> >				\
  { static const char  *name () { return NAME "<2,float>"; } };		\
  template<> struct traits<TYPE<3,float> >				\
  { static const char  *name () { return NAME "<2,float>"; } };		\
  template<> struct traits<TYPE<2,double> >				\
  { static const char  *name () { return NAME "<2,double>"; } };	\
  template<> struct traits<TYPE<3,double> >				\
  { static const char  *name () { return NAME "<2,double>"; } };

  Geometry_TRAITS(Geometry::cube,"Geometry::cube")
  Geometry_TRAITS(Geometry::sphere,"Geometry::sphere")
  Geometry_TRAITS(Geometry::PointPair,"Geometry::PointPair")
  Geometry_TRAITS(Geometry::cuboid,"Geometry::cuboid")
  Geometry_TRAITS(Geometry::SearchSphere,"Geometry::SearchSphere")
#undef Geometry_TRAITS
} // namespace WDutils
#undef noexcept
#undef always_inline
// inline implementations
#ifndef WDutils_included_geometry_tcc
#  include <geometry.tcc>
#endif
#undef GeoVec
//
#endif // WDutils_included_geometry_h
