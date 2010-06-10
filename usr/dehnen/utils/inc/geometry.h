// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/geometry.h
///
/// \brief  simple geometrical algorithms for 2D and 3D.
///
/// \author Walter Dehnen
///
/// \date   2010
///
/// \version 07-06-2010 WD  tested: octant,contains,dist_sq,outside,inside
/// \version 08-06-2010 WD  new template layout, struct Algorithm<al,sse>
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2010 Walter Dehnen
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
#ifndef WDutils_included_tupel_h
#  include <tupel.h>
#endif

namespace WDutils {
  namespace Geometry {
    ///
    /// a cubic box
    ///
    template<int Dim, typename real>
    struct cube {
      tupel<Dim,real> X;         ///< centre of cube
      real            H;         ///< half side length of cube
    };
    ///
    /// a sphere
    ///
    template<int Dim, typename real>
    struct sphere {
      tupel<Dim,real> X;         ///< centre of sphere
      real            Q;         ///< radius-squared of sphere
    };
    ///
    /// simple geometry algorithms implemented using explicit SSE
    ///
    /// \note If template parameter @a aligned_to_16_bytes is true, we assume
    ///       that all data are 16-byte aligned, which often results in
    ///       considerable speed-up (up to factor 3 when @a use_sse is
    ///       true). If data are unaligned but both @a aligned_to_16_bytes and
    ///       @a use_sse are true, a run-time error will result.
    ///
    /// \note The SSE algorithms are at least as fast as the non-SSE ones, and
    ///       usually faster in particular if data are 16-byte aligned. It is
    ///       recommended to use the default setting for @a use_sse (if @a
    ///       use_sse is set, but SSE is unavailable, the non-SSE instructions
    ///       are used automatically). The explicit non-SSE versions are used
    ///       mainly for validation of the SSE versions.
    ///
    /// \note Unless otherwise stated, the algorithms (static methods below)
    ///       are implemented only for 2D and 3D and for single (float) and
    ///       double precision. Using other parameters causes a compile-time
    ///       error.
    ///
    template<bool aligned_to_16_bytes, bool use_sse = true>
    struct Algorithms
    {
      /// copy a cube
      /// \param[in]  in  cube to copy
      /// \param[out] out cube copied
      template<int Dim, typename real> static inline
      void copy(cube<Dim,real> const&in, cube<Dim,real> &out);
      /// copy a sphere
      /// \param[in]  in  sphere to copy
      /// \param[out] out sphere copied
      template<int Dim, typename real> static inline
      void copy(sphere<Dim,real> const&in, sphere<Dim,real> &out);
      /// move centre position to octant
      /// \param[in,out] c cube
      /// \param[in]     i octant
      /// \param[in]     r amount to move by (= radius of shrunk cube)
      template<int Dim, typename real> static inline
      void move_to_octant(tupel<Dim,real>&c, int i, real r);
      /// shrink cube to its octant
      /// \param[in,out] c cube
      /// \param[in]     i octant
      template<int Dim, typename real> static inline
      void shrink(cube<Dim,real>&c, int i)
      { c.H *= real(0.5); move_to_octant(c.X,i,c.H); }
      /// octant of a point w.r.t. a centre.
      /// \param[in] c  centre
      /// \param[in] x  point
      /// \return integer: ith bit equals x[i]>c[i]; bits @a D and beyond are 0.
      template<int Dim, typename real> static inline
      int octant(tupel<Dim,real> const&c, tupel<Dim,real> const&x);
      /// octant of a point w.r.t. a cube's centre
      /// \param[in] c  cube
      /// \param[in] x  point
      /// \return integer: ith bit equals x[i]>c[i]; bits @a D and beyond are 0.
      template<int Dim, typename real> static inline
      int octant(cube<Dim,real> const&c, tupel<Dim,real> const&x)
      { return octant(c.X,x); }
      /// does a cubic box contain a given position
      /// \param[in] c    cube
      /// \param[in] x    position
      /// \return is @a x[i] in [c.X[i]-c.R[i], c.X[i]+c.R[i]) for i=0..D-1?
      template<int Dim, typename real> static inline
      bool contains(cube<Dim,real> const&c, tupel<Dim,real> const&x);
      /// distance^2 between two points
      /// \param[in] x    point
      /// \param[in] y    point
      /// \return |x-y|^2
      template<int Dim, typename real> static inline
      real dist_sq(tupel<Dim,real> const&x, tupel<Dim,real> const&y);
      /// distance^2 from given point to the nearest point on a cube
      /// \param[in] c    cube
      /// \param[in] x    point
      /// \return squared distance of @a x to @a c; zero if @a x is inside @a c.
      template<int Dim, typename real> static inline
      real dist_sq(cube<Dim,real> const&c, tupel<Dim,real> const&x);
      /// is a sphere completely outside of a cube?
      /// \param[in] c    cube
      /// \param[in] s    sphere
      /// \return is sphere outside cube?
      /// \note Equivalent to, but on average faster than, 
      ///       \code s.Q < outside_dist_sq(c,s.X) \endcode
      template<int Dim, typename real> static inline
      bool outside(cube<Dim,real> const&c, sphere<Dim,real> const&s);
      /// is a sphere completely inside of a cube?
      /// \param[in] c    cube
      /// \param[in] s    sphere
      /// \return is sphere completely inside cube?
      template<int Dim, typename real> static inline
      bool inside(cube<Dim,real> const&c, sphere<Dim,real> const&s);
    };// struct WDutils::Geometry::Algorithms<aligned,sse>

    ///
    /// a search sphere
    ///
    /// \note If SSE instructions for @a __X are available, we implement in
    ///       geometry_inl.h specialisation which exploit them and, to this
    ///       end, have a different data representation. Therefore, nothing
    ///       must be assumed about the private member layout.
    ///
    /// \note Only @a __D = 2 or 3 and @a __X = float or double are possible.
    template<int __D, typename __X> struct SearchSphere
    {
      static const int Dim = __D;  ///< number of spatial dimensions
      typedef __X      real;       ///< floating point type
      /// default ctor
      SearchSphere() {}
      /// ctor from centre and radius^2
      /// \param[in] x  centre of sphere
      /// \param[in] q  radius^2 of sphere
      SearchSphere(tupel<Dim,real> const&x, real q) { reset(x,q); }
      /// ctor from sphere
      /// \param[in] s  sphere, need not be aligned
      SearchSphere(sphere<Dim,real> const&s) { reset(s); }
      /// reset centre and radius^2
      /// \param[in] x  centre of sphere
      /// \param[in] q  radius^2 of sphere
      void reset(tupel<Dim,real> const&x, real q) { S.X=x; S.Q=q; }
      /// reset centre and radius^2
      /// \param[in] x  centre of sphere
      /// \param[in] q  radius^2 of sphere
      void reset(sphere<Dim,real> const&s) { S.X=s.X; S.Q=s.Q; }
      /// reset just the centre
      void reset(tupel<Dim,real> const&x) { S.X=x; }
      /// reset just the radius^2
      void reset(real q) { S.Q = q; }
      /// radius^2 of search sphere
      real const&RadSq() const { return S.Q; }
      /// \name geometric relations with cubic box
      //@{
      /// distance^2 from cube to centre of sphere (zero if inside cube)
      /// \param[in]  c  cubic box
      /// \note If @a c is 16-byte aligned, use @a aligned == true, otherwise
      ///       @a aligned == false. In case SSE instructions are supported for
      ///       @a real, we specialise this routine in geometry_inl.h such that
      ///       alignement gives slightly faster execution.
      template<bool aligned>
      real dist_sq(cube<Dim,real> const&c) const
      { return Algorithms<aligned>::dist_sq(c,S.X); }
      /// is search sphere outside of a cubic box?
      /// \param[in]  c  cubic box
      /// \note If @a c is 16-byte aligned, use @a aligned == true, otherwise
      ///       @a aligned == false. In case SSE instructions are supported for
      ///       @a real, we specialise this routine in geometry_inl.h such that
      ///       alignement gives slightly faster execution.
      template<bool aligned>
      bool outside(cube<Dim,real> const&c) const
      { return Algorithms<aligned>::outside(c,S); }
      /// is search sphere inside of a cubic box?
      /// \param[in]  c  cubic box
      /// \note If @a c is 16-byte aligned, use @a aligned == true, otherwise
      ///       @a aligned == false. In case SSE instructions are supported for
      ///       @a real, we specialise this routine in geometry_inl.h such that
      ///       alignement gives slightly faster execution.
      template<bool aligned>
      bool inside(cube<Dim,real> const&c) const
      { return Algorithms<aligned>::inside(c,S); }
      //@}
      /// \name geometric relations with position
      //@{
      /// distance^2 from centre of sphere to some point
      /// \param[in]  x  position
      /// \note If @a x is 16-byte aligned, use @a aligned == true, otherwise
      ///       @a aligned == false. In case SSE instructions are supported for
      ///       @a real, we specialise this routine in geometry_inl.h such that
      ///       alignement gives slightly faster execution.
      template<bool aligned>
      real dist_sq(tupel<Dim,real> const&x) const
      { return WDutils::dist_sq(x,S.X); }
      /// is a position contained within the search sphere?
      /// \param[in]  x  position
      /// \note If @a x is 16-byte aligned, use @a aligned == true, otherwise
      ///       @a aligned == false. In case SSE instructions are supported for
      ///       @a real, we specialise this routine in geometry_inl.h such that
      ///       alignement gives slightly faster execution.
      template<bool aligned>
      bool contains(tupel<Dim,real> const&x) const
      { return dist_sq(x) < S.Q; }
      //@}
    private:
      WDutilsStaticAssert( ( __D == 2 || __D == 3 )               &&
			   meta::TypeInfo<__X>::is_floating_point    );
      WDutils__align16 sphere<Dim,real> S;   ///< search sphere data
    };// struct SearchSphere
  } // namespace WDutils::Geometry
} // namespace WDutils
// inline implementations
#ifndef WDutils_included_geometry_inl_h
#  include <geometry_inl.h>
#endif
//
#endif // WDutils_included_geometry_h
