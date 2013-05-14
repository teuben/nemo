// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/geometry.tcc
///
/// \brief  implements inline methods for geometry.h
///
/// \author Walter Dehnen
///
/// \date   2010,2012
///
////////////////////////////////////////////////////////////////////////////////
///
/// \version Jun-2010 WD  first tested version
/// \version Jun-2012 WD  using vector instead of tupel if C++11
/// \version Nov-2012 WD  enable AVX for 3D double precision version 
///                       re-written using generic code from sse_vec.h
///                       some speed-up by avoiding branching
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2010,2012 Walter Dehnen
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
#ifndef WDutils_included_geometry_tcc
#define WDutils_included_geometry_tcc

#ifndef WDutils_included_geometry_h
#  include <geometry.h>
#endif

#if __cplusplus < 201103L
#  define noexcept
#endif

//
namespace WDutils {
  namespace Geometry {
    //
    // struct Algorithms
    //
    // NOTE that we cannot simply write a generic implementation for
    //      Algorithms<> and then provide specialisations for SSE and
    //      alignement, because it is not permittable in C++ to specialise
    //      member function templates --- only class templates, such as
    //      AlgorithmsHelper<>.
    //
    namespace details {
#if defined(__clang__) || (defined(__GNUC__) && !defined(__INTEL_COMPILER))
#  define always_inline __attribute__((__always_inline__))
#else
#  define always_inline
#endif
      //
      // 0  generic non-SSE code
      //
      template<typename real>
      class AlgorithmsHelperBase
      {
	friend struct Algorithms<0,0>;
	friend struct Algorithms<0,1>;
	friend struct Algorithms<1,0>;
	friend struct Algorithms<1,1>;
      protected:
	// is @a x in the interval [@a c-r, @a c+r) ?
	static bool always_inline in_interval(real c, real r, real x) noexcept
	{
	  return x<c? c<=x+r : x<c+r;
	}
	// is @a x in the interval [@a c-r, @a c+r] ?
	static bool always_inline in_open_interval(real c, real r, real x)
	  noexcept
	{
	  return x<c? c<=x+r : x<=c+r;
	}
	// is the interval [a-b, a+b] in the interval [c-r, c+r]?
	static bool always_inline in_interval(real c, real r, real a, real b)
	  noexcept
	{
	  return std::abs(c-a) <= r-b;
	}
	// copy(cube,cube)
	template<int Dim>
	static void always_inline copy(cube<Dim,real> const&in,
				       cube<Dim,real> &out) noexcept
	{
	  out.X = in.X;
	  out.H = in.H;
	}
	// copy(sphere,sphere)
	template<int Dim>
	static void always_inline copy(sphere<Dim,real> const&in,
				       sphere<Dim,real> &out) noexcept
	{
	  out.X = in.X;
	  out.Q = in.Q;
	}
	// convert2cuboid(PointPair)
	template<int Dim>
	static void always_inline convert2cuboid(PointPair<Dim,real> &c)
	  noexcept
	{
	  GeoVec<Dim,real> C = real(0.5)*(c.X+c.Y);
	  GeoVec<Dim,real> H = real(0.5)*(c.Y-c.X);
	  static_cast<GeoVec<Dim,real>&>(c.X) = C;
	  static_cast<GeoVec<Dim,real>&>(c.Y) = H;
	}
      };// class AlgorithmsHelperBase<>
      //
      // 1  generic cases: not using SSE, ignoring alignment
      //
      template<int Dim, typename real, bool aligned, bool sse>
      class AlgorithmsHelper;
      //
      // 1.1  2D non-SSE case
      //
      template<typename real, bool aligned, bool sse>
      class AlgorithmsHelper<2,real,aligned,sse> :
	protected AlgorithmsHelperBase<real>
      {
	friend struct Algorithms<0,0>;
	friend struct Algorithms<0,1>;
	friend struct Algorithms<1,0>;
	friend struct Algorithms<1,1>;
      protected:
	using AlgorithmsHelperBase<real>::copy;
	using AlgorithmsHelperBase<real>::in_interval;
	using AlgorithmsHelperBase<real>::in_open_interval;
	using AlgorithmsHelperBase<real>::convert2cuboid;
	//
	typedef GeoVec<2,real> vec2;
	typedef cube  <2,real> cub2;
	typedef cuboid<2,real> box2;
	typedef sphere<2,real> sph2;
	// octant(point,point)
	static int always_inline octant(vec2 const&c, vec2 const&x) noexcept
	{
	  int oct(0);
	  if(x[0] >= c[0]) oct |= 1;
	  if(x[1] >= c[1]) oct |= 2;
	  return oct;
	}
	// contains(cube,point)
	static bool always_inline contains(cub2 const&c, vec2 const&x) noexcept
	{
	  return in_interval(c.X[0],c.H,x[0])
	    &&   in_interval(c.X[1],c.H,x[1]);
	}
	// contains_open(cube,point)
	static bool always_inline contains_open(cub2 const&c,
						vec2 const&x) noexcept
	{
	  return in_open_interval(c.X[0],c.H,x[0])
	    &&   in_open_interval(c.X[1],c.H,x[1]);
	}
	// contains(cuboid,point)
	static bool always_inline contains(box2 const&c, vec2 const&x) noexcept
	{
	  return in_interval(c.X[0],c.Y[0],x[0])
	    &&   in_interval(c.X[1],c.Y[1],x[1]);
	}
	// contains_open(cuboid,point)
	static bool always_inline contains_open(box2 const&c,
						vec2 const&x) noexcept
	{
	  return in_open_interval(c.X[0],c.Y[0],x[0])
	    &&   in_open_interval(c.X[1],c.Y[1],x[1]);
	}
	// contains(cuboid,cuboid)
	static bool always_inline contains(box2 const&c, box2 const&b) noexcept
	{
	  return in_interval(c.X[0],c.Y[0], b.X[0],b.Y[0])
	    &&   in_interval(c.X[1],c.Y[1], b.X[1],b.Y[1]);
	}
	// dist_sq(point,point)
	static real always_inline dist_sq(vec2 const&x, vec2 const&y) noexcept
	{
	  return WDutils::dist_sq(x,y);
	}
	// dist_sq(cube,point)
	static real always_inline dist_sq(cub2 const&c, vec2 const&x) noexcept
	{
	  real q(0),D;
	  D = abs(c.X[0]-x[0]); if(D>c.H) q+=square(D-c.H);
	  D = abs(c.X[1]-x[1]); if(D>c.H) q+=square(D-c.H);
	  return q;
	}
	// dist_sq(cuboid,point)
	static real always_inline dist_sq(box2 const&c, vec2 const&x) noexcept
	{
	  real q(0),D;
	  D = abs(c.X[0]-x[0])-c.Y[0]; if(D>0) q+=D*D;
	  D = abs(c.X[1]-x[1])-c.Y[1]; if(D>0) q+=D*D;
	  return q;
	}
	// inside(cube,sphere)
	static bool always_inline inside(cub2 const&c, sph2 const&s) noexcept
	{
	  real D;
	  D=abs(c.X[0]-s.X[0])-c.H; if(D>0 || s.Q>D*D) return false;
	  D=abs(c.X[1]-s.X[1])-c.H; if(D>0 || s.Q>D*D) return false;
	  return true;
	}
	// inside(cuboid,sphere)
	static bool always_inline inside(box2 const&c, sph2 const&s) noexcept
	{
	  real D;
	  D=abs(c.X[0]-s.X[0])-c.Y[0]; if(D>0 || s.Q>D*D) return false;
	  D=abs(c.X[1]-s.X[1])-c.Y[1]; if(D>0 || s.Q>D*D) return false;
	  return true;
	}
	// move to octant
	static void always_inline move_to_octant(vec2&c, int i, real r) noexcept
	{
	  if(i&1) c[0] += r;  else  c[0] -= r;
	  if(i&2) c[1] += r;  else  c[1] -= r;
	}
      };

      //
      // 1.2  3D non-SSE case
      //
      template<typename real, bool aligned, bool sse>
      class AlgorithmsHelper<3,real,aligned,sse> :
	protected AlgorithmsHelperBase<real>
      {
	friend struct Algorithms<0,0>;
	friend struct Algorithms<0,1>;
	friend struct Algorithms<1,0>;
	friend struct Algorithms<1,1>;
      protected:
	using AlgorithmsHelperBase<real>::copy;
	using AlgorithmsHelperBase<real>::in_interval;
	using AlgorithmsHelperBase<real>::in_open_interval;
	using AlgorithmsHelperBase<real>::convert2cuboid;
	//
	typedef GeoVec<3,real> vec3;
	typedef cube  <3,real> cub3;
	typedef cuboid<3,real> box3;
	typedef sphere<3,real> sph3;
	// octant(point,point)
	static int always_inline octant(vec3 const&c, vec3 const&x) noexcept
	{
	  int oct(0);
	  if(x[0] >= c[0]) oct |= 1;
	  if(x[1] >= c[1]) oct |= 2;
	  if(x[2] >= c[2]) oct |= 4;
	  return oct;
	}
	// contains(cube,point)
	static bool always_inline contains(cub3 const&c, vec3 const&x) noexcept
	{
	  return in_interval(c.X[0],c.H,x[0])
	    &&   in_interval(c.X[1],c.H,x[1])
	    &&   in_interval(c.X[2],c.H,x[2]);
	}
	// contains_open(cube,point)
	static bool contains_open(cub3 const&c, vec3 const&x) noexcept
	{
	  return in_open_interval(c.X[0],c.H,x[0])
	    &&   in_open_interval(c.X[1],c.H,x[1])
	    &&   in_open_interval(c.X[2],c.H,x[2]);
	}
	// contains(cuboid,point)
	static bool always_inline contains(box3 const&c, vec3 const&x) noexcept
	{
	  return in_interval(c.X[0],c.Y[0],x[0])
	    &&   in_interval(c.X[1],c.Y[1],x[1])
	    &&   in_interval(c.X[2],c.Y[2],x[2]);
	}
	// contains_open(cuboid,point)
	static bool always_inline contains_open(box3 const&c,
						vec3 const&x) noexcept
	{
	  return in_open_interval(c.X[0],c.Y[0],x[0])
	    &&   in_open_interval(c.X[1],c.Y[1],x[1])
	    &&   in_open_interval(c.X[2],c.Y[2],x[2]);
	}
	// contains(cuboid,cuboid)
	static bool always_inline contains(box3 const&c, box3 const&b) noexcept
	{
	  return in_interval(c.X[0],c.Y[0], b.X[0],b.Y[0])
	    &&   in_interval(c.X[1],c.Y[1], b.X[1],b.Y[1])
	    &&   in_interval(c.X[2],c.Y[2], b.X[2],b.Y[2]);
	}
	// dist_sq(point,point)
	static real always_inline dist_sq(vec3 const&x, vec3 const&y) noexcept
	{
	  return WDutils::dist_sq(x,y);
	}
	// dist_sq(cube,point)
	static real always_inline dist_sq(cub3 const&c, vec3 const&x) noexcept
	{
	  real q(0),D;
	  D = abs(c.X[0]-x[0]); if(D>c.H) q+=square(D-c.H);
	  D = abs(c.X[1]-x[1]); if(D>c.H) q+=square(D-c.H);
	  D = abs(c.X[2]-x[2]); if(D>c.H) q+=square(D-c.H);
	  return q;
	}
	// dist_sq(cuboid,point)
	static real always_inline dist_sq(box3 const&c, vec3 const&x) noexcept
	{
	  real q(0),D;
	  D = abs(c.X[0]-x[0])-c.Y[0]; if(D>0) q+=D*D;
	  D = abs(c.X[1]-x[1])-c.Y[1]; if(D>0) q+=D*D;
	  D = abs(c.X[2]-x[2])-c.Y[2]; if(D>0) q+=D*D;
	  return q;
	}
	// inside(cube,sphere)
	static bool always_inline inside(cub3 const&c, sph3 const&s) noexcept
	{
	  real D;
	  D=abs(c.X[0]-s.X[0])-c.H; if(D>0 || s.Q>D*D) return false;
	  D=abs(c.X[1]-s.X[1])-c.H; if(D>0 || s.Q>D*D) return false;
	  D=abs(c.X[2]-s.X[2])-c.H; if(D>0 || s.Q>D*D) return false;
	  return true;
	}
	// inside(cuboid,sphere)
	static bool always_inline inside(box3 const&c, sph3 const&s) noexcept
	{
	  real D;
	  D=abs(c.X[0]-s.X[0])-c.Y[0]; if(D>0 || s.Q>D*D) return false;
	  D=abs(c.X[1]-s.X[1])-c.Y[1]; if(D>0 || s.Q>D*D) return false;
	  D=abs(c.X[2]-s.X[2])-c.Y[2]; if(D>0 || s.Q>D*D) return false;
	  return true;
	}
	// move to octant
	static void always_inline move_to_octant(vec3&c, int i, real r) noexcept
	{
	  if(i&1) c[0] += r;  else  c[0] -= r;
	  if(i&2) c[1] += r;  else  c[1] -= r;
	  if(i&4) c[2] += r;  else  c[2] -= r;
	}
      };// class AlgorithmsHelper<3,real,aligned,sse>

#ifdef __SSE__

      using namespace WDutils::SSE;
      //
      // 2    specialisations using SSE and AVX
      //

      //
      // 2.1  SSE code for 2D with real=float
      //
      template<bool aligned> class AlgorithmsHelper<2,float,aligned,1> {
	friend struct Algorithms<aligned,1>;
	typedef GeoVec<2,float> vec2;
	typedef cube  <2,float> cub2;
	typedef cuboid<2,float> box2;
	typedef sphere<2,float> sph2;
	// load a fvec4
	static fvec4 always_inline load(vec2 const&p) noexcept
	{
	  return fvec4::template load_t<aligned>(p.data());
	}
	// store a fvec4
	static void always_inline store(fvec4 x, vec2&p) noexcept
	{
	  x.template store_t<aligned>(p.data());
	}
#if(0)
	// load a fvec4
	template<typename A>
	static fvec4 always_inline load(A const&p) noexcept
	{
	  return fvec4::template pack_t<aligned>(p);
	}
	// store a fvec4
	template<typename A>
	static void always_inline store(fvec4 x, A&p) noexcept
	{
	  x.template unpack_t<aligned>(p);
	}
#endif
	// octant(point,point)
	static int always_inline octant(vec2 const&c, vec2 const&x) noexcept
	{
	  return 3 & signbits(load(c) <= load(x));
	}
	// contains(cube,point)
	static bool always_inline contains(cub2 const&c, vec2 const&x) noexcept
	{
	  fvec4 C = load(c.X);
	  fvec4 H = single<2>(C);
	  fvec4 X = load(x);
	  return 3 == (3 & signbits((C-H <= X) & (X < C+H)));
	}
	// contains_open(cube,point)
	static bool always_inline contains_open(cub2 const&c,
						vec2 const&x) noexcept
	{
	  fvec4 C = load(c.X);
	  fvec4 H = single<2>(C);
	  fvec4 X = load(x);
	  return 3 == (3 & signbits((C-H <= X) & (X <= C+H)));
	}
	// contains(cuboid,point)
	static bool always_inline contains(box2 const&c, vec2 const&x) noexcept
	{
	  fvec4 C = load(c.X);
	  fvec4 H = movehl(C,C);
	  fvec4 X = load(x);
	  return 3 == (3 & signbits((C-H <= X) & (X < C+H)));
	}
	// contains_open(cuboid,point)
	static bool always_inline contains_open(box2 const&c,
						vec2 const&x) noexcept
	{
	  fvec4 C = load(c.X);
	  fvec4 H = movehl(C,C);
	  fvec4 X = load(x);
	  return 3 == (3 & signbits((C-H <= X) & (X <= C+H)));
	}
	// contains(cuboid,cuboid)
	// for unaligned access, non-SSE code is faster
	static bool always_inline contains(box2 const&c, box2 const&b) noexcept
	{
	  WDutilsStaticAssert(aligned);
	  fvec4 C = load(c.X);
	  fvec4 R = movehl(C,C);
	  fvec4 A = load(b.X);
	  fvec4 B = movehl(A,A);
	  C.make_diff(A);
	  R-= B;
	  return !(3&signbits(C>R));
	}
	// dist_sq(point,point)
	// for unaligned access, non-SSE code is faster
	static float always_inline dist_sq(vec2 const&x, vec2 const&y) noexcept
	{
	  WDutilsStaticAssert(aligned);
	  fvec4 X;
	  X = load(x);
	  X-= load(y);
	  X*= X;
	  return sum2(X);
	}
	// dist_sq(cube,point)
	static float always_inline dist_sq(cub2 const&c, vec2 const&x) noexcept
	{
	  fvec4 C = load(c.X);
	  fvec4 H = single<2>(C);
	  C.make_diff(load(x));
	  C-= H;
	  C = max(C,fvec4::zero());
	  C*= C;
	  return sum2(C);
	}
	// dist_sq(cuboid,point)
	static float always_inline dist_sq(box2 const&c, vec2 const&x) noexcept
	{
	  fvec4 C = load(c.X);
	  fvec4 H = movehl(C,C);
	  C.make_diff(load(x));
	  C-= H;
	  C = max(C,fvec4::zero());
	  C*= C;
	  return sum2(C);
	}
	// inside(cube,sphere)
	static bool always_inline inside(cub2 const&c, sph2 const&s) noexcept
	{
	  fvec4 C = load(c.X);
	  fvec4 H = single<2>(C);
	  fvec4 X = load(s.X);
	  fvec4 Q = single<2>(X);
	  X.make_diff(C);
	  C = H-X;
	  C*= C;
	  return ! ( 3 & signbits( (Q>C) | (X>H) ) );
	}
	// inside(cuboid,sphere)
	static bool always_inline inside(box2 const&c, sph2 const&s) noexcept
	{
	  fvec4 C = load(c.X);
	  fvec4 H = movehl(C,C);
	  fvec4 X = load(s.X);
	  fvec4 Q = single<2>(X);
	  X.make_diff(C);
	  C = H-X;
	  C*= C;
	  return ! ( 3 & signbits( (Q>C) | (X>H) ) );
	}
	// convert2cuboid(PointPair)
	static void always_inline convert2cuboid(PointPair<2,float> &c) noexcept
	{
	  WDutilsStaticAssert(aligned);
	  fvec4 XY = load(c.X);                                     // [X  , Y ]
	  fvec4 YX = shuffle<2,3,0,1>(XY,XY);                       // [Y  , X ]
	  YX ^= fvec4(fvec4::nil_mask_elem(),fvec4::nil_mask_elem(),// [Y  ,-X ]
		      fvec4::sgn_mask_elem(),fvec4::sgn_mask_elem());
	  XY += YX;                                                 // [X+Y,Y-X]
	  XY *= fvec4(0.5f);
	  store(XY,c.X);
	}
      };// class AlgorithmsHelper<2,float,aligned,1>

# ifdef __AVX__
      //
      // 2.2  SSE code for 3D with real=float
      //      AVX code for 3D with real=double
      //
      template<typename real, bool aligned>
      class AlgorithmsHelper<3,real,aligned,1>
# else
      //
      // 2.2  SSE code for 3D with real=float
      //
      template<bool aligned> class AlgorithmsHelper<3,float,aligned,1>
# endif
      {
	friend struct Algorithms<aligned,1>;
# ifndef __AVX__
	typedef float real;
# else
	WDutilsStaticAssert((is_same<real,float>::value ||
			     is_same<real,double>::value));
# endif
	static const bool real_is_float = is_same<real,float>::value;
	typedef GeoVec<3,real> vec3;
	typedef cube  <3,real> cub3;
	typedef cuboid<3,real> box3;
	typedef sphere<3,real> sph3;
	typedef packed<4,real> rvec4;
	// load a rvec4
	static rvec4 always_inline load(vec3 const&p) noexcept
	{
	  return rvec4::template load_t<aligned>(p.data());
	}
	// store a rvec4
	static void always_inline store(rvec4 x, vec3&p) noexcept
	{
	  x.template store_t<aligned>(p.data());
	}
#if(0)
	// load a rvec4
	template<typename A>
	static rvec4 always_inline load(A const&p) noexcept
	{
	  return rvec4::template pack_t<aligned>(p);
	}
	// store a rvec4
	template<typename A>
	static void always_inline store(rvec4 x, A&p) noexcept
	{
	  x.template unpack_t<aligned>(p);
	}
#endif
	// copy(cube,cube)
	static void always_inline copy(cub3 const&in, cub3 &out) noexcept
	{
	  store(load(in.X),out.X);
	}
	// copy(sphere,sphere)
	static void always_inline copy(sph3 const&in, sph3 &out) noexcept
	{
	  store(load(in.X),out.X);
	}
	// octant(point,point)
	static int always_inline octant(vec3 const&c, vec3 const&x) noexcept
	{
	  return 7 & signbits(load(c) <= load(x));
	}
	// contains(cube,point)
	static bool always_inline contains(cub3 const&c, vec3 const&x) noexcept
	{
	  rvec4 C = load(c.X);
	  rvec4 H = single<3>(C);
	  rvec4 X = load(x);
	  return 7 == (7 & signbits((C-H <= X) & (X < C+H)));
	}
	// contains_open(cube,point)
	static bool always_inline contains_open(cub3 const&c,
						vec3 const&x) noexcept
	{
	  rvec4 C = load(c.X);
	  rvec4 H = single<3>(C);
	  rvec4 X = load(x);
	  return 7 == (7 & signbits((C-H <= X) & (X <= C+H)));
	}
	// contains(cuboid,point)
	static bool always_inline contains(box3 const&c, vec3 const&x) noexcept
	{
	  rvec4 C = load(c.X);
	  rvec4 H = load(c.Y);
	  rvec4 X = load(x);
	  return 7 == (7 & signbits((C-H <= X) & (X < C+H)));
	}
	// contains_open(cuboid,point)
	static bool always_inline contains_open(box3 const&c,
						vec3 const&x) noexcept
	{
	  rvec4 C = load(c.X);
	  rvec4 H = load(c.Y);
	  rvec4 X = load(x);
	  return 7 == (7 & signbits((C-H <= X) & (X <= C+H)));
	}
	// contains(cuboid,cuboid)
	// for unaligned access, non-SSE code is faster
	static bool always_inline contains(box3 const&c, box3 const&b) noexcept
	{
	  WDutilsStaticAssert(!real_is_float || aligned);
	  rvec4 C = load(c.X);  C.make_diff(load(b.X));
	  rvec4 D = load(c.Y);  D-= load(b.Y);
	  return !(7&signbits( C > D ));
	}
	// dist_sq(point,point)
	// for unaligned access, non-SSE code is faster
	static real always_inline dist_sq(vec3 const&x, vec3 const&y) noexcept
	{
	  WDutilsStaticAssert(!real_is_float || aligned);
	  rvec4 X = load(x);  X-= load(y);  X*= X;
	  return sum3(X);
	}
	// dist_sq(cube,point)
	static real always_inline dist_sq(cub3 const&c, vec3 const&x) noexcept
	{
	  rvec4 X = load(c.X);
	  rvec4 H = single<3>(X);
	  X.make_diff(load(x));
	  X-= H;
	  X = max(X,rvec4::zero());
	  X*= X;
	  return sum3(X);
	}
	// dist_sq(cuboid,point)
	static real always_inline dist_sq(box3 const&c, vec3 const&x) noexcept
	{
	  rvec4 X = load(c.X);
	  X.make_diff(load(x));
	  X-= load(c.Y);
	  X = max(X,rvec4::zero());
	  X*= X;
	  return sum3(X);
	}
	// inside(cube,sphere)
	static bool inside(cub3 const&c, sph3 const&s) noexcept
	{
	  rvec4 C = load(c.X);
	  rvec4 H = single<3>(C);
	  rvec4 X = load(s.X);
	  rvec4 Q = single<3>(X);
	  X-= C;  X =abs(X);
	  return ! (7&signbits( X > H )) 
	    &&   ! (7&signbits( Q > square(H-X) ));
	}
	// inside(cuboid,sphere)
	static bool always_inline inside(box3 const&c, sph3 const&s) noexcept
	{
	  rvec4 X = load(s.X);
	  rvec4 Q = single<3>(X);
	  X-= load(c.X);  X = abs(X);
	  rvec4 H = load(c.Y);
	  return ! (7&signbits( X > H )) 
	    &&   ! (7&signbits( Q > square(H-X) ));
	}
	// convert2cuboid(PointPair)
	static void always_inline convert2cuboid(PointPair<3,real>&c) noexcept
	{
	  rvec4 X = load(c.X);
	  rvec4 Y = load(c.Y);
	  rvec4 H(real(0.5));
	  rvec4 T;
	  T = Y+X;  T*= H;  store(T,c.X);
	  T = Y-X;  T*= H;  store(T,c.Y);
	}
      };// class AlgorithmsHelper<3,real,aligned,1>
# ifdef __SSE2__
      //
      // 2.3  SSE code for 2D with real=double
      //
      template<bool aligned> class AlgorithmsHelper<2,double,aligned,1> {
	friend struct Algorithms<aligned,1>;
	typedef GeoVec<2,double> vec2;
	typedef cube  <2,double> cub2;
	typedef cuboid<2,double> box2;
	typedef sphere<2,double> sph2;
	// load a dvec2
	static dvec2 always_inline load(vec2 const&p) noexcept
	{
	  return dvec2::template load_t<aligned>(p.data());
	}
	// store a dvec2
	static void always_inline store(dvec2 x, vec2&p) noexcept
	{
	  x.template store_t<aligned>(p.data());
	}
#if(0)
	// load a dvec2
	template<typename A>
	static dvec2 always_inline load(A const&p) noexcept
	{
	  return dvec2::template pack_t<aligned>(p);
	}
	// store a dvec2
	template<typename A>
	static void always_inline store(dvec2 x, A&p) noexcept
	{
	  x.template unpack_t<aligned>(p);
	}
#endif
	// octant(point,point)
	static int always_inline octant(vec2 const&c, vec2 const&x) noexcept
	{
	  return signbits(load(c) <= load(x));
	}
	// contains(cube,point)
	static bool always_inline contains(cub2 const&c, vec2 const&x) noexcept
	{
	  dvec2 C = load(c.X);
	  dvec2 H (c.H);
	  dvec2 X = load(x);
	  return ! signbits( (C-H > X) | (X >= C+H) );
	}
	// contains_open(cube,point)
	static bool always_inline contains_open(cub2 const&c,
						vec2 const&x) noexcept
	{
	  dvec2 C = load(c.X);
	  dvec2 H (c.H);
	  dvec2 X = load(x);
	  return ! signbits( (C-H > X) | (X > C+H) );
	}
	// contains(cuboid,point)
	static bool always_inline contains(box2 const&c, vec2 const&x) noexcept
	{
	  dvec2 C = load(c.X);
	  dvec2 H = load(c.Y);
	  dvec2 X = load(x);
	  return ! signbits( (C-H > X) | (X >= C+H) );
	}
	// contains_open(cuboid,point)
	static bool always_inline contains_open(box2 const&c,
						vec2 const&x) noexcept
	{
	  dvec2 C = load(c.X);
	  dvec2 H = load(c.Y);
	  dvec2 X = load(x);
	  return ! signbits( (C-H > X) | (X > C+H) );
	}
	// contains(cuboid,cuboid)
	static bool always_inline contains(box2 const&c, box2 const&b) noexcept
	{
	  dvec2 C = diff(load(c.X),load(b.X));
	  dvec2 D = load(c.Y) - load(b.Y);
	  return ! signbits( C > D );
	}
	// dist_sq(point,point)
	static double always_inline dist_sq(vec2 const&x, vec2 const&y) noexcept
	{
	  dvec2 X = load(x) - load(y);
	  X*= X;
	  return sum(X);
	}
	// dist_sq(cube,point)
	static double always_inline dist_sq(cub2 const&c, vec2 const&x) noexcept
	{
	  dvec2 C = load(c.X);
	  dvec2 H (c.H);
	  C.make_diff(load(x));
	  C-= H;
	  C = max(C,dvec2::zero());
	  C*= C;
	  return sum(C);
	}
	// dist_sq(cuboid,point)
	static double always_inline dist_sq(box2 const&c, vec2 const&x) noexcept
	{
	  dvec2 C = load(c.X);
	  dvec2 H = load(c.Y);
	  C.make_diff(load(x));
	  C-= H;
	  C = max(C,dvec2::zero());
	  C*= C;
	  return sum(C);
	}
	// inside(cube,sphere)
	static bool always_inline inside(cub2 const&c, sph2 const&s) noexcept
	{
	  dvec2 C = load(c.X);
	  dvec2 H (c.H);
	  dvec2 X = load(s.X);
	  dvec2 Q (s.Q);
	  X.make_diff(C);
	  C = H-X;
	  C*= C;
	  return ! signbits( (Q>C) | (X>H) );
	}
	// inside(cuboid,sphere)
	static bool always_inline inside(box2 const&c, sph2 const&s) noexcept
	{
	  WDutilsStaticAssert(aligned);
	  dvec2 C = load(c.X);
	  dvec2 H = load(c.Y);
	  dvec2 X = load(s.X);
	  dvec2 Q (s.Q);
	  X.make_diff(C);
	  C = H-X;
	  C*= C;
	  return ! signbits( (Q>C) | (X>H) );
	}
	// convert2cuboid(PointPair)
	static void always_inline convert2cuboid(PointPair<2,double>&c) noexcept
	{
	  WDutilsStaticAssert(aligned);
	  dvec2 X = load(c.X);
	  dvec2 Y = load(c.Y);
	  dvec2 H(0.5);
	  dvec2 T;
	  T = Y+X;  T*= H;  store(T,c.X);
	  T = Y-X;  T*= H;  store(T,c.Y);
	}
      };// class AlgorithmsHelper<double,aligned,1>
#  ifndef __AVX__
      //
      // 2.4  SSE2 code for 3D with real=double
      //
      template<bool aligned> class AlgorithmsHelper<3,double,aligned,1> {
	friend struct Algorithms<aligned,1>;
	typedef GeoVec<3,double> vec3;
	typedef cube  <3,double> cub3;
	typedef cuboid<3,double> box3;
	typedef sphere<3,double> sph3;
	// load a dvec2
	template<int K, int Dim>
	static dvec2 always_inline load(GeoVec<Dim,double> const&p) noexcept
	{
	  WDutilsStaticAssert(K==0 || K==2);
	  return dvec2::template
	    load_t<aligned>(p.data()+K);
	}
	// store a dvec2
	template<int K, int Dim>
	static void always_inline store(dvec2 x, GeoVec<Dim,double>&p) noexcept
	{
	  WDutilsStaticAssert(K==0 || K==2);
	  x.template store_t<aligned>(p.data()+K);
	}
#if(0)
	// load a dvec2
	template<int K, typename A>
	static dvec2 always_inline load(A const&p) noexcept
	{
	  WDutilsStaticAssert(K==0 || K==2);
	  return dvec2::template
	    load_t<aligned>(static_cast<const double*>(p)+K);
	}
	// store a dvec2
	template<int K, typename A>
	static void always_inline store(dvec2 x, A&p) noexcept
	{
	  WDutilsStaticAssert(K==0 || K==2);
	  x.template store_t<aligned>(static_cast<double*>(p)+K);
	}
#endif
	// copy(cube,cube)
	static void always_inline copy(cub3 const&in, cub3 &out) noexcept
	{
	  store<0>(load<0>(in.X),out.X);
	  store<2>(load<2>(in.X),out.X);
	}
	// copy(sphere,sphere)
	static void always_inline copy(sph3 const&in, sph3 &out) noexcept
	{
	  store<0>(load<0>(in.X),out.X);
	  store<2>(load<2>(in.X),out.X);
	}
	// octant(point,point) 3D
	static int always_inline octant(vec3 const&c, vec3 const&x) noexcept
	{
	  return signbits(load<0>(c) <= load<0>(x))
	    |((1&signbits(load<2>(c) <= load<2>(x)))<<2);
	}
	// contains(cube,point)
	static bool always_inline contains(cub3 const&c, vec3 const&x) noexcept
	{
	  dvec2 C0 = load<0>(c.X);
	  dvec2 X0 = load<0>(x);
	  dvec2 C2 = load<2>(c.X);
	  dvec2 H  = single<1>(C2);
	  dvec2 X2 = load<2>(x);
	  return ! signbits( (C0-H > X0) | (X0 >= C0+H) )
	    && !(1&signbits( (C2-H > X2) | (X2 >= C2+H) ));
	}
	// contains_open(cube,point)
	static bool always_inline contains_open(cub3 const&c,
						vec3 const&x) noexcept
	{
	  dvec2 C0 = load<0>(c.X);
	  dvec2 X0 = load<0>(x);
	  dvec2 C2 = load<2>(c.X);
	  dvec2 H  = single<1>(C2);
	  dvec2 X2 = load<2>(x);
	  return ! signbits( (C0-H > X0) | (X0 > C0+H) )
	    && !(1&signbits( (C2-H > X2) | (X2 > C2+H) ));
	}
	// contains(cuboid,point)
	static bool always_inline contains(box3 const&c, vec3 const&x) noexcept
	{
	  dvec2 C0 = load<0>(c.X);
	  dvec2 H0 = load<0>(c.Y);
	  dvec2 X0 = load<0>(x);
	  dvec2 C2 = load<2>(c.X);
	  dvec2 H2 = load<2>(c.Y);
	  dvec2 X2 = load<2>(x);
	  return ! signbits( (C0-H0 > X0) | (X0 >= C0+H0) )
	    && !(1&signbits( (C2-H2 > X2) | (X2 >= C2+H2) ));
	}
	// contains_open(cuboid,point)
	static bool always_inline contains_open(box3 const&c,
						vec3 const&x) noexcept
	{
	  dvec2 C0 = load<0>(c.X);
	  dvec2 H0 = load<0>(c.Y);
	  dvec2 X0 = load<0>(x);
	  dvec2 C2 = load<2>(c.X);
	  dvec2 H2 = load<2>(c.Y);
	  dvec2 X2 = load<2>(x);
	  return ! signbits( (C0-H0 > X0) | (X0 > C0+H0) )
	    && !(1&signbits( (C2-H2 > X2) | (X2 > C2+H2) ));
	}
	// contains(cuboid,cuboid)
	static bool always_inline contains(box3 const&c, box3 const&b) noexcept
	{
	  dvec2 D0 = diff(load<0>(c.X),load<0>(b.X));
	  dvec2 S0 = load<0>(c.Y) - load<0>(b.Y);
	  dvec2 D2 = diff(load<2>(c.X),load<2>(b.X));
	  dvec2 S2 = load<2>(c.Y) - load<2>(b.Y);
	  return ! signbits( D0 > S0)
	    && !(1&signbits( D2 > S2));
	}

	// dist_sq(point,point)
	static double always_inline dist_sq(vec3 const&x, vec3 const&y) noexcept
	{
	  dvec2 X0 = load<0>(x) - load<0>(y);
	  dvec2 X2 = load<2>(x) - load<2>(y);
	  X0*= X0;
	  X2*= X2;
	  return sum(X0) + sum1(X2);
	}
	// dist_sq(cube,point)
	static double always_inline dist_sq(cub3 const&c, vec3 const&x) noexcept
	{
	  dvec2 D0 = load<0>(c.X);
	  dvec2 D2 = load<2>(c.X);
	  dvec2 H  = single<1>(D2);
	  D0.make_diff(load<0>(x));  D0-=H;  D0=max(D0,dvec2::zero());  D0*=D0;
	  D2.make_diff(load<2>(x));  D2-=H;  D2=max(D2,dvec2::zero());  D2*=D2;
	  return sum(D0) + sum1(D2);
	}
	// dist_sq(cuboid,point)
	static double always_inline dist_sq(box3 const&c, vec3 const&x) noexcept
	{
	  dvec2 D0 = load<0>(c.X);  D0.make_diff(load<0>(x));
	  D0-=load<0>(c.Y); D0=max(D0,dvec2::zero()); D0*=D0;
	  dvec2 D2 = load<2>(c.X);  D2.make_diff(load<2>(x));
	  D2-=load<2>(c.Y); D2=max(D2,dvec2::zero()); D2*=D2;
	  return sum(D0) + sum1(D2);
	}
	// inside(cube,sphere)
	static bool always_inline inside(cub3 const&c, sph3 const&s) noexcept
	{
	  WDutilsStaticAssert(aligned);
	  dvec2 X0 = load<0>(s.X);
	  X0.make_diff(load<0>(c.X));
	  dvec2 C2 = load<2>(c.X);
	  dvec2 H  = single<1>(C2);
	  dvec2 X2 = load<2>(s.X);
	  dvec2 Q  = single<1>(X2);
	  X2.make_diff(C2);
	  return ! signbits( (Q>square(H-X0)) | (X0>H) )
	    &&!(1& signbits( (Q>square(H-X2)) | (X2>H) ) );
	}
	// inside(cuboid,sphere)
	static bool always_inline inside(box3 const&c, sph3 const&s) noexcept
	{
	  WDutilsStaticAssert(aligned);
	  dvec2 H0 = load<0>(c.Y);
	  dvec2 X0 = load<0>(s.X);
	  X0.make_diff(load<0>(c.X));
	  dvec2 H2 = load<2>(c.Y);
	  dvec2 X2 = load<2>(s.X);
	  dvec2 Q  = single<1>(X2);
	  X2.make_diff(load<2>(c.X));
	  return ! signbits( (Q>square(H0-X0)) | (X0>H0) )
	    &&!(1& signbits( (Q>square(H2-X2)) | (X2>H2) ) );
	}
	// convert2cuboid(PointPair)
	static void always_inline convert2cuboid(PointPair<3,double>&c) noexcept
	{
	  dvec2 X = load<0>(c.X);
	  dvec2 Y = load<0>(c.Y);
	  dvec2 H(0.5);
	  dvec2 T;
	  T = Y+X;  T*= H;  store<0>(T,c.X);
	  T = Y-X;  T*= H;  store<0>(T,c.Y);
	  X = load<2>(c.X);
	  Y = load<2>(c.Y);
	  T = Y+X;  T*= H;  store<2>(T,c.X);
	  T = Y-X;  T*= H;  store<2>(T,c.Y);
	}
      };// class AlgorithmsHelper<3,double,aligned,1>
#  endif // ! __AVX__
# endif  // __SSE2__
#endif   // __SSE__
#undef always_inline
    } // namespace WDutils::Geometry::details
    // copy(cube,cube):  for 2D: use non-SSE
    template<bool _A, bool _S> template<int _D, typename _X> inline
    void Algorithms<_A,_S>::copy(cube<_D,_X> const&in,
				 cube<_D,_X> &out) noexcept
    { return details::AlgorithmsHelper<_D,_X,_A,_D==3&&_S>::copy(in,out); }
    // copy(sphere,sphere):  for 2D: use non-SSE
    template<bool _A, bool _S> template<int _D, typename _X> inline
    void Algorithms<_A,_S>::copy(sphere<_D,_X> const&in,
				 sphere<_D,_X> &out) noexcept
    { return details::AlgorithmsHelper<_D,_X,_A,_D==3&&_S>::copy(in,out); }
    // move_to_octant(point,int,real)
    template<bool _A, bool _S> template<int _D, typename _X> inline
    void Algorithms<_A,_S>::move_to_octant(GeoVec<_D,_X>&c,
					   int i, _X r) noexcept
    { return details::AlgorithmsHelper<_D,_X,0,0>::move_to_octant(c,i,r); }
    // octant(point,point)
    template<bool _A, bool _S> template<int _D, typename _X> inline
    int Algorithms<_A,_S>::octant(GeoVec<_D,_X> const&c,
				  GeoVec<_D,_X> const&x) noexcept
    { return details::AlgorithmsHelper<_D,_X,_A,_S>::octant(c,x); }
    // dist_sq(point,point):   prefer non-SSE if unaligned
    template<bool _A, bool _S> template<int _D, typename _X> inline
    _X Algorithms<_A,_S>::dist_sq(GeoVec<_D,_X> const&x,
				  GeoVec<_D,_X> const&y) noexcept
    { return details::AlgorithmsHelper<_D,_X,_A,_A&&_S>::dist_sq(x,y); }
    // contains(cube,point)
    template<bool _A, bool _S> template<int _D, typename _X> inline
    bool Algorithms<_A,_S>::contains(cube  <_D,_X> const&c,
				     GeoVec<_D,_X> const&x) noexcept
    { return details::AlgorithmsHelper<_D,_X,_A,_S>::contains(c,x); }
    // contains_open(cube,point)
    template<bool _A, bool _S> template<int _D, typename _X> inline
    bool Algorithms<_A,_S>::contains_open(cube  <_D,_X> const&c,
					  GeoVec<_D,_X> const&x) noexcept
    { return details::AlgorithmsHelper<_D,_X,_A,_S>::contains_open(c,x); }
    // dist_sq(cube,point)
    template<bool _A, bool _S> template<int _D, typename _X> inline
    _X Algorithms<_A,_S>::dist_sq(cube  <_D,_X> const&c,
				  GeoVec<_D,_X> const&x) noexcept
    { return details::AlgorithmsHelper<_D,_X,_A,_S>::dist_sq(c,x); }
    // outside(cube,sphere)
    template<bool _A, bool _S> template<int _D, typename _X> inline
    bool Algorithms<_A,_S>::outside(cube<_D,_X> const&c,
				    sphere<_D,_X> const&s) noexcept
    { return dist_sq(c,s.X) > s.Q; }
    // inside(cube,sphere)):  SSE only if aligned or float or 2D
    template<bool _A, bool _S> template<int _D, typename _X> inline
    bool Algorithms<_A,_S>::inside(cube<_D,_X> const&c,
				   sphere<_D,_X> const&s) noexcept
    {
      return details::AlgorithmsHelper<_D,_X,_A,
#ifndef __AVX__
	(_D==2 || is_same<float,_X>::value || _A) &&
#endif
	_S>::inside(c,s);
    }
    // contains(cuboid,point)
    template<bool _A, bool _S> template<int _D, typename _X> inline
    bool Algorithms<_A,_S>::contains(cuboid<_D,_X> const&c,
				     GeoVec<_D,_X> const&x) noexcept
    { return details::AlgorithmsHelper<_D,_X,_A,_S>::contains(c,x); }
    // contains_open(cuboid,point)
    template<bool _A, bool _S> template<int _D, typename _X> inline
    bool Algorithms<_A,_S>::contains_open(cuboid<_D,_X> const&c,
					  GeoVec<_D,_X> const&x) noexcept
    { return details::AlgorithmsHelper<_D,_X,_A,_S>::contains_open(c,x); }
    // contains(cuboid,cuboid):   prefer non-SSE if unaligned
    template<bool _A, bool _S> template<int _D, typename _X> inline
    bool Algorithms<_A,_S>::contains(cuboid<_D,_X> const&c,
				     cuboid<_D,_X> const&b) noexcept
    { return details::AlgorithmsHelper<_D,_X,_A,_A&&_S>::contains(c,b); }
    // dist_sq(cuboid,point)
    template<bool _A, bool _S> template<int _D, typename _X> inline
    _X Algorithms<_A,_S>::dist_sq(cuboid<_D,_X> const&c,
				  GeoVec<_D,_X> const&x) noexcept
    { return details::AlgorithmsHelper<_D,_X,_A,_S>::dist_sq(c,x); }
    // outside(cuboid,point)
    template<bool _A, bool _S> template<int _D, typename _X> inline
    bool Algorithms<_A,_S>::outside(cuboid<_D,_X> const&c,
				    sphere<_D,_X> const&s) noexcept
    { return dist_sq(c,s.X) > s.Q; }
    // inside(cuboid,sphere):  SSE only if aligned or float
    template<bool _A, bool _S> template<int _D, typename _X> inline
    bool Algorithms<_A,_S>::inside(cuboid<_D,_X> const&c,
				   sphere<_D,_X> const&s) noexcept
    { 
      return details::AlgorithmsHelper<_D,_X,_A,
	(is_same<float,_X>::value || _A) && _S>::inside(c,s);
    }
    // convert2cuboid(point-pair):   for unaligned: prefer non-SSE version
    template<bool _A, bool _S> template<int _D, typename _X> inline
    void Algorithms<_A,_S>::convert2cuboid(PointPair<_D,_X> &p) noexcept
    { return details::AlgorithmsHelper<_D,_X,_A,_A&&_S>::convert2cuboid(p); }
  } // nameapce WDutils::Geometry
} // namespace WDutils
#undef noexcept
//
#endif // WDutils_included_geometry_tcc
