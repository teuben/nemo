// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/geometry_inl.h
///
/// \brief  implements inline methods for geometry.h
///
/// \author Walter Dehnen
///
/// \date   2010,2012
///
////////////////////////////////////////////////////////////////////////////////
///
/// \version June-2010 WD  first tested version
/// \version June-2012 WD  using vector instead of tupel if C++11
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
#ifndef WDutils_included_geometry_inl_h
#define WDutils_included_geometry_inl_h

#ifndef WDutils_included_geometry_h
#  include <geometry.h>
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
    namespace Meta {
      //
      // 1  generic case: not using SSE, ignoring alignment
      //
      template<typename real, bool aligned, bool sse>
      class AlgorithmsHelper
      {
	friend struct Algorithms<aligned,sse>;
      protected:
	//
	typedef GeoVec<2,real> vec2;
	typedef GeoVec<3,real> vec3;
	typedef cube  <2,real> cub2;
	typedef cube  <3,real> cub3;
	typedef cuboid<2,real> box2;
	typedef cuboid<3,real> box3;
	typedef sphere<2,real> sph2;
	typedef sphere<3,real> sph3;
	// is @a x in the interval [@a c-r, @a c+r) ?
	static bool in_interval(real c, real r, real x)
	{
	  return x<c? c<=x+r : x<c+r;
	}
	// is the interval [a-b, a+b] in the interval [c-r, c+r]?
	static bool in_interval(real c, real r, real a, real b)
	{
	  return std::abs(c-a) <= r-b;
	}
	// copy(cube,cube)
	template<int Dim>
	static void copy(cube<Dim,real> const&in, cube<Dim,real> &out)
	{
	  out.X = in.X;
	  out.H = in.H;
	}
	// copy(sphere,sphere)
	template<int Dim>
	static void copy(sphere<Dim,real> const&in, sphere<Dim,real> &out)
	{
	  out.X = in.X;
	  out.Q = in.Q;
	}
	// octant(point,point)
	static int octant(vec2 const&c, vec2 const&x)
	{
	  int oct(0);
	  if(x[0] > c[0]) oct |= 1;
	  if(x[1] > c[1]) oct |= 2;
	  return oct;
	}
	static int octant(vec3 const&c, vec3 const&x)
	{
	  int oct(0);
	  if(x[0] > c[0]) oct |= 1;
	  if(x[1] > c[1]) oct |= 2;
	  if(x[2] > c[2]) oct |= 4;
	  return oct;
	}
	// contains(cube,point)
	static bool contains(cub2 const&c, vec2 const&x)
	{
	  return in_interval(c.X[0],c.H,x[0])
	    &&   in_interval(c.X[1],c.H,x[1]);
	}
	static bool contains(cub3 const&c, vec3 const&x)
	{
	  return in_interval(c.X[0],c.H,x[0])
	    &&   in_interval(c.X[1],c.H,x[1])
	    &&   in_interval(c.X[2],c.H,x[2]);
	}
	// contains(cuboid,point)
	static bool contains(box2 const&c, vec2 const&x)
	{
	  return in_interval(c.X[0],c.Y[0],x[0])
	    &&   in_interval(c.X[1],c.Y[1],x[1]);
	}
	static bool contains(box3 const&c, vec3 const&x)
	{
	  return in_interval(c.X[0],c.Y[0],x[0])
	    &&   in_interval(c.X[1],c.Y[1],x[1])
	    &&   in_interval(c.X[2],c.Y[2],x[2]);
	}
	// contains(cuboid,cuboid)
	static bool contains(box2 const&c, box2 const&b)
	{
	  return in_interval(c.X[0],c.Y[0], b.X[0],b.Y[0])
	    &&   in_interval(c.X[1],c.Y[1], b.X[1],b.Y[1]);
	}
	static bool contains(box3 const&c, box3 const&b)
	{
	  return in_interval(c.X[0],c.Y[0], b.X[0],b.Y[0])
	    &&   in_interval(c.X[1],c.Y[1], b.X[1],b.Y[1])
	    &&   in_interval(c.X[2],c.Y[2], b.X[2],b.Y[2]);
	}
	// dist_sq(point,point)
	static real dist_sq(vec2 const&x, vec2 const&y)
	{
	  return WDutils::dist_sq(x,y);
	}
	static real dist_sq(vec3 const&x, vec3 const&y)
	{
	  return WDutils::dist_sq(x,y);
	}
	// dist_sq(cube,point)
	static real dist_sq(cub2 const&c, vec2 const&x)
	{
	  real q(0),D;
	  D = abs(c.X[0]-x[0]); if(D>c.H) q+=square(D-c.H);
	  D = abs(c.X[1]-x[1]); if(D>c.H) q+=square(D-c.H);
	  return q;
	}
	static real dist_sq(cub3 const&c, vec3 const&x)
	{
	  real q(0),D;
	  D = abs(c.X[0]-x[0]); if(D>c.H) q+=square(D-c.H);
	  D = abs(c.X[1]-x[1]); if(D>c.H) q+=square(D-c.H);
	  D = abs(c.X[2]-x[2]); if(D>c.H) q+=square(D-c.H);
	  return q;
	}
	// dist_sq(cuboid,point)
	static real dist_sq(box2 const&c, vec2 const&x)
	{
	  real q(0),D;
	  D = abs(c.X[0]-x[0])-c.Y[0]; if(D>0) q+=D*D;
	  D = abs(c.X[1]-x[1])-c.Y[1]; if(D>0) q+=D*D;
	  return q;
	}
	static real dist_sq(box3 const&c, vec3 const&x)
	{
	  real q(0),D;
	  D = abs(c.X[0]-x[0])-c.Y[0]; if(D>0) q+=D*D;
	  D = abs(c.X[1]-x[1])-c.Y[1]; if(D>0) q+=D*D;
	  D = abs(c.X[2]-x[2])-c.Y[2]; if(D>0) q+=D*D;
	  return q;
	}
	// outside(cube,sphere)
	static bool outside(cub2 const&c, sph2 const&s)
	{
	  real q(0),D;
	  D=abs(c.X[0]-s.X[0])-c.H; if(D>0 && s.Q<(q+=D*D)) return true;
	  D=abs(c.X[1]-s.X[1])-c.H; if(D>0 && s.Q<(q+=D*D)) return true;
	  return false;
	}
	static bool outside(cub3 const&c, sph3 const&s)
	{
	  real q(0),D;
	  D=abs(c.X[0]-s.X[0])-c.H; if(D>0 && s.Q<(q+=D*D)) return true;
	  D=abs(c.X[1]-s.X[1])-c.H; if(D>0 && s.Q<(q+=D*D)) return true;
	  D=abs(c.X[2]-s.X[2])-c.H; if(D>0 && s.Q<(q+=D*D)) return true;
	  return false;
	}
	// outside(cuboid,sphere)
	static bool outside(box2 const&c, sph2 const&s)
	{
	  real q(0),D;
	  D=abs(c.X[0]-s.X[0])-c.Y[0]; if(D>0 && s.Q<(q+=D*D)) return true;
	  D=abs(c.X[1]-s.X[1])-c.Y[1]; if(D>0 && s.Q<(q+=D*D)) return true;
	  return false;
	}
	static bool outside(box3 const&c, sph3 const&s)
	{
	  real q(0),D;
	  D=abs(c.X[0]-s.X[0])-c.Y[0]; if(D>0 && s.Q<(q+=D*D)) return true;
	  D=abs(c.X[1]-s.X[1])-c.Y[1]; if(D>0 && s.Q<(q+=D*D)) return true;
	  D=abs(c.X[2]-s.X[2])-c.Y[2]; if(D>0 && s.Q<(q+=D*D)) return true;
	  return false;
	}
	// inside(cube,sphere)
	static bool inside(cub2 const&c, sph2 const&s)
	{
	  real D;
	  D=abs(c.X[0]-s.X[0])-c.H; if(D>0 || s.Q>D*D) return false;
	  D=abs(c.X[1]-s.X[1])-c.H; if(D>0 || s.Q>D*D) return false;
	  return true;
	}
	static bool inside(cub3 const&c, sph3 const&s)
	{
	  real D;
	  D=abs(c.X[0]-s.X[0])-c.H; if(D>0 || s.Q>D*D) return false;
	  D=abs(c.X[1]-s.X[1])-c.H; if(D>0 || s.Q>D*D) return false;
	  D=abs(c.X[2]-s.X[2])-c.H; if(D>0 || s.Q>D*D) return false;
	  return true;
	}
	// inside(cuboid,sphere)
	static bool inside(box2 const&c, sph2 const&s)
	{
	  real D;
	  D=abs(c.X[0]-s.X[0])-c.Y[0]; if(D>0 || s.Q>D*D) return false;
	  D=abs(c.X[1]-s.X[1])-c.Y[1]; if(D>0 || s.Q>D*D) return false;
	  return true;
	}
	static bool inside(box3 const&c, sph3 const&s)
	{
	  real D;
	  D=abs(c.X[0]-s.X[0])-c.Y[0]; if(D>0 || s.Q>D*D) return false;
	  D=abs(c.X[1]-s.X[1])-c.Y[1]; if(D>0 || s.Q>D*D) return false;
	  D=abs(c.X[2]-s.X[2])-c.Y[2]; if(D>0 || s.Q>D*D) return false;
	  return true;
	}
      public:
	// convert2cuboid(PointPair)
	template<int Dim>
	static void convert2cuboid(PointPair<Dim,real> &c)
	{
	  GeoVec<Dim,real> C = real(0.5)*(c.X+c.Y);
	  GeoVec<Dim,real> H = real(0.5)*(c.Y-c.X);
	  static_cast<GeoVec<Dim,real>&>(c.X) = C;
	  static_cast<GeoVec<Dim,real>&>(c.Y) = H;
	}
	// move to octant
	static void move_to_octant(vec2&c, int i, real r)
	{
	  if(i&1) c[0] += r;  else  c[0] -= r;
	  if(i&2) c[1] += r;  else  c[1] -= r;
	}
	static void move_to_octant(vec3&c, int i, real r)
	{
	  if(i&1) c[0] += r;  else  c[0] -= r;
	  if(i&2) c[1] += r;  else  c[1] -= r;
	  if(i&4) c[2] += r;  else  c[2] -= r;
	}
      };
#ifdef __SSE__
      //
      // 2   specialisations using SSE
      //
#define  PF(__X) static_cast<float*>(__X)
#define cPF(__X) static_cast<const float*>(__X)
      //
      // 2.1 specialisation for unaligned real=float with SSE
      //
      template<> class AlgorithmsHelper<float,0,1> {
	friend struct Algorithms<0,1>;
	typedef GeoVec<2,float> vec2;
	typedef GeoVec<3,float> vec3;
	typedef cube  <2,float> cub2;
	typedef cube  <3,float> cub3;
	typedef cuboid<2,float> box2;
	typedef cuboid<3,float> box3;
	typedef sphere<2,float> sph2;
	typedef sphere<3,float> sph3;
	// copy(cube,cube)
	static void copy(cub2 const&in, cub2 &out)
	{
	  out.X = in.X;
	  out.H = in.H;
	}
	static void copy(cub3 const&in, cub3 &out)
	{
	  _mm_storeu_ps(PF(out.X),_mm_loadu_ps(cPF(in.X)));
	}
	// copy(sphere,sphere)
	static void copy(sph2 const&in, sph2 &out)
	{
	  out.X = in.X;
	  out.Q = in.Q;
	}
	static void copy(sph3 const&in, sph3 &out)
	{
	  _mm_storeu_ps(PF(out.X),_mm_loadu_ps(cPF(in.X)));
	}
	// octant(point,point)
	static int octant(vec2 const&c, vec2 const&x)
	{
	  return 3&_mm_movemask_ps(_mm_cmplt_ps(_mm_loadu_ps(cPF(c)),
						_mm_loadu_ps(cPF(x))));
	}
	static int octant(vec3 const&c, vec3 const&x)
	{
	  return 7&_mm_movemask_ps(_mm_cmplt_ps(_mm_loadu_ps(cPF(c)),
						_mm_loadu_ps(cPF(x))));
	}
	// contains(cube,point)
	static bool contains(cub2 const&c, vec2 const&x)
	{
	  __m128 C = _mm_loadu_ps(cPF(c.X));
	  __m128 X = _mm_loadu_ps(cPF(x));
	  __m128 H = _mm_shuffle_ps(C,C,_MM_SHUFFLE(2,2,2,2));
	  return
	    3==(3&_mm_movemask_ps(_mm_and_ps(_mm_cmple_ps(_mm_sub_ps(C,H),X),
					     _mm_cmpgt_ps(_mm_add_ps(C,H),X))));
	}
	static bool contains(cub3 const&c, vec3 const&x)
	{
	  __m128 C = _mm_loadu_ps(cPF(c.X));
	  __m128 X = _mm_loadu_ps(cPF(x));
	  __m128 H = _mm_shuffle_ps(C,C,_MM_SHUFFLE(3,3,3,3));
	  return
	    7==(7&_mm_movemask_ps(_mm_and_ps(_mm_cmple_ps(_mm_sub_ps(C,H),X),
					     _mm_cmpgt_ps(_mm_add_ps(C,H),X))));
	}
	// contains(cuboid,point)
	static bool contains(box2 const&c, vec2 const&x)
	{
	  __m128 C = _mm_loadu_ps(cPF(c.X));
	  __m128 H = _mm_movehl_ps(C,C);
	  __m128 X = _mm_loadu_ps(cPF(x));
	  return
	    3==(3&_mm_movemask_ps(_mm_and_ps(_mm_cmple_ps(_mm_sub_ps(C,H),X),
					     _mm_cmpgt_ps(_mm_add_ps(C,H),X))));
	}
	static bool contains(box3 const&c, vec3 const&x)
	{
	  __m128 C = _mm_loadu_ps(cPF(c.X));
	  __m128 H = _mm_loadu_ps(cPF(c.Y));
	  __m128 X = _mm_loadu_ps(cPF(x));
	  return
	    7==(7&_mm_movemask_ps(_mm_and_ps(_mm_cmple_ps(_mm_sub_ps(C,H),X),
					     _mm_cmpgt_ps(_mm_add_ps(C,H),X))));
	}
	// contains(cuboid,cuboid)
	static bool contains(box2 const&c, box2 const&b)
	{
	  __m128 C = _mm_loadu_ps(cPF(c.X));
	  __m128 R = _mm_movehl_ps(C,C);
	  __m128 A = _mm_loadu_ps(cPF(b.X));
	  __m128 B = _mm_movehl_ps(A,A);
	  return 3==(3&_mm_movemask_ps(_mm_cmple_ps(SSE::diff(C,A),
						    _mm_sub_ps(R,B))));
	}
	static bool contains(box3 const&c, box3 const&b)
	{
	  __m128 C = _mm_loadu_ps(cPF(c.X));
	  __m128 R = _mm_loadu_ps(cPF(c.Y));
	  __m128 A = _mm_loadu_ps(cPF(b.X));
	  __m128 B = _mm_loadu_ps(cPF(b.Y));
	  return 7==(7&_mm_movemask_ps(_mm_cmple_ps(SSE::diff(C,A),
						    _mm_sub_ps(R,B))));
	}
	// dist_sq(point,point)
	// for unaligned access, non-SSE code is faster (don't ask me why)
#if(0)
	static float dist_sq(vec2 const&x, vec2 const&y)
	{
	  WDutils__align16 float q[4];
	  __m128 D = _mm_sub_ps(_mm_loadu_ps(cPF(x)),_mm_loadu_ps(cPF(y)));
	  _mm_store_ps(q,_mm_mul_ps(D,D));
	  return q[0]+q[1];
	}
	static float dist_sq(vec3 const&x, vec3 const&y)
	{
	  WDutils__align16 float q[4];
	  __m128 D = _mm_sub_ps(_mm_loadu_ps(cPF(x)),_mm_loadu_ps(cPF(y)));
	  _mm_store_ps(q,_mm_mul_ps(D,D));
	  return q[0]+q[1]+q[2];
	}
#else
	static float dist_sq(vec2 const&x, vec2 const&y)
	{ return WDutils::dist_sq(x,y); }
	static float dist_sq(vec3 const&x, vec3 const&y)
	{ return WDutils::dist_sq(x,y); }
#endif
	// dist_sq(cube,point)
	static float dist_sq(cub2 const&c, vec2 const&x)
	{
	  __m128 C = _mm_loadu_ps(cPF(c.X));
	  __m128 X = _mm_loadu_ps(cPF(x));
	  __m128 H = _mm_shuffle_ps(C,C,_MM_SHUFFLE(2,2,2,2));
	  X = _mm_sub_ps(SSE::diff(C,X),H);
	  WDutils__align16 float q[4];
	  _mm_store_ps(q,_mm_and_ps(_mm_mul_ps(X,X),
				    _mm_cmpgt_ps(X,_mm_setzero_ps())));
	  return q[0]+q[1];
	}
	static float dist_sq(cub3 const&c, vec3 const&x)
	{
	  WDutils__align16 float q[4];
	  __m128 C = _mm_loadu_ps(cPF(c.X));
	  __m128 X = _mm_loadu_ps(cPF(x));
	  __m128 H = _mm_shuffle_ps(C,C,_MM_SHUFFLE(3,3,3,3));
	  X = _mm_sub_ps(SSE::diff(C,X),H);
	  _mm_store_ps(q,_mm_and_ps(_mm_mul_ps(X,X),
				    _mm_cmpgt_ps(X,_mm_setzero_ps())));
	  return q[0]+q[1]+q[2];
	}
	// dist_sq(cuboid,point)
	static float dist_sq(box2 const&c, vec2 const&x)
	{
	  __m128 C = _mm_loadu_ps(cPF(c.X));
	  __m128 H = _mm_movehl_ps(C,C);
	  __m128 X = _mm_loadu_ps(cPF(x));
	  X = _mm_sub_ps(SSE::diff(C,X),H);
	  WDutils__align16 float q[4];
	  _mm_store_ps(q,_mm_and_ps(_mm_mul_ps(X,X),
				    _mm_cmpgt_ps(X,_mm_setzero_ps())));
	  return q[0]+q[1];
	}
	static float dist_sq(box3 const&c, vec3 const&x)
	{
	  WDutils__align16 float q[4];
	  __m128 C = _mm_loadu_ps(cPF(c.X));
	  __m128 H = _mm_loadu_ps(cPF(c.Y));
	  __m128 X = _mm_loadu_ps(cPF(x));
	  X = _mm_sub_ps(SSE::diff(C,X),H);
	  _mm_store_ps(q,_mm_and_ps(_mm_mul_ps(X,X),
				    _mm_cmpgt_ps(X,_mm_setzero_ps())));
	  return q[0]+q[1]+q[2];
	}
	// outside(cube,sphere)
	static bool outside(cub2 const&c, sph2 const&s)
	{ return dist_sq(c,s.X) > s.Q; }
	static bool outside(cub3 const&c, sph3 const&s)
	{ return dist_sq(c,s.X) > s.Q; }
	// outside(cuboid,sphere)
	static bool outside(box2 const&c, sph2 const&s)
	{ return dist_sq(c,s.X) > s.Q; }
	static bool outside(box3 const&c, sph3 const&s)
	{ return dist_sq(c,s.X) > s.Q; }
	// inside(cube,sphere)
	static bool inside(cub2 const&c, sph2 const&s)
	{
	  __m128 C = _mm_loadu_ps(cPF(c.X));
	  __m128 H = _mm_shuffle_ps(C,C,_MM_SHUFFLE(2,2,2,2));
	  __m128 X = _mm_loadu_ps(cPF(s.X));
	  C = SSE::diff(C,X);
	  if(3&_mm_movemask_ps(_mm_cmpgt_ps(C,H)))
	    return false;
	  X = _mm_shuffle_ps(X,X,_MM_SHUFFLE(2,2,2,2));
	  C = _mm_sub_ps(C,H);
	  C = _mm_mul_ps(C,C);
	  return ! (3&_mm_movemask_ps(_mm_cmpgt_ps(X,C)));
	}
	static bool inside(cub3 const&c, sph3 const&s)
	{
	  __m128 C = _mm_loadu_ps(cPF(c.X));
	  __m128 H = _mm_shuffle_ps(C,C,_MM_SHUFFLE(3,3,3,3));
	  __m128 X = _mm_loadu_ps(cPF(s.X));
	  C = SSE::diff(C,X);
	  if(7&_mm_movemask_ps(_mm_cmpgt_ps(C,H)))
	    return false;
	  X = _mm_shuffle_ps(X,X,_MM_SHUFFLE(3,3,3,3));
	  C = _mm_sub_ps(C,H);
	  C = _mm_mul_ps(C,C);
	  return ! (7&_mm_movemask_ps(_mm_cmpgt_ps(X,C)));
	}
	// inside(cuboid,sphere)
	static bool inside(box2 const&c, sph2 const&s)
	{
	  __m128 C = _mm_loadu_ps(cPF(c.X));
	  __m128 H = _mm_movehl_ps(C,C);
	  __m128 X = _mm_loadu_ps(cPF(s.X));
	  C = SSE::diff(C,X);
	  if(3&_mm_movemask_ps(_mm_cmpgt_ps(C,H)))
	    return false;
	  X = _mm_shuffle_ps(X,X,_MM_SHUFFLE(2,2,2,2));
	  C = _mm_sub_ps(C,H);
	  C = _mm_mul_ps(C,C);
	  return ! (3&_mm_movemask_ps(_mm_cmpgt_ps(X,C)));
	}
	static bool inside(box3 const&c, sph3 const&s)
	{
	  __m128 C = _mm_loadu_ps(cPF(c.X));
	  __m128 H = _mm_loadu_ps(cPF(c.Y));
	  __m128 X = _mm_loadu_ps(cPF(s.X));
	  C = SSE::diff(C,X);
	  if(7&_mm_movemask_ps(_mm_cmpgt_ps(C,H)))
	    return false;
	  X = _mm_shuffle_ps(X,X,_MM_SHUFFLE(3,3,3,3));
	  C = _mm_sub_ps(C,H);
	  C = _mm_mul_ps(C,C);
	  return ! (7&_mm_movemask_ps(_mm_cmpgt_ps(X,C)));
	}
	// convert2cuboid(PointPair)
	static void convert2cuboid(PointPair<2,float> &c)
	{
	  __m128 XY,YX;
	  XY =_mm_loadu_ps(cPF(c.X));
	  YX =_mm_shuffle_ps(XY,XY,_MM_SHUFFLE(1,0,3,2));
#ifdef __SSE2__
	  YX =_mm_xor_ps(YX,reinterpret_cast<__m128>
			 (_mm_set_epi32(0x80000000,0x80000000,0,0)));
#else
	  XY =_mm_xor_ps(YX,reinterpret_cast<__m128>
			 (_mm_set_ps(meta::__neg_mask.__F,
				     meta::__neg_mask.__F,0,0)));
#endif
	  _mm_storeu_ps(PF(c.X),_mm_mul_ps(_mm_set1_ps(0.5f),
					   _mm_add_ps(XY,YX)));
	}
	static void convert2cuboid(PointPair<3,float>&c)
	{
	  __m128 Xi = _mm_loadu_ps(cPF(c.X));
	  __m128 Xa = _mm_loadu_ps(cPF(c.Y));
	  __m128 h  = _mm_set1_ps(0.5f);
	  _mm_storeu_ps(PF(c.X),_mm_mul_ps(h,_mm_add_ps(Xi,Xa)));
	  _mm_storeu_ps(PF(c.Y),_mm_mul_ps(h,_mm_sub_ps(Xa,Xi)));
	}
      };// class AlgorithmsHelper<float,0,1>
      //
      // 2.2 specialisation for aligned real=float with SSE
      //
      template<> class AlgorithmsHelper<float,1,1> {
	friend struct Algorithms<1,1>;
	typedef GeoVec<2,float> vec2;
	typedef GeoVec<3,float> vec3;
	typedef cube  <2,float> cub2;
	typedef cube  <3,float> cub3;
	typedef cuboid<2,float> box2;
	typedef cuboid<3,float> box3;
	typedef sphere<2,float> sph2;
	typedef sphere<3,float> sph3;
	// copy(cube,cube)
	static void copy(cub2 const&in, cub2 &out)
	{
	  out.X = in.X;
	  out.H = in.H;
	}
	static void copy(cub3 const&in, cub3 &out)
	{
	  _mm_store_ps(PF(out.X),_mm_load_ps(cPF(in.X)));
	}
	// copy(sphere,sphere)
	static void copy(sph2 const&in, sph2 &out)
	{
	  out.X = in.X;
	  out.Q = in.Q;
	}
	static void copy(sph3 const&in, sph3 &out)
	{
	  _mm_store_ps(PF(out.X),_mm_load_ps(cPF(in.X)));
	}
	// octant(point,point)
	static int octant(vec2 const&c, vec2 const&x)
	{
	  return 3&_mm_movemask_ps(_mm_cmplt_ps(_mm_load_ps(cPF(c)),
						_mm_load_ps(cPF(x))));
	}
	static int octant(vec3 const&c, vec3 const&x)
	{
	  return 7&_mm_movemask_ps(_mm_cmplt_ps(_mm_load_ps(cPF(c)),
						_mm_load_ps(cPF(x))));
	}
	// contains(cube,point)
	static bool contains(cub2 const&c, vec2 const&x)
	{
	  __m128 C = _mm_load_ps(cPF(c.X));
	  __m128 X = _mm_load_ps(cPF(x));
	  __m128 H = _mm_shuffle_ps(C,C,_MM_SHUFFLE(2,2,2,2));
	  return
	    3==(3&_mm_movemask_ps(_mm_and_ps(_mm_cmple_ps(_mm_sub_ps(C,H),X),
					     _mm_cmpgt_ps(_mm_add_ps(C,H),X))));
	}
	static bool contains(cub3 const&c, vec3 const&x)
	{
	  __m128 C = _mm_load_ps(cPF(c.X));
	  __m128 X = _mm_load_ps(cPF(x));
	  __m128 H = _mm_shuffle_ps(C,C,_MM_SHUFFLE(3,3,3,3));
	  return
	    7==(7&_mm_movemask_ps(_mm_and_ps(_mm_cmple_ps(_mm_sub_ps(C,H),X),
					     _mm_cmpgt_ps(_mm_add_ps(C,H),X))));
	}
	// contains(cuboid,point)
	static bool contains(box2 const&c, vec2 const&x)
	{
	  __m128 C = _mm_load_ps(cPF(c.X));
	  __m128 H = _mm_movehl_ps(C,C);
	  __m128 X = _mm_load_ps(cPF(x));
	  return
	    3==(3&_mm_movemask_ps(_mm_and_ps(_mm_cmple_ps(_mm_sub_ps(C,H),X),
					     _mm_cmpgt_ps(_mm_add_ps(C,H),X))));
	}
	static bool contains(box3 const&c, vec3 const&x)
	{
	  __m128 C = _mm_load_ps(cPF(c.X));
	  __m128 H = _mm_load_ps(cPF(c.Y));
	  __m128 X = _mm_load_ps(cPF(x));
	  return
	    7==(7&_mm_movemask_ps(_mm_and_ps(_mm_cmple_ps(_mm_sub_ps(C,H),X),
					     _mm_cmpgt_ps(_mm_add_ps(C,H),X))));
	}
	// contains(cuboid,cuboid)
	static bool contains(box2 const&c, box2 const&b)
	{
	  __m128 C = _mm_load_ps(cPF(c.X));
	  __m128 R = _mm_movehl_ps(C,C);
	  __m128 A = _mm_load_ps(cPF(b.X));
	  __m128 B = _mm_movehl_ps(A,A);
	  return 3==(3&_mm_movemask_ps(_mm_cmple_ps(SSE::diff(C,A),
						    _mm_sub_ps(R,B))));
	}
	static bool contains(box3 const&c, box3 const&b)
	{
	  __m128 C = _mm_load_ps(cPF(c.X));
	  __m128 R = _mm_load_ps(cPF(c.Y));
	  __m128 A = _mm_load_ps(cPF(b.X));
	  __m128 B = _mm_load_ps(cPF(b.Y));
	  return 7==(7&_mm_movemask_ps(_mm_cmple_ps(SSE::diff(C,A),
						    _mm_sub_ps(R,B))));
	}
	// dist_sq(point,point)
	static float dist_sq(vec2 const&x, vec2 const&y)
	{
	  WDutils__align16 float q[4];
	  __m128 D = _mm_sub_ps(_mm_load_ps(cPF(x)),_mm_load_ps(cPF(y)));
	  _mm_store_ps(q,_mm_mul_ps(D,D));
	  return q[0]+q[1];
	}
	static float dist_sq(vec3 const&x, vec3 const&y)
	{
	  WDutils__align16 float q[4];
	  __m128 D = _mm_sub_ps(_mm_load_ps(cPF(x)),_mm_load_ps(cPF(y)));
	  _mm_store_ps(q,_mm_mul_ps(D,D));
	  return q[0]+q[1]+q[2];
	}
	// dist_sq(cube,point)
	static float dist_sq(cub2 const&c, vec2 const&x)
	{
	  WDutils__align16 float q[4];
	  __m128 C = _mm_load_ps(cPF(c.X));
	  __m128 X = _mm_load_ps(cPF(x));
	  __m128 H = _mm_shuffle_ps(C,C,_MM_SHUFFLE(2,2,2,2));
	  X = _mm_sub_ps(SSE::diff(C,X),H);
	  _mm_store_ps(q,_mm_and_ps(_mm_mul_ps(X,X),
				    _mm_cmpgt_ps(X,_mm_setzero_ps())));
	  return q[0]+q[1];
	}
	static float dist_sq(cub3 const&c, vec3 const&x)
	{
	  WDutils__align16 float q[4];
	  __m128 C = _mm_load_ps(cPF(c.X));
	  __m128 X = _mm_load_ps(cPF(x));
	  __m128 H = _mm_shuffle_ps(C,C,_MM_SHUFFLE(3,3,3,3));
	  X = _mm_sub_ps(SSE::diff(C,X),H);
	  _mm_store_ps(q,_mm_and_ps(_mm_mul_ps(X,X),
				    _mm_cmpgt_ps(X,_mm_setzero_ps())));
	  return q[0]+q[1]+q[2];
	}
	// dist_sq(cuboid,point)
	static float dist_sq(box2 const&c, vec2 const&x)
	{
	  __m128 C = _mm_load_ps(cPF(c.X));
	  __m128 H = _mm_movehl_ps(C,C);
	  __m128 X = _mm_load_ps(cPF(x));
	  X = _mm_sub_ps(SSE::diff(C,X),H);
	  WDutils__align16 float q[4];
	  _mm_store_ps(q,_mm_and_ps(_mm_mul_ps(X,X),
				    _mm_cmpgt_ps(X,_mm_setzero_ps())));
	  return q[0]+q[1];
	}
	static float dist_sq(box3 const&c, vec3 const&x)
	{
	  WDutils__align16 float q[4];
	  __m128 C = _mm_load_ps(cPF(c.X));
	  __m128 H = _mm_load_ps(cPF(c.Y));
	  __m128 X = _mm_load_ps(cPF(x));
	  X = _mm_sub_ps(SSE::diff(C,X),H);
	  _mm_store_ps(q,_mm_and_ps(_mm_mul_ps(X,X),
				    _mm_cmpgt_ps(X,_mm_setzero_ps())));
	  return q[0]+q[1]+q[2];
	}
	// outside(cube,sphere)
	static bool outside(cub2 const&c, sph2 const&s)
	{ return dist_sq(c,s.X) > s.Q; }
	static bool outside(cub3 const&c, sph3 const&s)
	{ return dist_sq(c,s.X) > s.Q; }
	// outside(cuboid,sphere)
	static bool outside(box2 const&c, sph2 const&s)
	{ return dist_sq(c,s.X) > s.Q; }
	static bool outside(box3 const&c, sph3 const&s)
	{ return dist_sq(c,s.X) > s.Q; }
	// inside(cube,sphere)
	static bool inside(cub2 const&c, sph2 const&s)
	{
	  __m128 C = _mm_load_ps(cPF(c.X));
	  __m128 H = _mm_shuffle_ps(C,C,_MM_SHUFFLE(2,2,2,2));
	  __m128 X = _mm_load_ps(cPF(s.X));
	  C = SSE::diff(C,X);
	  if(3&_mm_movemask_ps(_mm_cmpgt_ps(C,H)))
	    return false;
	  X = _mm_shuffle_ps(X,X,_MM_SHUFFLE(2,2,2,2));
	  C = _mm_sub_ps(C,H);
	  C = _mm_mul_ps(C,C);
	  return ! (3&_mm_movemask_ps(_mm_cmpgt_ps(X,C)));
	}
	static bool inside(cub3 const&c, sph3 const&s)
	{
	  __m128 C = _mm_load_ps(cPF(c.X));
	  __m128 H = _mm_shuffle_ps(C,C,_MM_SHUFFLE(3,3,3,3));
	  __m128 X = _mm_load_ps(cPF(s.X));
	  C = SSE::diff(C,X);
	  if(7&_mm_movemask_ps(_mm_cmpgt_ps(C,H)))
	    return false;
	  X = _mm_shuffle_ps(X,X,_MM_SHUFFLE(3,3,3,3));
	  C = _mm_sub_ps(C,H);
	  C = _mm_mul_ps(C,C);
	  return ! (7&_mm_movemask_ps(_mm_cmpgt_ps(X,C)));
	}
	// inside(cuboid,sphere)
	static bool inside(box2 const&c, sph2 const&s)
	{
	  __m128 C = _mm_load_ps(cPF(c.X));
	  __m128 H = _mm_movehl_ps(C,C);
	  __m128 X = _mm_load_ps(cPF(s.X));
	  C = SSE::diff(C,X);
	  if(3&_mm_movemask_ps(_mm_cmpgt_ps(C,H)))
	    return false;
	  X = _mm_shuffle_ps(X,X,_MM_SHUFFLE(2,2,2,2));
	  C = _mm_sub_ps(C,H);
	  C = _mm_mul_ps(C,C);
	  return ! (3&_mm_movemask_ps(_mm_cmpgt_ps(X,C)));
	}
	static bool inside(box3 const&c, sph3 const&s)
	{
	  __m128 C = _mm_load_ps(cPF(c.X));
	  __m128 H = _mm_load_ps(cPF(c.Y));
	  __m128 X = _mm_load_ps(cPF(s.X));
	  C = SSE::diff(C,X);
	  if(7&_mm_movemask_ps(_mm_cmpgt_ps(C,H)))
	    return false;
	  X = _mm_shuffle_ps(X,X,_MM_SHUFFLE(3,3,3,3));
	  C = _mm_sub_ps(C,H);
	  C = _mm_mul_ps(C,C);
	  return ! (7&_mm_movemask_ps(_mm_cmpgt_ps(X,C)));
	}
	// convert2cuboid(PointPair)
	static void convert2cuboid(PointPair<2,float>&c)
	{
	  __m128 XY,YX;
	  XY =_mm_load_ps(cPF(c.X));
	  YX =_mm_shuffle_ps(XY,XY,_MM_SHUFFLE(1,0,3,2));
#ifdef __SSE2__
	  YX =_mm_xor_ps(YX,reinterpret_cast<__m128>
			 (_mm_set_epi32(0x80000000,0x80000000,0,0)));
#else
	  XY =_mm_xor_ps(YX,reinterpret_cast<__m128>
			 (_mm_set_ps(meta::__neg_mask.__F,
				     meta::__neg_mask.__F,0,0)));
#endif
	  _mm_store_ps(PF(c.X),_mm_mul_ps(_mm_set1_ps(0.5f),
					  _mm_add_ps(XY,YX)));
	}
	static void convert2cuboid(PointPair<3,float>&c)
	{
	  __m128 Xi = _mm_load_ps(cPF(c.X));
	  __m128 Xa = _mm_load_ps(cPF(c.Y));
	  __m128 h  = _mm_set1_ps(0.5f);
	  _mm_store_ps(PF(c.X),_mm_mul_ps(h,_mm_add_ps(Xi,Xa)));
	  _mm_store_ps(PF(c.Y),_mm_mul_ps(h,_mm_sub_ps(Xa,Xi)));
	}
      };// class AlgorithmsHelper<float,1,1>
#undef  PF
#undef cPF
#endif // __SSE__

#ifdef __SSE2__
#define  PD(__X) static_cast<double*>(__X)
#define cPD(__X) static_cast<const double*>(__X)
#define  PD2(__X) static_cast<double*>(__X)+2
#define cPD2(__X) static_cast<const double*>(__X)+2
      //
      // 2.3 specialisation for unaligned real=double with SSE
      //
      template<> class AlgorithmsHelper<double,0,1> {
	friend struct Algorithms<0,1>;
	typedef GeoVec<2,double> vec2;
	typedef GeoVec<3,double> vec3;
	typedef cube  <2,double> cub2;
	typedef cube  <3,double> cub3;
	typedef cuboid<2,double> box2;
	typedef cuboid<3,double> box3;
	typedef sphere<2,double> sph2;
	typedef sphere<3,double> sph3;
	// copy(cube,cube)
	static void copy(cub2 const&in, cub2 &out)
	{
	  out.X = in.X;
	  out.H = in.H;
	}
	static void copy(cub3 const&in, cub3 &out)
	{
	  _mm_store_pd(PD (out.X),_mm_load_pd(cPD (in.X)));
	  _mm_store_pd(PD2(out.X),_mm_load_pd(cPD2(in.X)));
	}
	// copy(sphere,sphere)
	static void copy(sph2 const&in, sph2 &out)
	{
	  out.X = in.X;
	  out.Q = in.Q;
	}
	static void copy(sph3 const&in, sph3 &out)
	{
	  _mm_store_pd(PD (out.X),_mm_load_pd(cPD (in.X)));
	  _mm_store_pd(PD2(out.X),_mm_load_pd(cPD2(in.X)));
	}
	// octant(point,point)
	static int octant(vec2 const&c, vec2 const&x)
	{
	  return _mm_movemask_pd(_mm_cmplt_pd(_mm_loadu_pd(cPD(c)),
					      _mm_loadu_pd(cPD(x))));
	}
	static int octant(vec3 const&c, vec3 const&x)
	{
	  return (  _mm_movemask_pd(_mm_cmplt_pd(_mm_loadu_pd(cPD(c)),
						 _mm_loadu_pd(cPD(x)))))
	    |   ((1&_mm_movemask_pd(_mm_cmplt_pd(_mm_loadu_pd(cPD2(c)),
						 _mm_loadu_pd(cPD2(x)))))<<2);
	}
	// contains(cube,point)
	static bool contains(cub2 const&c, vec2 const&x)
	{
	  __m128d C = _mm_loadu_pd(cPD(c.X));
	  __m128d X = _mm_loadu_pd(cPD(x));
	  __m128d H = _mm_set1_pd(c.H);
	  return
	    3 == _mm_movemask_pd(_mm_and_pd(_mm_cmple_pd(_mm_sub_pd(C,H),X),
					    _mm_cmpgt_pd(_mm_add_pd(C,H),X)));
	}
	static bool contains(cub3 const&c, vec3 const&x)
	{
	  __m128d C = _mm_loadu_pd(cPD(c.X));
	  __m128d X = _mm_loadu_pd(cPD(x));
	  __m128d H = _mm_set1_pd(c.H);
	  if(3!=_mm_movemask_pd(_mm_and_pd(_mm_cmple_pd(_mm_sub_pd(C,H),X),
					   _mm_cmpgt_pd(_mm_add_pd(C,H),X))))
	    return false;
	  C = _mm_loadu_pd(cPD2(c.X));
	  X = _mm_loadu_pd(cPD2(x));
	  return 1&_mm_movemask_pd(_mm_and_pd(_mm_cmple_pd(_mm_sub_pd(C,H),X),
					      _mm_cmpgt_pd(_mm_add_pd(C,H),X)));
	}
	// contains(cuboid,point)
	static bool contains(box2 const&c, vec2 const&x)
	{
	  __m128d C = _mm_loadu_pd(cPD(c.X));
	  __m128d X = _mm_loadu_pd(cPD(x));
	  __m128d H = _mm_loadu_pd(cPD(c.Y));
	  return
	    3 == _mm_movemask_pd(_mm_and_pd(_mm_cmple_pd(_mm_sub_pd(C,H),X),
					    _mm_cmpgt_pd(_mm_add_pd(C,H),X)));
	}
	static bool contains(box3 const&c, vec3 const&x)
	{
	  __m128d C = _mm_loadu_pd(cPD(c.X));
	  __m128d X = _mm_loadu_pd(cPD(x));
	  __m128d H = _mm_loadu_pd(cPD(c.Y));
	  if(3!=_mm_movemask_pd(_mm_and_pd(_mm_cmple_pd(_mm_sub_pd(C,H),X),
					   _mm_cmpgt_pd(_mm_add_pd(C,H),X))))
	    return false;
	  C = _mm_loadu_pd(cPD2(c.X));
	  X = _mm_loadu_pd(cPD2(x));
	  H = _mm_loadu_pd(cPD2(c.Y));
	  return 1&_mm_movemask_pd(_mm_and_pd(_mm_cmple_pd(_mm_sub_pd(C,H),X),
					      _mm_cmpgt_pd(_mm_add_pd(C,H),X)));
	}
	// contains(cuboid,cuboid)
	static bool contains(box2 const&c, box2 const&b)
	{
	  __m128d C = _mm_loadu_pd(cPD(c.X));
	  __m128d R = _mm_loadu_pd(cPD(c.Y));
	  __m128d A = _mm_loadu_pd(cPD(b.X));
	  __m128d B = _mm_loadu_pd(cPD(b.Y));
	  return
	    3==_mm_movemask_pd(_mm_cmple_pd(SSE::diff(C,A),_mm_sub_pd(R,B)));
	}
	static bool contains(box3 const&c, box3 const&b)
	{
	  __m128d C = _mm_loadu_pd(cPD(c.X));
	  __m128d R = _mm_loadu_pd(cPD(c.Y));
	  __m128d A = _mm_loadu_pd(cPD(b.X));
	  __m128d B = _mm_loadu_pd(cPD(b.Y));
	  if(3!=_mm_movemask_pd(_mm_cmple_pd(SSE::diff(C,A),_mm_sub_pd(R,B))))
	    return false;
	  C = _mm_loadu_pd(cPD2(c.X));
	  R = _mm_loadu_pd(cPD2(c.Y));
	  A = _mm_loadu_pd(cPD2(b.X));
	  B = _mm_loadu_pd(cPD2(b.Y));
	  return
	    1 &_mm_movemask_pd(_mm_cmple_pd(SSE::diff(C,A),_mm_sub_pd(R,B)));
	}
	// dist_sq(point,point)
	// for unaligned access, non-SSE code is faster
#if(0)
	static double dist_sq(vec2 const&x, vec2 const&y)
	{
	  WDutils__align16 double q[2];
	  __m128d D = _mm_sub_pd(_mm_loadu_pd(cPD(x)),_mm_loadu_pd(cPD(y)));
	  _mm_store_pd(q,_mm_mul_pd(D,D));
	  return q[0]+q[1];
	}
	static double dist_sq(vec3 const&x, vec3 const&y)
	{
	  WDutils__align16 double q[2],p[2];
	  __m128d D = _mm_sub_pd(_mm_loadu_pd(cPD(x)),_mm_loadu_pd(cPD(y)));
	  _mm_store_pd(q,_mm_mul_pd(D,D));
	  D = _mm_sub_pd(_mm_loadu_pd(cPD2(x)),_mm_loadu_pd(cPD2(y)));
	  _mm_store_pd(p,_mm_mul_pd(D,D));
	  return q[0]+q[1]+p[0];
	}
#else
	static double dist_sq(vec2 const&x, vec2 const&y)
	{ return WDutils::dist_sq(x,y); }
	static double dist_sq(vec3 const&x, vec3 const&y)
	{ return WDutils::dist_sq(x,y); }
#endif
	// dist_sq(cube,point)
	static double dist_sq(cub2 const&c, vec2 const&x)
	{
	  WDutils__align16 double q[2];
	  __m128d C = _mm_loadu_pd(cPD(c.X));
	  __m128d X = _mm_loadu_pd(cPD(x));
	  __m128d H = _mm_set1_pd(c.H);
	  X = _mm_sub_pd(SSE::diff(C,X),H);
	  _mm_store_pd(q,_mm_and_pd(_mm_mul_pd(X,X),
				    _mm_cmpgt_pd(X,_mm_setzero_pd())));
	  return q[0]+q[1];
	}
	static double dist_sq(cub3 const&c, vec3 const&x)
	{
	  __m128d C,X,H;
	  C = _mm_loadu_pd(cPD(c.X));
	  X = _mm_loadu_pd(cPD(x));
	  H = _mm_set1_pd(c.H);
	  C = _mm_sub_pd(SSE::diff(C,X),H);
	  WDutils__align16 double q[2];
	  _mm_store_pd(q,_mm_and_pd(_mm_mul_pd(C,C),
				    _mm_cmpgt_pd(C,_mm_setzero_pd())));
	  register double d=q[0]+q[1];
	  C = _mm_loadu_pd(cPD2(c.X));
	  X = _mm_loadu_pd(cPD2(x));
	  C = _mm_sub_pd(SSE::diff(C,X),H);
	  _mm_store_pd(q,_mm_and_pd(_mm_mul_pd(C,C),
				    _mm_cmpgt_pd(C,_mm_setzero_pd())));
	  return d+=q[0];
	}
	// dist_sq(cuboid,point)
	static double dist_sq(box2 const&c, vec2 const&x)
	{
	  WDutils__align16 double q[2];
	  __m128d C = _mm_loadu_pd(cPD(c.X));
	  __m128d X = _mm_loadu_pd(cPD(x));
	  __m128d H = _mm_loadu_pd(cPD(c.Y));
	  X = _mm_sub_pd(SSE::diff(C,X),H);
	  _mm_store_pd(q,_mm_and_pd(_mm_mul_pd(X,X),
				    _mm_cmpgt_pd(X,_mm_setzero_pd())));
	  return q[0]+q[1];
	}
	static double dist_sq(box3 const&c, vec3 const&x)
	{
	  __m128d C,X,H;
	  C = _mm_loadu_pd(cPD(c.X));
	  X = _mm_loadu_pd(cPD(x));
	  H = _mm_loadu_pd(cPD(c.Y));
	  C = _mm_sub_pd(SSE::diff(C,X),H);
	  WDutils__align16 double q[2];
	  _mm_store_pd(q,_mm_and_pd(_mm_mul_pd(C,C),
				    _mm_cmpgt_pd(C,_mm_setzero_pd())));
	  register double d=q[0]+q[1];
	  C = _mm_loadu_pd(cPD2(c.X));
	  X = _mm_loadu_pd(cPD2(x));
	  H = _mm_loadu_pd(cPD2(c.Y));
	  C = _mm_sub_pd(SSE::diff(C,X),H);
	  _mm_store_pd(q,_mm_and_pd(_mm_mul_pd(C,C),
				    _mm_cmpgt_pd(C,_mm_setzero_pd())));
	  return d+=q[0];
	}
	// outside(cube,sphere)
	static bool outside(cub2 const&c, sph2 const&s)
	{
	  return dist_sq(c,s.X) > s.Q;
	}
	static bool outside(cub3 const&c, sph3 const&s)
	{
	  __m128d C,X,H,Q;
	  C = _mm_loadu_pd(cPD(c.X));
	  X = _mm_loadu_pd(cPD(s.X));
	  H = _mm_set1_pd(c.H);
	  X = _mm_sub_pd(SSE::diff(C,X),H);
	  Q = _mm_and_pd(_mm_mul_pd(X,X),_mm_cmpgt_pd(X,_mm_setzero_pd()));
	  WDutils__align16 double q[2];
	  _mm_store_pd(q,Q);
	  double d=q[0]+q[1];
	  if(d>s.Q) return true;
	  C = _mm_loadu_pd(cPD2(c.X));
	  X = _mm_loadu_pd(cPD2(s.X));
	  X = _mm_sub_pd(SSE::diff(C,X),H);
	  Q = _mm_and_pd(_mm_mul_pd(X,X),_mm_cmpgt_pd(X,_mm_setzero_pd()));
	  _mm_store_pd(q,Q);
	  return (d+=q[0]) > s.Q;
	}
	// outside(cuboid,sphere)
	static bool outside(box2 const&c, sph2 const&s)
	{
	  return dist_sq(c,s.X) > s.Q;
	}
	static bool outside(box3 const&c, sph3 const&s)
	{
	  __m128d C,X,H,Q;
	  C = _mm_loadu_pd(cPD(c.X));
	  X = _mm_loadu_pd(cPD(s.X));
	  H = _mm_loadu_pd(cPD(c.Y));
	  X = _mm_sub_pd(SSE::diff(C,X),H);
	  Q = _mm_and_pd(_mm_mul_pd(X,X),_mm_cmpgt_pd(X,_mm_setzero_pd()));
	  WDutils__align16 double q[2];
	  _mm_store_pd(q,Q);
	  double d=q[0]+q[1];
	  if(d>s.Q) return true;
	  C = _mm_loadu_pd(cPD2(c.X));
	  X = _mm_loadu_pd(cPD2(s.X));
	  H = _mm_loadu_pd(cPD2(c.Y));
	  X = _mm_sub_pd(SSE::diff(C,X),H);
	  Q = _mm_and_pd(_mm_mul_pd(X,X),_mm_cmpgt_pd(X,_mm_setzero_pd()));
	  _mm_store_pd(q,Q);
	  return (d+=q[0]) > s.Q;
	}
	// inside(cube,sphere)
	static bool inside(cub2 const&c, sph2 const&s)
	{
	  __m128d C = _mm_loadu_pd(cPD(c.X));
	  __m128d H = _mm_set1_pd(c.H);
	  __m128d X = _mm_loadu_pd(cPD(s.X));
	  C = SSE::diff(C,X);
	  if(_mm_movemask_pd(_mm_cmpgt_pd(C,H))) return false;
	  X = _mm_set1_pd(s.Q);
	  C = _mm_sub_pd(C,H);
	  C = _mm_mul_pd(C,C);
	  return ! _mm_movemask_pd(_mm_cmpgt_pd(X,C));
	}
	static bool inside(cub3 const&c, sph3 const&s)
	{
	  __m128d C = _mm_loadu_pd(cPD(c.X));
	  __m128d H = _mm_set1_pd(c.H);
	  __m128d X = _mm_loadu_pd(cPD(s.X));
	  C = SSE::diff(C,X);
	  if(_mm_movemask_pd(_mm_cmpgt_pd(C,H))) return false;
	  __m128d Q = _mm_set1_pd(s.Q);
	  C = _mm_sub_pd(C,H);
	  C = _mm_mul_pd(C,C);
	  if(_mm_movemask_pd(_mm_cmpgt_pd(Q,C))) return false;
	  C = _mm_loadu_pd(cPD2(c.X));
	  X = _mm_loadu_pd(cPD2(s.X));
	  C = SSE::diff(C,X);
	  if(1&_mm_movemask_pd(_mm_cmpgt_pd(C,H))) return false;
	  C = _mm_sub_pd(C,H);
	  C = _mm_mul_pd(C,C);
	  return ! (1&_mm_movemask_pd(_mm_cmpgt_pd(Q,C)));
	}
	// inside(cuboid,sphere)
	static bool inside(box2 const&c, sph2 const&s)
	{
	  __m128d C = _mm_loadu_pd(cPD(c.X));
	  __m128d H = _mm_loadu_pd(cPD(c.Y));
	  __m128d X = _mm_loadu_pd(cPD(s.X));
	  C = SSE::diff(C,X);
	  if(_mm_movemask_pd(_mm_cmpgt_pd(C,H))) return false;
	  X = _mm_set1_pd(s.Q);
	  C = _mm_sub_pd(C,H);
	  C = _mm_mul_pd(C,C);
	  return ! _mm_movemask_pd(_mm_cmpgt_pd(X,C));
	}
	static bool inside(box3 const&c, sph3 const&s)
	{
	  __m128d C = _mm_loadu_pd(cPD(c.X));
	  __m128d H = _mm_loadu_pd(cPD(c.Y));
	  __m128d X = _mm_loadu_pd(cPD(s.X));
	  C = SSE::diff(C,X);
	  if(_mm_movemask_pd(_mm_cmpgt_pd(C,H))) return false;
	  __m128d Q = _mm_set1_pd(s.Q);
	  C = _mm_sub_pd(C,H);
	  C = _mm_mul_pd(C,C);
	  if(_mm_movemask_pd(_mm_cmpgt_pd(Q,C))) return false;
	  C = _mm_loadu_pd(cPD2(c.X));
	  X = _mm_loadu_pd(cPD2(s.X));
	  H = _mm_loadu_pd(cPD2(c.Y));
	  C = SSE::diff(C,X);
	  if(1&_mm_movemask_pd(_mm_cmpgt_pd(C,H))) return false;
	  C = _mm_sub_pd(C,H);
	  C = _mm_mul_pd(C,C);
	  return ! (1&_mm_movemask_pd(_mm_cmpgt_pd(Q,C)));
	}
	// convert2cuboid(PointPair)
#if(0)  // unaligned SSE code is slower than non-SSE code
	static void convert2cuboid(PointPair<2,double>&c)
	{
	  __m128d Xi = _mm_loadu_pd(cPD(c.X));
	  __m128d Xa = _mm_loadu_pd(cPD(c.Y));
	  __m128d h  = _mm_set1_pd(0.5);
	  _mm_storeu_pd(PD(c.X),_mm_mul_pd(h,_mm_add_pd(Xi,Xa)));
	  _mm_storeu_pd(PD(c.Y),_mm_mul_pd(h,_mm_sub_pd(Xa,Xi)));
	}
	static void convert2cuboid(PointPair<3,double>&c)
	{
	  
	  __m128d Xi = _mm_loadu_pd(cPD(c.X));
	  __m128d Xa = _mm_loadu_pd(cPD(c.Y));
	  __m128d h  = _mm_set1_pd(0.5);
	  _mm_storeu_pd(PD(c.X),_mm_mul_pd(h,_mm_add_pd(Xi,Xa)));
	  _mm_storeu_pd(PD(c.Y),_mm_mul_pd(h,_mm_sub_pd(Xa,Xi)));
	  Xi = _mm_loadu_pd(cPD2(c.X));
	  Xa = _mm_loadu_pd(cPD2(c.Y));
	  _mm_storeu_pd(PD2(c.X),_mm_mul_pd(h,_mm_add_pd(Xi,Xa)));
	  _mm_storeu_pd(PD2(c.Y),_mm_mul_pd(h,_mm_sub_pd(Xa,Xi)));
	}
#else
	static void convert2cuboid(PointPair<2,double>&c)
	{ return AlgorithmsHelper<double,0,0>::convert2cuboid(c); }
	static void convert2cuboid(PointPair<3,double>&c)
	{ return AlgorithmsHelper<double,0,0>::convert2cuboid(c); }
#endif
      }; // class AlgorithmsHelper<double,0,1>
      //
      // 2.4 specialisation for aligned real=double with SSE
      //
      template<> class AlgorithmsHelper<double,1,1> {
	friend struct Algorithms<1,1>;
	typedef GeoVec<2,double> vec2;
	typedef GeoVec<3,double> vec3;
	typedef cube  <2,double> cub2;
	typedef cube  <3,double> cub3;
	typedef cuboid<2,double> box2;
	typedef cuboid<3,double> box3;
	typedef sphere<2,double> sph2;
	typedef sphere<3,double> sph3;
	// copy(cube,cube)
	static void copy(cub2 const&in, cub2 &out)
	{
	  out.X = in.X;
	  out.H = in.H;
	}
	static void copy(cub3 const&in, cub3 &out)
	{
	  _mm_store_pd(PD (out.X),_mm_load_pd(cPD (in.X)));
	  _mm_store_pd(PD2(out.X),_mm_load_pd(cPD2(in.X)));
	}
	// copy(sphere,sphere)
	static void copy(sph2 const&in, sph2 &out)
	{
	  out.X = in.X;
	  out.Q = in.Q;
	}
	static void copy(sph3 const&in, sph3 &out)
	{
	  _mm_store_pd(PD (out.X),_mm_load_pd(cPD (in.X)));
	  _mm_store_pd(PD2(out.X),_mm_load_pd(cPD2(in.X)));
	}
	// octant(point,point)
	static int octant(vec2 const&c, vec2 const&x)
	{
	  return 3&_mm_movemask_pd(_mm_cmplt_pd(_mm_load_pd(cPD(c)),
						_mm_load_pd(cPD(x))));
	}
	static int octant(vec3 const&c, vec3 const&x)
	{
	  return (3&_mm_movemask_pd(_mm_cmplt_pd(_mm_load_pd(cPD(c)),
						 _mm_load_pd(cPD(x)))))
	    |   ((1&_mm_movemask_pd(_mm_cmplt_pd(_mm_load_pd(cPD2(c)),
						 _mm_load_pd(cPD2(x)))))<<2);
	}
	// contains(cube,point)
	static bool contains(cub2 const&c, vec2 const&x)
	{
	  __m128d C = _mm_load_pd(cPD(c.X));
	  __m128d X = _mm_load_pd(cPD(x));
	  __m128d H = _mm_set1_pd(c.H);
	  return
	    3 == _mm_movemask_pd(_mm_and_pd(_mm_cmple_pd(_mm_sub_pd(C,H),X),
					    _mm_cmpgt_pd(_mm_add_pd(C,H),X)));
	}
	static bool contains(cub3 const&c, vec3 const&x)
	{
	  __m128d C = _mm_load_pd(cPD(c.X));
	  __m128d X = _mm_load_pd(cPD(x));
	  __m128d H = _mm_set1_pd(c.H);
	  if(3!=_mm_movemask_pd(_mm_and_pd(_mm_cmple_pd(_mm_sub_pd(C,H),X),
					   _mm_cmpgt_pd(_mm_add_pd(C,H),X))))
	    return false;
	  C = _mm_load_pd(cPD2(c.X));
	  X = _mm_load_pd(cPD2(x));
	  return 1&_mm_movemask_pd(_mm_and_pd(_mm_cmple_pd(_mm_sub_pd(C,H),X),
					      _mm_cmpgt_pd(_mm_add_pd(C,H),X)));
	}
	// contains(cuboid,point)
	static bool contains(box2 const&c, vec2 const&x)
	{
	  __m128d C = _mm_load_pd(cPD(c.X));
	  __m128d X = _mm_load_pd(cPD(x));
	  __m128d H = _mm_load_pd(cPD(c.Y));
	  return
	    3 == _mm_movemask_pd(_mm_and_pd(_mm_cmple_pd(_mm_sub_pd(C,H),X),
					    _mm_cmpgt_pd(_mm_add_pd(C,H),X)));
	}
	static bool contains(box3 const&c, vec3 const&x)
	{
	  __m128d C = _mm_load_pd(cPD(c.X));
	  __m128d X = _mm_load_pd(cPD(x));
	  __m128d H = _mm_load_pd(cPD(c.Y));
	  if(3!=_mm_movemask_pd(_mm_and_pd(_mm_cmple_pd(_mm_sub_pd(C,H),X),
					   _mm_cmpgt_pd(_mm_add_pd(C,H),X))))
	    return false;
	  C = _mm_load_pd(cPD2(c.X));
	  X = _mm_load_pd(cPD2(x));
	  H = _mm_load_pd(cPD2(c.Y));
	  return 1&_mm_movemask_pd(_mm_and_pd(_mm_cmple_pd(_mm_sub_pd(C,H),X),
					      _mm_cmpgt_pd(_mm_add_pd(C,H),X)));
	}
	// contains(cuboid,cuboid)
	static bool contains(box2 const&c, box2 const&b)
	{
	  __m128d C = _mm_load_pd(cPD(c.X));
	  __m128d R = _mm_load_pd(cPD(c.Y));
	  __m128d A = _mm_load_pd(cPD(b.X));
	  __m128d B = _mm_load_pd(cPD(b.Y));
	  return
	    3==_mm_movemask_pd(_mm_cmple_pd(SSE::diff(C,A),_mm_sub_pd(R,B)));
	}
	static bool contains(box3 const&c, box3 const&b)
	{
	  __m128d C = _mm_load_pd(cPD(c.X));
	  __m128d R = _mm_load_pd(cPD(c.Y));
	  __m128d A = _mm_load_pd(cPD(b.X));
	  __m128d B = _mm_load_pd(cPD(b.Y));
	  if(3!=_mm_movemask_pd(_mm_cmple_pd(SSE::diff(C,A),_mm_sub_pd(R,B))))
	    return false;
	  C = _mm_load_pd(cPD2(c.X));
	  R = _mm_load_pd(cPD2(c.Y));
	  A = _mm_load_pd(cPD2(b.X));
	  B = _mm_load_pd(cPD2(b.Y));
	  return
	    1 &_mm_movemask_pd(_mm_cmple_pd(SSE::diff(C,A),_mm_sub_pd(R,B)));
	}
	// dist_sq(point,point)
	static double dist_sq(vec2 const&x, vec2 const&y)
	{
	  WDutils__align16 double q[2];
	  __m128d D = _mm_sub_pd(_mm_load_pd(cPD(x)),_mm_load_pd(cPD(y)));
	  _mm_store_pd(q,_mm_mul_pd(D,D));
	  return q[0]+q[1];
	}
	static double dist_sq(vec3 const&x, vec3 const&y)
	{
	  WDutils__align16 double q[2],p[2];
	  __m128d D = _mm_sub_pd(_mm_load_pd(cPD(x)),_mm_load_pd(cPD(y)));
	  _mm_store_pd(q,_mm_mul_pd(D,D));
	  D = _mm_sub_pd(_mm_load_pd(cPD2(x)),_mm_load_pd(cPD2(y)));
	  _mm_store_pd(p,_mm_mul_pd(D,D));
	  return q[0]+q[1]+p[0];
	}
	// dist_sq(cube,point)
	static double dist_sq(cub2 const&c, vec2 const&x)
	{
	  WDutils__align16 double q[2];
	  __m128d C = _mm_load_pd(cPD(c.X));
	  __m128d X = _mm_load_pd(cPD(x));
	  __m128d H = _mm_set1_pd(c.H);
	  X = _mm_sub_pd(SSE::diff(C,X),H);
	  _mm_store_pd(q,_mm_and_pd(_mm_mul_pd(X,X),
				    _mm_cmpgt_pd(X,_mm_setzero_pd())));
	  return q[0]+q[1];
	}
	static double dist_sq(cub3 const&c, vec3 const&x)
	{
	  __m128d C,X,H,Q;
	  C = _mm_load_pd(cPD(c.X));
	  X = _mm_load_pd(cPD(x));
	  H = _mm_set1_pd(c.H);
	  X = _mm_sub_pd(SSE::diff(C,X),H);
	  Q = _mm_and_pd(_mm_mul_pd(X,X),_mm_cmpgt_pd(X,_mm_setzero_pd()));
	  WDutils__align16 double q[2];
	  _mm_store_pd(q,Q);
	  double d=q[0]+q[1];
	  C = _mm_load_pd(cPD2(c.X));
	  X = _mm_load_pd(cPD2(x));
	  X = _mm_sub_pd(SSE::diff(C,X),H);
	  Q = _mm_and_pd(_mm_mul_pd(X,X),_mm_cmpgt_pd(X,_mm_setzero_pd()));
	  _mm_store_pd(q,Q);
	  return d+=q[0];
	}
	// dist_sq(cuboid,point)
	static double dist_sq(box2 const&c, vec2 const&x)
	{
	  WDutils__align16 double q[2];
	  __m128d C = _mm_load_pd(cPD(c.X));
	  __m128d X = _mm_load_pd(cPD(x));
	  __m128d H = _mm_load_pd(cPD(c.Y));
	  X = _mm_sub_pd(SSE::diff(C,X),H);
	  _mm_store_pd(q,_mm_and_pd(_mm_mul_pd(X,X),
				    _mm_cmpgt_pd(X,_mm_setzero_pd())));
	  return q[0]+q[1];
	}
	static double dist_sq(box3 const&c, vec3 const&x)
	{
	  __m128d C,X,H;
	  C = _mm_load_pd(cPD(c.X));
	  X = _mm_load_pd(cPD(x));
	  H = _mm_load_pd(cPD(c.Y));
	  C = _mm_sub_pd(SSE::diff(C,X),H);
	  WDutils__align16 double q[2];
	  _mm_store_pd(q,_mm_and_pd(_mm_mul_pd(C,C),
				    _mm_cmpgt_pd(C,_mm_setzero_pd())));
	  register double d=q[0]+q[1];
	  C = _mm_load_pd(cPD2(c.X));
	  X = _mm_load_pd(cPD2(x));
	  H = _mm_load_pd(cPD2(c.Y));
	  C = _mm_sub_pd(SSE::diff(C,X),H);
	  _mm_store_pd(q,_mm_and_pd(_mm_mul_pd(C,C),
				    _mm_cmpgt_pd(C,_mm_setzero_pd())));
	  return d+=q[0];
	}
	// outside(cube,sphere)
	static bool outside(cub2 const&c, sph2 const&s)
	{
	  return dist_sq(c,s.X) > s.Q;
	}
	static bool outside(cub3 const&c, sph3 const&s)
	{
	  __m128d C,X,H,Q;
	  C = _mm_load_pd(cPD(c.X));
	  X = _mm_load_pd(cPD(s.X));
	  H = _mm_set1_pd(c.H);
	  X = _mm_sub_pd(SSE::diff(C,X),H);
	  Q = _mm_and_pd(_mm_mul_pd(X,X),_mm_cmpgt_pd(X,_mm_setzero_pd()));
	  WDutils__align16 double q[2];
	  _mm_store_pd(q,Q);
	  double d=q[0]+q[1];
	  if(d>s.Q) return true;
	  C = _mm_load_pd(cPD2(c.X));
	  X = _mm_load_pd(cPD2(s.X));
	  X = _mm_sub_pd(SSE::diff(C,X),H);
	  Q = _mm_and_pd(_mm_mul_pd(X,X),_mm_cmpgt_pd(X,_mm_setzero_pd()));
	  _mm_store_pd(q,Q);
	  return (d+=q[0]) > s.Q;
	}
	// outside(cuboid,sphere)
	static bool outside(box2 const&c, sph2 const&s)
	{
	  return dist_sq(c,s.X) > s.Q;
	}
	static bool outside(box3 const&c, sph3 const&s)
	{
	  __m128d C,X,H,Q;
	  C = _mm_load_pd(cPD(c.X));
	  X = _mm_load_pd(cPD(s.X));
	  H = _mm_load_pd(cPD(c.Y));
	  X = _mm_sub_pd(SSE::diff(C,X),H);
	  Q = _mm_and_pd(_mm_mul_pd(X,X),_mm_cmpgt_pd(X,_mm_setzero_pd()));
	  WDutils__align16 double q[2];
	  _mm_store_pd(q,Q);
	  double d=q[0]+q[1];
	  if(d>s.Q) return true;
	  C = _mm_load_pd(cPD2(c.X));
	  X = _mm_load_pd(cPD2(s.X));
	  H = _mm_load_pd(cPD2(c.Y));
	  X = _mm_sub_pd(SSE::diff(C,X),H);
	  Q = _mm_and_pd(_mm_mul_pd(X,X),_mm_cmpgt_pd(X,_mm_setzero_pd()));
	  _mm_store_pd(q,Q);
	  return (d+=q[0]) > s.Q;
	}
	// inside(cube,sphere)
	static bool inside(cub2 const&c, sph2 const&s)
	{
	  __m128d C = _mm_load_pd(cPD(c.X));
	  __m128d H = _mm_set1_pd(c.H);
	  __m128d X = _mm_load_pd(cPD(s.X));
	  C = SSE::diff(C,X);
	  if(_mm_movemask_pd(_mm_cmpgt_pd(C,H))) return false;
	  X = _mm_set1_pd(s.Q);
	  C = _mm_sub_pd(C,H);
	  C = _mm_mul_pd(C,C);
	  return ! _mm_movemask_pd(_mm_cmpgt_pd(X,C));
	}
	static bool inside(cub3 const&c, sph3 const&s)
	{
	  __m128d C = _mm_load_pd(cPD(c.X));
	  __m128d H = _mm_set1_pd(c.H);
	  __m128d X = _mm_load_pd(cPD(s.X));
	  C = SSE::diff(C,X);
	  if( _mm_movemask_pd(_mm_cmpgt_pd(C,H))) return false;
	  __m128d Q = _mm_set1_pd(s.Q);
	  C = _mm_sub_pd(C,H);
	  C = _mm_mul_pd(C,C);
	  if( _mm_movemask_pd(_mm_cmpgt_pd(Q,C))) return false;
	  C = _mm_load_pd(cPD2(c.X));
	  X = _mm_load_pd(cPD2(s.X));
	  C = SSE::diff(C,X);
	  if(1 & _mm_movemask_pd(_mm_cmpgt_pd(C,H)) ) return false;
	  C = _mm_sub_pd(C,H);
	  C = _mm_mul_pd(C,C);
	  return ! (1&_mm_movemask_pd(_mm_cmpgt_pd(Q,C)));
	}
	// inside(cuboid,sphere)
	static bool inside(box2 const&c, sph2 const&s)
	{
	  __m128d C = _mm_load_pd(cPD(c.X));
	  __m128d H = _mm_load_pd(cPD(c.Y));
	  __m128d X = _mm_load_pd(cPD(s.X));
	  C = SSE::diff(C,X);
	  if(_mm_movemask_pd(_mm_cmpgt_pd(C,H))) return false;
	  X = _mm_set1_pd(s.Q);
	  C = _mm_sub_pd(C,H);
	  C = _mm_mul_pd(C,C);
	  return ! _mm_movemask_pd(_mm_cmpgt_pd(X,C));
	}
	static bool inside(box3 const&c, sph3 const&s)
	{
	  __m128d C = _mm_load_pd(cPD(c.X));
	  __m128d H = _mm_load_pd(cPD(c.Y));
	  __m128d X = _mm_load_pd(cPD(s.X));
	  C = SSE::diff(C,X);
	  if(_mm_movemask_pd(_mm_cmpgt_pd(C,H))) return false;
	  __m128d Q = _mm_set1_pd(s.Q);
	  C = _mm_sub_pd(C,H);
	  C = _mm_mul_pd(C,C);
	  if(_mm_movemask_pd(_mm_cmpgt_pd(Q,C))) return false;
	  C = _mm_load_pd(cPD2(c.X));
	  X = _mm_load_pd(cPD2(s.X));
	  H = _mm_load_pd(cPD2(c.Y));
	  C = SSE::diff(C,X);
	  if(1&_mm_movemask_pd(_mm_cmpgt_pd(C,H))) return false;
	  C = _mm_sub_pd(C,H);
	  C = _mm_mul_pd(C,C);
	  return ! (1&_mm_movemask_pd(_mm_cmpgt_pd(Q,C)));
	}
	// convert2cuboid(PointPair)
	static void convert2cuboid(PointPair<2,double>&c)
	{
	  __m128d Xi = _mm_load_pd(cPD(c.X));
	  __m128d Xa = _mm_load_pd(cPD(c.Y));
	  __m128d h  = _mm_set1_pd(0.5);
	  _mm_store_pd(PD(c.X),_mm_mul_pd(h,_mm_add_pd(Xi,Xa)));
	  _mm_store_pd(PD(c.Y),_mm_mul_pd(h,_mm_sub_pd(Xa,Xi)));
	}
	static void convert2cuboid(PointPair<3,double>&c)
	{
	  
	  __m128d Xi = _mm_load_pd(cPD(c.X));
	  __m128d Xa = _mm_load_pd(cPD(c.Y));
	  __m128d h  = _mm_set1_pd(0.5);
	  _mm_store_pd(PD(c.X),_mm_mul_pd(h,_mm_add_pd(Xi,Xa)));
	  _mm_store_pd(PD(c.Y),_mm_mul_pd(h,_mm_sub_pd(Xa,Xi)));
	  Xi = _mm_load_pd(cPD2(c.X));
	  Xa = _mm_load_pd(cPD2(c.Y));
	  _mm_store_pd(PD2(c.X),_mm_mul_pd(h,_mm_add_pd(Xi,Xa)));
	  _mm_store_pd(PD2(c.Y),_mm_mul_pd(h,_mm_sub_pd(Xa,Xi)));
	}
      }; // class AlgorithmsHelper<double,1,1>
#undef  PD
#undef cPD
#undef  PD2
#undef cPD2
#endif // __SSE2__
    } // namespace WDutils::Geometry::Meta
    //
    template<bool __A, bool __S> template<int __D, typename __X> inline
    void Algorithms<__A,__S>::copy(cube<__D,__X> const&in,
				   cube<__D,__X> &out)
    { return Meta::AlgorithmsHelper<__X,__A,__S>::copy(in,out); }
    //
    template<bool __A, bool __S> template<int __D, typename __X> inline
    void Algorithms<__A,__S>::copy(sphere<__D,__X> const&in,
				   sphere<__D,__X> &out)
    { return Meta::AlgorithmsHelper<__X,__A,__S>::copy(in,out); }
    //
    template<bool __A, bool __S> template<int __D, typename __X> inline
    void Algorithms<__A,__S>::move_to_octant(GeoVec<__D,__X>&c, int i, __X r)
    { return Meta::AlgorithmsHelper<__X,0,0>::move_to_octant(c,i,r); }
    //
    template<bool __A, bool __S> template<int __D, typename __X> inline
    int Algorithms<__A,__S>::octant(GeoVec<__D,__X> const&c,
				    GeoVec<__D,__X> const&x)
    { return Meta::AlgorithmsHelper<__X,__A,__S>::octant(c,x); }
    //
    template<bool __A, bool __S> template<int __D, typename __X> inline
    __X Algorithms<__A,__S>::dist_sq(GeoVec<__D,__X> const&x,
				     GeoVec<__D,__X> const&y)
    { return Meta::AlgorithmsHelper<__X,__A,__S>::dist_sq(x,y); }
    //
    template<bool __A, bool __S> template<int __D, typename __X> inline
    bool Algorithms<__A,__S>::contains(cube  <__D,__X> const&c,
				       GeoVec<__D,__X> const&x)
    { return Meta::AlgorithmsHelper<__X,__A,__S>::contains(c,x); }
    //
    template<bool __A, bool __S> template<int __D, typename __X> inline
    __X Algorithms<__A,__S>::dist_sq(cube  <__D,__X> const&c,
				     GeoVec<__D,__X> const&x)
    { return Meta::AlgorithmsHelper<__X,__A,__S>::dist_sq(c,x); }
    //
    template<bool __A, bool __S> template<int __D, typename __X> inline
    bool Algorithms<__A,__S>::outside(cube<__D,__X> const&c,
				      sphere<__D,__X> const&s)
    { return Meta::AlgorithmsHelper<__X,__A,__S>::outside(c,s); }
    //
    template<bool __A, bool __S> template<int __D, typename __X> inline
    bool Algorithms<__A,__S>::inside(cube<__D,__X> const&c,
				     sphere<__D,__X> const&s)
    { return Meta::AlgorithmsHelper<__X,__A,__S>::inside(c,s); }
    //
    template<bool __A, bool __S> template<int __D, typename __X> inline
    bool Algorithms<__A,__S>::contains(cuboid<__D,__X> const&c,
				       GeoVec<__D,__X> const&x)
    { return Meta::AlgorithmsHelper<__X,__A,__S>::contains(c,x); }
    //
    template<bool __A, bool __S> template<int __D, typename __X> inline
    bool Algorithms<__A,__S>::contains(cuboid<__D,__X> const&c,
				       cuboid<__D,__X> const&b)
    { return Meta::AlgorithmsHelper<__X,__A,__S>::contains(c,b); }
    //
    template<bool __A, bool __S> template<int __D, typename __X> inline
    __X Algorithms<__A,__S>::dist_sq(cuboid<__D,__X> const&c,
				     GeoVec<__D,__X> const&x)
    { return Meta::AlgorithmsHelper<__X,__A,__S>::dist_sq(c,x); }
    //
    template<bool __A, bool __S> template<int __D, typename __X> inline
    bool Algorithms<__A,__S>::outside(cuboid<__D,__X> const&c,
				      sphere<__D,__X> const&s)
    { return Meta::AlgorithmsHelper<__X,__A,__S>::outside(c,s); }
    //
    template<bool __A, bool __S> template<int __D, typename __X> inline
    bool Algorithms<__A,__S>::inside(cuboid<__D,__X> const&c,
				     sphere<__D,__X> const&s)
    { return Meta::AlgorithmsHelper<__X,__A,__S>::inside(c,s); }
    //
    template<bool __A, bool __S> template<int __D, typename __X> inline
    void Algorithms<__A,__S>::convert2cuboid(PointPair<__D,__X> &p)
    { return Meta::AlgorithmsHelper<__X,__A,__S>::convert2cuboid(p); }
  } // nameapce WDutils::Geometry
} // namespace WDutils
//
#endif
