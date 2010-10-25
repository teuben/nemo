/*
     PLIB - A Suite of Portable Game Libraries
     Copyright (C) 1998,2002  Steve Baker
 
     This library is free software; you can redistribute it and/or
     modify it under the terms of the GNU Library General Public
     License as published by the Free Software Foundation; either
     version 2 of the License, or (at your option) any later version.
 
     This library is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     Library General Public License for more details.
 
     You should have received a copy of the GNU Library General Public
     License along with this library; if not, write to the Free Software
     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
 
     For further information visit http://plib.sourceforge.net

     $Id$
*/


#include "sg.h"

sgdVec3 _sgdGravity = { 0.0f, 0.0f, -9.8f } ;

void sgdVectorProductVec3 ( sgdVec3 dst, const sgdVec3 a, const sgdVec3 b )
{
  dst[0] = a[1] * b[2] - a[2] * b[1] ;
  dst[1] = a[2] * b[0] - a[0] * b[2] ;
  dst[2] = a[0] * b[1] - a[1] * b[0] ;
}

inline SGDfloat _sgdClampToUnity ( const SGDfloat x )
{
  if ( x >  SGD_ONE ) return  SGD_ONE ;
  if ( x < -SGD_ONE ) return -SGD_ONE ;
  return x ;
}

int sgdCompare3DSqdDist( const sgdVec3 v1, const sgdVec3 v2, const SGDfloat sqd_dist )
{
  sgdVec3 tmp ;

  sgdSubVec3 ( tmp, v2, v1 ) ;

  SGDfloat sqdist = tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2] ;

  if ( sqdist > sqd_dist ) return  1 ;
  if ( sqdist < sqd_dist ) return -1 ;
  return 0 ;
}

void sgdMakeRotMat4( sgdMat4 mat, const SGDfloat angle, const sgdVec3 axis )
{
  sgdVec3 ax ;
  sgdNormalizeVec3 ( ax, axis ) ; 

  SGDfloat temp_angle = angle * SGD_DEGREES_TO_RADIANS ;
  SGDfloat s = sin ( temp_angle ) ;
  SGDfloat c = cos ( temp_angle ) ;
  SGDfloat t = SGD_ONE - c ;
   
  mat[0][0] = t * ax[0] * ax[0] + c ;
  mat[0][1] = t * ax[0] * ax[1] + s * ax[2] ;
  mat[0][2] = t * ax[0] * ax[2] - s * ax[1] ;
  mat[0][3] = SGD_ZERO ;

  mat[1][0] = t * ax[1] * ax[0] - s * ax[2] ;
  mat[1][1] = t * ax[1] * ax[1] + c ;
  mat[1][2] = t * ax[1] * ax[2] + s * ax[0] ;
  mat[1][3] = SGD_ZERO ;

  mat[2][0] = t * ax[2] * ax[0] + s * ax[1] ;
  mat[2][1] = t * ax[2] * ax[1] - s * ax[0] ;
  mat[2][2] = t * ax[2] * ax[2] + c ;
  mat[2][3] = SGD_ZERO ;

  mat[3][0] = SGD_ZERO ;
  mat[3][1] = SGD_ZERO ;
  mat[3][2] = SGD_ZERO ;
  mat[3][3] = SGD_ONE ;
}




void sgdMakePickMatrix( sgdMat4 mat, sgdFloat x, sgdFloat y,
                       sgdFloat width, sgdFloat height, sgdVec4 viewport )
{
   sgdFloat sx =   viewport[2] / width  ;
   sgdFloat sy =   viewport[3] / height ;
   sgdFloat tx = ( viewport[2] + SGD_TWO * (viewport[0] - x) ) / width  ;
   sgdFloat ty = ( viewport[3] + SGD_TWO * (viewport[1] - y) ) / height ;
 
   mat[0][0] =    sx   ;
   mat[0][1] = SGD_ZERO ;
   mat[0][2] = SGD_ZERO ;
   mat[0][3] = SGD_ZERO ;

   mat[1][0] = SGD_ZERO ;
   mat[1][1] =    sy   ;
   mat[1][2] = SGD_ZERO ;
   mat[1][3] = SGD_ZERO ;

   mat[2][0] = SGD_ZERO ;
   mat[2][1] = SGD_ZERO ;
   mat[2][2] = SGD_ONE  ;
   mat[2][3] = SGD_ZERO ;

   mat[3][0] =    tx   ;
   mat[3][1] =    ty   ;
   mat[3][2] = SGD_ZERO ;
   mat[3][3] = SGD_ONE  ;
}                                                                               



void sgdMakeLookAtMat4 ( sgdMat4 dst, const sgdVec3 eye,
                                    const sgdVec3 center,
                                    const sgdVec3 up )
{
  // Caveats:
  // 1) In order to compute the line of sight, the eye point must not be equal
  //    to the center point.
  // 2) The up vector must not be parallel to the line of sight from the eye
  //    to the center point.

  /* Compute the direction vectors */
  sgdVec3 x,y,z;

  /* Y vector = center - eye */
  sgdSubVec3 ( y, center, eye ) ;

  /* Z vector = up */
  sgdCopyVec3 ( z, up ) ;

  /* X vector = Y cross Z */
  sgdVectorProductVec3 ( x, y, z ) ;

  /* Recompute Z = X cross Y */
  sgdVectorProductVec3 ( z, x, y ) ;

  /* Normalize everything */
  sgdNormaliseVec3 ( x ) ;
  sgdNormaliseVec3 ( y ) ;
  sgdNormaliseVec3 ( z ) ;

  /* Build the matrix */
  sgdSetVec4 ( dst[0], x[0], x[1], x[2], SGD_ZERO ) ;
  sgdSetVec4 ( dst[1], y[0], y[1], y[2], SGD_ZERO ) ;
  sgdSetVec4 ( dst[2], z[0], z[1], z[2], SGD_ZERO ) ;
  sgdSetVec4 ( dst[3], eye[0], eye[1], eye[2], SGD_ONE ) ;
}

// -dw- inconsistent linkage!

sgdFloat sgdTriArea( sgdVec3 p0, sgdVec3 p1, sgdVec3 p2 )
{
  /* 
    From comp.graph.algorithms FAQ

	2A(P) = abs(N.(sum_{i=0}^{n-1}(v_i x v_{i+1})))
	This is an optimized version for a triangle
	but easily extended for planar polygon's with more sides
	by passing in the number of sides and the vv array
	sgdTriArea( int nsides, float **vv )
	and changing the normal calculation and the for loop appropriately
	sgdMakeNormal( norm, vv[0], vv[1], vv[2] )
	for( int i=0; i<n; i++ )
  */

  sgdVec3 sum;
  sgdZeroVec3( sum );

  sgdVec3 norm;
  sgdMakeNormal( norm, p0, p1, p2 );

  sgdFloat *vv[3];
  vv[0] = p0;
  vv[1] = p1;
  vv[2] = p2;

  for( int i=0; i<3; i++ )
  {
    int ii = (i+1) % 3;

    sum[0] += (vv[i][1] * vv[ii][2] - vv[i][2] * vv[ii][1]) ;
    sum[1] += (vv[i][2] * vv[ii][0] - vv[i][0] * vv[ii][2]) ;
    sum[2] += (vv[i][0] * vv[ii][1] - vv[i][1] * vv[ii][0]) ;
  }

  sgdFloat area = sgdAbs ( sgdScalarProductVec3 ( norm, sum ) ) ;

  return area / 2.0f ;
}

/***************************************************\
*   functions to get the angle between two vectors  *
\***************************************************/

SGDfloat sgdAngleBetweenVec3 ( sgdVec3 v1, sgdVec3 v2 )
{
  sgdVec3 nv1, nv2 ;

  sgdNormalizeVec3 ( nv1, v1 ) ;
  sgdNormalizeVec3 ( nv2, v2 ) ;
  return sgdAngleBetweenNormalizedVec3 ( nv1, nv2 ) ;
}

SGDfloat sgdAngleBetweenNormalizedVec3 (sgdVec3 first, sgdVec3 second, sgdVec3 normal)
{ 
 // result is in the range  0..360 degrees
 //
 // Attention: first and second have to be normalized
 // the normal is needed to decide between for example 0.123
 // looking "from one side" and -0.123 looking fomr the other

  SGDfloat myCos, abs1, abs2, SProduct, deltaAngle, myNorm;

  if((normal[0]==0) && (normal[1]==0) && (normal[2]==0))
  {	
    ulSetError ( UL_WARNING, "sgdGetAngleBetweenVectors: Normal is zero.");
    return 0.0 ;
  }

  sgdVec3 temp;

  sgdVectorProductVec3( temp, first, second);

  myNorm = sgdLengthVec3 ( temp );

  if ( (sgdScalarProductVec3(temp, normal))<0 )
    myNorm = -myNorm;

  if ( myNorm < -0.99999 )
    deltaAngle = -SGD_PI*0.5;
  else
  if ( myNorm > 0.99999 )
    deltaAngle = SGD_PI*0.5;
  else
    deltaAngle = (SGDfloat)asin((double)myNorm);

  // deltaAngle is in the range -SGD_PI*0.5 to +SGD_PI*0.5 here
  // However, the correct result could also be 
  // deltaAngleS := pi - deltaAngle
  // Please note that:
  // cos(deltaAngleS)=cos(pi-deltaAngle)=-cos(deltaAngle)
  // So, the question is whether + or - cos(deltaAngle)
  // is sgdScalarProductVec3(first, second)
  
  if ( deltaAngle < 0 ) 
    deltaAngle = deltaAngle + 2*SGD_PI; // unnessecary?
  
  SProduct = sgdScalarProductVec3(first, second);
  myCos = (SGDfloat) cos(deltaAngle);
  
  abs1 = SProduct - myCos;
  abs2 = SProduct + myCos;

  if ( abs1 < 0 ) abs1 = -abs1 ;
  if ( abs2 < 0 ) abs2 = -abs2 ;
  
  assert( (abs1 < 0.1) || (abs2 < 0.1) ) ;

  if ( abs2 < abs1 )
  {
    // deltaAngleS is the correct result

    if ( deltaAngle <= SGD_PI )
      deltaAngle = SGD_PI - deltaAngle ;
    else
      deltaAngle = 3*SGD_PI - deltaAngle ;
  }

  assert ( deltaAngle >= 0.0 ) ;
  assert ( deltaAngle <= 2.0*SGD_PI ) ;

  return deltaAngle * SGD_RADIANS_TO_DEGREES ;
}


SGDfloat sgdAngleBetweenVec3 ( sgdVec3 v1, sgdVec3 v2, sgdVec3 normal )
{
  // nornmal has to be normalized.
  sgdVec3 nv1, nv2 ;

  sgdNormalizeVec3 ( nv1, v1 ) ;
  sgdNormalizeVec3 ( nv2, v2 ) ;
  return sgdAngleBetweenNormalizedVec3 ( nv1, nv2, normal ) ;
}


/*********************\
*    sgdBox routines   *
\*********************/


void sgdBox::extend ( const sgdVec3 v )
{
  if ( isEmpty () )
  {
    sgdCopyVec3 ( min, v ) ;
    sgdCopyVec3 ( max, v ) ;
  }
  else
  {
    if ( v[0] < min[0] ) min[0] = v[0] ;
    if ( v[1] < min[1] ) min[1] = v[1] ;
    if ( v[2] < min[2] ) min[2] = v[2] ;
    if ( v[0] > max[0] ) max[0] = v[0] ;
    if ( v[1] > max[1] ) max[1] = v[1] ;
    if ( v[2] > max[2] ) max[2] = v[2] ;
  }
}


void sgdBox::extend ( const sgdBox *b )
{
  if ( b -> isEmpty () )
    return ;

  if ( isEmpty () )
  {
    sgdCopyVec3 ( min, b->getMin() ) ;
    sgdCopyVec3 ( max, b->getMax() ) ;
  }
  else
  {
    extend ( b->getMin() ) ;
    extend ( b->getMax() ) ;
  }
}


void sgdBox::extend ( const sgdSphere *s )
{
  if ( s -> isEmpty () ) 
    return ;

  /*
    In essence, this extends around a box around the sphere - which
    is still a perfect solution because both boxes are axially aligned.
  */

  sgdVec3 x ;

  sgdSetVec3 ( x, s->getCenter()[0]+s->getRadius(),
                 s->getCenter()[1]+s->getRadius(),
                 s->getCenter()[2]+s->getRadius() ) ;
  extend ( x ) ;

  sgdSetVec3 ( x, s->getCenter()[0]-s->getRadius(),
                 s->getCenter()[1]-s->getRadius(),
                 s->getCenter()[2]-s->getRadius() ) ;
  extend ( x ) ;
}


int sgdBox::intersects ( const sgdVec4 plane ) const
{
  /*
    Save multiplies by not redoing Ax+By+Cz+D for each point.
  */

  SGDfloat Ax_min        = plane[0] * min[0] ;
  SGDfloat By_min        = plane[1] * min[1] ;
  SGDfloat Cz_min_plus_D = plane[2] * min[2] + plane[3] ;

  SGDfloat Ax_max        = plane[0] * max[0] ;
  SGDfloat By_max        = plane[1] * max[1] ;
  SGDfloat Cz_max_plus_D = plane[2] * max[2] + plane[3] ;

  /*
    Count the number of vertices on the positive side of the plane.
  */

  int count = ( Ax_min + By_min + Cz_min_plus_D > SGD_ZERO ) +
              ( Ax_min + By_min + Cz_max_plus_D > SGD_ZERO ) +
              ( Ax_min + By_max + Cz_min_plus_D > SGD_ZERO ) +
              ( Ax_min + By_max + Cz_max_plus_D > SGD_ZERO ) +
              ( Ax_max + By_min + Cz_min_plus_D > SGD_ZERO ) +
              ( Ax_max + By_min + Cz_max_plus_D > SGD_ZERO ) +
              ( Ax_max + By_max + Cz_min_plus_D > SGD_ZERO ) +
              ( Ax_max + By_max + Cz_max_plus_D > SGD_ZERO ) ;

  /*
    The plane intersects the box unless all 8 are positive
    or none of them are positive.
  */
              
  return count != 0 && count != 8 ;
}



/**********************\
*  sgdSphere routines   *
\**********************/

void sgdSphere::extend ( const sgdVec3 v )
{
  if ( isEmpty () )
  {
    sgdCopyVec3 ( center, v ) ;
    radius = SGD_ZERO ;
    return ;
  }

  SGDfloat d = sgdDistanceVec3 ( center, v ) ;

  if ( d <= radius )  /* Point is already inside sphere */
    return ;

  SGDfloat new_radius = (radius + d) / SGD_TWO ;  /* Grow radius */

  SGDfloat ratio = (new_radius - radius) / d ;

  center[0] += (v[0]-center[0]) * ratio ;    /* Move center */
  center[1] += (v[1]-center[1]) * ratio ;
  center[2] += (v[2]-center[2]) * ratio ;

  radius = new_radius ;
}


void sgdSphere::extend ( const sgdBox *b )
{
  if ( b -> isEmpty () )
    return ;

  if ( isEmpty() )
  {
    sgdAddVec3   ( center, b->getMin(), b->getMax() ) ;
    sgdScaleVec3 ( center, SGD_HALF ) ;
    radius = sgdDistanceVec3 ( center, b->getMax() ) ;
    return ;
  }

  /*
    I can't think of a faster way to get an
    utterly minimal sphere.

    The tighter algorithm:- enclose each
    of eight vertices of the box in turn - it
    looks like being pretty costly.
    [8 sqrt()'s]

    The looser algorithm:- enclose the box
    with an empty sphere and then do a
    sphere-extend-sphere. This algorithm
    does well for close-to-cube boxes, but
    makes very poor spheres for long, thin
    boxes.
    [2 sqrt()'s]
  */

#ifdef DONT_REALLY_NEED_A_TIGHT_SPHERE_EXTEND_BOX

  /* LOOSER/FASTER sphere-around-sphere-around-box */
  sgdSphere s ;
  s.empty   ()    ;
  s.enclose ( b ) ;  /* Fast because s is empty */
    enclose ( s ) ;

#else

  /* TIGHTER/EXPENSIVE sphere-around-eight-points */
  sgdVec3 x ;
                                                        extend ( b->getMin() ) ;
  sgdSetVec3 ( x, b->getMin()[0],b->getMin()[1],b->getMax()[2] ) ; extend ( x ) ;
  sgdSetVec3 ( x, b->getMin()[0],b->getMax()[1],b->getMin()[2] ) ; extend ( x ) ;
  sgdSetVec3 ( x, b->getMin()[0],b->getMax()[1],b->getMax()[2] ) ; extend ( x ) ;
  sgdSetVec3 ( x, b->getMax()[0],b->getMin()[1],b->getMin()[2] ) ; extend ( x ) ;
  sgdSetVec3 ( x, b->getMax()[0],b->getMin()[1],b->getMax()[2] ) ; extend ( x ) ;
  sgdSetVec3 ( x, b->getMax()[0],b->getMax()[1],b->getMin()[2] ) ; extend ( x ) ;
                                                        extend ( b->getMax() ) ;
#endif
}


void sgdSphere::extend ( const sgdSphere *s )
{
  if ( s->isEmpty () )
    return ;

  if ( isEmpty () )
  {
    sgdCopyVec3 ( center, s->getCenter() ) ;
    radius = s->getRadius() ;
    return ;
  }

  /* 
    d == The distance between the sphere centers
  */

  SGDfloat d = sgdDistanceVec3 ( center, s->getCenter() ) ;

  if ( d + s->getRadius() <= radius )  /* New sphere is already inside this one */
    return ;

  if ( d + radius <= s->getRadius() )  /* New sphere completely contains this one */
  {
    sgdCopyVec3 ( center, s->getCenter() ) ;
    radius = s->getRadius() ;
    return ;
  } 

  /*
    Build a new sphere that completely contains the other two:

    The center point lies halfway along the line between
    the furthest points on the edges of the two spheres.
    Computing those two points is ugly - so we'll use similar
    triangles
  */

  SGDfloat new_radius = (radius + d + s->getRadius() ) / SGD_TWO ;

  SGDfloat ratio = ( new_radius - radius ) / d ;

  center[0] += ( s->getCenter()[0] - center[0] ) * ratio ;
  center[1] += ( s->getCenter()[1] - center[1] ) * ratio ;
  center[2] += ( s->getCenter()[2] - center[2] ) * ratio ;
  radius = new_radius ;
}


int sgdSphere::intersects ( const sgdBox *b ) const
{
  sgdVec3 closest ;

  if ( b->getMin()[0] > center[0] ) closest[0] = b->getMin()[0] ; else
  if ( b->getMax()[0] < center[0] ) closest[0] = b->getMax()[0] ; else
                                    closest[0] = center[0] ;

  if ( b->getMin()[1] > center[1] ) closest[1] = b->getMin()[1] ; else
  if ( b->getMax()[1] < center[1] ) closest[1] = b->getMax()[1] ; else
                                    closest[1] = center[1] ;

  if ( b->getMin()[2] > center[2] ) closest[2] = b->getMin()[2] ; else
  if ( b->getMax()[2] < center[2] ) closest[2] = b->getMax()[2] ; else
                                    closest[2] = center[2] ;

  return sgdCompare3DSqdDist ( closest, center, sgdSquare ( radius ) ) <= 0 ;
}


/************************\
*   sgdFrustum routines   *
\************************/

void sgdFrustum::update ()
{
  if ( fabs ( ffar - nnear ) < 0.1 )
  {
    ulSetError ( UL_WARNING, "sgdFrustum: Can't support depth of view <0.1 units.");
    return ;
  }

  if ( hfov != SGD_ZERO && vfov != SGD_ZERO )
  {
    if ( fabs ( hfov ) < 0.1 || fabs ( vfov ) < 0.1 )
    {
      ulSetError ( UL_WARNING, ortho ? 
		   "sgFrustum: Can't support width or height <0.1 units." : 
		   "sgFrustum: Can't support fields of view narrower than 0.1 degrees." ) ;
      return ;
    }

    if ( ortho )
    {
      right = SGD_HALF * hfov ;
      top   = SGD_HALF * vfov ;
    }
    else
    {
    right = nnear * tan ( hfov * SGD_DEGREES_TO_RADIANS / SGD_TWO ) ;
    top   = nnear * tan ( vfov * SGD_DEGREES_TO_RADIANS / SGD_TWO ) ;
    }
  
    left  = -right ;
    bot   = -top   ;
  }


  /* Compute the projection matrix */

  SGDfloat width  = right - left  ;
  SGDfloat height = top   - bot   ;
  SGDfloat depth  = ffar  - nnear ;

  if ( ortho )
  {
    /* orthographic */

    mat[0][0] =  SGD_TWO / width ;
    mat[0][1] =  SGD_ZERO ;
    mat[0][2] =  SGD_ZERO ;
    mat[0][3] =  SGD_ZERO ;

    mat[1][0] =  SGD_ZERO ;
    mat[1][1] =  SGD_TWO / height ;
    mat[1][2] =  SGD_ZERO ;
    mat[1][3] =  SGD_ZERO ;

    mat[2][0] =  SGD_ZERO ;
    mat[2][1] =  SGD_ZERO ;
    mat[2][2] = -SGD_TWO / depth ;
    mat[2][3] =  SGD_ZERO ;

    mat[3][0] = -( left  + right ) / width ;
    mat[3][1] = -( bot   + top   ) / height ;
    mat[3][2] = -( nnear + ffar  ) / depth ;
    mat[3][3] =  SGD_ONE ;
  }
  else
  {
    /* perspective */

    mat[0][0] =  SGD_TWO * nnear / width ;
    mat[0][1] =  SGD_ZERO ;
    mat[0][2] =  SGD_ZERO ;
    mat[0][3] =  SGD_ZERO ;

    mat[1][0] =  SGD_ZERO ;
    mat[1][1] =  SGD_TWO * nnear / height ;
    mat[1][2] =  SGD_ZERO ;
    mat[1][3] =  SGD_ZERO ;

    mat[2][0] =  ( right + left ) / width ;
    mat[2][1] =  ( top   + bot  ) / height ;
    mat[2][2] = -( ffar  + nnear ) / depth ;
    mat[2][3] = -SGD_ONE ;

    mat[3][0] =  SGD_ZERO ;
    mat[3][1] =  SGD_ZERO ;
    mat[3][2] = -SGD_TWO * nnear * ffar / depth ;
    mat[3][3] =  SGD_ZERO ;
  }


  /*
   * The clip planes are derived from the projection matrix.
   *
   * After projection (in clip coordinates), the clip planes are simply:
   * 
   *  left:    (  1,  0,  0,  1 )
   *  right:   ( -1,  0,  0,  1 )
   *  bottom:  (  0,  1,  0,  1 )
   *  top:     (  0, -1,  0,  1 )
   *  near:    (  0,  0,  1,  1 )
   *  far:     (  0,  0, -1,  1 )
   *
   * These can easily be transformed *backwards* by 
   * multiplying by the transposed projection matrix, i.e:
   *
   *  ( A )            ( A')
   *  ( B )  =  mat^T  ( B')
   *  ( C )            ( C')
   *  ( D )            ( D')
   *
   * where (A',B',C',D') represents a plane in clip coordinates,
   * and (A,B,C,D) is the same plane expressed in eye coordinates.
   */

  sgdSetVec4( plane[ SG_LEFT_PLANE  ],   SGD_ONE,  SGD_ZERO,  SGD_ZERO,  SGD_ONE );
  sgdSetVec4( plane[ SG_RIGHT_PLANE ],  -SGD_ONE,  SGD_ZERO,  SGD_ZERO,  SGD_ONE );
  sgdSetVec4( plane[ SG_BOT_PLANE   ],  SGD_ZERO,   SGD_ONE,  SGD_ZERO,  SGD_ONE );
  sgdSetVec4( plane[ SG_TOP_PLANE   ],  SGD_ZERO,  -SGD_ONE,  SGD_ZERO,  SGD_ONE );
  sgdSetVec4( plane[ SG_NEAR_PLANE  ],  SGD_ZERO,  SGD_ZERO,   SGD_ONE,  SGD_ONE );
  sgdSetVec4( plane[ SG_FAR_PLANE   ],  SGD_ZERO,  SGD_ZERO,  -SGD_ONE,  SGD_ONE );

  for ( int i = 0 ; i < 6 ; i++ )
  {
    sgdVec4 tmp ;

    for ( int j = 0 ; j < 4 ; j++ )
      tmp[j] = sgdScalarProductVec4 ( plane[i], mat[j] ) ;

    sgdScaleVec4 ( plane[i], tmp, SGD_ONE / sgdLengthVec3 ( tmp ) ) ;
  }
}



#define OC_LEFT_SHIFT   0
#define OC_RIGHT_SHIFT  1
#define OC_TOP_SHIFT    2
#define OC_BOT_SHIFT    3
#define OC_NEAR_SHIFT   4
#define OC_FAR_SHIFT    5

#define OC_ALL_ON_SCREEN 0x3F
#define OC_OFF_TRF      ((1<<OC_TOP_SHIFT)|(1<<OC_RIGHT_SHIFT)|(1<<OC_FAR_SHIFT))
#define OC_OFF_BLN      ((1<<OC_BOT_SHIFT)|(1<<OC_LEFT_SHIFT)|(1<<OC_NEAR_SHIFT))

int sgdFrustum::getOutcode ( const sgdVec3 pt ) const
{
  /* Transform the point by the Frustum's transform. */

  sgdVec4 tmp ;

  tmp [ 0 ] = pt [ 0 ] ;
  tmp [ 1 ] = pt [ 1 ] ;
  tmp [ 2 ] = pt [ 2 ] ;
  tmp [ 3 ] =  SGD_ONE  ;

  sgdXformPnt4 ( tmp, tmp, mat ) ;

  /*
    No need to divide by the 'w' component since we are only checking for
    results in the range 0..1
  */

  return (( tmp[0] <=  tmp[3] ) << OC_RIGHT_SHIFT ) |
         (( tmp[0] >= -tmp[3] ) << OC_LEFT_SHIFT  ) |
         (( tmp[1] <=  tmp[3] ) << OC_TOP_SHIFT   ) |
         (( tmp[1] >= -tmp[3] ) << OC_BOT_SHIFT   ) |
         (( tmp[2] <=  tmp[3] ) << OC_FAR_SHIFT   ) |
         (( tmp[2] >= -tmp[3] ) << OC_NEAR_SHIFT  ) ;
}

int sgdFrustum::contains ( const sgdVec3 pt ) const
{
  return getOutcode ( pt ) == OC_ALL_ON_SCREEN ;
}


int sgdFrustum::contains ( const sgdSphere *s ) const
{

  const SGDfloat *center = s->getCenter() ;
  const SGDfloat  radius = s->getRadius() ;

  /*
    Lop off half the database (roughly) with a quick near-plane test - and
    lop off a lot more with a quick far-plane test
  */

  if ( -center[2] + radius < nnear || -center[2] - radius > ffar )
    return SG_OUTSIDE ;  

  /*
    OK, so the sphere lies between near and far.

    Measure the distance of the center point from the four sides of the frustum,
    if it's outside by more than the radius then it's history.

    It's tempting to do a quick test to see if the center point is
    onscreen using sgdFrustumContainsPt - but that takes a matrix transform
    which is 16 multiplies and 12 adds - versus this test which does the
    whole task using only 12 multiplies and 8 adds.
  */

  /*
    A few operations are saved by observing that certain values in the plane 
    equations are zero or one. These are specific to orthographic and perspective 
    projections respectively.
  */

  SGDfloat sp1, sp2, sp3, sp4 ;
  
  if ( ortho )
  {
    /*
      left:    (  1,  0,  0,  x  )
      right:   ( -1,  0,  0,  x  )
      bottom:  (  0,  1,  0,  x  )
      top:     (  0, -1,  0,  x  )
    */
    sp1 = plane[ SG_LEFT_PLANE  ][3] + center[0] ;
    sp2 = plane[ SG_RIGHT_PLANE ][3] - center[0] ;
    sp3 = plane[ SG_BOT_PLANE   ][3] + center[1] ;
    sp4 = plane[ SG_TOP_PLANE   ][3] - center[1] ;
  }
  else
  {
    /*
      left:    (  x,  0,  x,  0  )
      right:   (  x,  0,  x,  0  )
      bottom:  (  0,  x,  x,  0  )
      top:     (  0,  x,  x,  0  )
    */
    sp1 = plane[ SG_LEFT_PLANE  ][0] * center[0] + plane[ SG_LEFT_PLANE  ][2] * center[2] ;
    sp2 = plane[ SG_RIGHT_PLANE ][0] * center[0] + plane[ SG_RIGHT_PLANE ][2] * center[2] ;
    sp3 = plane[ SG_BOT_PLANE   ][1] * center[1] + plane[ SG_BOT_PLANE   ][2] * center[2] ;
    sp4 = plane[ SG_TOP_PLANE   ][1] * center[1] + plane[ SG_TOP_PLANE   ][2] * center[2] ;
  }

  /* 
     Note: in the general case, we would have to do:

     sp1 = sgdScalarProductVec3 (  left_plane, center ) +  left_plane[3] ;
     sp2 = sgdScalarProductVec3 ( right_plane, center ) + right_plane[3] ;
     ...
     sp6 = sgdScalarProductVec3 (   far_plane, center ) +   far_plane[3] ;
  */
  

  if ( -sp1 > radius || -sp2 > radius || -sp3 > radius || -sp4 > radius )
    return SG_OUTSIDE ;
  
  /*
    If it's inside by more than the radius then it's *completely* inside
    and we can save time elsewhere if we know that for sure.
  */

  if ( sp1 >= radius && sp2 >= radius && sp3 >= radius && sp4 >= radius
       && -center[2] - radius >= nnear && -center[2] + radius <= ffar )
    return SG_INSIDE ;

  return SG_STRADDLE ;
}


int sgdFrustum::contains ( const sgdBox *b ) const 
{
  sgdVec3 p[8] = {
    { b->getMin()[0], b->getMin()[1], b->getMin()[2] },
    { b->getMax()[0], b->getMin()[1], b->getMin()[2] },
    { b->getMin()[0], b->getMax()[1], b->getMin()[2] },
    { b->getMax()[0], b->getMax()[1], b->getMin()[2] },
    { b->getMin()[0], b->getMin()[1], b->getMax()[2] },
    { b->getMax()[0], b->getMin()[1], b->getMax()[2] },
    { b->getMin()[0], b->getMax()[1], b->getMax()[2] },
    { b->getMax()[0], b->getMax()[1], b->getMax()[2] },
  } ;
  
  int all = -1 ;
  int one =  0 ;
  
  for (int i = 0 ; i < 8 ; i++ ) 
  {
    int tmp = ~ getOutcode ( p[i] ) ;
    all &= tmp ;
    one |= tmp ;
  }
  
  return ( all ? SG_OUTSIDE : one ? SG_STRADDLE : SG_INSIDE ) ;
}




SGDfloat sgdDistSquaredToLineVec3 ( const sgdLine3 line, const sgdVec3 pnt )
{
  sgdVec3 r ; sgdSubVec3 ( r, pnt, line.point_on_line ) ;
 
  return sgdScalarProductVec3 ( r, r ) -
         sgdScalarProductVec3 ( r, line.direction_vector ) ;
}



SGDfloat sgdDistSquaredToLineSegmentVec3 ( const sgdLineSegment3 line,
                                         const sgdVec3 pnt )
{
  sgdVec3 v ; sgdSubVec3 ( v, line.b, line.a ) ;
  sgdVec3 r1 ; sgdSubVec3 ( r1, pnt, line.a ) ;

  SGDfloat r1_dot_v = sgdScalarProductVec3 ( r1, v ) ;

  if ( r1_dot_v <= 0 )  /* Off the "A" end  */
    return sgdScalarProductVec3 ( r1, r1 ) ;

  sgdVec3 r2 ; sgdSubVec3 ( r2, pnt, line.b ) ;

  SGDfloat r2_dot_v = sgdScalarProductVec3 ( r2, v ) ;

  if ( r2_dot_v >= 0 )  /* Off the "B" end */
    return sgdScalarProductVec3 ( r2, r2 ) ;

  /* Closest point on line is on the line segment */
 
  return sgdScalarProductVec3 ( r1, r1 ) - r1_dot_v * r1_dot_v / sgdScalarProductVec3 ( v, v ) ;
}



void sgdMakeCoordMat4 ( sgdMat4 m, const SGDfloat x, const SGDfloat y, const SGDfloat z, const SGDfloat h, const SGDfloat p, const SGDfloat r )
{
  SGDfloat ch, sh, cp, sp, cr, sr, srsp, crsp, srcp ;

  if ( h == SGD_ZERO )
  {
    ch = SGD_ONE ;
    sh = SGD_ZERO ;
  }
  else
  {
    sh = sgdSin( h ) ;
    ch = sgdCos( h ) ;
  }

  if ( p == SGD_ZERO )
  {
    cp = SGD_ONE ;
    sp = SGD_ZERO ;
  }
  else
  {
    sp = sgdSin( p ) ;
    cp = sgdCos( p ) ;
  }

  if ( r == SGD_ZERO )
  {
    cr   = SGD_ONE ;
    sr   = SGD_ZERO ;
    srsp = SGD_ZERO ;
    srcp = SGD_ZERO ;
    crsp = sp ;
  }
  else
  {
    sr   = sgdSin( r ) ;
    cr   = sgdCos( r ) ;
    srsp = sr * sp ;
    crsp = cr * sp ;
    srcp = sr * cp ;
  }

  m[0][0] =  ch * cr - sh * srsp ;
  m[1][0] = -sh * cp ;
  m[2][0] =  sr * ch + sh * crsp ;
  m[3][0] =  x ;

  m[0][1] =  cr * sh + srsp * ch ;
  m[1][1] =  ch * cp ;
  m[2][1] =  sr * sh - crsp * ch ;
  m[3][1] =  y ;

  m[0][2] = -srcp ;
  m[1][2] =  sp ;
  m[2][2] =  cr * cp ;
  m[3][2] =  z ;

  m[0][3] =  SGD_ZERO ;
  m[1][3] =  SGD_ZERO ;
  m[2][3] =  SGD_ZERO ;
  m[3][3] =  SGD_ONE ;
}


void sgdMakeTransMat4 ( sgdMat4 m, const sgdVec3 xyz )
{
  m[0][1] = m[0][2] = m[0][3] =
  m[1][0] = m[1][2] = m[1][3] =
  m[2][0] = m[2][1] = m[2][3] = SGD_ZERO ;
  m[0][0] = m[1][1] = m[2][2] = m[3][3] = SGD_ONE ;
  sgdCopyVec3 ( m[3], xyz ) ;
}


void sgdMakeTransMat4 ( sgdMat4 m, const SGDfloat x, const SGDfloat y, const SGDfloat z )
{
  m[0][1] = m[0][2] = m[0][3] =
  m[1][0] = m[1][2] = m[1][3] =
  m[2][0] = m[2][1] = m[2][3] = SGD_ZERO ;
  m[0][0] = m[1][1] = m[2][2] = m[3][3] = SGD_ONE ;
  sgdSetVec3 ( m[3], x, y, z ) ;
}


void sgdSetCoord ( sgdCoord *dst, const sgdMat4 src )
{
  sgdCopyVec3 ( dst->xyz, src[3] ) ;
    
  sgdMat4 mat ;

  SGDfloat s = sgdLengthVec3 ( src[0] ) ;

  if ( s <= 0.00001 )
  {
    ulSetError ( UL_WARNING, "sgdMat4ToCoord: ERROR - Bad Matrix." ) ;
    sgdSetVec3 ( dst -> hpr, SGD_ZERO, SGD_ZERO, SGD_ZERO ) ;
    return ;
  }

  sgdScaleMat4 ( mat, src, SGD_ONE / s ) ;
    
  dst->hpr[1] = sgdASin ( _sgdClampToUnity ( mat[1][2] ) ) ;

  SGDfloat cp = sgdCos ( dst->hpr[1] ) ;
    
  /* If pointing nearly vertically up - then heading is ill-defined */

  if ( cp > -0.00001 && cp < 0.00001 )
  {
    SGDfloat cr = _sgdClampToUnity ( mat[0][1] ) ; 
    SGDfloat sr = _sgdClampToUnity (-mat[2][1] ) ;

    dst->hpr[0] = SGD_ZERO ;
    dst->hpr[2] = sgdATan2 ( sr, cr ) ;
  }
  else
  {
    cp = SGD_ONE / cp ;
    SGDfloat sr = _sgdClampToUnity ( -mat[0][2] * cp ) ;
    SGDfloat cr = _sgdClampToUnity (  mat[2][2] * cp ) ;
    SGDfloat sh = _sgdClampToUnity ( -mat[1][0] * cp ) ;
    SGDfloat ch = _sgdClampToUnity (  mat[1][1] * cp ) ;
  
    if ( (sh == SGD_ZERO && ch == SGD_ZERO) || (sr == SGD_ZERO && cr == SGD_ZERO) )
    {
      cr = _sgdClampToUnity ( mat[0][1] ) ;
      sr = _sgdClampToUnity (-mat[2][1] ) ;

      dst->hpr[0] = SGD_ZERO ;
    }
    else
      dst->hpr[0] = sgdATan2 ( sh, ch ) ;

    dst->hpr[2] = sgdATan2 ( sr, cr ) ;
  }
}


void sgdMakeNormal(sgdVec2 dst, const sgdVec2 a, const sgdVec2 b )
{
  sgdSubVec2 ( dst, b, a ) ;
  SGDfloat tmp = dst [ 0 ] ; dst [ 0 ] = -dst [ 1 ] ; dst [ 1 ] = tmp ;
  sgdNormaliseVec2 ( dst ) ;
}


void sgdMakeNormal(sgdVec3 dst, const sgdVec3 a, const sgdVec3 b, const sgdVec3 c )
{
  sgdVec3 ab ; sgdSubVec3 ( ab, b, a ) ;
  sgdVec3 ac ; sgdSubVec3 ( ac, c, a ) ;
  sgdVectorProductVec3 ( dst, ab,ac ) ; sgdNormaliseVec3 ( dst ) ;
}


void sgdPreMultMat4( sgdMat4 dst, const sgdMat4 src )
{
  sgdMat4 mat ;
  sgdMultMat4 ( mat, dst, src ) ;
  sgdCopyMat4 ( dst, mat ) ;
}

void sgdPostMultMat4( sgdMat4 dst, const sgdMat4 src )
{
  sgdMat4 mat ;
  sgdMultMat4 ( mat, src, dst ) ;
  sgdCopyMat4 ( dst, mat ) ;
}

void sgdMultMat4( sgdMat4 dst, const sgdMat4 m1, const sgdMat4 m2 )
{
  for ( int j = 0 ; j < 4 ; j++ )
  {
    dst[0][j] = m2[0][0] * m1[0][j] +
		m2[0][1] * m1[1][j] +
		m2[0][2] * m1[2][j] +
		m2[0][3] * m1[3][j] ;

    dst[1][j] = m2[1][0] * m1[0][j] +
		m2[1][1] * m1[1][j] +
		m2[1][2] * m1[2][j] +
		m2[1][3] * m1[3][j] ;

    dst[2][j] = m2[2][0] * m1[0][j] +
		m2[2][1] * m1[1][j] +
		m2[2][2] * m1[2][j] +
		m2[2][3] * m1[3][j] ;

    dst[3][j] = m2[3][0] * m1[0][j] +
		m2[3][1] * m1[1][j] +
		m2[3][2] * m1[2][j] +
		m2[3][3] * m1[3][j] ;
  }
}


void sgdTransposeNegateMat4 ( sgdMat4 dst, const sgdMat4 src )
{
  /* Poor man's invert - can be used when matrix is a simple rotate-translate */

  dst[0][0] = src[0][0] ;
  dst[1][0] = src[0][1] ;
  dst[2][0] = src[0][2] ;
  dst[3][0] = - sgdScalarProductVec3 ( src[3], src[0] ) ;

  dst[0][1] = src[1][0] ;
  dst[1][1] = src[1][1] ;
  dst[2][1] = src[1][2] ;
  dst[3][1] = - sgdScalarProductVec3 ( src[3], src[1] ) ;
                                                                               
  dst[0][2] = src[2][0] ;                                                      
  dst[1][2] = src[2][1] ;                                                      
  dst[2][2] = src[2][2] ;                                                      
  dst[3][2] = - sgdScalarProductVec3 ( src[3], src[2] ) ;
                                                                               
  dst[0][3] = SGD_ZERO ;
  dst[1][3] = SGD_ZERO ;                                                        
  dst[2][3] = SGD_ZERO ;                                                        
  dst[3][3] = SGD_ONE  ;
}


void sgdTransposeNegateMat4 ( sgdMat4 dst )
{
  sgdMat4 src ;
  sgdCopyMat4 ( src, dst ) ;
  sgdTransposeNegateMat4 ( dst, src ) ;
}



void sgdInvertMat4 ( sgdMat4 dst, const sgdMat4 src )
{
  sgdMat4 tmp ;

  sgdCopyMat4 ( tmp, src ) ;
  sgdMakeIdentMat4 ( dst ) ;

  for ( int i = 0 ; i != 4 ; i++ )
  {
    SGDfloat val = tmp[i][i] ;
    int ind = i ;
    int j ;

    for ( j = i + 1 ; j != 4 ; j++ )
    {
      if ( fabs ( tmp[i][j] ) > fabs(val) )
      {
        ind = j;
        val = tmp[i][j] ;
      }
    }

    if ( ind != i )
    {                   /* swap columns */
      for ( j = 0 ; j != 4 ; j++ )
      {
        SGDfloat t ;
        t = dst[j][i]; dst[j][i] = dst[j][ind]; dst[j][ind] = t ;
        t = tmp[j][i]; tmp[j][i] = tmp[j][ind]; tmp[j][ind] = t ;
      }
    }

    // if ( val == SG_ZERO)
    if ( fabs(val) <= DBL_EPSILON )
    {
      ulSetError ( UL_WARNING, "sg: ERROR - Singular matrix, no inverse!" ) ;
      sgdMakeIdentMat4 ( dst ) ;  /* Do *something* */
      return;
    }

    SGDfloat ival = SGD_ONE / val ;

    for ( j = 0 ; j != 4 ; j++ )
    {
      tmp[j][i] *= ival ;
      dst[j][i] *= ival ;
    }

    for (j = 0; j != 4; j++)
    {
      if ( j == i )
        continue ;

      val = tmp[i][j] ;

      for ( int k = 0 ; k != 4 ; k++ )
      {
        tmp[k][j] -= tmp[k][i] * val ;
        dst[k][j] -= dst[k][i] * val ;
      }
    }
  }
}



void sgdXformVec3 ( sgdVec3 dst, const sgdVec3 src, const sgdMat4 mat )
{
  SGDfloat t0 = src[ 0 ] ;
  SGDfloat t1 = src[ 1 ] ;
  SGDfloat t2 = src[ 2 ] ;

  dst[0] = ( t0 * mat[ 0 ][ 0 ] +
             t1 * mat[ 1 ][ 0 ] +
             t2 * mat[ 2 ][ 0 ] ) ;

  dst[1] = ( t0 * mat[ 0 ][ 1 ] +
             t1 * mat[ 1 ][ 1 ] +
             t2 * mat[ 2 ][ 1 ] ) ;

  dst[2] = ( t0 * mat[ 0 ][ 2 ] +
             t1 * mat[ 1 ][ 2 ] +
             t2 * mat[ 2 ][ 2 ] ) ;
}


void sgdXformPnt3 ( sgdVec3 dst, const sgdVec3 src, const sgdMat4 mat )
{
  SGDfloat t0 = src[ 0 ] ;
  SGDfloat t1 = src[ 1 ] ;
  SGDfloat t2 = src[ 2 ] ;

  dst[0] = ( t0 * mat[ 0 ][ 0 ] +
             t1 * mat[ 1 ][ 0 ] +
             t2 * mat[ 2 ][ 0 ] +
                  mat[ 3 ][ 0 ] ) ;

  dst[1] = ( t0 * mat[ 0 ][ 1 ] +
             t1 * mat[ 1 ][ 1 ] +
             t2 * mat[ 2 ][ 1 ] +
                  mat[ 3 ][ 1 ] ) ;

  dst[2] = ( t0 * mat[ 0 ][ 2 ] +
             t1 * mat[ 1 ][ 2 ] +
             t2 * mat[ 2 ][ 2 ] +
                  mat[ 3 ][ 2 ] ) ;
}


void sgdXformPnt4 ( sgdVec4 dst, const sgdVec4 src, const sgdMat4 mat )
{
  SGDfloat t0 = src[ 0 ] ;
  SGDfloat t1 = src[ 1 ] ;
  SGDfloat t2 = src[ 2 ] ;
  SGDfloat t3 = src[ 3 ] ;

  dst[0] = ( t0 * mat[ 0 ][ 0 ] +
             t1 * mat[ 1 ][ 0 ] +
             t2 * mat[ 2 ][ 0 ] +
             t3 * mat[ 3 ][ 0 ] ) ;

  dst[1] = ( t0 * mat[ 0 ][ 1 ] +
             t1 * mat[ 1 ][ 1 ] +
             t2 * mat[ 2 ][ 1 ] +
             t3 * mat[ 3 ][ 1 ] ) ;

  dst[2] = ( t0 * mat[ 0 ][ 2 ] +
             t1 * mat[ 1 ][ 2 ] +
             t2 * mat[ 2 ][ 2 ] +
             t3 * mat[ 3 ][ 2 ] ) ;

  dst[3] = ( t0 * mat[ 0 ][ 3 ] +
             t1 * mat[ 1 ][ 3 ] +
             t2 * mat[ 2 ][ 3 ] +
             t3 * mat[ 3 ][ 3 ] ) ;
}


void sgdFullXformPnt3 ( sgdVec3 dst, const sgdVec3 src, const sgdMat4 mat )
{
  sgdVec4 tmp ;

  tmp [ 0 ] = src [ 0 ] ;
  tmp [ 1 ] = src [ 1 ] ;
  tmp [ 2 ] = src [ 2 ] ;
  tmp [ 3 ] =   SGD_ONE  ;

  sgdXformPnt4 ( tmp, tmp, mat ) ;
  sgdScaleVec3 ( dst, tmp, SGD_ONE / tmp [ 3 ] ) ;
}

void sgdHPRfromVec3 ( sgdVec3 hpr, sgdVec3 src )
{
  sgdVec3 tmp ;
  sgdCopyVec3 ( tmp, src ) ;
  sgdNormaliseVec3 ( tmp ) ;
  hpr[0] = - sgdATan2 ( tmp [ 0 ], tmp [ 1 ] ) ;
  hpr[1] = - sgdATan2 ( tmp [ 2 ], sgdSqrt ( sgdSquare ( tmp [ 0 ] ) +
                                           sgdSquare ( tmp [ 1 ] ) ) ) ;
  hpr[2] = SGD_ZERO ;
}



/*
  Quaternion routines are Copyright (C) 1999
  Kevin B. Thompson <kevinbthompson@yahoo.com>
  Modified by Sylvan W. Clebsch <sylvan@stanford.edu>
  Largely rewritten by "Negative0" <negative0@earthlink.net>
*/


void sgdQuatToAngleAxis ( SGDfloat *angle,
                         SGDfloat *x, SGDfloat *y, SGDfloat *z,
                         const sgdQuat src )
{
  sgdVec3 axis ;

  sgdQuatToAngleAxis ( angle, axis, src ) ;

  *x = axis [ 0 ] ;
  *y = axis [ 1 ] ;
  *z = axis [ 2 ] ;
}


void sgdQuatToAngleAxis ( SGDfloat *angle, sgdVec3 axis, const sgdQuat src )
{
  SGDfloat a = (SGDfloat) acos ( src[SGD_W] ) ;
  SGDfloat s = (SGDfloat) sin  ( a ) ;

  *angle = a * SGD_RADIANS_TO_DEGREES * SGD_TWO ;

  if ( s == SGD_ZERO )
    sgdSetVec3 ( axis, SGD_ZERO, SGD_ZERO, SGD_ONE );
  else
  {
    sgdSetVec3   ( axis, src[SGD_X], src[SGD_Y], src[SGD_Z] ) ;
    sgdScaleVec3 ( axis, SGD_ONE / s ) ;
  }
}


void sgdAngleAxisToQuat ( sgdQuat dst,
                         const SGDfloat angle,
                         const SGDfloat x, const SGDfloat y, const SGDfloat z )
{
  sgdVec3 axis ; 
  sgdSetVec3 ( axis, x, y, z ) ;
  sgdAngleAxisToQuat ( dst, angle, axis ) ;
}


void sgdAngleAxisToQuat ( sgdQuat dst, const SGDfloat angle, const sgdVec3 axis )
{
  SGDfloat temp_angle = angle * SGD_DEGREES_TO_RADIANS / SGD_TWO ;

  sgdVec3 ax ;
  sgdNormaliseVec3 ( ax, axis ) ;

  SGDfloat s = - (SGDfloat) sin ( temp_angle ) ;

  dst[SGD_W] = (SGDfloat) cos ( temp_angle ) ;
  sgdScaleVec3 ( dst, ax, s ) ;
}


//from gamasutra.com
//by nb

void sgdMatrixToQuat( sgdQuat quat, const sgdMat4 m )
{
  SGDfloat tr, s, q[4] ;
  int   i, j, k ;

  int nxt[3] = {1, 2, 0};

  tr = m[0][0] + m[1][1] + m[2][2];

  // check the diagonal
  if (tr > SGD_ZERO )
  {
    s = (SGDfloat) sqrt (tr + SGD_ONE);
    quat[SGD_W] = s / SGD_TWO;
    s = SGD_HALF / s;
    quat[SGD_X] = (m[1][2] - m[2][1]) * s;
    quat[SGD_Y] = (m[2][0] - m[0][2]) * s;
    quat[SGD_Z] = (m[0][1] - m[1][0]) * s;
  }
  else
  {		
    // diagonal is negative
   	i = 0;
    if (m[1][1] > m[0][0]) i = 1;
    if (m[2][2] > m[i][i]) i = 2;
    j = nxt[i];
    k = nxt[j];
    s = sqrt ((m[i][i] - (m[j][j] + m[k][k])) + SGD_ONE);
    q[i] = s * SGD_HALF;
            
    if (s != SGD_ZERO) s = SGD_HALF / s;

    q[3] = (m[j][k] - m[k][j]) * s;
    q[j] = (m[i][j] + m[j][i]) * s;
    q[k] = (m[i][k] + m[k][i]) * s;

    quat[SGD_X] = q[0];
    quat[SGD_Y] = q[1];
    quat[SGD_Z] = q[2];
    quat[SGD_W] = q[3];
  }

  // seems to yield the inverse rotation, so:
  quat[SG_W] = - quat[SG_W];
}


void sgdMultQuat ( sgdQuat dst, const sgdQuat a, const sgdQuat b )
{
  /* [ ww' - v.v', vxv' + wv' + v'w ] */

  SGDfloat t[8];

  t[0] = (a[SGD_W] + a[SGD_X]) * (b[SGD_W] + b[SGD_X]);
  t[1] = (a[SGD_Z] - a[SGD_Y]) * (b[SGD_Y] - b[SGD_Z]);
  t[2] = (a[SGD_X] - a[SGD_W]) * (b[SGD_Y] + b[SGD_Z]);
  t[3] = (a[SGD_Y] + a[SGD_Z]) * (b[SGD_X] - b[SGD_W]);
  t[4] = (a[SGD_X] + a[SGD_Z]) * (b[SGD_X] + b[SGD_Y]);
  t[5] = (a[SGD_X] - a[SGD_Z]) * (b[SGD_X] - b[SGD_Y]);
  t[6] = (a[SGD_W] + a[SGD_Y]) * (b[SGD_W] - b[SGD_Z]);
  t[7] = (a[SGD_W] - a[SGD_Y]) * (b[SGD_W] + b[SGD_Z]);

  dst[SGD_W] =  t[1] + ((-t[4] - t[5] + t[6] + t[7]) * SGD_HALF);
  dst[SGD_X] =  t[0] - (( t[4] + t[5] + t[6] + t[7]) * SGD_HALF);
  dst[SGD_Y] = -t[2] + (( t[4] - t[5] + t[6] - t[7]) * SGD_HALF);
  dst[SGD_Z] = -t[3] + (( t[4] - t[5] - t[6] + t[7]) * SGD_HALF);
}

//from gamasutra.com
//by nb@netcom.ca 

void sgdMultQuat2 ( sgdQuat dst, const sgdQuat a, const sgdQuat b )
{
  SGDfloat A, B, C, D, E, F, G, H;

  A = (a[SGD_W] + a[SGD_X]) * (b[SGD_W] + b[SGD_X]) ;
  B = (a[SGD_Z] - a[SGD_Y]) * (b[SGD_Y] - b[SGD_Z]) ;
  C = (a[SGD_X] - a[SGD_W]) * (b[SGD_Y] + b[SGD_Z]) ;
  D = (a[SGD_Y] + a[SGD_Z]) * (b[SGD_X] - b[SGD_W]) ;
  E = (a[SGD_X] + a[SGD_Z]) * (b[SGD_X] + b[SGD_Y]) ;
  F = (a[SGD_X] - a[SGD_Z]) * (b[SGD_X] - b[SGD_Y]) ;
  G = (a[SGD_W] + a[SGD_Y]) * (b[SGD_W] - b[SGD_Z]) ;
  H = (a[SGD_W] - a[SGD_Y]) * (b[SGD_W] + b[SGD_Z]) ;


  dst[SGD_W] =  B + (-E - F + G + H) / SGD_TWO ;
  dst[SGD_X] =  A - ( E + F + G + H) / SGD_TWO ; 
  dst[SGD_Y] = -C + ( E - F + G - H) / SGD_TWO ;
  dst[SGD_Z] = -D + ( E - F - G + H) / SGD_TWO ;
}

//from gamasutra.com
//by nb@netcom.ca 

void sgdEulerToQuat(sgdQuat quat, const sgdVec3 hpr )
{
  SGDfloat cr, cp, cy, sr, sp, sy, cpcy, spsy;

// calculate trig identities
  cr = (SGDfloat) cos(hpr[2]*SGD_DEGREES_TO_RADIANS/SGD_TWO);
  cp = (SGDfloat) cos(hpr[1]*SGD_DEGREES_TO_RADIANS/SGD_TWO);
  cy = (SGDfloat) cos(hpr[0]*SGD_DEGREES_TO_RADIANS/SGD_TWO);

  sr = (SGDfloat) sin(hpr[2]*SGD_DEGREES_TO_RADIANS/SGD_TWO);
  sp = (SGDfloat) sin(hpr[1]*SGD_DEGREES_TO_RADIANS/SGD_TWO);
  sy = (SGDfloat) sin(hpr[0]*SGD_DEGREES_TO_RADIANS/SGD_TWO);
  
  cpcy = cp * cy;
  spsy = sp * sy;

  quat[SGD_W] = cr * cpcy + sr * spsy;
  quat[SGD_X] = sr * cpcy - cr * spsy;
  quat[SGD_Y] = cr * sp * cy + sr * cp * sy;
  quat[SGD_Z] = cr * cp * sy - sr * sp * cy;
}

//from darwin3d.com
// jeffl@darwin3d.com

void sgdQuatToEuler( sgdVec3 hpr, const sgdQuat quat )
{
  SGDfloat matrix[3][3];
  SGDfloat cx,sx;
  SGDfloat cy,sy;
  SGDfloat cz,sz;

  // CONVERT QUATERNION TO MATRIX - I DON'T REALLY NEED ALL OF IT

  matrix[0][0] = SGD_ONE - (SGD_TWO * quat[SGD_Y] * quat[SGD_Y])
                         - (SGD_TWO * quat[SGD_Z] * quat[SGD_Z]);
//matrix[0][1] = (SGD_TWO * quat->x * quat->y) - (SGD_TWO * quat->w * quat->z);
//matrix[0][2] = (SGD_TWO * quat->x * quat->z) + (SGD_TWO * quat->w * quat->y);

  matrix[1][0] = (SGD_TWO * quat[SGD_X] * quat[SGD_Y]) +
                          (SGD_TWO * quat[SGD_W] * quat[SGD_Z]);
//matrix[1][1] = SGD_ONE - (SGD_TWO * quat->x * quat->x)
//                      - (SGD_TWO * quat->z * quat->z);
//matrix[1][2] = (SGD_TWO * quat->y * quat->z) - (SGD_TWO * quat->w * quat->x);

  matrix[2][0] = (SGD_TWO * quat[SGD_X] * quat[SGD_Z]) -
                 (SGD_TWO * quat[SGD_W] * quat[SGD_Y]);
  matrix[2][1] = (SGD_TWO * quat[SGD_Y] * quat[SGD_Z]) +
                 (SGD_TWO * quat[SGD_W] * quat[SGD_X]);
  matrix[2][2] = SGD_ONE - (SGD_TWO * quat[SGD_X] * quat[SGD_X])
                        - (SGD_TWO * quat[SGD_Y] * quat[SGD_Y]);

  sy = -matrix[2][0];
  cy = sgdSqrt(SGD_ONE - (sy * sy));
  
  hpr[1] = sgdATan2( sy, cy );

  // AVOID DIVIDE BY ZERO ERROR ONLY WHERE Y= +-90 or +-270 
  // NOT CHECKING cy BECAUSE OF PRECISION ERRORS
  if (sy != SGD_ONE && sy != -SGD_ONE)	
  {
    cx = matrix[2][2] / cy;
    sx = matrix[2][1] / cy;
    hpr[0] = sgdATan2 ( sx, cx );

    cz = matrix[0][0] / cy;
    sz = matrix[1][0] / cy;
    hpr[2] = sgdATan2 ( sz, cz );
  }
  else
  {
    // SINCE Cos(Y) IS 0, I AM SCREWED.  ADOPT THE STANDARD Z = 0
    // I THINK THERE IS A WAY TO FIX THIS BUT I AM NOT SURE.  EULERS SUCK
    // NEED SOME MORE OF THE MATRIX TERMS NOW

    matrix[1][1] = SGD_ONE - (SGD_TWO * quat[SGD_X] * quat[SGD_X])
                          - (SGD_TWO * quat[SGD_Z] * quat[SGD_Z]);
    matrix[1][2] = (SGD_TWO * quat[SGD_Y] * quat[SGD_Z]) -
                   (SGD_TWO * quat[SGD_W] * quat[SGD_X]);

    cx =  matrix[1][1];
    sx = -matrix[1][2];
    hpr[0] = sgdATan2 ( sx, cx );

    cz = SGD_ONE ;
    sz = SGD_ZERO ;
    hpr[2] = sgdATan2 ( sz, cz );
  }
}


void sgdQuatToMatrix ( sgdMat4 dst, const sgdQuat q )
{
  SGDfloat two_xx = q[SGD_X] * (q[SGD_X] + q[SGD_X]) ;
  SGDfloat two_xy = q[SGD_X] * (q[SGD_Y] + q[SGD_Y]) ;
  SGDfloat two_xz = q[SGD_X] * (q[SGD_Z] + q[SGD_Z]) ;

  SGDfloat two_wx = q[SGD_W] * (q[SGD_X] + q[SGD_X]) ;
  SGDfloat two_wy = q[SGD_W] * (q[SGD_Y] + q[SGD_Y]) ;
  SGDfloat two_wz = q[SGD_W] * (q[SGD_Z] + q[SGD_Z]) ;

  SGDfloat two_yy = q[SGD_Y] * (q[SGD_Y] + q[SGD_Y]) ;
  SGDfloat two_yz = q[SGD_Y] * (q[SGD_Z] + q[SGD_Z]) ;

  SGDfloat two_zz = q[SGD_Z] * (q[SGD_Z] + q[SGD_Z]) ;

  sgdSetVec4 ( dst[0], SGD_ONE-(two_yy+two_zz), two_xy-two_wz, two_xz+two_wy, SGD_ZERO ) ;
  sgdSetVec4 ( dst[1], two_xy+two_wz, SGD_ONE-(two_xx+two_zz), two_yz-two_wx, SGD_ZERO ) ;
  sgdSetVec4 ( dst[2], two_xz-two_wy, two_yz+two_wx, SGD_ONE-(two_xx+two_yy), SGD_ZERO ) ;
  sgdSetVec4 ( dst[3], SGD_ZERO, SGD_ZERO, SGD_ZERO, SGD_ONE ) ;
}


//from gamasutra.com
//by nb@netcom.ca 

/************************************
 DEPRECATED - use sgdQuatToMatrix instead.
*************************************/

void sgdMakeRotMat42( sgdMat4 m, sgdQuat quat ){
  SGDfloat wx, wy, wz, xx, yy, yz, xy, xz, zz, x2, y2, z2;

  // calculate coefficients
  x2 = quat[SGD_X] + quat[SGD_X]; y2 = quat[SGD_Y] + quat[SGD_Y]; 
  z2 = quat[SGD_Z] + quat[SGD_Z];
  xx = quat[SGD_X] * x2;   xy = quat[SGD_X] * y2;   xz = quat[SGD_X] * z2;
  yy = quat[SGD_Y] * y2;   yz = quat[SGD_Y] * z2;   zz = quat[SGD_Z] * z2;
  wx = quat[SGD_W] * x2;   wy = quat[SGD_W] * y2;   wz = quat[SGD_W] * z2;

  m[0][0] = SGD_ONE- (yy + zz); 	m[0][1] = xy - wz;
  m[0][2] = xz + wy;		m[0][3] = SGD_ZERO ;
 
  m[1][0] = xy + wz;		m[1][1] = SGD_ONE- (xx + zz);
  m[1][2] = yz - wx;		m[1][3] = SGD_ZERO ;

  m[2][0] = xz - wy;		m[2][1] = yz + wx;
  m[2][2] = SGD_ONE- (xx + yy);		m[2][3] = SGD_ZERO ;

  m[3][0] = 0;			m[3][1] = 0;
  m[3][2] = 0;			m[3][3] = 1;
}

//from gamasutra.com
//by nb@netcom.ca 

void sgdSlerpQuat2( sgdQuat dst, const sgdQuat from, const sgdQuat to, const SGDfloat t )
{
  SGDfloat           to1[4];
  SGDfloat        omega, cosom, sinom, scale0, scale1;
  
  // calc cosine
  cosom = from[SGD_X] * to[SGD_X] +
          from[SGD_Y] * to[SGD_Y] +
          from[SGD_Z] * to[SGD_Z] +
          from[SGD_W] * to[SGD_W];

  // adjust signs (if necessary)

  if ( cosom < SG_ZERO )
  { 
    cosom = -cosom; 
    to1[0] = - to[SGD_X];
    to1[1] = - to[SGD_Y];
    to1[2] = - to[SGD_Z];
    to1[3] = - to[SGD_W];
  }
  else
  {
    to1[0] = to[SGD_X];
    to1[1] = to[SGD_Y];
    to1[2] = to[SGD_Z];
    to1[3] = to[SGD_W];
  }
  
  // calculate coefficients
#define DELTA SGD_ZERO 
  if ( (SGD_ONE- cosom) > DELTA ) {
    // standard case (slerp)
    omega = acos(cosom);
    sinom = sin(omega);
    scale0 = sin((SGD_ONE- t) * omega) / sinom;
    scale1 = sin(t * omega) / sinom;
    
  }
  else
  {
    // "from" and "to" quaternions are very close 
    //  ... so we can do a linear interpolation
    scale0 = SGD_ONE- t;
    scale1 = t;
  }
  
  // calculate final values
  dst[SGD_X] = scale0 * from[SGD_X] + scale1 * to1[0];
  dst[SGD_Y] = scale0 * from[SGD_Y] + scale1 * to1[1];
  dst[SGD_Z] = scale0 * from[SGD_Z] + scale1 * to1[2];
  dst[SGD_W] = scale0 * from[SGD_W] + scale1 * to1[3];
}

void sgdSlerpQuat( sgdQuat dst, const sgdQuat from, const sgdQuat to, const SGDfloat t )
{
  SGDfloat co, scale0, scale1;
  bool flip = false ;

  /* SWC - Interpolate between to quaternions */

  co = sgdScalarProductVec4 ( from, to ) ;

  if ( co < SGD_ZERO )
  {
    co = -co;
    flip = true ;
  }

  if ( co < SGD_ONE - (SGDfloat) 1e-6 )
  {
    SGDfloat o = (SGDfloat) acos ( co );
    SGDfloat so = SGD_ONE / (SGDfloat) sin ( o );

    scale0 = (SGDfloat) sin ( (SGD_ONE - t) * o ) * so;
    scale1 = (SGDfloat) sin ( t * o ) * so;
  }
  else
  {
    scale0 = SGD_ONE - t;
    scale1 = t;
  }

  if ( flip )
  {
    scale1 = -scale1 ;
  }

  dst[SGD_X] = scale0 * from[SGD_X] + scale1 * to[SGD_X] ;
  dst[SGD_Y] = scale0 * from[SGD_Y] + scale1 * to[SGD_Y] ;
  dst[SGD_Z] = scale0 * from[SGD_Z] + scale1 * to[SGD_Z] ;
  dst[SGD_W] = scale0 * from[SGD_W] + scale1 * to[SGD_W] ;
}




/* Function to rotate a vector through a given quaternion using the formula
 * R = Q r Q-1 -- this gives the components of a ROTATED vector in a STATIONARY
 * coordinate system.  We assume that Q is a unit quaternion.
 */
void sgdRotateVecQuat ( sgdVec3 vec, sgdQuat q )
{
  sgdVec3 rot ;
  sgdFloat qwqw = q[SG_W] * q[SG_W] ;
  sgdFloat qwqx = q[SG_W] * q[SG_X] ;
  sgdFloat qwqy = q[SG_W] * q[SG_Y] ;
  sgdFloat qwqz = q[SG_W] * q[SG_Z] ;
  sgdFloat qxqx = q[SG_X] * q[SG_X] ;
  sgdFloat qxqy = q[SG_X] * q[SG_Y] ;
  sgdFloat qxqz = q[SG_X] * q[SG_Z] ;
  sgdFloat qyqy = q[SG_Y] * q[SG_Y] ;
  sgdFloat qyqz = q[SG_Y] * q[SG_Z] ;
  sgdFloat qzqz = q[SG_Z] * q[SG_Z] ;
  rot[SG_X] = ( qwqw + qxqx - qyqy - qzqz ) * vec[SG_X] + 2.0f * ( qxqy - qwqz ) * vec[SG_Y] + 2.0f * ( qxqz + qwqy ) * vec[SG_Z] ;
  rot[SG_Y] = ( qwqw - qxqx + qyqy - qzqz ) * vec[SG_Y] + 2.0f * ( qyqz - qwqx ) * vec[SG_Z] + 2.0f * ( qxqy + qwqz ) * vec[SG_X] ;
  rot[SG_Z] = ( qwqw - qxqx - qyqy + qzqz ) * vec[SG_Z] + 2.0f * ( qxqz - qwqy ) * vec[SG_X] + 2.0f * ( qyqz + qwqx ) * vec[SG_Y] ;
  sgdCopyVec3 ( vec, rot ) ;
}

/* Function to rotate a vector through a given quaternion using the formula
 * R = Q-1 r Q -- this gives the components of a STATIONARY vector in a ROTATED
 * coordinate system.  We assume that Q is a unit quaternion.
 */
void sgdRotateCoordQuat ( sgdVec3 vec, sgdQuat q )
{
  sgdVec3 rot ;
  sgdFloat qwqw = q[SG_W] * q[SG_W] ;
  sgdFloat qwqx = q[SG_W] * q[SG_X] ;
  sgdFloat qwqy = q[SG_W] * q[SG_Y] ;
  sgdFloat qwqz = q[SG_W] * q[SG_Z] ;
  sgdFloat qxqx = q[SG_X] * q[SG_X] ;
  sgdFloat qxqy = q[SG_X] * q[SG_Y] ;
  sgdFloat qxqz = q[SG_X] * q[SG_Z] ;
  sgdFloat qyqy = q[SG_Y] * q[SG_Y] ;
  sgdFloat qyqz = q[SG_Y] * q[SG_Z] ;
  sgdFloat qzqz = q[SG_Z] * q[SG_Z] ;
  rot[SG_X] = ( qwqw + qxqx - qyqy - qzqz ) * vec[SG_X] + 2.0f * ( qxqy + qwqz ) * vec[SG_Y] + 2.0f * ( qxqz - qwqy ) * vec[SG_Z] ;
  rot[SG_Y] = ( qwqw - qxqx + qyqy - qzqz ) * vec[SG_Y] + 2.0f * ( qyqz + qwqx ) * vec[SG_Z] + 2.0f * ( qxqy - qwqz ) * vec[SG_X] ;
  rot[SG_Z] = ( qwqw - qxqx - qyqy + qzqz ) * vec[SG_Z] + 2.0f * ( qxqz + qwqy ) * vec[SG_X] + 2.0f * ( qyqz - qwqx ) * vec[SG_Y] ;
  sgdCopyVec3 ( vec, rot ) ;
}


sgdFloat sgdDistSquaredToLineLineSegment ( const sgdLineSegment3 seg, const sgdLine3 line )
{
  /* Convert the line segment to a line.  We will deal with the segment limits later. */
  sgdLine3 line2 ;
  sgdLineSegment3ToLine3 ( &line2, seg ) ;

  /* Get the dot product of the two direction vectors */
  sgdFloat t1_dot_t2 = sgdScalarProductVec3 ( line.direction_vector, line2.direction_vector ) ;

  /* If the lines are parallel, distance is just the distance from a point to the other line */
  if ( fabs ( t1_dot_t2 ) >= 1.0 )
    return sgdDistSquaredToLineVec3 ( line, seg.a ) ;

  /* Get the parametric coordinates of the closest points on the two lines.  The first line
   * is parameterized by r = r1 + t1 u while the second is parameterized by r = r2 + t2 v.
   * The square of the distance between them is [ ( r1 + t1 u ) - ( r2 + t2 v ) ] dot
   * [ ( r1 + t1 u ) - ( r2 + t2 v ) ].  Differentiating this dot product with respect to
   * u and v and setting the derivatives to zero gives a matrix equation:
   * [ 1         -(t1 dot t2) ] [ u ] = [  ( r1 - r2 ) dot t1 ]
   * [ -(t1 dot t2)         1 ] [ v ]   [ -( r1 - r2 ) dot t2 ]
   * We invert the matrix to get the equations below.
   */
  sgdVec3 r1_minus_r2 ;
  sgdSubVec3 ( r1_minus_r2, line.point_on_line, line2.point_on_line ) ;

  /* t1_t2_t2_minus_t1 = ( t1 dot t2 ) t2 - t1
   * t2_minus_t1_t2_t1 = t2 - ( t1 dot t2 ) t1
   */
  sgdVec3 t1_t2_t2_minus_t1, t2_minus_t1_t2_t1 ;
  sgdAddScaledVec3 ( t1_t2_t2_minus_t1, line.direction_vector, line2.direction_vector, -t1_dot_t2 ) ;
  sgdNegateVec3 ( t1_t2_t2_minus_t1 ) ;
  sgdAddScaledVec3 ( t2_minus_t1_t2_t1, line2.direction_vector, line.direction_vector, -t1_dot_t2 ) ;

  sgdFloat u = sgdScalarProductVec3 ( r1_minus_r2, t1_t2_t2_minus_t1 ) / ( 1.0f - t1_dot_t2 * t1_dot_t2 ) ;
  sgdFloat v = sgdScalarProductVec3 ( r1_minus_r2, t2_minus_t1_t2_t1 ) / ( 1.0f - t1_dot_t2 * t1_dot_t2 ) ;

  /* Since line 2 is a line segment, we limit "v" to between 0 and the distance between the points. */
  sgdFloat length = sgdDistanceVec3 ( seg.a, seg.b ) ;
  if ( v < 0.0 ) v = 0.0 ;
  if ( v > length ) v = length ;

  sgdVec3 point1, point2 ;
  sgdAddScaledVec3 ( point1, line.point_on_line, line.direction_vector, u ) ;
  sgdAddScaledVec3 ( point2, line2.point_on_line, line2.direction_vector, v ) ;
  return sgdDistanceSquaredVec3 ( point1, point2 ) ;
}


 
void sgdReflectInPlaneVec3 ( sgdVec3 dst, const sgdVec3 src, const sgdVec4 plane )
{
  SGDfloat src_dot_norm  = sgdScalarProductVec3 ( src, plane ) ;
 
  sgdVec3 tmp ;

  sgdScaleVec3 ( tmp, plane, SGD_TWO * src_dot_norm ) ;
  sgdSubVec3 ( dst, src, tmp ) ;
}



int sgdClassifyMat4 ( const sgdMat4 m )
{
  const SGDfloat epsilon = 1e-6 ;

  int flags = 0 ;


  SGDfloat sx, sy, sz ;

  if ( m[0][1] == SGD_ZERO && m[0][2] == SGD_ZERO &&
       m[1][0] == SGD_ZERO && m[1][2] == SGD_ZERO &&
       m[2][0] == SGD_ZERO && m[2][1] == SGD_ZERO )
  {

    int n = ( m[0][0] < 0 ) + ( m[1][1] < 0 ) + ( m[2][2] < 0 ) ;

    if ( n > 1 )
      flags |= SG_ROTATION ;

    if ( n % 2 != 0 )
      flags |= SG_MIRROR ;

    sx = m[0][0] * m[0][0] ;
    sy = m[1][1] * m[1][1] ;
    sz = m[2][2] * m[2][2] ;

  }
  else
  {

    flags |= SG_ROTATION ;

    if ( sgdAbs ( sgdScalarProductVec3 ( m[1], m[2] ) ) > epsilon ||
         sgdAbs ( sgdScalarProductVec3 ( m[2], m[0] ) ) > epsilon ||
         sgdAbs ( sgdScalarProductVec3 ( m[0], m[1] ) ) > epsilon )
    {
      flags |= SG_NONORTHO ;
    }

    sgdVec3 temp ;
    sgdVectorProductVec3 ( temp, m[0], m[1] ) ;
    SGDfloat det = sgdScalarProductVec3 ( temp, m[2] ) ;

    if ( det < 0 )
      flags |= SG_MIRROR ;

    sx = sgdScalarProductVec3 ( m[0], m[0] ) ;
    sy = sgdScalarProductVec3 ( m[1], m[1] ) ;
    sz = sgdScalarProductVec3 ( m[2], m[2] ) ;

  }


  if ( sgdAbs ( sx - sy ) > epsilon ||
       sgdAbs ( sx - sz ) > epsilon )
  {
    flags |= SG_NONORTHO ;
    flags |= SG_GENERAL_SCALE ; // also set general scale bit, though it may be deleted in the future
  }
  else
  {
    if ( sgdAbs ( sx - SGD_ONE ) > epsilon )
      flags |= SG_SCALE ;
  }


  if ( m[3][0] != SGD_ZERO || m[3][1] != SGD_ZERO || m[3][2] != SGD_ZERO )
  {
    flags |= SG_TRANSLATION ;
  }


  if ( m[0][3] != SGD_ZERO || m[1][3] != SGD_ZERO || m[2][3] != SGD_ZERO ||
       m[3][3] != SGD_ONE )
  {
    flags |= SG_PROJECTION ;
  }


  return flags ;
}

 
SGDfloat sgdTriangleSolver_ASAtoArea ( SGDfloat angA, SGDfloat lenB, SGDfloat angC )
{
  /* Get the third angle */

  SGDfloat angB = SGD_180 - (angA + angC) ;

  /* Use Sine Rule to get length of a second side - then use SAStoArea. */

  SGDfloat sinB = sgdSin ( angB ) ;

  if ( sinB == SGD_ZERO )
    return SGD_ZERO ;

  SGDfloat lenA = lenB * sgdSin(angA) / sinB ;

  return sgdTriangleSolver_SAStoArea ( lenA, angC, lenB ) ;
}


SGDfloat sgdTriangleSolver_SAStoArea ( SGDfloat lenA, SGDfloat angB, SGDfloat lenC )
{
  return SGD_HALF * lenC * lenA * sgdSin ( angB ) ;
}


SGDfloat sgdTriangleSolver_SSStoArea ( SGDfloat lenA, SGDfloat lenB, SGDfloat lenC )
{
  /* Heron's formula */

  SGDfloat s = ( lenA + lenB + lenC ) / SGD_TWO ;
  SGDfloat q = s * (s-lenA) * (s-lenB) * (s-lenC) ;

  /* Ikky illegal triangles generate zero areas. */

  return ( q <= SGD_ZERO ) ? SGD_ZERO : sgdSqrt ( q ) ;
}


SGDfloat sgdTriangleSolver_ASStoArea ( SGDfloat angB, SGDfloat lenA, SGDfloat lenB,
                                     int angA_is_obtuse )
{
  SGDfloat lenC ;

  sgdTriangleSolver_ASStoSAA ( angB, lenA, lenB, angA_is_obtuse,
                                                         &lenC, NULL, NULL ) ;

  return sgdTriangleSolver_SAStoArea ( lenA, angB, lenC ) ;
}

SGDfloat sgdTriangleSolver_SAAtoArea ( SGDfloat lenA, SGDfloat angB, SGDfloat angA )
{
  SGDfloat lenC ;

  sgdTriangleSolver_SAAtoASS ( lenA, angB, angA, NULL, NULL, &lenC ) ;

  return sgdTriangleSolver_SAStoArea ( lenA, angB, lenC ) ;
}

void sgdTriangleSolver_SSStoAAA ( SGDfloat  lenA, SGDfloat  lenB, SGDfloat  lenC,
                                 SGDfloat *angA, SGDfloat *angB, SGDfloat *angC )
{
  SGDfloat aa, bb, cc ;

  int flag =  ( lenA == SGD_ZERO )     |
             (( lenB == SGD_ZERO )<<1) |
             (( lenC == SGD_ZERO )<<2) ;

  /* Ikky zero-sized triangles generate zero/90 angles appropriately. */
  /* Ikky triangles with all lengths zero generate 60 degree angles. */
  /* Ikky impossible triangles generate all zero angles. */

  switch ( flag )
  {
    case 0 :  /* no zero-lengthed sides */
     /* Cosine law */
     aa = sgdACos (( lenB*lenB + lenC*lenC - lenA*lenA )/(SGD_TWO*lenB*lenC)) ;
     bb = sgdACos (( lenA*lenA + lenC*lenC - lenB*lenB )/(SGD_TWO*lenA*lenC)) ;
     cc = sgdACos (( lenA*lenA + lenB*lenB - lenC*lenC )/(SGD_TWO*lenA*lenB)) ;
     break ;

    case 1 :  /* lenA is zero */
     aa = SGD_ZERO ;
     bb = cc = SGD_90 ;
     break ;

    case 2 : /* lenB is zero */
     bb = SGD_ZERO ;
     aa = cc = SGD_90 ;
     break ;

    case 4 : /* lenC is zero */
      cc = SGD_ZERO ;
      aa = bb = SGD_90 ;
      break ;

    case 3 : /* Two lengths are zero and the third isn't?!? */
    case 5 :
    case 6 :
      aa = bb = cc = SGD_ZERO ;
      break ;

    default : /* All three sides are zero length */
      aa = bb = cc = SGD_60 ;
      break ;
  }

  if ( angA ) *angA = aa ;
  if ( angB ) *angB = bb ;
  if ( angC ) *angC = cc ;
}

void sgdTriangleSolver_SAStoASA ( SGDfloat  lenA, SGDfloat  angB, SGDfloat  lenC,
                                 SGDfloat *angA, SGDfloat *lenB, SGDfloat *angC )
{
  /* Get third side using Cosine Rule */

  SGDfloat s = lenC * lenC +
              lenA * lenA - SGD_TWO * lenC * lenA * sgdCos( angB ) ;

  SGDfloat lb = ( s <= SGD_ZERO ) ? SGD_ZERO : (SGDfloat) sqrt ( s ) ;

  if ( lenB ) *lenB = lb ;

  sgdTriangleSolver_SSStoAAA ( lenA, lb, lenC, angA, NULL, angC ) ;
}


void sgdTriangleSolver_ASAtoSAS ( SGDfloat  angA, SGDfloat  lenB, SGDfloat  angC,
                                 SGDfloat *lenA, SGDfloat *angB, SGDfloat *lenC )
{
  /* Find the missing angle */

  SGDfloat bb = SGD_180 - (angA + angC) ;

  if ( angB ) *angB = bb ;

  /* Use Sine Rule */

  SGDfloat sinB = sgdSin ( bb ) ;

  if ( sinB == SGD_ZERO )
  {
    if ( lenA ) *lenA = lenB / SGD_TWO ;  /* One valid interpretation */
    if ( lenC ) *lenC = lenB / SGD_TWO ;
  }
  else
  {
    if ( lenA ) *lenA = lenB * sgdSin(angA) / sinB ;
    if ( lenC ) *lenC = lenB * sgdSin(angC) / sinB ;
  }
}


void sgdTriangleSolver_ASStoSAA ( SGDfloat angB, SGDfloat lenA, SGDfloat lenB,
                                 int angA_is_obtuse,
                                 SGDfloat *lenC, SGDfloat *angA, SGDfloat *angC )
{
  /* Sine law */

  SGDfloat aa = (lenB == SGD_ZERO ) ? SGD_ZERO : sgdASin (lenA * sgdSin(angB)/lenB) ;

  if ( angA_is_obtuse )
    aa = SGD_180 - aa ;

  if ( angA ) *angA = aa ;

  /* Find the missing angle */

  SGDfloat cc = SGD_180 - (aa + angB) ;

  if ( angC ) *angC = cc ;

  /* Use SAStoASA to get the last length */

  sgdTriangleSolver_SAStoASA ( lenA, cc, lenB, NULL, lenC, NULL ) ;
}


void sgdTriangleSolver_SAAtoASS ( SGDfloat  lenA, SGDfloat  angB, SGDfloat  angA,
                                 SGDfloat *angC, SGDfloat *lenB, SGDfloat *lenC )
{
  /* Find the missing angle */

  SGDfloat cc = SGD_180 - (angB + angA) ;

  if ( angC ) *angC = cc ;

  sgdTriangleSolver_ASAtoSAS ( cc, lenA, angB, lenC, NULL, lenB ) ;
}



