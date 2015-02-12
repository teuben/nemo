// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2015                                  
// e-mail:   Jean-Charles.Lambert@lam.fr                                      
// address:  Centre de donneeS Astrophysique de Marseille (CeSAM)              
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#ifndef GLNEMOVEC3D_H
#define GLNEMOVEC3D_H
#include <math.h>
#include <iostream>
/**
	@author Jean-Charles Lambert <jean-charles.lambert@lam.fr>
*/
namespace glnemo {


class Vec3D{
public:

         ~Vec3D() {};
	// Data
	float x, y, z;

	// Ctors
	Vec3D( float InX, float InY, float InZ ) : x( InX ), y( InY ), z( InZ )
		{
		}
        Vec3D( const Vec3D& V) : x( V.x ), y( V.y ), z( V.z )
		{
		}
	Vec3D( ) : x(0), y(0), z(0)
		{
		}
        inline void set( const float InX, const float InY, const float InZ ) {
            x = InX; y = InY; z = InZ;
        }

	// Operator Overloads
	inline bool operator== (const Vec3D& V2) const
		{
		return (x == V2.x && y == V2.y && z == V2.z);
		}
        inline Vec3D& operator=(const Vec3D& V2)
		{
                  x=V2.x; y=V2.y; z=V2.z;
		  return *this;
		}
	inline Vec3D operator+ (const Vec3D& V2) const
		{
		return Vec3D( x + V2.x,  y + V2.y,  z + V2.z);
		}
	inline Vec3D operator- (const Vec3D& V2) const
		{
		return Vec3D( x - V2.x,  y - V2.y,  z - V2.z);
		}
	inline Vec3D operator- ( ) const
		{
		return Vec3D(-x, -y, -z);
		}

	inline Vec3D operator/ (float S ) const
		{
		float fInv = 1.0f / S;
		return Vec3D (x * fInv , y * fInv, z * fInv);
		}
	inline Vec3D operator/ (const Vec3D& V2) const
		{
		return Vec3D (x / V2.x,  y / V2.y,  z / V2.z);
		}
	inline Vec3D operator* (const Vec3D& V2) const
		{
		return Vec3D (x * V2.x,  y * V2.y,  z * V2.z);
		}
	inline Vec3D operator* (float S) const
		{
		return Vec3D (x * S,  y * S,  z * S);
		}

	inline void operator+= ( const Vec3D& V2 )
		{
		x += V2.x;
		y += V2.y;
		z += V2.z;
		}
	inline void operator-= ( const Vec3D& V2 )
		{
		x -= V2.x;
		y -= V2.y;
		z -= V2.z;
		}

	inline float operator[] ( int i )
		{
		if ( i == 0 ) return x;
		else if ( i == 1 ) return y;
		else return z;
		}

	// Functions
	inline float Dot( const Vec3D &V1 ) const
		{
		return V1.x*x + V1.y*y + V1.z*z;
		}

	inline Vec3D CrossProduct( const Vec3D &V2 ) const
		{
		return Vec3D(
			y * V2.z  -  z * V2.y,
			z * V2.x  -  x * V2.z,
			x * V2.y  -  y * V2.x 	);
		}

	// Return vector rotated by the 3x3 portion of matrix m
	// (provided because it's used by bbox.cpp in article 21)
	Vec3D RotByMatrix( const float m[16] ) const
	{
	return Vec3D(
		x*m[0] + y*m[4] + z*m[8],
		x*m[1] + y*m[5] + z*m[9],
		x*m[2] + y*m[6] + z*m[10] );
 	}

	// These require math.h for the sqrtf function
	inline float Magnitude( ) const
		{
		return sqrt( x*x + y*y + z*z );
		}

	inline float Distance( const Vec3D &V1 ) const
		{
		return ( *this - V1 ).Magnitude();	
		}

	inline void Normalize()
		{
		float fMag = ( x*x + y*y + z*z );
		if (fMag == 0) {return;}

		float fMult = 1.0f/sqrtf(fMag);            
		x *= fMult;
		y *= fMult;
		z *= fMult;
		return;
		}
        inline void Display(std::string s="") {
		std::cerr <<s<< x << " " << y << " " << z << "\n";
        }
};

}

#endif
