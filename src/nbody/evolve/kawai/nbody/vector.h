/*
 *  vector.h: 3D vector operations include file
 *.............................................................................
 *    version 1:  Dec 1992   Piet Hut, Steve McMillan, Jun Makino
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  1) definition of class vector
 *.............................................................................
 */

// This is slightly modified version of vector header file
// taken from STARLAB
// -- J. Makino

#ifndef  STARLAB_VECTOR_H
#  define  STARLAB_VECTOR_H

//#include "stdinc.h"

#ifdef PERIODIC
inline real readjust_r(real x)
{
    return fmod(x+2.5, 1.0)-0.5;
}
#else
inline real readjust_r(real x)
{
    return x;
}
#endif

/*-----------------------------------------------------------------------------
 *  vector  --  a class for 3-dimensional vectors
 *-----------------------------------------------------------------------------
 */

const int ndim = 3;

class vector
{
private:
    
    real element[3];
    
public:
    
    //	Default: initialize to zero.
    
    vector(real c = 0)
    {element[0] = element[1] = element[2] = c;}
    
    vector(real x, real y, real z)
    {element[0] = x; element[1] = y; element[2] = z;}
    
    //  []: the return type is declared as a reference (&), so that it can be used
    //  on the left-hand side of an asignment, as well as on the right-hand side,
    //  i.e.  v[1] = 3.14  and  x = v[2]  are both allowed and work as expected.
    
    real & operator [] (int i)       {return element[i];}
    
    inline void print() {cout << element[0] << " " << element[1] << " "
			      << element[2] << "\n";}
#ifdef PERIODIC
    // PERIODIC basic reajustment menber function
    vector readjust(){
	return vector(readjust_r(element[0]),
		      readjust_r(element[1]),
		      readjust_r(element[2]));
    }
#else
    vector  readjust(){
	return vector(*this);
    }
#endif    
	
	
	
//	Unary -
	
        vector operator - ()
	    {return vector(-element[0], -element[1], -element[2]);}
	
	//	Dot product.
	
        real operator * (const vector& b)
	    {return element[0]*b.element[0]
		  + element[1]*b.element[1]
		  + element[2]*b.element[2];}

//	Outer product.

        vector operator ^ (const vector &b)
	    {return vector(element[1]*b.element[2] - element[2]*b.element[1],
			   element[2]*b.element[0] - element[0]*b.element[2],
			   element[0]*b.element[1] - element[1]*b.element[0]);}

//	Vector +, -

        vector operator + (const vector &b)
	    {return vector(element[0]+b.element[0],
			   element[1]+b.element[1],
			   element[2]+b.element[2]);}
        vector operator - (const vector &b)
	    {return vector(element[0]-b.element[0],
		 	   element[1]-b.element[1],
			   element[2]-b.element[2]);}

        friend vector operator + (real, const vector & );
        friend vector operator + (const vector &, real);

//	Scalar *, /

        friend vector operator * (real, const vector & );
        friend vector operator * (const vector &, real);
        friend vector operator / (const vector &, real);

//	Vector +=, -=, *=, /=

        vector& operator += (const vector& b)
	    {element[0] += b.element[0];       
	     element[1] += b.element[1];
	     element[2] += b.element[2];
	     return *this;}

	vector& operator -= (const vector& b)
	    {element[0] -= b.element[0];
	     element[1] -= b.element[1];
	     element[2] -= b.element[2];
	     return *this;}

	vector& operator *= (const real b)
	    {element[0] *= b;
	     element[1] *= b;
	     element[2] *= b;
	     return *this;}

	vector& operator /= (const real b)
	    {element[0] /= b;
	     element[1] /= b;
	     element[2] /= b;
	     return *this;}

//      Input / Output

        friend ostream & operator << (ostream & , const vector & );

	friend istream & operator >> (istream & , vector & );
};

inline  ostream & operator << (ostream & s, const vector & v)
	    {return s << v.element[0] << "  " << v.element[1]
		      << "  " << v.element[2];}

inline  istream & operator >> (istream & s, vector & v)
	    {s >> v.element[0] >> v.element[1] >> v.element[2];
	     return s;}

inline  real square(vector v) {return v*v;}
inline  real abs(vector v)    {return sqrt(v*v);}

inline  vector operator + (real b, const vector & v)
	    {return vector(b+v.element[0],
			   b+v.element[1],
			   b+v.element[2]);}

inline  vector operator + (const vector & v, real b)
	    {return vector(b+v.element[0],
			   b+v.element[1],
			   b+v.element[2]);}

inline  vector operator * (real b, const vector & v)
	    {return vector(b*v.element[0],
			   b*v.element[1],
			   b*v.element[2]);}

inline  vector operator * (const vector & v, real b)
	    {return vector(b*v.element[0],
			   b*v.element[1],
			   b*v.element[2]);}

inline  vector operator / (const vector & v, real b)
	    {return vector(v.element[0]/b,
			   v.element[1]/b,
			   v.element[2]/b);}

#endif

//=======================================================================//
//  +---------------+        _\|/_        +------------------------------\\
//  |  the end of:  |         /|\         |  inc/vector.h
//  +---------------+                     +------------------------------//
//========================= STARLAB =====================================\\
 
