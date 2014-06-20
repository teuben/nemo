#ifndef CVECUTILS_H
#define CVECUTILS_H

namespace vectutils {
#define NDIM 3
template <class T> inline void clrv(T * v)			// CLeaR Vector
{
    register int _i;
    for (_i = 0; _i < NDIM; _i++)
        (v)[_i] = 0.0;
}

template <class T> inline void setvs(T * v,T s)		        // SET Vector to Scalar
{
    int _i;
    for (_i = 0; _i < NDIM; _i++)
        (v)[_i] = (s);
}

template <class T> inline void setv(T * v, T * u)		// SET Vector
{
    int _i;
    for (_i = 0; _i < NDIM; _i++)
        (v)[_i] = (u)[_i];
}

template <class T> inline void  mulvs(T * _vp, T * _up , T s)	// MULtiply Vector by Scalar
{
    *_vp++ = (*_up++) * (s);
    *_vp++ = (*_up++) * (s);
    *_vp   = (*_up  ) * (s);
}

template <class T> inline void  addvs(T * v, T * u, T s)	// ADD Vector and Scalar
{
    register int _i;
    for (_i = 0; _i < NDIM; _i++)
        (v)[_i] = (u)[_i] + (s);
}

template <class T> inline void  addmulvs(T * v,T * u,T s)       // MUL Vect by Scalar, ADD to vect
{
    (v)[0] += (u)[0] * (s);
    (v)[1] += (u)[1] * (s);
    (v)[2] += (u)[2] * (s);
}

template <class T> inline void  addv(T *_vp,T * _up,T *_wp)	// ADD Vector
{

    *_vp++ = (*_up++) + (*_wp++);
    *_vp++ = (*_up++) + (*_wp++);
    *_vp   = (*_up  ) + (*_wp  );
}
template <class T> inline void  subv(T * v, T * u, T * w)	// SUBtract Vector
{
    register int _i;
    for (_i = 0; _i < NDIM; _i++)
        (v)[_i] = (u)[_i] - (w)[_i];
}

template <class T> inline void  divvs(T * v,T * u, T s)		// DIVide Vector by Scalar
{
    register int _i;
    for (_i = 0; _i < NDIM; _i++)
        (v)[_i] = (u)[_i] / (s);
}

#undef dotvp
template <class T> inline void  dotvp(T &s, T * _vp, T * _up)	/* DOT Vector Product */
{
    s  = (*_vp++) * (*_up++);
    s += (*_vp++) * (*_up++);
    s += (*_vp  ) * (*_up  );
}

}

#endif // CVECUTILS_H
