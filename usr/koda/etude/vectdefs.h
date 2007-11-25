/****************************************************************************/
/* VECTDEFS.H: definitions from vectmath.h, separated for programs which    */
/* need to define vectors without loading the whole mess of definitions.    */
/* Copyright (c) 1999 by Joshua E. Barnes, Tokyo, JAPAN.                    */
/*    revised for SPH calculation by Jin Koda, Tokyo, JAPAN. 2000           */
/****************************************************************************/

#ifndef _vectdefs_h
#define _vectdefs_h

#define TWODIM

#if !defined(NDIM) && !defined(TWODIM) && !defined(THREEDIM)
#define THREEDIM                              /* supply default dimensions  */
#endif

#if defined(THREEDIM) || (NDIM==3)
#undef  TWODIM
#define THREEDIM
#define NDIM 3
#endif

#if defined(TWODIM) || (NDIM==2)
#undef  THREEDIM
#define TWODIM
#define NDIM 2
#endif

typedef real vector[NDIM];
typedef real matrix[NDIM][NDIM];

#endif  /* ! _vectdefs_h */
