/****************************************************************************/
/* MATHFNS.H: header file for system and zeno math functions; assumes role  */
/* of math.h.  Defines real-valued synonyms for system functions (eg, rsqrt */
/* for square root) and zeno functions (eg, seval), depending on precision  */
/* switch (MIXEDPREC, SINGLEPREC, or DOUBLEPREC).                           */
/* Copyright (c) 2001 by Joshua E. Barnes, Honolulu, Hawai`i.               */
/****************************************************************************/

#ifndef _mathfns_h
#define _mathfns_h

#include <math.h>

/*
 * System math functions.  Use double-precision versions in mixed or
 * double precision, and single-precisions versions otherwise.
 */

#if defined(MIXEDPREC) || defined(DOUBLEPREC)
#define rsqrt    sqrt
#define rcbrt    cbrt
#define rsin     sin
#define rcos     cos
#define rtan     tan
#define rasin    asin
#define racos    acos
#define ratan    atan
#define ratan2   atan2
#define rlog     log
#define rexp     exp
#define rlog10   log10
#define rsinh    sinh
#define rcosh    cosh
#define rtanh    tanh
#define rpow     pow
#define rabs     fabs
#define rfloor   floor
#define rceil    ceil
#endif

#if defined(SINGLEPREC)

#if !defined(LINUX)

#define rsqrt    fsqrt
#define rsin     fsin
#define rcos     fcos
#define rtan     ftan
#define rasin    fasin
#define racos    facos
#define ratan    fatan
#define ratan2   fatan2
#define rlog     flog
#define rexp     fexp
#define rlog10   log10f
#define rsinh    fsinh
#define rcosh    fcosh
#define rtanh    ftanh
#define rpow     powf
#define rabs     fabsf
#define rfloor   ffloor
#define rceil    fceil

#else

#define rsqrt    sqrtf
#define rsin     sinf
#define rcos     cosf
#define rtan     tanf
#define rasin    asinf
#define racos    acosf
#define ratan    atanf
#define ratan2   atan2f
#define rlog     logf
#define rexp     expf
#define rlog10   log10f
#define rsinh    sinhf
#define rcosh    coshf
#define rtanh    tanhf
#define rpow     powf
#define rabs     fabsf
#define rfloor   floorf
#define rceil    ceilf

#endif

#endif

/*
 * Functions in mathfns.c; invoked just like those above.
 */

#if defined(MIXEDPREC) || defined(DOUBLEPREC)
#define rsqr     sqr
#define rqbe     qbe
#define rlog2    log2
#define rexp2    exp2
#define rdex     dex
#endif

#if defined(SINGLEPREC)
#define rsqr     fsqr
#define rqbe     fqbe
#define rlog2    flog2
#define rexp2    fexp2
#define rdex     fdex
#define rcbrt    fcbrt
#endif

real rsqr(real);
real rqbe(real);
real rlog2(real);
real rexp2(real);
real rdex(real);

#if defined(SINGLEPREC)
float rcbrt(float);
#endif

/*
 * Random number functions available only in double precision.
 */

double xrandom(double, double);
double grandom(double, double);

/*
 * Functions which traffic in pointers to real values must be provided
 * in all three variants.
 */

#if defined(MIXEDPREC)
#define pickshell mpickshell
#define pickball  mpickball
#define pickbox   mpickbox
#endif

#if defined(DOUBLEPREC)
#define pickshell dpickshell
#define pickball  dpickball
#define pickbox   dpickbox
#endif

#if defined(SINGLEPREC)
#define pickshell fpickshell
#define pickball  fpickball
#define pickbox   fpickbox
#endif

/* Functions in pickpnt.c */

void pickshell(real *, int, real);
void pickball(real *, int, real);
void pickbox(real *, int, real);

#endif  /* ! _mathfns_h */
