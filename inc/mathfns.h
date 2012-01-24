/****************************************************************************/
/* MATHFNS.H: header file for system and zeno math functions; assumes role  */
/* of math.h.  Defines real-valued synonyms for system functions (eg, rsqrt */
/* for square root) and zeno functions (eg, seval), depending on precision  */
/* switch (MIXEDPREC, SINGLEPREC, or DOUBLEPREC).                           */
/* Copyright (c) 2001 by Joshua E. Barnes, Honolulu, Hawai`i.               */
/*                                                                          */
/* 22-jun-01     adapted for NEMO - included by stdinc.h                    */
/*               note, NEMO does not use MIXEDPREC yet                      */
/* 18-Sep-08  WD made sqr,qbe,dex inline                                    */
/* ??-???-??  WD adapted for MAC                                            */
/* 24-Jan-12  WD commented decl of rexp2() which clashes with math.h        */
/****************************************************************************/

#ifndef _mathfns_h
#define _mathfns_h

#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

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

#if defined(LINUX) || defined(linux) || defined (__CYGWIN__) || defined(darwin)

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

#else

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

inline static double sqr(double x) { return x*x; }
inline static double qbe(double x) { return x*x*x; }
/* inline static double dex(double x) { return pow(10.0,x); } */
inline static double dex(double x) { return exp(M_LN10*x); }
inline static float fsqr(float x) { return x*x; }
inline static float fqbe(float x) { return x*x*x; }
inline static float fdex(float x) { return rexp(M_LN10*x); }

/* real rsqr(real); */
/* real rqbe(real); */
/* real rdex(real); */

/* real rlog2(real);*/
/* real rexp2(real);*/

#if defined(SINGLEPREC)
float rcbrt(float);
#endif



/*
 * Random number functions available only in double precision.
 */

int init_xrandom(string);
int set_xrandom(int);
double xrandom(double, double);
double grandom(double, double);
double frandom(double, double, real_proc);

#ifdef NEMO

void pickshell(real *, int, double);
void pickball(real *, int, double);
void pickbox(real *, int, double);

#else
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

#endif /* !NEMO */

#ifdef __cplusplus
}
#endif

#endif  /* ! _mathfns_h */
