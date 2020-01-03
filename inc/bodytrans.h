/*
 * BODYTRANS.H: include file for body-transformation functions.
 *		
 *	12-aug-92	added TWODIM and THREEDIM security	PJT
 *			relies on CFLAGS environment variable
 *	 1-mar-94       <math.h> loaded - don't declare them again
 *	12-apr-95	no more ARGS  - defer math stuff to stdinc.h
 *      31-dec-02       gcc3/SINGLEPREC
 *      24-sep-04       added macro defining r as specified in man page  WD
 */

#ifndef _bodytrans_h
#define _bodytrans_h

#include <stdinc.h>
#include <vectmath.h>
#include <snapshot/body.h>

/*
 * Better prototypes for the 'proc's (C++ needs is, gcc3 also for -DSINGLEPREC)
 */

typedef real (*rproc_body)(Body *, real, int);
typedef int  (*iproc_body)(Body *, real, int);

extern rproc_body btrtrans(string expr);
extern iproc_body btitrans(string expr);

#ifndef _bodytransc_h
/*
 * Macros for standard components of a body b.-- only needed in true bodytrans
 * routines, never never in nemo_main() applications, unless you promise not
 * to use the variable names that we macro-fied below here
 */

#define m     Mass(b)

#define pos   Pos(b)
#define x     Pos(b)[0]
#define y     Pos(b)[1]
#if defined(THREEDIM)
#define z     Pos(b)[2]
#define r     sqrt(x*x+y*y+z*z)
#else
#define r     sqrt(x*x+y*y)
#endif

#define vel   Vel(b)
#define vx    Vel(b)[0]
#define vy    Vel(b)[1]
#if defined(THREEDIM)
#define vz    Vel(b)[2]
#endif

#define phi   Phi(b)

//#define acc   Acc(b)
#define ax    Acc(b)[0]
#define ay    Acc(b)[1]
#if defined(THREEDIM)
#define az    Acc(b)[2]
#endif

#define aux   Aux(b)

#define key   Key(b)

#define dens  Dens(b)

#define eps   Eps(b)

#endif /* _bodytransc_h */
#endif /* _bodytrans_h  */
