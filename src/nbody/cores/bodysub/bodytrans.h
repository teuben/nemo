/*
 * BODYTRANS.H: include file for body-transformation functions.
 *		
 *	12-aug-92	added TWODIM and THREEDIM security	PJT
 *			relies on CFLAGS environment variable
 *	 1-mar-94       <math.h> loaded - don't declare them again
 */

#include <stdinc.h>
#include <vectmath.h>
#include <snapshot/body.h>

/*
 * External functions available to body-transformations.
 * Most are provided by -lm (math.h) with additional functions
 * that we keep in NEMO's kernel:
 *     dex(x) = pow(10.0, x)
 *     qbe(x) = x*x*x
 *     sqr(x) = x*x
 * supplied by bodytrans.c
 */

extern double dex ARGS((double));
extern double qbe ARGS((double));
extern double sqr ARGS((double));

/*
 * Macros for standard components of a body b.
 */

#define m     Mass(b)

#define pos   Pos(b)
#define x     Pos(b)[0]
#define y     Pos(b)[1]
#if defined(THREEDIM)
#define z     Pos(b)[2]
#endif

#define vel   Vel(b)
#define vx    Vel(b)[0]
#define vy    Vel(b)[1]
#if defined(THREEDIM)
#define vz    Vel(b)[2]
#endif

#define phi   Phi(b)

#define acc   Acc(b)
#define ax    Acc(b)[0]
#define ay    Acc(b)[1]
#if defined(THREEDIM)
#define az    Acc(b)[2]
#endif

#define aux   Aux(b)

#define key   Key(b)
