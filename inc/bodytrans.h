/*
 * BODYTRANS.H: include file for body-transformation functions.
 *		
 *	12-aug-92	added TWODIM and THREEDIM security	PJT
 *			relies on CFLAGS environment variable
 *	 1-mar-94       <math.h> loaded - don't declare them again
 *	12-apr-95	no more ARGS  - defer math stuff to stdinc.h
 */

#include <stdinc.h>
#include <vectmath.h>
#include <snapshot/body.h>

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
