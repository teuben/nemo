/****************************************************************************/
/* MATHFNS.C: utility routines for various sorts of math operations. Most   */
/* these functions work with real values, meaning that they can handle      */
/* either floats or doubles, depending on compiler switches.                */
/* Copyright (c) 1999 by Joshua E. Barnes, Tokyo, JAPAN.                    */
/*    revised for SPH calculation by Jin Koda, Tokyo, JAPAN. 2000           */
/****************************************************************************/

#include "stdinc.h"
#include "mathfns.h"


/*
 * RSQR, RQBE: compute x*x and x*x*x.
 */

real rsqr(real x)
{
    return (x * x);
}

real rqbe(real x)
{
    return (x * x * x);
}

/*
 * RLOG2, REXP2: log, inverse log to base two.
 */

real rlog2(real x)
{
    return (rlog(x) / M_LN2);
}

real rexp2(real x)
{
    return (rexp(M_LN2 * x));
}

/*
 * RDEX: inverse log base ten.
 */

real rdex(real x)
{
    return (rexp(M_LN10 * x));
}

#if defined(SINGLEPREC)

/*
 * FCBRT: floating cube root.
 */

float fcbrt(float x)
{
    return ((float) cbrt((double) x));
}

#endif

/*
 * XRANDOM: floating-point random number routine.
 */

double xrandom(double xl, double xh)
{

    return (xl + (xh - xl) * random0());
}

/*
 * GRANDOM: normally distributed random number (polar method).
 * Reference: Knuth, vol. 2, p. 104.
 */

double grandom(double mean, double sdev)
{
    double v1, v2, s;

    do {
        v1 = xrandom(-1.0, 1.0);
        v2 = xrandom(-1.0, 1.0);
        s = v1*v1 + v2*v2;
    } while (s >= 1.0);
    return (mean + sdev * v1 * sqrt(-2.0 * log(s) / s));
}

/*
 * PICKSHELL: pick point on shell.
 */

void pickshell(real vec[], int ndim, real rad)
{
    real rsq, rscale;
    int i;

    do {
        rsq = 0.0;
        for (i = 0; i < ndim; i++) {
            vec[i] = xrandom(-1.0, 1.0);
            rsq = rsq + vec[i] * vec[i];
        }
    } while (rsq > 1.0);
    rscale = rad / rsqrt(rsq);
    for (i = 0; i < ndim; i++)
        vec[i] = vec[i] * rscale;
}

/*
 * PICKBALL: pick point within ball.
 */

void pickball(real vec[], int ndim, real rad)
{
    real rsq;
    int i;

    do {
        rsq = 0.0;
        for (i = 0; i < ndim; i++) {
            vec[i] = xrandom(-1.0, 1.0);
            rsq = rsq + vec[i] * vec[i];
        }
    } while (rsq > 1.0);
    for (i = 0; i < ndim; i++)
        vec[i] = vec[i] * rad;
}

/*
 * PICKBOX: pick point within box.
 */

void pickbox(real vec[], int ndim, real size)
{
    int i;

    for (i = 0; i < ndim; i++)
        vec[i] = xrandom(- size, size);
}

/*
 * SRANDOM0: Initialize portable random number generator. This functuion
 * and the following one are adapted from the RAN1 random number generator
 * (see Press, W. et al. 1989, Numerical Recipes, Cambridge Univ. Press).
 */

#define NTAB 32
static long rtab[NTAB];
static long seed;
static long last;

void srandom0(long seed0)
{
    long a=16807, m=2147483647, q=127773, r=2836, k;
    int i;

    k = seed0 / q;
    seed = a * (seed0 - k * q) - k * r;
    for (i=NTAB+7; i>=0;  i--) {
	k = seed / q;
	seed = a * (seed - k * q) - k * r;
	if (seed < 0) seed += m;
	if (i<NTAB) rtab[i] = seed;
    }
    last = seed;
}

/*
 * RANDOM0: Generate random number between 0 and 1.
 */

double random0(void)
{
    long a=16807, m=2147483647, q=127773, r=2836, k;
    int i;
    double rand0;

    k = seed/q;
    seed = a * (seed - k * q) - k * r;
    if (seed < 0) seed += m;
    i = last / (1+(m-1)/NTAB);
    if (i < 0 || i >= NTAB)
	error("random0: error in rand\n");
    last = rtab[i];
    rtab[i] = seed;
    rand0 = (double) ((double) last / (double) m);
    return(rand0);
}
