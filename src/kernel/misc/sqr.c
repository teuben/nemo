/*
 * (r)sqr  compute x*x
 * (r)qbe  compute x*x*x.
 * (r)dex  compute inverse log^10: 10**x
 *
 *	dark ages	created			JEB
 *	22-jan-95	ansi proto		PJT
 *      22-jun-01       ZENO stuff added for mixed precision
 */

#include <stdinc.h>
#include <math.h>
#include <mathfns.h>

double sqr(double x)
{
    return x*x;
}

double qbe(double x)
{
    return x*x*x;
}

double dex(double x)
{
    /* extern double pow(double,double); */
    return pow(10.0, x);
}

/* ---  mixed math versions --- */

#if !defined(DOUBLEPREC)

real rsqr(real x)
{
    return (x * x);
}

real rqbe(real x)
{
    return (x * x * x);
}

real rdex(real x)
{
    return (rexp(M_LN10 * x));
}

#endif
