/*
 * VECTMATH.C: source code for vector/matrix operations.
 * Joshua Barnes  1 December 1986  Princeton, NJ.
 *	29-jul-92	convex fix		PJT
 *	22-jan-95	ansi prototypes - functions return real ????
 *      23-jun-01       ZENOisms
 */

#include <stdinc.h>

#if defined(DOUBLEPREC)
#  define _dotvp	_dotvp_d
#  define _absv		_absv_d
#  define _distv	_distv_d
#  define _tracem	_tracem_d
#endif

#if defined(SINGLEPREC)
#  define _dotvp	_dotvp_f
#  define _absv		_absv_f
#  define _distv	_distv_f
#  define _tracem	_tracem_f
#endif

double _dotvp(real *v, real *u, int n)
{
    real s;

    s = 0.0;
    while (--n >= 0)
	s += (*v++) * (*u++);
    return (s);
}

double _absv(real *v, int n)
{
    real s;

    s = 0.0;
    while (--n >= 0) {
	s += (*v) * (*v);
	v++;
    }
    return (sqrt(s));
}

double _distv(real *v, real *u, int n)
{
    real d, s;

    s = 0.0;
    while (--n >= 0) {
	d = (*v++) - (*u++);
	s += d * d;
    }
    return (sqrt(s));
}

double _tracem(real *p, int n)
{
    double s;
    register int i;

    s = 0.0;
    for (i = n; --i >= 0; ) {
	s += (*p);
	p += (n + 1);			/* next diag. element */
    }
    return (s);
}
