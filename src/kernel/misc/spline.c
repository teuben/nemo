/*
 * SPLINE.C: functions to evaluate cubic spline approximations.
 *
 *        spline()	set up spline coeffs
 * real   seval()	spline evaluator
 * real   spldif()	1st derivative of spline
 * real   spldif2()     2nd derivative of spline
 *
 *        splsub()      (local) worker routine for spline
 *      interval()      (local) worker routine to find interval
 *
 * Reference: Forsythe, Malcolm & Moler, "Computer Methods for
 *	      Mathematical Computations", pp. 76-79.
 *
 *	sep-90   added spldif2() for second derivs			PJT
 *	oct-90	 reworked line 45 because of cc compiler bug on 3b1	PJT
 *               (it cannot handle double array assigments  a[]=b[]=c)
 *   25-feb-92   happy gcc2.0                                           PJT
 *   29-jul-92   remember last invocation accross seval, spldif, 
 *               and spldif2  using interval()                          PJT
 *    5-nov-93   TOOLBOX version to find derivatives                    pjt
 *   26-feb-94   1.1 : ansi headers TESTBED incomplete
 *    1-mar-94   1.2 : double is now real				pjt
 *   15-apr-98   2.0 : finishing TOOLBOX version                        pjt
 *   16-apr-98   2.1 : TOOLBOX in its own tabspline.c (see kernel/tab)  pjt
 *   20-jun-01   gcc3                                                   pjt
 */

#include <stdinc.h>

local int lasti = -1;       /* remember index of last invocation */

local void splsub(real*, real*, real*, real*, real*, int);
local int interval(real, real *, int);

/*
 * SPLINE: compute cubic spline coefs.
 *	real coef[];          storage for computed coefficients 
 * 	real x[], y[];        data points to interpolate between 
 *	int n;                number of data points 
 */

void spline(real *coef, real *x, real *y, int n)
{
    splsub(&coef[0], &coef[n], &coef[2*n], x, y, n);
    lasti = -1;         /* at least reset the last index */
}

/*
 * SPLSUB: local worker routine that sets up the spline coeffs
 */

local void splsub(
    real b[], real c[], real d[],   /* 1st to 3rd order polynomial coefs */
    real x[], real y[],
    int n)
{
    int i;
    real t1, tn, t;

    if (n < 3) {
        error("spline: n=%d lt 3",n);
        return;
    }
    d[0] = x[1] - x[0];
    c[1] = (y[1] - y[0]) / d[0];
    for (i = 1; i <= n-2; i++) {
        d[i] = x[i+1] - x[i];
        b[i] = 2.0 * (d[i-1] + d[i]);
        c[i+1] = (y[i+1] - y[i]) / d[i];
        c[i] = c[i+1] - c[i];
    }
    b[ 0 ] = - d[ 0 ];
    b[n-1] = - d[n-2];
    if (n == 3) {
        c[0] = 0.0; c[n-1] = 0.0;     /* 3b1 compiler chocked if one stmt. */
    } else {
        t1 = c[ 2 ] / (x[ 3 ] - x[ 1 ]) - c[ 1 ] / (x[ 2 ] - x[ 0 ]);
        tn = c[n-2] / (x[n-1] - x[n-3]) - c[n-3] / (x[n-2] - x[n-4]);
        c[ 0 ] =   t1 * d[ 0 ]*d[ 0 ] / (x[ 3 ] - x[ 0 ]);
        c[n-1] = - tn * d[n-2]*d[n-2] / (x[n-1] - x[n-4]);
    }
    for (i = 1; i <= n-1; i++) {
        t = d[i-1] / b[i-1];
        b[i] = b[i] - t * d[i-1];
        c[i] = c[i] - t * c[i-1];
    }
    c[n-1] = c[n-1] / b[n-1];
    for (i = n-2; i >= 0; i--)
        c[i] = (c[i] - d[i] * c[i+1]) / b[i];
    b[n-1] = (y[n-1] - y[n-2]) / d[n-2] + d[n-2] * (c[n-2] + 2 * c[n-1]);
    for (i = 0; i <= n-2; i++) {
        b[i] = (y[i+1] - y[i]) / d[i] - d[i] * (c[i+1] + 2 * c[i]);
        d[i] = (c[i+1] - c[i]) / d[i];
        c[i] = 3 * c[i];
    }
    c[n-1] = 3 * c[n-1];
    d[n-1] = d[n-2];
}

/*
 * SEVAL: evaluate cubic spline interpolation, using Horner's rule
 *	real x0;              point at which to compute function 
 *	real x[], y[];        data points to interpolate between 
 *	real coef[];          storage for computed coefficients
 *	int n;                number of data points 
 */

real seval(real x0, real *x, real *y, real *coef, int n)
{
    int i;
    real u;

    i = interval(x0,x,n);
    u = x0 - x[i];
    return y[i] + u * (coef[i] + u * (coef[n+i] + u * coef[2*n+i]));
}

/*
 * SPLDIF: evaluate derivative of cubic spline.
 *	real x0;              point at which to compute function
 *	real x[], y[];        data points to interpolate between
 *	real coef[];          storage for computed coefficients
 *	int n;                number of data points
 */

real spldif(real x0, real *x, real *y, real *coef, int n)
{
    int i;
    real u;

    i = interval(x0,x,n);
    u = x0 - x[i];
    return coef[i] + u * (2.0*coef[i+n] + u * 3.0*coef[i+2*n]);
}

/*
 * SPLDIF2: evaluate second derivative of cubic spline
 *	real x0;              point at which to compute function
 *	real x[], y[];        data points to interpolate between
 *	real coef[];          storage for computed coefficients
 *	int n;                number of data points
 */

real spldif2(real x0, real *x, real *y, real *coef, int n)
{
    int i;
    real u;

    i = interval(x0,x,n);
    u = x0 - x[i];
    return 2.0*coef[i+n] + u * 6.0*coef[i+2*n];
}

/*
 * INTERVAL: find the correct interval, optimized to first try
 *           the last interval tried. Values beyond the two edges
 *           use the edge data
 */
local int interval(real u, real *x, int n)
{
    int i, j, k;

    if (lasti>=0 && lasti<n-1 && x[lasti]<=u && u<x[lasti+1]) /* check last */
        return lasti;
    else if (u < x[0])             /* check extrapolation left */
        lasti = 0;
    else if (u >= x[n-1])          /* and right */
        lasti = n-1;
    else {                          /* hence u wasn't in last interval */
        i = 0;                      /* and a binary search is performed */
        k = n;
        while (i+1 < k) {
            j = (i + k) / 2;
            if (x[j] <= u)
                i = j;
            else
                k = j;
        }
        lasti = i;                  /* and remember this interval       */
    }
    return lasti;
}


#ifdef TESTBED

#if !defined(N)
#define N       11
#endif

#define fun(x0)  (sqrt(x0) * cos(PI*x0))
#define fpr(x0)  (0.5*cos(PI*x0) / sqrt(x0) - sqrt(x0) * PI*sin(PI*x0))
#define fpr2(x0) (-0.75*cos(PI*x0)/(sqrt(x0)*x0) - PI*cos(PI*x0)/(2*sqrt(x0)) \
    - PI*sin(PI*x0)/(2*sqrt(x0)) - PI*PI*sqrt(x0)*cos(PI*x0))

main()
{
    int i;
    real x[N], y[N], coef[3*N], x0;

    for (i = 0; i < N; i++) {
        x[i] = i / (N - 1.0);
        y[i] = fun(x[i]);
    }
    spline(coef, x, y, N);
    printf("\n%12s%12s%12s%12s%12s\n", "x", "y", "b", "c", "d");
    for (i = 0; i < N; i++) {
        printf("%12.6f%12.6f%12.6f%12.6f%12.6f\n",
               x[i], y[i], coef[i], coef[N+i], coef[2*N+i]);
    }
    printf("\n");
    for ( ; ; ) {
        printf("x0: ");
        scanf("%lf", &x0);
	printf("        %12s\t%12s\t%12s\n", 
                                    "f(x)", "f'(x)","f''(x)");
	printf("exact   %12.6f\t%12.6f\t%12.6f\n", 
                                    fun(x0), fpr(x0),fpr2(x0));
	printf("spline  %12.6f\t%12.6f\t%12.6f\n",
	                            seval(x0, x, y, coef, N),
	                            spldif(x0, x, y, coef, N),
                                    spldif2(x0, x, y, coef, N));
    }
}

#endif
