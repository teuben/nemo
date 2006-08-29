/*
 *  FRANDOM: Returns random numbers between
 *        for a given probability distribution f(x),
 *        taken between values a and b.
 *
 *
One normally uses acceptance-rejection methods, which are
in the literature, including Knuth, or acceptance-replacement
methods, which are not.  Even then, there are lots of choices.

One of the best sources of algorithms is the book by
Luc Devroye.  There are also programs in many of the
libraries.  I find most of them at best fair.
 *
 *
 *	25-feb-92  pjt   happy gccV2.0 
 *       7-sep-95  pjt   prototyping
 *	16-feb-97  pjt   fixed for SINGLEPREC 
 *	 8-sep-01  pjt   init_xrandom
 *      29-aug-06  pjt   fixing for prototypes in TESTBED
 *
 */
#include <stdinc.h>
#include <spline.h>

#ifndef MAXN
#define MAXN 1000
#endif

local int n = 0;
local real_proc lastfun=NULL;
local real coeff[MAXN*3], t[MAXN], f[MAXN], cf[MAXN];

double frandom(double a, double b, real_proc fun)
{
   double x, dx;
   int i;

   if (b<=a) return a;         /* do not allow zero or negative intervals */

   if (n==0 || (n>0 && lastfun!=fun)) {
      dprintf(1,"frandom: Setting up a new lookup table\n");
      n = MAXN;
      lastfun = fun;
      dx = (b-a) /(n-1);
      for (i=0; i<n; i++) {       /* get function values at regular interval */
	 t[i] = a + i*dx;
	 f[i] = (*fun)(t[i]);
      }
      spline(coeff,t,f,n);         /* put spline through it for f(t) */
      cf[0] = 0.0;
      for (i=0; i<n-1; i++)       /* get cumulative distribution function */
	 cf[i+1] = cf[i] + dx * (f[i] +
			     dx * (coeff[i] / 2.0 +
			       dx * (coeff[i+n] / 3.0 +
				 dx * (coeff[i+2*n] / 4.0))));
      for (i=0; i<n; i++)         /* normalize it */
	 cf[i] /= cf[n-1];
      spline(coeff,cf,t,n);        /* get inverse spline coef for t(cf) */
   }
   x = xrandom(0.0,1.0);
   x = seval(x, cf,t,coeff, n);
   if (x<a || x>b)
      dprintf(0,"Warning: frandom returns %f; not in [%f,%f]\n",
		x,a,b);
   return x ;
}


#ifdef TESTBED
#include <getparam.h>

string defv[] = {
   "a=-3\n     Left size of interval",
   "b=3\n      Right size of interval",
   "n=20\n     Number of experiments",
   "mode=1\n   Function mode:  1=gauss 2=power",
   "seed=0\n   Random seed",
   "VERSION=2\n 29-aug-06 PJT",
   NULL,
};

string usage="frandom test";

string cvsid="$Id$";

real gauss(real x)
{
   return exp(-0.5*x*x);  /* gauss: mean 0, sigma: 1 */
}

real power(real x)
{
  return pow(x,-2.5);
}

nemo_main()
{
   double a,b;
   int n;
   int mode=getiparam("mode");
   real_proc fun;

   n = getiparam("n");
   a = getdparam("a");
   b = getdparam("b");
   init_xrandom(getparam("seed"));

   if (mode==1)
     fun = (real_proc) gauss;
   else if (mode == 2)
     fun = (real_proc) power;
   else
     error("Bad mode %d",mode);

   while (n-- > 0)
     printf("%g\n", frandom(a,b,fun));

}

#endif
