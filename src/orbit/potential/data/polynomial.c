/*
 * polynomial.c:  (spherical) polynomial potential
 *
 *	jul-2002b	created for N4313 (Apj 437,162)
 *      jan-2003        MAXN -> MAXPOW
 *      sep-2004	float/double
 */

/*CTEX
 *	{\bf potname=polynomial
 *       potpars={\it $\Omega,a0,a1,a2,a3,....$}}
 *
 *  Polynomial potential
 *
 * $$
 *    \Phi =  a_0 + a_1 r + a_2 r^2 + ....  a_N r^N
 * $$
 *
 *  where any unused coefficients will be set to 0. Up to 16 (defined
 *  as MAXPOW) can be used.
 */                     
  
#include <stdinc.h>
#include <vectmath.h>
#include <potential_float.h>

#ifndef MAXPOW
#define MAXPOW  16
#endif

local double omega = 0.0;

local double a[MAXPOW];
local int na = 0;

void inipotential (int *npar, double *par, string name)
{
    int i,n;

    n = *npar;
    if (n>0) omega = par[0];
    na = n-1;
    if (na < 1) error("polynomial potential: needs coefficients");
    if (na > MAXPOW) error("polynomial potential: MAXPOW=%d, you requested %d",MAXPOW,na);
    for (i=0; i<na; i++)
      a[i] = par[i+1];
    dprintf(0,"polynomial: Found %d polynomial coefficients %g .. %g\n",
	    na,a[0],a[na-1]);
    par[0] = omega;
}


void potential_double (int *ndim,double *pos,double *acc,double *pot,double *time)
{
    double r2 =  pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
    double r = sqrt(r2);
    double p,f;
    int i;


    p = f = 0.0;
    for (i=na; i>0; i--) {
      p = p*r + a[i-1];
      if (i>0) f = f*r + i*a[i-1];
    }
    f /= r;
    dprintf(1,"f=%g p=%g\n",f,p);

    *pot = p;
    acc[0] = f*pos[0];
    acc[1] = f*pos[1];
    acc[2] = f*pos[2];
}

