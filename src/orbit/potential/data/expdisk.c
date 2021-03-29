/*
 *  expdisk.c:	thin exponential disk
 *
 */


/*CTEX
 *
 * {\bf potname=expdisk
 *	 potpars={\it $\Omega,M,a$}}
 *
 * Exponential disk (BT, pp.77 or BT2008 pp.100)
 * $$
 * \Phi = - {M \over r_d} x \left[ I_0(x)K_1(x) - I_1(x)K_0(x) \right]
 * $$
 *
 * $$
 *     M = 2 \pi \Sigma_0 r_d
 * $$
 * and
 * $$
 *     \Sigma = \Sigma_0 e^{-r/r_d}
 * $$
 */
 
 
#include <stdinc.h>             
#include <potential_float.h>

local double omega = 0.0;       /* just put to zero until implemented */
local double mass = 1.0;	/* total mass */
local double a = 1.0;		/* scale length */

local double alpha;

extern double bessi0(double), bessk0(double), bessi1(double), bessk1(double);

void inipotential (int *npar, double *par, char *name)
{
    if (*npar>0) omega = par[0];
    if (*npar>1) mass = par[1];
    if (*npar>2) a = par[2];
    if (*npar>3) warning("Expdisk: skipped potential parameters beyond 3");

    dprintf (1,"INI_POTENTIAL Exponential disk potential\n");
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  mass = %f   scalelength = %f\n",mass, a);

    alpha = 1.0/a;
    par[0] = omega;
}

void potential_double (int *ndim,double *pos,double *acc,double *pot,double *time)
{
    int    i;
    double r2, r, x, i0, k0, i1, k1, f;
        
    for (i=0, r2=0; i<*ndim; i++)              /* radius - squared */
        r2 += sqr(pos[i]);
    r = sqrt(r2);
    x = 0.5*alpha*r;

    if (r2==0.0) {
        *pot = -mass;
        for (i=0; i<*ndim; i++)
            acc[i] = 0.0;
    } else {
        i0=bessi0(x);
        k0=bessk0(x);
        i1=bessi1(x);
        k1=bessk1(x);
        *pot = -mass*x*(i0*k1-i1*k0);
        f = -0.5*sqr(alpha)*alpha*mass*(i0*k0-i1*k1);
        for (i=0; i<*ndim; i++)
            acc[i] = f*pos[i];
    }
}

