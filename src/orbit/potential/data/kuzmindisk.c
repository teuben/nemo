/*
 *  kuzmindisk.c:   an infinitesimally thin Kuzmin disk
 *
 *	sep-04	added float/double		PJT
 *
 */


/*CTEX
 *
 * {\bf potname=kuzmin
 *	 potpars={\it $\Omega,M,a$}}
 *
 * Kuzmin (1956) found a closed expression for the potential of 
 * an infinitesimally thin disk with a Plummer potential in the
 * plane of the disk (see also BT pp43, eq. 2-49a and 2-49b):
 * $$
 * \Phi = - {  G M \over {\sqrt{r^2 + (a+{|z|})^2}}}
 * $$
 * and corresponding surface brightness ({\it check units})
 * $$
 * \Sigma = {  {a M} \over {2 \pi {(a^2 + r^2)}^{-3/2}}}
 * $$
 * With $GMa^2 = V_0^2$.
 * This potential is also known as a Toomre n=1 disk, since it
 * was re-derived by Toomre (1963) as part of a series of disks
 * with index $n$, where this disk has $n=1$.
 */
 
 
#include <stdinc.h>
#include <potential_float.h>

local double omega = 0.0;       /* just put to zero until implemented */
local double mass = 1.0;	/* total mass */
local double a = 1.0;		/* scale length */

void inipotential(int *npar, double *par, string name)
{
    if (*npar>0) omega = par[0];
    if (*npar>1) mass = par[1];
    if (*npar>2) a = par[2];
    if (*npar>3) warning("Kuzmindisk: skipped potential parameters beyond 3");

    dprintf (1,"INI_POTENTIAL Kuzmin disk potential\n");
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  mass = %f   scalelength = %f\n",mass, a);

    par[0] = omega;
}

void potential_double (int *ndim,double *pos,double *acc, double *pot,double *time)
{
    int    i;
    double r2, r, f;

    r2 = sqr(pos[0])+sqr(pos[1]);
    r = sqrt(r2);

    if (r2==0.0) {
      *pot = mass/a;
      for (i=0; i<*ndim; i++)
	acc[i] = 0.0;
    } else {
      *pot = -mass/sqrt(r2+sqr(a+ABS(pos[2])));
      // fix the next expression for z <> 0   !!!
      f = (*pot)/(r2+a*a);      
      for (i=0; i<*ndim; i++)
	acc[i] = f*pos[i];
    }
}

