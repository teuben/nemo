/*
 *  kuzmindisk.c:   a thin Kuzmin disk
 *
 */


/*CTEX
 *
 * {\bf potname=kuzmin
 *	 potpars={\it $\Omega,M,a$}}
 *
 * Kuzmin (1956) found a closed expression for the potential of 
 * an infinitesimally thin disk with a Plummer potential in the
 * plane of the disk:
 * $$
 * \Phi = - {  G M \over {\sqrt{r^2 + (a+\abs{z})^2}}}
 * $$
 * and corresponding surface brightness (check units)
 * $$
 * \Sigma = - {  {G M} \over {2 \pi}} {\sqrt{1 + (r/a)^2}^{3/2}}
 * $$
 */
 
 
#include <stdinc.h>                     /* standard Nemo include */

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

void potential (int *ndim,double *pos,double *acc, double *pot,double *time)
{
    int    i;
    double r2, r, x, i0, k0, i1, k1, f;

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

