/*
 * HOM:  homogeneous sphere, with exponentially decaying mass
 *		M(t) = M(0) * exp(-time/tau)
 *
 *	feb-90  Created for open cluster project - no real boundary though
 *	mar-92	Happy gcc2.0 - re-introduced rmax boundary          pjt
 *	oct-93  get_pattern, and fixed center point bug             pjt
 *	sep-04  float/double					    pjt
 */

/*CTEX
 *  {\bf potname=hom
 *       potpars={\it $\Omega,M,R,\tau$}}
 */
 

#include <stdinc.h>
#include <potential_float.h>

local  double omega= 0.0;	     /* just put to zero until implemented */
local  double mtot=  1.0; 	     /* Mass of sphere */
local  double rmax=  1.0;	     /* Radius of sphere */
local  double tau=   0.0;            /* decay time; 0.0 means no decay */

local  double w2, rmax2, poffset;

void inipotential (int  *npar, double *par, char *name)
{
    if (*npar>0) omega = par[0];
    if (*npar>1) mtot = par[1];  
    if (*npar>2) rmax = par[2];
    if (*npar>3) tau = par[3];
    
    rmax2 = rmax*rmax;
    w2 = mtot/(rmax2*rmax);     /* harmonic coefficient: psi = w2*r^2/2 */
    poffset = 1.5 * w2*rmax2;   /* offset harmonic and keplerian */

    dprintf (1,"INI_POTENTIAL decaying homsph\n");
    dprintf (1,"  mtot=%f rmax=%f tau=%f\n",mtot,rmax,tau);
    par[0] = omega;
}
    
void potential_double (int *ndim,double *pos,double *acc,double *pot,double *time)
{
    double mass, r2, r, ft;
    int    i;

    for (r2=0, i=0; i<*ndim; i++)
        r2 += sqr(pos[i]);

    if (r2==0.0) {
    	*pot = -poffset;
	for (i=0; i<*ndim; i++) acc[i] = 0.0;
	return;
    }
	
    if (tau == 0.0)
        ft = w2;
    else
        ft = w2 * exp(-(*time)/tau);  /* time decay factor */


    if (r2<=rmax2) {                    /* inside: constant density */
        *pot = 0.0;
        for (i=0; i<*ndim; i++) {
            *pot += ft*sqr(pos[i]);
            acc[i] = -ft*pos[i];
        }
        *pot *= 0.5;
        *pot -= poffset;
    } else {                            /* outside: constant mass */
        r = sqrt(r2);
        mass = ft*rmax*rmax2;           
        ft = mass/(r2*r);               /* radial force / r */
        for (i=0; i<*ndim; i++)
            acc[i] = -ft*pos[i];
        *pot = -mass/r;
    }
}
