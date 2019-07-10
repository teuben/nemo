/*
 * rotcure.c: potential & forces as defined by a rotation curve
 *
 *	20-jun-2019	derived from rotcur1
 *       9-jul-2019     swap par order r0,v0
 */

/*CTEX
 *  {\bf potname=rotcure
 *       potpars={\it $\Omega,v_0,r_0}}
 *	 
 * The forces returned are the axisymmetric forces as defined by
 * a parameterized rotation curve as defined by the turnover point $r_0,v_0$
 * and exp shape.
 */
 

#include <stdinc.h>
#include <potential_float.h>

local double omega = 0.0;
local double v0    = 1.0;
local double r0    = 1.0;
local double m     = 2.0;    /* m=2: logarithmic potential */

void inipotential (int *npar, double *par, string name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) v0    = par[1];
    if (n>2) r0    = par[2];
    if (n>3) warning("Rotcure potential: only 3 parameters usable");
    
    dprintf (1,"INIPOTENTIAL rotcure Potential %s\n",name);
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  V0 = %g  R0 = %g  \n", v0, r0);

    par[0] = omega;     /* return pattern speed again */
}
    
void potential_double (int *ndim, double *pos,double *acc,double *pot,double *time)
{
    real r, r2, v, f;
    int    i;

    for (i=0, r2=0.0; i<2; i++)
        r2 += sqr(pos[i]);
    r=sqrt(r2);

    v = v0 * (1 - exp(-r/r0));

    if (r > 0)
        f = sqr(v/r);
    else
    	f = 0;
    dprintf(2,"r=%g v=%g f=%g\n",r,v,f);

    *pot = 0.0;             /* no potentials... for now */
    acc[0] = -f*pos[0]; 
    acc[1] = -f*pos[1]; 
    acc[2] = 0.0;
}
