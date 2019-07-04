/*
 * persic.c: potential & forces as defined by a rotation curve
 *
 *	3-jul-2019	testing plots
 */

/*CTEX
 *  {\bf potname=persic
 *       potpars={\it $\Omega,v_0,r_0,LBS}}
 *	 
 * The forces returned are the axisymmetric forces as defined by
 * a parameterized rotation curve as defined by the point $r_0,v_0$
 * and shape parameter LBS
 */
 

#include <stdinc.h>
#include <potential_float.h>

local double omega = 0.0;
local double r0    = 1.0;
local double v0    = 1.0;
local double lbs   = 1.0;

local double c1,c2,c3,c4;

void inipotential (int *npar, double *par, string name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) v0    = par[1];
    if (n>2) r0    = par[2];
    if (n>3) lbs   = par[3];
    if (n>4) warning("Persic potential: only 3 parameters usable");
    
    dprintf (1,"INIPOTENTIAL Persic Potential %s\n",name);
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  V0 = %g R0 = %g  L_B/L_B* = %g\n", v0, r0, lbs);

    // Persic uses log(L_B*) = 10.4
    c1 = 0.72 + 0.44 * log10(lbs);
    c2 = 1.6 * exp(-0.6*lbs);
    c3 = 1.5 * 1.5 * pow(lbs,0.40);
    c4 = 0.78*0.78;
    
    

    par[0] = omega;     /* return pattern speed again */
}
    
void potential_double (int *ndim, double *pos,double *acc,double *pot,double *time)
{
    real r, r2, v, f, x;
    int    i;

    for (i=0, r2=0.0; i<2; i++)
        r2 += sqr(pos[i]);
    r=sqrt(r2);

    x = r/r0;

    v = v0 * sqrt(c1 * 1.97*pow(x,1.22) / pow(x*x + c4,1.43)     + c2 * x*x / ( x*x + c3) );

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
