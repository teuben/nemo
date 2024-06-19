/*
 * brandt.c: potential & forces as defined by a Brandt rotation curve
 *
 *	20-jun-2019	derived from rotcurm
 *      19-jun-2024     add norm=1 as default now; switched order or (v,r)  v is first now
 */

/*CTEX
 *  {\bf potname=brand1,
 *       potpars={\it $\Omega,v_0,r_0,m$}}
 *	 
 * The forces returned are the axisymmetric forces as defined by
 * a parameterized rotation curve as defined by max $v_0$ around $r_0$
 *
 * See also Brandt, J.C. 1960, ApJ 131, 293. eq. (26)
 */
 

#include <stdinc.h>
#include <spline.h>
#include <table.h>
#include <potential_float.h>

local double omega = 0.0;
local double v0    = 1.0;   /* this is not the vmax if norm=0 (v0*r0 if norm=0) */
local double r0    = 1.0;
local double m     = 2.0;
local int    norm  =   1;  /* non-zero will normalize v0 so it is the vmax */

/* the 1/3 and 2/3 are from Greisen 2009, but Brandt 1960 eq.26 seems to use 1 and 1 */

#if 0
#define C1 1.0
#define C2 1.0
#else
#define C1 1.0/3.0
#define C2 2.0/3.0
#endif

void inipotential (int *npar, double *par, string name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) v0    = par[1];
    if (n>2) r0    = par[2];
    if (n>3) m     = par[3];
    if (n>4) norm  = (int) par[4];
    if (n>5) warning("Brandt potential: only 4 parameters usable");
    
    dprintf (1,"INIPOTENTIAL Brandt potential %s (%g,%g)\n",name,C1,C2);
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    if (norm != 0)
      v0 = v0 / r0;
    dprintf (1,"  V0 = %g  R0 = %g m=%g  norm=%d\n", v0, r0, m, norm);

    par[0] = omega;     /* return pattern speed again */
}
    
void potential_double (int *ndim, double *pos,double *acc,double *pot,double *time)
{
    real r, r2, v, f;
    int    i;

    for (i=0, r2=0.0; i<2; i++)
        r2 += sqr(pos[i]);
    r=sqrt(r2);

    v = v0 * r / pow(C1 + (C2)*pow(r/r0,m),1.5/m);

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
