/*
 * hubble.c: modified triaxial hubble potential
 *           See also BT, p.40
 *             
 *
 *	7-mar-92  happy gcc2.0 - pjt
 *	19-may-93 proper scaling by mass etc.; allow triaxial and
 *                default to prolate if short axis not given
 *		  -- some optimized version is still possible --
 *	   oct-93 get_pattern
 *	   mar-95 no output of NULL name
 */

/*CTEX
 *  {\bf potname=hubble
 *       potpars={\it $\Omega,M,R,b,c$}} 
 */
 
 
#include <stdinc.h>		    /* NEMO only */

local double omega = 0.0;           /* just put to zero until implemented */
local double hubble_mass = 1.0;	    /* core mass */
local double hubble_radius = 1.0;   /* length scale */
local double hubble_b = 1.0;        /* intermediate b/a axis ratio */
local double hubble_c = 1.0;        /* short c/a axis ratio */

local double eps[3], gr, gr2;       /* useful local's */

#define SQR(x) ((x)*(x))

void inipotential (int *npar, double *par, string name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) hubble_mass = par[1];
    if (n>2) hubble_radius = par[2];
    if (n>3) hubble_b = par[3];
    if (n>4) hubble_c = par[4];
    else hubble_c = hubble_b;		/* prolate if 'c' not given */
    if (n>5) warning("Hubble potential: only 5 parameters used");

    eps[0] = 1.0;               /* for: long axis */
    eps[1] = 1/SQR(hubble_b);	/* for: intermediate axis */
    eps[2] = 1/SQR(hubble_c);	/* for: short axis */
    gr = hubble_mass/hubble_radius;
    gr2 = gr / (hubble_radius);

    dprintf (1,"INIPOTENTIAL Prolate modified Hubble\n");
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  mass, radius, b/a, c/a = %f %f %f %f\n",
                        hubble_mass,hubble_radius,hubble_b,hubble_c);

    par[0] = omega;

}
    
void potential (int *ndim,double *pos,double *acc,double *pot,double *time)
{
    double tmp1, tmp2, tmp3, tmp4, r, r2 = 0.0;
    int    i;
    
    for (i=0; i<*ndim; i++)
        r2 += eps[i]*SQR(pos[i]);
    r2 /= hubble_radius;            /* make it dimensionless */
    
    r=sqrt(r2);
    tmp1 = sqrt(1+r2);
    tmp2 = r + tmp1;
    tmp3 = log(tmp2)/r;             /* ln(r+V1+r^2)/r  */

    *pot = -tmp3 * gr;

    tmp4 = (1/tmp1 - tmp3)/r2 * gr2;

    for (i=0; i<*ndim; i++)
        acc[i] = tmp4*eps[i]*pos[i];
}
