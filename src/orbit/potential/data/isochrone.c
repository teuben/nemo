/*
 * isochrone.c: isochrone potential
 *           See also BT, p.38
 *             
 *      nov-90	Created		PJT
 *	mar-92  happy gcc2.0	PJT
 *	oct-93  get_pattern()	PJT
 *      sep-04  float/double	PJT
 */

/*CTEX
 *  {\bf potname=isochrone
 *       potpars={\it $\Omega,M,R$}}
 */

 
#include  <stdinc.h>
#include <potential_float.h>

local double omega = 0.0;           /* just put to zero until implemented */
local double iso_mass = 1.0;
local double iso_radius = 1.0;

local double iso_radius2;

void inipotential (int  *npar, double *par, char *name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) iso_mass = par[1];
    if (n>2) iso_radius = par[2];
    if (n>3) warning("Isochrone potential: only 3 parameters used");

    iso_radius2 = sqr(iso_radius);

    dprintf (1,"INIPOTENTIAL Isochrone potential %s\n",name);
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  mass, radius, = %f %f\n", iso_mass,iso_radius);
    par[0] = omega;
}
    
void potential_double (int *ndim,double *pos,double *acc,double *pot,double *time)
{
    double a, tmp;
    int    i;

    for (i=0, a=iso_radius2; i<*ndim; i++)
        a += sqr(pos[i]);
    a=sqrt(a);

    *pot = -iso_mass / (iso_radius + a);

    tmp = *pot / (a * (iso_radius + a));

    for (i=0; i<*ndim; i++)
        acc[i] = tmp*pos[i];
	
}
