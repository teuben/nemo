/*
 * mestel.c: Mestel disk - only valid in the XY plane
 *
 *	nov-90  created		PJT
 *	mar-92  happy gcc2.0	PJT
 *	oct-93  get_pattern     PJT
 *      sep-04  float/double    PJT
 */

/*CTEX
 *  {\bf potname=mestel
 *       potpars={\it $\Omega,M,R$}}
 *
 *  \smallskip
 *  See also the description for the Kuzmin disk.
 *
 *  A mestel disk has a flat rotation curve, and up until the radius R,
 *  a surface brightness $ v_o^2 / (2\pi G R) $.
 */
 

#include <stdinc.h>
#include <potential_float.h>

local double omega = 0.0;           /* just put to zero until implemented */
local double mestel_mass = 1.0;
local double mestel_radius = 1.0;

local double vc2;


void inipotential (int  *npar, double *par, char *name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) mestel_mass = par[1];
    if (n>2) mestel_radius = par[2];
    if (n>3) warning("Mestel disk potential: only 3 parameters used");

    dprintf (1,"INIPOTENTIAL Mestel disk %s\n",name);
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  mass, radius, q = %f %f\n", mestel_mass,mestel_radius);
    vc2 = mestel_mass/mestel_radius;
    par[0] = omega;
}
    
void potential_double (int *ndim,double *pos,double *acc,double *pot,double *time)
{
    double tmp, r, r2;
    int    i, n;

    n = MIN(*ndim,2);
    for (i=0, r2=0.0; i<n; i++)
        r2 += sqr(pos[i]);
    if (r2==0.0) {
        *pot = 0.0;         /* a lie */
        for (i=0; i<n; i++)
            acc[i] = 0.0;
    } else {
        r=sqrt(r2);
        *pot = vc2 * log(r);        /*** ??? dimensionality ??? ***/
        tmp = -vc2 / r2;
        for (i=0; i<n; i++)
            acc[i] = tmp*pos[i];
    }
}
