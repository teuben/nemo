/*  triax.c - inipotential, potential, sqr    */

/*
 * triax.c:  procedures for intializing and calculating the forces and
 *             potential of a growing bi/triaxial potential
 *
 *	may 90 - Created for triaxial halo project - PJT
 *
 */

 /*CTEX
 *  {\bf potname=triax}
 * 
 * A growing bi/triaxial potential
 */


#include <stdinc.h>
#include <potential_float.h>
 
#define MPAR 5                       /* omega and remaining */
local double omega = 0.0;           /* just put to zero until implemented */
local double t_cenpot = 1.0;	/* central potential */
local double t_radius = 1.0;	/* length scale */
local double t_epsmin = 0.5;	/* lowest value of b/a */
local double t_tcrit = 10;	/* timescale that b/a varies from 1->epsmin */

local double depsdt;           /* scratch variable (derivative) */
/*------------------------------------------------------------------------------
 * INIPOTENTIAL: initializes the potential.
 *      input: npar, the number of parameters
 *             par[] an array of npar parameters
 *      If npar=0 defaults are taken (remember to initialize them as static
 *      variables in this file.
 *------------------------------------------------------------------------------
 */
void inipotential (int  *npar, double *par, char *name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) t_cenpot = par[1];
    if (n>2) t_radius = par[2];
    if (n>3) t_epsmin = par[3];
    if (n>4) t_tcrit = par[4];
    if (n>5) warning("Triax: too many parameters; only 5 accepted");

    dprintf (1,"INIPOTENTIAL Triax: name=%s\n",name);
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  Cenpot, Radius, Epsmin, Tcrit = %f %f %f %f\n",
		t_cenpot, t_radius, t_epsmin, t_tcrit);		

    if (t_tcrit > 0.0)
      depsdt = (1-t_epsmin)/t_tcrit;
    else
      depsdt = 0.0;

    par[0] = omega;
}
    
/*------------------------------------------------------------------------------
 *  POTENTIAL: the worker routine. Determines at any given point x,y,z) the
 *      forces and potential. Although the naming of parameters suggests
 *      cartesian coordinates, they need not necessarely be, as long as the
 *      the equations of motion can be written in the same way.
 *      Note that this routine is good for 2 and 3D
 *------------------------------------------------------------------------------
 */
void potential_double (int *ndim,double *pos,double *acc,double *pot,double *time)
{
    double tmp, eps;
        
    if (*time < t_tcrit)
        eps = 1.0/sqr(1.0 - depsdt* (*time));
    else
        eps = 1.0/sqr(t_epsmin);

    tmp = 1 + sqr(pos[0]) + eps*(sqr(pos[1])+sqr(pos[2]));

    *pot = t_cenpot*log(tmp);
    tmp = -2*t_cenpot/tmp;

    acc[0] = tmp*pos[0];
    acc[1] = tmp*pos[1]*eps;
    acc[2] = tmp*pos[2]*eps;
}
