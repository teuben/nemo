/*  halo.c - inipotential, potential, sqr	*/

/*
 *
 *	Mar 92:		happy gcc2.0		pjt
 *	Oct 93:	get_pattern() support		pjt
 *	Sep 04: float/double			pjt
 */

/*CTEX
 *  {\bf potname=halo
 *       potpars={\it $\Omega,v_0,r_c$}}
 *
 *   rotcurm with m=2
 *
 */


#include <stdinc.h>
#include <potential_float.h>

local double omega = 0.0;	     /* just put to zero until implemented */
local double v0 = 1.0;		/* rotation velocity at large radii */
local double rc = 1.0;		/* core radius */

local double vc, r2;		/* scratch variables */

/*------------------------------------------------------------------------------
 * INIPOTENTIAL: initializes the potential.
 *      input: npar, the number of parameters
 *             par[] an array of npar parameters
 *      If npar=0 defaults are taken (remember to initialize them as static
 *      variables in this file.
 *------------------------------------------------------------------------------
 */

void inipotential (npar, par, name)
int    *npar;
double par[];
char *name;
{
    if (*npar>0) omega = par[0];
    if (*npar>1) v0 = par[1];
    if (*npar>2) rc = par[2];
    
    dprintf (1,"INI_POTENTIAL 3D-Logarithmic potential; name=%s\n",name);
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  V_0, R_c, %f %f\n",v0,rc);
    
    vc = 0.5*v0*v0;
    r2 = rc*rc;
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
	double rad, tmp;
	
	rad = r2 + pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
	*pot = vc * log(rad);
	tmp = -2*vc/rad;
	acc[0] = tmp*pos[0];
	acc[1] = tmp*pos[1];
	acc[2] = tmp*pos[2];
}
