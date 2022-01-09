/*  2-body potential
 *
 *	 9-jan-2022 Peter Teuben - after a gradmap talk by Geme
 *
 */

/*CTEX
 *	{\bf potname=twobody
 *       potpars={\it $\Omega,m,r,\epsilon$}}
 *
 *
 * This is the well known 2-body potential in the rotating frame
 * of reference.
 */
 
#include <stdinc.h>
#include <potential_float.h>

static double omega = 0.0;                      /* will be re-computed if < 0 */

static double m1 = 1.0;                         /* particle 1 */
static double x1[3] = { -0.1, 0.0, 0.0 };

static double m2 = 0.1;                         /* particle 2 */
static double x2[3] = {  1.0, 0.0, 0.0 };

static double eps = 0, eps2;                    /* classic softening */

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
    if (*npar>0) omega = par[0];
    if (*npar>1) m2    = par[1];
    if (*npar>2) x2[0] = par[2];
    if (*npar>3) eps   = par[3];

    if (omega < 0) {   // special case in the matching rotating frame of reference
      omega = 1/sqrt(x2[0])/x2[0]/(1+m2);
      x1[0] = -x2[0] * m2/m1;
    }
    eps2 = eps*eps;

    dprintf (1,"INI_POTENTIAL Two Body problem name=%s\n",name);
    dprintf (1," omega=%g   eps=%g\n", omega,eps);
    dprintf (1," 1: m=%g  pos=%g %g %g\n",m1,x1[0],x1[1],x1[2]);
    dprintf (1," 2: m=%g  pos=%g %g %g\n",m2,x2[0],x2[1],x2[2]);
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
  double dr1, dr2;
	
    *pot = 0.0;
    dr1 = sqr(pos[0]-x1[0]) + sqr(pos[1]-x1[1]) + sqr(pos[2]-x1[2]) + eps2;
    dr2 = sqr(pos[0]-x2[0]) + sqr(pos[1]-x2[1]) + sqr(pos[2]-x2[2]) + eps2;
    *pot -= m1/sqrt(dr1);
    *pot -= m2/sqrt(dr2);
    acc[0] = -m1*(pos[0]-x1[0])/dr1/sqrt(dr1)
             -m2*(pos[0]-x2[0])/dr2/sqrt(dr2);
    acc[1] = -m1*(pos[1]-x1[1])/dr1/sqrt(dr1)
             -m2*(pos[1]-x2[1])/dr2/sqrt(dr2);
    acc[2] = -m1*(pos[2]-x1[2])/dr1/sqrt(dr1)
             -m2*(pos[2]-x2[2])/dr2/sqrt(dr2);
}
