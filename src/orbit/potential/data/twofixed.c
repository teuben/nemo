/*  'two fixed point' potential:
 *      two (fixed) particles, variable mass and position.
 *
 * Two well known limits:
 *	"stark problem"	- 1 of the fixed bodies is far away from the
 *			  inner (close to circular) orbit
 *
 *	"xxx"		- two fixed particles very close together
 *
 *	3-feb-98    Peter Teuben - after a suggestion by Kevin Rauch
 *	22-jun-01   compiler complaints
 *      19-sep-04   float/double
 *
 */

/*CTEX
 *	{\bf potname=twofixed
 *       potpars={\it $\Omega,M_1,x_1,y_1,z_1,M_2,x_2,y_2,z_2$}}
 *
 *
 * This potential is defined by two fixed points, with different masses
 * and positions. Orbits in this potential exhibit a number of interesting
 * properties. One well known limit is the {\tt stark problem}, where one
 * of the two bodies is far from the other and near-circular orbits near
 * the central particles are studied. Another is the limit or two particles
 * near to other and orbits that circumscribe both particles.
 */
 
#include <stdinc.h>
#include <potential_float.h>

static double omega = 0.0;

static double m1 = 1.0;                         /* particle 1 */
static double x1[3] = {  0.0, 0.0, 0.0 };

static double m2 = 0.1;                         /* particle 2 */
static double x2[3] = { -1.5, 0.0, 0.0 };

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
    if (*npar>1) m1    = par[1];
    if (*npar>2) x1[0] = par[2];
    if (*npar>3) x1[1] = par[3];
    if (*npar>4) x1[2] = par[4];
    if (*npar>5) m2    = par[5];
    if (*npar>6) x2[0] = par[6];
    if (*npar>7) x2[1] = par[7];
    if (*npar>8) x2[2] = par[8];
        

    dprintf (1,"INI_POTENTIAL Two Fixed Point problem name=%s\n",name);
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
	dr1 = sqr(pos[0]-x1[0]) + sqr(pos[1]-x1[1]) + sqr(pos[2]-x1[2]);
	dr2 = sqr(pos[0]-x2[0]) + sqr(pos[1]-x2[1]) + sqr(pos[2]-x2[2]);
	*pot -= m1/sqrt(dr1);
	*pot -= m2/sqrt(dr2);
	acc[0] = -m1*(pos[0]-x1[0])/dr1/sqrt(dr1)
                 -m2*(pos[0]-x2[0])/dr2/sqrt(dr2);
	acc[1] = -m1*(pos[1]-x1[1])/dr1/sqrt(dr1)
                 -m2*(pos[1]-x2[1])/dr2/sqrt(dr2);
	acc[2] = -m1*(pos[2]-x1[2])/dr1/sqrt(dr1)
                 -m2*(pos[2]-x2[2])/dr2/sqrt(dr2);
}
