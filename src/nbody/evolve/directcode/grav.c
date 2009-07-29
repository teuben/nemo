/*
 * GRAV.C: routines to compute gravity.
 *   
 *      16-feb-04    cloned from hackcode1 for DirectCode
 *      29-jul-09    eps < 0 allowed for pseudo-newtonian
 *
 */

#include "code.h"


/* must define either USE_LOCAL *or* USE_STACK here  */

/*
 * HACKGRAV: evaluate grav field at a given particle.
 *           for a Direct N-body code this is awfully simple:
 *           loop over all particles, but yourself, and accumulate
 *           the forces, and potential
 */



void hackgrav(bodyptr p)
{
  bodyptr p0;			/* body itself  */
  vector pos0;			/* point to evaluate field at */
  real phi0;			/* resulting potential at pos0 */
  vector acc0;			/* resulting acceleration at pos0 */
  vector dr;			/* between gravsub and subdivp */
  real drsq;
  real drabs, phii, mor3;
  vector ai;

  p0 = p;					/* exclude p from f.c.      */
  SETV(pos0, Pos(p));				/* set field point          */
  phi0 = 0.0;					/* init potential, etc      */
  CLRV(acc0);

  for (p=bodytab; p<bodytab+nbody; p++) {
    if (p == p0) continue;
    SUBV(dr,Pos(p),pos0);
    DOTVP(drsq,dr,dr);

    if (eps >= 0.0) {
      drsq += eps*eps;                         /* use standard softening   */
      drabs = sqrt(drsq);
      phii = Mass(p) / drabs;
      phi0 -= phii;                               /* add to grav. pot.        */
      mor3 = phii / drsq;
      MULVS(ai, dr, mor3);
      ADDV(acc0, acc0, ai);                       /* add to net accel.        */
    } else if (eps < 0) {
      drabs = sqrt(drsq) + eps;
      if (drabs < 0) error("PN violation at time=%g",tnow);
      drsq = drabs*drabs;                         
      phii = Mass(p) / drabs;
      phi0 -= phii;                               /* add to grav. pot.        */
      mor3 = phii / drsq;
      MULVS(ai, dr, mor3);
      ADDV(acc0, acc0, ai);                       /* add to net accel.        */
      
    }
  }

  Phi(p0) = phi0;				/* stash the pot.           */
  SETV(Acc(p0), acc0);        			/* and the acceleration     */
}


