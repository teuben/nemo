/*
 * GRAV.C: routines to compute gravity.
 *   
 *      16-feb-04    for DirectCode
 *
 */

#include "code.h"

/*
 * HACKGRAV: evaluate grav field at a given particle.
 */

local bodyptr p0;			/* body itself  */
local vector pos0;			/* point to evaluate field at */
local real phi0;			/* resulting potential at pos0 */
local vector acc0;			/* resulting acceleration at pos0 */

local vector dr;			/* between gravsub and subdivp */
local real drsq;
local real drabs, phii, mor3;
local vector ai;

void hackgrav(bodyptr p)
{
  p0 = p;					/* exclude p from f.c.      */
  SETV(pos0, Pos(p));				/* set field point          */
  phi0 = 0.0;					/* init potential, etc      */
  CLRV(acc0);
  
  for (p=bodytab; p<bodytab+nbody; p++) {
    if (p == p0) continue;
    SUBV(dr,Pos(p),pos0);
    DOTVP(drsq,dr,dr);

    drsq += eps*eps;                            /* use standard softening   */
    drabs = sqrt(drsq);
    phii = Mass(p) / drabs;
    phi0 -= phii;                               /* add to grav. pot.        */
    mor3 = phii / drsq;
    MULVS(ai, dr, mor3);
    ADDV(acc0, acc0, ai);                       /* add to net accel.        */
  }

  Phi(p0) = phi0;				/* stash the pot.           */
  SETV(Acc(p0), acc0);        			/* and the acceleration     */
}


