/*
 * op73.c:  Ostriker & Peebles halo, as used in their classic paper
 *		(1973 ApJ 186, 467)
 *
 *
 *      oct-90  created (though never made working)     PJT
 *	mar-92  happy gcc2.0 - actually made it work 	pjt
 *	oct-93  get_pattern
 *      jun-01  fixed potential un-initialized r=0 variables    PJT
 *      sep-04  float/double					PJT
 *
 */

/*CTEX
 *	{\bf potname=op73
 *       potpars={\it $\Omega,M_H,r_c,r_h$}}
 *
 * Ostriker-Peebles 1973 potential
 *	(1973, ApJ {\bf 186}, 467).
 * Their potential is given in the form of the radial force law in the disk
 * plane:
 * $$
 *    F = { M \over R_h^2 }
 *            {  {(R_h+R_c)}^2 \over {(r+R_c)}^2 }
 *            { r \over R_h }
 * $$
 */                     
 
#include <stdinc.h>
#include <potential_float.h>

static double omega = 0.0;      /* just put to zero until implemented */
static double mh = 1.0;         /* mass of halo, within cutoff radius R */
static double rc = 0.1;         /* core radius */
static double rh = 1.0;         /* Radius of halo (cutoff) */

static double rc2, rh2, rhc, mhc;	/* local useful things */

void inipotential (int *npar, double *par, string name)
{
    if (*npar>0) omega = par[0];
    if (*npar>1) mh = par[1];
    if (*npar>2) rc = par[2];
    if (*npar>3) rh = par[3];
    
    dprintf (1,"INI_POTENTIAL OP73 potential");
    dprintf (1,"  Parameters : Pattern Speed = %f \n",omega);
    dprintf (1,"  M_H,  r_c, R= %f %f %f\n",mh,rc,rh);

    rc2 = rc*rc;    /* useful constants */
    rh2 = rh*rh;
    rhc = sqr(1+rc/rh);
    mhc = mh * sqr(rh+rc) / (rh2*rh);
    par[0] = omega;
}
    
void potential_double (int *ndim, double *pos, double *acc, double *pot, double *time)
{
    double rad, f, r;
    int i;
        
    for (i=0, rad=0; i<*ndim; i++)  /* radius */
        rad += sqr(pos[i]);

    if (rad>rh2) {                  /* OUTSIDE the halo: Newton */
      r = sqrt(rad);
      f = mh/(rad*r);             /* radial force / r */
      for (i=0; i<*ndim; i++)     /* cartesian acc's */
	acc[i] = -f*pos[i];
      *pot = -mh/r;               /* potential */
    } else {                        /* INSIDE the halo:  O&P */
      if (rad==0.0) {
	r = f = 0.0;
	for (i=0; i<*ndim; i++)  /* cartesian acc's */            
	  acc[i] = 0.0;
      } else {
	r = sqrt(rad);
	f = mhc / sqr( r + rc );
	for (i=0; i<*ndim; i++)  /* cartesian acc's */            
	  acc[i] = -f*pos[i];
      }
      *pot = -mh/rh * ( 1 - rhc*(log((r+rc)/(rh+rc)) + 
				 rc*(rh-r)/(rh+rc)/(r+rc)) );
    }
}
