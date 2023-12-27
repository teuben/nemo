/*
 * log2.c:  logarithmic potential with m=2 in NDIM=2
 *     jun 2009:  one of the possible examples in NH06
 *     dec 2023:  optionally add q (as from log) only if f=0
 *                this speeds of the log potential by ~10%
 */

/*CTEX
 *    {\bf potname=log2
 *       potpars={\it $\Omega,M_c,r_c,f,q$}} 
 *
 * The Logarithmic Potential (BT, pp.45, eq. 2.54 and eq. 3.77) has
 * been often used in orbit calculations because of its flat rotation
 * curve. This modification adds a constant m=2 harmonic
 * distortion for f>0. The potential is now given by
 *
 * $$
 *    \Phi = {1\over 2} v_0^2
 *                     \ln{ \left( r_c^2 + r^2 \right) } (1 + f \cos(m \phi))
 * $$
 *
 * with $ M_c \equiv {1\over 2} r_c v_0^2 $ defined as the ``core mass''.
 *
 * Note this potential is only defined in the XY plane.
 */

#include <stdinc.h>
#include <potential_float.h>

/* default parameters */

static double omega = 0.0;	/* pattern speed */
static double mc = 1.0;		/* core mass */
static double rc = 1.0;		/* core radius */
static double f2 = 0.0;         /* amplitude of m=2 mode */
static double q  = 1.0;		/* asphericity (q<=1)	*/

static int    m;                /* mode  ; hardcoded for now */
static double mor, rc2, iq2;    /* scratch variables  */

void inipotential (int *npar, double *par, char *name)
{
    if (*npar>0) omega = par[0];
    if (*npar>1) mc = par[1];
    if (*npar>2) rc = par[2];
    if (*npar>3) f2 = par[3];
    if (*npar>4) q = par[4];
    if (f2==0.0)
      m = 0;
    else
      m = 2;                      /* activate the hardcoded m=2 distortion */
      

    dprintf(1,"INI_POTENTIAL Logarithmic potential %s [2D-%s]\n",name, m ? "f" : "q");
    dprintf(1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf(1,"  m_c, r_c, f2, q, m= %f %f %f %f %d\n",mc,rc,f2,q,m);

    if (rc>0)    
        mor = mc/rc;
    else
    	mor = mc;
    rc2 = rc*rc;
    iq2 = 1.0/sqr(q);
    par[0] = omega;

	
}


void potential_double1 (int *ndim, double *pos, double *acc, double *pot, double *time)
{
  double f;
  double rad = rc2;
  rad += pos[0]*pos[0] + pos[1]*pos[1]*iq2;

  *pot = mor * log(rad);
  f = -2.0*mor/rad;

  acc[0] = f*pos[0];
  acc[1] = f*pos[1]*iq2;
  acc[2] = 0.0;
}



void potential_double2 (int *ndim, double *pos, double *acc, double *pot, double *time)
{
  double rad, fr, fp, x2, y2, r2, phi, f1;
	
  x2 = pos[0]*pos[0];
  y2 = pos[1]*pos[1];

  r2 = x2 + y2;
  rad = r2 + rc2;

  phi = mor * log(rad);                  /* unperturbed potential */
  f1 = 1.0 + f2*(x2-y2)/r2;
  *pot = phi * f1;

  fr = -2.0*mor/rad*f1;                  /* f_rad / r */
  fp = 4.0*f2*phi/(r2*r2);               /* f_tan / r */

  acc[0] = pos[0] * (fr - y2*fp);
  acc[1] = pos[1] * (fr + x2*fp);
  acc[2] = 0.0;
}


void potential_double (int *ndim, double *pos, double *acc, double *pot, double *time)
{
  if (m==0)
    potential_double1(ndim,pos,acc,pot,time);
  else
    potential_double2(ndim,pos,acc,pot,time);
}


