/*
 * log2.c:  logarithmic potential with m=2 in NDIM=2
 *     jun 2009:  one of the possible examples in NH06 
 */

/*CTEX
 *    {\bf potname=log2
 *       potpars={\it $\Omega,M_c,r_c,f,m$}} 
 *
 * The Logarithmic Potential (BT, pp.45, eq. 2.54 and eq. 3.77) has
 * been often used in orbit calculations because of its flat rotation
 * curve. This modification adds a constant m=2 harmonic
 * distortion. The potential is now given by
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
static double f  = 0.0;		/* asphericity (q<=1)	*/
static int    m  = 2;           /* mode  ; hardcoded for now */

static double mor, rc2;         /* scratch variables  */

void inipotential (int *npar, double *par, char *name)
{
    if (*npar>0) omega = par[0];
    if (*npar>1) mc = par[1];
    if (*npar>2) rc = par[2];
    if (*npar>3) f = par[3];
    m = 2;                      /* hardcoded below */

    dprintf(1,"INI_POTENTIAL Logarithmic potential %s\n",name);
    dprintf(1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf(1,"  m_c, r_c, f, m= %f %f %f %d\n",mc,rc,f,m);

    if (rc>0)    
        mor = mc/rc;
    else
    	mor = mc;
    rc2 = rc*rc;
    par[0] = omega;
}
    

void potential_double (int *ndim, double *pos, double *acc, double *pot, double *time)
{
  double rad, fr, fp, x2, y2, r2, phi, f1;
  int    i;
  double theta;
	
  x2 = pos[0]*pos[0];
  y2 = pos[1]*pos[1];

  r2 = x2 + y2;
  rad = r2 + rc2;

  phi = mor * log(rad);                  /* unperturbed potential */
  f1 = 1.0 + f*(x2-y2)/r2;
  *pot = phi * f1;

  fr = -2.0*mor/rad*f1;                  /* f_rad / r */
  fp = 4.0*f*phi/(r2*r2);                /* f_tan / r */

  acc[0] = pos[0] * (fr - y2*fp);
  acc[1] = pos[1] * (fr + x2*fp);
  acc[2] = 0.0;
}
