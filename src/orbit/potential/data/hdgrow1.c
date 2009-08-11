/*
 * hdgrow1.c:  growing an isochrone disk in a logarithmic potential
 *
 *        11-aug-2009     Created     PJT/RPO
 */

/*CTEX
 *    {\bf potname=hdgrow1
 *       potpars={\it $\Omega,M_h,r_h,M_d,r_d,M_{ce},\tau_c,M_{de},\tau_d$}} 
 *
 * A Logarithmic Potential (BT, pp.45, eq. 2.54 and eq. 3.77) with
 * an embedded isochrone disk (BT).
 * Both can then grow on different timescales.
 *  
 */

#include <stdinc.h>
#include <potential_float.h>

/* default parameters */

static double omega = 0.0;	/* pattern speed */
static double m_h = 1.0;	/* 1 core mass */
static double r_h = 1.0;	/* 2 core radius */
static double m_d = 1.0;	/* 3 disk mass */
static double r_d = 2.0;	/* 4 disk scale length */
static double m_he = 0.0;	/* 5 extra mass to grow for halo */
static double t_h  = 0.0;	/* 6 timescale for halo to grow */
static double m_de = 0.0;	/* 7 extra mass to grow for disk */
static double t_d  = 0.0;	/* 8 timescale for disk to grow */


static double r2_h, r2_d;
static bool   Qh, Qd;
static double m1, m2;



void inipotential (int *npar, double *par, char *name)
{
  if (*npar>0) omega = par[0];
  if (*npar>1) m_h   = par[1];
  if (*npar>2) r_h   = par[2];
  if (*npar>3) m_d   = par[3];
  if (*npar>4) r_d   = par[4];
  if (*npar>5) m_he  = par[5];
  if (*npar>6) t_h   = par[6];
  if (*npar>7) m_de  = par[7];
  if (*npar>8) t_d   = par[8];
  if (*npar>9) error("hdgrow1: too many parameters given");


  dprintf(1,"INI_POTENTIAL hdgrow1 %s\n",name);
  dprintf(1,"  Parameters : Pattern Speed = %f\n",omega);
  dprintf(1,"  Parameters : Halo:  m,a,me,tau=%f %f %f %f\n",m_h,r_h,m_he,t_h);
  dprintf(1,"  Parameters : Disk:  m,a,me,tau=%f %f %f %f\n",m_d,r_d,m_de,t_d);

  /* only the mass changes, the shapes don't */
  r2_h = r_h*r_h;
  r2_d = r_d*r_d;

  /* for speedier mass computations */
  Qh = (m_he != 0.0);
  Qd = (m_de != 0.0);

  m1 = m_h;
  m2 = m_d;

  par[0] = omega;
}
    

void potential_double (int *ndim, double *pos, double *acc, double *pot, double *time)
{
  double rad, f1, f2, pot1, pot2, mor;
  int    i;

  if (Qh) m1 = m_h + m_he * (1-exp(-(*time/t_h)));
  if (Qd) m2 = m_d + m_de * (1-exp(-(*time/t_d)));

  mor = (r_h > 0) ? m1/r_h : m1;

  rad = sqr(pos[0]) + sqr(pos[1]) + sqr(pos[2]);
  pot1 = mor * log(rad+r2_h);
  f1   = -2.0*mor/(rad+r2_h);

  rad = sqrt(rad+r2_d);
  pot2 = -m2 / (r_d + rad);
  f2   = pot2 / (rad*(r_d+rad));
  
  *pot = pot1 + pot2;
  for (i=0; i<*ndim; i++)
    acc[i] = (f1+f2)*pos[i];
}

