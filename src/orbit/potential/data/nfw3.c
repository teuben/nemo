/*
 * the NFW (Navarro,Frank \& White, 1997 ApJ 490, 493) profile rotcurshape
 * since Dehnen's version doesn't use a concentration parameter ??
 * Note order of parameters a bit odd:
 *        v0 = peak rotation curve
 *        a = scale length
 *        c = concentration parameter (r200/a)
 *
 * this was a quick hack to compare to rotcurshape, and only returns
 * forces, no potentials.
 *
 * 24-may-2005   PJT    Hacked Up Q&D
 * 13-jul-2007   PJT    added a 'q' to check Rachel's thesis
 *
 * 
 */


#include <stdinc.h>

static double omega = 0.0;
static double v0=1;
static double a=1;
static double c=1;
static double q=1;  /* axis ratio in the XY plane */

static double ia, lnc, iq2[3];


void inipotential(int *npar, double *par, string file) 
{
  int n = *npar;
  if (n>0) omega = par[0];
  if (n>1) v0 = par[1];
  if (n>2) a = par[2];
  if (n>3) c = par[3];
  if (n>4) q = par[4];
  dprintf (0,"INI_POTENTIAL NFW3 potential v0=%g a=%g c=%g q=%g\n",v0,a,c,q);

  ia = 1/a;
  lnc = log(1+c) - c/(1+c);
  iq2[0] = 1.0;
  iq2[1] = 1.0/sqr(q);
  iq2[2] = 1.0/sqr(q);
}
//------------------------------------------------------------------------------
#define POTENTIAL(TYPE)							\
void potential_##TYPE(int*NDIM, TYPE*X, TYPE*F, TYPE*P, TYPE*T) {	\
  register double x,r,r2,fr,ir;					       	\
  r2 = X[0]*X[0]*iq2[0] + X[1]*X[1]*iq2[1];				\
  if(*NDIM > 2)   r2+= X[2]*X[2]*iq2[2];				\
  if (r2==0.0) {                \
    F[0] = F[1] = F[2] = 0.0;   \
    *P = 0.0;                   \
    return;                     \
  }                             \
  r  = sqrt(r2);							\
  ir = 1./r;								\
  x = c*ia*r;                                                           \
  fr = log(1+x) - x/(1+x);						\
  fr *= -(a*v0*v0)/(r*r2*lnc);                                          \
  F[0] = fr * X[0] * iq2[0];						\
  F[1] = fr * X[1] * iq2[1];						\
  if(*NDIM>2) F[2] = fr * X[2] * iq2[2];				\
  *P = 0.0;                                                             \
}

POTENTIAL(float)
POTENTIAL(double)

#undef POTENTIAL
//------------------------------------------------------------------------------
