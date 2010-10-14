/* -*- C -*-                                                                   |
 *-----------------------------------------------------------------------------+
 *                                                                             |
 * nfw.c                                                                       |
 *                                                                             |
 * C code                                                                      |
 *                                                                             |
 * Copyright Walter Dehnen, 2002                                               |
 * e-mail:   wdehnen@aip.de                                                    |
 * address:  Astrophysikalisches Institut Potsdam,                             |
 *           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
 *                                                                             |
 *-----------------------------------------------------------------------------+
 *                                                                             |
 * defines as NEMO potential                                                   |
 *                                                                             |
 *             M0                           ln(1+r/a)                          |
 * rho(r) = ---------    Phi(r) = - 4 Pi M0 ---------                          |
 *          r (r+a)^2                           r                              |
 *                                                                             |
 *-----------------------------------------------------------------------------+
 *                                                                             |
 * Versions                                                                    |
 * 0.0   14-aug-2002    created                                           WD   |
 * 0.1   18-nov-2002    converted from C++ to C                           WD   |
 * 0.2   24-may-2005    bit more dprintf() output                        PJT   |
 * 0.3    7-apr-2009    add shapes to play with non-spherical            PJT   |
 *                                                                             |
 *----------------------------------------------------------------------------*/

/*CTEX
 *
 * The NFW (Navarro,Frank \& White, 1997 ApJ 490, 493) density is given by
 *
 * $$
 *    \rho = {  M_0  \over { r (r+a)^2}}
 * $$
 *
 * and the potential by
 * $$
 *     \Phi = -4 \pi M_0 { \ln{(1+r/a)} \over r }
 * $$
 */
#include <stdinc.h>

static double a,ia,fac,qb,qc,iq2,iq3;

void inipotential(int *npar, double *par, string file) {
  double omega, scale, vcmax;
  omega = (*npar>0)? par[0] : 0.;
  scale = (*npar>1)? par[1] : 1.;
  vcmax = (*npar>2)? par[2] : 1.;
  qb    = (*npar>3)? par[3] : 1.;
  qc    = (*npar>4)? par[4] : qb;
  if (*npar>5) warning("Skipped potential parameters for potname=nfw beyond 5");
  a   = scale;
  ia  = 1./a;
  iq2 = 1./(qb*qb);
  iq3 = 1./(qc*qc);

  fac = a * vcmax*vcmax / 0.2162165954;
  dprintf (1,"INI_POTENTIAL NFW potential s\n");
}
//------------------------------------------------------------------------------
#define POTENTIAL(TYPE)							\
void potential_##TYPE(int*NDIM, TYPE*X, TYPE*F, TYPE*P, TYPE*T) {	\
  register double r,fr,ir;					       	\
  r = X[0]*X[0] + X[1]*X[1]*iq2;					\
  if(*NDIM > 2)   r+= X[2]*X[2]*iq3;					\
  r  = sqrt(r);								\
  ir = 1./r,								\
  fr = log(1+ia*r) * ir;						\
  *P =-fac*fr;								\
  fr*= ir;								\
  fr-= ir/(r+a);							\
  fr*=-fac*ir;								\
  F[0] = fr * X[0];							\
  F[1] = fr * X[1]*iq2;							\
  if(*NDIM>2) F[2] = fr * X[2]*iq3;					\
}

POTENTIAL(float)
POTENTIAL(double)

#undef POTENTIAL
//------------------------------------------------------------------------------
