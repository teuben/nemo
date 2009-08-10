/* -*- C -*-                                                                   |
 *-----------------------------------------------------------------------------+
 *                                                                             |
 * nfw2.c                                                                      |
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
 * but with the NH06 style f*cos(2.phi) harmonic distortion                    |
 * this version only works in the XY plane                                     |
 *-----------------------------------------------------------------------------+
 *                                                                             |
 * Versions                                                                    |
 * 0.0   14-aug-2002    created                                           WD   |
 * 0.1   18-nov-2002    converted from C++ to C                           WD   |
 * 0.2   24-may-2005    bit more dprintf() output                        PJT   |
 * 0.3    7-apr-2009    add shapes to play with non-spherical            PJT   |
 * 0.4   23-jun-2009    add HN06 style cos(2.theta) term                 PJT   |
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
 *
 * This version (nfw2) only works in the X-Y plane
 */
#include <stdinc.h>

static double a,ia,fac,f;

void inipotential(int *npar, double *par, string file) {
  double omega, scale, vcmax;
  omega = (*npar>0)? par[0] : 0.;
  scale = (*npar>1)? par[1] : 1.;
  vcmax = (*npar>2)? par[2] : 1.;
  f     = (*npar>3)? par[3] : 0.;
  if (*npar>3) warning("Skipped potential parameters for potname=nfw2 beyond 3");
  a   = scale;
  ia  = 1./a;

  fac = a * vcmax*vcmax / 0.2162165954;
  dprintf (1,"INI_POTENTIAL NFW2 potential\n");
}
//------------------------------------------------------------------------------
#define POTENTIAL(TYPE)							\
void potential_##TYPE(int*NDIM, TYPE*X, TYPE*F, TYPE*P, TYPE*T) {	\
  register double r,r2,fr,fp,ir,x2,y2,f1,phi;				\
  x2 = X[0]*X[0];          	     				        \
  y2 = X[1]*X[1];          	     				        \
  r2 = x2 + y2;          	     				        \
  r  = sqrt(r2);							\
  f1 = 1.0 + f*(x2-y2)/r2;  				        	\
  ir = 1./r,								\
  fr = log(1+ia*r) * ir;						\
  phi -fac*fr;  						       	\
  *P = phi*f1; 						        	\
  fr*= ir;								\
  fr-= ir/(r+a);							\
  fr*=-fac*ir;								\
  fp = 4.0*f*phi/(r2*r2);				        	\
  F[0] = X[0] * (fr - y2*fp);				        	\
  F[1] = X[1] * (fr - x2*fp);				        	\
}

POTENTIAL(float)
POTENTIAL(double)

#undef POTENTIAL
//------------------------------------------------------------------------------
