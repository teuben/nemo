/* -*- C -*-                                                                   |
 *-----------------------------------------------------------------------------+
 *                                                                             |
 * grow_plum.c                                                                 |
 *                                                                             |
 * C code                                                                      |
 *                                                                             |
 * Copyright Walter Dehnen, 2002-2003                                          |
 * e-mail:   wdehnen@aip.de                                                    |
 * address:  Astrophysikalisches Institut Potsdam,                             |
 *           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
 *                                                                             |
 *-----------------------------------------------------------------------------+
 *                                                                             |
 * a growing plummer potential:                                                |
 *                                                                             |
 *               GM A(t)                                                       |
 * Phi(r) = ---------------                                                    |
 *          sqrt[r^2 + a^2]                                                    |
 *                                                                             |
 * where A(t) may have various functional forms, see file timer.h              |
 *                                                                             |
 * parameters are:                                                             |
 *  par[0] = index for A(t)  default: 0        (note: no pattern speed)        |
 *  par[1] = GM              default: 1                                        |
 *  par[2] = a               default: 1                                        |
 *  par[3] = t0              default: 0                                        |
 *  par[4] = tau             default: 1                                        |
 *-----------------------------------------------------------------------------+
 *                                                                             |
 * Versions                                                                    |
 * 0.0   13-nov-2002    created (whilst in Marseille)                     WD   |
 * 0.1   18-nov-2002    debugged                                          WD   |
 * 0.2   08-jan-2003    using timer.h                                     WD   |
 *                                                                             |
 *----------------------------------------------------------------------------*/
#include <math.h>
#include <stdinc.h>
#include "timer.h"
/*
 * global data
 */
static double    Rqd, GMd;                        /* scale rad^2, -total mass */
static float     Rqf, GMf;                        /* scale rad^2, -total mass */

void inipotential(int*npar, double*par, string file) {

  if(*npar < 5) {
    warning("grow_plum potential: insufficient parameters\n"
	    " the following parameters are recognized:\n"
	    " par[0] = index (%s) [0]\n"
	    " par[1] = GM of Plummer sphere                          [1]\n"
	    " par[2] = scale radius of Plummer sphere                [1]\n"
	    " par[3] = starting time for growth                      [0]\n"
	    " par[4] = tau: time-scale of growth                     [1]\n",
	    timers);
    /*
     * Note, we don't bother here about Peter's idea that the first parameter
     *       shall be a pattern speed (a pattern speed makes no sense with a
     *       spherical model anyway ...)
     * Peter's reaction: ieck.... ok, but then you can't integrate orbits...
     */
  }
  init_timer((*npar>0)? (int)(par[0]) : 0,
	     (*npar>3)?       par[3]  : 0.,
	     (*npar>4)?       par[4]  : 1.);
  GMf  = GMd  = (*npar>1)? -fabs(par[1]) : -1.;
  Rqf  = Rqd  = (*npar>2)? par[2]*par[2] : 1.;

  if (*npar>5) warning("Skipped potential parameters for grow_plum beyond 5");
  dprintf (1,"INI_POTENTIAL growing Plummer sphere\n");
}

void potential_double(int*NDIM, double*X, double*F, double*P, double*T) {
  register double XX, P0;
  if(*T <= t0d) {
    *P   = 0.;
    F[0] = 0.;
    F[1] = 0.;
    if(*NDIM > 2) F[2] = 0.;
    return;
  }
  XX = Rqd + X[0]*X[0] + X[1]*X[1];
  if(*NDIM > 2)   XX += X[2]*X[2];                 /* R^2+X^2                */
  XX = 1./XX;                                      /* 1/(R^2+X^2)            */
  P0 = GMd * timer_double(*T) * sqrt(XX);          /* -A(t)*GM/sqrt(R^2+X^2) */
  *P   = P0;                                       /* potential              */
  P0  *= XX;                                       /* -A*GM/[R^2+X^2]^(3/2)  */
  F[0] = P0 * X[0];
  F[1] = P0 * X[1];
  if(*NDIM > 2) F[2] = P0 * X[2];
}

void potential_float(int*NDIM, float*X, float*F, float*P, float*T) {
  register float XX, P0;
  if(*T <= t0f) {
    *P   = 0.f;
    F[0] = 0.f;
    F[1] = 0.f;
    if(*NDIM > 2) F[2] = 0.f;
    return;
  }
  XX = Rqf + X[0]*X[0] + X[1]*X[1];
  if(*NDIM > 2)   XX += X[2]*X[2];                 /* R^2+X^2                */
  XX = 1.f/XX;                                     /* 1/(R^2+X^2)            */
  P0 = GMf * timer_float(*T) * sqrtf(XX);          /* -A(t)*GM/sqrt(R^2+X^2) */
  *P   = P0;                                       /* potential              */
  P0  *= XX;                                       /* -A*GM/[R^2+X^2]^(3/2)  */
  F[0] = P0 * X[0];
  F[1] = P0 * X[1];
  if(*NDIM > 2) F[2] = P0 * X[2];
}
