/* -*- C -*-                                                                   |
 *-----------------------------------------------------------------------------+
 *                                                                             |
 * Grow_plum.c                                                                 |
 *                                                                             |
 * C code                                                                      |
 *                                                                             |
 * Copyright Walter Dehnen, 2003                                               |
 * e-mail:   wdehnen@aip.de                                                    |
 * address:  Astrophysikalisches Institut Potsdam,                             |
 *           An der Sternwarte 16, D-14482 Potsdam, Germany                    |
 *                                                                             |
 *-----------------------------------------------------------------------------+
 *                                                                             |
 * a growing plummer potential:                                                |
 *                                                                             |
 *                    1-A(t)               A(t)                                |
 * Phi(r) = GM ( ---------------  +  --------------- )                         |
 *               sqrt[r^2 + a^2]     sqrt[r^2 + b^2]                           |
 *                                                                             |
 * where A(t) may have various functional forms, see file timer.h              |
 *                                                                             |
 * parameters are:                                                             |
 *  par[0] = index for A(t)  default: 0        (note: no pattern speed)        |
 *  par[1] = GM              default: 1                                        |
 *  par[2] = a               default: 100                                      |
 *  par[2] = b               default: 1                                        |
 *  par[3] = t0              default: 0                                        |
 *  par[4] = tau             default: 1                                        |
 *-----------------------------------------------------------------------------+
 *                                                                             |
 * Versions                                                                    |
 * 0.0   08-jan-2003    created                                           WD   |
 *                                                                             |
 *----------------------------------------------------------------------------*/
#include <math.h>
#include <stdinc.h>
#include "timer.h"
/*
 * global data
 */
static double aqd, bqd, GMd;                     /* scale rad^2, -total mass */
static float  aqf, bqf, GMf;                     /* scale rad^2, -total mass */

void inipotential(int*npar, double*par, string file) {

  if(*npar < 6) {
    warning("Grow_plum potential: insufficient parameters\n"
	    " the following parameters are recognized:\n"
	    " par[0] = index (%s) [0]\n"
	    " par[1] = total GM of Plummer spheres                   [1]\n"
	    " par[2] = scale radius of outer Plummer sphere          [100]\n"
	    " par[3] = scale radius of inner Plummer sphere          [1]\n"
	    " par[4] = starting time for growth                      [0]\n"
	    " par[5] = tau: time-scale of growth                     [1]\n",
	    timers);
    /*
     * Note, we don't bother here about Peter's idea that the first parameter
     *       shall be a pattern speed (a pattern speed makes no sense with a
     *       spherical model anyway ...)
     */
  }
  init_timer((*npar>0)? (int)(par[0]) : 0,
	     (*npar>3)?       par[4]  : 0.,
	     (*npar>4)?       par[5]  : 1.);
  GMf  = GMd  = (*npar>1)? -fabs(par[1]) : -1.;
  aqf  = aqd  = (*npar>2)? par[2]*par[2] : 100.;
  bqf  = bqd  = (*npar>3)? par[3]*par[3] : 1.;

  if (*npar>6) warning("Skipped potential parameters for Grow_plum beyond 6");
  dprintf (1,"INI_POTENTIAL Growing Plummer sphere\n");
}


void potential_double(int*NDIM, double*X, double*F, double*P, double*T) {
  register double Rq, At, XX, P0, F0;
  if(*T <= t0d) {
    *P   = 0.;
    F[0] = 0.;
    F[1] = 0.;
    if(*NDIM > 2) F[2] = 0.;
    return;
  }
  Rq = X[0]*X[0] + X[1]*X[1];
  if(*NDIM > 2)   Rq += X[2]*X[2];                 /* R^2                    */
  XX = 1./(bqd+Rq);                                /* 1/(R^2 + b^2)          */
  At = timer_double(*T);                           /* A(T)                   */
  P0 = GMd * At * sqrt(XX);                        /* -GM*A/SQRT(R^2+b^2)    */
  *P = P0;                                         /* Potential              */
  F0 = P0 * XX;                                    /* -GM*A/[R^2+b^2]^(3/2)  */
  XX = 1./(aqd+Rq);                                /* 1/(R^2 + a^2)          */
  P0 = GMd * (1-At) * sqrt(XX);                    /* -GM*(1-A)/SQRT(R^2+b^2)*/
  *P+= P0;                                         /* add to Potential       */
  F0+= P0 * XX;                                    /* add to acceleration    */
  F[0] = F0 * X[0];
  F[1] = F0 * X[1];
  if(*NDIM > 2) F[2] = F0 * X[2];
}

void potential_float(int*NDIM, float*X, float*F, float*P, float*T) {
  register float Rq, At, XX, P0, F0;
  if(*T <= t0f) {
    *P   = 0.f;
    F[0] = 0.f;
    F[1] = 0.f;
    if(*NDIM > 2) F[2] = 0.f;
    return;
  }
  Rq = X[0]*X[0] + X[1]*X[1];
  if(*NDIM > 2)   Rq += X[2]*X[2];                 /* R^2                    */
  XX = 1.f/(bqf+Rq);                               /* 1/(R^2 + b^2)          */
  At = timer_float(*T);                            /* A(T)                   */
  P0 = GMf * At * sqrt(XX);                        /* -GM*A/SQRT(R^2+b^2)    */
  *P = P0;                                         /* Potential              */
  F0 = P0 * XX;                                    /* -GM*A/[R^2+b^2]^(3/2)  */
  XX = 1.f/(aqf+Rq);                               /* 1/(R^2 + a^2)          */
  P0 = GMf * (1-At) * sqrt(XX);                    /* -GM*(1-A)/SQRT(R^2+b^2)*/
  *P+= P0;                                         /* add to Potential       */
  F0+= P0 * XX;                                    /* add to acceleration    */
  F[0] = F0 * X[0];
  F[1] = F0 * X[1];
  if(*NDIM > 2) F[2] = F0 * X[2];
}


