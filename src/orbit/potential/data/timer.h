/* -*- C -*-                                                                   |
 *-----------------------------------------------------------------------------+
 *                                                                             |
 * timer.h                                                                     |
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
 * timer functions for external potentials                                     |
 *                                                                             |
 * various functional forms are supported:                                     |
 *                                                                             |
 * index   name         functional form                                        |
 * -----+--------------+-----------------------------------------------------  |
 *  0   | adiabatic    |  a smooth transition from 0 at t<t0 to 1 at t>t0+tau  |
 *  1   | saturate     |  1-exp([t-t0]/tau)                                    |
 *  2   | quasi-linear |  sqrt(x^2+1) - 1, x=[t-t0]/tau for t>t0, 0 for t<t0   |
 *  3   | linear       |  [t-t0]/tau for t>t0, 0 for t<t0                      |
 *                                                                             |
 * initialization: void   init_timer(int index, double t0, double tau);        |
 * call:           double timer_double(double t);                              |
 * call:           float  timer_float (float  t);                              |
 *                                                                             |
 *-----------------------------------------------------------------------------+
 *                                                                             |
 * Versions                                                                    |
 * 0.0   08-jan-2003    created                                           WD   |
 *                                                                             |
 *----------------------------------------------------------------------------*/
#ifndef included_timer_h
#define included_timer_h
static double t0d=0., t1d=1., taud=1., itaud=1.;   /* t_0, t_1, tau, 1/tau    */
static float  t0f=0., t1f=1., tauf=1., itauf=1.;   /* t_0, t_1, tau, 1/tau    */
static int    timer=0;                             /* index of timer function */
static const char* timers = "adiabatic, saturate, quasi-linear, linear";

/*
 * timer functions (assume t > t0)
 */
inline double square_double(double x) { return x*x; }
inline float  square_float (float  x) { return x*x; }
/*
 *      3   5    5   3    15       1           t - t0
 * A = --- x  - --- x  + ---- x + ---;  x = 2 -------- - 1;
 *     16        8        16       2            tau
 *
 * for t in [t0,t0+tau]
 */
inline double adiabatic_double(double t) {
  register double xi, xq;
  if(t >= t1d) return 1.;
  xi  = 2*(t-t0d)*itaud-1;
  xq  = xi*xi;
  xi *= ((0.1875*xq - 0.625)*xq + 0.9375);
  return xi + 0.5;
}

inline float adiabatic_float(float t) {
  register float xi, xq;
  if(t >= t1f) return 1.f;
  xi  = 2*(t-t0f)*itauf-1;
  xq  = xi*xi;
  xi *= ((0.1875f*xq - 0.625f)*xq + 0.9375f);
  return xi + 0.5f;
}
/*
 *              t - t0
 * A = 1 - exp(--------);
 *               tau
 */
inline double saturate_double(double t) {
  return 1.  - exp((t0d-t)*itaud);
}

inline float saturate_float(float t) {
  return 1.f - expf((t0f-t)*itauf);
}
/*
 *            2                  t - t0
 * A = sqrt( x  + 1 ) - 1;  x = --------;
 *                                tau
 */
inline double quasilinear_double(double t) {
  return sqrt(square_double((t-t0d)*itaud)+1.)-1.;
}

inline float quasilinear_float(float t) {
  return sqrtf(square_float((t-t0f)*itauf)+1.f)-1.f;
}
/*
 *          t - t0
 * A = x = --------;
 *           tau
 */
inline double linear_double(double t) {
  return t<t0d? 0. : (t-t0d)*itaud;
}

inline float inear_float(float t) {
  return t<t0f? 0. : (t-t0f)*itauf;
}
/*
 * routines to be used in externally linkable potentials
 */
inline void init_timer(int    index,
		       double t0,
		       double tau)
{
  timer = index;
  t0f   = t0d   = t0;
  tauf  = taud  = tau;
  t1f   = t1d   = t0+tau;
  itauf = itaud = tau? 1./tau : 0.;
  if(index > 3) warning("init_timer(): timer index=%d is out of range [0,3]\n"
			" defaulting to index=0 [=adiabatic]\n",index);
}

inline double timer_double(double t)
{
  if (t <= t0d) return 0.;
  switch(timer) {
  case 3:  return linear_double(t);
  case 2:  return quasilinear_double(t);
  case 1:  return saturate_double(t);
  default: return adiabatic_double(t);
  }
}

inline float timer_float(float t)
{
  if (t <= t0f) return 0.f;
  switch(timer) {
  case 3:  return linear_float(t);
  case 2:  return quasilinear_float(t);
  case 1:  return saturate_float(t);
  default: return adiabatic_float(t);
  }
}
#endif

