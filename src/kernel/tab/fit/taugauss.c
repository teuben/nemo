/*
 * example fit function optical depth gaussian function
 *
 *	tau = b*exp(-(x-c)^2/(2*d^2))
 *      f = a*(1 - exp(-tau))
 *
 * 3-jan-2012	created 
 */

#include <stdinc.h>

static int debug_first = 1;

real func_taugauss(real *x, real *p, int np)
{
  real a,b,arg,tau;

  if (debug_first) {
    dprintf(0,"taugauss: p1*(1-exp(-p1*exp(-(x-p3)^2/(2*p4^2))))\n");
    debug_first = 0;
  }

  a = p[2]-x[0];
  b = p[3];
  arg = a*a/(2*b*b);
  tau =  p[1] * exp(-arg);

  return p[0]*(1-exp(-tau));
}

void derv_taugauss(real *x, real *p, real *e, int np)
{
  real a,b,arg,tau,tau1;
  a = x[0]-p[2];
  b = p[3];
  arg = a*a/(2*b*b);
  tau = p[1] * exp(-arg);
  tau1 = 1-exp(-tau);
  e[0] = tau1;
  e[1] = tau1*exp(-arg);
  e[2] = tau1*p[1]*e[1] *  a  /  (b*b);
  e[3] = tau1*p[1]*e[1] * a*a / (b*b*b);
}

