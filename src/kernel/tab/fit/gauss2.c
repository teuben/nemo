/*
 * example fit function for 2 gaussians
 *
 *	f = a + b1*exp(-(x-c1)^2/(2*d1^2)) + b2*exp(-(x-c2)^2/(2*d2^2))
 *
 * 15-may-2004	created 
 */

#include <stdinc.h>

real func_gauss2(real *x, real *p, int np)
{
  real a,b,arg1,arg2;
  a = p[2]-x[0];
  b = p[3];
  arg1 = a*a/(2*b*b);

  a = p[5]-x[0];
  b = p[6];
  arg2 = a*a/(2*b*b);
  return p[0] + p[1]*exp(-arg1) + p[4]*exp(-arg2);

}
 
void derv_gauss2(real *x, real *p, real *e, int np)
{
  real a1,b1,arg1,a2,b2,arg2;
  a1 = p[2]-x[0];
  b1 = p[3];
  arg1 = a1*a1/(2*b1*b1);
  a2 = p[5]-x[0];
  b2 = p[6];
  arg2 = a2*a2/(2*b2*b2);
  e[0] = 1.0;
  e[1] = exp(-arg1);
  e[2] = -p[1]*e[1] * a1 / (b1*b1);
  e[3] = p[1] * e[1] * a1 * a1 / (b1*b1*b1);
  e[4] = exp(-arg2);
  e[5] = -p[4]*e[4] * a2 / (b2*b2);
  e[6] = p[4] * e[4] * a2 * a2 / (b2*b2*b2);

}

