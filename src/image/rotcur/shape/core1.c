#include <nemo.h>

real rotcur_core(real r, int n, real *p, real *d)
{
  real x = r / p[1];
  d[0] = x/(1+x);
  d[1] = -p[0]*d[0]/(p[1]*(1+x));
  return p[0] * d[0];
}

