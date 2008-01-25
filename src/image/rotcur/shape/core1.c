/*
 * naming conventions:      rotcur_NAME where NAME is the name
 *
 * input:    r    radius
 *           n    number of parameters
 *           p    p[n], parameters to the rotation curve
 * returns:  the rotation curve
 *           d    d[n], partial derivatives dV/dp
 *
 */

#include <nemo.h>

static int first = 1;

real rotcur_core1(real r, int n, real *p, real *d)
{
  real x = r / p[1];
  d[0] = x/(1+x);
  d[1] = -p[0]*d[0]/(p[1]*(1+x));
#if 1
  if (first) {
    dprintf(0,"First time in rotcur_core\n");
    first = 0;
  }
#endif
  return p[0] * d[0];
}

