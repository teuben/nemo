#include <nemo.h>

static int first = 1;

real rotcur_core123(real r, int n, real *p, real *d)
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

