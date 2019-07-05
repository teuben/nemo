/*
 * example fit function for a rotcurm rotation curve
 *
 *	f = p0 * (x/p1) / (1+(x/p1)**m)**(1/m)
 *
 * 4-jul-2019	created 
 */

#include <stdinc.h>

static int debug_first = 1;

real func_rotcurm(real *x, real *p, int np)
{
  real r, arg1;

  if (debug_first) {
    dprintf(0,"rotcurm(%d): p0 * (x/p1) / (1+(x/p1)**p2)**(1/p2)\n",np);
    debug_first = 0;
  }

  r = x[0]/p[1];
  arg1 = 1+pow(r,p[2]);
  
  return p[0] * r / pow( arg1, 1/p[2]);
}
 
void derv_rotcurm(real *x, real *p, real *e, int np)
{
  real r, arg1;

  r = x[0]/p[1];
  arg1 = 1+pow(r,p[2]);

  e[0] = r/pow( arg1, 1/p[2]);
  e[1] = r*(pow(r,p[2])/arg1 - 1)/pow(arg1,1/p[2]);
  e[2] = p[0]*r*log(arg1)*log(r)*pow(r,p[2])/(p[2]*p[2]*pow(arg1,1/p[2]));
}

