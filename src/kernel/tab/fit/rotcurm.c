/*
 * fit function for a rotcurm rotation curve
 *
 *	f = p0 * (x/p1) / (1+(x/p1)**p2)**(1/p2)
 *
 * 4-jul-2019	created
 * 9-jul-2019   args..... this is the same as the "core" function in rotcurshape
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

  r    = x[0]/p[1];
  arg1 = 1 + pow(r,p[2]);
  
  return p[0] * r / pow( arg1, 1/p[2]);
}
 
void derv_rotcurm(real *x, real *p, real *e, int np)
{
  real r, arg1;

  r    = x[0]/p[1];
  arg1 = 1+pow(r,p[2]);

  e[0] = r/pow( arg1, 1/p[2]);
  e[1] = p[0]*x[0]*x[0]*pow(x[0]/p[1],p[2]-1)*pow(arg1,-1-1/p[2])/(p[1]*p[1]*p[1]) -
         p[0]*x[0]/(p[1]*p[1]*pow(arg1,1/p[2]));
  e[2] = p[0]*r/pow(arg1,1/p[2]) * (log(arg1)/(p[2]*p[2]) - pow(r,p[2])*log(r)/(p[2]*arg1));
}


