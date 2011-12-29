/*
 * example fit function for N gaussians
 *
 *	f = a + b1*exp(-(x-c1)^2/(2*d1^2)) + b2*exp(-(x-c2)^2/(2*d2^2)) ...
 *
 * 28-dec-2011  created from gauss2
 */

#include <stdinc.h>

static int debug_first = 1;

real func_gaussn(real *x, real *p, int np)
{
  real a,b,arg,val;
  int i,ng;

  ng = np - 1;
  if (ng%3) error("Number of parameters not  1 + 3N_g  (%d)",np);
  ng /= 3;

  if (debug_first) {
    dprintf(0,"gaussn[%d]: p1 + p2*exp(-(x-p3)^2/(2*p4^2)) + ...\n",ng);
    debug_first = 0;
  }

  val = p[0];
  for (i=0; i<ng; i++) {
    a = p[3*i+2] - x[0];
    b = p[3*i+3];
    arg = a*a/(2*b*b);
    val += p[3*i+1] * exp(-arg);
  }
  return val;
}
 
void derv_gaussn(real *x, real *p, real *e, int np)
{
  real a,b,arg;
  int i,ng;

  /* sanity check */
  ng = np - 1;
  if (ng%3) error("Number of parameters not  1 + 3N_g  (%d)",np);
  ng /= 3;

  e[0] = 1.0;
  for (i=0; i<ng; i++) {
    a = p[3*i+2]-x[0];
    b = p[3*i+3];
    arg = a*a/(2*b*b);
    e[3*i+1] = exp(-arg);
    e[3*i+2] = -p[3*i+1]*e[3*i+1] *  a  /  (b*b);
    e[3*i+3] =  p[3*i+1]*e[3*i+1] * a*a / (b*b*b);
  }
}

