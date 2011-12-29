/*
 * example fit function for N gaussians
 *
 *	f = a + b1*exp(-(x-c1)^2/(2*d1^2)) + b2*exp(-(x-c2)^2/(2*d2^2)) ...
 *
 * 28-dec-2011  created from gauss2
 * 29-dec-2011  added version where all dispersions the same
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



/* fixing all widths to be the same */


real func_gaussn_1(real *x, real *p, int np)
{
  real a,b,arg,val;
  int i,ng;

  ng = np - 2;
  if (ng%2) error("Number of parameters not  2 + 2N_g  (%d)",np);
  ng /= 2;

  if (debug_first) {
    dprintf(0,"gaussn_1[%d]: p1 + p3*exp(-(x-p4)^2/(2*p2^2)) + ...[p2 all same]\n",ng);
    debug_first = 0;
  }

  val = p[0];
  for (i=0; i<ng; i++) {
    a = p[2*i+3] - x[0];
    b = p[1];
    arg = a*a/(2*b*b);
    val += p[2*i+2] * exp(-arg);
  }
  return val;
}
 
void derv_gaussn_1(real *x, real *p, real *e, int np)
{
  real a,b,arg;
  int i,ng;

  /* sanity check */
  ng = np - 2;
  if (ng%2) error("Number of parameters not  2 + 2N_g  (%d)",np);
  ng /= 2;

  e[0] = 1.0;
  e[1] = 0.0;
  for (i=0; i<ng; i++) {
    a = p[2*i+3] - x[0];
    b = p[1];
    arg = a*a/(2*b*b);
    // @todo: better double check this, seems to work though
    /* amp, vel */
    e[2*i+2] = exp(-arg);
    e[2*i+3] = -p[2*i+2]*e[2*i+2] *  a  /  (b*b);    
    /* accum sig */
    e[1] += p[2*i+2] * e[2*i+2] *  (a*a) / (b*b*b);
  }
}

