/*
 * hh.c - Henon Heiles 1964 potential   ApJ 69, 73-79. - 
 *
 *       8-dec-2019  finally
 */

/*CTEX
 * {\bf potname=htt
 *       potpars={\it $\Omega,\lambda$}}
 *
 *
 * Henon Heiles 1964 potential
 * $$
 *        \Phi = {1 \over 2} ( x^2 + x^2 ) + \lambda ( x^2 y - {1\over 3} y^3 )
 * $$
 */
 

#include <stdinc.h>
#include <potential_float.h>

static double omega = 0.0;     
static double lambda = 1.0;


void inipotential (int *npar, double *par, string name)
{
  if (*npar>0) omega  = par[0];
  if (*npar>1) lambda = par[1];
  
  dprintf (1,"INI_POTENTIAL Henon Heiles 1964\n");
  dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
  dprintf (1,"               lambda = %f\n",lambda);
  par[0] = omega;
}

void potential_double (int *ndim,double *pos,double *acc,double *pot,double *time)
{
  int    i;
  double x = pos[0];
  double y = pos[1];
  double xx = x*x;
  double yy = y*y;
  
  *pot = 0.5 * (xx + yy) + lambda * (xx - yy/3) * y;
  acc[0] = -x - lambda*2*x*y;
  acc[1] = -y - lambda*(xx - yy);
  acc[2] = 0.0;
}
