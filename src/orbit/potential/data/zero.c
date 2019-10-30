/*
 * zero.c:  a true zero potential, no forces anywhere
 *
 */

/*CTEX
 *	{\bf potname=zero}
 *
 *  Zero potential, no forces. Most potentials have a scaling factor that can
 *  be set to 0 which results in the same, but at the cost of some useless
 *  operations.
 * 
 *
 * $$
 *    \Phi =  0
 * $$
 */                     
  
#include <stdinc.h>
 
void inipotential (int *npar, double *par, string name)
{
  int n;

  n = *npar;
  if (n>0) warning("zero: npar=%d no parameters accepted",n);
}

void potential_double (int *ndim,double *pos,double *acc,double *pot,double *time)
{
  *pot = acc[0] = acc[1] = acc[2] = 0.0;
}
void potential_float (int *ndim,float *pos,float *acc,float *pot,float *time)
{
  *pot = acc[0] = acc[1] = acc[2] = 0.0;
}

