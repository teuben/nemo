/*
 * zero.c:  a true zero potential
 *
 */

/*CTEX
 *	{\bf potname=zero}
 *
 *  Zero potential
 *
 * $$
 *    \Phi =  0
 * $$
 */                     
  
#include <stdinc.h>
#include <vectmath.h>	/* define DIMensionality */
 
void inipotential (int *npar, double *par, string name)
{
    int n;

    n = *npar;
    if (n>0) warning("zero: npar=%d no parameters accepted",n);
}

void potential (int *ndim,double *pos,double *acc,double *pot,double *time)
{
    *pot = 0.0;
    acc[0] = 0.0;
    acc[1] = 0.0;
    acc[2] = 0.0;
}
