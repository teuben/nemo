/*
 * gauss.c:  (spherical) gauss potential
 *
 *      19-mar-2021    spring break in covid time... while reading Cappelleri 2002
 *
 */

/*CTEX
 *	{\bf potname=gauss
 *       potpars={\it $\Omega,M,\sigma,p,q$}}
 *
 *  Gauss potential (e.g. Cappellari 2002 from the MGE formalism).  For the spherical case
 *  (cases with P and q < 1 will be considered another time)
 *
 * $$
 *    \rho = { M \over {(\sigma \sqrt{2\pi})}^3  }   e^{- {R^2}\over{2\sigma^2}}
 * $$
 *
 * $$
 *    \Phi = -  {  M  \over R } erf({ {R}\over{\sigma\sqrt{2}}})
 * $$
 */                     
  
#include <stdinc.h>
 
local double omega = 0.0;
local double g_mass = 1.0;
local double g_sigma = 1.0;
local double g_q = 1.0;         // discarded for now
local double g_p = 1.0;         // discarded for now

local double s2,s3;

void inipotential (int *npar, double *par, string name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) g_mass = par[1];
    if (n>2) g_sigma = par[2];
    if (n>3) warning("gauss: npar=%d only up to 3 parameters accepted",n);

    dprintf (1,"INIPOTENTIAL Gauss: \n");
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  mass, sigma = %f %f \n",g_mass,g_sigma);
	
    par[0] = omega;
    s2 = sqrt(2);
    s3 = sqrt(PI);
}


void potential_double (int *ndim,double *pos,double *acc,double *pot,double *time)
{
  double r2, r, tmp, x;

  r2 = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
  
  if (r2 == 0) {
    *pot = -g_mass / g_sigma * sqrt(2/PI);
    acc[0] = acc[1] = acc[2] = 0;
    return;
  }

  r = sqrt(r2);
  x = r/(s2*g_sigma);

  tmp = erf(x) / r;
    
  *pot = -g_mass * tmp;

  tmp = tmp - s2*exp(-x*x)/g_sigma/s3;
  tmp = -tmp * g_mass/r2;
  
  acc[0] = tmp*pos[0];
  acc[1] = tmp*pos[1];
  acc[2] = tmp*pos[2];
}

