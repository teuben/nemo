/*
 * point.c:  point mass potential, with optional pseudo-newtonian term(s)
 *
 *
 */

/*CTEX
 *	{\bf potname=point
 *       potpars={\it $\Omega,M,e$}}
 *
 *  Point mass potential with optional Pseudo-Newtonian term(s)
 *
 * $$
 *    \Phi = -  {  M  \over
 *                    {   {r - \eps}  }
 * $$
 */                     
  
#include <stdinc.h>
#include <vectmath.h>	/* define DIMensionality */
 
local double omega = 0.0;
local double mass = 1.0;
local double eps = 0.0;	    /* radius=0 gives point mass */


void inipotential (int *npar, double *par, string name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) mass = par[1];
    if (n>2) eps = par[2];
    if (n>3) warning("point: npar=%d only 3 parameters accepted",n);

    dprintf (1,"INIPOTENTIAL Point: [3d version]\n");
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  mass, eps = %f %f \n",mass,eps);
	
    par[0] = omega;
}


void potential_double (int *ndim,double *pos,double *acc,double *pot,double *time)
{
  double tmp, r;

  r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
  if (r<=eps) error("NP violation");
  tmp = mass/(r-eps);
  *pot = -tmp;
  tmp = -tmp/((r-eps)*r);
  acc[0] = tmp*pos[0];
  acc[1] = tmp*pos[1];
  acc[2] = tmp*pos[2];

}

void potential_float (int *ndim,float *pos,float *acc,float *pot,float *time)
{
  float tmp, r;

  r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
  if (r<=eps) error("NP violation");
  tmp = mass/(r-eps);
  *pot = -tmp;
  tmp = -tmp/((r-eps)*r);
  acc[0] = tmp*pos[0];
  acc[1] = tmp*pos[1];
  acc[2] = tmp*pos[2];

}

