/*
 * plummer2.c:  2D a-spherical plummer potential 
 *              (see plummer.c for a full description)
 */
  
#include <stdinc.h>
#include <potential_float.h>
 
local double omega = 0.0;
local double plummer_mass = 1.0;
local double plummer_radius = 1.0;
local double q = 1.0;

local double r2, iq2[3];      /* scratch variables */

void inipotential (int *npar, double *par, string name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) plummer_mass = par[1];
    if (n>2) plummer_radius = par[2];
    if (n>3) q = par[3];
    if (n>4) warning("plummer2: npar=%d only 3 parameters accepted",n);

    dprintf (1,"INIPOTENTIAL Plummer2:\n");
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  mass, radius, q = %f %f %f\n",plummer_mass,plummer_radius,q);
	
    r2 = sqr(plummer_radius);
    iq2[0] = 1.0;
    iq2[1] = 1.0/sqr(q);
    iq2[2] = 1.0/sqr(q);
    par[0] = omega;
}

void potential_double(int *ndim, double *pos, double *acc, double *pot, double *time)
{
  double tmp, rad;

  rad = sqr(pos[0])*iq2[0] + sqr(pos[1])*iq2[1];
  tmp = 1.0/(rad+r2);

  *pot = -sqrt(tmp);

  if (rad > 0) {
    tmp *= (*pot) * plummer_mass;
  } else
    tmp = 0;
  *pot *= plummer_mass;
  acc[0] = tmp*pos[0]*iq2[0];
  acc[1] = tmp*pos[1]*iq2[1];
  acc[2] = 0.0;
}
