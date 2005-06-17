/*
 * cp80.c:   Used by Contopoulos & Papayanopoulos (1980, AA 92, 33)
 *           in their study of orbits in barred galaxies
 *           (See also Barbanis & Woltjer, 1967)
 *
 *  V = V_0 + V_1
 *      where the axisymmetric part (V_0) is the Isochrone model
 *      and a cos(2*phi) perturbation V_1=eps.sqrt(r).(16-r)cos(2*phi)
 *
 *	7-mar-92   happy gcc2.0				pjt
 *	  oct-93   get_pattern
 *        feb-03   ANSI coding
 *        sep-04   double/float
 *
 */

/*CTEX
 *    {\bf potname=cp80
 *	 potpars={\it $\Omega,\epsilon$}}
 *
 *  Contopoulos \& Papayannopoulos (1980, A\&A, 92,33)
 *  used this potential
 *  in their study of orbits in barred galaxies, and has
 *  been used in other studies thereafter.
 *  Note that their
 *  ``bar'' is oriented along the Y-axis, an axis ratio is not
 *  well defined, and for larger values of $\epsilon$ the density
 *  can be negative. The potential used is given by adding an
 *  axisymmetric component to a m=2 fourier component:
 *  $$
 *     \Phi = \Phi_1 + \Phi_2
 *  $$
 *  where $\Phi_1$ is the Isochrone potential with unit scalelength and
 *  mass, and $\Phi_2$ the Barbanis \& Woltjer (1965) potential:
 *  $$
 *        \Phi_1 = - { 1 \over { (1 + \sqrt{1+r^2})}}
 *  $$
 *  and
 *  $$
 *         \Phi_2 = \epsilon  r (16-r) cos(2\phi)
 *  $$
 *
 *  A value of $\epsilon=0.00001$ is the default for a moderate bar,
 *  whereas 0.001 is a strong bar!
 */
 
#include <stdinc.h>
#include <potential_float.h>

local double omega = 0.05;         /* pattern speed */
local double eps = 0.00001;        /* 'bar' perturbation: 0.001 is strong! */

/*  Following two are not considered parameters in the model, but fixed */
/*  yet in the code they are treated as parameters                      */
local double iso_mass = 1.0;
local double iso_radius = 1.0;

local double iso_radius2;

void inipotential (int *npar, double *par, char *name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) eps = par[1];
    if (n>2) warning("CP80 potential: only 2 parameters used");

    iso_radius2 = sqr(iso_radius);

    dprintf (0,"INIPOTENTIAL CP80 [2D version]\n");
    dprintf (1,"  Parameters : Pattern Speed=%g  eps=%g\n",omega,eps);
    par[0] = omega;
}
    
void potential_double (int *ndim,double *pos,double *acc,double *pot,double *time)
{
    double a, tmp, rad2, rad, rads, ft, fr, cosp, sinp, costp, sintp;

    rad2 = sqr(pos[0]) + sqr(pos[1]);

/* Isochrone */
    a=sqrt(rad2+iso_radius2);
    *pot = -iso_mass / (iso_radius + a);
    tmp = *pot / (a * (iso_radius + a));
    acc[0] = tmp*pos[0];
    acc[1] = tmp*pos[1];
    if (*ndim>2) acc[2] = 0.0;
    if (rad2==0.0) return;

/* Barbanis & Woltjer's bar: */
    rad = sqrt(rad2);
    rads = sqrt(rad);
    costp=(sqr(pos[0])-sqr(pos[1]))/rad2;
    sintp=2*pos[0]*pos[1]/rad2;
    cosp=pos[0]/rad;
    sinp=pos[1]/rad;

    *pot += eps*rads*(16-rad)*costp;
    fr = eps*(16-3*rad)*costp/rads;	/* (neg of) radial and */
    ft = 2*eps*(rad-16)*sintp/rads;	/* tangential force */
    acc[0] -= fr*cosp - ft*sinp;
    acc[1] -= fr*sinp + ft*cosp;
}

