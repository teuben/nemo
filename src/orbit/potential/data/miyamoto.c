/*
 * miyamoto.c: procedures for intializing and calculating the forces and
 *             potential of a Miyamoto-Nagai (generalized 
 *             Kuzmin/Plummer) potential
 *          Refs: BT pp. 43-44; Miyamoto and Nagai PASJ 27, 533 (1975)
 *          Potential Phi_m (R, z) = -GM / sqrt (R^2+(a+sqrt(z^2+b^2))^2)
 *             Parameters: a, b (shape parameters), M (mass); G=1
 *             Names used: miya_ascal, miya_bscal, miya_mass
 *          If a=0 standard Plummer model; if b=0 standard Kuzmin disk
 *             Realistic choice: b = 0.2 a
 *
 * For a realistic galaxy model, see e.g. http://arxiv.org/abs/1502.00627
 *
 *  March 90 Stefano Casertano, University of Pittsburgh
 * 10-Nov-90 inserted omega as first parameter for standard Nemo  PJT
 *  6-oct-91 fixed bug - and made code accept both XYZ and XZY versions (pjt)
 *  7-mar-92 merged sun and 3b1 versions once more			 pjt
 *    oct-93 get_pattern
 *
 */

/*CTEX
 *  {\bf potname=miyamoto
 *       potpars={\it $\Omega,a,b,M$}}
 *
 *  $$
 *    \Phi = -  {  M  \over
 *                    {   ....  }
 *		}
 * $$
 */

#include <stdinc.h>
#include <potential_float.h>

local double omega = 0.0;		/* pattern speed */
local double miya_ascal = 0.0;
local double miya_bscal = 1.0;
local double miya_mass = 1.0;

#if !defined(Y) && !defined(Z)
/*                                          default: XZY setup (re-oriented) */
#define X 0
#define Y 2
#define Z 1
#else
/*                      the normal axisymmetric case, as defined in e.g. B&T */
#define X 0
#define Y 1
#define Z 2
#endif

void inipotential (int *npar, double *par, string name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) miya_ascal = par[1];
    if (n>2) miya_bscal = par[2];
    if (n>3) miya_mass  = par[3];
    if (n>4) warning("Miyamoto: only first 4 parameters recognized");

    dprintf (1,"INIPOTENTIAL Miyamoto-Nagai: X=%d Y=%d Z=%d\n",X,Y,Z);
    dprintf (1,"  Parameters : omega, ascale, bscale, mass = ");
    dprintf (1,"  %f %f %f %f \n", omega, miya_ascal, miya_bscal, miya_mass);
    par[0] = omega;
}
    
void potential_double(int *ndim,double *pos,double *acc,double *pot,double *time)
{
    double qpar, spar, rcyl, tmp;
    int i;

    rcyl = sqrt (sqr(pos[X])+sqr(pos[Y]));
    qpar = sqrt (sqr(pos[Z])+sqr(miya_bscal));
    spar = sqrt (sqr(rcyl)+sqr(miya_ascal+qpar));

    *pot = - miya_mass / spar;
    tmp = *pot / sqr(spar);
    for (i=0; i<(*ndim); i++)
        acc[i] = tmp*pos[i];

    if (miya_ascal > 0.0) 
        acc[Z] *= (miya_ascal+qpar)/qpar;

}
