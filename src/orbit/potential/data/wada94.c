/*
 * wada94.c:  2D Toomre potential + bar perturbation
 *	      as done in Sanders (1977, ApJ 217, 916)
 *	      See also :
 *              Wada & Habe (1992, MN258, 82)
 *              Wada (1994, PASJ 46, 165)
 *
 *	16-jun-96	fixed bug in setting wada_e (n>2 -> n>3)	pjt
 *      19-sep-04	float/double
 */
/*CTEX
 *    {\bf potname=wada94
 *	 potpars={\it $\Omega,c,a,\epsilon$}}
 *
 *  Wada (1994, PASJ 46, 165) and also
 *  Wada \& Have (1992, MN 258, 82)
 *  used this potential
 *  in the study of gaseous orbits in barred galaxies.
 *  $$
 *     \Phi = \Phi_0 + \Phi_b
 *  $$
 *  where $\Phi_1$ is the Toomre potential with scalelength $a$ 
 *  $$
 *        \Phi_0 = - { 1 \over \sqrt{R^2 + a^2}}
 *  $$
 *  and
 *  $$
 *         \Phi_b = -\epsilon  {   {a R^2} \over { {(R^2 + a^2)}^2 } }
 *  $$
 *  A relationship for the axisymmetric component is
 *  $$
 *      -\sqrt(27/4)
 *  $$
 */

#include <stdinc.h>
#include <vectmath.h>
#include <potential_float.h>

local double omega = 0.0;
local double wada_c = 1.0;    /* mass scaling: c = v_max * a * (27/4)**(1/4) */
local double wada_a = 1.0;    /* length scaling */
local double wada_e = 0.0;    /* eps: small coefficient */

local double a2, c2a;

void inipotential (int *npar, double *par, string name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) wada_c = par[1];
    if (n>2) wada_a = par[2];
    if (n>3) wada_e = par[3];
    if (n>4) warning("wada94: npar=%d only 4 parameters accepted",n);

    dprintf(0,"INIPOTENTIAL Wada: [2d version]\n");
    dprintf(1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf(1,"  c, a, eps_0 = %f %f %f\n",wada_c, wada_a, wada_e);
	
    a2 = sqr(wada_a);
    c2a = sqr(wada_c) / wada_a;     /* c^2/a : units: G*M */
    par[0] = omega;
}

void potential_double(int *ndim, double *pos, double *acc, double *pot, double *time)
{
    double r2, w, ws, ws3, fr, ft, cos2t, eps;

    r2  = pos[0]*pos[0] + pos[1]*pos[1];
    w = 1.0/(a2 + r2);
    ws = sqrt(w);
    ws3 = w*ws;

    *pot = 1.0;
    fr = -1.0;
    ft = 0.0;
    if (r2>0 && wada_e>0) {
        eps =  wada_e * wada_a *  r2 * ws3;
        cos2t = 2*pos[0]*pos[0]/r2-1;
        *pot +=  eps * cos2t;
        fr += 2 * wada_e * wada_a * cos2t * (a2-r2) * ws3;
        ft = -(c2a*ws) * 4 * wada_e * wada_a * ws3 * pos[0] * pos[1] / r2;
    }
    *pot *= -c2a*ws;
    fr *= c2a*ws3;

    acc[0] = fr * pos[0] - ft * pos[1];
    acc[1] = fr * pos[1] + ft * pos[0];
    acc[2] = 0.0;

}
