/*
 * harmonic.c - harmonic (constant density) potential
 *
 *	7-mar-92  happy gcc2.0		pjt
 *	  oct-93  get_pattern		pjt
 * 	  sep-04  double/float		PJT
 *        aug-08  get rid of sqr        WD  (cvs committed jul-2010)
 */

/*CTEX
 * {\bf potname=harmonic
 *       potpars={\it $\Omega,\omega_x^2,\omega_z^2,\omega_z^2$}}
 *
 *
 * Harmonic potential
 * $$
 *        \Phi = {1 \over 2} \omega_x^2 x^2
 *             + {1 \over 2} \omega_y^2 y^2
 *             + {1 \over 2} \omega_z^2 z^2
 * $$
 */
 

#include <stdinc.h>
#include <potential_float.h>

static double omega = 0.0;           /* just put to zero until implemented */
static double h[3] = {1.0,1.0,1.0};  /* default parameters harmonic potential */

/* have our own square (to avoid linkage problems) */
inline double square(double x) { return x*x; }

void inipotential (int *npar, double *par, string name)
{
    int i;

    if (*npar>0) omega = par[0];
    for (i=1; i<(*npar); i++)
       h[i-1] = square(par[i]);      /* WD; original: h[i-1] = sqr(par[i]); */

    dprintf (1,"INI_POTENTIAL Harmonic\n");
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  wx^2,wy^2,wz^2= %f %f %f\n",h[0],h[1],h[2]);
    par[0] = omega;
}

void potential_double (int *ndim,double *pos,double *acc,double *pot,double *time)
{
        int    i;
        
        *pot = 0.0;
        for (i=0; i<*ndim; i++) {
                (*pot) += h[i] * square(pos[i]);
                acc[i] = -h[i] * pos[i];
        }
        *pot *= 0.5;
}
