/*
 *  Hernquist:		approximate R^(1/4) law which is a bit easier
 *			to write down in potential/density form:
 *
 *
 *	7-mar-92	happy gcc2.0 - optimized acc[] a bit	    pjt
 *	  oct-93        get_pattern				    pjt
 *	  feb-94	pretty bad sign error 			    pjt
 *	  jun-97	documented dehnen
 *	may-2002	provide both a _double and _float 
 *                      version for the new potproc interface       wd
 *      sep-2004        replaced call to sqr(A) with A*A
 *                      sqr() is bullshit and should never be used! wd
 */

/*CTEX
 *  {\bf potname=hernquist
 *       potpars={\it $\Omega,M,r_c$}}
 *
 * The Hernquist potential (ApJ, 356, pp.359, 1990) is a special $\gamma=1$ case
 * of the Dehnen potential. The potential is given by:
 * $$
 *     \Phi = - { M \over  {(r_c+r)}}
 * $$
 * and mass
 * $$      
 *         M(r) = M { r^2 \over {(r+r_c)}^2 }
 * $$
 * and density
 * $$
 *       \rho = { M \over {2\pi}} {r_c \over r} { 1 \over {(r+r_c)}^3}
 * $$
 */

 
#include <stdinc.h>                     /* standard Nemo include */

local real omega = 0.0;         /* pattern speed */
local real hmass = 1.0;		/* total mass */
local real a = 1.0;		/* scale lenght */

void inipotential (int *npar, double *par, string name)
{
    if (*npar>0) omega = par[0];
    if (*npar>1) hmass = par[1];
    if (*npar>2) a = par[2];
    if (*npar>3) warning("Skipped potential parameters beyond 3");

    dprintf (1,"INI_POTENTIAL Hernquist potential\n");
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  mass = %f   scalelength = %f\n",hmass, a);
    par[0] = omega;
}

#define POT								\
{									\
    int    i;								\
    real   r2,r,f;							\
									\
    for (i=0, r2=0; i<*ndim; i++)            /* radius - squared */	\
        r2 += pos[i]*pos[i];						\
									\
    if (r2==0.0) {							\
        *pot = -hmass/a;						\
        for (i=0; i<*ndim; i++)						\
            acc[i] = 0.0;						\
    } else {								\
        r = sqrt(r2);                        /* radius */		\
        f = 1.0/(r+a);                       /* temporary storage */	\
        *pot = -hmass * f;                   /* returned potential */	\
        f = (*pot) * f / r;                  /* radial force / r  */	\
        for (i=0; i<*ndim; i++)						\
            acc[i] = pos[i] * f;             /* make cartesian forces*/ \
    }									\
}

void potential_double (int *ndim,
		       double *pos,
		       double *acc,
		       double *pot,
		       double *time) POT
void potential_float  (int *ndim,
		       float *pos,
		       float *acc,
		       float *pot,
		       float *time) POT
#undef POT
