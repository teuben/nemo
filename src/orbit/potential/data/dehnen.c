/*
 *	3-nov-93	created				                   pjt
 *     26-jun-96        finalized, special hernq and jaffe models builtin  PJT
 *     17-may-02        added potential_double, potential_float            WD 
 *                          
 */

/*CTEX
 *  {\bf potname=dehnen
 *       potpars={\it $\Omega,M,a,\gamma$}}
 *
 *  Walter Dehnen (1993, MN {\bf 265}, 250-256) introduced a
 *  family of potential-density pairs for spherical systems,
 *  which he often refers to as the ``gamma'' models.
 *
 * The potential is given by:
 * $$
 *	\Phi = { G M \over a } {1\over{2-\gamma}}
 *		 {\left[ 1 - {\left(r\over{r+a}\right)}^{2-\gamma}\right]}
 * $$
 * cumulative mass by
 * $$      
 *         M(r) = M { r \over {(r+a)}^{3-\gamma} }
 * $$
 * and density by
 * $$
 *       \rho = { {(3-\gamma)M} \over {4\pi}}
 *		{ a \over {r^{\gamma} (r+a)^{4-\gamma}}}
 * $$
 * with $0 <= \gamma < 3$.
 * Special cases are the Hernquist potential ($\gamma=1$), and the
 * Jaffe model ($\gamma=2$). The model with $\gamma=3/2$ seems to
 * give the best comparison withe de Vaucouleurs $R^{1/4}$ law.
 *
 * See also Tremaine et al. (1994, AJ, 107, 634) in which they describe
 * the same density models with $\eta=3-\gamma$ and call them
 * $\eta$-models.
 */

 
#include <stdinc.h>                     /* standard Nemo include */

local real omega = 0.0;
local real m = 1.0;		/* total mass */
local real a = 1.0;		/* scale lenght */
local real gam = 1.5;		/* exponent */

local real dmass, vc;
local bool Qjaffe, Qhernq;


void inipotential (int *npar, double *par, string name)
{
    if (*npar>0) omega = par[0];
    if (*npar>1) m = par[1];
    if (*npar>2) a = par[2];
    if (*npar>3) gam = par[3];
    if (*npar>4) warning("Skipped potential parameters beyond 3");

    dprintf (1,"INI_POTENTIAL Dehnen potential\n");
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  mass = %f scalelength = %f gamma = %f\n",m, a, gam);
    par[0] = omega;

    dmass = m / (2-gam);
    vc = m / a;

    Qjaffe = (gam == 2.0);
    Qhernq = (gam == 1.0);

    if (Qjaffe) warning("Dehnen: gamma=2 Jaffe model");
    if (Qhernq) warning("Dehnen: gamma=1 Hernquist model");
}

#define DEHNEN_POT								\
{										\
    int    i;									\
    real   r2,r,f,g;								\
										\
    for (i=0, r2=0; i<*ndim; i++)            /* radius - squared */		\
        r2 += sqr(pos[i]);							\
										\
    if (Qjaffe) {								\
        if (r2==0.0) {								\
            warning("dehnen/jaffe r=0");					\
            *pot = 0.0;								\
            for (i=0; i<*ndim; i++) acc[i] = 0.0;				\
            return;								\
        }									\
        r = sqrt(r2);                        /* radius */			\
        f = 1.0/(r+a);                       /* temporary storage */		\
        *pot = vc * log(r*f);                /* returned potential */		\
        f *= m/r2;                           /* radial_force / r */		\
        for (i=0; i<*ndim; i++)							\
            acc[i] = -pos[i]*f;              /* radial force to cartesian */	\
        return;        								\
    } else if (Qhernq) {							\
        if (r2==0.0) {								\
            *pot = -m / a;							\
            for (i=0; i<*ndim; i++) acc[i] = 0.0;				\
            return;								\
        }									\
        r = sqrt(r2);                        /* radius */			\
        f = 1.0/(r+a);                       /* temporary storage */		\
        *pot = -m * f;                       /* returned potential */		\
        f = (*pot) * f / r;                  /* radial force / r  */		\
        for (i=0; i<*ndim; i++)							\
            acc[i] = pos[i] * f;             /* make cartesian forces */	\
    } else {									\
        if (r2==0.0) {								\
            *pot = -dmass/a;							\
            for (i=0; i<*ndim; i++) acc[i] = 0.0;				\
            return;								\
        }									\
        r = sqrt(r2);                        /* radius */			\
        f = r/(r+a);                         /* temporary storage */		\
        g = pow(f,-gam);							\
        *pot = -dmass * (1-f*f*g) / a;						\
	f = -m/qbe(a+r)*g;	             /* radial_force/r	*/		\
        for (i=0; i<*ndim; i++)							\
            acc[i] = pos[i] * f;             /* make cartesian forces */       	\
    }										\
}

void potential_double(int   *ndim,            	    /* I: # dims, ignored     */
		      double*pos,             	    /* I: (x,y,z)             */
		      double*acc,            	    /* O: (ax,ay,az)          */
		      double*pot,            	    /* O: potential           */
		      double*time) DEHNEN_POT       /* I: time                */

void potential_float (int   *ndim,            	    /* I: # dims, ignored     */
		      float *pos,             	    /* I: (x,y,z)             */
		      float *acc,            	    /* O: (ax,ay,az)          */
		      float *pot,            	    /* O: potential           */
		      float *time) DEHNEN_POT       /* I: time                */

#undef DEHNEN_POT
