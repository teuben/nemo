/*
 *  Jaffe:		Jaffe model
 *		jun-97	documented Dehnen
 *
 */

/*CTEX
 *	{\bf potname=jaffe
 *       potpars={\it $\Omega,M,r_c$}}
 *
 * The Jaffe potential (BT, pp.237, see also MNRAS 202, 995 (1983))),
 * is another special $\gamma=2$ case of the Dehnen potential.
 *
 * $$
 *    \Phi = - {  M \over r_c}  \ln{
 *                         \left( { r \over {r_c + r} } \right) }
 * $$
 */
 
 
#include <stdinc.h>                     /* standard Nemo include */

local double omega = 0.0;           /* just put to zero until implemented */
local double hmass = 1.0;	/* total mass */
local double a = 1.0;		/* scale lenght */

local double vc;

void inipotential (int *npar, double *par, string name)
{
    if (*npar>0)   omega = par[0];
    if (*npar>1)   hmass = par[1];
    if (*npar>2)   a = par[2];
    if (*npar>3)   warning("Skipped potential parameters beyond 3");

    dprintf (1,"INI_POTENTIAL Jaffe potential\n");
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  mass = %f   scalelength = %f\n",hmass, a);

    vc = hmass/a;       /* velocity near center */

    par[0] = omega;

}

#define POT									\
{										\
    int    i;									\
    double r2,r,f;								\
        									\
    for (i=0, r2=0; i<*ndim; i++)              /* radius - squared */		\
        r2 += sqr(pos[i]);							\
										\
    if (r2==0.0) {								\
        *pot = 0.0;		       	       /* a lie though */		\
        for (i=0; i<*ndim; i++)							\
            acc[i] = 0.0;							\
    } else {									\
        r = sqrt(r2);                          /* radius */			\
        f = 1.0/(r+a);                         /* temporary storage */		\
        *pot = vc * log(r*f);                  /* returned potential */		\
        f *= hmass/r2;								\
        for (i=0; i<*ndim; i++)							\
            acc[i] = -pos[i]*f;                /* radial force to cartesian */	\
    }										\
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

