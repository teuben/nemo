/*
 * bulge1:	homogeneous bulge with certain axis ratio
 *
 *        18-may-04    made a float/double, for Alan Peyaud    PJT
 *
 */

/*CTEX
 *  {\bf potname=bulge1
 *	 potpars={\it $\Omega,M,R,c/a$}}
 *
 *  homogeneous oblate bulge with mass $M$, radius $R$, and axis ratio $c/a$
 */
 
 
#include <stdinc.h>			/* standard Nemo include */

static double omega = 0.0;	     /* just put to zero until implemented */
static double bumass = 1.0;
static double rbul = 1.0;
static double axirat = 1.0;

static double ecc, rbul2;

#define G  1.0

void inipotential (int *npar, double *par, char *name)
{
    if (*npar>0)
        omega = par[0];
    if (*npar>1)
        bumass = par[1];
    if (*npar>2)
        rbul = par[2];
    if (*npar>3)
        axirat = par[3];
    if (*npar>4)
        warning("Skipped potential parameters beyond 4");

    dprintf (1,"INI_POTENTIAL Homogenous oblate bulge %s\n",name);
    dprintf (1,"  Parameters : Pattern Speed = %f (should be forced 0.0)\n",omega);
    dprintf (1,"  mass = %f   length = %f   axirat = %f\n",bumass, rbul, axirat);

    if (axirat < 0.0) axirat = 1.0;
    if (axirat > 1.0) axirat = 1.0/axirat;
    if (axirat == 1.0)
        ecc = 0.0;
    else 
        ecc = sqrt(1.0 - sqr(axirat));
    rbul2 = sqr(rbul);

    dprintf(2,"  ecc = %f\n",ecc);
}


void potential_double (int *ndim, double *pos,double *acc,double *pot,double *time)
{
    int    i;
    double r,x,f;
        
    for (i=0, r=0; i<*ndim; i++)    /* find radius */
        r += sqr(pos[i]);
    dprintf(3,"r = %f",sqrt(r));
    if (ecc>0.0) {
        if (*ndim>2 && pos[2] != 0.0) 
    	    warning("bulge1_double: doesn't work outside xy-plane");
        r = sqrt(r);
        x = ecc / MAX(1.0, r/rbul);
        f = r * 3*G*bumass * (asin(x)-x*sqrt(1-x*x)) /
                    (2*rbul*rbul*rbul*ecc*ecc*ecc);
    } else {
        if (r>rbul2) {                 /* outside bulge */
            f = G * bumass / r;        /* Good old Newton */
            r = sqrt(r);
        } else {                       /* inside bulge */
            r = sqrt(r);
            f = G * bumass / (rbul*rbul*rbul) * r;
        }
    }
    if (r>0.0) {
        for (i=0; i<*ndim; i++) 
            acc[i] = -(pos[i]/r) * f;        /* radial force to cartesian */
    } else 
        for (i=0; i<*ndim; i++) 
            acc[i] = 0.0;
    *pot = 0.0;
    dprintf(3," f = %f acc=%f %f %f\n",f,acc[0],acc[1],acc[2]);
    
}
void potential_float (int *ndim, float *pos, float *acc, float *pot, float *time)
{
    int    i;
    float  r,x,f;
        
    for (i=0, r=0; i<*ndim; i++)    /* find radius */
        r += sqr(pos[i]);
    dprintf(3,"r = %f",sqrt(r));
    if (ecc>0.0) {
        if (*ndim>2 && pos[2] != 0.0) 
    	    warning("bulge1_float   doesn't work outside xy-plane");
        r = sqrt(r);
        x = ecc / MAX(1.0, r/rbul);
        f = r * 3*G*bumass * (asin(x)-x*sqrt(1-x*x)) /
                    (2*rbul*rbul*rbul*ecc*ecc*ecc);
    } else {
        if (r>rbul2) {                 /* outside bulge */
            f = G * bumass / r;        /* Good old Newton */
            r = sqrt(r);
        } else {                       /* inside bulge */
            r = sqrt(r);
            f = G * bumass / (rbul*rbul*rbul) * r;
        }
    }
    if (r>0.0) {
        for (i=0; i<*ndim; i++) 
            acc[i] = -(pos[i]/r) * f;        /* radial force to cartesian */
    } else 
        for (i=0; i<*ndim; i++) 
            acc[i] = 0.0;
    *pot = 0.0;
    dprintf(3," f = %f acc=%f %f %f\n",f,acc[0],acc[1],acc[2]);
}
