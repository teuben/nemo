/* pickshell, pickball, pickbox */
/*
 * PICKPNT.C: pick points at random in various distributions.
 *	 2-jun-88	some original version	JEB
 *	25-feb-92	happy gcc2.0
 *	24-mar-94	ansi
 *      22-jun-01       added ZENO compatibility
 */

#include <stdinc.h>

extern double xrandom(double, double);

void pickshell(			/* pick point on shell */
    real vec[],			/*   coordinates chosen */
    int ndim,			/*   number of dimensions */
    double rad)			/*   radius of shell */
{
    double rsq, rscale;
    int i;

    do {
	rsq = 0.0;
	for (i = 0; i < ndim; i++) {
	    vec[i] = xrandom(-1.0, 1.0);
	    rsq = rsq + vec[i] * vec[i];
	}
    } while (rsq > 1.0);
    rscale = rad / sqrt(rsq);
    for (i = 0; i < ndim; i++)
	vec[i] = vec[i] * rscale;
}

void pickball(			/* pick point within ball */
    real vec[],			/*   coordinates chosen */
    int ndim,			/*   number of dimensions */
    double rad)			/*   radius of ball */
{
    double rsq;
    int i;

    do {
	rsq = 0.0;
	for (i = 0; i < ndim; i++) {
	    vec[i] = xrandom(-1.0, 1.0);
	    rsq = rsq + vec[i] * vec[i];
	}
    } while (rsq > 1.0);
    for (i = 0; i < ndim; i++)
	vec[i] = vec[i] * rad;
}

void pickbox(			/* pick point within box */
    real vec[],			/*   coordinates chosen */
    int ndim,			/*   number of dimensions */
    double size)			/*   from -size to size */
{
    int i;

    for (i = 0; i < ndim; i++)
	vec[i] = xrandom(-size, size);
}

/* end of: pickpnt.c */
