/*
 * UTIL.C: various useful routines and functions.
 *
 *      17-feb-94 fixed bug with -DSINGLEPREC               pjt
 */

#include "code.h"

extern  double xrandom(double,double);

/*
 * PICKVEC: generate random coordinates within a unit sphere.
 */

void pickvec(
	     vector x,                               /* coord vector to generate */
	     bool cf)                                /* pick from 1/r^2 profile */
{

    if (debug)
        printf("pickvec: cf = %d\t", cf);
    if (cf)					/* cent. concentrated?      */
	pickshell(x, NDIM, xrandom(0.0, 1.0));	/*   pick from M(r) = r     */
    else
	pickball(x, NDIM, 1.0);		/*   use uniform distr.     */
    if (debug)
        printf("x = [%8.4f,%8.4f,%8.4f]\n", x[0], x[1], x[2]);
}
