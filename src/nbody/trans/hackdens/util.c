/*
 * UTIL.C: various useful routines and functions.
 *
 *      18-jul-92  PJT  replaced many if(debug)printf(...) by dprintf(1,...)
 */

#include "defs.h"

/*
 * PICKVEC: generate random coordinates within a unit sphere.
 */

pickvec(x, cf)
vector x;                               /* coord vector to generate */
bool cf;                                /* pick from 1/r^2 profile */
{
    double xrandom(), xd[NDIM];

    dprintf(1,"pickvec: cf = %d\t", cf);
    if (cf)					/* cent. concentrated?      */
	pickshell(xd, NDIM, xrandom(0.0, 1.0));	/*   pick from M(r) = r     */
    else
	pickball(xd, NDIM, 1.0);		/*   use uniform distr.     */
    SETV(x, xd);				/* copy, for SINGLEPREC     */
    dprintf(1,"x = [%8.4f,%8.4f,%8.4f]\n", x[0], x[1], x[2]);
}
