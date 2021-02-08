/****************************************************************************/
/* PICK.C: utility routines for various sorts of math operations. Most   */
/* these functions work with real values, meaning that they can handle      */
/* either floats or doubles, depending on compiler switches.                */
/* Copyright (c) 1999 by Joshua E. Barnes, Tokyo, JAPAN.                    */
/*    revised for SPH calculation by Jin Koda, Tokyo, JAPAN. 2000           */
/*    adapted for NEMO for Koda's SPH and Josh' treecode1                   */
/****************************************************************************/

// see src/nbody/cores/pickpnt.c for original

#include "stdinc.h"
#include "mathfns.h"

/*
 * PICKSHELL: pick point on shell.
 */

void pickshell1(real *vec, int ndim, double rad)
{
    real rsq, rscale;
    int i;

    do {
        rsq = 0.0;
        for (i = 0; i < ndim; i++) {
            vec[i] = xrandom(-1.0, 1.0);
            rsq = rsq + vec[i] * vec[i];
        }
    } while (rsq > 1.0);
    rscale = rad / rsqrt(rsq);
    for (i = 0; i < ndim; i++)
        vec[i] = vec[i] * rscale;
}

/*
 * PICKBALL: pick point within ball.
 */

void pickball1(real *vec, int ndim, double rad)
{
    real rsq;
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

/*
 * PICKBOX: pick point within box.
 */

void pickbox1(real *vec, int ndim, double size)
{
    int i;

    for (i = 0; i < ndim; i++)
        vec[i] = xrandom(- size, size);
}
