/*
 * ZEROCMS.C: routines to find and zero the center of mass.
 *
 *	?? 	  ??	 Created
 *	7-mar-92  pjt    changed all double's to real's - happy gcc2.0 
 *	5-mar-95  pjt    ansi prototypes
 */

#include <stdinc.h>

#define MDIM 12

void findcms(real *cms, real *space, int ndim, real *mass, int npnt)
{
    real mtot, *sp, *mp;
    int i, j;

    mtot = 0.0;
    for (j = 0; j < ndim; j++)
	cms[j] = 0.0;
    for (sp = space, mp = mass, i = 0; i < npnt; i++, mp++) {
	mtot = mtot + *mp;
	for (j = 0; j < ndim; j++, sp++)
	    cms[j] += (*mp) * (*sp);
    }
    for (j = 0; j < ndim; j++)
	cms[j] /= mtot;
}

void zerocms(real *space, int ndim, real *mass, int npnt, int nzer)
{
    real cms[MDIM], *sp;
    int i, j;

    if (ndim > MDIM) 
        error("zerocms: too many dimensions; ndim=%d MDIM=%d, ndim, MDIM");
    findcms(cms, space, ndim, mass, nzer);
    for (sp = space, i = 0; i < npnt; i++)
	for (j = 0; j < ndim; j++, sp++)
	    *sp = *sp - cms[j];
}

