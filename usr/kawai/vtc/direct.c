#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "vtc.h"
#include "vtclocal.h"

/* 'typical' params for direct-summation algorithm */
void
vtc_get_default_direct_params(Forceinfo *fi)
{
    fi->eps = 0.02;
    fi->calculator = HOST;
    fi->ninteraction = 0;
}

#define NIMAX (2200000)
#define NJMAX (200000) /* on MDM this need to be smaller than JMEMMAX. I don't know why... */

void vtc_get_force_direct(Forceinfo *fi, Nbodyinfo *nb)
{
    int i, k, nj, off;
    double epsinv;
    static double atmp[NIMAX][3], ptmp[NIMAX];

    assert(NIMAX > nb->n);

    if (fi->calculator != GRAPE_POTENTIALONLY && fi->calculator != HOST_POTENTIALONLY) {
	for (i = 0; i < nb->n; i++) {
	    for (k = 0; k < 3; k++) {
		nb->a[i][k] = 0.0;
	    }
	}
    }
    if (fi->calculator != GRAPE_FORCEONLY && fi->calculator != HOST_FORCEONLY) {
	for (i = 0; i < nb->n; i++) {
	    nb->p[i] = 0.0;
	}
    }
    off = 0;
    nj = NJMAX;
    if (fi->calculator == HOST) {
	nj = 100;
    }
    while (off < nb->n) {
	if (off + nj > nb->n) {
	    nj = nb->n - off;
	}
	fprintf(stderr, "off: %d n: %d nj: %d\n", off, nb->n, nj);

	(vtc_force_calculator[fi->calculator])(nb->n, nb->x, nj, nb->x+off, nb->m+off, fi->eps,
					       atmp, ptmp);
	if (fi->calculator != GRAPE_POTENTIALONLY && fi->calculator != HOST_POTENTIALONLY) {
	    for (i = 0; i < nb->n; i++) {
		for (k = 0; k < 3; k++) {
		    nb->a[i][k] += atmp[i][k];
		}
	    }
	}
	if (fi->calculator != GRAPE_FORCEONLY && fi->calculator != HOST_FORCEONLY) {
	    for (i = 0; i < nb->n; i++) {
		nb->p[i] += ptmp[i];
	    }
#if 0
	    /* no correction is necessary for MD2.
	       MD2 automatically eliminate self interaction. */
	    if (fi->eps != 0.0) { 
		epsinv = 1.0/fi->eps;
		for (i = 0; i < nb->n; i++) {
		    nb->p[i] += nb->m[i] * epsinv;
		}
	    }
#endif
	}
	off += nj;
    }
}


#if 0
    if (fi->eps != 0.0) { 
	epsinv = 1.0/fi->eps;
	for (i = 0; i < nb->n; i++) {
#if 1       /* no correction is necessary MD2 automatically eliminate self interaction */
	    nb->p[i] = nb->p[i];
#else       /* but on vpp5k this correction seems to be necessary. I don't know why... */
	    nb->p[i] = nb->p[i] + nb->m[i] * epsinv;
#endif
	}
    }
#endif
