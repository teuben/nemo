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

void vtc_get_force_direct(Forceinfo *fi, Nbodyinfo *nb)
{
    int i, k;
    double epsinv;

    assert(NJMAX > nb->n);
    (vtc_force_calculator[fi->calculator])(nb->n, nb->x, nb->n, nb->x, nb->m, fi->eps,
					   nb->a, nb->p);
    if (fi->eps != 0.0) { 
	epsinv = 1.0/fi->eps;
	for (i = 0; i < nb->n; i++) {
	    nb->p[i] = nb->p[i] + nb->m[i] * epsinv;
	}
    }
}
