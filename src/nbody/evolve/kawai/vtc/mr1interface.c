/*
 * interface to MR1 functions
 * necessary to link libvtc.a to AMBER_MDGRAPE2
 */

#include <stdio.h>
#include <stdlib.h>
#include "vtc.h"
#include "vtclocal.h"

static Forceinfo mr1fi;
static Nbodyinfo mr1nb;
static int params_initialized = FALSE;

#ifdef __linux__
#define FNAME(x) (x ## __)
#else
#define FNAME(x) (x ## _)
#endif

static void
mr1calccoulomb_init(int n, double eps, double theta, int ncrit,
		    int node_div_crit, int me_order,
		    int test_id)
{
    int k;

    mr1nb.n = n;
    mr1nb.m = (double *)malloc(sizeof(double)*n);
    mr1nb.x = (double (*)[3])malloc(sizeof(double)*3*n);
    mr1nb.v = (double (*)[3])malloc(sizeof(double)*3*n);
    mr1nb.a = (double (*)[3])malloc(sizeof(double)*3*n);
    mr1nb.p = (double *)malloc(sizeof(double)*n);
    if (NULL == mr1nb.p) {
	perror("mr1calccoulomb_init");
	exit(1);
    }
    if (!params_initialized) {
	params_initialized = TRUE;
	vtc_get_default_tree_params(&mr1fi);
	/* modify some of them */
	mr1fi.node_div_crit = node_div_crit;
	mr1fi.eps = eps;
	mr1fi.theta = theta;
	mr1fi.ncrit = ncrit;
	mr1fi.test_id = test_id;
	mr1fi.p = me_order;
	mr1fi.negativemass = TRUE;
	mr1fi.calculator = GRAPE_FORCEONLY;
    }
}

void
mr1calccoulomb_set_tree_param(double eps, double theta, int ncrit,
			      int node_div_crit, int me_order,
			      int test_id)
{
    if (!params_initialized) {
	params_initialized = TRUE;
	vtc_get_default_tree_params(&mr1fi);
    }
    mr1fi.node_div_crit = node_div_crit;
    mr1fi.eps = eps;
    mr1fi.theta = theta;
    mr1fi.ncrit = ncrit;
    mr1fi.p = me_order;
    mr1fi.negativemass = TRUE;
    mr1fi.test_id = test_id;
}

/* currently n is assumed to be a constant, periodicflag==0, and natchangeflag ==0.
   rscale is ignored. */
void
MR1calccoulomb_tree(double *x, int n, double *chg,
		    double rscale, int tblno, double xmax, int periodicflag, int natchangeflag,
		    double *force)
{
    int i, k;
    static int firstcall = 1;
    static int cnt = 0;

    if (firstcall) {
	firstcall = 0;
	mr1calccoulomb_init(n, 0.0, 0.7, 4000, 8, 1, -1);
	/* # of particle
	 * softening parameter
	 * opening parameter
	 * vectorization parameter
	 * test id (-1 for normal operation)
	 * node division parameter
	 */
    }
    for (i = 0; i < n; i++) {
	for (k = 0; k < 3; k++) {
	    mr1nb.x[i][k] = x[i*3+k];
	}
	mr1nb.m[i] = chg[i];
    }
    switch (tblno) {
    case 0:
	mr1fi.calculator = GRAPE_FORCEONLY;
	break;
    case 1:
	mr1fi.calculator = GRAPE_POTENTIALONLY;
	break;
    case 2:
	mr1fi.calculator = HOST_FORCEONLY;
	break;
    case 3:
	mr1fi.calculator = HOST_POTENTIALONLY;
	break;
    default:
	fprintf(stderr, "MR1calccoulomb_tree: unknown force mode\n");
	exit(1);
    }
    Cfprintf(stderr, "p: %d cnt: %d\n", mr1fi.p, cnt);
    vtc_get_force_tree(&mr1fi, &mr1nb);
    vtc_close_grape();
    switch (tblno) {
    case 0:
    case 2:
	for (i = 0; i < n; i++) {
	    for (k = 0; k < 3; k++) {
		force[i*3+k] = mr1nb.a[i][k];
	    }
	}
	break;
    case 1:
    case 3:
	for (i = 0; i < n; i++) {
	    force[i*3+0] = -mr1nb.p[i];
	}
	break;
    default:
	fprintf(stderr, "MR1calccoulomb_tree: unknown force mode\n");
	exit(1);
    }
    cnt++;
}

void
FNAME(mr1calccoulomb_tree)(double *x, int *n, double *chg,
		     double *rscale, int *tblno, double *xmax,
		     int *periodicflag, int *natchangeflag,
		     double *force)
{
    MR1calccoulomb_tree(x, *n, chg, *rscale, *tblno, *xmax,
			*periodicflag, *natchangeflag, force);
}

void
FNAME(mr1calccoulomb_set_tree_param)(double *eps, double *theta, int *ncrit,
			       int *node_div_crit, int *me_order,
			       int *test_id)
{
    mr1calccoulomb_set_tree_param(*eps, *theta, *ncrit,
				  *node_div_crit, *me_order,
				  *test_id);
}
