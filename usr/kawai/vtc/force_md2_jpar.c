#define _FORCE_C_ (1)
#include <stdio.h>
#include <math.h>
#include "vtc.h"
#include "vtclocal.h"

static void get_force_none(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
			      double eps, double (*a)[3], double *p);
static void get_force_host(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
			   double eps, double (*a)[3], double *p);
static void get_force_grape(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
			    double eps, double (*a)[3], double *p);
static void get_potential_grape(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
				double eps, double (*a)[3], double *p);
static void get_force_noneMB(int nboard, int *ni, double (**xi)[3],
			      int *nj, double (**xj)[3], double **mj,
			      double eps, double (**a)[3], double **p);
static void get_force_hostMB(int nboard, int *ni, double (**xi)[3],
			      int *nj, double (**xj)[3], double **mj,
			      double eps, double (**a)[3], double **p);
static void get_force_grapeMB(int nboard, int *ni, double (**xi)[3],
			      int *nj, double (**xj)[3], double **mj,
			      double eps, double (**a)[3], double **p);
static void get_potential_grapeMB(int nboard, int *ni, double (**xi)[3],
			      int *nj, double (**xj)[3], double **mj,
			      double eps, double (**a)[3], double **p);
void (*vtc_force_calculator[])(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
			       double eps, double (*a)[3], double *p) = {
				   get_force_none,          /* no force calculation */
				   get_force_host,          /* force & potential on host */
				   get_force_grape,         /* force & potential on GRAPE */
				   get_force_grape,         /* force on GRAPE */
				   get_potential_grape,     /* potential on GRAPE */
			       };
void (*vtc_force_calculatorMB[])(int nboard,
				 int *ni, double (**xi)[3], int *nj, double (**xj)[3], double **mj,
				 double eps, double (**a)[3], double **p) = {
				     get_force_noneMB,
				     get_force_hostMB,
				     get_force_grapeMB,
				     get_force_grapeMB,
				     get_potential_grapeMB,
				 };

void
get_force_none(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
	       double eps, double (*a)[3], double *p)
{
    int i, k;
    static int firstcall = 1;

    if (firstcall) {
	firstcall = 0;
	fprintf(stderr, "get_force_grape DOES NOT calculate force\n");
    }
    for (i = 0; i < ni; i++) {
	for (k = 0; k < 3; k++) {
	    a[i][k] = 0.0;
	}
	p[i] = 0.0;
    }
}

void
get_force_host(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
	       double eps, double (*a)[3], double *p)
{
    int i, j, k;
    double r, r2, r3;

    for (i = 0; i < ni; i++) {
	for (k = 0; k < 3; k++) {
	    a[i][k] = 0.0;
	}
	p[i] = 0.0;
	for (j = 0; j < nj; j++) {
	    r2 = eps*eps;
	    for (k = 0; k < 3; k++) {
		r2 += (xi[i][k]-xj[j][k])*(xi[i][k]-xj[j][k]);
	    }
	    if (r2 == 0.0) {
		r2 = 1e-10;
	    }
	    r3 = pow(r2, 1.5);
	    for (k = 0; k < 3; k++) {
		a[i][k] += mj[j]*(xj[j][k]-xi[i][k])/r3;
	    }
	    r = sqrt(r2);
	    p[i] -= mj[j]/r;
	}
    }
}

static void
get_force_noneMB(int nboard, int *ni, double (**xi)[3], int *nj, double (**xj)[3], double **mj,
		double eps, double (**a)[3], double **p)
{
    /* nop */
}

static void
get_force_hostMB(int nboard, int *ni, double (**xi)[3], int *nj, double (**xj)[3], double **mj,
		double eps, double (**a)[3], double **p)
{
    /* nop */
}


#ifdef NOGRAPE
void
vtc_set_scale(double xscale)
{
    /* nop */
}

void
get_force_grape(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
		double eps, double (*a)[3], double *p)
{
    /* nop */
}

void
get_potential_grape(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
		double eps, double (*a)[3], double *p)
{
    /* nop */
}

#else /* !NOGRAPE */

static double holdtime = 0.0;
static int is_grape_opened = 0;

#include "def.h"
#include "emulator_switch.h"
#include "mdgp2defs.h"
#include "m2_unit.h"
#include "m2_pipeline_allocator.h"

void md2_set_mode(int mode);
int md2_get_mode(void);
void md2_set_xmax(double xmax);
double md2_get_xmax(void);

static M2_UNIT *md2_unit = NULL;
static M2_UNIT *md2_units[NBOARDMAX];
static double md2_xmax = 32.0;
static int md2_mode = M2_FORCE;
static int md2_rescaled = TRUE;
static int md2_nboard = 0;

void
vtc_set_scale(double xscale)
{
    md2_set_xmax(xscale);
    md2_rescaled = TRUE;
    fprintf(stderr, ">>> rescaled. xscale: %f\n", xscale);
}

void
md2_set_mode(int mode)
{
    md2_mode = mode;
}

int
md2_get_mode(void)
{
    return (md2_mode);
}

void
md2_set_xmax(double xmax)
{
    md2_xmax = xmax;
}

double
md2_get_xmax(void)
{
    return (md2_xmax);
}

static double grav(double x){
	return pow(x, -1.5);
}

static double gravpot(double x){
	return pow(x, -0.5);
}

static M2_UNIT*
md2_alloc(char *tablename, int mode, double (*func)(double), double xmax)
{
    M2_UNIT *mu;

    mu = m2_allocate_unit(tablename, /* table file name */
			  mode, /* M2_FORCE or M2_POTENTIAL */
			  0, xmax, /* xmin and xmax */
			  NULL_INT); /* specify # of pipelines.
					set NULL_INT to use all pipes. */
    if (mu != NULL) {
	m2_set_host_function(mu, func);
	m2_set_rscale(mu, 1.0);
    }
    return (mu);
}

static M2_UNIT*
md2_alloc_onepb(char *tablename, int mode, double (*func)(double), double xmax)
{
    M2_UNIT *mu;

    mu = m2_allocate_unit(tablename, /* table file name */
			  mode, /* M2_FORCE or M2_POTENTIAL */
			  0, xmax, /* xmin and xmax */
			  MDGP2_NUM_PIPE*MDGP2_NUM_CHIP); /* get pipes on one board */
    if (mu != NULL) {
	m2_set_host_function(mu, func);
	m2_set_rscale(mu, 1.0);
    }
    return (mu);
}

double
grape_holdtime(void)
{
    return (holdtime);
}

int
grape_is_opened(void)
{
    return (is_grape_opened);
}

void
grape_openMB(void)
{
    int b = 0;

    if (is_grape_opened) {
	Cfprintf(stderr, "open_grape: already opened\n");
	return;
    }
    while (1) {
	if (md2_mode == M2_FORCE) {
	    md2_units[b] = md2_alloc_onepb("grav.table", md2_mode, grav, md2_xmax);
	}
	else {
	    md2_units[b] = md2_alloc_onepb("gravpot.table", md2_mode, gravpot, md2_xmax);
	}
	Cfprintf(stderr, "md2_units[%d]: %x\n", b, md2_units[b]);
	if (md2_units[b] != NULL) {
	    b++;
	}
	else {
	    break;
	}
    }
    holdtime = vtc_get_cputime();
    is_grape_opened = 1;
    md2_nboard = b;
}

void
grape_closeMB(void)
{
    int b;

    if (!is_grape_opened) {
	Cfprintf(stderr, "close_grape: not opened\n");
	return;
    }
    for (b = 0; b < md2_nboard; b++) {
	m2_free_unit(md2_units[b]);
	md2_units[b] = NULL;
    }
    is_grape_opened = 0;
    md2_nboard = b;
}

void
grape_open(void)
{
    if (is_grape_opened) {
	Cfprintf(stderr, "open_grape: already opened\n");
	return;
    }
    if (md2_mode == M2_FORCE) {
	md2_unit = md2_alloc("grav.table", md2_mode, grav, md2_xmax);
    }
    else {
	md2_unit = md2_alloc("gravpot.table", md2_mode, gravpot, md2_xmax);
    }
    holdtime = vtc_get_cputime();
    is_grape_opened = 1;
}

void
grape_close(void)
{
    if (!is_grape_opened) {
	Cfprintf(stderr, "close_grape: not opened\n");
	return;
    }
    m2_free_unit(md2_unit);
    md2_unit = NULL;
    is_grape_opened = 0;
}

int /* have a look at # of boards available */
vtc_get_grape(void)
{
    int ret;

    grape_openMB();
    ret = md2_nboard;
    grape_closeMB();
    return (ret);
}


#if 0
static void
get_force_grapeMB(int nboard, int *ni, double (**xi)[3], int *nj, double (**xj)[3], double **mj,
		double eps, double (**a)[3], double **p)
{
    int i, k, b;

    if ((grape_is_opened() && md2_get_mode() != M2_FORCE) || md2_rescaled){
	grape_close();
	md2_set_mode(M2_FORCE);
	md2_rescaled = FALSE;
    }
    if (!grape_is_opened()){
	md2_set_mode(M2_FORCE);
	grape_open();
	m2_set_softening(md2_unit, eps);
    }

    for (b = 0; b < nboard; b++) {
	m2_set_positions(md2_unit, xj[b], nj[b]);
	m2_set_charges(md2_unit, mj[b], nj[b]);
	{
	    M2_CELL c[1];
	    c[0].base = 0;
	    c[0].size = nj[b];
	    m2_set_cells(md2_unit, c, 1);
	}
	m2_calculate_forces(md2_unit, xi[b], ni[b], a[b]);
    }
}
#else
static void
get_force_grapeMB(int nboard, int *ni, double (**xi)[3], int *nj, double (**xj)[3], double **mj,
		double eps, double (**a)[3], double **p)
{
    int i, k, b;

    if ((grape_is_opened() && md2_get_mode() != M2_FORCE) || md2_rescaled){
	grape_closeMB();
	md2_set_mode(M2_FORCE);
	md2_rescaled = FALSE;
    }
    if (!grape_is_opened()){
	md2_set_mode(M2_FORCE);
	grape_openMB();
	fprintf(stderr, "allocated %d boards\n", md2_nboard);
	if (nboard > md2_nboard) {
	    fprintf(stderr, "md2_nboard: %d smaller than nboard: %d\n", md2_nboard, nboard);
	    exit(1);
	}
	for (b = 0; b < nboard; b++) {
	    m2_set_softening(md2_units[b], eps);
	}
    }
#if 1
    for (k = 0; k < 3; k++) {
	b = 0;
	m2_do_posconv(1);
	m2_set_positions_one_component(md2_units[b], xj[b], nj[b], k);
	m2_do_posconv(0);
	for (b = 1; b < md2_nboard; b++) {
	    m2_set_positions_one_component(md2_units[b], xj[b], nj[b], k);
	}
	m2_do_posconv(1);
    }
#else
    for (b = 0; b < md2_nboard; b++) {
	m2_set_positions(md2_units[b], xj[b], nj[b]);
    }
#endif
    for (b = 0; b < md2_nboard; b++) {
	m2_set_charges(md2_units[b], mj[b], nj[b]);
    }
    for (b = 0; b < nboard; b++) {
	m2_start_force_calculation(md2_units[b], xi[b], ni[b], a[b]);
    }
    for (b = 0; b < nboard; b++) {
	m2_wait_calculation(md2_units[b]);
    }
}
#endif

#if 1
extern void m2_set_positions_one_component(M2_UNIT *mu, double (*rj0)[3], long n, int k);

void
get_force_grape(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
		double eps, double (*a)[3], double *p)
{
    int i, k, b, nii, off = 0;

    if ((grape_is_opened() && md2_get_mode() != M2_FORCE) || md2_rescaled){
	grape_closeMB();
	md2_set_mode(M2_FORCE);
	md2_rescaled = FALSE;
    }
    if (!grape_is_opened()){
	md2_set_mode(M2_FORCE);
	grape_openMB();
	fprintf(stderr, "allocated %d boards\n", md2_nboard);
	for (b = 0; b < md2_nboard; b++) {
	    m2_set_softening(md2_units[b], eps);
	}
    }

    for (k = 0; k < 3; k++) {
	b = 0;
	m2_do_posconv(1);
	m2_set_positions_one_component(md2_units[b], (double (*)[3])xj, nj, k);
	m2_do_posconv(0);
	for (b = 1; b < md2_nboard; b++) {
	    m2_set_positions_one_component(md2_units[b], (double (*)[3])xj, nj, k);
	}
	m2_do_posconv(1);
    }
    for (b = 0; b < md2_nboard; b++) {
	m2_set_charges(md2_units[b], mj, nj);
    }

#if 1
    nii = ni/md2_nboard+1;
    for (b = 0; b < md2_nboard; b++) {
	m2_start_force_calculation(md2_units[b], xi+off, nii, a+off);
	off += nii;
	if (off+nii >= ni) {
	    nii = ni - off;
	}
    }
    for (b = 0; b < md2_nboard; b++) {
	m2_wait_calculation(md2_units[b]);
    }
#endif
}
#else
void
get_force_grape(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
		double eps, double (*a)[3], double *p)
{
    int i, k;

    if ((grape_is_opened() && md2_get_mode() != M2_FORCE) || md2_rescaled){
	grape_close();
	md2_set_mode(M2_FORCE);
	md2_rescaled = FALSE;
    }
    if (!grape_is_opened()){
	md2_set_mode(M2_FORCE);
	grape_open();
	m2_set_softening(md2_unit, eps);
    }
    m2_set_positions(md2_unit, (double (*)[3])xj, nj);
    m2_set_charges(md2_unit, mj, nj);
    {
	M2_CELL c[1];
	c[0].base = 0;
	c[0].size = nj;
	m2_set_cells(md2_unit, c, 1);
    }
/*
    m2_calculate_forces_test(md2_unit, (double (*)[3])xi, ni, (double (*)[3])a, 2);
*/
    m2_calculate_forces(md2_unit, (double (*)[3])xi, ni, (double (*)[3])a);
}
#endif

static void
get_potential_grapeMB(int nboard, int *ni, double (**xi)[3], int *nj, double (**xj)[3], double **mj,
		double eps, double (**a)[3], double **p)
{
}

void
get_potential_grape(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
		double eps, double (*a)[3], double *p)
{
    int i, k;

    if ((grape_is_opened() && md2_get_mode() != M2_POTENTIAL) || md2_rescaled) {
	grape_close();
	md2_set_mode(M2_POTENTIAL);
	md2_rescaled = FALSE;
    }
    if (!grape_is_opened()){
	md2_set_mode(M2_POTENTIAL);
	grape_open();
	m2_set_softening(md2_unit, eps);
    }
    m2_set_positions(md2_unit, (double (*)[3])xj, nj);
    m2_set_charges(md2_unit, mj, nj);
    {
	M2_CELL c[1];
	c[0].base = 0;
	c[0].size = nj;
	m2_set_cells(md2_unit, c, 1);
    }
    m2_calculate_potentials(md2_unit, (double (*)[3])xi, ni, p);
    for (i = 0; i < ni; i++) {
	p[i] *= -1.0;
    }
}
#endif /* NOGRAPE */
