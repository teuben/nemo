#define _FORCE_C_ (1)
#include <stdio.h>
#include <math.h>
#include "vtc.h"
#include "vtclocal.h"

static void get_force_none(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
			      double eps, double (*a)[3], double *p);
static void get_force_and_potential_host(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
			   double eps, double (*a)[3], double *p);
static void get_force_host(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
			   double eps, double (*a)[3], double *p);
static void get_potential_host(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
			   double eps, double (*a)[3], double *p);
static void get_force_grape(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
			    double eps, double (*a)[3], double *p);
static void get_potential_grape(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
				double eps, double (*a)[3], double *p);


void (*vtc_force_calculator[])(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
			       double eps, double (*a)[3], double *p) = {
				   get_force_none,          /* no force calculation */
				   get_force_and_potential_host,          /* force & potential on host */
				   get_force_host,          /* force on host */
				   get_potential_host,          /* potential on host */
				   get_force_grape,         /* force & potential on GRAPE */
				   get_force_grape,         /* force on GRAPE */
				   get_potential_grape,     /* potential on GRAPE */
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
	p[i] = 1e-5;
    }
}

void
get_force_and_potential_host(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
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
	}
    }
}

void
get_potential_host(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
		   double eps, double (*a)[3], double *p)
{
    int i, j, k;
    double r, r2;

    for (i = 0; i < ni; i++) {
        p[i] = 0.0;
        for (j = 0; j < nj; j++) {
            r2 = eps*eps;
            for (k = 0; k < 3; k++) {
                r2 += (xi[i][k]-xj[j][k])*(xi[i][k]-xj[j][k]);
            }
            if (r2 != 0.0) { /* avoid self interaction for zero-softening calculation */
                r = sqrt(r2);
                p[i] -= mj[j]/r;
            }
        }
    }
}

#ifdef NOGRAPE

int
vtc_does_include_self_interaction(int calculator)
{
    return FALSE;
}

void
vtc_set_scale(double xmax, double mmin)
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
#include "m2_unit.h"
#include "m2_pipeline_allocator.h"

void md2_set_mode(int mode);
int md2_get_mode(void);
void md2_set_xmax(double xmax);
double md2_get_xmax(void);

static M2_UNIT *md2_unit = NULL;
static double md2_xmax = 32.0;
static int md2_mode = M2_FORCE;
static int md2_rescaled = TRUE;

int
vtc_does_include_self_interaction(int calculator)
{
    int ret;

    switch (calculator) {
    case GRAPE:
    case GRAPE_FORCEONLY:
    case GRAPE_POTENTIALONLY:
	ret = TRUE;
	break;
    default:
	ret = FALSE;
    }
    return ret;
}

void
vtc_set_scale(double xmax, double mmin)
{
    md2_set_xmax(xmax);
    md2_rescaled = TRUE;
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
    M2_UNIT *mu = NULL;
    int ntry = 0;

    while (!mu) {
	mu = m2_allocate_unit(tablename, /* table file name */
			      mode, /* M2_FORCE or M2_POTENTIAL */
			      0, /* period */
			      xmax, /* xmax */
#if 0
			      NULL_INT); /* specify # of pipelines.
					    set NULL_INT to use all pipes. */
#else
	MD2_NPIPES*1);
#endif
	
	switch (ntry) {
	case 0:
	    break;
	case 1:
	    fprintf(stderr, "somebody is using MDGRAPE-2. wait");
	    break;
	default:
	    fprintf(stderr, ".");
	    break;
	}
	if (ntry != 0) {
	    sleep(5);
	}
	ntry++;
    }
    m2_set_host_function(mu, func);
    m2_set_rscale(mu, 1.0);
    m2_set_neighbor_radius(mu, 0);
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

void
vtc_close_grape(void)
{
    grape_close();
}

int
vtc_set_rscale(double rscale)
{
    if (!grape_is_opened()){
	return (-1);
    }
    m2_set_rscale(md2_unit, rscale);
    return (0);
}

void
get_force_grape(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
		double eps, double (*a)[3], double *p)
{
    int i, k, off, nj0;
    static double atmp[JMEMMAX][3];

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

    if (nj < JMEMMAX) { /* normal operation */
	m2_set_charges(md2_unit, mj, nj);
	m2_set_positions(md2_unit, (double (*)[3])xj, nj);
	m2_calculate_forces(md2_unit, (double (*)[3])xi, ni, (double (*)[3])a);
    }
    else { /* for very long interaction list */
	for (i = 0; i < ni; i++) {
	    for (k = 0; k < 3; k++) {
		a[i][k] = 0.0;
	    }
	}
	off = 0;
	nj0 = JMEMMAX;
	while (off < nj) {
	    if (off + nj0 > nj) {
		nj0 = nj - off;
	    }
	    m2_set_charges(md2_unit, mj+off, nj0);
	    m2_set_positions(md2_unit, (double (*)[3])xj+off, nj0);
	    m2_calculate_forces(md2_unit, (double (*)[3])xi, ni, (double (*)[3])atmp);
	    for (i = 0; i < ni; i++) {
		for (k = 0; k < 3; k++) {
		    a[i][k] += atmp[i][k];
		}
	    }
	    off += nj0;
	}
    }
}

void
get_potential_grape(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
		double eps, double (*a)[3], double *p)
{
    int i, off, nj0;
    static double ptmp[JMEMMAX];

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
    if (nj < JMEMMAX) { /* normal operation */
	m2_set_charges(md2_unit, mj, nj);
	m2_set_positions(md2_unit, (double (*)[3])xj, nj);
	m2_calculate_potentials(md2_unit, (double (*)[3])xi, ni, p);
    }
    else { /* for very long interaction list */
	for (i = 0; i < ni; i++) {
	    p[i] = 0.0;
	}
	off = 0;
	nj0 = JMEMMAX;
	while (off < nj) {
	    if (off + nj0 > nj) {
		nj0 = nj - off;
	    }
	    m2_set_charges(md2_unit, mj+off, nj0);
	    m2_set_positions(md2_unit, (double (*)[3])xj+off, nj0);
	    m2_calculate_potentials(md2_unit, (double (*)[3])xi, ni, ptmp);
	    for (i = 0; i < ni; i++) {
		p[i] += ptmp[i];
	    }
	    off += nj0;
	}
    }
    for (i = 0; i < ni; i++) {
	p[i] *= -1.0;
    }
}
#endif /* NOGRAPE */
