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
				   get_force_grape,         /* potential on GRAPE */
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
	    r = sqrt(r2);
	    p[i] -= mj[j]/r;
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

#include "gp5util.h"

int
vtc_does_include_self_interaction(int calculator)
{
    return FALSE;
}

void
vtc_set_scale(double xmax, double mmin)
{
    if (!grape_is_opened()) {
	grape_open();
    }
    g5_set_range(-xmax, xmax, mmin);
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
    fprintf(stderr, "open GRAPE-5\n");
    g5_open();
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
    fprintf(stderr, "close GRAPE-5\n");
    g5_close();
    is_grape_opened = 0;
}

void
vtc_close_grape(void)
{
    grape_close();
}

#if 0

void
get_force_grape(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
		double eps, double (*a)[3], double *p)
{
    int offs, offr, nii, c, c0, i, np, nc, k;

    if (JMEMSIZE < nj) {
	fprintf(stderr, "nj: %d exceeded GRAPE-5 JMEMSIZE (%d)\n", nj, JMEMSIZE);
	exit(1);
    }
    if (!grape_is_opened()) {
	grape_open();
    }
    g5_set_mj(0, nj, mj);
    g5_set_xj(0, nj, xj);
    g5_set_n(nj);
    g5_set_eps_to_all(eps);

    g5_calculate_force_on_x(xi, a, p, ni);

    for (i = 0; i < ni; i++) {
	p[i] *= -1;
    }
    if (vtc_get_cputime()-grape_holdtime() > 150.0) {
	grape_close();
    }

    i = 10;
    fprintf(stderr, "!!! nj: %d a: %f %f %f p: %f\n",
	    nj, a[i][0], a[i][1], a[i][2], p[i]);
}

#else

void
get_force_grape(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj,
		double eps, double (*a)[3], double *p)
{
    int offs, offr, nii, c, c0, i, np, nc, k;

    if (JMEMSIZE < nj) {
	fprintf(stderr, "nj: %d exceeded GRAPE-5 JMEMSIZE (%d)\n", nj, JMEMSIZE);
	exit(1);
    }
    if (!grape_is_opened()) {
	grape_open();
    }
    np = g5_get_number_of_pipelines_per_board();
    nc = g5_get_number_of_boards();
    c0 = g5_get_firstcluster();

    //    g5_set_range(-100.0, 100.0, 1.0/4096);

    g5_set_mj(0, nj, mj);
    g5_set_xj(0, nj, xj);
    g5_set_n(nj);
    g5_set_eps_to_all(eps);

    offs = 0;
    for (c = c0; c < nc+c0 && offs < ni; c++) {
	nii = np;
	if (offs+nii > ni) {
	    nii = ni - offs;
	}
	g5_set_xiMC(c, nii, (double (*)[3])xi[offs]);
	g5_runMC(c);
	offs += nii;
    }

    for (offr = 0; offr < ni;) {
	for (c = c0; c < nc+c0; c++) {
	    if (offr < ni) {
		nii = np;
		if (offr+nii > ni) {
		    nii = ni - offr;
		}
		g5_get_forceMC(c, nii, (double (*)[3])a[offr], &p[offr]);
		offr += nii;
	    }
	    if (offs < ni) {
		nii = np;
		if (offs+nii > ni) {
		    nii = ni - offs;
		}
		g5_set_xiMC(c, nii, (double (*)[3])xi[offs]);
		g5_runMC(c);
		offs += nii;
	    }
	}
    }
    for (i = 0; i < ni; i++) {
	p[i] *= -1;
    }
    if (vtc_get_cputime()-grape_holdtime() > 150.0) {
	grape_close();
    }
}
#endif 

#endif /* NOGRAPE */
