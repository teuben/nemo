#define _FORCE_C_ (1)
#include <stdio.h>
#include <math.h>
#include "vtc.h"
#include "vtclocal.h"

#include "force_host.c"

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

static double holdtime = 0.0;
static int is_grape_opened = 0;

#include "gp5util.h"

void
vtc_set_scale(double xmax, double mmin)
{
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
    np = g5_get_number_of_pipelines_per_board();
    nc = g5_get_number_of_boards();
    c0 = g5_get_firstcluster();

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
    if (vtc_get_cputime()-grape_holdtime() > 15.0) {
	grape_close();
    }
}

