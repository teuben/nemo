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
				 get_potential_grape,     /* potential on GRAPE */
			       };
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
static double md2_xmax = 32.0;
static int md2_mode = M2_FORCE;
static int md2_rescaled = TRUE;

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
			      0, xmax, /* xmin and xmax */
			      NULL_INT); /* specify # of pipelines.
					    set NULL_INT to use all pipes. */
	
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
#if 1
    for (i = 0; i < ni; i++) {
	p[i] *= -1.0;
    }
#endif
}

