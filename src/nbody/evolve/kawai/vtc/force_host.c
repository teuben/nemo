#ifndef _FORCE_C_
#define _FORCE_C_ (1)
#endif /* _FORCE_C_ */
#include <stdio.h>
#include <math.h>
#include "vtc.h"
#include "vtclocal.h"

#define PROTOTYPEOF(func) static void func(int ni, double (*xi)[3], int nj, double (*xj)[3], double *mj, double eps, double (*a)[3], double *p)

PROTOTYPEOF(get_force_none);
PROTOTYPEOF(get_force_and_potential_host);
PROTOTYPEOF(get_force_host);
PROTOTYPEOF(get_potential_host);
PROTOTYPEOF(get_force_grape);
PROTOTYPEOF(get_potential_grape);

#ifdef NOGRAPE
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
#endif

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

#if NOGRAPE
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
#endif /* NOGRAPE */
