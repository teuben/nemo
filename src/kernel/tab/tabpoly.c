/* TABPOLY:   roots of polynomials, evaluate polynomials
 *
 *
 *   7-sep-01    1.0	written , first example with GSL		pjt
 */

#include <stdinc.h> 
#include <getparam.h>

#ifdef HAVE_GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>
#else
#error GSL not enabled
#endif

string defv[] = {
    "coef=???\n     Polynomial coefficients, c_0 first",
    "x=\n           Evaluate polynial at these X's",
    "format=%g\n    Output format",
    "VERSION=1.0\n  7-sep-01 PJT",
    NULL,
};

string usage="Evaluate polynomials or find roots of polynomials";

#define MAXZERO     64
#define MAXDATA  16384
#define MAXCOEF    100

nemo_main()
{
    double x, y, det, ap, bp, cp, xr[MAXZERO], c[MAXCOEF], xp[MAXDATA];
    double z[2*MAXZERO];
    gsl_complex zr[MAXZERO];
    gsl_poly_complex_workspace *w;
    string fmt;
    char fmt1[100], fmt2[200];
    int i, j, n, nx, ny, nmax, nr, ncoef;

    ncoef = nemoinpd(getparam("coef"),c,MAXCOEF);
    if (ncoef < 0) error("Error %d parsing %s",ncoef,getparam("coef"));

    fmt = getparam("format");
    sprintf(fmt1,"%s\n",fmt,fmt);
    sprintf(fmt2,"%s %s\n",fmt,fmt);
    dprintf(1,"Using format=\"%s\"\n",fmt);

    if (hasvalue("x")) {

      nx = nemoinpr(getparam("x"),xp,MAXDATA);
      if (nx<0) error("Parsing x=%s",getparam("x"));
      for (j=0; j<nx; j++) {
        x = xp[j];
        y = gsl_poly_eval(c, ncoef, x);
        printf(fmt,x,y);
      }

    } else {

      if (ncoef < 2) {
        warning("No roots for f(x)=%g",c[0]);
      } else if (ncoef == 2) {
        x = -c[0]/c[1];
        printf(fmt1,x);
      } else if (ncoef == 3) {
        nr = gsl_poly_complex_solve_quadratic(c[2],c[1],c[0],&zr[0],&zr[1]);
        for (i=0; i<nr; i++)
            printf(fmt2,zr[i].dat[0],zr[i].dat[1]);
      } else if (ncoef == 4) {
        c[2] /= c[3];		/* renormalize coefficient !!! */
        c[1] /= c[3];
        c[0] /= c[3];
        nr = gsl_poly_complex_solve_cubic(c[2],c[1],c[0],&zr[0],&zr[1],&zr[2]);
        for (i=0; i<nr; i++)
            printf(fmt2,zr[i].dat[0],zr[i].dat[1]);
      } else {
        w = gsl_poly_complex_workspace_alloc(ncoef);
        nr = gsl_poly_complex_solve(c,ncoef,w,z);
	if (nr == GSL_EFAILED) error("QR reduction did not converge");
        for (i=0; i<ncoef-1; i++)
            printf(fmt2,z[2*i],z[2*i+1]);
      }
    }
}


