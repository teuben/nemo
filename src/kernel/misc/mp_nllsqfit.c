/*
 *  mp_nllsqfit:   wrapper function that acts like GIPSY's nllsqfit, but
 *                 calls MINPACK's mpfit() to do a Marquardt-Levenberg
 *                 minimization
 *
 *     9-oct-2013  Created                                  Peter Teuben
 */


#include <stdinc.h>
#include <mpfit.h>

#define MAXPAR 10

static rproc old_f;
static iproc old_df;

static int old_xdim = -1;

/* my_func: wrapper to our function */

int my_func(int m, int n, double *x, double *fvec, double **dvec, void *private_data)
{
  static real old_x, old_a[MAXPAR], old_d[MAXPAR];
  static int i;

#if 0

  /* to be called as follows: */
  // funct(m,n,x,fvec,dvec,priv);

  old_x = x;
  for (i=0; i<na; i++)
    old_a[i] = a[i+1];
  *y = (*old_f)(&old_x, old_a, na);
  (*old_df)(&old_x, old_a, old_d, na);
  for (i=0; i<na; i++)
    dyda[i+1] = old_d[i];
  dprintf(2,"my_func(%d): x=%g y=%g a[1]=%g a[2]=%g\n",na,x,*y,a[1],a[2]);
#endif
}


int mp_nllsqfit(
    real *xdat,       /*  x[ndat][xdim]   or   x(xdim,ndat)   */
    int xdim,         /*  */
    real *ydat,       /*  y[ndat] */
    real *wdat,       /*  w[ndat] */
    real *ddat,       /*  d[ndat] */
    int ndat,         /*  */
    real *fpar,       /*  f[npar] */
    real *epar,       /*  e[npar] */
    int *mpar,        /*  m[npar] */
    int npar,         /*  */
    real tol,         /* tolerance to convergence */
    int its,          /* # iterations */
    real lab,         /* (small) mixing parameter (0 for linear) */
    rproc f,          /*  f */
    iproc df)         /*  df/da */
{
  float *x, *y, *sig, **covar, **alpha;
  float *a, chisq, ochisq, lamda, fac;
  int i, k, itst, *ia;
  int itc = 0;
  mp_config  *config = NULL;
  mp_par *pars = NULL;
  void *private_data = NULL;
  mp_result result;

  if (npar > MAXPAR) error("MAXPAR too small for npar");
  if (xdim > MAXPAR) error("MAXPAR too small for xdim");


  if (xdim != 1) error("mp_nllsqfit: cannot deal with xdim=%d",xdim);

  old_f = f;   /* mpfit only uses f, doesn't need df */

  mpfit(my_func, ndat, npar,
	xdat, pars, config, private_data, 
	&result);


  /* to be called as follows: */
  // funct(m,n,x,fvec,dvec,priv);

  return 0;                       /* return number of iterations */
}

