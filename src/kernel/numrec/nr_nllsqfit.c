
/*
 *  nr_nllsqfit:   wrapper function that acts like GIPSY's nllsqfit, but
 *                 calls the NumRec routine mrqmin() and helpers instead
 *                 to do the Marquardt-Levenberg
 *
 *    16-jul-2002  Created                                  Peter Teuben
 *    17-jul-2002  added lfit() call for linear fits                 PJT
 *    18-jul-2000  added support for xdim > 1                        PJT
 *    18-apr-2004  weights were computed wrong                       PJT
 */


#include <stdinc.h>
#include <nr.h>
#include <nrutil.h>

#define MAXPAR 10

static rproc old_f;
static iproc old_df;

static int old_xdim = -1;

/* nr_func: for non-linear fits */

void nr_func(float x, float *a, float *y, float *dyda, int na)  
{
  static real old_x, old_a[MAXPAR], old_d[MAXPAR];
  static int i;

  old_x = x;
  for (i=0; i<na; i++)
    old_a[i] = a[i+1];
  *y = (*old_f)(&old_x, old_a, na);
  (*old_df)(&old_x, old_a, old_d, na);
  for (i=0; i<na; i++)
    dyda[i+1] = old_d[i];
  dprintf(2,"nr_func(%d): x=%g y=%g a[1]=%g a[2]=%g\n",na,x,*y,a[1],a[2]);
}

/* nr_funcx:  for non-linear fits that also allow xdim > 1        */
/*           needs patches NumRec routines mrqminx and mrqcofx   */

void nr_funcx(float *x, float *a, float *y, float *dyda, int na)  
{
  static real old_x[MAXPAR], old_a[MAXPAR], old_d[MAXPAR];
  static int i;

  for (i=0; i<old_xdim; i++)
    old_x[i] = x[i];
  for (i=0; i<na; i++)
    old_a[i] = a[i+1];
  *y = (*old_f)(old_x, old_a, na);
  (*old_df)(old_x, old_a, old_d, na);
  for (i=0; i<na; i++)
    dyda[i+1] = old_d[i];
  dprintf(2,"nr_funcx(%d): x=%g y=%g a[1]=%g a[2]=%g\n",na,*x,*y,a[1],a[2]);
}

/* nr_funcl: for linear fits using lfit(), just returns the basis functions */
/*          it has a cumbersome way to get the basis functions, the way  */
/*          lfit() wants them to our f/df functions that nllsqfit() wants them */

void nr_funcl(float xa, float *a, int na)
{
  static real x, p[MAXPAR];
  static int i;

  x = xa;
  for (i=0; i<na; i++)
    p[i] = 0;
  for (i=0; i<na; i++) {
    p[i] = 1;
    a[i+1] = (*old_f)(&x, p, na);
    p[i] = 0;
  }
  dprintf(2,"nr_funcl(%d): x=%g\n",na,x);
}

void nr_funclx(float *xa, float *a, int na)
{
  static real x[MAXPAR], p[MAXPAR];
  static int i;

  for (i=0; i<old_xdim; i++)
    x[i] = xa[i];
  for (i=0; i<na; i++)
    p[i] = 0;
  for (i=0; i<na; i++) {
    p[i] = 1;
    a[i+1] = (*old_f)(x, p, na);
    p[i] = 0;
  }
}

int nr_nllsqfit(
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

  if (npar > MAXPAR) error("MAXPAR too small for npar");
  if (xdim > MAXPAR) error("MAXPAR too small for xdim");

  old_f  = f;
  old_df = df;
  old_xdim = xdim;

  x   = fvector(1,ndat*xdim);        /* allocate NumRec arrays and matrices */
  y   = fvector(1,ndat);
  sig = fvector(1,ndat);
  covar = fmatrix(1,npar,1,npar);
  alpha = fmatrix(1,npar,1,npar);
  a = fvector(1,npar);
  ia = ivector(1,npar);

  for (i=0; i<ndat; i++) {          /* copy our input data into NumRec ones */
    y[i+1] = ydat[i];
    if (wdat)
      sig[i+1] = 1/sqrt(wdat[i]);           /* mrqmin uses errors, not weights .... */
    else
      sig[i+1] = 1.0;
  }
  for (i=0; i<ndat*xdim; i++)
    x[i+1] = xdat[i];
  for (i=0; i<npar; i++) {
    a[i+1] = fpar[i];
    ia[i+1] = mpar[i];
  }

  if (lab > 0.0) {
    lamda = -1;
    dprintf(1,"MRQMIN\n");
    mrqminx(x,xdim,y,sig,ndat, a,ia,npar, covar,alpha,&chisq,nr_funcx,&lamda);
  } else { 
    warning("nr_nllsqfit: testing linear fit using lfit");
    lfitx(x,xdim,y,sig,ndat, a,ia,npar, covar,&chisq,nr_funclx);
  }
  k=1;
  itst=0;
  while (lab > 0.0) {
    dprintf(1,"Iteration # %2d chi-squared: %10.4f lamda: %9.2e\n",k,chisq,lamda);
    dprintf(1,"a:::");
    for (i=0;i<npar;i++) dprintf(1,"%9.4f",a[i+1]);
    dprintf(1,"\n");
    k++;
    ochisq=chisq;
    mrqminx(x,xdim,y,sig,ndat, a,ia,npar, covar,alpha,&chisq,nr_funcx,&lamda);
    itc++;
    if (chisq > ochisq)
      itst=0;
    else if (fabs(ochisq-chisq) < 0.1)
      itst++;
    if (itst < 4) continue;
    /* reached convergence, reset lamda=0 to get covar */
    lamda=0.0;
    mrqminx(x,xdim,y,sig,ndat, a,ia,npar, covar,alpha,&chisq,nr_funcx,&lamda);
    itc++;
    dprintf(1,"Uncertainties:\n");
    for (i=0;i<npar;i++) dprintf(1,"%9.4f",sqrt(covar[i+1][i+1]));
    dprintf(1,"\n");
    break;
  }

  if (wdat)
    fac = 1;
  else {
    fac = sqrt(chisq/(ndat-npar));
  }

  for (i=0;i<npar;i++) {                /* copy info back to our arrays */
    epar[i] = fac * sqrt(covar[i+1][i+1]);
    fpar[i] = a[i+1];
  }

  free_fmatrix(alpha,1,npar,1,npar);      /* free up the NumRec arrays and matrices */
  free_fmatrix(covar,1,npar,1,npar);
  free_fvector(sig,1,ndat);
  free_fvector(y,1,ndat);
  free_fvector(x,1,ndat*xdim);
  free_ivector(ia,1,npar);

  dprintf(1,"chisq=%g\n",chisq);
  return itc;                       /* return number of iterations */
}

