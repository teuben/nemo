
/*
 *  nr_nllsqfit:   wrapper function that acts like GIPSY's nllsqfit, but
 *                 calls the NumRec routine mrqmin() and helpers instead
 *                 to do the Marquand-Levenbergh
 *
 *    16-jul-2002  Created                                  Peter Teuben
 *    17-jul-2002  added lfit() call for linear fits                 PJT
 */


#include <stdinc.h>
#include <nr.h>
#include <nrutil.h>

#define MAXPAR 10

static rproc old_f;
static iproc old_df;

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

/* nr_func: for linear fits using lfit(), just returns the basis functions */
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
  float *a, chisq, ochisq;
  float lamda;
  int i, k, itst, *ia;
  int itc = 0;

  if (xdim != 1) error("nr: Cannot handle xdim > 1 yet (ever?)");
  if (npar > MAXPAR) error("MAXPAR too small");

  old_f  = f;
  old_df = df;

  x   = vector(1,ndat);            /* allocate NumRec arrays and matrices */
  y   = vector(1,ndat);
  sig = vector(1,ndat);
  covar = matrix(1,npar,1,npar);
  alpha = matrix(1,npar,1,npar);
  a = vector(1,npar);
  ia = ivector(1,npar);

  for (i=0; i<ndat; i++) {          /* copy our input data into NumRec ones */
    x[i+1] = xdat[i];
    y[i+1] = ydat[i];
    if (wdat)
      sig[i+1] = wdat[i];
    else
      sig[i+1] = 1.0;
  }
  for (i=0; i<npar; i++) {
    a[i+1] = fpar[i];
    ia[i+1] = mpar[i];
  }

  if (lab > 0.0) {
    lamda = -1;
    dprintf(1,"MRQMIN\n");
    mrqmin(x,y,sig,ndat, a,ia,npar, covar,alpha,&chisq,nr_func,&lamda);
  } else { 
    warning("nr_nllsqfit: testing linear fit");
    lfit(x,y,sig,ndat, a,ia,npar, covar,&chisq,nr_funcl);
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
    mrqmin(x,y,sig,ndat, a,ia,npar, covar,alpha,&chisq,nr_func,&lamda);
    itc++;
    if (chisq > ochisq)
      itst=0;
    else if (fabs(ochisq-chisq) < 0.1)
      itst++;
    if (itst < 4) continue;
    /* reached convergence, reset lamda=0 to get covar */
    lamda=0.0;
    mrqmin(x,y,sig,ndat, a,ia,npar, covar,alpha,&chisq,nr_func,&lamda);
    itc++;
    dprintf(1,"Uncertainties:\n");
    for (i=0;i<npar;i++) dprintf(1,"%9.4f",sqrt(covar[i+1][i+1]));
    dprintf(1,"\n");
    break;
  }

  for (i=0;i<npar;i++) {                /* copy info back to our arrays */
    epar[i] = sqrt(covar[i+1][i+1]);
    fpar[i] = a[i+1];
  }

  free_matrix(alpha,1,npar,1,npar);      /* free up the NumRec arrays and matrices */
  free_matrix(covar,1,npar,1,npar);
  free_vector(sig,1,ndat);
  free_vector(y,1,ndat);
  free_vector(x,1,ndat);
  free_ivector(ia,1,npar);

  return itc;                       /* return number of iterations */

}

#if defined(TESTBED)
/*
 * For testing purposes only. We try to fit a one-dimensional Gaussian
 * to data with a uniform noise distribution. Although the algorithm
 * is developed for Gaussian noise patterns, it should not cause to
 * much problems. For Poison noise the algorithm is not suitable!
 */

#if     defined(LINEAR)
   int  fopt = 2;
#else
   int  fopt = 1;
#endif

real func_c( real *xdat, real *fpar, int *npar )
/*
 * if *fopt == 1:
 * f(x) = a + b * exp( -1/2 * (c - x)^2 / d^2 ) (one-dimension Gauss)
 * if *fopt == 2:
 * f(x) = a + b * x + c * x^2 + d * x^3         (polynomial)
 */
{
   if (fopt == 1) {
      real a;
      real b;
      real arg;
   
      a = fpar[2] - xdat[0];
      b = fpar[3];
      arg = 0.5 * a / b * a / b;
      return( fpar[0] + fpar[1] * exp( -arg ) );
   } else if (fopt == 2) {
      real x = *xdat;

      return( fpar[0] + fpar[1] * x + fpar[2] * x * x + fpar[3] * x * x * x );
   }
   return( 0.0 );
}

void derv_c( real *xdat, real *fpar, real *epar, int *npar )
{
   if (fopt == 1) {
      real a;
      real b;
      real arg;
   
      a = fpar[2] - xdat[0];
      b = fpar[3];
      arg = 0.5 * a / b * a / b;
      epar[0] = 1.0;
      epar[1] = exp( -arg );
      epar[2] = -fpar[1] * epar[1] * a / b / b;
      epar[3] = fpar[1] * epar[1] * a * a / b / b / b;
   } else if (fopt == 2) {
      real x = *xdat;

      epar[0] = 1.0;
      epar[1] = x;
      epar[2] = x * x;
      epar[3] = x * x * x;      
   }
}

#if defined(mc68k)
#define RAND_MAX 32767
#else
#define RAND_MAX 2147483647
#endif

void main( )
{
   int  m;
   int  n;
   int  r;
   int  xdim = 1;
   int  its = 50;
   int  ndat = 30;
   int  nfit;
   int  npar = 4;
   int  mpar[4];
   real tol = 0.0;
#if     defined(LINEAR)
   real lab = 0.0;
#else
   real lab = 0.01;
#endif
   real xdat[30];
   real ydat[30];
   real wdat[30];
   real fpar[4];
   real epar[4];


   mpar[0] = 1; mpar[1] = 0; mpar[2] = 1; mpar[3] = 1;
   for (nfit = 0; nfit < 10; nfit++) {
      fpar[0] = 0.0; fpar[1] = 10.0; fpar[2] = 15.0; fpar[3] = 2.0;
      for (n = 0; n < ndat; n++) {
         real rndm;
         xdat[n] = (real) n;
         wdat[n] = 1.0;
         ydat[n] = func_c( &xdat[n], fpar, &npar );
         rndm = ((real) ( rand( ) - RAND_MAX / 2 )) / (real) RAND_MAX * 1.0;
#if defined(SCALE)
	 rndm *= (real) nfit;
#endif
         ydat[n] += rndm;
      }
      for (m = 0; m < npar; m++) {
         real rndm;
         rndm = (real) ( rand( ) - RAND_MAX / 2 ) / (real) RAND_MAX * 3.0;
         fpar[n] += rndm;
      }
      printf( "nllsqfit: on entry\n" );
      printf( "fpar: %10f %10f %10f %10f\n", fpar[0], fpar[1], fpar[2], fpar[3] );
      r = nllsqfit( xdat, xdim, ydat, wdat, ndat, fpar, epar, mpar, npar, tol, 
                its, lab, func_c, derv_c );
      printf( " %ld", r );
      if (r < 0) {
         printf( ", 	error!\n" );
      } else {
         printf( ", 	success!\n" );
         printf( "fpar: %10f %10f %10f %10f\n", fpar[0], fpar[1], fpar[2], fpar[3] );
         printf( "epar: %10f %10f %10f %10f\n", (real) epar[0], (real) epar[1], (real) epar[2], (real) epar[3] );
      }
   }
}
#endif
