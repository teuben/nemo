/*
        Copyright (c) Kapteyn Laboratorium Groningen 1990
        All Rights Reserved.

Function:     NLLSQFIT

Purpose:      NLLSQFIT is a routine for making a least-squares fit of a
              function to a set of data points. The method used is
              described in: Marquardt, J.Soc.Ind.Appl.Math. 11, 431 (1963).
              This  method is a mixture of the steepest descent method and
              the Taylor method.

Category:     MATH

File:         nllsqfit.c

Author:       K.G. Begeman, P.J. Teuben (NEMO adaptations)

Use:          INTEGER NLLSQFIT( XDAT ,    Input      REAL ARRAY
                                XDIM ,    Input      INTEGER
                                YDAT ,    Input      REAL ARRAY
                                WDAT ,    Input      REAL ARRAY
                                NDAT ,    Input      INTEGER
                                FPAR ,   In/Output   REAL ARRAY
                                EPAR ,    Output     REAL ARRAY
                                MPAR ,    Input      INTEGER ARRAY
                                NPAR ,    Input      INTEGER
                                TOL  ,    Input      REAL
                                ITS  ,    Input      INTEGER
                                LAB  ,    Input      REAL)

            NLLSQFIT   Returns number of iterations needed to achieve
                       convergence according to TOL. When this
                       number is negative, the fitting was not
                       continued because a fatal error occurred:
                       -1 Too many free parameters, maximum is 32.
                       -2 No free parameters.
                       -3 Not enough degrees of freedom.
                       -4 Maximum number of iterations too small to
                          obtain a solution which satisfies TOL.
                       -5 Diagonal of matrix contains elements which
                          are zero.
                       -6 Determinant of the coefficient matrix is zero.
                       -7 Square root of negative number.
              XDAT     Contains coordinates of data points.
                       XDAT is two-dimensional: XDAT(XDIM,NDAT)
              XDIM     Dimension of fit.
              YDAT     Contains data points.
              WDAT     Contains weigths for data points.
              NDAT     Number of data points.
              FPAR     On input contains initial estimates of the
                       parameters for non-linear fits, on output the
                       fitted parameters.
              EPAR     Contains estimates of errors in fitted
                       parameters.
              MPAR     Logical mask telling which parameters are
                       free (MPAR(J)=non-zero) and which parameters
                       are fixed (MPAR(J)=0).
              NPAR     Number of parameters (free+fixed).
              TOL      Relative tolerance. NLLSQFIT stops when
                       successive iterations fail to produce a
                       decrement in reduced chi-squared less than
                       TOL. If TOL is less than the minimum tolerance
                       possible, TOL will be set to this value. This
                       means that maximum accuracy can be obtained by
                       setting TOL=0.0.
              ITS      Maximum number of iterations.
              LAB      Mixing parameter, LAB determines the initial
                       weight of steepest descent method relative to
                       the Taylor method. LAB should be a small
                       value (i.e. 0.01). LAB can only be zero when
                       the partial derivatives are independent of
                       the parameters. In fact in this case LAB
                       should be exactly equal to zero.
	      FUNC	External Function, see below
	      DERV	External Subroutine, see below

Notes:        The following routines have to be defined by the user:
              
              REAL FUNC( XDAT ,   Input    REAL ARRAY
                         FPAR ,   Input    REAL ARRAY
                         NPAR ,   Input    INTEGER)

              FUNC    Returns the function value of the function to
                      be fitted.
              XDAT    Coordinate(s) of data point.
              FPAR    Parameter list.
              NPAR    Number of parameters.
              
              CALL DERV( XDAT ,   Input    REAL ARRAY
                         FPAR ,   Input    REAL ARRAY
                         DPAR ,   Output   REAL ARRAY
                         NPAR ,   Input    INTEGER)

              XDAT    Coordinate(s) of data point.
              FPAR    Parameter list.
              EPAR    Partial derivatives to the parameters of the
                      function to be fitted.
              NPAR    Number of parameters.

              In FORTRAN applications you need to specify the
              following f2cvv syntax for the C to Fortran
              interface somewhere in your source code:
              
              C@ real function func( real, real, integer, integer )
              C@ subroutine derv( real, real, real, integer, integer )
              
Example:      Fitting y(x) = a + b * x 

              REAL FUNCTION FUNC( XDAT, FPAR, NPAR)
              ...
              FUNC = FPAR(1) + FPAR(2) * XDAT
              RETURN
              END
              
              SUBROUTINE DERV( XDAT, FPAR, DPAR, NPAR)
              ...
              DPAR(1) = 1.0
              DPAR(2) = XDAT
              RETURN
              END

Updates:      May  7, 1990: KGB, Document created.
              May 14, 1990: MXV, Document refereed.
              Apr 30, 1991: PJT, NEMO version, more like old 'fit'
              Oct 15, 1999: PJT  Added residual computations
              Jun 20, 2001: PJT  gcc3 prototpypes 
	      Jul 12, 2002: PJT  allow wdat to be NULL, in which case all weights = 1 (deja vu???)
              Apr 18, 2004: PJT  fixed wdat normalization error for chi2 computation

*/

/*
 *    Here are the steps that every nonlinear regression program follows: (or should?)

    1.Start with an initial estimated value for each variable in the equation. 
    2.Generate the curve defined by the initial values. Calculate the sum-of-squares (the
      sum of the squares of the vertical distances of the points from the curve). 
    3.Adjust the variables to make the curve come closer to the data points. There are
      several algorithms for adjusting the variables. The most commonly used method was
      derived by Levenberg and Marquardt (often called simply the Marquardt method). 
    4.Adjust the variables again so that the curve comes even closer to the points. 
    5.Keep adjusting the variables until the adjustments make virtually no difference in the
      sum-of-squares. 
    6.Report the best-fit results. The precise values you obtain will depend in part on the
      initial values chosen in step 1 and the stopping criteria of step 5. This means that
      repeat analyses of the same data will not always give exactly the same results. 

 */

#include <stdinc.h>
#if !defined(FLT_EPSILON)
#define FLT_EPSILON 1.0e-6
#endif


#define LABFAC  10.0                            /* labda step factor */
#define LABMAX  1.0e+10                         /* maximum value for labda */
#define LABMIN  1.0e-10                         /* minimum value for labda */
#define MAXPAR  32                              /* number of free parameters */

static  real  chi1;                           /* old reduced chi-squared */
static  real  chi2;                           /* new reduced chi-squared */
static  real  labda;                          /* mixing parameter */
static  real  tolerance;                      /* accuracy */
static  real  vector[MAXPAR];                 /* correction vector */
static  real  matrix1[MAXPAR][MAXPAR];        /* original matrix */
static  real  matrix2[MAXPAR][MAXPAR];        /* inverse of matrix1 */
static  int    itc;                            /* fate of fit */
static  int    found;                          /* solution found ? */
static  int    nfree;                          /* number of free parameters */
static  int    nuse;                           /* number of useable data points */
static  int    parptr[MAXPAR];                 /* parameter pointer */

typedef real (*my_proc1)(real *, real *, int);
typedef void (*my_proc2)(real *, real *, real *, int);

static my_proc1 fitfunc_c;
static my_proc2 fitderv_c;

static int invmat()
/*
 * invmat calculates the inverse of matrix2. The algorithm used is the
 * Gauss-Jordan algorithm described in Stoer, Numerische matematik, 1 Teil.
 * Did Tom Oosterloo write the original version of this routine?
 */
{
   real even;
   real hv[MAXPAR];
   real mjk;
   real rowmax;
   int   evin;
   int   i;
   int   j;
   int   k;
   int   per[MAXPAR];
   int   row;

   for (i = 0; i < nfree; i++) per[i] = i;      /* set permutation array */
   for (j = 0; j < nfree; j++) {                /* in j-th column, ... */
      rowmax = fabs( matrix2[j][j] );           /* determine row with ... */
      row = j;                                  /* largest element. */
      for (i = j + 1; i < nfree; i++) {
         if (fabs( matrix2[i][j] ) > rowmax) {
            rowmax = fabs( matrix2[i][j] );
            row = i;
         }
      }
      if (matrix2[row][j] == 0.0) return( -6 ); /* determinant is zero! */
      if (row > j) {                            /* if largest element not ... */
         for (k = 0; k < nfree; k++) {          /* on diagonal, then ... */
            even = matrix2[j][k];               /* permutate rows. */
            matrix2[j][k] = matrix2[row][k];
            matrix2[row][k] = even;
         }
         evin = per[j];                         /* keep track of permutation */
         per[j] = per[row];
         per[row] = evin;
      }
      even = 1.0 / matrix2[j][j];               /* modify column */
      for (i = 0; i < nfree; i++) matrix2[i][j] *= even;
      matrix2[j][j] = even;
      for (k = 0; k < j; k++) {
         mjk = matrix2[j][k];
         for (i = 0; i < j; i++) matrix2[i][k] -= matrix2[i][j] * mjk;
         for (i = j + 1; i < nfree; i++) matrix2[i][k] -= matrix2[i][j] * mjk;
         matrix2[j][k] = -even * mjk;
      }
      for (k = j + 1; k < nfree; k++) {
         mjk = matrix2[j][k];
         for (i = 0; i < j; i++) matrix2[i][k] -= matrix2[i][j] * mjk;
         for (i = j + 1; i < nfree; i++) matrix2[i][k] -= matrix2[i][j] * mjk;
         matrix2[j][k] = -even * mjk;
      }
   }
   for (i = 0; i < nfree; i++) {                /* finally, repermute the ... */
      for (k = 0; k < nfree; k++) {             /* columns. */
         hv[per[k]] = matrix2[i][k];
      }
      for (k = 0; k < nfree; k++) {
         matrix2[i][k] = hv[k];
      }
   }
   return( 0 );                                 /* all is well */
} /* invmat */


static void getmat(                             /* build up the matrix */
        real *xdat, int xdim, 
        real *ydat, real *wdat, real *ddat, int ndat,
        real *fpar, real *epar, int npar)
{
   real wd;
   real wn;
   real yd;
   int   i;
   int   j;
   int   n;

   for (j = 0; j < nfree; j++) {
      vector[j] = 0.0;                          /* zero vector ... */
      for (i = 0; i <= j; i++) {                /* and matrix ... */
         matrix1[j][i] = 0.0;                   /* only on and below diagonal */
      }
   }
   chi2 = 0.0;                                  /* reset reduced chi-squared */
   for (n = 0; n < ndat; n++) {              /* loop trough data points */
      wn = wdat ? wdat[n] : 1.0;
      if (wn > 0.0) {                           /* legal weight ? */
         (*fitderv_c)( &xdat[xdim * n], fpar, epar, npar );
         yd = ydat[n] - (*fitfunc_c)( &xdat[xdim * n], fpar, npar );
         if (ddat) ddat[n] = yd;
         chi2 += yd * yd * wn;                  /* add to chi-squared */
         for (j = 0; j < nfree; j++) {
            wd = epar[parptr[j]] * wn;          /* weighted derivative */
            vector[j] += yd * wd;               /* fill vector */
            for (i = 0; i <= j; i++) {          /* fill matrix */
               matrix1[j][i] += epar[parptr[i]] * wd;
            }
         }
      } 
   }
} /* getmat */

static int getvec(
    real *xdat, int xdim, 
    real *ydat, real *wdat, int ndat, 
    real *fpar, real *epar, int npar)
/*
 * getvec calculates the correction vector. The matrix has been built by
 * getmat, we only have to rescale it for the current value for labda.
 * The matrix is rescaled so that the diagonal gets the value 1 + labda.
 * Next we calculate the inverse of the matrix and then the correction
 * vector.
 */
{
   real dj, dy, mii, mjj, mji, wn;
   int   i, j, n, r;

   for (j = 0; j < nfree; j++) {                /* loop to modify and ... */
      mjj = matrix1[j][j];                      /* scale the matrix */
      if (mjj <= 0.0) return( -5 );             /* diagonal element wrong! */ 
      mjj = sqrt( mjj );
      for (i = 0; i < j; i++) {                 /* scale it */
         mji = matrix1[j][i] / mjj / sqrt( matrix1[i][i] );
         matrix2[i][j] = mji;
	 matrix2[j][i] = mji;
      }
      matrix2[j][j] = 1.0 + labda;              /* scaled value on diagonal */
   }
   if (r = invmat( )) return( r );              /* invert matrix inplace */
   for (i = 0; i < npar; i++) epar[i] = fpar[i];
   for (j = 0; j < nfree; j++) {                /* loop to calculate ... */
      dj = 0.0;                                 /* correction vector */
      mjj = matrix1[j][j];
      if (mjj <= 0.0) return( -7 );             /* not allowed! */
      mjj = sqrt( mjj );
      for (i = 0; i < nfree; i++) {
         mii = matrix1[i][i];
         if (mii <= 0.0) return( -7 );
         mii = sqrt( mii );
         dj += vector[i] * matrix2[j][i] / mjj / mii;
      }
      epar[parptr[j]] += dj;                    /* new parameters */
   }
   chi1 = 0.0;                                  /* reset reduced chi-squared */
   for (n = 0; n < ndat; n++) {                 /* loop through data points */
      wn = wdat ? wdat[n] : 1.0;                /* get weight */
      if (wn > 0.0) {                           /* legal weight */
         dy = ydat[n] - (*fitfunc_c)( &xdat[xdim * n], epar, npar );
         chi1 += wn * dy * dy;
      }
   }
   return( 0 );
} /* getvec */

int nllsqfit(
    real *xdat, 
    int xdim, 
    real *ydat, 
    real *wdat, 
    real *ddat,
    int ndat, 
    real *fpar, 
    real *epar,
    int *mpar, 
    int npar, 
    real tol, 
    int its, 
    real lab, 
    my_proc1 f, 
    my_proc2 df)
{
   int   i, n, r;

   fitfunc_c = f;                       /* save for local routines */
   fitderv_c = df;
   itc = 0;                             /* fate of fit */
   found = 0;                           /* reset */
   nfree = 0;                           /* number of free parameters */
   nuse = 0;                            /* number of legal data points */
   if (tol < (FLT_EPSILON * 10.0)) {
      tolerance = FLT_EPSILON * 10.0;   /* default tolerance */
   } else {
      tolerance = tol;                 /* tolerance */
   }
   labda = fabs( lab ) * LABFAC;       /* start value for mixing parameter */
   for (i = 0; i < npar; i++) {
      epar[i] = 0.0;
      if (mpar[i]) {
         if (nfree > MAXPAR) return( -1 );      /* too many free parameters */
         parptr[nfree++] = i;           /* a free parameter */
      }
   }
   if (nfree == 0) {
     if (labda == 0.0) 
       warning("Not computing differences properly");
     getmat( xdat, xdim, ydat, wdat, ddat, ndat, fpar, epar, npar ); /* get diff */
     return -2;           /* no free parameters */
   }
   for (n = 0; n < ndat; n++) {
     if (wdat && wdat[n] > 0.0) nuse++;        /* legal weight */
     else nuse++;
   }
   if (nfree >= nuse) return( -3 );     /* no degrees of freedom */

   if (labda == 0.0) {                  /* linear fit */

      for (i = 0; i < nfree; fpar[parptr[i++]] = 0.0);
      getmat( xdat, xdim, ydat, wdat, ddat, ndat, fpar, epar, npar );
      r = getvec( xdat, xdim, ydat, wdat, ndat, fpar, epar, npar );
      if (r) return( r );               /* error */
      for (i = 0; i < npar; i++) {
         fpar[i] = epar[i];             /* save new parameters */
         epar[i] = 0.0;                 /* and set errors to zero */
      }
      chi1 = sqrt( chi1 / (real) (nuse - nfree) );
      for (i = 0; i < nfree; i++) {
         if ((matrix1[i][i] <= 0.0) || (matrix2[i][i] <= 0.0)) return( -7 );
         epar[parptr[i]] = chi1 * sqrt( matrix2[i][i] ) / sqrt( matrix1[i][i] );
      }
      /* somehow ddat is not set in linear mode in getmat()..... */
      if (ddat) {
	for (n = 0; n < ndat; n++) {
	  ddat[n] = ydat[n] - (*fitfunc_c)( &xdat[xdim * n], fpar, npar );
	}
      }

   } else {                             /* Non-linear fit */

      /*
       * The non-linear fit uses the steepest descent method in combination
       * with the Taylor method. The mixing of these methods is controlled
       * by labda. In the outer loop (called the iteration loop) we build
       * the matrix and calculate the correction vector. In the inner loop
       * (called the interpolation loop) we check whether we have obtained
       * a better solution than the previous one. If so, we leave the
       * inner loop, else we increase labda (give more weight to the
       * steepest descent method), calculate the correction vector and check
       * again. After the inner loop we do a final check on the goodness of
       * the fit and if this satisfies the tolerance we calculate the
       * errors of the fitted parameters.
       */

      while (!found) {                          /* iteration loop */
         if (itc++ == its) return( -4 );     /* increase iteration counter */
         getmat( xdat, xdim, ydat, wdat, ddat, ndat, fpar, epar, npar );
         /*
          * here we decrease labda since we may assume that each iteration
          * brings us closer to the answer.
          */
         if (labda > LABMIN) labda /= LABFAC;   /* decrease labda */
         r = getvec( xdat, xdim, ydat, wdat, ndat, fpar, epar, npar );
         if (r) return( r );            /* error */
         while (chi1 >= chi2) {         /* interpolation loop */
            /*
             * The next statement is based on experience, not on the
             * mathematics of the problem although I (KGB) think that it
             * is correct to assume that we have reached convergence
             * when the pure steepest descent method does not produce
             * a better solution. Think about this somewhat more, anyway,
             * as already stated, the next statement is based on experience.
             */
            if (labda > LABMAX) break;  /* assume solution found */
            labda *= LABFAC;            /* Increase mixing parameter */
            r = getvec( xdat, xdim, ydat, wdat, ndat, fpar, epar, npar );
            if (r) return( r );         /* error */
         }
         if (labda <= LABMIN) {         /* save old parameters */
            for (i = 0; i < npar; i++) fpar[i] = epar[i];
         }
         if (fabs( chi2 - chi1 ) <= (tolerance * chi1) || (labda > LABMAX)) {
            /*
             * We have a satisfying solution, so now we need to calculate
             * the correct errors of the fitted parameters. This we do
             * by using the pure Taylor method because we are very close
             * to the real solution.
             */
            labda = 0.0;                /* for Taylor solution */
            getmat( xdat, xdim, ydat, wdat, ddat, ndat, fpar, epar, npar );
            r = getvec( xdat, xdim, ydat, wdat, ndat, fpar, epar, npar );
            if (r) return( r );         /* error */
            for (i = 0; i < npar; i++) {
               fpar[i] = epar[i];       /* save new parameters */
               epar[i] = 0.0;           /* and set error to zero */
            }
            chi1 = sqrt( chi1 / (real) (nuse - nfree) );
            for (i = 0; i < nfree; i++) {
               if ((matrix1[i][i] <= 0.0) || (matrix2[i][i] <= 0.0)) return( -7);
#if 1
	       /* original */
               epar[parptr[i]] = chi1 * sqrt( matrix2[i][i] ) / sqrt( matrix1[i][i] );
#else
	       /* somewhat like the nr_ version */
               epar[parptr[i]] = sqrt( matrix2[i][i] ) / sqrt( matrix1[i][i] );
	       if (wdat == NULL)
		 epar[parptr[i]] *= chi1;
#endif
            }
            found = 1;                  /* we found a solution */
         }
      }
   }
   if (ddat) {
     real chisq = 0.0;
     real w;
     for (n = 0; n < ndat; n++) {
       w = wdat ? wdat[i] : 1.0;
       chisq += sqr(ddat[i])*w;
     }
     dprintf(1,"chisq=%g chi1,2=%g %g\n",chisq,chi1,chi2);
   }
   return itc;                       /* return number of iterations (0 for linear) */
}

#if     defined(TESTBED)
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

real func_c( real *xdat, real *fpar, int npar )
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

void derv_c( real *xdat, real *fpar, real *epar, int npar )
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

int main( )
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
   real ddat[30];
   real fpar[4];
   real epar[4];


   mpar[0] = 1; mpar[1] = 0; mpar[2] = 1; mpar[3] = 1;
   for (nfit = 0; nfit < 10; nfit++) {
      fpar[0] = 0.0; fpar[1] = 10.0; fpar[2] = 15.0; fpar[3] = 2.0;
      for (n = 0; n < ndat; n++) {
         real rndm;
         xdat[n] = (real) n;
         wdat[n] = 1.0;
         ydat[n] = func_c( &xdat[n], fpar, npar );
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
      r = nllsqfit( xdat, xdim, ydat, wdat, ddat, ndat, fpar, epar, mpar, npar, tol, 
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
   return 0;
}
#endif
