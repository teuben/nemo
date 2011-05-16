#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include <limits.h>
#include <memory.h>
#include "dopri5.h"


static long      nfcn, nstep, naccpt, nrejct;
static double    hout, xold, xout;
static unsigned  nrds, *indir;
static double    *yy1, *k1, *k2, *k3, *k4, *k5, *k6, *ysti;
static double    *rcont1, *rcont2, *rcont3, *rcont4, *rcont5;


long nfcnRead (void)
{
  return nfcn;

} /* nfcnRead */


long nstepRead (void)
{
  return nstep;

} /* stepRead */


long naccptRead (void)
{
  return naccpt;

} /* naccptRead */


long nrejctRead (void)
{
  return nrejct;

} /* nrejct */


double hRead (void)
{
  return hout;

} /* hRead */


double xRead (void)
{
  return xout;

} /* xRead */


static double sign (double a, double b)
{
  return (b > 0.0) ? fabs(a) : -fabs(a);

} /* sign */


static double min_d (double a, double b)
{
  return (a < b)?a:b;

} /* min_d */


static double max_d (double a, double b)
{
  return (a > b)?a:b;

} /* max_d */


static double hinit (unsigned n, FcnEqDiff fcn, double x, double* y,
	      double posneg, double* f0, double* f1, double* yy1, int iord,
	      double hmax, double* atoler, double* rtoler, int itoler)
{
  double   dnf, dny, atoli, rtoli, sk, h, h1, der2, der12, sqr;
  unsigned i;

  dnf = 0.0;
  dny = 0.0;
  atoli = atoler[0];
  rtoli = rtoler[0];

  if (!itoler)
    for (i = 0; i < n; i++)
    {
      sk = atoli + rtoli * fabs(y[i]);
      sqr = f0[i] / sk;
      dnf += sqr*sqr;
      sqr = y[i] / sk;
      dny += sqr*sqr;
    }
  else
    for (i = 0; i < n; i++)
    {
      sk = atoler[i] + rtoler[i] * fabs(y[i]);
      sqr = f0[i] / sk;
      dnf += sqr*sqr;
      sqr = y[i] / sk;
      dny += sqr*sqr;
    }

  if ((dnf <= 1.0E-10) || (dny <= 1.0E-10))
    h = 1.0E-6;
  else
    h = sqrt (dny/dnf) * 0.01;

  h = min_d (h, hmax);
  h = sign (h, posneg);

  /* perform an explicit Euler step */
  for (i = 0; i < n; i++)
    yy1[i] = y[i] + h * f0[i];
  fcn (n, x+h, yy1, f1);

  /* estimate the second derivative of the solution */
  der2 = 0.0;
  if (!itoler)
    for (i = 0; i < n; i++)
    {
      sk = atoli + rtoli * fabs(y[i]);
      sqr = (f1[i] - f0[i]) / sk;
      der2 += sqr*sqr;
    }
  else
    for (i = 0; i < n; i++)
    {
      sk = atoler[i] + rtoler[i] * fabs(y[i]);
      sqr = (f1[i] - f0[i]) / sk;
      der2 += sqr*sqr;
    }
  der2 = sqrt (der2) / h;

  /* step size is computed such that h**iord * max_d(norm(f0),norm(der2)) = 0.01 */
  der12 = max_d (fabs(der2), sqrt(dnf));
  if (der12 <= 1.0E-15)
    h1 = max_d (1.0E-6, fabs(h)*1.0E-3);
  else
    h1 = pow (0.01/der12, 1.0/(double)iord);
  h = min_d (100.0 * h, min_d (h1, hmax));

  return sign (h, posneg);

} /* hinit */


/* core integrator */
static int dopcor (unsigned n, FcnEqDiff fcn, double x, double* y, double xend,
		   double hmax, double h, double* rtoler, double* atoler,
		   int itoler, FILE* fileout, SolTrait solout, int iout,
		   long nmax, double uround, int meth, long nstiff, double safe,
		   double beta, double fac1, double fac2, unsigned* icont)
{
  double   facold, expo1, fac, facc1, facc2, fac11, posneg, xph;
  double   atoli, rtoli, hlamb, err, sk, hnew, yd0, ydiff, bspl;
  double   stnum, stden, sqr;
  int      iasti, iord, irtrn, reject, last, nonsti;
  unsigned i, j;
  double   c2, c3, c4, c5, e1, e3, e4, e5, e6, e7, d1, d3, d4, d5, d6, d7;
  double   a21, a31, a32, a41, a42, a43, a51, a52, a53, a54;
  double   a61, a62, a63, a64, a65, a71, a73, a74, a75, a76;

  /* initialisations */
  switch (meth)
  {
    case 1:

      c2=0.2, c3=0.3, c4=0.8, c5=8.0/9.0;
      a21=0.2, a31=3.0/40.0, a32=9.0/40.0;
      a41=44.0/45.0, a42=-56.0/15.0; a43=32.0/9.0;
      a51=19372.0/6561.0, a52=-25360.0/2187.0;
      a53=64448.0/6561.0, a54=-212.0/729.0;
      a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0;
      a64=49.0/176.0, a65=-5103.0/18656.0;
      a71=35.0/384.0, a73=500.0/1113.0, a74=125.0/192.0;
      a75=-2187.0/6784.0, a76=11.0/84.0;
      e1=71.0/57600.0, e3=-71.0/16695.0, e4=71.0/1920.0;
      e5=-17253.0/339200.0, e6=22.0/525.0, e7=-1.0/40.0;
      d1=-12715105075.0/11282082432.0, d3=87487479700.0/32700410799.0;
      d4=-10690763975.0/1880347072.0, d5=701980252875.0/199316789632.0;
      d6=-1453857185.0/822651844.0, d7=69997945.0/29380423.0;

      break;
  }

  facold = 1.0E-4;
  expo1 = 0.2 - beta * 0.75;
  facc1 = 1.0 / fac1;
  facc2 = 1.0 / fac2;
  posneg = sign (1.0, xend-x);

  /* initial preparations */
  atoli = atoler[0];
  rtoli = rtoler[0];
  last  = 0;
  hlamb = 0.0;
  iasti = 0;
  fcn (n, x, y, k1);
  hmax = fabs (hmax);
  iord = 5;
  if (h == 0.0)
    h = hinit (n, fcn, x, y, posneg, k1, k2, k3, iord, hmax, atoler, rtoler, itoler);
  nfcn += 2;
  reject = 0;
  xold = x;
  if (iout)
  {
    irtrn = 1;
    hout = h;
    xout = x;
    solout (naccpt+1, xold, x, y, n, &irtrn);
    if (irtrn < 0)
    {
      if (fileout)
	fprintf (fileout, "Exit of dopri5 at x = %.16e\r\n", x);
      return 2;
    }
  }

  /* basic integration step */
  while (1)
  {
    if (nstep > nmax)
    {
      if (fileout)
	fprintf (fileout, "Exit of dopri5 at x = %.16e, more than nmax = %li are needed\r\n", x, nmax);
      xout = x;
      hout = h;
      return -2;
    }

    if (0.1 * fabs(h) <= fabs(x) * uround)
    {
      if (fileout)
	fprintf (fileout, "Exit of dopri5 at x = %.16e, step size too small h = %.16e\r\n", x, h);
      xout = x;
      hout = h;
      return -3;
    }

    if ((x + 1.01*h - xend) * posneg > 0.0)
    {
      h = xend - x;
      last = 1;
    }

    nstep++;

    /* the first 6 stages */
    for (i = 0; i < n; i++)
      yy1[i] = y[i] + h * a21 * k1[i];
    fcn (n, x+c2*h, yy1, k2);
    for (i = 0; i < n; i++)
      yy1[i] = y[i] + h * (a31*k1[i] + a32*k2[i]);
    fcn (n, x+c3*h, yy1, k3);
    for (i = 0; i < n; i++)
      yy1[i] = y[i] + h * (a41*k1[i] + a42*k2[i] + a43*k3[i]);
    fcn (n, x+c4*h, yy1, k4);
    for (i = 0; i <n; i++)
      yy1[i] = y[i] + h * (a51*k1[i] + a52*k2[i] + a53*k3[i] + a54*k4[i]);
    fcn (n, x+c5*h, yy1, k5);
    for (i = 0; i < n; i++)
      ysti[i] = y[i] + h * (a61*k1[i] + a62*k2[i] + a63*k3[i] + a64*k4[i] + a65*k5[i]);
    xph = x + h;
    fcn (n, xph, ysti, k6);
    for (i = 0; i < n; i++)
      yy1[i] = y[i] + h * (a71*k1[i] + a73*k3[i] + a74*k4[i] + a75*k5[i] + a76*k6[i]);
    fcn (n, xph, yy1, k2);
    if (iout == 2)
      if (nrds == n)
	for (i = 0; i < n; i++)
	{
	  rcont5[i] = h * (d1*k1[i] + d3*k3[i] + d4*k4[i] + d5*k5[i] + d6*k6[i] + d7*k2[i]);
	}
      else
	for (j = 0; j < nrds; j++)
	{
	  i = icont[j];
	  rcont5[j] = h * (d1*k1[i] + d3*k3[i] + d4*k4[i] + d5*k5[i] + d6*k6[i] + d7*k2[i]);
	}
    for (i = 0; i < n; i++)
      k4[i] = h * (e1*k1[i] + e3*k3[i] + e4*k4[i] + e5*k5[i] + e6*k6[i] + e7*k2[i]);
    nfcn += 6;

    /* error estimation */
    err = 0.0;
    if (!itoler)
      for (i = 0; i < n; i++)
      {
	sk = atoli + rtoli * max_d (fabs(y[i]), fabs(yy1[i]));
	sqr = k4[i] / sk;
	err += sqr*sqr;
      }
    else
      for (i = 0; i < n; i++)
      {
	sk = atoler[i] + rtoler[i] * max_d (fabs(y[i]), fabs(yy1[i]));
	sqr = k4[i] / sk;
	err += sqr*sqr;
      }
    err = sqrt (err / (double)n);

    /* computation of hnew */
    fac11 = pow (err, expo1);
    /* Lund-stabilization */
    fac = fac11 / pow(facold,beta);
    /* we require fac1 <= hnew/h <= fac2 */
    fac = max_d (facc2, min_d (facc1, fac/safe));
    hnew = h / fac;

    if (err <= 1.0)
    {
      /* step accepted */

      facold = max_d (err, 1.0E-4);
      naccpt++;

      /* stiffness detection */
      if (!(naccpt % nstiff) || (iasti > 0))
      {
	stnum = 0.0;
	stden = 0.0;
	for (i = 0; i < n; i++)
	{
	  sqr = k2[i] - k6[i];
	  stnum += sqr*sqr;
	  sqr = yy1[i] - ysti[i];
	  stden += sqr*sqr;
	}
	if (stden > 0.0)
	  hlamb = h * sqrt (stnum / stden);
	if (hlamb > 3.25)
	{
	  nonsti = 0;
	  iasti++;
	  if (iasti == 15)
	    if (fileout)
	      fprintf (fileout, "The problem seems to become stiff at x = %.16e\r\n", x);
	    else
	    {
	      xout = x;
	      hout = h;
	      return -4;
	    }
	}
	else
	{
	  nonsti++;
	  if (nonsti == 6)
	    iasti = 0;
	}
      }

      if (iout == 2)
	if (nrds == n)
	  for (i = 0; i < n; i++)
	  {
	    yd0 = y[i];
	    ydiff = yy1[i] - yd0;
	    bspl = h * k1[i] - ydiff;
	    rcont1[i] = y[i];
	    rcont2[i] = ydiff;
	    rcont3[i] = bspl;
	    rcont4[i] = -h * k2[i] + ydiff - bspl;
	  }
	else
	  for (j = 0; j < nrds; j++)
	  {
	    i = icont[j];
	    yd0 = y[i];
	    ydiff = yy1[i] - yd0;
	    bspl = h * k1[i] - ydiff;
	    rcont1[j] = y[i];
	    rcont2[j] = ydiff;
	    rcont3[j] = bspl;
	    rcont4[j] = -h * k2[i] + ydiff - bspl;
	  }

      memcpy (k1, k2, n * sizeof(double)); 
      memcpy (y, yy1, n * sizeof(double));
      xold = x;
      x = xph;

      if (iout)
      {
	hout = h;
	xout = x;
	solout (naccpt+1, xold, x, y, n, &irtrn);
	if (irtrn < 0)
	{
	  if (fileout)
	    fprintf (fileout, "Exit of dopri5 at x = %.16e\r\n", x);
	  return 2;
	}
      }

      /* normal exit */
      if (last)
      {
	hout=hnew;
	xout = x;
	return 1;
      }

      if (fabs(hnew) > hmax)
	hnew = posneg * hmax;
      if (reject)
	hnew = posneg * min_d (fabs(hnew), fabs(h));

      reject = 0;
    }
    else
    {
      /* step rejected */
      hnew = h / min_d (facc1, fac11/safe);
      reject = 1;
      if (naccpt >= 1)
	nrejct=nrejct + 1;
      last = 0;
    }

    h = hnew;
  }

} /* dopcor */


/* front-end */
int dopri5
 (unsigned n, FcnEqDiff fcn, double x, double* y, double xend, double* rtoler,
  double* atoler, int itoler, SolTrait solout, int iout, FILE* fileout, double uround,
  double safe, double fac1, double fac2, double beta, double hmax, double h,
  long nmax, int meth, long nstiff, unsigned nrdens, unsigned* icont, unsigned licont)
{
  int       arret, idid;
  unsigned  i;

  /* initialisations */
  nfcn = nstep = naccpt = nrejct = arret = 0;
  rcont1 = rcont2 = rcont3 = rcont4 = rcont5 = NULL;
  indir = NULL;

  /* n, the dimension of the system */
  if (n == UINT_MAX)
  {
    if (fileout)
      fprintf (fileout, "System too big, max. n = %u\r\n", UINT_MAX-1);
    arret = 1;
  }

  /* nmax, the maximal number of steps */
  if (!nmax)
    nmax = 100000;
  else if (nmax <= 0)
  {
    if (fileout)
      fprintf (fileout, "Wrong input, nmax = %li\r\n", nmax);
    arret = 1;
  }

  /* meth, coefficients of the method */
  if (!meth)
    meth = 1;
  else if ((meth <= 0) || (meth >= 2))
  {
    if (fileout)
      fprintf (fileout, "Curious input, meth = %i\r\n", meth);
    arret = 1;
  }

  /* nstiff, parameter for stiffness detection */
  if (!nstiff)
    nstiff = 1000;
  else if (nstiff < 0)
    nstiff = nmax + 10;

  /* iout, switch for calling solout */
  if ((iout < 0) || (iout > 2))
  {
    if (fileout)
      fprintf (fileout, "Wrong input, iout = %i\r\n", iout);
    arret = 1;
  }

  /* nrdens, number of dense output components */
  if (nrdens > n)
  {
    if (fileout)
      fprintf (fileout, "Curious input, nrdens = %u\r\n", nrdens);
    arret = 1;
  }
  else if (nrdens)
  {
    /* is there enough memory to allocate rcont12345&indir ? */
    rcont1 = (double*) malloc (nrdens*sizeof(double));
    rcont2 = (double*) malloc (nrdens*sizeof(double));
    rcont3 = (double*) malloc (nrdens*sizeof(double));
    rcont4 = (double*) malloc (nrdens*sizeof(double));
    rcont5 = (double*) malloc (nrdens*sizeof(double));
    if (nrdens < n)
      indir = (unsigned*) malloc (n*sizeof(unsigned));

    if (!rcont1 || !rcont2 || !rcont3 || !rcont4 || !rcont5 || (!indir && (nrdens < n)))
    {
      if (fileout)
	fprintf (fileout, "Not enough free memory for rcont12345&indir\r\n");
      arret = 1;
    }

    /* control of length of icont */
    if (nrdens == n)
    {
      if (icont && fileout)
	fprintf (fileout, "Warning : when nrdens = n there is no need allocating memory for icont\r\n");
      nrds = n;
    }
    else if (licont < nrdens)
    {
      if (fileout)
	fprintf (fileout, "Insufficient storage for icont, min. licont = %u\r\n", nrdens);
      arret = 1;
    }
    else
    {
      if ((iout < 2) && fileout)
	fprintf (fileout, "Warning : put iout = 2 for dense output\r\n");
      nrds = nrdens;
      for (i = 0; i < n; i++)
	indir[i] = UINT_MAX;
      for (i = 0; i < nrdens; i++)
	indir[icont[i]] = i;
    }
  }

  /* uround, smallest number satisfying 1.0+uround > 1.0 */
  if (uround == 0.0)
    uround = 2.3E-16;
  else if ((uround <= 1.0E-35) || (uround >= 1.0))
  {
    if (fileout)
      fprintf (fileout, "Which machine do you have ? Your uround was : %.16e\r\n", uround);
    arret = 1;
  }

  /* safety factor */
  if (safe == 0.0)
    safe = 0.9;
  else if ((safe >= 1.0) || (safe <= 1.0E-4))
  {
    if (fileout)
      fprintf (fileout, "Curious input for safety factor, safe = %.16e\r\n", safe);
    arret = 1;
  }

  /* fac1, fac2, parameters for step size selection */
  if (fac1 == 0.0)
    fac1 = 0.2;
  if (fac2 == 0.0)
    fac2 = 10.0;

  /* beta for step control stabilization */
  if (beta == 0.0)
    beta = 0.04;
  else if (beta < 0.0)
    beta = 0.0;
  else if (beta > 0.2)
  {
    if (fileout)
      fprintf (fileout, "Curious input for beta : beta = %.16e\r\n", beta);
    arret = 1;
  }

  /* maximal step size */
  if (hmax == 0.0)
    hmax = xend - x;

  /* is there enough free memory for the method ? */
  yy1 = (double*) malloc (n*sizeof(double));
  k1 = (double*) malloc (n*sizeof(double));
  k2 = (double*) malloc (n*sizeof(double));
  k3 = (double*) malloc (n*sizeof(double));
  k4 = (double*) malloc (n*sizeof(double));
  k5 = (double*) malloc (n*sizeof(double));
  k6 = (double*) malloc (n*sizeof(double));
  ysti = (double*) malloc (n*sizeof(double));

  if (!yy1 || !k1 || !k2 || !k3 || !k4 || !k5 || !k6 || !ysti)
  {
    if (fileout)
      fprintf (fileout, "Not enough free memory for the method\r\n");
    arret = 1;
  }

  /* when a failure has occured, we return -1 */
  if (arret)
  {
    if (ysti)
      free (ysti);
    if (k6)
      free (k6);
    if (k5)
      free (k5);
    if (k4)
      free (k4);
    if (k3)
      free (k3);
    if (k2)
      free (k2);
    if (k1)
      free (k1);
    if (yy1)
      free (yy1);
    if (indir)
      free (indir);
    if (rcont5)
      free (rcont5);
    if (rcont4)
      free (rcont4);
    if (rcont3)
      free (rcont3);
    if (rcont2)
      free (rcont2);
    if (rcont1)
      free (rcont1);

    return -1;
  }
  else
  {
    idid = dopcor (n, fcn, x, y, xend, hmax, h, rtoler, atoler, itoler, fileout,
		   solout, iout, nmax, uround, meth, nstiff, safe, beta, fac1, fac2, icont);
    free (ysti);
    free (k6);
    free (k5);    /* reverse order freeing too increase chances */
    free (k4);    /* of efficient dynamic memory managing       */
    free (k3);
    free (k2);
    free (k1);
    free (yy1);
    if (indir)
      free (indir);
    if (rcont5)
    {
      free (rcont5);
      free (rcont4);
      free (rcont3);
      free (rcont2);
      free (rcont1);
    }

    return idid;
  }

} /* dopri5 */


/* dense output function */
double contd5 (unsigned ii, double x)
{
  unsigned i, j;
  double   theta, theta1;

  i = UINT_MAX;

  if (!indir)
    i = ii;
  else
    i = indir[ii];

  if (i == UINT_MAX)
  {
    printf ("No dense output available for %uth component", ii);
    return 0.0;
  }

  theta = (x - xold) / hout;
  theta1 = 1.0 - theta;

  return rcont1[i] + theta*(rcont2[i] + theta1*(rcont3[i] + theta*(rcont4[i] + theta1*rcont5[i])));

} /* contd5 */

