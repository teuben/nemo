/*WIP: adding mode=plane */
/*TODO:adding mode=poly*/
/*BUG: large tables will fail in old table input method */
/*
 * TABLSQFIT: a general (linear) fitting program for tabular data 
 *
 *       5-aug-87  created using NumRec in double C     Peter Teuben
 *      10-nov-87  V1.2 need new 'drange' version skipping blanks
 *      11-Mar-88  V1.3 dprintf() output now used
 *      17-May-88  V1.4 optional output of data+residuals
 *      31-May-88  V2.0 for NumRec beta-release (float instead of double)
 *                      fixed small buf for interactive input (see also tabmath)
 *      30-aug-90  V2.1 fitting ellipses (still unfinished)
 *      13-nov-90  V2.2 drange->nemoinp; using lsq() not NumRec
 *       9-apr-91  V2.3 optional output with deleting high residuals data
 *      15-nov-91  V2.4 fixed ellipse fit routine (PJT)
 *	21-nov-91       optional plot ellipse and or residuals - hyperbola
 *	12-mar-92  V2.5 added the imageshift fit (complex y=ax+b)
 *	13-apr-92  V2.5a  fixed bad nemoinpi("x...   calls
 *	 6-nov-93  V2.6 added nmax= to allow dynamic input
 *		       a  handle no-data cases
 *      12-jun-94  V2.7 general overhaul to use a 'column' data structure
 *                      to make code extendible for different fit modes
 *	 8-jun-95  V2.8 added fit option peak
 *
 *	30-jan-97     a resurrected the xrange selection for nxcol==1 !!
 *	26-jan-98  V2.9 added gamma functions and goodness of fit for fit=line
 *       4-feb-98     a compute 'r' (corr.coeff.) for fit=line
 *      14-feb-98     c fixed bug in poly fitting column requirement
 *	17-apr-00     d fixed edge finding in do_peak()
 *      14-aug-00  V3.0 added tab=t|f for tabular short output
 *      23-mar-01     a fixed xrange=  for singular column files
 *      14-apr-01  V3.1 added convex hull computation as 'area'
 *      18-jun-01     b fixed boundary check for xrange, also use natof now in setxrange
 *       8-aug-01     c compute error in axis ratio for ellipses (w/ Mousumi Das)
 *      10-sep-01  V3.2 GSL enabled for linear fit
 *                      figuring out error bars?
 *      29-oct-02  V3.2c: ellipse fit cleanup (still bug in ellipse semi major axis?)
 *                 V3.3: add fourier mode (see also snapfour)
 *                     a: add out= for fourier
 *      24-feb-03  V3.4  add fit=zero
 *       4-oct-03      a fix nsigma>0 for fit=line
 *       3-may-05  V3.4b add x/x0+y/y0=1 variant for a linear fit
 *      21-nov-05  V3.4c added gauss2d
 *      16-feb-13  V3.5  added fit=slope from miriad::immerge
 *      28-may-13   4.0e fixed bug in fit=peak value
 *
 * TODO:   check 'r', wip gives slightly different numbers
 */

/*
 * literature:
 *
   W.Gander, G.H.Golub and RStrebel, "Least-Square Fitting of Circles and Ellipses",
    BIT, No,43 pp.558-578, 1998
or
    A.Fitzbibbon, M. Pilu and R.B.Fisher, "Direct Least Square Fitting of Ellipse",
    IEEE Tran. Pattern Analysis and Maschine Intelligence, Vol.21, No.5, pp.476-480, May 1999

    http://www.mai.liu.se/~akbjo/NMbook.html
	ch.11 has a section on fitting circles and ellipses

 */

#include <stdinc.h>  
#include <getparam.h>


/* pick (n)one TESTNR= numrec    TESTMP = mpfit */
//#define TESTNR 1
#define TESTMP 1


#ifdef TESTNR
#include "testnr.h"
#endif


#ifdef TESTMP
#include <mpfit.h>
#endif


/*    we don't allow this code with GSL yet, output is
 *    not unified with older version - compile tablsqfit_gsl
 */

#if HAVE_GSL
#include <gsl/gsl_fit.h>
#endif 

string defv[] = {
    "in=???\n           input (table) file name",
    "xcol=1\n           column(s) for x",
    "ycol=2\n           column(s) for y",
    "dxcol=\n           column for X errors (only used in fit=line)",
    "dycol=\n           column for Y errors",
    "xrange=\n          in case restricted range is used (1D only)",
    "fit=line\n         fitmode (line, ellipse, imageshift, plane, poly, peak, area, zero, gauss2d, slope)",
    "order=0\n		Order of plane/poly fit",
    "out=\n             optional output file for some fit modes",
    "nsigma=-1\n        delete points more than nsigma away?",
    "estimate=\n        optional estimates (e.g. for ellipse center)",
    "nmax=10000\n       Default max allocation",
    "mpfit=0\n          fit mode for mpfit",
    "tab=f\n            short one-line output?",
    "VERSION=4.0e\n     28-may-2013 PJT",
    NULL
};

string usage="a linear least square fitting program";

string cvsid="$Id$";


/**************** SOME GLOBAL VARIABLES ************************/

#if !defined(HUGE)
#define HUGE 1e20
#endif

#ifndef MAXCOL
#define MAXCOL 16
#endif

typedef struct column {
    int maxdat;     /* allocated length of data */          /* not used */
    int ndat;       /* actual length of data */             /* not used */
    real *dat;      /* pointer to data */
    int colnr;      /* column number this data came from */ /* not used */
} a_column;

int nxcol, nycol, xcolnr[MAXCOL], ycolnr[MAXCOL], dxcolnr, dycolnr; 
a_column            xcol[MAXCOL],   ycol[MAXCOL],   dxcol, dycol;

real xrange[MAXCOL*2];      /* ??? */

string method;              /* fit method (line, ellipse, ....) */
stream instr, outstr;       /* input / output file */


int    nmax;                /* allocated space */
int    npt;                 /* actual number of points from table */
real   nsigma;              /* fractional sigma removal */

real  a,b;                  /* fit parameters in: y=ax+b  */
int order;

int mpfit_mode;

bool Qtab;                  /* do table output ? */


/****************************** START OF PROGRAM **********************/

nemo_main()
{

    setparams();
    read_data();

    if (scanopt(method,"line")) {
        do_line();
    } else if (scanopt(method,"slope")) {
        do_slope();
    } else if (scanopt(method,"ellipse")) {
        do_ellipse();
    } else if (scanopt(method,"imageshift")) {
        do_imageshift();
    } else if (scanopt(method,"plane")) {
    	do_plane();
    } else if (scanopt(method,"gauss1d")) {
    	do_gauss1d();
    } else if (scanopt(method,"gauss2d")) {
    	do_gauss2d();
    } else if (scanopt(method,"poly")) {
    	do_poly();
    } else if (scanopt(method,"area")) {
        do_area();
    } else if (scanopt(method,"peak")) {
    	do_peak();
    } else if (scanopt(method,"zero")) {
    	do_zero();
    } else if (scanopt(method,"fourier")) {
    	do_fourier();
    } else
        error("fit=%s invalid; try [line,ellipse,imageshift,plane,gauss1d,gauss2d,poly,area,peak,zero,fourier]",
	      getparam("fit"));

    if (outstr) strclose(outstr);
}

setparams()
{
    string inname = getparam("in");
    nmax = nemo_file_lines(inname,getiparam("nmax"));
    if (nmax<0) error("Error opening %s",inname);
    if (nmax==0) error("No data?");
    instr = stropen (inname,"r");

    if (hasvalue("out"))
        outstr=stropen(getparam("out"),"w");
    else
        outstr=NULL;

    nxcol = nemoinpi(getparam("xcol"),xcolnr,MAXCOL);
    if (nxcol<0) error("Illegal xcol= nxcol=%d",nxcol);
    nycol = nemoinpi(getparam("ycol"),ycolnr,MAXCOL);
    if (nycol<0) error("Illegal ycol= nycol=%d",nycol);

    if (hasvalue("dxcol"))
        dxcolnr = getiparam("dxcol");
    else
        dxcolnr = 0;

    if (hasvalue("dycol"))
        dycolnr = getiparam("dycol");
    else
        dycolnr = 0;

    if (hasvalue("xrange"))
        setrange(xrange,getparam("xrange"));
    else {
        xrange[0] = -HUGE;
        xrange[1] = HUGE;
    } 
    
    method = getparam("fit");
    nsigma = getdparam("nsigma");
    order = getiparam("order");
    if (order<0) error("order=%d of %s cannot be negative",order,method);
    Qtab = getbparam("tab");

    mpfit_mode = getiparam("mpfit");
}

setrange(real *rval, string rexp)
{
    char *cptr;

    cptr = strchr(rexp, ':');
    if (cptr) {
        rval[0] = natof(rexp);
        rval[1] = natof(cptr+1);
    } else {
        rval[0] = 0.0;
        rval[1] = natof(rexp);
    	warning("Range taken from 0 - %g",rval[1]);
    }
}

read_data()
{
    real *coldat[2*MAXCOL+1];
    int colnr[2*MAXCOL+1], ncols = 0, i, j;

    dprintf(0,"%s: reading X column(s) %s and Y column(s) %s\n",
	    getparam("in"),getparam("xcol"),getparam("ycol"));

    for (i=0; i<nxcol; i++) {
        coldat[ncols] = xcol[i].dat = (real *) allocate(nmax * sizeof(real));
        colnr[ncols] = xcolnr[i];        
        ncols++;
    }
    for (i=0; i<nycol; i++) {
        coldat[ncols] = ycol[i].dat = (real *) allocate(nmax * sizeof(real));
        colnr[ncols] = ycolnr[i];        
        ncols++;
    }
    if (dxcolnr>0) {
        coldat[ncols] = dxcol.dat = (real *) allocate(nmax * sizeof(real));
        colnr[ncols] = dxcolnr;
        ncols++;
    }
    if (dycolnr>0) {
        coldat[ncols] = dycol.dat = (real *) allocate(nmax * sizeof(real));
        colnr[ncols] = dycolnr;
        ncols++;
    }
    
    npt = get_atable(instr,ncols,colnr,coldat,nmax);
    if (npt < 0) {
        npt = -npt;
       	warning("Could only read %d data",npt);
    }


    /* special case for nxcol=1  ... what to do for nxcol > 1 ??? */
    /* should also handle nycol > 1  but does not yet             */

    if (nxcol == 1 && nycol == 1) {
        for(i=0, j=0; i<npt; i++) {
          if(xrange[0] <= xcol[0].dat[i] && xcol[0].dat[i] <= xrange[1]) {    /* sub-select on X */
              xcol[0].dat[j] = xcol[0].dat[i];
              ycol[0].dat[j] = ycol[0].dat[i];
              if (dxcolnr>0) dxcol.dat[j] = dxcol.dat[i];
              if (dycolnr>0) dycol.dat[j] = dycol.dat[i];
              j++;
           }
        }
        dprintf(1,"Copied over %d/%d data within xrange's\n",j,npt);
	npt = j;
    }
       
    if (npt==0) error("No data");
}


write_data(stream outstr)              /* only for straight line fit */
{
#if 0    
    int i;
    real d;

    for (i=0; i<npt; i++) {
        d = y[i] - a*x[i] - b;
        fprintf (outstr,"%f %f %f\n",x[i],y[i],d);
    }
#endif    
}

/* helper stuff for CMPFIT */
struct vars_struct {
  double *x;
  double *y;
  double *ex;
  double *ey;
  int mode;
};

int linfitex(int m, int n, double *p, double *dy, double **dvec, void *vars)
{
  int i, mode;
  struct vars_struct *v = (struct vars_struct *) vars;
  double *x, *y, *ex, *ey, f;
  static int count=0;

  x = v->x;
  y = v->y;
  ex = v->ex;
  ey = v->ey;
  mode = v->mode;

  if (count==0) {
    count++;
    for (i=0; i<m; i++) {
      dprintf(1,"%d : %g %g %g %g\n",i,x[i],y[i],ex[i],ey[i]);
    }
  }

  for (i=0; i<m; i++) {
    f = p[0] + p[1]*x[i];     /* Linear fit function */
    if (mode==0)
      dy[i] = (y[i] - f);
    else if (mode==1)
      dy[i] = (y[i] - f)/(p[1]*ex[i]);
    else if (mode==2)
      dy[i] = (y[i] - f)/ey[i];
    else
      dy[i] = (y[i] - f)/sqrt(ey[i]*ey[i] + p[1]*p[1]*ex[i]*ex[i]);
  }
}


do_line()
{
    real *x, *y, *dx, *dy, *dz;
    int i,j, mwt;
    real chi2,q,siga,sigb, sigma, d, sa, sb;
    real cov00, cov01, cov11, sumsq;
    real r, prob, z;
    void fit(), pearsn();

    if (nxcol < 1) error("nxcol=%d",nxcol);
    if (nycol < 1) error("nycol=%d",nycol);
    x = xcol[0].dat;
    y = ycol[0].dat;
    dx = (dxcolnr>0 ? dxcol.dat : NULL);
    dy = (dycolnr>0 ? dycol.dat : NULL);

#if HAVE_GSL
    if (dx) error("Cannot use GSL with errors in X column");
    gsl_fit_linear (x, 1, y, 1, npt, &b, &a, 
		    &cov00, &cov01, &cov11, &sumsq);
    printf("y=ax+b:  a=%g b=%g  cov00,01,11=%g %g %g   sumsq=%g\n",
	     a,b, cov00,cov01,cov11, sumsq);

    if (dy) {
      for (i=0; i<npt; i++)
	dy[i] = 1/(dy[i]*dy[i]);
      gsl_fit_wlinear (x, 1, dy, 1, y, 1, npt, &b, &a, 
		    &cov00, &cov01, &cov11, &chi2);
      printf("y=ax+b:  a=%g b=%g  cov00,01,11=%g %g %g   chi^2=%g\n",
	     a,b, cov00,cov01,cov11, chi2);

    }
#else

    if (dx) {
#if defined(TESTNR)
      warning("new FITEXY method");
      dz = (real *) allocate(npt*sizeof(real));
      for (i=0; i<npt; i++) dz[i] = 0.0;

      fitexy(x,y,npt,dx,dy,&b,&a,&sigb,&siga,&chi2,&q);
      printf("fitexy(x,y,dx,dy)    a=%g   b=%g  %g %g\n",b,a,sigb,siga);
      fitexy(x,y,npt,dz,dy,&a,&b,&siga,&sigb,&chi2,&q);
      printf("fitexy(x,y, 0,dy)    a=%g   b=%g  %g %g\n",a,b,siga,sigb);
      fitexy(x,y,npt,dx,dz,&a,&b,&siga,&sigb,&chi2,&q);
      printf("fitexy(x,y,dx, 0)    a=%g   b=%g  %g %g\n",a,b,siga,sigb);
      fit   (x,y,npt,dy, 1,&a,&b,&siga,&sigb,&chi2,&q);
      printf("fit   (x,y, 0,dy)    a=%g   b=%g  %g %g\n",a,b,siga,sigb);
      fit   (y,x,npt,dx, 1,&a,&b,&siga,&sigb,&chi2,&q);
      printf("fit   (y,x, 0,dx)    a=%g   b=%g  %g %g\n",a,b,siga,sigb);

      sa=sqrt(siga*siga+sqr(sigb*(a/b)))/b;
      sb=sigb/(b*b);
      printf("FITEXY: %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f\n",
	     -a/b,1.0/b,sa,sb,chi2,q); 
#elif defined(TESTMP)
      struct vars_struct v;      /* private data for mpfit() */
      int status;
      mp_result result;
      double p[2] = {1.0, 1.0};
      double perror[2];

      warning("MPFIT method; mode mpfit=%d",mpfit_mode);

      bzero(&result,sizeof(result));
      result.xerror = perror;

      v.x = xcol[0].dat;
      v.y = ycol[0].dat;
      v.ex = dxcol.dat;
      v.ey = dycol.dat;
      v.mode = mpfit_mode;
      
      status = mpfit(linfitex, npt, 2, p, 0, 0, (void *) &v, &result);
      dprintf(1,"*** mpfit status = %d\n", status);

      printf("  CHI-SQUARE = %f    (%d DOF)\n", 
	     result.bestnorm, result.nfunc-result.nfree);
      printf("        NPAR = %d\n", result.npar);
      printf("       NFREE = %d\n", result.nfree);
      printf("     NPEGGED = %d\n", result.npegged);
      printf("       NITER = %d\n", result.niter);
      printf("        NFEV = %d\n", result.nfev);
      printf("\n");
      for (i=0; i<result.npar; i++) {
	printf("  P[%d] = %f +/- %f\n", 
	       i, p[i], result.xerror[i]);
      }
#else
      error("no dycol= implemented yet");
#endif
    } else {

      for(mwt=0;mwt<=1;mwt++) {

#if 1
	if (mwt>0 && nsigma>0) {
	  fit(x,y,npt,dy,0,&b,&a,&sigb,&siga,&chi2,&q);
	} else {
	  if (mwt>0 && dycolnr==0)
	    continue;             
	  fit(x,y,npt,dy,mwt,&b,&a,&sigb,&siga,&chi2,&q);
	}
#else
	if (mwt==0)
	  fit(x,y,npt,dy,mwt,&b,&a,&sigb,&siga,&chi2,&q);
	else
	  myfit(x,y,npt,dy,mwt,&b,&a,&sigb,&siga,&chi2,&q);
#endif
	dprintf(2,"\n");
	dprintf (1,"Fit:   y=ax+b\n");
	if (mwt == 0)
	  dprintf(1,"Ignoring standard deviations\n");
	else
	  dprintf(1,"Including standard deviation\n");
	printf("%12s %9.6f %18s %9.6f \n",
	       "a  =  ",a,"uncertainty:",siga);
	printf("%12s %9.6f %18s %9.6f \n",
	       "b  =  ",b,"uncertainty:",sigb);

	printf("%12s %9.6f %18s %9.6f \n",
	       "x0 =  ",-b/a,"uncertainty:",sqrt(sqr(sigb/a)+sqr(siga*b/(a*a))));
	printf("%12s %9.6f %18s %9.6f \n",
	       "y0 =  ",b,"uncertainty:",sigb);

	printf("%19s %14.6f \n","chi-squared: ",chi2);
	printf("%23s %10.6f %s\n","goodness-of-fit: ",q,
	       q==1 ? "(no Y errors supplied [dycol=])" : "");
      
	pearsn(x, y, npt, &r, &prob, &z);
      
	printf("%9s %g\n","r: ",r);
	printf("%12s %g\n","prob: ",prob);
	printf("%9s %g\n","z: ",z);
      
	if (mwt==0 && nsigma>0) {                
	  sigma = 0.0;
	  for(i=0; i<npt; i++)
	    sigma += sqr(y[i] - a*x[i] - b);
	  sigma /= (real) npt;
	  sigma = nsigma * sqrt(sigma);   /* critical distance */
	
	  for(i=0, j=0; i<npt; i++) {     /* loop over points */
	    d = ABS(y[i] - a*x[i] - b);
	    if (d > sigma) continue;
	    x[j] = x[i];                  /* shift them over */
	    y[j] = y[i];
	    if (dy) dy[j] = dy[i];
	    j++;
	  }
	  dprintf(0,"%d/%d points outside %g*sigma (%g)\n",
		  npt-j,npt,nsigma,sigma);
	  npt = j;
	}
      } /* mwt */
    } /* dxcol */
    
    if (outstr) write_data(outstr);
#endif
}




static inline real signo(real a, real b)
{
  if (b < 0) return -a;
  return a;
}
/*
 *  Determine the function needed to solve for the median.
 */

real medfunc(real a, real *x, real *y, int n)
{
  real rv = 0.0;
  int i;
  for (i=0; i<n; i++) {
    rv += y[i]*signo(1.0,x[i]-a*y[i]) - x[i]/a/a*signo(1.0,y[i]-x[i]/a);
  }
  return rv;
}


do_slope()
{
    real *x, *y, *dx, *dy, *dz;
    int i,j, n, mwt;
    real chi2,q,siga,sigb, sigma, d, sa, sb;
    real cov00, cov01, cov11, sumsq;
    real r, prob, z;
    real sumxx, sumxy, sumyy, a, f;
    real rms, a1, f1, a2, f2;

    warning("Some experimental code from miriad::immerge");
    n = npt;

    if (nxcol < 1) error("nxcol=%d",nxcol);
    if (nycol < 1) error("nycol=%d",nycol);
    x = ycol[0].dat;
    y = xcol[0].dat;

    sumxx = sumxy = sumyy = 0.0;
    for (i=0; i<n; i++) {
      sumxx += x[i]*x[i];
      sumxy += x[i]*y[i];
      sumyy += y[i]*y[i];
    }
    if (sumxx==0.0 || sumyy==0.0) error("zero data");
    a = sumxy/sumyy;
    rms = (sumxx + a*a*sumyy - 2*a*sumxy)/(n*sumyy);
    rms = sqrt(rms);

    a1 = a;
    f1 = medfunc(a1,x,y,n);
    a2 = a1 + signo(3*rms,f1);
    f2 = medfunc(a2,x,y,n);
    while (f1*f2 > 0.0) {
      a = 2*a2 - a1;
      a1 = a2;
      f1 = f2;
      a2 = a;
      f2 = medfunc(a2,x,y,n);
    }
    while ( ABS(a1-a2) > 0.001 * rms) {
      a = 0.5*(a1+a2);
      if (a==a1 || a==a2) break;
      f = medfunc(a,x,y,n);
      if (f*f1 > 0) {
	f1 = f;
	a1 = a;
      } else {
	f2 = f;
	a2 = a;
      }
    }
    printf("Fitting y=ax:\n");
    printf("a= %g\n",a);
}


char name[10] = "ABCDE";

do_ellipse()
{
    real *xdat, *ydat;
    real mat[5*5], vec[5], sol[5], a[6], cnt, xmean, ymean, x, y;
    real aa,bb,cc,dd,ee,pa,pp,cospp,sinpp,cospa,sinpa,ecc,al,r,ab,
         s1,s2,s3,y1,y2,y3,x0,y0, sum0, sum1, sum2, dx, dy, rr,
	 radmean, radsig, dr;
    real aaa, bbb, xp, yp, delta, sigfac;
    real siga, sigb, sigc, sigd, sige, fac1, fac2, fac3, dr1da, dr1db, dr1dc, sigr;
    int i, j;

    if (nxcol < 1) error("nxcol=%d",nxcol);
    if (nycol < 1) error("nycol=%d",nycol);
    xdat = xcol[0].dat;
    ydat = ycol[0].dat;

    if (npt < 5) {
        warning("Got %d; need minimum 5 points for an ellipse",npt);
        return 0;
    }
    xmean = ymean = 0.0;
    for (i=0; i<npt; i++) {             /* get mean of (x,y) as first guess */
      xmean += xdat[i];                 /* of center of ellips */
      ymean += ydat[i];                 /* to prevent all too much roundoff */
      dprintf(2,"Data: %f %f\n",xdat[i], ydat[i]);
    } 
    xmean /= npt;
    ymean /= npt;
    dprintf(1,"Estimate for center of ellips using %d points: %f %f\n",
					npt,xmean,ymean);
    if (hasvalue("estimate")) {
        if (nemoinpr(getparam("estimate"),vec,2) != 2) 
            error("estimate=%s needs two values for center of ellipse",
                   getparam("estimate"));
        xmean = vec[0];
        ymean = vec[1];
        dprintf(1,"Reset center of ellips: %f %f\n", xmean,ymean);
    }

    lsq_zero(5,mat,vec);

    for (i=0; i<npt; i++) {       /* gather all the stuff in matrix */
        x = xdat[i]-xmean;          /* treat (x,y) w.r.t. the mean center */
        y = ydat[i]-ymean;          /* of all points to prevent rounding err */
        a[0] = sqr(x);
        a[1] = 2*x*y;
        a[2] = sqr(y);
        a[3] = x;
        a[4] = y;
        a[5] = 1.0;
        lsq_accum(5,mat,vec,a,1.0); /* spread a[] into mat[] and vec[] */
    }

    for (i=0; i<5; i++)		/* print input matrix */
      dprintf(1,"( %9.3e %9.3e %9.3e %9.3e %9.3e  ) * ( %c ) = ( %9.3e )\n",
        mat[i*5+0], mat[i*5+1], mat[i*5+2], mat[i*5+3], mat[i*5+4],
        name[i], vec[i]);

    lsq_solve(5,mat,vec,sol);

    for (i=0; i<5; i++)		/* print inverse matrix & solution */
      dprintf(1,"( %9.3e %9.3e %9.3e %9.3e %9.3e  ) ; %c = %9.3e\n",
        mat[i*5+0], mat[i*5+1], mat[i*5+2], mat[i*5+3], mat[i*5+4],
        name[i], sol[i]);

    sigfac = 0.0;
    for (i=0; i<npt; i++) {
        x = xdat[i]-xmean;          /* treat (x,y) w.r.t. the mean center */
        y = ydat[i]-ymean;          /* of all points to prevent rounding err */
	sigfac += sqr(sol[0]*x*x+2*sol[1]*x*y+sol[2]*y*y+sol[3]*x+sol[4]*y-1);
    }
    sigfac /= (npt-5);
    sigfac = sqrt(sigfac);
    dprintf(1,"Sigma factor=%g\n",sigfac);
    siga = sigfac * sqrt(mat[0]);     /* trace elements are errors */
    sigb = sigfac * sqrt(mat[6]);
    sigc = sigfac * sqrt(mat[12]);
    sigd = sigfac * sqrt(mat[18]);
    sigd = sigfac * sqrt(mat[24]);

    /* Now that we have the coefficient a,b,c,d,e in:
     *      a.x^2 + 2bxy + cy^2 + dx + ey = 1
     * They have to be converted to human readable ones
     */

    aa = sol[0]; bb = sol[1]; cc = sol[2]; dd = sol[3]; ee = sol[4];
    dprintf(1,"Solutions  a..e: %g %g %g %g %g\n",aa,bb,cc,dd,ee);
    dprintf(1,"Errors  in a..e: %g %g %g %g %g\n",siga,sigb,sigc,sigd,sige);
    dprintf(1,"FracErr in a..e: %g %g %g %g %g\n",siga/aa,sigb/bb,sigc/cc,sigd/dd,sige/ee);

    delta = aa*cc-bb*bb;
    if (delta < 0)
        warning("You seem to have an hyperbola, not an ellipse");
    pp = atan(2.0*bb/(aa-cc));
    pa = 0.5*pp;        /* P.A. of undetermined axis */
    cospp = cos(pp);
    cospa = cos(pa);
    sinpp = sin(pp);
    sinpa = sin(pa);
    fac1 = sinpp*(aa+cc) + 2*bb;
    fac2 = sinpp*(aa+cc) - 2*bb;
    fac3 = sqr(aa-cc)+4*sqr(bb);
    r = fac1/fac2;
    r = sqrt(ABS(r));           /* r is now axis ratio b/a or a/b for ell/hyp */
    s1 = bb*ee-cc*dd;
    s2 = bb*dd-aa*ee;
    s3 = aa*cc-bb*bb;
    x0 = 0.5*s1/s3;
    y0 = 0.5*s2/s3;

    lsq_zero(3,mat,vec);
    vec[0] = aa*sqr(x0)+1; 
    vec[1] = bb*sqr(x0);
    vec[2] = cc*sqr(x0);
    lsq_cfill(3,mat,0,vec);
    vec[0] = 2*aa*x0*y0; 
    vec[1] = 2*bb*x0*y0+1;
    vec[2] = 2*cc*x0*y0; 
    lsq_cfill(3,mat,1,vec);
    vec[0] = aa*sqr(y0);
    vec[1] = bb*sqr(y0);
    vec[2] = cc*sqr(y0)+1;
    lsq_cfill(3,mat,2,vec);
    vec[0] = aa;
    vec[1] = bb;
    vec[2] = cc;
    lsq_solve(3,mat,vec,sol);
    dprintf(1,"Real A,B,C = %f %f %f \n",sol[0],sol[1],sol[2]);
#if 1
    if (ABS(cospp) > ABS(sinpp))    /* two ways to find ab: which is a^2 now */
        ab = (2/(sol[0]+sol[2] + (sol[0]-sol[2])/cospp));
    else                            /* use the biggest divisor */
        ab = (2/(sol[0]+sol[2] + 2*sol[1]/sinpp));
#else
    /* solve bug  that Kartik and I both identified Oct 2002 ? */
    ab = (2/(sol[0]+sol[2] + (sol[0]-sol[2])/cospp));
#endif
    /* compute error in axis ratio : for this one we only need 3 partial derivitives*/
    /* see also:  http://iraf.noao.edu/ADASS/adass_proc/adass_95/buskoi/buskoi.html */
    /* and the Jedrzejewski, R. 1987, MNRAS, 226, 747 reference                     */
#if 0
    /* first Mousumi version */
    dr1da = 4*bb*(sinpp-(2*bb*(aa+cc)*cospp)/fac3)/sqr(fac2);
    dr1dc = 4*bb*(sinpp+(2*bb*(aa+cc)*cospp)/fac3)/sqr(fac2);
#else
    /* new version, but it's really the same math, just refactorized:-) */
    dr1da = 4*bb*(sinpp-((aa+cc)*(aa-cc)*sinpp)/fac3)/sqr(fac2);
    dr1dc = 4*bb*(sinpp+((aa+cc)*(aa-cc)*sinpp)/fac3)/sqr(fac2);
#endif
    dr1db = (8*bb*cospp*(sqr(aa)-sqr(cc))+4*sinpp*(aa+cc)*fac3)/(fac3*sqr(fac2)); 
    sigr  = sqrt((sqr(siga*dr1da)+sqr(sigb*dr1db)+sqr(sigc*dr1dc))/(2*r));
    dprintf(1,"sigr=%g\n",sigr);
    dprintf(2,"aa,bb,cc,sinpp,cospp=%g %g %g %g %g\n",
  	       aa,bb,cc,sinpp,cospp);

    pa *= 180/PI;           /* PA is now in degrees */
    dprintf(1,"Before re-arranging: ab, r, pa = %f %f %f\n",ab,r,pa);
    if (delta > 0) {            /* ELLIPSE */
        if (r>1.0) { /* ellipse with r>1 means we've got a & b mixed up */
            r = 1/r;            /* so make r < 1 */
            ab = sqrt(ab)/r;    /* and ab is now the true semi major axis */
            pa -= 90.0;         /* and change PA by 90 degrees */
	    sigr *= sqr(r);     /* correct errors */
       }
       ecc = sqrt(1-r*r);       /* eccentricity of the ellipse ecc=sqrt(a^2-b^2)/a */
    } else {                    /* HYPERBOLA */
        if (ab < 0) {             /* if ab<0 means we've got a and b mixed up */
            r = 1/r;              /* so make r < 1 */
            ab = sqrt(ABS(ab))/r; /* so take abs value as semi major real axis */
            pa -= 90.0;           /* and change PA by 90 degrees */
        } else {
            ab = sqrt(ab); 
        }
        ecc = sqrt(1+r*r);      /* eccentricity of the hyperbola ecc=sqrt(a^2+b^2)/a */
    }

    x0 += xmean;        /* and finally correct to the real center */
    y0 += ymean;    

    if (nsigma<0) {         /* if no rejection of 'bad' points */
        if (delta>0) {
	  if(Qtab) {
	    printf("%f %f %f %f %f %f\n",ab,ecc,r,x0,y0,pa+90);
	  } else {
            printf("Ellips fit: \n");
            printf(" semi major axis: %g\n",ab);
            printf(" eccentricity:    %g axis ratio: %g  error: %g\n",ecc,r,sigr);
            printf(" x-center:        %g\n",x0);
            printf(" y-center:        %g\n",y0);
            printf(" P.A.:            %g\n",pa+90);
	  }
        } else {
            printf("Hyperbole fit: \n");
            printf(" semi real axis: %g\n",ab);
            printf(" eccentricity:   %g\n",ecc);
            printf(" x-center:       %g\n",x0);
            printf(" y-center:       %g\n",y0);
            printf(" P.A.:           %g\n",pa+90);
        }
    } else {
        dprintf(1,"Ellfit[%d]: a,b/a,x,y,pa=%g %g %g %g %g\n",npt,ab,r,x0,y0,pa+90);
			
        dprintf(0,"Now looking to delete points > %f sigma\n",nsigma);
        sum0 = sum1 = sum2 = 0.0;
        pa *= PI/180;           /* pa now in radians */
        for (i=0; i<npt; i++) {
            x = xdat[i] - x0;
            y = ydat[i] - y0;
            pp = atan2(y,x)-pa;
            rr = sqrt((sqr(x)+sqr(y))*(sqr(cos(pp))+sqr(sin(pp)/r)));
            sum0 += 1.0;
            sum1 += rr;
            sum2 += sqr(rr);
        }
        radmean = sum1/sum0;
        radsig = sqrt(sum2/sum0 - sqr(radmean));
	dprintf(0,"Deprojected radii:  mean: %f  sigma: %f\n",
			radmean, radsig);
	
        for (i=0; i<npt; i++) {
            x = xdat[i] - x0;
            y = ydat[i] - y0;
            pp = atan2(y,x)-pa;
            rr = sqrt((sqr(x)+sqr(y))*(sqr(cos(pp))+sqr(sin(pp)/r)));
            dr = (rr - radmean)/radsig;
            if (ABS(dr) < nsigma)
                printf("%g %g %g\n",xdat[i],ydat[i],dr);
            else
                dprintf(0,"Data Point %d deviates %f sigma from mean\n",
			i,dr);
        }
#if 0
        for (i=0; i<npt; i++) {
            pp = atan2(ydat[i]-y0,xdat[i]-x0)*180.0/PI;
            x = xdat[i] - xmean;
            y = ydat[i] - ymean;
            rr = aa*x*x + 2*bb*x*y + cc*y*y + dd*x + ee*y;
            printf("%d %g %g\n",i+1,(pp<0?pp+360:pp),rr);
        }
#endif    
    }

    if (outstr) {
       if (delta<0) error("Can't compute output hyperbola table yet");
       bbb = r*ab;        /* b = semi minor axis  */
       aaa = 1-r*r;            /* eps^2 */
       for (i=0; i<=360; i++) {
           cospp = cos(PI*i/180.0);
           rr = bbb/sqrt(1-aaa*cospp*cospp);
           xp = x0 + rr*cos(PI*(i+pa)/180.0);
           yp = y0 + rr*sin(PI*(i+pa)/180.0);
           fprintf(outstr,"%d %g %g\n", i+1, xp, yp);
       }
    }
}

do_imageshift()
{
    real *x, *y, *u, *v;
    /* this code was on 3b1 ... before CVS ... lost for now */
}


my_poly(bool Qpoly)
{ 
  real mat[(MAXCOL+1)*(MAXCOL+1)], vec[MAXCOL+1], sol[MAXCOL+1], a[MAXCOL+2], sum;
  int i, j;

  if (nycol<1) error("Need 1 value for ycol=");
  if (nxcol<order && !Qpoly) error("Need %d value(s) for xcol=",order);

  lsq_zero(order+1, mat, vec);
  for (i=0; i<npt; i++) {
    a[0] = 1.0;
    for (j=0; j<order; j++) {
      if (Qpoly)
	a[j+1] = a[j] * xcol[0].dat[i];     /* polynomial */
      else
	a[j+1] = xcol[j].dat[i];            /* plane */
    }
    a[order+1] = ycol[0].dat[i];
    lsq_accum(order+1,mat,vec,a,1.0);
  }
  if (order==0) printf("TEST = %g %g\n",mat[0], vec[0]);
  lsq_solve(order+1,mat,vec,sol);
  printf("%s fit of order %d:\n", Qpoly ? "Polynomial" : "Planar" , order);
  for (j=0; j<=order; j++) printf("%g ",sol[j]);
  printf("\n");

  if (outstr && Qpoly) {           /* output fitted values, if need be */
    for(i=0; i<npt; i++) {
      sum=sol[order];
      for (j=order-1; j>=0; j--)
	sum = (xcol[0].dat[i] * sum + sol[j]);
      fprintf(outstr,"%g %g %g\n",xcol[0].dat[i], ycol[0].dat[i], sum);
    }
  }
}

/*
 * PLANE:       y = b_0 + b_1 * x_1 + b_2 * x_2 + ... + b_n * x_n
 *
 *      used:   n = dimensionality of space in which hyper plane is fit
 */
 
do_plane()
{
    my_poly(FALSE);
}
/*
 * POLYNOMIAL:  y = b_0 + b_1 * x^1 + b_2 * x^2 + ... + b_n * x^n
 *
 *      used:   n = order of polynomial
 */
 
do_poly()
{
    my_poly(TRUE);
}


/*
 *   GAUSS2D:   y = A * exp( -[(x-x0)^2]/2b^2 )
 *        needs min. 3 points for a linear fit
 *            
 */

do_gauss1d()
{
  real mat[(MAXCOL+1)*(MAXCOL+1)], vec[MAXCOL+1], sol[MAXCOL+1], a[MAXCOL+2], sum;
  int i, j, gorder=2, neg=0;
  real A, b, x0, x,y;

  if (nycol<1) error("Need 1 value for ycol=");
  if (nxcol<1) error("Need 1 values for xcol=");
  if (npt<3) error("Need at least 3 data points for gauss1d fit");

  lsq_zero(gorder+1, mat, vec);
  for (i=0; i<npt; i++) {
    a[0] = 1.0;                     /* ln A - (x0^2+y0^2)/2b^2 */
    a[1] = xcol[0].dat[i];          /* x0/b^2  */
    a[2] = sqr(a[1]);               /* -1/2b^2 */
    if (ycol[0].dat[i] <= 0) {
      neg++;
      continue;
    }
    a[3] = log(ycol[0].dat[i]);
    lsq_accum(gorder+1,mat,vec,a,1.0);
  }
  if (neg > 0) {
    warning("Ignored %d negative data",neg);
    if (npt-neg < 3) error("Too many points rejected for a gauss1d fit");
  }
  lsq_solve(gorder+1,mat,vec,sol);
  printf("gauss1d fit:\n");
  for (j=0; j<=gorder; j++) printf("%g ",sol[j]);
  printf("\n\n");
  printf("  y = A * exp( -[(x-x0)^2]/2b^2 ):\n\n");
  if (sol[2] > 0) {
    warning("Bad gauss1d fit: 1/b^2 = %g\n",sol[2]);
    return 0;
  }

  x0 = -sol[1]/(2*sol[2]);
  b = sqrt(-1/(2*sol[2]));
  A = exp(sol[0] - sol[2]*(x0*x0));
  printf("     A  = %g\n",A);
  printf("     x0 = %g\n",x0);
  printf("     b  = %g\n",b);
  for (i=0; i<npt; i++) {
    x = xcol[0].dat[i];
    y = sol[0] + sol[1]*x + sol[2]*x*x;
    if (ycol[0].dat[i] <= 0) continue;
    dprintf(1,"%g %g %g => %g\n",x,log(ycol[0].dat[i]),y,log(ycol[0].dat[i])-y);
  }
}

/*
 *   GAUSS2D:   y = A * exp( -[(x-x0)^2 + (y-y0)^2]/2b^2 )
 *        needs min. 4 points for a linear fit
 *            
 */

do_gauss2d()
{
  real mat[(MAXCOL+1)*(MAXCOL+1)], vec[MAXCOL+1], sol[MAXCOL+1], a[MAXCOL+2], sum;
  int i, j, gorder=3, neg=0;
  real A, b, x0, y0, x,y,z;

  if (nycol<1) error("Need 1 value for ycol=");
  if (nxcol<2) error("Need 2 values for xcol=");
  if (npt<4) error("Need at least 4 data points for gauss2d fit");

  lsq_zero(gorder+1, mat, vec);
  for (i=0; i<npt; i++) {
    a[0] = 1.0;                     /* ln A - (x0^2+y0^2)/2b^2 */
    a[1] = xcol[0].dat[i];          /* x0/b^2  */
    a[2] = xcol[1].dat[i];          /* y0/b^2  */
    a[3] = sqr(a[1]) + sqr(a[2]);   /* -1/2b^2 */
    if (ycol[0].dat[i] <= 0) {
      neg++;
      continue;
    }
    a[4] = log(ycol[0].dat[i]);
    lsq_accum(gorder+1,mat,vec,a,1.0);
  }
  if (neg > 0) {
    warning("Ignored %d negative data",neg);
    if (npt-neg < 4) error("Too many points rejected for a gauss2d fit");
  }
  lsq_solve(gorder+1,mat,vec,sol);
  printf("gauss2d fit:\n");
  for (j=0; j<=gorder; j++) printf("%g ",sol[j]);
  printf("\n\n");
  printf("  y = A * exp( -[(x-x0)^2 + (y-y0)^2]/2b^2 ):\n\n");
  if (sol[3] > 0) {
    warning("Bad gauss2d fit: 1/b^2 = %g\n",sol[3]);
    return 0;
  }

  x0 = -sol[1]/(2*sol[3]);
  y0 = -sol[2]/(2*sol[3]);
  b = sqrt(-1/(2*sol[3]));
  A = exp(sol[0] - sol[3]*(x0*x0+y0*y0));
  printf("     A  = %g\n",A);
  printf("     x0 = %g\n",x0);
  printf("     y0 = %g\n",y0);
  printf("     b  = %g\n",b);
  for (i=0; i<npt; i++) {
    x = xcol[0].dat[i];
    y = xcol[1].dat[i];
    z = sol[0] + sol[1]*x + sol[2]*y + sol[3]*(x*x+y*y);
    if (ycol[0].dat[i] <= 0) continue;
    dprintf(1,"%g %g %g %g => %g\n",x,y,log(ycol[0].dat[i]),z,log(ycol[0].dat[i])-z);
  }
}




do_fourier()
{
  real mat[(MAXCOL+1)*(MAXCOL+1)], vec[MAXCOL+1], sol[MAXCOL+1], a[MAXCOL+2];
  real sum, theta, amp, pha, sigma;
  int i, j, dim = 2*order+1;

  if (dim > MAXCOL) error("order=%d too high",order);

  lsq_zero(dim, mat, vec);

  for (i=0; i<npt; i++) {
    a[0] = 1.0;
    for (j=0; j<order; j++) {
      theta = xcol[0].dat[i] * PI/180;
      a[2*j+1] = cos((j+1)*theta);
      a[2*j+2] = sin((j+1)*theta);
    }
    a[dim] = ycol[0].dat[i];
    lsq_accum(dim,mat,vec,a,1.0);
  }
  lsq_solve(dim,mat,vec,sol);
  printf("fourier fit of order %d:\n", order);
  printf("\ncos/sin amplitudes:\n");
  for (j=0; j<=2*order; j++) printf("%g ",sol[j]);
  printf("\ncos amp/phase solutions:\n");
  printf("%g ",sol[0]);
  for (j=0; j<order; j++) {
    amp = sqrt(sqr(sol[2*j+1])+sqr(sol[2*j+2]));
    pha = atan2( sol[2*j+2] , sol[2*j+1]) * 180 / PI;
    printf("%g %g ",amp,pha);
  }
  printf("\n");


  sigma = 0.0;
  for(i=0; i<npt; i++) {
    sum=sol[0];
    for (j=0; j<order; j++) {
      theta =  xcol[0].dat[i] * PI/180;
      sum += sol[2*j+1]*cos((j+1)*theta) + sol[2*j+2]*sin((j+1)*theta);
    }
    if (outstr) fprintf(outstr,"%g %g %g\n",xcol[0].dat[i], ycol[0].dat[i], sum);
    sigma += sqr(ycol[0].dat[i]-sum);
  }
  sigma = sqrt(sigma/(npt-2*order-1));
  printf("chisq=%g\n",sigma);

}



/* fit a peak, 
 * for now this is the my_poly() code, though forced with order=2 
 * first find the peak, then take the two points on either side
 * to find an exact solution
 */

do_peak()
{
    real mat[(MAXCOL+1)*(MAXCOL+1)], vec[MAXCOL+1], sol[MAXCOL+1], a[MAXCOL+2];
    real *x, *y, ymax;
    int i, j, k, range=1;

    order = 2;

    if (nycol<1) error("Need 1 value for ycol=");

    x = xcol[0].dat;
    y = ycol[0].dat;
    for (i=1, j=0; i<npt; i++)			/* find the peak, at j */
        if (y[j] < y[i]) j = i;
    if (j==0 || j==npt-1) {			/* handle edge cases */
        warning("Peak at the edge");
        printf("%g %g\n",x[j],y[j]);
        return 0;
    }
    if (range==2) {
    	if (j==1 || j==npt-2) {
    	    warning("Peak too close to edge");
    	    return 0;
    	}
    }

    lsq_zero(order+1, mat, vec);
    for (i=j-range; i<=j+range; i++) {
        a[0] = 1.0;
        for (k=0; k<order; k++) {
            a[k+1] = a[k] * (x[i]-x[j]);
        }
        a[order+1] = y[i];
        lsq_accum(order+1,mat,vec,a,1.0);
    }
    lsq_solve(order+1,mat,vec,sol);
    dprintf(1,"Poly2 fit near j=%d (%g,%g)\n",j+1,x[j],y[j]);
    printf("Peak:x,y= %g %g\n",
            x[j] - sol[1]/(2*sol[2]),
	   sol[0]-sol[1]*sol[1]/(4*sol[2]));
    return 0;
}

/* find a zero point
 * for now this is the my_poly() code, though forced with order=2 
 * first find the peak, then take the two points on either side
 * to find an exact solution
 */

do_zero()
{
    real *x, *y, zero;
    int i, j;
    int nzero = 1, izero = 0;

    order = 2;

    if (nycol<1) error("Need 1 value for ycol=");

    x = xcol[0].dat;
    y = ycol[0].dat;
    j = -1;
    for (i=1; i<npt; i++) {	       /* find the zero point */
      if (y[i]*y[i-1] < 0) {
	j = i;
	izero++;
	if (izero == nzero) break;
      }
    }
    if (j < 0) error("Could not find a zero");
    zero = x[j-1] - y[j-1]*(x[j]-x[j-1])/(y[j]-y[j-1]);
    printf("Zero:x= %g\n",zero);
}



do_area()
{
    real *x, *y, ymax, xmean, ymean;
    int i, j, k;

    order = 2;

    if (nycol<1) error("Need 1 value for ycol=");

    x = xcol[0].dat;
    y = ycol[0].dat;
    for (i=0; i<npt; i++) {
      xmean += x[i];
      ymean += y[i];
    }
    xmean /= npt;
    ymean /= npt;

    error("ah, not done coding here yet");

    /* allocate temp array, 
       compute angles 
       sort an index array by angles
       compute half the sum of x[i+1]*y[i]-x[i]*y[i+1]
       that's the area
    */
}


    /* NumRec emulator - no chi^2 */
myfit(real *x,real *y,int npt,real *dy,int mwt,
      real *b, real *a,real *sigb,real *siga,real *chi2, real *q)
{
  real mat[4], vec[2], sol[2], aa[3];
  int i;

  dprintf(0,"local fit\n");
  lsq_zero(2,mat,vec);
  for (i=0; i<npt; i++) {
    aa[0] = 1.0;       /* mat */
    aa[1] = x[i];
    aa[2] = y[i];       /* rsh */
    lsq_accum(2, mat,  vec, aa, 1.0);    /* 1.0 -> dy */
  }
  dprintf(0,"fit: mat = %f %f %f %f\n", mat[0], mat[1], mat[2], mat[3]);
  dprintf(0,"fit: vec = %f %f\n", vec[0], vec[1]);
  lsq_solve(2, mat, vec, sol);
  *a = sol[1];
  *b = sol[0];
}

/* NR version of fit(), using gamma functions */

extern real gammq(real a, real x);

void fit(real *x, real *y,int ndata,real *sig,int mwt,
	 real *a, real *b,real *siga,real *sigb,real *chi2,real *q)
{
  int i;
  real wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;

  *b=0.0;
  if (mwt) {
    ss=0.0;
    for (i=0;i<ndata;i++) {
      wt=1.0/sqr(sig[i]);
      ss += wt;
      sx += x[i]*wt;
      sy += y[i]*wt;
    }
  } else {
    for (i=0;i<ndata;i++) {
      sx += x[i];
      sy += y[i];
    }
    ss=ndata;
  }
  sxoss=sx/ss;
  if (mwt) {
    for (i=0;i<ndata;i++) {
      t=(x[i]-sxoss)/sig[i];
      st2 += t*t;
      *b += t*y[i]/sig[i];
    }
  } else {
    for (i=0;i<ndata;i++) {
      t=x[i]-sxoss;
      st2 += t*t;
      *b += t*y[i];
    }
  }
  *b /= st2;
  *a=(sy-sx*(*b))/ss;
  *siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
  *sigb=sqrt(1.0/st2);
  *chi2=0.0;
  if (mwt == 0) {
    for (i=0;i<ndata;i++)
      *chi2 += sqr(y[i]-(*a)-(*b)*x[i]);
    *q=1.0;
    sigdat=sqrt((*chi2)/(ndata-2));
    *siga *= sigdat;
    *sigb *= sigdat;
  } else {
    for (i=0;i<ndata;i++)
      *chi2 += sqr((y[i]-(*a)-(*b)*x[i])/sig[i]);
    *q=gammq(0.5*(ndata-2),0.5*(*chi2));	
  }
}


void fit1(x,y,ndata,sig,mwt,a,b,siga,sigb,chi2,q)
real x[],y[],sig[],*a,*b,*siga,*sigb,*chi2,*q;
int ndata,mwt;
{
  real mat[4], vec[2], sol[2], aa[3];
  int i;

  warning("testing fit using lsq_ routines");

  lsq_zero(2,mat,vec);
  for (i=0; i<ndata; i++) {
    aa[0] = 1.0;
    aa[1] = x[i];
    aa[2] = y[i];
    lsq_accum(2,mat,vec,aa,1.0);
  }
  printf("mat   : %g %g => %g \n",mat[0],mat[1],vec[0]);
  printf("      : %g %g => %g \n",mat[2],mat[3],vec[1]);
  lsq_solve(2,mat,vec,a);
  printf("a+bx  :  a=%g  b=%g\n", aa[0],aa[1]);
  printf("mat^-1: %g %g \n       %g %g\n", mat[0],mat[1],mat[2],mat[3]);


  stop(0);       
}

