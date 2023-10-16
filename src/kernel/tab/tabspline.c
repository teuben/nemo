/* TABSPLINE:   interpolation into a table: spline and linear
 *
 *
 *   16-apr-98   2.0 : split off from TOOLBOX in spline.c  (kernel/misc)    pjt
 *   22-jan-01   2.1 : allow x= and y= to be arrays, for Mousumi Das	    pjt
 *   5-apr-01	 2.2 : added format=					    pjt
 *   8-apr-01       a: fixed SINGLEPREC operation
 *   9-sep-01    3.0   GSL enabled
 *  23-sep-01       a  ->nemo_file_lines
 *  26-may-02    3.1   add derivatives (old TOOLBOX from spline.c) for non-GSL
 *  27-jul-02    3.2   implement nder= using spline.c if GSL not present    PJT
 *  28-aug-02    3.3   allow to use 'all' for  x= or y= ?
 *   9-aug-22    3.4   disabled pipe reading until table v2 implemented 
 *                     also: remove data when either X or Y values repeat   PJT
 *
 *   TODO:
 *      - spline vs. linear
 *      - error bars, via derivatives
 *      - (GSL) integrate function
 */

#include <stdinc.h> 
#include <getparam.h>
#include <table.h>

#if HAVE_GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#else
extern void spline(real *coef, real *x, real *y, int n);
extern real seval(real x0, real *x, real *y, real *coef, int n);
extern real spldif(real x0, real *x, real *y, real *coef, int n);
extern real spldif2(real x0, real *x, real *y, real *coef, int n);
#endif

string defv[] = {
    "in=???\n       Input table file",
    "xcol=1\n       X column (1=first)",
    "ycol=2\n       Y Column",
    "x=\n           For these X's, find Y; 'all' uses all xcol's",
    "y=\n           For these Y's, find X; 'all' uses all ycol's",
    "n=0\n          Which roots one to find (0=all) - not used in GSL",
    "format=%g\n    Output format",
#if HAVE_GSL
    "type=cspline\n Spline interpolation type (only for GSL)",
#endif
    "nder=0\n       Number of derivates to show (0,1,2)",
    "VERSION=3.4\n  9-aug-2022 PJT",
    NULL,

};

string usage="interpolation and first two derivatives of a function table";

#define MAXZERO     64
#define MAXDATA  16384

void nemo_main()
{
    int colnr[2];
    real *coldat[2], *xdat, *ydat, xmin, xmax, ymin, ymax;
    int *mdat;
    real fy[MAXZERO], fx[MAXZERO], xp[MAXDATA], yp[MAXDATA], x, y, xd, yd;
    stream instr;
    string fmt, stype;
    char fmt1[100], fmt2[100];
    int i, j, n, nx, ny, nmax, izero;
    int nzero = getiparam("n");
    int nder = getiparam("nder");
    int ndup;
    bool Qx, Qy, Qxall, Qyall;
#if HAVE_GSL
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline;
#else
    real *sdat;
#endif

    if (nzero > MAXZERO)
        error("MAXZERO=%d: Too many zero's to search for",MAXZERO);

#if 0
    nmax = nemo_file_lines(getparam("in"),MAXLINES);
#else    
    nmax = nemo_file_lines(getparam("in"),0);  // until we use table V2
#endif    
    xdat = coldat[0] = (real *) allocate(nmax*sizeof(real));
    ydat = coldat[1] = (real *) allocate(nmax*sizeof(real));
    mdat = (int *) allocate(nmax*sizeof(int));
    colnr[0] = getiparam("xcol");
    colnr[1] = getiparam("ycol");
    Qx = hasvalue("x");
    Qy = hasvalue("y");
    fmt = getparam("format");
    sprintf(fmt2,"%s %s",fmt,fmt);
    sprintf(fmt1,"%s %s",fmt,fmt);
    dprintf(1,"Using format=\"%s\"\n",fmt2);
    Qxall = (Qx && streq("all",getparam("x")));
    Qyall = (Qy && streq("all",getparam("y")));

    instr = stropen(getparam("in"),"r");
    n = get_atable(instr,2,colnr,coldat,nmax);

    ndup = 0;
    for(i=0; i<n; i++) {
        dprintf(2,fmt,xdat[i],ydat[i]);
	mdat[i] = 1;
        if (i==0) {
            xmin = xmax = xdat[0];
            ymin = ymax = ydat[0];
        } else {
            xmax = MAX(xmax,xdat[i]);
            ymax = MAX(ymax,ydat[i]);
            xmin = MIN(xmin,xdat[i]);
            ymin = MIN(ymin,ydat[i]);
	    if (xdat[i] == xdat[i-1]) {
	      mdat[i] = 0;
	      ndup++;
	      dprintf(1,"X=%f at %d duplicated\n", xdat[i], i);
	    }
	    if (ydat[i] == ydat[i-1]) {
	      mdat[i] = 0;
	      ndup++;
	      dprintf(1,"Y=%f at %d duplicated\n", ydat[i], i);
	    }
        }
    }
    dprintf(1,"Xrange: %g : %g  Yrange: %g : %g\n",xmin,xmax,ymin,ymax);
    if (ndup>0) {
      warning("There were %d/%d X or Y duplications - removing from table (a hack)",ndup,n);
      for(i=0, j=0; i<n; i++) {
	if (mdat[i]) {
	  if (j<i) {
	    xdat[j] = xdat[i];
	    ydat[j] = ydat[i];
	  }
	  j++;
	}
      }
      n=j;
      warning("New n=%d", n);
    }

#if HAVE_GSL
    stype = getparam("type");
    dprintf(1,"Using interpolation type=%s\n",stype);
    if (streq(stype,"linear"))
      spline = gsl_spline_alloc(gsl_interp_linear, n);
    else if (streq(stype,"cspline"))
      spline = gsl_spline_alloc(gsl_interp_cspline, n);
    else if (streq(stype,"cspline_periodic"))
      spline = gsl_spline_alloc(gsl_interp_cspline_periodic, n);
    else if (streq(stype,"akima"))
      spline = gsl_spline_alloc(gsl_interp_akima, n);
    else if (streq(stype,"akima_periodic"))
      spline = gsl_spline_alloc(gsl_interp_akima_periodic, n);
    else
      error("Illegal spline interpolation type %s,\n"
	    "use: linear,cspline[_periodic],akima[_periodic]",stype);
    if (sizeof(real) != sizeof(double))
      error("Program not compiled with real==double, cannot use GSL");

    if (Qx) {
      dprintf(1,"Evaluating Y(x) using GSL\n");
      gsl_spline_init(spline, xdat, ydat, n);
      if (Qxall ) {
	nx = n;
	for (i=0; i<nx; i++) xp[i] = xdat[i];
      } else {
	nx = nemoinpr(getparam("x"),xp,MAXDATA);
	if (nx<0) error("Parsing x=%s",getparam("x"));
      }
      for (j=0; j<nx; j++) {
        x = xp[j];
        if (x<xmin || x>xmax) error("x=%g out of range %g : %g",x,xmin,xmax);
	y = gsl_spline_eval(spline, x, acc);
        printf(fmt2,x,y);
	if (nder > 0) {
	  yd = gsl_spline_eval_deriv(spline, x, acc);
	  printf(fmt1,yd);
	  if (nder > 1) {
	    yd = gsl_spline_eval_deriv2(spline, x, acc);
	    printf(fmt1,yd);
	  }
	}
	printf("\n");
      }
      gsl_spline_free(spline);
    }

    if (Qy) {
      dprintf(1,"Evaluating X(y) using GSL\n");
      gsl_spline_init(spline, ydat, xdat, n);
      if (Qyall ) {
	ny = n;
	for (i=0; i<ny; i++) yp[i] = ydat[i];
      } else {
	ny = nemoinpr(getparam("y"),yp,MAXDATA);
	if (ny<0) error("Parsing y=%s",getparam("y"));
      }
      for (j=0; j<ny; j++) {
        y = yp[j];
        if (y<ymin || y>ymax) error("y=%g out of range %g : %g",y,ymin,ymax);
	x = gsl_spline_eval(spline, y, acc);
        printf(fmt2,x,y);
	if (nder > 0) {
	  xd = gsl_spline_eval_deriv(spline, y, acc);
	  printf(fmt1,xd);
	  if (nder > 1) {
	    xd = gsl_spline_eval_deriv2(spline, y, acc);
	    printf(fmt1,xd);
	  }
	}
	printf("\n");
      }
      gsl_spline_free(spline);
    }
    gsl_interp_accel_free(acc);
#else

    if (Qx) {
      dprintf(1,"Evaluating Y(x) using spline\n");
      if (Qxall ) {
	nx = n;
	for (i=0; i<nx; i++) xp[i] = xdat[i];
      } else {
	nx = nemoinpr(getparam("x"),xp,MAXDATA);
	if (nx<0) error("Parsing x=%s",getparam("x"));
      }
      dprintf(1,"[Using spline.c] n=%d\n",n);
      sdat = (real *) allocate(sizeof(real)*n*3);
      spline(sdat,xdat,ydat,n);
      for (j=0; j<nx; j++) {
        x = xp[j];
        if (x<xmin || x>xmax) error("x=%g out of range %g : %g",x,xmin,xmax);
        izero = 0;
#if 0
        for (i=1; i<n; i++) {
            if (xdat[i-1] <= x && x < xdat[i] ||
                xdat[i] <= x && x < xdat[i-1]) {
                if (nzero==0 || izero < nzero) {
                    fx[izero] = ydat[i-1] +
                        (x-xdat[i-1])*(ydat[i]-ydat[i-1])/(xdat[i]-xdat[i-1]);
                    printf(fmt2,x,fx[izero]);
                }
                izero++;
            }
        }
        dprintf(1,"# Found %d zero's\n",izero);
#else
	y = seval(x,xdat,ydat,sdat,n);
        printf(fmt2,x,y);
	if (nder > 0) {
	  yd = spldif(x,xdat,ydat,sdat,n);
	  printf(fmt1,yd);
	  if (nder > 1) {
	    yd = spldif2(x,xdat,ydat,sdat,n);
	    printf(fmt1,yd);
	  }
	}
	printf("\n");
#endif
      }
    }

    if (Qy) {
      dprintf(1,"Evaluating X(y)\n");
      if (Qyall ) {
	ny = n;
	for (i=0; i<ny; i++) yp[i] = ydat[i];
      } else {
	ny = nemoinpr(getparam("y"),yp,MAXDATA);
	if (ny<0) error("Parsing y=%s",getparam("y"));
      }
      for (j=0; j<ny; j++) {
        y = yp[j];
        if (y<ymin || y>ymax) error("y=%g out of range %g : %g",y,ymin,ymax);
        izero = 0;
        for (i=1; i<n; i++) {
            if (ydat[i-1] <= y && y < ydat[i] ||
                ydat[i] <= y && y < ydat[i-1]) {
                if (nzero==0 || izero < nzero) {
                    fy[izero] = xdat[i-1] +
                        (y-ydat[i-1])*(xdat[i]-xdat[i-1])/(ydat[i]-ydat[i-1]);
                    printf(fmt2,fy[izero],y);
                }
                izero++;
            }
        }
        dprintf(1,"# Found %d zero's\n",izero);
      }
    }
#endif
}
