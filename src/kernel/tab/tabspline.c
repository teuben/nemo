/* TABSPLINE:   interpolation into a table: spline and linear
 *
 *
 *   16-apr-98   2.0 : split off from TOOLBOX in spline.c  (kernel/misc)    pjt
 *   22-jan-01   2.1 : allow x= and y= to be arrays, for Mousumi Das	    pjt
 *   5-apr-01	 2.2 : added format=					    pjt
 *   8-apr-01       a: fixed SINGLEPREC operation
 *   9-sep-01    3.0   GSL enabled
 *  23-sep-01       a  ->nemo_file_lines
 *
 *   TODO:
 *      - spline vs. linear
 *      - error bars, via derivatives
 *      - (GSL) integrate function
 */

#include <stdinc.h> 
#include <getparam.h>

#if HAVE_GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#endif

string defv[] = {
    "in=???\n       Input table file",
    "xcol=1\n       X column",
    "ycol=2\n       Y Column",
    "x=\n           For these X's, find Y",
    "y=\n           For these Y's, find X",
    "n=0\n          Which one to find (0=all) - not used in GSL",
    "format=%g\n    Output format",
#if HAVE_GSL
    "type=cspline\n Spline interpolation type (only for GSL)",
    "nder=0\n       Number of derivates to show (0,1,2)",
#endif
    "VERSION=3.0\n  9-sep-01 PJT",
    NULL,

};

string usage="first and second derivatives of a table";

#define MAXZERO     64
#define MAXDATA  16384

nemo_main()
{
    int colnr[2];
    real *coldat[2], *xdat, *ydat, xmin, xmax, ymin, ymax;
    real fy[MAXZERO], fx[MAXZERO], xp[MAXDATA], yp[MAXDATA], x, y, xd, yd;
    stream instr;
    string fmt, stype;
    char fmt1[100], fmt2[100];
    int i, j, n, nx, ny, nmax, izero, nder, nzero = getiparam("n");
    bool Qx, Qy;
#if HAVE_GSL
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline;
#endif

    if (nzero > MAXZERO)
        error("MAXZERO=%d: Too many zero's to search for",MAXZERO);

    nmax = nemo_file_lines(getparam("in"),MAXLINES);
    xdat = coldat[0] = (real *) allocate(nmax*sizeof(real));
    ydat = coldat[1] = (real *) allocate(nmax*sizeof(real));
    colnr[0] = getiparam("xcol");
    colnr[1] = getiparam("ycol");
    Qx = hasvalue("x");
    Qy = hasvalue("y");
    fmt = getparam("format");
    sprintf(fmt2,"%s %s",fmt,fmt);
    sprintf(fmt1," %s",fmt,fmt);
    dprintf(1,"Using format=\"%s\"\n",fmt2);


    instr = stropen(getparam("in"),"r");
    n = get_atable(instr,2,colnr,coldat,nmax);

    for(i=0; i<n; i++) {
        dprintf(2,fmt,xdat[i],ydat[i]);
        if (i==0) {
            xmin = xmax = xdat[0];
            ymin = ymax = ydat[0];
        } else {
            xmax = MAX(xmax,xdat[i]);
            ymax = MAX(ymax,ydat[i]);
            xmin = MIN(xmin,xdat[i]);
            ymin = MIN(ymin,ydat[i]);
        }
    }
    dprintf(1,"Xrange: %g : %g  Yrange: %g : %g\n",xmin,xmax,ymin,ymax);

#if HAVE_GSL
    stype = getparam("type");
    nder = getiparam("nder");
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
      dprintf(0,"Evaluating Y(x) using GSL\n");
      gsl_spline_init(spline, xdat, ydat, n);
      nx = nemoinpr(getparam("x"),xp,MAXDATA);
      if (nx<0) error("Parsing x=%s",getparam("x"));
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
      dprintf(0,"Evaluating X(y) using GSL\n");
      gsl_spline_init(spline, ydat, xdat, n);
      ny = nemoinpr(getparam("y"),yp,MAXDATA);
      if (ny<0) error("Parsing y=%s",getparam("y"));
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
    strcat(fmt2,"\n");

    if (Qx) {
      nx = nemoinpr(getparam("x"),xp,MAXDATA);
      if (nx<0) error("Parsing x=%s",getparam("x"));
      for (j=0; j<nx; j++) {
        x = xp[j];
        if (x<xmin || x>xmax) error("x=%g out of range %g : %g",x,xmin,xmax);
        izero = 0;
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
      }
    }

    if (Qy) {
      ny = nemoinpr(getparam("y"),yp,MAXDATA);
      if (ny<0) error("Parsing y=%s",getparam("y"));
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
