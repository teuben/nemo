/* TABSPLINE:   interpolation into a table: spline and linear
 *
 *
 *   16-apr-98   2.0 : split off from TOOLBOX in spline.c  (kernel/misc)    pjt
 *   22-jan-01   2.1 : allow x= and y= to be arrays, for Mousumi Das	    pjt
 *
 *   TODO:
 *      - spline vs. linear
 *      - error bars, via derivatives
 */

#include <stdinc.h> 
#include <getparam.h>

string defv[] = {
    "in=???\n       Input table file",
    "xcol=1\n       X column",
    "ycol=2\n       Y Column",
    "x=\n           For these X's, find Y",
    "y=\n           For these Y's, find X",
    "n=0\n          Which one to find (0=all)",
    "VERSION=2.1\n  24-jan-00 PJT",
    NULL,

};

string usage="first and second derivatives of a table";

#define MAXZERO     64
#define MAXDATA  16384

nemo_main()
{
    int colnr[2];
    real *coldat[2], *xdat, *ydat, xmin, xmax, ymin, ymax;
    real fy[MAXZERO], fx[MAXZERO], xp[MAXDATA], yp[MAXDATA], x, y;
    stream instr;
    int i, j, n, nx, ny, nmax, izero, nzero = getiparam("n");
    bool Qx, Qy;

    if (nzero > MAXZERO)
        error("MAXZERO=%d: Too many zero's to search for",MAXZERO);

    nmax = file_lines(getparam("in"),MAXLINES);
    xdat = coldat[0] = (real *) allocate(nmax*sizeof(real));
    ydat = coldat[1] = (real *) allocate(nmax*sizeof(real));
    colnr[0] = getiparam("xcol");
    colnr[1] = getiparam("ycol");
    Qx = hasvalue("x");
    Qy = hasvalue("y");

    instr = stropen(getparam("in"),"r");
    n = get_atable(instr,2,colnr,coldat,nmax);

    for(i=0; i<n; i++) {
        dprintf(2,"%g %g\n",xdat[i],ydat[i]);
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

    if (Qx) {
      nx = nemoinpd(getparam("x"),xp,MAXDATA);
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
                    printf("%g %g\n",x,fx[izero]);
                }
                izero++;
            }
        }
        dprintf(1,"# Found %d zero's\n",izero);
      }
    }

    if (Qy) {
      ny = nemoinpd(getparam("y"),yp,MAXDATA);
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
                    printf("%g %g\n",fy[izero],y);
                }
                izero++;
            }
        }
        dprintf(1,"# Found %d zero's\n",izero);
      }
    }
}


