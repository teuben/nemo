#ifdef TOOLBOX
#include <getparam.h>

string defv[] = {
    "in=???\n       Input table file",
    "xcol=1\n       X column",
    "ycol=2\n       Y Column",
    "x=\n           For this X, find Y",
    "y=\n           For this Y, find X",
    "n=0\n          Which one to find (0=all)",
    "VERSION=2.0\n  15-apr-98 PJT",
    NULL,

};

string usage="first and second derivatives of a table";

#define MAXZERO  64

nemo_main()
{
    int colnr[2];
    real *coldat[2], *xdat, *ydat, xmin, xmax, ymin, ymax;
    real x,y, fy[MAXZERO], fx[MAXZERO];
    stream instr;
    int i, n, nmax, izero, nzero = getiparam("n");
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
        x = getdparam("x");
        if (x<xmin || x>xmax) error("x=%g out of range %g : %g",x,xmin,xmax);
        izero = 0;
        for (i=1; i<n; i++) {
            if (xdat[i-1] <= x && x < xdat[i]) {
                if (nzero==0 || izero < nzero) {
                    fx[izero] = ydat[i-1] +
                        (x-xdat[i-1])*(ydat[i]-ydat[i-1])/(xdat[i]-xdat[i-1]);
                    printf("%g %g\n",x,fx[izero]);
                }
                izero++;
            }
        }
        dprintf(0,"# Found %d zero's\n",izero);
    }

    if (Qy) {
        y = getdparam("y");
        if (y<ymin || y>ymax) error("y=%g out of range %g : %g",y,ymin,ymax);
        izero = 0;
        for (i=1; i<n; i++) {
            if (ydat[i-1] <= y && y < ydat[i]) {
                if (nzero==0 || izero < nzero) {
                    fy[izero] = xdat[i-1] +
                        (y-ydat[i-1])*(ydat[i]-ydat[i-1])/(xdat[i]-xdat[i-1]);
                    printf("%g %g\n",fy[izero],y);
                }
                izero++;
            }
        }
        dprintf(0,"# Found %d zero's\n",izero);
    }

}

#endif

