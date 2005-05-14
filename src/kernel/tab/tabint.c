/* TABINT:   integrate a sorted table
 *
 *
 *      13-may-05    Q&D first version, cloned off tabfilter
 */

#include <stdinc.h> 
#include <getparam.h>
#include <spline.h>
#include <extstring.h>
#include <table.h>

string defv[] = {
  "in=???\n         Input table file",
  "xcol=1\n         Column with X coordinate (better be sorted in that column",
  "ycol=2\n         Column with Y coordinate of function",
  "step=\n          Integration step if resampling used",
  "normalize=f\n    Normalize integral",
  "VERSION=0.1\n    13-may-05 PJT",
  NULL,

};

string usage="integrate a sorted table";

string cvsid="$Id$";


#define MAXDATA  16384

extern string *burststring(string, string);
extern void freestrings(string *);

void nemo_main()
{
  int colnr[2], i, n, nmax;
  real *coldat[2], *xdat, *ydat, xmin, xmax, ymin, ymax;
  real x, y, xold, yold, dx, sum, sum0, *sdat;
  stream instr;
  string spectrum = getparam("in");
  bool Qnorm = getbparam("normalize");
  
  nmax = nemo_file_lines(spectrum,MAXLINES);
  xdat = coldat[0] = (real *) allocate(nmax*sizeof(real));
  ydat = coldat[1] = (real *) allocate(nmax*sizeof(real));
  colnr[0] = getiparam("xcol");
  colnr[1] = getiparam("ycol");
  instr = stropen(spectrum,"r");
  n = get_atable(instr,2,colnr,coldat,nmax);
  strclose(instr);
  
  for(i=0; i<n; i++) {
    dprintf(2,"%g %g\n",xdat[i],ydat[i]);
    if (i==0) {
      xmin = xmax = xdat[0];
      ymin = ymax = ydat[0];
    } else {
      if (xdat[i] <= xdat[i-1]) 
	error("Column %d must be sorted",colnr[0]);
      xmax = MAX(xmax,xdat[i]);
      ymax = MAX(ymax,ydat[i]);
      xmin = MIN(xmin,xdat[i]);
      ymin = MIN(ymin,ydat[i]);
    }
  }
  dprintf(1,"X range: %g : %g\n",xmin,xmax);
  dprintf(1,"Y range: %g : %g\n",ymin,ymax);

  sum = sum0 = 0.0;

  if (hasvalue("step")) {     /* integration via summation of a resampled function using splines */
    dx = getdparam("step");

    /* setup a spline interpolation table into the function */
    sdat = (real *) allocate(sizeof(real)*n*3);
    spline(sdat,xdat,ydat,n);

    yold = ydat[0];
    xold = xmin;
    for (x=xmin+dx; x<xmax; x += dx) {
      y = seval(x,xdat,ydat,sdat,n);
      dprintf(1,"%g -> %g\n",x,y);

      sum  += 0.5*(yold+y)*(x-xold);
      sum0 += (x-xold);
      yold = y;
      xold = x;
    }
    y = ydat[n-1];
    sum  += 0.5*(yold+y)*(xmax-xold);
    sum0 += (xmax-xold);
  } else {              /* simple summation via the datapoints itself */
    for (i=1; i<n; i++) {
      sum += 0.5*(ydat[i]+ydat[i-1])*(xdat[i]-xdat[i-1]);
      sum0 += (xdat[i]-xdat[i-1]);
      dprintf(1,"%d %g %g %g\n",i,ydat[i],sum,sum0);
    }
  }
  if (Qnorm)
    sum /= sum0;
  printf("%g\n",sum);
}

