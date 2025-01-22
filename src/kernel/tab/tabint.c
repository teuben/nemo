/* TABINT:   integrate a sorted table
 *
 *
 *      13-may-05    Q&D first version, cloned off tabfilter
 *      26-may-05    allow step<0 for logsteps (xmin must be > 0)
 *      25-apr-22    convert to use the new table V2
 *      22-feb-23    add xmin,xmax,scale keywords
 *      31-jul-23    add mom=
 *       6-apr-24    cleanup
 *      17-jan-25    add higher moments 3,4
 */

#include <stdinc.h> 
#include <getparam.h>
#include <spline.h>
#include <extstring.h>
#include <table.h>
#include <mdarray.h>

string defv[] = {
  "in=???\n         Input table file",
  "xcol=1\n         Column with X coordinate (must be sorted in that column)",
  "ycol=2\n         Column with Y coordinate of function",
  "xmin=\n          value if data below xmin to be discarded",
  "xmax=\n          value if data above xmax to be discarded",
  "step=\n          (override) Integration step, or if resampling used (<0 for logarithmic steps)",
  "normalize=f\n    Normalize integral",
  "cumulative=f\n   Show accumulation of integral",
  "scale=1\n        Scale factor to apply to integral",
  "mom=0\n          0=flux 1=weighted mean 2=dispersion 3=skewness 4=kurtosis",
  "VERSION=1.0\n    17-jan-2025 PJT",
  NULL,

};

string usage="integrate or compute moments of a (sorted) table";


extern int minmax(int, real *, real *, real *);
local void reverse(int n, real *x, real *y);
local void sortxy(int n, real *x, real *y);

void nemo_main()
{
  int colnr[2], i, n, nmax, nsteps;
  real *xdat, *ydat, xmin, xmax, ymin, ymax, zmin, zmax;
  real x, y, z, s, xold, yold, dx, dz, sum, sum0, *sdat;
  string spectrum = getparam("in");
  bool Qnorm = getbparam("normalize");
  bool Qcum = getbparam("cumulative");
  bool Qmin = hasvalue("xmin");
  bool Qmax = hasvalue("xmax");
  real scale = getrparam("scale");
  int mom = getiparam("mom");
  
  /* read the data */

  tableptr t = table_open(stropen(spectrum,"r"),0);
  n = nmax = table_nrows(t);
  int ncols = table_ncols(t);
  dprintf(1,"%s has %d x %d table\n", spectrum, n, ncols);
  colnr[0] = getiparam("xcol");
  colnr[1] = getiparam("ycol");
  mdarray2 d = table_md2cr(t,2,colnr,0,0);
  xdat = &d[0][0];
  ydat = &d[1][0];

  /* reverse arrays if not sorted properly */
  if (xdat[0] > xdat[1])
    reverse(n, xdat, ydat);

  /* check if array is now properly sorted */
  int nbad = 0;
  for (i=1; i<n; i++)
    if (xdat[i] < xdat[i-1]) nbad++;
  if (nbad > 0) sortxy(n, xdat, ydat);
  
  /* figure out an appropriate min/max */
  minmax(n,xdat,&xmin,&xmax);
  minmax(n,ydat,&ymin,&ymax);
  dprintf(1,"X range: %g : %g\n",xmin,xmax);
  dprintf(1,"Y range: %g : %g\n",ymin,ymax);
  if (Qmin || Qmax) {
    if (Qmin) xmin = getrparam("xmin");
    if (Qmax) xmax = getrparam("xmax");
    dprintf(1,"X range: %g : %g   (reset)\n",xmin,xmax) ;
    int imin = n;
    int imax = 0;
    for (i=0; i<n; i++) {
      if (Qmin && xdat[i] < xmin) imin = i;
      if (Qmax && xdat[i] < xmax) imax = i;
    }
    dprintf(1,"New range %d - %d   %g - %g\n", imin, imax, xdat[imin], xdat[imax]);
    /* new data */
    n = imax - imin + 1;
    xdat = &xdat[imin];
    ydat = &ydat[imin];
  }

  sum = sum0 = 0.0;
  nsteps = 0;
  if (hasvalue("step")) {     /* resample the function using splines ? */
    dx = getdparam("step");

    /* setup a spline interpolation table into the function */
    sdat = (real *) allocate(sizeof(real)*n*3);
    spline(sdat,xdat,ydat,n);
    if (dx > 0) {                        /* dx > 0 : linear steps from xmin->xmax */
        yold = ydat[0];
      xold = xmin;
      for (x=xmin+dx; x<xmax; x += dx) { /* loop over the sampled interval */
	nsteps++;
	y = seval(x,xdat,ydat,sdat,n);
	sum  += 0.5*(yold+y)*(x-xold);
	sum0 += (x-xold);
	yold = y;
	xold = x;
      }
      y = ydat[n-1];
      sum  += 0.5*(yold+y)*(xmax-xold);     /* it's possible we never */
      sum0 += (xmax-xold);                  /* closed the interval */
    } else if (xmin > 0 || xmax < 0) {      /* dx < 0: log steps */
      dz = -dx;
      s = xmin > 0 ? 1.0 : -1.0;
      zmin = log10(s*xmin);
      zmax = log10(s*xmax);
      yold = ydat[0];
      xold = xmin;
      for (z=zmin+dz; z<zmax; z += dz) {
	nsteps++;
	x = s*pow(10.0,z);
	y = seval(x,xdat,ydat,sdat,n);
	sum  += 0.5*(yold+y)*(x-xold);
	sum0 += (x-xold);
	yold = y;
	xold = x;
      }
      y = ydat[n-1];
      sum  += 0.5*(yold+y)*(xmax-xold);     /* it's possible we never */
      sum0 += (xmax-xold);                  /* closed the interval */
    } else
      error("Cannot handle log steps with dx=%g xmin=%g xmax=%g",dx,xmin,xmax);
  } else {  /* use the datapoints itself */
    dx = xdat[1]-xdat[0];
    real dxmin = dx;       // also keep track of any variable step size
    real dxmax = dx;
    if (mom == 0) {
      if (Qcum) printf("# X SUM SUM/SUM0\n");
      for (i=1; i<n; i++) {   // integration loop using trapezoid rule
	dx = xdat[i]-xdat[i-1];
	sum  += 0.5*(ydat[i]+ydat[i-1])*dx;
	sum0 += dx;
	if (Qcum) printf("%g %g %g\n",xdat[i],sum,sum/sum0);
	dxmin = MIN(dxmin,dx);
	dxmax = MAX(dxmax,dx);
      }
      dprintf(1,"dx min/max: %g %g\n",dxmin,dxmax);
    } else {
      // this is for mom > 0
      double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0, retval=0.0;
      
      for (i=0; i<n; i++) {
	sum0 +=  ydat[i];
	sum1 +=  ydat[i]*xdat[i];
	sum2 +=  ydat[i]*xdat[i]*xdat[i];
      }
      sum1 /= sum0;
      sum2 /= sum0;
      dprintf(1,"dx=%g\n", xdat[1]-xdat[0]);
      dprintf(1,"sum0=%g\n",      sum0);
      dprintf(1,"sum1/sum0=%g\n", sum1);
      dprintf(1,"sum2/sum0=%g\n", sum2);
      sum2 = sqrt(sum2 - sum1*sum1);
      dprintf(1,"dispersion=%g\n", sum2);
      if (mom==1) retval=sum1;
      if (mom==2) retval=sum2;
      if (mom>2) warning("new unchecked feature");
      if (mom==3) {
	// skewness is 0 for symmetric
	sum3 = 0.0;
	for (i=0; i<n; i++)
	  sum3 += ydat[i]*qbe(xdat[i]-sum1);
	retval = sum3/qbe(sum2)/sum0;
      }
      if (mom==4) {
	// kurtosis is 0 for gaussian,   < 0 
	sum4 = 0.0;
	for (i=0; i<n; i++)
	  sum4 += ydat[i]*sqr(sqr(xdat[i]-sum1));
	retval = -3.0 + sum4/sqr(sqr(sum2))/sum0;
      }
      printf("%g\n",retval);
      return;
    }
  }
  dprintf(1,"xmin=%g xmax=%g dx=%g sum=%g sum0=%g nsteps=%d\n",
	  xmin,xmax,dx,sum,sum0,nsteps);
  if (Qnorm)
    sum /= sum0;
  printf("%s%g\n",  Qcum ? "# " : "",  sum * scale);
}

local void reverse(int n, real *x, real *y)
{
  real t;
  for (int i = 0; i<n/2; i++) {
    t = x[i];  x[i] = x[n-i-1];  x[n-i-1] = t;
    t = y[i];  y[i] = y[n-i-1];  y[n-i-1] = t;
  }
}

local void sortxy(int n, real *x, real *y)
{
  error("Your data were not sorted, and sortxy is not implemented yet");
}
