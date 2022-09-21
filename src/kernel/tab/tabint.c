/* TABINT:   integrate a sorted table
 *
 *
 *      13-may-05    Q&D first version, cloned off tabfilter
 *      26-may-05    allow step<0 for logsteps (xmin must be > 0)
 *      25-apr-22    convert to use the new table V2
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
  "step=\n          Integration step if resampling used (<0 for logarithmic steps)",
  "normalize=f\n    Normalize integral",
  "cumulative=f\n   Show accumulation of integral",
  "VERSION=0.5\n    20-sep-2022 PJT",
  NULL,

};

string usage="integrate a sorted table";



extern int minmax(int, real *, real *, real *);


void nemo_main()
{
  int colnr[2], i, n, nmax, nsteps;
  real *xdat, *ydat, xmin, xmax, ymin, ymax, zmin, zmax;
  real x, y, z, s, xold, yold, dx, dz, sum, sum0, *sdat;
  string spectrum = getparam("in");
  bool Qnorm = getbparam("normalize");
  bool Qcum = getbparam("cumulative");

  
  /* read the data */

#if 1
  // #ifdef TABLE2
  tableptr t = table_open(stropen(spectrum,"r"),0);
  n = nmax = table_nrows(t);
  int ncols = table_ncols(t);
  dprintf(1,"%s has %d x %d table\n", spectrum, n, ncols);
  colnr[0] = getiparam("xcol");
  colnr[1] = getiparam("ycol");
  mdarray2 d = table_md2cr(t,2,colnr,0,0);
  xdat = &d[0][0];
  ydat = &d[1][0];
#else
  real *coldat[2];
  nmax = nemo_file_lines(spectrum,MAXLINES);
  xdat = coldat[0] = (real *) allocate(nmax*sizeof(real));
  ydat = coldat[1] = (real *) allocate(nmax*sizeof(real));
  colnr[0] = getiparam("xcol");
  colnr[1] = getiparam("ycol");
  stream instr = stropen(spectrum,"r");
  n = get_atable(instr,2,colnr,coldat,nmax);
  strclose(instr);
#endif
  
  /* figure out some min/max */
  minmax(n,xdat,&xmin,&xmax);
  minmax(n,ydat,&ymin,&ymax);
  dprintf(1,"X range: %g : %g\n",xmin,xmax);
  dprintf(1,"Y range: %g : %g\n",ymin,ymax);

  sum = sum0 = 0.0;
  nsteps = 0;
  if (hasvalue("step")) {     /* resample the function using splines */
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
    } else if (xmin > 0 || xmax < 0) {     /* dx < 0: log steps */
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
    dx = 0.0;
    for (i=1; i<n; i++) {
      sum += 0.5*(ydat[i]+ydat[i-1])*(xdat[i]-xdat[i-1]);
      sum0 += (xdat[i]-xdat[i-1]);
      if (Qcum) printf("%g %g %g\n",xdat[i],sum,sum/sum0);
    }
  }
  dprintf(1,"xmin=%g xmax=%g dx=%g sum=%g sum0=%g nsteps=%d\n",
	  xmin,xmax,dx,sum,sum0,nsteps);
  if (Qnorm)
    sum /= sum0;
  printf("%s%g\n",  Qcum ? "# " : "",  sum);
}

