/*
 * TABCLIP:  clips points from a table based on some smoothness criteria
 *
 *    23-nov-2002  0.1 Simple version, with just beta=
 *    27-nov-2002  0.3 merged 0.2 with Rahul/Stuart's afternoon hacking
 *    22-dec-2025  0.4 ANSI C
 */

/**************** INCLUDE FILES ********************************/ 

#include <stdinc.h>	
#include <getparam.h>
#include <table.h>

#define MAXCOL  256

/**************** COMMAND LINE PARAMETERS **********************/

string defv[] = {                /* DEFAULT INPUT PARAMETERS */
    "in=???\n            Input file name (an ascii table)",
    "out=???\n           Output table",
    "xcol=1\n		 x coordinate column",
    "ycol=2\n		 y coordinate column",
    "clip=\n             Clipping if abs(derivative) exceeds this value",
    "eta=\n              Max deviation of local 2nd order polynomial allowed",
    "logic=and\n         Use AND or OR",
    "nmax=100000\n       Hardcoded allocation space, if needed for piped data",
    "nf=3\n              Max number of channels of contiguous outliers",
    "comment=t\n         Keep clipped points as commented lines?",
    "VERSION=0.5\n	 28-dec-2025 PJT",
    NULL
};

string usage = "clip points from a table based on some criteria";

/**************** GLOBAL VARIABLES ************************/

local string input, output;		/* filename */
local stream instr, outstr;		/* input file */

local int xcol, ycol;                   /* columns for X and Y */

local real *x, *y; 			/* data from file */
local real  xmin,xmax,ymin,ymax;        /* min and max from the table */
local bool *ok;                         /* masking array */
local real *d;                          /* difference array */

local int   npt;			/* actual number of data points */
local int   nmax;			/* lines to allocate */
local int   nf;

local bool Qand;                        /* logic and/or dealing with both sides */


/****************************** START OF PROGRAM **********************/

void setparams()
{
  string  logic = getparam("logic");

  input = getparam("in");               /* input table file */
  output = getparam("out");             /* output table file */
  xcol = getiparam("xcol");
  ycol = getiparam("ycol");
  if (streq(logic,"and"))
    Qand = TRUE;
  else if (streq(logic,"or"))
    Qand = FALSE;
  else 
    error("logic %s must be 'and' or 'or'",logic);
  nf = getiparam("nf");
  if (nf > 3) error("Code cannot handle nf > 3");

  nmax = nemo_file_lines(input,getiparam("nmax"));
  dprintf(1,"Allocated %d lines for table\n",nmax);
  x = (real *) allocate(sizeof(real) * (nmax+1));   /* X data */
  y = (real *) allocate(sizeof(real) * (nmax+1));   /* Y data */
}

void read_data()
{
  real *coldat[1+MAXCOL];
  int i, j, k, colnr[1+MAXCOL];
		
  dprintf (2,"Reading datafile, xcol,ycol=%d,%d,\n",xcol,ycol);
  colnr[0]  = xcol;
  coldat[0] = x;
  colnr[1]  = ycol;
  coldat[1] = y;
  
  npt = get_atable(instr,2,colnr,coldat,nmax);    /* get data */
  if (npt < 0) {
    npt = -npt;
    warning("Could only read first set of %d data",npt);
  }
  ok = (bool *) allocate(npt*sizeof(bool));     /* boolean array if we're keeping this point */
  for (i=0; i<npt; i++) ok[i] = TRUE;
  d = (real *) allocate(npt * sizeof(real));
  
  xmin = ymin =  HUGE;
  xmax = ymax = -HUGE;
  for (i=0; i<npt; i++) {     /* loop to find global min and max in X and Y */
      xmax=MAX(x[i],xmax);
      xmin=MIN(x[i],xmin);
      ymax=MAX(y[i],ymax);
      ymin=MIN(y[i],ymin);
      if (i>0 && x[i-1]>x[i]) 
	error("Data in xcol must be sorted.....");
  }
  dprintf(0,"X: min/max = %g %g\n",xmin,xmax);
  dprintf(0,"Y: min/max = %g %g\n",ymin,ymax);
}


void deriv_data()
{
  int    i;
  real d1, d2, d1a, d2a, dmin = getdparam("clip");
  real d_min, d_max;
  bool first = TRUE;

  for (i=0; i<npt; i++) {                  /* loop over all points */
    if (i>1 ) {
      d1 = (y[i]-y[i-1])/(x[i]-x[i-1]);
      if (first) {
	d_min = d_max = ABS(d1);
	first = FALSE;
      }
    } else
      d1 = dmin+1;
    if (i>2 )
      d1a = (y[i]-y[i-2]);
    else
      d1a = 0;
    if (i<npt-2)
      d2a = y[i]-y[i+2];
    else
      d2a = 0.;
    if (i<npt-1) {
      d2 = (y[i]-y[i+1])/(x[i]-x[i+1]);
      if (first) {
	d_min = d_max = ABS(d1);
	first = FALSE;
      }
    } else
      d2 = dmin+1;
    d1 = ABS(d1);  
    d2 = ABS(d2);  
    d1a = ABS(d1a);
    d2a = ABS(d2a);
    if (d1 < dmin) {
      d_min = MIN(d_min,d1);
      d_min = MIN(d_min,d2);
    }
    if (d2 < dmin) {
      d_max = MAX(d_max,d1);
      d_max = MAX(d_max,d2);}
    if (d1 > dmin || d2 > dmin || d1a > dmin || d2a > dmin) ok[i] = FALSE;
    dprintf(1,"DERIV: %d %g %g %g %g %g %g %g %g %d\n",i,y[i],y[i-1],dmin,dmin,d1,d2,d1a,d2a,ok[i]);
  }
  dprintf(0,"Min/Max for deriv = %g %g (dmin=%g)\n",d_min,d_max,dmin);
}

void deriv_data2(int nf)
{
  int  i;
  real d1, d2, d1a, d2a, clip = getdparam("clip");
  real d_min, d_max;
  bool first = TRUE;
  real sum0, sum1, sum2, sigma;

  sum0 = sum1 = sum2 = 0.0;

  for (i=1; i<npt; i++) {                  /* loop over all points */
    d1 = y[i]-y[i-1];
    d[i] = fabs(d1);
    if (d[i] > clip) continue;             /* ignore obvious outliers */
    if (first) {
      d_min = d_max = d1;
      first = FALSE;
    } else {
      d_min = MIN(d_min, d1);
      d_max = MAX(d_max, d1);
    }
    sum0 += 1.0;
    sum1 += d1;
    sum2 += d1*d1;
  }
  sigma = sum2/sum0 - sum1*sum1/(sum0*sum0);
  if (sigma <= 0.0) error("Bad sigma^2 = %g\n", sigma);
  sigma = sqrt(sigma);
  dprintf(0,"Min/Max for deriv = %g %g (with clip=%g)\n",d_min,d_max,clip);
  dprintf(0,"Sigma = %g    clip/sigma=%g\n",sigma,clip/sigma);

  i = 1;
  while (i < npt-1-nf) {
    if (d[i] > clip && d[i+2] < clip) {
      ok[i] = FALSE;
      i += 2;
      continue;
    }
    if (d[i] > clip && d[i+1] > clip && d[i+3] < clip) {
      ok[i] = FALSE;
      ok[i+1] = FALSE;
      i += 3;
      continue;
    }
    if (d[i] > clip && d[i+1] > clip && d[i+2] > clip && d[i+4] < clip) {
      ok[i] = FALSE;
      ok[i+1] = FALSE;
      ok[i+2] = FALSE;
      i += 4;
      continue;
    }
    i += 1;
  }
  warning("This code can handle up to 3 contiguous outlier channels");
}

void delta_data()
{
  int    i;
  real d1, d2, d1a, d2a, dmin = getdparam("clip");
  real d_min, d_max;
  bool first = TRUE;

  for (i=0; i<npt; i++) {                  /* loop over all points */
    if (i>1 ) {
      d1 = (y[i]-y[i-1])/(x[i]-x[i-1]);
      if (first) {
	d_min = d_max = fabs(d1);
	first = FALSE;
      }
    } else
      d1 = dmin+1;
    if (i>2 )
      d1a = (y[i]-y[i-2]);
    else
      d1a = 0;
    if (i<npt-2)
      d2a = y[i]-y[i+2];
    else
      d2a = 0.;
    if (i<npt-1) {
      d2 = (y[i]-y[i+1])/(x[i]-x[i+1]);
      if (first) {
	d_min = d_max = fabs(d1);
	first = FALSE;
      }
    } else
      d2 = dmin+1;
    d1 = ABS(d1);  
    d2 = ABS(d2);  
    d1a = ABS(d1a);
    d2a = ABS(d2a);
    if (d1 < dmin) {
      d_min = MIN(d_min,d1);
      d_min = MIN(d_min,d2);
    }
    if (d2 < dmin) {
      d_max = MAX(d_max,d1);
      d_max = MAX(d_max,d2);}
    if (d1 > dmin || d2 > dmin || d1a > dmin || d2a > dmin) ok[i] = FALSE;
    dprintf(1,"DERIV: %d %g %g %g %g %g %g %g %g %d\n",i,y[i],y[i-1],dmin,dmin,d1,d2,d1a,d2a,ok[i]);
  }
  dprintf(0,"Min/Max for deriv = %g %g (dmin=%g)\n",d_min,d_max,dmin);
}

/*
 * fit a local polynomial to the 3 points defined by 0,1,2
 * and compute the 'velocity' and 'acceleration'
 * then predict the 'position' where it should be at point 3
 */
real get_av(int i0, int i1, int i2, int i3)
{
  real t0 = x[i0];
  real t1 = x[i1]-t0;
  real t2 = x[i2]-t0;
  real t3 = x[i3]-t0;
  real x0 = y[i0];
  real x1 = y[i1]-x0;
  real x2 = y[i2]-x0;
  real x,a,v;

  a =  2*(t1*x2-t2*x1)/(t1*t2*(t2-t1));
  v = (t1*t1*x2-t2*t2*x1)/(t1*t2*(t1-t2));
  x = x0 + v*t3 + 0.5*a*t3*t3;
  return x;
}

void eta_data()
{
  int    i;
  real d1, d2, d1a, d2a, dmin = getdparam("clip");
  real eta = getdparam("eta");
  real d_min, d_max;
  bool first = TRUE;
  dprintf(0,"ETA: Using eta=%g\n",eta);

  /* can't do much about the first 3 and last 3 points .... */

  for (i=3; i<npt-3; i++) {                  /* loop over all points */
    d1 = get_av(i-3,i-2,i-1,i);
    d2 = get_av(i+3,i+2,i+1,i);
    d1 = ABS(d1-y[i]);
    d2 = ABS(d2-y[i]);
    if (Qand) {
      if (d1 > eta && d2 > eta) ok[i] = FALSE;
    } else {
      if (d1 > eta || d2 > eta) ok[i] = FALSE;
    }
    dprintf(1,"ETA: %d %g %g %g %g\n",i,x[i],y[i],d1,d2);
  }
}

void write_data()
{
  int i, nok=0;
  bool Qcomment = getbparam("comment");

  for(i=0; i<npt; i++) {
    if (ok[i]) {
      fprintf(outstr,"%g %g\n",x[i],y[i]);
      nok++;
    } else if (Qcomment)
      fprintf(outstr,"# %g %g\n",x[i],y[i]);
  }
  dprintf(0,"Wrote %d/%d points\n",nok,npt);
}



void nemo_main()
{
    setparams();

    instr = stropen(input,"r");
    outstr = stropen(output,"w");
    read_data();
    if (hasvalue("clip")) deriv_data2(nf);
    if (hasvalue("eta"))     eta_data();
    write_data();
}

