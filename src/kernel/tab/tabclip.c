/*
 * TABCLIP:  clips points from a table based on some smoothness criteria
 *
 *    23-nov-2002  0.1 Simple version, with just beta=
 *    27-nov-2002  0.3 merged 0.2 with Rahul/Stuart's afternoon hacking
 */

/**************** INCLUDE FILES ********************************/ 

#include <stdinc.h>	
#include <getparam.h>

#define MAXCOL  256

/**************** COMMAND LINE PARAMETERS **********************/

string defv[] = {                /* DEFAULT INPUT PARAMETERS */
    "in=???\n            Input file name (an ascii table)",
    "out=???\n           Output table",
    "xcol=1\n		 x coordinate column",
    "ycol=2\n		 y coordinate column",
    "deriv=\n            Clipping if abs(derivative) exceeds this value",
    "eta=\n              Max deviation of local 2nd order polynomial allowed",
    "logic=and\n         Use AND or OR",
    "nmax=100000\n       Hardcoded allocation space, if needed for piped data",
    "comment=f\n         Keep clipped points as commented lines?",
    "VERSION=0.3\n	 27-nov-02 PJT",
    NULL
};

string usage = "clip points from a table based on some criteria";

/**************** GLOBAL VARIABLES ************************/

local string input, output;		/* filename */
local stream instr, outstr;		/* input file */

local int xcol, ycol;                   /* columns for X and Y */

local real  *x, *y; 			/* data from file */
local real  xmin,xmax,ymin,ymax;        /* min and max from the table */
local bool *ok;                         /* masking array */

local int    npt;			/* actual number of data points */
local int    nmax;			/* lines to allocate */

local bool Qand;


/****************************** START OF PROGRAM **********************/

nemo_main()
{
    setparams();

    instr = stropen(input,"r");
    outstr = stropen(output,"w");
    read_data();
    if (hasvalue("deriv")) deriv_data();
    if (hasvalue("eta"))     eta_data();
    write_data();
}

setparams()
{
  string  logic = getparam("logic");

  input = getparam("in");             /* input table file */
  output = getparam("out");             /* input table file */
  xcol = getiparam("xcol");
  ycol = getiparam("ycol");
  if (streq(logic,"and"))
    Qand = TRUE;
  else if (streq(logic,"or"))
    Qand = FALSE;
  else 
    error("logic %s must be 'and' or 'or'",logic);

  nmax = nemo_file_lines(input,getiparam("nmax"));
  dprintf(1,"Allocated %d lines for table\n",nmax);
  x = (real *) allocate(sizeof(real) * (nmax+1));   /* X data */
  y = (real *) allocate(sizeof(real) * (nmax+1));   /* Y data */
}

read_data()
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


deriv_data()
{
  int    i;
  real d1, d2, d1a, d2a, dmin = getdparam("deriv");
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

delta_data()
{
  int    i;
  real d1, d2, d1a, d2a, dmin = getdparam("deriv");
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

eta_data()
{
  int    i;
  real d1, d2, d1a, d2a, dmin = getdparam("deriv");
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

write_data()
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
