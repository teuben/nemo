/*
 * TABCLIP:  clips points from a table based on some smoothness criteria
 *
 *    23-nov-2002  Simple version, with just beta=
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
    "nmax=100000\n       Hardcoded allocation space, if needed for piped data",
    "VERSION=0.1\n	 23-nov-02 PJT",
    NULL
};

string usage = "clip a table";

/**************** GLOBAL VARIABLES ************************/

local string input, output;		/* filename */
local stream instr, outstr;		/* input file */

local int xcol, ycol;                   /* columns for X and Y */

local real  *x, *y; 			/* data from file */
local real  xmin,xmax,ymin,ymax;        /* min and max from the table */
local bool *ok;                         /* masking array */

local int    npt;			/* actual number of data points */
local int    nmax;			/* lines to allocate */


/****************************** START OF PROGRAM **********************/

nemo_main()
{
    setparams();

    instr = stropen(input,"r");
    outstr = stropen(output,"w");
    read_data();
    if (hasvalue("deriv")) deriv_data();
    write_data();
}

setparams()
{
    input = getparam("in");             /* input table file */
    output = getparam("out");             /* input table file */
    xcol = getiparam("xcol");
    ycol = getiparam("ycol");

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
  real d1, d2, dmin = getdparam("deriv");

  for (i=0; i<npt; i++) {
    if (i>1 ) 
      d1 = (y[i]-y[i-1])/(x[i]-x[i-1]);
    else
      d1 = 0.0;
    if (i<npt-1)
      d2 = (y[i]-y[i+1])/(x[i]-x[i+1]);
    else
      d2 = 0.0;
    d1 = ABS(d1);  
    d2 = ABS(d2);  
    if (d1 > dmin && d2 > dmin) ok[i] = FALSE;
  }
}

write_data()
{
  int i, nok=0;

  for(i=0; i<npt; i++) {
    if (ok[i]) {
      fprintf(outstr,"%g %g\n",x[i],y[i]);
      nok++;
    }
  }
  dprintf(0,"Wrote %d/%d points\n",nok,npt);
}
