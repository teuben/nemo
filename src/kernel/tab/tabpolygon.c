/*
 * TABPOLYGON: select points that are either inside or outside a polygon
 *
 *       5-mar-2013   v0.1 : created for Erin -                     peter
 */

/* TODO:
 *     - lots?
 */

/**************** INCLUDE FILES ********************************/ 
 
#include <stdinc.h>	
#include <getparam.h>

#define MAXCOL  256
#define MAXCOORD 16

/**************** COMMAND LINE PARAMETERS **********************/

string defv[] = {                /* DEFAULT INPUT PARAMETERS */
    "in=???\n            Input file name (table)",
    "xcol=1\n		 x coordinate column(s)",
    "ycol=2\n		 y coordinate column(s)",
    "polygon=\n          Give pairs of X,Y that define the polygon",
    "inside=t\n          Select inside of region to be output",
    "comment=t\n         Also copy comment lines?",
    "nmax=100000\n       Hardcoded allocation space, if needed for piped data",
    "VERSION=0.1\n	 5-mar-2013 PJT",
    NULL
};

string usage = "select points that are either inside or outside a polygon";

string cvsid="$Id$";


/**************** GLOBAL VARIABLES ************************/

#define MAXPOLY    128

local string input;				/* filename */
local stream instr;				/* input file */

local int xcol, ycol;

local real  *x[MAXCOL], *y[MAXCOL]; 			/* data from file */
local int    npt;				/* actual number of data points */
local int    np;                              /* actual number of rebinned data */

local int    nmax;				/* lines to allocate */
local int    np;
local real   polygon[2*MAXPOLY];
local real   xp[MAXPOLY], yp[MAXPOLY];

local bool Qinside;
local bool Qcomment;

void setparams(void);
int pairs(string p, real *pairs, real *x, real *y, int maxp);
int inpolygon (int n, real *x, real *y, real x0, real y0);

#define MVAL            128

#ifndef MAX_LINELEN
#define MAX_LINELEN  2048
#endif




/****************************** START OF PROGRAM **********************/

nemo_main()
{
  char line[MAX_LINELEN];
  real vals[MVAL], x0,y0;
  int nv, in;
  setparams();

  instr = stropen (input,"r");
  while(1) {
    if (fgets(line,MAX_LINELEN,instr) == NULL) break;
    if (line[0] == '#')  {
      if (Qcomment)
	printf("%s",line);
      else
	continue;
    }
    line[strlen(line)-1] = '\0';
    nv = nemoinpr(line,vals,MVAL);
    if (nv<1) error("parsing error line: %s",line);
    x0 = vals[xcol-1];
    y0 = vals[ycol-1];
    in = inpolygon(np,xp,yp,x0,y0);
    if ((Qinside && in) || (!Qinside && !in))
      printf("%s\n",line);
    else
      printf("# %s\n",line);
  }
}

void setparams()
{
  input = getparam("in");             /* input table file */
  xcol = getiparam("xcol");
  ycol = getiparam("ycol");
  np = pairs(getparam("polygon"),polygon,xp,yp,2*MAXPOLY);
  
  Qinside = getbparam("inside");
  Qcomment = getbparam("comment");
}



int pairs(string p, real *pairs, real *xp, real *yp, int maxp)
{
  int n, i;

  n = nemoinpr(p,pairs,maxp);
  if (n>=0 && n<2) error("%s= needs at least 2 values",p);
  if (n%2) error("Need even number of points for pairs");
  if (n<0) error("Parsing error %s",p);
  for (i=0; i<n; i++) {
    if (i%2)
      yp[i/2] = pairs[i];
    else
      xp[i/2] = pairs[i];
  }
  n = n/2;
  for (i=0; i<n; i++) 
    dprintf(0,"polygon-%d: %g %g\n",i,xp[i],yp[i]);
  return n;
}


/*	Early version, for updates see snapplot.c or some library */
/*  The current algorithm has some flaws, may not be the fastest,
 *  (I have seen a faster one somewhere, but can't recall where),
 *  works. It's from CACM aug 1962, algorithm 112 by M. Shimrat (p434)
 *  but see also remark in CACM dec 1962, p606)
 */

/*		Code for this is now in editplot/snapplot */

int inpolygon (int n, real *x, real *y, real x0, real y0)
{
  int i,b;              /* b=0:outside b=1:inside */
        
  x[n] = x[0];        /* make sure polygon is closed */
  y[n] = y[0];

  b = 0;              /* set initially to false */
  for (i=0; i<n; i++) {
    if ( ((y0<y[i])==(y0>y[i+1])) &&
	 ((x0-x[i]-(y0-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i])) < 0 ))
      b = !b;
  }
  return b;
}
