/*
 *  LINEID:   draft
 *
 *    13-dec-2024 0.1 draft
 */

/**************** INCLUDE FILES ********************************/ 
 
#include <stdinc.h>	
#include <getparam.h>
#include <axis.h>
#include <table.h>
#include <mdarray.h>
#include <mks.h>

#define MAXCOL    2

/**************** COMMAND LINE PARAMETERS **********************/

string defv[] = {                /* DEFAULT INPUT PARAMETERS */
  "in=???\n            Input (table) file name",
  "xcol=1\n		 x coordinate column", 
  "ycol=2\n		 y coordinate column",
  "xunit=GHz\n           X axis unit (GHz, MHz, A, nm)",
  "minchan=3\n           Minimum channels for a peak",
  "maxchan=\n            Maximum channels for a peak",
  "vlsr=0\n              VLSR of object",
  "restfreq=\n           If line is know, give restfreq in xunits",
  "VERSION=0.1\n	 13-dec-2024 PJT",
  NULL
};

string usage = "lineid draft";


/**************** GLOBAL VARIABLES ************************/

local string input;				/* filename */
local stream instr;				/* input file */
local table *tptr;                              /* table */
local mdarray2 d2;                              /* data[col][row] */

local int xcol, ycol;

local real  *x, *y;    			/* data from file */
local int    npt;				/* actual number of data points */

local real vlsr;
local real restfreq;

void setparams(void);
void read_data(void);
void peak_data(void);

void nemo_main()
{
  setparams();
  read_data();
  peak_data();
}

void setparams()
{
  input = getparam("in");             /* input table file */
  instr = stropen (input,"r");
  tptr = table_open(instr,0);
  
  xcol = getiparam("xcol");
  ycol = getiparam("ycol");

  vlsr = getdparam("vlsr");
  if (hasvalue("restfreq"))
    restfreq = getdparam("restfreq");
  else
    restfreq = -1.0;
}

#define MVAL 		 64
#define MLINELEN	512

void read_data(void)
{
  int colnr[1+MAXCOL];
		
  dprintf(2,"Reading datafile, xcol,ycol=%d..,%d,...\n",xcol,ycol);
    
  colnr[0]  = xcol;
  colnr[1]  = ycol;
  d2 = table_md2cr(tptr, 2, colnr,0,0);

  x = &d2[0][0];
  y = &d2[1][0];

  npt = table_nrows(tptr);
  dprintf(0,"Found %d lines\n", npt);
}

void peak_data(void)
{
  dprintf(0,"C(mks)=%g\n",c_MKS);
  
  for (int i=0; i<npt; i++)
    printf("%d %g %g\n", i+1, x[i], y[i]);
}

