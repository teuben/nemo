/*
 * TABBIN: bin two hack
 *          
 */

/**************** INCLUDE FILES ********************************/ 

#include <stdinc.h>
#include <strlib.h>
#include <getparam.h>
#include <table.h>
#include <axis.h>

/**************** COMMAND LINE PARAMETERS **********************/

string defv[] = {
    "in=???\n                     Input file name",
    "xcol=1\n			  Column(s) to use",
#if 0    
    "nmax=100000\n                max size if a pipe",
#endif
    "VERSION=0.1\n		  2-apr-2026 PJT",
    NULL
};

string usage = "binning (by 2) hack";


/**************** SOME GLOBAL VARIABLES ************************/

#ifndef MAXHIST
#define MAXHIST	1024
#endif

#ifndef MAXCOL
#define MAXCOL 256
#endif

#define MAXCOORD 16

local string input;			/* filename */
local stream instr;			/* input file */
local table *tptr;                      /* table */

local int ncol;                         /* number of columns used */
local int col[MAXCOL];			/* column number(s) */

local int    npt;			/* actual number of points */
local mdarray2 coldat;                  /* the table data */


local void setparams(void);
local void read_data(void); 
local void bin_data(void);



/****************************** START OF PROGRAM **********************/

void nemo_main(void)
{
    setparams();			/* read the parameters */
    read_data();
    bin_data();
}

local void setparams()
{
    input = getparam("in");
    ncol = nemoinpi(getparam("xcol"),col,MAXCOL);
    if (ncol < 0) error("parsing error col=%s",getparam("col"));
    
    instr = stropen (input,"r");
    tptr = table_open(instr,0);

}



local void read_data()
{
    coldat = table_md2cr(tptr, ncol, col, 0, 0);
    dprintf(1,"Reading %d column(s)\n",ncol);
    npt = table_nrows(tptr);
}



local void bin_data(void)
{
  int i, n;
  real aver;

  n = npt%2  ? npt-1 : npt;
  
  for (i=0; i<n; i+=2) {
    aver = (coldat[0][i] + coldat[0][i+1])/2;
    printf("%d %g\n",i,aver);
  }
}

