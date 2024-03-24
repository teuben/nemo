/*
 * TABRANGE: only keep rows in which a column is in some range
 *          
 *   23-mar-2023 Q&D    - PJT
 */

/**************** INCLUDE FILES ********************************/ 

#include <stdinc.h>
#include <strlib.h>
#include <getparam.h>
#include <moment.h>
#include <yapp.h>
#include <table.h>
#include <axis.h>
#include <mdarray.h>

/**************** COMMAND LINE PARAMETERS **********************/

string defv[] = {
    "in=???\n                     Input file name",
    "xcol=1\n                     Column for X to use for range",
    "xrange=0,1\n                 Min and Max value to keep in range",
    "VERSION=0.1\n		  23-mar-2024 PJT",
    NULL
};

string usage = "keep rows where given column is in range";


/**************** GLOBAL VARIABLES *****************************/


#define MAXCOL 2

local string input;			/* filename */
local stream instr;			/* input file */
local table *tptr;                      /* table */

local int ncol;                         /* number of columns used */
local int col[MAXCOL];			/* column number(s) */

local real xrange[2];                   /* min and max of the selected column */

local int    npt;			/* actual number of points */
local mdarray2 coldat;                  /* the table data */


local void setparams(void);
local void read_data(void); 
local void range_data(void);



/****************************** START OF PROGRAM **********************/

void nemo_main(void)
{
    setparams();
    read_data();
    if (npt > 0)
      range_data();
}

local void setparams()
{
    input = getparam("in");
    ncol = nemoinpi(getparam("xcol"),col,MAXCOL);
    if (ncol != 1) error("parsing error col=%s",getparam("xcol"));
    int nx = nemoinpr(getparam("xrange"),xrange,2);
    if (nx != 2) error("parsing error xrange=%s",getparam("xrange"));

    dprintf(1,"xcol=%d xrange=%g,%g\n",col[0],xrange[0],xrange[1]);
}

local void read_data()
{
    instr = stropen (input,"r");
    tptr = table_open(instr,0);
    npt = table_nrows(tptr);
    if (npt > 0) {
      coldat = table_md2cr(tptr, ncol, col, 0, 0);
      dprintf(1,"Reading %d column(s)\n",ncol);
    }
}

local void range_data(void)
{
  int i;
  real *xdat;

  for (i=0, xdat=coldat[0]; i<npt; i++, xdat++) {
    if (*xdat < xrange[0] || *xdat > xrange[1]) continue;
    printf("%s\n",table_row(tptr, i));
  }
}
