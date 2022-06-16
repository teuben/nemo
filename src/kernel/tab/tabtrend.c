/*
 * TABTREND: difference rows from previous rows
 *          
 *   12-may-05   Created quick & dirty for CARMA               PJT
 *   15-jun-2016 cumul= option
 *   23-may-2020 orig=T/F option
 *   17-mar-2021 add first= option
 *    6-jun-2022 allow delta > 1
 *                        
 * 
 * TODO:
 *     option to do dual-pass to subtract the mean before computing
 *     the higher order moments - needed for accuracy
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
    "xcol=1\n			  Column(s) to use",
#if 0    
    "nmax=100000\n                max size if a pipe",
#endif
    "cumul=f\n                    cumulative instead?",
    "delta=1\n                    How far away the neighbor",
    "orig=f\n                     show original column as well?",
    "first=f\n                    add first row?",
    "VERSION=0.4\n		  17-mar-2021 PJT",
    NULL
};

string usage = "difference rows from previous rows, or cumulate them";


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
local int delta;                        /* neighbor distance */

local int    npt;			/* actual number of points */
local mdarray2 coldat;                  /* the table data */

local bool   Qcumul;
local bool   Qorig;
local bool   Qfirst;


local void setparams(void);
local void read_data(void); 
local void trend_data(void);
local void cumul_data(void);



/****************************** START OF PROGRAM **********************/

void nemo_main(void)
{
    setparams();			/* read the parameters */
    read_data();
    if (Qcumul)
      cumul_data();
    else
      trend_data();
}

local void setparams()
{
    input = getparam("in");
    ncol = nemoinpi(getparam("xcol"),col,MAXCOL);
    if (ncol < 0) error("parsing error col=%s",getparam("col"));
    delta = getiparam("delta");
    
    instr = stropen (input,"r");
    tptr = table_open(instr,0);

    Qcumul = getbparam("cumul");
    Qorig = getbparam("orig");
    Qfirst = getbparam("first");
}



local void read_data()
{
    coldat = table_md2cr(tptr, ncol, col, 0, 0);
    dprintf(1,"Reading %d column(s)\n",ncol);
    npt = table_nrows(tptr);
}


local void trend_data(void)
{
  int i,j;

  for (j=0; j<ncol;  j++) {  
    if (Qfirst) {
      if (Qorig)
	printf("%g %g ",coldat[j][0],coldat[j][0]);
      else
	printf("%g ",coldat[j][0]);
      printf("\n");       
    }
  }

  for (i=delta; i<npt; i++) {
    for (j=0; j<ncol;  j++) {
      if (Qorig)
	printf("%g %g ",coldat[j][i]-coldat[j][i-delta],coldat[j][i-delta]);
      else
	printf("%g ",coldat[j][i]-coldat[j][i-delta]);
      printf("\n");
    }
  }
}

local void cumul_data(void)
{
  int i,j;
  real sum[MAXCOL];

  for (i=0; i<npt; i++) {
    if (i==0)
      for (j=0; j<ncol;  j++)
	sum[j] = -coldat[j][0];
	  
    for (j=0; j<ncol;  j++) {
	sum[j] += coldat[j][i];
	printf("%g ",sum[j]);
    }
    printf("\n");
  }
}

