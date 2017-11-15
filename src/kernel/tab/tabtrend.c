/*
 * TABTREND: difference rows from previous rows
 *          
 *   12-may-05   Created quick & dirty for CARMA               PJT
 *   15-jun-2016 cumul= option
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
#include <axis.h>
#include <mdarray.h>

/**************** COMMAND LINE PARAMETERS **********************/

string defv[] = {
    "in=???\n                     Input file name",
    "xcol=1\n			  Column(s) to use",
    "nmax=100000\n                max size if a pipe",
    "cumul=f\n                    cumulative instead?",
    "VERSION=0.2\n		  15-jun-2016 PJT",
    NULL
};

string usage = "difference rows from previous rows, or cumulate them";

string cvsid = "$Id$";

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

local int ncol;                         /* number of columns used */
local int col[MAXCOL];			/* histogram column number(s) */

real *coldat[MAXCOL];
local int    nmax;			/* lines to use at all */
local int    npt;			/* actual number of points */

local bool   Qcumul;


local void setparams(void);
local void read_data(void); 
local void trend_data(void);
local void cumul_data(void);



/****************************** START OF PROGRAM **********************/

nemo_main()
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
    
    nmax = nemo_file_lines(input,getiparam("nmax"));
    if (nmax<1) error("Problem reading from %s",input);

    instr = stropen (input,"r");

    Qcumul = getbparam("cumul");
}



local void read_data()
{
    int   i,j,k;
    
    dprintf(0,"Reading %d column(s)\n",ncol);
    for (i=0; i<ncol; i++)
      coldat[i] = (real *) allocate(sizeof(real)*nmax);
    npt = get_atable(instr,ncol,col,coldat,nmax);        /* read it */
    if (npt == -nmax) {
    	warning("Could only read %d data",nmax);
    	npt = nmax;
    }
}


local void trend_data(void)
{
  int i,j;

  for (i=1; i<npt; i++) {
    for (j=0; j<ncol;  j++)
      printf("%g ",coldat[j][i]-coldat[j][i-1]);
    printf("\n");
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

