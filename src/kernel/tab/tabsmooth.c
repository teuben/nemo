/*
 * TABSMOOTH: smooth a table column - hanning for now 
 *          
 *   24-oct-07   Created quick & dirty               PJT
 *                        
 * 
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
    "VERSION=0.2\n		  26-jun-2012 PJT",
    NULL
};

string usage = "(hanning) smooth columns of a table";

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


local void setparams(void);
local void read_data(void); 
local void smooth_data(void);



/****************************** START OF PROGRAM **********************/

nemo_main()
{
    setparams();			/* read the parameters */
    read_data();
    smooth_data();
}

local void setparams()
{
    input = getparam("in");
    ncol = nemoinpi(getparam("xcol"),col,MAXCOL);
    if (ncol < 0) error("parsing error xcol=%s",getparam("xcol"));
    
    nmax = nemo_file_lines(input,getiparam("nmax"));
    if (nmax<1) error("Problem reading from %s",input);

    instr = stropen (input,"r");
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


local void smooth_data(void)
{
  int i,j;

  for (j=0; j<ncol;  j++)
    printf("%g ",0.667*coldat[j][1] + 0.333*coldat[j][0]);
  printf("\n");

  for (i=1; i<npt-1; i++) {
    for (j=0; j<ncol;  j++)
      printf("%g ",0.25*coldat[j][i-1] + 0.5*coldat[j][i] + 0.25*coldat[j][i+1]);
    printf("\n");
  }
  for (j=0; j<ncol;  j++)
    printf("%g ",0.667*coldat[j][npt-2] + 0.333*coldat[j][npt-1]);
  printf("\n");
}

