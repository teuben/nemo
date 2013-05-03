/*
 * TABSMOOTH: smooth a table column - hanning for now 
 *          
 *   24-oct-07   Created quick & dirty               PJT
 *          12   keyword error
 *      may-13   smooth= added
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
    "smooth=0.25,0.5,0.25\n       Smoothing array",
    "VERSION=0.4\n		  2-may-2013 PJT",
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
#define MAXSM    128

local string input;			/* filename */
local stream instr;			/* input file */

local int ncol;                         /* number of columns used */
local int col[MAXCOL];			/* histogram column number(s) */

real *coldat[MAXCOL];
local int    nmax;			/* lines to use at all */
local int    npt;			/* actual number of points */

local int  nsm;
local real sm[MAXSM];                   /* smoothing array */


local void setparams(void);
local void read_data(void); 
local void smooth_data(void);
local void smooth_data_old(void);



/****************************** START OF PROGRAM **********************/

nemo_main()
{
    setparams();			/* read the parameters */
    read_data();
    smooth_data();
}

local void setparams()
{
    real sum;
    int i;

    input = getparam("in");
    ncol = nemoinpi(getparam("xcol"),col,MAXCOL);
    if (ncol < 0) error("parsing error xcol=%s",getparam("xcol"));
    
    nmax = nemo_file_lines(input,getiparam("nmax"));
    if (nmax<1) error("Problem reading from %s",input);

    instr = stropen (input,"r");

    nsm = nemoinpr(getparam("smooth"),sm,MAXSM);
    if (nsm < 0) error("smooth=%s parsing error",getparam("smooth"));
    if (nsm % 2 != 1) error("smooth=%s needs an odd number of values",getparam("smooth"));
    for (i=0, sum=0.0; i<nsm; i++) sum += sm[i];
    dprintf(0,"Smooth sum= %g\n",sum);
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


local void smooth_data_old(void)
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

local void smooth_data(void)
{
  int i,j,k,ipk;
  real sum[MAXCOL];
  int nsmh = (nsm-1)/2;  /* nsm better be odd */

  for (i=0; i<npt; i++) {
    for (j=0; j<ncol;  j++) {
      sum[j] = 0.0;
      for (k=-nsmh; k<=nsmh; k++) {
	ipk = i+k;
	if (ipk<0 || ipk>=npt) continue;
	sum[j] += coldat[j][ipk] * sm[k+nsmh];
      }
      printf("%g ",sum[j]);
    }
    printf("\n");
  }
}

