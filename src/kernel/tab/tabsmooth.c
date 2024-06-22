/*
 * TABSMOOTH: smooth a table column - hanning for now 
 *          
 *   24-oct-07   Created quick & dirty               PJT
 *          12   keyword error
 *      may-13   smooth= added
 *   28-sep-23   0.5 filter= added
 *   29-sep-23   0.6 converted to table V2 format
 *   21-jan-24   0.8 added show= and pars=
 *
 */

/**************** INCLUDE FILES ********************************/ 

#include <stdinc.h>
#include <strlib.h>
#include <getparam.h>
#include <moment.h>
#include <yapp.h>
#include <axis.h>
#include <table.h>
#include <mdarray.h>

/**************** COMMAND LINE PARAMETERS **********************/

string defv[] = {
    "in=???\n                     Input file name",
    "xcol=1\n			  Column(s) to use (1=first)",
    "filter=0\n                   Select one of the test filters (-1,0,1,2,...)",
    "pars=\n                      Optional parameters for the smoothing filter",
    "smooth=\n                    Optional explicit smoothing array",
    "tcol=\n                      Optional independant variable ('time')",
    "show=f\n                     Show the kernel",
    "edge=0\n                     How to deal with edge effects (not implemented yet)",
    "VERSION=0.8\n		  21-jun-2024 PJT",
    NULL
};

string usage = "smooth columns of a table (Hanning, Savitzky-Golay)";


/**************** LOCAL VARIABLES *****************************/

#ifndef MAXHIST
#define MAXHIST	1024
#endif

#ifndef MAXCOL
#define MAXCOL 256
#endif

#define MAXSM    100

local string input;			/* filename */
local stream instr;			/* input file */
local table *tptr;                      /* table pointer */

local int xcol[MAXCOL+1];		/* column number(s) ; one extra is tcol= used */
local int nxcol;                        /* actual number of columns used */
local mdarray2  x;                      /* x[col][row] */
local int    npt;			/* actual number of points (rows) in table */
 
local real sm[MAXSM];                   /* smoothing array (the kernel) */
local int  nsm;                         /* actual length of smoothing array */
local bool Qt;                          /* is tcol= used ? */
local bool Qshow;                       /* show the smoothing kernel */

local real pars[MAXSM];                 /* (optional) smoothing kernel parameters */
local int  npars;                       /* actual number of parameters given */


/* Savitzky-Golay tables for fixed stepsize */

typedef struct cst {
  int n;              // number of points in the filter (needs to be odd)
  real norm;          // normalization for coeff
  real coeff[MAXSM];  // coefficients
} cst, *cstptr;
  

local cst csts[] = {
  {3,  4, { 1,  2,  1}},                  // 0:  Hanning (3pt)
  {5, 35, {-3, 12, 17, 12, -3}},          // 1:  SK2 - smooth (5pt)
  {5, 12, { 1, -8,  0,  8, -1}},          // 2:  SK2 - 1st der
  {5,  7, { 2, -1, -2, -1,  2}},          // 3:  SK2 - 2nd der
  {7,231, { 5,-30, 75,131, 75,-30, 5}},   // 4:  SK4 - smooth (7pt)
  {7,252, {22,-67,-58,  0, 58, 67,-22}},  // 5:  SK4
  // 11:   box         pars=width
  // 12:   gauss       pars=stdev,nsigma
  // 13:   trapezoidal pars=width
};

/* forward references */

local void setparams(void);
local void read_data(void);
local void build_filter(int);
local void smooth_data(void);
local void show_kernel(void);



/****************************** START OF PROGRAM **********************/

void nemo_main()
{
    setparams();
    if (Qshow)
      show_kernel();
    else {
      read_data();
      smooth_data();
    }
}

local void setparams()
{
    real sum;
    int i;
    int filter = getiparam("filter");
    int nrows, ncols;

    Qshow = getbparam("show");
     
    if (!Qshow) {
      input = getparam("in");
      nxcol = nemoinpi(getparam("xcol"),xcol,MAXCOL+1);
      if (nxcol < 0) error("parsing error xcol=%s",getparam("xcol"));
    
      Qt = hasvalue("tcol");
      if (Qt)
	xcol[nxcol] = getiparam("tcol");

      instr = stropen (input,"r");
      tptr  = table_open(instr, 0);
      nrows = table_nrows(tptr);
      ncols = table_ncols(tptr);
      dprintf(1,"Table: %d x %d\n", nrows, ncols);
    }
    npars = nemoinpr(getparam("pars"), pars, MAXSM);
 
    nsm = nemoinpr(getparam("smooth"),sm,MAXSM);
    if (nsm < 0) error("smooth=%s parsing error",getparam("smooth"));
    if (nsm == 0) 
      build_filter(filter);
    else {
      if (nsm % 2 != 1) error("smooth=%s needs an odd number of values",getparam("smooth"));
      for (i=0, sum=0.0; i<nsm; i++) sum += sm[i];
      for (i=0; i<nsm; i++) sm[i] /= sum;
    }

    for (i=0, sum=0.0; i<nsm; i++) sum += sm[i];
    dprintf(1,"Smooth sum= %g\n",sum);
}



local void read_data()
{
  if (Qt)
    x = table_md2cr(tptr, nxcol+1, xcol, 0, 0);
  else
    x = table_md2cr(tptr, nxcol,   xcol, 0, 0);
  npt = tptr->nr;
  if (nxcol == 0) nxcol = tptr->nc;  
}

local void build_filter(int filter)
{
  int i;
  
  if (filter < 0) {
    warning("passthrough");
    nsm = 1;
    sm[0] = 1;
    return;
  }
  int maxfilter = sizeof(csts)/sizeof(cst);
  if (filter < maxfilter) {
    cst *c = &csts[filter];
    nsm = c->n;
    for (i=0; i<nsm; i++)
      sm[i] = c->coeff[i]/c->norm;
  } else if (filter == 11) {              // box
    if (npars < 1) error("%d need one pars=",npars);
    nsm = (int) pars[0];
    if (nsm >  MAXSM) error("MAXSM=%d too small",MAXSM);
    for (i=0; i<nsm; i++)
      sm[i] = 1.0/nsm;
  } else if (filter == 12) {              // gauss
    if (npars < 1) error("%d need one pars=",npars);
    int stdev = (int) pars[0];
    int nsigma = (npars < 2 ? 4 : (int) pars[1]);
    nsm = 1+2*nsigma*stdev;
    if (nsm >  MAXSM) error("MAXSM=%d too small",MAXSM);
    dprintf(1,"gaussian kernel: stdev=%d nsigma=%d  nsm=%d\n", stdev, nsigma, nsm);
    double sum = 0.0;
    for (i=0; i<nsm; i++) {
      sm[i] = exp(-0.5*sqr((i-nsigma*stdev)/(1.0*stdev)));
      sum += sm[i];
    }
    for (i=0; i<nsm; i++)
      sm[i] /= sum;
  } else if (filter == 13) {              // trapezoidal
    error("trapezoidal not implemented yet");
  } else if (filter == 14) {              // moffat
    error("moffat not implemented yet");    
  } else
    error("Unknown filter=%d", filter);
}


local void show_kernel(void)
{
  for (int i=0; i<nsm; i++)
    printf("%g\n",sm[i]);
}

local void smooth_data(void)
{
  int i,j,k,ipk;
  real sum[MAXCOL];
  int nsmh = (nsm-1)/2;  /* nsm better be odd */

  for (i=0; i<npt; i++) {
    if (Qt) printf("%g ", x[nxcol][i]);
    for (j=0; j<nxcol;  j++) {
      sum[j] = 0.0;
      for (k=-nsmh; k<=nsmh; k++) {
	ipk = i+k;
	if (ipk<0 || ipk>=npt) continue;   // @todo deal with edge=
	sum[j] += x[j][ipk] * sm[k+nsmh];
      }
      printf("%g ",sum[j]);
    }
    printf("\n");
  }
}

