/*
 * TABPEAK: find peaks - see also tablsqfit fit=peak
 *          
 *   28-may-2013   0.1 Created quick & dirty for ASTUTE               PJT
 *   30-may-2013   0.2 Also search for valleys
 *    1-jun-2013   0.3 Allow intensity weighted mean
 *   23-mar-2022   0.4 3pt vs. 5 pt
 *
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
#include <lsq.h>
#include <table.h>

/**************** COMMAND LINE PARAMETERS **********************/

string defv[] = {
    "in=???\n                     Input file name",
    "xcol=1\n			  X-Column",
    "ycol=2\n                     Y-column",
    "clip=0\n                     Only consider points above this",
    "pmin=3\n                     Min of points part of the peak [3=exact fit]",
    "edge=1\n                     Edge to ignore (1 or higher)",
    "valley=f\n                   Also find the valleys?",
    "mean=f\n                     Intensity weighted mean",
    "nmax=100000\n                max size if a pipe",
    "VERSION=0.5\n		  25-mar-2012 PJT",
    NULL
};

string usage = "peaks in a table";


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

local int col[2], ncol;

real *xcol, *ycol, *coldat[2];
local int    nmax;			/* lines to use at all */
local int    npt;			/* actual number of points */
local int    pmin;
local int    edge;
real clip;
bool  Qvalley, Qmean;

local void setparams(void);
local void read_data(void); 
local void peak_data(void);
local void peak_fit(void);
local void mean_data(void);



/****************************** START OF PROGRAM **********************/

void nemo_main()
{
    setparams();			/* read the parameters */
    read_data();
    if (Qmean)
      mean_data();
    else if (pmin == 3)
      peak_data();
    else
      peak_fit();
}

local void setparams()
{
    input = getparam("in");
    col[0] = getiparam("xcol");
    col[1] = getiparam("ycol");
    clip = getrparam("clip");
    pmin = getiparam("pmin");
    if (pmin < 3) error("pmin=%d too small, 3 or higher",pmin);
    edge = getiparam("edge");
    if (edge < 1) error("edge=%d too small, 1 or higher",edge);
    Qvalley = getbparam("valley");
    Qmean = getbparam("mean");
    if (Qmean && Qvalley) warning("Valley fitting not supported in mean mode");
    
    nmax = nemo_file_lines(input,getiparam("nmax"));
    if (nmax<1) error("Problem reading from %s",input);

    instr = stropen (input,"r");
}



local void read_data()
{
    int   i,j,k;
    
    ncol = 2;
    dprintf(0,"Reading %d column(s)\n",ncol);
    xcol = (real *) allocate(sizeof(real)*nmax);
    ycol = (real *) allocate(sizeof(real)*nmax);
    coldat[0] = xcol;
    coldat[1] = ycol;

    npt = get_atable(instr,ncol,col,coldat,nmax);        /* read it */
    if (npt == -nmax) {
    	warning("Could only read %d data",nmax);
    	npt = nmax;
    }
}


local void peak_data(void)
{
  int i,j;
  real mat[9], vec[3], sol[3], a[4];

  /* loop over all interior points and find peaks or valleys, fit local polynomial */

  for (i=edge; i<npt-edge; i++) {
    if (            (ycol[i]> clip && ycol[i]>ycol[i-1] && ycol[i]>ycol[i+1]) ||
         (Qvalley && ycol[i]<-clip && ycol[i]<ycol[i-1] && ycol[i]<ycol[i+1]) ) {
      lsq_zero(3,mat,vec);
      for (j=i-1; j<=i+1; j++) {
	a[0] = 1.0;
	a[1] = (xcol[j]-xcol[i]);
	a[2] = (xcol[j]-xcol[i]) * a[1];
	a[3] = ycol[j];
	lsq_accum(3,mat,vec,a,1.0);
      }
      lsq_solve(3,mat,vec,sol);
      dprintf(1,"Poly2 fit near i=%d (%g,%g)   %g %g %g\n",i+1,xcol[i],ycol[i],sol[0],sol[1],sol[2]);
      printf("%f %f \n",xcol[i] - sol[1]/(2*sol[2]),
	     sol[0]-sol[1]*sol[1]/(4*sol[2]));
    } 
  }
}

local void peak_fit(void)
{
  int i,j,p;
  real mat[9], vec[3], sol[3], a[4];
  int pedge = (pmin-1)/2;

  dprintf(0,"pedge=%d\n",pedge);

  /* loop over all interior points and find peaks or valleys, fit local polynomial */

  for (i=pedge; i<npt-pedge; i++) {
    if (            (ycol[i]> clip && ycol[i]>ycol[i-1] && ycol[i]>ycol[i+1]) ||
         (Qvalley && ycol[i]<-clip && ycol[i]<ycol[i-1] && ycol[i]<ycol[i+1]) ) {

      lsq_zero(3,mat,vec);
      for (j=i-1; j<=i+1; j++) {
	a[0] = 1.0;
	a[1] = (xcol[j]-xcol[i]);
	a[2] = (xcol[j]-xcol[i]) * a[1];
	a[3] = ycol[j];
	lsq_accum(3,mat,vec,a,1.0);
      }
      
      // found a peak, or valley, now ensure we get at least pmin points,
      // constantly requiring to go downhill/uphill
      p = 3;
      if (ycol[i-2] < ycol[i-1]) {
	j = i-2;
	a[0] = 1.0;
	a[1] = (xcol[j]-xcol[i]);
	a[2] = (xcol[j]-xcol[i]) * a[1];
	a[3] = ycol[j];
	lsq_accum(3,mat,vec,a,1.0);	
	p++;
      }
      if (ycol[i+2] < ycol[i+1]) {
	j = i+2;
	a[0] = 1.0;
	a[1] = (xcol[j]-xcol[i]);
	a[2] = (xcol[j]-xcol[i]) * a[1];
	a[3] = ycol[j];
	lsq_accum(3,mat,vec,a,1.0);	
	p++;
      }
      if (p < pmin) continue;
      
      lsq_solve(3,mat,vec,sol);
      dprintf(1,"Poly2 fit near i=%d (%g,%g)   %g %g %g\n",i+1,xcol[i],ycol[i],sol[0],sol[1],sol[2]);
      printf("%f %f %d\n",xcol[i] - sol[1]/(2*sol[2]),
	     sol[0]-sol[1]*sol[1]/(4*sol[2]),p);
    } 
  }
}

local void mean_data(void)
{
  int i,i0,i1,ipeak;
  real peak, sum0, sum1, sum2, xmean, xsig;


  /* find first occurence > clip */
  peak = ycol[0];
  i0   = -1;
  for (i=0; i<npt; i++) {
    if (ycol[i]>peak) { 
      peak = ycol[i];
      ipeak = i;
    }
    if (ycol[i]>clip) {
      i0 = i;
      break;
    }
  }
  dprintf(1,"First peak %g at %d\n",peak,i0);
  if (i0 < 0) error("No data above clip=%g, peak %g at %d",clip,peak,ipeak);

  while (1) {                         /* enter loop searching for sections > clip */

    sum0 = sum1 = sum2 = 0.0;
    peak = ycol[i0];                 /* first point is guarenteed above clip */
    for (i=i0; i<npt; i++) {         
      if (ycol[i] < clip) {
	i0 = i;
	break;
      }
      if (ycol[i]>peak) peak = ycol[i];
      sum0 += ycol[i];
      sum1 += ycol[i]*xcol[i];
      sum2 += ycol[i]*xcol[i]*xcol[i];
    }
    xmean = sum1/sum0;
    xsig = sum2/sum0 - xmean*xmean;
    if (xsig>0) xsig=sqrt(xsig);
    printf("%f %f %f %d\n",xmean,xsig,peak,i0);

    /* search for next peak , i0 is known to have < clip */
    for (i=i0; i<npt; i++) {
      if (ycol[i] > clip) {
	i0 = i;
	break;
      }
    }
    if (i >= npt-1) break;
  }
  
  
}


