/*
 *  LINEID:   draft
 *
 *    13-dec-2024 0.2 simple brightest peak finder.
 */

/**************** INCLUDE FILES ********************************/ 
 
#include <stdinc.h>
#include <getparam.h>
#include <axis.h>
#include <table.h>
#include <mdarray.h>
#include <moment.h>
#include <mks.h>
#include <lsq.h>

#define MAXCOL    2

/**************** COMMAND LINE PARAMETERS **********************/

string defv[] = {                /* DEFAULT INPUT PARAMETERS */
  "in=???\n            Input (table) file name",
  "xcol=1\n		 x coordinate column", 
  "ycol=2\n		 y coordinate column",
  "xunit=GHz\n           X axis unit (GHz, MHz, A, nm)",
  "minchan=3\n           Minimum channels for a peak",
  "maxchan=\n            Maximum channels for a peak",
  "vlsr=\n               VLSR of object, if known (km/s)",
  "restfreq=\n           RESTFREQ in xunits, if known",
  "linelist=\n           ASCII linelist (freq,label)",
  "rms=\n                Do not fit peaks when this RMS is reached",
  "VERSION=0.2\n	 13-dec-2024 PJT",
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

local bool Qvlsr, Qrest;
local real vlsr;
local real restfreq;
local string linelist = NULL;

#define UNIT_MHZ 1
#define UNIT_GHZ 2
#define UNIT_A   3
#define UNIT_nm  4

void setparams(void);
void read_data(void);
void peak_data(void);
void do_peak(real *xpeak, real *xerr, real *ypeak);
  
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

  Qvlsr = hasvalue("vlsr");
  if (Qvlsr) vlsr = getdparam("vlsr");
  
  Qrest = hasvalue("restfreq");
  if (Qrest) restfreq = getdparam("restfreq");

  if (hasvalue("linelist"))
    linelist = getparam("linelist");
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
  npt = table_nrows(tptr);
  dprintf(0,"Found %d channels\n", npt);
  
  x = &d2[0][0];
  y = &d2[1][0];



  // get minmax in X and Y
}

void peak_data(void)
{
  real xpeak, xerr, ypeak;
  real z;
  
  dprintf(1,"C(mks)=%g\n",c_MKS);

  // simple brightest peak finder
  do_peak(&xpeak, &xerr, &ypeak);

  if (Qrest) {
    z = 1 - xpeak/restfreq;   // @todo fix
    dprintf(0,"Line at %g has z=%g or vlsr=...\n",xpeak,z);
  } else if (Qvlsr) {
    z = vlsr/c_MKS*1000.0;
    restfreq = xpeak / (1-z);  // @todo fix
    dprintf(0,"Line at %g, Look near freq = %g for a line\n", xpeak, restfreq);
  } else {
    dprintf(0,"Peak: %g %g\n", xpeak, ypeak);
  }

}


/* from: tablsqfit */
void do_peak(real *xpeak, real *xerr, real *ypeak)
{
    real mat[(MAXCOL+1)*(MAXCOL+1)], vec[MAXCOL+1], sol[MAXCOL+1], a[MAXCOL+2];
    int i, j, k, range=1;
    int order = 2;

    for (i=1, j=0; i<npt; i++)			/* find the peak, at j */
        if (y[j] < y[i]) j = i;
    if (j==0 || j==npt-1) {			/* handle edge cases */
        warning("Peak at the edge");
        printf("%g %g\n",x[j],y[j]);
    }
    if (range==2) {
    	if (j==1 || j==npt-2) {
    	    warning("Peak too close to edge");
    	}
    }

    lsq_zero(order+1, mat, vec);
    for (i=j-range; i<=j+range; i++) {
        a[0] = 1.0;
        for (k=0; k<order; k++) {
            a[k+1] = a[k] * (x[i]-x[j]);
        }
        a[order+1] = y[i];
        lsq_accum(order+1,mat,vec,a,1.0);
    }
    lsq_solve(order+1,mat,vec,sol);
    dprintf(1,"Poly2 fit near j=%d (%g,%g)\n",j+1,x[j],y[j]);
    dprintf(1,"Peak:x,y= %g %g\n",
            x[j] - sol[1]/(2*sol[2]),
	    sol[0]-sol[1]*sol[1]/(4*sol[2]));
    *xpeak = x[j] - sol[1]/(2*sol[2]);
    *xerr = 0.0;   // TBD
    *ypeak = sol[0]-sol[1]*sol[1]/(4*sol[2]);
}

// below here routines from ccdmom which we need for multiple peaks


/*
 * peak_spectrum:
 * return location of peak for (-1,y1) (0,y2) (1,y3)
 * as determined from the 2nd order polynomial going through
 * these 3 points.
 * (also finds valleys)
 *
 * Returns a number between -0.5 and 0.5 from the center
 * because it is guarenteed y2 is a local extremum
 * (though this is NOT checked here)
 *
 * It returns an exact fit, because of the 2nd order poly,
 * could try 5 points and do a lsqfit or gaussfit....
 */

local real peak_spectrum(int n, real *spec, int p)
{
  real y1, y2, y3;

  y1 = spec[p-1];
  y2 = spec[p];
  y3 = spec[p+1];
  if (y1+y3 == 2*y2) return 0.0;
  return 0.5*(y1-y3)/(y1+y3-2*y2);
}

/*
 * @todo:    have an option to use abs() for mom2 calculations,
 *           but this assumes Qcontsub has been applied
 */

local real peak_mom(int n, real *spec, int *smask, int peak, int mom, bool Qcontsub, bool Qabs, bool Qzero)
{
  int i;
  Moment m;
  bool Qfirst = TRUE;
  real cont = 0.0;

  if (Qzero) {
    for (i=0; i<n; i++) {
      if (smask[i] && spec[i] < 0.0) smask[i] = 0;
    }
  }

  ini_moment(&m, mom, n);
  if (Qcontsub) {
    for (i=0; i<n; i++) {
      if (smask[i]==peak) {
	if (Qfirst) {
	  cont = spec[i];
	  Qfirst = FALSE;
	}
	if (spec[i] < cont) cont = spec[i];
      }
    }
  }
  dprintf(1,"peak_mom %d %d %g\n",peak,mom,cont);
  for (i=0; i<n; i++)
    if (smask[i]==peak)
      accum_moment(&m, i, spec[i]-cont);
  if (mom==0)
    return sum_moment(&m);
  else if (mom==1)
    return mean_moment(&m);
  else if (mom==2)
    return sigma_moment(&m);
  else if (mom==3)
    //return skewness_moment(&m); 
    return h3_moment(&m); 
  else if (mom==4)
    //return kurtosis_moment(&m); 
    return h4_moment(&m); 
  else
    return 0.0;
}

/* 
 * this routine can be called multiple times
 * each time it will find a peak, and then walk down the peak
 * mask[] contains 0 if not part of a peak, 1 for the first peak
 * 2 for the 2nd etc.
 *
 * During the first call, the mask[] array is set to all 0's.
 * i.e. unassigned to a peak, you will need npeak=0
 *
 * Returns: peak location (an integer)  and mask[] was modified
 *
 * TODO:    after all peaks have been find, need to resolve the issue
 *          of possibly re-assigning membership to the neighbor peak
 */


local int peak_find(int n, real *data, int *mask, int npeak)
{
  int i, ipeak=0, apeak=-1;
  real peakvalue=0, oldvalue=0;  // fool compiler

  dprintf(1,"peak_find %d\n",n);
  if (npeak==0) {               /* initialize by resetting the mask */
    for(i=0; i<n; i++)
      mask[i] = 0;              /* 0 means no peak assigned */
    return -1;
  }

  while (1) {                   /* iterate until a good peak found ? */
    ipeak++;
    dprintf(1,"iter ipeak=%d npeak=%d\n",ipeak,npeak);
    if (ipeak > npeak+3) break;

    for(i=0; i<n; i++) {        /* go over all good points and find the peak */
      if (mask[i]==0) {         /* for all good data not yet assigned to a peak */
	if (apeak<0) {          /* make sure we initialize first good value as peak */
	  peakvalue = data[i];   
	  apeak = i; 
	  continue;
	}
	if (data[i] > peakvalue) {   /* find largest values and remember where it was */
	  peakvalue = data[i];
	  apeak = i;
	}
      }
    }
    dprintf(1,"apeak=%d\n",apeak);
    if (apeak==0) {             /* never allow first point to be the peak ... */
      mask[apeak] = 0;
      apeak = -1;
      for (i=1, oldvalue=peakvalue; mask[i]==0 && i<n; i++) {   /* walk down, and mask out */
	if (data[i] > oldvalue) break;                         /* util data increase again */
	oldvalue = data[i];
	mask[i] = npeak;
      }
      continue;
    }
    if (apeak==(n-1)) {         /* ...or last point, since they have no neighbors */
      mask[apeak]= 0;
      apeak = -1;
      // compiler: peakvalue may be used uninitialized
      for (i=n-2; oldvalue=peakvalue, mask[i]==0 && i>0; i--) {  /* walk down, and mask out */
	if (data[i] > oldvalue) break;                        /* until data increase again */
	oldvalue = data[i];
	mask[i] = npeak;
      }
      // break;
      continue;
    }
    if (apeak > 0) {            /* done, found a peak, but mask down the peak */
      mask[apeak] = npeak;
      for (i=apeak+1, oldvalue=peakvalue; mask[i]==0 && i<n;  i++) {
	if (data[i] > oldvalue) break;
	oldvalue = data[i];
	mask[i] = npeak;
      }
      for (i=apeak-1, oldvalue=peakvalue; mask[i]==0 && i>=0; i-- ) {
	if (data[i] > oldvalue) break;
	oldvalue = data[i];
	mask[i] = npeak;
      }
      break;
    }

  } /* while() */
  return apeak;
}

static inline int outside(real x,int i0,int i1) {
  if(x<i0) return 0;
  if(x>i1) return 0;
  return 1;
}

local void peak_assign(int n, real *d, int *s)
{
  int i;
  real p;

  for (i=1; i<n-1; i++) {
    if (s[i] && s[i-1] && s[i]!=s[i-1]) {      /* look for neighbor peaks */
      if (d[i] < d[i-1]) {
	p = peak_spectrum(n,d,i);
	//dprintf(0,"assign-1: %d %d %d   %g   %g %g\n",i,s[i],s[i-1],p,d[i-1],d[i]);
	if (p>0) {
	  dprintf(1,"re-assign-1 @%d from %d to %d\n",i,s[i],s[i-1]);
	  s[i] = s[i-1];                       /* re-assign */
	}
      } else {
	p = peak_spectrum(n,d,i-1);
	//dprintf(0,"assign-2: %d %d %d   %g   %g %g\n",i,s[i],s[i-1],p,d[i-1],d[i]);	
	if (p<0) {
	  dprintf(1,"re-assign-2 @%d from %d to %d\n",i,s[i-1],s[i]);
	  s[i-1] = s[i];                       /* re-assign */
	}
      } 
    }
  }
}


