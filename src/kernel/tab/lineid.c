/*
 *  LINEID:   spectral line ID tool
 *
 *    13-dec-2024   0.2 simple brightest peak finder.
 *     3-jul-2025   0.3 add mode=
 *    15-jan-206    0.4 all modes now work, added simple nearest neighbor line_id
 *
 *
 *  Here's a contrived example of taking some lines, rounding them to 2 digits, and fitting them to the same linelist
 
 tabcols $NEMODAT/z_lines.list 1 | tabmath - - '%1/(1+0.001)' all format=%.2f | lineid - mode=0 linelist=$NEMODAT/z_lines.list vel=300 dv=2.5

 */

/**************** INCLUDE FILES ********************************/ 
 
#include <stdinc.h>
#include <getparam.h>
#include <extstring.h>
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
  "xunit=GHz\n           x-axis unit (GHz, km/s)",
  "minchan=3\n           Minimum channels for a peak (not implemented)",
  "maxchan=\n            Maximum channels for a peak (not implemented)",
  "vel=\n                VEL of object, if known (km/s)",
  "restfreq=\n           RESTFREQ (or RESTWAVE?) in xunits, if known",
  "linelist=\n           ASCII linelist (freq,label)",
  "dv=10\n               Slop in velocity to allow in line_id (km/s)",
  "clip=\n               Do not fit peaks below this clip level",
  "mode=1\n              0: peaks given in xcol  1: fit peak(s) in (xcol,ycol)",
  "velmode=OPT\n         OPT, RAD or REL (not implemented)",
  "VERSION=0.4\n	 15-jan-2026 PJT",
  NULL
};

string usage = "spectral line ID tool";


/**************** GLOBAL VARIABLES ************************/

local string input;				/* filename */
local stream instr;				/* input file */
local table *tptr;                              /* table */
local mdarray2 d2;                              /* data[col][row] */

local int xcol, ycol;

local real *rfreq;                              /* restfreq estimated go here */
local int nrfreq;                               /* active length of this array */

local real  *x, *y;    			        /* data from file */
local int    npt;				/* actual number of data points */

local int mode;

local bool Qvlsr, Qrest;
local real vlsr;
local real restfreq;
local real dv;
local string linelist = NULL;

#define MAXL 100
local real   lfreq[MAXL];
local string lname[MAXL];
local int    lnum = 0;

local string xunit;
local int xunit_mode = 0;

#define UNIT_MHZ 1
#define UNIT_GHZ 2
#define UNIT_A   3
#define UNIT_NM  4
#define UNIT_KMS 6
#define UNIT_Z   7

void set_params(void);
void read_data(void); 
void peak_data(void); 
void line_id(void);
void do_peak(real *xpeak, real *xerr, real *ypeak);
  
void nemo_main()
{
  set_params();
  read_data();
  peak_data();
  line_id();
}

void set_params()
{
  mode = getiparam("mode");
  input = getparam("in");             /* input table file */
  instr = stropen (input,"r");
  tptr = table_open(instr,0);
  
  xcol = getiparam("xcol");
  if (mode==0)
    ycol = 0;
  else
    ycol = getiparam("ycol");

  dv = getdparam("dv");

  if (hasvalue("linelist"))
    linelist = getparam("linelist");

  xunit = getparam("xunit");
  if (streq(xunit,"GHz"))  xunit_mode = UNIT_GHZ;
  if (streq(xunit,"km/s")) xunit_mode = UNIT_KMS;
  if (xunit_mode == 0) error("xunit %s not supported yet",xunit);

  Qvlsr = hasvalue("vel");
  Qrest = hasvalue("restfreq");
  if (Qvlsr) vlsr = getdparam("vel");
  if (Qrest) restfreq = getdparam("restfreq");
  if (Qvlsr && Qrest) {
    if (xunit_mode == UNIT_KMS) {
      warning("Mode not implemented yet");
    } else
      error("Cannot give both vel= and restfreq=");
  }

  
}

void read_data()
{
  int colnr[1+MAXCOL];
		
  dprintf(2,"Reading datafile, xcol,ycol=%d..,%d,...\n",xcol,ycol);
    
  colnr[0]  = xcol;
  colnr[1]  = ycol;
  d2 = table_md2cr(tptr, 2, colnr,0,0);
  npt = table_nrows(tptr);
  dprintf(1,"Found %d rows\n", npt);
  
  x = &d2[0][0];
  y = &d2[1][0];

  // get minmax in X and Y

  rfreq = (real *) allocate(npt * sizeof(real));
  nrfreq = 0;

  // manually parse the linelist, only look at first and second column
  // (freq,name)

  if (linelist != NULL) {
    char line[MAX_LINELEN];
    string *words;
    int nw;

    instr = stropen(linelist,"r");
    while (fgets(line, MAX_LINELEN, instr) != NULL) {
      if (line[0] == '#') continue;
      words = burststring(line," \n");
      nw = xstrlen(words,sizeof(string))-1;
      dprintf(1,"%d: %g %s\n", nw, atof(words[0]), words[1]);
      lfreq[lnum] = atof(words[0]);
      lname[lnum] = strdup(words[1]);
      lnum++;
    }
    strclose(instr);
    dprintf(0,"Found %d lines in the listlist %s\n",lnum,linelist);
  }
}

void peak_data()
{
  real xpeak, xerr, ypeak;
  real v, z, z2, rf2, sf2;
  int i;
  
  dprintf(1,"C(mks)=%g\n",c_MKS);
  dprintf(1,"xunit=%s\n",xunit);

  if (mode == 1 ) {
    //warning("Lines from peak fits in table:");    
    // for now: simple brightest peak finder
    do_peak(&xpeak, &xerr, &ypeak);
  } else if (mode == 0) {
    //warning("Direct line from table:");
    xpeak = x[0];
    ypeak = y[0];
  } else
    error("mode=%d not implemented", mode);

  if (mode == 1) {  // peak fit to (xcol,ycol)
    if (Qrest) {
      z = restfreq/xpeak - 1;
      v = z*c_MKS/1000.0;
      dprintf(1,"1a: Line at %f %s has z=%g or vel=%g km/s (cz)\n",xpeak,xunit,z,v);
    } else if (Qvlsr) {
      z = vlsr/c_MKS*1000.0;
      restfreq = xpeak * (1+z);
      dprintf(1,"1b: Line at %f, Look near freq = %f %s for a lineid; vel=%g\n", xpeak, restfreq, xunit,vlsr);
      // add
      rfreq[nrfreq++] = restfreq;
    } else {
      dprintf(0,"Peak: %g %g\n", xpeak, ypeak);
    }
  } else { // mode=0:    (xcol has freq/vel; ycol not used)
    if (xunit_mode == UNIT_KMS) {
      if (Qrest) {
	for (i=0; i<npt; i++) {
	  z2 = x[i]*1000/c_MKS;
	  sf2 = restfreq/(1+z2);
	  dprintf(1,"0a: Line at %f %s has skyfreq %f   (z=%f)\n",x[i],xunit,sf2,z2);
	}
      } else if (Qvlsr) {
	error("not implementable");
      }
    } else if (xunit_mode == UNIT_GHZ) {
      if (Qrest) {
	for (i=0; i<npt; i++) {
	  z2 = restfreq/x[i] - 1.0;
	  dprintf(1,"0b: Line at %f %s has z=%f or vel=%f km/s (cz)\n",x[i],xunit,z2,z2*c_MKS/1000.0);
	}
      } else if (Qvlsr) {
	for (i=0; i<npt; i++) {
	  z2 = vlsr*1000/c_MKS;
	  rf2 = x[i] * (1+z2);
	  dprintf(1,"0c: Line at %f %s has restfreq %f %s z=%g\n",x[i],xunit,rf2,xunit,z2);
	  // add
	  rfreq[nrfreq++] = rf2;
	}
      }
    }
  }
}

void line_id()
{
  int i, j, lmin=0;
  real dfmin, df, delta;
  char comment[3];
  
  if (nrfreq==0) {
    printf("# There are no lines to identify:\n");
    return;
  }
  
  printf("# There are %d restfreq's to identify:\n",nrfreq);

  if (lnum > 0) {
    printf("# Commented lines do not fall within dv=%g km/s\n",dv);
    printf("#    guess    ID    LINE     dFreq    deltaV  NAME\n");
  }

  for (i=0; i<nrfreq; i++) {   // loop over all frequencies found
    if (lnum == 0) {
      printf("%f\n",rfreq[i]);
      continue;
    }
    dfmin = -1.0;
    for (j=0; j<lnum; j++) {   // loop over linelist and find nearest match
      if (dfmin < 0) {
	dfmin = ABS(lfreq[j] - rfreq[i]);
	lmin = j;
	continue;
      }
      df = ABS(lfreq[j] - rfreq[i]);
      if (df < dfmin) {
	dfmin = df;
	lmin = j;
      }
    }
    delta = (1-rfreq[i]/lfreq[lmin])*c_MKS/1000.0;
    if (ABS(delta) > dv)
      strcpy(comment,"#");
    else
      strcpy(comment," ");
    
    printf("%s %f   %d %f %f %g %s\n",comment,rfreq[i], lmin, lfreq[lmin], dfmin, delta, lname[lmin]);
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

// code below is not used yet
// below here routines from ccdmom which we will eventually need for multiple peaks


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


