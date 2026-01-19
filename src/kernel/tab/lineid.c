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


/**************** COMMAND LINE PARAMETERS **********************/

string defv[] = {                /* DEFAULT INPUT PARAMETERS */
  "in=???\n              Input (table) file name",
  "xcol=1\n		 x coordinate column of spectrum or linesq", 
  "ycol=\n		 y coordinate column of spectrum",
  "dxcol=\n              If used, columns with errors in the x coordinate [not implemented]",
  "xunit=GHz\n           x-axis unit (GHz, km/s)",
  "vel=\n                VEL of object, if known (km/s)",
  "restfreq=\n           RESTFREQ (or RESTWAVE?) in xunits, if known",
  "linelist=\n           ASCII linelist (freq,label), e.g. $NEMODAT/z_lines.list",
  "dv=10\n               Slop in velocity to allow in line_id (km/s)",
  "clip=\n               Do not fit peaks below this clip level",
  "veldef=OPT\n          doppler_convention: OPT, RAD or REL",
  "VERSION=0.7\n	 18-jan-2026 PJT",
  NULL
};

string usage = "spectral line ID tool";


/**************** GLOBAL VARIABLES ************************/

local string input;				/* filename */
local stream instr;				/* input file */
local table *tptr;                              /* table */
local mdarray2 d2;                              /* data[col][row] */

local int xcol, ycol, dxcol;
local int mode = -1;                            /* 0: lines known, just need to ID   1: spectrum given */

local real *rfreq;                              /* restfreq estimated go here */
local int nrfreq;                               /* active length of this array */
local real *dvel ;                              /* doppler velocity estimates */
local int ndvel;                                /* active length of this array */

local real  *x, *y;    			        /* data from file */
local int    npt;				/* actual number of data points */

local bool Qvlsr, Qrest;
local real vlsr;
local real restfreq;
local real dv;
local string linelist = NULL;

#define MAXL 1000
local real   lfreq[MAXL];
local string lname[MAXL];
local int    lnum = 0;

#define MAXU 64
local string xunit;
local int xunit_mode = 0;
local char vlsr_unit[MAXU];
local string veldef_s;
local int veldef;

#define MAXCOL    2
#define ckms 299792.458


#define UNIT_MHZ 1
#define UNIT_GHZ 2
#define UNIT_A   3
#define UNIT_NM  4
#define UNIT_KMS 6
#define UNIT_Z   7

#define VEL_RAD  1
#define VEL_REL  2
#define VEL_OPT  3


void set_params(void);
void read_data(void); 
void peak_data(void); 
void line_id(void);
void do_peak(real *xpeak, real *xerr, real *ypeak);
void get_nu(string kv, real *value, string unit, string defunit);
  
void nemo_main()
{
  set_params();
  read_data();
  peak_data();
  line_id();
}

/*
 *     get_nu :   parse a   NUMBER,UNIT    string
 *     input:   kv         e.g.   "1,pc"
 *              defunit    e.g.   "pc"
 *     output:  value      
 *              unit
 *
 */

void get_nu(string kv, real *value, string unit, string defunit)
{
  string *sp = burststring(kv,",");
  int nsp = xstrlen(sp,sizeof(string)) - 1;
  int nv;

  *value = 1.0;
  strcpy(unit,defunit);
  if (nsp == 0) return;

  nv = nemoinpr(sp[0],value,1);
  if (nv != 1) error("error parsing %s",sp[0]);
  if (nsp == 1) return;

  strcpy(unit,sp[1]);
}

// radio to optical
inline real vrad2vopt(real vrad) {    return vrad/(1-vrad/ckms); }
// optical to radio
inline real vopt2vrad(real vopt) {    return vopt/(1+vopt/ckms); }
// optical to z
inline real vopt2z(real vopt) {       return vopt/ckms; }
// z to optical
inline real z2vopt(real z) {          return ckms*z; }
// z to radio
inline real z2vrad(real z) {          return vopt2vrad( z2vopt(z) ); }
// radio to z
inline real vrad2z(real vrad) {       return vopt2z( vrad2vopt(vrad) ); }
// relativistic to z
inline real vrel2z(real vrel) {       real b = vrel/ckms;          return sqrt((1+b)/(1-b)) - 1.0; }
// z to beta
inline real z2beta(real z) {          real z1=1.0+z;   z1=z1*z1;   return (z1-1)/(z1+1); }
// z to relativistic
inline real z2vrel(real z) {          return z2beta(z) * ckms; }



void set_params()
{
  input = getparam("in");             /* input table file */
  instr = stropen (input,"r");
  tptr = table_open(instr,0);

  xcol = getiparam("xcol");
  if (hasvalue("ycol")) {
    ycol = getiparam("ycol");
    mode = 1;
  } else {
    mode = 0;
  }

  if (hasvalue("dxcol"))
    warning("dxcol= not supported yet");

  dv = getdparam("dv");

  if (hasvalue("linelist"))
    linelist = getparam("linelist");

  xunit = getparam("xunit");
  if (streq(xunit,"GHz"))  xunit_mode = UNIT_GHZ;
  if (streq(xunit,"km/s")) xunit_mode = UNIT_KMS;
  if (xunit_mode == 0) error("xunit %s not supported yet",xunit);

  veldef_s = getparam("veldef");
  veldef = 0;
  if (strncmp(veldef_s, "opt", 1) == 0) veldef = VEL_OPT;
  if (strncmp(veldef_s, "rad", 2) == 0) veldef = VEL_RAD;
  if (strncmp(veldef_s, "rel", 2) == 0) veldef = VEL_REL;
  if (strncmp(veldef_s, "OPT", 1) == 0) veldef = VEL_OPT;
  if (strncmp(veldef_s, "RAD", 2) == 0) veldef = VEL_RAD;
  if (strncmp(veldef_s, "REL", 2) == 0) veldef = VEL_REL;
  if (veldef==0) error("veldef=%s not a valid option",veldef_s);
  dprintf(1,"veldef=%d\n",veldef);

  Qvlsr = hasvalue("vel");
  Qrest = hasvalue("restfreq");
  if (Qvlsr) {
    get_nu(getparam("vel"), &vlsr, vlsr_unit, "km/s");
    if (streq(vlsr_unit,"z")) {
      if (veldef != VEL_OPT) warning("You used z=%g, but did not specify veldef=optical convention",vlsr);
      vlsr *= c_MKS/1000;
    }
    if (veldef == VEL_REL) {
      if (vlsr > ckms) error("Bad vlsr for veldef=RELATIVISTIC");
      vlsr = z2vopt(vrel2z(vlsr));
    }
    if (veldef == VEL_RAD) {
      if (vlsr > ckms) error("Bad vlsr for veldef=RADIO");      
      vlsr = vrad2vopt(vlsr);
    }
    dprintf(0,"Vrad=%g Vrel=%g Vopt=%g z=%g\n",vopt2vrad(vlsr),z2vrel(vopt2z(vlsr)),vlsr,vopt2z(vlsr));
  }
  if (Qrest) restfreq = getdparam("restfreq");
  if (Qvlsr && Qrest) {
    if (xunit_mode == UNIT_KMS) {
      warning("Mode not implemented yet");
    } else
      error("Cannot give both vel= and restfreq= yet");
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

  // code will populate either of these two arrays
  rfreq  = (real *) allocate(npt * sizeof(real));
  dvel   = (real *) allocate(npt * sizeof(real));
  nrfreq = 0;
  ndvel  = 0;

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
      if (lnum==MAXL) {
	warning("Could only read %d lines, increase MAXL in the code", lnum);
	break;
      }
    }
    strclose(instr);
    dprintf(0,"Found %d lines in the linelist=%s\n",lnum,linelist);
  }
}

void peak_data()
{
  real xpeak, xerr, ypeak;
  real v, z, z1, z2, rf2, sf2;
  int i;
  
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
      // add vel
      dvel[ndvel++] = v;
      // add skyfreq in the restfreq slot
      rfreq[nrfreq++] = xpeak;
    } else if (Qvlsr) {
      z = vlsr/c_MKS*1000.0;
      restfreq = xpeak * (1+z);
      dprintf(1,"1b: Line at %f, Look near freq = %f %s for a lineid; vel=%g\n", xpeak, restfreq, xunit,vlsr);
      // add freq
      rfreq[nrfreq++] = restfreq;
    } else {
      dprintf(0,"Peak: %g %g\n", xpeak, ypeak);
      //rfreq[nrfreq++] = xpeak;
    }
  } else { // mode=0:    (xcol has freq/vel; ycol not used)
    if (xunit_mode == UNIT_KMS) {
      if (Qrest) {
	for (i=0; i<npt; i++) {
	  if (veldef == VEL_OPT)
	    z2 = x[i]*1000/c_MKS;
	  else if (veldef == VEL_RAD)
	    z2 =  vrad2z(x[i]);
	  else if (veldef == VEL_REL)
	    z2 =  vrel2z(x[i]);
	  sf2 = restfreq/(1+z2);
	  if (i==0) // first one is reference 
	    z1 = z2;
	  rf2 = sf2 * (1+z1);
	  dprintf(1,"0a: Line at %f %s has skyfreq %f  restfreq %f\n",x[i],xunit,sf2,rf2);
	  // add freq
	  rfreq[nrfreq++] = rf2;
	}
      } else if (Qvlsr) {
	//warning("not implementable");
	return;
      }
    } else if (xunit_mode == UNIT_GHZ) {
      if (Qrest) {
	for (i=0; i<npt; i++) {
	  z2 = restfreq/x[i] - 1.0;
	  v = z2*c_MKS/1000.0;
	  dprintf(1,"0b: Line at %f %s has z=%f or vel=%f km/s (cz)\n",x[i],xunit,z2,v);
	  // add vel
	  dvel[ndvel++] = v;
	}
      } else if (Qvlsr) {
	for (i=0; i<npt; i++) {
	  // convert to opt if not opt
	  if (veldef == VEL_OPT)
	    z2 = vlsr*1000/c_MKS;
	  else if (veldef == VEL_RAD)
	    z2 = vrad2z(vlsr);
	  else if (veldef == VEL_REL)
	    z2 = vrel2z(vlsr);
	  rf2 = x[i] * (1+z2);
	  dprintf(1,"0c: Line at %f %s has restfreq %f %s z=%g\n",x[i],xunit,rf2,xunit,z2);
	  // add freq
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

  if (ndvel > 0) {
    // need to check ndvel first, since skyfreq is hidden in rfreq for spectra
    printf("# restfreq=%g\n",restfreq);
    printf("# SkyFreq Vrad  Vrel   Vopt   z\n");
    printf("# GHz     km/s  km/s   km/s\n");
    for (i=0; i<ndvel; i++) {
      real v = dvel[i];
      printf("%g  %g %g %g  %g\n", mode==0? x[i] : rfreq[i], vopt2vrad(v), z2vrel(vopt2z(v)), v, vopt2z(v));
    }
    return;
  }
  
  if (nrfreq==0) {
    printf("# There are no lines to identify:\n");
    return;
  }
  
  printf("# There are %d restfreq's to identify:%s\n",
	 nrfreq, lnum==0 ? " (use linelist= to find their names)" : "");

  if (lnum > 0) {
    printf("# Commented lines do not fall within dv=%g km/s\n",dv);
    printf("#    guess    ID    LINE     dFreq    deltaV  NAME\n");
    printf("#    GHz            Ghz      GHz      km/s\n");
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

// Possible routines from other codes if fancier peak detection is needed
// see tabpeak and ccdmom
//  peak_spectrum
//  peak_mom
//  peak_find
//  peak_assign
