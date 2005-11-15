/*
 * SNAPSHELL.C: compute various diagnostic properties in a set of shells
 *              
 * 
 *     13-nov-01     V1.0   derived from snapkinem                             PJT
 *     17-nov-01     V1.1   implemented two options for svar= (in UA 1020 !!)  PJT
 *     19-nov-02     V1.2   process all snapshots in input if requested        PJT
 *     14-nov-05     V2.0   changed svar= to rvar=, no more sort=              PJT
 *                   V2.1   added mvar= cumulative=                            PJT
 *
 * TODO: use constant number (or mass?) fraction shells as option
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <moment.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <bodytransc.h>


string defv[] = {	
    "in=???\n			 Input file name (snapshot)",
    "radii=???\n                 (normalized) radii for shell boundaries (see also cumulative=)",
    "pvar=vt\n                   Variables to print statistics of in each shell",
    "rvar=r\n                    shell radius variable (snapshot needs sorted in this)",
    "mvar=m\n                    Mass variable if cumulative= is selected",
    "weight=1\n			 weighting for particles",
    "axes=1,1,1\n                X,Y,Z axes for spatial spheroidal normalization",
    "stats=mean,disp,npt\n       Statistics to print (mean,disp,skew,kurt,min,max,median,npt)",
    "format=%g\n                 Format used for output columns",
    "normalized=f\n              Use normalized rvar/mvar radii?",
    "cumulative=f\n              Use mvar= as cumulative in radii=",
    "first=t\n                   Process only first snapshot?",
    "rstat=f\n                   Add stats in 'r' also ?",
    "VERSION=2.1\n		 14-nov-05 PJT",
    NULL,
};

string usage="compute statistics of bodyvariables in a set of shells";

string cvsid="$Id$";

extern int match(string, string, int *);
extern string *burststring(string,string);


#ifndef MAXRAD
#define MAXRAD 10000
#endif

local int nrad;                        /* number of radii (#shells is nrad-1) */
local real radii[MAXRAD];              /* radii of shells */

local Body *btab = NULL;		/* pointer to array of bodies		    */
local int nbody;			/* number of bodies in array		    */
local real tsnap;		        /* time associated with data		    */

local rproc weight;			/* weighting function for bodies	    */
local rproc pvar;
local rproc rvar;
local rproc mvar;

local vector axes;                    /* normalization radii for shells           */
local bool Qaxes;
local bool Qnorm;                     /* rvar in normalized space ? */
local bool Qrstat;
local bool Qcumul;


local string p_format;

local string stat_options = "npt,mean,dispersion,skewness,kurtosis,min,max,median,sigma";

/*   careful: the order of the following MACRO's need to reflect those in stat_options */
#define STAT_NPT   (1<<0)
#define STAT_MEA   (1<<1)
#define STAT_DIS   (1<<2)
#define STAT_SKE   (1<<3)
#define STAT_KUR   (1<<4)
#define STAT_MIN   (1<<5)
#define STAT_MAX   (1<<6)
#define STAT_MED   (1<<7)
#define STAT_SIG   (1<<8)

local string *sel_options;
local int n_sel, *n_mask;


local void print_stat(Moment *m, bool Qhead, string name);
local void printvec(string name, vector vec);

nemo_main()
{
    stream instr;
    int i, bits, ndim;
    bool Qfirst = getbparam("first");
    Qrstat = getbparam("rstat");
    Qcumul = getbparam("cumulative");

    instr = stropen(getparam("in"), "r");
    nrad = nemoinpd(getparam("radii"),radii,MAXRAD);
    if (nrad<0) error("Parsing rad=");
    get_history(instr);
    weight = btrtrans(getparam("weight"));
    pvar = btrtrans(getparam("pvar"));
    rvar = btrtrans(getparam("rvar"));
    mvar = btrtrans(getparam("mvar"));
    Qnorm = getbparam("normalized");
    p_format = getparam("format");
    ndim = nemoinpd(getparam("axes"),axes,3);
    if (ndim != NDIM) error("Not enough values for axes=");
    if (axes[1] == 1) {
      Qaxes = TRUE;
      for (i=1; i<ndim; i++)
	if (axes[i] != axes[0]) Qaxes = FALSE;
    } else
      Qaxes = FALSE;
    sel_options = burststring(getparam("stats"),",");
    n_sel = xstrlen(sel_options,sizeof(string))-1;
    if (n_sel <= 0) error("bad stats=%g",getparam("stats"));
    n_mask = (int *) allocate(sizeof(int)*n_sel);
    for (i=0; i<n_sel; i++) {
      if (match(sel_options[i],stat_options,&n_mask[i]) < 0)
	error("No match for %s in %s",sel_options[i],stat_options);
      dprintf(1,"match %d -> %d\n",i,n_mask[i]);
    }
    dprintf(1,"Masking range %d - %d\n",STAT_NPT,STAT_MED);
    dprintf(1,"Axes: %g %g %g\n",axes[0],axes[1],axes[2]);
    while (get_snap(instr, &btab, &nbody, &tsnap, &bits)) {
      if (bits & PhaseSpaceBit) {
	reshape(1);
	shells();
      }
      if (Qfirst) break;
    }
}


int n_tot;			/* number of bodies with positive weight    */

vector cm_pos;			/* rough center of mass position	    */

int n_sum;			/* number of bodies contributing	    */

real w_sum;			/* sum of body weights			    */

vector w_pos;			/* weighted center of mass position	    */
vector w_vel;			/* weighted center of mass velocity	    */
vector w_jvec;			/* specific angular momentum		    */

matrix w_qpole;			/* weighted quadrupole moment		    */
matrix w_keten;			/* weighted kinetic energy tensor	    */

reshape(int dir)
{
    int i;
    Body *b;
    vector tmpv, iaxes;

    if (!Qaxes) return;          /* if all axes == 1, no work needed here */

    if (dir == 1) {
      for (i = 0, b = btab; i < nbody; i++, b++) {
	MULVV(tmpv, Pos(b), axes);
	SETV(Pos(b),tmpv);
      }
    } else if (dir == -1) {
      for (i=0; i<NDIM; i++)
	iaxes[i] = 1/axes[i];
      for (i = 0, b = btab; i < nbody; i++, b++) {
	MULVV(tmpv, Pos(b), iaxes);
	SETV(Pos(b),tmpv);
      }
    } else
      error("bad dir=%d",dir);


}


shells()
{
  int i, j, irad, nviol;
  Body *b;
  real rad, r, rmin, rmax, wt, unew, uold;
  vector tmpv, pos_b, vel_b;
  Moment mq, mr, ms;
  bool Qhead = TRUE;
  real *cmass, sum, totmas;

  ini_moment(&mq,4,0);
  ini_moment(&mr,4,0);
  ini_moment(&ms,4,0);

  if (Qcumul) {
    warning("new option: cumulative=t not well tested");
    cmass = (real *) allocate(nbody*sizeof(real));
    for (i=0, b=btab, sum=0.0; i<nbody; i++, b++) {
      sum += (mvar)(b, tsnap, i);
      cmass[i] = sum;
    }
    totmas = cmass[nbody-1];
    if (Qnorm) {                    /* normalize cmass to (0,1] */
      for (i=0; i<nbody; i++)
	cmass[i] /= cmass[nbody-1];
      cmass[nbody-1] = 1.0;         /* ensure it's 1 if rounding error */
    } else {
      dprintf(0,"Total for cmas = %g\n",totmas);
    }
    dprintf(0,"Cmass range: %g .. %g\n",cmass[0],cmass[nbody-1]);
  } else {
    rmin = (rvar)(btab, tsnap, 0);
    rmax = (rvar)(btab+nbody-1, tsnap, nbody-1);
    if (Qnorm) {
      for (j=0; j<nrad; j++) {  /* rescale to 0..nbody for easy comparisons */
	if (radii[j] < 0.0 || radii[j] > 1.0)
	  error("Normalized radii need to be in range 0..1: %d->%g",
		j+1,radii[j]);
	radii[j] *= nbody;
      }
    } else {      
      dprintf(0,"Range rvar=%s  from %g to %g\n",
	      getparam("rvar"),rmin,rmax);
      if (rmin==rmax)
	error("Cannot normalize, all values for svar=%s are %g",
	      getparam("rvar"),rmin); 
    }
  }
  dprintf(0,"Shell range %g .. %g\n",radii[0],radii[nrad-1]);

  /*
   * notice the quircky double i,j loops, looping over the particles
   * the j-loop is needed to detect particles before radii[0]
   * this way is probably faster than trying to index into the
   * radii[] array, but requires the particles to be sorted.
   * 
   *
   */

  for (i = 0, b = btab, irad=0, nviol=0, sum=0.0; i < nbody; ) {
    reset_moment(&mq);
    reset_moment(&mr);
    reset_moment(&ms);

    if (Qcumul) {
      for (j=0; i<nbody; b++, i++, j++) {
	rad = absv(Pos(b));
	r = (rvar)(b,tsnap,i);
	if (cmass[i] >= radii[irad] && cmass[i] < radii[irad+1]) {
	  accum_moment(&mr, rad, 1.0);
	  accum_moment(&ms, r, 1.0);
	  accum_moment(&mq, (pvar)(b, tsnap, i), (weight)(b, tsnap, i));
	} else if (cmass[i] < radii[irad]) {
	  dprintf(3,"Skipping %d\n",j);
	} else {
	  irad++;
	  break;
	}
      }
    } else if (Qnorm) {                      /* select if rvar normalized */
      for (j=0; i<nbody; b++, i++, j++) {
	rad = absv(Pos(b));
	r = (rvar)(b,tsnap,i);
	if (i>0 && r < uold) nviol++;
	uold = r;
	dprintf(2,"%g %d checking %d[%g %g]\n",
		rad,i,irad,radii[irad],radii[irad+1]);
	if (i >= radii[irad] && i < radii[irad+1]) {       /* radii are floats 0..nbody */
	  accum_moment(&mr, rad, 1.0);
	  accum_moment(&ms, r, 1.0);
	  accum_moment(&mq, (pvar)(b, tsnap, i), (weight)(b, tsnap, i));
	} else if (i < radii[irad]) {
	  dprintf(3,"Skipping %d\n",j);
	} else {
	  irad++;
	  break;
	}
      }
    } else {                                /* select if rvar space */
      for (j=0; i<nbody; b++, i++, j++) {
	rad = absv(Pos(b));
	r = (rvar)(b,tsnap,i);
	if (i>0 && r < uold) nviol++;
	uold = r;
	dprintf(2,"%g checking %d[%g %g]\n",
		rad,irad,radii[irad],radii[irad+1]);
	if (r >= radii[irad] && r < radii[irad+1]) {      /* radii are in svar units */
	  accum_moment(&mr, rad, 1.0);
	  accum_moment(&ms, r, 1.0);
	  accum_moment(&mq, (pvar)(b, tsnap, i), (weight)(b, tsnap, i));
	} else if (r < radii[irad]) {
	  dprintf(3,"Skipping %d\n",j);
	} else {
	  irad++;
	  break;
	}
      }
    }
    if (n_moment(&mr)) {       /* only print shells that have data */
      if (Qhead) {
	print_stat(&ms,Qhead,"rvar");
	print_stat(&mq,Qhead,"pvar");
	if (Qrstat) print_stat(&mr,Qhead,"r");
	print_stat(0,Qhead,"");
	Qhead = FALSE;
      }
      print_stat(&ms,Qhead,"");
      print_stat(&mq,Qhead,"");
      if (Qrstat) print_stat(&mr,Qhead,"");
      print_stat(0,Qhead,"");
    }
    if (irad >= nrad) break;
  } /* for (i=0, irad=0; ; i < nbody */
  if (nviol)
    warning("There were %d/%d particles in the snapshot out of sort order",
	    nviol,nbody);
  if (Qcumul) {
    free(cmass);
  }
}


local void print_stat(Moment *m, bool Qhead, string name)
{
  int i;

  if (m) {
    if (Qhead) {
      printf("#[%s] ",name);
      for (i=0; i<n_sel;i++) {
	switch (n_mask[i]) {
	case STAT_NPT:  
	  printf("npt ");
	  break;
	case STAT_MEA:  
	  printf("mea ");
	  break;
	case STAT_DIS:
	case STAT_SIG:    
	  printf("dis ");
	  break;
	case STAT_SKE:  
	  printf("ske ");
	  break;
	case STAT_KUR:  
	  printf("kur ");
	  break;
	case STAT_MIN:  
	  printf("min ");
	  break;
	case STAT_MAX:  
	  printf("max ");
	  break;
	case STAT_MED:  
	  printf("med ");
	  break;
	default: 	      error("Bad stats %d selected",n_mask[i]);
	}
	printf(" ");
      }
      return;
    }
    for (i=0; i<n_sel;i++) {
      switch (n_mask[i]) {
      case STAT_NPT:  
	printf("%d", n_moment(m));            
	break;
      case STAT_MEA:  
	printf(p_format, mean_moment(m));     
	break;
      case STAT_DIS:  
      case STAT_SIG:  
	printf(p_format, sigma_moment(m));    
	break;
      case STAT_SKE:  
	printf(p_format, skewness_moment(m)); 
	break;
      case STAT_KUR:  
	printf(p_format, kurtosis_moment(m)); 
	break;
      case STAT_MIN:  
	printf(p_format, min_moment(m));      
	break;
      case STAT_MAX:  
	printf(p_format, max_moment(m));      
	break;
      case STAT_MED:  
	error("cannot do median yet");        
	break;
      default: 	      error("Bad stats %d selected",n_mask[i]);

      }
      printf(" ");
    }
    printf("  ");
  } else
    printf("\n");
}

local void printvec(string name, vector vec)
{
    printf("%12s  %10.5f  %10.5f  %10.5f  %10.5f\n",
	   name, absv(vec), vec[0], vec[1], vec[2]);
}

