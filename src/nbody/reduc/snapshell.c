/*
 * SNAPSHELL.C: compute various diagnostic properties in a set of shells
 *              
 * 
 *     13-nov-01     V1.0   derived from snapkinem                             PJT
 *     17-nov-01     V1.1   implemented two options for svar= (in UA 1020 !!)  PJT
 *     19-nov-02     V1.2   process all snapshots in input if requested        PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <moment.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>

string defv[] = {	
    "in=???\n			 Input file name (snapshot)",
    "radii=???\n                 (normalized) radii",
    "pvar=vt\n                   Variables to print statistics of",
    "svar=\n                     sorting variable if shells sorted by an expression",
    "weight=1\n			 weighting for particles",
    "axes=1,1,1\n                X,Y,Z axes for spatial spheroidal normalization",
    "stats=mean,disp,n\n         Statistics to print (mean,disp,skew,kurt,min,max,median,n)",
    "format=%g\n                 Format used for output columns",
    "normalized=t\n              Use normalized radii is svar= is used?",
    "sort=t\n                    Sort snapshot in svar, if needed",
    "first=t\n                   Process only first snapshot?",
    "VERSION=1.2\n		 19-nov-02 PJT",
    NULL,
};

string usage="compute statistics of bodyvariables in a set of shells";


Body *btab = NULL;		/* pointer to array of bodies		    */
int nbody;			/* number of bodies in array		    */
real tsnap;		        /* time associated with data		    */

rproc weight;			/* weighting function for bodies	    */
rproc pvar;
rproc svar;

vector axes;                    /* normalization radii for shells           */
bool Qrad;                      /* if radii are true radii (or governed by svar) */
bool Qnorm;                     /* svar in normalized space ? */
bool Qsort;                     /* presort snapshot in svar ? */

#define MAXRAD 10000
int nrad;
real radii[MAXRAD];

string p_format;

extern int match(string, string, int *);


string stat_options = "npt,mean,dispersion,skewness,kurtosis,min,max,median,sigma";

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

string *sel_options;
int n_sel, *n_mask;


extern string *burststring(string,string);

local void print_stat(Moment *m, bool Qhead, string name);

nemo_main()
{
    stream instr;
    rproc btrtrans();
    int i, bits, ndim;
    bool Qfirst = getbparam("first");

    instr = stropen(getparam("in"), "r");
    nrad = nemoinpd(getparam("radii"),radii,MAXRAD);
    if (nrad<0) error("Parsing rad=");
    get_history(instr);
    weight = btrtrans(getparam("weight"));
    pvar = btrtrans(getparam("pvar"));
    if (hasvalue("svar")) {
      svar = btrtrans(getparam("svar"));
      Qrad = FALSE;
    } else {
      Qrad = TRUE;
    }
    Qnorm = getbparam("normalized");
    Qsort = getbparam("sort");
    p_format = getparam("format");
    ndim = nemoinpd(getparam("axes"),axes,3);
    if (ndim != NDIM) error("Not enough values for axes=");
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
            findmoment();
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


findmoment()
{
  int i, j, irad, nviol=0;
  Body *b;
  real rad, s, smin, smax, wt, unew, uold;
  vector tmpv, pos_b, vel_b;
  Moment mq, mr, ms;
  bool Qhead = TRUE;

  ini_moment(&mq,4);
  ini_moment(&mr,4);
  ini_moment(&ms,4);

  if (!Qrad) {
    smin = (svar)(btab, tsnap, 0);
    smax = (svar)(btab+nbody-1, tsnap, nbody-1);
    if (Qnorm) {
      for (j=0; j<nrad; j++) {  /* rescale to 0..nbody for easy comparisons */
	if (radii[j] < 0.0 || radii[j] > 1.0)
	  error("Normalized radii need to be in range 0..1: %d->%g",
		j+1,radii[j]);
	radii[j] *= nbody;
      }
    } else {      
      dprintf(0,"Range svar=%s  from %g to %g\n",
	    getparam("svar"),smin,smax);
      if (smin==smax)
	error("Cannot normalize, all values for svar=%s are %g",
	      getparam("svar"),smin); 
    }

  }

  for (i = 0, b = btab, irad=0; i < nbody; ) {
    reset_moment(&mq);
    reset_moment(&mr);
    reset_moment(&ms);

    if (Qrad) {                              /* select in configuration space */
      for (j=0; i < nbody ; b++, i++, j++) {
	rad = absv(Pos(b));
	uold = rad;
	if (i) {
	  if (unew < uold) nviol++;
	  uold = unew;
	} 
	  
	dprintf(2,"%g checking %d[%g %g]\n",
		rad,irad,radii[irad],radii[irad+1]);
	if (rad >= radii[irad] && rad < radii[irad+1]) {    /* radii are ellipsoidal space */
	  accum_moment(&mr, rad, 1.0);
	  accum_moment(&mq, (pvar)(b, tsnap, i), (weight)(b,tsnap,i));
	} else if (rad < radii[irad]) {
	  dprintf(3,"Skipping %d\n",j);
	} else {
	  irad++;
	  break;
	}
      }
    } else if (Qnorm) {                      /* select in direct svar space */
      for (j=0; i<nbody; b++, i++, j++) {
	rad = absv(Pos(b));
	s = (svar)(b,tsnap,i);
	uold = s;
	if (i) {
	  if (unew < uold) nviol++;
	  uold = unew;
	} 
	dprintf(2,"%g %d checking %d[%g %g]\n",
		rad,i,irad,radii[irad],radii[irad+1]);
	if (i >= radii[irad] && i < radii[irad+1]) {       /* radii are floats 0..nbody */
	  accum_moment(&mr, rad, 1.0);
	  accum_moment(&ms, s, 1.0);
	  accum_moment(&mq, (pvar)(b, tsnap, i), (weight)(b, tsnap, i));
	} else if (i < radii[irad]) {
	  dprintf(3,"Skipping %d\n",j);
	} else {
	  irad++;
	  break;
	}
      }
    } else {                                /* select in normalized svar space */
      for (j=0; i<nbody; b++, i++, j++) {
	rad = absv(Pos(b));
	s = (svar)(b,tsnap,i);
	uold = s;
	if (i) {
	  if (unew < uold) nviol++;
	  uold = unew;
	} 
	dprintf(2,"%g checking %d[%g %g]\n",
		rad,irad,radii[irad],radii[irad+1]);
	if (s >= radii[irad] && s < radii[irad+1]) {      /* radii are in svar units */
	  accum_moment(&mr, rad, 1.0);
	  accum_moment(&ms, s, 1.0);
	  accum_moment(&mq, (pvar)(b, tsnap, i), (weight)(b, tsnap, i));
	} else if (s < radii[irad]) {
	  dprintf(3,"Skipping %d\n",j);
	} else {
	  irad++;
	  break;
	}
      }
    }
    if (n_moment(&mr)) {
      if (Qhead) {
	print_stat(&mr,Qhead,"r");
	if (!Qrad) print_stat(&ms,Qhead,"svar");
	print_stat(&mq,Qhead,"pvar");
	print_stat(0,Qhead,"");
	Qhead = FALSE;
      }
      print_stat(&mr,Qhead,"");
      if (!Qrad) print_stat(&ms,Qhead,"");
      print_stat(&mq,Qhead,"");
      print_stat(0,Qhead,"");
      Qhead = FALSE;
    }
    if (irad >= nrad) break;
  }
  if (nviol)
    warning("There were %d/%d particles in the snapshot out of sort order",
	    nviol,nbody);
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

printvec(name, vec)
string name;
vector vec;
{
    printf("%12s  %10.5f  %10.5f  %10.5f  %10.5f\n",
	   name, absv(vec), vec[0], vec[1], vec[2]);
}

