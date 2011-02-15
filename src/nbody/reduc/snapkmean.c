/*
 *  SNAPKMEAN: find kmean in a selected phase space 
 *
 *	24-sep-07	V1.0 created, at ADASS     		PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <moment.h>
#include <vectmath.h>		/* otherwise NDIM undefined */
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>

#include <mdarray.h>

string defv[] = {		
  "in=???\n	              Input file (snapshot)",
  "var=x\n                    Variables to use for coordinates",
  "k=2\n                      Number of means to find",
  "mean=-1,1\n                Initial estimates of the means",
  "times=all\n                Times of snapshot",
  "VERSION=1.0\n	      26-sep-07 pjt",
  NULL,
};

string usage = "find kmean in selected phase space of a snapshot";

string cvsid = "$Id$";

#define MAXOPT    6
#define MAXK      10

void do_kmean(int k,int ndim,int nbody,real **x, real **xmean, int *idx);

nemo_main()
{
  stream instr, tabstr;
  real   tsnap, ekin, etot, dr, r, rv, v, vr, vt, aux;
  real   varmin[MAXOPT], varmax[MAXOPT];
  real   var0[MAXOPT], var1[MAXOPT], var2[MAXOPT];
  real   mean[MAXOPT*MAXK];
  mdarray2 xmean, x;
  string headline=NULL, options, times, mnmxmode;
  Body *btab = NULL, *bp, *bq;
  bool   Qmin, Qmax, Qmean, Qsig, Qtime, scanopt();
  int i, k, n, nbody, bits, nsep, isep, ndim, ParticlesBit, *idx, nmean;
  char fmt[20],*pfmt;
  string *burststring(), *opt;
  rproc btrtrans(), fopt[MAXOPT], faux;
  
  ParticlesBit = (MassBit | PhaseSpaceBit | PotentialBit | AccelerationBit |
		  AuxBit | KeyBit);
  instr = stropen(getparam("in"), "r");	/* open input file */
  k = getiparam("k");


  opt = burststring(getparam("var"),", ");
  ndim = 0;					/* count options */
  while (opt[ndim]) {				/* scan through options */
    fopt[ndim] = btrtrans(opt[ndim]);
    ndim++;
    if (ndim==MAXOPT) {
      dprintf(0,"\n\nMaximum number of var's = %d exhausted\n",MAXOPT);
      break;
    }
  }
  for (i=0; i<ndim; i++)
    dprintf(0,"%s ",opt[i]);
  dprintf(0,"\n");

  nmean = nemoinpr(getparam("mean"),mean,MAXOPT*MAXK);
  if (nmean != ndim*k) error("not enough means given (found %d, need %d)",nmean,ndim*k);
  
  tabstr = stdout;

  times = getparam("times");

  xmean = allocate_mdarray2(k,ndim);


  get_history(instr);                 /* read history */

  for(;;) {                /* repeating until first or all times are read */
    get_history(instr);
    if (!get_tag_ok(instr, SnapShotTag))
      break;                                  /* done with work */
    get_snap(instr, &btab, &nbody, &tsnap, &bits);
    if (!streq(times,"all") && !within(tsnap,times,0.0001))
      continue;                   /* skip work on this snapshot */
    if ( (bits & ParticlesBit) == 0)
      continue;                   /* skip work, only diagnostics here */

    x = allocate_mdarray2(nbody,ndim);
    idx = (int *) allocate(nbody*sizeof(int));       /* idx[nbody] */
    
    
    for (bp = btab, i=0; bp < btab+nbody; bp++, i++) {
      idx[i] = -1;
      for (n=0; n<ndim; n++) {
	x[i][n]= fopt[n](bp,tsnap,i);
      }
    }
    for (i=0; i<k; i++)
      for (n=0; n<ndim; n++)
	xmean[i][n] = i;

    do_kmean(k,ndim,nbody,x,xmean,idx);

    free_mdarray2(x,nbody,ndim);
    free(idx);
  }
  
strclose(instr);
}

real distance(int ndim, real *x1, real *x2)
{
  int i;
  real d;

  if (ndim==1) {
    d = *x1-*x2;
    if (d<0) d = -d;
  } else {
    for (i=0, d=0.0; i<ndim; i++)
      d += (x1[i]-x2[i])*(x1[i]-x2[i]);
  }
  return d;
}

/*
 *  k               number of means
 *  ndim            dimension of space
 *  nbody           number of points
 *  x[nbody][ndim]  coordinates of points
 *  xmean[k][ndim]  coordinates of means
 *  idx[nbody]      membershop to mean (0..k-1)
 */

void do_kmean(int k, int ndim, int nbody, mdarray2 x, mdarray2 xmean, int *idx)
{
  int i,j,m,nq,iter=0;
  real d,dmin,q;
  bool again;

  dprintf(0,"iter0: ndim=%d xmean[0][0] = %g\n",ndim,xmean[0][0]);

  while(1) {
    again = FALSE;
    iter++;
    for (i=0; i<nbody; i++) {
      dmin = distance(ndim,x[i],xmean[0]);
      dprintf(1,"%d: %d %g     %g %g\n",i,idx[i],dmin,x[i][0],xmean[0][0]);
      m = 0;
      for (j=1; j<k; j++) {
	d =  distance(ndim,x[i],xmean[j]);
	if (d<dmin) {
	  dmin = d;
	  m = j;
	}
      }
      if (m != idx[i]) {
	again = TRUE;
	idx[i] = m;
      }
      dprintf(1,"%d: %d %g\n",i,idx[i],dmin);
    }

    for (j=0; j<k; j++) {
      q=0.0;
      nq = 0;
      for (i=0; i<nbody; i++) {
	if (idx[i] != j) continue;
	q += x[i][0];
	nq++;
      }
      if (nq==0) error("odd nq=0 for j=%d",j);
      xmean[j][0] = q/nq;
    }
    dprintf(0,"iterating %d  q=%g\n",iter,xmean[0][0]);    

    if (!again) break;
  } 
  

  
}
