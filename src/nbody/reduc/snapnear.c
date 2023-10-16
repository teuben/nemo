/*
 *  SNAPNEAR: find nearest point in a snapshot
 *
 *   22-jul-2020    V0.1   drafted
 */

#include <stdinc.h>
#include <getparam.h>
#include <math.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <bodytransc.h>

string defv[] = {
    "in=???\n			Input file (snapshot)",
    "vals=\n                    Values of the things to compare",
    "options=x,y,z\n	        Things to compare",
    "times=all\n		Times to select snapshot",
    "VERSION=0.1\n		22-jul-2020 PJT",
    NULL,
};

string usage="find a near point in a snapshot";

#define MAXOPT    50

extern string *burststring(string,string);

void nemo_main()
{
    stream instr, tabstr;
    real   tsnap, dr, aux, d, dmin;
    real   vars[MAXOPT], bvars[MAXOPT];
    string times;
    Body *btab = NULL, *bp, *bq;
    int i, n, nbody, bits, nsep, isep, nopt, ParticlesBit, nvals;
    int imin;
    string *opt;
    rproc_body fopt[MAXOPT];

    ParticlesBit = (MassBit | PhaseSpaceBit | PotentialBit | AccelerationBit |
            AuxBit | KeyBit | DensBit | EpsBit);
    instr = stropen(getparam("in"), "r");	/* open input file */

    opt = burststring(getparam("options"),", ");
    nopt = 0;					/* count options */
    while (opt[nopt]) {				/* scan through options */
      fopt[nopt] = btrtrans(opt[nopt]);
      nopt++;
      if (nopt==MAXOPT) {
	dprintf(0,"\n\nMaximum number of options = %d exhausted\n",MAXOPT);
	break;
      }
    }
    
    nvals = nemoinpr(getparam("vals"),vars,nopt);
    if (nopt != nvals)
      error("Need %d values",nopt);
    times = getparam("times");

    get_history(instr);                 /* read history */
    for(;;) {                /* repeating until first or all times are read */
      get_history(instr);
      if (!get_tag_ok(instr, SnapShotTag))
	break;                                  /* done with work */
      get_snap_by_t(instr, &btab, &nbody, &tsnap, &bits, times);	
      if ( (bits & ParticlesBit) == 0)
	continue;                   /* skip work, only diagnostics here */

      for (bp = btab, i=0; bp < btab+nbody; bp++, i++) {     
	d = 0.0;
	for (n=0; n<nopt; n++) {
	  bvars[n] = fopt[n](bp,tsnap,i);
	  d += sqr(bvars[n]-vars[n]);
	}
	d = sqrt(d);
	if (i==0) {
	  dmin = d;
	  imin = 0;
	} else if (d < dmin) {
	  dmin = d;
	  imin = i;
	}
      }	
      bp = btab+imin;
      printf("%g %g %g ",Pos(bp)[0],Pos(bp)[1],Pos(bp)[2]);
      printf("%g %g %g ",Vel(bp)[0],Vel(bp)[1],Vel(bp)[2]);
      printf("%d %g",imin,dmin);
      printf("\n");
    }
    strclose(instr);
}
