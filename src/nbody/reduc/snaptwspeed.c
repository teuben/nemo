/*
 * SNAPTWSPEED:  TW pattern speed using particle densities (see twspeed)
 *
 *	12-mar-90  V0.1  quick hack for JCL 11:24 - 11:54         PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h> 
#include <snapshot/get_snap.c>

#include <grid.h>

string defv[] = {
  "in=???\n        Input file name (snapshot)",
  "y=-2:2:0.1\n    Sample slices in Y to take to average",
  "mass=f\n        Use mass instead of density?",
  "VERSION=0.1\n   11-mar-2011 PJT",
  NULL,
};

string usage="snapshot tw speed";

#define MAXY 1024



nemo_main()
{
  stream instr;
  real   tsnap, x, y, v, d;
  int    i, nbody, bits, iy, ny, nout;
  Body *btab = NULL, *bp;
  Grid g;
  real yrange[MAXY], sumd[MAXY], sumvd[MAXY], sumxd[MAXY];
  bool Qmass = getbparam("mass");
  
  instr = stropen(getparam("in"), "r");           /* open input file */
  ny = nemoinpr(getparam("y"),yrange,MAXY);
  if (ny<2) error("yrange syntax error or not enuf slices");
#if 0  
  inia_grid(&g, ny, yrange);
#else
  inil_grid(&g, ny-1, yrange[0], yrange[ny-1]);
#endif

  get_history(instr);			    /* accumulate data history */
  for(;;) {				 	 /* loop for all times */
    get_history(instr);                         /* for paranoidici */
    if (!get_tag_ok(instr, SnapShotTag))          /* check if done */
      break;
    get_snap(instr, &btab, &nbody, &tsnap, &bits);      /* get one */
    if ((bits & PhaseSpaceBit) == 0)
      continue;     /* if no positions -  skip */
    
    for (i=0; i<ny; i++)
      sumd[i] = sumxd[i] = sumvd[i] = 0.0;

    for (bp = btab; bp < btab+nbody; bp++) {    /* loop all bodies */
      y = Pos(bp)[1];
      iy = index_grid(&g, y);
      if (iy < 0) continue;
      d = Qmass ? Mass(bp) : Dens(bp);
      x = Pos(bp)[0];
      v = Vel(bp)[1];
      sumd[iy]  += d;
      sumvd[iy] += v*d;
      sumxd[iy] += x*d;
    }
    break; /* for now just first snapshot */
  }
  nout = 0;
  for (i=0; i<ny-1; i++) {
    if (sumd[i] > 0.0) {
      nout++;
      printf("%g %g %g %g\n",(yrange[i-1]+yrange[i])/2.0,
	     sumxd[i]/sumd[i], sumvd[i]/sumd[i], sumd[i]);
    }
  }
  if (nout==0)
    warning("No densities found in slices %s",getparam("y"));
} /* nemo_main() */
