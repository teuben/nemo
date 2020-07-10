/*
 *  SNAPIPDIST: some stats on interparticle distance (see also snapstat and snapkmean)
 *
 *	24-jun-2020	V0.1 Q&D     PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <moment.h>
#include <vectmath.h>
#include <history.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>

#include <mdarray.h>

string defv[] = {		
  "in=???\n	              Input file (snapshot)",
  "times=all\n                Times of snapshot",
  "pos=f\n                    Show x,y,z as well?",
  "key=f\n                    Use key instead of ordinal ID (0...)",
  "fmt=%g\n                   Output format of distance",
  "VERSION=0.5\n	      10-jul-2020 PJT",
  NULL,
};

string usage = "some stats on interparticle distance";

string cvsid = "$Id$";

#define MAXOPT    6
#define MAXK      10

real distance(int ndim, real *x1, real *x2);

void nemo_main(void)
{
  stream instr;
  real   tsnap, ekin, etot, dr, r, rv, v, vr, vt, aux, d, d0;
  string headline=NULL, options, times;
  string fmt = getparam("fmt");
  Body *btab = NULL, *bp1, *bp2;
  bool   Qtime, Qpos, Qkey;
  int i, j, j0, k, n, nbody, bits, ParticlesBit, ndim;

  instr = stropen(getparam("in"), "r");	/* open input file */
  times = getparam("times");
  Qpos = getbparam("pos");
  Qkey = getbparam("key");
  ParticlesBit = PhaseSpaceBit ;
  
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

    for (bp1 = btab, i=0; bp1 < btab+nbody; bp1++, i++) {
      d0 = -1;
      j0 = -1;
      for (bp2 = btab, j=0; bp2 < btab+nbody; bp2++, j++) {
	if (i==j) continue;
	d = distance(NDIM, Pos(bp1), Pos(bp2));
	if (j0<0 || d < d0) {
	  d0 = d;
	  j0 = j;
	}
	dprintf(1,"%d %d %g %g\n",j,i,d,d0);
      }
      if (Qkey)
	printf("%d %d ",Key(btab+i),Key(btab+j0));
      else
	printf("%d %d ",i,j0);
      printf(fmt,sqrt(d0));
      if (Qpos)
	printf(" %g %g %g\n",Pos(btab+i)[0], Pos(btab+i)[1], Pos(btab+i)[2]);
      else
	printf("\n");
    }
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

