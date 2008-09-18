/*
 *  SNAPPOT:    add a potential force/acc to a snapshot
 *
 *  25-mar-05   Created         Peter Teuben
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

#include <potential.h>


string defv[] = {
  "in=???\n             Input file (snapshot)",
  "out=???\n            Output file (snapshot)",
  "potname=???\n        name of potential(5)",
  "potpars=\n           parameters to potential",
  "potfile=\n           optional filename to potential",
  "times=all\n          Which times to work on",
  "VERSION=1.0\n	4-apr-97 pjt",
  NULL,
};

string usage="add analytical potentials/forces to an N-body system";

string cvsid="$Id$";

#define TIMEFUZZ	0.0001	/* tolerance in time comparisons */

void nemo_main(void)
{
  stream instr, outstr;
  real   tsnap;
  string times;
  Body *btab = NULL, *bp;
  int nbody, bits;
  potproc_real pot;
  real ome,ome2,half_ome2,two_ome;

  vector lacc,lpos;
  real   lphi;
  int    ndim=NDIM;


  times = getparam("times");

  pot = get_potential (getparam("potname"),
		       getparam("potpars"), 
		       getparam("potfile"));
  ome = get_pattern();     /* pattern speed first par of potential */
  ome2 = ome*ome;
  half_ome2 = 0.5 * ome2;
  two_ome = 2.0 * ome;

					/* open files */
  instr = stropen(getparam("in"), "r");
  outstr = stropen(getparam("out"), "w");

  get_history(instr);		/* get data history from input file */
  put_history(outstr);	/* write "    "     to output file */
  for(;;) {			/* do the work in an infinite loop */
    get_history(instr);     /* skip any unneeded history/headline stuff */
    if (!get_tag_ok(instr, SnapShotTag))
      break;			/* done with work in loop */
    get_snap(instr, &btab, &nbody, &tsnap, &bits);
    if ((bits & MassBit) == 0 && (bits & PhaseSpaceBit) == 0) continue;

    if ((bits & TimeBit) == 0)
      tsnap = 0.0;
    else if (!streq(times,"all") && !within(tsnap, times, TIMEFUZZ))
      continue;		/* however skip this snapshot */
    dprintf (1,"Snapshot time=%f shifting\n",tsnap);
    for (bp = btab; bp < btab+nbody; bp++) {
      SETV(lpos,Pos(bp));
      (*pot)(&ndim,lpos,lacc,&lphi,&tsnap);

      if (ome!=0.0) {
	lphi -= half_ome2*(sqr(lpos[0])+sqr(lpos[1]));
	lacc[0] += ome2*lpos[0] + two_ome*Vel(bp)[1];
	lacc[1] += ome2*lpos[1] - two_ome*Vel(bp)[0];
      }

      Phi(bp) = lphi;
      SETV(Acc(bp),lacc);
    }
    bits |= (PotentialBit|AccelerationBit|TimeBit);
    put_snap(outstr, &btab, &nbody, &tsnap, &bits);
  }
  strclose(instr);
  strclose(outstr);
}


