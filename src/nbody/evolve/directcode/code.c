/*
 * CODE.C: simple direct N-body code.
 *	
 * updates:
 *     16-feb-04  V1.0   example, prove a point to Jim       PJT
 */

#include "code.h"

string defv[] = {		/* DEFAULT PARAMETER VALUES */

    /* file names for structured binary input/output */
    "in=\n			  snapshot of initial conditions ",
    "out=\n			  stream of output snapshots ",

    /* params used only if "in", "restart" and "continue" not given */
    "nbody=128\n		  number of particles to generate ",
    "seed=123\n			  random number generator seed ",
    "cencon=false\n		  centrally concentrated system ",

    /* params to control N-body integration */
    "freq=32.0\n		  fundamental integration frequency ",
    "eps=0.05\n			  usual potential softening ",
    "options=mass,phase\n	  misc. control options ",

    "tstop=2.0\n		  time to stop integration ",
    "freqout=4.0\n		  major data-output frequency ",
    "minor_freqout=32.0\n	  minor data-output frequency ",

    "VERSION=1.0\n		  16-feb-04 PJT",
    NULL,
};

string usage = "simple direct N-body code";

string headline = "DirectCode";

extern  bool scanopt(string, string);

void nemo_main(void)
{
  startrun();					/* set params, input data   */
  initoutput();				/* begin system output      */
  while (tnow < tstop + 0.1/freq)		/* while not past tstop     */
    stepsystem();				/*   advance N-body system  */
  stopoutput();				/* finish up output         */
}

/*
 * STARTRUN: startup hierarchical N-body code.
 */

void startrun(void)
{
  infile = getparam("in");			/* set I/O file names       */
  outfile = getparam("out");
  options = getparam("options");		/* set control options      */
  
  if (hasvalue("in"))
    inputdata(infile);			/*     read inital data     */
  else {					/*   make initial conds?    */
    nbody = getiparam("nbody");		/*     get nbody parameter  */
    if (nbody < 1)			/*     is value absurd?     */
      error("startrun: absurd nbody");
    init_xrandom(getparam("seed"));	/*     set random generator */
    testdata(getbparam("cencon"));	/*     make test model      */
  }
  freq = getdparam("freq");		/*   get various parameters */
  eps = getdparam("eps");
  tstop = getdparam("tstop");
  freqout = getdparam("freqout");
  minor_freqout = getdparam("minor_freqout");
  nstep = 0;				/*   start counting steps   */
  minor_tout = tout = tnow;		/*   schedule first output  */
}

/*
 * TESTDATA: generate initial conditions for test runs.
 * NOTE: Should really make a Plummer Model! 
 */

void testdata(bool cencon)
{
  vector cmr, cmv;
  bodyptr p;
  
  headline = "Direct code: test data";	/* supply default headline  */
  tnow = 0.0;					/* pos, vel set at t = 0    */
  bodytab = (bodyptr) allocate(nbody * sizeof(body));      /* MEMORY LEAK */
  CLRV(cmr);					/* init cm pos, vel         */
  CLRV(cmv);
  for (p = bodytab; p < bodytab+nbody; p++) {	/* loop over particles      */
    Mass(p) = 1.0 / nbody;			/*   set masses equal       */
    pickvec(Pos(p), cencon);		/*   pick position          */
    ADDV(cmr, cmr, Pos(p));
    pickvec(Vel(p), FALSE);			/*   pick velocity          */
    ADDV(cmv, cmv, Vel(p));
  }
  DIVVS(cmr, cmr, (real) nbody);		/* normalize cm coords      */
  DIVVS(cmv, cmv, (real) nbody);
  for (p = bodytab; p < bodytab+nbody; p++) {	/* loop over particles      */
    SUBV(Pos(p), Pos(p), cmr);		/*   offset by cm coords    */
    SUBV(Vel(p), Vel(p), cmv);
  }
}

/*
 * STEPSYSTEM: advance N-body system one time-step.
 */

void stepsystem(void)
{
  real dthf, dt;
  bodyptr p;
  vector acc1, dacc, dvel, vel1, dpos;

  dt = 1.0 / freq;				/* get basic time-step      */
  dthf = 0.5 * dt;				/* and basic half-step      */
  for (p = bodytab; p < bodytab+nbody; p++) { /* loop over particles      */
    SETV(acc1, Acc(p));			/*   save old acceleration  */
    hackgrav(p);				/*   compute new acc for p  */
    if (nstep > 0) {			/*   if past first step?    */
      SUBV(dacc, Acc(p), acc1);		/*     use change in accel  */
      MULVS(dvel, dacc, dthf);		/*     to make 2nd order    */
      ADDV(Vel(p), Vel(p), dvel);		/*     correction to vel    */
    }
  }
  output();					/* do major or minor output */
  for (p = bodytab; p < bodytab+nbody; p++) {	/* loop advancing bodies    */
    MULVS(dvel, Acc(p), dthf);		/*   use current accel'n    */
    ADDV(vel1, Vel(p), dvel);		/*   find vel at midpoint   */
    MULVS(dpos, vel1, dt);			/*   find pos at endpoint   */
    ADDV(Pos(p), Pos(p), dpos);		/*   advance position       */
    ADDV(Vel(p), vel1, dvel);		/*   advance velocity       */
  }
  nstep++;					/* count another mu-step    */
  tnow = tnow + dt;				/* finally, advance time    */
}
