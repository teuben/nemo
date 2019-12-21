/*
 * CODE.C: simple direct N-body code (directcode)
 *	
 * updates:
 *     16-feb-04  V1.0   example, prove a point to Jim       PJT
 *     22-feb-04  V1.1   leapfrog corrected                  PJT
 *     25-feb-04   1.1a  corrected the last correction       PJT
 *     24-dec-04   1.1b  use global for MacOSX               PJT
 *     21-jul-09   1.1c  added code to check euler steps at PiTP09
 *     29-jul-09   1.2   added option eps < 0 for PN force   PJT
 *     30-jul-09   1.3   added option gravc=                 PJT
 */

#define global
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
    "eps=0.05\n			  if > 0 usual potential softening, if < 0, pseudo-newtonion ",
    "options=mass,phase\n	  misc. control options ",

    "tstop=2.0\n		  time to stop integration ",
    "freqout=4.0\n		  major data-output frequency ",
    "minor_freqout=32.0\n	  minor data-output frequency ",

    /* constants */

    "gravc=1\n                    Gravitatonal constant",

    "VERSION=1.3b\n		  6-dec-2016 PJT",
    NULL,
};

string usage = "simple direct N-body code";

string headline = "DirectCode";

string cvsid="$Id$";

extern  bool scanopt(string, string);


/* apple lvvm could not deal with this?
  #define stepsystem stepsystem_leapfrog
 */

void nemo_main(void)
{
  startrun();				/* set params, input data   */
  initoutput();				/* begin system output      */
  while (tnow < tstop + 0.1/freq)	/* while not past tstop     */
    stepsystem_leapfrog();		/*   advance N-body system  */
  stopoutput();				/* finish up output         */
}

/*
 * STARTRUN: startup hierarchical N-body code.
 */

void startrun(void)
{
  int seed = init_xrandom(getparam("seed")); /*  set random generator */

  infile = getparam("in");		/* set I/O file names       */
  outfile = getparam("out");
  options = getparam("options");	/* set control options      */
  
  if (hasvalue("in"))
    inputdata(infile);			/*     read inital data     */
  else {				/*   make initial conds?    */
    nbody = getiparam("nbody");		/*     get nbody parameter  */
    if (nbody < 1)			/*     is value absurd?     */
      error("invalid nbody=%d",nbody);
    testdata(getbparam("cencon"));	/*     make test model      */
  }
  freq = getdparam("freq");		/*   get various parameters */
  eps = getdparam("eps");               /*   softening length       */
  gravc = getdparam("gravc");           /*   grav constant          */  
  tstop = getdparam("tstop");           /*   stop time              */
  freqout = getdparam("freqout");       /*   output frequency       */
  minor_freqout = getdparam("minor_freqout");
  nstep = 0;				/*   start counting steps   */
  minor_tout = tout = tnow;		/*   schedule first output  */
}

/*
 * TESTDATA: generate initial conditions for test runs.
 *           this is exactly the same code as in hackcode1
 */

void testdata(bool cencon)
{
  vector cmr, cmv;
  bodyptr p;
  
  headline = "Direct code: test data";	/* supply default headline  */
  tnow = 0.0;					/* pos, vel set at t = 0    */
  bodytab = (bodyptr) allocate(nbody * sizeof(body));      /* MEMORY LEAK   */
  CLRV(cmr);					/* init cm pos, vel         */
  CLRV(cmv);
  for (p = bodytab; p < bodytab+nbody; p++) {	/* loop over particles      */
    Mass(p) = 1.0 / nbody;			/*   set masses equal       */
    pickvec(Pos(p), cencon);		        /*   pick position          */
    ADDV(cmr, cmr, Pos(p));
    pickvec(Vel(p), FALSE);			/*   pick velocity          */
    ADDV(cmv, cmv, Vel(p));
  }
  DIVVS(cmr, cmr, (real) nbody);		/* normalize cm coords      */
  DIVVS(cmv, cmv, (real) nbody);
  for (p = bodytab; p < bodytab+nbody; p++) {	/* loop over particles      */
    SUBV(Pos(p), Pos(p), cmr);	         	/*   offset by cm coords    */
    SUBV(Vel(p), Vel(p), cmv);
  }
}

/*
 * STEPSYSTEM: advance N-body system one time-step using a leapfrog stepper
 */

void stepsystem_leapfrog(void)
{
  real dthf, dt;
  bodyptr p;

  dt = 1.0 / freq;				/* get basic time-step      */
  dthf = 0.5 * dt;				/* and basic half-step      */
  if (nstep==0) {
    for (p = bodytab; p < bodytab+nbody; p++) {
      hackgrav(p);				/*   compute new acc for p  */
    }
  }
  output();					/* do major or minor output */
  for (p = bodytab; p < bodytab+nbody; p++) {	/* loop advancing bodies    */
    ADDMULVS(Vel(p), Acc(p), dthf);             /* advance v by 1/2 step    */
    ADDMULVS(Pos(p), Vel(p), dt);               /* advance r by 1 step      */
  }
  for (p = bodytab; p < bodytab+nbody; p++) {	/* loop getting new forces  */
    hackgrav(p);				/*   compute new acc for p  */
  }
  for (p = bodytab; p < bodytab+nbody; p++) {   /* loop over all bodies     */
    ADDMULVS(Vel(p), Acc(p), dthf);             /* advance v by 1/2 step    */
  }
  nstep++;					/* count another mu-step    */
  tnow = tnow + dt;				/* finally, advance time    */
}

/*
 * STEPSYSTEM: advance N-body system one time-step using a forward euler stepper
 */

void stepsystem_euler(void)
{
  real dthf, dt;
  bodyptr p;

  dt = 1.0 / freq;				/* get basic time-step      */

  if (nstep==0) {
    for (p = bodytab; p < bodytab+nbody; p++) {
      hackgrav(p);				/*   compute new acc for p  */
    }
  }
  output();					/* do major or minor output */
  for (p = bodytab; p < bodytab+nbody; p++) {	/* loop advancing bodies    */
    ADDMULVS(Pos(p), Vel(p), dt);               /* advance r by 1 step      */
    ADDMULVS(Vel(p), Acc(p), dt);               /* advance v by 1 step      */
  }
  for (p = bodytab; p < bodytab+nbody; p++) {	/* loop getting new forces  */
    hackgrav(p);				/*   compute new acc for p  */
  }

  nstep++;					/* count another mu-step    */
  tnow = tnow + dt;				/* finally, advance time    */
}

/*
 * here is the very old leapfrogger we used in 1986 .... see also hackcode1
 */

void stepsystem_old(void)
{
  real dthf, dt;
  bodyptr p;
  vector acc1, dacc, dvel, vel1, dpos;

  dt = 1.0 / freq;				/* get basic time-step      */
  dthf = 0.5 * dt;				/* and basic half-step      */
  for (p = bodytab; p < bodytab+nbody; p++) {   /* loop over particles      */
    SETV(acc1, Acc(p));			        /*   save old acceleration  */
    hackgrav(p);				/*   compute new acc for p  */
    if (nstep > 0) {			        /*   if past first step?    */
      SUBV(dacc, Acc(p), acc1);		        /*     use change in accel  */
      MULVS(dvel, dacc, dthf);		        /*     to make 2nd order    */
      ADDV(Vel(p), Vel(p), dvel);		/*     correction to vel    */
    }
  }
  output();					/* do major or minor output */
  for (p = bodytab; p < bodytab+nbody; p++) {	/* loop advancing bodies    */
    MULVS(dvel, Acc(p), dthf);		        /*   use current accel'n    */
    ADDV(vel1, Vel(p), dvel);		        /*   find vel at midpoint   */
    MULVS(dpos, vel1, dt);			/*   find pos at endpoint   */
    ADDV(Pos(p), Pos(p), dpos);		        /*   advance position       */
    ADDV(Vel(p), vel1, dvel);		        /*   advance velocity       */
  }
  nstep++;					/* count another mu-step    */
  tnow = tnow + dt;				/* finally, advance time    */
}
