/*
 * CODE.C: hierarchical N-body code.
 *	
 * updates:
 * 	7-jul-89  V1.1   changed names of 'options'         (PJT)
 *     26-may-91  V1.2   fixed IMAX in defs.h for Cray      (PJT)
 *     18-nov-91  V1.2a  malloc -> allocate                 (PJT)
 *     15-aug-92  V1.2b  nemo_main such to return 0 to shell(PJT)
 *     18-jan-94  V1.3   srandom -> set_xrandom()	     pjt 
 *     28-nov-00  documented a memory leak
 *      8-sep-01  init_xrandom
 *     29-mar-04  MacOS forcing us to use global/extern      pjt
 *                plus LOTS of prototype cleanup
 */

#define global                                  /* don't default to extern  */
#include "code.h"

string defv[] = {		/* DEFAULT PARAMETER VALUES */

    /* file names for structured binary input/output */
    "in=\n			  snapshot of initial conditions ",
    "out=\n			  stream of output snapshots ",

    /* file names for saving/restoring program state */
    "restart=\n			  input state and set controls ",
    "continue=\n		  input state and continue run ",
    "save=\n			  output state as code runs ",

    /* params used only if "in", "restart" and "continue" not given */
    "nbody=128\n		  number of particles to generate ",
    "seed=123\n			  random number generator seed ",
    "cencon=false\n		  centrally concentrated system ",

    /* params to control N-body integration */
    "freq=32.0\n		  fundamental integration frequency ",
    "eps=0.05\n			  usual potential softening ",
    "tol=1.0\n			  cell subdivision tolerence ",
    "fcells=1.0\n		  cell allocation parameter ",
    "options=mass,phase\n	  misc. control options ",

    "tstop=2.0\n		  time to stop integration ",
    "freqout=4.0\n		  major data-output frequency ",
    "minor_freqout=32.0\n	  minor data-output frequency ",

    "debug=false\n		  turn on debugging messages ",
    "VERSION=1.4\n		  29-mar-04 PJT",
    NULL,
};

string usage = "hierarchical N-body code";

string headline = "Hack code";	/* default id for run */

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
    string restfile, contfile;
    bool scanopt();

    infile = getparam("in");			/* set I/O file names       */
    outfile = getparam("out");
    restfile = getparam("restart");
    contfile = getparam("continue");
    savefile = getparam("save");
    options = getparam("options");		/* set control options      */
    debug = getbparam("debug");
    if (*contfile)				/* resume interrupted run   */
	restorestate(contfile);
    else if (*restfile) {			/* resume w/ new parameters */
	restorestate(restfile);
	/* NOTE: someday, I will have a way to tell which, if any, of these *
	 * parameters are actually input from the command line, and only    *
	 * change them.  ANY NON-DEFAULT ARGS MUST BE SPECIFIED AT RESTART. */
	eps = getdparam("eps");			/*   get modified params    */
	tol = getdparam("tol");
	options = getparam("options");		/*   restorestate overwrite */
	fcells = getdparam("fcells");
	tstop = getdparam("tstop");
	freqout = getdparam("freqout");
	minor_freqout = getdparam("minor_freqout");
	if (scanopt(options, "new_tout")) {	/*   reset output times?    */
	    tout = tnow + 1 / freqout;		/*     offset from present  */
	    minor_tout = tnow + 1 / minor_freqout;
	}
    } else {					/* start new calculation    */
	if (*infile)				/*   was data file given?   */
	    inputdata(infile);			/*     read inital data     */
	else {					/*   make initial conds?    */
	    nbody = getiparam("nbody");		/*     get nbody parameter  */
	    if (nbody < 1)			/*     is value absurd?     */
		error("startrun: absurd nbody\n");
	    init_xrandom(getparam("seed"));	/*     set random generator */
	    make_testdata(getbparam("cencon"));	/*     make test model      */
	}
	freq = getdparam("freq");		/*   get various parameters */
	eps = getdparam("eps");
	tol = getdparam("tol");
	fcells = getdparam("fcells");
	tstop = getdparam("tstop");
	freqout = getdparam("freqout");
	minor_freqout = getdparam("minor_freqout");
	nstep = 0;				/*   start counting steps   */
	minor_tout = tout = tnow;		/*   schedule first output  */
	SETVS(rmin, -2.0);			/*   init box scaling       */
	rsize = -2.0 * rmin[0];
    }
}

/*
 * MAKE_TESTDATA: generate initial conditions for test runs.
 * NOTE: Should really make a Plummer Model! 
 */

void make_testdata(bool cencon)			/* make concentrated system */
{
    vector cmr, cmv;
    register bodyptr p;

    headline = "Hack code: test data";		/* supply default headline  */
    tnow = 0.0;					/* pos, vel set at t = 0    */
    bodytab = (bodyptr) allocate(nbody * sizeof(body));      /* MEMORY LEAK */
    CLRV(cmr);					/* init cm pos, vel         */
    CLRV(cmv);
    for (p = bodytab; p < bodytab+nbody; p++) {	/* loop over particles      */
	Type(p) = BODY;				/*   tag as a body          */
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
    register bodyptr p;
    vector acc1, dacc, dvel, vel1, dpos;

    dt = 1.0 / freq;				/* get basic time-step      */
    dthf = 0.5 * dt;				/* and basic half-step      */
    maketree(bodytab, nbody);			/* load bodies into tree    */
    nfcalc = n2bcalc = nbccalc = 0;		/* zero interaction counts  */
    for (p = bodytab; p < bodytab+nbody; p++) { /* loop over particles      */
	SETV(acc1, Acc(p));			/*   save old acceleration  */
	hackgrav(p);				/*   compute new acc for p  */
	nfcalc++;				/*   count force calcs      */
	n2bcalc += n2bterm;			/*   and 2-body terms       */
	nbccalc += nbcterm;			/*   and body-cell terms    */
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
