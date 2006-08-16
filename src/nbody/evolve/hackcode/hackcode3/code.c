/*
 * CODE.C: hierarchical N-body code with various kludges...
 *	
 * updates:
 *     xx-xxx-xx  V1.0  cloned off hackcode1		JEB/PJT
 * 	7-jul-89  V1.1  changed names of 'options'	PJT
 *     26-jan-90  V1.2  added optional rigid particles, PJT
 *				and external potential
 *	3-apr-90  V1.2a bug removed			PJT
 *     28-jan-94  V1.2b srandom -> set_xrandom for Solaris  PJT
 *                      also use nemo_main now
 *     24-mar-94  V1.2c ansi + using allocate()
 *	8-sep-01	init_xrandom
 *     11-jan-02      d debug= renamed to hdebug=
 *     15-aug-06  V1.3  prototype and global/extern as in hackcode1	PJT
 */

#define global
#include "code.h"

string defv[] = {		/* DEFAULT PARAMETER VALUES */

    /* file names for structured binary input/output */
    "in=\n			Snapshot of initial conditions",
    "out=\n			Stream of output snapshots",

    /* file names for saving/restoring program state */
    "restart=\n			Input state and set controls",
    "continue=\n		Input state and continue run",
    "save=\n			Output state as code runs",

    /* params used only if "in", "restart" and "continue" not given */
    "nbody=128\n		Number of particles to generate",
    "seed=123\n			Random number generator seed",
    "cencon=false\n		Centrally concentrated system",

    /* params to control N-body integration */
    "freq=32.0\n		Fundamental integration frequency",
    "eps=0.05\n			Usual potential softening",
    "tol=1.0\n			Cell subdivision tolerence",
    "fcells=1.0\n		Cell allocation parameter",
    "options=mass,phase\n	Misc. control options",

    "nrigid=0\n                 Number of rigid particles",

    "tstop=2.0\n		Time to stop integration",
    "freqout=4.0\n		Major data-output frequency",
    "minor_freqout=32.0\n	Minor data-output frequency",

    "potname=\n			Name for potential(5)",
    "potpars=\n			Parameters for potential(5)",
    "potfile=\n			Filename for potential(5)",

    "hdebug=false\n		Turn on debugging messages",
    "VERSION=1.3a\n		15-aug-06 PJT",
    NULL,
};

string headline = "Hack code3";	/* default id for run */

string usage = "hierarchical N-body code, with potential(5NEMO) descriptors";

string cvsid="$Id$";



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

startrun()
{
    string restfile, contfile;
    bool scanopt();
    proc get_potential();
    infile = getparam("in");			/* set I/O file names       */
    outfile = getparam("out");
    restfile = getparam("restart");
    contfile = getparam("continue");
    savefile = getparam("save");
    options = getparam("options");		/* set control options      */
    debug = getbparam("hdebug");
    if (debug)
	dprintf(0,"hdebug is turned on");
    nrigid = getiparam("nrigid");
    if (*contfile)	         		/* resume interrupted run   */
	restorestate(contfile);
    else if (*restfile) {	        	 /* resume w/ new parameters */
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
	if (*infile)	                 	/*   was data file given?   */
	    inputdata(infile);			/*     read inital data     */
	else {					/*   make initial conds?    */
	    nbody = getiparam("nbody");		/*     get nbody parameter  */
	    if (nbody < 1)			/*     is value absurd?     */
		error("startrun: absurd nbody\n");
	    init_xrandom(getparam("seed"));	/*     set random generator */
	    testdata(getbparam("cencon"));	/*     make test model      */
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
    contfile = getparam("potname");
    if (*contfile) {
        extpot = get_potential(contfile, 
                    getparam("potpars"),getparam("potfile"));
    }
}

/*
 * TESTDATA: generate initial conditions for test runs.
 * NOTE: Should really make a Plummer Model! 
 */

testdata(cencon)
bool cencon;				/* make concentrated system */
{
    vector cmr, cmv;
    register bodyptr p;

    headline = "Hack code: test data";		/* supply default headline  */
    tnow = 0.0;					/* pos, vel set at t = 0    */
    bodytab = (bodyptr) allocate(nbody * sizeof(body));
    if (bodytab == NULL)
	error("testdata: not enuf memory\n");
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

stepsystem()
{
    real dthf, dt;
    register bodyptr p;
    vector acc1, dacc, dvel, vel1, dpos;
    int i;

    dt = 1.0 / freq;				/* get basic time-step      */
    dthf = 0.5 * dt;				/* and basic half-step      */
    maketree(bodytab, nbody);			/* load bodies into tree    */
    nfcalc = n2bcalc = nbccalc = 0;		/* zero interaction counts  */
    for (p = bodytab,i=0;i<nbody; p++,i++) {    /* loop over particles      */
      SETV(acc1, Acc(p));                       /*   save old acceleration  */
      hackgrav(p);                              /*   compute new acc for p  */
      nfcalc++;                                 /*   count force calcs      */
      n2bcalc += n2bterm;                       /*   and 2-body terms       */
      nbccalc += nbcterm;                       /*   and body-cell terms    */
      if (i>=nrigid) {				/* move them if not fixed   */
	if (nstep > 0) {			/*   if past first step?    */
	    SUBV(dacc, Acc(p), acc1);		/*     use change in accel  */
	    MULVS(dvel, dacc, dthf);		/*     to make 2nd order    */
	    ADDV(Vel(p), Vel(p), dvel);		/*     correction to vel    */
	}
      }
    }
    output();					/* do major or minor output */
    for (p = bodytab, i=0; i<nbody; p++,i++) {  /* loop advancing bodies    */
        if (i>=nrigid) {                        /* only for non-rigid part's*/
	    MULVS(dvel, Acc(p), dthf);		/*   use current accel'n    */
	    ADDV(vel1, Vel(p), dvel);		/*   find vel at midpoint   */
	    MULVS(dpos, vel1, dt);		/*   find pos at endpoint   */
	    ADDV(Pos(p), Pos(p), dpos);		/*   advance position       */
	    ADDV(Vel(p), vel1, dvel);		/*   advance velocity       */
        }
    }
    nstep++;					/* count another mu-step    */
    tnow = tnow + dt;				/* finally, advance time    */
}
