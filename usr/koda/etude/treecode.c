/****************************************************************************/
/* TREECODE.C: new hierarchical N-body/SPH code.                            */
/* Copyright (c) 1999 by Joshua E. Barnes, Tokyo, JAPAN.                    */
/*    revised for SPH calculation by Jin Koda, Tokyo, JAPAN. 2000           */
/****************************************************************************/

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#define global                                  /* don't default to extern  */
#include "treecode.h"

/*
 * Default values for input parameters.
 */

#if defined(QUICKSCAN)
#  define SCANNER  "(quick scan)"
#else
#  define SCANNER  "(full scan)"
#endif

string defv[] = {               ";Hierarchical N-body/SPH code " SCANNER,
    "in=",                      ";Input file with initial conditions",
    "out=",                     ";Output file of N-body frames",
    "param=",                   ";Parameter file of calculation",
    "freq=32.0",                ";Fundamental integration frequency",
    "eps=0.025",                ";Density smoothing length",
    "selfgrav=false",           ";Compute self-gravity",
#if !defined(QUICKSCAN)
    "theta=1.0",                ";Force accuracy parameter",
#endif
    "usequad=false",            ";If true, use quad moments",
    "options=",                 ";Various control options",
    "tstop=3.0",                ";Time to stop integration",
    "freqout=30.0",             ";Data output quency",
    "nbody=5000",               ";Number of bodies for test run",
    "seed=123",                 ";Random number seed for test run",
    "save=",                    ";Write state file as code runs",
    "restore=",                 ";Continue run from state file",
    "mode=2.0",                 ";Mode of spiral pattern",  
    "pitch=10.0",               ";Pitch angle of spirals [deg]",
    "rinit=3.0",                ";Radius of initial disk [kpc]",
    "rcore=1.0",                ";Core radius [kpc]",
    "vmax=220.0",               ";Maximum disk rotation [km/s]",
    "omgb=26.1",                ";Bar pattern speed [km/s/kpc]",
    "fgas=0.05",                ";Mass fraction of gas in unit radius",   
    "fbar=0.1",                 ";Strength of bar potential",
    "rcorebh=0.01",             ";Core radius of central black hole [kpc]",
    "fbh=0.0",                  ";Mass ratio of central BH to the galaxy",
    "alpha=1.0",                ";Coefficient of viscousity",
    "beta=2.0",                 ";Coefficient of viscousity",
    "nnbr=32",                  ";Requested number of neighbors",
    "nmax=42",                  ";Maximum number of neighbors",
    "nmin=22",                  ";Minimum number of neighbors",
    "VERSION=beta1.5",          ";Jin Koda    Dec 12 2002",
    NULL,
};

/* Prototypes for local procedures. */

local void calcforce(void);                     /* do force calculation     */
local void stepsystem(void);                    /* advance by one time-step */
local real settimestep(void);                   /* set time step 'dt'       */
local void startrun(void);                      /* initialize system state  */

/*
 * MAIN: toplevel routine for hierarchical N-body code.
 */

void main(int argc, string argv[])
{
    initparam(argv, defv);                      /* initialize param access  */
    headline = defv[0] + 1;                     /* skip ";" in headline     */
    startrun();                                 /* get params & input data  */
    startoutput();                              /* activate output code     */
    if (nstep == 0) {                           /* if data just initialized */
#if defined(QLOOK)
	plot(bodytab, nbody);
#endif
        calcforce();                            /* do complete calculation  */
        output();                               /* and report diagnostics   */
    }

    if (freq != 0.0)                            /* if time steps requested  */
        while (tstop - tnow > 0.01/freq) {      /* while not past tstop     */
#if defined(QLOOK)
	    if (nstep%10 == 0) plot(bodytab, nbody);
#endif
            stepsystem();                       /* advance step by step     */
            output();                           /* and output results       */
        }
}

/*
 * CALCFORCE: common parts of force calculation.
 */

local void calcforce(void)
{
    bodyptr p;

    for (p = bodytab; p < bodytab+nbody; p++) { /* loop over all bodies     */
        Update(p) = TRUE;                       /* update all forces        */
	CLRV(Acc(p));                           /* clear acceleration       */ 
    }
    planttree(bodytab, nbody);                  /* construct tree structure */
    sphcalc(bodytab, nbody);                    /* compute SPH force        */
    if (selfgrav) {                             /* if compute self-gravity  */
	managetree();                           /* calc com and branck link */
	gravcalc();                             /* then compute them        */
    }
    extforce(bodytab, nbody, tnow);             /* compute external forces  */
    forcereport();                              /* print force statistics   */
}

/*
 * STEPSYSTEM: advance N-body system using simple leap-frog.
 */

local void stepsystem(void)
{
    real dt;

    bodyptr p;

    dt = settimestep();

    for (p = bodytab; p < bodytab+nbody; p++) { /* loop over all bodies     */
        ADDMULVS(Pos(p), Vel(p), 0.5 * dt);     /* advance r by 1/2 step    */
	ADDVMULVS(Velm(p), Vel(p), Acc(p), 0.5 * dt);
                                                /* predict v by 1/2 step    */
    }
    stephknl(bodytab, nbody, dt);               /* update smoothing length  */
    calcforce();                                /* perform force calc.      */
    for (p = bodytab; p < bodytab+nbody; p++) { /* loop over all bodies     */
        ADDMULVS(Vel(p), Acc(p), dt);           /* advance v by 1/2 step    */
        ADDMULVS(Pos(p), Vel(p), 0.5 * dt);     /* advance r by 1 step      */
    }
    nstep++;                                    /* count another time step  */
    tnow = tnow + dt;                           /* finally, advance time    */
}

/*
 * SETTIMESTEP: set time increment to advance.
 */

local real settimestep(void)
{
    bodyptr p;
    real dt, v, a;

    dtvel = 1.0e9;
    dtacc = 1.0e9;
    for (p = bodytab; p < bodytab+nbody; p++) { /* loop over all bodies     */
	ABSV(v, Vel(p));                        /* absolute of velocity     */
	ABSV(a, Acc(p));                        /* absolute of acceleration */
	if (v > 0.0)                            /* compute minimum time     */
	    dtvel = MIN(dtvel, Hknl(p) / v);    /* increment from velocity  */
	if (a > 0.0)                            /* criteria, and that       */
	    dtacc = MIN(dtacc, rsqrt(Hknl(p) / a)); /* from acceleration    */
    }
    dt = 0.3 * MIN(dtvel, dtacc);               /* find minimum of all the  */
    dt = MIN((0.4 * dtcfl), dt);                /* time criteria            */
    return(dt);
}

/*
 * STARTRUN: startup hierarchical N-body code.
 */

local void startrun(void)
{
    real temp;

    infile = getparam("in");                    /* set I/O file names       */
    outfile = getparam("out");
    savefile = getparam("save");
    if (strnull(getparam("restore"))) {         /* if starting a new run    */
        eps = getdparam("eps");                 /* get input parameters     */
        freq = getdparam("freq");
	selfgrav=getbparam("selfgrav");
#if !defined(QUICKSCAN)
        theta = getdparam("theta");
#endif
        usequad = getbparam("usequad");
        tstop = getdparam("tstop");
        freqout = getdparam("freqout");
        options = getparam("options");
	mode = getdparam("mode");
	pitch = getdparam("pitch");
	rinit = getdparam("rinit");
	vmax = getdparam("vmax");
	rcore = getdparam("rcore");
	omgb = getdparam("omgb");
	fgas = getdparam("fgas");
	fbar = getdparam("fbar");
	rcorebh = getdparam("rcorebh");
	fbh = getdparam("fbh");
	alpha = getdparam("alpha");
	beta = getdparam("beta");
	nnbr = getiparam("nnbr");
	nmax = getiparam("nmax");
	nmin = getiparam("nmin");

	rscale = 8.0;                           /* length unit              */
	mscale = massext(rscale);               /* mass within unit length  */
	tscale = 4.71477e11 * rsqrt(rpow(rscale,3.0) / mscale);
                                                /* time unit, from G=1      */
	massdk = massext(rinit) / mscale;       /* tot. disk mass in sys uni*/

	printf("\n\tSystem Unit\n");
	printf("\t\trscale = %.2f [kpc]\n", rscale);
	printf("\t\tmscale = %.2f [Msun]\n", mscale);
	printf("\t\ttscale = %.2f [year]\n", tscale);

        rsize = 1.0;                            /* start root w/ unit cube  */
        nstep = 0;                              /* begin counting steps     */
        tout = tnow;                            /* schedule first output    */
	initextf();                             /* init. ext. force params. */
        if (! strnull(infile))                  /* if data file was given   */
            inputdata();                        /* then read inital data    */
        else {                                  /* else make initial data   */
            tnow = 0.0;                         /* reset elapsed model time */
            nbody = getiparam("nbody");         /* get number of bodies     */
            if (nbody < 1)                      /* check for silly values   */
                error("startrun: absurd value for nbody\n");
            srandom0(getiparam("seed"));        /* set random number gen.   */
            testdata();                         /* and make plummer model   */
        }
    } else {                                    /* else restart old run     */
        restorestate(getparam("restore"));      /* read in state file       */
        if (getparamstat("eps") & ARGPARAM)     /* if given, set new params */
            eps = getdparam("eps");
        if (getparamstat("selfgrav") & ARGPARAM)
            usequad = getbparam("selfgrav");
#if !defined(QUICKSCAN)
        if (getparamstat("theta") & ARGPARAM)
            theta = getdparam("theta");
#endif
        if (getparamstat("usequad") & ARGPARAM)
            usequad = getbparam("usequad");
        if (getparamstat("options") & ARGPARAM)
            options = getparam("options");
        if (getparamstat("tstop") & ARGPARAM)
            tstop = getdparam("tstop");
        if (getparamstat("freqout") & ARGPARAM)
            freqout = getdparam("freqout");
        if (scanopt(options, "new-tout"))       /* if output time reset     */
            tout = tnow + 1 / freqout;          /* then offset from now     */
	temp = rsqrt(27.0/4.0) * rcore * vmax * vmax;
	initextf();                             /* init. ext. force params. */
    }
}
