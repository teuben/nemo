/****************************************************************************/
/* TREECODE.C: new hierarchical N-body code.                                */
/* Copyright (c) 2001 by Joshua E. Barnes, Honolulu, Hawai`i.               */
/* 22-jun-01   NEMOfied for NEMO V3                                         */
/* 22-feb-04   dtime->dtimes                                                */
/* 25-apr-04   implemented USE_NEMO_IO                                      */ 
/****************************************************************************/

#include <stdinc.h>
#include <mathfns.h>
#include <vectmath.h>	
#include <getparam.h>
#define global                                  /* don't default to extern  */
#include "treecode.h"

/*
 * Default values for input parameters.
 */

string defv[] = {
#if defined(USE_NEMO_IO)
    "in=\n                       Input file with initial conditions (snapshot format)",
    "out=\n                      Output file of N-body frames (snapshot format)",
#else
    "in=\n                       Input file with initial conditions",
    "out=\n                      Output file of N-body frames",
#endif
#if defined(USEFREQ)
    "freq=32.0\n                 Leapfrog integration frequency",
#else
    "dtime=1/32\n                Leapfrog integration timestep",
#endif
    "eps=0.025\n                 Density smoothing length",
#if !defined(QUICKSCAN)
    "theta=1.0\n                 Force accuracy parameter",
#endif
    "usequad=false\n             If true, use quad moments",
    "options=\n                  Various control options",
    "tstop=2.0\n                 Time to stop integration",
#if defined(USEFREQ)
    "freqout=4.0\n               Data output frequency",
#else
    "dtout=1/4\n                 Data output timestep",
#endif
    "nbody=4096\n                Number of bodies for test run",
    "seed=123\n                  Random number seed for test run",
    "save=\n                     Write state file as code runs",
    "restore=\n                  Continue run from state file",
    "VERSION=1.4.2\n             25-apr-04 PJT",
    NULL,
};

#if defined(QUICKSCAN)
string usage = "Hierarchical N-body code (quick scan)";
#else
string usage = "Hierarchical N-body code (theta scan)";
#endif


/* Prototypes for local procedures. */

local void treeforce(void);                     /* do force calculation     */
local void stepsystem(void);                    /* advance by one time-step */
local void startrun(void);                      /* initialize system state  */
local void testdata(void);                      /* generate test data       */

/*
 * MAIN: toplevel routine for hierarchical N-body code.
 */

int nemo_main()
{
    headline = usage;                           /* headline                 */
    startrun();                                 /* get params & input data  */
    startoutput();                              /* activate output code     */
    if (nstep == 0) {                           /* if data just initialized */
        treeforce();                            /* do complete calculation  */
        output();                               /* and report diagnostics   */
    }
#if defined(USEFREQ)
    if (freq != 0.0)                            /* if time steps requested  */
        while (tstop - tnow > 0.01/freq) {      /* while not past tstop     */
            stepsystem();                       /* advance step by step     */
            output();                           /* and output results       */
        }
#else
    if (dtimes != 0.0)                          /* if time steps requested  */
        while (tstop - tnow > 0.01 * dtimes) {  /* while not past tstop     */
            stepsystem();                       /* advance step by step     */
            output();                           /* and output results       */
        }
#endif
    return 0;                                   /* end with proper status   */
}

/*
 * TREEFORCE: common parts of force calculation.
 */

local void treeforce(void)
{
    bodyptr p;

    for (p = bodytab; p < bodytab+nbody; p++)   /* loop over all bodies     */
        Update(p) = TRUE;                       /* update all forces        */
    maketree(bodytab, nbody);                   /* construct tree structure */
    gravcalc();                                 /* compute initial forces   */
    forcereport();                              /* print force statistics   */
}

/*
 * STEPSYSTEM: advance N-body system using simple leap-frog.
 */

local void stepsystem(void)
{
#if defined(USEFREQ)
    real dtimes = 1.0 / freq;                   /* set basic time-step      */
#endif
    bodyptr p;

    for (p = bodytab; p < bodytab+nbody; p++) { /* loop over all bodies     */
        ADDMULVS(Vel(p), Acc(p), 0.5 * dtimes); /* advance v by 1/2 step    */
        ADDMULVS(Pos(p), Vel(p), dtimes);       /* advance r by 1 step      */
    }
    treeforce();                                /* perform force calc.      */
    for (p = bodytab; p < bodytab+nbody; p++) { /* loop over all bodies     */
        ADDMULVS(Vel(p), Acc(p), 0.5 * dtimes); /* advance v by 1/2 step    */
    }
    nstep++;                                    /* count another time step  */
    tnow = tnow + dtimes;                       /* finally, advance time    */
}

/*
 * STARTRUN: startup hierarchical N-body code.
 */

local void startrun(void)
{
    infile = getparam("in");                    /* set I/O file names       */
    outfile = getparam("out");
    savefile = getparam("save");
    if (strnull(getparam("restore"))) {         /* if starting a new run    */
        eps = getdparam("eps");                 /* get input parameters     */
#if defined(USEFREQ)
        freq = getdparam("freq");
#else
        dtimes = getdparam("dtime");
#endif
#if !defined(QUICKSCAN)
        theta = getdparam("theta");
#endif
        usequad = getbparam("usequad");
        tstop = getdparam("tstop");
#if defined(USEFREQ)
        freqout = getdparam("freqout");
#else
        dtout = getdparam("dtout");
#endif
        options = getparam("options");
        if (! strnull(infile))                  /* if data file was given   */
            inputdata();                        /* then read inital data    */
        else {                                  /* else make initial data   */
            nbody = getiparam("nbody");         /* get number of bodies     */
            if (nbody < 1)                      /* check for silly values   */
                error("startrun: %d is an absurd value for nbody",nbody);
            srandom(getiparam("seed"));         /* set random number gen.   */
            testdata();                         /* and make plummer model   */
            tnow = 0.0;                         /* reset elapsed model time */
        }
        rsize = 1.0;                            /* start root w/ unit cube  */
        nstep = 0;                              /* begin counting steps     */
        tout = tnow;                            /* schedule first output    */
    } else {                                    /* else restart old run     */
#if 1
      error("restore= not yet implemented for NEMO version");
#else
        restorestate(getparam("restore"));      /* read in state file       */
        if (getparamstat("eps") & ARGPARAM)     /* if given, set new params */
            eps = getdparam("eps");
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
#if defined(USEFREQ)
        if (getparamstat("freqout") & ARGPARAM)
            freqout = getdparam("freqout");
        if (scanopt(options, "new-tout"))       /* if output time reset     */
            tout = tnow + 1 / freqout;          /* then offset from now     */
#else
	dtout = getdparam("dtout");
        if (scanopt(options, "new-tout"))       /* if output time reset     */
            tout = tnow + dtout;                /* then offset from now     */
#endif
#endif
    }
}

/*
 * TESTDATA: generate Plummer model initial conditions for test runs,
 * scaled to units such that M = -4E = G = 1 (Henon, Heggie, etc).
 * See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37, 183.
 */

#define MFRAC  0.999                            /* cut off 1-MFRAC of mass  */

local void testdata(void)
{
    real rsc, vsc, r, v, x, y;
    vector rcm, vcm;
    bodyptr p;

    bodytab = (bodyptr) allocate(nbody * sizeof(body));
                                                /* alloc space for bodies   */
    rsc = (3 * PI) / 16;                        /* and length scale factor  */
    vsc = rsqrt(1.0 / rsc);                     /* find speed scale factor  */
    CLRV(rcm);                                  /* zero out cm position     */
    CLRV(vcm);                                  /* zero out cm velocity     */
    for (p = bodytab; p < bodytab+nbody; p++) { /* loop over bodies         */
        Type(p) = BODY;                         /* tag as a body            */
        Mass(p) = 1.0 / nbody;                  /* set masses equal         */
        x = xrandom(0.0, MFRAC);                /* pick enclosed mass       */
        r = 1.0 / rsqrt(rpow(x, -2.0/3.0) - 1); /* find enclosing radius    */
        pickshell(Pos(p), NDIM, rsc * r);       /* pick position vector     */
        do {                                    /* select from fn g(x)      */
            x = xrandom(0.0, 1.0);              /* for x in range 0:1       */
            y = xrandom(0.0, 0.1);              /* max of g(x) is 0.092     */
        } while (y > x*x * rpow(1 - x*x, 3.5)); /* using von Neumann tech   */
        v = x * rsqrt(2.0 / rsqrt(1 + r*r));    /* find resulting speed     */
        pickshell(Vel(p), NDIM, vsc * v);       /* pick velocity vector     */
        ADDMULVS(rcm, Pos(p), 1.0 / nbody);     /* accumulate cm position   */
        ADDMULVS(vcm, Vel(p), 1.0 / nbody);     /* accumulate cm velocity   */
    }
    for (p = bodytab; p < bodytab+nbody; p++) { /* loop over bodies again   */
        SUBV(Pos(p), Pos(p), rcm);              /* subtract cm position     */
        SUBV(Vel(p), Vel(p), vcm);              /* subtract cm velocity     */
    }
}
