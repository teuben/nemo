/*
 * CODE_IO.C: I/O routines for potential + non-interacting particles N-body
 * Public routines: inputdata(), initoutput(), stopoutput(), output(),
 *		    savestate(), restorestate().
 *
 *    6-oct-92  repaired the long_changed convention out_mass -> mass 
 *              as was done in hackcode1 ages ago                       PJT
 *   10-apr-01  gcc warnings
 *    7-feb-07  gcc4 fix for prototype
 */

#include "defs.h"

#include <filestruct.h>
#include <history.h>
#include <snapshot/snapshot.h>
#include <snapshot/get_snap.c>

local void _put_snap_diagnostics(stream outstr, int *ofptr);

#define put_snap_diagnostics  _put_snap_diagnostics
#include <snapshot/put_snap.c>

local void diagnostics(void);

local bool Qkey, Qaux;

/*
 * INPUTDATA: read initial conditions from input file.
 */

void inputdata(void)
{
    stream instr;
    bodyptr btab, bp;
    int bits;

    instr = stropen(infile, "r");		/* open input stream        */
    get_history(instr);				/* read history */
    btab = &bodytab[0];				/* prepare input pointer    */
    nbody = MBODY;				/* set maximum size         */
    get_snap(instr, &btab, &nbody, &tnow, &bits);
    						/* invoke generic input     */
    if ((bits & PhaseSpaceBit) == 0)
	error("inputdata: essential data missing\tbits = %o", bits);

    if ((bits & MassBit) == 0) {                /* if no masses present */
	warning("No masses present: setting all to 1.0");
        for(bp=btab; bp<btab+nbody;bp++)
            Mass(bp)=1.0;                       /* set all masses to dummy 1*/
    }
    Qkey = bits & KeyBit;
    Qaux = bits & AuxBit;   /* will never be used, since we use this otherwise */

    if (scanopt(options, "reset_time") || (bits & TimeBit) == 0)
						/* no time specified?       */
    tnow = 0.0;				/*   then supply default    */
    strclose(instr);				/* close input stream       */
}

/*
 * INITOUTPUT: initialize output routines.
 */

local stream outstr = NULL;		/* output stream pointer */

void initoutput(void)
{
    printf("\n%s\n\n", headline);               /* headline log stream      */
    if (*outfile != 0) {                        /* output file given?       */
        outstr = stropen(outfile, "w");         /*   setup out. stream      */
	put_history(outstr);			/* output history */
        put_string(outstr, HeadlineTag, headline);
    }

    printf("%12s%12s%12s%12s%12s%12s%12s\n",
           "nbody", "freq", "eta", "sigma", "ome  ", "mode", "tstop");
    printf("%12d%12.4f%12.4f%12.4f%12.4f%12d%12.2f\n\n",
	   nbody, freq, eta, sigma, ome, mode, tstop);
    if (*options != 0)
	printf("\toptions: %s\n\n", options);
    minor_tout = tout = tnow;			/* schedule 1st outputs     */
}

/*
 * STOPOUTPUT: finish up after a run.
 */

void stopoutput(void)
{
  if (outstr) strclose(outstr);
}

/*
 * Diagnostics computed for output routines.
 */

local real mtot;                /* total mass of N-body system */
local real etot[3];             /* binding, kinetic, potential energy */
local matrix keten;		/* kinetic energy tensor */
local matrix peten;		/* potential energy tensor */
local matrix amten;		/* antisymmetric ang. mom. tensor */
local vector cmphase[2];	/* center of mass coordinates */

/*
 * OUTPUT: compute diagnostics and output binary data.
 */

local bool firstmass = TRUE;	/* set if masses have yet to be output */

void output(void)
{
    int k, bits;
    bodyptr btab;

    diagnostics();
    printf("\n  %12s%12s%12s%12s\n",
	   "tnow", "T+U", "T/U", "cputime");
    printf("  %12.3f%12.6f%12.4f%12.2f\n\n",
	   tnow, etot[0], etot[1]/etot[2], cputime());
    printf("\t    %10s", "cm pos");
    for (k = 0; k < NDIM; k++)
	printf("%10.4f", cmphase[0][k]);
    printf("\n\t    %10s", "cm vel");
    for (k = 0; k < NDIM; k++)
	printf("%10.4f", cmphase[1][k]);
    printf("\n");
    bits = 0;
    if (minor_freqout > 0.0 && (minor_tout - 0.01/freq) <= tnow) {
	minor_tout += 1.0 / minor_freqout;
	bits |= TimeBit;
    }
    if (freqout > 0.0 && (tout - 0.01/freq) <= tnow) {
	tout += 1.0 / freqout;
	bits |= TimeBit | PhaseSpaceBit;
	if (scanopt(options, "mass") || firstmass) {
	    bits |= MassBit;
	    if (Qkey) bits |= KeyBit;
	    firstmass = FALSE;
	}
	if (scanopt(options, "phi"))
	    bits |= PotentialBit;
	if (scanopt(options, "acc"))
	    bits |= AccelerationBit;
        if (scanopt(options, "aux"))
            bits |= AuxBit;
    }
    if (bits != 0 && outstr != NULL) {
	btab = &bodytab[0];
	put_snap(outstr, &btab, &nbody, &tnow, &bits);
	if (bits & PhaseSpaceBit)
	    printf("\n\tparticle data written\n");
    }
}

/*
 * DIAGNOSTICS: compute set of dynamical diagnostics.
 */

local void diagnostics(void)
{
    register bodyptr p;
    real velsq;
    vector tmpv;
    matrix tmpt;

    mtot = 0.0;					/* zero total mass */
    etot[1] = etot[2] = 0.0;			/* zero total KE and PE */
    CLRM(keten);				/* zero ke tensor */
    CLRM(peten);				/* zero pe tensor */
    CLRM(amten);				/* zero am tensor */
    CLRV(cmphase[0]);				/* zero c. of m. position */
    CLRV(cmphase[1]);				/* zero c. of m. velocity */
    for (p = bodytab; p < bodytab+nbody; p++) {	/* loop over all bodies */
	mtot += Mass(p);                        /*   sum body masses */
	DOTVP(velsq, Vel(p), Vel(p));		/*   square vel vector */
	etot[1] += 0.5 * Mass(p) * velsq;	/*   sum current KE */
	etot[2] +=       Mass(p) * Phi(p);	/*   and current PE (non-int.) */
	MULVS(tmpv, Vel(p), 0.5 * Mass(p));	/*   sum 0.5 m v_i v_j */
	OUTVP(tmpt, tmpv, Vel(p));
	ADDM(keten, keten, tmpt);
	MULVS(tmpv, Pos(p), Mass(p));		/*   sum m r_i a_j */
	OUTVP(tmpt, tmpv, Acc(p));
	ADDM(peten, peten, tmpt);
	OUTVP(tmpt, tmpv, Vel(p));		/*   sum m r_i v_j */
	ADDM(amten, amten, tmpt);
	MULVS(tmpv, Pos(p), Mass(p));		/*   sum cm position */
	ADDV(cmphase[0], cmphase[0], tmpv);
	MULVS(tmpv, Vel(p), Mass(p));		/*   sum cm momentum */
	ADDV(cmphase[1], cmphase[1], tmpv);
    }
    etot[0] = etot[1] + etot[2];                /* sum KE and PE */
    TRANM(tmpt, amten);				/* anti-sym. AM tensor */
    SUBM(amten, amten, tmpt);
    DIVVS(cmphase[0], cmphase[0], mtot);        /* normalize cm coords */
    DIVVS(cmphase[1], cmphase[1], mtot);
}

/*
 * PUT_SNAP_DIAGNOSTICS: output various N-body diagnostics.
 */

local void _put_snap_diagnostics(stream outstr, int *ofptr)
{
    real cput;

    cput = cputime();
    put_set(outstr, DiagnosticsTag);
    put_data(outstr, EnergyTag, RealType, etot, 3, 0);
    put_data(outstr, KETensorTag, RealType, keten, NDIM, NDIM, 0);
    put_data(outstr, PETensorTag, RealType, peten, NDIM, NDIM, 0);
    put_data(outstr, AMTensorTag, RealType, amten, NDIM, NDIM, 0);
    put_data(outstr, CMPhaseSpaceTag, RealType, cmphase, 2, NDIM, 0);
    put_data(outstr, "cputime", RealType, &cput, 0);
    put_tes(outstr, DiagnosticsTag);
}


/*
 * SAVESTATE: write current state to disk file.
 */

void savestate(string file)
{
    stream str;

    str = stropen(file, "a");			/* open state output file   */
    fseek(str, 0L, 0);				/* rewind stream to origin  */
    put_string(str, "program", getargv0());
    put_string(str, "version", getparam("VERSION"));
    put_string(str, "headline", headline);	/* save control parameters  */
    put_data(str, "freq", RealType, &freq, 0);
    put_data(str, "mode", IntType, &mode, 0);
    put_data(str, "eta", RealType, &eta, 0);
    put_data(str, "sigma", RealType, &sigma, 0);
    put_data(str, "dr", RealType, dr, NDIM, 0);
    put_string(str, "options", options);
    put_data(str, "tstop", RealType, &tstop, 0);
    put_data(str, "minor_freqout", RealType, &minor_freqout, 0);
    put_data(str, "freqout", RealType, &freqout, 0);
    put_data(str, "tnow", RealType, &tnow, 0);	/* save state variables     */
    put_data(str, "minor_tout", RealType, &minor_tout, 0);
    put_data(str, "tout", RealType, &tout, 0);
    put_data(str, "nbody", IntType, &nbody, 0);
    put_data(str, "bodytab", AnyType, bodytab, nbody, sizeof(body), 0);
    strclose(str);
}

/*
 * RESTORESTATE: restore state from disk file.
 */

void restorestate(string file)
{
    stream str;
    string program, version;

    str = stropen(file, "r");			/* open state input file    */
    program = get_string(str, "program");
    version = get_string(str, "version");
    if (! streq(program, getargv0()) ||		/* check program, version   */
	  ! streq(version, getparam("VERSION")))
	warning("state file may be outdated\n\n");
    headline = get_string(str, "headline");	/* read control parameters  */
    get_data(str, "freq", RealType, &freq, 0);
    get_data(str, "mode", IntType, &mode, 0);
    get_data(str, "eta", RealType, &eta, 0);
    get_data(str, "sigma", RealType, &sigma, 0);
    get_data(str, "dr", RealType, dr, NDIM, 0);
    options = get_string(str, "options");
    get_data(str, "tstop", RealType, &tstop, 0);
    get_data(str, "minor_freqout", RealType, &minor_freqout, 0);
    get_data(str, "freqout", RealType, &freqout, 0);
    get_data(str, "tnow", RealType, &tnow, 0);	/* read system state        */
    get_data(str, "minor_tout", RealType, &minor_tout, 0);
    get_data(str, "tout", RealType, &tout, 0);
    get_data(str, "nbody", IntType, &nbody, 0);
    get_data(str, "bodytab", AnyType, bodytab, nbody, sizeof(body), 0);
    strclose(str);
}
