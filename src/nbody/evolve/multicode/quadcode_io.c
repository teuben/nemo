/*
 * QUADCODE_IO.C: I/O routines for quadrupole N-body code.
 * Public routines: inputdata(), initoutput(), stopoutput(), output(),
 *		    savestate(), restorestate().
 *
 *	xx-xxx-xx	JEB	original version
 *	16-mar-90	PJT	made GCC happy
 *	12-nov-91	PJT	made options compatible with hackcode1
 *	21-may-92	PJT	sgi needs some forward decl.
 *      20-may-94       pjt     allocate() decl. into header
 *      15-aug-06       pjt     prototype fixes
 *
 *    'BUG': headline is taken from input file if it was present
 */

#include "quaddefs.h"

#include <filestruct.h>
#include <history.h>

local diagnostics(), put_quadfield();

#include <snapshot/snapshot.h>
#include <snapshot/get_snap.c>
#define put_snap_diagnostics  my_put_snap_diagnostics
local void put_snap_diagnostics(stream, int*);

#include <snapshot/put_snap.c>

/*
 * INPUTDATA: read initial conditions from input file.
 */

void inputdata()
{
    stream instr;
    int bits;

    instr = stropen(infile, "r");		/* open input stream        */
    get_history(instr);
#if 0
    if (ask_headline() != NULL && streq(headline, ""))
	headline = ask_headline();
#endif
    bodytab = NULL;				/* prepare input pointer    */
    nbody = 0;
    get_snap(instr, &bodytab, &nbody, &tnow, &bits);
    						/* invoke generic input     */
    strclose(instr);				/* close input stream       */
    if ((bits & MassBit) == 0 || (bits & PhaseSpaceBit) == 0)
	error("inputdata: essential data missing\tbits = 0x%x", bits);
    if (nbody > MBODY)
	error("inputdata: nbody(%d) > MBODY(%d) (recompile)", nbody,MBODY);
    if (scanopt(options, "reset_time") || (bits & TimeBit) == 0)
						/* no time specified?       */
	tnow = 0.0;				/*   then supply default    */
}

/*
 * INITOUTPUT: initialize output routines.
 */

local stream outstr = NULL;		/* output stream pointer */

local stream quadstr = NULL;		/* quadrupole field output */

void initoutput()
{
    if (! streq(headline, "")) {		/* nontrivial headline?     */
	printf("\n%s\n", headline);		/* headline log stream      */
	set_headline(headline);			/* and save with history    */
    }
    if (*outfile) {                             /* output file given?       */
        outstr = stropen(outfile, "w");         /*   setup out. stream      */
	put_history(outstr);
    }
    if (*quadfile) {                            /* field file given?        */
        quadstr = stropen(quadfile, "w");       /*   setup out. stream      */
        put_history(quadstr);
    }
    printf("\n%12s%12s%12s%12s%12s%12s\n",
           "nbody", "freq", "eps_r", "eps_t", "mode", "tstop");
    printf("%12d%12.4f%12.4f%12.4f%12d%12.2f\n\n",
	   nbody, freq, eps1, eps2, mode, tstop);
    if (*options)
	printf("\toptions: %s\n\n", options);
    minor_tout = tout = tnow;			/* schedule 1st outputs     */
    if (*savefile)       			/* restart file enabled?    */
	savestate(savefile);			/*   save inital state      */
}

/*
 * STOPOUTPUT: finish up after a run.
 */

void stopoutput()
{
    if (outstr != NULL)
        strclose(outstr);
    if (quadstr != NULL)
        strclose(quadstr);
}

/*
 * OUTPUT: compute diagnostics and output binary data.
 */

local real mtot;                /* total mass of N-body system              */
local real etot[3];             /* binding, kinetic, potential energy       */
local matrix keten;		/* kinetic energy tensor                    */
local matrix peten;		/* potential energy tensor                  */
local matrix amten;		/* antisymmetric ang. mom. tensor           */
local vector cmphase[2];	/* center of mass coordinates               */

local bool firstmass = TRUE;	/* set if masses have yet to be output      */

void output()
{
    int k, bits;

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
	if (firstmass || scanopt(options, "mass"))
	    bits |= MassBit;
	firstmass = FALSE;
	if (scanopt(options, "phi"))
	    bits |= PotentialBit;
	if (scanopt(options, "acc"))
	    bits |= AccelerationBit;
    }
    if (bits != 0 && outstr != NULL) {
	put_snap(outstr, &bodytab, &nbody, &tnow, &bits);
	if (bits & PhaseSpaceBit)
	    printf("\n\tparticle data written\n");
    }
    if ((bits & PhaseSpaceBit) != 0 && quadstr != NULL)
	put_quadfield(quadstr);
    if (*savefile)        			/* state file specified?    */
	savestate(savefile);			/*   save system data       */
}

/*
 * DIAGNOSTICS: compute set of dynamical diagnostics.
 */

local diagnostics()
{
    int i;
    Body *p;
    real velsq;
    vector tmpv;
    matrix tmpt;

    mtot = 0.0;					/* zero total mass          */
    etot[1] = etot[2] = 0.0;			/* zero total KE and PE     */
    CLRM(keten);				/* zero KE tensor           */
    CLRM(peten);				/* zero PE tensor           */
    CLRM(amten);				/* zero AM tensor           */
    CLRV(cmphase[0]);				/* zero c. of m. position   */
    CLRV(cmphase[1]);				/* zero c. of m. velocity   */
    for (p = bodytab; p < bodytab+nbody; p++) {	/* loop over all bodies     */
	mtot += Mass(p);                        /*   sum body masses        */
	DOTVP(velsq, Vel(p), Vel(p));		/*   square vel vector      */
	etot[1] += 0.5 * Mass(p) * velsq;	/*   sum current KE         */
	etot[2] += 0.5 * Mass(p) * Phi(p);	/*   and current PE         */
	MULVS(tmpv, Vel(p), 0.5 * Mass(p));	/*   sum 0.5 m v_i v_j      */
	OUTVP(tmpt, tmpv, Vel(p));
	ADDM(keten, keten, tmpt);
	MULVS(tmpv, Pos(p), Mass(p));		/*   sum m r_i a_j          */
	OUTVP(tmpt, tmpv, Acc(p));
	ADDM(peten, peten, tmpt);
	OUTVP(tmpt, tmpv, Vel(p));		/*   sum m r_i v_j          */
	ADDM(amten, amten, tmpt);
	MULVS(tmpv, Pos(p), Mass(p));		/*   sum cm position        */
	ADDV(cmphase[0], cmphase[0], tmpv);
	MULVS(tmpv, Vel(p), Mass(p));		/*   sum cm momentum        */
	ADDV(cmphase[1], cmphase[1], tmpv);
    }
    etot[0] = etot[1] + etot[2];                /* sum KE and PE            */
    TRANM(tmpt, amten);				/* antisymmetrize AM tensor */
    SUBM(amten, amten, tmpt);
    DIVVS(cmphase[0], cmphase[0], mtot);        /* normalize cm coords      */
    DIVVS(cmphase[1], cmphase[1], mtot);
}

/*
 * MY_PUT_SNAP_DIAGNOSTICS: output various N-body diagnostics.
 */

local my_put_snap_diagnostics(stream outstr, int *ofptr)
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
 * PUT_QUADFIELD: output tabulation of quadrupole field.
 */

local put_quadfield(stream quadstr)
{
    put_set(quadstr, "QuadField");
    put_data(quadstr, TimeTag, RealType, &tnow, 0);
    put_data(quadstr, "nqtab", IntType, &qfld.nqtab, 0);
    put_data(quadstr, "radtab", RealType, qfld.radtab, qfld.nqtab, 0);
    put_data(quadstr, "Q00tab", RealType, qfld.Q00tab, qfld.nqtab, 0);
    put_data(quadstr, "P00tab", RealType, qfld.P00tab, qfld.nqtab, 0);
    put_data(quadstr, "Q10tab", RealType, qfld.Q10tab, qfld.nqtab, 0);
    put_data(quadstr, "P10tab", RealType, qfld.P10tab, qfld.nqtab, 0);
    put_data(quadstr, "Q11tab", RealType, qfld.Q11tab, NDIM, qfld.nqtab, 0);
    put_data(quadstr, "P11tab", RealType, qfld.P11tab, NDIM, qfld.nqtab, 0);
    put_data(quadstr, "Q22tab", RealType, qfld.Q22tab,
	     NDIM, NDIM, qfld.nqtab, 0);
    put_data(quadstr, "P22tab", RealType, qfld.P22tab,
	     NDIM, NDIM, qfld.nqtab, 0);
    put_data(quadstr, "eps1", RealType, &eps1, 0);
    put_data(quadstr, "eps2", RealType, &eps2, 0);
    put_tes(quadstr, "QuadField");
    fflush(quadstr);
}

/*
 * SAVESTATE: write current state to disk file.
 */

savestate(string file)
{
    stream str;

    str = stropen(file, "a");			/* open state output file   */
    fseek(str, 0L, 0);				/* rewind stream to origin  */
    put_string(str, "program", getargv0());
    put_string(str, "version", getparam("VERSION"));
    put_string(str, "headline", headline);	/* save control parameters  */
    put_data(str, "freq", RealType, &freq, 0);
    put_data(str, "mode", IntType, &mode, 0);
    put_data(str, "eps1", RealType, &eps1, 0);
    put_data(str, "eps2", RealType, &eps2, 0);
    put_string(str, "options", options);
    put_data(str, "tstop", RealType, &tstop, 0);
    put_data(str, "minor_freqout", RealType, &minor_freqout, 0);
    put_data(str, "freqout", RealType, &freqout, 0);
    put_data(str, "tnow", RealType, &tnow, 0);	/* save state variables     */
    put_data(str, "minor_tout", RealType, &minor_tout, 0);
    put_data(str, "tout", RealType, &tout, 0);
    put_data(str, "nbody", IntType, &nbody, 0);
    put_data(str, "bodytab", AnyType, bodytab, nbody, sizeof(body), 0);
    put_data(str, "qfld", AnyType, &qfld, sizeof(qfld), 0);
    strclose(str);
}

/*
 * RESTORESTATE: restore state from disk file.
 */

restorestate(string file)
{
    stream str;
    string program, version;

    str = stropen(file, "r");			/* open state input file    */
    program = get_string(str, "program");
    version = get_string(str, "version");
    if (! streq(version, getparam("VERSION")))
	printf("warning: state file may be outdated\n\n");
    headline = get_string(str, "headline");	/* read control parameters  */
    get_data(str, "freq", RealType, &freq, 0);
    get_data(str, "mode", IntType, &mode, 0);
    get_data(str, "eps1", RealType, &eps1, 0);
    get_data(str, "eps2", RealType, &eps2, 0);
    options = get_string(str, "options");
    get_data(str, "tstop", RealType, &tstop, 0);
    get_data(str, "minor_freqout", RealType, &minor_freqout, 0);
    get_data(str, "freqout", RealType, &freqout, 0);
    get_data(str, "tnow", RealType, &tnow, 0);	/* read system state        */
    get_data(str, "minor_tout", RealType, &minor_tout, 0);
    get_data(str, "tout", RealType, &tout, 0);
    get_data(str, "nbody", IntType, &nbody, 0);
    bodytab = (Body *) allocate(nbody * sizeof(body));
    get_data(str, "bodytab", AnyType, bodytab, nbody, sizeof(body), 0);
    get_data(str, "qfld", AnyType, &qfld, sizeof(qfld), 0);
    strclose(str);
}
