/*
 * CODE_IO.C: I/O routines for a simple direct N-body code.
 * Public routines: inputdata(), initoutput(), stopoutput(), output(),
 *
 *   16-feb-04   cleanup while cloning from hackcode1
 */

#include "code.h"
#include <filestruct.h>
#include <history.h>

/*	Snapshot I/O routines - special local one for diagnostics */
#include <snapshot/snapshot.h>
#include <snapshot/get_snap.c>
#define put_snap_diagnostics  my_put_snap_diagnostics
local void put_snap_diagnostics(stream, int *);
#include <snapshot/put_snap.c>


local void diagnostics(void);

extern bool scanopt(string, string);
extern double cputime(void);

/*
 * INPUTDATA: read initial conditions from input file.
 */

void inputdata(string file)
{
    stream instr;
    int bits;
    bodyptr p;

    instr = stropen(file, "r");			/* open input stream        */
    get_history(instr);				/* read file history data   */
    if (ask_headline() != NULL)			/* if headline was present  */
	headline = ask_headline();		/*   set headline for run   */
    bodytab = NULL;				/* request new input data   */
    get_snap(instr, &bodytab, &nbody, &tnow, &bits);
    						/* invoke generic input     */
    strclose(instr);				/* close input stream       */
    if ((bits & MassBit) == 0 || (bits & PhaseSpaceBit) == 0)
	error("inputdata: essential data missing\tbits = 0x%x", bits);
    if ((bits & TimeBit) == 0 || scanopt(options, "reset_time"))
						/* time missing or reset?   */
	tnow = 0.0;				/*   then supply default    */
}

/*
 * INITOUTPUT: initialize output routines.
 */

local stream outstr;                  /* output stream pointer */

void initoutput(void)
{
    printf("\n%s\n\n", headline);               /* print headline, params   */
    printf("%12s%12s%12s\n",
           "nbody", "freq", "eps");
    printf("%12d%12.2f%12.4f\n\n",
           nbody, freq, eps);
    if (*options)
        printf("\toptions: %s\n", options);
    if (*outfile) { 		                /* output file specified?   */
        outstr = stropen(outfile, "w");         /*   setup output stream    */
	put_history(outstr);			/*   write file history     */
    } else
        outstr = NULL;				/*   prevent binary output  */
}

/*
 * STOPOUTPUT: finish up after a run.
 */

void stopoutput(void)
{
  if (outstr) strclose(outstr);
}

/*
 * Counters and accumulators for output routines.
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

local bool firstmass = TRUE;	/* if true, output mass data */

void output(void)
{
    int k, bits;

    diagnostics();				/* compute std diagnostics  */
    printf("\n  %10s%10s%10s%10s\n",
           "tnow", "T+U", "T/U", "cputime");
    printf("  %10.3f%10.4f%10.4f%10.2f\n\n",
           tnow, etot[0], etot[1]/etot[2], cputime());
    printf("\t    %10s", "cm pos");
    for (k = 0; k < NDIM; k++)
        printf("%10.4f", cmphase[0][k]);
    printf("\n\t    %10s", "cm vel");
    for (k = 0; k < NDIM; k++)
        printf("%10.4f", cmphase[1][k]);
    printf("\n");
    bits = 0;					/* collect output bit flags */
    if (minor_freqout > 0.0 && (minor_tout - 0.01/freq) <= tnow) {
	minor_tout += 1.0 / minor_freqout;
	bits |= TimeBit;
    }
    if (freqout > 0.0 && (tout - 0.01/freq) <= tnow) {
	tout += 1.0 / freqout;
	bits |= TimeBit | PhaseSpaceBit;
	if (scanopt(options, "mass") || firstmass) {
	    bits |= MassBit;
	    firstmass = FALSE;
	}
	if (scanopt(options, "phi"))
	    bits |= PotentialBit;
	if (scanopt(options, "acc"))
	    bits |= AccelerationBit;
    }
    if (bits != 0 && outstr != NULL) {		/* output ready and able?   */
	put_snap(outstr, &bodytab, &nbody, &tnow, &bits);
	if (bits & PhaseSpaceBit)
	    printf("\n\tparticle data written\n");
    }
}

/*
 * DIAGNOSTICS: compute set of dynamical diagnostics.
 */

local void diagnostics(void) 
{
    bodyptr p;
    real velsq;
    vector tmpv;
    matrix tmpt;

    mtot = 0.0;					/* zero total mass          */
    etot[1] = etot[2] = 0.0;			/* zero total KE and PE     */
    CLRM(keten);				/* zero ke tensor           */
    CLRM(peten);				/* zero pe tensor           */
    CLRM(amten);				/* zero am tensor           */
    CLRV(cmphase[0]);				/* zero c. of m. position   */
    CLRV(cmphase[1]);				/* zero c. of m. velocity   */
    for (p = bodytab; p < bodytab+nbody; p++) {	/* loop over all particles  */
	mtot += Mass(p);                        /*   sum particle masses    */
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
    TRANM(tmpt, amten);				/* anti-sym. AM tensor      */
    SUBM(amten, amten, tmpt);
    DIVVS(cmphase[0], cmphase[0], mtot);        /* normalize cm coords      */
    DIVVS(cmphase[1], cmphase[1], mtot);
}

/*
 * MY_PUT_SNAP_DIAGNOSTICS: output various N-body diagnostics.
 */

local void my_put_snap_diagnostics(stream outstr, int *ofptr)
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
