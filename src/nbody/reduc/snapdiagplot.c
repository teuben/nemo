/*
 * SNAPDIAGPLOT.C: read a snapshot output file and plot diagnostics.
 *	24-sep-87  V1.0a JEB  changed number of tickmarks
 *	22-jan-89  V1.1  JEB  Some formal version
 *	14-nov-91  V1.2  PJT  NEMO V2.  - try new trange[] 
 *      17-jun-92  V1.3  PJT  malloc -> allocate; process times=
 *                            and report fractional energy loss
 *	11-nov-92  V1.4  PJT  now using get_data_coerced in 4 places
 *	 5-aug-96  V1.4b PJT fixed local 'times' decl.
 *      16-feb-97  V1.4c PJT  SINGLEPREC
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <vectmath.h>
#include <snapshot/snapshot.h>
#include <yapp.h>
#include <axis.h>


string defv[] = {               /* DEFAULT INPUT PARAMETERS */
    "in=???\n                   input file name",
    "times=all\n		times to plot",
    "exactpot=false\n		if true, compute O(N^2) PE",
    "eps=0.05\n			using this softening",
    "formal=false\n		publication-style plot",
    "VERSION=1.4e\n		29-jul-09 PJT",
    NULL,
};

string usage="plot diagnostics of a snapshot file";

local string input;
local string times;
local bool exactpot;
local real eps;
local bool formal;

#define TIMEFUZZ 0.00001
    
/*
 * Arrays used to store diagnostic info.
 */

local int ndiagfr = 0;			/* number of observations */

#define MAXDIAGFR 100000		/* max. no. of observations */

local real Time[MAXDIAGFR];			/* time of observation */
local real deltar[MAXDIAGFR];		/* offset in position */
local real deltav[MAXDIAGFR];		/* offset in velocity */
local real KEnergy[MAXDIAGFR];		/* total kinetic energy */
local real PEnergy[MAXDIAGFR];		/* total potential energy */

local int npartfr = 0;			/* number of exact PE values */

#define MAXPARTFR 64			/* max. no. of exact PE values */

local int dfcount[MAXPARTFR];			/* observation number */
local real ExactPE[MAXPARTFR];		/* and computed PE value */

/*
 * Plotting ranges, to be expanded as needed.
 */

local real trange[2] =  {  0.00, 1.50 };	/* ### changed from 1.0  24/9/87 */
local real drrange[2] = {  0.00, 0.04 };
local real dvrange[2] = {  0.00, 0.04 };
local real dErange[2] = { -0.02, 0.02 };
local real Erange[2]  = { -1.00, 1.00 };

real ttrans(real), drtrans(real), Etrans(real), dEtrans(real);


#define NXTICK  5			/* ### changed from 7  24/9/87 */
#define NYTICK  3

nemo_main()
{
    stream instr;

    setparams();
    instr = stropen(input, "r");
    get_history(instr);
    readinput(instr);
    diagplot();
}

setparams()
{
    input = getparam("in");
    times = getparam("times");
    exactpot = getbparam("exactpot");
    eps = getdparam("eps");
    formal = getbparam("formal");
}


real poten(phsp, mass, nobj)
real phsp[][2][NDIM];
real mass[];
int nobj;
{
    int i, j;
    real epssq, petot, *posi, pei, dpos[NDIM], dpossq;

    epssq = eps * eps;
    petot = 0.0;
    for (i = 1; i < nobj; i++) {
	posi = phsp[i][0];
	pei = 0.0;
	for (j = 0; j < i; j++) {
	    SUBV(dpos, posi, phsp[j][0]);
	    DOTVP(dpossq, dpos, dpos);
	    pei = pei - mass[j] / sqrt(dpossq + epssq);
	}
	petot += mass[i] * pei;
    }
    return (petot);
}

readinput(stream instr)
{
    int nobj;
    real etot[3], cmphase[2][3];
    real *mptr, *pptr;

    mptr = NULL;
    while (get_tag_ok(instr, SnapShotTag)) {
	if (ndiagfr >= MAXDIAGFR)
	    error("readinput: diagnostics table overflow");
	get_set(instr, SnapShotTag);
	get_set(instr, ParametersTag);
	get_data_coerced(instr, TimeTag, RealType, &Time[ndiagfr], 0);
	if (!streq(times, "all") && !within(Time[ndiagfr], times, TIMEFUZZ))
            continue;
	get_data(instr, NobjTag, IntType, &nobj, 0);
	get_tes(instr, ParametersTag);
	if (exactpot && get_tag_ok(instr, ParticlesTag)) {
	    if (npartfr >= MAXPARTFR)
		error("particle table overflow");
	    get_set(instr, ParticlesTag);
	    if (mptr == NULL) {
		mptr = (real *) allocate(nobj * sizeof(real));
		pptr = (real *) allocate(nobj * 2 * NDIM * sizeof(real));
		get_data_coerced(instr, MassTag, RealType, mptr, nobj, 0);
	    }
	    get_data_coerced(instr, PhaseSpaceTag, RealType, pptr, nobj, 2, NDIM, 0);
	    get_tes(instr, ParticlesTag);
	    dfcount[npartfr] = ndiagfr;			/* save frame no. */
	    ExactPE[npartfr] = poten(pptr, mptr, nobj);
	    npartfr++;
	}
	get_set(instr, DiagnosticsTag);
	get_data_coerced(instr, EnergyTag, RealType, etot, 3, 0);
	get_data_coerced(instr, CMPhaseSpaceTag, RealType, cmphase, 2, NDIM, 0);
	get_tes(instr, DiagnosticsTag);
	KEnergy[ndiagfr] = etot[1];
	PEnergy[ndiagfr] = etot[2];
	deltar[ndiagfr] = absv(cmphase[0]);
	deltav[ndiagfr] = absv(cmphase[1]);
	get_tes(instr, SnapShotTag);
	ndiagfr++;
    }
    printf("%d diagnostic frames read\n", ndiagfr);
}


diagplot()
{
    real Etot, deltaE;
    int i, j;
    char msg[128];

    setlimits();
    if (! formal) {
	plinit("", 0.0, 20.0, 0.0, 20.0);
	xaxis( 2.0,  2.0, 16.0, trange, - NXTICK, ttrans, "time");
	xaxis( 2.0, 10.0, 16.0, trange, - NXTICK, ttrans, NULL);
	yaxis( 2.0,  2.0,  8.0, dErange, - NYTICK, dEtrans, "delta E / E");
	yaxis(18.0,  2.0,  8.0, dErange, - NYTICK, dEtrans, NULL);
	xaxis( 2.0, 18.0, 16.0, trange, - NXTICK, ttrans, NULL);
	xaxvar(-1, -1.0, 0.2, -1.0, -1.0);
	xaxis( 2.0, 10.0, 16.0, trange, - NXTICK, ttrans, NULL);
	yaxis( 2.0, 10.0,  8.0, drrange, - NYTICK, drtrans, "delta r");
	yaxis(18.0, 10.0,  8.0, drrange, - NYTICK, drtrans, NULL);
	sprintf(msg, "File: %s", input);
	pltext(msg, 2.0, 18.4, 0.32, 0.0);
    } else {
	plinit("", -1.0, 19.0, 0.0, 20.0);
	formalaxis = TRUE;
	xaxisvar.labdn = 0.44;
	xaxisvar.szlab = 0.40;
	yaxisvar.numdn = 0.24;
	yaxisvar.labdn = 0.24;
	yaxisvar.szlab = 0.40;
	xaxis( 2.0,  2.0, 16.0, trange, - NXTICK, ttrans, "t");
	xaxis( 2.0, 10.0, 16.0, trange, - NXTICK, ttrans, NULL);
	yaxis( 2.0,  2.0,  8.0, dErange, - NYTICK, dEtrans, "dE/E");
	yaxis(18.0,  2.0,  8.0, dErange, - NYTICK, dEtrans, NULL);
    }
    Etot = KEnergy[0] + PEnergy[0];
    for (i = 0; i < ndiagfr; i++) {
	if (! formal) {
	    plpoint(ttrans(Time[i]), drtrans(deltar[i]));
	    if (i % 4 == 0) {
		plpoint(ttrans(Time[i]), Etrans(KEnergy[i]));
		plpoint(ttrans(Time[i]), Etrans(PEnergy[i]));
	    }
	}
	deltaE = (KEnergy[i] + PEnergy[i]) / Etot - 1.0;
	plpoint(ttrans(Time[i]), dEtrans(deltaE));
    }
    for (j = 0; j < npartfr; j++) {
	i = dfcount[j];
	deltaE = (KEnergy[i] + ExactPE[j]) / Etot - 1.0;
	plcircle(ttrans(Time[i]), dEtrans(deltaE), 0.1);
    }
    plstop();
}

setlimits()
{
    real Etot, deltaE, WorstTime, WorstDeltaE=0.0;
    int i;

    if (0.99*Time[ndiagfr-1] + 0.01*Time[ndiagfr-2] > trange[1] ) {   /* up */
        while (0.99*Time[ndiagfr-1] + 0.01*Time[ndiagfr-2] > trange[1])
            trange[1] = 2 * trange[1];  
    } else {                                                 /* kludge down */
        trange[0] = Time[0] - 0.05*(Time[ndiagfr-1]-Time[0]);
        trange[1] = Time[ndiagfr-1] + 0.05*(Time[ndiagfr-1]-Time[0]);
	warning("Autoscaling time. MinMax=%g %g",trange[0],trange[1]);
    }
    Etot = KEnergy[0] + PEnergy[0];
    for (i = 0; i < ndiagfr; i++) {
        while (deltar[i] > drrange[1])
	    drrange[1] = 2 * drrange[1];
        while (deltav[i] > dvrange[1])
	    dvrange[1] = 2 * dvrange[1];
	deltaE = (KEnergy[i] + PEnergy[i]) / Etot - 1.0;
        while (deltaE < dErange[0] || dErange[1] < deltaE) {
	    dErange[0] = 2 * dErange[0];
	    dErange[1] = 2 * dErange[1];
	}
	while (PEnergy[i] < Erange[0] || KEnergy[i] > Erange[1]) {
	    Erange[0] = 2 * Erange[0];
	    Erange[1] = 2 * Erange[1];
	}
        if (ABS(deltaE) > ABS(WorstDeltaE))  {
            WorstTime = Time[i];
            WorstDeltaE = deltaE;
        }
    }
    printf("Worst fractional energy loss dE/E = (E_t-E_0)/E_0 = %g at T = %g\n",
           WorstDeltaE,WorstTime);
}

real ttrans(real t)
{
    return (2.0 + 16.0 * (t - trange[0]) / (trange[1] - trange[0]));
}

real drtrans(real dr)
{
    return (10.0 + 8.0 * (dr - drrange[0]) / (drrange[1] - drrange[0]));
}

real dvtrans(real dv)
{
    return (10.0 + 8.0 * (dv - dvrange[0]) / (dvrange[1] - dvrange[0]));
}

real dEtrans(real dE)
{
    return (2.0 + 8.0 * (dE - dErange[0]) / (dErange[1] - dErange[0]));
}

real Etrans(real E)
{
    return (2.0 + 8.0 * (E - Erange[0]) / (Erange[1] - Erange[0]));
}
