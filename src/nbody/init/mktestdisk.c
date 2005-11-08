/*
 * MKTESTDISK.C: set up a uniform-density test disk in a 
 *	spherical N-body galaxy.
 *	22-jan-89  V1.2	Josh
 *	15-nov-90  V1.3 NEMO 2.x	PJT
 *      28-mar-97  V1.4 SINGLEPREC/proto's      PJT
 *                      but doesn't work in SINGLEPREC yet
 *       8-sep-01       a   init_xrandom
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>
#include <spline.h>

string defv[] = {	/* DEFAULT INPUT PARAMETERS */
    "in=???\n		  spheroid snapshot file ",
    "out=???\n		  output file name ",
    "rmin=0.1\n		  inner disk radius ",
    "rmax=0.4\n		  outer cutoff radius ",
    "eps=0.0\n		  softening (use w/ central point) ",
    "ndisk=2048\n	  number of disk particles ",
    "ncenter=0\n	  particles used to center disk ",
    "seed=0\n     	  usual random number seed ",
    "headline=\n	  text headline for output ",
    "VERSION=1.4b\n	  8-aug-05 PJT",
    NULL,
};

string usage="uniform density test-disk inside a spherical N-body";

string cvsid="$Id$";

local real rmin, rmax, eps;

local int ndisk, ncenter, nspheroid, ngalaxy;

local Body *spheroid, *galaxy;

void nemo_main()
{
    rmin = getdparam("rmin");
    rmax = getdparam("rmax");
    eps = getdparam("eps");
    ndisk = getiparam("ndisk");
    ncenter = getiparam("ncenter");
    init_xrandom(getparam("seed"));
    readspheroid(getparam("in"));
    if (ncenter > 0)
	centersnap(spheroid, ncenter, nspheroid);
    setsphprof();
    testdisk();
    centersnap(galaxy, ngalaxy, ngalaxy);
    writegalaxy(getparam("out"), getparam("headline"));
}

/*
 * READSPHEROID: read spheroid model from input.
 */

readspheroid(string name)
{
    stream instr;
    real tsnap;
    int bits;

    instr = stropen(name, "r");
    get_history(instr);
    get_snap(instr, &spheroid, &nspheroid, &tsnap, &bits);
    strclose(instr);
    if ((bits & MassBit) == 0 || (bits & PhaseSpaceBit) == 0)
	error("readspheroid: required data missing: bits = %#o", bits);
}

/*
 * WRITEGALAXY: write galaxy model to output.
 */

writegalaxy(string name, string headline)
{
    stream outstr;
    real tsnap = 0.0;
    int bits = MassBit | PhaseSpaceBit;

    if (! streq(headline, ""))
	set_headline(headline);
    outstr = stropen(name, "w");
    put_history(outstr);
    put_snap(outstr, &galaxy, &ngalaxy, &tsnap, &bits);
    strclose(outstr);
}

/*
 * SETSPHPROF: process spheroid to generate mass profile table [rsph, msph].
 */

#define NTAB  (256 + 1)				/* allow one for center     */

real rsph[NTAB];				/* radii of selected bodies */
real msph[4*NTAB];				/* mass in rsph, inclusive */
						/* spline coefs follow mass */

setsphprof()
{
    Body **rsort;
    int i, j, skip, i1, rankrad();

    rsort = (Body **) allocate(nspheroid * sizeof(Body *));
    for (i = 0; i < nspheroid; i++)
	rsort[i] = spheroid + i;
    qsort(rsort, nspheroid, sizeof(Body *), rankrad);
    j = 0;
    msph[j] = rsph[j] = 0.0;
    skip = (int) ceil(MAX((real) nspheroid / (NTAB-1.0), 1.0));
    for (i = 0; i < nspheroid; i++)
	if (absv(Pos(rsort[i])) <= rsph[j])
	    msph[j] = msph[j] + Mass(rsort[i]);
	else {
	    if (++j >= NTAB)
		error("setsphprof: table ovf nspheroid= %d  skip= %d NTAB=%d",
		      nspheroid, skip,NTAB);
	    msph[j] = msph[j-1] + Mass(rsort[i]);
	    i1 = MIN(i + skip, nspheroid) - 1;
	    rsph[j] = absv(Pos(rsort[i1]));
        }
    while (++j < NTAB) {
	msph[j] = msph[j-1];
	rsph[j] = rsph[j-1] + 1.0;
    }
    spline(&msph[NTAB], &rsph[0], &msph[0], NTAB);
}

int rankrad(Body **x, Body **y)
{
    real rxsq = dotvp(Pos(*x), Pos(*x));
    real rysq = dotvp(Pos(*y), Pos(*y));

    return (rxsq < rysq ? -1 : rxsq > rysq ? 1 : 0);
}

/*
 * TESTDISK: use tabulated spheroid mass profile to make a uniform
 * density test disk.  Note that the softening correction is only
 * exact for a central point mass.
 */

testdisk()
{
    Body *gp, *sp;
    real rmin2, rmax2, r_i, theta_i, msph_i, vcir_i;
    real cost, sint;
    int i;

    ngalaxy = nspheroid + ndisk;
    galaxy = (Body *) allocate(ngalaxy * sizeof(Body));
    rmin2 = rmin * rmin;
    rmax2 = rmax * rmax;
    for (i = 0, gp=galaxy; i < ndisk; i++, gp++) {
	Mass(gp) = 0.0;					/* zero mass        */
	r_i = sqrt(rmin2 + i * (rmax2 - rmin2) / (ndisk - 1.0));
	theta_i = xrandom(0.0, TWO_PI);
        cost = cos(theta_i);
        sint = sin(theta_i);
	Phase(gp)[0][0] = r_i * sint;		/* set positions    */
	Phase(gp)[0][1] = r_i * cost;
	Phase(gp)[0][2] = 0.0;
	if (r_i < rsph[NTAB-1])
	    msph_i = seval(r_i, &rsph[0], &msph[0], &msph[NTAB], NTAB);
	else
	    msph_i = msph[NTAB-1];
	vcir_i = sqrt(msph_i * r_i*r_i /
		        pow(r_i*r_i + eps*eps, 1.5));
	Phase(gp)[1][0] =   vcir_i * cost;	/* set velocities   */
	Phase(gp)[1][1] = - vcir_i * sint;
	Phase(gp)[1][2] = 0.0;
    }
    sp = spheroid;
    for (i = 0; i < nspheroid; i++)
	*gp++ = *sp++;
}

/*
 * CENTERSNAP: transform coordinates to the center defined
 * by a subset of the particles.
 */

centersnap(btab, nzero, nbody)
Body *btab;			/* array of bodies to translate		    */
int nzero;			/* first nzero bodies used to find center   */
int nbody;			/* number of bodies in table		    */
{
    real mtot;
    vector cmphase[2], tmp;
    Body *bp;

    if (nzero > nbody)
	error("centersnap: nzero = %d > nbody = %d", nzero, nbody);
    mtot = 0.0;
    CLRV(cmphase[0]);
    CLRV(cmphase[1]);
    for (bp = btab; bp < btab + nzero; bp++) {
	mtot = mtot + Mass(bp);
	MULVS(tmp, Phase(bp)[0], Mass(bp));
	ADDV(cmphase[0], cmphase[0], tmp);
	MULVS(tmp, Phase(bp)[1], Mass(bp));
	ADDV(cmphase[1], cmphase[1], tmp);
    }
    MULVS(cmphase[0], cmphase[0], 1.0/mtot);
    MULVS(cmphase[1], cmphase[1], 1.0/mtot);
    for (bp = btab; bp < btab + nbody; bp++) {
	SUBV(Phase(bp)[0], Phase(bp)[0], cmphase[0]);
	SUBV(Phase(bp)[1], Phase(bp)[1], cmphase[1]);
    }
}
