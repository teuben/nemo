/*
 * MKCUBE.C: set up a uniform cube
 *	
 *	25-aug-93  V1.0	Created (cloned off mkhomsph)	PJT
 *	23-mar-96  V1.0a proto cleanup, free memory 	PJT
 *       9-sep-01      b    gsl/xrandom
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/put_snap.c>

string defv[] = {	/* DEFAULT INPUT PARAMETERS */
    "out=???\n      Output file name",
    "nbody=256\n    Number of particles",
    "size=1.0\n     Size of cube",
    "sigma=0.0\n    Isotropic velocity dispersion",
    "seed=0\n       Random number seed",
    "zerocm=t\n     Center c.o.m. ?",
    "headline=\n    Text headline for output",
    "VERSION=1.0b\n 9-sep-01 PJT",
    NULL,
};

string usage = "create a uniform cube of equal massive stars";

local real rmin, rmax;
local real sigma;
local bool zerocm;

local Body *btab;
local int nbody;

extern double xrandom(double,double), grandom(double,double);


void nemo_main()
{
    int seed;

    rmin = -0.5 * getdparam("size");
    rmax = -rmin;
    nbody = getiparam("nbody");
    sigma = getdparam("sigma");
    seed = init_xrandom(getparam("seed"));
    zerocm = getbparam("zerocm");
    mkcube();
    writegalaxy(getparam("out"), getparam("headline"));
    free(btab);
}

/*
 * WRITEGALAXY: write galaxy model to output.
 */

writegalaxy(name, headline)
string name;
string headline;
{
    stream outstr;
    real tsnap = 0.0;
    int bits = MassBit | PhaseSpaceBit;

    if (! streq(headline, ""))
	set_headline(headline);
    outstr = stropen(name, "w");
    put_history(outstr);
    put_snap(outstr, &btab, &nbody, &tsnap, &bits);
    strclose(outstr);
}

/*
 * MKCUBE: homogeneous cube
 */

mkcube()
{
    Body *bp;
    real rmin3, rmax3, r_i, theta_i, phi_i, mass_i;
    int i, ndim=NDIM;

    btab = (Body *) allocate(nbody * sizeof(Body));
    mass_i = 1.0/nbody;
    for (bp=btab, i = 0; i < nbody; bp++, i++) {
	Mass(bp) = mass_i;
	Phase(bp)[0][0] = xrandom(rmin,rmax);
	Phase(bp)[0][1] = xrandom(rmin,rmax);
	Phase(bp)[0][2] = xrandom(rmin,rmax);
	Phase(bp)[1][0] = (sigma > 0) ? grandom(0.0,sigma) : 0.0;
	Phase(bp)[1][1] = (sigma > 0) ? grandom(0.0,sigma) : 0.0;
	Phase(bp)[1][2] = (sigma > 0) ? grandom(0.0,sigma) : 0.0;
    }
    if (zerocm)
        centersnap(btab,nbody);
}

centersnap(btab, nb)
Body *btab;
int nb;
{
    real mtot;
    vector cmphase[2], tmp;
    Body *bp;

    mtot = 0.0;
    CLRV(cmphase[0]);
    CLRV(cmphase[1]);
    for (bp = btab; bp < btab + nb; bp++) {
	mtot = mtot + Mass(bp);
	MULVS(tmp, Phase(bp)[0], Mass(bp));
	ADDV(cmphase[0], cmphase[0], tmp);
	MULVS(tmp, Phase(bp)[1], Mass(bp));
	ADDV(cmphase[1], cmphase[1], tmp);
    }
    MULVS(cmphase[0], cmphase[0], 1.0/mtot);
    MULVS(cmphase[1], cmphase[1], 1.0/mtot);
    for (bp = btab; bp < btab + nb; bp++) {
	SUBV(Phase(bp)[0], Phase(bp)[0], cmphase[0]);
	SUBV(Phase(bp)[1], Phase(bp)[1], cmphase[1]);
    }
}

