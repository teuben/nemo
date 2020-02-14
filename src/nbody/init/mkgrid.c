/*
 * MKGRID.C: set up a cubic grid to apply for getting potential.
 *	
 *	20-dec-05  V1.0 created,		 					JJF
 *             entirely based on mkcube V1.0b by PJT
 *      13-feb-2020 V1.1     prototypes + typo fix      PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/put_snap.c>

string defv[] = {	/* DEFAULT INPUT PARAMETERS */
    "out=???\n      Output file name",
    "ngrid=2\n      Number of points on one axis of the grid",
    "size=1.0\n     Size of cube",
    "zerocm=f\n     Center c.o.m. ?",
    "headline=\n    Text headline for output",
    "VERSION=1.1\n  13-feb-2020 JJF",
    NULL,
};

string usage = "create a cube of equal massive stars equally spaced in space";

string cvsid="$Id$";


local real rmin, rmax;
local bool zerocm;

local Body *btab;
local int ngrid, nbody;

void writegalaxy(string name, string headline);
void mkgrid(void);
void centersnap(body *btab, int nb);


void nemo_main(void)
{
    int seed;

    rmin = -0.5 * getdparam("size");
    rmax = -rmin;
    ngrid = getiparam("ngrid");
    nbody = ngrid * ngrid * ngrid ;
    zerocm = getbparam("zerocm");
    mkgrid();
    writegalaxy(getparam("out"), getparam("headline"));
    free(btab);
}

/*
 * WRITEGALAXY: write galaxy model to output.
 */

void writegalaxy(string name, string headline)
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
 * MKGRID: Cubic grid
 */

void mkgrid()
{
    Body *bp;
    real rmin3, rmax3, r_i, theta_i, phi_i, mass_i, delta;
    int i, j, k, l, ndim=NDIM;

    btab = (Body *) allocate(nbody * sizeof(Body));
    mass_i = 2.0/nbody;
    j = k = l = 0; 
    delta = (rmax-rmin)/(ngrid-1) ;
    for (bp=btab, i = 0; i < nbody; bp++, i++) {
                Mass(bp) = mass_i;
		Phase(bp)[0][0] =  rmin + j* delta;
		Phase(bp)[0][1] =  rmin + k* delta;
		Phase(bp)[0][2] =  rmin + l* delta;
		Phase(bp)[1][0] =  0.0;
		Phase(bp)[1][1] =  0.0;
		Phase(bp)[1][2] =  0.0;
		j++;
		if (j==ngrid) {j=0;k++;}
		if (k==ngrid) {k=0;l++;}
    }
    if (zerocm)
        centersnap(btab,nbody);
}

void centersnap(Body *btab, int nb)
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

