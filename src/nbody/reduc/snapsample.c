/*
 * SNAPSAMPLE.C: read a snapshot, optionally sort by binding energy,
 * and write a uniformly sampled subset of the particles.
 *
 * V1.1 22-jan-89    JEB
 * V1.2 19-nov-93    PJT  NEMO V2.x
 */

#include <stdinc.h>
#include <vectmath.h>
#include <getparam.h>
#include <filestruct.h>
#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

string defv[] = {
    "in=???\n			input file name",
    "out=???\n			output file name",
    "nsamp=32\n			number of particles to sample",
    "offset=0\n			select sampling to use",
    "besort=false\n		if true, sort by binding energy",
    "VERSION=1.2\n		 22 Jan 1989",
    NULL,
};

string usage="sub-sample of a snapshot";


nemo_main()
{
    stream instr, outstr;
    string headline = NULL;
    Body *btab = NULL;
    int nsamp, offset, nbody, bits, rank_be();
    real tsnap;

    instr = stropen(getparam("in"), "r");
    get_history(instr);
    outstr = stropen(getparam("out"), "w");
    put_history(outstr);
    nsamp = getiparam("nsamp");
    offset = getiparam("offset");
    get_snap(instr, &btab, &nbody, &tsnap, &bits);
    if (getbparam("besort")) {
	if ((bits & PhaseSpaceBit) == 0 || (bits & PotentialBit) == 0)
	    error("%s: essential data missing\tbits = %#o\n",
		  getargv0(), bits);
	qsort(btab, nbody, sizeof(Body), rank_be);
    }
    snapsample(btab, nbody, nsamp, offset);
    put_snap(outstr, &btab, &nsamp, &tsnap, &bits);
}

int rank_be(a, b)
Body *a, *b;
{
    real be_a, be_b;

    be_a = Phi(a) + 0.5 * dotvp(Phase(a)[1], Phase(a)[1]);
    be_b = Phi(b) + 0.5 * dotvp(Phase(b)[1], Phase(b)[1]);
    return (be_a < be_b ? -1 : 1);
}

snapsample(btab, nbody, nsamp, joff)
Body btab[];
int nbody;
int nsamp;
int joff;
{
    int jskip, i, j;

    if (nbody % nsamp != 0)
	error("snapsample: nbody%nsamp not zero\n");
    jskip = nbody / nsamp;
    if (joff >= jskip)
	error("snapsample: offset too big\n");
    for (i = 0, j = 0; i < nsamp; i++, j += jskip) {
	btab[i] = btab[j + joff];
	Mass(&btab[i]) = jskip * Mass(&btab[i]);
    }
}
